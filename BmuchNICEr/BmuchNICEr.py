#! usr/bin/env python

"""
This script initiates a BNICE run and waits for completion. It then reads BNICE output (flatfiles and molfiles) and
loads the information into a Mongo Database. Duplicate compounds and reactions are excluded so this script may be used
to merge reaction networks

usage: %BmuchNICEr.py [options] run_name
Use BmuchNICE.py -h to see options
"""

#Created by James Jeffryes

import sys
import shutil
import time
import os
import openbabel as ob
import hashlib
import platform
import re
import json
from subprocess import call
from pymongo import MongoClient
from optparse import OptionParser


def parse_options():
    #This block handles user flags and arguments. For more information see the optparse API documentation
    usage = "usage: %prog [options] run_name"
    parser = OptionParser(usage)
    parser.add_option("-f", "--start-compound-file", dest="start_compound_path", default="Input/StartCompounds.dat",
                      help="The path to the desired start compound list .dat file")
    parser.add_option("-b", "--build-start-compound-file", dest="build_start_compounds", default="null",
                      help="Specify a folder containing molfiles (Must be in the Molfiles directory) to use as start"
                           " compounds w/ standard cofactors")
    parser.add_option("-k", "--kegg-compound", dest="kegg_compound", default="null",
                      help="Specify a single KEGG compound to use as a start compound w/ standard cofactors")
    parser.add_option("-o", "--operator-list", dest="operator_path", default="Input/AllOperators(Forward)_NoNonG.dat",
                      help="The path to the desired operator list .dat file FROM CORE DIRECTROY")
    parser.add_option("-s", "--single-operator", dest="single_operator", default="null",
                      help="specify a path FROM CORE FOLDER to a single operator to use")
    parser.add_option("-d", "--database", dest="DB_name", default="null",
                      help="The name of the mongo database into which the BNICE data will be loaded. If unspecified,"
                           " loading does not occur")
    parser.add_option("-g", "--generations", dest="generation_cap", type="int", default="1",
                      help="The number of generations that NetGen will attempt to create. Set to 1 by default.")
    parser.add_option("-r", "--retro", dest="retro_flag", action="store_true", default=False,
                      help="Tag as a retro run. Reaction directions will be swapped when loaded into the database.")
    parser.add_option("-t", "--trash-collection", dest="trash_collection", type="int", default="1",
                      help="If set to 0, no files are removed after database loading. If set to 1, the input files the"
                           " script creates as well as sapphire log files. If set to 2, the script will also remove the"
                           " molfile directory. if set to 3, all BNICE output is removed. Set to 1 by default.")
    parser.add_option("-T", "--Thermo", dest="thermo_flag", action="store_true", default=False,
                      help="Run thermo on the netgen output")
    parser.add_option("-S", "--SimIndex", dest="sim_target", default="null",
                      help="Specify a path target compound for a SimIndex run FROM THE MAIN DIRECTORY(i.e. one level above "
                           "BmuchNICEr) If unspecified a normal NetGen run takes place.")
    parser.add_option("--SimIndex-tol", dest="sim_tol", type="int", default="1",
                      help="specify the tolerance for a SimIndex run. 1 by default.")
    (options, args) = parser.parse_args()
    #Ensure that the right number of arguments are passed to the script.
    if len(args) == 0:
        sys.exit("Please supply a name for the netgen run as an argument")
    if len(args) > 1:
        sys.exit("Too many arguments submitted. Usage: [options] run_name")
    run_name = args[0]
    return run_name, options


def generate_parameters_file(options, run_name):
    #if the KEGG compound option is active we generate a new compound file with that compound,
    # otherwise just use the specified start compound file.
    if options.kegg_compound != "null":
        start_compound_path = kegg_start_compound(options.kegg_compound, run_name)
    elif options.build_start_compounds != "null":
        start_compound_path = build_start_compounds_list(options.build_start_compounds, run_name)
    else:
        start_compound_path = options.start_compound_path

    #if the single operator option is active we create a operators list file, otherwise we simply use the specified
    # operate list file.
    if options.single_operator != "null":
        operator_path = generate_operators_file(options.single_operator, run_name)
    else:
        operator_path = options.operator_path

    #open the template and copy the data into memory as a list where each line is a string
    try:
        infile = open("./NetGen_template.txt")
    except IOError:
        sys.exit("NetGen template not found")
    NetGen_data = infile.readlines()
    infile.close()

    #replace the strings of parameters that get tweaked
    NetGen_data[1] = "start compounds filename|%s\n" % start_compound_path
    NetGen_data[2] = "operator list filename|%s\n" % operator_path
    NetGen_data[3] = "output folder|NetGen/%s\n" % run_name
    NetGen_data[7] = "rank limit|%s\n" % options.generation_cap

    #write the parameter info in a new file in the parameters folder for NetGen to read
    parameters_file = "../Core/Parameters/NetGen-%s.txt" % run_name
    outfile = open(parameters_file, 'w')
    outfile.writelines(NetGen_data)
    outfile.close()

    return parameters_file


def build_start_compounds_list(start_compound_folder, run_name):
    #open the template and copy the data into memory as a list where each line is a string
    try:
        infile = open("./StartCompounds_template.dat")
    except IOError:
        sys.exit("StartCompounds template not found")
    StartCompounds_data = infile.readlines()
    infile.close()

    if not os.path.exists("../Molfiles/%s" % start_compound_folder):
        sys.exit("Folder of start compounds not found. Ensure compounds folder is present in Molfiles folder")

    #Walk the indicated directory and add the files to the list of starting compounds
    for x in os.listdir("../Molfiles/%s" % start_compound_folder):
        if '.mol' in x:
            StartCompounds_data.append("%s/%s\n" % (start_compound_folder, x))
    compound_count = len(StartCompounds_data) - 1

    #replace the strings to update the starting compound count and clip off the new line character at the end of the
    # file (which pisses of BNICE if left in.)
    StartCompounds_data[0] = "%s\n" % compound_count
    last_compound = str(StartCompounds_data[compound_count])
    StartCompounds_data[compound_count] = last_compound.rstrip("\n")

    #write the start compound info in a new file in the input folder for NetGen to read
    StartCompounds_file = "../Core/Input/StartCompounds-%s.dat" % run_name
    outfile = open(StartCompounds_file, 'w')
    outfile.writelines(StartCompounds_data)
    outfile.close()

    return StartCompounds_file


def kegg_start_compound(kegg_compound, run_name):
    #open the template and copy the data into memory as a list where each line is a string
    try:
        infile = open("./StartCompounds_template.dat")
    except IOError:
        sys.exit("StartCompounds template not found")
    StartCompounds_data = infile.readlines()
    infile.close()

    #replace the strings to add the compound and update the starting count
    StartCompounds_data[0] = "35\n"
    StartCompounds_data.append("Kegg/%s.mol" % kegg_compound)

    #write the start compound info in a new file in the input folder for NetGen to read
    StartCompounds_file = "../Core/Input/StartCompounds-%s.dat" % run_name
    outfile = open(StartCompounds_file, 'w')
    outfile.writelines(StartCompounds_data)
    outfile.close()

    return StartCompounds_file


def generate_operators_file(single_operator, run_name):

    Operator_data = ["1\n", single_operator]

    #write the operator info in a new file in the input folder for NetGen to read
    Operator_file = "../Core/Input/Operator-%s.dat" % run_name
    outfile = open(Operator_file, 'w')
    outfile.writelines(Operator_data)
    outfile.close()

    return Operator_file


def call_SimIndex(options, run_name):

    #Check if the target specified is a full path, KEEG compound or bullshit.
    if os.path.exists("../"+options.sim_target):
        sim_target = options.sim_target
    elif os.path.exists("../Molfiles/Kegg/%s.mol" % options.sim_target):
        sim_target = "./Molfiles/Kegg/%s.mol" % options.sim_target
    else:
        sys.exit("Could not locate the target compound specified for SimIndex")

    os.chdir('..')
    #Calls a run of BNICE through the shell.
    rc = call(['python SIModule.py %s NetGen-%s.txt %d' %
             (sim_target, run_name, options.sim_tol)], shell=True)

    #If the python script doesn't finish successfully, quit.
    if rc != 0:
        sys.exit("SimIndex error")

    os.chdir("BmuchNICEr")


def call_netgen(run_name):
    #Goes to the core folder to execute (required for paths in parameters files to work correctly)
    execution_dir = os.getcwd()
    os.chdir("../Core")

    #Creates a file to feed bnice its parameters as standard input one at a time
    standard_input = open("%s_input.txt" % run_name, "w+")
    standard_input.write("NetGen-%s\n1\n%s\n" % (run_name, run_name))
    standard_input.seek(0, 0)

    #Calls a run of BNICE through the shell. The return code is collected and ignored
    call(["."+str(os.sep)+"bnice.exe"], stdin=standard_input, shell=True)

    #Removes the standard input file and goes back to the script's home directory
    os.remove("%s_input.txt" % run_name)
    os.chdir(execution_dir)


def run_thermo(run_name):
    #open the template and copy the data into memory as a list where each line is a string
    try:
        infile = open("./Thermo_template.txt")
    except IOError:
        sys.exit("Thermo template not found")
    Thermo_data = infile.readlines()
    infile.close()

    #replace the strings of parameters that get tweaked
    Thermo_data[0] = "molfiles|../Output/NetGen/%smolfiles/\n" % run_name

    #write the parameter info in a new file in the parameters folder for BNICE to read
    thermo_file = "../Core/Parameters/Thermo-%s.txt" % run_name
    outfile = open(thermo_file, 'w')
    outfile.writelines(Thermo_data)
    outfile.close()

    call_thermo(run_name)

    #replace the strings of parameters that get tweaked
    Thermo_data[18] = "balance reactions|0\n"
    Thermo_data[38] = "automatically add electrons|0\n"
    Thermo_data[39] = "automatically add H|0\n"

    #overwrite the old file with changes.
    thermo_file = "../Core/Parameters/Thermo-%s.txt" % run_name
    outfile = open(thermo_file, 'w')
    outfile.writelines(Thermo_data)
    outfile.close()

    call_thermo(run_name)


def call_thermo(run_name):
    #Goes to the core folder to execute (required for paths in parameters files to work correctly)
    execution_dir = os.getcwd()
    os.chdir("../Core")

    #Creates a file to feed bnice its parameters as standard input one at a time
    standard_input = open("%s_input.txt" % run_name, "w+")
    standard_input.write("Thermo-%s\n2\n%s\n" % (run_name, run_name))
    standard_input.seek(0, 0)
    standard_output = open("../Output/ThermoOutput.log", "a")

    #Calls a run of BNICE through the shell. The return code is collected and ignored
    # noinspection PyUnusedLocal
    rc = call(["."+str(os.sep)+"bnice.exe"], stdin=standard_input, stdout=standard_output, shell=True)

    #Removes the standard input file and the old flatfile
    os.remove("%s_input.txt" % run_name)
    os.remove("./Systems/%s.txt" % run_name)

    #Renames the new flatfile and goes back to the script's home directory
    os.rename("./Systems/%sNew.txt" % run_name, "./Systems/%s.txt" % run_name)
    os.chdir(execution_dir)


class Network:
    """This class makes passing around information from function to function when loading the database cleaner"""
    def __init__(self, run_name, DB_name, retro_flag):
        self.run_name = run_name
        self.DB_name = DB_name
        self.retro = retro_flag
        self.smiles = []
        self.smile_hashes = []
        self.comp_dict_list = []
        self.rxn_dict_list = []
        self.compound_generation = []
        self.rxn_generation = []
        self.add_db_links = False

    def __str__(self):
        pass

    def __repr__(self):
        pass

    def create_smiles(self):

        #I am lazy so we just move to the molfile directory to execute the function
        execution_dir = os.getcwd()
        os.chdir("../Output/NetGen/%smolfiles" % self.run_name)

        #Count the number of molfiles to convert.
        mol_number = len([x for x in os.listdir('.') if '.mol' in x])
        print("Converting %d mol files" % mol_number)

        #Suppress openbabel/pyble error output.
        error_log_check = hasattr(ob, "obErrorLog")
        if error_log_check:
            ob.obErrorLog.ClearLog()
            ob.obErrorLog.SetOutputLevel(-1)

        #initalize conversion object and set it to output canonical smiles w/o stereochemistry
        obConv = ob.OBConversion()
        obConv.SetInAndOutFormats("mol", "can")
        obConv.AddOption('i')

        #Read each molfile in the directory, convert and append to a list
        #Watch out! If there are files in the molfile directory other than molfiles, this code will throw an index error
        molecule = ob.OBMol()
        for molfile in range(1, mol_number + 1):
            obConv.ReadFile(molecule, "%s.mol" % molfile)
            line = obConv.WriteString(molecule)

            # Edit the smile stringcode to correct pseudoatoms that BNICE uses but open babel can't understand.
            if "CoA" in line:
                line = line.replace("*", "CoA")
            if "R" in line:
                line = line.replace("*", "R")
            if "E" in line:
                line = line.replace("*", "E")
            if line == '\tH\n':
                line = "H+"

            #Discard the chemical formula (it's already in the flatfile)
            line = line.split("\t")[0]

            #collect the smile string code and generate a unique hash from this string code that will serve as a unique
            #compound ID
            smile_hash = hashlib.sha1(line).hexdigest()

            #build lists of smiles and hashes for later use
            self.smiles.append(line)
            self.smile_hashes.append(smile_hash)

        #Return back to where we started
        os.chdir(execution_dir)

    def parse_text_file(self):

        #read the file, select the compound information, and split it line by line into a list with 1 compound per entry
        flat_data = open("../Core/Systems/%s.txt" % self.run_name).read()
        flat_comps = flat_data.split('REACTIONS\n')[0].split('\n')[1:-1]

        #further split the data vertically to create a lists of lists. if we hit the operator section, stop
        split_comp_info = []
        for x in flat_comps:
            if "OPERATORS" == x:
                break
            split_comp_info.append(x.split(';'))

        #this code converts the internal lists into a dictionary with keys as the column headers
        comp_dict_list = [{key.strip(): value for key, value in zip(split_comp_info[0], row)} for row in split_comp_info[1:]]

        #read reactions section of the file and split it into lines. if no reactions present exit the method
        try:
            rxn_info = flat_data.split('REACTIONS\n')[1].splitlines()

            if not 'ENTRY;' in rxn_info[0]:
                print "No Reactions present"
                return comp_dict_list, []

            #Split each reaction entry vertically into fields. Stop if we hit OPERATORS or ferredoxin (for thermo).
            split_rxn_info = []
            for x in rxn_info:
                if ("OPERATORS" == x) or ("oxidoreductase" in x):
                    break
                split_rxn_info.append(x.split(';'))

            #this code turns the internal lists into dictionaries with the headers as keys
            rxn_dict_list = [{key.strip(): value for key, value in zip(split_rxn_info[0], row)} for row in split_rxn_info[1:]]

        except IndexError:
            print "No Reactions present"
            rxn_dict_list = []

        return comp_dict_list, rxn_dict_list

    def parse_compound_information(self, comp_dict_list):

        #this checks to make sure we don't have missing compounds from the flat files or molfiles
        if len(comp_dict_list) == len(self.smile_hashes) and len(comp_dict_list) == len(self.smiles):
                for i, x in enumerate(comp_dict_list):
                    #replace entry number with a unique hash id which explicitly diferentiates cofactors
                    if x['COFACTOR'] == '0':
                        self.smile_hashes[i] = 'C' + self.smile_hashes[i]
                        #rather than delete keys from the existing dictionary, just create a new dictionary with only
                        # the fields we care about
                        self.comp_dict_list.append(dict(_id=self.smile_hashes[i]))
                    if x['COFACTOR'] == '1':
                        self.smile_hashes[i] = 'X' + self.smile_hashes[i]
                        self.comp_dict_list.append(dict(_id=self.smile_hashes[i]))
                        self.comp_dict_list[i]['Name'] = x['NAME']
                    #replace BNICE's arbitrary stringcode with a cannonical smile string
                    self.comp_dict_list[i]['Stringcode'] = self.smiles[i]

                    #stores information about what generation the compound appears for use in a dedicated collection
                    self.compound_generation.append({'_id': self.smile_hashes[i], 'Generation': x['GENERATION']})

                    #transfer over information from old to new dictionary, including thermo information if it exists
                    self.comp_dict_list[i]['Formula'] = x['FORMULA']
                    if 'ENERGY' in x:
                        if x['ENERGY'] != '':
                            self.comp_dict_list[i]['Delta G'] = x['ENERGY']
                            self.comp_dict_list[i]['Error'] = x['ERROR']

        else:
            sys.exit("Missmatch in flatfile and molfile compound counts")

    def parse_reaction_information(self, rxn_dict_list):

        #the following code parses the equation field for the reaction data to extract the juicy morsels of
        #information from it then save that information as a tuple of the stoichometric coefficient and the compound hash
        for i, x in enumerate(rxn_dict_list):
            eqn_data = x['EQUATION']
            
            #take the either the left or right hand of the equation based on the retro flag
            if self.retro:
                raw_reactants = eqn_data.split('<==>')[1]
            else:
                raw_reactants = eqn_data.split('<==>')[0]
            reactant_hashes = []

            #iterate through reactants one compound at a time
            for y in raw_reactants.split('+'):
                #pull the numbers into a list and ignore the rest
                nums = re.findall('\d+', y)
                if len(nums) == 1:
                    #if only one number is present it must be the compound number
                    reactant_hashes.append(("1", self.smile_hashes[int(nums[0])-1]))
                elif len(nums) == 2:
                    #if 2 numbers are present the first is the stoic and second the compound number
                    reactant_hashes.append((nums[0], self.smile_hashes[int(nums[1])-1]))
                else:
                    raise ValueError("Failed to parse reaction reactants")

            #as with the compounds we create a new dictionary with only the information from BNICE we care about
            self.rxn_dict_list.append(dict(Reactants=reactant_hashes))

            #parseing the products proceeds exactly like reactants so reference the comments in the above section
            if self.retro:
                raw_products = eqn_data.split('<==>')[0]
            else:
                raw_products = eqn_data.split('<==>')[1]
            product_hashes = []
            for y in raw_products.split('+'):
                nums = re.findall('\d+', y)
                if len(nums) == 1:
                    product_hashes.append(("1", self.smile_hashes[int(nums[0])-1]))
                elif len(nums) == 2:
                    product_hashes.append((nums[0], self.smile_hashes[int(nums[1])-1]))
                else:
                    raise ValueError("Failed to parse reaction products")

            self.rxn_dict_list[i]['Products'] = product_hashes

            #like the hash of the stringcode forms a unique id for a compound, a hash of the reactants and products
            #uniquely identifies a reaction. The hashes are sorted first to ensure uniqueness
            reactant_hashes.sort()
            product_hashes.sort()
            sequenced_rxn = "+".join(str(z) for z in reactant_hashes)+'<>' + "+".join(str(z) for z in product_hashes)
            self.rxn_dict_list[i]['_id'] = 'R' + hashlib.sha1(sequenced_rxn).hexdigest()

            #store generation information in a separate list
            self.rxn_generation.append({'_id': self.rxn_dict_list[i]['_id'], 'Generation': x['GENERATION']})

            #if there is thermo data in the flatfile we copy that over (ignoring duplicate information)
            try:
                self.rxn_dict_list[i]['Energy'] = x['No_O2_ENERGY']
                self.rxn_dict_list[i]['Error'] = x['ERROR_O2']
                self.rxn_dict_list[i]['mM delta G'] = x['MMDELTAG']
            except KeyError:
                pass
            self.rxn_dict_list[i]['Operators'] = x['OPERATORS']

    def load_MongoDB(self):

        #this load the mongo client and quits if unable to establish a connection
        try:
            if 'node' in platform.node():
                client = MongoClient(host='master')
            else:
                client = MongoClient()
        except:
            sys.exit("Failed to load database client. Please verify that mongod is running")

        #this loads the database with the user's name of choice and two collections for compounds and reactions
        #note: if the user specifies a new database name, the database and collections will not be created until the
        # next step when data is actually inserted
        db = client[self.DB_name]

        start_comp_count = db.compounds.count()
        start_rxn_count = db.reactions.count()

        #This loads the dictionaries into the database into their respective collections.
        #Mongo will use the _id field as the primary index for products and reactants.
        #Mongo checks the the _id value to make sure that it does not already match a compound in the database.
        #If a duplicate is detected it overwrites the duplicate (not very efficent but robust)

        for x in self.comp_dict_list:
            db.compounds.update({'_id': x['_id']}, x, upsert=True)
        for x in self.rxn_dict_list:
            db.reactions.update({'_id': x['_id']}, x, upsert=True)

        #store the generation information in separate collections
        for x in self.compound_generation:
            db.compounds_generation.update({'_id': x['_id']}, {'$set': {'%s Generation' % self.run_name: x['Generation']}}, upsert=True)
        for x in self.rxn_generation:
            db.reactions_generation.update({'_id': x['_id']}, {'$set': {'%s Generation' % self.run_name: x['Generation']}}, upsert=True)

        #record the number of documents that we attempted to insert into the database 
        number_comp_attempted = len(self.comp_dict_list)
        number_rxns_attempted = len(self.rxn_dict_list)

        print "%s compounds and %s reactions were submitted to the database" \
              % (number_comp_attempted, number_rxns_attempted)

        #A little math to determine how many compounds were skipped over
        #WARNING: This will give nonsensical results if multiple scripts are writing to the database simultaneously
        end_comp_count = db.compounds.count()
        end_rxn_count = db.reactions.count()
        duplicate_comps = number_comp_attempted - (end_comp_count - start_comp_count)
        duplicate_rxns = number_rxns_attempted - (end_rxn_count - start_rxn_count)

        print "%s compounds were determined to be duplicates and not inserted. The database now holds %s compounds." \
              % (duplicate_comps, end_comp_count)

        print "%s reactions were determined to be duplicates and not inserted. The database now holds %s reactions." \
              % (duplicate_rxns, end_rxn_count)


def trash_collection(trash_flag, run_name):
    if trash_flag == 0:
        print "No files will be removed"
        return

    if trash_flag >= 1:
        print "Removing script generated input, parameter, and log files"
        
        #Attempt to remove uninformative input, parameter and log files
        os.remove("../Core/Parameters/NetGen-%s.txt" % run_name)

        if os.path.exists("../Core/Input/StartCompounds-%s.dat" % run_name):
            os.remove("../Core/Input/StartCompounds-%s.dat" % run_name)

        if os.path.exists("../Core/Input/Operator-%s.dat" % run_name):
            os.remove("../Core/Input/Operator-%s.dat" % run_name)

        if os.path.exists("../Output/NetGen/%sOutput.log" % run_name):
            os.remove("../Output/NetGen/%sOutput.log" % run_name)

        if os.path.exists("../Output/ThermoOutput.log"):
            os.remove("../Output/ThermoOutput.log")

    if trash_flag >= 2:
        print "Removing molfile directory"
        shutil.rmtree("../Output/NetGen/%smolfiles" % run_name)

    if trash_flag >= 3:
        #This should remove any remining BNICE output
        print "Removing all other BNICE outputs"

        if os.path.exists("../Core/Systems/%s.json" % run_name):
            os.remove("../Core/Systems/%s.json" % run_name)

        if os.path.exists("../Core/Systems/%s.txt" % run_name):
            os.remove("../Core/Systems/%s.txt" % run_name)

        if os.path.exists("../Output/NetGen/%sStats%s.txt" % (run_name, run_name)):
            os.remove("../Output/NetGen/%sStats%s.txt" % (run_name, run_name))

        if os.path.exists("../Output/ThermoError.log"):
            os.remove("../Output/ThermoError.log")

        if os.path.exists("../Output/NetGen/%sError.log" % run_name):
            os.remove("../Output/NetGen/%sError.log" % run_name)

    if trash_flag == 11:
        print "IT GOES TO 11!?!"


######################################################################################
#                                   Main code                                        #
######################################################################################
if __name__ == '__main__':
    time_1 = time.time()
    run_name, options = parse_options()

    if os.path.exists('../Output/%smolfiles' % run_name):
        sys.exit('Molfiles from previous run of the same name are in the Output folder. Please delete molfiles or rename run.')

    parameters_file = generate_parameters_file(options, run_name)
    time_2 = time.time()
    print "Generating necessary files took %s s" % (time_2-time_1)

    #Run SimIndex if target specified
    if options.sim_target != "null":
        call_SimIndex(options, run_name)
        time_3 = time.time()
        print "SimIndex run took %s s" % (time_3-time_2)

    #Run Netgen
    else:
        call_netgen(run_name)
        time_3 = time.time()
        print "NetGen run took %s s" % (time_3-time_2)

    #Run Thermo
    if options.thermo_flag:
        time_tstart = time.time()
        run_thermo(run_name)
        time_tend = time.time()
        print "Thermo runs took %s s" % (time_tend-time_tstart)

    #Parse files and load database
    if options.DB_name != "null":
        time_3 = time.time()
        new_net = Network(run_name, options.DB_name, options.retro_flag)

        new_net.create_smiles()
        time_4 = time.time()
        print "Converting mol files to smiles took %s s" % (time_4-time_3)

        #if the output is a json, parsing is straightforward. if it's a txt file we call a special parseing function
        if os.path.exists("../Core/Systems/%s.json" % run_name):
            with open("../Core/Systems/%s.json" % run_name) as sys_json:
                flat_dict = json.load(sys_json)
            comp_dict_list = flat_dict['COMPOUNDS']
            rxn_dict_list = flat_dict['REACTIONS']
        else:
            #parse .txt systems files
            comp_dict_list, rxn_dict_list = new_net.parse_text_file()

        new_net.parse_compound_information(comp_dict_list)
        new_net.parse_reaction_information(rxn_dict_list)
        time_5 = time.time()
        print "Parsing reaction and compound information took %s s" % (time_5-time_4)

        new_net.load_MongoDB()
        time_6 = time.time()
        print "Loading reaction network into the database took %s s" % (time_6-time_5)

    time_7 = time.time()
    trash_collection(options.trash_collection, run_name)
    time_final = time.time()
    print "Trash collection took %s s" % (time_final-time_7)
    print "BmuchNICEr run complete. Script executed in %s s." % (time_final-time_1)