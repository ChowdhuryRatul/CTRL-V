# Author: Yee Chuen Teoh
'''
update 4/8/2024: version 7
    - (NOT TESTED YET) update to MPNN to be a list of singlepoint mutation.

update 3/4/2024: version 6
    - (TESTED) Error handling if config file does not have enough parameter.

update 2/29/2024: version 6
    - (TESTED) add option to function createPreviousQueue() to do continuation or report.
    - (TESTED) update printExampleConfigFile() with most updated config file
    - (TESTED) add function to get IRL mutation from csv file. csv format refer to 
                "/work/ratul1/chuen/viral_escape/pyRosetta_script/CoV_IRL_mutation.csv"

update 2/28/2024: version 6
    - (TESTED) add option to sort Ag-Ab energy rank in ascending or descending.
    - (TESTED) add option to sort Ag-R energy rank in ascending or descending.

update 2/23/2024: version 6
    - incoorporate MPNN prediction

update 2/21/2024: version 5
    - Disable repeated mutation on the same index.
    
update 2/20/2024: version 4
    - Allow a new configuration parameter to disable fastrelax

update (2/19/2024):
    - added new option in configuration file to only use given epitope window.
'''

#_____
# import
import os
import shutil

#from VE2_protein_helper import getEpitopeWindow

#_____
# variable
aa_list = ['G','I','F','A','P','L','V','W','Y','M','S','T','D','E','C','N','Q','R','H','K']
aa_mapping = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}
str_list = [
    "JobName:TempTestJob",
    "AbComplex:pdbs/LY-CoV1404_clean.pdb",
    "AbChain:AB_C",
    "A2Complex:pdbs/6M0J_clean.pdb",
    "A2Chain:A_E",
    "EpitopeWindow:444,371,368,446",
    "Mutations:3",
    "MaxGeneration:3",
    "useGivenWindow:Yes",
    "relaxPose:No",
    "AllowRepeatMutation:No",
    "UseMPNN:No",
    "EnergyWindow:-60,-10",
    "SortRankAg-Ab:Ascending",
    "SortRankAg-R:Ascending",
]

#_____
# functions

def getNewGenerationPath(base_job_path, path_list, mutations):
    new_gen = int(path_list[0].split('_')[0]) + 1
    new_var = int(path_list[1].split('_')[0]) * mutations

    #print(new_gen, new_var) # <-- development and debug
    new_gen_var_path_list = []

    for i in range(mutations):
        tmp_var = new_var + i
        new_gen_var_path_list.append(f"{base_job_path}/{str(new_gen)}_gen/{str(tmp_var)}_var")

    return new_gen_var_path_list

'''
if 'task' is 'continue'
    Given a CTRL-V job directroy, 
    create queue and past mutation to continue from previous job.

if 'task' is 'report'
    Given a CTRL-V job directory, 
    report the total mutation made.
'''
def createPreviousQueue(base_job_path, task = "continue"):
    queue = []
    past_mutation = []

    # v-- update (2/29/2024) get all mutation.
    all_simulation_mutations = []
    #done_simulation_mutations = [] # <-- redundant, past_mutation have the same element.

    dir_list = os.listdir(base_job_path)
    dir_list.sort(key=lambda x: int(x.split("_")[0]))
    #print(dir_list) # <-- development and debug
    for gen_dir in dir_list:
        tmp_gen_path = base_job_path + "/" + gen_dir
        tmp_gen_dir_list = os.listdir(tmp_gen_path)
        tmp_gen_dir_list.sort(key=lambda x: int(x.split("_")[0]))
        #print(tmp_gen_dir_list) # <-- development and debug
        for var_dir in tmp_gen_dir_list:
            tmp_var_path = tmp_gen_path + "/" + var_dir
            check_list = os.listdir(tmp_var_path)
            tmp_var_info_path = f'{tmp_var_path}/info'
            generation_info_dict = getGenerationInfo(tmp_var_info_path) # <-- update (2/29/2024) move outside of else statement
            #print(check_list) # <-- development and debug
            if 'done' not in check_list:
                queue.append(tmp_var_info_path) # <-- the variant file is not done, hence add to queue.

            else:
                #generation_info_dict = getGenerationInfo(tmp_var_info_path) # <-- update (2/21/2024, about past mutation mutation)
                past_mutation = past_mutation + generation_info_dict['MutateInfo'].split(",")
            all_simulation_mutations = all_simulation_mutations + generation_info_dict['MutateInfo'].split(",")
    
    if task == "continue":
        return queue, past_mutation
    if task ==  "report":
        return all_simulation_mutations, past_mutation

    raise ValueError(f"\nVariable 'task' needs to be one of the following:\n[continue, report]\nGiven value: {task}")

def getGenerationInfo(info_path, by_pass_space = False):
    rtn_dict = {}
    with open(info_path, 'r') as f:
        file_content = f.readlines()
        for line in file_content:
            if by_pass_space and " " in line:
                raise ValueError("\nIssue with reading Info file from generation.")
            line = line.replace('\n', '')
            rtn_dict[line.split(":")[0]] = line.split(":")[-1]
    return rtn_dict

def createGenerationInfo(gen_0_path, info_list):
    path_list = gen_0_path.split("/")
    if not os.path.exists('/'.join(path_list[-2:-1])): os.mkdir('/'.join(path_list[-2:-1]))
    if os.path.exists('/'.join(path_list[-2:])): raise ValueError(f"\nVariant directory already Exists, should not happen.\n{gen_0_path}")
    os.mkdir('/'.join(path_list[-2:]))
    info_file_path = f"{gen_0_path}/info"
    name_list = [
        'ParentAbComplex',
        'AbChain',
        'ParentA2Complex',
        'A2Chain',
        'Window',
        'Mutations',
        'MutateInfo'
    ]
    if len(name_list) != len(info_list): raise ValueError("\nMatching error creating generation info.")
    with open(info_file_path, 'w') as f:
        for i in range(len(name_list)):
            if type(info_list[i]) is list:  # <-- update (2/21/2024, automatically convert list of Mutation info to str, ie "1A,6B,19M")
                tmp = ','.join(info_list[i])
                f.write(f"{name_list[i]}:{tmp}\n")
            else:
                f.write(f"{name_list[i]}:{info_list[i]}\n")

    return info_file_path

'''
Print an example configuration file for ve workflow.

Make sure to update this function whenever configuration is updated.
'''
def printExampleConfigFile(): 
    # v-- depreciated, minor update (2/29/2024)
    #str = "JobName:TempTestJob\nAbComplex:pdbs/LY-CoV1404_clean.pdb\nAbChain:HL_C\nA2Complex:pdbs/6M0J_clean.pdb\nA2Chain:A_E\nEpitopeWindow:345,346,437,438,439,440\nMutations:5\nMaxGeneration:3\n"
    return "\n".join(str_list) + "\n"

'''
Get Information from the given config file.
'''       
def getInfoFromConfig(config_file):
    rtn_list = []
    with open(config_file, 'r') as f:
        file_content = f.readlines()
        for i, line in enumerate(file_content):
            if line == "\n": continue
            if " " in line:
                number = i + 1
                show = ""
                for x in line:
                    if x == " ": show += "^"
                    if x != " ": show += " "
                line = line + show
                raise ValueError(f"\nDetected Spacing in line {number}:\n{line}\nPlease remove any space in config file, including file name.\n\nExample config file:\n\n{printExampleConfigFile()}")
            line = line.replace('\n', '')
            rtn_list.append(line.split(":")[-1])

    # v-- update 3/4/2024 error handling if given paramter is lesser than required.
    if len(rtn_list) != len(str_list): raise ValueError(f"\nMissing parameter in config file.\nMake sure config file have the following format.\n\n{printExampleConfigFile()}")

    # 2, 4, 5, 7 need to be checked.
    if "_" not in rtn_list[2]: raise ValueError("\nMissing receptor chain or viral chain for Ab-Ag complex.\nRemember to add \'_\' between Ab chain and Ag chain, for example: HL_C.")
    if len(rtn_list[2].split("_")) != 2: raise ValueError("\nMake sure there is only two seperate protein.\n<Ab_chain>_<Ag_chain>, for example: HL_C.")
    if "_" not in rtn_list[4]: raise ValueError("\nMissing receptor chain or viral chain for ACE2-Ag complex.\nRemember to add \'_\' between ACE2 chain and Ag chain, for example: A_E.")
    if len(rtn_list[4].split("_")) != 2: raise ValueError("\nMake sure there is only two seperate protein.\n<ACE2_chain>_<Ag_chain>, for example: A_E.")
    if not rtn_list[6].isnumeric() : raise ValueError("\nMutations need to be numeric.")
    if not rtn_list[7].isnumeric() : raise ValueError("\nMax Generation needs to be numeric.")
    if rtn_list[8] not in ['Yes', 'No']: raise ValueError(f"\nUse given window needs to be [Yes, No].\ngiven value: {rtn_list[8]}")
    if rtn_list[9] not in ['Yes', 'No']: raise ValueError(f"\nRelax pose needs to be [Yes, No].\ngiven value: {rtn_list[9]}")
    # v-- update (2/21/2024), to allow repeat mutation on the same position
    if rtn_list[10] not in ['Yes', 'No']: raise ValueError(f"\nAllow repeat mutation needs to be [Yes, No].\ngiven value: {rtn_list[10]}")
    # v-- update (2/23/2024), use Window from ProteinMPNN, automate this in the future
    if rtn_list[11] not in ['', 'No']: print(f"WARNING: \nMake sure UseMPNN needs to be list of single point mutation,\ni.e. '338A,444T,444G' or 'No' or ''.\ngiven value: {rtn_list[11]}")
    
    # v-- update (2/24/2024), add energy Window.
    if len(rtn_list[12].split(",")) != 2: raise ValueError(f"\nEnergy window parameter needs to be only two value separated by comma.\nExample: '-100,100' in ascending.\ngiven value: {rtn_list[12]}")

    # v-- update (2/28/2024), add option to sort Ag-Ab energy rank and Ag-R energy rank
    if rtn_list[13] not in ['Ascending', 'Descending']: raise ValueError(f"\nSortRankAg-Ab needs to be ['Ascending', 'Descending'].\ngiven value: {rtn_list[13]}")
    if rtn_list[14] not in ['Ascending', 'Descending']: raise ValueError(f"\nSortRankAg-R needs to be ['Ascending', 'Descending'].\ngiven value: {rtn_list[14]}")

    return rtn_list

'''
Get Full information after parsing config file. 
'''
def getInfoAfterConfig(args_config, epitopeWindow):
    # v-- update (2/28/2024), future update to allow repeat mutation on the same position
    [job_name, ab_pdb, ab_chain, a2_pdb, a2_chain, window, mutations, max_gen, use_given_window, given_relax_pose, given_allow_repeat, use_MPNN, energy_window, sort_AgAb, sort_AgR] = getInfoFromConfig(args_config)

    use_window = None
    if use_given_window == "Yes": use_window = True
    if use_given_window == "No": use_window = False
    relax_pose = None
    if given_relax_pose == "Yes": relax_pose = True
    if given_relax_pose == "No": relax_pose = False

    # v-- update (2/21/2024), future update to allow repeat mutation on the same position
    allow_repeat = None
    if given_allow_repeat == "Yes": allow_repeat = True
    if given_allow_repeat == "No": allow_repeat = False

    # v-- update (2/23/2024), added MPNN parameter.
    MPNN = None
    # v-- Future, automate MPNN is wanted to,
    # This parameter is after using MPNN without removing known aa, then add all IRL mutation into it.
    # ** NOTE: MANUAL MPNN CHANGE HERE **
    # Run protein_mpnn.py, then it will produce a result file named: protein_mpnn_result.txt, replace below with "Guesses" from MPNN.
    
    # v-- original MPNN prediction
    #if use_MPNN == "Yes": MPNN = {'339': {'G', 'K', 'N', 'A', 'T', 'D', 'S', 'H', 'Q', 'E'}, '346': {'K', 'N', 'T', 'D', 'S', 'E'}, '368': {'V', 'I', 'L'}, '371': {'N', 'F', 'L', 'D'}, '373': {'K', 'R', 'N', 'A', 'T', 'D', 'S', 'P'}, '375': {'K', 'N', 'F', 'T', 'D', 'S', 'I', 'V', 'E'}, '376': {'K', 'R', 'T', 'A', 'H', 'S', 'M', 'I', 'L', 'V', 'Q', 'E'}, '405': {'K', 'R', 'N', 'F', 'D', 'Y', 'S', 'H', 'C', 'Q', 'E'}, '408': {'G', 'K', 'R', 'A', 'T', 'D', 'S', 'V', 'Q', 'E', 'P'}, '417': {'K', 'R', 'N', 'T', 'D', 'Y', 'S', 'M', 'I', 'L', 'V', 'Q', 'E', 'P'}, '439': {'K', 'R', 'A', 'T', 'S'}, '440': {'G', 'K', 'N', 'A', 'D', 'H', 'S', 'E'}, '444': {'K', 'R', 'A', 'T', 'Q'}, '445': {'K', 'N', 'A', 'T', 'D', 'S', 'I', 'L', 'V', 'Q', 'E', 'P'}, '446': {'G', 'S', 'N'}, '452': {'K', 'R', 'T', 'L', 'Q', 'E'}, '456': {'Y', 'F', 'H', 'L'}, '460': {'K', 'A', 'P'}, '477': {'K', 'R', 'N', 'A', 'T', 'D', 'S', 'V', 'E'}, '478': {'K', 'R', 'A', 'T', 'I', 'L', 'V', 'Q', 'E'}, '484': {'K', 'T', 'A', 'I', 'L', 'V', 'Q', 'E'}, '486': {'K', 'R', 'A', 'T', 'D', 'S', 'V', 'E', 'P'}, '490': {'K', 'R', 'A', 'T', 'D', 'Y', 'H', 'S', 'L', 'Q', 'E', 'P'}, '493': {'K', 'R', 'A', 'T', 'S', 'V', 'Q', 'E'}, '496': {'K', 'R', 'N', 'A', 'T', 'H', 'S', 'C'}, '498': {'K', 'N', 'R', 'D', 'H', 'S', 'L', 'Q'}, '501': {'A', 'D', 'Y', 'S', 'Q', 'E'}, '505': {'K', 'R', 'N', 'D', 'H', 'L', 'Q', 'E'}}
    
    # 3/26/2024, remove this one and use the original MPNN prediction instead. 
    # this is for checking difference between MPNN and without
    #if use_MPNN == "Yes": MPNN = {'484': {'F', 'W'}, '439': {'G', 'T', 'Q', 'A', 'S', 'C'}, '498': {'I', 'M', 'F'}, '493': {'T', 'N', 'H'}, '501': {'P', 'Q', 'S', 'A', 'V'}, '505': {'N', 'W'}, '444': {'T', 'A', 'P', 'N'}, '440': {'G', 'R', 'K', 'A', 'S', 'D', 'H'}, '346': {'I', 'D', 'H', 'E'}, '445': {'C'}, '405': {'E', 'Y'}, '376': {'V', 'Q'}, '375': {'V', 'I'}, '371': {'L'}, '490': {'Y'}, '452': {'R'}, '368': {'I'}, '373': {'P', 'N', 'R'}}
    if use_MPNN == "No" or use_MPNN == "": 
        MPNN = None
    else:
        MPNN = use_MPNN.split(",")

    print(f"\n--- Starting job: {job_name} ---")

    print("Antibody Complex Path:\n  ", ab_pdb)
    print("Antibody chain:\n  ", ab_chain)
    print("ACE2 Complex Path:\n  ", a2_pdb)
    print("ACE2 chain:\n  ", a2_chain)
    print("Number of viral mutations:\n  ", mutations)
    print("Only use given epitope window:\n  ", use_window)

    if not window: 
        window1 = epitopeWindow(ab_pdb, ab_chain)
        #print(window1)
        window2 = epitopeWindow(a2_pdb, a2_chain)
        #print(window2)
        tmp = window1 + window2
        #print(tmp)
        tmp = set(tmp)
        #print(tmp)
        window = list(tmp)
        window.sort()
        #print(window)
        #print(",".join(window))
    else: # <-- always use given window for first iteration.
        window = window.split(",")
        if "" in window: window.remove("")

    print("Epitope window:\n  ", ",".join(window))
    print("Max generations:\n  ", max_gen)
    
    print("Allow relax pose:\n  ", relax_pose)
    # v-- update (2/21/2024), to allow repeat mutation on the same position
    print("Allow repeat mutation:\n  ", allow_repeat)
    # v-- update (2/21/2024), determine if we are using MPNN prediction
    print("Allow use MPNN:\n  ", use_MPNN)
    # v-- update (2/24/2024), energy window
    print("Energy window information:\n  ", energy_window)

    # v-- update (2/28/2024), sorting option on energy ranking (sort_AgAb, sort_AgR)
    print("Energy sorting for Antigen-Antibody Energy Ranking:\n  ", sort_AgAb)
    print("Energy sorting for Antigen-Receptor Energy Ranking:\n  ", sort_AgR)
    if sort_AgAb == "Ascending": sort_AgAb = False # 'Ascending', 'Descending'
    if sort_AgAb == "Descending": sort_AgAb = True # 'Ascending', 'Descending'
    if sort_AgR == "Ascending": sort_AgR = False # 'Ascending', 'Descending'
    if sort_AgR == "Descending": sort_AgR = True # 'Ascending', 'Descending'

    print("--- ---")
    # v-- update (2/21/2024), future update to allow repeat mutation on the same position
    return job_name, ab_pdb, ab_chain, a2_pdb, a2_chain, window, mutations, int(max_gen), use_window, relax_pose, allow_repeat, MPNN, energy_window, sort_AgAb, sort_AgR

'''
Resolve job directory issue where the job_name exists.
'''
def resolveJobDirectory(job_name, status):
    continue_status = False
    # Resolve job directory issue.
    if os.path.exists(job_name):
        print("Detected previous job with the same directory name")
        if status == 'restart':
            print("Attempting to remove the previous job and start a new job.")
            shutil.rmtree(job_name)
            os.mkdir(job_name)
            #raise ValueError("\n--status restart is not implemented yet.")

        elif status == 'continue':
            print("Attempting to continue from the previous job.")
            continue_status = True
            #raise ValueError("\n--status continue is not implemented yet.")
        else:
            raise ValueError("\nJob directory exists,\nuse --status restart, to remove previous job,\nuse --status continue, to continue previous job")
    else:
        os.mkdir(job_name)

    base_job_path = os.getcwd() + '/' + job_name
    return base_job_path, continue_status

'''
Get Mutation only from csv file

for example:
    return set --> {444T, 399I, ...}
'''
def getMutationOnlyInfo(csv_file = "/work/ratul1/chuen/viral_escape/pyRosetta_script/CoV_IRL_mutation.csv", 
                        task = "none"):
    tasks = ["mutations", "irlDictionary"]

    contents = []
    with open(csv_file, 'r') as f:
        contents = f.readlines()
        for i in range(len(contents)):
            row = contents[i].replace('\n', '')
            row_list = row.split(',')
            contents[i] = row_list

    mutant_set = set()
    mutant_dict = {}
    for i in range(len(contents[0])):
        mutant_name = ''
        for j in range(len(contents)):
            # v-- this is for full irl info dictionary 
            if j == 0:
                mutant_dict[contents[j][i]] = []
                mutant_name = contents[j][i]
            elif contents[j][i] != '':
                mutant_dict[mutant_name].append(contents[j][i])

            # v-- this is for mutantion only    
            if contents[j][i] == "" or j == 0: continue
            mutant_set.add(contents[j][i][3:])

    if task == tasks[0]: return mutant_set 
    if task == tasks[1]: return mutant_dict

    raise ValueError(f"\nTask options needs to be in {tasks},\nGiven: {task}")
