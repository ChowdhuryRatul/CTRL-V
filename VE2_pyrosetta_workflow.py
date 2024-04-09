# Author: Yee Chuen Teoh
'''
curr version: 7.0 (4/8/2024)

requirement:
    - make sure that PyRosetta module is loaded
    - make sure to clean PDB before using pdb, 
        --> just use cleanATOM(<file.pdb>) as example in VE2_pyrosetta_test.py
    - make sure to enter the chain information correctly.

required helpers:
    - VE2_directory_handling
    - VE2_rosetta_helper

known issue(irrelevant to CTRL-V):
    - sometimes when running sbatch job, it will abort due to the following error:
        FATAL:   aborting failed to close file descriptor: bad file descriptor
        this seems to be error within the sbatch itself, resubmitting the sbatch script would be alright.

update 2/28/2024:
    - (FIXED) issue where if a variant has no potential mutation, it does not get marked as done in the directory.

Usage:
** to restart job **
python VE2_pyrosetta_workflow.py --config ve_config.ve --status restart
python VE2_pyrosetta_workflow.py --config test.ve --status restart

** to continue previous job **
python VE2_pyrosetta_workflow.py --config ve_config.ve --status continue
python VE2_pyrosetta_workflow.py --config test.ve --status continue


python VE2_pyrosetta_workflow.py --config ve_configIRL_v6_report.ve --status continue
'''

#_____
# import
import os
import shutil
import sys
import random
import argparse

# from VE2_protein_helper import performMutation
from VE2_directory_handling import getInfoAfterConfig, resolveJobDirectory, createPreviousQueue, createGenerationInfo, getGenerationInfo, getNewGenerationPath
from VE2_rosetta_helper import performMutationPyRosetta, getEpitopeWindowPyRosetta

#_____
# variable
aa_list = ['G','I','F','A','P','L','V','W','Y','M','S','T','D','E','C','N','Q','R','H','K']

#_____
# functions (directory handler here.)
def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=str, required=True)
    parser.add_argument('--status', type=str, required=False)
    '''
    config file that contain
    AbComplex:<path/to/AntibodyComplex>
    AbChain:<AbChain_AgChain>  
    A2Complex:<path/to/Ace2Complex>
    A2Chain:<A2Chain_AgChain> 
    Window:<1,2,6,7,10,12,...,n>
    NOTE**: top example is depreciated, check example below.

    v-- updated example.
    example config file (C is viral chain, E is viral chain):
    JobName:TempTestJob
    AbComplex:pdbs/LY-CoV1404_clean.pdb
    AbChain:AB_C
    A2Complex:pdbs/6M0J_clean.pdb
    A2Chain:A_E
    EpitopeWindow:444,371,368,446
    Mutations:3
    MaxGeneration:3
    useGivenWindow:Yes
    relaxPose:No
    AllowRepeatMutation:No
    UseMPNN:No
    EnergyWindow:-60,-10
    SortRankAg-Ab:Ascending
    SortRankAg-R:Ascending
    '''
    args = parser.parse_args()
    return args

#_____
# main
if __name__ == "__main__":
    args = parser()
    job_name, ab_pdb, ab_chain, a2_pdb, a2_chain, window, mutations, max_gen, use_window, relax_pose, allow_repeat, MPNN, energy_window, sort_AgAb, sort_AgR = getInfoAfterConfig(args.config, getEpitopeWindowPyRosetta)
    ab_pdb = os.getcwd() + '/' + ab_pdb
    a2_pdb = os.getcwd() + '/' + a2_pdb
    #print(ab_pdb, ab_chain, a2_pdb, a2_chain, window)

    base_job_path, continue_status = resolveJobDirectory(job_name, args.status)
    os.chdir(base_job_path) # <-- always reset to base directory

    queue = []
    past_mutation = []
    if continue_status: # <-- handle continue from previous job.
        queue, past_mutation = createPreviousQueue(base_job_path)
    else:
        gen_0_path = f"{base_job_path}/0_gen/0_var"
        gen_0_info_file = createGenerationInfo(gen_0_path, [ab_pdb, ab_chain, a2_pdb, a2_chain, window, mutations, "WildType"])

        queue = [gen_0_info_file] # <-- here set the default queue if the job is restart or new.
        
    print("") # <-- spacing to start loop.

    max_gen_limit = max_gen + 1
    #max_gen = 5 # <-- development and debug
    total_mutation = set(past_mutation)
    while queue:
        os.chdir(base_job_path) #<-- always reset to base directory
        to_mutate_info_file_path = queue.pop(0)

        path_list = to_mutate_info_file_path.split('/')
        to_mutate_gen_path_list = path_list[-3:-1]
        to_mutate_gen_path = '/'.join(to_mutate_gen_path_list)

        if max_gen_limit <= int(to_mutate_gen_path.split("_")[0]): # <-- as long as directory generation is larger than max_gen_limit
            print(f"\nReached generation limit of: {max_gen} generations\n")
            break

        print(f"----- Working on generation {to_mutate_gen_path} ------")
        os.chdir(to_mutate_gen_path)

        to_mutate_info_file_dict = getGenerationInfo(to_mutate_info_file_path)
        # v-- perform the mutation here, Ratul work flow here.
        new_mutations_info_list = performMutationPyRosetta(to_mutate_info_file_dict, use_given_window = use_window, relax_pose = relax_pose, allow_repeat = allow_repeat, MPNN = MPNN, energy_window = energy_window, sort_AgAb = sort_AgAb, sort_AgR = sort_AgR) 
        # v-- this is for directory creation.
        new_gen_var_path_list = getNewGenerationPath(base_job_path, to_mutate_gen_path_list, int(mutations))
        # v-- update (2/21/2024) to remember new mutation.
        total_mutation.add(to_mutate_info_file_dict['MutateInfo'].split(",")[-1])
        print("Current total mutation:", total_mutation) # <-- development and debug purposes.

        #print(new_gen_var_path_list) # <-- development and debug

        # create generation directory.
        os.chdir(base_job_path) #<-- always reset to base directory
        # v-- remove this testing because sometimes it generates less than top-k when we restrict repeating aa
        #if len(new_mutations_info_list) != len(new_gen_var_path_list):
        #    raise ValueError("\nNumber of mutations does NOT match number of path generated.")

        # create generation info file for new mutation
        for i in range(len(new_mutations_info_list)):
            gen_0_info_file = createGenerationInfo(new_gen_var_path_list[i], new_mutations_info_list[i])

        # v-- added a 'done' flag as long as this variant has gone through the workflow.
        #       regardless of if it produces any mutation suggestion.
            queue.append(new_gen_var_path_list[i] + '/info')
        with open(f"{to_mutate_gen_path}/done", 'w') as f: # <-- tracker to track if a variant has gone through the workflow
            f.write("mutation workflow complete.")
        print("")

    total_mutation.remove('WildType') # <-- WildType is not considered mutation, remove from the reporting set.
    print(total_mutation)
