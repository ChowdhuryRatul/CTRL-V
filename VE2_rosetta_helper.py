# Author: Yee Chuen Teoh
# Description: helper functions to be used in VE2_rosettaworkflow.py
# NOTE: This script can only run with PyRosetta Module in HPC.
'''
update 4/8/2024: version 7
    - (NOT TESTED YET) update to MPNN to be a list of singlepoint mutation.

update 3/26/2024
    - (Error) check if it key in MPNN

update 2/28/2024: version 6
    - (TESTED) add option to sort Ag-Ab energy rank in ascending or descending.
    - (TESTED) add option to sort Ag-R energy rank in ascending or descending.
    - (TESTED) add function to get all possible mutation, used in VE2_report.py.

update 2/24/2024: version 6
    - (TESTED) added energy window, previously is a constant, now it is a parameter in the config file.

update 2/22/2024: version 6
    - (TESTED) added ProteinMPNN prediction. before ranking, 
        ** NOTE: MPNN is not automated, needed to be changed manually in VE2_directory_handling **

update 2/21/2024: version 5
    - (TESTED) Disable repeated mutation on the same index.

update 2/20/2024: version 4
    - Not to reverse ACE2 ranking after seeing scatterplot.
    - (TESTED) Since mutations energy are always within a treshold, we use disregard any energy outside of it.
    - (TESTED) Allow a new configuration parameter to disable fastrelax

update 2/19/2024: version 3
    - (TESTED) mutate residue only mutate when PyRosetta actually changes it.
        sometimes PyRosetta doesn't do certain mutation on certain residue, reasoning is not investigated
    - (TESTED) Ranking residue now ignores differences in Ab and ACE2, 
        only looks at residue ranking in Ab, and ACE2.
    - (TESTED) Binding Energy now uses IAM get seperated interface energy.
    - (FIXED) Fix issue where splitPdb takes in partner_chain, it only takes in viral chain list. 
    - (TESTED) Using relax pose.
            -> Relax pose during in mainMutation.
            -> Relax pose if getAllMutation uses non-relax pose.

issue:
    - (FIXED) unable to relax poses. PyRosetta issue.

'''

#____________________________________________________________________________________________________
# imports 
import argparse
import os
import sys
import math
import random # <-- development and debug purposes to randomize mutation.
import statistics
#import numpy as np

# PyRosetta
from pyrosetta import *
init()
from pyrosetta.toolbox import *
from pyrosetta.teaching import * # <-- for score function
from pyrosetta.rosetta.protocols import *
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.core.select import *
from pyrosetta.rosetta.protocols.relax import FastRelax

#____
# helper variable
aa_list = ['G','I','F','A','P','L','V','W','Y','M','S','T','D','E','C','N','Q','R','H','K']
aa_mapping = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

#_____
# Function helper
'''
Get the distance information from pdb_file.
'''
def calculateDistancePyRosetta(pdb_file):
    pose = pose_from_pdb(pdb_file)
    print(pose.pdb_info())

    distance_dict = {}

    for i in range(1, pose.total_residue() + 1):
        a_chain = pose.pdb_info().chain(i)
        a_number = pose.pdb_info().number(i)
        # v-- currently uses Carbon Atom, 
        # unable to find Omega atom in pose.residue(i), N = nitrogen, C = carbon, O = oxygen, unsure what is the omega atom.
        a_coor = pose.residue(i).atom("CA").xyz() # <-- change ATOM to compare.

        a_name = f"{a_chain}{a_number}"

        for j in range(1, pose.total_residue() + 1):
            b_chain = pose.pdb_info().chain(j)
            b_number = pose.pdb_info().number(j)
            b_coor = pose.residue(j).atom("CA").xyz() # <-- change ATOM to compare.


            b_name = f"{b_chain}{b_number}"

            dist = a_coor.distance(b_coor)

            distance_dict[(a_name, b_name)] = dist

    #print(distance_dict) # <-- development and debug purposes.
    return distance_dict


'''
Create epitope window from complex file

input:
    pdb_file --> path to the complex pdb file
    chain --> "AbChain_AgChain" for example: 'HL_C'
'''
def getEpitopeWindowPyRosetta(pdb_file, chain, angstrom = 10, add_window = 1):
    [_, viral_chain] = chain.split("_")

    dist_dict = calculateDistancePyRosetta(pdb_file)
    viral_chain_list = list(viral_chain)

    #print(dist_dict)

    epitope_window = set()

    for tuple in dist_dict:
        if any(chain in tuple[0] for chain in viral_chain_list) and not any(chain in tuple[1] for chain in viral_chain_list):
            if float(dist_dict[tuple]) <= angstrom:
                epitope_window.add(tuple[0][1:])
                number = int(tuple[0][1:])
                for i in range(1, add_window + 1):
                    epitope_window.add(str(number + i))
                    epitope_window.add(str(number - i))

    epitope_window_list = list(epitope_window)
    epitope_window_list.sort()

    return epitope_window_list


'''
(Not used)
split pdb complex file into two pdb file based on chain.
'''
def splitPdb(pdbFile, chain2, name1, name2):
    fileContent = ''
    with open(pdbFile, 'r') as f:
        fileContent = f.read()
    fileContentList = fileContent.split('\n')
    newList = []

    saveName_1 = name1 + '.pdb'
    pdbList_1 = []
    saveName_2 = name2 + '.pdb'
    pdbList_2 = []

    for string in fileContentList:
        if string[0:4] == 'ATOM':
            chain = string[21]
            if chain in chain2:
                pdbList_2.append(string)
            else:
                pdbList_1.append(string)
    
    with open(saveName_1, 'w') as f:
        for line in pdbList_1:
            f.write(line + '\n')
            
    with open(saveName_2, 'w') as f:
        for line in pdbList_2:
            f.write(line + '\n')

    return saveName_1, saveName_2


'''
(Not used)
unbind the Pose complex.

partner_chain is the ab_chain/a2_chain information, for example HL_C.
'''
def unbindComplex(pose, partner_chain):
    STEP_SIZE = 100
    JUMP = 2
    docking.setup_foldtree(pose, partner_chain, Vector1([-1,-1,-1]))
    trans_mover = rigid.RigidBodyTransMover(pose,JUMP)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose)


'''
(Not used)
Get the binding energy of the complex, score(complex) - (score(anti) + score(viral))

3 functions below uses different implementation.
'''
def getBindingEnergyUnbind(pose_complex, partner_chain):
    sfxn = get_score_function(True)  

    bound_binding_energy = sfxn(pose_complex)
    #print(bound_binding_energy) # <-- development and debug purposes

    unbound_pose_complex = Pose()
    unbound_pose_complex.assign(pose_complex)
    unbindComplex(unbound_pose_complex, partner_chain)

    unbound_binding_energy = sfxn(unbound_pose_complex)
    #print(unbound_binding_energy) # <-- development and debug purposes

    return bound_binding_energy - unbound_binding_energy

'''
(Not used)
'''
def getBindingEnergy(pose_complex, pose_viral, pose_recep):
    sfxn = get_score_function(True)  

    bound_binding_energy = sfxn(pose_complex)
    #print(bound_binding_energy) # <-- development and debug purposes

    unbound_binding_energy = sfxn(pose_viral) + sfxn(pose_recep)
    #print(unbound_binding_energy) # <-- development and debug purposes

    return bound_binding_energy - unbound_binding_energy

def getBindingEnergyIAM(pose, partner_chain):
    # Create an instance of InterfaceAnalyzerMover
    interface_analyzer = rosetta.protocols.analysis.InterfaceAnalyzerMover()

    # Set the pose for analysis
    interface_analyzer.set_input_pose(pose)

    interface_analyzer.set_interface(partner_chain)

    # Analyze the interface
    interface_analyzer.apply(pose)

    # Get the interface energy# Get the dG_separated value
    energy = interface_analyzer.get_separated_interface_energy()
    # get_interface_dG
    # get_separated_interface_energy()
    #energy_crossterm = interface_analyzer.get_crossterm_interface_energy()

    # v-- developement and debug purposes.
    #sfxn = get_score_function(True)  
    #bound_binding_energy = sfxn(pose)
    # Print the interface energy
    #print("Energy IAM dG: ", energy_dG)
    #print("Energy IAM crossterm: ", energy_crossterm)
    #print("Energy from sfxn: ", bound_binding_energy)

    return energy

'''
Relax the given pose.
'''
def relaxPose(pose):
    # Create a FastRelax protocol
    fast_relax = rosetta.protocols.relax.FastRelax()
    # Set up the FastRelax options
    fast_relax.set_scorefxn(rosetta.core.scoring.get_score_function(True))
    fast_relax.constrain_relax_to_start_coords(True)  # Optional: Constrain the relaxation to the starting coordinates
    # Relax the pose
    fast_relax.apply(pose)
    return pose

'''
(Not used)
Get poses for complex, antibody and viral
'''
def getPoses(pdb_file, partner_chain):
    # (Depreciated: Viral and Recep file) updated 2/17/2024
    name = ''.join(pdb_file.split(".")[:-1])


    recep_file = f"{name}_anti"
    viral_file = f"{name}_viral"
    recep_file, viral_file = splitPdb(pdb_file, partner_chain.split("_")[1], recep_file, viral_file)

    pose_complex = pose_from_pdb(pdb_file)  
    pose_recep = pose_from_pdb(recep_file)  
    pose_viral = pose_from_pdb(viral_file)

    #print(pose_complex.pdb_info()) # <-- development and debug
    #print(pose_recep.pdb_info()) # <-- development and debug
    #print(pose_viral.pdb_info()) # <-- development and debug

    #relax.apply(pose_complex)
    #relax.apply(pose_recep)
    #relax.apply(pose_viral)

    return pose_complex, pose_viral, pose_recep

'''
Get one pose for input pdb.
'''
def getPose(pdb_file):
    # updated 2/19/2024
    pose_complex = pose_from_pdb(pdb_file) 
    return pose_complex

'''
given pdb, with chain to find
the binding score for each mutation at each position.

parameter:
    --> pack_range: this parameter is used for packing during mutation, default to 0, 
                    indicate 0 angstrom to be repack with PyRosetta after mutate_residue
    --> toadd: Depreciated. can be deleted.
    --> debug: only set to true during development and debug, can be deleted.
    --> MPNN: ProteinMPNN prediction window. 

return a dictionary where 
    --> key is the residue index with mutated aa.
    --> value is the binding score.
'''
def getAllMutationScore(pdb_name, partner_chain, epitope_window, 
                        pack_range = 0, toadd = 0, debug = False, MPNN = None, energy_window = [-100,100]):
    [_, viral_chain] = partner_chain.split("_")
    viral_chain_list = list(viral_chain)

    # v-- update 2/24/2024, energy window
    lower_E = energy_window[0]
    higher_E = energy_window[1]

    # v-- updated version 2/19/2024
    pose_complex = getPose(pdb_name)

    # v-- (Depreciated) older version
    #pose_complex, pose_viral, pose_recep = getPoses(pdb_name, viral_chain_list)
    #ori_seq_complex = pose_complex.sequence() # <-- development and debug
    #ori_seq_viral = pose_viral.sequence() # <-- development and debug

    #binding_score = getBindingEnergy(pose_complex)
    #print(pose_complex.pdb_info()) # <-- development and debug
    #print(pose_viral.pdb_info()) # <-- development and debug
    #print('pdb length: ', pose_viral.total_residue()) # <-- development and debug

    viral_idx_in_complex = -1 # <-- the viral pose numbering in the complex pose. number in pose =/= number in pdb.
    viral_size = 0 # <-- get the viral sequence size.
    for i in range(pose_complex.total_residue()):
        pose_idx = i + 1
        # get the size of viral sequence in antigen
        if pose_complex.pdb_info().chain(pose_idx) in viral_chain_list:
            if viral_idx_in_complex == -1:
                viral_idx_in_complex = pose_idx
            #print(pose_idx) # <-- development and debug
            viral_size += 1
            
    #print(pose_complex.pdb_info().number(viral_idx_in_complex)) # <-- development and debug

    # v-- new updated 2/19/2024 comment/remove all binding score related lines.
    #all_binding_score = []

    scoring_dict = {}
    for i in range(1, viral_size + 1):
        if debug and i > 5: break # <-- development and debug purposes.
        pose_idx = i
        for aa in aa_list:
            # (Depreciated, updated 2/19/2024)
            #if pose_viral.pdb_info().number(pose_idx) != pose_complex.pdb_info().number(viral_idx_in_complex):
            #    print(pose_idx, viral_idx_in_complex) # <-- development and debug purposes.
            #    print(pose_viral.pdb_info().chain(pose_idx)) # <-- development and debug purposes.
            #    print(pose_viral.pdb_info().number(pose_idx)) # <-- development and debug purposes.
            #    print(pose_complex.pdb_info().number(viral_idx_in_complex)) # <-- development and debug purposes.
            #    raise ValueError("\nIncorrect index mapping between viral and complex.")
            
            # get the pdb numbering.
            curr_numbering = str(pose_complex.pdb_info().number(viral_idx_in_complex))

            energy = None
            key = curr_numbering + aa

            # v-- update (2/23/2024) 
            # adding MPNN window
            # key[-1] not in MPNN[key[:-1]]
            # v-- update (3/26/2024)
            # key[:-1] is not a key in MPNN and

            # v-- update (4/8/2024)
            # update MPNN to be a list of single point mutation
            # i.e. MPNN = [339A, 444T, ...], where the user will input, 339A,444T,... in the config file.
            if MPNN:
                if key not in MPNN:
                    continue # <-- skip if there is no such key in MPNN.
                #if key[-1] not in MPNN[key[:-1]]: 
                #    continue # <-- skip if the mutation is not in MPNN's prediction.
            elif curr_numbering not in epitope_window:
                #print("Not in window issue", curr_numbering) # <-- development and debug purposes
                continue # <-- skip if the pdb number is not part of the window.

            # v-- chain checker
            if pose_complex.pdb_info().chain(viral_idx_in_complex) not in viral_chain_list:
                print("Given pose number is not residue from viral chain.")
                continue

            # create pose for mutation.
            #mutate_viral = Pose() # <-- Depreciated
            mutate_complex = Pose()

            #print("Performing mutation:", key) # <-- development and debug
            try:
                #mutate_viral.assign(pose_viral) # <-- Depreciated <-- copy poses to do mutation on.
                mutate_complex.assign(pose_complex) # <-- copy pose
                #print(mutate_viral.pdb_info()) # <-- development and debug
                #print(mutate_complex.pdb_info()) # <-- development and debug

                # general score function for mutation
                sfxn = get_score_function(True)
                # v-- low resolution mutation
                #mutate_residue(mutate_viral, pose_idx, aa, pack_range, sfxn)

                # v-- older version
                mutate_residue(mutate_complex, viral_idx_in_complex, aa, pack_range, sfxn)
                
                # v-- new updated 2/19/2024
                # found PyRosetta perform successfully one certain position. 
                # for none successful ones, we ignore.
                # check to make sure it is a different res from original
                if mutate_complex.residue(viral_idx_in_complex).name1() == pose_complex.residue(viral_idx_in_complex).name1():
                    print(f"Fail mutation from PyRosetta, we ignore this mutation: {key}")
                    continue

                # v-- new version as of 2/19/2024, uses IAM to calculated energy.
                energy = getBindingEnergyIAM(mutate_complex, partner_chain)
                
                # v-- new verion as of 2/17/2024 base of
                # https://nbviewer.org/github/RosettaCommons/PyRosetta.notebooks/blob/master/notebooks/06.08-Point-Mutation-Scan.ipynb
                # doesn't seem to work
                #packMutation(mutate_complex, viral_idx_in_complex, aa, sfxn)
                
                # v-- new version
                #energy = getBindingEnergyUnbind(mutate_complex, partner_chain)

                # v-- updated (2/20/2024) 
                # Ab energy threshold: -10 - -60
                # ACE2 energy threshold: -20 - -60
                # use treshold -10 - -60
                # since most SARS COV 2 are seen to be within this window.
                if not (lower_E < energy and energy < higher_E): 
                    print(f"Mutation out of energy window [{lower_E},{higher_E}]:", key, energy) # <-- development and debug
                    continue # <-- we skip energy outside of this threshold.

                #v-- debug to test error
                if debug and i == 5:
                    raise ValueError("Random error")
                #all_binding_score.append(energy)

                # v-- updated (2/19/2024)
                # v-- only successful mutation gets added to scoring dict.
                scoring_dict[key] = energy
                #print("Successful mutation on:", key) # <-- development and debug

            except Exception as error:
                print("Erorr message:", error)
                energy = "Error"
                # v-- new updated 2/19/2024 to ignore error mutation
                continue

            # (depreciated.)
            #scoring_dict[key] = energy

            #print(scoring_dict) # <-- development and debug
            # (depreciated) reset old to original amino acid. 
            #mutate_residue(pose_viral, pose_idx, ori_aa_viral)
            #mutate_residue(pose_complex, viral_idx_in_complex, ori_aa_complex)
            #if ori_seq_complex != pose_complex.sequence() or ori_seq_viral != pose_viral.sequence():
            #    raise ValueError("\nIncorrect mutate back to original strand.")
        viral_idx_in_complex += 1

    #all_binding_score.sort()
    #print(all_binding_score) # <-- development and debug
    #for key in scoring_dict:
    #    if scoring_dict[key] == "Error": # v-- assign median of all binding energy to error sequence.
    #        scoring_dict[key] = statistics.median(all_binding_score)

    return scoring_dict

'''
Get the IRL mutation information from csv file,

Return a dictionary where key = mutation_name : value = [S:{ori_aa}{index}{mutate_aa}, S:{ori_aa}{index}{mutate_aa}, ...]

(2/19/2024) Not used in PyRosetta workflow.
'''
def getMutationInfo(csv_file):
    contents = []
    with open(csv_file, 'r') as f:
        contents = f.readlines()
        for i in range(len(contents)):
            row = contents[i].replace('\n', '')
            row_list = row.split(',')
            contents[i] = row_list

    mutant_dict = {}
    for i in range(len(contents[0])):
        mutant_name = ''
        for j in range(len(contents)):
            if j == 0:
                mutant_dict[contents[j][i]] = []
                mutant_name = contents[j][i]
            elif contents[j][i] != '':
                mutant_dict[mutant_name].append(contents[j][i])

    return mutant_dict

'''
Save the ranking dictionary generated by getAllMutationScore

(2/19/2024) Used to save the energy ranking. 
            Used to plot scatter plot as well.
'''
def saveRankDict(dict1, savename):
    sorted_dict = sorted(dict1.items(), key=lambda item: item[0])
    with open(savename, 'w') as f:
        for tuple in sorted_dict:
            mutation = tuple[0]
            energy_aa = tuple[1]
            to_write = f'{mutation}:{energy_aa}\n'
            f.write(to_write)


'''
Main mutation used at the start of performMutationPyRosetta to have high resolution complex.

mutate_info(structure) = "478D" for example.
'''
def mainMutationPyRosetta(pdb_file, mutate_info, new_pdb_tosave, chain, dist = 15, relax_pose = True):
    [_, viral_chain] = chain.split("_")
    viral_chain_list = list(viral_chain)

    pose = pose_from_pdb(pdb_file) 
    print(pose.pdb_info())

    print(mutate_info)
    mutation_idx = int(mutate_info[:-1]) # <-- this is pdb numbering not pose numbering.
    mutation_res = mutate_info[-1]
    #print(mutation_idx, mutation_res) # <-- development and debug purposes.

    pose_numbering = -1
    for i in range(1, pose.total_residue() + 1):
        if int(pose.pdb_info().number(i)) == mutation_idx:
            if pose.pdb_info().chain(i) not in viral_chain_list:
                continue
                #raise ValueError("\nAttemping to mutate non viral residue.\nError chain matching during main mutation.")
            pose_numbering = i # <-- this is the index of the of the viral_chain

    # v-- development and debug purposes
    #print("residue numbering", pose.pdb_info().number(pose_numbering))
    #print("residue chain", pose.pdb_info().chain(pose_numbering))
    #print("residue info", pose.residue(pose_numbering).name1(), pose.residue(pose_numbering).name3())

    sfxn = get_score_function(True)
    
    # v-- older version.
    mutate_residue(pose, pose_numbering, mutation_res, dist, sfxn)
    # v-- new version as of 2/17/2024
    #packMutation(pose, pose_numbering, mutation_res, sfxn)

    # v-- development and debug purposes
    #print("residue numbering after Mutate", pose.pdb_info().number(pose_numbering))
    #print("residue chain after Mutate", pose.pdb_info().chain(pose_numbering))
    #print("residue info after Mutate", pose.residue(pose_numbering).name1(), pose.residue(pose_numbering).name3())

    # v-- update: 2/20/2024, only relax the pose if we allow it.
    # v-- update: 2/19/2024, relax the pose before dump_pdb
    if relax_pose: pose = relaxPose(pose)

    dump_pdb(pose, new_pdb_tosave)
    
    # v-- development and debug purposes
    #return pose.residue(pose_numbering).name1()

'''
Get ranking after mutation done by PyRosetta on Ab and ACE2

Ranking based on the binding energy rank in antibody and ACE2. 

(depreciated)
older version does this: both ascending
    reverse_state = True
inverse = True
if reverse_state: inverse = False
sorted_ab_list = sorted(... reverse = inverse)
sorted_a2_list = sorted(... reverse = inverse)
'''
# v-- update option to sort ascending or descending. False = ascending, True = descending
def getRankingFromScoringDict(ab_dict, a2_dict,
                                rankdiff = False, sort_AgAb = False, sort_AgR = False, sort_diff = False): 
    sorted_ab_list = sorted(ab_dict.items(), key=lambda item: item[1], reverse = sort_AgAb)
                                                    # v-- updated (2/20/2024) to not reverse ACE2 due to scatterplot.
    sorted_a2_list = sorted(a2_dict.items(), key=lambda item: item[1], reverse = sort_AgR)

    diff_dict = {}
    for key in ab_dict:
        if key in a2_dict:
            diff_dict[key] = abs(ab_dict[key] - a2_dict[key])

    sorted_diff_list = sorted(diff_dict.items(), key=lambda item: item[1], reverse = sort_diff)

    ab_number_res_rank_list = [tuple[0] for tuple in sorted_ab_list]
    a2_number_res_rank_list = [tuple[0] for tuple in sorted_a2_list]
    diff_number_res_rank_list = [tuple[0] for tuple in sorted_diff_list]

    # v-- development and debug purposes
    #print(sorted_ab_list)
    #print(ab_number_res_rank_list)
    #print(sorted_a2_list)
    #print(a2_number_res_rank_list)
    #print(sorted_diff_list)
    #print(diff_number_res_rank_list)

    # v-- updated 2/19/2024 to ignore diff rank.
    average_rank = []
    for elem in diff_number_res_rank_list:
        rank_score = 0
        if rankdiff:
            rank_score = sum([ab_number_res_rank_list.index(elem), a2_number_res_rank_list.index(elem), diff_number_res_rank_list.index(elem)]) / 3
        else:
            rank_score = sum([ab_number_res_rank_list.index(elem), a2_number_res_rank_list.index(elem)]) / 2
        average_rank.append((elem, rank_score))
    
    average_rank.sort(key = lambda x: x[1])
    average_ranking_list = [tuple[0] for tuple in average_rank]

    # v-- development and debug purposes
    print("Ranking...")
    print(average_ranking_list)

    return average_ranking_list

'''
Perform Mutation using PyRosetta

VE workflow here.

1. get/make the latest pdb/pose.
    --> in this implementation as of 2/19/2024, 
        --> we rely on PyRosetta mutate_residue, and its packing tool.
        --> perform relax on the new mutated pdb.
        --> reference the new mutated pdb in this workflow. 
    --> Optional: few ways to do (1), could do structure, then docking or MDsim.

2. get new epitope window based on latest pdb.
    --> here we use PyRosetta to get window.
    --> Optional: can use Biopython as well, but its the same. 
        only issue is that we are unable to load other module when loading PyRosetta

3. get the new Ab and ACE2 energy ranking from the latest pdb using the new window.
    --> here we use PyRosetta IAM binding energy.
    --> Optional: can use Interger optimization model, or other different ways to rank mutation.

4. select the top-n ranking mutation.

5. return the list.
'''
def performMutationPyRosetta(to_mutate_info_file_dict, 
                            use_given_window = False, relax_pose = True, allow_repeat = True, MPNN = None, 
                            energy_window = "-100,100", sort_AgAb = False, sort_AgR = False):
                                                      # sort_AgAb = False, sort_AgR = False
    parent_ab_complex_path = to_mutate_info_file_dict['ParentAbComplex']
    ab_chain = to_mutate_info_file_dict['AbChain']
    parent_a2_complex_path = to_mutate_info_file_dict['ParentA2Complex']
    a2_chain = to_mutate_info_file_dict['A2Chain']
    window_str = to_mutate_info_file_dict['Window']
    mutations = to_mutate_info_file_dict['Mutations']
    mutate_info = to_mutate_info_file_dict['MutateInfo']

    # v-- update (2/24/2024) checker for energy window data structure
    lower_E = int(energy_window.split(",")[0])
    higher_E = int(energy_window.split(",")[1])
    if len(energy_window.split(",")) != 2: raise ValueError("\nEnergy window in performMutationPyRosetta() has length != 2.")
    if lower_E >= higher_E: 
        tmp = higher_E
        higher_E = lower_E
        lower_E = tmp

    # v-- updated (2/21/2024) remember past mutation
    past_mutate_list = mutate_info.split(",")
    mutate_info = past_mutate_list[-1]
    for i in range(len(past_mutate_list)):
        if past_mutate_list[i] == "WildType": past_mutate_list[i] = past_mutate_list[i]
        if past_mutate_list[i] != "WildType": past_mutate_list[i] = past_mutate_list[i][:-1]
    #print("Debug for 2/21/2024 update") # <-- development and debug purposes.
    #print(past_mutate_list) # <-- development and debug purposes.
    #print(mutate_info) # <-- development and debug purposes.
    #print(to_mutate_info_file_dict['MutateInfo'].split(",")) # <-- development and debug purposes.

    # check info was passed.
    for checker in [parent_ab_complex_path, ab_chain, parent_a2_complex_path, a2_chain, window_str, mutations, mutate_info]:
        if not checker or checker == "": raise ValueError(f"\nEmpty generation information.")

    '''
    Here we perform mutation to the viral complex based on mutate_info. to get most updated protein

    Using PyRosetta + Docking here to get the new complex for selected mutation.
    '''
    # below is development only.
    new_ab_pdb = os.getcwd() + "/Ab_V.pdb" # <-- perform viral mutation from previous
    new_a2_pdb = os.getcwd() + "/ACE2_V.pdb" # <-- perform viral mutation from previous
    # write new pdb to file.
    if mutate_info == "WildType":
        print(f"0th generation, copying wild type complex. No mutations")
        # v-- update (2/19/2024) use PyRosetta to do the copy.
        copy_ab_pose = pose_from_pdb(parent_ab_complex_path)
        if relax_pose: copy_ab_pose = relaxPose(copy_ab_pose)
        copy_ab_pose.dump_pdb(new_ab_pdb)

        copy_a2_pose = pose_from_pdb(parent_a2_complex_path)
        if relax_pose: copy_a2_pose = relaxPose(copy_a2_pose)
        copy_a2_pose.dump_pdb(new_a2_pdb)

    else:
        print(f"Perform mutation with: {mutate_info}") # <-- all generation other than 0th,
        # parent_ab_complex_path, ab_chain, parent_a2_complex_path
        mainMutationPyRosetta(parent_ab_complex_path, mutate_info, new_ab_pdb, ab_chain, relax_pose = relax_pose)
        mainMutationPyRosetta(parent_a2_complex_path, mutate_info, new_a2_pdb, a2_chain, relax_pose = relax_pose)

    '''
    Here get new epitope window as requested from Ratul
    '''
    window_list = []
    if use_given_window:
        window_list = window_str.split(',')
    else:
        window1 = getEpitopeWindowPyRosetta(new_ab_pdb, ab_chain)
        window2 = getEpitopeWindowPyRosetta(new_a2_pdb, a2_chain)
        tmp = window1 + window2
        tmp = set(tmp)
        window_list = list(tmp)
        window_list.sort()
    #print(window_list) # <-- development and debug
    
    
    # v-- updated (2/21/2024) remember past mutation
    '''
    Here check if we allow same position mutation from past generation. 
    for example, if we dont allow, gen 0th mutate 1A, and not we do not allow mutation at index 1.
    '''
    #print("Debug for 2/21/2024 update") # <-- development and debug purposes.
    #print(window_list) # <-- development and debug purposes.
    if not allow_repeat:
        removed_past = [x for x in window_list if str(x) not in past_mutate_list]
        window_list = removed_past
    #print("New Window:", window_list) # <-- development and debug purposes.
    
    '''
    Try all possible mutation from window. 

    Possible perform Integer Optimization Model or any energy ranking model.
    '''
    # (below is Depreciated, used during development and debug only.)
    #new_mutations_score_dict = {}
    #for residue_position in window_list:
    #    for aa in aa_list:
    #        key = residue_position + aa
    #        '''
    #        Get Binding Energy here. Use PyRosetta.
    #
    #        Currently: below is for development and debug purposes. using random.
    #        '''
    #        new_mutations_score_dict[key] = random.uniform(-100, 100)
    #        # need a way to rank score from antibody and ACE2
    #
    #        '''
    #        Perform Docking here / MDsimulation, to relax and pack the complex.
    #        '''
    #new_mutations_score_tuple_list = sorted(new_mutations_score_dict.items(), key=lambda item: item[1])
    #print("MPNN info:\n", MPNN) # <-- development and debug for 2/23/2024 update

    #print("ACE 2 all mutation score") # <-- development and debug
    ace2_scoring_dict = getAllMutationScore(new_a2_pdb, a2_chain, window_list, MPNN = MPNN, energy_window = [lower_E, higher_E])
    #print(ace2_scoring_dict) # <-- development and debug

    #print("Antibody all mutation score") # <-- development and debug
    anti_scoring_dict = getAllMutationScore(new_ab_pdb, ab_chain, window_list, MPNN = MPNN, energy_window = [lower_E, higher_E])
    #print(anti_scoring_dict) # <-- development and debug

    new_mutations_score_list = getRankingFromScoringDict(anti_scoring_dict, ace2_scoring_dict, sort_AgAb = sort_AgAb, sort_AgR = sort_AgR)
    #print(new_mutations_score_list) # <-- development and debug

    if mutate_info == "WildType": print(f"Get binding energy ranking for: {mutate_info}")
    if mutate_info != "WildType": print(f"Get binding energy ranking for mutant: {mutate_info}")
  
    new_mutations_info_list = [] # <-- contain all info 
    top_mutation_list = [] # <-- for printing purposes
    memory_mutation = set()
    for i in range(len(new_mutations_score_list)):
        # v-- updated (2/21/2024) to remember past mutation.
        if len(memory_mutation) >= int(mutations): break
        curr_mutation = new_mutations_score_list[i] # <-- updated
        if curr_mutation[:-1] in memory_mutation:
            continue
        new_mutation = to_mutate_info_file_dict['MutateInfo'].split(",") + [curr_mutation] # structure: [4A, 8M, ...]

        #print("Debug for 2/21/2024 update") # <-- development and debug purposes.
        #print(new_mutation) # <-- development and debug purposes.

        '''
        Possible Perform Docking here / MDsimulation for chosen Mutation.
        '''
        # Resolved: below is fixed.
        top_mutation_list.append(curr_mutation)
        new_mutations_info_list.append([new_ab_pdb, ab_chain, new_a2_pdb, a2_chain, window_list, mutations, new_mutation])
        
        # v-- updated (2/21/2024) force different position mutation
        memory_mutation.add(curr_mutation[:-1])

    print(f"top-{mutations} diff position binding energy mutation: {top_mutation_list}") # <-- development and debug.
    #print(len(memory_mutation), memory_mutation) # <-- development and debug.
    return new_mutations_info_list


'''
currently using mutate_residue function from PyRosetta library. (2/19/2024)
not sure if this funciton works.

(Depreciated) Not using this anymore, was testing this function 
                as one of the potential mutation that does packing.
                since PyRosetta's mutate_residue does packing, 
                this function is obselete.
'''
def packMutation(pose, posi, amino, scorefxn):
    # Select Mutate Position
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(mut_posi.apply(pose)))

    # Select Neighbor Position
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose)))

    # Select No Design Area
    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)
    #print(pyrosetta.rosetta.core.select.get_residues_from_subset(not_design.apply(pose)))

    # The task factory accepts all the task operations
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    # These are pretty standard
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
    tf.push_back(prevent_subset_repacking)

    # Disable design
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
        pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

    # Enable design
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))
    
    # Create Packer
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
    packer.task_factory(tf)

    #Perform The Move
    if not os.getenv("DEBUG"):
      packer.apply(pose)

'''
update (2/29/2024)
Get all possible mutation from two complex file. 

this script is used mainly in VE2_report.py
'''
def getAllPossibleMutations(generation_info_dict):
    #print(generation_info_dict) # <-- development and debug purposes.
    parent_ab_complex_path = generation_info_dict['ParentAbComplex']
    ab_chain = generation_info_dict['AbChain']
    parent_a2_complex_path = generation_info_dict['ParentA2Complex']
    r_chain = generation_info_dict['A2Chain'] 

    # v-- get the complex pose.
    ab_complex_pose = getPose(parent_ab_complex_path)
    r_complex_pose = getPose(parent_a2_complex_path)

    # v-- get all possible mutation
    ab_all_mutations = getAllMutations(ab_complex_pose, ab_chain)
    r_all_mutations = getAllMutations(r_complex_pose, r_chain)
 
    #print(ab_all_mutations) # <-- development and debug purposes
    #print(r_all_mutations) # <-- development and debug purposes

    all_possible_mutations = ab_all_mutations | r_all_mutations # <-- | operation is union, can also do a.union(b)
    #print(len(all_possible_mutations)) # <-- development and debug. expected: viral.length * 20, in this case, 195 * 20

    return all_possible_mutations

'''
helper script for getAllPossibleMutations(), 
to return all possible mutation of the viral chain given complex and chain info
'''
def getAllMutations(complex_pose, partner_chain):
    viral_chain = partner_chain.split("_")[1]
    viral_chain_list = list(viral_chain)

    all_mutations = set()

    for i in range(1, complex_pose.total_residue() + 1):
        # v-- check to make sure that residue we are looking at is part of viral chain.
        if complex_pose.pdb_info().chain(i) in viral_chain_list:
            for aa in aa_list: # <-- for all possible amino acid
                all_mutations.add(str(complex_pose.pdb_info().number(i)) + aa)
                # v-- below is depreciated
                #if aa != complex_pose.residue(i).name1(): # <-- it is a mutation if it is not the wildtype aa.
                #    all_mutations.add(str(complex_pose.pdb_info().number(i)) + aa)
                #else:
                #    #print(str(complex_pose.pdb_info().number(i)) + complex_pose.residue(i).name1()) # <-- development and debug purposes
                #    continue
    
    return all_mutations

