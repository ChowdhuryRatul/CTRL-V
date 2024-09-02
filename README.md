![CTRL-V logo](https://github.com/ChowdhuryRatul/CTRL-V/blob/main/CTRL-V_logo_2.png?raw=true)

# CTRL-V: Computational TRacking of Likely Variants

CTRL-V is a viral-escape-inspired, modular, protein biosensor design platform. Current implementation is limited to designing protein biosensors that iteratively evolves to discriminate between two protein targets. The biosensor molecule is called the design molecule (DM), and the preferred binder is called target molecule (TM) and the molecule to avoid is called non-target molecule (non-TM). If you want to use CTRL-V to discriminate between non protein molecules (ions and small molecules), reach out to Prof. Chowdhury (ratul@iastate.edu)

Official github repo for CTRL-V software.
- [Github repo](https://github.com/ChowdhuryRatul/CTRL-V)

Official webserver on Google Collab:
- [CTRL-V](https://colab.research.google.com/drive/1Pkw8MNW4uGnqp5wfW6AIyfM1nORsOyRE?usp=sharing)

Official publication for CTRL-V:
- To be updated. (ProQuest)
- Bioinformatics Journal.

## Abstract

A generalizable computational platform, CTRL-V (Computational TRacking of Likely Variants), is introduced to design selective binding biosensor proteins. Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) has been construed as a lethal, evolving biosensor which has iteratively evolved to distinguish and selectively bind to human entry receptors over neutralizing antibodies. CTRL-V prioritizes mutations that reduce antibody binding while enhancing/ retaining binding with the host entry receptor. CTRL-V categorized 20% (of the 39) reported SARS-CoV-2 point mutations across 30 circulating, infective strains as responsible for immune escape from commercial antibody LY-CoV1404. Specifically, CTRL-V successfully identifies ~70% (five out of seven) single point mutations (371F, 373P, 440K, 445H, 456L) in the latest circulating KP.2 variant and offers detailed structural insights to the escape mechanism. While other viral escape variant predictor tools have been proposed which aim to identify a biochemical objective that justifies the observed genomic transitions in escaping viruses, CTRL-V is a protein-based biosensor design platform at the core – which uses viral escape as an appropriate paradigm to assess its efficacy. We demonstrate three versions of CTRL-V and the usage of integer optimization, stochastic sampling using PyRosetta, and deep learning-based ProteinMPNN for structure-guided biosensor design. Fidelity in capturing preferred amino acid transitions at multiple loci of viral protein demonstrates the utility of CTRL-V as a promising protein-based biosensor design platform.

## Notes

In the [CTRL-V web server](https://colab.research.google.com/drive/1Pkw8MNW4uGnqp5wfW6AIyfM1nORsOyRE?usp=sharing), default values are configured for a CTRL-V job that performs a small viral escape simulation, estimated to take about 20-30 minutes to complete. For any custom CTRL-V jobs or to replicate jobs from the CTRL-V study, it is recommended to clone the [github repo](https://github.com/ChowdhuryRatul/CTRL-V), and run the CTRL-V job on a high-performance computer.

## Installation & Usage

1. Download PyRosetta software as directed on [PyRosetta](https://www.pyrosetta.org/downloads).

2. Clone the repository to directory location,

```
git clone https://github.com/ChowdhuryRatul/CTRL-V
```

3. Add pdbs to the pdb files directory,
```
cd ./CTRL-V/pdbs
```

4. Create a config file, following the example and parameters in ```test.ve```.

5. Run CTRL-V using,
```
python VE2_pyrosetta_workflow.py --config ve_config.ve
```

## Config file parameter.

Checkout more example config file in ```./CTRL-V/examples```.

| Config parameter | description |
| :---:   | :---: |
| JobName | The CTRL-V job name, results from CTRL-V will be saved in a new directory with <JobName>. |
| AbComplex |  The pdb path to antigen-antibody complex file, for example ```./pdbs/LY-CoV1404_clean.pdb```.  |
| AbChain |  The partner chain information for antigen-antibody complex given, use format ```<antibody chain>_<antigen chain>```, for example ```AB_C```.  |
| A2Complex |  The pdb path to antigen-receptor complex file, for example ```./pdbs/6M0J_clean.pdb```.  |
| A2Chain |  The partner chain information for antigen-receptor complex given, use format ```<receptor chain>_<antigen chain>```, for example ```A_E```.  |
| EpitopeWindow |  A list of epitope position for antigen protein in both complex, use format ```<position>,<position>,...```, for example ```444,371,368,446```. Empty input will causes CTRL-V to calculate epitope windows on its own. |
| Mutations |  The number of mutations(single point mutation) for each parent antigen protein, for example ```3``` will allow for 3 mutations per wildtype/variants in each generation. |
| MaxGeneration |  The number of generation to simulate, for example ```3``` will simulate mutation for up to the 3rd generation. |
| useGivenWindow |  Input value: ```[Yes/No]```, Yes: indicate wildtype/variants will always use given ```EpitopeWindow``` value, No: allows CTRL-V to recalculate Epitope Window. |
| relaxPose |  Input value: ```[Yes/No]```, Yes: indicate variants will be relax using PyRosetta, No: disable relax pose.  |
| AllowRepeatMutation |  Input value: ```[Yes/No]```, Yes: enable mutation to the same position, No: disable mutation to the same position  |
| UseMPNN |  Accepts empty value, ```No``` or a list of single point mutation, for example ```444R,444K,444T,444A,444Q,371D,371F,371L,371N,368I,368V,368L,446G,446S,446N```.  |
| EnergyWindow |  Accepts a list of two number, for example ```-100,100```, use ```-100,100``` as default. |
| SortRankAg-Ab |  Input value: ```[Ascending/Descending]```, indicate ranking single point mutations for antigen against antibody in either ascending or descending based on energy value. use ```Ascending``` as default  |
| SortRankAg-R |  Input value: ```[Ascending/Descending]```, indicate ranking single point mutations for antigen against receptor in either ascending or descending based on energy value. use ```Ascending``` as default  |

