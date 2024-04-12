# CTRL-V
Computationally tracking of likely viral escape variants by iterative optimization. Official github repo for CTRL-V software.

Official webserver on Google Collab:
- [CTRL-V](https://colab.research.google.com/drive/1YWHpwb-8hn7VPOSCAl4jMeRXOJn1zt8m?usp=sharing)

## Abstract

SARS-CoV-2 emerged in China In late 2019 and has since had a profound global impact, infecting over 670 million individuals and resulting in more than 6.8 million fatalities by 2023. Although vaccines have demonstrated effectiveness against the native SARS-CoV-2 strain and its variants, the Omicron variant has emerged as the predominant variant due to critical mutations in its spike protein. These mutations have led to an increase in binding affinity to human cell receptor ACE2 (angiotensin-converting enzyme 2), allowing for viral entry and antigen replication within the host. 

In response, the human immune system generates B-cells and T-cells to protect against the virus. Opsonization is a crucial immune mechanism for antigen elimination, involving binding antibodies produced by B-cells. The survival of the viral antigen relies on its ability to escape this binding process while maintaining its binding affinity with entry receptors like ACE2, often achieved through mutations in the spike protein region. A priori knowledge of a virus's future infective and transmissive strains will enable humans to be therapeutically prepared with escape-proof antibody formulations. 

To this end, we introduce a simulation platform, CTRL-V (Computational Tracking of Likely Viral Escape variants), that iteratively model the viral escape process of a viral antigenic protein when confronted with human antibodies. We show a test case demonstration in correctly recovering the infective strains of SARS-CoV-2 starting with the wildtype spike receptor-binding domain against known commercial neutralizing antibodies. CTRL-V unveils a putative viral mutational landscape by leveraging only the information about the antibody-antigen complex structure and amino acid interaction preferences in proteins. This information is critical to surveillance and antibody design strategies for preventing viral diseases in humans and livestock.

## Installation

1. Download PyRosetta software as directed on [PyRosetta](https://www.pyrosetta.org/downloads)

2. Clone the repository to directory location,

```
git clone https://github.com/YeeChuen/CTRL-V
```

3. Add pdbs to the pdb files directory.
```
cd ./CTRL-V/pdbs
```

4. Create a config file, following the example and parameters in ```test.ve```

5. Run CTRL-V using
```
python VE2_pyrosetta_workflow.py --config ve_config.ve
```

## Config file parameter.

Checkout example config file in ```./CTRL-V/examples```

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

