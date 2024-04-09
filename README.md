# CTRL-V
Computationally tracking of likely viral escape variants by iterative optimization. Official github repo for CTRL-V software.

## Abstract

SARS-CoV-2 emerged in China In late 2019 and has since had a profound global impact, infecting over 670 million individuals and resulting in more than 6.8 million fatalities by 2023. Although vaccines have demonstrated effectiveness against the native SARS-CoV-2 strain and its variants, the Omicron variant has emerged as the predominant variant due to critical mutations in its spike protein. These mutations have led to an increase in binding affinity to human cell receptor ACE2 (angiotensin-converting enzyme 2), allowing for viral entry and antigen replication within the host. In response, the human immune system generates B-cells and T-cells to protect against the virus. Opsonization is a crucial immune mechanism for antigen elimination, involving binding antibodies produced by B-cells. The survival of the viral antigen relies on its ability to escape this binding process while maintaining its binding affinity with entry receptors like ACE2, often achieved through mutations in the spike protein region. A priori knowledge of a virus's future infective and transmissive strains will enable humans to be therapeutically prepared with escape-proof antibody formulations. To this end, we introduce a simulation platform, CTRL-V (Computational Tracking of Likely Viral Escape variants), that iteratively model the viral escape process of a viral antigenic protein when confronted with human antibodies. We show a test case demonstration in correctly recovering the infective strains of SARS-CoV-2 starting with the wildtype spike receptor-binding domain against known commercial neutralizing antibodies. CTRL-V unveils a putative viral mutational landscape by leveraging only the information about the antibody-antigen complex structure and amino acid interaction preferences in proteins. This information is critical to surveillance and antibody design strategies for preventing viral diseases in humans and livestock.

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

