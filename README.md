# Histo.fyi Pipeline

A CLI based pipeline for creating the dataset for histo.fyi

## Python libraries used 

To install the required python libraries

`pip install -r requirements.txt`


## Set up localpdb


## Pipeline steps

All commands should be run from this folder.


### Update localpdb

This should be performed before running the pipeline. 

`python steps/pdb_update.py` 

This step will update the localpdb instance to the latest version. It creates a specific version file beforehand so that no error occurs in one specific localpdb step

`python steps/pdb_version.py`

This step simply updates a version file within the 'warehouse' directory

### Initial querying

`python steps/pdb_query.py --mhc_class=class_i`

This step runs a set of queries against the localpdb instance. It is designed to retrieve a wide variety of types of MHC Class I complexes, namely alpha1/2/3 soluble constructs, truncated alpha1/2 constructus, erroneously deposited full length MHC Class I sequences and single chain MHC Class I/beta-2m/peptide constructs.

This generates a collection of pdb_code sets which can be combined and deduplicated to form the basis of the dataset

## Pipeline actions (and order to run them)

`python steps/initialise_records.py --mhc_class=class_i`

This step initialises the core records.

`python steps/fetch_stcrdab_dataset.py --mhc_class=class_i`

This step downloads data on TCRs

`python steps/fetch_sabdab_dataset.py --mhc_class=class_i`

This step downloads data on antibodies

`python steps/find_alike_chains.py --mhc_class=class_i`

This step finds chains similar in sequence to each other to reduce the complexity in chain assignment

`python steps/classify_peptide.py --mhc_class=class_i`

This step identifies the peptide and then compares the peptide sequence in the PDB file with the sequence of the peptide structure to look for gaps and numbering errors. 

TODO: This is quite slow and could do with taking note of previously done work. At some point the tmp files created should be internal sets.

`python steps/assign_chain_types.py --mhc_class=class_i`

This step attempts to assign specific chain types using fuzzy matching

`python steps/assign_species.py --mhc_class=class_i`

This step assigns the species automatically based on the alpha chain. There are sometimes errors in the data and these are fixed by adding a pdb file to the species_overrides constants file.

`python steps/assign_alpha_chain.py --mhc_class=class_i`

This step assigns the alpha chain of the molecule against a collection of canonical sequences











