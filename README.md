# IEDB data for IMMREP'25


![IEDB](iedb_immrep/iedb_logo.png)

This repository contains the code used to produce the positive alpha-beta TCR-pMHC pairs visible on Kaggle.

## Instructions
If you'd like to reproduce the analysis, install this repo:

```bash

git clone https://github.com/erichardson97/iedb_immrep25/
cd iedb_immrep25
pip install .
stitchrdl -s human
python -m iedb_immrep --export_file iedb_immrep/dat/export_2_21_25.parquet --output_dir iedb_immrep/dat/results
```

This reproduces the paired dataset. Like everyone else, we have many more unpaired beta chains.
If you'd like those instead, add "--chains beta" to the above.


## Input and output

As input, we selected all TCRs in the IEDB from **Humans** which show recognition of epitopes presented in the context of **MHC Class I**.
We subset for non-modified peptides and assays which have a fully-specified MHC Class I alpha allele (this is the parquet file in iedb_immrep/dat/results).

Because our data has been manually curated from >300 papers with our earliest TCR reference from 1994, our curated TCR and BCR fields can contain values that are difficult to work with for
bioinformaticians as they are author-reported, and may no longer reflect naming conventions. We are standardizing these values in the future so you'll just be able to use the Calculated fields; in this release, we make use of the package
[**tidytcells**](https://github.com/yutanagano/tidytcells) to standardize our values. We also use [**Stitchr**](https://github.com/JamieHeather/stitchr) to produce full-length alpha and beta sequences from our data.

We provide two outputs:

1. A CSV file (immrep_IEDB.csv). This includes the fields: TRAV, TRAJ, TRBV, TRBJ, and each of the resultant CDRs and TRA/TRB sequence output by Stitchr (with leader+constant removed). In addition, we have the field "receptor_ids" which traces back to our unique receptor ID identifiers, "references" which maps back to our reference IRIs, and finally a field "just_10X". This field is "True" where a given TCR-pHLA pair has *only* been seen in 10X experiments, a filter which was used to exclude data in the training of [NetTCR2.2](https://elifesciences.org/articles/93934) demonstrated to improve performance. This is with the exception of [iTRAP-corrected data](https://elifesciences.org/articles/81810). 
2. An AIRR-formatted TSV file; the additional fields are "epitope" and "mhc".

## Other datasets

VDJdb did a special release for the occasion! https://github.com/antigenomics/vdjdb-db/releases/tag/pyvdjdb-2025-02-21.
If you'd like to apply similar processing as we did here, you can run process_vdjdb.py.

## Contact
Please contact us if you'd like to learn more about our data and tools at help@iedb.org.




