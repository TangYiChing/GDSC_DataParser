# GDSC_DataParser


Data parser for GDSC datasets

    * filename  source/link
     
    * baseline exoression: Cell_line_RMA_proc_basalExp.txt, source link: https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html
    * baseline expression: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-3983/Downloads?ref=aebrowse
    * mutation: mutations_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
    * copy number variation: cnv_gistic_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
    * cell line: Cell_Lines_Details.xlsx, source link: https://www.cancerrxgene.org/downloads/bulk_download
    * compound: screened_compounds_rel_8.1.csv, source link: https://www.cancerrxgene.org/downloads/bulk_download
    * drug sensitivity: GDSC2_fitted_dose_response_15Oct19.xlsx, source link: https://www.cancerrxgene.org/downloads/bulk_download
    * model list: https://cellmodelpassports.sanger.ac.uk/downloads
    * gene list: gene_identifiers_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
    * cancer gene list: cancer_genes_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads

# How to use it?

Download files listed above and put it in a folder (e.g., /data/DR/db/GDSC/):

```{python}
DB_PATH = '/data/DR/db/GDSC/'
DB_FILE = {'MODEL':'model_list_20200204.csv',
           'EXP_rma':'Cell_line_RMA_proc_basalExp.txt',    # E-MTAB-3610 in ArrayExpress
           'EXP_tpm':'E-MTAB-3983-query-results.tpms.tsv', # E-MTAB-3983 in ArrayExpress
           'CNV_gistic':'cnv_gistic_20191101.csv',              # log2 cnv
           'CNV_abs':'cnv_abs_copy_number_picnic_20191101.csv', # absolute cnv
           'MUT':'mutations_20191101.csv',
           'RESP':'GDSC2_fitted_dose_response_15Oct19.xlsx',
           'CELL': 'Cell_Lines_Details.xlsx',
           'DRUG':'screened_compounds_rel_8.1.csv',
           'GENE':'gene_identifiers_20191101.csv'}

```

type $python useGDSC.py to generate the following files in a tidy format ready for analysis

```{python}
gene expression: cell line by gene matrix
somatic mutation: cell line by gene matrix
copy number variation: cell line by gene matrix
drug sensitivity: compound by cell line matrix
```
