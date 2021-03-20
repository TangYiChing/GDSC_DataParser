"""
Data parser for GDSC
    filename  source/link
    baseline exoression: Cell_line_RMA_proc_basalExp.txt, source link: https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Home.html
    baseline expression: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-3983/Downloads?ref=aebrowse
    mutation: mutations_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
    copy number variation: cnv_gistic_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
    cell line: Cell_Lines_Details.xlsx, source link: https://www.cancerrxgene.org/downloads/bulk_download
    compound: screened_compounds_rel_8.1.csv, source link: https://www.cancerrxgene.org/downloads/bulk_download
    drug sensitivity: GDSC2_fitted_dose_response_15Oct19.xlsx, source link: https://www.cancerrxgene.org/downloads/bulk_download
    model list: https://cellmodelpassports.sanger.ac.uk/downloads
    gene list: gene_identifiers_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
    cancer gene list: cancer_genes_20191101.csv, source link: https://cellmodelpassports.sanger.ac.uk/downloads
"""


#import built-in pkgs
import os
import sys
import numpy as np
import pandas as pd
import itertools as itrs
#import customized pkgs
import util as ut
import plot_util as pt
import externaldb as exdb

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

class UseGDSC:
    # initialize
    def __init__(self, dbPathStr=DB_PATH, dbFileDict=DB_FILE):
        """
        Data parser for GDSC dataset: https://www.cancerrxgene.org/downloads/anova
        
        :param dbPathStr: string representing path to the database
        :param dbFileDict: dict representing file location to each data        
        """
        self.folder_str = dbPathStr
        self.file_dict = dbFileDict

    # parse raw data
    def parseDRUG(self):
        """
        Read file: screened_compounds_rel_8.1.csv
        
        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['DRUG']
        df = pd.read_csv(f_str, header=0, sep=",")
        return df

    def parseRESP(self):
        """
        Read file: GDSC2_fitted_dose_response_15Oct19.xlsx

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['RESP']
        df = pd.read_excel(f_str)
        return df

    def parseGENE(self):
        """
        Read file: gene_identifiers_20191101.csv

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['GENE']
        df = pd.read_csv(f_str, header=0, sep=",")
        return df

    def parseCELL(self):
        """
        Read file: Cell_Lines_Details.xlsx

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['CELL']
        df = pd.read_excel(f_str)
        # clean data
        df = df.iloc[0:-1,:] # remove the last line
        return df

    def parseCNV(self, use='gistic'):
        """
        Read file: 
            cnv_gistic_20191101.csv                   # log2 cnv
            cnv_abs_copy_number_picnic_20191101.csv   # absolute cnv

        :param use: string representing data type, options=[ gistic| abs ], default:gistic
        :return df: dataframe with headers
        """
        if use == 'gistic':
            f_str = self.folder_str + self.file_dict['CNV_gistic']
            df = pd.read_csv(f_str, header=0, sep=',', low_memory=False)
            df = df.drop('Unnamed: 1', axis=1).drop(df.index[[0,1]])
        elif use == 'abs':
            f_str = self.folder_str + self.file_dict['CNV_abs']
            df = pd.read_csv(f_str, header=0, sep=',', low_memory=False)
            df = df.drop('Unnamed: 1', axis=1).drop(df.index[[0,1]])
        else:
            print( 'ERROR: expecting [abs | gistic], got {:}'.format(use) )
            sys.exit(1)
        df.set_index('model_id', inplace=True)
        df = df.T # model_id*gene_id
        return df

    def parseMUT(self):
        """
        Read file: mutations_20191101.cs

        :return df: dataframe with headers
        """
        f_str = self.folder_str + self.file_dict['MUT']
        df = pd.read_csv(f_str, header=0, sep=',')
        return df

    def parseEXP(self, use='TPM'):
        """
        Read file:
           Cell_line_RMA_proc_basalExp.txt    # E-MTAB-3610 in ArrayExpress
           E-MTAB-3983-query-results.tpms.tsv # E-MTAB-3983 in ArrayExpress

        :param use: string representing data type, options=[ TPM | RMA ], default:TPM
        :return df: dataframe with headers
        """
        if use == 'TPM': 
            f_str = self.folder_str + self.file_dict['EXP_tpm']
            df = pd.read_csv(f_str, header=0, index_col=1, sep="\t", skiprows=4)

        elif use == 'RMA':
            f_str = self.folder_str + self.file_dict['EXP_rma']
            df = pd.read_csv(f_str, header=0, index_col=0, sep="\t")

        else:
            print( 'ERROR: expecting [TPM | RMA], got {:}'.format(use) )
            sys.exit(1)
        return df

    # retrieve processed data
    def getMODEL(self):
        """
        Retrieve model information.

        :return df: dataframe with headers

        Note:
        =====
        Return dataframe was merged from two files:
            1. CCLE's sample_info.csv
            2. Cell Model Passport's model_list_20200204.csv
        This method calls the retrieveMergedMODEL() method of externaldb.py
        """
        df = exdb.retrieveMergedMODEL(db='GDSC')
        return df

    def getDRUG(self, use=['target', 'smile', 'resp']):
        """
        Retrieve drug that has data in given data types

        :param use: list containing data types, options=[target | smile], default: [target | smile]
        :return inter_drugname_list: list containing drug name that has data in given data types
        """
        # load data
        target = self.getCompoundTARGET(use_exdb=True)
        smile = self.getCompoundSMILE(use_exdb=True)
        resp = self.parseRESP()
        drugname_dict = { 'target': list(target['drug'].unique()), 'smile': list(smile['drug'].unique()), 
                          'resp': list(resp['DRUG_NAME'].unique()) }
        # get intersection drug names
        use_drugname_list = []
        for data_str in use:
            use_drugname_list.append( set(drugname_dict[data_str]) )
            #print(data_str, len(drugname_dict[data_str]))
        inter_drugname_list = sorted(list( set.intersection( *use_drugname_list ) ))
        #print( 'data used={:}, shared drugs={:}'.format(use, len(inter_drugname_list)) )
        return inter_drugname_list

    def getCELL(self, use=['exp', 'cnv', 'mut', 'resp']):
        """
        Retrieve cell that has data in given data types

        :param use: list containing data types, options=[exp | cnv | mut | resp], default: [exp, cnv, mut, resp]
        :return inter_modelID_list: list containing cell modelIDs that has data in given data types
        """
        # load data
        exp = self.getEXP(use='TPM')
        cnv = self.getCNV(use='gistic')
        mut = self.parseMUT()
        resp = self.parseRESP()
        modelID_dict = { 'exp': exp.index.tolist(), 'cnv': cnv.index.tolist(), 
                         'mut': list(mut['model_id'].unique()), 'resp': list(resp['SANGER_MODEL_ID'].unique()) }
        for key, value in modelID_dict.items():
            print('data={:}, #cells={:}, {:}'.format(key, len(value), value[:5]))
        # get intersection cell lines
        use_modelID_list = []
        for data_str in use:
            use_modelID_list.append( set(modelID_dict[data_str]) )
        inter_modelID_list = sorted(list( set.intersection( *use_modelID_list ) ))
        print( 'data used={:}, shared cells={:}'.format(use, len(inter_modelID_list)) )
        return inter_modelID_list

    def getCELLMap(self, cellList, colStr):
        """
        Retrieve cell mapping for the given cells

        :param cellList: list containing a list of cells (e.g., model_id, model_name)
        :param colStr: string representing column name of the model where itemList can be retrieved. 
                       options=[model_id, model_name, BROAD_ID, RRID, COSMIC_ID]
        :return df: dataframe contains headers 

        Note:
        =====
        This method calls the retrieveMergedMODEL() method of externaldb.py
        """
        # load data
        model = self.getMODEL(db='GDSC')
        # check if colStr in model
        if not colStr in model.columns:
            print( 'ERROR: {:} not in model, options={:}'.format(colStr, '[model_id, model_name, BROAD_ID, RRID, COSMIC_ID]') )
            sys.exit(1)
        # return
        found_in_model = list(set(cellList) & set(model[colStr]))
        if len(found_in_model) != len(cellList):
            not_found_in_model = list(set(cellList) - set(model[colStr]))
            print( 'WARNING: Not Found={:},\n    {:}'.format(len(not_found_in_model), not_found_in_model) )
            df = model.loc[model[colStr].isin(cellList)]
        else:
            df = model.loc[model[colStr].isin(cellList)]
        return df

    def getMUT(self, cellList=None, to_dict=False):
        """
        Retrieve cell:mutation dictionary for the given cells

        :param cellList: list containing a list of cells (i.e., model_id)
        :param to_dict: boolean indicating to output dictionary or dataframe
        :return cell_mutList_dict: dictionary containing cell: mutation_list pairs
        
        Note:
        =====
        if cellList=None, the program will return all cells available in parseMUT()
        """
        # load data
        df = self.parseMUT()
        # clean data
        use_cols = ['model_id', 'gene_symbol']
        df = df[use_cols].copy()
        df = df.drop_duplicates(keep='first')
        # get cellList
        if cellList == None:
            cellList = df['model_id'].values.tolist()
        else:
            cellList = sorted(list( set(cellList) & set(df['model_id']) ))
        # subsetting to include only cellList
        use_df = df.loc[df['model_id'].isin(cellList)]
        # create dict
        if to_dict == True:
            cell_mutList_dict = { cid: gnm['gene_symbol'].values.tolist() for cid, gnm in use_df.groupby('model_id') }
            outputs = cell_mutList_dict
        else:
            outputs = use_df
        # return
        return outputs

    def getCNV(self, use='gistic', cellList=None, categorical=True):
        """
        Retrieve copy number for the given cells

        :param cellList: list containing a list of cells (i.e., model_id)
        :param categorical: boolean indicating to whether to return gistic score or one-hot CNV
        :return df: dataframe with model_id by Gene

        Note:
        =====
        if cellList=None, the program will return all cells available in parseCNV()
        if categorical is True, gistic score != 0 will set to 1, otherwise remains 0       

        Referebce:
        ==========
        the GISTIC (Genomic Identification of Significant Targets in Cancer) calls comprising:
        -2 (deletion), -1 (loss), 0 (diploid), 1 (gain), and 2 (amplification) were made using GISTIC2.0 [20].
        https://www.thno.org/v10p3939.htm
        """
        # load data
        if use == 'gistic':
            raw_cnv = self.parseCNV(use='gistic')
        elif use == 'abs':
            raw_cnv = self.parseCNV(use='abs')
        else:
            print( 'ERROR: expecting [abs | gistic], got {:}'.format(use) )
            sys.exit(1)
        # select by cell ids
        if cellList == None:
            cellList = raw_cnv.index.tolist()
        else:
            cellList = sorted(list( set(cellList) & set(raw_cnv.index.tolist()) ))
        # subsetting
        df = raw_cnv.loc[cellList]
        # convert string to float
        df = df.astype(float, errors='ignore')
        # categorical arg only applys to gistic score
        if use == 'gistic' and categorical == True:
            df[(df!=0)] = 1     # loose criteria
        return df

    def getCompoundSMILE(self, use_exdb=True, drugList=None):
        """
        Retrieve drug chemical structure: SMILE

        :param drugList: list containing  a list of drug name
        :param use_exdb: boolean indicating whether use external database or not
        :return df: dataframe contains headers=[drug, smile]

        Note:
        =====
        This method calls the retrieveSMILE() method of externaldb.py
        If drugList is None, this program will use drug name available
        """
        # load data
        resp = self.parseRESP()
        # replace name
        new_name_dict = { 'Cisplatin': 'cis-Platin', 'Nutlin-3a (-)': 'nutlin-3A',
                          'Mirin': 'MRN-ATM Pathway Inhibitor, Mirin', 'Oxaliplatin': 'Oxalitin',
                          'Picolinici-acid': 'Picolinate', 'KRAS (G12C) Inhibitor-12': 'K-Ras(G12C) inhibitor 12',
                          'GDC0810': 'Brilanestrant', 'BPD-00008900': 'Bpd-MA'} # oldName: newName
        resp['DRUG_NAME'].replace(to_replace=new_name_dict, inplace=True)
        drug_list = list( resp['DRUG_NAME'].unique() )
        # use external database
        drug_df = exdb.retrieveSMILE(drug_list, db='PubChem')
        # change drug name back to original ones
        old_name_dict = { value:key for key,value in new_name_dict.items() }
        drug_df['drug'].replace(to_replace=old_name_dict, inplace=True)
        if drug_df.shape[0] < len(drug_list):
            print( '{:}/{:} have SMILE string'.format(len(drug_df), len(drug_list)) )
        # select by drug name
        if drugList == None:
            drugList = drug_df['drug'].values.tolist()
        df = drug_df.loc[drug_df['drug'].isin(drugList)]
        # return 
        if df.shape[0] < len(drugList):
            print( 'these drugs do not have SMILE string: {:}'.format(list(set(drugList)-set(df['drug']))) )
        return df #headers=[drug, smile]
        
    def getCompoundTARGET(self, use_exdb=True, drugList=None):
        """
        Retrieve drug target gene/protein 
  
        :param drugList: list containing  a list of drug name
        :param use_exdb: boolean indicating whether use external database or not
        :return df: dataframe contains headers=[drug, target]

        Note:
        =====
        This method calls the retrieveTARGET() method of externaldb.py
        If drugList is None, this program will use drug name available
        """
        # load data
        df = self.parseDRUG()
        use_cols = ['DRUG_NAME',  'PUTATIVE_TARGET']
        df = df[use_cols].replace( to_replace={np.nan:'not defined'} )
        target_df = ut.wide2long(df, 'DRUG_NAME', 'PUTATIVE_TARGET', ',')
        value_dict = { val: val.strip() for val in target_df['PUTATIVE_TARGET'].values.tolist() }
        target_df['PUTATIVE_TARGET'] = target_df['PUTATIVE_TARGET'].replace(to_replace=value_dict)# remove space 
        ######target_df = target_df.loc[target_df['PUTATIVE_TARGET']!='not defined']# clean data by removing na
        target_df.columns = ['drug', 'target']
        # check condition
        if use_exdb:
            # retrieve targets from external database
            map_dict = {'Nelarabine':'NELARABINE', 'TPCA-1':'IKK-2 INHIBITOR IV', 'Dactinomycin':'DACTINOMYCIN', 'SN-38':'7-ETHYL-10-HYDROXY-CAMPTOTHECIN',
                        'CZC24832':'CHEMBL2064571', 'JNK Inhibitor VIII':'CHEMBL210618', 'Zibotentan':'ZIBOTENTAN', 'Linifanib':'LINIFANIB',
                        'Dactolisib':'DACTOLISIB', 'Tozasertib':'TOZASERTIB', 'Avagacestat':'AVAGACESTAT', 'Leflunomide':'LEFLUNOMIDE',
                        'Masitinib':'MASITINIB', 'Voxtalisib':'VOXTALISIB', 'Uprosertib':'UPROSERTIB', 'Apitolisib':'APITOLISIB', 'Pevonedistat':'PEVONEDISTAT',
                        'Alisertib':'ALISERTIB', 'Entospletinib':'ENTOSPLETINIB', 'Tozasertib':'TOZASERTIB', 'Brivanib, BMS-540215':'BRIVANIB',
                        'Obatoclax Mesylate':'OBATOCLAX MESYLATE', 'Sepantronium bromide':'SEPANTRONIUM BROMIDE', 'Pictilisib':'PICTILISIB',
                        'AZD8186':'AZD-8186', 'PRT062607':'CHEMBL2177736', 'MK-8776':'SCH-900776', 'EPZ5676':'CHEMBL3087499', 'I-BET-762':'CHEMBL1232461',
                        'AZD5363':'AZD-5363',  'ZM447439':'CHEMBL202721', '5-Fluorouracil':'Fluorouracil',
                        'Doramapimod':'DORAMAPIMOD', 'Luminespib':'NVP-AUY922', '(5Z)-7-Oxozeaenol':'5Z-7-OXOZEAENOL'}
            target_df['drug'] = target_df['drug'].replace(to_replace=map_dict) #modify drug name to get more matches
            ex_target_df = exdb.retrieveTARGET(target_df['drug'].values.tolist(), db='DGIdb')
            # return
            drug_df = pd.concat([target_df, ex_target_df], axis=0, sort='True') # merge
            map2_dict = { value:key for key, value in map_dict.items() } #change back to original drug name
            drug_df['drug'] = drug_df['drug'].replace(to_replace=map2_dict) # replace name
            drug_df.drop_duplicates(inplace=True) # remove duplicated drug-target pair
            # select by gene name
            if drugList == None:
                drugList = drug_df['drug'].values.tolist()
            df = drug_df.loc[drug_df['drug'].isin(drugList)]
        else:
            # return 
            drug_df = target_df
            drug_df.drop_duplicates(inplace=True) # remove duplicated drug-target pair
            # select by gene name
            if drugList == None:
                drugList = drug_df['drug'].values.tolist()
            df = drug_df.loc[drug_df['drug'].isin(drugList)]

        df = df.dropna(axis=0)
        return df # headers = ['drug', 'target']

    def _replace_cell_name(self, use=['TPM', 'RMA']):
        """
        Return model_id for cell name in given data
     
        :param use: list containing data type
        :return cid_dict: dictionary containing model_id for each given data
        """
        # cheating set
        tpm_dict = { 'DIFI': 'SIDM00049', 'GEO': 'SIDM00068', '143b':'SIDM00804',
                     'COLO-320-HSR':'SIDM00842', 'H2373':'SIDM00103', 'H2461':'SIDM00102',
                     'H2591':'SIDM00101', 'H2595':'SIDM00100', 'H2722':'SIDM00099', 'H513':'SIDM00114',
                     'H2731':'SIDM00098', 'H2795':'SIDM00154', 'H2803':'SIDM00309',
                     'H2804':'SIDM00310', 'H2810':'SIDM00311', 'H2818':'SIDM00520',
                     'H2869':'SIDM00519', 'Hs633T':'SIDM00667', 'KMH-2':'SIDM00619',
                     'KO52':'SIDM00018', 'MMAC-SF':'SIDM01242', 'NB-TU-1-10':'SIDM00579',
                     'NCI-SNU-1':'SIDM01146', 'NCI-SNU-16':'SIDM01145', 'NTERA-2cl-D1':'SIDM01203',
                     'UWB1-289':'SIDM00815', 'U031':'SIDM00112', 'SR':'SIDM00096' }

        rma_dict = { '1240139': 'SIDM00518', '1290808': 'SIDM00104', '907796':'SIDM01219', '910569':'SIDM00842',
                     '1290907': 'SIDM00597', '1240140': 'SIDM00517', '1503362.1':'SIDM00445',
                     '906815':'SIDM01483', '907284':'SIDM00018', '1247873':'SIDM00046', '1331031':'SIDM00400',
                     '1330983':'SIDM00461', '1240156':'SIDM00037', '1479987':'SIDM00969'}

        # program start
        cid_dict = {'TPM':{}, 'RMA':{}}
        # load data
        tpm = self.parseEXP(use='TPM')
        rma = self.parseEXP(use='RMA')
        model = self.getMODEL() # model_id, model_name
        ## Cell name in TPM is a string containing (cname, cancer_type, tissue)
        ## use cname to find corresponding model_id 
        not_found_cname_list = []
        for cname_str in tpm.columns.tolist()[1:]:
            cname, cancer_type, tissue = cname_str.split(',')
            if cname in model['model_name'].values.tolist():
                cid_dict['TPM'].update( {cname_str: model.loc[model['model_name']==cname]['model_id'].values[0]} )
            else:
                if cname in tpm_dict.keys():
                    cid_dict['TPM'].update( {cname_str: tpm_dict[cname]} )
                else:
                    not_found_cname_list.append(cname)
        # send warning message
        if len(not_found_cname_list) > 0:
            print( 'WARNING: cells in EXP_TPM not_found model_id={:}'.format(not_found_cname_list) )
            
        ## Cell name in RMA is a string containg DATA.COSMIC_ID
        ## use COSMIC_ID to find corresponding model_id
        not_found_cname_list = []
        for cname_str in rma.columns.tolist()[1:]:
            cosmic = cname_str.split('DATA.')[1]
            if cosmic in model['COSMIC_ID'].values.tolist():
                cid_dict['RMA'].update( {cname_str: model.loc[model['COSMIC_ID']==cosmic]['model_id'].values[0]} )
            else:
                if cosmic in rma_dict.keys():
                    cid_dict['RMA'].update( {cname_str: rma_dict[cosmic]} )
                else:
                    not_found_cname_list.append(cosmic)
        # send warning message
        if len(not_found_cname_list) > 0:
            print( 'WARNING: cells in EXP_RMA not_found model_id={:}'.format(not_found_cname_list) )
    
        # return
        return cid_dict #keys=[TPM, RMA]


    def getEXP(self, use='TPM', cellList=None):
        """
        Retrieve processed expression data
     
        :param cellList: list containing a list of cells (i.e., DepMap_ID)
        :param use: string representing data type to be retrieved, options=[TPM | RMA]
        :return mat: matrix of model_id by gene

        Note:
        =====
        if use='TPM', then perform log2(TPM+1) for the expression values
        Due to name:id replacement, output from getEXP() may have fewer cells than the output from parseEXP()
        Due to missing gene name, output from getEXP() may have fewer genes than the output from parseEXP()
        """
        # load data dictionary to replace cname with cid
        cid_dict = self._replace_cell_name(use=['TPM', 'RMA'])

        # program start
        if use == 'TPM':
            # load data
            raw_exp = self.parseEXP(use='TPM')
            # clean data
            use_cols = raw_exp.columns.tolist()[1:] # remove col: Gene ID
            cname_list = sorted(list( set(raw_exp[use_cols].columns) & set(cid_dict['TPM'].keys()) )) # retrieve cnames that have cid
            exp = raw_exp[use_cols][cname_list] # select columns by cnames that have cid
            exp.rename(columns=cid_dict['TPM'], inplace=True) # replace cname with cid
            # convert
            exp = exp.T # conver to cell*gene
            exp = np.log2( exp+1 )  # Log2 transformed, using a pseudo-count of 1.

        elif use == 'RMA':
            # load data
            raw_exp = self.parseEXP(use='RMA')
            # clean data
            use_cols = raw_exp.columns.tolist()[1:] # remove col: GENE_title
            cname_list = sorted(list( set(raw_exp[use_cols].columns) & set(cid_dict['RMA'].keys()) )) # retrieve cnames that have cid
            exp = raw_exp[use_cols][cname_list] # select columns by cnames that have cid
            exp.rename(columns=cid_dict['RMA'], inplace=True) # replace cname with cid
            # convert
            exp = exp.T # conver to cell*gene

        else:
            print( 'ERROR: expecting [TPM | RMA], got {:}'.format(use) )
            sys.exit(1)
        # return
        exp.columns.name = 'gene'
        exp.index.name = 'model_id'
        raw_exp = exp.loc[:, exp.columns.notnull()] # remove gene is na
        # select by cell ids
        if cellList == None:
            cellList = raw_exp.index.tolist()
        else:
            cellList = sorted(list( set(cellList) & set(raw_exp.index.tolist()) ))
        # return
        df = raw_exp.loc[cellList]
        return df

    def getCIDMAPDICT(self, keyStr, valueStr):
        # load data
        model = self.getMODEL()
        # check columns
        for col in [keyStr, valueStr]:
            if not col in model.columns:
                print( 'ERROR: {:} not in model'.format(col) )
                sys.exit(1)
        # get mapping dict
        map_dict = dict( zip(model[keyStr].values.tolist(), model[valueStr].values.tolist()) )
        return map_dict

    def getGIDMAPDICT(self, keyStr, valueStr):
        gene = self.parseGENE()
        # get mapping dict
        map_dict = dict( zip(gene[keyStr].values.tolist(), gene[valueStr].values.tolist()) )
        return map_dict

    def getContRESPDICT(self, use='IC50', drugList=None, cellList=None):
        """
        Retrieve neumeric drug response

        :param use: options = [IC50 | LN_IC50 | Z_SCORE]
        :param drugList: list containing  a list of drug name
        :param cellList: list containing a list of cells (e.g., model_id, model_name)
        :return respMat_dict: dictionary containing respMat of drug by mode_id
        """
        # setting
        resp = self.parseRESP()
        use_dict = { 'IC50': 'IC50', 'LN_IC50':'LN_IC50', 'Z_SCORE': 'Z_SCORE', 'AUC':'AUC'}
        use_cols =  ['SANGER_MODEL_ID', 'DRUG_NAME'] + [ use_dict[use] ]
        print(use_cols)
        if use == 'IC50':
            resp['IC50'] = np.exp(resp['LN_IC50'] ) # convert LN_IC50 back to IC50
        # indexing
        df1 = resp.set_index(['DRUG_NAME']).sort_index()
        df2 = resp.set_index(['DRUG_NAME', 'MIN_CONC', 'MAX_CONC']).sort_index()
        # collect wanted drugs by checking replication
        EXP1_list, EXP2_list = self.checkRESP(use=use)
        # retrieve data and check duplication by averaging them
        exp_list = [EXP1_list, EXP2_list]
        exp_df_dict = {}
        for i in range(len(exp_list)):
            exp_int = i+1
            if len(exp_list[i]) != 0:
                # load data
                exp_df = df2.loc[exp_list[i]].reset_index()[use_cols]
                exp_df = exp_df.groupby(['DRUG_NAME', 'SANGER_MODEL_ID']).mean().reset_index()
                exp_df_dict.update( {'exp'+str(exp_int):exp_df} )
        # conver to respMat
        respMat_dict = {}
        for key, df in exp_df_dict.items():
            # select by drug name
            if drugList != None:
                #drugList = df.index.tolist()
                df = df.loc[df['DRUG_NAME'].isin(drugList)]
            # select by cell name
            if cellList != None:
                #cellList = df.columns.tolist()
                df = df.loc[df['SANGER_MODEL_ID'].isin(cellList)]

            # long to wide
            respMat = df.pivot(index='DRUG_NAME', columns='SANGER_MODEL_ID', values=use)
            respMat.fillna(np.nan, inplace=True)
            # save to dict
            respMat.index.name = 'drug'
            respMat.columns.name = 'model_id'
            respMat_dict.update( {key:respMat} )
        return respMat_dict # key:value = exp1:respMat (i.e., drug-cell response matrix)

    def _getDeltaResponse(self, respMat, by='cell'): # deprecated due to speed, use respMat2respDf methond.
        """
        Return delta response of pairs

        :param respMat: matrix with Compound by DepMap_ID
        :param by: string representing key of return delta_dict, options=['cell', 'drug'], default=cell
        :return delta_dict: dictionary containing delta response of pairs

        Note:
        =====
        if by='cell', the program will calculate delta response of two drugs for each cell
        if by='drug', the program will calculate delta response of two cells for each drug
        """
        # create result dict
        delta_dict = {}
        # condition
        if by == 'cell':
            respMat = respMat.T
        elif by == 'drug':
            respMat = respMat
        else:
            print( 'ERROR: by={:}, options=[cell|drug]'.format(by) )

        # all-pair list
        pair_list = [ pair for pair in itrs.combinations(respMat.columns, 2) ]
        # looping rows and calculate delta response
        for idx in respMat.index: # if by=cell, idx will be cell id, if by=drug, idx will be drug name
            # load data
            idx_df = respMat.loc[idx] # Series
            idx_df = idx_df.to_frame(name='resp') # dataframe
            idx_df.dropna(axis=0, inplace=True)
            # calculate pairwise difference among possible column pair (i.e., cellpair, drugpair)
            delta_difference_df = pd.DataFrame(np.abs( idx_df['resp'].values - idx_df['resp'].values[:, None] ),
                                               columns = idx_df.index.tolist(), index = idx_df.index.tolist())
            # convert to long-form df
            delta_df = ut.sym2half(delta_difference_df, keepSelfPairs=False) # columns = [idxpair, similarity]
            # change column name
            if by == 'cell':
                delta_df.columns = ['drugpair', 'delta response']
            else:
                delta_df.columns = ['cellpair', 'delta response']
            # update to result
            delta_dict.update( {idx: delta_df} )
        # return
        #print(delta_dict)
        return delta_dict

    def respMat2respDf(self, respMat):
        """
        Return delta response of cell-pairs again certain drug

        :param respMat: matrix with Compound by DepMap_ID
        :return delta_dict: dictionary containing delta response of cell-pairs for each compund in respMat

        Note:
        =====
        1. delta response of cell-pairs is defined by response value of one cell - response value of the other cell
                for example: ACH-000001 has response value of 0.245932 against Compound 17-AAG, and
                             ACH-000005 has response value of 1.616860 against Compound 17-AAG, therefore
                             delta response of ACH-000001-ACH-000005 pair is 0.245932 - 1.616860 = -1.370928
        2. If any of cell-pair has one na, the cell pair will not be saved into result
        """
        # create result dict
        delta_dict = {}
        # create cell-pair list
        cellpair_list = [ pair for pair in itrs.combinations(respMat.columns, 2) ]
        # loop through drug and calculate delta response
        for drug in respMat.index:
            # load data
            drug_df = respMat.loc[drug] # Series
            drug_df = drug_df.to_frame(name='resp') # dataframe
            drug_df.dropna(axis=0, inplace=True)
            # calculate pairwise difference among possible column pair (i.e., cellpair)
            delta_difference_df = pd.DataFrame(np.abs( drug_df['resp'].values - drug_df['resp'].values[:, None] ),
                                               columns = drug_df.index.tolist(), index = drug_df.index.tolist())
            # convert to long-form df
            delta_df = ut.sym2half(delta_difference_df, keepSelfPairs=False) # columns = [cellpair, similarity]
            delta_df.columns = ['cellpair', 'delta response']
            # update to result
            delta_dict.update( {drug: delta_df} )
        # return
        return delta_dict

    def checkRESP(self, use='IC50'):
        """
        :param use: options = [IC50 | LN_IC50 | Zscore | AUC]
        """
        # setting
        resp = self.parseRESP()
        use_dict = { 'IC50': 'IC50', 'LN_IC50':'LN_IC50', 'Z_SCORE': 'Z_SCORE', 'AUC':'AUC'}
        use_cols =  ['SANGER_MODEL_ID', 'DRUG_NAME'] + [ use_dict[use] ]
        #print(use_cols) 
        if use == 'IC50':
            resp['IC50'] = np.exp (resp['LN_IC50'] ) # convert LN_IC50 back to IC50
        # indexing
        df1 = resp.set_index(['DRUG_NAME']).sort_index()
        df2 = resp.set_index(['DRUG_NAME', 'MIN_CONC', 'MAX_CONC']).sort_index()
        # orgaining data
        data_list = [] # list of dict for each dose range: [{'DRUG_NAME':, 'MIN_CONC':, 'MAX_CONC':, 'CELL':[]}, {'', '', '', ''}....]
        drug_idx_dict = { drug:[] for drug in set(df1.index) } # {drug:[0, 1]} index of data list
        for idx in set(df2.index):
            cell = df2.loc[idx, ['SANGER_MODEL_ID',use_dict[use]]]
            drug_dict = {'DRUG_NAME': idx[0], 'MIN_CONC': idx[1], 'MAX_CONC': idx[2],
                         'CELL': df2.loc[idx, ['SANGER_MODEL_ID',use_dict[use]]]}
            data_list.append( drug_dict )
        for i in range(len(data_list)):
            drug = data_list[i]['DRUG_NAME']
            drug_idx_dict[drug].append(i)
        # Branching based on dose range
        BASE_list = [] # If a drug has no multiple dose range
        EXP1_list = [] # If a drug has multiple dose range had different max conc, 
        EXP2_list = [] # then put into 2 lists: EXP1, EXP2
        # get some stats also
        n_drug_has_multiple_dose_range = 0
        n_drug_has_diff_max_dose = 0
        # looping to inspect drug by drug
        for drug, idxList in drug_idx_dict.items():
           
            if len(idxList) > 1: # has multiple dose range
                n_drug_has_multiple_dose_range+=1
                drug_exp_dict = {drug: {'maxc':None, 'exp1':[], 'exp2':[]} } # branching
                for idx in idxList: 
                    drug_dict = data_list[idx]
                    drug_df = drug_dict['CELL']
                    maxc = list(set(drug_df.index))[0][2] 
                    if drug_exp_dict[drug]['maxc'] == None:
                        drug_exp_dict[drug]['maxc'] = maxc
                        EXP1_list.append(list(set(drug_df.index))[0])
                    else:
                        if  maxc == drug_exp_dict[drug]['maxc'] and not list(set(drug_df.index))[0] in EXP1_list:
                            EXP1_list.append(list(set(drug_df.index))[0])
                            drug_exp_dict[drug]['exp1'].append(list(set(drug_df.index))[0])
                        else:
                            EXP2_list.append(list(set(drug_df.index))[0])
                            drug_exp_dict[drug]['exp2'].append(list(set(drug_df.index))[0])
                # summary
                if len(drug_exp_dict[drug]['exp2']) > 0:
                    n_drug_has_diff_max_dose+=1

            else:
                drug_dict = data_list[idxList[0]]
                drug_df = drug_dict['CELL']
                BASE_list.append( list(set(drug_df.index))[0] ) # e.g., ('Alpelisib', 0.005003, 10.0)
                
        #print( 'total unique drugs={:}'.format(len(drug_idx_dict)) )
        #print( '    {:} has 2+ dose range'.format(n_drug_has_multiple_dose_range) )
        #print( '        {:} has different max dose'.format(n_drug_has_diff_max_dose) )
        print( 'BASE_list={:}, EXP1_list={:}, EXP2_list={:}'.format(len(BASE_list), len(EXP1_list), len(EXP2_list)) )
        #print( 'EXP1={:}\nEXP2={:}'.format(EXP1_list, EXP2_list) )
        return list(set(BASE_list+EXP1_list)), list(set(BASE_list+EXP2_list))

    def queryRESP(self, drugList, cellList):
        """
        Return dataframe containing drug sensitivity data for drugList and cellList

        :param drugList: a list contatining drug names
        :param cellList: a list containing cell id
        :return df: dataframe containing drug sensitivity data for drugList and cellList
        """
        # load data
        resp = self.parseRESP()
        # subsetting resp
        df = resp[ resp['DRUG_NAME'].isin(drugList) & resp['SANGER_MODEL_ID'].isin(cellList) ]
        # convert back to IC50
        df['IC50'] = np.exp( df['LN_IC50'] )
        # return
        use_cols = ['DRUG_NAME', 'SANGER_MODEL_ID', 'MIN_CONC', 'MAX_CONC', 'IC50', 'LN_IC50', 'AUC', 'RMSE', 'Z_SCORE']
        df = df[use_cols]
        print(df)
        return df

# execute from script for debugging
if __name__ == "__main__":
    #  TESTING
    gdsc = UseGDSC() # initiate an instance
    print(__doc__)
    # set choices
    test = False
    save = False

    if test:
        # TESTING
        #save raw data: exp, mut, cnv, target, resp
        #gene = gdsc.parseGENE()
        #gid_gnm_dict = dict(zip(gene['gene_id'], gene['cosmic_gene_symbol']))
        #exp = gdsc.getEXP(use='TPM', cellList=None)
        #print(exp)
        #exp = exp.dropna(axis=1, how='all')
        #trans_exp = exp.fillna(exp.mean())
        #print(trans_exp)
        #trans_exp.to_csv('./GDSC.RAW.EXP.TPM.Mat.txt',header=True, index=True, sep="\t")

        #exp = gdsc.getEXP(use='RMA', cellList=None)
        #print(exp.isnull().sum().sum())
        #exp.to_csv('./GDSC.RAW.EXP.RMA.Mat.txt',header=True, index=True, sep="\t")
        #cnv_c = gdsc.getCNV(use='gistic', cellList=None, categorical=True)
        #mut = gdsc.getMUT(cellList=None, to_dict=False)
        #smile = gdsc.getCompoundSMILE(use_exdb=True, drugList=None)
        #smile.to_csv('/repo4/ytang4/PHD/pathwayNet/data/GDSC/GDSC.DRUG.isoSMILE.txt', header=True, index=False, sep="\t")
        #target = gdsc.getCompoundTARGET(use_exdb=True, drugList=None)
        resp_dict1 = gdsc.getContRESPDICT(use='IC50', drugList=None, cellList=None)
        resp_drug = resp_dict1['exp1'].index.tolist()
        #model = gdsc.getMODEL()
        #print(resp_drug[:5])
        #print('missing={:}'.format(set(resp_drug)-set(target['drug'])))
        #resp_dict2 = gdsc.getContRESPDICT(use='LN_IC50', drugList=None, cellList=None)
        #resp_dict3 = gdsc.getContRESPDICT(use='Z_SCORE', drugList=None, cellList=None)
        #resp_dict4 = gdsc.getContRESPDICT(use='AUC', drugList=None, cellList=None)
        # df_list = []
        #for idx in cnv_c.index:
        #    data = cnv_c.loc[[idx]].T 
        #    ones = data[data[idx]==1]
        #    genes = ones.index.tolist()
        #    df = pd.DataFrame({'cell':[idx]*len(genes), 'gene':genes})
        #    df_list.append(df)
        #cnv = pd.concat(df_list, axis=0)
        #cnv['gene'] = cnv['gene'].replace(gid_gnm_dict)
        DB_PATH = '/data/DR/db/GDSC/processed/' #'/repo4/ytang4/PHD/perturbNet/data/GDSC/'
        #target.to_csv(DB_PATH+'GDSC.RAW.TARGET.GeneName.Mat.txt',header=True, index=False, sep="\t")
        #exp.to_csv(DB_PATH+'GDSC.RAW.EXP.TPM.Mat.txt',header=True, index=True, sep="\t")
        #cnv_c.to_csv(DB_PATH+'GDSC.CNV.GisticScore2.OneHot.Mat.txt',header=True, index=True, sep="\t")
        #mut.to_csv(DB_PATH+'GDSC.MUT.GeneName.Mat.txt',header=True, index=False, sep="\t")
        #cnv.to_csv(DB_PATH+'GDSC.CNV.GeneName.Mat.txt',header=True, index=False, sep="\t")
        resp_dict1['exp1'].to_csv(DB_PATH+'GDSC.RAW.RESP.IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict2['exp1'].to_csv(DB_PATH+'GDSC.RAW.RESP.LN_IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict3['exp1'].to_csv(DB_PATH+'GDSC.RAW.RESP.Z_SCORE.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict4['exp1'].to_csv(DB_PATH+'GDSC.RAW.RESP.AUC.exp1.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict1['exp2'].to_csv(DB_PATH+'GDSC.RAW.RESP.IC50.exp2.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict2['exp2'].to_csv(DB_PATH+'GDSC.RAW.RESP.LN_IC50.exp2.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict3['exp2'].to_csv(DB_PATH+'GDSC.RAW.RESP.Z_SCORE.exp2.Mat.txt',header=True, index=True, sep="\t")
        #resp_dict4['exp2'].to_csv(DB_PATH+'GDSC.RAW.RESP.AUC.exp2.Mat.txt',header=True, index=True, sep="\t")
        #model.to_csv(DB_PATH+'GDSC.MODEL.Annotation.Mat.txt',header=True, index=False, sep="\t")

          
    if save:
        # get cell, drug
        cell_list = gdsc.getCELL(use=['exp', 'cnv', 'mut', 'resp'])
        drug_list = gdsc.getDRUG(use=['target', 'smile', 'resp'])
        # generate processed data
        smile = gdsc.getCompoundSMILE(use_exdb=True, drugList=drug_list)
        target = gdsc.getCompoundTARGET(use_exdb=True, drugList=drug_list)
        exp = gdsc.getEXP(use='TPM', cellList=cell_list)
        cnv_c = gdsc.getCNV(use='gistic', cellList=cell_list, categorical=True)
        cnv = gdsc.getCNV(use='gistic', cellList=cell_list, categorical=False)
        mut = gdsc.getMUT(cellList=cell_list, to_dict=False)
        resp_dict1 = gdsc.getContRESPDICT(use='IC50', drugList=drug_list, cellList=cell_list)
        resp_dict2 = gdsc.getContRESPDICT(use='LN_IC50', drugList=drug_list, cellList=cell_list)
        resp_dict3 = gdsc.getContRESPDICT(use='Z_SCORE', drugList=drug_list, cellList=cell_list)
        resp_dict4 = gdsc.getContRESPDICT(use='AUC', drugList=drug_list, cellList=cell_list)
        model = gdsc.getMODEL()
        gene = gdsc.parseGENE()
        # convert gene id to gene name
        gid_gnm_dict = dict( zip(gene['gene_id'], gene['hgnc_symbol']) )
        cnv_c = cnv_c.rename(columns=gid_gnm_dict)
        cnv = cnv.rename(columns=gid_gnm_dict)
        # save files to DP_PATH/processed/ folder
        smile.to_csv(DB_PATH+'/processed/GDSC.TARGET.SMILE.Mat.txt',header=True, index=True, sep="\t")
        target.to_csv(DB_PATH+'/processed/GDSC.TARGET.GeneName.Mat.txt',header=True, index=False, sep="\t")
        exp.to_csv(DB_PATH+'/processed/GDSC.EXP.TPM.Mat.txt',header=True, index=True, sep="\t")
        cnv_c.to_csv(DB_PATH+'/processed/GDSC.CNV.OneHot.Mat.txt',header=True, index=True, sep="\t")
        cnv.to_csv(DB_PATH+'/processed/GDSC.CNV.Log2.Mat.txt',header=True, index=True, sep="\t")
        mut.to_csv(DB_PATH+'/processed/GDSC.MUT.GeneName.Mat.txt',header=True, index=False, sep="\t")
        resp_dict1['exp1'].to_csv(DB_PATH+'/processed/GDSC.RESP.IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        resp_dict2['exp1'].to_csv(DB_PATH+'/processed/GDSC.RESP.LN_IC50.exp1.Mat.txt',header=True, index=True, sep="\t")
        resp_dict3['exp1'].to_csv(DB_PATH+'/processed/GDSC.RESP.Z_SCORE.exp1.Mat.txt',header=True, index=True, sep="\t")
        resp_dict4['exp1'].to_csv(DB_PATH+'/processed/GDSC.RESP.AUC.exp1.Mat.txt',header=True, index=True, sep="\t")
        resp_dict1['exp2'].to_csv(DB_PATH+'/processed/GDSC.RESP.IC50.exp2.Mat.txt',header=True, index=True, sep="\t")
        resp_dict2['exp2'].to_csv(DB_PATH+'/processed/GDSC.RESP.LN_IC50.exp2.Mat.txt',header=True, index=True, sep="\t")
        resp_dict3['exp2'].to_csv(DB_PATH+'/processed/GDSC.RESP.Z_SCORE.exp2.Mat.txt',header=True, index=True, sep="\t")
        resp_dict4['exp2'].to_csv(DB_PATH+'/processed/GDSC.RESP.AUC.exp2.Mat.txt',header=True, index=True, sep="\t")
        model.to_csv(DB_PATH+'/processed/GDSC.MODEL.Annotation.Mat.txt',header=True, index=False, sep="\t")
        gene.to_csv(DB_PATH+'/processed/GDSC.GENE.Annotation.Mat.txt',header=True, index=False, sep="\t")
        # display print message
        print('save processed files to {:}'.format(DB_PATH+'/processed/'))
        print('    #common drugs (combo removed) ={:}, common cells={:}'.format(len(drug_list), len(cell_list)))
        print('    resp (drug,cell)={:}, target (drug,target)={:}, smile (drug, smile)={:}'.format(resp_dict1['exp1'].shape, target.shape, smile.shape))
        print('    exp (cell,gene)={:}, cnv (cell,gene)={:}, mut(cell,gene)={:}'.format(exp.shape, cnv.shape, mut.shape))
