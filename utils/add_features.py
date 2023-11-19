#!/usr/bin/python3.9

# Annotate the table using Predictions:
# 1. Flagellins/Flagellum proteins directly remapped using metadata files
# 2. Secreted proteins predicted using signalp6 in slow mode
# 3. Membrane proteins predicted using topgraph annotation tool
# 4. Toxins predicted running diamond on uniprot swisprot reported Toxins DB

# NB: for Signalp
# The prediction score for signalp is nested inside a column composed by strings
# containing different information: take only the Prediction score and add it as
# another column, to give the possibility to filter it based on value
import pandas as pd
import numpy as np

def replace_vals(row):
        if type(row['uniref_org']) != float:
            if ' ' not in row['uniref_org'] and type(row['IEDB_org']) != float:
                return row['IEDB_org']
            elif ' ' not in row['uniref_org'] and type(row['VFDB_org']) != float:
                return row['VFDB_org']
            elif ' ' not in row['uniref_org'] and type(row['bac_source_org']) != float:
                return row['bac_source_org']
            else:
                return row['uniref_org']
        else:
            return row['uniref_org']


def add_pred_score(df_signal):
    df_signal = df_signal[df_signal['CS Position'].isna() == False].copy()
    df_signal.loc[:,'Prediction Score'] = [float(x.split('Pr: ')[1]) for x in df_signal['CS Position']]
    return df_signal


def add_signalp6(df_agilent, threshold = 0.85):
    # Extract info from signalp6 file (TODO: add it in data folder)
    signalp6_slow = pd.read_csv('signalp6_slow.txt', sep='\t')
    signalp6_slow = add_pred_score(signalp6_slow)
    secreted_slow = signalp6_slow[(signalp6_slow['Prediction'] != 'OTHER') & 
                              (signalp6_slow['Prediction Score'] >= threshold)]['# ID']
    df_agilent['signalp6_slow'] = np.where(df_agilent['prot_num'].isin(secreted_slow), 1, 0)
    return df_agilent

def add_flagellum(df_agilent, all_columns=True):
    if all_columns:
        # Take Flagellum name from all columns + diamond (old)
        flagel_containing_rows = df_agilent[df_agilent.applymap(lambda cell: isinstance(cell, str) and 'Flagel' in cell).any(axis=1)|\
            df_agilent.applymap(lambda cell: isinstance(cell, str) and 'flagel' in cell).any(axis=1)]
    else:
        # Take it from full name column + diamond (old)
        flagel_containing_rows = df_agilent[((df_agilent['full name'].str.contains('Flagel')) | 
        (df_agilent['full name'].str.contains('flagel'))) &\
            (df_agilent['full name'].isnull() == False)]
    df_agilent['is_flagellum'] = np.where((df_agilent['peptide_name'].isin(flagel_containing_rows['peptide_name'])) | (df_agilent['diamond_flagellum'] == 1), 1,0)
    return df_agilent


def add_toxin(df_agilent):
    # Toxins: diamond+mmseqs2 predictions
    diamond_res_tox = pd.read_table('diamond_mmseqs/toxin/matches_toxin_SV.tsv', sep='\t',
                                header=None)
    mmseqs_res_tox = pd.read_table('diamond_mmseqs/toxin/mmseqs_res.tsv', sep='\t',
                                header=None)
    mmseqs_res_tox.columns = diamond_res_tox.columns = ['query_accession', 'target_accession',
                       'seq_identity', 'length',
                       'mismatches', 'gap_openings',
                       'query_start', 'query_end',
                       'target_start', 'target_end',
                       'E_value', 'Bit_score']
    # Go for inner merge
    merged_mmseqs_diamond = pd.merge(mmseqs_res_tox, diamond_res_tox, on='query_accession', how='inner')
    # Calculate E-value mean between the 2 methods, then
    # group by name and mean and take the minimum E-value
    merged_mmseqs_diamond['Mean_E_value'] = merged_mmseqs_diamond[['E_value_x', 'E_value_y']].mean(axis=1)
    merged_mmseqs_diamond = merged_mmseqs_diamond.loc[merged_mmseqs_diamond.groupby('query_accession')['Mean_E_value'].idxmin()]
    # Add in the dataframe as presence/absen column
    df_agilent['diamond_mmseqs_intersec_toxin'] = np.where(df_agilent['prot_num'].isin(merged_mmseqs_diamond['query_accession']), 1,0)
    return df_agilent

def add_membrane_prot(df_agilent):

    # Topgraph membrane protein prediction tool
    topgraph_res = pd.read_csv('membrane_prediction/final_out_topgraph.tsv', sep='\t',header=None)
    topgraph_res.columns = ['protein_name', 'is_topgraph']
    topgraph_pos_prot = topgraph_res[topgraph_res['is_topgraph'] == 1] #['protein_name']
    df_agilent['is_topgraph_new'] = np.where(df_agilent['prot_num'].isin(topgraph_pos_prot['protein_name']), 1,0)
    # Sum this to the old result (old and new)
    # ADD also the case in which is_topgraph_new and is_topgraph_old are together
    df_agilent['is_topgraph_new_&_old'] = np.where(((df_agilent['is_topgraph'] == 1) | (df_agilent['is_topgraph_new'] == 1)), 1,0)
    return df_agilent

def map_predictions(df_agilent):
    from functools import reduce
    # Flagellum: map it on the metadata and get the rows in which it appears 'flagel' or 'Flagel'
    #   then add it as a binary info in 'is_flagellum' column
    # Secreted proteins : use prediction threshold of 0.85
    # Toxins: using diamond + mmseqs2 E_value average
    # Topgraph: prediction tool
    transformation_functions = [add_flagellum, add_signalp6, add_membrane_prot, add_toxin, add_complete_organism_info] 
    df_agilent = reduce(lambda df,func: func(df), transformation_functions, df_agilent)
    return df_agilent

def add_complete_organism_info(df_agilent):
    organism_info_df = df_agilent[['uniref_func', 'IEDB_organism_name', 'toxin_prot_name', 'bac_src']]
    organism_info_df_index = organism_info_df[((organism_info_df['uniref_func'].str.contains('Tax=')) & \
        (organism_info_df['uniref_func'].isnull() == False))].index
    organism_tox_index = organism_info_df[(organism_info_df['toxin_prot_name'].isnull() == False)].index

    organism_IEDB_index = organism_info_df[(organism_info_df['IEDB_organism_name'].isnull() == False)].index
    organism_microb_index = organism_info_df[(organism_info_df['bac_src'].isnull() == False)].index

    # Create a new df for easier calculations
    df_agilent_organism = df_agilent.copy()
    df_agilent_organism['uniref_org'] = df_agilent_organism.iloc[organism_info_df_index,:]['uniref_func'].apply(lambda x: x.split(' Tax=')[1].split(' T')[0])
    df_agilent_organism['IEDB_org'] = df_agilent_organism.iloc[organism_IEDB_index,:]['IEDB_organism_name'].apply(lambda x: x.split('\'')[1].strip())
    df_agilent_organism['VFDB_org'] = df_agilent_organism.iloc[organism_tox_index,:]['toxin_prot_name'].apply(lambda x: x.split('[')[2].split(']')[0].strip())
    df_agilent_organism['bac_source_org'] = df_agilent_organism.iloc[organism_microb_index,:]['bac_src'].apply(lambda x: x.split(' &')[0].strip())
    # Add to agilent total df the complete name of the organism
    df_agilent_organism['Organism_complete_name'] = df_agilent_organism.apply(replace_vals, axis=1)



    df_agilent['Organism_complete_name'] = df_agilent_organism['Organism_complete_name']
    return df_agilent