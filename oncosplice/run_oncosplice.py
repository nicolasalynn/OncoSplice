import shutil
import os
import pandas as pd
import numpy as np
from geney import unload_json, parse_in_args, get_correct_gene_file
from oncosplice import oncosplice_setup
from oncosplice.spliceai_utils import find_missplicing_spliceai_adaptor
from oncosplice.Gene import AnnotatedGene
from oncosplice.variant_utils import EpistaticSet, Mutation
from oncosplice.protein_scoring_utils import generate_report
from pandas.errors import EmptyDataError
import json
from geney.general.performance_utils import check_ti

fake_mut_id = 'KRAS:12:25227343:G:T'

def main(mut_id, sai_threshold=0.25, force=False, save_flag=True):
    tic = check_ti()

    if '|' in mut_id:
        input = EpistaticSet(mut_id)
    else:
        input = Mutation(mut_id)

    annot_file = get_correct_gene_file(input.gene, target_directory=oncosplice_setup['MRNA_PATH'])
    if not annot_file:
        print(f'No annotations for gene: {input.gene}...')
        return pd.DataFrame()

    print(f'>> Processing: {input}')

    tic = check_ti(tic, 'set up space')

    ################### MISSPLICING
    missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=sai_threshold, force=force, save_flag=save_flag)
    print(f'\tMissplicing: {missplicing}')
    tic = check_ti(tic, 'run oncosplice')

    ################### VARIANT ANNOTATIONS
    reference_gene = AnnotatedGene(annot_file)          # the next step is to have spliceai not have to repopen annotation data and rather to use the sequence build in the reference seq object... also we can check if a mutation is implemented before running splice ai.
    variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
    ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome()

    ################### GENERATE VARIANT REPORT
    report = generate_report(ref_proteome, var_proteome, missplicing, input)
    if report.empty:
        return report

    report = pd.merge(report, reference_gene.tranex_tpm, on=['ensembl_transcript_id'], how='left')
    tic = check_ti(tic, 'generate transcriptome predictions')

    return report

def converter(instr, s):
    return np.fromstring(instr[1:-1], count=s, sep=' ')
def calculate_final_score(file='', df=None):
    if file:
        try:
            df = pd.read_csv(file)
        except EmptyDataError:
            return pd.Series(dtype='float64')
        if df.empty:
            return pd.Series(dtype='float64')

    if 'transcipt_id' in df.columns:
        df.rename(columns={'transcipt_id': 'transcript_id'}, inplace=True)

    df = df.loc[~df.gene.isna()]

    missplicing = json.loads(df.iloc[0].full_missplicing)
    highest_ms = 0
    for i, k in missplicing.items():
        for j, k2 in k.items():
            if abs(k2['delta']) > highest_ms:
                highest_ms = abs(k2['delta'])

    df = df[df.isoform_prevalence >= 0.05]
    df['oncosplice_preservation1'] = df.oncosplice_score / (1 - df.preservation)

    tracker = {}
    tracker['mut_id'] = df.iloc[0].mut_id
    tracker['highest_splicing_penetrance'] = highest_ms
    tracker['new_oncosplice_score'] = df.groupby('transcript_id').oncosplice_score.mean().max()
    tracker['new_oncosplice_score2'] = df.groupby('transcript_id').oncosplice_preservation1.mean().max()
    tracker['new_oncosplice_score3'] = df.groupby('transcript_id').oncosplice_window.mean().max()
    tracker['preservation'] = df.groupby('transcript_id').preservation.mean().min()
    tracker['legacy_oncosplice_score'] = df.groupby('transcript_id').legacy_oncosplice_score.mean().max()
    temp = pd.Series(np.array([np.array(v) for v in list(tracker.values())]), index=list(tracker.keys()))
    temp.name = temp.mut_id
    temp.index.name = 'mut_id'
    return temp

def short_results(mut_id, sai_threshold=0.25):
    return calculate_final_score(main(mut_id, sai_threshold=sai_threshold))


#
# def get_higehest_ms(row):
#     missplicing = json.loads(row.full_missplicing)
#     highest_ms = 0
#     for i, k in missplicing.items():
#         for j, k2 in k.items():
#             if abs(k2['delta']) > highest_ms:
#                 highest_ms = abs(k2['delta'])
#     return highest_ms
# def apply_to_mut(g):
#     return g.groupby('transcript_id').oncosplice_score.mean().max()
# def apply_to_mut2(g):
#     return g.groupby('transcript_id').legacy_oncosplice_score.mean().max()
# def calculate_final_score(df):
#     if 'transcipt_id' in df.columns:
#         df.rename(columns={'transcipt_id': 'transcript_id'}, inplace=True)
#     df = df.loc[~df.gene.isna()]
#     df = df[df.isoform_prevalence >= 0.25]
#     df.oncosplice_score /= (1 - df.preservation)
#     temp = df[['mut_id', 'full_missplicing']].drop_duplicates(keep='first')
#     temp['highest_ms'] = temp.apply(get_higehest_ms, axis=1)
#     temp.drop(columns=['full_missplicing'], inplace=True)
#     temp2 = df.groupby('mut_id').apply(apply_to_mut)
#     temp3 = df.groupby('mut_id').apply(apply_to_mut2)
#     temp2.name='oncosplice_score'
#     temp3.name='legacy_oncosplice_score'
#     return pd.merge(temp, pd.concat([temp2, temp3], axis=1), on=['mut_id'])
#
# def split(a, n):
#     k, m = divmod(len(a), n)
#     return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))
#
# all_files = os.listdir()
# all_data, all_data_in = [], []
# for g in split(all_files, 100):
#     for file in tqdm(g):
#         all_data_in.append(pd.read_csv(file))
#     all_data_in = pd.concat(all_data_in)
#     all_data.append(calculate_final_score(all_data_in))
#     if isinstance(all_data, list) and len(all_data) > 1:
#         all_data = pd.concat(all_data)
#         all_data.to_csv('/tamir2/nicolaslynn/temp/clinvar_scores.csv')
#         all_data = [all_data]
#     all_data_in = []
