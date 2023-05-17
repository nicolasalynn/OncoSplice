import shutil
import os
import pandas as pd
import numpy as np

def main(mut_id, sai_threshold=0.25, force=False, save_flag=True):
    from geney import unload_json, parse_in_args, get_correct_gene_file
    from oncosplice import oncosplice_setup
    from oncosplice.spliceai_utils import find_missplicing_spliceai_adaptor
    from oncosplice.Gene import AnnotatedGene
    from oncosplice.variant_utils import EpistaticSet, Mutation
    from oncosplice.protein_scoring_utils import generate_report

    if '|' in mut_id:
        input = EpistaticSet(mut_id)
    else:
        input = Mutation(mut_id)

    annot_file = get_correct_gene_file(input.gene, target_directory=oncosplice_setup['MRNA_PATH'])
    if not annot_file:
        print(f'No annotations for gene: {input.gene}...')
        return pd.DataFrame()

    print(f'>> Processing: {input}')

    ################### MISSPLICING
    missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=sai_threshold, force=force, save_flag=save_flag)
    print(f'\tMissplicing: {missplicing}')

    ################### VARIANT ANNOTATIONS
    reference_gene = AnnotatedGene(annot_file)          # the next step is to have spliceai not have to repopen annotation data and rather to use the sequence build in the reference seq object... also we can check if a mutation is implemented before running splice ai.
    variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
    ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome()

    ################### GENERATE VARIANT REPORT
    report = generate_report(ref_proteome, var_proteome, missplicing, input)
    if report.empty:
        return report

    report = pd.merge(report, reference_gene.tranex_tpm, on=['ensembl_transcript_id'], how='left')
    return report

def converter(instr, s):
    return np.fromstring(instr[1:-1], count=s, sep=' ')
def calculate_final_score(file='', df=None):
    from pandas.errors import EmptyDataError
    import json

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
    df['weighted_score'] = df.oncosplice_score * df.isoform_prevalence

    tracker = {}
    tracker['mut_id'] = df.iloc[0].mut_id
    tracker['interesting'] = df.iloc[0].missplicing_flag
    tracker['highest_splicing_penetrance'] = highest_ms

    tracker['oncosplice_score_gof'] = max(0, df.groupby('transcript_id').oncosplice_score.mean().max())
    tracker['oncosplice_score_weighted_gof'] = max(0, df.groupby('transcript_id').weighted_score.mean().max())
    tracker['oncosplice_score_lof'] = min(0, df.groupby('transcript_id').oncosplice_score.mean().min())
    tracker['oncosplice_score_weighted_lof'] = min(0, df.groupby('transcript_id').weighted_score.mean().min())

    tracker['legacy_oncosplice_score'] = df.groupby('transcript_id').legacy_oncosplice_score.mean().max()
    # tracker['legacy_oncosplice_score_cons'] = df[df.cons_available].groupby('transcript_id').legacy_oncosplice_score.mean().max()

    temp = pd.Series(np.array([np.array(v) for v in list(tracker.values())]), index=list(tracker.keys()))
    temp.name = temp.mut_id
    temp.index.name = 'mut_id'
    return temp

def short_results(mut_id, sai_threshold=0.25):
    return calculate_final_score(main(mut_id, sai_threshold=sai_threshold))

