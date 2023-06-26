import shutil
import os
import pandas as pd
import numpy as np
from geney import unload_json, parse_in_args, get_correct_gene_file, append_line
from oncosplice import oncosplice_setup
from oncosplice.spliceai_utils import find_missplicing_spliceai_adaptor, check_splicing_difference
from oncosplice.Gene import AnnotatedGene
from oncosplice.variant_utils import EpistaticSet, Mutation
from oncosplice.protein_scoring_utils import generate_report
from pandas.errors import EmptyDataError
import json
from geney.general.performance_utils import check_ti
sample_mut_id = 'KRAS:12:25227343:G:T'

def main(mut_id, sai_threshold=0.25, min_coverage=2500, force=False, save_flag=True, show_output=False):
    oncosplice_setup['show_output'] = show_output
    if '|' in mut_id:
        input = EpistaticSet(mut_id)
    else:
        input = Mutation(mut_id)

    annot_file = get_correct_gene_file(input.gene, target_directory=oncosplice_setup['MRNA_PATH'])
    if not annot_file:
        print(f'No annotations for gene: {input.gene}...')
        return pd.DataFrame()

    if len(input.file_identifier) > 70:
        append_line(oncosplice_setup['failed_mut_path'], mut_id)
        return pd.DataFrame()

    missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=sai_threshold, min_coverage=min_coverage, force=force, save_flag=save_flag)
    print(f'>> Processing: {input}')
    print(f'\tMissplicing: {missplicing}')

    reference_gene = AnnotatedGene(annot_file)
    variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
    ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome()
    report = generate_report(ref_proteome, var_proteome, missplicing, input)

    if report.empty:
        return report

    report = pd.merge(report, reference_gene.tranex_tpm, on=['ensembl_transcript_id'], how='left')
    return report, missplicing

def run_pairwise_and_constituents(epistasis):
    m1, m2 = epistasis.split('|')
    data = {}
    data[epistasis], me = main(epistasis)
    data[m1], m1 = main(m1)
    data[m2], m2 = main(m2)
    data['m1_vs_me'], _ = check_splicing_difference(me, m1, 0.4)
    data['m2_vs_me'], _ = check_splicing_difference(me, m2, 0.4)
    return data

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
    tracker = {}
    tracker['mut_id'] = df.iloc[0].mut_id
    tracker['highest_splicing_penetrance'] = highest_ms
    tracker['missplicing'] = json.dumps(missplicing)
    tracker['legacy_oncosplice_score'] = df.groupby('transcript_id').legacy_oncosplice_score.mean().max()
    df['weighted_lof'] = df.oncosplice_score_lof * df.isoform_prevalence
    df['weighted_gof'] = df.oncosplice_score_gof * df.isoform_prevalence
    tracker['lof'] = df.groupby('transcript_id').weighted_lof.sum().mean()
    tracker['gof'] = df.groupby('transcript_id').weighted_gof.sum().mean()
    temp = pd.Series(np.array([np.array(v) for v in list(tracker.values())]), index=list(tracker.keys()))
    temp.name = temp.mut_id
    temp.index.name = 'mut_id'
    return temp

def short_results(mut_id, sai_threshold=0.25):
    return calculate_final_score(main(mut_id, sai_threshold=sai_threshold))


# def converter(instr, s):
#     return np.fromstring(instr[1:-1], count=s, sep=' ')