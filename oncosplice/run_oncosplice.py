import shutil
import os
from geney import unload_json, parse_in_args, get_correct_gene_file
import pandas as pd
import numpy as np
from pandas.errors import EmptyDataError

from oncosplice import oncosplice_setup
from oncosplice.spliceai_utils import find_missplicing_spliceai_adaptor
from oncosplice.Gene import AnnotatedGene
from oncosplice.variant_utils import EpistaticSet, Mutation
from oncosplice.protein_scoring_utils import generate_report

def main(mut_id, sai_threshold=25):
    if '|' in mut_id:
        input = EpistaticSet(mut_id)
    else:
        input = Mutation(mut_id)

    annot_file = get_correct_gene_file(input.gene, target_directory=oncosplice_setup['MRNA_PATH'])
    if not annot_file:
        print(f'No annotations for gene: {input.gene}...')
        return pd.DataFrame()

    print(f'>> Processing: {input}')

    try:
        ################### MISSPLICING
        missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=round(sai_threshold/100, 3), force=False)
        print(f'\tMissplicing: {missplicing}')

        ################### VARIANT ANNOTATIONS
        reference_gene = AnnotatedGene(annot_file)          # the next step is to have spliceai not have to repopen annotation data and rather to use the sequence build in the reference seq object... also we can check if a mutation is implemented before running splice ai.
        variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
        ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome()

        ################### GENERATE VARIANT REPORT
        report = generate_report(ref_proteome, var_proteome, missplicing, input)
        if report.empty:
            return report

        report = pd.merge(report, reference_gene.tranex_tpm, on=['ensembl_transcript_id'])
        return report

    except:
        return pd.DataFrame()

def converter(instr, s):
    return np.fromstring(instr[1:-1], count=s, sep=' ')
def calculate_final_score(file):
    try:
        df = pd.read_csv(file)
    except EmptyDataError:
        return pd.Series()
    if df.empty:
        return pd.Series()

    tracker = {}

    if 'tpm_sum' not in df.columns:
        df['tpm_sum'] = df.apply(lambda row: 1, axis=1)

    unique_tmps = sum(df['tpm_sum'].unique())
    df['tpm_weight'] = df['tpm_sum'] / unique_tmps
    tracker['mut_id'] = df.iloc[0].mut_id

    tracker['oncosplice_sum'] = df.oncosplice_score.sum()

    temp = df.oncosplice_score * df.isoform_prevalence
    tracker['oncosplice_isoform_sum'] = sum(temp)

    temp = df.oncosplice_score * df.isoform_prevalence * df.tpm_weight * 100
    tracker['oncosplice_isoform_tpm_sum'] = sum(temp)

    temp = df.gof_score * df.isoform_prevalence * df.tpm_weight * 100
    tracker['gof_score'] = sum(temp)

    temp = df.lof_score * df.isoform_prevalence * df.tpm_weight * 100
    tracker['lof_score'] = sum(temp)

    temp = pd.Series(np.array(list(tracker.values())), index=list(tracker.keys()))
    temp.name = temp.mut_id
    temp.index.name = 'mut_id'
    return temp

def short_results(mut_id, sai_threshold=25):
    return calculate_final_score(main(mut_id, sai_threshold=sai_threshold))


