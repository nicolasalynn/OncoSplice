import shutil
import os
import pandas as pd
import numpy as np



def main(mut_id, sai_threshold=0.25, force=False):
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
    missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=sai_threshold, force=force)
    print(f'\tMissplicing: {missplicing}')

    # try:
    ################### VARIANT ANNOTATIONS
    reference_gene = AnnotatedGene(annot_file)          # the next step is to have spliceai not have to repopen annotation data and rather to use the sequence build in the reference seq object... also we can check if a mutation is implemented before running splice ai.
    variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
    ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome()

    ################### GENERATE VARIANT REPORT
    report = generate_report(ref_proteome, var_proteome, missplicing, input)
    if report.empty:
        return report

    report = pd.merge(report, reference_gene.tranex_tpm, on=['ensembl_transcript_id'], how='outer')
    return report

    # except:
    #     return pd.DataFrame()

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

    elif df == None:
        print(f'Must define a file or a dataframe.')
        return pd.Series(dtype='float64')

    tracker = {}
    tracker['mut_id'] = df.iloc[0].mut_id
    missplicing = json.loads(df.iloc[0].full_missplicing)
    tracker['interesting'] = any([missplicing.get('missed_donors', {}), missplicing.get('missed_acceptors', {}), missplicing.get('discovered_donors', {}), missplicing.get('discovered_acceptors', {})])

    lof_score_min = 0
    gof_score_max = 0
    for tid, g in df.groupby('transcript_id'):
        g = g[g.isoform_prevalence == max(g.isoform_prevalence)]
        lof_score = min(g.lof_score)
        gof_score = max(g.gof_score)
        if lof_score < lof_score_min:
            lof_score_min = lof_score
        if gof_score > gof_score_max:
            gof_score_max = gof_score

    tracker['gof_score'] = gof_score_max
    tracker['lof_score'] = lof_score_min
    tracker['oncosplice_score'] = np.mean(df.oncosplice_score * df.isoform_prevalence)

    tracker['legacy_oncosplice_score'] = df.groupby('transcript_id').mean().max()

    temp = pd.Series(np.array(list(tracker.values())), index=list(tracker.keys()))
    temp.name = temp.mut_id
    temp.index.name = 'mut_id'
    return temp

def short_results(mut_id, sai_threshold=0.25):
    return calculate_final_score(main(mut_id, sai_threshold=sai_threshold))

