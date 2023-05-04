import shutil
import os
from geney import unload_json, parse_in_args, get_correct_gene_file
import pandas as pd

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
        return None

    print(f'>> Processing: {input}')

    # try:
    ################### MISSPLICING
    missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=round(sai_threshold/100, 3), force=False)
    print(f'\tMissplicing: {missplicing}')

    ################### VARIANT ANNOTATIONS
    reference_gene = AnnotatedGene(annot_file)          # the next step is to have spliceai not have to repopen annotation data and rather to use the sequence build in the reference seq object... also we can check if a mutation is implemented before running splice ai.
    variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
    ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome()

    ################### GENERATE VARIANT REPORT
    report = generate_report(ref_proteome, var_proteome, missplicing, input)
    report = pd.merge(report, reference_gene.tranex_tpm, on=['ensembl_transcript_id'])

    return report
    #
    # except:
    #     return None


