import shutil
import os
from geney import unload_json, parse_in_args, get_correct_gene_file

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
    print(f'>> Processing: {input}')

    ################### MISSPLICING
    missplicing = find_missplicing_spliceai_adaptor(input=input, sai_threshold=round(sai_threshold/100, 3))

    ################### VARIANT ANNOTATIONS
    annot_file = get_correct_gene_file(input.gene, target_dir=oncosplice_setup['MRNA_PATH'])
    if not annot_file:
        return f"No data for gene {input.gene}."

    reference_gene = AnnotatedGene(annot_file)
    variant_gene = reference_gene.create_gene_isoform(mut_ids=mut_id, aberrant_splicing=missplicing)
    ref_proteome, var_proteome = reference_gene.develop_proteome(), variant_gene.develop_proteome(experimental=True)

    ################### GENERATE VARIANT REPORT
    report = generate_report(ref_proteome, var_proteome, missplicing, mut_id.split('|'))
    return variant_gene.__dict__, report

if __name__ == '__main__':
    # prs, _ = parse_in_args()
    # args_data = unload_json(prs.args_file)
    # out_dir = args_data['results']
    # mut_data = args_data['mutations']
    # annotations, summary = main(mut_data)
    # sh_file = args_data['sh_file']
    # args_file = args_data['args_file']
    # done_path = os.path.join(os.path.dirname(sh_file), 'done')
    # shutil.move(args_file, os.path.join(done_path, os.path.basename(args_file)))
    # shutil.move(sh_file, os.path.join(done_path, os.path.basename(sh_file)))
    # print("Done.")
    print("Not executable.")


else:
    print('Importing functions.')

