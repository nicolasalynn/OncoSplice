__version__ = '0.1.0'

from pathlib import Path
import os, sys


def get_data_dir():
    data_dir = os.path.join(os.path.expanduser('~'), '.oncosplice', 'config.txt')
    with open(data_dir, 'r') as config_file:
        return config_file.readline().split('=')[1].strip()


# temp_loc = '/'.join(os.path.dirname(__file__).split('/')[:-1])
#
# oncosplice_setup = {}
# if sys.platform == 'linux':
#     print('OncoSplice, meet Power.')
#     oncosplice_setup.update({
#         'CHROM_SOURCE': Path('/tamir1/lab_resources/Genomes/Human/human_hg38/Chromosome'),
#         'MRNA_PATH': Path('/tamir2/nicolaslynn/data/ensembl/mrna_database_ensembl_v110/protein_coding'),
#         'CONS_PATH': Path('/tamir2/nicolaslynn/data/Conservation/data2'),
#         'MISSPLICING_PATH': Path('/tamir2/nicolaslynn/experimental_data/variant_missplicing'),
#         'TRANEX_PATH': Path('/tamir2/nicolaslynn/data/HumanProteinAtlas/tranex_sum_tpm_database/'),
#     })
#     from oncosplice.oncosplice import oncosplice
#
# else:
#     print('OncoSplice, meet Mac.')
#     oncosplice_setup.update({
#         'CHROM_SOURCE': Path('/Users/nl/Documents/data/Genomes/Human/human_hg38/Chromosome'),
#         'MRNA_PATH': Path('/Users/nl/Documents/phd/projects/parse_genome_annotations/data/mrna_database_ensembl_v110/protein_coding'),
#         'MISSPLICING_PATH': Path('/Users/nl/Downloads/'),
#
#     })
