__version__ = '0.1.0'

from pathlib import Path
import os, sys

temp_loc = '/'.join(os.path.dirname(__file__).split('/')[:-1])

oncosplice_setup = {}
if sys.platform != 'darwin':
    print('OncoSplice, meet Power.')
    MACHINE = 'POWER'
    oncosplice_setup.update({
        'MACHINE': 'HOME',
        'HOME': True,
        'CHROM_SOURCE': Path('/tamir1/lab_resources/Genomes/Human/human_hg38/Chromosome'),
        'MRNA_PATH': Path('/tamir2/nicolaslynn/data/ensembl/mrna_database_ensembl_v110/protein_coding'),
        'CONS_PATH': Path('/tamir1/nicolaslynn/data/Conservation/data'),
        'MISSPLICING_PATH': Path('/tamir2/nicolaslynn/experimental_data/variant_missplicing'),
        'TRANEX_PATH': Path('/tamir2/nicolaslynn/data/HumanProteinAtlas/tranex_sum_tpm_database/'),
        'flag': '0'

    })
    from oncosplice.oncosplice import oncosplice

else:
    print('OncoSplice, meet Mac.')
    MACHINE = 'MAC'
    oncosplice_setup.update({
        'MACHINE': 'REMOTE',
        'HOME': False,
        'CHROM_SOURCE': Path('/Users/nl/Documents/data/Genomes/Human/human_hg38/Chromosome'),
        'MRNA_PATH': Path('/Users/nl/Documents/phd/projects/parse_genome_annotations/data/mrna_database_ensembl_v110/protein_coding'),
        'MISSPLICING_PATH': Path('/Users/nl/Downloads/'),
    })
