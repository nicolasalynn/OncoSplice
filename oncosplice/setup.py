from pathlib import Path
import os
from gtfparse import read_gtf
from geney import dump_json, dump_pickle, unload_pickle
import pandas as pd
from tqdm import tqdm
import requests
import shutil
import argparse
import gzip

# script_dir = os.path.dirname(__file__)
# resources_dir = os.path.join(script_dir, 'resources')
resources_dir = '/tamir2/nicolaslynn/tools/OncoSplice/oncosplice/resources'
CONS_DATA = unload_pickle(os.path.join(resources_dir, 'conservation.pkl'))

class Download:
    def __init__(self, external_url, local_path):
        self.external_url = external_url
        self.local_path = Path(local_path)
        self.local_file = self.get_local_file_path()
        self.download_and_process_file()

    def get_local_file_path(self):
        file_name = Path(self.external_url).name
        local_file_path = self.local_path / file_name
        # Check if a non-gzipped version of the file already exists
        if local_file_path.suffix == '.gz':
            potential_unzipped_file = Path(local_file_path.as_posix().rstrip('.gz'))
            if potential_unzipped_file.exists():
                return potential_unzipped_file
        return local_file_path

    def download_and_process_file(self):
        if not self.local_file.exists():
            try:
                response = requests.get(self.external_url, stream=True)
                response.raise_for_status()  # Raises a HTTPError if the HTTP request returned an unsuccessful status code
                with open(self.local_file, 'wb') as f:
                    f.write(response.content)
                self.process_file()
            except Exception as e:
                print(f"Error during download: {e}")

    def process_file(self):
        if self.local_file.suffix == '.gz':
            try:
                shutil.unpack_archive(str(self.local_file), extract_dir=self.local_path)
                self.local_file = Path(self.local_file.as_posix().rstrip('.gz'))
            except Exception as e:
                print(f"Error during file processing: {e}")


def process_transcript(transcript_df, rev):
    if transcript_df.empty:
        return None

    transcript = transcript_df[transcript_df.feature == 'transcript'].squeeze()
    if transcript.empty:
        return None

    exon_df = transcript_df[transcript_df.feature == 'exon']
    cds_df = transcript_df[transcript_df.feature == 'CDS']

    # Simplifying start and end assignments
    transcript_start, transcript_end = (transcript.end, transcript.start) if rev else (transcript.start, transcript.end)

    # Handling exons
    exon_starts, exon_ends = (exon_df.end, exon_df.start) if rev else (exon_df.start, exon_df.end)
    exon_starts, exon_ends = exon_starts.tolist(), exon_ends.tolist()

    if transcript_start not in exon_starts or transcript_end not in exon_ends:
        raise ValueError('Transcript start or end not in exons')

    acceptors, donors = list(exon_starts), list(exon_ends)
    acceptors.remove(transcript_start)
    donors.remove(transcript_end)

    if len(acceptors) != len(donors):
        raise ValueError('Different number of acceptors and donors')

    data = {
        'transcript_id': transcript.transcript_id,
        'transcript_biotype': transcript.transcript_biotype,
        'transcript_start': int(transcript_start),
        'transcript_end': int(transcript_end),
        'tag': transcript.tag,
    }

    if acceptors and donors:
        data.update({'donors': donors, 'acceptors': acceptors})

    # Handling CDS
    if not cds_df.empty:
        cds_start, cds_end = (cds_df.end, cds_df.start) if rev else (cds_df.start, cds_df.end)
        cds_start, cds_end = [c for c in cds_start.tolist() if c not in acceptors], [c for c in cds_end.tolist() if
                                                                                     c not in donors]
        if len(cds_start) != 1 or len(cds_end) != 1:
            return None
        cds_start, cds_end = cds_start[0], cds_end[0]
        data.update({'TIS': cds_start, 'TTS': cds_end, 'protein_id': transcript.protein_id})

    if transcript.transcript_id in CONS_DATA:
        data.update({'cons_available': True, 'cons_vector': CONS_DATA[transcript.transcript_id]})
    else:
        data.update({'cons_available': False})

    return data


def retrieve_and_parse_ensembl_annotations(local_path, valid_biotypes=['protein_coding']):
    gtex_url = 'https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'
    gtex_address = Download(gtex_url, local_path)
    gtex_df = pd.read_csv(gtex_address.local_file, delimiter='\t', header=2)
    gtex_df.Name = gtex_df.apply(lambda row: row.Name.split('.')[0], axis=1)
    gtex_df = gtex_df.set_index('Name').drop(columns=['Description'])

    ensembl_url = 'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz'
    ensembl_address = Download(ensembl_url, local_path)
    ensembl_version = ensembl_address.local_file.name.split('.')[-2]
    db_dir = ensembl_address.local_path / f'grch38_ensembl_v{ensembl_version}'
    if not db_dir.exists():
        db_dir.mkdir()

    annotations = read_gtf(ensembl_address.local_file)
    for gene_id, gene_df in tqdm(annotations.groupby('gene_id')):
        biotype = gene_df.gene_biotype.unique().tolist()
        chrm = gene_df.seqname.unique().tolist()
        strand = gene_df.strand.unique().tolist()
        gene_attribute = gene_df[(gene_df.feature == 'gene')]
        assert len(biotype) == 1, f'Multiple Biotypes: {biotype}'
        assert len(chrm) == 1, f'Multiple Chromosomes: {chrm}'
        assert len(strand) == 1, f'Multiple Strands: {strand}'
        assert len(gene_attribute) == 1, f"Multiple gene attributes: {gene_attribute.size}"

        if biotype[0] not in valid_biotypes:
            continue

        gene_attribute = gene_attribute.squeeze()
        rev = True if gene_attribute.strand == '-' else False
        json_data = {
            'gene_name': gene_attribute.gene_name,
            'chrm': gene_attribute.seqname.replace('chr', ''),
            'gene_id': gene_attribute.gene_id,
            'gene_start': gene_attribute.start,
            'gene_end': gene_attribute.end,
            'rev': rev,
            'tag': gene_attribute.tag.split(','),
            'primary_transcript': True if 'Ensembl' in gene_attribute.tag else False,
            'biotype': gene_attribute.gene_biotype,
            'transcripts': {transcript_id: process_transcript(transcript_df, rev) for transcript_id, transcript_df in
                            gene_df.groupby('transcript_id') if transcript_id},
            'tissue_expression': gtex_df.loc[gene_id].squeeze().to_dict() if gene_id in gtex_df.index else {},
        }

        json_data['transcripts'] = {tid: v for tid, v in json_data['transcripts'].items() if v is not None}

        if not json_data['transcripts']:
            continue

        biotype_path = db_dir / biotype[0]
        if not biotype_path.exists():
            biotype_path.mkdir()

        if gene_attribute.gene_name == '' or gene_id == '':
            continue

        file_name = biotype_path / f'mrnas_{gene_id}_{gene_attribute.gene_name}.pkl'
        dump_pickle(file_name, json_data)


def load_grch38_fasta(local_path):
    gtex_address = Download('https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz', local_path)
    split_gzipped_fasta(gtex_address.local_file, gtex_address.local_path / 'grch38')
    return None


def split_gzipped_fasta(input_file, output_directory):
    """
    Splits a gzipped FASTA file into multiple files, each containing a single sequence.
    File names are derived from the sequence identifier in the FASTA header.

    :param input_file: Path to the input gzipped FASTA file.
    :param output_directory: Directory where output FASTA files will be saved.
    """
    with gzip.open(input_file, 'rt') as file:  # 'rt' for reading text from a gzipped file
        sequence = ''
        header = ''
        for line in file:
            if line.startswith('>'):
                if sequence:
                    write_sequence(output_directory, header, sequence)
                    sequence = ''
                header = line[1:].split()[0]  # Assumes first word after '>' is the sequence identifier
            else:
                sequence += line.strip()

        if sequence:
            write_sequence(output_directory, header, sequence)


def write_sequence(output_directory, header, sequence):
    if not output_directory.exists():
        output_directory.mkdir()

    if '_' in header:
        return None

    print(f"{header}.fasta")
    output_file = Path(output_directory) / f"{header}.fasta"
    with open(output_file, 'w') as out:
        out.write(f'>{header}\n{sequence}\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Oncosplice Data Setup")
    parser.add_argument("basedir", help="The location of the directory where you want to establish a Grch38 database.")
    args = parser.parse_args()

    _ = retrieve_and_parse_ensembl_annotations(args.basedir)
    _ = load_grch38_fasta(args.basedir)

    import oncosplice
    init_path = os.path.join(oncosplice.__file__, '__init__.py')

    # Write the content
    with open(init_path, 'a') as file:
        file.write(f'database_path = \'{args.basedir}\'')

    print(f"Finished mounding database in {args.basedir}.")

