from pathlib import Path
import os
from gtfparse import read_gtf
from geney import dump_json, dump_pickle, unload_pickle
import pandas as pd
from tqdm import tqdm
import requests
import argparse
from sh import gunzip


def download_and_ungzip(external_url, local_path):
    local_file = Path(external_url).name
    local_file_path = Path(local_path) / local_file
    try:
        response = requests.get(external_url, stream=True)
        response.raise_for_status()  # Raises a HTTPError if the HTTP request returned an unsuccessful status code
        with open(local_file_path, 'wb') as f:
            f.write(response.content)

    except Exception as e:
        print(f"Error during download: {e}")

    gunzip(str(local_file_path))
    local_file_path = Path(local_file_path.as_posix().rstrip('.gz'))
    return local_file_path


def process_transcript(transcript_df, rev, cons_data):
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

    if transcript.transcript_id in cons_data:
        data.update({'cons_available': True, 'cons_vector': cons_data[transcript.transcript_id]['scores']})

    else:
        data.update({'cons_available': False})

    return data


def retrieve_and_parse_ensembl_annotations(local_path, annotations_file, gtex_file, cons_file, valid_biotypes=('protein_coding')):
    cons_data = unload_pickle(cons_file)

    gtex_df = pd.read_csv(gtex_file, delimiter='\t', header=2)
    gtex_df.Name = gtex_df.apply(lambda row: row.Name.split('.')[0], axis=1)
    gtex_df = gtex_df.set_index('Name').drop(columns=['Description'])

    annotations = read_gtf(annotations_file)

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

        biotype_path = local_path / biotype[0]
        if not biotype_path.exists():
            biotype_path.mkdir()

        file_name = biotype_path / f'mrnas_{gene_id}_{gene_attribute.gene_name}.pkl'
        if file_name.exists():
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
            'transcripts': {transcript_id: process_transcript(transcript_df, rev, cons_data) for transcript_id, transcript_df in
                            gene_df.groupby('transcript_id') if transcript_id},
            'tissue_expression': gtex_df.loc[gene_id].squeeze().to_dict() if gene_id in gtex_df.index else {},
        }

        json_data['transcripts'] = {tid: v for tid, v in json_data['transcripts'].items() if v is not None}
        if not json_data['transcripts']:
            continue

        if gene_attribute.gene_name == '' or gene_id == '':
            continue

        dump_pickle(file_name, json_data)


def split_fasta(input_file, output_directory):
    """
    Splits a gzipped FASTA file into multiple files, each containing a single sequence.
    File names are derived from the sequence identifier in the FASTA header.

    :param input_file: Path to the input gzipped FASTA file.
    :param output_directory: Directory where output FASTA files will be saved.
    """

    with open(input_file, 'r') as file:
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


def main():
    config_dir = Path(os.path.join(os.path.expanduser('~'), '.oncosplice_setup'))
    if config_dir.exists():
        for file in config_dir.glob('*'):
            file.unlink()
        config_dir.rmdir()
    config_dir.mkdir()

    parser = argparse.ArgumentParser(description="Conservation file location")
    parser.add_argument("consfile", help="The location of the conservation data file.", required=True)
    parser.add_argument("basepath", help="The location of the data we are mounting.", required=True)
    parser.add_argument("splicepath", help="The location of the data we are mounting.", required=False, default=None)
    args = parser.parse_args()

    config_file = config_dir / 'config.json'
    config_paths = {
        'CHROM_SOURCE': os.path.join(args.basepath, 'chromosomes'),
        'MRNA_PATH': os.path.join(args.basepath, 'annotations'),
        'MISSPLICING_PATH': args.splicepath if args.splicepath is not None else os.path.join(args.basepath, 'missplicing')
    }
    dump_json(config_file, config_paths)

    cons_file = Path(args.consfile)
    assert cons_file.exists(), f"{cons_file} does not exist. Please provide a path to the conservation data."

    base_path = Path(args.basepath)
    if base_path.exists() and len(os.listdir(base_path)) > 0:
        raise FileExistsError(f"Directory {base_path} not empty.")
        # pass
    #
    elif not base_path.exists():
        print(f"Initializing data folder at {base_path}.")
        base_path.mkdir()

    gtex_url = 'https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'
    gtex_file = download_and_ungzip(gtex_url, base_path)

    fasta_url = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
    fasta_file = download_and_ungzip(fasta_url, base_path)
    fasta_build_path = base_path / f'chromosomes'
    fasta_build_path.mkdir()
    split_fasta(fasta_file, fasta_build_path)

    ensembl_url = 'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz'
    ensembl_file = download_and_ungzip(ensembl_url, base_path)
    ensembl_annotation_path = base_path / f'annotations'
    ensembl_annotation_path.mkdir()
    retrieve_and_parse_ensembl_annotations(ensembl_annotation_path, ensembl_file, gtex_file, cons_file)

    fasta_file.unlink()
    gtex_file.unlink()
    ensembl_file.unlink()
    print(f"Finished mounding database in {args.basepath}.")

if __name__ == '__main__':
    main()