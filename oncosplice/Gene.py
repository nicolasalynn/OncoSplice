
import copy
from geney import dump_json, unload_json

from oncosplice.pre_mRNA import pre_mRNA
from oncosplice.mature_mRNA import mature_mRNA
from oncosplice.Protein import Protein

class EmptyGene:
    def __init__(self):
        self.gene_name = ''
        self.gene_id = ''
        self.rev = None
        self.chrm = ''
        self.gene_start = 0
        self.gene_end = 0
        self.transcripts = {}
        self.mutations = []

    def __repr__(self):
        return 'Gene(gene_name={gname})'.format(gname=self.gene_name)

    def __len__(self):
        return len(self.transcripts)

    def __getitem__(self, position):
        return list(self.transcripts.keys())[position]

    def __str__(self):
        return '{gname}, {ntranscripts} transcripts'.format(gname=self.gene_name, ntranscripts=self.__len__())

    def __copy__(self):
        return copy.deepcopy(self)

    def write_annotation_file(self, file_name):
        dump_json(file_name, self.__dict__)

    def add_transcript(self, transcript):
        tid = transcript.get('transcript_id', '')
        if not tid:
            print('Cannot add a transcript without a meaningful transcript_id.')
        elif tid in self.transcripts.keys():
            print('Re-writing annotations for transcript {tid}.')

        self.transcripts.update({tid: transcript})

    def develop_pre_mrna(self, tid):
        return pre_mRNA(gene_name=self.gene_name, rev=self.rev, chrm=self.chrm, transcript_id=tid,
                        transcript_start=self.transcripts[tid]['transcript_start'],
                        transcript_end=self.transcripts[tid]['transcript_end'],
                        transcript_type=self.transcripts[tid]['transcript_type'])

    def develop_mature_mrna(self, tid):
        return mature_mRNA(gene_name=self.gene_name,
                           rev=self.rev,
                           chrm=self.chrm,
                           transcript_id=tid,
                           transcript_start=self.transcripts[tid]['transcript_start'],
                           transcript_end=self.transcripts[tid]['transcript_end'],
                           transcript_type=self.transcripts[tid]['transcript_type'],
                           donors=self.transcripts[tid].get('donors', []),
                           acceptors=self.transcripts[tid].get('acceptors', []))

    def develop_protein(self, tid):
        return Protein(gene_name=self.gene_name,
                       rev=self.rev,
                       chrm=self.chrm,
                       transcript_id=tid,
                       transcript_start=self.transcripts[tid]['transcript_start'],
                       transcript_end=self.transcripts[tid]['transcript_end'],
                       transcript_type=self.transcripts[tid]['transcript_type'],
                       donors=self.transcripts[tid].get('donors', []),
                       acceptors=self.transcripts[tid].get('acceptors', []),
                       used_tis=self.transcripts[tid]['TIS'],
                       used_tts=self.transcripts[tid]['TTS'])

    def develop_proteome(self, experimental=False):
        proteome = {}
        for tid in self.transcripts.keys():
            if 'TIS' in self.transcripts[tid].keys():
                proteome[tid] = self.develop_protein(tid)
                proteome[tid].generate_protein(mutations=self.mutations, regenerate_mature_mrna=True, regenerate_pre_mrna=True, experimental=experimental)
        return proteome

    def ntranscripts(self, filter=None):
        if filter:
            count = 0
            for t in self.transcripts.values():
                if t['transcript_type'] == filter:
                    count += 1
            return count

        else:
            return len(self.transcripts)

    def nexons(self):
        acceptors, donors = [], []
        for t in self.transcripts.values():
            acceptors.append(t['transcript_start'])
            acceptors.extend(t['acceptors'])
            donors.extend(t['donors'])
            donors.append(t['transcript_end'])

        nunique_acceptors = len(set(acceptors))
        nunique_donors = len(set(donors))
        exons = list(zip(acceptors, donors))
        nexons = len(set(exons))
        return {'exon_count': nexons, 'acceptor_count': nunique_acceptors, 'donor_count': nunique_donors}

class AnnotatedGene(EmptyGene):
    def __init__(self, file):
        # self.file = file
        super().__init__()
        if not file.exists():
            raise ValueError(f"Gene {file} annotations not available.")

        self.data = unload_json(file)
        self.gene_name = self.data['gene_name']
        if not all(['gene_name' in self.data.keys(), 'rev' in self.data.keys(), 'chrm' in self.data.keys(),
                    'gene_start' in self.data.keys(), 'gene_end' in self.data.keys()]):
            raise ValueError(f"Gene {self.gene_name} annotations not complete. Please ensure annotations include gene "
                             f"name, rev, chrm, gene_start, and gene_end")
        self.gene_name = self.data['gene_name']
        self.rev = self.data['rev']
        self.chrm = self.data['chrm']
        self.gene_start = self.data['gene_start']
        self.gene_end = self.data['gene_end']

        self.gene_id = self.data.get('gene_id', 'undefinied')
        self.tag = self.data.get('tag', 'undefined')
        self.transcripts = self.data.get('transcripts', {})
        self.transcripts = {k: v for k, v in self.transcripts.items() if 'donors' in v.keys() and 'acceptors' in v.keys() and 'TIS' in v.keys() and 'TTS' in v.keys()}


    def __copy__(self):
        return copy.deepcopy(self)

    def create_gene_isoform(self, mut_ids, aberrant_splicing=None):

        variant_gene = self.__copy__()
        variant_gene.mutations = mut_ids.split('|')
        for tid in self.transcripts.keys():
            ref_transcript = self.develop_mature_mrna(tid)

            for tid_counter, new_blueprints in enumerate(
                    ref_transcript.develop_aberrant_splicing(
                                              aberrant_splicing=aberrant_splicing,
                                              rev=self.rev)):

                var_transcript = copy.deepcopy(self.transcripts[tid])
                var_transcript['protein_seq'] = var_transcript['transcript_seq'] = ''
                var_transcript['transcript_id'] += f'-{tid_counter}'
                var_transcript['donors'] = new_blueprints['donors']
                var_transcript['acceptors'] = new_blueprints['acceptors']
                var_transcript['penetrance_weight'] = new_blueprints['path_weight']
                variant_gene.add_transcript(var_transcript)

        return variant_gene




