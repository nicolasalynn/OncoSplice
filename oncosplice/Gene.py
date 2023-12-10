from copy import copy
from geney import *
# from geney import get_correct_gene_file, reverse_complement, pull_fasta_seq_endpoints
from Bio.Seq import Seq
# import json
from oncosplice.variant_utils import generate_mut_variant, Mutation, find_new_tis, find_new_tts
from pathlib import Path
from oncosplice import oncosplice_setup
file = Path('/Users/nl/Documents/phd/data/ensembl/mRNAs/protein_coding/mrnas_ENSG00000006283.18_CACNA1G.json')

class Gene:

    def __init__(self, gene_name=None, file=None, dict_data=None, target_directory=oncosplice_setup['MRNA_PATH']):
        self.gene_name = gene_name
        self.gene_id = ''
        self.rev = None
        self.chrm = ''
        self.gene_start = 0
        self.gene_end = 0
        self.transcripts = {}
        if gene_name:
            file = get_correct_gene_file(gene_name, target_directory=target_directory)
            if file.exists():
                self.load_from_file(file)
        elif dict_data:
            self.load_from_dict(dict_data)
        elif file:
            self.load_from_file(file)

    def __repr__(self):
        return 'Gene(gene_name={gname})'.format(gname=self.gene_name)

    def __len__(self):
        return len(self.transcripts)

    def __str__(self):
        return '{gname}, {ntranscripts} transcripts'.format(gname=self.gene_name, ntranscripts=self.__len__())

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __getitem__(self, index):
        return Transcript(list(self.transcripts.values())[index])

    def load_from_file(self, file_name):
        if not file_name.exists():
            raise ValueError(f"'{str(file_name)}' does not exist.")

        data = json.load(open(file_name))
        self.load_from_dict(dict_data=data)
        return self

    def load_from_dict(self, dict_data=None):
        for k, v in dict_data.items():
            setattr(self, k, v)
        return self

    def write_annotation_file(self, file_name):
        with open(file_name, 'w') as out:
            json.dump(out, self.__dict__)

    def generate_transcript(self, tid):
        return Transcript(self.transcripts[tid])

    @property
    def primary_transcript(self):
        temp = [k for k, v in self.transcripts.items() if v.get('primary_transcript', False) or ('Ensembl_canonical' in v['tag'] and ('CCDS' in v['tag']))]
        return self.generate_transcript(temp[0])

    def transcript(self, tid):
        temp = [k for k, v in self.transcripts.items() if tid in k]
        return self.generate_transcript(temp[0])

class Transcript:
    def __init__(self, d=None):
        self.transcript_id = None
        self.transcript_start = None
        self.transcript_end = None
        self.transcript_type = None
        self.acceptors, self.donors = [], []
        self.TIS, self.TTS = None, None
        self.transcript_seq, self.transcript_indices = '', []
        self.rev = None
        self.chrm = ''
        self.pre_mrna = ''
        self.orf = ''
        self.protein = ''
        if d:
            self.load_from_dict(d)

    def __repr__(self):
        return 'pre_mRNA(transcript_id={tid})'.format(tid=self.transcript_id)

    def __len__(self):
        return len(self.transcript_seq)

    def __str__(self):
        return 'Transcript {tid}, Transcript Type: ' \
               '{protein_coding}'.format(
                tid=self.transcript_id, protein_coding=self.transcript_type)

    def __eq__(self, other):
        return self.transcript_seq == other.transcript_seq

    def __contains__(self, subvalue):
        if isinstance(subvalue, str):
            return subvalue in self.transcript_seq
        elif isinstance(subvalue, int):
            return subvalue in self.transcript_indices
        else:
            print(
                "Pass an integer to check against the span of the gene's coordinates or a string to check against the "
                "pre-mRNA sequence.")
            return False

    def __copy__(self, other):
        return copy(self)

    @property
    def constructor(self):
        core_attributes = ['transcript_id', 'transcript_start', 'transcript_end', 'transcript_type', 'acceptors', 'donors', 'TIS', 'TTS', 'rev', 'chrm']
        return {k: v for k, v in self.__dict__.items() if k in core_attributes}

    def load_from_dict(self, data):
        for k, v in data.items():
            setattr(self, k, v)
        self.__arrange_boundaries()
        self.generate_mature_mrna(inplace=True)
        return self

    @property
    def exons(self):
        return list(zip(self.acceptors, self.donors))

    def set_exons(self, boundaries):
        self.acceptors, self.donors = boundaries['acceptors'], boundaries['donors']
        self.__arrange_boundaries()
        return self

    @property
    def introns(self):
        return list(zip([v for v in self.donors if v != self.transcript_end], [v for v in self.acceptors if v != self.transcript_start]))

    def __exon_coverage_check(self):
        if sum([abs(a-b) + 1 for a, b in self.exons]) == len(self):
            return True
        else:
            return False
    @property
    def exons_pos(self):
        temp = self.exons
        if self.rev:
            temp = [(b, a) for a, b in temp[::-1]]
        return temp
    @property
    def mrna_indices(self):
        temp = [lst for lsts in [list(range(a, b+1)) for a, b in self.exons_pos] for lst in lsts]
        return sorted(temp, reverse=self.rev)

    def __arrange_boundaries(self):
        self.acceptors.append(self.transcript_start)
        self.donors.append(self.transcript_end)
        self.acceptors = list(set(self.acceptors))
        self.donors = list(set(self.donors))
        self.acceptors.sort(reverse=self.rev)
        self.donors.sort(reverse=self.rev)
        return self

    def positive_strand(self):
        if self.rev:
            return reverse_complement(self.transcript_seq)
        else:
            return self.transcript_seq

    def __pos2sense(self, mrna, indices):
        if self.rev:
            mrna = reverse_complement(mrna)
            indices = indices[::-1]
        return mrna, indices

    def pull_pre_mrna_pos(self):
        if self.rev:
            return pull_fasta_seq_endpoints(self.chrm, self.transcript_end,
                                                                   self.transcript_start)
        else:
            return pull_fasta_seq_endpoints(self.chrm, self.transcript_start,
                                                                   self.transcript_end)

    def generate_pre_mrna_pos(self, mutations=[]):
        seq, indices = self.pull_pre_mrna_pos()
        for mutation in mutations:
            mutation = Mutation(mutation)
            seq, indices, _, _ = generate_mut_variant(seq, indices, mut=mutation)

        self.pre_mrna, _ = self.__pos2sense(seq, indices)
        return seq, indices

    def generate_pre_mrna(self, mutations=[], inplace=True):
        pre_mrna, pre_indices = self.__pos2sense(*self.generate_pre_mrna_pos(mutations))
        self.pre_mrna = pre_mrna
        if inplace:
            return self
        return pre_mrna

    def generate_mature_mrna_pos(self, mutations=[]):
        mature_mrna, mature_indices = '', []
        pre_seq, pre_indices = self.generate_pre_mrna_pos(mutations)
        for i, j in self.exons_pos:
            rel_start, rel_end = pre_indices.index(i), pre_indices.index(j)
            mature_mrna += pre_seq[rel_start:rel_end + 1]
            mature_indices.extend(pre_indices[rel_start:rel_end + 1])
        return mature_mrna, mature_indices

    def generate_mature_mrna(self, mutations=[], inplace=True):
        if inplace:
            self.transcript_seq, self.transcript_indices = self.__pos2sense(*self.generate_mature_mrna_pos(mutations))
            return self
        return self.__pos2sense(*self.generate_mature_mrna_pos(mutations))

    def generate_protein(self, inplace=True):
        rel_start = self.transcript_indices.index(self.TIS)
        rel_end = self.transcript_indices.index(self.TTS)
        orf = self.transcript_seq[rel_start:rel_end + 1 + 3]
        protein = str(Seq(orf).translate())
        if inplace:
            self.orf = orf
            self.protein = protein
            return self
        return protein

    def generate_translational_boundaries(self):
        if self.TIS not in self.transcript_indices:
            if 'ATG' in self.transcript_seq:
                self.TIS = self.transcript_indices[self.transcript_seq.index('ATG')]
            else:
                self.TIS = self.transcript_indices[0]
            # self.TIS = find_new_tis(self.transcript_seq, self.transcript_indices, self.TIS, self.TTS)
        self.TTS = find_new_tts(self.transcript_seq, self.transcript_indices, self.TIS)
        return self

