import copy

from geney import pull_fasta_seq_endpoints, reverse_complement

from oncosplice.variant_utils import generate_mut_variant, Mutation

class pre_mRNA:

    def __init__(self, transcript_start: int, transcript_end: int, rev: bool, chrm: str, gene_name='undefined',
                 transcript_id='undefined', transcript_type='undefined', mutations=[]):
        self.transcript_start = transcript_start
        self.transcript_end = transcript_end
        self.rev = rev
        self.chrm = chrm
        self.gene_name = gene_name
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type

        self.transcript_upper_bound = max([self.transcript_start, self.transcript_end])
        self.transcript_lower_bound = min([self.transcript_start, self.transcript_end])

        # Features of greater classes
        self.pre_mrna, self.pre_indices = '', []
        self.applied_mutations = []
        self.strand_state_rev = False
        self.__generate_pre_mrna()
        if mutations:
            self.apply_mutations(mutations)

    def __repr__(self):
        return 'pre_mRNA(transcript_id={tid})'.format(tid=self.transcript_id)

    def __len__(self):
        return self.transcript_upper_bound - self.transcript_lower_bound

    def __str__(self):
        return 'Gene {gname}, Transcript {tid}, Applied Mutations: {mutations}, Transcript Type: ' \
               '{protein_coding}'.format(
                gname=self.gene_name, tid=self.transcript_id, mutations=bool(self.applied_mutations),
                protein_coding=self.transcript_type)

    def __eq__(self, other):
        return self.pre_mrna == other.pre_mrna

    def __contains__(self, subvalue):
        if isinstance(subvalue, str):
            return subvalue in self.pre_mrna
        elif isinstance(subvalue, int):
            return subvalue in self.pre_indices
        else:
            print(
                "Pass an integer to check against the span of the gene's coordinates or a string to check against the "
                "pre-mRNA sequence.")
            return False

    def __copy__(self, other):
        return copy.deepcopy(self)

    def __valid_pre_mrna(self):
        # Description: in order for a pre-mRNA to be valid we need at least least one exon start and one exon end.
        if isinstance(self.transcript_start, int) and isinstance(self.transcript_end, int):
            return True
        else:
            print("Invalid pre mRNA. pre_mRNA requires transcript start and end genomic coordiantes.")
            return False

    def __to_positive_strand(self):
        if self.strand_state_rev:  #self.rev and self.pre_indices[0] > self.pre_indices[-1]:
            self.pre_mrna, self.pre_indices = reverse_complement(self.pre_mrna), self.pre_indices[::-1]
            self.strand_state_rev = False
        return self
    def __to_true_strand(self):
        if self.rev and not self.strand_state_rev: #self.pre_indices[0] < self.pre_indices[-1]:
            self.pre_mrna, self.pre_indices = reverse_complement(self.pre_mrna), self.pre_indices[::-1]
            self.strand_state_rev = True
        return self

    def __generate_pre_mrna(self):
        if self.__valid_pre_mrna():
            self.applied_mutations = []
            self.pre_mrna, self.pre_indices = pull_fasta_seq_endpoints(self.chrm, self.transcript_lower_bound,
                                                                       self.transcript_upper_bound)
            self.__to_true_strand()
        return self

    def apply_mutations(self, ms):
        ms = ms if isinstance(ms, list) else [ms]
        self.__to_positive_strand()

        mrna_seq, mrna_indices = self.pre_mrna, self.pre_indices
        for mutation in ms:
            mutation = Mutation(mutation)
            mrna_seq, mrna_indices, successfully_applied = \
                generate_mut_variant(
                    mrna_seq,
                    mrna_indices,
                    mut=mutation)

            if successfully_applied:
                self.applied_mutations.append(mutation.mut_id)

        self.pre_mrna, self.pre_indices = mrna_seq, mrna_indices
        self.__to_true_strand()
        return self

    def is_mutated(self):
        return True if self.applied_mutations else False

