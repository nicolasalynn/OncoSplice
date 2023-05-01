import copy

from geney import pull_fasta_seq_endpoints, reverse_complement

from oncosplice.variant_utils import generate_mut_variant, Mutation

class pre_mRNA:
    """ Description: a class dedicated to representing a pre_mrna sequence.
        Parameters:
            - rev (required): whether the transcript is encoded on the positive or negative strand
            - chrm: the chromosome (required) (str) that the transcript/gene is encoded in
            - transcript_start (required) (relative to the sense strand)
            - transcript_end (required) (relative to the sense strand)
            - gene_name (optional): the name of the gene that the transcript corresponds to
            - transcript_type (optional): a description of the role of the transcript (eg. protein coding)
            - transcript_id (optional): the transcript's identifying id with version

        Output: an instance of this object stores all the information required to construct the transcript sequence and
                chromososomal coverage. Depending on whether the transcript is encoded on the + or - strand, the correct
                start and end positions of the transcript are acquired and the pre-mRNA can be obtained using reference
                build genome fasta files.

                The two items that are obtained using the data provided are a sequence (with mutations) and chromsomal
                coverage as a list (positioned so that can be overlaid on sequence). "

        Functionality: - __valid_pre_mRNA: checks to make sure that the both the transcript start and end are integer
                         values and makes sure that the chromosome defined is valid and has associated fasta data
                         to_positive_strand: converts inplace the mRNA sequence and genomic indices to the positive
                                             strand
                         to_true_strand: converts inplace the mRNA sequence and genomic indices to the true encoded
                                         strand
                         generate_pre_mRNA: if __valid_pre_mRNA is satisfied, this funciton will generate the pre_mRNA
                                            sequence and corresponding indices and will apply the passed mutations to
                                            the sequence and coordinates. Once the mutations are applied, if necessary,
                                            the sequence is then converted to the sense strand.
                         apply_mutations: takes a list of mutation ids as an argument and will

    """

    def __init__(self, transcript_start: int, transcript_end: int, rev: bool, chrm: str, gene_name='undefined',
                 transcript_id='undefined', transcript_type='undefined'):
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
        self.generation_report = ''

        # self.generate_pre_mRNA()

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

    def to_positive_strand(self):
        if self.rev and self.pre_indices[0] > self.pre_indices[-1]:
            self.pre_mrna, self.pre_indices = reverse_complement(self.pre_mrna), self.pre_indices[::-1]

    def to_true_strand(self):
        if self.rev and self.pre_indices[0] < self.pre_indices[-1]:
            self.pre_mrna, self.pre_indices = reverse_complement(self.pre_mrna), self.pre_indices[::-1]

    def generate_pre_mrna(self, mutations=None) -> None:
        # Generates the pre-mRNA sequence. chrm_path must point to a
        # chromosome fasta file with equal number of bps per row
        if self.__valid_pre_mrna():
            self.generation_report = ''
            self.applied_mutations = []
            self.pre_mrna, self.pre_indices = pull_fasta_seq_endpoints(self.chrm, self.transcript_lower_bound,
                                                                       self.transcript_upper_bound)

            if mutations:
                self.apply_mutations(mutations)

            self.to_true_strand()

    def apply_mutations(self, ms):
        if ms:
            ms = ms if isinstance(ms, list) else [ms]
            if ms:
                self.to_positive_strand()

            mrna_seq, mrna_indices = self.pre_mrna, self.pre_indices
            for m in ms:
                expanded_mut = Mutation(m)
                mrna_seq, mrna_indices, successfully_applied, affected_indices, var_inclusion_report = \
                    generate_mut_variant(
                        mrna_seq,
                        mrna_indices,
                        start_pos=expanded_mut['Start_Position'],
                        end_pos=expanded_mut['End_Position'],
                        var_type=expanded_mut['Variant_Type'],
                        ref=expanded_mut['Reference_Allele'],
                        mut=expanded_mut['Tumor_Seq_Allele2'])

                if successfully_applied:
                    self.applied_mutations.append(m)
                    self.generation_report += f'Mutation {m} incorporated.'

                else:
                    self.generation_report += f'Mutation {m} not incorporated.)'
                    if self.rev:
                        if expanded_mut['Start_Position'] > self.transcript_start:
                            distance_before = expanded_mut['Start_Position'] - self.transcript_start
                            self.generation_report += f'Mutation occurs {distance_before} ' \
                                                      f'nucleotides before the transcript start site.'
                        elif expanded_mut['Start_Position'] < self.transcript_end:
                            distance_after = self.transcript_end - expanded_mut['Start_Position']
                            self.generation_report += f'Mutation occurs {distance_after} ' \
                                                      f'nucleotides after the transcript end site.'
                    else:
                        if expanded_mut['Start_Position'] < self.transcript_start:
                            distance_before = self.transcript_start - expanded_mut['Start_Position']
                            self.generation_report += f'Mutation occurs {distance_before} ' \
                                                      f'nucleotides before the transcript start site.'
                        elif expanded_mut['Start_Position'] > self.transcript_end:
                            distance_after = expanded_mut['Start_Position'] - self.transcript_end
                            self.generation_report += f'Mutation occurs {distance_after} ' \
                                                      f'nucleotides after the transcript end site.'

                if self.transcript_start in affected_indices:
                    self.generation_report += f'Mutation {m} affects the Transcription Start Site.\n'

                if self.transcript_end in affected_indices:
                    self.generation_report += f'Mutation {m} affects the Transcription End Site.\n'

            self.pre_mrna, self.pre_indices = mrna_seq, mrna_indices
            # self.to_true_strand()

    def is_mutated(self):
        return True if self.applied_mutations else False
