from Bio.Seq import Seq
from geney import find_end_codon, unload_pickle

from oncosplice import oncosplice_setup
from trash.mature_mRNA import mature_mRNA
from oncosplice.tis_utils import run_through_titer, build_titer_model

titer_model = build_titer_model()

class Protein(mature_mRNA):
    def __init__(self, transcript_start: int, transcript_end: int, rev: bool, chrm: str, donors: list, acceptors: list,
                 gene_name='undefined', transcript_id='undefined', transcript_type='undefined',
                 used_tis: int = None, used_tts: int = None, penetrance: float = 1.0, mutations=[]):

        mature_mRNA.__init__(self, transcript_start=transcript_start, transcript_end=transcript_end, rev=rev,
                             chrm=chrm, donors=donors, acceptors=acceptors, gene_name=gene_name,
                             transcript_id=transcript_id, transcript_type=transcript_type, penetrance=penetrance,
                             mutations=mutations)

        self.used_tis = used_tis
        self.used_tts = used_tts
        self.orf = ''
        self.protein = ''
        self.conservation_vector = []
        self.__generate_protein()
        self.__access_conservation_data()

    def __repr__(self):
        return 'Protein(transcript_id={tid})'.format(tid=self.transcript_id)

    def __len__(self):
        return len(self.protein)

    def __eq__(self, other):
        return self.protein == other.protein

    def __contains__(self, subvalue):
        if isinstance(subvalue, str):
            return subvalue in self.protein
        elif isinstance(subvalue, int):
            return subvalue in self.mature_indices
        else:
            print("Pass an integer to check against the span of the protein's coordinates or a string to "
                  "check against the protein sequence.")
            return False

    def __valid_protein(self):
        if self.used_tis in self.mature_indices:
            return True
        else:
            return False

    def __access_conservation_data(self):
        files = [file for file in oncosplice_setup['CONS_PATH'].glob(f"*{self.transcript_id.split('.')[0]}*")]
        if len(files) > 0:
            cons_data = unload_pickle(files[0])
            if cons_data['seq'] == self.protein:
                self.cons_seq = cons_data['seq']
                self.conservation_vector = cons_data['scores']
                self.cons_available = True
                return self

        self.cons_seq = self.protein
        self.conservation_vector = [1] * len(self.cons_seq)
        self.cons_available = False
        return self
    def __generate_protein(self):
        if self.used_tis not in self.mature_indices:
            self.used_tis, _, _, _, _ = run_through_titer(mut_seq=self.mature_mrna,
                                               mut_coords=self.mature_indices,
                                               ref_sc_coord=self.used_tis,
                                               ref_tts_coord=self.used_tts,
                                               ref_id=self.transcript_id,
                                               data_path='/tamir1/shaicohen1/share/titer-master',
                                               titer_model=titer_model,
                                               all_test_seqs={})

        if self.__valid_protein():
            rel_start = self.mature_indices.index(self.used_tis)
            self.used_tts = find_end_codon(self.mature_mrna[rel_start:], self.mature_indices[rel_start:])
            rel_end = self.mature_indices.index(self.used_tts)
            self.orf = self.mature_mrna[rel_start:rel_end + 1]
            self.protein = str(Seq(self.orf).translate())

        self.__access_conservation_data()


