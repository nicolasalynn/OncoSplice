


def generate_mut_variant(seq: str, indices: list, start_pos: int, end_pos: int, var_type: str,
                         ref: str, mut: str, suppress=False):

    if (var_type == 'SNP' and start_pos not in indices) or (var_type == 'INS' and (start_pos not in indices or start_pos + 1 not in indices)) or (var_type == 'DEL' and not any([v in indices for v in range(start_pos, start_pos + len(ref))])):
        real_indices = [v for v in indices if v > 0]
        if suppress:
            print(
                f"Mutation ({start_pos}:{end_pos}:{ref}:{mut}) not within context of sequence. ie, gene covers ({min(real_indices[0], real_indices[-1])}, {max(real_indices[0], real_indices[-1])}) while mutation occurs at {start_pos}")
        # print(indices)
        return seq, indices, False, [], 'Mutation not within context of sequence.'

    if var_type in ['SNP', 'DNP', 'TNP', 'ONP']:
        '''
        In case of SNP, all characters from start_pos to end_pos are altered.
        Thus, both start_pos and end_pos in seq are altered.
        '''

        assert seq[indices.index(
            start_pos)] == ref, f'Reference allele does not match position in SNP. {seq[indices.index(start_pos)]}, {start_pos}, {ref}'

        return seq[:indices.index(start_pos)] + mut + seq[indices.index(start_pos) + 1:], indices, True, [start_pos], ''


    elif var_type == 'INS':
        '''
        In case of INS, the inserted mut is between start_pos and end_pos.
        Thus, both start_pos and end_pos positions in seq are not altered.
        '''
        rel_pos = indices.index(start_pos)
        new_seq, new_indices = seq[:rel_pos + 1] + mut + seq[rel_pos + 1:], indices[:rel_pos + 1] + [
            indices[rel_pos] + v / 1000 for v in list(range(1, len(mut) + 1))] + indices[rel_pos + 1:]

        return new_seq, new_indices, True, [], ''

    elif var_type == 'INDEL':
        '''
        In case of INDEL, all characters from start_pos to end_pos are altered.
        Thus, both start_pos and end_pos in seq are altered.
        '''

        assert seq[indices.index(
            start_pos):indices.index(end_pos)+1] == ref, f'Reference allele does not match position in SNP. {seq[indices.index(start_pos):indices.index(end_pos)+1]}, {start_pos}, {ref}'

        rel_start, rel_end = indices.index(start_pos), indices.index(end_pos) + 1
        return seq[:rel_start] + mut + seq[rel_end:], indices[:rel_start] + [start_pos + i/1000 for i in range(1, len(mut))+1] + indices[rel_end:], True, [start_pos], ''


    elif var_type == 'DEL':
        '''
        In case of DEL, all charaters from start_pos to end_pos are removed.
        '''
        # new_indices, new_seq = indices, seq
        for sp, ref_val in list(zip(range(start_pos, start_pos+len(ref)), ref)):
            if ref_val not in indices:
                continue
            rel_pos = indices.index(sp)
            assert seq[rel_pos] == ref_val, f'Reference allele does not match position in DEL. {ref_val} not {seq[rel_pos]}'
            indices = indices[:rel_pos] + indices[rel_pos + 1:]
            seq = seq[:rel_pos] + seq[rel_pos + 1:]

        return seq, indices, True, list(range(start_pos, start_pos+len(ref))), ''

    else:
        print(f"generate_mut_variant: var_type {var_type} not supported!!")
        return None, None, False, [], f'var_type {var_type} not supported.'

class EpistaticSet:
    def __init__(self, epistatic_set):
        self.epistatic_id = epistatic_set
        self.variants = [Mutation(m) for m in self.epistatic_id.split('|')]
        self.gene = self.variants[0].gene
        self.chrom = self.variants[0].chrom.strip('chr')
        self.file_identifier = ','.join([v.file_identifier for v in self.variants])

    def __str__(self):
        return '|'.join([m.mut_id for m in self.variants])

    def __repr__(self):
        return '|'.join([m.mut_id for m in self.variants])


class Mutation:
    def __init__(self, mid):
        self.mut_id = mid
        gene, chrom, pos, ref, alt = mid.split(':')
        self.gene = gene
        self.chrom = chrom.strip('chr')
        self.start = int(pos)
        self.ref = ref
        self.alt = alt

        if self.ref != '-' and self.alt != '-' and (len(self.ref) > 1 or len(self.alt) > 1):
            self.vartype = 'INDEL'
        elif self.ref != '-' and self.alt != '-':
            self.vartype = 'SNP'
        elif self.ref == '-':
            self.vartype = 'INS'
        elif self.alt == '-':
            self.vartype = 'DEL'
        elif len(self.alt) > 1 and len(self.ref) > 1:
            self.vartype = 'INDEL'

        if self.vartype == 'DEL' or self.vartype == 'INDEL':
           self.end = self.start + len(self.ref) - 1
        else:
            self.end = self.start

        self.file_identifier = f'{self.gene}_{self.chrom}_{self.start}_{self.ref}_{self.alt}'

    def __str__(self):
        return self.mut_id

    def __repr__(self):
        return self.mut_id
