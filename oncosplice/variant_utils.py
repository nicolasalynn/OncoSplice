
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
        self.file_identifier = self.mut_id.replace(':', '_')

        gene, chrom, pos, ref, alt = mid.split(':')
        self.gene = gene
        self.chrom = chrom.strip('chr')
        self.start = int(pos)
        self.ref = ref if ref != '-' else ''
        self.alt = alt if alt != '-' else ''

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

    def __str__(self):
        return self.mut_id

    def __repr__(self):
        return self.mut_id


def generate_mut_variant(seq: str, indices: list, mut: Mutation, suppress=False):
    offset = 1 if not mut.ref else 0
    check_indices = list(range(mut.start, mut.start + len((mut.ref) + offset)))
    if any([m not in indices for m in check_indices]):
        print(f"Mutation {mut} not in indices: {min(indices)} - {max(indices)}.")
        raise IndexError
        return seq, indices
    # if (mut.vartype == 'SNP' and mut.start not in indices) or (mut.vartype == 'INS' and (mut.start not in indices or mut.start + 1 not in indices)) or (mut.vartype == 'DEL' and not any([v in indices for v in range(mut.start, mut.start + len(mut.ref))])):
    #     real_indices = [v for v in indices if v > 0]
    #     if suppress:
    #         print(
    #             f"Mutation ({mut}) not within context of sequence. ie, gene covers ({min(real_indices[0], real_indices[-1])}, {max(real_indices[0], real_indices[-1])})")
    #     return seq, indices

    # if mut.vartype in ['SNP', 'DNP', 'TNP', 'ONP', 'DEL', 'IN']:
    '''
    In case of SNP, all characters from start_pos to end_pos are altered.
    Thus, both start_pos and end_pos in seq are altered.
    '''

    assert seq[indices.index(
        mut.start)] == mut.ref, f'Reference allele does not match position in SNP. {seq[indices.index(mut.start)]}, {mut.start}, {mut.ref}'

    new_seq = seq[:indices.index(mut.start) + offset] + mut.alt + seq[indices.index(mut.start) + len(mut.ref) + offset:]

    if len(mut.ref) == len(mut.alt) > 0:
        new_indices = list(range(mut.start, mut.start + len(mut.ref)))
    else:
        new_indices = [indices[indices.index(mut.start)] + v / 1000 for v in list(range(len(1, mut.alt+1)))]

    new_indices = indices[:indices.index(mut.start) + offset] + new_indices + indices[indices.index(mut.start) + len(mut.ref) + offset:]

    assert len(new_seq) == len(new_indices), f'Error in variant modification: {mut}, {len(new_seq)}, {len(new_indices)}'
    return new_seq, new_indices


    # elif mut.vartype == 'INS':
    #     '''
    #     In case of INS, the inserted mut is between start_pos and end_pos.
    #     Thus, both start_pos and end_pos positions in seq are not altered.
    #     '''
    #     rel_pos = indices.index(start_pos)
    #     new_seq, new_indices = seq[:rel_pos + 1] + mut + seq[rel_pos + 1:], indices[:rel_pos + 1] + [
    #         indices[rel_pos] + v / 1000 for v in list(range(1, len(mut) + 1))] + indices[rel_pos + 1:]
    #
    #     return new_seq, new_indices, True, [], ''

    # elif var_type == 'INDEL':
    #     '''
    #     In case of INDEL, all characters from start_pos to end_pos are altered.
    #     Thus, both start_pos and end_pos in seq are altered.
    #     '''
    #
    #     # assert seq[indices.index(
    #     #     start_pos):indices.index(end_pos)] == ref, f'Reference allele does not match position in INDEL. {seq[indices.index(start_pos):indices.index(end_pos)]}, {start_pos}, {ref}'
    #     #
    #     # rel_start, rel_end = indices.index(start_pos), indices.index(end_pos)
    #     # return seq[:rel_start] + mut + seq[rel_end+1:], indices[:rel_start] + [start_pos + i/1000 for i in list(range(1, len(mut)+1))] + indices[rel_end:], True, [start_pos], ''
    #
    #     for sp, ref_val in list(zip(range(start_pos, start_pos+len(ref)), ref)):
    #         # if ref_val not in indices:
    #         #
    #         #     continue
    #         assert ref_val in indices, f'Ref index not in indices, {ref_val}, {indices}'
    #         rel_pos = indices.index(sp)
    #         assert seq[rel_pos] == ref_val, f'Reference allele does not match position in DEL. {ref_val} not {seq[rel_pos]}'
    #         indices = indices[:rel_pos] + indices[rel_pos + 1:]
    #         seq = seq[:rel_pos] + seq[rel_pos + 1:]
    #     indices = indices[:indices.index(start_pos)]
    #     seq = seq[:]
    #
    # elif var_type == 'DEL':
    #     '''
    #     In case of DEL, all charaters from start_pos to end_pos are removed.
    #     '''
    #     # new_indices, new_seq = indices, seq
    #     for sp, ref_val in list(zip(range(start_pos, start_pos+len(ref)), ref)):
    #         # if ref_val not in indices:
    #         #     continue
    #         assert ref_val in indices, f'Ref index not in indices, {ref_val}, {indices}'
    #         rel_pos = indices.index(sp)
    #         assert seq[rel_pos] == ref_val, f'Reference allele does not match position in DEL. {ref_val} not {seq[rel_pos]}'
    #         indices = indices[:rel_pos] + indices[rel_pos + 1:]
    #         seq = seq[:rel_pos] + seq[rel_pos + 1:]
    #
    #     return seq, indices, True, list(range(start_pos, start_pos+len(ref))), ''
    #
    # else:
    #     print(f"generate_mut_variant: var_type {var_type} not supported!!")
    #     return None, None, False, [], f'var_type {var_type} not supported.'

