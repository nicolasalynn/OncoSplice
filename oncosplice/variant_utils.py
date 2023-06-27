from geney import is_monotonic
class EpistaticSet:
    def __init__(self, epistatic_set):
        epistatic_set = {int(m.split(':')[2]): m for m in epistatic_set.split('|')}
        epistatic_set = {k: v for k, v in sorted(epistatic_set.items())}
        epistatic_set = '|'.join(epistatic_set.values())
        self.mut_id = epistatic_set
        self.variants = [Mutation(m) for m in self.mut_id.split('|')]
        self.start = self.variants[0].start
        self.positions = [v.start for v in self.variants]
        self.ref = ','.join([m.ref for m in self.variants])
        self.alt = ','.join([m.alt for m in self.variants])
        self.gene = self.variants[0].gene
        self.chrom = self.variants[0].chrom.strip('chr')
        self.file_identifier = f'{self.gene}_{self.chrom}' + '_' + '_'.join([v.file_identifier_short for v in self.variants])

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

        self.file_identifier = self.mut_id.replace(':', '_')
        self.file_identifier_short = f'{self.start}_{ref}_{alt}'

        self.ref = ref if ref != '-' else ''
        self.alt = alt if alt != '-' else ''

        if len(self.ref) == len(self.alt) == 1:
            self.vartype = 'SNP'
        elif len(self.ref) == len(self.alt) > 1:
            self.vartype = 'SUB'
        elif self.ref and not self.alt:
            self.vartype = 'DEL'
        elif self.alt and not self.ref:
            self.vartype = 'INS'
        else:
            self.vartype = 'INDEL'

    def __str__(self):
        return self.mut_id

    def __repr__(self):
        return self.mut_id


def generate_mut_variant(seq: str, indices: list, mut: Mutation):
    offset = 1 if not mut.ref else 0

    check_indices = list(range(mut.start, mut.start + len(mut.ref) + offset))
    check1 = all([m in indices for m in check_indices])
    if not check1:
        print(f"Mutation {mut} not within transcript bounds: {min(indices)} - {max(indices)}.")
        return seq, indices, False, False

    rel_start, rel_end = indices.index(mut.start)+offset, indices.index(mut.start)+offset+len(mut.ref)
    acquired_seq = seq[rel_start:rel_end]
    check2 = acquired_seq == mut.ref
    if not check2:
        print(f'Reference allele does not match genome_build allele. {acquired_seq}, {mut.ref}, {mut.start}')
        consensus_allele = False
    else:
        consensus_allele = True
    if len(mut.ref) == len(mut.alt) > 0:
        temp_indices = list(range(mut.start, mut.start + len(mut.ref)))
    else:
        temp_indices = [indices[indices.index(mut.start)] + v / 1000 for v in list(range(1, len(mut.alt)+1))]

    print(f'mutation: {mut}')
    print(f"New Indiced: {temp_indices}")
    new_indices = indices[:rel_start] + temp_indices + indices[rel_end:]
    new_seq = seq[:rel_start] + mut.alt + seq[rel_end:]

    assert len(new_seq) == len(new_indices), f'Error in variant modification: {mut}, {len(new_seq)}, {len(new_indices)}'
    assert is_monotonic(list(filter((-1).__ne__, new_indices))), f'Mut indices are not monotonic. {new_indices}'
    return new_seq, new_indices, True, consensus_allele



# def generate_mut_variant(seq: str, indices: list, mut: Mutation):
#     offset = 1 if not mut.ref else 0
#
#     check_indices = list(range(mut.start, mut.start + len(mut.ref) + offset))
#     check1 = all([m in indices for m in check_indices])
#     if not check1:
#         print(f"Mutation {mut} not within transcript bounds: {min(indices)} - {max(indices)}.")
#         return seq, indices, False, False
#
#     rel_start, rel_end = indices.index(mut.start)+offset, indices.index(mut.start)+offset+len(mut.ref)
#     acquired_seq = seq[rel_start:rel_end]
#     check2 = acquired_seq == mut.ref
#     if not check2:
#         print(f'Reference allele does not match genome_build allele. {acquired_seq}, {mut.ref}, {mut.start}')
#         consensus_allele = False
#     else:
#         consensus_allele = True
#     if len(mut.ref) == len(mut.alt) > 0:
#         temp_indices = list(range(mut.start, mut.start + len(mut.ref)))
#     else:
#         temp_indices = [indices[indices.index(mut.start)] + v / 1000 for v in list(range(1, len(mut.alt)+1))]
#
#     new_indices = indices[:rel_start] + temp_indices + indices[rel_end:]
#     new_seq = seq[:rel_start] + mut.alt + seq[rel_end:]
#
#     assert len(new_seq) == len(new_indices), f'Error in variant modification: {mut}, {len(new_seq)}, {len(new_indices)}'
#     return new_seq, new_indices, True, consensus_allele
