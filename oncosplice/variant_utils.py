import numpy as np
import networkx as nx
import random
from pathlib import Path
from dataclasses import dataclass
from geney import is_monotonic
from geney.genomic.translation import END_CODONS

def check_pairwise_variant_compatibility(mut1, mut2):
    if mut1 == mut2:
        return True
    if mut1.vartype == 'DEL':
        deleted_positions1 = list(range(mut1.start, mut1.start + len(mut1.ref)))
    else:
        deleted_positions1 = []
    if mut2.vartype == 'DEL':
        deleted_positions2 = list(range(mut2.start, mut2.start + len(mut2.ref)))
    else:
        deleted_positions2 = []

    if mut1.vartype == 'SNP':
        snp_pos1 = mut1.start
    else:
        snp_pos1 = 0
    if mut2.vartype == 'SNP':
        snp_pos2 = mut2.start
    else:
        snp_pos2 = 0

    if mut1.vartype == 'INS':
        ins_pos1 = mut1.start
    else:
        ins_pos1 = 0
    if mut2.vartype == 'INS':
        ins_pos2 = mut2.start
    else:
        ins_pos2 = 0

    if snp_pos1 == snp_pos2 != 0:
        return False

    if snp_pos1 in deleted_positions2 or snp_pos2 in deleted_positions1:
        return False

    if ins_pos1 == ins_pos2 != 0:
        return False

    if len(np.intersect1d(deleted_positions1, deleted_positions2)) > 0:
        return False

    return True
def check_variant_compatibility(mutations):
    return all([check_pairwise_variant_compatibility(mut1, mut2) for mut1 in mutations for mut2 in mutations])


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

    def __eq__(self, other):
        return all([self.chrom == other.chrom, self.start == other.start, self.ref == other.ref, self.alt == other.alt])
    def __lt__(self, other):
        return self.start < other.start


class Variations:
    def __init__(self, epistatic_set):
        self.variants = sorted([Mutation(m) for m in epistatic_set.split('|')])
        self.mut_id = epistatic_set
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

    def __iter__(self):
        self.current_index = 0
        return self
    def __next__(self):
        if self.current_index < len(self.variants):
            x = self.variants[self.current_index]
            self.current_index += 1
            return x
        raise StopIteration

    @property
    def congruent(self):
        return check_variant_compatibility(self.variants)

    @property
    def file_identifier_json(self):
        return Path(self.file_identifier + '.json')

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

    # print(f'mutation: {mut}')
    # print(f"New Indiced: {temp_indices}")
    new_indices = indices[:rel_start] + temp_indices + indices[rel_end:]
    new_seq = seq[:rel_start] + mut.alt + seq[rel_end:]

    assert len(new_seq) == len(new_indices), f'Error in variant modification: {mut}, {len(new_seq)}, {len(new_indices)}'
    assert is_monotonic(list(filter((-1).__ne__, new_indices))), f'Mut indices are not monotonic.'
    return new_seq, new_indices, True, consensus_allele

def find_new_tis(seq, indices, tis, tts, tid='', data_path='/tamir2/shaicohen1/share/titer-master'):
    from oncosplice.titer_utils import run_through_titer
    new_tis, _, _, _, _ = run_through_titer(mut_seq=seq,
                                           mut_coords=indices,
                                           ref_sc_coord=tis,
                                           ref_tts_coord=tts,
                                           ref_id=tid,
                                           data_path=data_path,
                                           titer_model='',
                                           all_test_seqs={})
    return new_tis

def find_new_tts(seq, indices, tis):
    seq, indices = seq[indices.index(tis):], indices[indices.index(tis):]
    pos_options = [i for i in list(range(0, len(seq), 3)) if seq[i:i + 3] in END_CODONS and i+3 <= len(seq)]
    if len(pos_options) == 0:
        return indices[0] #[len(seq) - (len(seq) % 3) - 1]
    assert min(pos_options) % 3 == 0, f'{min(pos_options)} not divisible by three.'
    return indices[min(pos_options)+2]


def develop_aberrant_splicing(exons, aberrant_splicing):
    boundaries = [lst for lsts in [[a, b] for a, b in exons] for lst in lsts]
    exon_starts, exon_ends = list(zip(*exons))
    transcript_start, transcript_end = exon_starts[0], exon_ends[-1]

    if exon_starts[0] > exon_ends[0]:
        rev = True
    else:
        rev = False

    upper_range, lower_range = max(boundaries), min(boundaries)
    exon_starts = {v: 1 for v in exon_starts}
    exon_ends = {v: 1 for v in exon_ends}

    for k, v in aberrant_splicing.get('missed_donors', {}).items():
        if k in exon_ends.keys():
            exon_ends[k] = max(v['absolute'], 0.001)

    exon_ends.update(
        {k: v['absolute'] for k, v in aberrant_splicing.get('discovered_donors', {}).items() if lower_range <= k <= upper_range})

    for k, v in aberrant_splicing.get('missed_acceptors', {}).items():
        if k in exon_starts.keys():
            exon_starts[k] = max(v['absolute'], 0.001)

    exon_starts.update(
        {k: v['absolute'] for k, v in aberrant_splicing.get('discovered_acceptors', {}).items() if lower_range <= k <= upper_range})

    nodes = [SpliceSite(pos=pos, ss_type=0, prob=prob) for pos, prob in exon_ends.items() if
             lower_range <= pos <= upper_range] + \
            [SpliceSite(pos=pos, ss_type=1, prob=prob) for pos, prob in exon_starts.items() if
             lower_range <= pos <= upper_range]

    nodes = [s for s in nodes if s.prob > 0]
    nodes.sort(key=lambda x: x.pos, reverse=rev)

    # while nodes[0].ss_type == 0:
    #     nodes = nodes[1:]
    #
    # while nodes[-1].ss_type == 1:
    #     nodes = nodes[:-1]

    G = nx.DiGraph()
    G.add_nodes_from([n.pos for n in nodes])
    for i in range(len(nodes)):
        trailing_prob, in_between = 0, []
        for j in range(i + 1, len(nodes)):
            curr_node, next_node = nodes[i], nodes[j]
            spread = curr_node.ss_type in in_between
            in_between.append(next_node.ss_type)

            if curr_node.ss_type != next_node.ss_type:
                if spread:
                    new_prob = next_node.prob - trailing_prob
                    if new_prob <= 0:
                        break

                    G.add_edge(curr_node.pos, next_node.pos)
                    G.edges[curr_node.pos, next_node.pos]['weight'] = new_prob
                    trailing_prob += next_node.prob

                else:
                    G.add_edge(curr_node.pos, next_node.pos)
                    G.edges[curr_node.pos, next_node.pos]['weight'] = next_node.prob
                    trailing_prob += next_node.prob

    new_paths, prob_sum = {}, 0
    for i, path in enumerate(nx.all_simple_paths(G, transcript_start, transcript_end)):
        curr_prob = path_weight_mult(G, path, 'weight')
        prob_sum += curr_prob
        new_paths[i] = {'acceptors': sorted([p for p in path if p in exon_starts.keys() and p != transcript_start], reverse=rev),
                        'donors': sorted([p for p in path if p in exon_ends.keys() and p != transcript_end], reverse=rev),
                        'path_weight': curr_prob}

    for i, d in new_paths.items():
        d['path_weight'] = round(d['path_weight'] / prob_sum, 2)

    new_paths = {k: v for k, v in new_paths.items() if v['path_weight'] > 0.01}
    return list(new_paths.values())


def path_weight_mult(G, path, weight):
    multigraph = G.is_multigraph()
    cost = 1
    if not nx.is_path(G, path):
        raise nx.NetworkXNoPath("path does not exist")
    for node, nbr in nx.utils.pairwise(path):
        if multigraph:
            cost *= min(v[weight] for v in G[node][nbr].values())
        else:
            cost *= G[node][nbr][weight]
    return cost


def generate_random_as(t):
    ma = random.sample(t.acceptors, 1)[0]
    md = random.sample(t.donors, 1)[0]
    da = random.sample(list(range(min(t.acceptors), max(t.acceptors))), 1)[0]
    dd = random.sample(list(range(min(t.donors), max(t.donors))), 1)[0]
    return {
        'discovered_acceptors': {da: {'absolute': 0.9}},
        'discovered_donors': {dd: {'absolute': 0.6}},
        'missed_donors': {ma: {'absolute': 0.2}},
        'missed_acceptors': {md: {'absolute': 0.1}},
    }


@dataclass
class SpliceSite(object):
    pos: int
    ss_type: int
    prob: float

    def __post_init__(self):
        pass

    def __lt__(self, other):
        return self.pos < other.pos

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
