import networkx as nx
import numpy as np

from oncosplice.pre_mRNA import pre_mRNA
from oncosplice.SpliceSite import SpliceSite

class mature_mRNA(pre_mRNA):

    def __init__(self, transcript_start: int, transcript_end: int, rev: bool, chrm: str, donors, acceptors,
                 gene_name='undefined', transcript_id='undefined', transcript_type='undefined'):
        pre_mRNA.__init__(self, transcript_start=transcript_start, transcript_end=transcript_end,
                          rev=rev, chrm=chrm, gene_name=gene_name, transcript_id=transcript_id,
                          transcript_type=transcript_type)

        self.donors = donors
        self.acceptors = acceptors

        self.donors.sort(reverse=self.rev)
        self.acceptors.sort(reverse=self.rev)

        self.mature_mrna, self.mature_indices = '', []

    def __len__(self):
        return len(self.mature_mrna)

    def __repr__(self):
        return 'mature_mRNA(transcript_id={tid})'.format(tid=self.transcript_id)

    def __eq__(self, other):
        return self.mature_mrna == other.mature_mrna

    def __contains__(self, pattern):
        if isinstance(pattern, str):
            return pattern in self.mature_mrna
        elif isinstance(pattern, int):
            return pattern in self.mature_indices
        else:
            print("Pass an integer to check against the span of the mature mRNAs's coordinates or a string to check "
                  "against the mature mRNA sequence.")
            return False

    def __valid_mature_mrna(self):
        """
            Description: a mature mRNA is valid if
                1. the length of acceptors equals the length of the donors
                2. the donors and acceptors are sorted in a staggering manner
                3. the first donor is larger (or smaller) than the first acceptor depending on rev or not.
        """

        check1 = len(self.donors) == len(self.acceptors)  # there must be an equal number of acceptors and donors

        # print(f'Acceptors: {self.acceptors}')
        # print(f'Donors: {self.donors}')

        if self.rev:
            check2 = all([self.acceptors[i - 1] > self.donors[i] > self.acceptors[i]
                          for i in range(1, len(self.donors))])
        else:
            check2 = all([self.acceptors[i - 1] < self.donors[i] < self.acceptors[i]
                          for i in range(1, len(self.donors))])

        check3 = all([d in self.pre_indices for d in self.donors + self.acceptors])
        check4 = self.transcript_start in self.pre_indices and self.transcript_end in self.pre_indices
        check5 = len([v for v in self.donors if v in self.acceptors]) == 0
        if all([check1, check2, check3, check4, check5]):
            return True

        else:
            print(f'Cannot generate mature mRNA.\n')
            if not check1:
                print(f"No one-to-one pairing. Donors: {len(self.donors)}, Acceptors: {len(self.acceptors)}")
            if not check2:
                print(f"Donors/acceptor arrangement nonsensical: {self.donors}, {self.acceptors}")
            if not check3:
                print(f"Donors({[v for v in self.donors if v not in self.pre_indices]}),"
                      f" Acceptors({[v for v in self.acceptors if v not in self.pre_indices]}) "
                      f"missing from pre_mRNA transcript indices.")
            if not check4:
                print(f"Transcript Start ({self.transcript_start}) "
                      f"available: {self.transcript_start in self.pre_indices}")
                print(f"Transcript End ({self.transcript_end}) "
                      f"available: {self.transcript_end in self.pre_indices}")
            if not check5:
                print(f"Acceptor/Donor overlap: {[v for v in self.donors if v in self.acceptors]}")
            return False

    def generate_mature_mrna(self, mutations=None, regenerate_pre_mrna=False):
        if not self.pre_mrna or regenerate_pre_mrna:
            self.generate_pre_mrna(mutations=mutations)

        if self.__valid_mature_mrna():
            # Description: using the state's exon start and end indices along with the developed pre mRNA
            # indices and sequence, will cut using the splicing blueprints
            # Modifies the states mature mrna, indices, and exon description in place.

            mature_mrna, mature_indices = '', []
            for i, j in self.exon_boundaries():
                if i == self.transcript_start and i not in self.pre_indices:
                    print("ERROR 1 mature_mRNA")
                    i = self.pre_indices[0]
                if j == self.transcript_end and i not in self.pre_indices:
                    print("ERROR 2 mature_mRNA")
                    j = self.pre_indices[-1]

                rel_start, rel_end = self.pre_indices.index(i), self.pre_indices.index(j)
                mature_mrna += self.pre_mrna[rel_start:rel_end + 1]
                mature_indices.extend(self.pre_indices[rel_start:rel_end + 1])

            self.mature_mrna, self.mature_indices = mature_mrna, mature_indices
        return self

    def exon_boundaries(self):
        exon_starts = [self.transcript_start] + self.acceptors
        exon_ends = self.donors + [self.transcript_end]
        exon_starts.sort(reverse=self.rev)
        exon_ends.sort(reverse=self.rev)
        return list(zip(exon_starts, exon_ends))


    def develop_aberrant_splicing(self, aberrant_splicing, rev):
        upper_range, lower_range = max(self.transcript_start, self.transcript_end), min(self.transcript_start, self.transcript_end)
        exon_starts = {v: 1 for v in self.acceptors + [self.transcript_start]}
        exon_ends = {v: 1 for v in self.donors + [self.transcript_end]}

        for k, v in aberrant_splicing['missed_donors'].items():
            if k in exon_ends.keys():
                exon_ends[k] = v['absolute']

        exon_ends.update(
            {k: v['absolute'] for k, v in aberrant_splicing['discovered_donors'].items() if lower_range <= k <= upper_range})

        for k, v in aberrant_splicing['missed_acceptors'].items():
            if k in exon_starts.keys():
                exon_starts[k] = v['absolute']

        exon_starts.update(
            {k: v['absolute'] for k, v in aberrant_splicing['discovered_acceptors'].items() if lower_range <= k <= upper_range})

        nodes = [SpliceSite(pos=pos, ss_type=0, prob=prob) for pos, prob in exon_ends.items() if
                 lower_range <= pos <= upper_range] + \
                [SpliceSite(pos=pos, ss_type=1, prob=prob) for pos, prob in exon_starts.items() if
                 lower_range <= pos <= upper_range]

        nodes.sort(key=lambda x: x.pos, reverse=rev)

        while nodes[0].ss_type == 0:
            nodes = nodes[1:]

        while nodes[-1].ss_type == 1:
            nodes = nodes[:-1]

        G = nx.DiGraph()
        G.add_nodes_from([n.pos for n in nodes])
        for n in [1, 2]:
            for i in range(len(nodes)):
                trailing_prob, in_between = 0, []
                for j in range(i + 1, len(nodes)):
                    curr_node, next_node = nodes[i], nodes[j]
                    spread = curr_node.ss_type in in_between
                    in_between.append(next_node.ss_type)

                    if curr_node.ss_type != next_node.ss_type:
                        if spread:
                            new_prob = next_node.prob - trailing_prob if n == 1 else curr_node.prob - trailing_prob
                        else:
                            new_prob = next_node.prob if n == 1 else curr_node.prob

                        trailing_prob += next_node.prob

                        if new_prob <= 0:
                            break

                        if n == 1:
                            G.add_edge(curr_node.pos, next_node.pos)
                            G.edges[curr_node.pos, next_node.pos]['weight'] = new_prob

                        else:
                            if G.has_edge(next_node.pos, curr_node.pos):
                                G.edges[next_node.pos, curr_node.pos]['weight'] = max(
                                    G.edges[next_node.pos, curr_node.pos]['weight'], new_prob)
                            else:
                                G.add_edge(next_node.pos, curr_node.pos)
                                G.edges[next_node.pos, curr_node.pos]['weight'] = new_prob

            nodes.reverse()

        new_paths, prob_sum = {}, 0
        for i, path in enumerate(nx.all_simple_paths(G, self.transcript_start, self.transcript_end)):
            curr_prob = path_weight_mult(G, path, 'weight')
            prob_sum += curr_prob
            new_paths[i] = {'acceptors': sorted([p for p in path if p in exon_starts.keys() and p != self.transcript_start], reverse=rev),
                            'donors': sorted([p for p in path if p in exon_ends.keys() and p != self.transcript_end], reverse=rev),
                            'path_weight': curr_prob}


        for i, d in new_paths.items():
            d['path_weight'] = round(d['path_weight'] / prob_sum, 2)

        new_paths = {k: v for k, v in new_paths.items() if v['path_weight'] > 0.05}
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