import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
import pandas as pd
import re
import json
from oncosplice import oncosplice_setup
from oncosplice.variant_utils import Mutation, EpistaticSet
from oncosplice.oncosplice_score import calculate_oncosplice_scores, calculate_legacy_oncosplice_score
from oncosplice.spliceai_utils import missplicing_bool

aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -3
aligner.extend_gap_score = 0
aligner.target_end_gap_score = 0
aligner.query_end_gap_score = 0

def generate_report(ref_proteome, var_proteome, missplicing, mutation):
    full_report = []
    for (ref_id, var_id) in [(ref_id, var_id) for ref_id in ref_proteome.keys() for var_id in var_proteome.keys() if
                             ref_id == var_id.split('-')[0]]:

        ### Compare and Score Ref protein and Isoform protein
        ref_prot, var_prot = ref_proteome[ref_id], var_proteome[var_id]
        if len(ref_prot.protein) < 20:
            continue

        # if not ref_prot.cons_available:
        #     continue

        no_start_codon = False
        if not var_prot.protein:
            no_start_codon = True
            var_prot.protein = '*'

        if oncosplice_setup.get('show_output', False):
            print(f"Transcript: {var_prot.transcript_id.split('-')[-1]}")

        alignment, num_ins, num_del = get_logical_alignment(ref_prot.protein, var_prot.protein)
        deleted, inserted, aligned = get_insertions_and_deletions(alignment)

        window_length = min(76, len(ref_prot.protein))
        cons_vector = np.array(ref_prot.conservation_vector, dtype=float)
        scores = calculate_oncosplice_scores(deleted, inserted, cons_vector, window_length//2)
        legacy_scores = calculate_legacy_oncosplice_score(deleted, inserted, cons_vector,
                                                          window_length)

        pes, pir, es, ne, ir = define_missplicing_events(ref_prot.exon_boundaries(), var_prot.exon_boundaries(),
                                                         ref_prot.rev)
        description = '|'.join([v for v in [pes, pir, es, ne, ir] if v])

        ref_exons = ref_prot.exon_boundaries()
        ref_introns = [(ref_exons[i][1], ref_exons[i + 1][0]) for i in range(len(ref_exons) - 1)]
        affected_exon, affected_intron, closest_donor, closest_acceptor = '-', '-', '-', '-'
        for ex_num, (ex1, ex2) in enumerate(ref_exons):
            if (not ref_prot.rev and ex1 <= mutation.start <= ex2) or (ref_prot.rev and ex1 >= mutation.start >= ex2):
                affected_exon = ex_num
                closest_donor = abs(ex2 - mutation.start)
                closest_acceptor = abs(ex1 - mutation.start)
                break

        for int_num, (in1, in2) in enumerate(ref_introns):
            if (not ref_prot.rev and in1 < mutation.start < in2) or (ref_prot.rev and in1 > mutation.start > in2):
                affected_intron = int_num
                closest_donor = abs(in1 - mutation.start)
                closest_acceptor = abs(in2 - mutation.start)
                break

        report = {}
        report['gene'] = ref_prot.gene_name
        report['chrom'] = ref_prot.chrm
        report['mut_id'] = mutation.mut_id
        report['pos'] = mutation.start if isinstance(mutation, Mutation) else ','.join(
            [str(m.start) for m in mutation.variants])
        report['ref'] = mutation.ref if isinstance(mutation, Mutation) else ','.join([m.ref for m in mutation.variants])
        report['alt'] = mutation.alt if isinstance(mutation, Mutation) else ','.join([m.alt for m in mutation.variants])
        report['strand'] = '+' if not ref_prot.rev else '-'
        report['transcript_id'] = ref_prot.transcript_id
        report['ensembl_transcript_id'] = ref_prot.transcript_id.split('.')[0]
        report['isoform_id'] = var_prot.transcript_id.split('-')[-1]
        report['full_missplicing'] = json.dumps(missplicing)
        # report['missplicing_flag'] = missplicing_bool(missplicing)
        report['isoform_prevalence'] = var_prot.penetrance
        report['no_start_codon_found'] = no_start_codon
        report['missplicing_event'] = description
        report['missplicing_event_summary'] = summarize_missplicing_event(pes, pir, es, ne, ir)
        report['reference_protein_length'] = len(ref_prot.protein)
        report['variant_protein_length'] = len(var_prot.protein)
        report['preservation'] = aligned/len(ref_prot.protein)
        report['num_insertions'] = num_ins
        report['num_deletions'] = num_del
        report['insertions'] = json.dumps(inserted)
        report['deletions'] = json.dumps(deleted)
        report['mut_exon_residence'] = affected_exon
        report['mut_intron_residence'] = affected_intron
        report['mutation_distance_from_5'] = closest_acceptor
        report['mutation_distance_from_3'] = closest_donor
        report['consensus_allele_match'] = var_prot.consensus_allele_match
        report['cons_available'] = ref_prot.cons_available
        report.update(legacy_scores)
        report.update(scores)
        report = pd.Series(report)
        full_report.append(report)

    if not full_report:
        full_report = pd.DataFrame()
    else:
        full_report = pd.concat(full_report, axis=1).transpose()
    return full_report


def define_missplicing_events(ref_exons, var_exons, rev):
    ref_introns = [(ref_exons[i][1], ref_exons[i + 1][0]) for i in range(len(ref_exons) - 1)]
    var_introns = [(var_exons[i][1], var_exons[i + 1][0]) for i in range(len(var_exons) - 1)]
    num_ref_exons = len(ref_exons)
    num_ref_introns = len(ref_introns)
    if not rev:
        partial_exon_skipping = ','.join(
            [f'Exon {exon_count + 1}/{num_ref_exons} truncated: {(t1, t2)} --> {(s1, s2)}' for (s1, s2) in var_exons for
             exon_count, (t1, t2) in enumerate(ref_exons) if (s1 == t1 and s2 < t2) or (s1 > t1 and s2 == t2)])
        partial_intron_retention = ','.join(
            [f'Intron {intron_count + 1}/{num_ref_introns} partially retained: {(t1, t2)} --> {(s1, s2)}' for (s1, s2)
             in var_introns for intron_count, (t1, t2) in enumerate(ref_introns) if
             (s1 == t1 and s2 < t2) or (s1 > t1 and s2 == t2)])

    else:
        partial_exon_skipping = ','.join(
            [f'Exon {exon_count + 1}/{num_ref_exons} truncated: {(t1, t2)} --> {(s1, s2)}' for (s1, s2) in var_exons for
             exon_count, (t1, t2) in enumerate(ref_exons) if (s1 == t1 and s2 > t2) or (s1 < t1 and s2 == t2)])
        partial_intron_retention = ','.join(
            [f'Intron {intron_count + 1}/{num_ref_introns} partially retained: {(t1, t2)} --> {(s1, s2)}' for (s1, s2)
             in var_introns for intron_count, (t1, t2) in enumerate(ref_introns) if
             (s1 == t1 and s2 > t2) or (s1 < t1 and s2 == t2)])

    exon_skipping = ','.join(
        [f'Exon {exon_count + 1}/{num_ref_exons} skipped: {(t1, t2)}' for exon_count, (t1, t2) in enumerate(ref_exons)
         if
         t1 not in [s1 for s1, s2 in var_exons] and t2 not in [s2 for s1, s2 in var_exons]])
    novel_exons = ','.join([f'Novel Exon: {(t1, t2)}' for (t1, t2) in var_exons if
                            t1 not in [s1 for s1, s2 in ref_exons] and t2 not in [s2 for s1, s2 in ref_exons]])
    intron_retention = ','.join(
        [f'Intron {intron_count + 1}/{num_ref_introns} retained: {(t1, t2)}' for intron_count, (t1, t2) in
         enumerate(ref_introns) if
         t1 not in [s1 for s1, s2 in var_introns] and t2 not in [s2 for s1, s2 in var_introns]])

    return partial_exon_skipping, partial_intron_retention, exon_skipping, novel_exons, intron_retention


def summarize_missplicing_event(pes, pir, es, ne, ir):
    event = []
    if pes:
        event.append('PES')
    if es:
        event.append('ES')
    if pir:
        event.append('PIR')
    if ir:
        event.append('IR')
    if ne:
        event.append('NE')

    if len(event) > 1:
        return event
    elif len(event) == 1:
        return event[0]
    else:
        return '-'


def get_insertions_and_deletions(alignment):
    def find_last_real_pos(pos, pos_map):
        if pos <= min(list(pos_map.keys())):
            return pos_map[min(list(pos_map.keys()))]

        while pos not in pos_map.keys() and pos > min(list(pos_map.keys())):
            pos -= 1

        return pos_map[pos]

    ref_map = {j: p for p, j in enumerate([i for i, nt in enumerate(list(alignment.seqA)) if nt != '-'])}
    aligned_pos, deleted_pos, inserted_pos, missmatched_pos = [], [], [], []
    insertions, deletions = {}, {}
    last_event = 'ALIGN'
    del_start_pos = ''
    for rel_pos, (ch_a, ch_b) in enumerate(list(zip(alignment.seqA, alignment.seqB))):
        if ch_a == ch_b != '-':
            aligned_pos.append(rel_pos)
            last_event = 'ALIGN'

        elif (ch_a != ch_b == '-') or (ch_a != ch_b and ch_a != '-' and ch_b != '-'):
            deleted_pos.append(rel_pos)
            if last_event == 'DEL':
                deletions[del_start_pos] += ch_a
            else:
                last_event = 'DEL'
                del_start_pos = ref_map[rel_pos]
                deletions[del_start_pos] = ch_a

        elif ch_b != ch_a == '-':
            point_of_insertion = find_last_real_pos(rel_pos, ref_map)
            inserted_pos.append(point_of_insertion)
            if last_event == 'INS':
                insertions[point_of_insertion] += ch_b
            else:
                last_event = 'INS'
                insertions[point_of_insertion] = ch_b

    return deletions, insertions, len(aligned_pos)


class optimal_alignment_class:
    def __init__(self, a, b):
        self.seqA = a
        self.seqB = b


def get_logical_alignment(r, v):
    '''
        get_logical_alignment aims to reduce the number of unique gaps. For the context of dealing with aberrant
        splicing, it is most common that blocks are inserted or deleted and therefore the most likely comparison is
        one in which gaps are minimalized and correspond to those alternative splicing blocks.
    '''
    if len(r) * len(v) > 1:
        alignments = aligner.align(r, v)
        if len(alignments) == 1:
            optimal_alignment = alignments[0]
        else:
            # This calculates the number of gaps in each alignment.
            number_of_gaps = [re.sub('-+', '-', al[0, :]).count('-') + re.sub('-+', '-', al[1, :]).count('-') for al in
                              alignments]
            optimal_alignment = alignments[number_of_gaps.index(min(number_of_gaps))]

        num_insertions = re.sub('-+', '-', optimal_alignment[0, :]).count('-')
        num_deletions = re.sub('-+', '-', optimal_alignment[1, :]).count('-')

        print(optimal_alignment)
        optimal_alignment = optimal_alignment_class(optimal_alignment[0, :], optimal_alignment[1, :])

    else:
        alignments = pairwise2.align.globalms(r, v, 1, -1, -3, 0, penalize_end_gaps=(True, False))

        if len(alignments) == 1:
            optimal_alignment = alignments[0]
        else:
            # This calculates the number of gaps in each alignment.
            number_of_gaps = [re.sub('-+', '-', al.seqA).count('-') + re.sub('-+', '-', al.seqB).count('-') for al in
                              alignments]
            optimal_alignment = alignments[number_of_gaps.index(min(number_of_gaps))]

        num_insertions = re.sub('-+', '-', optimal_alignment.seqA).count('-')
        num_deletions = re.sub('-+', '-', optimal_alignment.seqB).count('-')

        print(format_alignment(*optimal_alignment))
    # We return the alignment with the smallest number of gaps.
    return optimal_alignment, num_insertions, num_deletions


# def smooth_cons_scores(cons_scores, W):
#     PADDING = W // 2
#     new_scores = []
#     for i in range(len(cons_scores)):
#         if i < PADDING:
#             temp_vec = cons_scores[:i + PADDING + 1]
#         elif i > len(cons_scores) - PADDING:
#             temp_vec = cons_scores[i - PADDING:]
#         else:
#             temp_vec = cons_scores[i - PADDING:i + PADDING + 1]
#         len_temp_vec = len(temp_vec)
#         new_scores.append(sum(temp_vec) / len_temp_vec)
#     new_scores = np.array(new_scores) / max(new_scores)
#     assert len(new_scores) == len(
#         cons_scores), f'Smoothed scores are not same length... {len(cons_scores)}, {len(new_scores)}'
#     # assert round(sum(new_scores), 10) == 1, f'New score sum != to 1: {sum(new_scores)}'
#     # fig = tpl.figure()
#     # x = np.arange(0, len(new_scores))
#     # fig.plot(x, new_scores, width=60, height=20)
#     # fig.show()
#     return new_scores

#
# def calculate_del_penalty(deleted_domains, cons_scores, W):
#     penalty = np.zeros(cons_scores.size)
#     for dp_pos, dp_seq in deleted_domains.items():
#         dw = max(1.0, len(dp_seq) / W)
#         penalty[dp_pos:dp_pos + len(dp_seq) + 1] = cons_scores[dp_pos:dp_pos + len(dp_seq) + 1] * dw
#     return penalty
#
#
# def calculate_ins_penalty(inserted_domains, cons_scores, W):
#     penalty = np.zeros(cons_scores.size)
#     for ip_pos, ip_seq in inserted_domains.items():
#         # reach = min(W, len(ip_seq))
#         iw = max(1.0, len(ip_seq) / W)
#         # penalty[ip_pos] = iw*cons_scores[ip_pos]
#         # for ipv in list(range(ip_pos-reach, ip_pos + reach+1)):
#         for ipv in [ip_pos, ip_pos + 1]:
#             if ipv > len(cons_scores) - 1:
#                 pass
#             elif ipv < 0:
#                 pass
#             else:
#                 penalty[ipv] = iw * cons_scores[ipv] + penalty[ipv]
#     return penalty
#
# def combine_ins_and_del_scores(d_cons_scores, i_cons_scores, W):
#     PADDING = W // 2
#     combined_scores = [a + b for a, b in list(zip(d_cons_scores, i_cons_scores))]
#     penalty = []
#
#     for i in range(PADDING, len(combined_scores) - PADDING):
#         penalty.append(sum(combined_scores[i - PADDING:i + PADDING + 1]))
#
#     return max(penalty)
# #

# def find_unmodified_positions(lp, deletions, insertions, W):
#     unmodified_positions = list(range(lp))
#     modified_pos = []
#     for pos, deletion in deletions.items():
#         modified_pos.extend(list(range(pos, pos + len(deletion))))
#
#     max_reach = W // 2
#     for pos, insertion in insertions.items():
#         reach = min(len(insertion) // 2, max_reach)
#         modified_pos.extend(list(range(pos - reach, pos + reach + 1)))
#
#     return [v for v in unmodified_positions if v not in list(set(modified_pos))]
#
# def window_matching(unmodified_positions, l_p, W):
#     alignment_vector = [1 if i in unmodified_positions else 0 for i in range(l_p)]
#     convolver = np.ones(W)
#     convolving_length = np.array([min(l_p + W - i, W, i) for i in range(W // 2, l_p + W // 2)])
#     match_ratios = np.convolve(alignment_vector, convolver, mode='same') / (convolving_length // 2) - 1
#     return match_ratios
#
# def mask_matched_positions(cons_vec, unmodified_positions):
#     masked_positions = [val for val in list(range(len(cons_vec))) if val in unmodified_positions]
#     cons_vec[masked_positions] = 0
#     return cons_vec
#
# def window_conv(cons_vec, W):
#     return np.convole(cons_vec, np.zeros(W), mode='same')
#
# def transform_conservation_vector(c, W):
#     temp_W = W // 4
#     convolver = np.ones(temp_W)
#     convolving_length = np.array([min(len(c) + temp_W - i, temp_W, i) for i in range(temp_W // 2, len(c) + temp_W // 2)])
#     c1 = np.convolve(c, convolver, mode='same') / convolving_length
#     c2 = -c1 / max(abs(c1))
#     c3 = 100 ** c2
#     c_sum = sum(c3)
#     c_len = len(c3)
#     c4 = c3 * c_len / c_sum
#     assert len(c4) == len(c), f'Length of modified conservation not equal.'
#     # W_max = max(np.convolve(c4, np.ones(W), mode='same'))
#     return c4, 1 #, W_max
#
# def new_oncosplice_scoring(unmodified_positions, cons_vec, W):
#     cons_vec, _ = transform_conservation_vector(cons_vec, W)
#     alignment_ratio = window_matching(unmodified_positions, len(cons_vec), W)
#     functional_loss = mask_matched_positions(cons_vec.copy(), unmodified_positions)
#     s = alignment_ratio * functional_loss / len(cons_vec)
#     s = np.convolve(s, np.ones(W), 'same')
#     return cons_vec, s.min(), s.max(), sum(s)
#
#
