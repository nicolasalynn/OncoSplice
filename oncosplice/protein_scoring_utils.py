import numpy as np
from Bio import pairwise2
# from Bio.pairwise2 import format_alignment

import pandas as pd
import re


def generate_report_single_mut(ref_proteome, var_proteome, missplicing, mutation):
    full_report = []

    for (ref_id, var_id) in [(ref_id, var_id) for ref_id in ref_proteome.keys() for var_id in var_proteome.keys() if ref_id == var_id.split('-')[0]]:

        ### Compare and Score Ref protein and Isoform protein
        ref_prot, var_prot = ref_proteome[ref_id], var_proteome[var_id]

        alignment = get_logical_alignment(ref_prot.protein, var_prot.protein)
        deleted, inserted = get_insertions_and_deletions(alignment)
        W = 76
        if W >= len(ref_prot.protein):
            W = len(ref_prot.protein)-1

        smoothed_conservation_vector = smooth_cons_scores(ref_prot.conservation_vector, W, W // 2)

        deconv_del = calculate_del_penalty(deleted, smoothed_conservation_vector, W, W // 2)
        deconv_ins = calculate_ins_penalty(inserted, smoothed_conservation_vector, W, W // 2)
        oncosplice_score = combine_ins_and_del_scores(deconv_del, deconv_ins, W, W // 2)

        pes, pir, es, ne, ir = define_missplicing_events(ref_prot.exon_boundaries, var_prot.exon_boundaries, ref_prot.rev)
        description = '|'.join([v for v in [pes, pir, es, ne, ir] if v])

        ### Record Data
        report = pd.Series()
        report.gene = ref_prot.gene_name
        report.chrom = ref_prot.chrm
        report.pos = mutation.start
        report.ref = mutation.ref
        report.alt = mutation.alt
        report.transcipt_id = ref_prot.transcript_id
        report.isoform_id = var_prot.transcript_id.split('-')[-1]
        report.missed_acceptors = [pos for pos in missplicing.get('missed_acceptors', {}).keys() if pos in ref_prot.acceptors]
        report.missed_donors = [pos for pos in missplicing.get('missed_donors', {}).keys() if pos in ref_prot.donors]
        report.discovered_acceptors = list(missplicing.get('discovered_acceptors', {}).keys())
        report.discovered_donors = list(missplicing.get('discovered_donors', {}).keys())
        report.isoform_prevalence = var_prot.penetrance_weight
        report.missplicing_event = description
        report.missplicing_event_summary = summarize_missplicing_event(pes, pir, es, ne, ir)
        report.reference_protein_length = len(ref_prot.protein)
        report.variant_protein_length = len(var_prot.protein)
        report.insertions = ','.join(list(inserted.values()))
        report.deletions = ','.join(list(deleted.values()))
        report.oncosplice_score = oncosplice_score
        report.deconv_ins = deconv_ins
        report.deconv_del = deconv_del
        report.mutation_distance_from_5 = min([abs(pos - mutation.start) for pos in ref_prot.acceptors + [ref_prot.transcript_start]])
        report.mutation_distance_from_3 = min([abs(pos - mutation.start) for pos in ref_prot + [ref_prot.transcript_end]])

        full_report.append(report)

    return pd.concat(full_report, axis=1).transpose()


def define_missplicing_events(ref_exons, var_exons, rev):

    ref_introns = [(ref_exons[i][1], ref_exons[i+1][0]) for i in range(len(ref_exons) - 1)]
    var_introns = [(var_exons[i][1], var_exons[i+1][0]) for i in range(len(var_exons) - 1)]

    if not rev:
        partial_exon_skipping = ','.join([f'Exon {exon_count+1} truncated: {(t1, t2)} --> {(s1, s2)}' for (s1, s2) in var_exons for exon_count, (t1, t2) in enumerate(ref_exons) if (s1 == t1 and s2 < t2) or (s1 > t1 and s2 == t2)])
        partial_intron_retention = ','.join([f'Intron {intron_count + 1} partially retained: {(t1, t2)} --> {(s1, s2)}' for (s1, s2) in var_introns for intron_count, (t1, t2) in enumerate(ref_introns) if (s1 == t1 and s2 < t2) or (s1 > t1 and s2 == t2)])

    else:
        partial_exon_skipping = ','.join([f'Exon {exon_count+1} truncated: {(t1, t2)} --> {(s1, s2)}' for (s1, s2) in var_exons for exon_count, (t1, t2) in enumerate(ref_exons) if (s1 == t1 and s2 > t2) or (s1 < t1 and s2 == t2)])
        partial_intron_retention = ','.join([f'Intron {intron_count + 1} partially retained: {(t1, t2)} --> {(s1, s2)}' for (s1, s2) in var_introns for intron_count, (t1, t2) in enumerate(ref_introns) if (s1 == t1 and s2 > t2) or (s1 < t1 and s2 == t2)])


    exon_skipping = ','.join([f'Exon {exon_count + 1} skipped: {(t1, t2)}' for exon_count, (t1, t2) in enumerate(ref_exons) if
                     t1 not in [s1 for s1, s2 in var_exons] and t2 not in [s2 for s1, s2 in var_exons]])
    novel_exons = ','.join([f'Novel Exon: {(t1, t2)}' for (t1, t2) in var_exons if
                   t1 not in [s1 for s1, s2 in ref_exons] and t2 not in [s2 for s1, s2 in ref_exons]])
    intron_retention = ','.join([f'Intron {intron_count + 1} retained: {(t1, t2)}' for intron_count, (t1, t2) in enumerate(ref_introns) if
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

    return event if len(event) > 1 else event[0]


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

    return deletions, insertions


def get_logical_alignment(r, v):
    '''
        get_logical_alignment aims to reduce the number of unique gaps. For the context of dealing with aberrant
        splicing, it is most common that blocks are inserted or deleted and therefore the most likely comparison is
        one in which gaps are minimalized and correspond to those alternative splicing blocks.
    '''
    alignments = pairwise2.align.globalms(r, v, 1, -1, -3, -0.5)

    if len(alignments) == 1:
        seqa, seqb = re.sub('-+', '-', alignments[0].seqA), re.sub('-+', '-', alignments[0].seqB)
        return alignments[0], seqa.count('-') + seqb.count('-')

    # This calculates the number of gaps in each alignment.
    block_counts = [re.sub('-+', '-', al.seqA).count('-') + re.sub('-+', '-', al.seqB).count('-') for al in
                    alignments]

    # We return the alignment with the smallest number of gaps.
    try:
        return alignments[block_counts.index(min(block_counts))]
    except AttributeError:
        print("Error get_logical_alignment")
        print(block_counts, alignments, r, v)


def smooth_cons_scores(cons_scores, W, PADDING):
    new_scores = []
    for i in range(len(cons_scores)):
        if i < PADDING:
            temp_vec = cons_scores[:i + PADDING + 1]
        elif i > len(cons_scores) - PADDING:
            temp_vec = cons_scores[i - PADDING:]
        else:
            temp_vec = cons_scores[i - PADDING:i + PADDING + 1]
        len_temp_vec = len(temp_vec)
        new_scores.append(sum(temp_vec) / len_temp_vec)
    new_scores = np.array(new_scores) / max(new_scores)
    assert len(new_scores) == len(
        cons_scores), f'Smoothed scores are not same length... {len(cons_scores)}, {len(new_scores)}'
    # assert round(sum(new_scores), 10) == 1, f'New score sum != to 1: {sum(new_scores)}'
    # fig = tpl.figure()
    # x = np.arange(0, len(new_scores))
    # fig.plot(x, new_scores, width=60, height=20)
    # fig.show()
    return new_scores


def calculate_del_penalty(deleted_domains, cons_scores, W, PADDING):
    penalty = np.zeros(cons_scores.size)
    for dp_pos, dp_seq in deleted_domains.items():
        dw = max(1.0, len(dp_seq) / W)
        penalty[dp_pos:dp_pos + len(dp_seq) + 1] = cons_scores[dp_pos:dp_pos + len(dp_seq) + 1] * dw
    return penalty


def calculate_ins_penalty(inserted_domains, cons_scores, W, PADDING):
    penalty = np.zeros(cons_scores.size)
    for ip_pos, ip_seq in inserted_domains.items():
        reach = min(W, len(ip_seq))
        iw = max(1.0, len(ip_seq) / W)
        # penalty[ip_pos] = iw*cons_scores[ip_pos]
        # for ipv in list(range(ip_pos-reach, ip_pos + reach+1)):
        for ipv in [ip_pos, ip_pos + 1]:
            if ipv > len(cons_scores) - 1:
                pass
            elif ipv < 0:
                pass
            else:
                penalty[ipv] = iw * cons_scores[ipv] + penalty[ipv]
    return penalty


def combine_ins_and_del_scores(d_cons_scores, i_cons_scores, W, PADDING):
    combined_scores = [a + b for a, b in list(zip(d_cons_scores, i_cons_scores))]
    penalty = []
    # for i in range(PADDING, len(d_cons_scores) - PADDING):
    for i in range(PADDING, len(combined_scores) - PADDING):
        penalty.append(sum(combined_scores[i - PADDING:i + PADDING + 1]))

    # print(f"Penalty ({sum(penalty)}: {penalty}")
    # median_p = np.median(penalty)
    # penalty = [p - median_p for p in penalty]
    return max(penalty)



