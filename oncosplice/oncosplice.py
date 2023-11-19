import numpy as np
import pandas as pd
from Bio import Align
import re
import json
from copy import deepcopy
from geney import query_rate4site_db
from oncosplice.spliceai_utils import PredictSpliceAI
from oncosplice.Gene import Gene
from oncosplice.variant_utils import Variations, develop_aberrant_splicing

sample_mut_id = 'KRAS:12:25227343:G:T'


aligner = Align.PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = -1
aligner.open_gap_score = -3
aligner.extend_gap_score = 0
aligner.target_end_gap_score = 0
aligner.query_end_gap_score = 0

def oncosplice(mutation, sai_threshold=0.25, explicit=False):
    print(f'>> Processing: {mutation}')
    mutation = Variations(mutation)
    aberrant_splicing = PredictSpliceAI(mutation, sai_threshold)

    reports = []
    reference_transcript = Gene(mutation.gene).primary_transcript
    for i, new_boundaries in enumerate(develop_aberrant_splicing(reference_transcript.exons, aberrant_splicing.aberrant_splicing)):
        variant_transcript = deepcopy(reference_transcript)
        variant_transcript.transcript_seq, variant_transcript.indices = variant_transcript.generate_mature_mrna(mutations=mutation.mut_id.split('|'))
        variant_transcript.set_exons(new_boundaries).generate_translation_boundaries()

        # Generating data
        report = compare_transcripts(reference_transcript, variant_transcript)
        report['missplicing'] = bool(aberrant_splicing)
        report['aberrant_splicing'] = aberrant_splicing.aberrant_splicing
        report['isoform_prevalence'] = new_boundaries['path_weight']
        reports.append(pd.Series(report))

    reports = pd.concat(reports)
    reports['weighted_oncosplice'] = reports.oncosplice_score * reports.isoform_prevalence
    if explicit:
        return explicit
    else:
        return reports.groupby('transcript_id').weighted_oncosplice.mean().max()


def compare_transcripts(reference_transcript, variant_transcript, mut):
    reference_protein, variant_protein = reference_transcript.generate_protein(), variant_transcript.generate_protein()
    cons_seq, cons_vector = query_rate4site_db(reference_transcript.transcript_id)
    if cons_seq == reference_protein:
        cons_available = True
        cons_vector = reference_protein.conservation_vector
    else:
        cons_available = False
        cons_vector = [1] * len(reference_protein)

    alignment, num_ins, num_del = get_logical_alignment(reference_protein, variant_protein)
    deleted, inserted, aligned = get_insertions_and_deletions(alignment)
    window_length = min(76, len(reference_protein.protein))
    cons_vector = np.array(cons_vector, dtype=float)
    # scores = calculate_oncosplice_scores(deleted, inserted, cons_vector, window_length//2)

    affected_exon, affected_intron, distance_from_5, distance_from_3 = None, None, None, None
    for i, (ex_start, ex_end) in enumerate(reference_transcript.exons()):
        if min(ex_start, ex_end) <= mut.start <= max(ex_start, ex_end):
            affected_exon = i + 1
            distance_from_5 = abs(mut.start - ex_start)
            distance_from_3 = abs(mut.start - ex_end)

    for i, (in_start, in_end) in enumerate(reference_transcript.exons()):
        if min(in_start, in_end) <= mut.start <= max(in_start, in_end):
            affected_exon = i + 1
            distance_from_5 = abs(mut.start - in_end)
            distance_from_3 = abs(mut.start - in_start)

    report = reference_transcript.__dict__
    report['isoform_id'] = variant_transcript.transcript_id.split('-')[-1]
    report['var_TIS'] = variant_transcript.TIS
    report['var_TIS'] = variant_transcript.TTS
    report['exon_changes'] = '|'.join([v for v in define_missplicing_events(reference_protein.exons(), variant_protein.exons(),
                              reference_protein.rev)])
    report['ref_prot_length'] = len(reference_protein)
    report['var_prot_length'] = len(variant_protein)
    report['preservation'] = aligned/len(reference_protein)
    report['num_insertions'] = num_ins
    report['num_deletions'] = num_del
    report['insertions'] = json.dumps(inserted)
    report['deletions'] = json.dumps(deleted)
    report['affected_exon'] = affected_exon
    report['affected_intron'] = affected_intron
    report['mutation_distance_from_5'] = distance_from_5
    report['mutation_distance_from_3'] = distance_from_3
    report['cons_available'] = cons_available
    report['legacy_oncosplice_score'] = calculate_legacy_oncosplice_score(deleted, inserted, cons_vector,
                                                      window_length)
    report.update(calculate_oncosplice_scores(deleted, inserted, cons_vector))

    return pd.Series(report)


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
    # if len(r) * len(v) > 1:
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
    #
    # else:
    #     alignments = pairwise2.align.globalms(r, v, 1, -1, -3, 0, penalize_end_gaps=(True, False))
    #
    #     if len(alignments) == 1:
    #         optimal_alignment = alignments[0]
    #     else:
    #         # This calculates the number of gaps in each alignment.
    #         number_of_gaps = [re.sub('-+', '-', al.seqA).count('-') + re.sub('-+', '-', al.seqB).count('-') for al in
    #                           alignments]
    #         optimal_alignment = alignments[number_of_gaps.index(min(number_of_gaps))]
    #
    #     num_insertions = re.sub('-+', '-', optimal_alignment.seqA).count('-')
    #     num_deletions = re.sub('-+', '-', optimal_alignment.seqB).count('-')
    #
    #     print(format_alignment(*optimal_alignment))
    # We return the alignment with the smallest number of gaps.
    return optimal_alignment, num_insertions, num_deletions


def moving_average_conv(vector, W):
    convolving_length = np.array([min(len(vector) + W - i, W, i) for i in range(W // 2, len(vector) + W // 2)], dtype=float)
    return sum_conv(vector, W) / convolving_length
def moving_average_conv_modified(vector, W):
    convolving_length = np.array([min(len(vector) + W - i, W, i) for i in range(W // 2, len(vector) + W // 2)], dtype=float)
    return sum_conv(vector, W) / (convolving_length / 2)
def sum_conv(vector, W):
    return np.convolve(vector, np.ones(W), mode='same')

def transform_conservation_vector(c, W=9):
    factor = (100/W) // 1
    c = np.exp(factor*np.negative(moving_average_conv(c, W)))
    return c*100/max(c)
def find_unmodified_positions(lp, deletions, insertions):
    unmodified_positions = np.ones(lp, dtype=float)
    for pos, deletion in deletions.items():
        unmodified_positions[pos:pos+len(deletion)] = 0

    # max_reach = 32 #W // 2
    for pos, insertion in insertions.items():
        reach = min(len(insertion) // 2, 38)
        # temp = [reach / i for i in range(1, reach//2)]
        unmodified_positions[pos-reach:pos+reach+1] = 0

    return unmodified_positions
def calculate_oncosplice_scores(deletions, insertions, cons_vec):
    unmodified_positions = find_unmodified_positions(len(cons_vec), deletions=deletions, insertions=insertions)

    functional_loss_vector_5 = transform_conservation_vector(cons_vec, W=5) * (1 - unmodified_positions)
    functional_loss_vector_5 = sum_conv(functional_loss_vector_5, W=5)

    if len(cons_vec) < 76:
        W=len(cons_vec)-1
    else:
        W = 76

    functional_loss_vector_76 = transform_conservation_vector(cons_vec, W=W) * (1 - unmodified_positions)
    functional_loss_vector_76 = sum_conv(functional_loss_vector_76, W=W)

    return {'oncosplice_score_lof': max(functional_loss_vector_76), 'oncosplice_score_gof': max(functional_loss_vector_5)}


    # return {'cons_vec': np.array2string(np.around(cons_vec), 3), 'lof_score': abs(min(0, s.min())), 'gof_score': max(0, s.max()), 'oncosplice_score': sum(stemp)/len(cons_vec)}


##### LEGACY ONCOSPLICE CALCS
def legacy_smooth_cons_scores(cons_scores, W):
    return np.exp(np.negative(moving_average_conv_modified(cons_scores, W)))

def calculate_del_penalty(deleted_domains, cons_scores, W):
    penalty = np.zeros(len(cons_scores))
    for dp_pos, dp_seq in deleted_domains.items():
        dw = max(1.0, len(dp_seq) / W)
        penalty[dp_pos:dp_pos + len(dp_seq)] = cons_scores[dp_pos:dp_pos + len(dp_seq)] * dw
    return penalty

def calculate_ins_penalty(inserted_domains, cons_scores, W):
    penalty = np.zeros(len(cons_scores)) #cons_scores.copy()
    for ip_pos, ip_seq in inserted_domains.items():
        reach = min(W//2, len(ip_seq)//2)
        iw = max(1.0, len(ip_seq) / W)
        penalty[ip_pos-reach:ip_pos+reach] = iw*cons_scores[ip_pos-reach:ip_pos+reach]
    return penalty
def combine_ins_and_del_scores(d_cons_scores, i_cons_scores, W):
    combined_scores = d_cons_scores + i_cons_scores
    penalty = sum_conv(combined_scores, W)
    return max(penalty)
def calculate_legacy_oncosplice_score(deletions, insertions, cons_vec, W):
    smoothed_conservation_vector = legacy_smooth_cons_scores(cons_vec, W)
    deconv_del = calculate_del_penalty(deletions, smoothed_conservation_vector, W)
    deconv_ins = calculate_ins_penalty(insertions, smoothed_conservation_vector, W)
    return combine_ins_and_del_scores(deconv_del, deconv_ins, W)

