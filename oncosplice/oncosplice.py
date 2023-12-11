import numpy as np
import pandas as pd
from Bio import pairwise2
import re
from copy import deepcopy
from pathlib import Path
from geney import access_conservation_data, get_correct_gene_file
from oncosplice.spliceai_utils import PredictSpliceAI
from oncosplice.Gene import Gene, Transcript
from oncosplice.variant_utils import Variations, develop_aberrant_splicing

sample_mut_id = 'KRAS:12:25227343:G:T'

def oncosplice(mutation, sai_threshold=0.25, prevalence_threshold=0.25, target_transcripts=None, primary_transcript=False, target_directory=Path('/tamir2/nicolaslynn/projects/mutation_colabs/mrna_database')):
    print(f'>> Processing: {mutation}')
    mutation = Variations(mutation)
    # file = get_correct_gene_file(mutation.gene, target_directory=target_directory)
    file = [f for f in target_directory.glob(f'mrna_*_{mutation.gene}.json')]
    if len(file) == 1:
        file = file[0]
    else:
        return None

    gene = Gene(file=file)
    aberrant_splicing = PredictSpliceAI(mutation, threshold=sai_threshold)

    if target_transcripts:
        reports = pd.concat([oncosplice_transcript(gene.transcripts[transcript_id].generate_protein(), mutation, aberrant_splicing, prevalence_threshold) for
                             transcript_id in target_transcripts if gene.transcripts[transcript_id]['transcript_type'] == 'protein_coding'], axis=1)
    elif primary_transcript:
        reports = oncosplice_transcript(gene.primary_transcript.generate_protein(), mutation, aberrant_splicing, prevalence_threshold)
    else:
        reports = pd.concat([oncosplice_transcript(reference_transcript.generate_protein(), mutation, aberrant_splicing, prevalence_threshold) for
                             reference_transcript in gene if reference_transcript.transcript_type == 'protein_coding'], axis=1)
    return reports

def oncosplice_transcript(reference_transcript, mutation, aberrant_splicing, prevalence_threshold=0.0):
    reports = []
    for i, new_boundaries in enumerate(develop_aberrant_splicing(reference_transcript.exons, aberrant_splicing.aberrant_splicing)):
        variant_transcript = Transcript(deepcopy(reference_transcript).__dict__).set_exons(new_boundaries).generate_mature_mrna(mutations=mutation.mut_id.split('|'), inplace=True).generate_translational_boundaries().generate_protein()
        report = compare_transcripts(reference_transcript, variant_transcript, mutation)
        report['missplicing'] = bool(aberrant_splicing)
        report['aberrant_splicing'] = aberrant_splicing.aberrant_splicing
        report['isoform_prevalence'] = new_boundaries['path_weight']
        report['mutation'] = mutation.mut_id
        report['isoform'] = i
        reports.append(report)

    reports = pd.concat(reports, axis=1).transpose()
    reports['weighted_oncosplice'] = reports.legacy_oncosplice_score * reports.isoform_prevalence
    reports = reports[reports.isoform_prevalence >= prevalence_threshold]
    return reports

def compare_transcripts(reference_transcript, variant_transcript, mut):
    cons_seq, cons_vector = access_conservation_data(reference_transcript.transcript_id)
    if cons_seq == reference_transcript.protein:
        cons_available = True
        cons_vector = cons_vector
    else:
        cons_available = False
        cons_vector = [1] * len(reference_transcript.protein)

    alignment, num_ins, num_del = get_logical_alignment(reference_transcript.protein, variant_transcript.protein)
    deleted, inserted, aligned, unified_seq = get_insertions_and_deletions(alignment)
    window_length = min(76, len(reference_transcript.protein))
    cons_vector = np.array(cons_vector, dtype=float)

    affected_exon, affected_intron, distance_from_5, distance_from_3 = None, None, None, None
    for i, (ex_start, ex_end) in enumerate(reference_transcript.exons):
        if min(ex_start, ex_end) <= mut.start <= max(ex_start, ex_end):
            affected_exon = i + 1
            distance_from_5 = abs(mut.start - ex_start)
            distance_from_3 = abs(mut.start - ex_end)

    for i, (in_start, in_end) in enumerate(reference_transcript.exons):
        if min(in_start, in_end) <= mut.start <= max(in_start, in_end):
            affected_exon = i + 1
            distance_from_5 = abs(mut.start - in_end)
            distance_from_3 = abs(mut.start - in_start)

    report = {f'reference_{k}': v for k, v in reference_transcript.constructor.items()}
    report['reference_mRNA'] = reference_transcript.transcript_seq
    report['reference_CDS_start'] = reference_transcript.transcript_indices.index(reference_transcript.TIS)
    report['reference_pre_mrna'] = reference_transcript.pre_mrna
    report['reference_ORF'] = reference_transcript.pre_mrna
    report['reference_protein'] = reference_transcript.protein

    report.update({f'variant_{k}': v for k, v in variant_transcript.constructor.items()})
    report['variant_mRNA'] = variant_transcript.transcript_seq
    report['variant_CDS_start'] = variant_transcript.transcript_indices.index(variant_transcript.TIS)
    report['variant_pre_mrna'] = variant_transcript.pre_mrna
    report['variant_ORF'] = variant_transcript.pre_mrna
    report['variant_protein'] = variant_transcript.protein

    report['protein_view'] = unified_seq
    descriptions = define_missplicing_events(reference_transcript.exons, variant_transcript.exons,
                              reference_transcript.rev)
    report['exon_changes'] = '|'.join([v for v in descriptions if v])
    report['splicing_codes'] = summarize_missplicing_event(*descriptions)
    report['ref_prot_length'] = len(reference_transcript.protein)
    report['var_prot_length'] = len(variant_transcript.protein)
    report['preservation'] = aligned/len(reference_transcript.protein)
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
    unified_string = []
    for rel_pos, (ch_a, ch_b) in enumerate(list(zip(alignment.seqA, alignment.seqB))):
        if ch_a == ch_b != '-':
            if last_event == 'ALIGN':
                pass
            elif last_event == 'DEL':
                unified_string.append('<DEL_END>')
            elif last_event == 'INS':
                unified_string.append('<INS_END>')
            aligned_pos.append(rel_pos)
            last_event = 'ALIGN'
            unified_string.append(ch_a)

        elif (ch_a != ch_b and ch_a != '-' and ch_b != '-'):
            deleted_pos.append(rel_pos)
            if last_event == 'DEL':
                deletions[del_start_pos] += ch_a
            elif last_event == 'ALIGN':
                last_event = 'MM'
                del_start_pos = ref_map[rel_pos]
                deletions[del_start_pos] = ch_a
            elif last_event == 'INS':
                last_event = 'MM'
                del_start_pos = ref_map[rel_pos]
                deletions[del_start_pos] = ch_a

            unified_string.append(f'({ch_a}>{ch_b})')

        elif (ch_a != ch_b == '-'):
            deleted_pos.append(rel_pos)
            if last_event == 'DEL':
                deletions[del_start_pos] += ch_a
            elif last_event == 'ALIGN':
                last_event = 'DEL'
                del_start_pos = ref_map[rel_pos]
                deletions[del_start_pos] = ch_a
                unified_string.append('<DEL_START>')
            elif last_event == 'INS':
                last_event = 'DEL'
                del_start_pos = ref_map[rel_pos]
                deletions[del_start_pos] = ch_a
                unified_string.append('<INS_END>')
                unified_string.append('<DEL_START>')
            unified_string.append(ch_a)

        elif ch_b != ch_a == '-':
            point_of_insertion = find_last_real_pos(rel_pos, ref_map)
            inserted_pos.append(point_of_insertion)
            if last_event == 'INS':
                insertions[point_of_insertion] += ch_b
            elif last_event == 'ALIGN':
                last_event = 'INS'
                insertions[point_of_insertion] = ch_b
                unified_string.append('<INS_START>')
            elif last_event == 'DEL':
                last_event = 'INS'
                insertions[point_of_insertion] = ch_b
                unified_string.append('<DEL_END>')
                unified_string.append('<INS_START>')
            unified_string.append(ch_b)

    return deletions, insertions, len(aligned_pos), ''.join(unified_string)


def get_logical_alignment(r, v):
    '''
        get_logical_alignment aims to reduce the number of unique gaps. For the context of dealing with aberrant
        splicing, it is most common that blocks are inserted or deleted and therefore the most likely comparison is
        one in which gaps are minimalized and correspond to those alternative splicing blocks.
    '''
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

    for pos, insertion in insertions.items():
        reach = min(len(insertion) // 2, 38)
        end = np.linspace(0, 1, reach).round(2)
        start = end[::-1]
        # filler = np.concatenate([start, np.zeros(1), end])
        front_end = pos - reach

        if front_end < 0:
            front_end = 0

        back_end = pos + reach + 1
        if back_end > len(unmodified_positions):
            back_end = len(unmodified_positions)

        # rel_start_pos, rel_end_pos = 0, 0
        unmodified_positions[front_end:back_end] = np.zeros(back_end-front_end, dtype=float) #filler[rel_start_pos:rel_end_pos]

        # if back_end > end_pos:
        #     fill_pos = end_pos - front_end
        #     unmodified_positions[front_end:] = filler[:fill_pos]
        # else:
        #     print(front_end, back_end, len(filler))
        #     unmodified_positions[front_end:back_end] = filler

    return unmodified_positions


def calculate_oncosplice_scores(deletions, insertions, cons_vector, W=10):
    modified_positions = 1 - find_unmodified_positions(len(cons_vector), deletions=deletions, insertions=insertions)
    cons_vec = transform_conservation_vector(cons_vector, W=W)
    modified_cons_vector = cons_vec * modified_positions
    modified_cons_vector = sum_conv(modified_cons_vector, W=W) / W
    tenth_largest_score = sorted(list(modified_cons_vector.flatten().tolist()))[-W*2]
    max_score = max(modified_cons_vector)

    if max_score > 0:
        gof_prob = (max_score - tenth_largest_score) / max_score
        lof_prob = 1 - gof_prob
        return {'gof': gof_prob, 'lof': lof_prob, 'pof': 0, 'oncosplice_score': max_score}
    else:
        return {'gof': 0, 'lof': 0, 'pof': 1, 'oncosplice_score': max_score}


# def calculate_oncosplice_scores(deletions, insertions, cons_vec):
#     unmodified_positions = find_unmodified_positions(len(cons_vec), deletions=deletions, insertions=insertions)
#     functional_loss_vector_5 = transform_conservation_vector(cons_vec, W=5) * (1 - unmodified_positions)
#     functional_loss_vector_5 = sum_conv(functional_loss_vector_5, W=5)
#
#     if len(cons_vec) < 76:
#         W=len(cons_vec)-1
#     else:
#         W = 76
#
#     functional_loss_vector_76 = transform_conservation_vector(cons_vec, W=W) * (1 - unmodified_positions)
#     functional_loss_vector_76 = sum_conv(functional_loss_vector_76, W=W)
#
#     return {'oncosplice_score_lof': max(functional_loss_vector_76), 'oncosplice_score_gof': max(functional_loss_vector_5)}


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

