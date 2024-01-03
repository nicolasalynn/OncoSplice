import numpy as np
import pandas as pd
from Bio import pairwise2
import re
from copy import deepcopy
from pathlib import Path
from geney import access_conservation_data
from oncosplice.spliceai_utils import PredictSpliceAI
from oncosplice.Gene import Gene, Transcript
from oncosplice.variant_utils import Variations, develop_aberrant_splicing
sample_mut_id = 'KRAS:12:25227343:G:T'

def oncosplice(mutation: str, sai_threshold=0.25, prevalence_threshold=0.25, target_directory=Path('/tamir2/nicolaslynn/projects/mutation_colabs/mrna_database')) -> pd.DataFrame:
    '''
        :param mutation: str
                        the genomic variation
        :param sai_threshold: float
                        the threshold for including missplicing predictions in gene builds
        :param prevalence_threshold: float
                        the minimum threshold needed to consider a predicted isoform as valid
        :param target_directory: pathlib.Path
                        the directory on the machine where the mrna annotation files are stored
        :return: a dataframe object
                will contain columns pertinant to assessing mutation pathogenicity including oncosplice score, GOF score, legacy oncosplice score, missplicing,
    '''

    print(f'>> Processing: {mutation}')
    mutation = Variations(mutation)                                             # Generate mutation object
    file = target_directory / f'mrna_{mutation.gene}.json'                      # Gene annotations should be available in the target directory under the file name mrna_gene.json
    if not file.exists():                                                       # Ensure that the annotation file exists
        print(f'Missing annotation files for: {mutation.gene}')
        return pd.DataFrame()

    gene = Gene(file=file)                                                      # We obtain the annotation file and convert it into a Gene object
    aberrant_splicing = PredictSpliceAI(mutation, threshold=sai_threshold)      # SpliceAI predictions are processed and obtained for each mutation
    # Oncosplice obtains predictions for each transcript in the annotation file
    results = pd.concat([oncosplice_transcript(reference_transcript.generate_protein(), mutation, aberrant_splicing, prevalence_threshold) for
                         reference_transcript in gene])

    # Append some additional, uniform information to the results dataframe
    results['mut_id'] = mutation.mut_id
    results['missplicing'] = bool(aberrant_splicing)
    return results

def oncosplice_transcript(reference_transcript: Transcript, mutation: Variations, aberrant_splicing:PredictSpliceAI, prevalence_threshold=0.0, full_output=False) -> pd.DataFrame:
    '''
    :param reference_transcript:
    :param mutation:
    :param aberrant_splicing:
    :param prevalence_threshold:
    :param full_output:
    :return:
    '''
    reports = []
    cons_seq, cons_vector, cons_available = access_conservation_data(reference_transcript.transcript_id)        # access conservation data
    if cons_seq != reference_transcript.protein:                                                                # in the case that conservation is not available, we assign equal importances to all amino acids
        cons_available, cons_vector = False, np.ones(len(reference_transcript.protein), dtype=float)

    # For each transcript, we generate a series of isoforms based on the splice site predictions; each isoform is assigned a prevalence score
    # obtained using simple graph theory where the probability of the edges taken to generate the isoform are multiplied together
    for i, new_boundaries in enumerate(develop_aberrant_splicing(reference_transcript, aberrant_splicing.aberrant_splicing)):

        # The variant transcript is duplicated from the reference transcript and all needed modifications are performed
        variant_transcript = Transcript(deepcopy(reference_transcript).__dict__).set_exons(new_boundaries).generate_mature_mrna(mutations=mutation.mut_id.split('|'), inplace=True).generate_translational_boundaries().generate_protein()

        # The optimal alignment that minimizes gaps between the trnascripts is obtained
        alignment, _, _ = get_logical_alignment(reference_transcript.protein, variant_transcript.protein)

        # Based on the optimal alignment, we can generate the relative locations of insertions and deletions
        deleted, inserted = find_indels_with_mismatches_as_deletions(alignment)

        report = {
            'isoform': i,
            'isoform_prevalence': new_boundaries['path_weight'],
            'legacy_oncosplice_score': calculate_legacy_oncosplice_score(deleted, inserted, cons_vector,
                                                      min(76, len(reference_transcript.protein))),
        }
        report.update(calculate_oncosplice_scores(deleted, inserted, cons_vector))

        if full_output:
            report.update(compare_transcripts(reference_transcript, variant_transcript, mutation))

        reports.append(report)

    reports = pd.DataFrame(reports)
    reports['cons_available'] = cons_available
    return reports[reports.isoform_prevalence >= prevalence_threshold]

def get_logical_alignment(ref_prot, var_prot):
    '''
    :param ref_prot:
    :param var_prot:
    :return:
    '''

    alignments = pairwise2.align.globalms(ref_prot, var_prot, 1, -1, -3, 0, penalize_end_gaps=(True, False))
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


def find_indels_with_mismatches_as_deletions(seqA, seqB):
    # Convert sequences to numpy arrays for element-wise comparison
    ta, tb = np.array(list(seqA)), np.array(list(seqB))

    # Find mismatch positions
    mismatch_positions = (ta != tb) & (ta != '-') & (tb != '-')

    # Replace mismatch positions in seqB with '-'
    tb[mismatch_positions] = '-'
    modified_seqB = ''.join(tb)

    # Function to find continuous gaps using regex
    def find_continuous_gaps(sequence):
        return [(m.start(), m.end()) for m in re.finditer(r'-+', sequence)]

    # Find gaps in both sequences
    gaps_in_A = find_continuous_gaps(seqA)
    gaps_in_B = find_continuous_gaps(modified_seqB)

    # Identify insertions and deletions
    insertions = {start: modified_seqB[start:end].replace('-', '') for start, end in gaps_in_A if
                  seqB[start:end].strip('-')}
    deletions = {start: seqA[start:end].replace('-', '') for start, end in gaps_in_B if seqA[start:end].strip('-')}

    return deletions, insertions

##### NEW ONCOSPLICE SCORING ########################
def moving_average_conv(vector, window_size, factor=1):
    """
    Calculate the moving average convolution of a vector.

    :param vector: Input vector.
    :param window_size: Size of the convolution window.
    :return: Convolved vector as a numpy array.
    """
    convolving_length = np.array([min(len(vector) + window_size - i, window_size, i)
                                  for i in range(window_size // 2, len(vector) + window_size // 2)], dtype=float)
    return np.convolve(vector, np.ones(window_size), mode='same') / (convolving_length / factor)

def transform_conservation_vector(conservation_vector, window_size=10):
    """
    Transforms a conservation vector by applying a moving average convolution and scaling.

    :param conservation_vector: Array of conservation scores.
    :param window_size: Window size for the moving average convolution. Defaults to 10, the average binding site length.
    :return: Transformed conservation vector.
    """
    factor = (100 / window_size) // 1
    moving_avg = moving_average_conv(conservation_vector, window_size)
    transformed_vector = np.exp(factor * np.negative(moving_avg))
    return transformed_vector * 100 / sum(transformed_vector) # or max(transformed_vector)


def find_unmodified_positions(sequence_length, deletions, insertions, reach_limit=38):
    """
    Identify unmodified positions in a sequence given deletions and insertions.

    :param sequence_length: Length of the sequence.
    :param deletions: Dictionary of deletions.
    :param insertions: Dictionary of insertions.
    :param reach_limit: Limit for considering the effect of insertions/deletions.
    :return: Array indicating unmodified positions.
    """
    unmodified_positions = np.ones(sequence_length, dtype=float)

    for pos, deletion in deletions.items():
        deletion_length = len(deletion)
        unmodified_positions[pos:pos + deletion_length] = 0

    for pos, insertion in insertions.items():
        reach = min(len(insertion) // 2, reach_limit)
        front_end, back_end = max(0, pos - reach), min(sequence_length, pos + reach + 1)
        len_start, len_end = pos - front_end, back_end - pos
        unmodified_positions[front_end:back_end + 1] = np.concatenate([np.linspace(0, 1, len_start)[::-1], [0], np.linspace(0, 1, len_end)])

    return unmodified_positions



def calculate_oncosplice_scores(deletions, insertions, cons_vector, window_size=10):
    """
    Calculate oncosplice scores based on conservation vectors and detected sequence modifications.

    :param deletions: Dictionary of deletions in the sequence.
    :param insertions: Dictionary of insertions in the sequence.
    :param cons_vector: Conservation vector.
    :param window_size: Window size for calculations.
    :return: Dictionary of oncosplice scores.
    """
    modified_positions = 1 - find_unmodified_positions(len(cons_vector), deletions, insertions)
    cons_vec = transform_conservation_vector(cons_vector, window_size)
    modified_cons_vector = np.convolve(cons_vec * modified_positions, np.ones(window_size), mode='same') / window_size

    max_score = np.max(modified_cons_vector)
    max_score_indices = np.where(modified_cons_vector == max_score)[0]

    # Exclude windows within one window_size of the max scoring window
    exclusion_zone = set().union(*(range(max(i - window_size, 0), min(i + window_size, len(modified_cons_vector))) for i in max_score_indices))
    second_highest_score = np.max([score for i, score in enumerate(modified_cons_vector) if i not in exclusion_zone])

    gof_prob = (max_score - second_highest_score) / max_score
    return {'gof': gof_prob, 'oncosplice_score': max_score}



def calculate_penalty(domains, cons_scores, W, is_insertion=False):
    """
    Calculate the penalty for mutations (either insertions or deletions) on conservation scores.

    :param domains: Dictionary of mutations (inserted or deleted domains).
    :param cons_scores: Conservation scores.
    :param W: Window size.
    :param is_insertion: Boolean flag to indicate if the mutation is an insertion.
    :return: Penalty array.
    """
    penalty = np.zeros(len(cons_scores))
    for pos, seq in domains.items():
        mutation_length = len(seq)
        weight = max(1.0, mutation_length / W)

        if is_insertion:
            reach = min(W // 2, mutation_length // 2)
            penalty[pos - reach:pos + reach] = weight * cons_scores[pos - reach:pos + reach]
        else:  # For deletion
            penalty[pos:pos + mutation_length] = cons_scores[pos:pos + mutation_length] * weight

    return penalty


def calculate_legacy_oncosplice_score(deletions, insertions, cons_vec, W):
    """
    Calculate the legacy Oncosplice score based on deletions, insertions, and conservation vector.

    :param deletions: Dictionary of deletions.
    :param insertions: Dictionary of insertions.
    :param cons_vec: Conservation vector.
    :param W: Window size.
    :return: Legacy Oncosplice score.
    """
    smoothed_conservation_vector = np.exp(np.negative(moving_average_conv(cons_vec, W, 2)))
    del_penalty = calculate_penalty(deletions, smoothed_conservation_vector, W, is_insertion=False)
    ins_penalty = calculate_penalty(insertions, smoothed_conservation_vector, W, is_insertion=True)
    combined_scores = del_penalty + ins_penalty
    return np.max(np.convolve(combined_scores, np.ones(W), mode='same'))



##### Additional Insight

def compare_transcripts(reference_transcript, variant_transcript, mut):
    return {}
    # cons_seq, cons_vector = access_conservation_data(reference_transcript.transcript_id)
    # if cons_seq == reference_transcript.protein:
    #     cons_available = True
    #     cons_vector = cons_vector
    # else:
    #     cons_available = False
    #     cons_vector = [1] * len(reference_transcript.protein)

    # alignment, num_ins, num_del = get_logical_alignment(reference_transcript.protein, variant_transcript.protein)
    # deleted, inserted, aligned, unified_seq = get_insertions_and_deletions(alignment)
    # window_length = min(76, len(reference_transcript.protein))
    # cons_vector = np.array(cons_vector, dtype=float)

    # report['cons_available'] = cons_available
    # report['legacy_oncosplice_score'] = calculate_legacy_oncosplice_score(deleted, inserted, cons_vector,
    #                                                   window_length)
    # report.update(calculate_oncosplice_scores(deleted, inserted, cons_vector))
    # return pd.Series(report)

    # affected_exon, affected_intron, distance_from_5, distance_from_3 = None, None, None, None
    # for i, (ex_start, ex_end) in enumerate(reference_transcript.exons):
    #     if min(ex_start, ex_end) <= mut.start <= max(ex_start, ex_end):
    #         affected_exon = i + 1
    #         distance_from_5 = abs(mut.start - ex_start)
    #         distance_from_3 = abs(mut.start - ex_end)
    #
    # for i, (in_start, in_end) in enumerate(reference_transcript.exons):
    #     if min(in_start, in_end) <= mut.start <= max(in_start, in_end):
    #         affected_exon = i + 1
    #         distance_from_5 = abs(mut.start - in_end)
    #         distance_from_3 = abs(mut.start - in_start)

    # report = {f'reference_{k}': v for k, v in reference_transcript.constructor.items()}
    # report['reference_mRNA'] = reference_transcript.transcript_seq
    # report['reference_CDS_start'] = reference_transcript.transcript_indices.index(reference_transcript.TIS)
    # report['reference_pre_mrna'] = reference_transcript.pre_mrna
    # report['reference_ORF'] = reference_transcript.pre_mrna
    # report['reference_protein'] = reference_transcript.protein

    # report.update({f'variant_{k}': v for k, v in variant_transcript.constructor.items()})
    # report['variant_mRNA'] = variant_transcript.transcript_seq
    # report['variant_CDS_start'] = variant_transcript.transcript_indices.index(variant_transcript.TIS)
    # report['variant_pre_mrna'] = variant_transcript.pre_mrna
    # report['variant_ORF'] = variant_transcript.pre_mrna
    # report['variant_protein'] = variant_transcript.protein

    # report['protein_view'] = unified_seq
    # descriptions = define_missplicing_events(reference_transcript.exons, variant_transcript.exons,
    #                           reference_transcript.rev)
    # report['exon_changes'] = '|'.join([v for v in descriptions if v])
    # report['splicing_codes'] = summarize_missplicing_event(*descriptions)
    # report['ref_prot_length'] = len(reference_transcript.protein)
    # report['var_prot_length'] = len(variant_transcript.protein)
    # report['preservation'] = aligned/len(reference_transcript.protein)
    # report['affected_exon'] = affected_exon
    # report['affected_intron'] = affected_intron
    # report['mutation_distance_from_5'] = distance_from_5
    # report['mutation_distance_from_3'] = distance_from_3
    # report['cons_available'] = cons_available
    # report['legacy_oncosplice_score'] = calculate_legacy_oncosplice_score(deleted, inserted, cons_vector,
    #                                                   window_length)
    # report.update(calculate_oncosplice_scores(deleted, inserted, cons_vector))
    # return pd.Series(report)


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


# def get_insertions_and_deletions(alignment):
#     def find_last_real_pos(pos, pos_map):
#         if pos <= min(list(pos_map.keys())):
#             return pos_map[min(list(pos_map.keys()))]
#
#         while pos not in pos_map.keys() and pos > min(list(pos_map.keys())):
#             pos -= 1
#
#         return pos_map[pos]
#
#     ref_map = {j: p for p, j in enumerate([i for i, nt in enumerate(list(alignment.seqA)) if nt != '-'])}
#     aligned_pos, deleted_pos, inserted_pos, missmatched_pos = [], [], [], []
#     insertions, deletions = {}, {}
#     last_event = 'ALIGN'
#     del_start_pos = ''
#     unified_string = []
#     for rel_pos, (ch_a, ch_b) in enumerate(list(zip(alignment.seqA, alignment.seqB))):
#         if ch_a == ch_b != '-':
#             if last_event == 'ALIGN':
#                 pass
#             elif last_event == 'DEL':
#                 unified_string.append('<DEL_END>')
#             elif last_event == 'INS':
#                 unified_string.append('<INS_END>')
#             aligned_pos.append(rel_pos)
#             last_event = 'ALIGN'
#             unified_string.append(ch_a)
#
#         elif (ch_a != ch_b and ch_a != '-' and ch_b != '-'):
#             deleted_pos.append(rel_pos)
#             if last_event == 'DEL':
#                 deletions[del_start_pos] += ch_a
#             elif last_event == 'ALIGN':
#                 last_event = 'MM'
#                 del_start_pos = ref_map[rel_pos]
#                 deletions[del_start_pos] = ch_a
#             elif last_event == 'INS':
#                 last_event = 'MM'
#                 del_start_pos = ref_map[rel_pos]
#                 deletions[del_start_pos] = ch_a
#
#             unified_string.append(f'({ch_a}>{ch_b})')
#
#         elif (ch_a != ch_b == '-'):
#             deleted_pos.append(rel_pos)
#             if last_event == 'DEL':
#                 deletions[del_start_pos] += ch_a
#             elif last_event == 'ALIGN':
#                 last_event = 'DEL'
#                 del_start_pos = ref_map[rel_pos]
#                 deletions[del_start_pos] = ch_a
#                 unified_string.append('<DEL_START>')
#             elif last_event == 'INS':
#                 last_event = 'DEL'
#                 del_start_pos = ref_map[rel_pos]
#                 deletions[del_start_pos] = ch_a
#                 unified_string.append('<INS_END>')
#                 unified_string.append('<DEL_START>')
#             unified_string.append(ch_a)
#
#         elif ch_b != ch_a == '-':
#             point_of_insertion = find_last_real_pos(rel_pos, ref_map)
#             inserted_pos.append(point_of_insertion)
#             if last_event == 'INS':
#                 insertions[point_of_insertion] += ch_b
#             elif last_event == 'ALIGN':
#                 last_event = 'INS'
#                 insertions[point_of_insertion] = ch_b
#                 unified_string.append('<INS_START>')
#             elif last_event == 'DEL':
#                 last_event = 'INS'
#                 insertions[point_of_insertion] = ch_b
#                 unified_string.append('<DEL_END>')
#                 unified_string.append('<INS_START>')
#             unified_string.append(ch_b)
#
#     return deletions, insertions, len(aligned_pos), ''.join(unified_string)
#

# def find_unmodified_positions(lp, deletions, insertions):
#     unmodified_positions = np.ones(lp, dtype=float)
#     for pos, deletion in deletions.items():
#         unmodified_positions[pos:pos+len(deletion)] = 0
#
#     for pos, insertion in insertions.items():
#         reach = min(len(insertion) // 2, 38)
#         front_end = pos - reach
#         if front_end < 0:
#             front_end = 0
#         back_end = pos + reach + 1
#         if back_end >= len(unmodified_positions):
#             back_end = len(unmodified_positions)-1
#         len_start = pos - front_end
#         len_end = back_end - pos
#         unmodified_positions[front_end:back_end+1] = np.concatenate([np.linspace(0, 1, len_start)[::-1], np.zeros(1), np.linspace(0, 1, len_end)])
#
#     return unmodified_positions


# def calculate_oncosplice_scores(deletions, insertions, cons_vector, W=10):
#     modified_positions = 1 - find_unmodified_positions(len(cons_vector), deletions=deletions, insertions=insertions)
#     cons_vec = transform_conservation_vector(cons_vector, W=W)
#     modified_cons_vector = cons_vec * modified_positions
#     modified_cons_vector = sum_conv(modified_cons_vector, W=W) / W
#     tenth_largest_score = sorted(list(modified_cons_vector.flatten().tolist()))[-W*2]
#     max_score = max(modified_cons_vector)
#
#     if max_score > 0:
#         gof_prob = (max_score - tenth_largest_score) / max_score
#         lof_prob = 1 - gof_prob
#         return {'gof': gof_prob, 'lof': lof_prob, 'pof': 0, 'oncosplice_score': max_score}
#     else:
#         return {'gof': 0, 'lof': 0, 'pof': 1, 'oncosplice_score': max_score}
## def calculate_oncosplice_scores(deletions, insertions, cons_vec):
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
# def legacy_smooth_cons_scores(cons_scores, W):
#     return np.exp(np.negative(moving_average_conv(cons_scores, W, 2)))

# def calculate_del_penalty(deleted_domains, cons_scores, W):
#     penalty = np.zeros(len(cons_scores))
#     for dp_pos, dp_seq in deleted_domains.items():
#         dw = max(1.0, len(dp_seq) / W)
#         penalty[dp_pos:dp_pos + len(dp_seq)] = cons_scores[dp_pos:dp_pos + len(dp_seq)] * dw
#     return penalty

# def calculate_ins_penalty(inserted_domains, cons_scores, W):
    # penalty = np.zeros(len(cons_scores)) #cons_scores.copy()
    # for ip_pos, ip_seq in inserted_domains.items():
    #     reach = min(W//2, len(ip_seq)//2)
    #     iw = max(1.0, len(ip_seq) / W)
    #     penalty[ip_pos-reach:ip_pos+reach] = iw*cons_scores[ip_pos-reach:ip_pos+reach]
    # return penalty

# def combine_ins_and_del_scores(d_cons_scores, i_cons_scores, W):
#     combined_scores = d_cons_scores + i_cons_scores
#     penalty = np.convolve(combined_scores, np.ones(W), mode='same')
#     return max(penalty)
# def calculate_legacy_oncosplice_score(deletions, insertions, cons_vec, W):
#     smoothed_conservation_vector = legacy_smooth_cons_scores(cons_vec, W)
#     deconv_del = calculate_del_penalty(deletions, smoothed_conservation_vector, W)
#     deconv_ins = calculate_ins_penalty(insertions, smoothed_conservation_vector, W)
#     return combine_ins_and_del_scores(deconv_del, deconv_ins, W)