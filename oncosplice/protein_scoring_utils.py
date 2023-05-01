import numpy as np
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re


# import termplotlib as tpl


def generate_report(ref_prot, var_prot, missplicing, mutations):
    '''
    Report 1: The oncosplice score
    Report 2: The insertion score
    Report 3: The deletion score
    Report 4: General Information
            a. gene function
            b. gene association to oncogene or TSG
            c. Number of affected transcripts
            d. Verbal splicing description
            e. Presence of mutation in dbSNP, ClinVar, and 1K Genome Project

    Report 4: Per transcript isoform report -
            a. Protein alignment
            b. Novel peptides
            c. Number of positions conserved
            d. Number of added amino acids
            f. Number of deleted amino acids
            e. Verbal description of event (can take one of several values such as
                        - frame shift and early termination where first {x} / {n} positions of the protein are conserved.
                        - generally conserved with a novel {x} amino acid long peptide interrputing position {x} of {k}.
                        - protein unchanged
            f. Verbal description at level of splice sites/exons
            g. If available, include tissue expression data of transcript
    '''
    sample = list(ref_prot.values())[0]
    full_report, transcript_comprisons = {'gene_name': sample.gene_name}, {}
    del_score, ins_score, full_score = [], [], []
    for rtid, ref_protein in ref_prot.items():
        transcript_comprisons_isoforms = {}
        temp_del, temp_ins, temp_full = [], [], []
        for vtid, var_protein in var_prot.items():
            if var_protein.transcript_id.split('-')[0] == rtid:
                transcript_comprisons_isoforms[vtid] = protein_change_report(ref_protein.protein, var_protein.protein,
                                                                             ref_protein.conservation_vector, 76)
                transcript_comprisons_isoforms[vtid]['verbal_description'] = describe_general_change(
                    transcript_comprisons_isoforms[vtid])
                transcript_comprisons_isoforms[vtid]['splice_site_difference'] = describe_splice_site_difference(
                    ref_protein, var_protein)
                temp_del.append(transcript_comprisons_isoforms[vtid]['deletion_score'])
                temp_ins.append(transcript_comprisons_isoforms[vtid]['insertion_score'])
                temp_full.append(transcript_comprisons_isoforms[vtid]['score'])
        del_score.append(sum(temp_del) / len(temp_del))
        ins_score.append(sum(temp_ins) / len(temp_ins))
        full_score.append(sum(temp_full) / len(temp_full))

        transcript_comprisons[rtid] = transcript_comprisons_isoforms

    # full_report['gene_description'] = query_gene_description(ref_prot.gene_name)
    # full_report['gene_cancer_role'] = query_gene_cancer_associations(ref_prot.gene_name)
    # full_report['mutation_crossref'] = {m: query_mutation_db(m) for m in mutations}
    full_report['mutations'] = ', '.join(mutations)
    full_report['missplicing'] = missplicing
    full_report['scores'] = {
        'score': max(full_score),
        'deletion_score': max(del_score),
        'insertion_score': max(ins_score)
    }

    full_report['transcript_comparison'] = transcript_comprisons

    return full_report


def query_gene_description(gene):
    pass


def query_gene_cancer_associations(gene):
    pass


def query_mutation_db(mut_id):
    pass


def describe_general_change(data):
    print(data)
    log = ''
    if data['score'] == 0:
        return 'Protein is conserved. No changes observed.'

    if data['deletion_score'] > 100 or data['aligned_positions'] / data['reference_protein_length'] < 0.75:
        log += 'Deletion score indicates some loss of function. '

    if data['insertion_score'] > 0 and data['aligned_positions'] / data['reference_protein_length'] > 0.90:
        log += 'Insertion score and general conservation indicates a novel peptide that may modify protein function.'

    return log


def describe_splice_site_difference(ref_prot, var_prot):
    zipped_reference = list(
        zip([ref_prot.transcript_start] + ref_prot.acceptors, ref_prot.donors + [ref_prot.transcript_end]))
    zipped_variant = list(
        zip([var_prot.transcript_start] + var_prot.acceptors, var_prot.donors + [var_prot.transcript_end]))
    i = 0
    for i, ref_exon in enumerate(zipped_reference):
        var_exon = zipped_variant[i]
        if ref_exon[0] != var_exon[0] and ref_exon[1] != var_exon[1]:
            if ref_exon in zipped_variant[i:] and var_exon not in zipped_reference:
                return i + 1, f'Novel exon in exon {i}/{len(zipped_reference)}: {var_exon}'
            else:
                return i + 1, f'Exon skipped in exon {i}/{len(zipped_reference)}: {ref_exon}'

        elif ref_exon[0] != var_exon[0] and ref_exon[1] == var_exon[1]:
            return i + 1, f'Using alternative acceptor in exon {i}/{len(zipped_reference)}: {ref_exon[0]} --> {var_exon[0]}'

        elif ref_exon[0] == var_exon[0] and ref_exon[1] != var_exon[1]:
            if var_exon[1] in [val[1] for val in zipped_reference[i + 1:]]:
                return i + 1, f'Intron retention in exon {i}/{len(zipped_reference)}: {ref_exon} --> {var_exon}'
            else:
                return i + 1, f'Using alternative donor in exon {i}/{len(zipped_reference)}: {ref_exon[1]} --> {var_exon[1]}'

    return 'No splicing difference detected.'


def protein_change_report(r, v, cons_scores, W=76):
    """
        Purpose: a function that allows us to get a ngenome_general comparison in protein functionality. The method can
                 operate simply via alignment, but it is able to use conservation estimates (conservation values for
                 each nt in some reference sequence) to make its quantitative measurement more effective.

        Input:   1. r(eference_seq) - some reference sequence (we call it reference in the case that we have
                                       conservation data for it, otherwise it is simply an arbitary sequence in the
                                       comparison).
                 2. v(ariant_seq)    - some variant sequence(we call it variant in the case that we think of this
                                       sequence as some isoform of the reference sequence, otherwise it can be an
                                       arbitary sequence).
                 3. cons_scores      - an array of conservation estimates for each nucleotide in the reference sequence.
                                       In the case that conservation values are unavailable or unknown, it need not be
                                       specified and the function will treat each nt position as equally weighted.

        Output:  1. score - a value in the range of [0, 1] where 1 refers to complete preservation between the r and v
                            and 0 refers to complete uniqueness between the two sequences.
                 2. description - some mechanical descriptor of the differences of the two sequences.
                 3. reference_protein_alignment -
                 4. variant_protein_alignment -
    """

    def find_last_real_pos(pos, pos_map):
        if pos <= min(list(pos_map.keys())):
            return pos_map[min(list(pos_map.keys()))]

        while pos not in pos_map.keys() and pos > min(list(pos_map.keys())):
            pos -= 1

        return pos_map[pos]

    assert len(cons_scores) == len(
        r), 'Length of conservation vector provided does not match the length of the reference protein.'
    assert r, 'Reference protein must consist of amino acids. Protein passed is empty.'
    if not v:
        report = {
            'score': sum(cons_scores),
            'deletion_score': sum(cons_scores),
            'insertion_score': 0,
            'reference_protein_length': len(r),
            'variant_protein_length': 0,
            'reference_protein_alignment': '-',
            'variant_protein_alignment': '-',
            'insertions': {},
            'deletions': {0: r},
            'aligned_positions': 0,
        }
        return report

    if W > len(r) // 2:
        W = len(r) // 2

    if W % 2 != 0:
        W -= 1

    PADDING = W // 2

    best_alignment, gaps = get_logical_alignment(r, v)
    print(format_alignment(*best_alignment))

    ref_map = {j: p for p, j in enumerate([i for i, nt in enumerate(list(best_alignment.seqA)) if nt != '-'])}
    # ref_map_inverted = {p: j for j, p in ref_map.items()}
    aligned_pos, deleted_pos, inserted_pos = [], [], []
    insertions, deletions = {}, {}
    last_event = 'ALIGN'
    del_start_pos = ''
    for rel_pos, (ch_a, ch_b) in enumerate(list(zip(best_alignment.seqA, best_alignment.seqB))):
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

    bar, bav = str(best_alignment.seqA), str(best_alignment.seqB)
    score, del_score, ins_score = score_alignment(deletions, insertions, cons_scores, W, PADDING)
    if ins_score > del_score:
        print(f"CAPTURED INTERESTING EVENT!!!!\n{bar}\n{bav}")

    report = {
        'score': score,
        'deletion_score': del_score,
        'insertion_score': ins_score,
        'reference_protein_length': len(r),
        'variant_protein_length': len(v),
        'reference_protein_alignment': bar,
        'variant_protein_alignment': bav,
        'insertions': {f'{i}/{len(r)}': seq for i, seq in insertions.items()},
        'deletions': deletions,
        'aligned_positions': best_alignment.score,
    }

    return report


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
        return alignments[block_counts.index(min(block_counts))], min(block_counts)
    except AttributeError:
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


def score_alignment(deletions, insertions, cons_scores, W, PADDING):
    cons_scores = smooth_cons_scores(cons_scores, W, PADDING)
    del_penalty = calculate_del_penalty(deletions, cons_scores, W, PADDING)
    ins_penalty = calculate_ins_penalty(insertions, cons_scores, W, PADDING)
    return combine_ins_and_del_scores(del_penalty, ins_penalty, W, PADDING), combine_ins_and_del_scores(del_penalty,
                                                                                                        del_penalty, W,
                                                                                                        PADDING), combine_ins_and_del_scores(
        ins_penalty, ins_penalty, W, PADDING)



