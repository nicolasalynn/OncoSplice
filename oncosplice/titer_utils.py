# should these be not in this file? or in a flag of whether to import?

# from keras.constraints import maxnorm
# from keras.layers import Conv1D, MaxPool1D, LSTM, Dropout, Flatten, Dense, Activation
### from keras.layers import GRU
# from keras import Sequential, Input
### # from sklearn.metrics import average_precision_score, roc_auc_score
import numpy as np
import os
import re
import pickle
from scipy.stats import percentileofscore
from Bio import pairwise2
from keras.constraints import maxnorm
from keras.layers import Conv1D, MaxPool1D, LSTM, Dropout, Flatten, Dense, Activation
from keras import Sequential, Input

def run_through_titer_adaptor(ref_cds, var_cds):
    new_cds_start, titer_score, _, _, _ = run_through_titer(mut_seq=var_cds.mature_mrna,
                                                            mut_coords=var_cds.mature_indices,
                                                            ref_seq=ref_cds.mature_mrna,
                                                            ref_coords=ref_cds.mature_indices,
                                                            ref_sc_coord=ref_cds.used_tis,
                                                            ref_tts_coord=ref_cds.used_tts,
                                                            ref_id=ref_cds.transcript_id,
                                                            data_path='/tamir1/shaicohen1/share/titer-master',
                                                            titer_model=titer_model,
                                                            all_test_seqs={})
    if var_cds.use_tis != new_cds_start:
        var_cds.use_tis = new_cds_start
        var_cds.generation_report += f'TITER predicts change in TIS site: {new_cds_start}, score: {titer_score}'

    return None

def run_through_titer(mut_seq, mut_coords, ref_sc_coord, ref_tts_coord, ref_id, data_path,
                      titer_model, all_test_seqs):
    '''

    :param mut_seq: str - nts of mutated mRNA
    :param mut_coords: list of ints - coords to each nt in mut_seq
    :param ref_coords: str - nts of ref mRNA
    :param ref_seq: list of ints - coords to each nt in ref_seq
    :param ref_sc_coord: int - coord of ref start codon
    :param ref_sc_coord: int - coord of ref stop codon
    :param ref_id: str - id of transcript ("ENST....")
    :param data_path: str - path to titer-master folder (i.e. where the folders codes,data and model exist)
    :return:
        - suitable start codon coord (can be the same as the original)
        - TITER score of suitable start codon
    '''
    suitable_score_coord, suitible_score, suitable_score_percentile = -1, -np.Inf, 0.0
    data_file = os.path.join(data_path, 'codes', 'titer_canon_sc_scores_sorted.pickle')
    # print(data_file)
    with open(data_file, 'rb') as f:
        sc_table = pickle.load(f)

    if ~np.any(sc_table['transcript_id'] == ref_id):
        print(f'Did not find {ref_id} in start codon database!')
        return -1, -np.Inf, 0.0, 0.0, all_test_seqs

    ref_sc_table = sc_table.loc[sc_table['transcript_id'] == ref_id]
    ref_sc_context_seq = ref_sc_table.loc[ref_sc_table.index[0], 'sc_context_seq']
    ref_sc_context_score = ref_sc_table.loc[ref_sc_table.index[0], 'can_sc_score']

    if ref_sc_coord in mut_coords:
        sc_idx = mut_coords.index(ref_sc_coord)
        sc_context_seq, _ = get_sc_context_seq(sc_idx, mut_seq,
                                               mut_coords)  # idx refers to the 0:len(seq) of the sequence
        if sc_context_seq == ref_sc_context_seq:
            print('in titer utils - outputting start codon context is conserved')
            return ref_sc_coord, ref_sc_context_score, percentileofscore(sc_table['can_sc_score'].to_list(),
                                                                         ref_sc_context_score), percentileofscore(
                sc_table['can_sc_score'].to_list(),
                ref_sc_context_score), all_test_seqs  # start codon and context were preserved

    # if we got here it means start codon\context weren't preserved

    sc_align_context = 25
    # seq_to_align = ''.join([x for x in ref_sc_context_seq[100 - sc_align_context:103 + sc_align_context].split('$') if x!=''])
    best_al = pairwise2.align.localms(ref_sc_context_seq[100 - sc_align_context:103 + sc_align_context], mut_seq, 2, -1,
                                      -.5, -.1, one_alignment_only=True)
    mut_near_sc_idx = best_al[0].start + sc_align_context

    # start looking at all sc candidates in a window of 20 centered around mut_near_sc_pos (check 21 positions)
    # and expand the search until finding a suitable start codon or until reaching 400 nts around the start codon
    ref_score_percentile = percentileofscore(sc_table['can_sc_score'].to_list(), ref_sc_context_score)
    window_400 = [[max(0, mut_near_sc_idx - 400), min(mut_near_sc_idx + 400, len(mut_seq))]]
    sc_coords, cand_sc_idxs, scores, all_test_seqs = pred_trans_init_TITER([mut_seq], [mut_coords], [ref_tts_coord],
                                                                           window_400, data_path,
                                                                           titer_model=titer_model,
                                                                           all_test_seqs=all_test_seqs)
    windows = [10, 50, 100, 200, 400]
    for w_dist in windows:
        curr_sc_idxs = [i for i in range(len(scores)) if abs(mut_near_sc_idx - cand_sc_idxs[i]) <= w_dist]
        if curr_sc_idxs != []:
            curr_sc_scores = [scores[scores_idx][0] for scores_idx in curr_sc_idxs]
            best_sc_score_in_curr_window = max(curr_sc_scores)
            best_sc_coord_in_curr_window = sc_coords[curr_sc_idxs[curr_sc_scores.index(best_sc_score_in_curr_window)]]
            if best_sc_score_in_curr_window > suitible_score:
                suitable_score_coord = best_sc_coord_in_curr_window
                suitible_score = best_sc_score_in_curr_window
                suitable_score_percentile = percentileofscore(sc_table['can_sc_score'].to_list(), suitible_score)

            if (suitible_score >= ref_sc_context_score) or ((ref_score_percentile - suitable_score_percentile) <= 5):
                return suitable_score_coord, suitible_score, suitable_score_percentile, ref_score_percentile, all_test_seqs
    print(f"Used titer to find alternative start site: {suitable_score_coord}")
    # no suitable start codon was found (only if we searched 400 nts up\downstream from original start codon and didn't find any suitable start codon)  - so return the start codon that got the best score..
    return suitable_score_coord, suitible_score, suitable_score_percentile, ref_score_percentile, all_test_seqs


def pred_trans_init_TITER(mRNA_seqs, genomic_positions, ref_stop_codon_coord, windows_to_analyze, data_path,
                          titer_model=[], NNN_to_analyze='NTG_ACG_only', use_prior=True, all_test_seqs={}):
    # '''
    #    # INPUT:
    #    # :param mRNA_seqs: list of strings - each an mRNA to be analyzed (to find start codons)
    #    # :param genomic_positions: list of lists of ints - for each mRNA a corresponding list of genomic positions
    #    # :param ref_stop_codon_coord: list of ints - for each mRNA a corresponding stop codon position
    #    # :param windows_to_analyze: list of lists of 2 ints - for each mRNA a corresponding list of two integers - one denoting the start of the window, and one the end of it
    #####IMPORTANT!!! the ints in the lists in windows_to_analyze are the indices of the seq not in the genome but relative to the mRNA (i.e. from 0 to len(mRNA))
    #####IMPORTANT!!! the 2nd int in the lists in windows_to_analyze will be analyzed i.e. the function will take an extra 2 nts after it so that it can be analyzed
    #    # :param data_path: path to titer-master folder (i.e. where the folders codes,data and model exist)
    #    # :param scan_NTG_ACG_only: bool for whether to only scan for NTG (A\C\T\G GT) or ACG. Can also be ATG_only and All_Codons
    #    # :param use_prior: bool for whether to use prior distribiution of start codons (they did in the article)
    #    # :return:
    #    # scores - list of lists of floats - scores of each candidate start codon in each mRNA - [[m1sc1, m1sc2],[m2sc1],...]
    #    # start_codon_positions - list of lists of ints - genomic positions of each candidate start codon in each mRNA - [[m1pos1, m1pos2],[m2pos1],...]
    # '''
    start_codon_positions, sc_idxs = [], []
    if titer_model == []:
        titer_model = build_titer_model()
    codon_tis_prior = np.load(os.path.join(data_path, r'codes/dict_piror_front_Gaotrain.npy'), allow_pickle=True)
    codon_tis_prior = codon_tis_prior.item()
    test_seqs = []
    test_priors = []
    if use_prior and NNN_to_analyze == 'All_Codons':
        for c in codon_tis_prior.keys():
            if codon_tis_prior[c] == 'never':
                codon_tis_prior[c] = -np.Inf

    if NNN_to_analyze == 'ATG_only':
        codon_list = ['ATG']
    elif NNN_to_analyze == 'NTG_ACG_only':
        ## how they found the codons to anaylze in the example they attached..
        # for c in codon_tis_prior.keys():
        #     if codon_tis_prior[c] != 'never' and codon_tis_prior[c] >= 1:
        #         codon_list.append(c)
        codon_list = ['ATG', 'CTG', 'ACG', 'TTG', 'GTG']

    elif NNN_to_analyze == 'All_Codons':
        codon_list = list(codon_tis_prior.keys())  # all codons
    else:
        print(r'NNN_to_analyze takes only ATG_only\NTG_ACG_only\All_Codons')
        return [], [], [], []

    # get all seqs surrounding codons from codon_list to be analyzed from each mRNA
    for s_mRNA_idx in range(len(mRNA_seqs)):
        s_mRNA = mRNA_seqs[s_mRNA_idx]
        # start_codon_positions.append([])
        # sc_idxs.append([])
        curr_window = windows_to_analyze[s_mRNA_idx]
        if len(curr_window) != 2:
            s_mRNA_subseq = s_mRNA
            curr_window = [0, 0]
        else:
            s_mRNA_subseq = s_mRNA[curr_window[0]:min(curr_window[1] + 3,
                                                      len(s_mRNA))]  # curr_window[1]+3 because we want to analyze every nt in the window including the last 2, so need to pad

        for c in codon_list:
            if use_prior:
                curr_codon_prior = codon_tis_prior[c]
            else:
                curr_codon_prior = 1  # the prediction for a start codon is the prior for the codon times the prob from the Neural Network

            c_idxs = [m.start() + curr_window[0] for m in re.finditer(c, s_mRNA_subseq)]
            # take only codons that are in-frame with the original protein
            c_idxs = [x for x in c_idxs if
                      (genomic_positions[s_mRNA_idx].index(ref_stop_codon_coord[s_mRNA_idx]) + 3 - x) % 3 == 0]
            for curr_sc_idx in c_idxs:
                curr_sc_tss_seq = s_mRNA[
                                  curr_sc_idx:genomic_positions[s_mRNA_idx].index(ref_stop_codon_coord[s_mRNA_idx]) + 3]

                curr_stop_codons = [[i, curr_sc_tss_seq[i:i + 3]] for i in range(0, len(curr_sc_tss_seq), 3) if
                                    curr_sc_tss_seq[i:i + 3] in ['TGA', 'TAG', 'TAA']]
                if len(curr_stop_codons) > 1:  # if there's another stop codon upstream to the original one then this candidate is excluded because we want proteins that have a chance of being similar to the original one
                    continue
                curr_sc_context_start = max(curr_sc_idx - 100, 0)
                curr_sc_context_stop = min(curr_sc_idx + 103, len(s_mRNA))
                curr_test_seq = s_mRNA[
                                curr_sc_context_start:curr_sc_context_stop]  # start codon needs to be in 100:103 and the length should be 203
                if len(curr_test_seq) < 203:
                    curr_test_seq = '$' * (max(100 - curr_sc_idx, 0)) + curr_test_seq + '$' * max(
                        (curr_sc_idx + 103) - len(s_mRNA), 0)
                test_seqs.append(curr_test_seq)
                test_priors.append(curr_codon_prior)
                curr_genomic_pos = genomic_positions[s_mRNA_idx][curr_sc_idx]
                assert (len(curr_test_seq) == 203)
                assert (curr_test_seq[100:103] == c)
                start_codon_positions.append(int(curr_genomic_pos))
                sc_idxs.append(curr_sc_idx)

    scores = np.zeros((len(test_seqs), 1))
    if len(test_seqs) == 0:
        print(f'no possible start codons to analyze in this window! window: {windows_to_analyze}')
        return start_codon_positions, sc_idxs, scores, all_test_seqs

    test_priors_to_analyze = []
    test_seqs_to_analyze = []
    idxs_analyzed_scores = []
    for i in range(len(test_seqs)):
        if test_seqs[i] in all_test_seqs.keys():
            scores[i] = all_test_seqs[test_seqs[i]]
        else:
            idxs_analyzed_scores.append(i)
            test_seqs_to_analyze.append(test_seqs[i])
            test_priors_to_analyze.append(test_priors[i])

    # breakpoint()
    print(f'{len(test_seqs)} test seqs found in window...')
    if len(test_seqs_to_analyze) == 0:
        print('Predictions for all test seqs for this window have already been acquired previously!\n not running NN!')
        return start_codon_positions, sc_idxs, scores, all_test_seqs

    print(
        f'Predicting on {len(test_seqs_to_analyze)} test seqs ({len(test_seqs) - len(test_seqs_to_analyze)} already analyzed in other isoforms or transcripts)...')

    processed_test_seqs = seq_matrix(test_seqs_to_analyze)
    test_priors_to_analyze = np.array(test_priors_to_analyze).reshape(len(test_priors_to_analyze), 1)
    analyzed_scores = np.zeros((len(test_seqs_to_analyze), 1))
    for i in range(32):
        titer_model.load_weights(os.path.join(data_path, r"model/bestmodel_" + str(i) + ".hdf5"))
        y_test_pred = titer_model.predict(processed_test_seqs, verbose=0)
        analyzed_scores += y_test_pred * test_priors_to_analyze

    for i in range(len(analyzed_scores)):
        scores[idxs_analyzed_scores[i]] = analyzed_scores[i]
        all_test_seqs.update({test_seqs_to_analyze[i]: analyzed_scores[i]})

    return start_codon_positions, sc_idxs, scores, all_test_seqs

    # count = 0
    # scores_organized = [[] for x in start_codon_positions]
    # for imRNA in range(len(start_codon_positions)):
    #    if len(start_codon_positions[imRNA]) == 0:
    #        continue
    #    scores_organized[imRNA] = [x[0] for x in scores[count:count + len(start_codon_positions[imRNA])]]
    #    count += len(start_codon_positions[imRNA])
    # return start_codon_positions, sc_idxs, scores_organized, test_seqs


def get_sc_context_seq(idx, seq, seq_coords):
    sc_context_start_idx = max(idx - 100, 0)
    sc_context_stop_idx = min(idx + 103, len(seq))
    sc_context_coords = seq_coords[sc_context_start_idx:sc_context_stop_idx]
    sc_context_seq = seq[
                     sc_context_start_idx:sc_context_stop_idx]  # start codon needs to be in 100:103 and the length should be 203
    if len(sc_context_seq) < 203:
        sc_context_seq = '$' * (max(100 - idx, 0)) + sc_context_seq + '$' * max(
            (idx + 103) - len(sc_context_seq), 0)
    return sc_context_seq.upper(), sc_context_coords


def seq_matrix(seq_list):
    tensor = np.zeros((len(seq_list), 203, 8))
    for i in range(len(seq_list)):
        seq = seq_list[i]
        j = 0
        for s in seq:
            if s == 'A' and (j < 100 or j > 102):
                tensor[i][j] = [1, 0, 0, 0, 0, 0, 0, 0]
            if s == 'T' and (j < 100 or j > 102):
                tensor[i][j] = [0, 1, 0, 0, 0, 0, 0, 0]
            if s == 'C' and (j < 100 or j > 102):
                tensor[i][j] = [0, 0, 1, 0, 0, 0, 0, 0]
            if s == 'G' and (j < 100 or j > 102):
                tensor[i][j] = [0, 0, 0, 1, 0, 0, 0, 0]
            if s == '$':
                tensor[i][j] = [0, 0, 0, 0, 0, 0, 0, 0]
            if s == 'A' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 1, 0, 0, 0]
            if s == 'T' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 0, 1, 0, 0]
            if s == 'C' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 0, 0, 1, 0]
            if s == 'G' and (j >= 100 and j <= 102):
                tensor[i][j] = [0, 0, 0, 0, 0, 0, 0, 1]
            j += 1
    return tensor


def build_titer_model():
    print('Building model...')
    # from keras.constraints import maxnorm
    # from keras.layers import Conv1D, MaxPool1D, LSTM, Dropout, Flatten, Dense, Activation
    # # from keras.layers import GRU
    # from keras import Sequential, Input
    # from keras.constraints import maxnorm
    # from keras.layers import Conv1D, MaxPool1D, LSTM, Dropout, Flatten, Dense, Activation
    # # from keras.layers import GRU
    # from keras import Sequential, Input
    model = Sequential()
    model.add(Input(shape=(203, 8)))
    model.add(Conv1D(filters=128,
                     kernel_size=3,
                     padding='valid',
                     kernel_constraint=maxnorm(3),
                     activation='relu'))  # ,
    # subsample_length=1)) #this was in the older implementation in the article with an older version of keras..How do we replace this? do we need to?

    model.add(MaxPool1D(3))
    model.add(Dropout(rate=0.21370950078747658))
    model.add(LSTM(units=256,
                   return_sequences=True))
    model.add(Dropout(rate=0.7238091317104384))
    model.add(Flatten())
    model.add(Dense(1))
    model.add(Activation('sigmoid'))

    model.compile(loss='binary_crossentropy',
                  optimizer='nadam',
                  metrics=['accuracy'])
    return model


if __name__ == '__main__':
    print("only utils")
    titer_model = build_titer_model()
else:
    # from time import time
    titer_model = build_titer_model()

################### DEBUGGING #############################################################
# def get_mRNA(exon_starts, exon_ends, strand, chr_num):
#     mRNA = ''
#     genomic_idxs = []
#     for idx in range(len(exon_starts)):
#         curr_exon_end = max(exon_starts[idx], exon_ends[idx])
#         curr_exon_start = min(exon_starts[idx], exon_ends[idx])
#         curr_exon = chrs['chr' + chr_num][curr_exon_start - 1:curr_exon_end].seq.upper()
#         curr_genomic_idxs = list(range(curr_exon_start, curr_exon_end + 1))  # 1 based - so we add +1
#         if strand == '-':
#             curr_exon = reverse_complement(curr_exon)
#             curr_genomic_idxs.reverse()
#         mRNA = mRNA + curr_exon
#         genomic_idxs.extend(curr_genomic_idxs)
#     return mRNA, genomic_idxs
#
# from pyfaidx import Fasta
# import pandas as pd
# import numpy as np
# import sys
# sys.path.append(r"C:\Users\shaic\PycharmProjects\EXPosition- for shaked\Utils")
# from translation import reverse_complement
# import pickle
# import json
# with open('./titer_can_sc_scores_0_7500.pickle','rb') as f:
#     sc_table = pickle.load(f)
# file = r"C:\Users\shaic\PycharmProjects\titer-master\codes\CYB5RL_annotation.json"
# with open(file) as f:
#     gene = json.load(f)
# ref_id = gene['transcripts'][0]['transcript_id']
# ref_sc_coord = gene['transcripts'][0]['CDS_start'][0]
# exon_starts = gene['transcripts'][0]['exon_start']
# exon_ends = gene['transcripts'][0]['exon_end']
# strand = ['-' if gene['rev'] == True else '+'][0]
# chr_num = gene['chrm']
# fasta_name = r"./hg38.fa"
# chrs = Fasta(fasta_name)
# ref_seq, ref_coords = get_mRNA(exon_starts, exon_ends, strand, chr_num)
# mut_seq, mut_coords = get_mRNA(exon_starts, exon_ends, strand, chr_num)
# SNP
# mut_seq= mut_seq[:ref_coords.index(ref_sc_coord)] + 'G' + mut_seq[ref_coords.index(ref_sc_coord)+1:]

# INS
# mut_seq= mut_seq[:ref_coords.index(ref_sc_coord)+1] + 'G' + mut_seq[ref_coords.index(ref_sc_coord)+1:]
# mut_coords.insert(ref_coords.index(ref_sc_coord)+1,ref_sc_coord+.0001)

# DEL
# mut_seq= mut_seq[:ref_coords.index(ref_sc_coord)-10] + mut_seq[ref_coords.index(ref_sc_coord)+10:] + 'G'*21
# mut_coords[mut_coords.index(ref_sc_coord)-10:mut_coords.index(ref_sc_coord)+9] = [0]*20

# suit_sc_coord, suit_sc_score = run_through_titer(mut_seq, mut_coords,
#                                                  ref_coords, ref_seq,
#                                                  ref_sc_coord, ref_id)
# breakpoint()

