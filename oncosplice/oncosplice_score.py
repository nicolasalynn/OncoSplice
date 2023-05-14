import numpy as np

def moving_average_conv(vector, W):
    convolving_length = np.array([min(len(vector) + W - i, W, i) for i in range(W // 2, len(vector) + W // 2)], dtype=float)
    return sum_conv(vector, W) / convolving_length
def moving_average_conv_modified(vector, W):
    convolving_length = np.array([min(len(vector) + W - i, W, i) for i in range(W // 2, len(vector) + W // 2)], dtype=float)
    return sum_conv(vector, W) / (convolving_length / 2)
def sum_conv(vector, W):
    return np.convolve(vector, np.ones(W), mode='same')

def transform_conservation_vector(c, W=9):
    c[c<0] /= abs(min(c))                           # normalizings
    c[c>0] /= max(c)
    c = np.exp(np.negative(sum_conv(c, W)))                      # smoothed and inverted evolutionary rate values; exponential polarizes values
    return c * len(c) / sum(c)                      # normalize so sum of conservation values is equal to the length of the protein, this can also be done with some arbitrary value such as 1000

def find_unmodified_positions(lp, deletions, insertions, W):
    unmodified_positions = np.ones(lp, dtype=float)
    for pos, deletion in deletions.items():
        unmodified_positions[pos:pos+len(deletion)] = 0

    max_reach = W // 2
    for pos, insertion in insertions.items():
        reach = min(len(insertion) // 2, max_reach)
        unmodified_positions[pos-reach:pos+reach+1] = 0

    return unmodified_positions
def calculate_oncosplice_scores(deletions, insertions, cons_vec, W):
    cons_vec = transform_conservation_vector(cons_vec)
    unmodified_positions = find_unmodified_positions(len(cons_vec), deletions=deletions, insertions=insertions, W=W)
    alignment_ratio_vector = moving_average_conv_modified(unmodified_positions, W) - 1
    functional_loss_vector = cons_vec * (1 - unmodified_positions)
    s = alignment_ratio_vector * functional_loss_vector / len(cons_vec)
    # s = moving_average_conv_modified(s, W)
    stemp = abs(s)
    return {'cons_vec': np.array2string(np.around(cons_vec), 3), 'lof_score': abs(min(0, s.min())), 'gof_score': max(0, s.max()), 'oncosplice_score': sum(stemp)/len(cons_vec)}


##### LEGACY ONCOSPLICE CALCS
def legacy_smooth_cons_scores(cons_scores, W):
    c = np.exp(np.negative(moving_average_conv_modified(cons_scores, W)))
    return c - c.min()

def calculate_del_penalty(deleted_domains, cons_scores, W):
    penalty = np.zeros(len(cons_scores)) #cons_scores.copy()
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
    return {'legacy_oncosplice_score': combine_ins_and_del_scores(deconv_del, deconv_ins, W)}
