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
def calculate_oncosplice_scores(deletions, insertions, cons_vec, W):
    unmodified_positions = find_unmodified_positions(len(cons_vec), deletions=deletions, insertions=insertions)

    functional_loss_vector_5 = transform_conservation_vector(cons_vec, W=5) * (1 - unmodified_positions)
    functional_loss_vector_5 = sum_conv(functional_loss_vector_5, W=5)

    functional_loss_vector_76 = transform_conservation_vector(cons_vec, W=76) * (1 - unmodified_positions)
    functional_loss_vector_76 = sum_conv(functional_loss_vector_76, W=76)

    return {'cons_vec': np.array2string(np.around(cons_vec), 3), 'oncosplice_score_lof': max(functional_loss_vector_76), 'oncosplice_score_gof': max(functional_loss_vector_5)}


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
    return {'legacy_oncosplice_score': combine_ins_and_del_scores(deconv_del, deconv_ins, W)}
