from param import *
from Bayes import bayes_evidence, merge_mean
from models import C3F2
import numpy as np
import torch
import os
import math


def load_models(models_path, priors_file, dev='cpu'):
    mdl_list = []
    priors_list = []
    for cv in range(4):
        model_state_dict = torch.load(f'{models_path}MM{cv}.cnn', weights_only=True)
        mdl = C3F2()
        mdl.load_state_dict(model_state_dict)
        mdl_list.append(mdl.to(dev))
    with open(priors_file, 'r') as fin:
        for line in fin:
            if len(line) < 3:
                continue
            lst = line.strip().split()
            priors_list.append(float(lst[1]))
    return mdl_list, priors_list


def create_out_dir():
    if not os.path.exists(output_path):
        os.system(f"mkdir {output_path}")


def read_p_matrix_dict1():
    ret = {}
    with open(p_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) <= 1:
                continue
            lst = line.split()
            if lst[0] == '.':
                aa_list = lst[1:]
                for aa in aa_list:
                    ret[aa] = []
            else:
                aa = lst[0]
                ret[aa] = np.array(lst[1:], dtype='float')
    return ret


def read_f5_dict():
    ret = {}
    with open(c_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            lst = line.split()
            if lst[0] not in ret:
                ret[lst[0]] = {}
            ret[lst[0]][lst[1]] = float(lst[2]) / 2.0
    return ret


def get_features_p(seq, p_dict, w1, w2, skip=0):
    w = w1 + w2 + skip
    aa_map = {'-': 20}
    features_p = {'w1': np.zeros(len(seq), dtype='float'), 'w2': np.zeros(len(seq), dtype='float')}
    for ii, aa in enumerate(p_dict):
        aa_map[aa] = ii
    _pad = '-' * w
    p2_seq = list(_pad + seq + _pad)
    for ii, aa in enumerate(p2_seq):
        if aa not in aa_map:
            p2_seq[ii] = '-'
    aa_count = {'w1': np.zeros(21, dtype='int32'), 'w2': np.zeros(21, dtype='int32')}
    aa_count['w1'][20] = w1
    aa_count['w2'][20] = w2
    for ii in range(w1):
        aai = aa_map[p2_seq[w + ii + 1 + skip]]
        aa_count['w1'][aai] += 1
    for ii in range(w1, w1 + w2):
        aai = aa_map[p2_seq[w + ii + 1 + skip]]
        aa_count['w2'][aai] += 1

    features_p['w1'][0] = np.dot(aa_count['w1'][:20], p_dict[seq[0]])
    features_p['w2'][0] = np.dot(aa_count['w2'][:20], p_dict[seq[0]])

    for i in range(1, len(seq)):
        ii = i + w
        aai = aa_map[p2_seq[ii - w1 - 1 - skip]]
        aa_count['w1'][aai] -= 1
        aai = aa_map[p2_seq[ii + skip]]
        aa_count['w1'][aai] -= 1
        aai = aa_map[p2_seq[ii - 1 - skip]]
        aa_count['w1'][aai] += 1
        aai = aa_map[p2_seq[ii + w1 + skip]]
        aa_count['w1'][aai] += 1

        aai = aa_map[p2_seq[ii - w - skip]]
        aa_count['w2'][aai] -= 1
        aai = aa_map[p2_seq[ii + w1 + skip]]
        aa_count['w2'][aai] -= 1

        aai = aa_map[p2_seq[ii - w1 - 1 - skip]]
        aa_count['w2'][aai] += 1

        aai = aa_map[p2_seq[ii + w]]
        aa_count['w2'][aai] += 1

        if p2_seq[i + w] in p_dict:
            features_p['w1'][i] = np.dot(aa_count['w1'][:20], p_dict[p2_seq[i + w]])
            features_p['w2'][i] = np.dot(aa_count['w2'][:20], p_dict[p2_seq[i + w]])
        else:
            features_p['w1'][i] = 0
            features_p['w2'][i] = 0
    return features_p


def gen_features_counts(seq, c_dict):
    features_counts = {'PDB': [], 'IDR': [], 'Linker': [], 'P_bind': [], 'N_bind': []}
    for ii, aa in enumerate(seq):
        for tg in ['PDB', 'IDR', 'Linker', 'P_bind', 'N_bind']:
            if aa in c_dict[tg]:
                features_counts[tg].append(c_dict[tg][aa])
            else:
                features_counts[tg].append(50.0)
    return features_counts


def code_one(seq, p_matrix, c_dict, verbose=False):
    # print(seq)
    features_counts = gen_features_counts(seq, c_dict)
    features_p = get_features_p(seq, p_dict=p_matrix, w1=10, w2=25, skip=0)
    features = features_counts
    for ftr in features_p:
        features[ftr] = features_p[ftr]
    if verbose:
        print(f">{len(seq)}\t{seq}")
        for ky in features_counts:
            print(f"{ky}({len(features_counts[ky])}):\t{features_counts[ky]}")
        for ky in features_p:
            print(f"{ky}({len(features_p[ky])}):\t{features_p[ky]}")
    return features


def assemble_features_one(ac, seq, features_dict, _w_out=10, _pad=250):
    w_in = _w_out + (2 * _pad)
    r_ftrs = []
    r_seq = []
    r_ac = []
    r_start = []
    sz = len(seq)
    features_list = list(features_dict.keys())
    wind_count = math.ceil(sz / _w_out)

    tmp_dict = {}
    zeros_left = [0] * _pad
    zeros_right = [0] * (_pad + (wind_count * _w_out) - sz)
    tmp_dict['seq'] = list(seq) + ['x'] * (_pad + (wind_count * _w_out) - sz)
    for ftr in features_list:
        tmp_dict[ftr] = list(zeros_left) + list(features_dict[ftr]) + list(zeros_right)

    for xst in range(wind_count):
        st = xst * _w_out
        r_ac.append(ac)
        r_start.append(st)
        r_seq.append(''.join(tmp_dict['seq'][st:st + _w_out]))
        tmp_lst = []
        for ftr in features_used:
            if ftr in features_list:
                tmp_lst.append(tmp_dict[ftr][st:st + w_in])
        r_ftrs.append(tmp_lst)
    return {'features': torch.tensor(np.array(r_ftrs, dtype='float32'), dtype=torch.float32),
            'seq': r_seq,
            'ac': r_ac,
            'start': r_start}


def smooth(s2):
    s1 = np.concatenate((np.array([s2[0]], dtype='float32'), s2))
    s0 = np.concatenate((np.array([s1[0]], dtype='float32'), s1))
    s3 = np.concatenate((s2, np.array([s2[-1]], dtype='float32')))
    s4 = np.concatenate((s3, np.array([s3[-1]], dtype='float32')))
    return (s0[:-2] + s1[:-1] + s2 + s3[1:] + s4[2:]) / 5.0


def score_ensemble(in_ac, in_seq, p_matrix, c_dict, models_list, priors_list, dev='cpu'):
    sigmoid = torch.nn.Sigmoid()
    features_dict = code_one(in_seq, p_matrix, c_dict, verbose=False)
    ftr_one = assemble_features_one(in_ac, in_seq, features_dict=features_dict, _w_out=w_out, _pad=pad)
    # print(ftr_one['features'].shape)
    output = []
    for ii, mm in enumerate(models_list):
        sc = sigmoid(mm(ftr_one['features'].to(dev))).reshape(-1).detach().numpy()
        output.append(bayes_evidence(sc, prior=priors_list[ii]))
    s2 = (output[0] + output[1] + output[2] + output[3]) / 4.0
    scores = smooth(s2)
    save_results(ac=in_ac, scores=scores, seq=in_seq)


def save_results(ac, scores, seq):
    for mdl in scores:
        out_file = f'{output_path}/{ac}.caid'
        with open(out_file, 'w') as fout:
            print(f'>{ac}', file=fout)
            for ii, aa in enumerate(seq):
                print(f"{ii+1}\t{aa}\t{scores[ii]:.5}", file=fout)
