import numpy as np
from annotated_fasta import *
import os


def aff_load_prd_merged_caid_scores(af, caid_scores_file, prd):
    with open(caid_scores_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue
            if line[0] == '>':
                ac = line[1:]
                if ac in af['data']:
                    if 'scores' not in af['data'][ac]:
                        af['data'][ac]['scores'] = {}
                    af['data'][ac]['scores'][prd] = []
                else:
                    ac = ''
                continue
            if len(ac) < 2:
                continue
            sc = float(line.split()[2])
            af['data'][ac]['scores'][prd].append(sc)
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        # foe each predictor
        if prd in af['data'][ac]['scores']:
            if len(af['data'][ac]['scores'][prd]) != len(af['data'][ac]['seq']):
                # scores do not match this seq
                print(f"{prd}\t{ac}\t{len(af['data'][ac]['scores'][prd])}\t{len(af['data'][ac]['seq'])}")
                del af['data'][ac]['scores'][prd]


def aff_load_protein_caid_scores(in_file):
    sc_list = []
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                ac = line[1:]
            else:
                sc_list.append(float(line.split()[2]))
    return np.array(sc_list, dtype='float32')


def aff_save_protein_caid_scores(af, ac, prd, out_file):
    with open(out_file, 'w') as fout:
        print(f">{ac}", file=fout)
        for ii, sc in enumerate(af['data'][ac]['scores'][prd]):
            print(f"{ii+1}\t{af['data'][ac]['seq'][ii]}\t{sc}", file=fout)
    return


def aff_load_caid_scores(af, scores_path, prd_list, merged=True, remove_missing_scores=False):

    path_files = os.listdir(scores_path)
    if merged:
        for prd in prd_list:
            if f"{prd}.caid" not in path_files:
                print(f"{prd}.caid not found")
                continue
            aff_load_prd_merged_caid_scores(af, caid_scores_file=f'{scores_path}{prd}.caid', prd=prd)
    else:
        for prd in prd_list:
            _p = f"{scores_path}{prd}/"
            ac_list = [x[:-5] for x in os.listdir(_p) if x[-5:] == '.caid']
            for ac in ac_list:
                if ac not in af['data']:
                    continue
                if 'scores' not in af['data'][ac]:
                    af['data'][ac]['scores'] = {}
                af['data'][ac]['scores'][prd] = aff_load_protein_caid_scores(in_file=f'{_p}{ac}.caid')
                if len(af['data'][ac]['scores'][prd]) != len(af['data'][ac]['seq']):
                    print(prd, ac, len(af['data'][ac]['scores'][prd]), len(af['data'][ac]['seq']))

    used_prd_set = set()
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        for prd in af['data'][ac]['scores']:
            used_prd_set.add(prd)

    if remove_missing_scores:
        aff_remove_missing_scores(af)
    return used_prd_set


