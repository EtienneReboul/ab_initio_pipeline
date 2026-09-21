from annotated_fasta import *


def aff_load_cdhit_clusters(af, cdhit_clstr_file):
    cluster_dict = {}
    with open(cdhit_clstr_file, 'r') as fin:
        cluster = ''
        p_identity = ''
        for line in fin:
            line = line.strip()
            if len(line) < 2:
                continue
            lst = line.split()
            if lst[0][0] == '>':
                cluster = lst[1]
                continue
            ac = lst[2].split('.')[0][1:]
            if lst[3] == '*':
                cluster_dict[cluster] = ac
                p_identity = '100.00%'
            elif lst[3] == 'at':
                p_identity = lst[4]
            else:
                print(f"ERROR ===============================================\t{line}", flush=True)
            if ac not in af['data']:
                print(ac, flush=True)
                continue
            af['data'][ac]['databases']['cluster'] = [cluster]
            af['data'][ac]['databases']['p_identity'] = [p_identity]
    ac_list = list(af['data'].keys())
    # print("----", len(ac_list), flush=True)
    # exit(0)
    for ac in ac_list:
        if 'cluster' not in af['data'][ac]['databases']:
            del af['data'][ac]
        # af['data'][ac]['c_center'] = cluster_dict[str(af['data'][ac]['cluster'])]
    # print("----", len(ac_list), flush=True)
    af['metadata']['database_list'].append('cluster')
    af['metadata']['database_list'].append('p_identity')
    # af['metadata']['database_list'].append('c_center')


def assemble():
    _p = '/home/nmalhis/Papers_data/25_Padua_1/sequences/af/'
    _f = '/home/nmalhis/Tools/CD-HIT/CLIP/all.fasta'
    af = aff_load_fasta(f"{_f}")
    print(len(af['data']))
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if 'CLIP_' not in ac:
            del af['data'][ac]
    print(len(af['data']))
    af2 = aff_load_fasta(f"{_p}/merged2_DBs_cleaned.fasta")
    print(len(af2['data']))
    for ac in af2['data']:
        af['data'][ac] = af2['data'][ac]
    print(len(af['data']))
    aff_save_fasta(af=af, f_name='/home/nmalhis/Tools/CD-HIT/CLIP2/all.fasta')


def filter():
    _p = '/home/nmalhis/Tools/CD-HIT/CLIP2/'
    af = aff_load_fasta(f"{_p}all40.fasta")
    cdhit_clstr_file = f"{_p}all_30.fasta.clstr"
    aff_load_cdhit_clusters(af, cdhit_clstr_file)
    ac_list = list(af['data'].keys())
    clusters_out = set()
    for ac in af['data']:
        if 'CLIP' in ac:
            clusters_out.add(af['data'][ac]['databases']['cluster'][0])
    print(len(clusters_out))
    for ac in ac_list:
        # if 'CLIP' not in ac:
        if af['data'][ac]['databases']['cluster'][0] in clusters_out:
            del af['data'][ac]
    print(len(af['data']))
    af['metadata']['database_list'] = []
    aff_save_fasta(af=af, f_name=f"{_p}all30.fasta")
