from annotated_fasta import *
import json


def aff_mdb_json_to_af(mdb_file_list: list = None, tag_list: list = None, quality_list: list = None, verbose=False):
    af = annotated_fasta(data_name='MobiDB from JSON')
    tag = 'lip'  # 'binding_mode_disorder_to_disorder'  # 'disorder'  # , 'disorder'
    for fl in mdb_file_list:
        with open(fl, 'r') as fin:
            data = json.load(fin)
            for xx in data:
                if 'sequence' not in xx:
                    # print(xx['acc'], list(xx.keys()))
                    # exit(0)
                    continue
                if f'curated-{tag}-merge' not in xx:
                    lst = [z for z in list(xx.keys()) if f'curated-{tag}' in z]
                    if len(lst) > 0:
                        print(xx['acc'], lst)

    return


def aff_mdb_fasta_to_af(mdb_file_list: list = None, tag_list: list = None, quality_list: list = None):
    if mdb_file_list is None:
        print(f"Bad MobiDB fasta files:\t{mdb_file_list}", flush=True)
        return None
    if tag_list is None:
        tag_list = ['disorder', 'lip', 'binding_mode_disorder_to_disorder']
    if quality_list is None:
        quality_list = ['curated', 'derived', 'homology']
    af = annotated_fasta(data_name="MobiDB-LIP")
    for fl in mdb_file_list:
        with open(fl, 'r') as fin:
            seq_next = False
            tag = ''
            ac_seq = ''
            ac = ''
            for line in fin:
                line = line.strip()
                if len(line) == 0:
                    continue
                if line[0] == '>':  # ==================================================
                    lst = line[1:].split('|')
                    if len(lst) >= 2:
                        # print(f"{line}\t*{lst[1].strip().split()[0]}*", flush=True)
                        if lst[1].strip().split()[0] == 'sequence':
                            if lst[0] not in af['data']:
                                ac_seq = lst[0]
                                af['data'][ac_seq] = {'seq': '', 'tags': {}, 'databases': {}, 'scores': {}}
                                seq_next = True
                            else:
                                ac_seq = ''
                                seq_next = False
                            continue
                        else:
                            seq_next = False
                            if lst[0] != ac_seq:
                                continue
                            t_lst = lst[1].split('-')
                            if t_lst[0] in quality_list and t_lst[1] in tag_list:
                                tag = lst[1]
                                continue
                            else:
                                tag = ''
                else:
                    if seq_next:
                        af['data'][ac_seq]['seq'] = line
                        continue
                    if len(tag) > 0 and len(ac_seq) > 0:
                        af['data'][ac_seq]['tags'][tag] = line
                        if tag not in af['metadata']['tags_dict']:
                            af['metadata']['tags_dict'][tag] = 'MobiDB tag'
                            af['metadata']['tags_list'].append(tag)
    for ac in af['data']:
        for tg in af['metadata']['tags_list']:
            if tg not in af['data'][ac]['tags']:
                af['data'][ac]['tags'][tg] = '-' * len(af['data'][ac]['seq'])

    af['metadata']['tags_list'].sort()
    # print(af['metadata']['tags_list'], flush=True)
    return af


# tags_quality_dict = {'HQ': ['curated'], 'LQ': ['derived', 'homology']}
def aff_mobidb_refine(af, tags_quality_dict):
    raf = annotated_fasta()
    raf['metadata']['database_list'] = af['metadata']['database_list']
    for ac in af['data']:
        # print(ac, flush=True)
        raf['data'][ac] = {'seq': af['data'][ac]['seq'], 'tags': {}, 'tmp': {},
                           'databases': af['data'][ac]['databases'],
                           'scores': {}}
        ## High Quality first
        for tag in af['data'][ac]['tags']:
            lst = tag.split('-')
            if len(lst) != 3:
                print(f"Error:\t{ac}\t{tag}", flush=True)
                continue
            t_qual = lst[0]
            tg = lst[1]
            tg_prd = lst[2]
            if tg_prd in ['merge', 'priority']:
                continue
            if t_qual in tags_quality_dict['HQ']:
                if tg not in raf['metadata']['tags_dict']:
                    raf['metadata']['tags_dict'][tg] = 'MobiDB Tag'
                    raf['metadata']['tags_list'].append(tg)
                if tg not in raf['data'][ac]['tmp']:
                    raf['data'][ac]['tmp'][tg] = ['-'] * len(af['data'][ac]['seq'])
                if len(af['data'][ac]['seq']) != len(af['data'][ac]['tags'][tag]):
                    print(ac, flush=True)
                    continue
                for ii in range(len(af['data'][ac]['seq'])):
                    if af['data'][ac]['tags'][tag][ii] == '1':
                        raf['data'][ac]['tmp'][tg][ii] = '1'
                    elif af['data'][ac]['tags'][tag][ii] == '0' and raf['data'][ac]['tmp'][tg][ii] != '1':
                        raf['data'][ac]['tmp'][tg][ii] = '0'

        ## Low quality second
        for tag in af['data'][ac]['tags']:
            lst = tag.split('-')
            if len(lst) != 3:
                print(f"Error:\t{ac}\t{tag}", flush=True)
                continue
            t_qual = lst[0]
            tg = lst[1]
            tg_prd = lst[2]
            if tg_prd in ['merge', 'priority']:
                continue
            if t_qual in tags_quality_dict['LQ']:
                if tg not in raf['metadata']['tags_dict']:
                    continue
                if tg not in raf['data'][ac]['tmp']:
                    continue
                if len(af['data'][ac]['seq']) != len(af['data'][ac]['tags'][tag]):
                    print(ac, flush=True)
                    continue
                for ii in range(len(af['data'][ac]['seq'])):
                    if af['data'][ac]['tags'][tag][ii] == '1' and raf['data'][ac]['tmp'][tg][ii] != '1':
                        raf['data'][ac]['tmp'][tg][ii] = '-'
        for tg in raf['data'][ac]['tmp']:
            raf['data'][ac]['tags'][tg] = ''.join(raf['data'][ac]['tmp'][tg])
    return raf
