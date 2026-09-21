import copy
# from crc64iso.crc64iso import crc64
from miscellaneous import get_url_response, get_uniparc_id
# import requests


# 3_updated
def annotated_fasta(data_name: str='Data has no Name', database_list: list=None, accession: str=None,
                    tags_dict: dict=None, tags_list: list=None):
    if database_list is None:
        database_list = []
    if tags_dict is None:
        tags_dict = {}
    if tags_list is None:
        tags_list = list(tags_dict.keys())
    if accession is not None:
        if accession not in database_list:
            database_list.append(accession)
    return {'data': {}, 'metadata': {'tags_dict': tags_dict, 'tags_list': tags_list, 'database_list': database_list,
                                     'counts': None, 'accession': accession, 'data_name': data_name}}


def aff_load0(in_file: str, data_name: str='Data has no Name', accession: str=None):  # , _mark=None
    af_sequences = {}
    tags_list = []
    tags_dict = {}
    database_list = []
    _more_tags = True
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                if _more_tags:
                    _lst = line.split()
                    if len(_lst) > 1:
                        if _lst[1] == 'TAG':
                            tags_list.append(_lst[2])
                            tags_dict[_lst[2]] = 'tag from an old af format'
                continue
            _more_tags = False

            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': '', 'tags': {}, 'databases': {}, 'scores': {}}
                # for extra in ac_lst[1:]:
                #     ex_lst = extra.split('=')
                #     if len(ex_lst) < 2:
                #         continue
                #     af_sequences[ac][ex_lst[0]] = ex_lst[1]
                #     if ex_lst[0] not in database_list:
                #         database_list.append(ex_lst[0])
                continue
            # af_sz = len(af_sequences[ac])
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
                ii = 0
                # tg_i0 = af_sz  # len(af_sequences[ac])
            else:
                af_sequences[ac]['tags'][tags_list[ii]] = line.replace('x', '-')
                ii += 1
            continue
    af = {'data': af_sequences, 'metadata': {'tags_dict': tags_dict, 'tags_list': tags_list,
                                             'database_list': database_list, 'statistics': None,
                                             'accession': '', 'data_name': 'No Name'}}
    # -------------------------------------------------------------------------------------------------------
    # for ac in af0['data']:
    #     for ntg in af0['metadata']['name_tags']:
    #
    # af = annotated_fasta(data_name=data_name)
    # for tg in af0['metadata']['tags']:
    #     af['metadata']['tags_dict'][tg] = 'xxx'
    # if accession is None:
    #     accession = 'Fasta'
    # af['metadata']['accession'] = accession
    # af['metadata']['database_list'].append(accession)
    # for ntg in af0['metadata']['name_tags']:
    #     af['metadata']['database_list'].append(ntg.split(';'))
    # for ac in af0['data']:
    #     af['data'][ac] = af0['data'][ac]
    #     # af['data'][ac][accession] = [ac]
    #     for ntg in af0['metadata']['database_list']:
    #         af['data'][ac]['databases'][ntg] = af0['data'][ac][ntg].split(';')
    #
    # _gen_database_counts(af)
    return af


def aff_load2(in_file: str, data_name: str='Data has no name'):  # , _mark=None
    af_sequences = {}
    tags_dict = {}
    tags_list = []
    names_list = []
    _more_tags = False
    _id_counts = False
    accession = None
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                if 'Data Name:' in line:
                    data_name = line.split('\t')[1]
                if 'Amino acid sequence' in line:
                    _more_tags = True
                    continue
                if _more_tags:
                    _lst = line.split('\t')
                    if len(_lst) > 1:
                        if len(_lst) > 2:
                            info = _lst[2]
                        else:
                            info = 'Annotation'
                        tags_dict[_lst[1]] = info
                        tags_list.append(_lst[1])
                    if 'ID Counts:' in line:
                        # print("ID counts", flush=True)
                        _more_tags = False
                        _id_counts = True
                    continue
                if _id_counts:
                    _lst = line.split()
                    if 'Unique#' in line:
                        continue
                    if 'Tags Counts:' in line:
                        _id_counts = False
                        continue
                    if len(_lst) < 4:
                        continue
                    if len(_lst) == 5:
                        accession = _lst[1]
                    names_list.append(_lst[1])
                continue
            if len(tags_dict) == 0:
                print(f"Error: No tags are found in {in_file}")
                return None

            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': '', 'scores': {}}
                for extra in ac_lst[1:]:
                    ex_lst = extra.split('=')
                    if len(ex_lst) != 2:
                        continue
                    lst2 = ex_lst[1].split(';')
                    af_sequences[ac][ex_lst[0]] = lst2
                    if ex_lst[0] not in names_list:
                        names_list.append(ex_lst[0])
                continue
            af_sz = len(af_sequences[ac])
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
                tg_i0 = af_sz
            else:
                af_sequences[ac][tags_list[af_sz - tg_i0]] = line.replace('x', '-')
            continue
    af = {'data': af_sequences, 'metadata': {'tags_dict': tags_dict, 'names_list': names_list, 'counts': None,
                                             'accession': accession, 'data_name': data_name}}
    for ac in af['data']:
        for ntg in af['metadata']['names_list']:
            if ntg not in af['data'][ac]:
                af['data'][ac][ntg] = []
    # _gen_name_counts(af)
    return af


# new func
def aff_clear_databases(af, clear_db_list=None, add_acc=False, clear_all_uniprot=True):
    _remove_up = False
    if 'UniProt' in af['metadata']['database_list']:
        if clear_db_list is None:
            _remove_up = True
        elif 'UniProt' in clear_db_list:
            _remove_up = True

    if clear_db_list is None:
        af['metadata']['database_list'] = []
    else:
        for db in clear_db_list:
            if db in af['metadata']['database_list']:
                af['metadata']['database_list'].remove(db)

    if add_acc:
        if 'Fasta' not in af['metadata']['database_list']:
            af['metadata']['database_list'].append('Fasta')
    if _remove_up:
        af['metadata']['database_list'].append('UniProt')

    for ac in af['data']:
        up0 = None
        if _remove_up:
            if len(af['data'][ac]['databases']['UniProt']) > 0:
                up0 = af['data'][ac]['databases']['UniProt'][0]
        af['data'][ac]['databases'] = {}
        if add_acc:
            af['data'][ac]['databases']['Fasta'] = [ac]
        if _remove_up:
            af['data'][ac]['databases']['UniProt'] = []
            if up0:
                af['data'][ac]['databases']['UniProt'].append(up0)


def _verify_file_keys(in_file: str):
    # 'ID Counts:'
    key_dict = {'Data Name:': 0,
                'Amino acid sequence': 0,
                'ID Counts:': 0,
                'Tag Counts:': 0
                }
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            # print(line, flush=True)
            if len(line) == 0:
                continue
            if line[0] == '#':
                for ky in key_dict:
                    if ky in line:
                        # print('===', flush=True)
                        key_dict[ky] += 1
            else:
                break
    for ky in key_dict:
        if key_dict[ky] < 1:
            print(f"Bad key: {ky} in {in_file}", flush=True)
            return False
    return True


# new func
def aff_load3(in_file: str, data_name: str='Data has no name'):  # , _mark=None
    if not _verify_file_keys(in_file):
        exit(0)
    af_sequences = {}
    tags_dict = {}
    tag_default_dict = {}
    tags_list = []
    database_list = []
    info_list = []
    _more_tags = False
    _id_counts = False
    accession = None
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                if 'Data Name:' in line:
                    # if data_name == 'Data has no name':
                    data_name = line.split('\t')[1]
                    continue
                if "# Info ....:" in line:
                    _lst = line.split('\t')
                    if len(_lst) > 1:
                        info_list.append(_lst[1].strip())
                    continue
                if 'Amino acid sequence' in line:
                    _more_tags = True
                    continue
                if _more_tags:
                    _lst = line.split('\t')
                    if len(_lst) > 1:
                        if len(_lst) > 2:
                            info = _lst[2]
                        else:
                            info = 'Annotation'
                        _tg = _lst[1]
                        if _tg[-1] == ':':
                            _tg = _tg[:-1]
                        tags_dict[_tg] = info
                        tags_list.append(_tg)
                        if 'Default:' in info:
                            tag_default_dict[_tg] = info.split()[-1]
                        else:
                            tag_default_dict[_tg] = '0'
                    if 'ID Counts:' in line:
                        # print("ID counts", flush=True)
                        _more_tags = False
                        _id_counts = True
                    continue
                if _id_counts:
                    _lst = line.split()
                    if 'Unique#' in line:
                        continue
                    if 'Tag Counts:' in line:
                        _id_counts = False
                        continue
                    if len(_lst) < 5:
                        continue
                    if len(_lst) == 6:
                        accession = _lst[1]
                    database_list.append(_lst[1])
                continue
            if len(tags_dict) == 0:
                print(f"Error: No tags are found in {in_file}")
                return None

            if line[0] == '>':
                ac_lst = line[1:].split('|')
                ac = ac_lst[0]
                af_sequences[ac] = {'seq': '', 'tags': {}, 'databases': {}, 'scores': {}}
                for extra in ac_lst[1:]:
                    ex_lst = extra.split('=')
                    if len(ex_lst) != 2:
                        continue
                    lst2 = ex_lst[1].split(';')
                    _db = ex_lst[0].split('(')[0]
                    af_sequences[ac]['databases'][_db] = lst2
                    if _db not in database_list:
                        database_list.append(_db)
                continue
            if af_sequences[ac]['seq'] == '':
                af_sequences[ac]['seq'] = line
            else:
                tag_idx = len(af_sequences[ac]['tags'])
                af_sequences[ac]['tags'][tags_list[tag_idx]] = line.replace('x', '-')
            continue
    af = {'data': af_sequences, 'metadata': {'tags_dict': tags_dict, 'tags_list': tags_list,
                                             'database_list': database_list, 'counts': None,
                                             'accession': accession, 'data_name': data_name,
                                             'tag_defaults': tag_default_dict, 'info_list': info_list}}
    for ac in af['data']:
        for ntg in af['metadata']['database_list']:
            if ntg not in af['data'][ac]['databases']:
                af['data'][ac]['databases'][ntg] = []
    return af


def aff_load_simple(in_file: str, data_name: str='Data has no name', tag: str= 'ANN',
                    tag_description: str= 'Source annotation'):
    af_sequences = {}
    tags_dict = {tag: tag_description}
    database_list = ['Fasta']
    _more_tags = False
    _id_counts = False
    accession = 'Fasta'
    l_num = 0
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                l_num = 1
                hd_lst = line[1:].split('|')
                ac = hd_lst[0]
                af_sequences[ac] = {'seq': '', 'tags': {}, 'databases': {'Fasta': [ac]}, 'scores': {}}  # accession: [ac]
                for hdn in hd_lst:
                    two_parts = hdn.replace(':', '=').split('=')
                    if len(two_parts) == 2:
                        af_sequences[ac]['databases'][two_parts[0]] = two_parts[1].split(';')
                        if two_parts[0] not in database_list:
                            database_list.append(two_parts[0])
                continue
            if l_num == 1:
                af_sequences[ac]['seq'] = line
                l_num = 2
                continue
            if l_num == 2:
                af_sequences[ac]['tags'][tag] = line
                l_num = 0
                continue
    af = {'data': af_sequences, 'metadata': {'tags_dict': tags_dict, 'tags_list': [tag], 'database_list': database_list,
                                             'counts': None, 'accession': accession, 'data_name': data_name}}
    for ac in af['data']:
        for ntg in af['metadata']['tags_list']:
            if ntg not in af['data'][ac]['tags']:
                af['data'][ac]['tags'][ntg] = ''
    _gen_database_counts(af)
    return af


def aff_load_fasta(in_file: str, data_name: str='Data has no name'):
    af = annotated_fasta(data_name=data_name)
    with open(in_file, 'r') as fin:
        ac = ''
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                ac = line[1:]
                af['data'][ac] = {'seq': '', 'tags': {}, 'databases': {}, 'scores': {}}
            else:
                af['data'][ac]['seq'] = af['data'][ac]['seq'] + line
    return af

def _validate_header_extra(he: str=None):
    _ret = True
    if he is None:
        return True
    lst = he.split('\n')
    for line in lst:
        if len(line) == 0:
            return False
        if line[0] != '#':
            return False
    return True


def aff_save2(af, f_name: str, header_top: str =None, header_bottom: str =None):
    with open(f_name, 'w') as fout:
        print(f"# Data Name:\t{af['metadata']['data_name']}\n#", file=fout)
        if header_top:
            if _validate_header_extra(header_top):
                print(header_top, file=fout)
            else:
                print(f"BAD header_top:\t{header_top}", flush=True)
        print(f"# Sequences:\t{len(af['data']):,}", file=fout)
        print("#", file=fout)
        print("# Format:", file=fout)
        print("#\t>accession", file=fout, end='')
        for n_tg in af['metadata']['names_list']:
            print(f"|{n_tg}={n_tg}_ID", file=fout, end='')
        print(file=fout)
        print("#\tAmino acid sequence", file=fout)
        for tg in af['metadata']['tags_dict']:
            print(f"#\t{tg}\t{af['metadata']['tags_dict'][tg]}", file=fout)
        print("#", file=fout)
        aff_gen_counts(af)
        str_counts = _get_string_counts(af)
        print(str_counts, file=fout)
        if header_bottom:
            if _validate_header_extra(header_bottom):
                print('#', file=fout)
                print(header_bottom, file=fout)
            else:
                print(f"BAD header_bottom:\t{header_bottom}", flush=True)
        print('#', file=fout)
        for ac in af['data']:
            ac_o = ac
            for tg in af['data'][ac]:
                if tg in af['metadata']['names_list']:
                    if len(af['data'][ac][tg]) > 0:
                        ac_o = f"{ac_o}|{tg}={af['data'][ac][tg][0]}"
                        if len(af['data'][ac][tg]) > 1:
                            for xx in af['data'][ac][tg][1:]:
                                ac_o = f"{ac_o};{xx}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)
            for tg in af['metadata']['tags_dict']:
                print(f"{af['data'][ac][tg]}", file=fout)
    return

# new func
def aff_save3(af, f_name: str, header_top: str =None, header_bottom: str =None):
    with open(f_name, 'w') as fout:
        print(f"# Data Name:\t{af['metadata']['data_name']}\n#", file=fout)
        for info in af['metadata']['info_list']:
            print(f"# Info ....:\t{info}", file=fout)
        print("#", file=fout)
        if header_top:
            if _validate_header_extra(header_top):
                print(header_top, file=fout)
            else:
                print(f"BAD header_top:\t{header_top}", flush=True)
        print(f"# Sequences:\t{len(af['data']):,}", file=fout)
        print("#", file=fout)
        print("# Format:", file=fout)
        print("#\t>accession", file=fout, end='')
        for _db in af['metadata']['database_list']:
            print(f"|{_db}={_db}_ID(s)", file=fout, end='')
        print("\n#\tAmino acid sequence", file=fout)
        for tg in af['metadata']['tags_list']:
            print(f"#\t{tg}:\t{af['metadata']['tags_dict'][tg]}", file=fout)
        print("#", file=fout)
        aff_gen_counts(af)
        str_counts = _get_string_counts(af)
        print(str_counts, file=fout)
        if header_bottom:
            if _validate_header_extra(header_bottom):
                print('#', file=fout)
                print(header_bottom, file=fout)
            else:
                print(f"BAD header_bottom:\t{header_bottom}", flush=True)
        print('#', file=fout)
        for ac in af['data']:
            # print(ac, list(af['data'][ac]['tags'].keys()))
            ac_o = ac
            for tg in af['metadata']['database_list']:
                if tg in af['data'][ac]['databases']:
                    db_sz = len(af['data'][ac]['databases'][tg])
                    if db_sz > 0:
                        ac_o = f"{ac_o}|{tg}({db_sz})={af['data'][ac]['databases'][tg][0]}"
                        if db_sz > 1:
                            for xx in af['data'][ac]['databases'][tg][1:]:
                                ac_o = f"{ac_o};{xx}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}", file=fout)
            for tg in af['metadata']['tags_list']:
                print(f"{af['data'][ac]['tags'][tg]}", file=fout)
    return


def aff_save_simple(af, f_name: str, tag: str):
    if tag not in af['metadata']['tags_dict']:
        print(f"Error in aff_save_simple {f_name}:\t{tag} not found")
        return
    with open(f_name, 'w') as fout:
        for ac in af['data']:
            ac_o = ac
            for ntg in af['metadata']['database_list']:
                if ntg in af['data'][ac]['databases']:
                    if len(af['data'][ac]['databases'][ntg]) > 0:
                        ac_o = f"{ac_o}|{ntg}={af['data'][ac]['databases'][ntg][0]}"
                        if len(af['data'][ac]['databases'][ntg]) > 1:
                            for xx in af['data'][ac]['databases'][ntg][1:]:
                                ac_o = f"{ac_o};{xx}"
            print(f">{ac_o}\n{af['data'][ac]['seq']}\n{af['data'][ac]['tags'][tag]}", file=fout)
    return


def aff_save_fasta(af, f_name: str, databases=False, mark: str =''):
    with open(f_name, 'w') as fout:
        for ac in af['data']:
            ac_o = ac
            if databases:
                for ntg in af['metadata']['database_list']:
                    if ntg in af['data'][ac]['databases']:
                        if len(af['data'][ac]['databases'][ntg]) > 0:
                            ac_o = f"{ac_o}|{ntg}={af['data'][ac]['databases'][ntg][0]}"
                            if len(af['data'][ac]['databases'][ntg]) > 1:
                                for xx in af['data'][ac]['databases'][ntg][1:]:
                                    ac_o = f"{ac_o};{xx}"
            print(f">{mark}{ac_o}\n{af['data'][ac]['seq']}", file=fout)


def aff_merge_simple(af_to, af_from):
    for ac in af_from['data']:
        if ac not in af_to['data']:
            af_to['data'][ac] = af_from['data'][ac]


def merge_annotations(ann1: str, ann2: str):
    if len(ann1) != len(ann2):
        return None
    sz = len(ann1)
    if sz <= 3:
        return None
    lst = ['-'] * sz
    for ii in range(sz):
        if ann1[ii] == '1' or ann2[ii] == '1':
            lst[ii] = '1'
        elif ann1[ii] == '0' or ann2[ii] == '0':
            lst[ii] = '0'
    return ''.join(lst)


# OK
def aff_remove_tags_list(af, tags_list_out: list):
    tags_list = list(af['metadata']['tags_list'])
    for tg in tags_list:
        if tg in tags_list_out:
            del af['metadata']['tags_dict'][tg]
            af['metadata']['tags_list'].remove(tg)
    for ac in af['data']:
        for tg in tags_list:
            if tg in tags_list_out:
                del af['data'][ac]['tags'][tg]


def aff_rename_tag(af, old_tag: str, new_tag: str, new_info: str):
    if old_tag in af['metadata']['tags_dict']:
        af['metadata']['tags_dict'][new_tag] = new_info
        ii = af['metadata']['tags_list'].index(old_tag)
        af['metadata']['tags_list'][ii] = new_tag
        del af['metadata']['tags_dict'][old_tag]

        for ac in af['data']:
            af['data'][ac]['tags'][new_tag] = af['data'][ac]['tags'][old_tag]
            del af['data'][ac]['tags'][old_tag]
    else:
        print(f"Tag {old_tag} not found", flush=True)


# OK
def aff_remove_missing_scores(af):
    used_scores_set = set()
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        for prd_sc in af['data'][ac]['scores']:
            used_scores_set.add(prd_sc)
    for ac in ac_list:
        for prd_sc in used_scores_set:
            if prd_sc not in af['data'][ac]['scores']:
                del af['data'][ac]
                break


# OK
def aff_remove_no_info_tag(af, tag: str):
    if tag not in af['metadata']['tags_dict']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if af['data'][ac]['tags'][tag].count('-') == len(af['data'][ac]['tags'][tag]):
            del af['data'][ac]


# OK
def aff_remove_no_class_tag(af, tag: str, cl: str):
    if tag not in af['metadata']['tags_dict']:
        return
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if cl not in af['data'][ac]['tags'][tag]:
            print('Removed:\t', ac)
            del af['data'][ac]


# OK
def aff_remove_no_class_any(af, cl: str):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        rmv = True
        for tg in af['metadata']['tags_dict']:
            if cl in af['data'][ac]['tags'][tg]:
                rmv = False
                break
        if rmv:
            del af['data'][ac]


# OK
def aff_remove_short(af, cut=15):
    ac_list = list(af['data'].keys())
    for ac in ac_list:
        if len(af['data'][ac]['seq']) < cut:
            del af['data'][ac]


# OK
def aff_add_tag(af, tag: str, info: str='New tag'):
    if tag not in af['metadata']['tags_dict']:
        af['metadata']['tags_list'].append(tag)
        af['metadata']['tags_dict'][tag] = info
        for ac in af['data']:
            af['data'][ac]['tags'][tag] = '-' * len(af['data'][ac]['seq'])
    else:
        print(f"{tag} exist in af", flush=True)


# new func
def aff_add_databases(af, requested_databases=None, max_id_count=10, d_file=None, uo_set=None,
                      source_id='Source_ID', verbose=False):
    id_dict = {}
    fout = None
    # print(uo_set)
    if d_file:
        print(f"Opening {d_file}", flush=True)
        fout = open(d_file, 'w')
    if 'UniProt' not in af['metadata']['database_list']:
        af['metadata']['database_list'].append('UniProt')
    if 'OX' not in af['metadata']['database_list']:
        af['metadata']['database_list'].append('OX')

    af2 = annotated_fasta()
    af2['metadata'] = af['metadata']
    if source_id not in af2['metadata']['database_list']:
        af2['metadata']['database_list'] = [source_id] + af2['metadata']['database_list']

    for ii, ac in enumerate(af['data']):
        if uo_set is not None:  # use only
            if ac not in uo_set:
                # print(ac)
                continue
        if verbose:
            print(f"{ii:,}\t{ac}", end='\t', flush=True)  # \t{list(af['data'][ac].keys())}
        else:
            print(f"{ii:,}\t{ac}", flush=True)
        up = aff_get_seq_databases(af=af, ac=ac, requested_databases=requested_databases,
                                   max_id_count=max_id_count, fout=fout, verbose=verbose)
        if up:
            id_dict[ac] = up
            if up not in af2['data']:
                af2['data'][up] = af['data'][ac]
            else:
                merge_entries(af2['data'][up], af['data'][ac])
            if source_id not in af2['data'][up]['databases']:
                af2['data'][up]['databases'][source_id] = [ac]
            elif ac not in af2['data'][up]['databases'][source_id]:
                af2['data'][up]['databases'][source_id].append(ac)
    return af2, id_dict


def merge_entries(into, entry):
    if len(into['seq']) != len(entry['seq']):
        print(f"Error len(seq) unequal:\t{len(into['seq'])}\t{len(entry['seq'])}")
        exit(1)
    if into['seq'] != entry['seq']:
        print(f"Error seq")
    for tg in entry['tags']:
        if tg not in into['tags']:
            into['tags'][tg] = entry['tags'][tg]
        else:
            into['tags'][tg] = merge_annotations(ann1=into['tags'][tg], ann2=entry['tags'][tg])


def _get_string_counts(af):
    _msg = ''
    if af['metadata']['counts'] is None:
        return _msg
    # if len(af['metadata']['database_list']) > 0:
    _msg = _msg + f"# ID Counts:\n#\tID \tSeq#\tTotal#\tUnique#"
    for ntg in af['metadata']['database_list']:
        _msg = _msg + f"\n#\t{ntg}\t{af['metadata']['counts']['database_dict'][ntg]['sequences']:,}"
        _msg = _msg + f"\t{af['metadata']['counts']['database_dict'][ntg]['total']:,}"
        _msg = _msg + f"\t{af['metadata']['counts']['database_dict'][ntg]['unique']:,}"
        # if ntg == af['metadata']['accession']:
        #     _msg = _msg + "\tAC"
    _msg = _msg + "\n#\n"
    _msg = _msg + "# Tag Counts:\n#\ttag\tSeq#\tSeg#\t'0'\t'1'\t'-'"
    for tg in af['metadata']['tags_list']:
        _msg = _msg + f"\n#\t{tg}"
        for cc in ['seq', 'seg', '0', '1', '-']:  # af['metadata']['counts']['tags_dict'][tg]:
            _msg = _msg + f"\t{af['metadata']['counts']['tags_dict'][tg][cc]:,}"
    return _msg


def aff_gen_counts(af):
    _gen_tag_counts(af)
    _gen_database_counts(af)


def _gen_tag_counts(af):
    if af['metadata']['counts'] is None:
        af['metadata']['counts'] = {}
    af['metadata']['counts']['tags_dict'] = {}
    for tg in af['metadata']['tags_dict']:
        af['metadata']['counts']['tags_dict'][tg] = {'seq': 0, 'seg': 0, '0': 0, '1': 0, '-': 0}
        for ac in af['data']:
            mask = str(af['data'][ac]['tags'][tg])
            c_set = set(list(mask))
            if '1' in c_set:
                c_set.remove('1')
                for oc in c_set:
                    if oc != '0':
                        mask = mask.replace(oc, '0')
                # =========================================================
                cnt = len([xx for xx in mask.split('0') if xx])
                if cnt > 0:
                    af['metadata']['counts']['tags_dict'][tg]['seq'] += 1
                    af['metadata']['counts']['tags_dict'][tg]['seg'] += cnt
                    for cc in ['0', '1', '-']:
                        af['metadata']['counts']['tags_dict'][tg][cc] += af['data'][ac]['tags'][tg].count(cc)


def _gen_database_counts(af):
    if af['metadata']['counts'] is None:
        af['metadata']['counts'] = {}
    af['metadata']['counts']['database_dict'] = {}
    ntg_set_dict = {}
    for ntg in af['metadata']['database_list']:
        cnt = 0
        af['metadata']['counts']['database_dict'][ntg] = {'sequences': 0,'total': 0, 'unique': 0}
        ntg_set_dict[ntg] = set()
        for ac in af['data']:
            if ntg in af['data'][ac]['databases']:
                if len(af['data'][ac]['databases'][ntg]) > 0:
                    af['metadata']['counts']['database_dict'][ntg]['sequences'] += 1
                    af['metadata']['counts']['database_dict'][ntg]['total'] += len(af['data'][ac]['databases'][ntg])
                    ntg_set_dict[ntg] = ntg_set_dict[ntg].union(set(af['data'][ac]['databases'][ntg]))
                    cnt += len(af['data'][ac]['databases'][ntg])
        # print(ntg, cnt)
    for ntg in af['metadata']['database_list']:
        af['metadata']['counts']['database_dict'][ntg]['unique'] = len(ntg_set_dict[ntg])


def _process_uniprot_obsolete(in_up_list, seq, max_up_count=10, verbose=False):
    ac_list = []
    for o_ac in in_up_list:
        _vv = o_ac.split('.')[1]
        _ac = o_ac.split('.')[0]
        url = f"https://rest.uniprot.org/unisave/{_ac}?format=fasta&versions={_vv}"
        response = get_url_response(url)
        if response is not None:
            r_seq = ''
            for ss in response.text.split('\n')[1:]:
                if len(ss) > 0:
                    r_seq = r_seq + ss
            if r_seq == seq:
                ac_list.append(o_ac)
            else:
                if verbose:
                    print('#', end='')
            if len(ac_list) >= max_up_count:
                break
    return ac_list


# new func
def _process_uniprot_list(in_up_list, seq, db_dict, md_db_list, requested_databases=None, max_up_count=10,
                          fout=None, verbose=False):
    up2 = None
    for ac in in_up_list:
        if '.' in ac:
            continue
        url = f'https://rest.uniprot.org/uniprotkb/{ac}.json'
        response = get_url_response(url)
        if response is None:
            print("_process_uniprot_list, response is None")
            continue
        data = response.json()

        # verify seq
        ss = data.get('sequence')
        if ss:
            d_seq = ss['value']
            if seq != d_seq:
                if verbose:
                    print('#', end='', flush=True)
                continue
        else:
            continue

        # ADD ac
        if ac not in db_dict['UniProt']:
            db_dict['UniProt'].append(ac)
        if up2 is None:
            up2 = get_uniparc_id(ac)
        tsv_out = {}

        # ADD ox
        ox = ''
        _taxa = data.get('organism')
        if _taxa:
            if 'taxonId' in _taxa:
                ox = str(_taxa['taxonId'])
                if ox not in db_dict['OX']:
                    db_dict['OX'].append(ox)

        # add databases
        for xref in data.get("uniProtKBCrossReferences", []):
            _db = xref["database"]
            if requested_databases is None or _db in requested_databases:
                if fout:
                    if _db not in tsv_out:
                        tsv_out[_db] = []
                    tsv_out[_db].append(xref["id"])

                if _db not in db_dict:
                    db_dict[_db] = []

                if xref["id"] not in db_dict[_db]:
                    db_dict[_db].append(xref["id"])
                    if _db not in md_db_list:
                        md_db_list.append(_db)
        if len(db_dict['UniProt']) >= max_up_count:
            break
        if fout:
            print(f">{ac}|OX={ox}", end='', file=fout)
            for _db in tsv_out:
                print(f"|{_db}={tsv_out[_db][0]}", end='', file=fout)
                for xx in tsv_out[_db][1:]:
                    print(f";{xx}", end='', file=fout)
            print(flush=True, file=fout)
    return up2


# new func
# def aff_get_databases(af, requested_other_ids=None, max_id_count=10, verbose=False):
#     used_set = set()
#     for ii, ac0 in enumerate(af['data']):
#         aff_get_seq_databases(af, ac=ac0, max_id_count=max_id_count, verbose=False)
#         if 'UniProt' in af['data'][ac0]['databases']:
#             for jj, ac in enumerate(af['data'][ac0]['UniProt']):
#                 if verbose:
#                     print(f"{ii}/{len(af['data']):,}\t{ac0}\t{jj}/{len(af['data'][ac0]['UniProt']):,}\t{ac}\tother_ids:",
#                           end='\t', flush=True)
#                 else:
#                     print(f"{ii}/{len(af['data']):,}\t{ac0}\t{jj}/{len(af['data'][ac0]['UniProt']):,}\t{ac}")
#
#                 ret_dict = _get_all_ids(ac)
#                 for _db in ret_dict:
#                     if requested_other_ids is not None:
#                         if _db not in requested_other_ids and _db != 'OX':
#                             continue
#                     if verbose:
#                         print(f"{_db}", end='\t', flush=True)
#                     if len(ret_dict[_db]) < 1:
#                         continue
#                     if _db in af['data'][ac0]:
#                         af['data'][ac0][_db] = list(set(ret_dict[_db] + af['data'][ac0][_db]))
#                     else:
#                         af['data'][ac0][_db] = ret_dict[_db]
#                     used_set.add(_db)
#                 if verbose:
#                     print(flush=True)
#     for _db in used_set:
#         if _db not in af['metadata']['names_list']:
#             af['metadata']['names_list'].append(_db)
#         for ac in af['data']:
#             if _db not in af['data'][ac]:
#                 af['data'][ac][_db] = []
#


def aff_get_seq_databases(af, ac, requested_databases=None, max_id_count=10, fout=None, verbose=False):
    seq = af['data'][ac]['seq']
    ac_list2 = [[], []]

    checksum = crc64(seq)
    response = get_url_response(f"https://rest.uniprot.org/uniparc/search?query=checksum: {checksum}")
    if response is not None:
        data = response.json()
        if verbose:
            print(f"{len(data['results'][0]['uniProtKBAccessions']):,}", end='', flush=True)
        for _ac in data["results"][0]['uniProtKBAccessions']:
            if '.' in _ac:
                ac_list2[1].append(_ac)
            else:
                ac_list2[0].append(_ac)
        if verbose:
            print(f"\t{len(ac_list2[0])}\t{len(ac_list2[1])}", end='', flush=True)
    else:
        # print somthing
        return None

    if 'UniProt' not in af['data'][ac]['databases']:
        af['data'][ac]['databases']['UniProt'] = []
    if 'OX' not in af['data'][ac]['databases']:
        af['data'][ac]['databases']['OX'] = []

    up2 = _process_uniprot_list(in_up_list=ac_list2[0], db_dict=af['data'][ac]['databases'],
                                md_db_list=af['metadata']['database_list'], max_up_count=max_id_count,
                                requested_databases=requested_databases, seq=seq, fout=fout, verbose=verbose)
    if len(af['data'][ac]['databases']['UniProt']) == 0:
        af['data'][ac]['databases']['UniProt'] = _process_uniprot_obsolete(in_up_list=ac_list2[1], seq=seq,
                                                                           max_up_count=max_id_count, verbose=verbose)
    if len(af['data'][ac]['databases']['UniProt']) > 1:
        up_list = af['data'][ac]['databases']['UniProt']
        for up in up_list:
            if '-' in up:
                up0 = up.split('-')[0]
                if up0 in af['data'][ac]['databases']['UniProt']:
                    af['data'][ac]['databases']['UniProt'].remove(up)

    if verbose:
        print(f"\t{len(af['data'][ac]['databases']['UniProt'])}", flush=True)
    return up2


def aff_tag_size(af, tag, sz_range=None):
    if sz_range is None:
        return
    if tag not in af['metadata']['tags_list']:
        return
    print(f"Size:\t{sz_range}")
    _from = sz_range[0]
    _to = sz_range[1]
    for ac in af['data']:
        lst = list(af['data'][ac]['tags'][tag])
        ss = af['data'][ac]['tags'][tag].replace('0', '-')
        lst2 = [xx for xx in ss.split('-') if xx != '']
        st = 0
        for xx1 in lst2:
            st = ss.find(xx1, st)
            if sz_range[0] <= len(xx1) < sz_range[1]:
                # print(ac, st+1, ss)
                pass
            else:
                for ii in range(st, st+len(xx1)):
                    lst[ii] = '-'
            st = st + len(xx1)
        af['data'][ac]['tags'][tag] = ''.join(lst)


def aff_mark_site(af, ac, tag, start, end, m_with='1', verbose=True):
    if ac not in af['data']:
        if verbose:
            print(f"AC not found", flush=True)
        return
    if tag not in af['data'][ac]['tags']:
        if verbose:
            print(f"Tag {tag} not found", flush=True)
        return
    if end > len(af['data'][ac]['seq']):
        if verbose:
            print(f"End {end:,} outside seq {len(af['data'][ac]['seq']):,}", flush=True)
        return
    lst = list(af['data'][ac]['tags'][tag])
    for ii in range(start-1, end):
        lst[ii] = m_with
    af['data'][ac]['tags'][tag] = ''.join(lst)


def aff_build_seq_map(af):
    af['seq_map'] = {}
    for ky in af['data']:
        seq = af['data'][ky]['seq']
        if seq in af['seq_map']:
            af['seq_map'][seq].append(ky)
        else:
            af['seq_map'][seq] = [ky]

def aff_seq_to_key_list(af, seq, candidate_list, verbose=False):
    key_list = []
    if 'seq_map' not in af:
        aff_build_seq_map(af)
    if seq in af['seq_map']:
        return af['seq_map'][seq]
    else:
        af['seq_map'][seq] = candidate_list
        for ky in candidate_list:
            if ky in af['data']:
                if verbose:
                    print(f"aff_seq_to_key_list: {ky} in af['data']", flush=True)
                continue
            key_list.append(ky)
            af['data'][ky] = {'seq': seq, 'tags': {}, 'scores': {}, 'databases': {}}
            for tg in af['metadata']['tags_list']:
                cc = '0'
                if tg in ['CAID-IDR_PDB', 'IDEAL-IDR']:
                    cc = '-'
                af['data'][ky]['tags'][tg] = cc * len(seq)
            for db in af['metadata']['database_list']:
                af['data'][ky]['databases'][db] = []
    return key_list


def add_databases_and_tags(_af, database_list, tag_dict, verbose=False):
    for db in database_list:
        if db not in _af['metadata']['database_list']:
            _af['metadata']['database_list'].append(db)
            for ac in _af['data']:
                _af['data'][ac]['databases'][db] = []
    for tg in tag_dict:
        aff_add_tag(_af, tag=tg, info=tag_dict[tg])

