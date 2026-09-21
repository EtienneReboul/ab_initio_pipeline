import requests
import xml.etree.ElementTree as et
import json
import os


dbs_database_list = ['srcUniProt']

def get_dbs_tags_dict(source):
    return {'IDR': f'{source} disorder', 'DtoO': f'{source} disorder to order transition',
            'binding_protein': f'{source} protein bind', 'binding_nucleic': f'{source} nucleic bind',
            'IDR_partner': f'{source} partner disorder', 'binding_partner': f'{source} partner protein bind'
            }

def get_dbs_ac_tags(sz):
    return {'IDR': '-' * sz, 'DtoO': '-' * sz, 'binding_protein': '-' * sz, 'binding_nucleic': '-' * sz,
            'IDR_partner': '-' * sz, 'binding_partner': '-' * sz,
            'list': {'IDR': ['-'] * sz, 'DtoO': ['0'] * sz, 'binding_protein': ['0'] * sz,
                     'binding_nucleic': ['0'] * sz, 'IDR_partner': ['-'] * sz, 'binding_partner': ['0'] * sz
                     }
            }


def get_url_response(url, **kwargs):
    response = None
    for attempt in range(1, 11):
        try:
            response = requests.get(url, **kwargs)
            if response is not None:
                break
            else: print(f"Attempt {attempt} failed, response is None:\t{url}")
        except Exception as e:
            print(f"Attempt {attempt} failed: {url}")

    if not response.ok:
        print(response.text)
        return None
    return response


def is_float(text):
    try:
        float(text)
        return True
    except ValueError:
        return False


def get_go_term_lineage(go_id, verbose=False):
    url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}/ancestors?relations=is_a,part_of"
    response = get_url_response(url)  # requests.get(url, headers={"Accept": "application/json"})
    if response is None:
        return None
    # get_url_response
    if response.status_code == 200:
        data = response.json()
        if 'ancestors' not in data['results'][0]:
            if verbose:
                print("KeyError: 'ancestors'")
            return None
        ancestors = data['results'][0]['ancestors']
        if verbose:
            print(f"GO Term: {go_id}")
            print("Ancestors:", ancestors)
        return ancestors
    else:
        print(f"Failed to retrieve data for {go_id}")
    return None


def get_xml_root(xml_file):
    root = None
    try:
        tree = et.parse(xml_file)
        root = tree.getroot()
    except FileNotFoundError:
        print(f"{xml_file}\tFile not found.")
        exit(1)
    except et.ParseError:
        print(f"{xml_file}\tInvalid XML format.")
        exit(1)
    return root

def load_json(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
        return data
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        return None
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return None


def get_uniprot_seq(ac):
    seq = None
    ox = None
    url = f'https://rest.uniprot.org/uniprotkb/{ac}.json'
    response = get_url_response(url)
    if response is None:
        print("_process_uniprot_list, response is None")
        return seq, ox
    data = response.json()

    ss = data.get('sequence')
    if not ss:
        return ss, ox
    seq = ss['value']

    _taxa = data.get('organism')
    if not _taxa:
        return ss, ox
    if 'taxonId' not in _taxa:
        return ss, ox

    ox = str(_taxa['taxonId'])

    return seq, ox


def run_iup(fasta_path, iup_path, out_path):
    f_list = os.listdir(fasta_path)
    for ff in f_list:
        cmd = f"python3 {iup_path}iupred3.py -a {fasta_path}{ff} short > {out_path}{ff.split('.')[0]}.iup3"
        os.system(cmd)


def split_fasta(in_file, out_path):
    from annotated_fasta import aff_load_fasta
    af = aff_load_fasta(in_file)
    for ac in af['data']:
        o_file = f"{out_path}{ac}.fasta"
        with open(o_file, 'w') as fout:
            print(f">{ac}\n{af['data'][ac]['seq']}", file=fout)

def split_morf_chibi(in_file, out_path):
    prd_list = ['MCW', 'MCL', 'MC']
    # prd_list = ['MCL', 'MC']
    # fout = {'MCW': None, 'MCL': None, 'MC': None}  # , 'IDP': None}
    ac = ''
    fout = {}
    for prd in prd_list:
        fout[prd] = None
        if not os.path.isdir(f"{out_path}/{prd}"):
            os.system(f"mkdir {out_path}/{prd}")
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '#':
                continue
            if line[0] == '>':
                ac = line[1:]
                for prd in fout:
                    if fout[prd] is not None:
                        if not fout[prd].closed:
                            fout[prd].close()
                    fout[prd] = open(f"{out_path}/{prd}/{ac}.caid", 'w')
                    print(f">{ac}", file = fout[prd])
                continue
            lst = line.split()
            for ii, prd in enumerate(prd_list):
                print(f"{lst[0]}\t{lst[1]}\t{lst[2 + ii]}", file = fout[prd])
        for prd in prd_list:
            if not fout[prd].closed:
                fout[prd].close()


def split_f_morf_pred(in_file, out_path):
    with open(in_file, 'r') as fin:
        ac = ''
        seq = ''
        cnt = 1
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '>':
                ac = line[1:]
                seq = ''
                cnt = 1
                continue
            if cnt == 1:
                seq = line.upper()
                cnt += 1
                continue
            if cnt == 2:
                sc_list = line.split(',')
                if len(ac) < 2 or len(seq) < 5:
                    continue
                with open(f"{out_path}{ac}.caid", 'w') as fout:
                    print(f">{ac}", file=fout)
                    for ii in range(len(seq)):
                        print(f"{ii+1}\t{seq[ii]}\t{sc_list[ii]}", file=fout)
                ac = ''
                seq = ''
            cnt += 1


def split_diso_rdp_bind(in_file, out_path):
    with open(in_file, 'r') as fin:
        ac = ''
        seq = ''
        cnt = 1
        for line in fin:
            line = line.strip()
            if len(line) == 0:
                continue
            if line[0] == '>':
                ac = line[1:]
                seq = ''
                cnt = 1
                continue
            if cnt == 1:
                seq = line.upper()
                cnt += 1
                continue
            if cnt == 7:
                sc_list = line.split(':')[1].split(',')
                if len(ac) < 2 or len(seq) < 5:
                    continue
                with open(f"{out_path}{ac}.caid", 'w') as fout:
                    print(f">{ac}", file=fout)
                    for ii in range(len(seq)):
                        print(f"{ii+1}\t{seq[ii]}\t{sc_list[ii]}", file=fout)
                ac = ''
                seq = ''
            cnt += 1


def split_portal(af, in_file, prd_dict, out_path):
    print(in_file, len(af['data']))
    dta = []
    seq_lst = []
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) < 3:
                continue
            if 'Position' in line:
                continue
            lst = line.split()
            dd = {}
            seq_lst.append(lst[1])
            dta.append({})
            for prd in prd_dict:
                dta[-1][prd] = float(lst[prd_dict[prd]])
    seq = ''.join(seq_lst)
    ac = ''
    for _ac in af['data']:
        if af['data'][_ac]['seq'] == seq:
            ac = _ac
            break
    for prd in prd_dict:
        out_file = f"{out_path}{prd}/{ac}.caid"
        with open(out_file, 'w') as fout:
            print(f">{ac}", file=fout)
            for ii in range(len(seq)):
                print(f"{ii+1}\t{seq[ii]}\t{dta[ii][prd]}", file=fout)


def get_uniparc_id(ac):
    url = f'https://rest.uniprot.org/uniprotkb/{ac}.json'
    response = get_url_response(url)
    if response is None:
        print(f"Error, response is None for {ac}", flush=True)
        return None
    data = response.json()
    if "extraAttributes" not in data:
        return None
    return data["extraAttributes"]['uniParcId']


def retrieve_uniprot(af, skip_db_list, verbose=False):
    for ac in af['data']:
        skp = False
        for skip_db in skip_db_list:
            if len(af['data'][ac]['databases'][skip_db]) > 0:
                skp = True
        if skp:
            continue
        seq = af['data'][ac]['seq']
        upi = seq_to_uniparc(seq)
        upkd_list, ox_list = uniparc_to_uniprotkb_ids(upi=upi, active_only=True, include_isoforms=True)
        if upi is None:
            upi = 'UPI not found'
        if len(upkd_list) > 0:
            af['data'][ac]['databases']['UniProt'] = upkd_list
            af['data'][ac]['databases']['UniParc'] = [upi]
        else:
            af['data'][ac]['databases']['UniProt'] = []
            af['data'][ac]['databases']['UniParc'] = [upi]
        if verbose:
            print(af['data'][ac]['databases']['UniParc'], af['data'][ac]['databases']['UniProt'], flush=True)

    return af

def get_uniprot_ids_from_uniparc(upi: str, active_only: bool = True):

    """
    Retrieve UniProtKB accessions + OX (taxonId) linked to a UniParc ID.

    Returns a list of dictionaries:
        [
            {"accession": "P04637", "ox": 9606, "organism": "Homo sapiens", "database": "UniProtKB/Swiss-Prot", "active": True},
            ...
        ]
    """
    url = f"https://rest.uniprot.org/uniparc/{upi}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers, timeout=30)
    response.raise_for_status()

    data = response.json()
    uniprot_ids = []

    for xref in data.get("uniParcCrossReferences", []):
        database = xref.get("database", "")
        # UniProtKB appears as "UniProtKB/Swiss-Prot" or "UniProtKB/TrEMBL"
        if "UniProtKB" in database:
            if active_only and not xref.get("active", False):
                continue
            acc = xref.get("id")
            if acc:
                uniprot_ids.append(acc)

    # Remove duplicates while preserving order
    return list(dict.fromkeys(uniprot_ids))


def uniparc_to_uniprotkb_ids(upi: str, active_only: bool = True, include_isoforms: bool = True):
    """
    Retrieve all UniProtKB accessions linked to a UniParc ID.
    ----------
    Parameters:
    upi : str
        UniParc identifier (e.g. "UPI000002ED67")
    active_only : bool
        If True, return only currently active cross-references
    include_isoforms : bool
        If False, exclude isoform accessions (those containing '-')
    ------------------
    Returns: a List of unique UniProtKB accessions (sorted)
    """
    url = f"https://rest.uniprot.org/uniparc/{upi}/databases/stream"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers, timeout=60)
    response.raise_for_status()

    results = response.json().get("results", [])

    entries = []
    # seen = set()  # to avoid exact duplicates (accession + ox)

    for entry in results:
        db = entry.get("database", "")
        if "UniProtKB" not in db:
            continue
        if active_only and not entry.get("active", False):
            continue

        acc = entry.get("id")
        if not acc:
            continue
        if not include_isoforms and "-" in acc:
            continue

        organism = entry.get("organism") or {}
        ox = organism.get("taxonId")

        entries.append({
            "accession": acc,
            "ox": ox,
        })

    ac_set = set()
    ox_set = set()
    for entry in entries:
        ac_set.add(entry.get("accession"))
        ox_set.add(entry.get("ox"))
    return list(ac_set), list(ox_set)

def seq_to_uniparc(sequence: str):
    """
    Return the UniParc ID (UPI...) for an exact protein sequence,
    or None if not found.
    """
    url = "https://www.ebi.ac.uk/proteins/api/uniparc/sequence"
    headers = {
        "Content-Type": "text/plain",
        "Accept": "application/json",
    }
    response = requests.post(url, data=sequence.strip().upper(), headers=headers, timeout=30)

    if response.status_code == 200:
        data = response.json()
        return data.get("accession")  # e.g. "UPI000002ED67"
    elif response.status_code == 404:
        return None  # sequence not present in UniParc
    else:
        response.raise_for_status()
