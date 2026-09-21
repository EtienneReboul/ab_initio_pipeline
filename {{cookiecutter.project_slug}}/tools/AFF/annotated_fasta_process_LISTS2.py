from annotated_fasta import *
# from miscellaneous import is_float, get_xml_root


def list_s2_to_caid(in_file, list_path, features_path):
    print(in_file, list_path)
    ac = ''
    fout = open('junk.junk', 'w')
    fout_f = open('junk.junk2', 'w')
    with open(in_file, 'r') as fin:
        for line in fin:
            line = line.strip()
            if len(line) < 5:
                continue
            if 'Submitted' in line:
                ac = line.split()[0]
                fout = open(f"{list_path}{ac}.caid", 'w')
                fout_f = open(f"{features_path}{ac}.caid", 'w')
                print(f">{ac}", file=fout)
                print(f">{ac}", file=fout_f)
                continue
            if 'Conservation' in line or 'Failed' in line:
                continue
            lst = line.split()
            print(f"{lst[0]}\t{lst[1]}\t{lst[3]}", file=fout)
            print(f"{lst[0]}\t{lst[1]}", end='', file=fout_f)
            # print(lst, flush=True)
            for ii in range(2, 23):
                print(f"\t{lst[ii]}", end='', file=fout_f)
            print(file=fout_f)