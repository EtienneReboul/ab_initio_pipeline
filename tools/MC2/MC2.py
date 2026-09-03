import sys
from param import aff_path
if aff_path not in sys.path:
    sys.path.append(aff_path)

from lib import *
from annotated_fasta import *
from datetime import datetime


if __name__ == '__main__':
    create_out_dir()
    mdl_list, priors_list = load_models(models_path=f'Models/', priors_file='Stuff/priors.tsv')
    fasta = aff_load_fasta(input_fasta)
    p_matrix = read_p_matrix_dict1()
    c_dict = read_f5_dict()
    with open(f"{output_path}timing.csv", 'w') as fout:
        for ac in fasta['data']:
            start_time = datetime.now()
            print(f"MoRFchibi 2.0 processing: {ac}", flush=True)
            score_ensemble(in_ac=ac, in_seq=fasta['data'][ac]['seq'], c_dict=c_dict, p_matrix=p_matrix,
                           models_list=mdl_list, priors_list=priors_list)
            end_time = datetime.now()
            td = end_time - start_time
            print(f"{ac} , {round(td.total_seconds() * 1000)}", file=fout, flush=True)
            # break

