The idea of the new repo is to generalize the current repos for ab initio predictions that share the same backbone/skeleton to a cookie cutter template that is easily customisable and runnable on a new data , the overall idea remains the same a 3 stage pipeline : 

1\) preprocessing  
2\) processing   
3\) postprocessing  

The overall split is : 

1\) preprocessing : 

1) Offload the MSA \+ template search to mmseq websever   
2) Make a manual or automatic domain annotation based on interproscan and prosite   
3) Use disopred AIupred to make disordered prediction  
4) Use Anchor2 MorfCHibi for Morf annotation   
5) Generate input for abcfold   
     
6) processing :   
1) Switch from sbatch submitting script to snakemake workflow with snakemake plugin   
2) Run abcfold with a customisable set of models and randomized seed strategy \+ parametrisable recycle samples etc…  
3) Metadata compression 

3\) postprocessing 

1) Alignment and kebesh on Calpha for with customizable list of anchor to determine rigid cores   
2) Chimerax Minimisation procedure with early  stop criterion for numeric explosion/garbage  
3) Clustering using hdbscan/optuna on Ca aligned  
4) Computing the IPSAE , ilis and Pinc metrics   
5) I want a plip run on minimized structure with pliparser for conversion into csv

For each stage I want to use the report option  from snakemake to replace mostly what visualisation the notebook would do  so that everything is nicely put in an html file at the end. 

These are the ideas I have for the report, but first I want to make clear that the goal is to have svg pdf  not png or jpeg.  

MMseq2 server :  i want sequence coverage like colabfold to show with the depth , i want also table with [https://datavzrd.github.io/](https://datavzrd.github.io/) to show pdb entries for selected templates and coverage statistics  of sequence in template 

For the domain annotation I want a figure with the traditional bar which is splitted in different sections corresponding to the different domain/linker.   
For the disordered prediction the graph i want supperposed line plot where i can see the different curves  and section where it has a strong signal

For Generate input for abcfold  : i want  a dumb statistic table showed where you break down the size of input files the number etc etc. use [https://datavzrd.github.io/](https://datavzrd.github.io/) for visualisation 

For the run of abcfold i think the automatic timer is enough but if you have a clever proposition i would be interested , one think i have in mind is monitor cpu/GPU usage but that is too dependent on the cluster architecture to be reliably assessed so i don’t know if it is a good idea or not, and the polling system needs to run at the same time as the actual jobs so i don’t how feasible it would be

For the metadata compression i would like to have table with [https://datavzrd.github.io/](https://datavzrd.github.io/) space saved more than anything else. For kebsh and alignment i just want some violin plot to see the RMSD spread an lineplot that underline which par are considered rigid core and which are not 

For chimerax minimisation I want to report a failure rate ratio in the table per backend, and I also want a line plot from seaborn where you can see the energy minimisation per backend. 

For clustering and rescoring metrics i want a customisable set of Dimensionality reduction (PCA,t-sne,umap, MDS , etc…) with different mode  coloration mode : Backend per rescoring metrics  and Clustering annotation

For plip i want different heatmaps , basically with a customisable set : 

- Per backend   
- Per cluster  
- Total 

For prot-prot interaction i want per domain, for nucleic acid i want domain with protein and index for nucleic acid partner , for nucleic acids and protein with small ligand i want heatmap per atom number and index or domain  heatmap