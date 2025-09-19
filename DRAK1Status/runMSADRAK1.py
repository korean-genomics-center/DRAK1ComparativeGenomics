# %%
dir_drak1_status = "/BiO/Research/ComparativeGenomics_DRAK1/Results/DRAK1_Status_Pep"

# %%
import subprocess

bin_mafft = "/BiO/Access/kyungwhan1998/miniconda3/envs/mafft/bin/mafft"
path_msa_input_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.faa"
fasta_query = path_msa_input_total
path_msa_output_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.aln.faa"
fasta_aln = path_msa_output_total

cmd = f"{bin_mafft} {fasta_query} > {fasta_aln}"
# subprocess.run(cmd, shell=True)

# %%
bin_prank = "/BiO/Share/Tool/prank/bin/prank"
fasta_querys = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.faa"
fasta_aln = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.prank.aln.faa"

cmd = f"{bin_prank} -d={fasta_query} -o={fasta_aln} -F -protein"
print(cmd)

# %%
path_msa_viz_total = f"{dir_drak1_status}/drak1_present_blastout_tophit_query_orf_longest.aln.png"
bin_pymsaviz = "/BiO/Access/kyungwhan1998/miniconda3/envs/mafft/bin/pymsaviz"
vis_aln = path_msa_viz_total
cmd = f"{bin_pymsaviz} -i {fasta_aln} -o {vis_aln} --color_scheme Taylor --show_consensus --show_count"
print(cmd)

