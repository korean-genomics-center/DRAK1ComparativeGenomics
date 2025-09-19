#!/usr/bin/env python3
# %%
"""
Comparative Genomics Analysis Pipeline
- Presence-Absence Variation Detection
- Phylogenetic Mapping
- Evolutionary Constraint Analysis
"""

# %%
import json
import os
from io import StringIO
from itertools import chain

import requests
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.PAML import codeml


# %%
class ComparativeGenomics:
    def __init__(self):
        self.ensembl_server = "https://rest.ensembl.org"
        self.mafft = "/BiO/Share/Tool/mafft-7.525-without-extensions/bin/mafft"

    def get_species_info(self):
        """
        Retrieve all species available from Ensembl.
        """
        
        ext = "/info/species?"
        url = self.ensembl_server + ext
        response = requests.get(url, headers={"Content-Type": "application/json"})

        if response.ok:
            decoded = response.json()
            species_dicts = list(map(lambda x: {x["name"]: x}, decoded["species"]))
            species_info = dict(chain(*map(dict.items, species_dicts)))
            
            return species_info
             
        else:
            print(f"[Error] Failed to retrieve species info: {response.status_code}")
            print(f"URL: {url}")       
    
    def get_gene_info(self, gene_id_or_symbol, species="homo_sapiens"):
        """
        Retrieve gene information from Ensembl.
        If species is provided, assumes gene_id_or_symbol is a gene symbol.
        If species is None, assumes it's an Ensembl gene ID.
        """
        if not str(gene_id_or_symbol).startswith("ENS"):
            # Symbol-based lookup
            ext = f"/lookup/symbol/{species}/{gene_id_or_symbol}?expand=1"
        else:
            # ID-based lookup
            ext = f"/lookup/id/{gene_id_or_symbol}?expand=1"

        url = self.ensembl_server + ext
        response = requests.get(url, headers={"Content-Type": "application/json"})

        if response.ok:
            return response.json()
        
        else:
            print(f"[Error] Failed to retrieve gene info: {response.status_code}")
            print(f"URL: {url}")
            
            return None

    def get_orthologs(self, gene_id_or_symbol, species="homo_sapiens"):
        """
        Get orthologs of a gene (by Ensembl ID) in specified species.
        If target_species is None, queries all species in self.species_list.
        """
        if not str(gene_id_or_symbol).startswith("ENS"):
            # Symbol-based lookup
            ext = f"/homology/symbol/{species}/{gene_id_or_symbol}?type=orthologues"
        else:
            ext = f"/homology/id/{species}/{gene_id_or_symbol}?type=orthologues"
     
        url = self.ensembl_server + ext
        response = requests.get(url, headers={"Content-Type": "application/json"})

        if response.ok:
            return response.json()
        else:
            print(f"[Error] Failed to fetch orthologs for {gene_id}: {response.status_code}")
            print(f"URL: {url}")
            return None
        
    def fetch_cds_sequences(self, dict_all_species_gene_ids):
        """Batch retrieve CDS sequences from Ensembl"""
        dict_all_species_cdna_seq = dict()
        for species, gene_id in dict_all_species_gene_ids.items():
            try:
                response = requests.get(
                    f"https://rest.ensembl.org/sequence/id/{gene_id}?",
                    headers={"Content-Type" : "application/json"}
                )
                if response.ok:
                    dict_all_species_cdna_seq[species] = response.json()
                    response.raise_for_status()
                    
            except requests.exceptions.RequestException:
                return 0

        return dict_all_species_cdna_seq
    
    def detect_pav_patterns(self, target_gene, species_list):
        """Detect presence-absence patterns for target gene"""
        pav_matrix = {}
        ortholog_data = self.get_orthologs(target_gene)
        
        if not ortholog_data:
            return None
        
        # Initialize presence matrix
        for species in species_list:
            pav_matrix[species] = {
                'present': False,
                'gene_id': None,
                'confidence': 0,
                'synteny_score': 0
            }
        
        # Process orthology relationships
        for homolog in ortholog_data.get('data', [{}])[0].get('homologies', []):
            target_species = homolog['target']['species']
            if target_species in pav_matrix:
                pav_matrix[target_species]['present'] = True
                pav_matrix[target_species]['gene_id'] = homolog['target']['id']
                pav_matrix[target_species]['confidence'] = homolog.get('dn_ds', 0)
        
        return pav_matrix

    def calculate_synteny_conservation(self, gene_id, species, window_size=10):
        """Calculate syntenic conservation around target gene"""
        gene_info = self.get_gene_info(gene_id, species)
        if not gene_info:
            return 0
        
        chromosome = gene_info['seq_region_name']
        start = gene_info['start'] - (window_size * 50000)
        end = gene_info['end'] + (window_size * 50000)
        
        ext = f"/overlap/region/{species}/{chromosome}:{start}-{end}?feature=gene"
        try:
            response = requests.get(
                self.ensembl_server + ext,
                headers={"Content-Type": "application/json"},
                timeout=10
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException:
            return 0
    
    def align_sequences(self, input_fa, num_threads=20):
        """Perform multiple sequence alignment"""
        from Bio.Align.Applications import MafftCommandline

        mafft_cline = MafftCommandline(self.mafft, input=input_fa, thread=num_threads)
        stdout, stderr = mafft_cline()
        
        return AlignIO.read(StringIO(stdout), "fa")

    def construct_phylogeny(self, sequences_file, model="GTR"):
        """Construct phylogenetic tree from sequence alignments"""
        try:
            phyml_cmd = PhymlCommandline(
                input=sequences_file,
                datatype='nt',
                model=model,
                alpha='e',
                categories=4,
                bootstrap=100
            )
            phyml_cmd()
            tree_file = sequences_file + "_phyml_tree.txt"
            return Phylo.read(tree_file, "newick")
        
        except Exception as e:
            print(f"Phylogeny construction failed: {e}")
            return None

    def calculate_dnds_ratios(self, alignment_file, tree_file):
        """Calculate dN/dS ratios using PAML codeml"""
        try:
            cml = codeml.Codeml(
                alignment=alignment_file,
                tree=tree_file,
                out_file="codeml_output.txt"
            )
            cml.set_options(
                seqtype=1,
                model=0,
                NSsites=[0, 1, 2],
                fix_omega=0,
                omega=0.4
            )
            results = cml.run()
            return {model: results[model]['omega'] for model in results}
        except Exception as e:
            print(f"dN/dS calculation failed: {e}")
            return None

# %%
# Initialize pipeline
cg = ComparativeGenomics()
workdir = "/BiO/Access/kyungwhan1998/comparativegenomics/Results/REST_ENSEMBL"
target_gene = "STK17A"
outdir = os.path.join(workdir, target_gene)
os.makedirs(outdir, exist_ok=True)
species_info = cg.get_species_info()
species_list = list(species_info.keys())

# %%
# Step 0: Get Homolog Data
ortholog_data = cg.get_orthologs(target_gene)

def fetch_homolog_info_all_species(ortholog_data):
    dict_all_species_homolog = dict()
    for homolog in ortholog_data.get('data', [{}])[0].get('homologies', []):
        target_species = homolog['target']['species']
        dict_all_species_homolog[target_species] = homolog['target']
    
    return dict_all_species_homolog

file_homolog = os.path.join(outdir, "homolog_info.json")
if not os.path.exists(file_homolog):
    dict_all_species_homolog = fetch_homolog_info_all_species(ortholog_data)
    with open(file_homolog, "w") as fw:
        json.dump(dict_all_species_homolog, fw)
else:
    with open(file_homolog, "r") as fr:
        dict_all_species_homolog = json.load(fr)

# %%
dict_all_species_gene_ids = {k: v["id"] for k, v in dict_all_species_homolog.items()}

file_sequence = os.path.join(outdir, "sequence_info.json")
if not os.path.exists(file_sequence):
    dict_all_species_cdna_seq = cg.fetch_cds_sequences(dict_all_species_gene_ids)
    with open(file_sequence, "w") as fw:
        json.dump(dict_all_species_cdna_seq, fw)
else:
    with open(file_sequence, "r") as fr:
        dict_all_species_cdna_seq = json.load(fr)

file_fasta = os.path.join(outdir, "sequences.fa")
if not os.path.exists(file_fasta):
    with open(file_fasta, "w") as fw:
        for species, info in dict_all_species_cdna_seq.items():
            fw.write(f">{species}\n")
            fw.write(f"{info['seq']}\n")

# %%
aligned = cg.align_sequences(file_fasta)
file_aligned_fasta = os.path.join(outdir, "mafft_aligned_sequences.fa")
AlignIO.write(aligned, file_aligned_fasta, "fasta")

# %%
# Step 1: Detect presence-absence patterns
# pav_matrix = cg.detect_pav_patterns(target_gene, species_list)
# pav_matrix
pav_matrix = {}
# Initialize presence matrix
for species in species_list:
    pav_matrix[species] = {
        'present': False,
        'gene_id': None,
        'confidence': 0,
        'synteny_score': 0
    }

# Process orthology relationships
for homolog in ortholog_data.get('data', [{}])[0].get('homologies', []):
    target_species = homolog['target']['species']
    if target_species in pav_matrix:
        pav_matrix[target_species]['present'] = True
        pav_matrix[target_species]['gene_id'] = homolog['target']['id']
        pav_matrix[target_species]['confidence'] = homolog.get('dn_ds', 0)

# %%
# Step 2: Calculate synteny conservation
for species in pav_matrix:
    if pav_matrix[species]['present']:
        gene_id = pav_matrix[species]['gene_id']
        pav_matrix[species]['synteny_score'] = cg.calculate_synteny_conservation(gene_id, species)

# %%
file_pav = os.path.join(outdir, "ortholog_info.json")
with open(file_pav, "w") as fw:
    json.dump(pav_matrix, fw)

# # %%
# Step 3: Phylogenetic analysis (example using dummy alignment)
# Note: Replace with actual alignment file path
phylogeny = cg.construct_phylogeny(file_aligned_fasta)

# # %%
# # Step 4: Evolutionary constraint analysis
# dnds_results = cg.calculate_dnds_ratios("alignment.phy", "species_tree.nwk")

# # Generate report
# report = {
#     "pav_matrix": pav_matrix,
#     "dnds_ratios": dnds_results,
#     "synteny_scores": {sp: pav_matrix[sp]['synteny_score'] for sp in pav_matrix}
# }
# # %%
# # Save results
# with open("comparative_genomics_report.json", "w") as f:
#     json.dump(report, f, indent=2)

# print("Analysis completed. Results saved to comparative_genomics_report.json")
# %%
