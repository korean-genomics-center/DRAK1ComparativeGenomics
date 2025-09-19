# DRAK1 Comparative Genomics

### Content
* [Analysis](#analysis)
  - [Genome Assembly Analysis](#genome-assembly-analysis)
  - [Endogenous Virus Element (EVE) Analysis](#endogenous-virus-element-eve-analysis)
  - [DRAK1/STK17A ORF Detection Analysis](#drak1/stk17a-orf-detection-analysis)
  - [Phylogenetic Tree Analysis](#phylogenetic-tree-analysis)
* [Data and Figures in the Manuscript](#data-and-figures-in-the-manuscript)
  - [Data](#data)
  - [Figures](#figures)

## Analysis
### Genome Assembly Analysis
* Download Genome Assemblies From NCBI GenBank By Taxon ID ([ncbi-datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download))
```bash
# chiroptera: taxon_id=9397
# rodentia: taxon_id=9989
datasets download genome taxon {taxon_id} --include genome,rna,protein,cds,gff3,gtf,gbff,seq-report --assembly-level chromosome,complete,contig,scaffold --filename {taxon_id}.zip --dehydrated
unzip {taxon_id}.zip -d {taxon_id}
datasets rehydrate --directory {taxon_id}
datasets summary genome taxon {taxon_id} > ncbi_datasets_summary_genome_taxon_{taxon_id}.json
```

* Check Assembly Quality from NCBI MetaData ([busco](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download))
```bash

```

* Check Assembly Completeness Using BUSCO ([busco](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/large-download))
```bash
# mode = genome
busco \
  -i {genome_fna} \
  -m {mode} \
  --lineage_dataset {mammalia_odb12} \
  --cpu 40 \
  --out_path {outdir} \
  --out {outfilename} \
  --metaeuk \
  --offline \
  --force"
```

### Endogenous Virus Element (EVE) Analysis
* Softmasking repeats with Repeatmasker
```bash
repeatmasker \
  -pa {threads} \
  -xsmall \
  quick \
  -species mammals \
  {genome_fna} \
  -dir {outdir}
```

* Download Viral Protein Sequences From NCBI Virus ([EEfinder](https://github.com/WallauBioinfo/EEfinder/wiki/Viral-Datasets))

* Download NCBI ID mapping files ([NCBI taxonomy accesion2taxid](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/))
download prot.accession2taxid.gz 

* Download NCBI Taxonomy Database ([NCBI taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy))
```bash
mkdir ncbi-taxdump
cd ncbi-taxdump
tar -xvzf taxdump.tar.gz
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
```

* create virus DB for MMSeq run ([mmseq:createtaxdb](https://github.com/soedinglab/mmseqs2/wiki#create-a-seqtaxdb-by-manual-annotation-of-a-sequence-database))
```bash
mmseqs \
  createtaxdb \
  virusSeqsDB \
  tmp \
  --ncbi-tax-dump ncbi-taxdump \
  --tax-mapping-file taxidmapping \
  2>&1 | tee log_taxon.txt
```

* Detect EVE by aligning with virus DB using MMSeq ([mmseq:easy-taxonomy](https://github.com/soedinglab/mmseqs2/wiki#taxonomy-top-hit-report))
```bash
mmseqs\
  easy-taxonomy \
  {genome_fna} \
  virusSeqsDB \
  {alnRes} \
  {tmpdir} \
  --threads {threads} \
  --mask-lower-case 1
```

### DRAK1/STK17A ORF Detection Analysis
* Make Index File for Genome Fasta using Samtools ([samtools](http://www.htslib.org))
```bash
samtools faidx {genome_fna}"
```

* First-pass Exonerate Run For Exon-Intron-Aware ORF Detection ([exonerate](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-31))
```bash
exonerate \
  --model protein2genome \
  --query {genome_fna} \
  --target {drak1_faa} \
  --showtargetgff yes \
  --showvulgar yes \
  --verbose 3 \
  > {outprefix}"
```

* Parse All Exonerate Outputs Into Bed File
```bash
# hit option: --top5 (for top5 only) or --all (for every) or both
python parse_exonerate_cds_to_bed.py \
  {file_exonerate} \
  {outprefix}.exo.out \
  {hit_option}
```

* Extract Flanking regions From DRAK1/STK17A Protein Sequneces as Bed File
```bash
# flanking regions (e.g. 1000)bp upstream and downstream to bed file
bedtools slop \
  -i {outprefix}.bed \
  -g {genome_fna}.fai \
  -b {flank_region} \
  > {outprefix}_{flank_region}bp_flank.bed

bedtools getfasta \
  -fi {genome_fna} \
  -bed {outprefix}_{flank_region}bp_flank.bed \
  > {outprefix}_{flank_region}bp_flank.fa
```

* Second-pass Exonerate Run For ORF Detection
```bash
exonerate \
  --model protein2genome:bestfit \
  --query {outprefix}_{flank_region}bp_flank.fa \
  --target {drak1_faa} \
  --showtargetgff yes \
  --showvulgar yes \
  --exhaustive yes \
  --bestn 1 \
  --verbose 3 \
  > {file_exonerate}.exo.bestfit.allhit.out
```

* Extract CDS Fasta from exonerate results
```bash
python extractExonerateBestFitCds.py \
  {file_exonerate}.exo.bestfit.allhit.out \
  {outprefix}_{flank_region}bp_flank.fa \
  {outprefix}.cds.fa"
```

* Genewise Run For Codon-aware ORF Detection ([genewise](https://www.ebi.ac.uk/jdispatcher/psa/genewise))
```bash
# genestats, BLOSUM matrix, and codon table found within wise2/wisecfg
genewise \
  {genome_fna} \
  {drak1_faa} \
  -genestats {gene.stat} \
  -matrix {BLOSUM62.bla} \
  -codon {codon.table} \
  -genesf -gff -cdna -pep -pseudo -pretty -ace -both -gener \
  > {outprefix}.genewise.out"

python extract_seq_from_genewise.py \
  --input {outprefix}.genewise.out
```

* Search ORF Best-aligned to DRAK1/STK7A on NR DB using Diamond ([diamond](https://github.com/bbuchfink/diamond))
```bash
diamond makedb --in {nr.fasta} -d {db}

# e.g. sensitivity = sensitive
diamond \
  blastp \
  -d {nr.fasta} \
  -q {fasta_query}.genewise.out.pep.faa \
  -o {blast_out}.genewise.out.pep.faa.blastout \
  --{sensitivity} \
  --threads 3
```

### Phylogenetic Tree Analysis
*Codon-aware Multiple Sequence Alignment Using MACSE ([macse](https://www.agap-ge2pop.org/macse/))
```bash
# ortholog = busco output
macse \
  -prog alignSequences \
  -seq {ortholog}.fna \
  -out_NT {ortholog}.aln.fna \
  -out_AA {ortholog}.aln.faa"
```

*Trim Alignments Using trimAl ([trimal](https://trimal.readthedocs.io/en/latest/))
```bash
trimal \
  -in {ortholog}.aln.trimmed.fna \
  -out {ortholog}.aln.trimmed.fna \
  -htmlout {ortholog}.html \
  -keepheader \
  -automated1
```

*Build Gene Tress using RAXML ([raxml-ng](https://github.com/amkozlov/raxml-ng))
```bash
# boostrap_num = 100
raxml \
  --threads 30 \
  --msa {ortholog}.aln.trimmed.cleaned.fna \
  --model GTR+G \
  --prefix {ortholog} \
  --all \
  --bs-trees {boostrap_num}

cat {dir_raxml}/*.bestTree >> concat.all.raxml.bestTree
```

*ASTRAL4 ([aster:astral4](https://github.com/chaoszhang/ASTER/blob/master/tutorial/astral4.md))
```bash
astral4 \
  -r 4 \
  -s 4 \
  -i concat.all.raxml.bestTree \
  -o astral.out.speciesTree \
  -g mammalian_265_species.pruned.nwk \
  -t 100 \
  2> astral.log
```

*ggtree ([ggtree](https://yulab-smu.top/treedata-book/chapter4.html))
```bash
Rscript 
```

### Positive Selection Test Analysis
*HyPhy ([hyphy](https://github.com/veg/hyphy))
```bash
# path_bf = path to batch file (i.e., REALX.bf)
# test = LOSS (i.e., complete loss)
# ref = INTACT (i.e., intact + partial loss)
# model = standard

hyphy {path_bf}\
  --alignment {ortholog}.aln.cleaned.matched.fna \
  --tree astral.out.pruned.annotated.speciesTree \
  --test {test} \
  --reference {ref} \
  --models {model} \
  --output {ortholog}.RELAX.out \
  | tee -a {ortholog}.RELAX.log
```