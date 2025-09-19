def write_first_fasta_entry(input_path, output_path):
    header = None
    sequence = []
    
    with open(input_path, 'r') as infile:
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header is None:
                    header = line
                else:
                    break  # Stop reading after first entry
            elif header:
                sequence.append(line)
    
    if header:
        with open(output_path, 'w') as outfile:
            outfile.write(header + '\n')
            # Optionally wrap sequence at 60 or 80 characters per line
            seq = ''.join(sequence)
            for i in range(0, len(seq), 60):
                outfile.write(seq[i:i+60] + '\n')
    else:
        print("No FASTA header found.")

input_fa = "/BiO/Share/GenomeAssemblies/taxon_9397/9397/ncbi_dataset/data/GCF_004115265.2/GCF_004115265.2_mRhiFer1_v1.p_genomic.fna"
output_fa = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/genome/GCF_004115265.2_mRhiFer1_v1.p_genomic_first_line_only.fna"
write_first_fasta_entry(input_fa, output_fa)
