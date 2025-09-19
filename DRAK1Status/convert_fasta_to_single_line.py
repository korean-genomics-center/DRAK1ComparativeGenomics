# %%
def convert_fasta_to_single_line(input_file, output_file, num_lines=None):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        seq = ''
        header = None
        cnt = 0
        for line in infile:
            if cnt > 10:
                break
            line = line.strip()
            
            if line.startswith('>'):
                if header is not None:
                    outfile.write(f"{header}\n{seq}\n")
                    cnt += 1
                    if cnt >= num_lines:
                        break

                header = line
                seq = ''
            
            else:
                seq += line
            

        if header is not None:
            outfile.write(f"{header}\n{seq}\n")

# %%
input = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/genome/protein.faa"
output = "/BiO/Access/kyungwhan1998/comparativegenomics/Resources/Data/genome/protein_modified.faa"
num_lines = 10
convert_fasta_to_single_line(input, output, num_lines=num_lines)

# %%
