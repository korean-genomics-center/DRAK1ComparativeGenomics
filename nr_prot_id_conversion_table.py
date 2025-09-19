import json
from collections import OrderedDict


def parse_fasta_headers_from_file(fasta_path):
    fasta_dict = OrderedDict()

    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().split(">")[0].strip()
                parts = header.split(maxsplit=1)
                if len(parts) == 2:
                    acc, desc = parts
                    fasta_dict[acc] = desc
                elif len(parts) == 1:
                    fasta_dict[parts[0]] = ""
    return fasta_dict

def save_dict_to_json(fasta_path, json_output_path):
    fasta_dict = parse_fasta_headers_from_file(fasta_path)
    with open(json_output_path, 'w', encoding='utf-8') as out_json:
        json.dump(fasta_dict, out_json, indent=2, ensure_ascii=False)
    print(f"✅ Saved to: {json_output_path}")



# ==== USAGE ====
fasta_input = "/BiO/Share/Databases/nr_db/nr.fasta"
output_json = "/BiO/Share/Databases/nr_db/nr_db_prot_id_conversion_table.json"
save_dict_to_json(fasta_input, output_json)