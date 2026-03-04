import argparse
import re

# Collect all protein IDs from the BLAST-like file
def collect_query_ids(input_file):
    query_ids = set()
    with open(input_file, "r") as infile:
        for line in infile:
            if line.strip():
                cols = re.split(r'\s+', line.strip())
                if len(cols) >= 2:
                    query_ids.add(cols[1])
    return query_ids

# Load only relevant accession->taxid mappings
def load_accession2taxid(mapping_file, query_ids):
    acc2taxid = {}
    with open(mapping_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                acc_ver, taxid = parts[1], parts[2]
                if acc_ver in query_ids:
                    acc2taxid[acc_ver] = taxid
    return acc2taxid

# Load taxonomy dump (nodes and names)
def load_taxdump(nodes_file, names_file):
    taxid_to_parent = {}
    taxid_to_name = {}

    with open(nodes_file, "r") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2:
                taxid, parent_taxid = parts[0], parts[1]
                taxid_to_parent[taxid] = parent_taxid

    with open(names_file, "r") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                taxid_to_name[parts[0]] = parts[1]

    return taxid_to_parent, taxid_to_name

# Build lineage from taxid
def get_lineage(taxid, taxid_to_parent, taxid_to_name):
    lineage = []
    while taxid in taxid_to_parent and taxid != "1":
        name = taxid_to_name.get(taxid, "")
        if name:
            lineage.append(name)
        taxid = taxid_to_parent[taxid]
    return "; ".join(reversed(lineage)) if lineage else "NA"

# Annotate BLAST-like file with taxonomy lineage
def process_file(input_file, output_file, acc2taxid, taxid_to_parent, taxid_to_name):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.strip():
                cols = re.split(r'\s+', line.strip())
                if len(cols) >= 2:
                    protein_id = cols[1]
                    taxid = acc2taxid.get(protein_id, "NA")
                    if taxid != "NA":
                        lineage = get_lineage(taxid, taxid_to_parent, taxid_to_name)
                    else:
                        lineage = "NA"
                    outfile.write(f"{line.strip()}\t{lineage}\n")
                else:
                    outfile.write(f"{line.strip()}\tNA\n")
    print(f"Annotated file saved to {output_file}")

# Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate BLAST-like file with taxonomy lineage using local NCBI taxdump.")
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    parser.add_argument("--acc2taxid", default="/N/scratch/lcaesar/prot.accession2taxid")
    parser.add_argument("--nodes", default="/N/scratch/lcaesar/ncbi_taxdump/nodes.dmp")
    parser.add_argument("--names", default="/N/scratch/lcaesar/ncbi_taxdump/names.dmp")
    args = parser.parse_args()

    print("Collecting query IDs...")
    query_ids = collect_query_ids(args.input_file)
    print(f"   Found {len(query_ids)} unique protein IDs in input")

    print("Loading accession→taxid mappings...")
    acc2taxid = load_accession2taxid(args.acc2taxid, query_ids)
    print(f"   Loaded {len(acc2taxid)} matching IDs from {args.acc2taxid}")

    print("Loading taxonomy dump...")
    taxid_to_parent, taxid_to_name = load_taxdump(args.nodes, args.names)
    print("   Taxdump loaded")

    print("Annotating input file...")
    process_file(args.input_file, args.output_file, acc2taxid, taxid_to_parent, taxid_to_name)

