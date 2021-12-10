# Author: Austin Seamann
# Version: 1.0
# Last Updated: December 10th, 2021
import argparse
import subprocess
import pandas as pd


# Method: grab_all_seq
# Goal: Collect all protein chains from FASTA file submitted.
# Input: Protein sequence FASTA file
# Output: Return dictionary of protein ID and sequence
def grab_all_seq(fasta_in):
    all_seq = {}  # ID:Seq
    with open(fasta_in, "r") as f1:
        current_id = ""
        for line in f1:
            if line[0] == ">":
                current_id = line.split("|")[1]
                all_seq[current_id] = ""
            elif current_id != "":
                all_seq[current_id] += line[:-1]
    return all_seq


# Method: fragment
# Goal: Upon each submitted protein chain, fragment into all possible x-mers.
# Input: Dictionary of protein sequences and fragment size
# Output: Dictionary of protein ID and fragment list
def fragment(seqs, size):
    fragment_dic = {}  # ID: Fragment List
    for current_id in seqs:
        if current_id not in fragment_dic.keys():
            fragment_dic[current_id] = []
        for position in range(0, len(seqs[current_id]) - size + 1):
            fragment_dic[current_id].append(seqs[current_id][position:position + size])
    return fragment_dic


# Method: grab_mhc
# Goal: Create a dictionary of HLAs based on population frequency rank
# Input: XLSX file from Be The Match Registry Haplotype Frequencies Tables
# Output: Dictionary containing Population as Key and HLA IDs for HLAs ranked in the top 50 per population
def grab_mhc(xlsx_files):
    hlas_info = {}  # Population:IDs > 50 rank
    first_file = True
    rank_headers = []
    for file_in in xlsx_files:
        df = pd.read_excel(file_in, engine='openpyxl')
        if first_file:
            headers = df.columns  # All headers in file
            for header in headers:
                if header.endswith("_rank"):
                    if header not in rank_headers:
                        rank_headers.append(header)
            first_file = False
        for header in rank_headers:
            if header not in hlas_info.keys():
                hlas_info[header] = []
            # Cutoff for rank of HLAs submitted
            # May change from rank to a different measure at some point
            df_top = df[df[header] <= 50.0]
            # Send in first letter of file
            hlas_info[header] = hlas_info[header] + df_top[file_in.split('/')[-1].split(".")[0]].values.tolist()
    return hlas_info


# Method: run_mhcflurry
# Goal: Automated submission of HLA and peptide pairs in batch submission to MHCFlurry
# Input:
#   fragments: dictionary of proteome fragments
#   hlas_info: dictionary of HLAs for each populations to be ran
#   supported_alleles: Text file containing supported alleles for MHCFlurry
# Output:
#   output_location: Location of resulting MHCFlurry output file
#   unsupported_alleles: List of alleles that were ranked in population but not in MHCFlurry supported list
def run_mhcflurry(fragments, hlas_info, supported_alleles):
    # Program
    commands = ["mhcflurry-predict"]
    # Keep track of unsupported alleles
    unsupported_alleles = []
    # Append alleles here
    commands.append("--alleles")
    for population in hlas_info:
        for allele in hlas_info[population]:
            if not allele[-1].isnumeric():
                allele = allele[:-1]
            full_allele = "HLA-" + allele
            if full_allele in supported_alleles:
                commands.append(full_allele)
            else:
                print("Unsupported allele: " + full_allele)
                unsupported_alleles.append(allele)  # Remove after for comparison at end
    # Append peptides
    commands.append("--peptides")
    for current_id in fragments:
        for fragment_ in fragments[current_id]:
            commands.append(fragment_)
    # commands.append(peptide_line[:-1])
    # Append output
    output_location = "predictions.csv"
    commands.append("--out")
    commands.append(output_location)
    subprocess.run(commands)
    return output_location, unsupported_alleles


# Method: parse_args
# Goal: Collect command line arguments from user
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Fasta file in containing protein sequences", type=str)
    parser.add_argument("-s", "--size", help="Sizes of fragment", type=int)
    parser.add_argument("-a", help="HLA-A", type=str)
    parser.add_argument("-b", help="HLA-B", type=str)
    parser.add_argument("-c", help="HLA-C", type=str)
    return parser.parse_args()


# Method: main
# Goal: Control the operation of program
def main():
    args = parse_args()
    all_seqs = grab_all_seq(args.fasta)  # Grab all protein sequences and ids
    fragments = fragment(all_seqs, args.size + 1)  # produce fragment size selected
    hla_files = []
    if args.a:
        hla_files.append(args.a)
    if args.b:
        hla_files.append(args.b)
    if args.c:
        hla_files.append(args.c)
    hlas = grab_mhc(hla_files)
    # Supported HLA List
    supported_alleles = []
    with open("supported_alleles.txt", "r") as f:
        for line in f:
            supported_alleles.append(line[:-1])
    # Remove redundant HLAs
    hlas_submit = hlas  # Copy of hlas but with removed redundancy
    long_hla_list = []
    for population in hlas_submit:
        for allele in hlas_submit[population]:
            if allele in long_hla_list:
                hlas_submit[population].remove(allele)
            else:
                long_hla_list.append(allele)
    print("Running MHCflurry...")
    results_final, unsupported_alleles = run_mhcflurry(fragments, hlas_submit, supported_alleles)
    for population in hlas_submit:
        for allele in hlas_submit[population]:
            if allele[-1].isnumeric():
                allele = allele[:-1]
            if allele in unsupported_alleles:
                print("Allele not found: " + allele + " Population: " + population[:-5])
    print("Results: " + results_final)
    print("Done!")


if __name__ == "__main__":
    main()
