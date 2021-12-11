# Author: Austin Seamann
# Version: 1.0
# Last Updated: December 10th, 2021
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import argparse
import pandas as pd
import os


# GLOBAL #
codes = {"Race": ["AAFA", "AFB", "AINDI", "AISC", "ALANAM", "AMIND", "CARB", "CARHIS", "CARIBI",
                  "EURCAU", "FILII", "HAWI", "JAPI", "KORI", "MENAFC", "MSWHIS", "NCHI", "SCAHIS",
                  "SCAMB", "SCSEAI", "VIET"],
         "Broad_race": ["AFA", "API", "CAU", "HIS", "NAM"]
         }

# Change to adjust what category of race codes are used
code_used = codes["Broad_race"]


# Method: process_file
# Goal: Generate working file for generating plots - Incorporate additional information to MHCFlurry output
# Input:
#   prediction_file: MHCFlurry raw output file
#   new_name: Name of the newly processed file
#   hla_inverse: Allele:Population List
#   hlas_freq: Allele:{Population:freq}
# Output: Produces processed file for later analysis
def process_file(prediction_file, new_name, hla_inverse, hlas_freq="..."):
    global code_used
    first_line = True
    # Updated csv to include Population column
    with open(new_name, "w") as out_file:
        with open(prediction_file, "r") as in_file:
            for line in in_file:
                if first_line:
                    # New header
                    out_file.write(line[:-1] + ",Population,Freq_in_Population,Allele_ID\n")
                    first_line = False
                else:
                    # Process each row
                    hla = line.split(",")[0].split("-")[1]
                    if hla in hla_inverse.keys():
                        new_line = line[:-1]
                        # Adding multiple lines if HLA is in several populations
                        for hla_population in hla_inverse[hla]:
                            if hla_population in code_used:
                                out_file.write(new_line + "," + hla_population + ","
                                               + str(hlas_freq[hla][hla_population]) + ","
                                               + str(hla + "_" + hla_population) + "\n")


# Method: box_plot
# Goal: Generate box plot for non-weighted number of binders per each HLA
def box_plot(prediction_file, hla_inverse, save, width_, height_):
    new_name = prediction_file[:-4] + "_tmp.csv"
    # Process file to add back repeated HLAs per population removed prior
    process_file(prediction_file, new_name, hla_inverse)
    # Grab data needed for plot
    csv_in = pd.read_csv(new_name).sort_values("Population")
    csv_in = csv_in[csv_in["mhcflurry_affinity"] <= 500].sort_values("Population")
    # Create boxplot
    sns.set(rc={"figure.figsize": (width_, height_)})
    ax = sns.boxplot(x="Population", y="mhcflurry_affinity", data=csv_in)
    # Collect data for number of high binders
    medians = csv_in.groupby(['Population'])['mhcflurry_affinity'].median().values
    count_values = csv_in['Population'].value_counts().values
    count_values = [str(x) for x in count_values.tolist()]
    count_values = ["n: " + i for i in count_values]
    pos = range(len(count_values))
    # Add number of high binders to each population
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick], medians[tick] + 0.03, count_values[tick], horizontalalignment='center',
                size='x-small', color='b', weight='semibold')
    if save != "...":
        plt.savefig(save, dpi=300)
    else:
        plt.show()


# Method: box_plot_freq
# Goal: Generate box plot based on weighted HLA frequency binder score
def box_plot_freq(prediction_file, hla_inverse, hlas_freq, save, width_, height_):
    global code_used
    new_name = prediction_file[:-4] + "_tmp.csv"
    # Process file to add back repeated HLAs per population removed prior
    process_file(prediction_file, new_name, hla_inverse, hlas_freq)
    # Grab data needed for plot
    csv_in = pd.read_csv(new_name)
    affinity = 500
    csv_in = csv_in[csv_in["mhcflurry_affinity"] <= affinity].sort_values("Allele_ID")
    with open("tmp.csv", "w") as f:
        f.write("Allele_ID,Allele,Population,Weighted_Value\n")
        for allele_id in csv_in['Allele_ID'].unique():
            # Get allele_id specific dataframe
            allele_id_df = csv_in[csv_in["Allele_ID"] == allele_id]
            # Collect count of allele detection of peptides
            # Number of unique peptides in list
            count_allele_id = allele_id_df["peptide"].nunique()
            # Collect freq information
            single_row = allele_id_df.iloc[0]
            freq_allele = single_row['Freq_in_Population']
            # Calculated weighted score
            weighted_score = count_allele_id * freq_allele
            # Other info
            population = allele_id.split("_")[1]
            allele = allele_id.split("_")[0]
            f.write(allele_id + "," + allele + "," + population + "," + str(weighted_score) + "\n")
    csv_in = pd.read_csv("tmp.csv")
    csv_in = csv_in.sort_values("Population")
    # Create boxplot
    sns.set(rc={"figure.figsize": (width_, height_)})
    title_ = "HLA binders in " + str(affinity) + "nM or better"
    # Create plot
    sns.set_style("white")
    ax = sns.boxplot(data=csv_in, x="Population", y="Weighted_Value")
    plt.title(title_, fontsize=20)
    plt.ylabel('Weighted Score', fontsize=20)
    plt.xlabel('Allele Binders by Population', fontsize=20)
    ax.set_yticklabels(ax.get_yticks(), size=15)
    if save != "...":
        plt.savefig(save, dpi=300)
    else:
        plt.show()
    os.remove("tmp.csv")


# Method: Histogram - Test method
def histogram(prediction_file, hla_inverse, save, width_, height_):
    global code_used
    new_name = prediction_file[:-4] + "_tmp.csv"
    # Process file to add back repeated HLAs per population removed prior
    process_file(prediction_file, new_name, hla_inverse)
    # Grab data needed for plot
    csv_in = pd.read_csv(new_name)
    affinity = 500
    csv_in = csv_in[csv_in["mhcflurry_affinity"] <= affinity].sort_values("Population")
    # Collect count
    count_values = csv_in['Population'].value_counts().values
    population_info = {}  # Population: Count
    count = 0
    for population in code_used:
        population_info[population] = count_values[count]
        count += 1
    # Create boxplot
    sns.set(rc={"figure.figsize": (width_, height_)})
    title_ = "HLA binders in " + str(affinity) + "nM and below group"
    plt.bar(population_info.keys(), population_info.values())
    plt.title(title_)
    plt.ylabel('Count')
    plt.xlabel('Population')
    if save != "...":
        plt.savefig(save, dpi=300)
    else:
        plt.show()


# Method: scatter_freq
# Goal: Generate scatter plot based on weighted HLA frequency binder score
def scatter_freq(prediction_file, hla_inverse, hlas_freq, save, width_, height_):
    global code_used
    new_name = prediction_file[:-4] + "_tmp.csv"
    # Process file to add back repeated HLAs per population removed prior
    process_file(prediction_file, new_name, hla_inverse, hlas_freq)
    # Grab data needed for plot
    csv_in = pd.read_csv(new_name)
    affinity = 500
    csv_in = csv_in[csv_in["mhcflurry_affinity"] <= affinity].sort_values("Allele_ID")
    with open("tmp.csv", "w") as f:
        f.write("Allele_ID,Allele,Population,Weighted_Value\n")
        for allele_id in csv_in['Allele_ID'].unique():
            # Get allele_id specific dataframe
            allele_id_df = csv_in[csv_in["Allele_ID"] == allele_id]
            # Collect count of allele detection of peptides
            # Number of unique peptides in list
            count_allele_id = allele_id_df["peptide"].nunique()
            # Collect freq information
            single_row = allele_id_df.iloc[0]
            freq_allele = single_row['Freq_in_Population']
            # Calculated weighted score
            weighted_score = count_allele_id * freq_allele
            # Other info
            population = allele_id.split("_")[1]
            allele = allele_id.split("_")[0]
            f.write(allele_id + "," + allele + "," + population + "," + str(weighted_score) + "\n")
    csv_in = pd.read_csv("tmp.csv")
    csv_in = csv_in.sort_values("Population")
    # Create boxplot
    sns.set(rc={"figure.figsize": (width_, height_)})
    title_ = "HLA binders in " + str(affinity) + "nM or better"
    # Determine positions for tick marks
    total_values = len(csv_in.index)
    section = total_values / len(code_used)
    positions = [section/2]  # First position in middle of section
    for each in range(1, len(code_used)):
        positions.append(positions[each - 1] + section)
    print(code_used)
    with sns.axes_style("white"):
        ax = sns.scatterplot(data=csv_in, x="Allele_ID", hue="Population",
                             y="Weighted_Value", legend=False)
        ax.set_yticklabels(ax.get_yticks(), size=15)
        ax.set_xticklabels(ax.get_xticks(), size=15)
        ax.xaxis.set_major_locator(ticker.FixedLocator(positions))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(code_used))
        ax.set_ylabel("Weighted Score", fontsize=20)
        ax.set_xlabel("Allele", fontsize=20)
        ax.set_title(title_, fontsize=20)
    if save != "...":
        plt.savefig(save, dpi=300)
    else:
        plt.show()
    os.remove("tmp.csv")


# Method: grab_mhc
# Goal: Generate dictionaries needed for analysis
# Input: XLSX file from Be The Match Registry Haplotype Frequencies Tables
# Output:
#   hlas_info: Population:IDs > rank
#   hlas_inverse: Allele:Population List
#   hlas_freq: Allele:{Population:freq}  ex. hlas_freq["C*07:01"]["NAM"]["0.0000003"]
def grab_mhc(xlsx_files):
    hlas_info = {}  # Population:IDs > 50 rank
    hlas_inverse = {}  # Allele:Population List
    hlas_freq = {}  # Allele:{Population:freq}  ex. hlas_freq["C*07:01"]["NAM"]["0.0000003"]
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
            df_top = df[df[header] <= 25.0]
            # Collect list of HLA IDs
            hlas_info[header] = hlas_info[header] + df_top[file_in.split('/')[-1].split(".")[0]].values.tolist()
            # Collect HLA frequency for Population
            alleles_list = df_top[file_in.split('/')[-1].split(".")[0]].values.tolist()
            freq_list = df_top[header.split("_")[0] + "_freq"].values.tolist()
            for position in range(len(alleles_list)):
                allele = alleles_list[position]
                # Collect allele listed if not already key
                if not allele[-1].isnumeric():  # Remove letters at end of HLA id
                    allele = allele[:-1]
                if allele not in hlas_freq.keys():
                    hlas_freq[allele] = {}
                hlas_freq[allele][header.split("_")[0]] = freq_list[position]
    # Create hlas_inverse
    for population in hlas_info:
        for hla in hlas_info[population]:
            if not hla[-1].isnumeric():  # Remove letters at end of HLA id
                hla = hla[:-1]
            if hla not in hlas_inverse.keys():  # Start new key in dictionary
                hlas_inverse[hla] = []
            hlas_inverse[hla].append(population.split("_")[0])
    return hlas_info, hlas_inverse, hlas_freq


# Method: parse_args
# Goal: Collect command line arguments from user
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("mhcflurry_csv", help="MHCflurry output csv", type=str)
    parser.add_argument("-a", help="HLA-A XLSX file from Be The Match", type=str)
    parser.add_argument("-b", help="HLA-B XLSX file from Be The Match", type=str)
    parser.add_argument("-c", help="HLA-C XLSX file from Be The Match", type=str)
    parser.add_argument("--code", help="Change from Broad_race codes to Race codes", action="store_true",
                        default=False)
    parser.add_argument("-s", help="Save plot", type=str, default="...")
    parser.add_argument("--width", help="Width of plot area, if saved", type=float, default=11)
    parser.add_argument("--height", help="Height of plot area, if saved", type=float, default=8.5)
    parser.add_argument("--box", help="Box plot of prediction percentile vs population", action="store_true",
                        default=False)
    parser.add_argument("--histo", help="Histogram plot of count of 90 percentile vs population", action="store_true",
                        default=False)
    parser.add_argument("--scatterF", help="Scatter plot of weighted HLA binder scores", action="store_true",
                        default=False)
    parser.add_argument("--boxF", help="Box plot of weighted HLA binder scores", action="store_true",
                        default=False)
    return parser.parse_args()


# Method: main
# Goal: Control the operation of program
def main():
    args = parse_args()
    # Change race code
    global code_used
    if args.code:
        code_used = codes["Race"]
    # Collect HLA information
    hla_files = []
    if args.a:
        hla_files.append(args.a)
    if args.b:
        hla_files.append(args.b)
    if args.c:
        hla_files.append(args.c)
    # hla - Population:IDs > rank
    # hlas_inverse - Allele:Population List
    # hlas_frq - Allele:{Population:freq}
    hlas, hlas_inverse, hlas_freq = grab_mhc(hla_files)
    # Supported HLA List
    supported_alleles = []
    with open("supported_alleles.txt", "r") as f:
        for line in f:
            supported_alleles.append(line[:-1])
    if args.box:
        box_plot(args.mhcflurry_csv, hlas_inverse, args.s, args.width, args.height)
    if args.histo:
        histogram(args.mhcflurry_csv, hlas_inverse, args.s, args.width, args.height)
    if args.scatterF:
        scatter_freq(args.mhcflurry_csv, hlas_inverse, hlas_freq, args.s, args.width, args.height)
    if args.boxF:
        box_plot_freq(args.mhcflurry_csv, hlas_inverse, hlas_freq, args.s, args.width, args.height)


if __name__ == "__main__":
    main()
