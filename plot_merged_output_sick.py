
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#This file plot the data from the arriba output merged files

def process_file(file_path):
    df = pd.read_csv(file_path, delimiter="\t")
    # Perform operations or analysis on the data frame
    df.rename(columns={'#gene1': 'gene1'}, inplace=True)
    df['gene1'] = df['gene1'].str.split(',')
    df['gene2'] = df['gene2'].str.split(',')
    df['multiple_genes'] = (df['gene1'].str.len()) * (df['gene2'].str.len())
    df = df.explode('gene1').explode('gene2')

    df = df.reset_index(drop=True)

    df['gene1'] = df['gene1'].str.replace(r'\(\d+\)', '')

    df['gene2'] = df['gene2'].str.replace(r'\(\d+\)', '')

    df['fusion_name'] = df['gene1'] + "--" + df['gene2']
    file_name = os.path.basename(file_path)
    if "discarded" in file_name:
        df['discarded'] = "yes"
    else:
        df['discarded'] = "no"
    df['sample_id'] = file_name.split(".fusions")[0]
    # filter out where site1 or site 2 are intron
    df = df[~df['site1'].isin(['intron'])
            & ~df['site2'].isin(['intron'])]
    chr1 = df['breakpoint1'].str.extract(r'(chr(\w+))')[0]
    chr2 = df['breakpoint2'].str.extract(r'(chr(\w+))')[0]
    df['chr'] = chr1
    # filter out the fusions of diffrenet chomosomes.
    filtered_df = df[chr1 == chr2]
    # filter rows where gene1 and gene2 the same and both sites are intergenic
    filtered_df = filtered_df[(filtered_df['gene1'] != filtered_df['gene2']) | (
        (filtered_df['site1'] != 'intergenic') & (filtered_df['site2'] == 'intergenic'))]
    # define the start as breakpint1 take what after the  :
    filtered_df['start'] = filtered_df['breakpoint1'].str.split(
        ':').str[1].astype(int)
    filtered_df['end'] = filtered_df['breakpoint2'].str.split(
        ':').str[1].astype(int)
    filtered_df['strand'] = filtered_df['strand1(gene/fusion)'].str.split(
        '/').str[1]
    filtered_df['start'], filtered_df['end'] = np.where(filtered_df['start'] > filtered_df['end'], [
        filtered_df['end'], filtered_df['start']], [filtered_df['start'], filtered_df['end']])
    filtered_df['diff'] = filtered_df['end'] - filtered_df['start']
    filtered_df['diff'] = filtered_df['diff'].abs()
    # filter out all the diff that are grader than 0.5MB
    filtered_df = filtered_df[filtered_df['diff'] < 500000]
    # filter out all the fusions that are translocation or inversion
    filtered_df = filtered_df[~filtered_df['type'].str.contains(
        'translocation|inversion')]
    if 'discarded' in file_name:
        filtered_df = filtered_df[~filtered_df['filters'].str.contains(
            'blacklist')]
    filtered_df = filtered_df.drop_duplicates(subset=['chr','start','end','strand'])
    # drop rowd where chr is M
    filtered_df = filtered_df[~filtered_df['chr'].str.contains('M')]
    filtered_df = filter_out_gtex_Star(filtered_df)
    return filtered_df


def filter_out_gtex_Star(data):
    # filter out the fusions that are in gtex
    gtex_fusion = []
    with open('/private/projects/kidney_rtt/bed_arriba_output/fusion_star.txt') as f:
       # add each fusion name to the list
        for line in f:
            gtex_fusion.append(line.split("\t")[0])

    filtered_df = data[~data['fusion_name'].isin(gtex_fusion)]
    return filtered_df


def plot_conf_all_samples(all_samples):
    fig, ax = plt.subplots()
    diff_low = []
    diff_medium = []
    diff_high = []
    for sample in all_samples:
        diff_low.append(sample[sample['confidence'] == 'low']['diff'])
        diff_medium.append(sample[sample['confidence'] == 'medium']['diff'])
        diff_high.append(sample[sample['confidence'] == 'high']['diff'])

    low = pd.concat(diff_low)
    medium = pd.concat(diff_medium)
    high = pd.concat(diff_high)

    if not low.empty:
        plt.hist(low, bins=50, alpha=0.3, color='r', label='Low Confidence')
    if not medium.empty:
        plt.hist(medium, color='b', bins=50,
                 alpha=0.3, label='Medium Confidence')
    if not high.empty:
        plt.hist(high, bins=50, color='g', alpha=0.3, label='High Confidence')

    ax.set_xlabel("Length of Fusion")
    ax.set_ylabel("Number of Fusions")
    ax.set_title("Frequency Distribution of Fusion Lengths")

    # Add legend with labels for each color
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    fig.savefig('/private/projects/kidney_rtt/all_plots/only_sick_plot/conf.png')


def count_fusion_for_top_samples(data_list):
    fusion_counts = {}

    for df in data_list:
        if len(df) == 0:
            continue
        else:
            sample_id = df['sample_id'].iloc[0]
            fusion_counts[sample_id] = fusion_counts.get(
                sample_id, 0) + len(df)

    sorted_fusion_counts = sorted(fusion_counts.items(),
                                  key=lambda x: x[1], reverse=True)

    top_genes = [gene for gene, count in sorted_fusion_counts[:20]]
    top_counts = [count for gene, count in sorted_fusion_counts[:20]]

    fig, ax = plt.subplots(figsize=(10, 12))

    # Create bar plot with sorted values
    ax.bar(top_genes, top_counts, color='g')
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Number of Fusions")
    ax.set_title("Number of Fusions for Each Sample")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()
    fig.savefig(
        '/private/projects/kidney_rtt/all_plots/only_sick_plot/count_fusion_for_top_sample.png')


def count_fusion_for_sample(data_list):
    fusion_counts = {}

    for df in data_list:
        if len(df) == 0:
            continue
        else:
            sample_id = df['sample_id'].iloc[0]
            fusion_counts[sample_id] = fusion_counts.get(
                sample_id, 0) + len(df)

    sorted_fusion_counts = sorted(fusion_counts.items(),
                                  key=lambda x: x[1], reverse=True)

    top_genes = [gene for gene, count in sorted_fusion_counts]
    top_counts = [count for gene, count in sorted_fusion_counts]

    fig, ax = plt.subplots(figsize=(10, 12))

    # Create bar plot with sorted values
    ax.bar(top_genes, top_counts, color='g')
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Number of Fusions")
    ax.set_title("Number of Fusions for Each Sample")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()
    fig.savefig(
        '/private/projects/kidney_rtt/all_plots/only_sick_plot/count_fusion_for_sample.png')


def count_fusion_for_gene(data_list):
    fig, ax = plt.subplots(figsize=(10, 12))
    fusion_counts = {}
    # Initialize an empty dictionary to store the gene counts
    for df in data_list:
        # Extract the fusion names from the 'fusion-name' column
        fusion_names = df['fusion_name']

        # Count the occurrences of each fusion name
        for fusion_name in fusion_names:
            fusion_counts[fusion_name] = fusion_counts.get(fusion_name, 0) + 1

    # Sort the gene counts in descending order
    sorted_counts = sorted(fusion_counts.items(),
                           key=lambda x: x[1], reverse=True)

    # Get the top 300 genes and their counts
    top_genes = [gene for gene, count in sorted_counts]
    top_counts = [count for gene, count in sorted_counts]

    ax.bar(top_genes, top_counts, color='g')
    ax.set_xlabel("Fusion Name")
    ax.set_ylabel("Count")
    ax.set_title("Frequency Distribution of Fusions")
    ax.set_xticklabels(top_genes, fontsize=5)
    plt.xticks(rotation=90)
    fig.savefig(
        '/private/projects/kidney_rtt/all_plots/only_sick_plot/count_fusion_for_gene.png')


def count_fusion_for_gene_top(data_list):
    fig, ax = plt.subplots(figsize=(10, 12))
    fusion_counts = {}
    # Initialize an empty dictionary to store the gene counts
    for df in data_list:
        # Extract the fusion names from the 'fusion-name' column
        fusion_names = df['fusion_name']

        # Count the occurrences of each fusion name
        for fusion_name in fusion_names:
            fusion_counts[fusion_name] = fusion_counts.get(fusion_name, 0) + 1

    # Sort the gene counts in descending order
    sorted_counts = sorted(fusion_counts.items(),
                           key=lambda x: x[1], reverse=True)

    # Get the top 300 genes and their counts
    top_genes = [gene for gene, count in sorted_counts[:20]]
    top_counts = [count for gene, count in sorted_counts[:20]]

    ax.bar(top_genes, top_counts, color='g')
    ax.set_xlabel("Fusion Name")
    ax.set_ylabel("Count")
    ax.set_title(
        "Frequency Distribution of Top Fusions")
    ax.set_xticklabels(top_genes, fontsize=5)
    plt.xticks(rotation=90)
    fig.savefig(
        '/private/projects/kidney_rtt/all_plots/only_sick_plot/count_fusion_for_top_gene.png')
    

def site1_site2_filters_read_through(data_list):
    all_read_through = []
    for df in data_list:
        if len(df) == 0:
            continue
        else:
            filters_read_throgh= df[(df['filters'].str.contains('read_through'))]
            all_read_through.append(filters_read_throgh)

    combined_df = pd.concat(all_read_through)
    #in the combined_df i have site1 column and site2 column i want to create a new column called site1_site2
    combined_df['site1_site2'] = combined_df['site1'] + '_' + combined_df['site2']
    #now i want to count the number of times each site1_site2 appears in the combined_df
    site1_site2_counts = combined_df['site1_site2'].value_counts()
    #now i want to plot the site1_site2_counts
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.bar(site1_site2_counts.index, site1_site2_counts.values, color='g')
    ax.set_xlabel("Site1_Site2")
    ax.set_ylabel("Count")
    ax.set_title(
        "Frequency Distribution of Site1_Site2 in read_through Filters")
    ax.set_xticklabels(site1_site2_counts.index, fontsize=5)
    plt.xticks(rotation=90)
    fig.savefig(
        '/private/projects/kidney_rtt/all_plots/only_sick_plot/site1_site2_read_through.png')
    

def main():
    gf_call_list = []
    # Directory path containing the files
    directory = "/private/projects/kidney_rtt/sick_samples"
    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        # Check if the item in the directory is a file
        if os.path.isfile(file_path):
            gf_call_one = process_file(file_path)
        gf_call_list.append(gf_call_one)

    combined_df_all = []
    done = []
    for one_df in gf_call_list:
        for other_df in gf_call_list:
            # check ig one_df and other_df are not empty
            if (one_df.empty or other_df.empty):
                continue
            if (one_df['sample_id'].iloc[0] == other_df['sample_id'].iloc[0] and one_df['discarded'].iloc[0] != other_df['discarded'].iloc[0]):
                if (one_df['sample_id'].iloc[0] not in done):
                    list_df = [one_df, other_df]
                    combined_df_one = pd.concat(list_df, ignore_index=True)
                    combined_df_all.append(combined_df_one)
                    done.append(one_df['sample_id'].iloc[0])

    

        
    for one_df in combined_df_all:
        bed_name = (one_df['sample_id'].iloc[1]) + ".bed"
        if one_df.empty:
            continue
        fusion_name = one_df['fusion_name']
        chr_fusion = one_df['chr']
        strand = one_df["strand"]
        start = one_df["start"]
        end = one_df["end"]
        score_series = one_df["multiple_genes"]
        # score_series = pd.Series([0] * len(chr_fusion))
        directory = "/private/projects/kidney_rtt/filters_arriba_bed/"
        bed_path = directory + str(bed_name)
        with open(bed_path, "w") as file:
            # Write the BED data to the file
            for row in zip(chr_fusion, start, end, fusion_name, score_series, strand):
                file.write("\t".join(str(value) for value in row) + "\n")
            

    count_fusion_for_gene(combined_df_all)
    count_fusion_for_sample(combined_df_all)
    plot_conf_all_samples(combined_df_all)
    count_fusion_for_gene_top(combined_df_all)
    count_fusion_for_top_samples(combined_df_all)


if __name__ == "__main__":
    main()
