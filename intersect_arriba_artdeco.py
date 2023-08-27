import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import glob
import numpy as np
import os



def plot():
    dir = "/private/projects/kidney_rtt/intersect_all_samples/summary_intersect.txt"

    # Define the data dictionary
    data = {}

    # Read the data from the text file
    with open(dir, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if lines[i].startswith('sample:'):
                sample = lines[i].split(':')[1].strip()
                percentage = float(lines[i+3].split(':')[1].strip().strip('%'))
                data[sample] = percentage

# Sort the data dictionary based on percentages in descending order
    sorted_data = dict(sorted(data.items(), key=lambda x: x[1], reverse=True))

    # Extract the sample names and percentages
    samples = list(sorted_data.keys())
    percentages = list(sorted_data.values())
    fig, ax = plt.subplots(figsize=(10, 12))
    ax.bar(samples, percentages)
    ax.set_xlabel('Samples')
    ax.set_ylabel('Percentage')
    ax.set_title('For Each Sample Percent of DoGs Found in Arriba')
    ax.set_xticks([])
    # Add percentage symbol to y-axis labels
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())
    plt.savefig(

        '/private/projects/kidney_rtt/filters_arriba_plot/intersect.png')


def discarded_or_not_plot():
    data = {}
    dir = "/private/projects/kidney_rtt/intersect_all_samples/"
    for input_file in os.listdir(dir):
        file_path = os.path.join(dir, input_file)
        #check if file end with ".bed"
        if not input_file.endswith(".bed"):
            continue
        # read each line in the file
        file_name = input_file.split("_", 1)[1].split(".dogs")[0]
        with open(file_path, 'r') as input:
            input_lines = input.readlines()
        directory = "/private/projects/kidney_rtt/arriba_outputs/"
        # Search for fusion files and discarded files in the directory
        fusion_file = glob.glob(directory + file_name + ".fusions.tsv")[0]
        discarded_file = glob.glob(
            directory + file_name + ".fusions.discarded.tsv")[0]

        fusions_lines = []
        discarded_lines = []
        for line in input_lines:
            line = line.strip()
            chr = line.split("\t")[0]
            start =line.split("\t")[1]
            end = line.split("\t")[2]
            strand = line.split("\t")[5]
            
            # split the line by tab and 6 column
            found = line.split("\t")[6]
            if found != 0:

                with open(fusion_file, "r") as file:
                    for line2 in file:
                        if line2.startswith('#'):
                           continue
                        # take chr start end and strand from the line in the file
                        chr2 = line2.split("\t")[4].split(
                            ':')[0]
                        start2 = line2.split("\t")[4].split(
                            ':')[1]
                        end2 = line2.split("\t")[4].split(
                            ':')[1]
                        strand2 = line2.split("\t")[2].split(
                            '/')[1]
                        #swap if start is bigger than end
                        if start > end:
                            start, end = end, start
                        if (chr==chr2 and strand==strand2 and start2==start and end2==end):
                            fusions_lines.append(line)
                            break


                with open(discarded_file, "r") as file:
                    for line2 in file:
                        if line2.startswith('#'):
                           continue
                        # take chr start end and strand from the line in the file
                        chr2 = line2.split("\t")[4].split(
                            ':')[0]
                        start2 = line2.split("\t")[4].split(
                            ':')[1]
                        end2 = line2.split("\t")[4].split(
                            ':')[1]
                        strand2 = line2.split("\t")[2].split(
                            '/')[1]
                        #swap if start is bigger than end
                        if start > end:
                            start, end = end, start
                        if (chr==chr2 and strand==strand2 and start2==start and end2==end):
                            fusions_lines.append(line)
                            break
        # save in dicationay how mant fusions and discarded for each sample
        data[file_name] = [len(fusions_lines)/len(input_lines)
                           * 100, len(discarded_lines)/len(input_lines)*100]
    samples = list(data.keys())
    fusion_percentages = [x[0] for x in data.values()]
    discarded_percentages = [x[1] for x in data.values()]

    fig, ax = plt.subplots(figsize=(10, 12))
    fusion_color = 'blue'
    discarded_color = 'red'
    # Plot the data
# Plot the data
    bar_width = 0.4
    index = np.arange(len(samples))

    ax.bar(index, fusion_percentages, width=bar_width,
           color=fusion_color, label='Fusions')
    ax.bar(index, discarded_percentages, width=bar_width,
           bottom=fusion_percentages, color=discarded_color, label='Discarded')
    ax.set_xlabel('Samples')
    ax.set_ylabel('Percentage')
    ax.set_title('For Each Sample Percent of DoGs Found in Arriba')
    ax.set_xticks([])
    # Add percentage symbol to y-axis labels
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())
    fig.savefig(

        '/private/projects/kidney_rtt/filters_arriba_plot/intersect2.png')
    # read fusion_lines and divide into 3 groups by site1 site2 column (non-coding gene1/CDS gene1 only/CDS gene1+CDS gene2)
    # cds_cds = []
    # sit1_cds = []
    # non_coding = []
    # for line in fusions_lines:
    #     site1 = line.split("\t")[6]
    #     site2 = line.split("\t")[7]
    #     if site1 == "CDS" and site2 == "CDS":
    #         cds_cds.append(line)
    #     elif site1 == "CDS":
    #         sit1_cds.append(line)
    #     else:
    #         non_coding.append(line)

def types_or_not_plot():
    data = {}
    dir = "/private/projects/kidney_rtt/intersect_all_samples/"
    for input_file in os.listdir(dir):
        file_path = os.path.join(dir, input_file)
        #check if file end with ".bed"
        if not input_file.endswith(".bed"):
            continue
        # read each line in the file
        file_name = input_file.split("_", 1)[1].split(".dogs")[0]
        with open(file_path, 'r') as input:
            input_lines = input.readlines()
        directory = "/private/projects/kidney_rtt/arriba_outputs/"
        # Search for fusion files and discarded files in the directory
        fusion_file = glob.glob(directory + file_name + ".fusions.tsv")[0]
        discarded_file = glob.glob(
            directory + file_name + ".fusions.discarded.tsv")[0]

        cds_cds = []
        sit1_cds = []
        non_coding = []
        for line in input_lines:
            line = line.strip()
            chr = line.split("\t")[0]
            start =line.split("\t")[1]
            end = line.split("\t")[2]
            strand = line.split("\t")[5]
            
            # split the line by tab and 6 column
            found = line.split("\t")[6]
            if found != 0:

                with open(fusion_file, "r") as file:
                    for line2 in file:
                        if line2.startswith('#'):
                           continue
                        # take chr start end and strand from the line in the file
                        chr2 = line2.split("\t")[4].split(
                            ':')[0]
                        start2 = line2.split("\t")[4].split(
                            ':')[1]
                        end2 = line2.split("\t")[4].split(
                            ':')[1]
                        strand2 = line2.split("\t")[2].split(
                            '/')[1]
                        site1=line2.split("\t")[6]
                        site2=line2.split("\t")[7]

                        #swap if start is bigger than end
                        if start > end:
                            start, end = end, start
                        if (chr==chr2 and strand==strand2 and start2==start and end2==end):
                            if site1 == "CDS" and site2 == "CDS":
                                cds_cds.append(line)
                            elif site1 == "CDS":
                                sit1_cds.append(line)
                            else:
                                non_coding.append(line)
                                                    


                with open(discarded_file, "r") as file:
                    for line2 in file:
                        if line2.startswith('#'):
                           continue
                        # take chr start end and strand from the line in the file
                        chr2 = line2.split("\t")[4].split(
                            ':')[0]
                        start2 = line2.split("\t")[4].split(
                            ':')[1]
                        end2 = line2.split("\t")[4].split(
                            ':')[1]
                        strand2 = line2.split("\t")[2].split(
                            '/')[1]
                        site1=line2.split("\t")[6]
                        site2=line2.split("\t")[7]

                        #swap if start is bigger than end
                        if start > end:
                            start, end = end, start
                        if (chr==chr2 and strand==strand2 and start2==start and end2==end):
                            if site1 == "CDS" and site2 == "CDS":
                                cds_cds.append(line)
                            elif site1 == "CDS":
                                sit1_cds.append(line)
                            else:
                                non_coding.append(line)
        # save in dicationay how mant fusions and discarded for each sample
        data[file_name] = [len(cds_cds)/len(input_lines)
                           * 100, len(sit1_cds)/len(input_lines)*100, len(non_coding)/len(input_lines)*100]
    samples = list(data.keys())
    cds_cds_percentages = [x[0] for x in data.values()]
    sit1_cds_percentages = [x[1] for x in data.values()]
    non_coding_percentages = [x[2] for x in data.values()]

    fig, plt2 = plt.subplots(figsize=(10, 12))
    # Define colors for the sections
    cds_cds_color = 'blue'
    sit1_cds_color = 'orange'
    non_coding_color = 'green'

    # Plot the data
    bar_width = 0.4
    index = np.arange(len(samples))

    plt2.bar(index, cds_cds_percentages, width=bar_width, color=cds_cds_color, label='CDS-CDS')
    plt2.bar(index, sit1_cds_percentages, width=bar_width, bottom=cds_cds_percentages, color=sit1_cds_color, label='CDS-non_CDS')
    plt2.bar(index, non_coding_percentages, width=bar_width, bottom=np.add(cds_cds_percentages, sit1_cds_percentages), color=non_coding_color, label='non_CDS-non_CDS')
    plt2.set_xlabel('Samples')
    plt2.set_ylabel('Percentage')
    plt2.title('For Each Sample Percent of DoGs Found in Arriba divided by Coding Types for each Sample')
    ax.set_xticks([])
    # Add percentage symbol to y-axis labels
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mticker.PercentFormatter())
    fig.savefig(
        '/private/projects/kidney_rtt/filters_arriba_plot/intersect3.png')
if __name__ == "__main__":
    plot()
    types_or_not_plot()
