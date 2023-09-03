from pybedtools import BedTool
# Create a BedTool object
#bed_file = "/home/ls/parshas/all_dogs_withnames.bed"
#bed_file = "/home/ls/parshas/vscode/5464_processing/elife4565_hg38.bed"
bed_file = "/home/ls/parshas/vscode/5464_processing/dogs_only_in_RTdis0.bed"
bedtool = BedTool(bed_file)

# Add the -slop option with a 80 kb value
#bedtool = bedtool.slop(b=80000)

# Generate IGV file
#igv_file = "screanshoots_igv.igv"
#igv_file = "screanshoots_igv_elife.igv"
igv_file = "screenshoots_igv_dogs_only_in_RTdis0"
bedtool.igv(output=igv_file, slop=80000)