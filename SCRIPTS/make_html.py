import json
import argparse
import os
import glob
import itertools
import pandas as pd
import sys
import time
from jinja2 import Environment, FileSystemLoader


########
# Main #
########

# Parser
argParser = argparse.ArgumentParser()
argParser.add_argument("-c", "--configfile", required=True, help="Jsonfile's path", dest="configFile")
argParser.add_argument("-t", "--templateFile", required=True, help="HTML Template's path", dest="templateFile")
argParser.add_argument("-p", "--projectName", required=True, help="Name of the project used to select samples in config file", dest="projectName")
argParser.add_argument("-o", "--output", required=True, help="Path of the output file", dest="output")

def getGenomeForProject(project):
	for s in config["samples"]:
		if(s["project"]==project):
			return s["species"]

# Recovery of arguments
args = argParser.parse_args()

# Load config file
config = json.load(open(args.configFile, "r"))
project = args.projectName

# Recovery json files' parts
multiqcFolder = os.path.join(config["maindir"],project,config["multiqc_folder"])
expressionFolder = os.path.join(config["maindir"],project,config["expression_folder"])

# Dict() creation for feeding source of the template
templateVars = dict()

# Load configuration file in template
templateVars["config"] = config

########################
# Project informations #
########################

templateVars["project"] = dict();
templateVars["project"]["projectName"] = project
templateVars["project"]["date"] = time.strftime('%B %d, %Y',time.localtime())
templateVars["project"]["genome"] = getGenomeForProject(project)

#####################
# FASTQC HighCharts #
#####################

# Opening of multiqc_general_stats.txt file
csvFile = pd.read_csv(os.path.join(multiqcFolder,"multiqc_data","multiqc_general_stats.txt"), header=0, sep="\t")

# Recovery of the samples' names
sampleNames = list(csvFile['Sample'])
	
# Recovery of total number sequences in every sample
totseq = list(csvFile['FastQC_mqc-generalstats-fastqc-total_sequences'])

# Feeding of the templateVars dict() with fastqc infos
fastqcChart = dict()
fastqcChart["samples"] = sampleNames
fastqcChart["totseq"] = totseq
templateVars["fastqcChart"] = fastqcChart


############################
# Barplot unq.well_summary #
############################

# Opening of DGE.unq.well_summary.dat file
csvFile2 = pd.read_csv(os.path.join(expressionFolder,project+".unq.well_summary.dat"), header=None, sep="\t").T

# Recovery of samples' names
Samples = csvFile2[csvFile2.columns[0]]
Samples = list(Samples)
del (Samples[0])

# Recovery of Assigned reads
Assigned = csvFile2[csvFile2.columns[1]]
Assigned = list(Assigned)
del (Assigned[0])
Assigned = ('[' + ', '.join(Assigned) + ']')

# Recovery of Aligned reads
Aligned = csvFile2[csvFile2.columns[2]]
Aligned = list(Aligned)
del (Aligned[0])
Aligned = ('[' + ', '.join(Aligned) + ']')

# Recovery of Refseq_Total
Refseq_Total = csvFile2[csvFile2.columns[3]]
Refseq_Total = list(Refseq_Total)
del (Refseq_Total[0])
Refseq_Total = ('[' + ', '.join(Refseq_Total) + ']')

# Recovery of Refseq_UMI
Refseq_UMI = csvFile2[csvFile2.columns[4]]
Refseq_UMI = list(Refseq_UMI)
del (Refseq_UMI[0])
Refseq_UMI = ('[' + ', '.join(Refseq_UMI) + ']')

# Recovery of Genes_Detected
Genes_Detected = csvFile2[csvFile2.columns[7]]
Genes_Detected = list(Genes_Detected)
del (Genes_Detected[0])
Genes_Detected = ('[' + ', '.join(Genes_Detected) + ']')

# Feeding of the templateVars dict() with Barplot infos
barplotChart = dict()
barplotChart["Samples"] = Samples
barplotChart["Assigned"] = Assigned
barplotChart["Aligned"] = Aligned
barplotChart["Refseq_Total"] = Refseq_Total
barplotChart["Refseq_UMI"] = Refseq_UMI
barplotChart["Genes_Detected"] = Genes_Detected
templateVars["barplotChart"] = barplotChart

###############
# Comparisons #
###############

if(project in config["comparisons"]):
	templateVars["project"]["minRep"] = config["comparisons"][project]["minRep"]
	templateVars["project"]["minGenes"] = config["comparisons"][project]["minGenes"]
	templateVars["project"]["minReads"] = config["comparisons"][project]["minReads"]
	sList = list()
	with open(os.path.join(config["maindir"],project,config["de_folder"],"samplesToRemove.txt")) as f:
		next(f)
		for line in f:
			s = dict()
			ls = line.strip('\n').split("\t")
			s["sample"] = ls[0]
			s["totalGenEx"] = int(ls[1])
			s["totalCounts"] = int(ls[2])
			sList.append(s)

	templateVars["project"]["filteredSamples"] = sList
	if(sList):
		templateVars["project"]["filteredExists"] = True

# Selection of versus conditions
conditionVS = list()
if(project in config["comparisons"]):
	conditionVS.extend(config["comparisons"][project]["comps"])
templateVars["conditionVS"] = conditionVS

templateVars["deseqres"] = dict()
templateVars["deseqres"]["nbdeg"] = list()
for c in conditionVS:
	# Create list for comparison
	templateVars["deseqres"][c["condition1"]+"__vs__"+c["condition2"]] = dict()
	templateVars["deseqres"][c["condition1"]+"__vs__"+c["condition2"]]["genes"] = list()
	# Initialize counter for number of differentially expressed genes
	count = 0
	# Create object for number of DEG table
	nbGenesObject = dict()
	nbGenesObject["cond1"] = c["condition1"]
	nbGenesObject["cond2"] = c["condition2"]
	# Opening of deseqRes file
	with open(os.path.join(config["maindir"],project,config["de_folder"],c["condition1"]+"__vs__"+c["condition2"],"DEseqResFiltered.tsv"), "r") as f:
		next(f)
		for line in f:
			count += 1
			ls = line.strip('\n').split("\t")
			g = dict()
			g["Gene"] = ls[0]
			g["baseMean"] = str("%.2f" % float(ls[1]))
			g["log2FC"] = str("%.4f" % float(ls[2]))
			g["lfcSE"] = str("%.4f" % float(ls[3]))
			g["stat"] = str("%.4f" % float(ls[4]))
			g["pvalue"] = ls[5]
			g["padj"] = ls[6]
			g["meanInComp"] = str("%.4f" % float(ls[7]))
			templateVars["deseqres"][c["condition1"]+"__vs__"+c["condition2"]]["genes"].append(g)
		
	# Insert counts of DEG in a variable in order to print ma-plot, volcano-plot and deseq table in report if necessary
	templateVars["deseqres"][c["condition1"]+"__vs__"+c["condition2"]]["counts"] = count
	# Insert counts in number of DEG table
	nbGenesObject["count"] = count
	templateVars["deseqres"]["nbdeg"].append(nbGenesObject)

########################
# FASTQC files listing #
########################

# Recovery of samples for project
samplelist = list()
for s in config["samples"]:
	if(s["project"]==project):
		samplelist.append(s)
templateVars["samples"] = samplelist


#################
# Tools listing #
#################

# Recovery of versions.sh results
pack = os.popen("bash "+os.path.join(sys.path[0],"versions.sh"),'r')
packages = pack.read()
pack.close()

packages = packages.replace("\n"," ; ")
templateVars["packages"] = packages


############
# TEMPLATE #
############ 

# Write in the template
template = args.templateFile
output = Environment(loader=FileSystemLoader(os.path.dirname(template))).get_template(os.path.basename(template)).render(templateVars)
with open(args.output, "wb") as f:
	f.write(output.encode("utf-8"))
