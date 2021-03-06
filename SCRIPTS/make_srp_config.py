from __future__ import print_function
import os
import json
import argparse
import csv
import re
import sys
import urllib.request
from collections import OrderedDict

## Main
parser = argparse.ArgumentParser()

parser.add_argument('-s', '--samplesheet', metavar='FILE', help='Tab delimited file with no header describing samples. Columns must be: "well index name project condition species". Only characters "A-Z","0-9","-" and "_" allowed. All columns are mandatory. (REQUIRED)', type=argparse.FileType('rt'), required=True, dest='samplesheet')
parser.add_argument('-w', '--workdir', metavar='DIR', help='Analysis working directory. Default: current directory', default='.', dest='workdir')
parser.add_argument('-i', '--illumina-dir', metavar='DIR', help='Directory containing the fastq input files generated by Illumina. The fastq files should be in paired-end mode. If your fastq files are not coming from an Illumina sequencer, please use option -f to specify a file listing the fastq input files. (REQUIRED if no "-f")', dest='illuminadir')
parser.add_argument('-f', '--fastq-file', metavar='FILE', help='File describing the fastq input file. This file should be tab delimited. First column: full path of Forward file; second column: full path of Reverse file. The fastq files should be in paired-end mode. If your fastq files were generated by an Illumina sequencer, you can use option "-i" to specify the directory containing the fastq input files. (REQUIRED if no "-i")', type=argparse.FileType('rt'), dest='fastqfile')
parser.add_argument('-r', '--reference-dir', metavar='DIR', help='Directory containing the reference files. It is recommended that you use this option if you have already used this pipeline and downloaded genome files.', dest='referencedir')
parser.add_argument('-c', '--conditions', metavar='FILE', help='Tab delimited file with no headers indicating which conditions to compare during differential expression analysis. Columns must be "project condition1 condition2". If not specified, only primary analysis will be performed', type=argparse.FileType('rt'), required=False, dest='conditions')
parser.add_argument('--minGenes', metavar='N', type=int, help='Minimum genes detected necessary for a sample to pass the filtering step in secondary analysis. (Default 5000)', default=5000, required=False, dest='minGenes')
parser.add_argument('--minReads', metavar='N', type=int, help='Minimum reads assigned necessary for a sample to pass the filtering step in secondary analysis. (Default 200000)', default=200000, required=False, dest='minReads')
parser.add_argument('--minLogFC', metavar='N', type=int, help='Minimum log Fold-Change threshold for differentially expessed gene. (Default 0.58 (1.5 FC))', default=0.58, required=False, dest='minLogFC')

args = parser.parse_args()
d = OrderedDict()

# Write working directory
workdir = os.path.abspath(args.workdir)
if not os.path.exists(workdir):
    try:
        os.makedirs(workdir)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise

d["maindir"] = workdir
d["fastq_folder"] = "FASTQ"
d["cutadapt_folder"] = "CUTADAPT"
d["fastqc_folder"] = "FASTQC"
d["multiqc_folder"] = "MULTIQC"
d["align_folder"] = "ALIGNMENT"
d["expression_folder"] = "EXPRESSION"
d["de_folder"] = "DE"
d["report_folder"] = "REPORT"

# Write reference directory
if(args.referencedir is None):
    referencedir = os.path.join(workdir, "REFERENCE")
else:
    referencedir = os.path.abspath(args.referencedir)

d["ref_folder"] = referencedir

# Write samples
d["samples"] = list()
alphanum = re.compile('^[\w-]+$')
alphanumStartLetter = re.compile('^[a-zA-Z][\w-]*$')
count = 1
try:
    for well,index,name,projet,condition,espece in csv.reader(args.samplesheet, delimiter='\t'):
        if(alphanum.match(well) and alphanum.match(name) and alphanum.match(index) and alphanum.match(projet) and alphanum.match(condition) and alphanum.match(espece)) :
            d["samples"].append(OrderedDict([("well",well),("name",name),("index",index),("project",projet),("condition",condition), ("species",espece)]))
        else:
            print("Invalid character or empty field in samplesheet at line "+str(count)+" :\""+well+"\t"+index+"\t"+name+"\t"+projet+"\t"+condition+"\t"+espece+"\"")
            sys.exit()
        if(not alphanumStartLetter.match(projet)):
            print("Invalid sequence for project name \""+projet+"\" in samplesheet at line "+str(count)+". Project name must start with a letter !")
            sys.exit()
        if(not alphanumStartLetter.match(condition)):
            print("Invalid sequence for condition name \""+condition+"\" in samplesheet at line "+str(count)+". Condition name must start with a letter !")
            sys.exit()
        count += 1

except:
    print("Invalid entry in samplesheet at line "+str(count))
    sys.exit()

# Verify if indexes are unique
indexes = list()
for sample in d["samples"]:
    if(sample["index"] in indexes):
        print("Invalid samplesheet. Index seen more than once: "+sample["index"])
        sys.exit()
    else:
        indexes.append(sample["index"])

# Verify if all samples in a project are from same species
projects = dict()
for sample in d["samples"]:
    if(not sample["project"] in projects):
        projects[sample["project"]] = sample["species"]
    if(sample["species"] != projects[sample["project"]]):
        print("Invalid samplesheet. More than one species specified for project: "+sample["project"])
        sys.exit()

# Verify that all sample identifiers are unique for a project
projects = dict()
for sample in d["samples"]:
    if(not sample["project"] in projects):
        projects[sample["project"]] = list()
    if(not sample["name"] in projects[sample["project"]]):
        projects[sample["project"]].append(sample["name"])
    else:
        print("Invalid samplesheet. "+sample["name"]+" found twice in project "+sample["project"])
        sys.exit() 

def calculateMinRep(project):
    cond2number = dict()
    for s in d["samples"]:
        if(s["project"]==project):
            if(not s["condition"] in cond2number):
                cond2number[s["condition"]] = 0
            cond2number[s["condition"]] += 1
    minRep = 999999999
    for cond in cond2number:
        if(cond2number[cond]<minRep):
            minRep = cond2number[cond]
    return minRep

def calculateMinRepForOneCondition(project):
    nbSamples = len([s for s in d["samples"] if s["project"]==project])
    return int(nbSamples*0.2)

def getGenomeForProject(project):
    for s in d["samples"]:
        if(s["project"]==project):
            return s["species"]

# Write conditions to analyse
d["comparisons"] = dict()
if (args.conditions is not None):
    for project,cond1,cond2 in csv.reader(args.conditions, delimiter='\t'):
        if(not project in d["comparisons"]):
            d["comparisons"][project] = OrderedDict()
            d["comparisons"][project]["species"] = getGenomeForProject(project)
            d["comparisons"][project]["minRep"] = calculateMinRep(project)
            d["comparisons"][project]["minGenes"] = args.minGenes
            d["comparisons"][project]["minReads"] = args.minReads
            d["comparisons"][project]["minLogFC"] = args.minLogFC
            d["comparisons"][project]["performComps"] = None
            d["comparisons"][project]["comps"] = list()
        if(cond1!=cond2):
            if(d["comparisons"][project]["performComps"]!=False):
                d["comparisons"][project]["comps"].append(OrderedDict([("condition1",cond1),("condition2",cond2)]))
                d["comparisons"][project]["performComps"] = True
            else:
                print("Invalid condition file. In project: \""+project+"\" both comparisons and first part secondary analysis defined.\nPlease choose if you want to perform comparisons (different conditions specified) or if you only want to perform first part of secondary analysis (identical conditions specified)")
                sys.exit()
        else:
            if(d["comparisons"][project]["performComps"]==None):
                d["comparisons"][project]["performComps"] = False
                d["comparisons"][project]["minRep"] = 1
                continue
            if(d["comparisons"][project]["performComps"]):
                print("Invalid condition file. In project: \""+project+"\" both comparisons and first part secondary analysis defined.\nPlease choose if you want to perform comparisons (different conditions specified) or if you only want to perform first part of secondary analysis (identical conditions specified)")
                sys.exit()
            


# Verify if all conditions specified for project exist in samplesheet
projects = dict()
for sample in d["samples"]:
    if(not sample["project"] in projects):
        projects[sample["project"]] = set()
    projects[sample["project"]].add(sample["condition"])
for project in d["comparisons"]:
    for cond in d["comparisons"][project]["comps"]:
        if(not project in projects):
            print("Invalid project found in conditions table. Project: \""+project+"\" not found in samplesheet")
            sys.exit()
        if(not cond["condition1"] in projects[project]):
            print("Invalid condition found in conditions table. Condition: \""+cond["condition1"]+"\" not found for project: \""+project+"\" in samplesheet")
            sys.exit()
        if(not cond["condition2"] in projects[project]):
            print("Invalid condition found in conditions table. Condition: \""+cond["condition2"]+"\" not found for project: \""+project+"\" in samplesheet")
            sys.exit()

# Write fastq input files
d["fastq_pairs"] = list()
if((args.illuminadir is None) and (args.fastqfile is None)):
    print("You must specify either an illumina directory containing fastq input files with option \"-i\" or a file describing your fastq input files with option \"-f\"")
    sys.exit()
elif((args.illuminadir is not None) and (args.fastqfile is not None)):
    print("You must specify EITHER an illumina directory containing fastq input files with option \"-i\" OR a file describing your fastq input files with option \"-f\", not both")
    sys.exit()
elif((args.illuminadir is None) and (args.fastqfile is not None)):
    count = 0
    try:
        for f,r in csv.reader(args.fastqfile, delimiter='\t'):
            count += 1
            if(os.path.isfile(f) and os.path.isfile(r)):
                 d["fastq_pairs"].append(OrderedDict({"read1":f, "read2":r}))
            elif(not os.path.isfile(f)):
                 print("File cannot be read: "+f)
                 sys.exit()
            elif(not os.path.isfile(r)):
                 print("File cannot be read: "+r)
                 sys.exit()
    except:
        print("Invalid entry in fastq file at line "+str(count+1))
        sys.exit()
elif((args.illuminadir is not None) and (args.fastqfile is None)):
    for fastqF in [f for f in os.listdir(args.illuminadir) if(re.match("^[\w-]+_R1_\d{3}.fastq(.gz)?", f))]:
        s = re.search("(^[\w-]+)_R1_(\d{3}.fastq(.gz)?)", fastqF)
        fastqFfile = os.path.join(args.illuminadir, fastqF)
        fastqRfile = os.path.join(args.illuminadir, s.group(1) + "_R2_" + s.group(2))
        if(os.path.isfile(fastqRfile)):
            d["fastq_pairs"].append(OrderedDict([("read1",os.path.abspath(fastqFfile)), ("read2",os.path.abspath(fastqRfile))]))

# Verify that genome can be downloaded if it does not exist
genomes = set()
for sample in d["samples"]:
    genomes.add(sample["species"])

for genome in genomes:
    try: 
        if ((not os.path.isfile(os.path.join(referencedir, genome, genome+"_chrM.fa.gz"))) and (not os.path.isfile(os.path.join(referencedir, genome, genome+"_chrM.fa")))):
            urllib.request.urlopen("http://hgdownload.cse.ucsc.edu/goldenPath/"+genome+"/chromosomes/chrM.fa.gz")
        if ((not os.path.isfile(os.path.join(referencedir, genome, genome+"_refGene.txt.gz"))) and (not os.path.isfile(os.path.join(referencedir, genome, genome+"_refGene.txt")))):
            urllib.request.urlopen("http://hgdownload.cse.ucsc.edu/goldenPath/"+genome+"/database/refGene.txt.gz")
        if ((not os.path.isfile(os.path.join(referencedir, genome, genome+"_refMrna.fa.gz"))) and (not os.path.isfile(os.path.join(referencedir, genome, genome+"_refMrna.fa")))):
            urllib.request.urlopen("http://hgdownload.cse.ucsc.edu/goldenPath/"+genome+"/bigZips/refMrna.fa.gz")
    except urllib.error.URLError as e:
        print("Cannot download data on UCSC for the assembly: "+genome+"\n"+e.reason+"\nYou may have to download the needed files manually (chrM.fa.gz, refGene.txt.gz, refMrna.fa.gz)!")
        sys.exit()
            
print(json.dumps(d , indent=4))
