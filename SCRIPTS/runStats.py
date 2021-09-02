import sys
import json
import os
import gzip

# Read config file
def read_config_file(config_path):
    with open(config_path, "r") as f:
    	config = json.load(f)
    return config

# Get number of reads in unknown.fastq
def count_unassigned_reads():
    unknown_fastq_path = os.path.join(config["maindir"], "unknown.fastq.gz")
    with gzip.open(unknown_fastq_path) as f:
        count = 0
        for _ in f:
            count += 1
    return (count+1)//4

def try_int(x):
    try:
        return int(x)
    except ValueError:
        return x

def getLogStats(exp):
    d = dict()
    d["WHOLERUN"] = dict()
    d["WHOLERUN"]["Total"] = 0
    d["WHOLERUN"]["Assigned"] = 0
    d["WHOLERUN"]["Aligned"] = 0
    d["WHOLERUN"]["Spike_Total"] = 0
    d["WHOLERUN"]["Spike_UMI"] = 0
    d["WHOLERUN"]["PolyA_Total"] = 0
    d["WHOLERUN"]["Mito_Total"] = 0
    d["WHOLERUN"]["Mito_UMI"] = 0
    d["WHOLERUN"]["Refseq_Total"] = 0
    d["WHOLERUN"]["Refseq_UMI"] = 0
    d["WHOLERUN"]["Unknown_Total"] = 0
    d["WHOLERUN"]["Unknown_UMI"] = 0

    for project in PROJECTS:
        with open(os.path.join(config["maindir"],project,config["expression_folder"],project+"."+exp+".log.dat")) as f:
            line = f.read().split('\n')[1]
            run_name,Assigned,Aligned,Spike_Total,Spike_UMI,PolyA_Total,Mito_Total,Mito_UMI,Refseq_Total,Refseq_UMI,Unknown_Total,Unknown_UMI = [try_int(x) for x in line.split("\t")]
            d[run_name] = dict()
            d[run_name]["Assigned"] = Assigned
            d[run_name]["Aligned"] = Aligned
            d[run_name]["Spike_Total"] = Spike_Total
            d[run_name]["Spike_UMI"] = Spike_UMI
            d[run_name]["PolyA_Total"] = PolyA_Total
            d[run_name]["Mito_Total"] = Mito_Total
            d[run_name]["Mito_UMI"] = Mito_UMI
            d[run_name]["Refseq_Total"] = Refseq_Total
            d[run_name]["Refseq_UMI"] = Refseq_UMI
            d[run_name]["Unknown_Total"] = Unknown_Total
            d[run_name]["Unknown_UMI"] = Unknown_UMI
            
            d["WHOLERUN"]["Total"] += Assigned
            d["WHOLERUN"]["Assigned"] += Assigned
            d["WHOLERUN"]["Aligned"] += Aligned
            d["WHOLERUN"]["Spike_Total"] += Spike_Total
            d["WHOLERUN"]["Spike_UMI"] += Spike_UMI
            d["WHOLERUN"]["PolyA_Total"] += PolyA_Total
            d["WHOLERUN"]["Mito_Total"] += Mito_Total
            d["WHOLERUN"]["Mito_UMI"] += Mito_UMI
            d["WHOLERUN"]["Refseq_Total"] += Refseq_Total
            d["WHOLERUN"]["Refseq_UMI"] += Refseq_UMI
            d["WHOLERUN"]["Unknown_Total"] += Unknown_Total
            d["WHOLERUN"]["Unknown_UMI"] += Unknown_UMI

    d["WHOLERUN"]["Total"] += UNKNOWNREADS
    return d

def writeLogStats(d, exp):
    with open(os.path.join(outdir, "run."+exp+".log.dat"), "w") as f:
        print("\t".join(["Run_Name",
                         "Total",
                         "Assigned",
                         "Aligned",
                         "Spike_Total",
                         "Spike_UMI",
                         "PolyA_Total",
                         "Mito_Total",
                         "Mito_UMI",
                         "Refseq_Total",
                         "Refseq_UMI",
                         "Unknown_Total",
                         "Unknown_UMI"]),
              file=f)
        for project in PROJECTS:
            print("\t".join([project] +
                         [str(x) for x in [ d[project]["Assigned"],
                                            d[project]["Assigned"],
                                            d[project]["Aligned"],
                                            d[project]["Spike_Total"],
                                            d[project]["Spike_UMI"],
                                            d[project]["PolyA_Total"],
                                            d[project]["Mito_Total"],
                                            d[project]["Mito_UMI"],
                                            d[project]["Refseq_Total"],
                                            d[project]["Refseq_UMI"],
                                            d[project]["Unknown_Total"],
                                            d[project]["Unknown_UMI"]]]),
                  file=f)

        print("\t".join(["WHOLERUN"] +
                         [str(x) for x in [ d["WHOLERUN"]["Total"],
                                            d["WHOLERUN"]["Assigned"],
                                            d["WHOLERUN"]["Aligned"],
                                            d["WHOLERUN"]["Spike_Total"],
                                            d["WHOLERUN"]["Spike_UMI"],
                                            d["WHOLERUN"]["PolyA_Total"],
                                            d["WHOLERUN"]["Mito_Total"],
                                            d["WHOLERUN"]["Mito_UMI"],
                                            d["WHOLERUN"]["Refseq_Total"],
                                            d["WHOLERUN"]["Refseq_UMI"],
                                            d["WHOLERUN"]["Unknown_Total"],
                                            d["WHOLERUN"]["Unknown_UMI"]]]),
                  file=f)

def getWellSummaryStats(exp):
    d = dict()
    for project in PROJECTS:
        with open(os.path.join(config["maindir"],project,config["expression_folder"],project+"."+exp+".well_summary.dat")) as f:
            samples = f.readline().strip('\n').split("\t")
            idx2sample = {idx: sample for idx, sample in enumerate(samples)}
            #print(idx2sample)
            for line in f:
                ls = line.strip('\n').split("\t")
                for i in range(1,len(ls)):
                    if(not idx2sample[i] in d):
                        d[idx2sample[i]] = dict()
                    d[idx2sample[i]][ls[0]] = ls[i]

    return d

def writeWellSummaryDat(d, exp):
    with open(os.path.join(outdir, "run."+exp+".well_summary.dat"), "w") as f:
        samples = sorted(d.keys())
        print("\t".join([""] + samples), file=f)
        print("\t".join(["Assigned"] +[d[sample]["Assigned"] for sample in samples]), file=f)
        print("\t".join(["Aligned"] +[d[sample]["Aligned"] for sample in samples]), file=f)
        print("\t".join(["Refseq_Total"] +[d[sample]["Refseq_Total"] for sample in samples]), file=f)
        print("\t".join(["Refseq_UMI"] +[d[sample]["Refseq_UMI"] for sample in samples]), file=f)
        print("\t".join(["Mito_Total"] +[d[sample]["Mito_Total"] for sample in samples]), file=f)
        print("\t".join(["Mito_UMI"] +[d[sample]["Mito_UMI"] for sample in samples]), file=f)
        print("\t".join(["Genes_Detected"] +[d[sample]["Genes_Detected"] for sample in samples]), file=f)
        print("\t".join(["Spike_Total"] +[d[sample]["Spike_Total"] for sample in samples]), file=f)
        print("\t".join(["Spike_UMI"] +[d[sample]["Spike_UMI"] for sample in samples]), file=f)

## Main
if len(sys.argv) != 3:
    print("Usage: python " + sys.argv[0] + " <config file> <outdir>",
            file=sys.stderr)
    sys.exit()

config_file = sys.argv[1]
outdir = sys.argv[2]
config = read_config_file(config_file)
PROJECTS = {s["project"] for s in config["samples"]}

UNKNOWNREADS = count_unassigned_reads()

stats = getLogStats("all")
writeLogStats(stats, "all")

stats = getLogStats("unq")
writeLogStats(stats, "unq")

ws = getWellSummaryStats("all")
writeWellSummaryDat(ws, "all")

ws = getWellSummaryStats("unq")
writeWellSummaryDat(ws, "unq")
