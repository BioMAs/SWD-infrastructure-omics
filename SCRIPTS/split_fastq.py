from __future__ import print_function
import os
import os.path
import sys
import gzip
import itertools
import json
import multiprocessing
import glob

## Parameters
r1_length = 16

# Open .fastq or .fastq.gz files for reading
def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "rU")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rt")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rt")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "rU")
    raise IOError("Unknown file: " + filename)

# Read config file
def read_config_file(config_path):
    config_file = open(config_path, "r")
    config = json.load(config_file)
    return config

# Write fastq file for each sample
def write_fastq(alignment_dir, sample_well, sample_name, sample_index, count, sample_name_seq_qual_list):
    if sample_well=="well_unknown":
        fastq_path = os.path.join(alignment_dir, ".".join(["unknown.fastq", str(count)]))
    else:
        fastq_path = os.path.join(alignment_dir, ".".join([sample_name,"fastq", str(count)]))
    if not os.path.exists(os.path.dirname(fastq_path)):
        try:
            os.makedirs(os.path.dirname(fastq_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with open(fastq_path, "w") as out:
    	for name, seq, qual in sample_name_seq_qual_list:
    	    print("\n".join(["@" + name, seq, "+", qual]), file=out)
    out.close

# Mask sequence by quality score
def mask(seq, qual, min_qual=10):
    return "".join((b if (ord(q) - 33) >= min_qual else "N") for b, q in itertools.zip_longest(seq, qual))

# Read fastq and split reads according to barcode
def split_reads(fastq_pair, count):
    index2sample = dict()
    for sample in config["samples"]:
    	sample["buf"] = list()
    	index2sample[sample["index"]] = sample

    index2sample["undefined"] = {"name":"sample_unknown","index":"index_unknown","well":"well_unknown"}
    index2sample["undefined"]["buf"] = list()
    readF = fastq_pair["read1"]
    readR = fastq_pair["read2"]
    with open_fastq_or_gz(readF) as r1_file, open_fastq_or_gz(readR) as r2_file:
    	r1_r2 = zip(r1_file, r2_file)
    	for header1, header2 in r1_r2:
            seq1, seq2 = next(r1_r2)
            plus1, plus2 = next(r1_r2)
            qual1, qual2 = next(r1_r2)

            read_name1, read_name2 = header1.split()[0][1:], header2.split()[0][1:]
            assert read_name1 == read_name2
            seq2, qual2 = seq2.rstrip(), qual2.rstrip()
            index, umi, seq, qual = mask(seq1[0:6], qual1[0:6], min_qual=10), mask(seq1[6:r1_length], qual1[6:r1_length]), seq2, qual2
            barcode = index+umi        
            barcoded_name = ":".join([read_name2, barcode])

            if index not in index2sample:
                index = "undefined"
    		
            index2sample[index]["buf"].append((barcoded_name, seq, qual))
    
    for index, sample in index2sample.items():
        if sample["well"]=="well_unknown":
            write_fastq(main_path, sample["well"], sample["name"], sample["index"], count, sample["buf"])
        else:
            write_fastq(os.path.join(main_path, sample["project"], config["fastq_folder"]), sample["well"], sample["name"], sample["index"], count, sample["buf"])

# Concatenate mutliple fastq files corresponding to sample
def concatenate_fastq():
    index2fastqs = dict()
    PROJECTS = {s["project"] for s in config["samples"]}
    for project in PROJECTS:
        #files = os.listdir(os.path.join(main_path, project, config["fastq_folder"]))
        files = glob.glob(os.path.join(main_path, project, config["fastq_folder"])+"/*.fastq.*")
        for f in files:
    	    filename, file_extension = os.path.splitext(f)
    	    if os.path.join(main_path, project, config["fastq_folder"], filename) not in index2fastqs:
    	        index2fastqs[os.path.join(main_path, project, config["fastq_folder"], filename)] = list()

    	    index2fastqs[os.path.join(main_path, project, config["fastq_folder"], filename)].append(os.path.join(main_path, project, config["fastq_folder"], f))

    #files = os.listdir(main_path)
    files = glob.glob(main_path+"/*.fastq.*")
    if files:
        for f in files:
    	    filename, file_extension = os.path.splitext(f)
    	    if os.path.join(main_path, filename) not in index2fastqs:
    	        index2fastqs[os.path.join(main_path, filename)] = list()

    	    index2fastqs[os.path.join(main_path, filename)].append(os.path.join(main_path, f))
        
    jobs = list()
    for f in index2fastqs:
    	jobs.append(multiprocessing.Process(target=concat_and_remove, args=(index2fastqs[f], f)))

    for j in jobs:
    	j.start()

    for j in jobs:
    	j.join()

# Remove multiple fastq files of same sample after concatenation
def concat_and_remove(flist, fout):
    print("Writing "+os.path.abspath(fout))
    with open(fout, "w") as outfile:
        for f in flist:
            with open(f, "r") as infile:
                for line in infile:
                    outfile.write(line)
                os.remove(f)

## Main
if len(sys.argv) != 2:
    print("Usage: python " + sys.argv[0] + " <config file>",
            file=sys.stderr)
    sys.exit()

config_file = sys.argv[1]
config = read_config_file(config_file)
main_path = config["maindir"]

count = 0
jobs = list()
for fastq_pair in config["fastq_pairs"]:
    count += 1
    jobs.append(multiprocessing.Process(target=split_reads, args=(fastq_pair, count)))

for j in jobs:
    j.start()

for j in jobs:
    j.join()

concatenate_fastq()
