from __future__ import print_function
import os
import sys
import re
import numpy as np
import json
import multiprocessing
import pysam

## Functions
def get_mismatches(seq1,seq2):
    mismatches = 0    
    for i in range(0,len(seq1)):
        if seq1[i] != seq2[i] or seq1[i] == "N" or seq2[i] == "N":
            mismatches += 1
    return mismatches

def getSym2refForProject(project):
	for s in config["samples"]:
		if(s["project"]==project):
			return config["ref_folder"]+"/"+s["species"]+"/"+s["species"]+"_sym2ref.dat"

# Read config file
def read_config_file(config_path):
    with open(config_path, "r") as f:
    	config = json.load(f)
    return config

# Count
def count_reads(multimappings):
    if multimappings:
        multi_file_suffix = ".all"
    else:
        multi_file_suffix = ".unq"

    def count_reads_in_samfile(samfile, out_q):
        # Initialize counters
        count_assigned = 0
        count_assigned_aligned = 0

        count_spike_total = np.zeros(len(ercc_list), np.int)
        count_spike_umi = np.zeros(len(ercc_list), np.int)

        count_total = np.zeros(len(gene_list), np.int)
        count_umi = np.zeros(len(gene_list), np.int)

        count_assigned_mito_reads = 0
        count_assigned_mito_umi = 0

        count_assigned_unknown_reads = 0
        count_assigned_unknown_umi = 0

        count_assigned_polyA = 0

        seen_umi = set()
        unknown = set()

        # Open file with pysam
        samfile = pysam.AlignmentFile(samfile)
        for read in samfile.fetch(until_eof=True):
            read_id = read.query_name
            barcode = read_id.split(":")[-1]
            # If umi contains 'N', skip the read
            if (re.search(amb_umi_pattern, barcode[7:21])):
                continue
            # Add to count assigned (if the read is in this file, it's already assigned)
            count_assigned += 1
            # If read does not map on anything, skip
            if read.reference_id == -1:
                continue
            aligned_id = read.reference_name
            
            seq = read.query_sequence
            edit_dist = read.get_tag("NM")
            best_hits = read.get_tag("X0")
            # If read has too many mismatchs with ref, or maps at too many positions, or contains polyA pattern, skip
            if (edit_dist > max_edit_dist) or (best_hits > max_best) or (re.search(polyA_tail_pattern, seq)):
                if (re.search(polyA_tail_pattern, seq)):
                    count_assigned_polyA += 1
                continue
            best_hits_list = list()
            try:
                read.get_tag("XA")
                best_hits_loc = read.get_tag("XA")
                split_loc = best_hits_loc.split(";")
                for loc in split_loc:
                    new_loc = loc.split(",")[0]
                    if new_loc:
                        best_hits_list.append(new_loc)
            except:
                pass
            count_assigned_aligned += 1
            # If read maps on spike
            if aligned_id.startswith("ERCC"):
                if (not multimappings) and best_hits_list:
                    continue
                umi = "_".join([barcode[7:21], aligned_id])
                count_spike_total[ercc_to_idx[aligned_id]] += 1
                if umi not in seen_umi:
                    count_spike_umi[ercc_to_idx[aligned_id]] += 1
                    seen_umi.add(umi)
            # If read maps on chromosome M
            elif aligned_id.startswith("chrM"):
                if (not multimappings) and best_hits_list:
                    continue
                umi = "_".join([barcode[7:21], aligned_id])
                count_assigned_mito_reads += 1
                if umi not in seen_umi:
                    count_assigned_mito_umi += 1
                    seen_umi.add(umi)
            # If read maps on a gene
            elif aligned_id in refseq_to_gene:
                gene = refseq_to_gene[aligned_id]
                multi_gene_hits = False
                if best_hits_list:
                    # Iterate over best hits to ckeck if they all belong to the same gene
                    for loc in best_hits_list:
                        try: 
                            refseq_to_gene[loc]
                            if refseq_to_gene[loc] != gene:
                                multi_gene_hits = True
                                break
                        except:
                            multi_gene_hits = True
                if (not multimappings) and multi_gene_hits:
                    continue                
                umi = "_".join([barcode[7:21], gene])
                count_total[gene_to_idx[gene]] += 1
                if umi not in seen_umi:
                    count_umi[gene_to_idx[gene]] += 1
                    seen_umi.add(umi)
            # If read does not map
            else:
                if (not multimappings) and best_hits_list:
                    continue
                unknown.add(aligned_id)
                umi = "_".join([barcode[7:21], aligned_id])
                count_assigned_unknown_reads += 1
                if umi not in seen_umi:
                    count_assigned_unknown_umi += 1
                    seen_umi.add(umi)
            
        # return a dict for sample with all counts
        name = os.path.splitext(os.path.basename(samfile.filename))[0].decode("utf-8")
        outdict = {}
        outdict[name] = {}
        outdict[name]["count_assigned"] = count_assigned
        outdict[name]["count_assigned_aligned"] = count_assigned_aligned

        outdict[name]["count_spike_total"] = count_spike_total
        outdict[name]["count_spike_umi"] = count_spike_umi

        outdict[name]["count_total"] = count_total
        outdict[name]["count_umi"] = count_umi

        outdict[name]["count_assigned_mito_reads"] = count_assigned_mito_reads
        outdict[name]["count_assigned_mito_umi"] = count_assigned_mito_umi

        outdict[name]["count_assigned_polyA"] = count_assigned_polyA

        outdict[name]["count_assigned_unknown_reads"] = count_assigned_unknown_reads
        outdict[name]["count_assigned_unknown_umi"] = count_assigned_unknown_umi

        outdict[name]["seen_umi"] = seen_umi
        outdict[name]["unknown"] = unknown

        out_q.put(outdict)

    out_q = multiprocessing.Queue()
    procs = []
    nprocs = 0
    for aligned_filename in [f for f in os.listdir(aligned_dir)
                            if f.endswith(".bam")]:
        nprocs += 1
        samfile = os.path.join(aligned_dir, aligned_filename)
        p = multiprocessing.Process(target=count_reads_in_samfile, args=(samfile, out_q))
        procs.append(p)
        p.start()
   
    results = {}
    for i in range(nprocs):
        results.update(out_q.get())

    for p in procs:
        p.join()

    return results

# Write summary of sample contents
def write_reports(res, multi_file_suffix):
    count_assigned_reads = 0
    count_assigned_aligned_reads = 0
    count_spike_total = 0
    count_spike_umi = 0
    count_assigned_polyA = 0
    count_assigned_mito_reads = 0
    count_assigned_mito_umi = 0
    count_total = 0
    count_umi = 0
    count_assigned_unknown_reads = 0
    count_assigned_unknown_umi = 0
    seen_umi = set()
    unknown = set()
    
    for name, sample in res.items():
        count_assigned_reads += sample["count_assigned"]
        count_assigned_aligned_reads += sample["count_assigned_aligned"]
        count_spike_total += sample["count_spike_total"].sum()
        count_spike_umi += sample["count_spike_umi"].sum()
        count_assigned_polyA += sample["count_assigned_polyA"]
        count_assigned_mito_reads += sample["count_assigned_mito_reads"]
        count_assigned_mito_umi += sample["count_assigned_mito_umi"]
        count_total += sample["count_total"].sum()
        count_umi += sample["count_umi"].sum()
        count_assigned_unknown_reads += sample["count_assigned_unknown_reads"]
        count_assigned_unknown_umi += sample["count_assigned_unknown_umi"]
        seen_umi.update(sample["seen_umi"])
        unknown.update(sample["unknown"])

    with open(os.path.join(dge_dir, project + multi_file_suffix + ".log.dat"), "w") as f:
        print("\t".join(["Run_Name",
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
        print("\t".join([project] +
                         [str(x) for x in [count_assigned_reads,
                                            count_assigned_aligned_reads,
                                            count_spike_total,
                                            count_spike_umi,
                                            count_assigned_polyA,
                                            count_assigned_mito_reads,
                                            count_assigned_mito_umi,
                                            count_total,
                                            count_umi,
                                            count_assigned_unknown_reads,
                                            count_assigned_unknown_umi]]),
                file=f)

    # Write list of RefSeq IDs that were not mapped to genes by UCSC
    with open(os.path.join(dge_dir, project + multi_file_suffix + ".unknown_list"), "w") as f:
        for name in unknown:
            print(name, file=f)

    # Write gene x well total alignment counts
    with open(os.path.join(dge_dir, project + multi_file_suffix + ".refseq.total.dat"), "w") as f:
        print("\t".join([""] + well_list), file=f)
        for gene in gene_list:
            print("\t".join([gene]+[str(res[name]["count_total"][gene_to_idx[gene]]) for name in well_list]), file=f)

    # Write gene x well UMI counts
    with open(os.path.join(dge_dir, project + multi_file_suffix + ".refseq.umi.dat"), "w") as f:
        print("\t".join([""] + well_list), file=f)
        for gene in gene_list:
            print("\t".join([gene]+[str(res[name]["count_umi"][gene_to_idx[gene]]) for name in well_list]), file=f)

    # Write ERCC x well total alignment counts
    with open(os.path.join(dge_dir, project + multi_file_suffix + ".spike.total.dat"), "w") as f:
        print("\t".join([""] + well_list), file=f)
        for ercc in ercc_list:
            print("\t".join([ercc]+[str(res[name]["count_spike_total"][ercc_to_idx[ercc]]) for name in well_list]), file=f)

    # Write ERCC x well UMI counts
    with open(os.path.join(dge_dir, project + multi_file_suffix + ".spike.umi.dat"), "w") as f:
        print("\t".join([""] + well_list), file=f)
        for ercc in ercc_list:
            print("\t".join([ercc]+[str(res[name]["count_spike_umi"][ercc_to_idx[ercc]]) for name in well_list]), file=f)

    # Write by-well total and UMI counts
    well_summary_file = os.path.join(dge_dir, project + multi_file_suffix + ".well_summary.dat")
    with open(well_summary_file, "w") as f:
        print("\t".join([""] + well_list), file=f)
        print("\t".join(["Assigned"] +[str(res[name]["count_assigned"]) for name in well_list]), file=f)
        print("\t".join(["Aligned"] +[str(res[name]["count_assigned_aligned"]) for name in well_list]), file=f)
        print("\t".join(["Refseq_Total"] +[str(res[name]["count_total"].sum()) for name in well_list]), file=f)
        print("\t".join(["Refseq_UMI"] +[str(res[name]["count_umi"].sum()) for name in well_list]), file=f)
        print("\t".join(["Mito_Total"] +[str(res[name]["count_assigned_mito_reads"]) for name in well_list]), file=f)
        print("\t".join(["Mito_UMI"] +[str(res[name]["count_assigned_mito_umi"]) for name in well_list]), file=f)
        print("\t".join(["Genes_Detected"] +[str((res[name]["count_umi"] != 0).sum(axis=0)) for name in well_list]), file=f)
        print("\t".join(["Spike_Total"] +[str(res[name]["count_spike_total"].sum(axis=0)) for name in well_list]), file=f)
        print("\t".join(["Spike_UMI"] +[str(res[name]["count_spike_umi"].sum(axis=0)) for name in well_list]), file=f)

## Main
if len(sys.argv) != 3:
    print("Usage: python " + sys.argv[0] + " <config file> <project_name>",
            file=sys.stderr)
    sys.exit()

config_file = sys.argv[1]
project = sys.argv[2]
config = read_config_file(config_file)
sym2ref = getSym2refForProject(project)
ercc_fasta=os.path.abspath(os.path.join(os.path.dirname(__file__), "ERCC92_polyAstrip.fa"))
aligned_dir=os.path.join(config["maindir"],project,config["align_folder"])
dge_dir=os.path.join(config["maindir"],project,config["expression_folder"])
if not os.path.exists(dge_dir):
    os.makedirs(dge_dir)

# Read gene symbol->RefSeq ID mapping
gene_list = list()
refseq_to_gene = dict()
max_best = 0
with open(sym2ref, "rU") as f:
    for line in f:
        items = line.split()
        gene, refseqs = items[0], items[1].split(",")
        if len(refseqs) > max_best:
            max_best=len(refseqs)
        gene_list.append(gene)
        refseq_to_gene.update([(rs, gene) for rs in refseqs])
gene_to_idx = {gene: idx for idx, gene in enumerate(gene_list)}

# Read ERCC control IDs
ercc_list = list()
with open(ercc_fasta, "rU") as f:
    for line in f:
        if line[0] == ">":
            items = line.split()
            ercc_list.append(items[0][1:])
ercc_to_idx = {ercc: idx for idx, ercc in enumerate(ercc_list)}

# Read barcode->well mapping
well_list = [s["name"] for s in config["samples"] if s["project"]==project]

polyA_tail_pattern = re.compile(r"A{10,}$")
amb_umi_pattern = re.compile(r"[A-Z]*N[A-Z]*")
max_edit_dist = 3

res_multi = count_reads(True)
res_uniq = count_reads(False)

write_reports(res_multi, ".all")
write_reports(res_uniq, ".unq")
    
