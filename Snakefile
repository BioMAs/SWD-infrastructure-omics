import shutil
import os
import datetime
import json
import sys
import subprocess

configfile: config["conf"]
wfbasedir = workflow.basedir
erccFile = workflow.basedir+"/SCRIPTS/ERCC92_polyAstrip.fa"

MAINDIR = config["maindir"]
PROJECTS = {s["project"] for s in config["samples"]}
SECONDARY_ANALYSIS_PROJECTS = {c for c in config["comparisons"] if config["comparisons"][c]["performComps"]}
NO_COMPARISONS_PROJECTS = {c for c in config["comparisons"] if not config["comparisons"][c]["performComps"]}
DATE = str(datetime.date.today())

# final file outputs suffixes for primary analysis
finalSuffixes = ["log.dat","refseq.total.dat","refseq.umi.dat","spike.total.dat","spike.umi.dat","unknown_list","well_summary.dat"]

def getSplitFastqOutputs():
	files = list()
	for s in config["samples"]:
		files.append(os.path.join(MAINDIR,s["project"],config["fastq_folder"],s["name"]+".fastq"))
	files.append(MAINDIR+"/unknown.fastq")
	return files

SPLITFASTQFILES = getSplitFastqOutputs()

def getBamFilesForProject(wildcards):
	files = list()
	for s in config["samples"]:
		if(s["project"]==wildcards.project):
			files.append(os.path.join(MAINDIR,wildcards.project,config["align_folder"],s["name"]+".bam"))
	return files

def getSym2refForProject(wildcards):
	for s in config["samples"]:
		if(s["project"]==wildcards.project):
			return config["ref_folder"]+"/"+s["species"]+"/"+s["species"]+"_sym2ref.dat"

def getFastaForProject(wildcards):
	for s in config["samples"]:
		if(s["project"]==wildcards.project):
			return config["ref_folder"]+"/"+s["species"]+"/"+s["species"]+"_ERCC_chrm_polyAstrip.fa"

def getFastaIndexForProject(wildcards):
	for s in config["samples"]:
		if(s["project"]==wildcards.project):
			return config["ref_folder"]+"/"+s["species"]+"/"+s["species"]+"_ERCC_chrm_polyAstrip.fa.bwt"

def getFastqcFilesForProject(wildcards):
	files = list()
	for s in config["samples"]:
		if(s["project"]==wildcards.project):
			files.append(os.path.join(MAINDIR,wildcards.project,config["fastqc_folder"],s["name"]+"_fastqc.html"))
	return files

def getFastqcFilesForRun(wildcards):
	files = list()
	for s in config["samples"]:
		files.append(os.path.join(MAINDIR,s["project"],config["fastqc_folder"],s["name"]+"_fastqc.html"))
	return files

def getAllFilesForReport(wildcards):
	files = dict()
	files["multiqc"] = MAINDIR+"/"+wildcards.project+"/"+config["multiqc_folder"]+"/multiqc_report.html"
	files["templateCopied"] = MAINDIR+"/"+wildcards.project+"/"+config["report_folder"]+"/templateCopied.txt"
	files["well_summary"] = MAINDIR+"/"+wildcards.project+"/"+config["expression_folder"]+"/"+wildcards.project+".unq.well_summary.dat"
	if(wildcards.project in SECONDARY_ANALYSIS_PROJECTS or wildcards.project in NO_COMPARISONS_PROJECTS):
		files["PCA"] = MAINDIR+"/"+wildcards.project+"/"+config["de_folder"]+"/PCA.png"
		files["heatmap"] = MAINDIR+"/"+wildcards.project+"/"+config["de_folder"]+"/HeatmapCorPearson.png"
	if(wildcards.project in SECONDARY_ANALYSIS_PROJECTS):
		files["deseqTables"] = list()
		files["annotation"] = list()
		for comp in config["comparisons"][wildcards.project]["comps"]:
			files["deseqTables"].append(os.path.join(MAINDIR,wildcards.project,config["de_folder"],comp["condition1"]+"__vs__"+comp["condition2"],"DEseqResFiltered.tsv"))
			files["annotation"].append(os.path.join(MAINDIR,wildcards.project,config["de_folder"],comp["condition1"]+"__vs__"+comp["condition2"],"gseGo.txt"))
			files["annotation"].append(os.path.join(MAINDIR,wildcards.project,config["de_folder"],comp["condition1"]+"__vs__"+comp["condition2"],"gseKegg.txt"))
		files["index"] = os.path.join(wfbasedir,"TEMPLATE","index_SA.html")
	if(wildcards.project in NO_COMPARISONS_PROJECTS):
		files["index"] = os.path.join(wfbasedir,"TEMPLATE","index_SAwoComp.html")
	if((not wildcards.project in NO_COMPARISONS_PROJECTS) and (not wildcards.project in SECONDARY_ANALYSIS_PROJECTS)):
		files["index"] = os.path.join(wfbasedir,"TEMPLATE","index_PA.html")
	return files

# get final files to drive analysis in rule all
def getTargetFiles():
	targets = list()

	# Targets for all projects
	for p in PROJECTS:
		targets.append(os.path.join(MAINDIR,p,config["multiqc_folder"],"multiqc_report.html"))
		targets.extend(expand(os.path.join(MAINDIR,p,config["expression_folder"],p+".{exp}.well_summary.pdf"), exp=["unq","all"]))
		for s in config["samples"]:
			if(s["project"]==p):
				targets.append(os.path.join(MAINDIR,p,config["align_folder"],s["name"]+".bai"))
		targets.append(os.path.join(MAINDIR,p,"report.html"))

	# Targets for projects with secondary analysis
	for project in SECONDARY_ANALYSIS_PROJECTS.union(NO_COMPARISONS_PROJECTS):
		targets.append(os.path.join(MAINDIR,project,config["de_folder"],"exprDatUPM.tsv"))
		
	# Targets for whole run
	targets.append(MAINDIR+"/"+config["multiqc_folder"]+"/multiqc_report.html")
	targets.extend(expand(os.path.join(MAINDIR,config["expression_folder"],"run.{exp}.well_summary.pdf"),exp=["unq","all"]))
	targets.append(MAINDIR+"/config_used_in_analysis.json")

	return targets

localrules: barplot

onstart:
	if not os.path.exists(MAINDIR):
		os.makedirs(MAINDIR)
	with open(os.path.join(MAINDIR,"RUNINFO.txt"), 'a+', buffering=1) as out:
		out.write("Run launched on :\n")
		out.write(DATE+"\n\n")
		out.write("Run launched by :\n")
		out.write(os.environ['USER']+"\n\n")
		out.write("Commit ID used in this analysis:\n")
		subprocess.run("git --git-dir="+wfbasedir+"/.git rev-parse HEAD", shell=True, stdout=out)
		out.write("If you want to use DGE-pipeline at the same stage as it was run for this analysis:\n")
		out.write("git clone \"https://gitlab.univ-nantes.fr/bird_pipeline_registry/DGE-pipeline.git\"\n")
		out.write("cd DGE-pipeline\n")
		out.write("git checkout <Commit_ID_above>\n")
		out.write("-----------------------------------------------------\n")

rule all:
	input:	
		getTargetFiles()

rule copyConfig:
	input:
		config["conf"]
	output:
		MAINDIR+"/config_used_in_analysis.json"
	shell:
		"""
		cp {input} {output}
		"""

rule buildReport:
	input:
		unpack(getAllFilesForReport)
	output:
		report = MAINDIR+"/{project}/report.html",
	shell:
		"""
		python {wfbasedir}/SCRIPTS/make_html.py -c {config[conf]} -t {input.index} -p {wildcards.project} -o {output.report}
		"""

rule copyTemplateFolder:
	output:
		temp(touch(MAINDIR+"/{project}/"+config["report_folder"]+"/templateCopied.txt"))
	params:
		outputFolder = MAINDIR+"/{project}/"+config["report_folder"]
	shell:
		"""
		if test -d {params.outputFolder}; then rm -rf {params.outputFolder};fi
		mkdir {params.outputFolder}
		cp -r {wfbasedir}/TEMPLATE/assets {wfbasedir}/TEMPLATE/images {params.outputFolder}
		"""

rule deAnnot:
	input:
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/DEseqRes.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/DEseqResFiltered.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/annotDbInstalled.tmp"
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/gseGo.txt",
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/gseKegg.txt"
	params:
		inDir = MAINDIR+"/{project}/"+config["de_folder"],
		build = lambda wildcards: config["comparisons"][wildcards.project]["species"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/annot.R {params.inDir} {wildcards.cond1} {wildcards.cond2} {params.build} {wfbasedir}/SCRIPTS/DE/corresIDorg.txt
		"""

rule downloadAnnotRef:
	output:
		temp(touch(MAINDIR+"/{project}/"+config["de_folder"]+"/annotDbInstalled.tmp"))
	params:
		build = lambda wildcards: config["comparisons"][wildcards.project]["species"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/installAnnotDatabase.R {params.build} {wfbasedir}/SCRIPTS/DE/corresIDorg.txt
		"""

rule deAnalysis:
	input:
		sampleTable = MAINDIR+"/{project}/"+config["de_folder"]+"/sampletable.tsv",
		dds = MAINDIR+"/{project}/"+config["de_folder"]+"/dds.RData",
		heatmap = MAINDIR+"/{project}/"+config["de_folder"]+"/heatmap.RData"
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/DEseqRes.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/DEseqResFiltered.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/clustDEgene.png",
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/MA-plot.png",
		MAINDIR+"/{project}/"+config["de_folder"]+"/{cond1}__vs__{cond2}/Volcano-plot.png"
	params:
		inDir = MAINDIR+"/{project}/"+config["de_folder"],
		minLogFC = lambda wildcards: config["comparisons"][wildcards.project]["minLogFC"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/deAnalysis.R {params.inDir} {wildcards.cond1} {wildcards.cond2} {params.minLogFC}
		"""

rule deClusterHeatmap:
	input:
		sampleTable = MAINDIR+"/{project}/"+config["de_folder"]+"/sampletable.tsv",
		expressionTable = MAINDIR+"/{project}/"+config["de_folder"]+"/exprTransformed.tsv",
		sampleAbstract = MAINDIR+"/{project}/"+config["de_folder"]+"/SamplesAbstract.tsv"
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/HeatmapCorPearson.png",
		MAINDIR+"/{project}/"+config["de_folder"]+"/heatmap.RData"
	params:
		outdir = MAINDIR+"/{project}/"+config["de_folder"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/clusterHeatmap.R {input.expressionTable} {input.sampleTable} {input.sampleAbstract} {params.outdir}
		"""

rule deQualityControls:
	input:
		sampleTable = MAINDIR+"/{project}/"+config["de_folder"]+"/sampletable.tsv",
		expressionTable = MAINDIR+"/{project}/"+config["de_folder"]+"/exprTransformed.tsv"
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/PCA.png"
	params:
		outdir = MAINDIR+"/{project}/"+config["de_folder"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/qualityControl.R {input.expressionTable} {input.sampleTable} {params.outdir}
		"""

rule deNormalization:
	input:
		sampleTable = MAINDIR+"/{project}/"+config["de_folder"]+"/sampletable.tsv",
		expressionTable = MAINDIR+"/{project}/"+config["de_folder"]+"/exprFiltered.tsv"
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/exprNormalized.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/exprTransformed.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/dds.RData"
	params:
		outdir = MAINDIR+"/{project}/"+config["de_folder"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/normalization.R {input.expressionTable} {input.sampleTable} {params.outdir}
		"""

rule deGenerateSampleTable:
	input:
		sk = MAINDIR+"/{project}/"+config["de_folder"]+"/samplesToKeep.txt"
	output:
		st = MAINDIR+"/{project}/"+config["de_folder"]+"/sampletable.tsv"
	run:
		with open(config["conf"], "r") as f:
			c = json.load(f)
		with open(output.st, "w") as fout:
			print("\t".join(["SampleName","Condition"]),file=fout)
			with open(input.sk) as f:
				for line in f:
					for s in c["samples"]:
						if(s["name"]==line.strip("\n") and wildcards.project==s["project"]):
							print("\t".join([line.strip("\n"),s["condition"]]),file=fout)

rule deUPM:
	input:
		MAINDIR+"/{project}/"+config["de_folder"]+"/exprFiltered.tsv"
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/exprDatUPM.tsv"
	params:
		outdir = MAINDIR+"/{project}/"+config["de_folder"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/UPM.R {input} {params.outdir}
		"""	

rule deFilter:
	input:
		exprDat = MAINDIR+"/{project}/"+config["expression_folder"]+"/{project}.unq.refseq.umi.dat",
		conf = config["conf"]
	output:
		MAINDIR+"/{project}/"+config["de_folder"]+"/SamplesAbstract.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/exprFiltered.tsv",
		MAINDIR+"/{project}/"+config["de_folder"]+"/samplesToKeep.txt",
		MAINDIR+"/{project}/"+config["de_folder"]+"/samplesToRemove.txt"
	params:
		outdir = MAINDIR+"/{project}/"+config["de_folder"],
		minRep = lambda wildcards: config["comparisons"][wildcards.project]["minRep"],
		minGenes = lambda wildcards: config["comparisons"][wildcards.project]["minGenes"],
		minReads = lambda wildcards: config["comparisons"][wildcards.project]["minReads"]
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/DE/filter.R {input.exprDat} {params.outdir} {params.minRep} {params.minGenes} {params.minReads}
		"""

rule fastqc:
	input:
		MAINDIR+"/{project}/"+config["cutadapt_folder"]+"/{sample}.fastq.gz"
	output:
		MAINDIR+"/{project}/"+config["fastqc_folder"]+"/{sample}_fastqc.html"
	params:
		outdir = MAINDIR+"/{project}/"+config["fastqc_folder"]
	shell:
		"""
		fastqc -o {params.outdir} {input}
		"""

rule runStats:
	input:	
		MAINDIR+"/unknown.fastq.gz",
		expand(MAINDIR+"/{project}/"+config["expression_folder"]+"/{project}.{exp}.{suffix}", project=PROJECTS, exp=["unq","all"], suffix=finalSuffixes)
	output:
		expand(MAINDIR+"/"+config["expression_folder"]+"/run.{exp}.log.dat", exp=["unq","all"]),
		expand(MAINDIR+"/"+config["expression_folder"]+"/run.{exp}.well_summary.dat", exp=["unq","all"])
	params:
		outdir = MAINDIR+"/"+config["expression_folder"]
	shell:
		"""
		python {wfbasedir}/SCRIPTS/runStats.py {config[conf]} {params.outdir}
		"""

rule runMultiqc:
	input:	
		getFastqcFilesForRun
	output:
		MAINDIR+"/"+config["multiqc_folder"]+"/multiqc_report.html"
	params:
		fastqcdirs = " ".join(expand(MAINDIR+"/{project}/"+config["fastqc_folder"],project=PROJECTS)),
		outdir = MAINDIR+"/"+config["multiqc_folder"]
	shell:
		"""
		multiqc -f -e general_stats -e tophat -e bowtie2 {params.fastqcdirs} -o {params.outdir}
		"""

rule multiqc:
	input:
		getFastqcFilesForProject
	output:
		MAINDIR+"/{project}/"+config["multiqc_folder"]+"/multiqc_report.html"
	params:
		fastqcdir = MAINDIR+"/{project}/"+config["fastqc_folder"],
		outdir = MAINDIR+"/{project}/"+config["multiqc_folder"]
	shell:
		"""
		multiqc -f -e general_stats -e tophat -e bowtie2 {params.fastqcdir} -o {params.outdir}
		"""

rule barplot:
	input:
		"{path}.well_summary.dat"
	output:
		"{path}.well_summary.pdf"
	shell:
		"""
		Rscript {wfbasedir}/SCRIPTS/barplot.R {input} {output}
		"""

rule unknown_index:
	input:
		unknown = MAINDIR+"/unknown.fastq.gz"
	output:
		count = MAINDIR+"/unknown_indexes.count"
	params:
		tmpdir = MAINDIR
	shell:
		"gunzip -c {input.unknown} | awk 'NR%4==1' | cut -d\":\" -f8 | cut -c 1-6 | LC_ALL=C sort -T {params.tmpdir} | uniq -c | sed -e 's/^[ ]*//' | sort -k1,1nr > {output.count}"

rule count:
	input:
		getBamFilesForProject,
		getSym2refForProject,
		ercc = erccFile
	output:
		expand("{maindir}/{{project}}/{expression_folder}/{{project}}.{exp}.{suffix}",maindir=MAINDIR, expression_folder=config["expression_folder"], exp=["all","unq"], suffix=finalSuffixes)
	shell:
		"""
		python {wfbasedir}/SCRIPTS/count.py {config[conf]} {wildcards.project}
		"""

rule download_refMrna:
	input:

	output:
		mrna = expand("{refdir}/{{ref}}/{{ref}}_refMrna.fa.gz", refdir=config["ref_folder"])
	shell:
		"wget \"http://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.ref}/bigZips/refMrna.fa.gz\" -O \"{output.mrna}\""

rule download_refGene:
	input:

	output:
		genes = expand("{refdir}/{{ref}}/{{ref}}_refGene.txt.gz", refdir=config["ref_folder"])
	shell:
		"wget \"http://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.ref}/database/refGene.txt.gz\" -O \"{output.genes}\""

rule download_chrM:
	input:

	output:
		chrm = expand("{refdir}/{{ref}}/{{ref}}_chrM.fa.gz", refdir=config["ref_folder"])
	shell:
		"wget \"http://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.ref}/chromosomes/chrM.fa.gz\" -O \"{output.chrm}\""

rule unzip_refMrna:
	input:
		mrna = expand("{refdir}/{{ref}}/{{ref}}_refMrna.fa.gz", refdir=config["ref_folder"])
	output:
		expand("{refdir}/{{ref}}/{{ref}}_refMrna.fa", refdir=config["ref_folder"])
	shell:
		"gunzip {input.mrna}"

rule unzip_refGenes:
	input:
		genes = expand("{refdir}/{{ref}}/{{ref}}_refGene.txt.gz", refdir=config["ref_folder"])
	output:
		expand("{refdir}/{{ref}}/{{ref}}_refGene.txt", refdir=config["ref_folder"])
	shell:
		"gunzip {input.genes}"

rule unzip_chrM:
	input:
		chrm = expand("{refdir}/{{ref}}/{{ref}}_chrM.fa.gz", refdir=config["ref_folder"])
	output:
		expand("{refdir}/{{ref}}/{{ref}}_chrM.fa", refdir=config["ref_folder"])
	shell:
		"gunzip {input.chrm}"

rule sym2ref:
	input:
		expand("{refdir}/{{ref}}/{{ref}}_refGene.txt", refdir=config["ref_folder"])
	output:
		expand("{refdir}/{{ref}}/{{ref}}_sym2ref.dat", refdir=config["ref_folder"])
	shell:
		"python {wfbasedir}/SCRIPTS/preprocess_refgene.py {input} > {output}"

rule polyAstrip:
	input:
		genes = expand("{refdir}/{{ref}}/{{ref}}_refMrna.fa", refdir=config["ref_folder"]),
		chrm = expand("{refdir}/{{ref}}/{{ref}}_chrM.fa", refdir=config["ref_folder"])
	output:
		genes_pas = expand("{refdir}/{{ref}}/{{ref}}_polyAstrip.fa", refdir=config["ref_folder"]),
		chrm_pas = expand("{refdir}/{{ref}}/{{ref}}_chrM_polyAstrip.fa", refdir=config["ref_folder"])
	shell:
		"cat {input.genes} | awk '/^>/ {{printf(\"%s%s\\t\",(N>0?\"\\n\":\"\"),$0);N++;next;}} {{printf(\"%s\",toupper($0));}} END {{printf(\"\\n\");}}' | sed -r 's/A{{5,}}$/AAAAA/' | awk -F'\\t' 'length($2)>300' | tr \"\\t\" \"\\n\" | fold -w 60 > {output.genes_pas} && cat {input.chrm} | awk '/^>/ {{printf(\"%s%s\\t\",(N>0?\"\\n\":\"\"),$0);N++;next;}} {{printf(\"%s\",toupper($0));}} END {{printf(\"\\n\");}}' | sed -r 's/A{{5,}}$/AAAAA/' | awk -F'\\t' 'length($2)>300' | tr \"\\t\" \"\\n\" | fold -w 60 > {output.chrm_pas}"

rule add_ercc_and_chrm:
	input:
		stripgenome = expand("{refdir}/{{ref}}/{{ref}}_polyAstrip.fa", refdir=config["ref_folder"]),
		chrm = expand("{refdir}/{{ref}}/{{ref}}_chrM_polyAstrip.fa", refdir=config["ref_folder"])
	output:
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa", refdir=config["ref_folder"])
	shell:
		"cat {input.stripgenome} {erccFile} {input.chrm} > {output}"

rule index_ref:
	input:
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa", refdir=config["ref_folder"])
	output:
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa.bwt", refdir=config["ref_folder"]),
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa.amb", refdir=config["ref_folder"]),
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa.ann", refdir=config["ref_folder"]),
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa.pac", refdir=config["ref_folder"]),
		expand("{refdir}/{{ref}}/{{ref}}_ERCC_chrm_polyAstrip.fa.sa", refdir=config["ref_folder"])
	shell:
		"bwa index {input}"

rule index_bam:
	input:
		"{sample}.bam"
	output:
		"{sample}.bai"
	shell:	
		"""
		samtools index {input} {output}
		"""

rule align:
	input: 
		bwaindex = getFastaIndexForProject,
		fastq = MAINDIR+"/{project}/"+config["cutadapt_folder"]+"/{sample}.fastq.gz",
		ref = getFastaForProject,
		symtoref = getSym2refForProject
	output: 
		MAINDIR+"/{project}/"+config["align_folder"]+"/{sample}.bam"
	shell:
		"""
		bwa aln -M 0 {input.ref} {input.fastq} | bwa samse -n `cat {input.symtoref} | awk -F "," ' NF > max {{max = NF}} END{{print max}}'` {input.ref} - {input.fastq} | samtools sort -T {output}.tmp -o {output} -
		"""

rule zip_fastq:
	input:
		"{sample}.fastq"
	output:
		"{sample}.fastq.gz"
	shell:
		"""
		[ -s {input} ] || echo -e "@DUMMYSEQUENCE:NNNNNNNNNNNNNNNN\nN\n+\nA" > {input}
		gzip -f {input}
		"""

rule trim_polyA:
	input:
		MAINDIR+"/{project}/"+config["fastq_folder"]+"/{sample}.fastq"
	output:
		MAINDIR+"/{project}/"+config["cutadapt_folder"]+"/{sample}.fastq"
	shell:
		"""
		cutadapt -a "A{{10}}" -o {output} {input}
		"""

rule split_fastq:
	input: 
		expand("{fastq_pair[read1]} {fastq_pair[read2]}".split(), fastq_pair=config["fastq_pairs"])
	output: 
		temp(SPLITFASTQFILES)
	shell: 
		"python {wfbasedir}/SCRIPTS/split_fastq.py {config[conf]}"

