ifeq ($(R1),)
$(error R1 undefined)
endif

ifeq ($(R2),)
$(error R2 undefined)
endif

ifeq ($(NJOBS),)
$(error NJOBS undefined)
endif

ifeq ($(OUTDIR),)
$(error OUTDIR undefined)
endif

ABSOUTDIR=$(abspath $(OUTDIR))

ifeq ($(ABSOUTDIR),)
$(error cannot resolve OUTDIR to full path)
endif

PREFIX=split.

$(ABSOUTDIR)/manifest.txt : $(ABSOUTDIR)/list1.txt $(ABSOUTDIR)/list2.txt
	cat $^ | xargs -l1 -t -r -P $(NJOBS) gzip -f --best
	paste $^| sed 's/_R1.fastq	/_R1.fastq.gz	/' | sed 's/_R2.fastq$$/_R2.fastq.gz/' > $@


$(ABSOUTDIR)/list1.txt : $(R1)
	rm -f "$(ABSOUTDIR)/$(PREFIX)*_R1.fastq"
	mkdir -p "$(ABSOUTDIR)"
	gunzip -c "$<" | split -l 16000000 --additional-suffix=_R1.fastq - "$(ABSOUTDIR)/$(PREFIX)"
	ls $(ABSOUTDIR)/$(PREFIX)*_R1.fastq | sort -V > $@

$(ABSOUTDIR)/list2.txt : $(R2)
	rm -f "$(ABSOUTDIR)/$(PREFIX)*_R2.fastq"
	mkdir -p "$(ABSOUTDIR)"
	gunzip -c "$<" | split -l 16000000 --additional-suffix=_R2.fastq - "$(ABSOUTDIR)/$(PREFIX)"
	ls $(ABSOUTDIR)/$(PREFIX)*_R2.fastq | sort -V > $@

