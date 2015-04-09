# Define some basic locations 
# Edit these or pass in your own values when calling make
FTP_ROOT:=ftp://ftp.arb-silva.de
REL?=$(shell curl -s -l $(FTP_ROOT)/ | grep release_ | perl -pe 's/release_//' | sort -n -r | head -1)
REL:=$(REL)

DB_SCRIPT_DIR?=.
BUILD_ROOT?=./Silva
BUILDDIR:=$(BUILD_ROOT)/Silva-$(REL)

LASTDB_ROOT?=/minilims/galaxy-data/tool-data/sequencedbs/lastdb
LASTDB_DIR:=$(LASTDB_ROOT)/Silva/$(REL)
LASTDBCHUNK?=100G

# The following is a hack to do doa  greaterthan conditional in make
PARSE_TAXONOMY:=$(shell if [ "$(REL)" -gt "118" ]; then echo True; else echo False; fi)

# The file names have changed over time, this works for releases 111 to 119
FTP_BASE:=$(FTP_ROOT)/release_$(REL)/Exports
SSUREF_FILE:=$(shell curl -s -l $(FTP_BASE)/ | grep $(REL) | grep "gz$$" | grep "SSURef_" | grep "_N[Rr]" | grep tax_silva_trunc)
SSUREF_URL:=$(FTP_BASE)/$(SSUREF_FILE)
SSU_FASTA:=$(BUILDDIR)/Silva_$(REL)_SSURef_NR99_tax_silva_trunc.fasta
LSUREF_FILE:=$(shell curl -s -l $(FTP_BASE)/ | grep $(REL) | grep "gz$$" | grep "LSURef_" | grep tax_silva_trunc)
LSUREF_URL:=$(FTP_BASE)/$(LSUREF_FILE)
LSU_FASTA:=$(BUILDDIR)/Silva_$(REL)_LSURef_tax_silva_trunc.fasta
#COMBINED_FASTA:=$(BUILDDIR)/Silva_$(REL)_SSULSURef_tagged_DNA.fasta

# The taxonomy files were new in 115 and changed locations in 119
# For now, I'm just going to support the new way in this makefile
SSU_TAXFILE_URL:=$(FTP_BASE)/taxonomy/tax_slv_ssu_nr_$(REL).txt
SSU_TAXFILE=$(BUILDDIR)//Silva_$(REL)_SSURef_NR99_tax_silva_trunc.tax
LSU_TAXFILE_URL:=$(FTP_BASE)/taxonomy/tax_slv_lsu_$(REL).txt
LSU_TAXFILE=$(BUILDDIR)/Silva_$(REL)_LSURef_tax_silva_trunc.tax

LSU_LASTDB_DIR:=$(LASTDB_DIR)/Silva_$(REL)_LSURef
LSU_LASTP:=$(LSU_LASTDB_DIR)/lastdb
LSU_LASTFILE=$(LSU_LASTP).prj
LSU_TAXIDMAP=$(LSU_LASTP).tax
LSU_HITIDMAP:=$(LSU_LASTP).ids
LSU_DBTAX:=$(LSU_LASTDB_DIR)/names.dmp

SSU_LASTDB_DIR:=$(LASTDB_DIR)/Silva_$(REL)_SSURef_NR99
SSU_LASTP:=$(SSU_LASTDB_DIR)/lastdb
SSU_LASTFILE=$(SSU_LASTP).prj
SSU_TAXIDMAP=$(SSU_LASTP).tax
SSU_HITIDMAP:=$(SSU_LASTP).ids
SSU_DBTAX:=$(SSU_LASTDB_DIR)/names.dmp

TAXSCRIPT:=$(DB_SCRIPT_DIR)/buildSilvaTaxFiles.py

ifeq ($(PARSE_TAXONOMY),True)
	ALL_TARGETS:=$(ALL_TARGETS) taxonomy
endif

all: report lastdb

taxfiles: $(SSU_TAXFILE) $(LSU_TAXFILE)

$(SSU_TAXFILE):
	curl $(SSU_TAXFILE_URL) > $(SSU_TAXFILE)

$(LSU_TAXFILE):
	curl $(LSU_TAXFILE_URL) > $(LSU_TAXFILE)

dbtax: $(LSU_DBTAX) $(SSU_DBTAX) 

$(LSU_DBTAX): $(LSU_TAXFILE) $(LSU_FASTA) | $(LSU_LASTDB_DIR)
	export PYTHONPATH=$(DB_SCRIPT_DIR)/..; python $(TAXSCRIPT) $^ $|

$(SSU_DBTAX): $(SSU_TAXFILE) $(SSU_FASTA) | $(SSU_LASTDB_DIR)
	export PYTHONPATH=$(DB_SCRIPT_DIR)/..; python $(TAXSCRIPT) $^ $|

report:
	@echo Release: $(REL)
	@echo FTP: $(FTP_BASE)
	@echo SSU URL: $(SSUREF_URL)
	@echo LSU URL: $(LSUREF_URL)

fasta: report $(SSU_FASTA) $(LSU_FASTA) taxfiles
lastdb: $(LSU_HITIDMAP) $(LSU_LASTFILE) $(SSU_HITIDMAP) $(SSU_LASTFILE) dbtax

$(LSU_LASTDB_DIR):
	mkdir -p $@

$(LSU_LASTFILE): $(LSU_FASTA) | $(LSU_LASTDB_DIR)
	@echo Formatting database for lastal searches
	lastdb -v -c -s $(LASTDBCHUNK) $(LSU_LASTP) $<

$(LSU_HITIDMAP): $(LSU_FASTA) | $(LSU_LASTDB_DIR)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

$(SSU_LASTDB_DIR):
	mkdir -p $@

$(SSU_LASTFILE): $(SSU_FASTA) | $(SSU_LASTDB_DIR)
	@echo Formatting database for lastal searches
	lastdb -v -c -s $(LASTDBCHUNK) $(SSU_LASTP) $<

$(SSU_HITIDMAP): $(SSU_FASTA) | $(SSU_LASTDB_DIR)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

$(SSU_FASTA):
	@echo downloading SSU from $(SSUREF_URL)
	@mkdir -p $(BUILDDIR)
	curl $(SSUREF_URL) | gunzip -c | perl -ne 'if (m/^>/) { print; } else {tr/uU/tT/; print;}' > $(SSU_FASTA)

$(LSU_FASTA):
	@echo downloading LSU from $(LSUREF_URL)
	@mkdir -p $(BUILDDIR)
	curl $(LSUREF_URL) | gunzip -c | perl -ne 'if (m/^>/) { print; } else {tr/uU/tT/; print;}' > $(LSU_FASTA)
