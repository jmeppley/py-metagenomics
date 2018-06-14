##########################
#
# download_silva.make
#
# Retrieves latest protein database and heirarchy from silva.
#
# Files are downloaded to ./Silva (set with BUILD_ROOT=)
# If requests, sequence databases for searching are compiled to ./seqdbs
#  change the location with SEQDB_ROOT=
#
# if BUILD_LASTDB is set to "True", 
#  a version is compiled for lastal
# if BUILD_BT2DB is set to "True", 
#  a version is compiled for bowtie2
# if BUILD_BWADB is set to "True", 
#  a version is compiled for BWA
#
###########################

# Define some basic locations 
# Edit these or pass in your own values when calling make
FTP_ROOT:=ftp://ftp.arb-silva.de
REL?=$(shell curl -s -l $(FTP_ROOT)/ | grep release_ | perl -pe 's/release_//' | sort -n -r | head -1)
REL:=$(REL)
REL_ROOT=$(basename $(word 1,$(subst _, ,$(REL))))
REL_TAX=$(subst _,.,$(REL))

DB_SCRIPT_DIR?=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
BUILD_ROOT?=./Silva
BUILDDIR:=$(BUILD_ROOT)/Silva-$(REL)

SEQDB_ROT?=./seqdbs
SEQDB_DIR:=$(SEQDB_ROOT)/Silva/$(REL)

BUILD_BT2DB:=False
BUILD_BWADB:=True
BUILD_LASTDB:=False
LASTDBCHUNK?=

ifeq ($(LASTDBCHUNK),)
	LASTDBCHUNK_OPTION:=
else
	LASTDBCHUNK_OPTION:= -s $(LASTDB_CHUNK)
endif

# The following is a hack to do a greaterthan conditional in make
PARSE_TAXONOMY:=$(shell if [ "$(REL_ROOT)" -gt "118" ]; then echo True; else echo False; fi)
#PARSE_TAXONOMY:=True

# The file names have changed over time, this works for releases 111 to 119
FTP_BASE:=$(FTP_ROOT)/release_$(subst .,_,$(REL))/Exports
SSUREF_FILE:=$(shell curl -s -l $(FTP_BASE)/ | grep $(REL) | grep "gz$$" | grep "SSURef_" | grep "_N[Rr]" | grep tax_silva_trunc)
SSU_FASTA:=$(BUILDDIR)/Silva_$(REL)_SSURef_NR99_tax_silva_trunc.fasta
ifeq ($(SSUREF_FILE),)
	SSU_REL=$(REL_ROOT)
	SSUREF_FILE:=$(shell curl -s -l $(FTP_BASE)/ | grep $(SSU_REL) | grep "gz$$" | grep "SSURef_" | grep "_N[Rr]" | grep tax_silva_trunc)
	SSU_FASTA:=$(BUILDDIR)/Silva_$(SSU_REL)_SSURef_NR99_tax_silva_trunc.fasta
else
	SSU_REL=$(REL)
endif
SSUREF_URL:=$(FTP_BASE)/$(SSUREF_FILE)
LSUREF_FILE:=$(shell curl -s -l $(FTP_BASE)/ | grep $(REL) | grep "gz$$" | grep "LSURef_" | grep tax_silva_trunc)
LSU_FASTA:=$(BUILDDIR)/Silva_$(REL)_LSURef_tax_silva_trunc.fasta
ifeq ($(LSUREF_FILE),)
	LSU_REL=$(REL_ROOT)
	LSUREF_FILE:=$(shell curl -s -l $(FTP_BASE)/ | grep $(LSU_REL) | grep "gz$$" | grep "LSURef_" | grep tax_silva_trunc)
	LSU_FASTA:=$(BUILDDIR)/Silva_$(LSU_REL)_LSURef_tax_silva_trunc.fasta
else
	LSU_REL=$(REL)
endif
LSUREF_URL:=$(FTP_BASE)/$(LSUREF_FILE)

# The taxonomy files were new in 115 and changed locations in 119
# For now, I'm only going to support the new way in this makefile
SSU_TAXFILE_URL:=$(FTP_BASE)/taxonomy/tax_slv_ssu_$(REL_TAX).txt
SSU_TAXFILE=$(BUILDDIR)/Silva_$(SSU_REL)_SSURef_NR99_tax_silva_trunc.tax
LSU_TAXFILE_URL:=$(FTP_BASE)/taxonomy/tax_slv_lsu_$(REL_TAX).txt
LSU_TAXFILE=$(BUILDDIR)/Silva_$(LSU_REL)_LSURef_tax_silva_trunc.tax

# LOCATIONS for compiled sequence dbs and associated files
LSU_SEQDB_DIR:=$(SEQDB_DIR)/Silva_$(LSU_REL)_LSURef
TAXIDMAP=acc_to_taxid.tsv
LSU_TAXIDMAP=$(LSU_SEQDB_DIR)/$(TAXIDMAP)
LSU_HITIDMAP:=$(LSU_SEQDB_DIR)/seqids_to_descriptions.tsv
LSU_DBTAX:=$(LSU_SEQDB_DIR)/names.dmp
SSU_SEQDB_DIR:=$(SEQDB_DIR)/Silva_$(SSU_REL)_SSURef_NR99
SSU_TAXIDMAP=$(SSU_SEQDB_DIR)/$(TAXIDMAP)
SSU_HITIDMAP:=$(SSU_SEQDB_DIR)/seqids_to_descriptions.tsv
SSU_DBTAX:=$(SSU_SEQDB_DIR)/names.dmp

# ..for lastdb
LSU_LASTP:=$(LSU_SEQDB_DIR)/lastdb
LSU_LASTDB_FILE=$(LSU_LASTP).prj
LSU_LASTDB_TAXIDMAP=$(LSU_LASTP).tax
LSU_LASTDB_HITIDMAP:=$(LSU_LASTP).ids
SSU_LASTP:=$(SSU_SEQDB_DIR)/lastdb
SSU_LASTDB_FILE=$(SSU_LASTP).prj
SSU_LASTDB_TAXIDMAP=$(SSU_LASTP).tax
SSU_LASTDB_HITIDMAP:=$(SSU_LASTP).ids

LASTDB_TARGETS=$(LSU_LASTDB_FILE) $(LSU_LASTDB_HITIDMAP) $(SSU_LASTDB_FILE) $(SSU_LASTDB_HITIDMAP)
ifeq ($(PARSE_TAXONOMY),True)
	LASTDB_TARGETS:=$(LASTDB_TARGETS) $(LSU_LASTDB_TAXIDMAP) $(LSU_TAXIDMAP) $(SSU_LASTDB_TAXIDMAP) $(SSU_TAXIDMAP)
endif

# ..for BWA
LSU_BWADB:=$(LSU_SEQDB_DIR)/bwadb
LSU_BWADB_TAXIDMAP=$(LSU_BWADB).tax
LSU_BWADB_HITIDMAP:=$(LSU_BWADB).ids
SSU_BWADB:=$(SSU_SEQDB_DIR)/bwadb
SSU_BWADB_TAXIDMAP=$(SSU_BWADB).tax
SSU_BWADB_HITIDMAP:=$(SSU_BWADB).ids

BWADB_TARGETS=$(LSU_BWADB).bwt $(LSU_BWADB_HITIDMAP) $(SSU_BWADB).bwt $(SSU_BWADB_HITIDMAP)
ifeq ($(PARSE_TAXONOMY),True)
	BWADB_TARGETS:=$(BWADB_TARGETS) $(LSU_BWADB_TAXIDMAP) $(LSU_TAXIDMAP) $(SSU_BWADB_TAXIDMAP) $(SSU_TAXIDMAP)
endif

# ..for bowtie2
LSU_BT2DB:=$(LSU_SEQDB_DIR)/bowtie2db
LSU_BT2DB_FILE=$(LSU_BT2DB).1.bt2
LSU_BT2TDB_TAXIDMAP=$(LSU_BT2DB).tax
LSU_BT2DB_HITIDMAP:=$(LSU_BT2DB).ids
SSU_BT2DB:=$(SSU_SEQDB_DIR)/bowtie2db
SSU_BT2DB_FILE=$(SSU_BT2DB).1.bt2
SSU_BT2DB_TAXIDMAP=$(SSU_BT2DB).tax
SSU_BT2DB_HITIDMAP:=$(SSU_BT2DB).ids

BT2DB_TARGETS=$(LSU_BT2DB_FILE) $(LSU_BT2DB_HITIDMAP) $(SSU_BT2DB_FILE) $(SSU_BT2DB_HITIDMAP)
ifeq ($(PARSE_TAXONOMY),True)
	BT2DB_TARGETS:=$(BT2DB_TARGETS) $(LSU_BT2DB_TAXIDMAP) $(LSU_TAXIDMAP) $(SSU_BT2DB_TAXIDMAP) $(SSU_TAXIDMAP)
endif

TAXSCRIPT:=$(DB_SCRIPT_DIR)/buildSilvaTaxFiles.py

ALL_TARGETS:=report fasta
ifeq ($(PARSE_TAXONOMY),True)
	ALL_TARGETS:=$(ALL_TARGETS) taxfiles
endif

DB_TARGETS:=
ifeq ($(BUILD_LASTDB),True)
	DB_TARGETS:=$(LASTDB_TARGETS)
endif
ifeq ($(BUILD_BWADB),True)
	DB_TARGETS:=$(DB_TARGETS) $(BWADB_TARGETS)
endif
ifeq ($(BUILD_BT2DB),True)
	DB_TARGETS:=$(DB_TARGETS) $(BT2DB_TARGETS)
endif

all: $(ALL_TARGETS) $(DB_TARGETS)

report:
	@echo Release: $(REL), rel_root: $(REL_ROOT), parse tax: $(PARSE_TAXONOMY)
	@echo FTP: $(FTP_BASE)
	@echo SSU URL: $(SSUREF_URL)
	@echo SSU release string: $(SSU_REL)
	@echo LSU URL: $(LSUREF_URL)
	@echo LSU release string: $(LSU_REL)
	@echo Parse Taxonomy: $(PARSE_TAXONOMY)
	@echo Build BWA: $(BUILD_BWADB)
	@echo Build Bowtie: $(BUILD_BT2DB)
	@echo Build lastdb: $(BUILD_LASTDB)

# Download fasta and supporting files
fasta: report $(SSU_FASTA).stats $(LSU_FASTA).stats

$(BUILDDIR):
	@mkdir -p $(BUILDDIR)

$(SSU_FASTA): | $(BUILDDIR)
	@echo downloading SSU from $(SSUREF_URL)
	curl $(SSUREF_URL) | gunzip -c | perl -ne 'if (m/^>/) { print; } else {tr/uU/tT/; print;}' > $(SSU_FASTA)

$(LSU_FASTA): | $(BUILDDIR)
	@echo downloading LSU from $(LSUREF_URL)
	curl $(LSUREF_URL) | gunzip -c | perl -ne 'if (m/^>/) { print; } else {tr/uU/tT/; print;}' > $(LSU_FASTA)

%.fasta.stats: %.fasta
	cat $^ | prinseq-lite.pl -fasta stdin -stats_len -stats_info > $@

# Download tax files and build custom taxdump for DB searches
taxfiles: $(SSU_TAXFILE) $(LSU_TAXFILE)

$(SSU_TAXFILE): | $(BUILDDIR)
	curl $(SSU_TAXFILE_URL) > $(SSU_TAXFILE)

$(LSU_TAXFILE): | $(BUILDDIR)
	curl $(LSU_TAXFILE_URL) > $(LSU_TAXFILE)

$(LSU_TAXIDMAP): $(LSU_TAXFILE) $(LSU_FASTA) | $(LSU_SEQDB_DIR)
	$(TAXSCRIPT) -o $(TAXIDMAP) -t $(LSU_TAXFILE) $(LSU_FASTA) $(LSU_SEQDB_DIR)

$(SSU_TAXIDMAP): $(SSU_TAXFILE) $(SSU_FASTA) | $(SSU_SEQDB_DIR)
	$(TAXSCRIPT) -o $(TAXIDMAP) -t $(SSU_TAXFILE) $(SSU_FASTA) $(SSU_SEQDB_DIR)

# Initialize DB dirs and create supporting files

$(LSU_SEQDB_DIR):
	mkdir -p $@

$(SSU_SEQDB_DIR):
	mkdir -p $@

$(LSU_HITIDMAP): $(LSU_FASTA) | $(LSU_SEQDB_DIR)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

$(SSU_HITIDMAP): $(SSU_FASTA) | $(SSU_SEQDB_DIR)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

# Lastdb

$(LSU_LASTDB_FILE): $(LSU_FASTA) | $(LSU_SEQDB_DIR)
	@echo Formatting database for lastal searches
	lastdb -v -c $(LASTDBCHUNK_OPTION) $(LSU_LASTP) $<

$(SSU_LASTDB_FILE): $(SSU_FASTA) | $(SSU_SEQDB_DIR)
	@echo Formatting database for lastal searches
	lastdb -v -c $(LASTDBCHUNK_OPTION) $(SSU_LASTP) $<

$(LSU_LASTDB_HITIDMAP): $(LSU_HITIDMAP)
	ln -s $^ $@

$(SSU_LASTDB_HITIDMAP): $(SSU_HITIDMAP)
	ln -s $^ $@

$(LSU_LASTDB_TAXIDMAP): $(LSU_TAXIDMAP)
	ln -s $^ $@

$(SSU_LASTDB_TAXIDMAP): $(SSU_TAXIDMAP)
	ln -s $^ $@

# BWA

$(LSU_BWADB).bwt: $(LSU_FASTA) | $(LSU_SEQDB_DIR)
	@echo Formatting database for BWA searches
	bwa index -p $(LSU_BWADB) $(LSU_FASTA)

$(SSU_BWADB).bwt: $(SSU_FASTA) | $(SSU_SEQDB_DIR)
	@echo Formatting database for BWA searches
	bwa index -p $(SSU_BWADB) $(SSU_FASTA)

$(LSU_BWADB_HITIDMAP): $(LSU_HITIDMAP)
	ln -s $^ $@

$(SSU_BWADB_HITIDMAP): $(SSU_HITIDMAP)
	ln -s $^ $@

$(LSU_BWADB_TAXIDMAP): $(LSU_TAXIDMAP)
	ln -s $^ $@

$(SSU_BWADB_TAXIDMAP): $(SSU_TAXIDMAP)
	ln -s $^ $@

# Bowtie2

$(LSU_BT2DB_FILE): $(LSU_FASTA) | $(LSU_SEQDB_DIR)
	@echo Formatting database for bowtie2 searches
	bowtie-build $(LSU_FASTA) $(LSU_BT2DB)

$(SSU_BT2DB_FILE): $(SSU_FASTA) | $(SSU_SEQDB_DIR)
	@echo Formatting database for bowtie2 searches
	bowtie-build $(SSU_FASTA) $(SSU_BT2DB)

$(LSU_BT2DB_HITIDMAP): $(LSU_HITIDMAP)
	ln -s $^ $@

$(SSU_BT2DB_HITIDMAP): $(SSU_HITIDMAP)
	ln -s $^ $@

$(LSU_BT2DB_TAXIDMAP): $(LSU_TAXIDMAP)
	ln -s $^ $@

$(SSU_BT2DB_TAXIDMAP): $(SSU_TAXIDMAP)
	ln -s $^ $@

