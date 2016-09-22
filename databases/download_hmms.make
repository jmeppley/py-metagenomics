##########################
#
# download_hmms.make
#
# Retrieves latest protein HMMs from a few different sources:
#  PFAM
#  COG
#  TIGRFAM
#
# Files are downloaded to ./seqdbs
#  change the location with SEQDB_ROOT=
#
# By default, only PFAM is downloaded. Control which DBs are downloaded with
# the GET_DB variables. EG:
#  GET_PFAM=False will skip PFAM
#  GET_TIGR=True will get TIGR
#  GET_COG=False will skip COG
#
###########################
SEQDB_ROOT?=./seqdbs
CDD_HMM_SCRIPT=./cdd-to-aln.prokka.pl

DBLIST=

# PFAM
GET_PFAM?=True
PFAM_FTP=ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
ifeq ($(GET_PFAM),True)
	DBLIST:=$(DBLIST) pfam
	PFAM_RELEASE?=$(shell curl -s $(PFAM_FTP)/relnotes.txt | perl -ne 'if (m/^\s*RELEASE\s+(\d[-_0-9.]+)/) { print $$1; }' | head -1)
endif

PFAM_DIR=$(SEQDB_ROOT)/PFAM/$(PFAM_RELEASE)
PFAM_HMM_NAME=Pfam-A.hmm
PFAM_HMM_DAT_NAME=$(PFAM_HMM_NAME).dat

# TIGR FAM
GET_TIGR?=False
TIGRFAM_FTP=ftp://ftp.jcvi.org/pub/data/TIGRFAMs
ifeq ($(GET_TIGR),True)
	DBLIST:=$(DBLIST) tigrfam
	TIGRFAM_README=$(shell curl -s -l $(TIGRFAM_FTP)/ | grep "RELEASE_NOTE_[0-9]" | sort -r | head -1)
	TIGRFAM_RELEASE=$(shell curl -s $(TIGRFAM_FTP)/$(TIGRFAM_README) | perl -ne 'if (m/\s+release\s+(\d[-_0-9.]+)/) { print "$$1\n"; }' | head -1)
	TIGRFAM_HMM_FILE=$(shell curl -s -l $(TIGRFAM_FTP)/ | grep "HMM.tar.gz" | sort | tail -n 1)
endif

TIGRFAM_DIR=$(SEQDB_ROOT)/TIGRFAM/$(TIGRFAM_RELEASE)
TIGRFAM_HMM_NAME=TIGRFAMS_HMM

# CDD files for COG (And maybe others in the future?)
CDD_README_URL=ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/README
CDD_FASTA_URL=ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/fasta.tar.gz
CDD_TABLE_URL=ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid_all.tbl.gz
CDD_TABLE=cddid_all.tbl.gz
CDD_FASTA=fasta.tar.gz
CDD_DATE=

# COG
GET_COG?=False
ifeq ($(GET_COG),True)
	DBLIST:=$(DBLIST) cog
	ifeq ($(CDD_DATE),)
		CDD_DATE=$(shell curl -s $(CDD_README_URL) | grep -m 1 revised | perl -pe 's/^.+revised\s+//; s/ /_/g;')
	endif
endif
CDD_DIR=$(SEQDB_ROOT)/cdd/$(CDD_DATE)
COG_DIR=$(CDD_DIR)/COG

all: $(DBLIST)

# universal rule
%.hmm_counts: %
	grep -c "^HMMER" $^ > $@

#PFAM
pfam: pfam_version $(PFAM_DIR)/relnotes.txt $(PFAM_DIR)/$(PFAM_HMM_NAME).h3p $(PFAM_DIR)/$(PFAM_HMM_DAT_NAME) $(PFAM_DIR)/$(PFAM_HMM_NAME).hmm_counts

pfam_version:
	@echo PFAM Release: $(PFAM_RELEASE)

$(PFAM_DIR):
	mkdir -p $@

$(PFAM_DIR)/relnotes.txt: | $(PFAM_DIR)
	curl -s $(PFAM_FTP)/relnotes.txt > $@

$(PFAM_DIR)/$(PFAM_HMM_NAME): | $(PFAM_DIR)
	curl -s $(PFAM_FTP)/$(PFAM_HMM_NAME).gz | gunzip -c > $@

$(PFAM_DIR)/$(PFAM_HMM_NAME).h3p: $(PFAM_DIR)/$(PFAM_HMM_NAME) | $(PFAM_DIR)
	hmmpress $^

$(PFAM_DIR)/$(PFAM_HMM_DAT_NAME): | $(PFAM_DIR)
	curl -s $(PFAM_FTP)/$(PFAM_HMM_DAT_NAME).gz | gunzip -c > $@

#TIGRFAM
tigrfam: tigrfam_version $(TIGRFAM_DIR)/relnotes.txt $(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME).h3p $(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME).hmm_counts

tigrfam_version:
	@echo TIGRFAM Release: $(TIGRFAM_RELEASE)

$(TIGRFAM_DIR):
	mkdir -p $@

$(TIGRFAM_DIR)/relnotes.txt: | $(TIGRFAM_DIR)
	curl -s $(TIGRFAM_FTP)/$(TIGRFAM_README) > $@

$(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME): $(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME).tar.gz | $(TIGRFAM_DIR)
	rm -f $@
	tar -tzvf $^ | sed -r 's/^.+\s(\S+)\s*$$/\1/' |  xargs tar -zxOf $^ > $@

$(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME).tar.gz: | $(TIGRFAM_DIR)
	curl -s $(TIGRFAM_FTP)/$(TIGRFAM_HMM_FILE) > $@

$(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME).h3p: $(TIGRFAM_DIR)/$(TIGRFAM_HMM_NAME) | $(TIGRFAM_DIR)
	hmmpress $^

# Generic CDD stuff
cdd_version:
	@echo CDD date: $(CDD_DATE)

$(CDD_DIR):
	mkdir -p $@

$(CDD_DIR)/README: | $(CDD_DIR)
	curl -s $(CDD_README_URL) > $@

$(CDD_DIR)/$(CDD_TABLE): | $(CDD_DIR)
	curl -s $(CDD_TABLE_URL) > $@

$(CDD_DIR)/$(CDD_FASTA): | $(CDD_DIR)
	curl -s $(CDD_FASTA_URL) > $@

# COG
cog: cdd_version $(CDD_DIR)/README $(COG_DIR)/COG.hmm.h3p $(COG_DIR)/COG.hmm.ascii.hmm_counts

$(COG_DIR):
	mkdir -p $@

$(COG_DIR)/COG.aln: $(CDD_DIR)/$(CDD_FASTA) $(CDD_DIR)/$(CDD_TABLE) | $(COG_DIR)
	$(CDD_HMM_SCRIPT) --srcdir $(CDD_DIR) --lib COG --dstdir $(COG_DIR) 
	sed -i -r 's/([0-9])AC /\1 AC /' $(COG_DIR)/COG.aln

$(COG_DIR)/COG.hmm.ascii: $(COG_DIR)/COG.aln
	hmmbuild $@ $^

$(COG_DIR)/COG.hmm: $(COG_DIR)/COG.hmm.ascii
	hmmconvert -b $^ > $@

$(COG_DIR)/COG.hmm.h3p: $(COG_DIR)/COG.hmm
	hmmpress $^

