##########################
#
# download_refseq.make
#
# Retrieves latest nonredundant protein database from refseq.
#
# to get an older version set the relese number with, eg: REL=75
#
# Files are downloaded to ./RefSeq (set with BUILD_ROOT=)
# if BUILD_LASTDB is set to "True", 
#  a version is compiled for lastal in ./lastdb (set with LASTDB_ROOT=)
#
###########################

# Define some basic locations 
# Edit these or pass in your own values when calling make

# what to download
FTP_ROOT=ftp://ftp.ncbi.nlm.nih.gov/refseq/release
REL?=$(shell curl $(FTP_ROOT)/RELEASE_NUMBER)
REL:=$(REL)

# where to put it
DB_SCRIPT_DIR?=.
SCRIPT_DIR?=..
BUILD_ROOT?=./RefSeq
RSDIR:=$(BUILD_ROOT)/RefSeq-$(REL)

# compile for lastal?
BUILD_LASTDB:=True
LASTDB_ROOT?=./lastdb
LASTDB_DIR:=$(LASTDB_ROOT)/RefSeq/$(REL)
LASTDBCHUNK?=

# Most folks won't need to edit below this line

# command line uption for running lastdb
ifeq ($(LASTDBCHUNK),)
	LASTDBCHUNK_OPTION:=
else
	LASTDBCHUNK_OPTION:= -s $(LASTDB_CHUNK)
endif

# Define the layout of the build directory
MDDIR:=$(RSDIR)/metadata
TAXDUMP_SOURCE:=$(RSDIR)/taxdump

COMPLETEFAA=$(RSDIR)/complete.protein.fasta

# Final protein outputs
FAA_NAME:=RefSeq-$(REL).AllProteins.faa
FAA_PREREQS=$(COMPLETEFAA)
FAA:=$(RSDIR)/$(FAA_NAME)

# LAST database will be packaged in it's own directory
LASTDIR=$(LASTDB_DIR)/$(FAA_NAME).ldb
LASTP=$(LASTDIR)/lastdb
LASTFILE=$(LASTP).prj
TAXDUMP_DB=$(LASTDIR)/nodes.dmp

ACCPREFREL=acc.to.taxid
ACCPREF=$(RSDIR)/$(ACCPREFREL)
ACCPREFP=$(ACCPREF).protein
ACCTAXMAPDB:=$(LASTP).tax
HITIDMAP:=$(LASTP).ids
ACCTAXMAP:=$(ACCPREFP)

TAXMAPSCRIPT=$(DB_SCRIPT_DIR)/buildRefSeqAccToTaxidMap.py

##
# Build the arguments for all
ifeq ($(BUILD_LASTDB),False)
	ALL_TARGETS:=fasta $(ACCTAXMAP)
	TAXDUMP=$(TAXDUMP_SOURCE)
else
	ALL_TARGETS:=lastdb $(ACCTAXMAPDB)
	TAXDUMP=$(TAXDUMP_DB)
endif
ALL_TARGETS:=$(ALL_TARGETS) $(TAXDUMP)

all: report $(MDDIR) $(ALL_TARGETS)

lastdb: $(LASTFILE) $(HITIDMAP)
fasta: $(FAA)

report:
	@echo RefSeq release number is: $(REL)
	@echo Building database in: $(RSDIR)
	@echo "Output fasta is $(FAA)"
	@echo "BUILD_LASTDB is $(BUILD_LASTDB)"
	@echo "all target list is $(ALL_TARGETS)"
	@if [ "$(BUILD_LASTDB)" != "False" ]; then echo Final database written to $(LASTDB_ROOT); else echo "Lastdb formatting will be skipped"; fi

$(LASTDIR):
	mkdir -p $(LASTDIR)

$(LASTFILE): $(FAA) | $(LASTDIR)
	@echo "==Formating last: $@"
	lastdb -v -c -p $(LASTDBCHUNK_OPTION) $(LASTP) $(FAA)

$(FAA): $(FAA_PREREQS)
	@echo "==Masking low complexity with tantan"
	tantan -p $^ | perl -ne 'if (m/^>(?!gi\|\d+)(.*)$$/) { if (defined $$n) { $$n++; } else { $$n=10000000000; } print ">gi|$$n|loc|$$1\n"; } else { print; }' > $@

%.protein.fasta: %/.download.complete.aa
	@echo "==Compiling $@ from gz archives"
	for FILE in $(RSDIR)/$(*F)/complete.[0-9]*.protein.gpff.gz; do gunzip -c $$FILE; done | python $(SCRIPT_DIR)/getSequencesFromGbk.py -F fasta -r > $@

$(RSDIR)/complete/.download.complete.aa:
	@echo "==Dowloading complete RefSeq proteins"
	mkdir -p $(RSDIR)/complete
	#cd $(RSDIR)/complete && wget -c $(FTP_ROOT)/complete/complete.nonredundant_protein.*.protein.gpff.gz
	cd $(RSDIR)/complete && wget -c $(FTP_ROOT)/complete/complete.[0-9]*.protein.gpff.gz
	touch $@

$(ACCTAXMAP).oldway: $(MDDIR) $(TAXDUMP_SOURCE) $(TAXMAPSCRIPT) 
	# The catalog only has one taxid even for multispecies entries, so
	# now we get the taxid maps from the gpff files
	export PYTHONPATH=$(DB_SCRIPT_DIR)/.. && gunzip -c $(MDDIR)/RefSeq-release$(REL).catalog.gz | python $(TAXMAPSCRIPT) $(TAXDUMP_SOURCE) > $@

$(ACCTAXMAP): $(RSDIR)/complete/.download.complete.aa
	# For some multispecies entries, you'll get multiple lines in the tax map
	gunzip -c $(RSDIR)/complete/complete.[0-9]*.protein.gpff.gz | perl -ne 'if (m/^ACCESSION\s+(\S+)\b/) { $$acc=$$1; } elsif (m/db_xref="taxon:(\d+)"/) { print "$$acc\t$$1\n"; }' > $@

$(ACCTAXMAPDB): $(ACCTAXMAP) | $(LASTDIR)
	cp $< $@

$(MDDIR):
	@echo "==Downloading metadata"
	mkdir -p $(MDDIR)
	cd $(MDDIR) && wget -c $(FTP_ROOT)/release-notes/RefSeq*.txt $(FTP_ROOT)/release-statistics/RefSeq-release*.stats.txt $(FTP_ROOT)/release-catalog/RefSeq-release$(REL).catalog.gz $(FTP_ROOT)/release-catalog/release$(REL)*

$(TAXDUMP_SOURCE):
	@echo "==Downloading taxonomy"
	mkdir -p $(TAXDUMP_SOURCE)
	cd $(TAXDUMP_SOURCE) && wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -zxvf taxdump.tar.gz

$(TAXDUMP_DB): $(TAXDUMP_SOURCE) | $(LASTDIR)
	cp $(TAXDUMP_SOURCE)/n??es.dmp $(LASTDIR)/

$(HITIDMAP): $(FAA) | $(LASTDIR)
	@echo "==Building map from hit ids to descriptions"
	grep "^>" $< | perl -pe 's/>(\S+)\s+(.*)$$/\1\t\2/' > $@
