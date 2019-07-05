##########################
#
# download_refseq.make
#
# Retrieves latest nonredundant protein database from refseq.
#
# to update an older version set the relese number with, eg: REL=75
#  I don't think this can download an old version, but you can use it to 
#  re-run the post-download processing.
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
MAX_ATTEMPTS=10

# where to put it
SEQDB_ROOT?=./seqdbs
BUILD_ROOT?=$(SEQDB_ROOT)/RefSeq
RSDIR:=$(BUILD_ROOT)/$(REL)

# compile for lastal?
BUILD_LASTDB:=True
LASTDB_ROOT?=$(SEQDB_ROOT)
LASTDB_DIR:=$(LASTDB_ROOT)/RefSeq/$(REL)
LASTDBCHUNK?=
LASTDBTHREADS?=10

# compile for Diamond?
BUILD_DMNDDB:=True
DMNDDB_FILE:=$(LASTDB_DIR)/$(FAA_NAME).dmndb
DMNDDBTHREADS?=10

# Most folks won't need to edit below this line

# command line options for running lastdb
ifeq ($(LASTDBCHUNK),)
    LASTDBCHUNK_OPTION:=
else
    LASTDBCHUNK_OPTION:= -s $(LASTDB_CHUNK)
endif
ifeq ($(LASTDBTHREADS),)
    LASTDBTHREADS_OPTION:=
else
    LASTDBTHREADS_OPTION:= -P $(LASTDBTHREADS)
endif

# Define the layout of the build directory
MDDIR:=$(RSDIR)/metadata
TAXDUMP_SOURCE:=$(RSDIR)/taxdump

# Final protein outputs
FAA_NAME:=RefSeq-$(REL).AllProteins.faa
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

##
# Build the arguments for all
ifeq ($(BUILD_LASTDB),False)
	TAXDUMP=$(TAXDUMP_SOURCE)
	TARGETS=fasta $(ACCTAXMAP) $(TAXDUMP)
else
	TAXDUMP=$(TAXDUMP_DB)
	TARGETS=lastdb $(ACCTAXMAPDB) $(TAXDUMP)
endif

ALL_TARGETS=report $(MDDIR) $(TARGETS)
all: $(ALL_TARGETS)
	@echo all deps were $^

lastdb: $(LASTFILE) $(HITIDMAP) $(FAA).stats
fasta: $(FAA).stats

report:
	@echo RefSeq release number is: $(REL)
	@echo Downloading database in: $(RSDIR)
	@echo "Output fasta is $(FAA)"
	@echo "BUILD_LASTDB is $(BUILD_LASTDB)"
	@if [ "$(BUILD_LASTDB)" != "False" ]; then echo Final database written to $(LASTDB_DIR); else echo "Lastdb formatting will be skipped"; fi
	@echo "target list is $(TARGETS)"

$(LASTDIR):
	mkdir -p $(LASTDIR)

$(LASTFILE): $(FAA) | $(LASTDIR)
	@echo "==Formating last: $@"
	#lastdb8 -v -c -p $(LASTDBCHUNK_OPTION) $(LASTDBTHREADS_OPTION) $(LASTP) $(FAA)
	lastdb -v -c -p $(LASTDBCHUNK_OPTION) $(LASTDBTHREADS_OPTION) $(LASTP) $(FAA)

%.stats: %
	cat $^ | prinseq-lite.pl -fasta stdin -aa -stats_len -stats_info > $@

$(FAA): $(RSDIR)/complete/.download.complete.aa
	@echo "==Compiling $@ from gz archives"
	@echo "... and masking low complexity with tantan"
	for FILE in $(RSDIR)/complete/complete.*.protein.gpff.gz; do gunzip -c $$FILE; done | get_sequences_from_gb.py -F fasta -r | tantan -p > $@

$(RSDIR)/complete:
	mkdir -p $@

# Dowload directory listing tto file so we can make sure we got everything
$(RSDIR)/complete/complete.ftp.ls: | $(RSDIR)/complete
	curl -s -l ${FTP_ROOT}/complete/ | grep "protein\.gpff" > $@

# This while loop tries up to (MAX_ATTEMPTS) times to download all files
$(RSDIR)/complete/.download.complete.aa: $(RSDIR)/complete/complete.ftp.ls
	@echo "==Dowloading complete RefSeq proteins"
	ATTEMPTS=0; ERF=$(RSDIR)/.tmp.errs; CTF=$(RSDIR)/.tmp.cnts; for F in $$ERF $$CTF; do rm -f $$F; touch $$F; done; echo "start" > $$ERF; while [ -s $$ERF ]; do echo starting attempt $$ATTEMPTS; for F in $$ERF $$CTF; do rm -f $$F; touch $$F; done; if [ "$$ATTEMPTS" -gt "0" ]; then echo Retrying RefSeq protein download; fi; if [ "$$ATTEMPTS" -gt "10" ]; then    echo ERROR: Exceeded 10 tries when downloading RefSeq; else cat $(RSDIR)/complete/complete.ftp.ls | while read CPGZ; do  CPGZ_PATH=$(RSDIR)/complete/$${CPGZ};  if [ ! -s $$CPGZ_PATH ]; then   echo $$CPGZ_PATH>>$$CTF;   curl -s ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/$${CPGZ} > $${CPGZ_PATH};  if [ ! -s $$CPGZ_PATH ]; then   echo $$CPGZ_PATH;    echo $$CPGZ_PATH>$$ERF;  fi; fi; done; echo downloaded `grep -c . $$CTF` files with `grep -c . $$ERF` errors; echo finished attempt $$ATTEMPTS; let ATTEMPTS=ATTEMPTS+1; fi; done
	touch $@

$(ACCPREFP): $(RSDIR)/complete/.download.complete.aa
	# For some multispecies entries, you'll get multiple lines in the tax map
	gunzip -c $(RSDIR)/complete/complete.*.protein.gpff.gz | perl -ne 'if (m/^ACCESSION\s+(\S+)\b/) { $$acc=$$1; } elsif (m/db_xref="taxon:(\d+)"/) { print "$$acc\t$$1\n"; }' > $@

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
