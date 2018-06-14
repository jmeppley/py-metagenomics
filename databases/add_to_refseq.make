# Define some basic locations 
# Edit these or pass in your own values when calling make
FTP_ROOT=ftp://ftp.ncbi.nlm.nih.gov/refseq/release
REL?=$(shell curl $(FTP_ROOT)/RELEASE_NUMBER)
REL:=$(REL)

ROOT_DIR?=.
BUILD_ROOT?=$(ROOT_DIR)/RefSeq
RSDIR:=$(BUILD_ROOT)/$(REL)

LASTDB_ROOT?=$(ROOT_DIR)
LASTDB_DIR:=$(LASTDB_ROOT)/RefSeq/$(REL)
LASTDBCHUNK?=
LASTDBTHREADS?=

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

ADDITIONS_SOURCE:=$(BUILD_ROOT)/additions
ADDITIONS_FAA:=$(ADDITIONS_SOURCE)/additions.protein.fasta
ADDITIONS_TAXIDS:=$(ADDITIONS_SOURCE)/acc.to.taxid.protein.additions
# If filter file is not empty, listed taxids will be removed from additions
ADDITIONS_FILTER:=$(ADDITIONS_SOURCE)/taxids.in.RefSeq.$(REL)
ADDITIONS_KOMAP:=$(ADDITIONS_SOURCE)/acc.to.ko.protein.additions

# Most folks won't need to edit below this line

# Define the layout of the build directory
MDDIR:=$(RSDIR)/metadata
COMPLETEFAA=$(RSDIR)/RefSeq-$(REL).AllProteins.faa

# Final protein outputs
FAA_NAME:=RefSeq-$(REL).AllProteinsPlus.faa
FAA:=$(RSDIR)/$(FAA_NAME)

# LAST database will be packaged in it's own directory
LASTDIR=$(LASTDB_DIR)/$(FAA_NAME).ldb
LASTP=$(LASTDIR)/lastdb
LASTFILE=$(LASTP).prj

ADDFAA=$(RSDIR)/additions.protein.fasta

ACCPREFREL=acc.to.taxid
ACCPREF=$(RSDIR)/$(ACCPREFREL)
ACCPREFP=$(ACCPREF).protein
ACCMAPP=$(ACCPREFP)
ADDACCMAPP=$(ACCPREFP).additions
PLUSACCMAPP=$(ACCPREFP).plus
ACCTAXMAPDB:=$(LASTP).tax
HITIDMAP:=$(LASTP).ids

##
# Build the arguments for all
ALL_TARGETS:=lastdb $(ACCTAXMAPDB)

all: report $(MDDIR) $(ALL_TARGETS)

lastdb: $(LASTFILE) $(HITIDMAP) $(FAA).stats
fasta: $(FAA).stats

report:
	@echo RefSeq release number is: $(REL)
	@echo Building database in: $(RSDIR)
	@echo "Output fasta is $(FAA)"
	@echo "all target list is $(ALL_TARGETS)"
	@echo Final database written to $(LASTDB_ROOT)
	@echo Adding sequences from $(ADDITIONS_SOURCE)

$(LASTDIR):
	mkdir -p $(LASTDIR)

$(LASTFILE): $(FAA) | $(LASTDIR)
	@echo "==Formating last: $@"
<<<<<<< HEAD
	lastdb8 -v -c -p $(LASTDBCHUNK_OPTION) $(LASTP) $(FAA)
=======
    lastdb -v -c -p $(LASTDBCHUNK_OPTION) $(LASTDBTHREADS_OPTION) $(LASTP) $(FAA)
>>>>>>> b77cd14f8dc993c4f34a2def576482e3c6d3bf54

$(FAA): $(ADDFAA) $(COMPLETEFAA)
	cat $^ > $@

$(ADDFAA): $(ADDACCMAPP) $(ADDITIONS_FAA)
	@echo "==Copying records from $@ that are included in $(ADDACCMAPP)"
	@echo "==Also masking low complexity with tantan"
	screen_list.py -a -k -C 0 $(ADDITIONS_FAA) -l $(ADDACCMAPP) | tantan -p > $@

$(ADDITIONS_FILTER): $(ADDITIONS_TAXIDS) $(MDDIR)/release$(REL).taxon.new
	#touch $(ADDITIONS_FILTER)
	cut -f 2 $(ADDITIONS_TAXIDS) | uniq | screen_table.py -l $(MDDIR)/release$(REL).taxon.new -k > $@

$(ADDACCMAPP): $(ADDITIONS_TAXIDS) $(ADDITIONS_FILTER)
	@echo "==Importing taxid map for additions"
	if [ -s $(ADDITIONS_FILTER) ]; then screen_table.py $(ADDITIONS_TAXIDS) -l $(ADDITIONS_FILTER) -c 1 -o $@; else cp $< $@; fi

$(PLUSACCMAPP): $(ADDACCMAPP) $(ACCMAPP)
	cat $(ADDACCMAPP) $(ACCMAPP) > $@

$(ACCTAXMAPDB): $(PLUSACCMAPP) | $(LASTDIR)
	cp $< $@

$(HITIDMAP): $(FAA) | $(LASTDIR)
	@echo "==Building map from hit ids to descriptions"
	grep "^>" $< | perl -pe 's/>(\S+)\s+(.*)$$/\1\t\2/' > $@

%.stats: %
	cat $^ | prinseq-lite.pl -fasta stdin -aa -stats_len -stats_info > $@

