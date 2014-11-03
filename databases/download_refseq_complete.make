# Define some basic locations 
# Edit these or pass in your own values when calling make
FTP_ROOT=ftp://ftp.ncbi.nlm.nih.gov/refseq/release
REL?=$(shell curl $(FTP_ROOT)/RELEASE_NUMBER)
REL:=$(REL)

DB_SCRIPT_DIR?=.
BUILD_ROOT?=./RefSeq
RSDIR:=$(BUILD_ROOT)/RefSeq-$(REL)

BUILD_LASTDB:=False
LASTDB_ROOT?=/minilims/galaxy-data/tool-data/sequencedbs/lastdb/RefSeq/$(REL)
LASTDBCHUNK?=100G

ADD_CUSTOM_SEQS:=False
ADDITIONS_SOURCE:=$(BUILD_ROOT)/additions
ADDITIONS_FAA:=$(ADDITIONS_SOURCE)/additions.protein.fasta
ADDITIONS_TAXIDS:=$(ADDITIONS_SOURCE)/acc.to.taxid.protein.additions
# If filter file is not empty, listed taxids will be removed from additions
ADDITIONS_FILTER:=$(ADDITIONS_SOURCE)/taxids.in.RefSeq.$(REL)

BUILD_KO_MAP:=False
KEGG_ROOT?=./KEGG
ifeq ($(BUILD_KO_MAP),False)
	KEGGDATE=LATEST	# Placeholder
else
	# finds the most recent DATE folder (can be overridden in `make` command)
	KEGGDATE:=$(shell ls -1rt $(KEGG_ROOT) | grep "20[0-9][0-9][0-9][0-9][0-9][0-9]$$" | tail -1)
endif

# Most folks won't need to edit below this line

# Define the layout of the build directory
MDDIR:=$(RSDIR)/metadata
TAXDUMP:=$(RSDIR)/taxdump

COMPLETEFAA=$(RSDIR)/complete.protein.fasta

# Final protein outputs
ifeq ($(ADD_CUSTOM_SEQS),False)
	FAA_NAME:=RefSeq-$(REL).AllProteins.faa
	FAA_PREREQS=$(COMPLETEFAA)
else
	FAA_NAME:=RefSeq-$(REL).AllProteinsPlus.faa
	FAA_PREREQS=$(ADDFAA) $(COMPLETEFAA)
endif
FAA:=$(RSDIR)/$(FAA_NAME)
# LAST database will be packaged in it's own directory
LASTDIR=$(LASTDB_ROOT)/$(FAA_NAME).ldb
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
ifeq ($(ADD_CUSTOM_SEQS),False)
	ACCTAXMAP:=$(ACCMAPP)
else
	ACCTAXMAP:=$(PLUSACCMAPP)
endif

# KEGG locations
KOMAP=$(RSDIR)/acc.to.ko.protein
KOMAP_PLUS=$(RSDIR)/acc.to.ko.protein.plus
KOMAP_ADD=$(RSDIR)/KEGG/acc.to.ko.additions
ifeq ($(ADD_CUSTOM_SEQS),False)
	KOMAP_ALL:=$(KOMAP)
else
	KOMAP_ALL:=$(KOMAP_PLUS)
endif
KOMAP_DB:=$(LASTDIR)/lastdb.kos
ifeq ($(BUILD_LASTDB),False)
	KOMAP_DEP:=$(KOMAP_ALL)
else
	KOMAP_DEP:=$(KOMAP_DB)
endif

KEGGLINKDIR:=$(KEGG_ROOT)/$(KEGGDATE)/links
KEGGGENE_GI_MAP=$(KEGGLINKDIR)/genes_ncbi-gi.list
KEGGGENE_KO_MAP=$(KEGGLINKDIR)/genes_ko.list
KOMAPSCRIPT=$(DB_SCRIPT_DIR)/buildAccKOMapping.py
TAXMAPSCRIPT=$(DB_SCRIPT_DIR)/buildRefSeqAccToTaxidMap.py

##
# Build the arguments for all
ifeq ($(BUILD_LASTDB),False)
	ALL_TARGETS:=fasta $(ACCTAXMAP)
else
	ALL_TARGETS:=lastdb  $(ACCTAXMAPDB)
endif
ALL_TARGETS:=$(ALL_TARGETS)
ifneq ($(BUILD_KO_MAP),False)
	ALL_TARGETS:=$(ALL_TARGETS) keggmap
endif

all: report $(ALL_TARGETS)

lastdb: $(LASTFILE) $(HITIDMAP)
fasta: $(FAA)
keggmap: $(KOMAP_DEP) 

report:
	@echo RefSeq release number is: $(REL)
	@echo Building database in: $(RSDIR)
	@echo "Output fasta is $(FAA)"
	@echo "BUILD_LASTDB is $(BUILD_LASTDB)"
	@echo "all target list is $(ALL_TARGETS)"
	@if [ "$(BUILD_LASTDB)" != "False" ]; then echo Final database written to $(LASTDB_ROOT); else echo "Lastdb formatting will be skipped"; fi
	@if [ "$(ADD_CUSTOM_SEQS)" != "False" ]; then echo Adding sequences from  $(ADDITIONS_SOURCE); fi
	@if [ "$(BUILD_KO_MAP)" != "False" ]; then echo Building KO map from link files in $(KEGG_ROOT); fi

$(LASTFILE): $(FAA)
	mkdir -p $(LASTDIR)
	@echo "==Formating last: $@"
	lastdb -v -c -p -s $(LASTDBCHUNK) $(LASTP) $(FAA)

$(FAA): $(FAA_PREREQS)
	@echo "==Masking low complexity with tantan"
	tantan -p $^ | perl -ne 'if (m/^>(?!gi\|\d+)(.*)$$/) { if (defined $$n) { $$n++; } else { $$n=10000000000; } print ">gi|$$n|loc|$$1\n"; } else { print; }' > $@
    
$(ADDFAA): $(ADDACCMAPP) $(ADDITIONS_FAA)
	@echo "==Copying $@ from previous installation"
	screen_list.py -a -k -C 0 $(ADDITIONS_FAA) -l $(ADDACCMAPP) -o $@

$(ADDITIONS_FILTER):
	touch $(ADDITIONS_FILTER)

$(ADDACCMAPP): $(ADDITIONS_TAXIDS) $(ADDITIONS_FILTER)
	@echo "==Importing taxid map for additions"
	if [ -s $(ADDITIONS_FILTER) ]; then screen_table.py $^ -l $(ADDITIONS_FILTER) -c 1 -o $@; else cp $< $@; fi
    
%.protein.fasta: %/.download.complete.aa
	@echo "==Compiling $@ from gz archives"
	for FILE in $(RSDIR)/$(*F)/*nonredundant_protein*gpff.gz; do gunzip -c $$FILE; done | getSequencesFromGbk.py -F fasta -r > $@

$(RSDIR)/complete/.download.complete.aa:
	@echo "==Dowloading complete RefSeq proteins"
	mkdir -p $(RSDIR)/complete
	cd $(RSDIR)/complete && wget -c $(FTP_ROOT)/complete/complete.nonredundant_protein.*.protein.gpff.gz
	touch $@

$(ACCMAPP): $(MDDIR) $(TAXDUMP) $(TAXMAPSCRIPT) 
	gunzip -c complete.nonredundant_protein.38.protein.gpff.gz | perl -ne 'if (m/^ACCESSION\s+(\S+)\b/) { $$acc=$$1; } elsif (m/db_xref="taxon:(\d+)"/) { print "$$acc\t$$1\n"; }' | sort > $@
	gunzip -c $(MDDIR)/RefSeq-release$(REL).catalog.gz | $(TAXMAPSCRIPT) $(TAXDUMP) | sort > $@.oldway

$(ACCTAXMAPDB): $(ACCTAXMAP)
	mkdir -p $(LASTDIR)
	cp $< $@

$(MDDIR):
	@echo "==Downloading metadata"
	mkdir -p $(MDDIR)
	cd $(MDDIR) && wget -c $(FTP_ROOT)/release-notes/RefSeq*.txt $(FTP_ROOT)/release-statistics/RefSeq-release*.stats.txt $(FTP_ROOT)/release-catalog/RefSeq-release$(REL).catalog.gz $(FTP_ROOT)/release-catalog/release$(REL)*

$(TAXDUMP):
	@echo "==Downloading taxonomy"
	mkdir -p $(TAXDUMP)
	cd $(TAXDUMP) && wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -zxvf taxdump.tar.gz
	if [ "$(BUILD_LASTDB)" != "False" ]; then mkdir -p $(LASTDIR); cp $(TAXDUMP)/n??es.dmp $(LASTDIR)/; fi

$(KOMAP_DB): $(KOMAP_ALL)
	mkdir -p $(LASTDIR)
	cp $^ $@

$(KOMAP_PLUS): $(KOMAP) $(KOMAP_ADD)
	@echo "==Combine kegg ko map with ko map from additions"
	cat $^ > $@

$(KOMAP_ADD):
	@echo "==Copy provided acc.to.ko map"
	mkdir -p $(RSDIR)/KEGG
	cp $(ADDITIONS_SOURCE)/acc.to.ko.additions $@

$(KOMAP): $(COMPLETEFAA) $(KEGGGENE_KO_MAP) $(KEGGGENE_GI_MAP) $(KOMAPSCRIPT)
	@echo "==Building map from accessions to kos"
	$(KOMAPSCRIPT) -v $(COMPLETEFAA) -l $(KEGGLINKDIR) | sort > $@

$(HITIDMAP): $(FAA)
	mkdir -p $(LASTDIR)
	@echo "==Building map from hit ids to descriptions"
	grep "^>" $^ | perl -pe 's/>(\S+)\s+(.*)$$/\1\t\2/' > $@
