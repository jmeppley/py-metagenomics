#DATE?=$(shell date +%Y%m%d)
#KGVER:=$(DATE)
FTP_ROOT:=ftp://ftp.bioinformatics.jp/kegg
KGVER:=$(shell curl --user $(KEGG_USER):$(KEGG_PASSWORD) $(FTP_ROOT)/RELEASE | head -1 | perl -pe 's/[^0-9]//g')
KEGG_USER?=
KEGG_PASSWORD?=

DB_SCRIPT_DIR?=.
BUILD_ROOT?=./KEGG
KGDIR:=$(BUILD_ROOT)/$(KGVER)

BUILD_LASTDB:=True
LASTDB_ROOT?=./lastdb
LASTDBCHUNK?=

ifeq ($(LASTDBCHUNK),)
	LASTDBCHUNK_OPTION:=
else
	LASTDBCHUNK_OPTION:= -s $(LASTDB_CHUNK)
endif

KOKEG=$(KGDIR)/ko00001.keg
KOKO=$(KGDIR)/ko/ko
PRKGENES=$(KGDIR)/fasta/prokaryotes.pep.gz
EUKGENES=$(KGDIR)/fasta/eukaryotes.pep.gz
GENESMETA=README.genes
KEGGGENES=$(KGDIR)/KeggGene.pep.$(KGVER).faa
LINKS=$(KGDIR)/links
GENOME=$(KGDIR)/genome
BRITE=$(KGDIR)/brite
TAX=$(KGDIR)/taxonomy

ifeq ($(BUILD_LASTDB),False)
	LASTDB_DIR:=$(KGDIR)
	KOMAP:=$(LASTDB_DIR)/KeggGene.pep.$(KGVER).kos
	HITIDMAP:=$(LASTDB_DIR)/KeggGene.pep.$(KGVER).ids
else
	LASTDB_DIR:=$(LASTDB_ROOT)/KEGG/KeggGene.pep.$(KGVER)
	LASTP:=$(LASTDB_DIR)/lastdb
	LASTFILE=$(LASTP).prj
	KOMAP=$(LASTP).kos
	HITIDMAP:=$(LASTP).ids
endif

WGET:=wget -c --user=$(KEGG_USER) --password=$(KEGG_PASSWORD)
FGENESPRK=$(FTP_ROOT)/genes/fasta/prokaryotes.pep.gz
FGENESEUK=$(FTP_ROOT)/genes/fasta/eukaryotes.pep.gz
FGENESMETA=$(FTP_ROOT)/genes/MD5.genes $(FTP_ROOT)/genes/README.genes $(FTP_ROOT)/RELEASE
FLINKS=$(FTP_ROOT)/genes/links/*gz
FGENOME=$(FTP_ROOT)/genes/genome.tar.gz
FKO=$(FTP_ROOT)/genes/ko.tar.gz
FTAX=$(FTP_ROOT)/genes/taxonomy
FBRITE=$(FTP_ROOT)/brite/*.tar.gz $(FTP_ROOT)/brite/brite* $(FTP_ROOT)/brite/*brite

ALL_TARGETS:=$(KOKO) $(KOKEG) $(LINKS) $(GENOME) $(TAX)
ifeq ($(KEGG_USER),)
	ALL_TARGETS:=nouser
else
	ifeq ($(BUILD_LASTDB),False)
		ALL_TARGETS:=$(ALL_TARGETS) $(KEGGGENES)
	else
		ALL_TARGETS:=$(ALL_TARGETS) $(LASTFILE) maps
	endif
endif

all: version $(ALL_TARGETS)

nouser:
	@echo You must set the KEGG_USER and KEGG_PASSWORD variables!

version:
	@echo KEGG version is $(KGVER)
	@echo KOMAP is $(KOMAP)
	@echo $(ALL_TARGETS)

maps: $(KOMAP) $(HITIDMAP)
links: $(LINKS)

$(KOMAP): $(LINKS) | $(LASTDB_DIR)
	cp $(LINKS)/genes_ko.list $@

$(HITIDMAP): $(KEGGGENES) | $(LASTDB_DIR)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

$(LASTFILE): $(KEGGGENES) | $(LASTDB_DIR)
	lastdb -v -c -p $(LASTDBCHUNK_OPTION) $(LASTP) $(KEGGGENES)

$(LASTDB_DIR):
	mkdir -p $(LASTDB_DIR)

$(KEGGGENES): $(GENESMETA) $(PRKGENES) $(EUKGENES)
	rm -f $(KEGGGENES)
	gunzip -c $(PRKGENES) $(EUKGENES) | tantan -p > $(KEGGGENES)

$(KGDIR):
	mkdir -p $(KGDIR)

$(KGDIR)/fasta:
	mkdir -p $(KGDIR)/fasta

$(GENESMETA): | $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FGENESMETA)

$(PRKGENES): | $(KGDIR)/fasta
	cd $(KGDIR)/fasta && $(WGET) $(FGENESPRK)

$(EUKGENES): | $(KGDIR)/fasta
	cd $(KGDIR)/fasta && $(WGET) $(FGENESEUK)
	
$(KOKO): | $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FKO)
	cd $(KGDIR) && tar -zxvf ko.tar.gz
	rm $(KGDIR)/ko.tar.gz

$(LINKS): | $(KGDIR)
	mkdir -p $(LINKS)
	cd $(LINKS) && $(WGET) $(FLINKS)
	cd $(LINKS) && gunzip *gz

$(GENOME): | $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FGENOME)
	cd $(KGDIR) && tar -zxvf genome.tar.gz
	rm $(KGDIR)/genome.tar.gz

$(TAX): | $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FTAX)

$(KOKEG): $(BRITE)
	cd $(BRITE) && tar -zxvf ko.tar.gz
	rm -f $(KOKEG)
	cp $(BRITE)/ko/ko00001.keg $(KOKEG)

$(BRITE):
	mkdir -p $(BRITE)
	cd $(BRITE) && $(WGET) $(FBRITE)
