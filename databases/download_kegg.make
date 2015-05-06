DATE?=$(shell date +%Y%m%d)
KGVER:=$(DATE)
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
GENES=$(KGDIR)/fasta/genes.pep.gz
KEGGGENES=$(KGDIR)/KeggGene.pep.$(DATE).faa
LASTDB_DIR:=$(LASTDB_ROOT)/KEGG/KeggGene.pep.$(DATE)
LASTP:=$(LASTDB_DIR)/lastdb
LASTFILE=$(LASTP).prj
KOMAP=$(LASTP).kos
HITIDMAP:=$(LASTP).ids
LINKS=$(KGDIR)/links
GENOME=$(KGDIR)/genome
BRITE=$(KGDIR)/brite
TAX=$(KGDIR)/taxonomy

WGET:=wget -c --user=$(KEGG_USER) --password=$(KEGG_PASSWORD)
FGENES=ftp://ftp.bioinformatics.jp/kegg/genes/MD5.genes ftp://ftp.bioinformatics.jp/kegg/genes/fasta/genes.pep.gz ftp://ftp.bioinformatics.jp/kegg/genes/README.genes
FLINKS=ftp://ftp.bioinformatics.jp/kegg/genes/links/*gz
FGENOME=ftp://ftp.bioinformatics.jp/kegg/genes/genome.tar.gz
FKO=ftp://ftp.bioinformatics.jp/kegg/genes/ko.tar.gz
FTAX=ftp://ftp.bioinformatics.jp/kegg/genes/taxonomy
FBRITE=ftp://ftp.bioinformatics.jp/kegg/brite/*.tar.gz ftp://ftp.bioinformatics.jp/kegg/brite/brite* ftp://ftp.bioinformatics.jp/kegg/brite/*brite

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

all: $(ALL_TARGETS)

nouser:
	@echo You must set the KEGG_USER and KEGG_PASSWORD variables!

maps: $(KOMAP) $(HITIDMAP)
links: $(LINKS)

$(KOMAP): $(LINKS)
	cp $(LINKS)/genes_ko.list $@

$(HITIDMAP): $(KEGGGENES)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

$(LASTFILE): $(KEGGGENES)
	mkdir -p $(LASTDB_DIR)
	lastdb -v -c -p $(LASTDBCHUNK_OPTION) $(LASTP) $(KEGGGENES)

$(KEGGGENES): $(GENES)
	rm -f $(KEGGGENES)
	gunzip -c $(GENES) | tantan -p > $(KEGGGENES)

$(GENES):
	mkdir -p $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FGENES)
	mkdir -p $(KGDIR)/fasta
	mv $(KGDIR)/genes.pep.gz $(GENES)
	
$(KOKO):
	mkdir -p $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FKO)
	cd $(KGDIR) && tar -zxvf ko.tar.gz
	rm $(KGDIR)/ko.tar.gz

$(LINKS):
	mkdir -p $(KGDIR)
	mkdir -p $(LINKS)
	cd $(LINKS) && $(WGET) $(FLINKS)
	cd $(LINKS) && gunzip *gz

$(GENOME):
	mkdir -p $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FGENOME)
	cd $(KGDIR) && tar -zxvf genome.tar.gz
	rm $(KGDIR)/genome.tar.gz

$(TAX):
	mkdir -p $(KGDIR)
	cd $(KGDIR) && $(WGET) $(FTAX)

$(KOKEG): $(BRITE)
	cd $(BRITE) && tar -zxvf ko.tar.gz
	rm -f $(KOKEG)
	cp $(BRITE)/ko/ko00001.keg $(KOKEG)

$(BRITE):
	mkdir -p $(BRITE)
	cd $(BRITE) && $(WGET) $(FBRITE)
