DATE=$(shell date +%Y%m%d)
KGDIR=$(DATE)

KOKEG=$(KGDIR)/ko00001.keg
KOKO=$(KGDIR)/ko/ko
GENES=$(KGDIR)/fasta/genes.pep.gz
KEGGGENES=KeggGene.pep.$(DATE).faa
BLASTDBDIR=/common/data/KEGG
LASTDBDIR=/common/ldb/KEGG
LASTDBCHUNK=100G
KEGGBLAST=$(BLASTDBDIR)/KeggGene.pep.$(DATE).faa
KEGGLAST=$(LASTDBDIR)/KeggGene.pep.$(DATE).ldb
KOMAP=$(LASTDBDIR)/KeggGene.pep.$(DATE).genes_ko.list
DESCMAP=$(LASTDBDIR)/KeggGene.pep.$(DATE).genes_desc.list
KEGGLAST_PRJ=$(KEGGLAST).prj
LINKS=$(KGDIR)/links
GENOME=$(KGDIR)/genome
BRITE=$(KGDIR)/brite
TAX=$(KGDIR)/taxonomy

WGET=wget -c --user=delong@mit.edu --password=mbari1
FGENES=ftp://ftp.bioinformatics.jp/kegg/genes/MD5.genes ftp://ftp.bioinformatics.jp/kegg/genes/fasta/genes.pep.gz ftp://ftp.bioinformatics.jp/kegg/genes/README.genes
FLINKS=ftp://ftp.bioinformatics.jp/kegg/genes/links/*gz
FGENOME=ftp://ftp.bioinformatics.jp/kegg/genes/genome.tar.gz
FKO=ftp://ftp.bioinformatics.jp/kegg/genes/ko.tar.gz
FTAX=ftp://ftp.bioinformatics.jp/kegg/genes/taxonomy
FBRITE=ftp://ftp.bioinformatics.jp/kegg/brite/*.tar.gz ftp://ftp.bioinformatics.jp/kegg/brite/brite* ftp://ftp.bioinformatics.jp/kegg/brite/*brite

all: $(KEGGBLAST) noblast
noblast: $(KOKO) $(KOKEG) $(LINKS) $(GENOME) $(TAX) $(KEGGLAST_PRJ) maps
maps: $(KOMAP) $(DESCMAP)

$(KOMAP): $(LINKS)
	cp $(LINKS)/genes_ko.list $@

$(DESCMAP): $(KEGGGENES)
	grep ">" $< | perl -pe 's/^>(\S+)\s+(\S.*)/\1\t\2/' > $@

$(KEGGLAST_PRJ): $(KEGGGENES)
	lastdb -v -c -p -s $(LASTDBCHUNK) $(KEGGLAST) $(KEGGGENES)

$(KEGGBLAST): $(KEGGGENES)
	edlformatdb $(KEGGGENES)

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
