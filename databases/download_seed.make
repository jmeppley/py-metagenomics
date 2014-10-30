FTPROOT=ftp://ftp.theseed.org
GENOMES_VER:=$(shell curl -l $(FTPROOT)/genomes/ | grep "data\.[0-9]" | sort | tail -1 | cut -d . -f 2)
SDDIR=SEED.$(GENOMES_VER)

#ROLEMAP=$(SDDIR)/maps/subsystems2role
ROLEMAP=$(SDDIR)/maps/role2subsystems
GENEMAP=$(SDDIR)/maps/gene2role
GENOMES=$(SDDIR)/genomes
GENOMESFAA=$(SDDIR)/seed.genomes.tt.faa
GENOMESHITMAP:=$(SDDIR)/maps/seed.genomes.hitd.to.desc
#GENOMESHITMAPFILTERED:=$(SDDIR)/maps/seed.genomes.hitd.to.desc.hits
#SCRIPTDIR:=/opt/edl/scripts
SCRIPTDIR:=~/work/delong/projects/scripts
LASTDBDIR:=/minilims/galaxy-data/tool-data/sequencedbs/seed/$(GENOMES_VER)
REFSEQ:=/minilims/galaxy-data/tool-data/sequencedbs/lastdb/RefSeq/61/RefSeq-61.AllProteinsPlus.faa.ldb
LASTDBCHUNK=100G
GENOMESDB=$(LASTDBDIR)/SEED.genomes
ROLEMAPSCRIPT=$(SCRIPTDIR)/scrapeSeedSubsystems.py

REF2SUBSYSTEMMAP=$(SDDIR)/maps/refseq2subsystem
REFSEQVSSEEDHITS=$(SDDIR)/RefSeq.vs.Seed.lastp
SEEDVSREFSEQHITS=$(REFSEQVSSEEDHITS).inverted

#ROLEMAPURL=$(FTPROOT)/subsystems/subsystems2role.gz
SEEDVIEWERURL=http://seed-viewer.theseed.org/seedviewer.cgi
GENOMESURL=$(FTPROOT)/genomes/data.$(GENOMES_VER)/SEED/*faa
WGET=wget -c 

all: $(ROLEMAP) $(REF2SUBSYSTEMMAP) $(GENOMESDB).prj

#$(ROLEMAP).gz:
#	mkdir -p $(SDDIR)/maps
#	cd $(SDDIR)/maps && $(WGET) $(ROLEMAPURL)
#	touch $(ROLEMAP).gz

#$(ROLEMAP): $(ROLEMAP).gz
#	gunzip -c $< > $@

#$(ROLEMAP): $(ROLEMAPSCRIPT)
#	$(ROLEMAPSCRIPT) -o $@ -u $(SEEDVIEWERURL)

$(ROLEMAP) $(GENEMAP): $(ROLEMAPSCRIPT)
	$(ROLEMAPSCRIPT) -o $(ROLEMAP) -g $(GENEMAP) -u $(SEEDVIEWERURL)

$(SEEDVSREFSEQHITS): $(REFSEQVSSEEDHITS)
	cat $(REFSEQVSSEEDHITS) | perl -ne '@cells=split("\t"); $$tmpref=$$cells[1]; $$cells[1]=$$cells[6]; $$cells[6]=$$tmpref; print join("\t", @cells);' | sort -T /minilims/tmp/ -t $$'\t' -k 7,7 -k 1rn,1 | $(SCRIPTDIR)/filter_blast_m8.py -f last -H 1 > $@

#$(GENOMESHITMAPFILTERED): $(GENOMESHITMAP) $(SEEDVSREFSEQHITS)
#	$(SCRIPTDIR)/screen_table.py -k -l $(SEEDVSREFSEQHITS) -C 1 $(GENOMESHITMAP) -o $@ 
	
$(REF2SUBSYSTEMMAP): $(SEEDVSREFSEQHITS) $(GENOMESHITMAP)
	cat $(SEEDVSREFSEQHITS) | cut -f 2,7 | perl -lane 'print "$$F[1]\t$$F[0]"' | $(SCRIPTDIR)/translateColumn.py -m $(GENOMESHITMAP) -c 2 -D 2 | perl -ne 's/\s+\[[^\[\]]+\]\s+/\n/; s/^\S+\|ref\|([-_A-Z0-9a-z]+)(:?\.\d+)?\|?\S*/\1/; if (m/\t\S+/) { print ; }' > $@

$(REFSEQVSSEEDHITS): $(REFSEQ).prj $(GENOMESFAA)
	$(SCRIPTDIR)/lastWrapper.py -N 8 -f 0 -u 2 -o $@ $(REFSEQ) $(GENOMESFAA)

$(GENOMESHITMAP): $(GENOMESFAA)
	grep ">" $(GENOMESFAA) | perl -pe 's/^>(\S+)\s+(\S.+\S)\s*$$/\1\t\2\n/' > $@

$(GENOMESDB).prj: $(GENOMESFAA)
	mkdir -p $(LASTDBDIR)
	lastdb -c -p -v -s $(LASTDBCHUNK) $(GENOMESDB) $(GENOMESFAA)

$(GENOMESFAA): $(GENOMES)
	cat $(GENOMES)/*faa | tantan -p > $@

$(GENOMES):
	mkdir -p $@
	cd $@ && $(WGET) $(GENOMESURL)
