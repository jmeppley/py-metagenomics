FTP_ROOT=ftp://ftp.ncbi.nlm.nih.gov/refseq/release
REL?=$(shell curl $(FTP_ROOT)/RELEASE_NUMBER)
REL:=$(REL)

BUILD_ROOT?=databases/RefSeq
RSDIR:=$(BUILD_ROOT)/RefSeq-$(REL)

report:
	echo RefSeq release number is: $(REL)
    echo Building database in: $(RSDIR)

