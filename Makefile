CC=g++

SEQLIB_DIR=$(shell pwd)

SYS := $(shell ${CC} -dumpmachine)

CFLAGS=-Wall -std=c++11 -DDEBUG=0 -march=native -Ofast 

INCLUDES=-I${SEQLIB_DIR}/SeqLib -I${SEQLIB_DIR}/SeqLib/htslib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa -L${SEQLIB_DIR}/SeqLib/htslib 

#LIBS=-lhts -lseqlib -lhts -lbwa -lz -lpthread
LIBS=-lseqlib -lbwa -lhts -llzma -lz -lbz2 -lpthread -ljsoncpp 

all: rtn.c
	${CC} ${CFLAGS} -Wl,-rpath,${SEQLIB_DIR}/SeqLib/htslib  ${INCLUDES}  -Wall -o rtn rtn.c ${LIBS}

static: rtn.c
	${CC} ${CFLAGS} -static -Wl,-rpath,${SEQLIB_DIR}/SeqLib/htslib  ${INCLUDES}  -Wall -o rtn rtn.c ${LIBS}

f2f: fasta2fastq.cpp
	${CC} ${CFLAGS} ${INCLUDES} -Wall -o f2f fasta2fastq.cpp ${LIBS}

