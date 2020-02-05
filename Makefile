CC=g++

# /usr/local/src/SeqLib
SEQLIB_DIR=.


SYS := $(shell ${CC} -dumpmachine)

CFLAGS=-Wall -std=c++11 -DDEBUG=0 -march=native -Ofast


#INCLUDES=  -I${SEQLIB_DIR}/SeqLib/htslib  -I${SEQLIB_DIR}/SeqLib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa -L${SEQLIB_DIR}/SeqLib/htslib
INCLUDES=-Wl,-rpath,${PWD}/SeqLib/htslib -I${SEQLIB_DIR}/SeqLib/htslib  -I${SEQLIB_DIR}/SeqLib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa -L${SEQLIB_DIR}/SeqLib/htslib

LIBS=-lhts -lseqlib -lhts -lbwa -lz -lpthread

all: rtn.c
	${CC} ${CFLAGS} ${INCLUDES} -Wall -o rtn rtn.c ${LIBS}

f2f: fasta2fastq.cpp
	${CC} ${CFLAGS} ${INCLUDES} -Wall -o f2f fasta2fastq.cpp ${LIBS}

