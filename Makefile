CC=g++

SEQLIB_DIR=/usr/local/src/SeqLib

SYS := $(shell ${CC} -dumpmachine)

#CFLAGS=-Wall -std=c++11 -DDEBUG=0 -march=native -Ofast
CFLAGS= -Wall -std=c++11 -static -static-libstdc++ -static-libgcc -DDEBUG=0 -Ofast

#INCLUDES=-Wl,-rpath,${PWD}/SeqLib/htslib -I${SEQLIB_DIR}/SeqLib/htslib  -I${SEQLIB_DIR}/SeqLib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa -L${SEQLIB_DIR}/SeqLib/htslib
INCLUDES=-I${SEQLIB_DIR}/SeqLib -I${SEQLIB_DIR}/SeqLib/htslib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa -L${SEQLIB_DIR}/SeqLib/htslib

#LIBS=-lhts -lseqlib -lhts -lbwa -lz -lpthread
LIBS=-lseqlib -lbwa -lhts -llzma -lz -lbz2 -lpthread -ljsoncpp -ldl

all: rtn.c
	${CC} ${CFLAGS} ${INCLUDES} -Wall -o rtn rtn.c ${LIBS}

f2f: fasta2fastq.cpp
	${CC} ${CFLAGS} ${INCLUDES} -Wall -o f2f fasta2fastq.cpp ${LIBS}

