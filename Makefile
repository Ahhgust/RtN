CC=g++

SEQLIB_DIR= /usr/local/src/SeqLib
#CFLAGS= -I/usr/local/src/SeqLib/SeqLib/htslib  -I/usr/local/src/SeqLib/SeqLib -L/usr/local/src/SeqLib/SeqLib/htslib  -L/usr/local/src/SeqLib/SeqLib/src -L/usr/local/src/SeqLib/SeqLib/bwa

CFLAGS= -I${SEQLIB_DIR}/SeqLib/htslib  -I${SEQLIB_DIR}/SeqLib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa
LIBS=-lhts -lseqlib -lhts -lbwa -lz -lpthread

all: rtn.c
	${CC} ${CFLAGS} -Wall -o rtn rtn.c ${LIBS}



