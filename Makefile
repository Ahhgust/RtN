CC=g++

SEQLIB_DIR= /usr/local/src/SeqLib

SYS := $(shell ${CC} -dumpmachine)

CFLAGS=-Wall -std=c++11 -Ofast -DDEBUG=0 -march=native -fomit-frame-pointer # -pg

INCLUDES= -I${SEQLIB_DIR}/SeqLib/htslib  -I${SEQLIB_DIR}/SeqLib -L${SEQLIB_DIR}/SeqLib/src -L${SEQLIB_DIR}/SeqLib/bwa
LIBS=-lhts -lseqlib -lhts -lbwa -lz -lpthread

all: rtn.c
	${CC} ${CFLAGS} ${INCLUDES} -Wall -o rtn rtn.c ${LIBS}



