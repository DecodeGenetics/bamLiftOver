# Use version 4.8.2 of g++
CXX=g++
CC=$(CXX)

CXXFLAGS+=-std=c++11

# Include SeqAnHTS (provided as submodule)
SEQAN_LIB=./SeqAnHTS/include
#SeqAnHTS requires htslib (available here:https://github.com/samtools/htslib)
HTS_LIB=./htslib/

CXXFLAGS+=-I$(SEQAN_LIB) -I$(HTS_LIB)/include -L$(HTS_LIB)/lib -Wl,-rpath,$(HTS_LIB)/lib -DSEQAN_HAS_ZLIB=1
CXXFLAGS+=-I.
LDLIBS=-lz -lpthread -lhts
DATE=on $(shell git log --pretty=format:"%cd" --date=iso | cut -f 1,2 -d " " | head -n 1)
CXXFLAGS+=-DDATE=\""$(DATE)"\"

# Enable warnings
CXXFLAGS+=-W -Wall -Wextra -Wno-long-long -pedantic -Wno-variadic-macros -Wno-unused-result -Wfatal-errors

# RELEASE build
CXXFLAGS+=-O3 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0


all: bamLiftOver

clean:
	rm -f *.o bamLiftOver
