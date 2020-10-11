HEAAN=./HEAAN/HEAAN
HEAAN_LIBS=-L$(HEAAN)/lib -lHEAAN -Wl,-rpath=$(HEAAN)/lib -lntl -Wl,-rpath=/usr/local/lib -lgmp -lm
HEAAN_INCLUDE=-I$(HEAAN)/src

FRNSHEAAN=./FullRNS-HEAAN
FRNSHEAAN_LIBS=-L$(FRNSHEAAN)/lib -lFRNSHEAAN -Wl,-rpath=$(FRNSHEAAN)/lib -lntl -Wl,-rpath=/usr/local/lib -lgmp -lm
FRNSHEAAN_INCLUDE=-I$(FRNSHEAAN)/src

PALISADE=/usr/local/include/palisade
PALISADE_INCLUDE=-I$(PALISADE) -I$(PALISADE)/core -I$(PALISADE)/pke
PALISADE_LIBS=-L/usr/local/lib -lPALISADEcore -lPALISADEpke -Wl,-rpath=/usr/local/lib -lntl

SEAL_INCLUDE=-I/usr/local/include/SEAL-3.5
SEAL_LIBS=-L/usr/local/lib -lseal-3.5 -lz -Wl,-rpath=/usr/local/lib -lntl

HELIB_INCLUDE=-I/usr/local/include
HELIB_LIBS=-Wl,-rpath=/usr/local/lib -L/usr/local/lib -lhelib -lntl -lgmp

all: attack palisade_attack seal_attack helib_attack

attack: attack.cpp eval.h eval.cpp ntl_utils.h ntl_utils.cpp
	g++ attack.cpp eval.cpp ntl_utils.cpp -o attack  -std=c++11 -O3 -pthread $(HEAAN_LIBS) $(HEAAN_INCLUDE)

rns_attack: rns_attack.cpp eval.h eval.cpp ntl_utils.h ntl_utils.cpp
	g++ rns_attack.cpp eval.cpp ntl_utils.cpp -o rns_attack  -std=c++11 -fopenmp -O3 -pthread $(FRNSHEAAN_LIBS) $(FRNSHEAAN_INCLUDE)

palisade_attack: palisade_attack.cpp palisade_utils.cpp eval.h eval.cpp
	g++ palisade_attack.cpp palisade_utils.cpp eval.cpp -o palisade_attack  -std=c++11 -O3 $(PALISADE_LIBS) $(PALISADE_INCLUDE)

seal_attack: seal_attack.cpp eval.h eval.cpp seal_utils.h seal_utils.cpp
	g++ seal_attack.cpp seal_utils.cpp eval.cpp -o seal_attack  -std=c++17 -fopenmp -O3 -pthread $(SEAL_LIBS) $(SEAL_INCLUDE)

helib_attack: helib_attack.cpp helib_utils.h helib_utils.cpp ntl_utils.h ntl_utils.cpp eval.h eval.cpp
	g++ helib_attack.cpp helib_utils.cpp ntl_utils.cpp eval.cpp -o helib_attack -std=c++14 -O3 -pthread $(HELIB_LIBS) $(HELIB_INCLUDE) 

clean:
	rm -f attack palisade_attack seal_attack helib_attack

.PHONY: all clean

