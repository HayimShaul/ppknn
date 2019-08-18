CCFLAGS = -g --std=c++11 -Wall 

NTLINCDIR = ../../ntl-10.5.0-multithread/include
NTLLIBDIR = ../../ntl-10.5.0-multithread/src

FHEINCDIR = ../../HElib-multithread/src
FHELIBDIR = ../../HElib-multithread/src

LIPHEINCDIR = ../../liphe/include
LIPHELIBDIR = ../../liphe/src

JSONDIR = -I/home/hayim/lib/json/src -I../json/src

LIBS = -L$(LIPHELIBDIR) -lliphe -lHelibLiphe  $(FHELIBDIR)/fhe.a $(NTLLIBDIR)/ntl.a  -lgmp  -lpthread
INCS = $(JSONDIR) -I$(NTLINCDIR) -I$(FHEINCDIR) -I$(LIPHEINCDIR)

#all: test_zp test_helib
all: test_zp test_helib test_kish_vs_k

#all: test min fast_min min2 min3

test_folding: test_folding.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_kish_vs_k: test_kish_vs_k.o get_percentile.o point.o configuration.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_zp: test_zp.o get_percentile.o point.o configuration.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_helib: test_helib.o get_percentile.o point.o configuration.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)


%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o test_zp test_helib build_sqrt_polynomial

