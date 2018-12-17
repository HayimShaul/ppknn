CCFLAGS = -g --std=c++11 -Wall

NTLINCDIR = -I../../ntl-10.5.0-multithread/include
NTLLIBDIR = ../../ntl-10.5.0-multithread/src

FHEINCDIR = -I../../HElib-multithread/src
FHELIBDIR = ../../HElib-multithread/src

HEADSUPINCDIR = -I../../liphe/include
HEADSUPLIBDIR = -L../../liphe/src

JSONDIR = -I/home/hayim/lib/json/src -I../json/src

LIBS = $(HEADSUPLIBDIR) -lliphe  $(FHELIBDIR)/fhe.a $(NTLLIBDIR)/ntl.a  -lgmp  -lpthread
INCS = $(JSONDIR) $(NTLINCDIR) $(FHEINCDIR) $(HEADSUPINCDIR)

#all: test_zp test_helib
all: test_zp test_helib

#all: test min fast_min min2 min3

test_folding: test_folding.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_zp: test_zp.o get_percentile.o point2d.o mem.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_helib: test_helib.o get_percentile.o point2d.o mem.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)


%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o test_zp test_helib build_sqrt_polynomial

