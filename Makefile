CCFLAGS = -g --std=c++11

#NTLINCDIR = -I/home/hayim/academic/fhe/ntl-9.6.2/include
#NTLLIBDIR = -L/home/hayim/academic/fhe/ntl-9.6.2/src
#
#FHEINCDIR = -I/home/hayim/academic/fhe/HElib-master/src
#FHELIBDIR = -L/home/hayim/academic/fhe/HElib-master/src

NTLINCDIR = -I/home/hayim/academic/fhe/ntl-10.5.0/include
NTLLIBDIR = -L/home/hayim/academic/fhe/ntl-10.5.0/src

FHEINCDIR = -I/home/hayim/academic/fhe/HElib/src
FHELIBDIR = -L/home/hayim/academic/fhe/HElib/src

HEADSUPINCDIR = -I../liphe/include
HEADSUPLIBDIR = -L../liphe/src

LIBS = $(HEADSUPLIBDIR) -lliphe $(FHELIBDIR) -lfhe $(NTLLIBDIR) -lntl -lgmp
INCS = $(NTLINCDIR) $(FHEINCDIR) $(HEADSUPINCDIR)

all: test_zp test_helib build_sqrt_polynomial

#all: test min fast_min min2 min3

build_sqrt_polynomial: build_sqrt_polynomial.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_zp: test_zp.o get_percentile.o mem.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_helib: test_helib.o get_percentile.o mem.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)


%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o test_zp test_helib build_sqrt_polynomial

