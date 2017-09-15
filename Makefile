CCFLAGS = -g --std=c++11

NTLINCDIR = -I/home/hayim/academic/fhe/ntl-9.6.2/include
NTLLIBDIR = -L/home/hayim/academic/fhe/ntl-9.6.2/src

FHEINCDIR = -I/home/hayim/academic/fhe/HElib-master/src
FHELIBDIR = -L/home/hayim/academic/fhe/HElib-master/src

HEADSUPINCDIR = -I../liphe/include
HEADSUPLIBDIR = -L../liphe/src

LIBS = $(HEADSUPLIBDIR) -lliphe $(FHELIBDIR) -lfhe $(NTLLIBDIR) -lntl
INCS = $(NTLINCDIR) $(FHEINCDIR) $(HEADSUPINCDIR)

all: test_zp 

#all: test min fast_min min2 min3

test_zp: test_zp.o get_percentile.o mem.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)

test_helib: test_helib.o one_mean.o
	g++ $(LDFLAGS) -o $@ $^ $(LIBS)


%.o: %.cc
	g++ $(CCFLAGS) -c  $(INCS) $<

clean:
	rm -f *.o test_zp test_helib

