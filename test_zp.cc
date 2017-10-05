#include <zp.h>
#include <helib_number.h>

#include <eq.h>
#include "polynomial.h" // Before it gets into Liphe
#include <cmp.h>
#include <binomial_tournament.h>
#include <first_non_zero.h>
#include <unsigned_word.h>
#include <average.h>

#include "get_percentile_framework.h"
#include "get_percentile.h"
#include "special_polynomials.h"

#define SIMD_FACTOR 2
typedef ZP<SIMD_FACTOR> MyZP;
typedef UnsignedWord<16, MyZP> MyZPBits;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::square_msd_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::sqrt_polynomial;

template<class Number>
std::vector<Polynomial<Number> > SpecialPolynomials<Number>::convert_to_bit;

int main(int argc, char **argv) {
	initialize(argc, argv);
	MyZP::set_global_p(p, r);

	SpecialPolynomials<MyZP>::init_polynomials(p);

	get_percentile<MyZP, MyZPBits>(rawData, 15.8);
}

