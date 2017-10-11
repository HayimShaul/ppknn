#include <iostream>
#include <string.h>

#include <carray_iterator.h>
#include <helib_number.h>
#include <eq.h>
#include "polynomial.h"
#include <cmp.h>
#include <binomial_tournament.h>
#include <unsigned_word.h>
#include <average.h>

#include "get_percentile_framework.h"
#include "get_percentile.h"
#include "special_polynomials.h"

typedef UnsignedWord<18, HelibNumber> HelibBits;


template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::square_msd_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::sqrt_polynomial;

template<class Number>
std::vector<Polynomial<Number> > SpecialPolynomials<Number>::convert_to_bit;



int main(int argc, char**argv) {
	measureAccuracy = false;

	initialize(argc, argv);

	HelibKeys keys;

	long R = 1;
	long d = 1;
	long c = 2;
	long k = 80;
	long s = 0;
	long chosen_m = 0;
	Vec<long> gens;
	Vec<long> ords;

	keys.initKeys(s, R, p, r, d, c, k, 64, L, chosen_m, gens, ords);
	HelibNumber::set_global_keys(&keys);

	SpecialPolynomials<HelibNumber>::init_polynomials(p);

	get_percentile<HelibNumber, HelibBits>(rawData, 15.8);
	return 0;
}

int main2(int argc, char**argv) {
	measureAccuracy = false;


	HelibKeys keys;

	p = 2;
	r = 4;

	long R = 1;
	long d = 1;
	long c = 2;
	long k = 80;
	long s = 0;
	long chosen_m = 0;
	Vec<long> gens;
	Vec<long> ords;

	keys.initKeys(s, R, p, r, d, c, k, 64, L, chosen_m, gens, ords);
	HelibNumber::set_global_keys(&keys);

	HelibNumber n(7);
	HelibBits nBits;

	nBits  = n.template to_digits<HelibBits>();

	std::cout << "bit[0] = " << nBits[0].to_int() << std::endl;
	std::cout << "bit[1] = " << nBits[1].to_int() << std::endl;
	std::cout << "bit[2] = " << nBits[2].to_int() << std::endl;
	std::cout << "bit[3] = " << nBits[3].to_int() << std::endl;

	HelibNumber bit = nBits[2];
	std::cout << "bit = " << bit.to_int() << std::endl;

	bit *= 8;
	std::cout << "2bit = " << bit.to_int() << std::endl;

	return 0;
}
