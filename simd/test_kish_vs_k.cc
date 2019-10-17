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

#define SIMD_FACTOR 600
typedef ZP MyZP;
typedef UnsignedWord<18, MyZP> MyZPBits;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::square_msd_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::sqrt_msd_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::sqrt_polynomial;

template<class Number>
Polynomial<Number> SpecialPolynomials<Number>::is_positive_polynomial;

//template<class Number>
//std::vector<Polynomial<Number> > SpecialPolynomials<Number>::convert_to_bit;

int main(int argc, char **argv) {
	ZP::set_global_simd_factor(SIMD_FACTOR);
	
	RETRIES = 5;
	MAX_CANDIDATES = -1;

	srand(time(NULL));

	initialize(argc, argv);
	MyZP::set_global_p(p, r);

//	SpecialPolynomials<MyZP>::init_polynomials(p);

//	secure_geo_search<MyZP, MyZPBits>(rawDiscreteData, discreteQuery);

//	secure_knn_classifier<MyZP, MyZPBits>(rawDiscreteData, rawDataClasses, discreteQuery);

	test_kish_classifier(rawDiscreteData, rawDataClasses);
	return 0;
}

int main2(int argc, char **argv) {
//	initialize(argc, argv);

	p = 23;
	r = 1;

	MyZP::set_global_p(p, r);

	SpecialPolynomials<MyZP>::init_polynomials(p);

	for (int i = 0; i < 1; ++i) {
		MyZP n(i);
		MyZP m = SpecialPolynomials<MyZP>::sqrt_polynomial.compute(n);

		std::cerr << "sqrt(" << n.to_int() << ") = " << m.to_int() << std::endl;
	}
	return 0;
}

