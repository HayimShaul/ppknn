#ifndef ___SPECIAL_POLYNOMIALS___
#define ___SPECIAL_POLYNOMIALS___

#include <vector>
#include "liphe/polynomial.h"

template<class Number>
class SpecialPolynomials {
public:
	static Polynomial<Number> square_msd_polynomial;

	static Polynomial<Number> sqrt_polynomial;

	static Polynomial<Number> sqrt_msd_polynomial;

	static Polynomial<Number> is_positive_polynomial;

//	static std::vector<Polynomial<Number> > convert_to_bit;

	static void init_polynomials(int p) {
		ThreadPool pool;

		pool.submit_job(std::function<void(void)>([p](){
			SpecialPolynomials::square_msd_polynomial = Polynomial<Number>::build_polynomial(p, p, [p](int x)->int{ return x*x/p; } );
		}));

		pool.submit_job(std::function<void(void)>([p](){
			SpecialPolynomials::sqrt_polynomial = Polynomial<Number>::build_polynomial(p, p, [p](int x)->int{ return sqrt(x); } );
		}));

		pool.submit_job(std::function<void(void)>([p](){
			SpecialPolynomials::sqrt_msd_polynomial = Polynomial<Number>::build_polynomial(p, p, [p](int x)->int{ return sqrt(x*p); } );
		}));

		pool.submit_job(std::function<void(void)>([p](){
			SpecialPolynomials::is_positive_polynomial = Polynomial<Number>::build_polynomial(p, p, [p](int x)->int{ return (x < p/2) ? 1 : 0; } );
		}));

//		int log_p = 1;
//		while ((1 << log_p) < p)
//			++log_p;
//		convert_to_bit.resize(log_p);
//		for (int bit = 0; bit < log_p; ++bit) {
//			pool.submit_job(std::function<void(void)>([bit, p](){
//				SpecialPolynomials::convert_to_bit[bit] = Polynomial<Number>::build_polynomial(p, p, [bit](int x)->int{ return (x >> bit) & 1; });
//			}));
//		}

		pool.process_jobs();
	}
};

//template<class NumberBits, class Number>
//void convert_to_bits(NumberBits &out, const Number &x, ThreadPool *threads) {
//	int bits = 0;
//	while ((1 << bits) < p)
//		++bits;
//
//	out.set_bit_length(bits);
//
//	int batch_size = 0;
//	Number *powers = SpecialPolynomials<Number>::convert_to_bit[0].compute_powers(x, batch_size);
//	for (int i = 0; i < bits; ++i) {
//		Number b = SpecialPolynomials<Number>::convert_to_bit[i].compute(x, powers, batch_size, threads);
//		out.set_bit(i, b);
//	}
//
//	delete[] powers;
//}




#endif
