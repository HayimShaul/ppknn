#ifndef __POLYNOMIAL__
#define __POLYNOMIAL__

#include <string.h>
#include <algorithm>
#include "assert.h"


template<class Number>
class Polynomial {
private:
	std::vector<int> _coef;
	int _mod;

	void set(int n, const char *x = NULL) {
		_coef.resize(2);
		_coef[0] = n;

		if (x == NULL) {
			_coef.resize(1);
		} else if (strcmp(x, "-x") == 0) {
			_coef[1] = -1;
		} else if ((strcmp(x, "+x") == 0) || (strcmp(x, "x") == 0)) {
			_coef[1] = 1;
		} else {
			fprintf(stderr, "unknown term %s\n", x);
			exit(1);
		}
	}
public:
	Polynomial(const Polynomial<Number> &p) : _coef(p._coef), _mod(p._mod) {}
	Polynomial(const std::vector<Number> &coef) : _coef(coef), _mod(0) {}

	Polynomial(const char *x, int n) : _mod(0) { set(n, x); }
	Polynomial(int n, const char *x = NULL) : _mod(0) { set(n, x); }

	Polynomial &set_mod(int n) { _mod = n; return *this; }

	Polynomial operator*(const Polynomial<Number> &p) const { Polynomial<Number> q(*this); q *= p; return q; }
	void operator*=(const Polynomial<Number> &p) {
		if (_mod == 0)
			_mod = p._mod;
		assert((_mod == p._mod) || (p._mod == 0));

		std::vector<int> c;
		c.resize(_coef.size() + p._coef.size());
		for (int i = 0; i < c.size(); ++i)
			c[i] = 0;

		for (int this_i = 0; this_i < _coef.size(); ++this_i) {
			for (int p_i = 0; p_i < p._coef.size(); ++p_i) {
				c[this_i + p_i] += _coef[this_i] * p._coef[p_i];
			}
		}

		_coef = c;

		apply_mod();
	}


	Polynomial<Number> operator+(const Polynomial<Number> &p) const { Polynomial<Number> q(*this); q += p; return q; }
	void operator+=(const Polynomial<Number> &p) {
		if (_mod == 0)
			_mod = p._mod;
		assert((_mod == p._mod) || (p._mod == 0));

		int prev_size = _coef.size();
		if (p._coef.size() > _coef.size())
			_coef.resize(p._coef.size());
		for (int i = prev_size - 1; i < _coef.size(); ++i)
			_coef[i] = 0;

		for (int i = 0; i < p._coef.size(); ++i) {
			_coef[i] += p._coef[i];
		}

		apply_mod();
	}

	Polynomial<Number> operator-(const Polynomial<Number> &p) const { Polynomial<Number> q(*this); q -= p; return q; }
	void operator-=(const Polynomial<Number> &p) { operator+=(-p); }

	Polynomial<Number> operator+(int p) const { return *this + Polynomial<Number>(p); }
	void operator+=(int p) { *this += Polynomial<Number>(p); }

	Polynomial<Number> operator-(int p) const { return *this - Polynomial<Number>(p); }
	void operator-=(int p) { *this -= Polynomial<Number>(p); }

	Polynomial<Number> operator-() const {
		Polynomial<Number> ret(*this);
		ret.set_mod(_mod);
		for (int i = 0; i < ret._coef.size(); ++i)
			ret._coef[i] = -ret._coef[i];
		ret.apply_mod();
		return ret;
	}

	// exponent of polynomials
	// no need to optimize this because this is happening in plain text anyway
	void operator^=(int p) { *this = (*this) ^ p; }
	Polynomial<Number> operator^(int p) const {

		if (p == 0)
			return Polynomial<Number>(1).set_mod(_mod);

		if (p == 1)
			return *this;

		Polynomial<Number> ret(*this);

		if ((p % 2) == 0) {
			ret *= ret;
			ret ^= p/2; 
			return ret;
//			return (ret*ret)^(p/2);
		}

		ret = ret * ((ret*ret)^(p/2));
		ret.apply_mod();
		return ret;
	}


	// composition
	Polynomial<Number> operator()(const Polynomial<Number> &p) {
		int mod = _mod;
		if (mod == 0)
			mod = p._mod;
		assert((mod == p._mod) || (p._mod == 0));

		Polynomial<Number> ret(0);
		ret.set_mod(mod);
		for (int i = 0; i < _coef.size(); ++i)
			ret += p^i * _coef[i];
		ret.apply_mod();
		return ret;
	}

	// Modulo, i.e. make it a polynomial in Z_p(x)
	// all coefficients are modulo p
	// and x^{p-1} = 1
	Polynomial<Number> apply_mod() { *this %= _mod; }
	Polynomial<Number> operator%(int p) const { Polynomial<Number> poly(*this);  poly %= p; }
	void operator%=(int p) {
		if (p == 0)
			return;

		std::vector<int> coef;

		int phi = ::phi(p);
		if (phi + 1 < _coef.size()) {
			coef.resize(phi + 1);
			for (int i = 0; i < coef.size(); ++i)
				coef[i] = _coef[i];
		} else {
			coef = _coef;
		}

		for (int i = coef.size(); i < _coef.size(); ++i) {
			// x^{phi+1} = x^1
			if ((i % phi) == 0)
				coef[phi] = _coef[i];
			else
				coef[i % phi] += _coef[i];
		}

		for (int i = 0; i < coef.size(); ++i)
			coef[i] %= p;

		_coef = coef;
	}

	int deg() const {
		int d = _coef.size() - 1;
		while ((d > 0) && (_coef[d] == 0))
			--d;
		return d;
	}


	void set_coef(int n, int c) {
		if (n >= _coef.size())
			_coef.resize(n + 1);
		_coef[n] = c;
	}

	Number batch_coefficient(const std::vector<int> &coef, int start, int end, const Number *powers) const {
		Number ret = coef[start];
		for (int i = start+1; i < end; ++i)
			ret += powers[i-start] * coef[i];
		return ret;
	}

	Number compute(const Number &x) const {
		Polynomial<Number> p = (*this) % x.p();

//		for (int i = 0; i < coef.size(); ++i)
//			std::cout << "    coef[" << i << "] = " << coef[i];
//		std::cout << std::endl;

		int batch_size = sqrt(p.deg());
		if (batch_size * batch_size < p.deg())
			++batch_size;

		int batch_number = p.deg() / batch_size;

		Number *powers = new Number[batch_size];
		// powers holds the values of x^1, x^2, ..., x^{sqrt_d - 1}
		powers[1] = x;
		for (int i = 2; i < batch_size; ++i) {
			int a = i/2;
			int b = (i+1)/2;
			powers[i] = powers[a] * powers[b];
		}

		Number *batch_multiplier = new Number[batch_number];
		// batch_multiplier holds the values of x^{sqrt_d}, x^{2sqrt_d}, ..., x^{d-1}
		batch_multiplier[1] = powers[batch_size/2] * powers[(batch_size+1)/2];

		Number ret = batch_coefficient(p._coef, /*start=*/ 0, /*end=*/batch_size, powers);
		ret += batch_coefficient(p._coef, /*start=*/ batch_size, /*end=*/ 2*batch_size, powers) * batch_multiplier[1];

		for (int i = 2; i < batch_number; ++i) {
			int a = i/2;
			int b = (i+1)/2;
			batch_multiplier[i] = batch_multiplier[a] * batch_multiplier[b];
			ret += batch_coefficient(p._coef, /*start=*/ i*batch_size, (i+1)*batch_size, powers) * batch_multiplier[i];
		}

		if (batch_number * batch_size < p.deg()) {
			int a = batch_number/2;
			int b = (batch_number+1)/2;
			Number batch_mult = batch_multiplier[a] * batch_multiplier[b];
			ret += batch_coefficient(p._coef, /*start=*/ batch_number*batch_size, p.deg(), powers) * batch_mult;
		}

		delete[] powers;
		delete[] batch_multiplier;

		return ret;
	}

};

#endif
