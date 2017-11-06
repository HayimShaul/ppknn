#ifndef ___INPUT_ITERATOR___
#define ___INPUT_ITERATOR___

#include <map>
#include <fstream>

#include "special_polynomials.h"

#define cache_prefix "ctxt_cache-"

template<class Number, class NumberBits>
class Distances {
public:
	typedef std::vector<Point2D<int> >  Sites;

private:
	const Sites &_sites;
	Point2D<int> _query;
	std::map<int, std::string> _cache;
public:
	Distances(const Sites &sites, const Point2D<int> &q) : _sites(sites), _query(q) {}

	unsigned int size() const { return _sites.size(); }
	NumberBits operator[] (unsigned int at) {

		std::string fname = _cache[at];
		if (fname != std::string("")) {
			std::ifstream f( fname );
			NumberBits ret;
			f >> ret;
			return ret;
		}

		std::vector<long int> xSite(Number::simd_factor());
		std::vector<long int> ySite(Number::simd_factor());

		std::vector<long int> xQuery(Number::simd_factor());
		std::vector<long int> yQuery(Number::simd_factor());

		unsigned int i = 0;
		while ((i < Number::simd_factor()) && (i + at < _sites.size())) {
			xSite[i] = _sites[i + at].x;
			ySite[i] = _sites[i + at].y;
			xQuery[i] = _query.x;
			yQuery[i] = _query.y;
			++i;
		}
		while (i < Number::simd_factor()) {
			xSite[i] = 0;
			ySite[i] = 0;
			xQuery[i] = 0;
			yQuery[i] = 0;
			++i;
		}

		Number xSiteEnc(xSite);
		Number ySiteEnc(ySite);

		Number xQueryEnc(xQuery);
		Number yQueryEnc(yQuery);

		Number ret = SpecialPolynomials<Number>::abs_polynomial.compute(xSiteEnc - xQueryEnc) + SpecialPolynomials<Number>::abs_polynomial.compute(ySiteEnc - yQueryEnc);

		NumberBits retBits;
		convert_to_bits<NumberBits, Number>(retBits, ret, NULL);

		std::stringstream sfname;
		sfname << cache_prefix << at << ".ctxt";

		_cache[at] = sfname.str();
		ofstream f(sfname.str());
		f << retBits;

		return retBits;
	}

	class Iterator {
	private:
		Distances &_db;
		unsigned int _loc;
	public:
		Iterator(Distances &db) : _db(db), _loc(0) {}
		Iterator(Distances &db, int i) : _db(db), _loc(i) {}

		bool operator==(const Iterator &i) const { return _loc == i._loc; }
		bool operator!=(const Iterator &i) const { return !operator==(i); }


		void operator++() {
			if (_loc == _db.size())
				return;
			for (unsigned int i = 0; i < Number::simd_factor(); ++i) {
				++_loc;
				if (_loc == _db.size())
					return;
			}
		}
		void operator++(int) { operator++(); }

		void operator+=(int a) { for (int i = 0; i < a; ++i) operator++(); }

		NumberBits operator*() { return _db[_loc]; }
	};


	Iterator begin() { return Iterator(*this, 0); }
	Iterator end() { return Iterator(*this, _sites.size()); }
};

#endif
