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
	const std::vector<int> _classes;
	Point2D<int> _query;
	std::map<int, std::string> _cache;
public:
	Distances(const Sites &sites, const Point2D<int> &q, ThreadPool *threads = NULL) : _sites(sites), _query(q) {}
	Distances(const Sites &sites, const std::vector<int> &classes, const Point2D<int> &q, ThreadPool *threads = NULL) : _sites(sites), _classes(classes), _query(q) {}

	unsigned int size() const { return _sites.size(); }
	NumberBits operator[] (unsigned int at) {
		return getDistances(at);
	}

	std::vector<long int> getClass(unsigned int at, int cls) {
		std::vector<long int> ret;
		ret.resize(Number::simd_factor());
		for (unsigned int i = 0; i < Number::simd_factor(); ++i)
			ret[i] = (_classes[at + i] == cls) ? 1 : 0;
		return ret;
	}

	int getPlaintextDistance(int i) {
		Point2D<int> dist = _query - _sites[i];
		return abs(dist.x) + abs(dist.y);
	}

	NumberBits getDistances(unsigned int at, ThreadPool *threads = NULL) {
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


		Number ret;

		{
			Number xSiteEnc(xSite);
			NumberBits xSiteBitsEnc(xSite);

			Number xQueryEnc(xQuery);
			NumberBits xQueryBitsEnc(xQuery);

			ret = (xSiteEnc - xQueryEnc) * ((xQueryBitsEnc < xSiteBitsEnc)*2 - 1);
		}


		{
			Number ySiteEnc(ySite);
			NumberBits ySiteBitsEnc(ySite);

			Number yQueryEnc(yQuery);
			NumberBits yQueryBitsEnc(yQuery);

			ret += (ySiteEnc - yQueryEnc) * ((yQueryBitsEnc < ySiteBitsEnc)*2 - 1);
		}



		NumberBits retBits;
		convert_to_bits<NumberBits, Number>(retBits, ret, threads);

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
		ThreadPool *_threads;
	public:
		Iterator(Distances &db, ThreadPool *thr = NULL) : _db(db), _loc(0), _threads(thr) {}
		Iterator(Distances &db, int i, ThreadPool *thr = NULL) : _db(db), _loc(i), _threads(thr) {}
		void set_thread_pool(ThreadPool *t) { _threads = t; }

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

		std::vector<long int> getClass(int cls) { return _db.getClass(_loc, cls); }
		NumberBits getDistances() { return _db.getDistances(_loc, _threads); }
		NumberBits operator*() { return getDistances(); }
	};


	Iterator begin() { return Iterator(*this, 0); }
	Iterator end() { return Iterator(*this, _sites.size()); }
};

#endif
