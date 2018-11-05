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

	static unsigned int simd_factor() {
#		ifdef USE_SIMD
		return Number::simd_factor();
#		elif defined DONT_USE_SIMD
		return 1;
#		else
#		error "must define either USE_SIMD or DONT_USE_SIMD"
#		endif
	}

public:
	Distances(const Sites &sites, const Point2D<int> &q, ThreadPool *threads = NULL) : _sites(sites), _query(q) {}
	Distances(const Sites &sites, const std::vector<int> &classes, const Point2D<int> &q, ThreadPool *threads = NULL) : _sites(sites), _classes(classes), _query(q) {}

	unsigned int size() const { return _sites.size(); }
	Number operator[] (unsigned int at) {
		return getDistances(at);
	}

	std::vector<long int> getClass(unsigned int at, int cls) {
		std::vector<long int> ret;
		ret.resize(Number::simd_factor());
		for (unsigned int i = 0; i < Number::simd_factor(); ++i)
			if ((i < simd_factor()) && (_classes.size() > at + i) && (_classes[at + i] == cls))
				ret[i] = 1;
			else
				ret[i] = 0;
		return ret;
	}

	std::vector<long int> getPlaintextClasses(unsigned int at) {
		std::vector<long int> ret;
		ret.resize(Number::simd_factor());
		for (unsigned int i = 0; i < Number::simd_factor(); ++i)
			if ((i < simd_factor()) && (_classes.size() > at + i))
				ret[i] = _classes[at + i];
			else
				ret[i] = 0;
		return ret;
	}

	std::vector<long int> getPlaintextDistances(int at) {
		std::vector<long int> ret;
		ret.resize(Number::simd_factor());
		for (unsigned int i = 0; i < Number::simd_factor(); ++i)
			if ((i < simd_factor()) && (_sites.size() > at + i)) {
				Point2D<int> dist = _query - _sites[at + i];
				ret[i] = abs(dist.x) + abs(dist.y);
			} else
				ret[i] = 0;
		return ret;
	}


	Number getDistances(unsigned int at, ThreadPool *threads = NULL) {
		std::string fname = _cache[at];
		if (fname != std::string("")) {
			std::ifstream f( fname );
			Number ret;
			f >> ret;
			return ret;
		}

		std::vector<long int> xSite(Number::simd_factor());
		std::vector<long int> ySite(Number::simd_factor());

		std::vector<long int> xQuery(Number::simd_factor());
		std::vector<long int> yQuery(Number::simd_factor());

		unsigned int i = 0;
		while ((i < simd_factor()) && (i + at < _sites.size())) {
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



//		NumberBits retBits;
//		convert_to_bits<NumberBits, Number>(retBits, ret, threads);

#		ifdef CACHE_CTEXT
		std::stringstream sfname;
		sfname << cache_prefix << at << ".ctxt";

		_cache[at] = sfname.str();
		ofstream f(sfname.str());
		f << ret;
#		endif

		return ret;
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
			for (unsigned int i = 0; i < simd_factor(); ++i) {
				++_loc;
				if (_loc == _db.size())
					return;
			}
		}
		void operator++(int) { operator++(); }

		void operator+=(int a) { for (int i = 0; i < a; ++i) operator++(); }

		std::vector<long int> getClass(int cls) { return _db.getClass(_loc, cls); }
		Number getDistances() { return _db.getDistances(_loc, _threads); }
		Number operator*() { return getDistances(); }

		std::vector<long int> getPlaintextDistances() { return _db.getPlaintextDistances(_loc); }
		std::vector<long int> getPlaintextClasses() { return _db.getPlaintextClasses(_loc); }

		unsigned int loc() const { return _loc; }
	};


	Iterator begin() { return Iterator(*this, 0); }
	Iterator end() { return Iterator(*this, _sites.size()); }
};

#endif
