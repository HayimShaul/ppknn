#ifndef ___POINT___
#define ___POINT___

#include <assert.h>
#include <ostream>
#include <vector>

template<class Number>
class Point {
public:
	std::vector<Number> _coordinate;

	Point() {}
	Point(const std::vector<Number> &a) : _coordinate(a) {}

	unsigned int dim() const { return _coordinate.size(); }
	void dim(unsigned int d) { _coordinate.resize(d); }
	Number &operator[](unsigned int i) { assert(i < dim()); return _coordinate[i]; }
	const Number &operator[](unsigned int i) const { assert(i < dim()); return _coordinate[i]; }

	bool operator==(const Point &a) const {
		assert(_coordinate.size() == a._coordinate.size());
		for (unsigned int i = 0; i < _coordinate.size(); ++i)
			if (_coordinate[i] != a._coordinate[i])
				return false;
		return true;
	}
	bool operator!=(const Point &a) const { return !operator==(a); }

	bool operator<(const Point &a) const {
		assert(_coordinate.size() == a._coordinate.size());
		for (unsigned int i = 0; i < _coordinate.size(); ++i) {
			if (_coordinate[i] < a._coordinate[i])
				return true;
			if (_coordinate[i] > a._coordinate[i])
				return false;
		}
		return false;
	}
	bool operator>(const Point &a) const { return a < (*this); }
	bool operator<=(const Point &a) const { return !operator>(a); }
	bool operator>=(const Point &a) const { return !operator<(a); }

	Point operator-(const Point &a) const {
		assert(_coordinate.size() == a._coordinate.size());
		Point ret;
		ret._coordinate.resize(_coordinate.size());
		for (unsigned int i = 0; i < _coordinate.size(); ++i)
			ret._coordinate[i] = _coordinate[i] - a._coordinate[i];
		return ret;
	}
	Point operator+(const Point &a) const {
		assert(_coordinate.size() == a._coordinate.size());
		Point ret;
		ret._coordinate.resize(_coordinate.size());
		for (unsigned int i = 0; i < _coordinate.size(); ++i)
			ret._coordinate[i] = _coordinate[i] + a._coordinate[i];
		return ret;
	}
	Point operator*(const Point &a) const {
		assert(_coordinate.size() == a._coordinate.size());
		Point ret;
		ret._coordinate.resize(_coordinate.size());
		for (unsigned int i = 0; i < _coordinate.size(); ++i)
			ret._coordinate[i] = _coordinate[i] * a._coordinate[i];
		return ret;
	}
	Point operator/(const Point &a) const {
		assert(_coordinate.size() == a._coordinate.size());
		Point ret;
		ret._coordinate.resize(_coordinate.size());
		for (unsigned int i = 0; i < _coordinate.size(); ++i)
			ret._coordinate[i] = _coordinate[i] / a._coordinate[i];
		return ret;
	}

	Number normL2sqr() const {
		Number ret = 0;
		for (unsigned int i = 0; i < dim(); ++i) {
			ret += _coordinate[i] * _coordinate[i];
		}
		return ret;
	}

	Number normL1() const {
		Number ret = 0;
		for (unsigned int i = 0; i < dim(); ++i) {
			ret += abs(_coordinate[i]);
		}
		return ret;
	}
};

template<class Number>
inline Point<Number> operator/(float a, const Point<Number> &b) {
	std::vector<Number> v(b.dim());
	for (unsigned int i = 0; i < b.di(); ++i)
		v[i] = a/b[i];
	return Point<Number>(v);
}


template<class Number>
inline Point<Number> min(const Point<Number> &a, const Point<Number> &b) {
	std::vector<Number> v(b.dim());
	for (unsigned int i = 0; i < b.dim(); ++i)
		v[i] = min(a[i], b[i]);
	return Point<Number>(v);
}

template<class Number>
inline Point<Number> max(const Point<Number> &a, const Point<Number> &b) {
	std::vector<Number> v(b.dim());
	for (unsigned int i = 0; i < b.dim(); ++i)
		v[i] = max(a[i], b[i]);
	return Point<Number>(v);
}

template<class Number>
inline std::ostream &operator<<(std::ostream &out, const Point<Number> &p) {
	out << "(" << p[0];
	for (unsigned int i = 1; i < p.dim(); ++i) 
		out << ", " << p[i];
	out << ")";
	return out;
}

extern Point<float> discreteBase;
extern Point<float> discreteResolution;

inline int myround(float f) { return (int)(f+0.5); }

inline Point<int> discretify(const Point<float> &a) {
	Point<float> _a = (a - discreteBase) * discreteResolution;
	std::vector<int> v(a.dim());
	for (unsigned int i = 0; i < a.dim(); ++i)
		v[i] = myround(_a[i]);
	return Point<int>(v);
}

inline Point<float> undiscretify(const Point<int> &a) {
	std::vector<float> v(a.dim());
	for (unsigned int i = 0; i < a.dim(); ++i)
		v[i] = a[i];
	Point<float> _a(v);
	return _a * discreteResolution + discreteBase;
}

#endif
