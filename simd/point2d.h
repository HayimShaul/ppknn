#ifndef ___POINT_2D__
#define ___POINT_2D__

#include <ostream>

template<class Number>
class Point2D {
public:
	Number x, y;

	Point2D() : x(0), y(0) {}
	Point2D(Number a, Number b) : x(a), y(b) {}

	bool operator==(const Point2D &a) const { return (x == a.x) && (y == a.y); }
	bool operator!=(const Point2D &a) const { return !operator==(a); }

	bool operator<(const Point2D &a) const { return (x < a.x) || ((x == a.x) && (y < a.y)); }
	bool operator>(const Point2D &a) const { return a < (*this); }
	bool operator<=(const Point2D &a) const { return !operator>(a); }
	bool operator>=(const Point2D &a) const { return !operator<(a); }

	Point2D operator-(const Point2D &a) const { return Point2D(x-a.x, y-a.y); }
	Point2D operator+(const Point2D &a) const { return Point2D(x+a.x, y+a.y); }
	Point2D operator*(const Point2D &a) const { return Point2D(x*a.x, y*a.y); }
	Point2D operator/(const Point2D &a) const { return Point2D(x/a.x, y/a.y); }
};

template<class Number>
inline Point2D<Number> operator/(float a, const Point2D<Number> &b) { return Point2D<Number>(a/b.x, a/b.y); }

template<class Number>
inline Point2D<Number> min(const Point2D<Number> &a, const Point2D<Number> &b) { return Point2D<Number>(min(a.x, b.x), min(a.y, b.y)); }

template<class Number>
inline Point2D<Number> max(const Point2D<Number> &a, const Point2D<Number> &b) { return Point2D<Number>(max(a.x, b.x), max(a.y, b.y)); }

template<class Number>
inline std::ostream &operator<<(std::ostream &out, const Point2D<Number> &p) {
	out << "(" << p.x << ", " << p.y << ")";
	return out;
}

extern Point2D<float> discreteBase;
extern Point2D<float> discreteResolution;

inline int myround(float f) { return (int)(f+0.5); }

inline Point2D<int> discretify(const Point2D<float> &a) {
	Point2D<float> _a = (a - discreteBase) * discreteResolution;
//	return Point2D<int>((int)_a.x, (int)_a.y);
	return Point2D<int>(myround(_a.x), myround(_a.y));
}

inline Point2D<float> undiscretify(const Point2D<int> &a) {
	Point2D<float> _a = Point2D<float>(a.x, a.y);
	return _a * discreteResolution + discreteBase;
}

#endif
