#ifndef ___ONE_MEAN_FRAMEWORK___
#define ___ONE_MEAN_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <vector>

#include "point.h"
#include "one_mean.h"

//template<class NUMBER>
//TempPoint sum_of_points(std::vector< Point<Number> > &points) {
//	BinomialSearch< Point<NUMBER> > sum( BinomialSearch< Point<NUMBER> >::add );
//
//	for (int i = 0; i < points.size(); ++i)
//		sum.add_to_tournament( points[i] );
//
//	TempPoint ret = sum.unite_all();
//	return ret;
//}
//
//template<class Number>
//TempPoint average_of_points(std::vector< Point<Number> > &points) {
//	TempPoint ret = sum_of_points(points);
//	ret = ret * (points[0][0].p() / points.size());
//	return ret;
//}

//template<class Number>
//Point<Number> average_of_points_by_probability(Point<Number> *points, int size) {
//	int dim = points[0].dim();
//
//	Point<Number> avg(dim);
//
//	for (int d = 0; d < dim; ++d) {
//		PointCoordinateIterator<Number> begin(points, size, d);
//		PointCoordinateIterator<Number> end(points, size, d, size);
//
//		avg[d] = average(begin, end);
//	}
//
//	return avg;
//}

Point<int> round_point(const Point<double> &p) {
	Point<int> ret;
	ret.dim(p.dim());
	for (int i = 0; i < p.dim(); ++i) {
		ret[i] = (int)(p[i]);
	}
	return ret;
}

template<class NUMBER>
Point<double> to_double(const Point<NUMBER> &p) {
	Point<double> ret;
	ret.dim(p.dim());
	for (int i = 0; i < p.dim(); ++i) {
		ret[i] = p[i].to_int();
	}
	return ret;
}

template<>
Point<double> to_double(const Point<int> &p) {
	Point<double> ret;
	ret.dim(p.dim());
	for (int i = 0; i < p.dim(); ++i) {
		ret[i] = p[i];
	}
	return ret;
}


template<class NUMBER>
Point<NUMBER> encode_point(const Point<double> &p, const Point<double> &min, int resolution) {
	Point<int> roundPoint = round_point((p - min) * resolution);
	Point<NUMBER> ret;
	ret.dim(roundPoint.dim());
	for (int j = 0; j < ret.dim(); ++j) {
//std::cerr << "encoding " << roundPoint[j] << std::endl;
		ret[j] = NUMBER::static_from_int(roundPoint[j]);
//std::cerr << ret[j] << std::endl;
	}
	return ret;
}

template<class NUMBER>
void encode_points(std::vector<Point<NUMBER> > &ret, const std::vector<Point<double> > &p, const Point<double> &min, int resolution) {
std::cout << "there are " << p.size() << " points " << std::endl;
	for (int i = 0; i < p.size(); ++i) {
		Point<int> roundPoint = round_point((p[i] - min) * resolution);
		Point<NUMBER> encPoint;
		encPoint.dim(roundPoint.dim());
		for (int j = 0; j < encPoint.dim(); ++j) {
//std::cerr << "encoding " << roundPoint[j] << std::endl;
			encPoint[j] = NUMBER::static_from_int(roundPoint[j]);
//std::cerr << encPoint[j] << std::endl;
		}
		ret.push_back(encPoint);
	}
}

//template<class NUMBER>
//void encode_points2(const std::vector<Point<double> > &p, const Point<double> &min, int resolution) {
//	std::vector<Point<NUMBER> > ret;
//
//	for (int i = 0; i < p.size(); ++i) {
//		Point<int> roundPoint = round_point((p[i] - min) * resolution);
//		Point<NUMBER> encPoint;
//		encPoint.dim(roundPoint.dim());
//		for (int j = 0; j < encPoint.dim(); ++j)
//			encPoint[j] = NUMBER::static_from_int(roundPoint[j]);
//		ret.push_back(encPoint);
//	}
//	return ret;
//}

template<class NUMBER>
Point<double> decode_point(const Point<NUMBER> &p, const Point<double> &min, int resolution) {
	return (to_double(p) * ((double)1/resolution)) + min;
}

//template<class OutNumber, class InNumber>
//Point<OutNumber> delme_average_by_probability(std::vector<Point<InNumber> > &p) {
//	Point<OutNumber> ret;
//	ret.dim(p[0].dim());
//
//	clock_t start = clock();
//	for (int d = 0; d < p[0].dim(); ++d) {
//		PointVectorCoordinateIterator<InNumber> begin(p, d, 0);
//		PointVectorCoordinateIterator<InNumber> end(p, d, p.size());
//
//		ret[d] = average_by_probability<OutNumber, CompareNative<InNumber>, TruncConversion, PointVectorCoordinateIterator<InNumber> >(begin, end);
//	}
//	clock_t end = clock();
//	return ret;
//}
//
//template<class OutNumber, class InNumber>
//void delme_add_to_average(Point<OutNumber> &avg, const Point<InNumber> &p, int n) {
//	for (int d = 0; d < p.dim(); ++d) {
//		add_to_average<OutNumber, CompareNative<InNumber>, TruncConversion, PointVectorCoordinateIterator<InNumber> >(avg[d], p[d], n);
//	}
//}





///////////////////////////////
// There be dragons above
///////////////////////////////

template<class OutNumber, class TempNumber>
Point<OutNumber> average(const std::vector<Point<double> > &rawPoints, const Point<double> &min, int resolution) {
	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgX(rawPoints.size());
	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgY(rawPoints.size());

	std::cout << std::endl << std::endl;

	for (int i = 0; i < rawPoints.size(); ++i) {
		Point<TempNumber> encPoint = encode_point<TempNumber>(rawPoints[i], min, resolution);
		std::cout << i << ":  " << rawPoints[i] << "  =>  " << convert_to_int_point<TempNumber>(encPoint) << std::endl;
		avgX.add(encPoint[0]);
		avgY.add(encPoint[1]);
	}

	Point<OutNumber> avg;
	avg.dim(2);
	avg[0] = avgX.getAverage();
	avg[1] = avgY.getAverage();

	return avg;
}


template<class OutNumber, class TempNumber>
void old_one_mean(const std::vector<Point<double> > &rawPoints, const Point<double> &minPointManhattan, int resolution) {

	clock_t start = clock();
	Point<OutNumber> encAvg = average<OutNumber, TempNumber>(rawPoints, minPointManhattan, resolution);
	clock_t end = clock();

	std::cout << "computing average took " << (end - start) << " microsecs" << std::endl;

	std::cout << "encAvg = " << encAvg << std::endl;

	Point<double> myAvg = decode_point(encAvg, minPointManhattan, resolution);

	std::cout << "HE average: " << myAvg << std::endl;
}





// helper function to compute 1/n \sum x_i^2
template<class Number>
Number sqr_Xi(const Number &x) {
	return x;
}




template<class OutNumber, class TempNumber, class InNumber, class Compare, class Convert>
void add_point_to_1_mean(
		AverageLiphe<OutNumber, TempNumber, Compare, Convert> &avgX,
		AverageLiphe<OutNumber, TempNumber, Compare, Convert> &avgY,
		AverageLiphe<OutNumber, TempNumber, Compare, Convert> &avgNorm,
		Point<InNumber> p1,
		Point<OutNumber> p2) {

	avgX.add(p1[0]);
	avgY.add(p1[1]);
	avgNorm.add_with_cost(p1[0], p2[0]);
	avgNorm.add_with_cost(p1[1], p2[1]);
}


time_t start = 0;
time_t total = 0;

void start_stat_interval() {
	start = clock();
}

void end_stat_interval() {
	total += clock() - start;
}

void print_stat_prefix(int i) {
	std::cout << "After ";

	if (i == -1)
		std::cout << "all";
	else
		std::cout << i;

	std::cout << " ";
}

int get_mem() {
	struct rusage rus;
	int res;

	res = getrusage(RUSAGE_SELF, &rus);
	if (res == -1) {
		perror("rusage");
		return 0;
	}

	return rus.ru_maxrss / 1024;
}

std::vector<Point<double> > roundedPoints;
Point<int> roundedAvgPoint;
int roundedAvgNorm;
int n;

void print_stat(int i) {
	print_stat_prefix(i);
	std::cout << "Took " << total << " micro" << std::endl;

	print_stat_prefix(i);
	std::cout << "Used " << get_mem() << " MegaBytes" << std::endl;
}

template<class NUMBER>
void print_detailed_stats(int step, Point<NUMBER> &avgTempPoint, int avgTempNorm) {
	int x = avgTempPoint[0].to_int();
	int y = avgTempPoint[1].to_int();

	print_stat_prefix(step);
	std::cout << "AvgX " << x << " / " << (roundedAvgPoint[0] / n) << " = " << ((double)x * n / roundedAvgPoint[0]) << std::endl;

	print_stat_prefix(step);
	std::cout << "AvgY " << y << " / " << (roundedAvgPoint[1] / n) << " = " << ((double)y * n / roundedAvgPoint[1]) << std::endl;

	print_stat_prefix(step);
	std::cout << "AvgNorm " << avgTempNorm << " / " << (roundedAvgNorm / n) << " = " << ((double)avgTempNorm * n / roundedAvgNorm) << std::endl;

}

inline int sqr(int x) { return x*x; }

int statistic_rate = 10;

template<class OutNumber, class TempNumber>
void one_mean(Point<double> &avgPoint, double &avgNorm, const std::vector<Point<double> > &rawPoints, const Point<double> &minPointManhattan, int resolution) {
	avgPoint.dim(minPointManhattan.dim());
	n = rawPoints.size();

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgX(rawPoints.size());
	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgY(rawPoints.size());

//	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgN(rawPoints.size(), sqr_Xi<OutNumber>);
	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgN(rawPoints.size());

	roundedAvgPoint.dim(minPointManhattan.dim());
	roundedAvgPoint[0] = roundedAvgPoint[1] = 0;
	roundedAvgNorm = 0;

	for (int i = 0; i < rawPoints.size(); ++i) {
		Point<OutNumber> roundedPoint = encode_point<OutNumber>(rawPoints[i], minPointManhattan, resolution);
		int x = roundedPoint[0].to_int();
		int y = roundedPoint[1].to_int();

		roundedPoints.push_back(decode_point(roundedPoint, minPointManhattan, resolution));

		roundedAvgPoint[0] += x;
		roundedAvgPoint[1] += y;
		roundedAvgNorm += x*x + y*y;

		Point<TempNumber> encBin = encode_point<TempNumber>(rawPoints[i], minPointManhattan, resolution);
		Point<OutNumber> encNative = encode_point<OutNumber>(rawPoints[i], minPointManhattan, resolution);

		std::cout << i << ":  " << rawPoints[i] << "  =>  " << convert_to_int_point<TempNumber>(encBin) << "  => " << convert_to_int_point<OutNumber>(encNative) << std::endl;

		start_stat_interval();
		add_point_to_1_mean(avgX, avgY, avgN, encBin, encNative);
		end_stat_interval();

		if ((i % statistic_rate) == statistic_rate - 1) {
			print_stat(i);

			if (measureAccuracy) {
				Point<OutNumber> avgTempPointEnc;
				avgTempPointEnc.dim(2);
				avgTempPointEnc[0] = avgX.getAverage();
				avgTempPointEnc[1] = avgY.getAverage();

//				Point<double> avgTempPoint = decode_point(avgTempPointEnc, minPointManhattan, resolution);

				OutNumber avgTempNormEnc;
				avgTempNormEnc = avgN.getAverage();

				print_detailed_stats(i, avgTempPointEnc, avgTempNormEnc.to_int());
			}
		}
	}


	Point<OutNumber> avgPointEnc;
	avgPointEnc.dim(2);
	avgPointEnc[0] = avgX.getAverage();
	avgPointEnc[1] = avgY.getAverage();

	OutNumber avgNormEnc;
	avgNormEnc = avgN.getAverage();

	print_stat(-1);
	if (measureAccuracy) {
		print_detailed_stats(-1, avgPointEnc, avgNormEnc.to_int());
	}


	avgPoint[0] = avgPointEnc[0].to_int();
	avgPoint[1] = avgPointEnc[1].to_int();

	avgNorm = avgNormEnc.to_int() / sqr(resolution)   -   sqr(minPointManhattan[0]) - sqr(minPointManhattan[1])  + 2*minPointManhattan[0]*avgPoint[0] + 2*minPointManhattan[1]*avgPoint[1];





	// the distance to (100,100)
	int testX = 100;
	int testY = 100;

	int realDist = 0;
	for (int i = 0; i < roundedPoints.size(); ++i) {
		Point<OutNumber> tempPoint;
		tempPoint.dim(2);
		tempPoint[0].from_int(roundedPoints[i][0]);
		tempPoint[1].from_int(roundedPoints[i][1]);
		Point<double> p = decode_point(tempPoint, minPointManhattan, resolution);
		realDist += sqr(p[0] - testX) + sqr(p[1] - testY);
	}
	realDist /= roundedPoints.size();

	int coreSetDistance =
		sqr(testX) + sqr(testY)
		+ avgNorm
		- 2*(testX*avgPoint[0] + testY*avgPoint[1]);

	std::cout << "Coreset gives: " << coreSetDistance << " / " << realDist << " = " << ((double)coreSetDistance / realDist) << std::endl;
}


// what we wan to test:
// 1. time
// 2. accuracy, once every 100 points check how accurate the average is (only in ZP)
// 3. accuracy, once evey 100 poinst, check how accurate the coreset is (only in ZP)


#endif
