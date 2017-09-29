#ifndef ___GET_PERCENTILE_FRAMEWORK___
#define ___GET_PERCENTILE_FRAMEWORK___

#include <sys/resource.h>
#include <sys/types.h>
#include <vector>

#include "polynomial.h"
#include "get_percentile.h"







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

int accuracyEnhancementBits = 3;
int realAvg;
int realAvgSqr;
int n;

void print_stat(int i) {
	print_stat_prefix(i);
	std::cout << "Took " << total << " micro" << std::endl;

	print_stat_prefix(i);
	std::cout << "Used " << get_mem() << " MegaBytes" << std::endl;
}

template<class NUMBER>
void print_detailed_stats(int step, NUMBER &avgTemp, NUMBER &avgSqrMsdTemp, NUMBER &avgSqrLsdTemp) {
	int avg = avgTemp.to_int();
	int avgSqrLsd = avgSqrLsdTemp.to_int();   // this is x mod p
	int avgSqrMsd = avgSqrMsdTemp.to_int();   // this is x/2^a

	int approx = avgSqrMsd << accuracyEnhancementBits;
	int approxMod = approx % p;

	int avgSqr = approx + avgSqrLsd - approxMod;

	print_stat_prefix(step);
	std::cout << std::endl;
	std::cout << "Avg " << avg << " / " << (realAvg / n) << " = " << ((double)avg * n / realAvg) << std::endl;
	std::cout << "AvgSqrLsd = " << avgSqrLsd << std::endl;
	std::cout << "AvgSqrMsd = " << avgSqrMsd << std::endl;
	std::cout << "AvgSqr " << avgSqr << " / " << (realAvgSqr / n) << " = " << ((double)avgSqr * n / realAvgSqr) << std::endl;

}

inline int sqr(int x) { return x*x; }

int statistic_rate = 10;







template<class Number, class Bits>
class PercentileIterator {
private:
	const std::vector<int> *_raw_data;
	Bits _threshold;
	int _location;

	Number _last_value;
	int _last_value_location;

	void copy(const PercentileIterator &a) { _raw_data = a._raw_data; _threshold = a._threshold; _location = a._location; _last_value_location = a._last_value_location; _last_value = a._last_value; }
public:
	PercentileIterator(const std::vector<int> *arr, const Bits &thr) : _raw_data(arr), _threshold(thr), _location(0), _last_value(-1) {}
	PercentileIterator(const PercentileIterator &a) { copy(a); }

	PercentileIterator &operator=(const PercentileIterator &a) { copy(a); return *this; }

	void begin() { _location = 0; }
	void end() { _location = _raw_data->size(); }

	bool is_end() { return _location == _raw_data->size(); }
	bool is_begin() { return _location == 0; }

	void operator++() { ++_location; }

	bool operator==(const PercentileIterator &a) const {
		return _location == a._location;
	}
	bool operator!=(const PercentileIterator &a) const { return !operator==(a); }

	const Number &operator*() {
		if (_location != _last_value_location) {
			Bits xi = Bits::static_from_int((*_raw_data)[_location]);
			_last_value = Number(1) - (xi < _threshold);

			if (((_last_value.to_int() == 1) && ((xi.to_int() < _threshold.to_int()))) ||
				((_last_value.to_int() == 0) && (!(xi.to_int() < _threshold.to_int())))) {
					std::cout << "Error comparing " << xi.to_int() << " and " << _threshold.to_int() << std::endl;
					_last_value = Number(1) - (xi < _threshold);
			}
			_last_value_location = _location;
		}
		return _last_value;
	}
};


template<class Number>
Polynomial<Number> build_sqrt_polynomial(int p) {
	ZZ_p::init(ZZ(p));
	Polynomial<Number> poly(0);

//	for (int i = 0; i < p; ++i) {
//		std::cerr << "Building sqrt polynomial " << i << " / " << p << "\r";
//		std::cerr << std::endl;
//		Polynomial<Number> i_x = Polynomial<Number>(i, "-x").set_mod(p);
//		i_x ^= phi(p);
//		
////		poly += (Polynomial<Number>(1) - i_x) * sqrt(i);
//////		poly += (Polynomial<Number>(1) - (Polynomial<Number>(i, "-x").set_mod(p))^phi(p)) * sqrt(i);
//	}


	int x = 1;
	int x2 = 1;
	while (x2 < p) {
		Polynomial<Number> val(0);

		while ((0.5+x)*(0.5+x) > x2) {
			val *= Polynomial<Number>(x2, "-x");
			++x2;
		}
		val ^= phi(p);
		val = (-val + 1) * x;

		++x;
	}

	return poly;
}


template<class Number>
Number get_bit(const Number &x, int bit) {
	Polynomial<Number> p(0);

	for (int i = 0; i < x.p(); ++i) {
		if ((i & (1 << bit)) == 1) {
			Polynomial<Number> i_x = Polynomial<Number>(i, "-x");
			i_x ^= phi(x.p());
			p += Polynomial<Number>(1) - i_x;

//			p += Polynomial<Number>(1) - Polynomial<Number>(i, "-x")^phi(x.p());
		}
	}

	return p.compute(x);
}


template<class Number, class BitNumber>
BitNumber get_bits(const Number &x) {
	BitNumber ret;

	int bit_size = 1;
	while ((1 << bit_size) < x.p())
		++bit_size;

	ret.set_bit_length(bit_size);

	for (int i = 0; i < bit_size; ++i) {
		ret.set_bit(i, get_bit(x, i));
	}

	return ret;
}


template<class Number, class NumberBits>
Number fhe_sqrt_approximation(const Number &_x) {
	assert(_x.p() == 2);

	NumberBits x = _x.template to_digits<NumberBits>();

std::cout << "converting to bits " << _x.to_int() << " turned into " << x.to_int() << std::endl;

std::cout << "x.bitlength = " << x.bitLength() << std::endl;

	NumberBits ret;
	ret.set_bit_length((x.bitLength() + 1) / 2);

	for (int i = 0; i < x.bitLength() / 2; ++i) {
		ret[i] = x[2*i] + x[2*i+1] - ( x[2*i] * x[2*i+1] );
	}

	if ((x.bitLength() % 2) == 1) {
		ret[x.bitLength() / 2] = x[x.bitLength() - 1];
	}

	Number _ret = ret.bits_to_number();
std::cout << "ret[0] = " << ret[0].to_int() << std::endl;
std::cout << "ret[1] = " << ret[1].to_int() << std::endl;
std::cout << "_ret = " << _ret.to_int() << std::endl;

	return ret.bits_to_number();
}

template<class OutNumber, class TempNumber>
void get_percentile(const std::vector<int> &_x, float percentile) {

	OutNumber avgValue;
	OutNumber avgSqrValue;

//	Polynomial<OutNumber> sqrtPoly = build_sqrt_polynomial<OutNumber>(OutNumber::global_p());

	n = _x.size();

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avg(_x.size());
	avg.compute_resample_constant(0.1, 0.05, *(std::max_element(_x.begin(), _x.end())));

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgSqrMsd(_x.size());
	avgSqrMsd.compute_resample_constant(0.1, 0.05, *(std::max_element(_x.begin(), _x.end())));
//	avgSqrMsd.set_m((n * avgValue.p()) >> accuracyEnhancementBits);

	AverageLiphe<OutNumber, TempNumber, CompareNative<TempNumber>, TruncConversion> avgSqrLsd(_x.size());
	avgSqrLsd.compute_resample_constant(0.1, 0.05, *(std::max_element(_x.begin(), _x.end())));


	for (int i = 0; i < _x.size(); ++i) {

		TempNumber xi(_x[i]);

		start_stat_interval();

		realAvg += _x[i];
		realAvgSqr += _x[i] * _x[i];

		avg.add(xi);
		OutNumber cost = xi.bits_to_number();
		avgSqrLsd.add_with_cost(xi, cost);
		
		OutNumber costShiftet = (xi >> accuracyEnhancementBits).bits_to_number();
		avgSqrMsd.add_with_cost(xi, costShiftet);

//		std::cout << "\n\n===================\n";
//		std::cout << "avg (raw) = " << avg.get_bits().to_bit_stream() << std::endl;
//		std::cout << "avg (bits) = " << avg.getAverage().to_bit_stream() << std::endl;
//		std::cout << "avg  = " << avg.getAverage().to_int() << std::endl;
//		std::cout << std::endl;
//		std::cout << "avgSqr (raw) = " << avgSqr.get_bits().to_bit_stream() << std::endl;
//		std::cout << "avgSqr (bits) = " << avgSqr.getAverage().to_bit_stream() << std::endl;

		end_stat_interval();

		if ((i % statistic_rate) == statistic_rate - 1) {
			print_stat(i);

			if (measureAccuracy) {
				OutNumber avgTempEnc = avg.getAverage();
				OutNumber avgSqrMsdTempEnc = avgSqrMsd.getAverage();
				OutNumber avgSqrLsdTempEnc = avgSqrLsd.getAverage();

				print_detailed_stats(i, avgTempEnc, avgSqrMsdTempEnc, avgSqrLsdTempEnc);
			}
		}
	}

	OutNumber avgEnc = avg.getAverage();
	OutNumber avgSqrMsdEnc = avgSqrMsd.getAverage();
	OutNumber avgSqrLsdEnc = avgSqrLsd.getAverage();

	print_stat(-1);
	if (measureAccuracy) {
		print_detailed_stats(-1, avgEnc, avgSqrMsdEnc, avgSqrLsdEnc);
	}

//	OutNumber sigma = sqrt_poly.compute(avgSqrEnc - avgEnc);

	std::cout << "avg = " << avgEnc.to_int() << std::endl;
	std::cout << "avgSqr = " << (avgSqrMsdEnc.to_int() * avgSqrMsdEnc.p() + avgSqrLsdEnc.to_int()) << std::endl;

//	OutNumber sigmaSqr = avgSqrEnc - avgEnc*avgEnc;
//	std::cout << "sigmaSqr = " << sigmaSqr.to_int() << std::endl;
//	OutNumber sigma = fhe_sqrt_approximation<OutNumber, TempNumber>(sigmaSqr);
	OutNumber sigma( sqrt( (avgSqrMsdEnc.to_int() * avgSqrMsdEnc.p() + avgSqrLsdEnc.to_int()) - avgEnc.to_int()*avgEnc.to_int() ) );

	std::cout << "sigma = " << sigma.to_int() << std::endl;


	OutNumber threshold = avgEnc;
	threshold += sigma;
	threshold += sigma;

	std::cout << "Threshold = " << threshold.to_int() << std::endl;

	TempNumber thresholdBits = threshold.template to_digits<TempNumber>();
std::cout << "converting to bits " << threshold.to_int() << " turned into " << thresholdBits.to_int() << std::endl;

//	if (percentile == 15.8) {
//		// take 1 sigma away from avg
//		threshold -= sigma;
//	} else if (percentile == 2.2) {
//		// take 2 sigma away from avg
//		threshold -= sigma * 2;
//	} else if (percentile == 0.1) {
//		// take 3 sigma away from avg
//		threshold -= sigma * 3;
//	}

//	TempNumber thresholdBits = get_bits<OutNumber, TempNumber>(threshold);

	PercentileIterator<OutNumber, TempNumber> sigmaDist(&rawData, thresholdBits);


	int i = 0;
	sigmaDist.begin();
	while (!sigmaDist.is_end()) {
		OutNumber x = *sigmaDist;
		std::cout << i << ") ";
		std::cout << "x = " << x.to_int() << "  " << std::endl;
		++sigmaDist;
		++i;
	}
	std::cout << std::endl;
//	report(sigmaDist, compactRepresentation)
}


#endif
