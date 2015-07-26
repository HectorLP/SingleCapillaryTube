#ifndef TOLERANCE_H_
#define TOLERANCE_H_

#include <cmath>

class tolerance {
public:
	tolerance(double eps) :
		_eps(eps) {
	}
	bool operator()(double a, double b) {
		return (fabs(b - a) <= _eps);
	}
private:
	double _eps;
};
#endif