#ifndef __UTILLITES_H__
#define __UTILLITES_H__

#include <assert.h>

typedef unsigned int uint32_t;


class Data
{
	typedef const double cdouble;

public:
	Data() :dt(0.0), v_ref(0.0), Lf(0.0), N(0.0) {}

	Data(const std::vector<uint32_t>&_starts, const std::vector<double>& _costs, cdouble &_dt, cdouble&_v_ref, cdouble &_Lf, cdouble &_N) :
		 starts(_starts),costs(_costs), dt(_dt), v_ref(_v_ref), Lf(_Lf), N(_N) {}

	Data& operator=(const Data& data)
	{
		costs	= data.costs;
		starts	= data.starts;
		dt		= data.dt;
		v_ref	= data.v_ref;
		Lf		= data.Lf;
		N		= data.N;

		return *this;
	}

	std::vector<double>		costs;
	std::vector<uint32_t>	starts;
	double					dt;
	double					v_ref;
	double					Lf;
	uint32_t 				N;
};


static double deg2rad(const double &deg)
{
	return deg * 0.01745329251994329576923690768489;
}
static double polyeval(const Eigen::VectorXd &coeffs, const double &x)
{
	double result = 0.0;
	for (uint32_t i = 0; i < coeffs.size(); ++i)
		result += coeffs[i] * std::pow(x, i);

	return result;
}
static Eigen::VectorXd polyfit(const Eigen::VectorXd &xvals, const Eigen::VectorXd& yvals, const uint32_t &order)
{
	assert(xvals.size() == yvals.size());
	assert(order >= 1 && order <= xvals.size() - 1);

	Eigen::MatrixXd A(xvals.size(), order + 1);

	for (uint32_t i = 0; i < xvals.size(); ++i)
		A(i, 0) = 1.0;

	for (uint32_t j = 0; j < xvals.size(); ++j)
		for (uint32_t i = 0; i < order; ++i)
			A(j, i + 1) = A(j, i) * xvals(j);


	return A.householderQr().solve(yvals);;
}
#endif