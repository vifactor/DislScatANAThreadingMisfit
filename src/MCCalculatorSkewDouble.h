/*
 * MCCalculatorSkewDouble.h
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#ifndef MCCALCULATORSKEWDOUBLE_H_
#define MCCALCULATORSKEWDOUBLE_H_

#include "MCSampleHex.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

class MCCalculatorSkewDouble
{
public:
	MCCalculatorSkewDouble(const MCSampleHex * sample, double alpha, double Phi,
			gsl_rng * rng, size_t nbCallsTotal = 10000);
	virtual ~MCCalculatorSkewDouble();

	void setResolution(double fwhm_qx, double fwhm_qz);
	virtual double I(const double q);
	virtual double I(const double q, double & error);

	friend double mc_skew_double_integrand(double * r, size_t dim, void * params);
protected:
	virtual double T(double z1, double z2) const;

	const MCSampleHex * m_sample;

	double m_q;
	double m_cotanPhi, m_sinPhi, m_alpha;
	double m_resol2_x, m_resol2_z;
	static const size_t m_integral_dimension = 2;

	double m_rl[m_integral_dimension];
	double m_ru[m_integral_dimension];

	gsl_rng * m_rng;
	gsl_monte_function m_function;
	size_t m_nbCallsTotal;

	gsl_monte_vegas_state * m_mc_vegas_state;
};

#endif /* MCCALCULATORSKEWDOUBLE_H_ */
