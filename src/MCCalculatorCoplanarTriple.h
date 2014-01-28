/*
 * MCCalculatorCoplanarTriple.h
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#ifndef MCCALCULATORCOPLANARTRIPLE_H_
#define MCCALCULATORCOPLANARTRIPLE_H_

#include "MCSampleHex.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

class MCCalculatorCoplanarTriple
{
public:
	MCCalculatorCoplanarTriple(const MCSampleHex * sample, gsl_rng * rng, size_t nbCallsTotal = 10000);
	virtual ~MCCalculatorCoplanarTriple();

	virtual void setResolution(double fwhm_qx, double fwhm_qz);
	virtual double I(const double qx, const double qz);
	virtual double I(const double qx, const double qz, double & error);

	friend double mc_coplanar_triple_integrand(double * r, size_t dim, void * params);
protected:
	virtual double T(double x, double z1, double z2) const;

	const MCSampleHex * m_sample;

	double m_qx, m_qz;
	double m_resol2_x, m_resol2_z;
	static const size_t m_integral_dimension = 3;

	double m_rl[m_integral_dimension];
	double m_ru[m_integral_dimension];

	gsl_rng * m_rng;
	gsl_monte_function m_function;
	size_t m_nbCallsTotal;

	gsl_monte_vegas_state * m_mc_vegas_state;
};

#endif /* MCCALCULATORCOPLANARTRIPLE_H_ */
