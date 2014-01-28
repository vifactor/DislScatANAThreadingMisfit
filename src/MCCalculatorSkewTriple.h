/*
 * MCCalculatorSkewTriple.h
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#ifndef MCCALCULATORSKEWTRIPLE_H_
#define MCCALCULATORSKEWTRIPLE_H_

#include <fstream>
#include "MCSampleHex.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

class MCCalculatorSkewTriple
{
public:
	MCCalculatorSkewTriple(const MCSampleHex * sample, double Psi, gsl_rng * rng, size_t nbCallsTotal = 10000);
	virtual ~MCCalculatorSkewTriple();

	void setResolution(double fwhm_qx, double fwhm_qz);
	double I(double qx, double qz);
	double I(double qx, double qz, double & error);

	friend double mc_skew_triple_integrand(double * r, size_t dim, void * params);
protected:
	double T(double y, double z1, double z2) const;

	const MCSampleHex * m_sample;

	double m_qx, m_qz;
	double m_cotanPsi, m_sinPsi;
	double m_resol2_x, m_resol2_z;
	static const size_t m_integral_dimension = 3;

	double m_rl[m_integral_dimension];
	double m_ru[m_integral_dimension];

	gsl_rng * m_rng;
	gsl_monte_function m_function;
	size_t m_nbCallsTotal;

	gsl_monte_vegas_state * m_mc_vegas_state;
};

#endif /* MCCALCULATORSKEWTRIPLE_H_ */
