/*
 * MCCalculatorSkewDouble.cpp
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#include "MCCalculatorSkewDouble.h"

double mc_skew_double_integrand(double * r, size_t dim, void * params)
{
	MCCalculatorSkewDouble * calculator;
	static double result, z1, z2;

	calculator = static_cast<MCCalculatorSkewDouble*> (params);

	z1 = r[0];
	z2 = r[1];

	/*z2 = z1 - z;*/

	result = exp(-calculator->T(z1, z2))
			* cos(calculator->m_q * (z1 - z2));

	return result;
}

MCCalculatorSkewDouble::MCCalculatorSkewDouble(const MCSampleHex * sample,
		double alpha, double Phi, gsl_rng * rng, size_t nbCallsTotal)
{
	m_sample = sample;
	m_rng = rng;
	m_nbCallsTotal 	= nbCallsTotal;

	m_function.f = mc_skew_double_integrand;
	m_function.params = this;
	m_function.dim = m_integral_dimension;

	m_mc_vegas_state = gsl_monte_vegas_alloc(m_integral_dimension);

	m_q = 0.0;

	m_resol2_x = 0.0;
	m_resol2_z = 0.0;

	m_rl[0] = 0.0;
	m_ru[0] = sample->m_thickness;

	m_rl[1] = 0.0;
	m_ru[1] = sample->m_thickness;

	/*m_rl[1] = -sample->m_size;
	m_ru[1] = sample->m_size;*/

	gsl_monte_vegas_init(m_mc_vegas_state);

	m_alpha = alpha;
	m_sinPhi = sin(Phi);
	m_cotanPhi = cos(Phi) / sin(Phi);
}

MCCalculatorSkewDouble::~MCCalculatorSkewDouble()
{
	gsl_monte_vegas_free (m_mc_vegas_state);
}

void MCCalculatorSkewDouble::setResolution(double fwhm_qx, double fwhm_qz)
{
	/* FWHM = 2 * sqrt(2 * log(2)) / sigma */
	/* resol2 = 1 / (2 * sigma**2) */
	m_resol2_x = gsl_pow_2(fwhm_qx / 4) / log(2.0);
	m_resol2_z = gsl_pow_2(fwhm_qz / 4) / log(2.0);
}

double MCCalculatorSkewDouble::I(const double q, double & err)
{
	static double result, abserr;

	m_q = q / m_sinPhi;

	result = 0.0;
	//gsl_monte_vegas_init(m_mc_vegas_state);

	/*int gsl_monte_vegas_integrate (gsl_monte_function * f, const double xl[], const double xu[], size_t dim,
				size_t calls, gsl_rng * r, gsl_monte_plain_state * s, double * result, double * abserr)*/
	gsl_monte_vegas_integrate(&m_function, m_rl, m_ru, m_integral_dimension,
			m_nbCallsTotal, m_rng, m_mc_vegas_state, &result, &abserr);

	err = abserr;
	return result;
}

double MCCalculatorSkewDouble::I(const double q)
{
	static double result, abserr;

	result = I(q, abserr);

	return result;
}

double MCCalculatorSkewDouble::T(double z1, double z2) const
{
	static double x, z;
	static double T_threading, T_misfit, T_resolution;
	static double result;

	z = z1 - z2;
	x = m_cotanPhi * z;
	T_resolution = m_resol2_x * gsl_pow_2(x) + m_resol2_z * gsl_pow_2(z);
	T_threading = m_sample->T_threading(fabs(x), m_alpha);
	/*initially x is along proj(Kout), whereas x in T_misfit is supposed to be along proj(Q)*/
	T_misfit = m_sample->T_misfit(x * cos(m_alpha), 0.0, z1, z2);

	result = T_threading + T_misfit + T_resolution;
	return result;
}

