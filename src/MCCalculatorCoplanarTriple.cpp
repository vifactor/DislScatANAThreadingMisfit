/*
 * MCCalculatorCoplanarTriple.cpp
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#include "MCCalculatorCoplanarTriple.h"

double mc_coplanar_triple_integrand(double * r, size_t dim, void * params)
{
	MCCalculatorCoplanarTriple * calculator;
	static double result, x, z1, z2;

	calculator = static_cast<MCCalculatorCoplanarTriple *> (params);

	x = r[0];
	z1 = r[1];
	z2 = r[2];

	result = exp(-calculator->T(x, z1, z2))
			* cos(calculator->m_qx * x + calculator->m_qz * (z1 - z2));

	return result;
}

MCCalculatorCoplanarTriple::MCCalculatorCoplanarTriple(
		const MCSampleHex * sample, gsl_rng * rng, size_t nbCallsTotal)
{

	m_sample = sample;
	m_rng = rng;
	m_nbCallsTotal 	= nbCallsTotal;


	m_function.f = mc_coplanar_triple_integrand;
	m_function.params = this;
	m_function.dim = m_integral_dimension;

	m_mc_vegas_state = gsl_monte_vegas_alloc(m_integral_dimension);

	m_qx = 0.0;
	m_qz = 0.0;

	m_resol2_x = 0.0;
	m_resol2_z = 0.0;

	m_rl[0] = -sample->m_size / 2;
	m_ru[0] = sample->m_size / 2;

	m_rl[1] = 0.0;
	m_ru[1] = sample->m_thickness;

	m_rl[2] = 0.0;
	m_ru[2] = sample->m_thickness;

	gsl_monte_vegas_init(m_mc_vegas_state);
}

MCCalculatorCoplanarTriple::~MCCalculatorCoplanarTriple()
{
	gsl_monte_vegas_free (m_mc_vegas_state);
}

double MCCalculatorCoplanarTriple::T(double x, double z1, double z2) const
{
	static double T_threading, T_misfit, T_resolution;
	static double result;

	/*arguments are correct for coplanar geometry*/
	T_threading = m_sample->T_threading(fabs(x), 0.0);

	T_misfit = m_sample->T_misfit(x, 0.0, z1, z2);

	T_resolution = m_resol2_x * x * x + m_resol2_z * gsl_pow_2(z1 - z2);

	result = T_threading + T_misfit + T_resolution;

	return result;
}

void MCCalculatorCoplanarTriple::setResolution(double fwhm_qx, double fwhm_qz)
{
	/* FWHM = 2 * sqrt(2 * log(2)) / sigma */
	/* resol2 = 1 / (2 * sigma**2) */
	m_resol2_x = gsl_pow_2(fwhm_qx / 4) / log(2.0);
	m_resol2_z = gsl_pow_2(fwhm_qz / 4) / log(2.0);
}

double MCCalculatorCoplanarTriple::I(const double qx, const double qz, double & err)
{
	static double result, abserr;

	m_qx = qx;
	m_qz = qz;

	result = 0.0;
	//gsl_monte_vegas_init(m_mc_vegas_state);

	/*int gsl_monte_vegas_integrate (gsl_monte_function * f, const double xl[], const double xu[], size_t dim,
				size_t calls, gsl_rng * r, gsl_monte_plain_state * s, double * result, double * abserr)*/
	gsl_monte_vegas_integrate(&m_function, m_rl, m_ru, m_integral_dimension,
			m_nbCallsTotal, m_rng, m_mc_vegas_state, &result, &abserr);

	err = abserr;
	return result;
}

double MCCalculatorCoplanarTriple::I(const double qx, const double qz)
{
	static double result, abserr;

	result = I(qx, qz, abserr);

	return result;
}

