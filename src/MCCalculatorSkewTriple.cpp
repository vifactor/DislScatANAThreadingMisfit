/*
 * MCCalculatorSkewTriple.cpp
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#include "MCCalculatorSkewTriple.h"

double mc_skew_triple_integrand(double * r, size_t dim, void * params)
{
	MCCalculatorSkewTriple * calculator;
	static double result, y, z1, z2;

	calculator = static_cast<MCCalculatorSkewTriple *> (params);

	y = r[0];
	z1 = r[1];
	z2 = r[2];

	result = exp(-calculator->T(y, z1, z2))
			* cos(calculator->m_qx * y + calculator->m_qz * (z1 - z2));

	return result;
}

MCCalculatorSkewTriple::MCCalculatorSkewTriple(const MCSampleHex * sample, double Psi, gsl_rng * rng, size_t nbCallsTotal)
{
	m_sample = sample;
	m_rng = rng;
	m_nbCallsTotal 	= nbCallsTotal;


	m_function.f = mc_skew_triple_integrand;
	m_function.params = this;
	m_function.dim = m_integral_dimension;

	m_mc_vegas_state = gsl_monte_vegas_alloc(m_integral_dimension);
	gsl_monte_vegas_init(m_mc_vegas_state);

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

	m_sinPsi = sin(Psi);
	m_cotanPsi = cos(Psi) / sin(Psi);

	//std::cout << "sinPsi:\t" << m_sinPsi << std::endl;
	//std::cout << "cotanPsi:\t" << m_cotanPsi << std::endl;
}

MCCalculatorSkewTriple::~MCCalculatorSkewTriple()
{
	gsl_monte_vegas_free (m_mc_vegas_state);
}

double MCCalculatorSkewTriple::T(double y, double z1, double z2) const
{
	static double x, z, r, phi;
	static double T_threading, T_misfit, T_resolution;
	static double result;

	/*arguments are correct for triple skew geometry*/
	z = z1 - z2;

	if(fabs(m_cotanPsi) < 1e-10)//symmetric reflection
	{
		x = 0;
		r = fabs(y);
		phi = M_PI / 2;
	}
	else
	{
		x = m_cotanPsi * z;
		r = sqrt(x * x + y * y);
		phi = atan2(y, x);
	}

	T_resolution = m_resol2_x * gsl_pow_2(r) + m_resol2_z * gsl_pow_2(z);
	T_threading = m_sample->T_threading(r, phi);
	T_misfit = m_sample->T_misfit(x, y, z1, z2);

	/*std::cout << "x:\t" << x << std::endl;
	std::cout << "y:\t" << y << std::endl;
	std::cout << "z:\t" << z << std::endl;
	std::cout << "z1:\t" << z1 << std::endl;
	std::cout << "z2:\t" << z2 << std::endl;*/
	/*std::cout << "phi:\t" << phi << std::endl;*/
	/*std::cout << "threading:\t" << T_threading << std::endl;
	std::cout << "misfit:\t" << T_misfit << std::endl;
	std::cout << "resolution:\t" << T_resolution << std::endl;*/


	result = T_threading + T_misfit + T_resolution;
	/*std::cout << "G:\t" << exp(- result) << std::endl;*/

	return result;
}

void MCCalculatorSkewTriple::setResolution(double fwhm_qx, double fwhm_qz)
{
	/* FWHM = 2 * sqrt(2 * log(2)) / sigma */
	/* resol2 = 1 / (2 * sigma**2) */
	m_resol2_x = gsl_pow_2(fwhm_qx / 4) / log(2.0);
	m_resol2_z = gsl_pow_2(fwhm_qz / 4) / log(2.0);

	/*std::ofstream fout("test_T.txt");
	double stepy, stepz;
	size_t ny, nz;
	double z1, z2;

	ny = 100; nz = 100;
	stepy = m_sample->m_size / (ny - 1);
	stepz = m_sample->m_thickness / (nz - 1);
	z1 = 10;
	for(double z = -m_sample->m_thickness; z < m_sample->m_thickness; z += stepz)
	{
		for(double y = -m_sample->m_size; y < m_sample->m_size; y += stepy)
		{
			z2 = z1 - z;
			fout << y << "\t" << z << "\t" << exp(-T(y, z1, z2)) << std::endl;
		}
	}
	fout.close();*/
}

double MCCalculatorSkewTriple::I(const double qx, const double qz, double & err)
{
	static double result, abserr;

	m_qx = qx;
	m_qz = qz / m_sinPsi;

	result = 0.0;
	gsl_monte_vegas_init(m_mc_vegas_state);

	/*int gsl_monte_vegas_integrate (gsl_monte_function * f, const double xl[], const double xu[], size_t dim,
				size_t calls, gsl_rng * r, gsl_monte_vegas_state * s, double * result, double * abserr)*/
	gsl_monte_vegas_integrate(&m_function, m_rl, m_ru, m_integral_dimension,
			m_nbCallsTotal, m_rng, m_mc_vegas_state, &result, &abserr);

	err = abserr;
	return result;
}

double MCCalculatorSkewTriple::I(double qx, double qz)
{
	static double result, abserr;

	result = I(qx, qz, abserr);

	return result;
}
