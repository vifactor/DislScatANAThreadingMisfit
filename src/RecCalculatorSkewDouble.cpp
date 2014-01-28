/*
 * RecCalculatorSkewDouble.cpp
 *
 *  Created on: 17 черв. 2013
 *      Author: kopp
 */

#include "RecCalculatorSkewDouble.h"

RecCalculatorSkewDouble::RecCalculatorSkewDouble(const MCSampleHex * sample, double alpha, double Phi,
		size_t z_sampling)
{
	m_sample = sample;

	m_resol2_x = 0.0;
	m_resol2_z = 0.0;

	m_alpha = alpha;
	m_sinPhi = sin(Phi);
	m_cotanPhi = cos(Phi) / sin(Phi);

	m_z_sampling = z_sampling;
	m_z1 = new double[m_z_sampling];
	m_z2 = new double[m_z_sampling];
	m_G = new double * [m_z_sampling];
	for(size_t i = 0; i < m_z_sampling; ++i)
	{
		m_G[i] = new double[m_z_sampling];
	}
}

RecCalculatorSkewDouble::~RecCalculatorSkewDouble()
{
	delete[] m_z1;
	delete[]  m_z2;
	for(size_t i = 0; i < m_z_sampling; ++i)
	{
		delete[] m_G[i];
	}
	delete[] m_G;
}

void RecCalculatorSkewDouble::setResolution(double fwhm_qx, double fwhm_qz)
{
	/* FWHM = 2 * sqrt(2 * log(2)) / sigma */
	/* resol2 = 1 / (2 * sigma**2) */
	m_resol2_x = gsl_pow_2(fwhm_qx / 4) / log(2.0);
	m_resol2_z = gsl_pow_2(fwhm_qz / 4) / log(2.0);

	setup();
}

double RecCalculatorSkewDouble::T(double z1, double z2) const
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

void RecCalculatorSkewDouble::setup()
{
	double dz;

	dz = m_sample->m_thickness / m_z_sampling;
	for(size_t i = 0; i < m_z_sampling; ++i)
	{
		m_z1[i] = dz / 2 + i * dz;
		for(size_t j = 0; j < m_z_sampling; ++j)
		{
			m_z2[j] = dz / 2 + j * dz;

			/*important that z1 != z2, otherwise T(z1, z2) == NaN*/
			if(m_z2[j] == m_z1[i])
			{
				m_z2[j] += dz / m_z_sampling;
			}
			m_G[i][j] = exp(-T(m_z1[i], m_z2[j])) * dz * dz;

			//std::cout << "G("<< m_z1[i] <<", " << m_z2[j] <<") = " << m_G[i][j] << std::endl;
		}
	}
}

double RecCalculatorSkewDouble::I(const double q)
{
	static double result, z;
	result = 0.0;

	for(size_t i = 0.0; i < m_z_sampling; ++i)
	{
		for(size_t j = 0.0; j < m_z_sampling; ++j)
		{
			z = m_z1[i] - m_z2[j];

			result += m_G[i][j] * cos(z / m_sinPhi * q);
		}
	}
	return result;
}
