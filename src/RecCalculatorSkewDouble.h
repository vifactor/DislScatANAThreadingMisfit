/*
 * RecCalculatorSkewDouble.h
 *
 *  Created on: 17 june 2013
 *      Author: kopp
 */

#ifndef RECCALCULATORSKEWDOUBLE_H_
#define RECCALCULATORSKEWDOUBLE_H_

#include "MCSampleHex.h"

class RecCalculatorSkewDouble
{
public:
	RecCalculatorSkewDouble(const MCSampleHex * sample, double alpha, double Phi,
			size_t z_sampling = 100);
	virtual ~RecCalculatorSkewDouble();

	void setResolution(double fwhm_qx, double fwhm_qz);
	virtual double I(const double q);
protected:
	virtual double T(double z1, double z2) const;

	const MCSampleHex * m_sample;

	double m_cotanPhi, m_sinPhi, m_alpha;
	double m_resol2_x, m_resol2_z;

	/*values of correlation function in corresponding points*/
	size_t m_z_sampling;
	double ** m_G;
	double *m_z1, *m_z2;
	void setup();
};

#endif /* RECCALCULATORSKEWDOUBLE_H_ */
