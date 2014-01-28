/*
 * MCSampleHex.h
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#ifndef MCSAMPLEHEX_H_
#define MCSAMPLEHEX_H_

#include "MisfitInterfaceHex.h"
#include "ThreadingLayerHex.h"
#include <vector>

class MCSampleHex
{
public:
	MCSampleHex(double thickness, double size);
	~MCSampleHex();
	void setQ(double Qx, double Qy, double Qz);

	void addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
			double Qx, double Qz, double nu);
	void addMisfitInterface(double rho, double bx, double bz, double nu, double d);
	double T_threading(double r, double phi) const;
	double T_misfit(double x, double y, double z1, double z2) const;
	double m_thickness;
	double m_size;
protected:
	std::vector<MisfitInterfaceHex * > m_interfaces;
	std::vector<ThreadingLayerHex * > m_layers;
};

#endif /* MCSAMPLEHEX_H_ */
