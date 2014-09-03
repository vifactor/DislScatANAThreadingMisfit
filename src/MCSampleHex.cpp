/*
 * MCSampleHex.cpp
 *
 *  Created on: 12 черв. 2013
 *      Author: kopp
 */

#include "MCSampleHex.h"

MCSampleHex::MCSampleHex(double thickness, double size)
{
	m_thickness = thickness; m_size = size;
}

MCSampleHex::~MCSampleHex()
{
	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		delete m_interfaces[i];
	}
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		delete m_layers[i];
	}
}

double MCSampleHex::T_threading(double r, double phi) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_layers.size(); ++i)
	{
		result += m_layers[i]->T(r, phi);
	}
	return result;
}

double MCSampleHex::T_misfit(double x, double y, double z1, double z2) const
{
	static double result;

	result = 0.0;
	for(size_t i = 0; i < m_interfaces.size(); ++i)
	{
		result += m_interfaces[i]->T(x, y, z1, z2);
	}
	return result;
}

void MCSampleHex::addThreadingLayer(double rho, double b_edge, double b_screw, double rc,
		double Qx, double Qz, double nu)
{
	m_layers.push_back(new ThreadingLayerHex(rho, b_edge, b_screw, rc, Qx, Qz, nu));
}

void MCSampleHex::addMisfitInterface(double rho, double bx, double by, double bz,
                                                 double Qx, double Qy, double Qz,
                                                 double nu, double d)
{
	m_interfaces.push_back(new AnalyticalMisfitInterfaceHex(rho, bx, by, bz, 
	                                                            Qx, Qy, Qz,
	                                                            nu, d));
}
