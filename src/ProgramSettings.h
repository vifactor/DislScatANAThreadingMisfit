/*
 * ProgramSettings.h
 *
 *  Created on: 11 june 2013
 *      Author: kopp
 */

#ifndef PROGRAMSETTINGS_H_
#define PROGRAMSETTINGS_H_

#include "StringTools.h"
#include <libconfig.h++>

struct Range
{
	double m_min;
	double m_max;
	size_t m_sampling;

	double getStep() const
	{
		return (m_sampling > 1) ? (m_max - m_min) / (m_sampling - 1) : 0.0;
	}
	void toVector(std::vector<double>& vec) const
	{
		double step = getStep();
		for (size_t i = 0; i < m_sampling; ++i)
		{
			vec.push_back(m_min + step * i);
		}
	}
	bool good() const
	{
		if ((m_min <= m_max) && (m_sampling > 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};

class ProgramSettings
{
public:
	class Exception: public std::exception
	{
	public:
		Exception(std::string m)
		{
			msg = "ProgramSettings::" + m;
		}
		~Exception() throw ()
		{
		}
		const char* what() const throw ()
		{
			return msg.c_str();
		}
	private:
		std::string msg;
	};
	struct SampleSettings
	{
		double nu;
		double thickness;
		double width;
		/*hexagonal lattice parameters*/
		double a0, c0;

		struct MisfitDislocationType
		{
			/*burgers components*/
			double b_x, b_z;
			double rho;
		};
		struct ThreadingDislocationType
		{
			double b_edge, b_screw;
			double rho;
			double rc;
		};

		std::vector<MisfitDislocationType> misfit;
		std::vector<ThreadingDislocationType> threading;
	};
	struct CalculatorSettings
	{
		enum {HEXDIM = 4};
		/*hexagonal index*/
		int Q[HEXDIM];
		/*X-ray wavelength*/
		double lambda;

		double qresolX, qresolZ;
		double scale;
		double background;

		union
		{
			long int nbMCCalls;
			long int sampling;
		};

		enum {methodMC, methodREC} method;
	};
	struct EngineSettings
	{
		std::string outfile;
		Range qxRange, qzRange;

		enum Geometry {geomSKEW, geomCOPLANAR, geomUNKNOWN} m_geometry;
		enum Diffractometry {diffDOUBLE, diffTRIPLE, diffUNKNOWN} m_diffractometry;
	};
	const SampleSettings& getSampleSettings() const
	{
		return m_sampleSettings;
	}
	const CalculatorSettings& getCalculatorSettings() const
	{
		return m_calculatorSettings;
	}
	const EngineSettings& getEngineSettings() const
	{
		return m_engineSettings;
	}
	const std::string& getConfigfile() const
	{
		return m_cfgfile;
	}
	ProgramSettings();
	virtual ~ProgramSettings();

	void read(const std::string& cfg);
	void print() const;
protected:
	void readCalculatorSettings(const libconfig::Setting& root);
	void readSampleSettings(const libconfig::Setting& root);
	void readEngineSettings(const libconfig::Setting& root);

	void readCoplanarSettings(const libconfig::Setting& stg);
	void readSkewSettings(const libconfig::Setting& stg);
	void readSkewDoubleSettings(const libconfig::Setting& stg);
	void readSkewTripleSettings(const libconfig::Setting& stg);

	void readMisfitDislocations(const libconfig::Setting& stg);
	void readThreadingDislocations(const libconfig::Setting& stg);

	EngineSettings::Diffractometry defineDiffractometry(const libconfig::Setting& stg);
	EngineSettings::Geometry defineGeometry(const libconfig::Setting& stg);

	void printSampleSettings() const;
	void printCalculatorSettings() const;
	void printEngineSettings() const;

	void printCoplanarSettings() const;
	void printSkewSettings() const;
	void printSkewDoubleSettings() const;
	void printSkewTripleSettings() const;

	void printMisfitDislocations() const;
	void printThreadingDislocations() const;

	SampleSettings m_sampleSettings;
	CalculatorSettings m_calculatorSettings;
	EngineSettings m_engineSettings;
	std::string m_cfgfile;
};

#endif /* PROGRAMSETTINGS_H_ */
