/*
 * ProgramSettings.cpp
 *
 *  Created on: 11 черв. 2013
 *      Author: kopp
 */

#include "ProgramSettings.h"

Range readRange(const libconfig::Setting& stg)
{
	Range range;

	range.m_min = stg[0][0];
	range.m_max = stg[0][1];
	range.m_sampling = stg[1];

	return range;
}

std::ostream& operator<<(std::ostream& out, const Range& range)
{
	out << "[" << range.m_min << ", " << range.m_max << "]:" << range.m_sampling;
	return out;
}

ProgramSettings::ProgramSettings()
{
}

ProgramSettings::~ProgramSettings()
{
}

void ProgramSettings::read(const std::string& cfgfile)
{
	libconfig::Config cfg;

	m_cfgfile = cfgfile;
	// Read the file. If there is an error, report it
	try
	{
		cfg.readFile(m_cfgfile.c_str());
		cfg.setAutoConvert(true);
		const libconfig::Setting& root = cfg.getRoot();

		readSampleSettings(root);
		readCalculatorSettings(root);
		readEngineSettings(root);
	} catch (const libconfig::FileIOException &fioex)
	{
		throw Exception(toString(fioex.what()) + " in\t" + cfgfile);
	} catch (const libconfig::ParseException &pex)
	{
		throw Exception(
				toString(pex.what()) + " in\t" + cfgfile + ":"
						+ toString(pex.getLine()) + " - "
						+ toString(pex.getError()));
	} catch (const libconfig::SettingNotFoundException &nfex)
	{
		throw Exception(
				toString(nfex.what()) + "\t" + toString(nfex.getPath())
						+ " in\t" + cfgfile);
	} catch (libconfig::SettingTypeException& tex)
	{
		throw Exception(
				toString(tex.what()) + "\t" + toString(tex.getPath()) + " in\t"
						+ cfgfile);
	}
}

void ProgramSettings::readCalculatorSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &calculator = root["Calculator"];

	/*reflection*/
	if(calculator["Q"].isArray() && calculator["Q"].getLength() == CalculatorSettings::HEXDIM)
	{
		for(int i = 0; i < CalculatorSettings::HEXDIM; ++i)
		{
			m_calculatorSettings.Q[i] = calculator["Q"][i];
		}
	}
	else
	{
		throw ProgramSettings::Exception(toString(calculator["Q"].getPath()));
	}
	/*check the property of hexagonal Miller indices*/
	if((m_calculatorSettings.Q[0] + m_calculatorSettings.Q[1] + m_calculatorSettings.Q[2]) != 0)
	{
		throw ProgramSettings::Exception(toString(calculator["Q"].getPath()));
	}

	/*X-ray wavelength*/
	m_calculatorSettings.lambda = calculator["lambda"];

	/*m_calculatorSettings.scale = calculator["scale"];
	m_calculatorSettings.background = calculator["background"];*/

	m_calculatorSettings.qresolX = calculator["resolution"]["x"];
	m_calculatorSettings.qresolZ = calculator["resolution"]["z"];
}

void ProgramSettings::readSampleSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &sample = root["Sample"];

	/*lattice parameters*/
	m_sampleSettings.a0 = sample["a0"];
	m_sampleSettings.c0 = sample["c0"];

	/*Poisson ratio*/
	m_sampleSettings.nu = sample["nu"];

	/*Sample sizes*/
	m_sampleSettings.thickness = sample["thickness"];
	m_sampleSettings.width = sample["width"];

	/*dislocation settings*/
	const libconfig::Setting &dislocations = sample["dislocations"];

	/*potential fit parameters*/
	if(dislocations.exists("misfit"))
	{
		const libconfig::Setting &misfit = dislocations["misfit"];
		readMisfitDislocations(misfit);
	}

	if(dislocations.exists("threading"))
	{
		const libconfig::Setting &threading = dislocations["threading"];
		readThreadingDislocations(threading);
	}
}

void ProgramSettings::readEngineSettings(const libconfig::Setting& root)
{
	const libconfig::Setting &engine = root["Engine"];

	m_engineSettings.m_geometry = defineGeometry(engine["geometry"]);

	if(m_engineSettings.m_geometry == EngineSettings::geomCOPLANAR)
	{
		const libconfig::Setting &settings = engine["coplanar_settings"];
		readCoplanarSettings(settings);
	}
	else if(m_engineSettings.m_geometry == EngineSettings::geomSKEW)
	{
		const libconfig::Setting &settings = engine["skew_settings"];
		readSkewSettings(settings);
	}
	else
	{
		throw Exception("Unknown geometry setup:\t" + toString(engine["geometry"].c_str()));
	}

	m_engineSettings.outfile = engine["outfile"].c_str();
}

void ProgramSettings::readCoplanarSettings(const libconfig::Setting& stg)
{
	m_engineSettings.qxRange = readRange(stg["qxrange"]);
	m_engineSettings.qzRange = readRange(stg["qzrange"]);
	m_calculatorSettings.nbMCCalls = stg["nbMCCalls"];
}

void ProgramSettings::readSkewSettings(const libconfig::Setting& stg)
{
	m_engineSettings.m_diffractometry = defineDiffractometry(stg["diffractometry"]);

	if(m_engineSettings.m_diffractometry == EngineSettings::diffDOUBLE)
	{
		readSkewDoubleSettings(stg["double_settings"]);
	}else if(m_engineSettings.m_diffractometry == EngineSettings::diffTRIPLE)
	{
		readSkewTripleSettings(stg["triple_settings"]);
	}
	else
	{
		throw Exception("Unknown diffractometry setup:\t" + toString(stg["diffractometry"].c_str()));
	}
}

ProgramSettings::EngineSettings::Diffractometry ProgramSettings::defineDiffractometry(const libconfig::Setting& stg)
{
	std::string diffractometry;

	diffractometry = stg.c_str();
	if (diffractometry.compare("DOUBLE") == 0)
	{
		return EngineSettings::diffDOUBLE;
	}
	else if (diffractometry.compare("TRIPLE") == 0)
	{
		return EngineSettings::diffTRIPLE;
	}
	else
	{
		return EngineSettings::diffUNKNOWN;
	}
}

ProgramSettings::EngineSettings::Geometry ProgramSettings::defineGeometry(const libconfig::Setting& stg)
{
	std::string diffractometry;

	diffractometry = stg.c_str();
	if (diffractometry.compare("COPLANAR") == 0)
	{
		return EngineSettings::geomCOPLANAR;
	}
	else if (diffractometry.compare("SKEW") == 0)
	{
		return EngineSettings::geomSKEW;
	}
	else
	{
		return EngineSettings::geomUNKNOWN;
	}
}

void ProgramSettings::print() const
{
	printEngineSettings();
	printSampleSettings();
	printCalculatorSettings();
}

void ProgramSettings::printCalculatorSettings() const
{
	std::cout << "---Calculator settings---" << std::endl;
	std::cout << "Reflection:\t[" << m_calculatorSettings.Q[0] << ", "
			<< m_calculatorSettings.Q[1] << ", " << m_calculatorSettings.Q[2]
			<< ", " << m_calculatorSettings.Q[3] << "]" << std::endl;
	std::cout << "X-ray wavelength:\t" << m_calculatorSettings.lambda << std::endl;
	std::cout << "Resolutions (dqx, dqz):\t" << m_calculatorSettings.qresolX
			<< "\t" << m_calculatorSettings.qresolZ << std::endl;
	//std::cout << "Intensity scale coefficient:\t" << m_calculatorSettings.scale << std::endl;
	//std::cout << "Intensity background:\t" << m_calculatorSettings.background << std::endl;
}

void ProgramSettings::printSampleSettings() const
{
	std::cout << "---Sample settings---" << std::endl;
	std::cout << "Lattice parameters: (a0, c0)\t" << m_sampleSettings.a0 << ", "
			<< m_sampleSettings.c0 << std::endl;
	std::cout << "Sample sizes (thickness width):\t" << m_sampleSettings.thickness << "\t"
			<< m_sampleSettings.width << std::endl;
	std::cout << "Poisson ratio:\t" << m_sampleSettings.nu << std::endl;

	std::cout << "Misfit dislocations:" << std::endl;
	printMisfitDislocations();

	std::cout << "Threading dislocations:" << std::endl;
	printThreadingDislocations();

}

void ProgramSettings::printEngineSettings() const
{
	std::cout << "---Engine settings---" << std::endl;

	switch(m_engineSettings.m_geometry)
	{
	case EngineSettings::geomCOPLANAR:
		std::cout << "Intensity is calculated in co-planar geometry" << std::endl;
		printCoplanarSettings();
		break;
	case EngineSettings::geomSKEW:
		std::cout << "Intensity is calculated in skew geometry" << std::endl;
		printSkewSettings();
		break;
	default:
		//never happens
		break;
	}

	std::cout << "Output basename:\t" << m_engineSettings.outfile << std::endl;
}

void ProgramSettings::printCoplanarSettings() const
{
	std::cout << "Qx range:\t" << m_engineSettings.qxRange << std::endl;
	std::cout << "Qz range:\t" << m_engineSettings.qzRange << std::endl;
	std::cout << "Max nb Monte-Carlo steps:\t" << m_calculatorSettings.nbMCCalls << std::endl;
}

void ProgramSettings::printSkewSettings() const
{
	switch(m_engineSettings.m_diffractometry)
	{
	case EngineSettings::diffDOUBLE:
		std::cout << "Double crystal diffractometry is used" << std::endl;
		std::cout << "Q range:\t" << m_engineSettings.qxRange << std::endl;
		if(m_calculatorSettings.method == CalculatorSettings::methodMC)
		{
			std::cout << "Max nb Monte-Carlo steps:\t" << m_calculatorSettings.nbMCCalls << std::endl;
		}
		else if(m_calculatorSettings.method == CalculatorSettings::methodREC)
		{
			std::cout << "Sampling of z range:\t" << m_calculatorSettings.sampling << std::endl;
		}
		break;
	case EngineSettings::diffTRIPLE:
		std::cout << "Triple crystal diffractometry is used" << std::endl;
		std::cout << "Qx range:\t" << m_engineSettings.qxRange << std::endl;
		std::cout << "Qz range:\t" << m_engineSettings.qzRange << std::endl;
		std::cout << "Max nb Monte-Carlo steps:\t" << m_calculatorSettings.nbMCCalls << std::endl;
		break;
	default:
		//never happens
		break;
	}
}

void ProgramSettings::readSkewDoubleSettings(const libconfig::Setting& stg)
{
	m_engineSettings.qxRange = readRange(stg["qrange"]);

	if(stg.exists("sampling"))
	{
		m_calculatorSettings.sampling = stg["sampling"];
		m_calculatorSettings.method = CalculatorSettings::methodREC;
	}else if(stg.exists("nbMCCalls"))
	{
		m_calculatorSettings.nbMCCalls = stg["nbMCCalls"];
		m_calculatorSettings.method = CalculatorSettings::methodMC;
	}
	else
	{
		m_calculatorSettings.sampling = 2500;
		m_calculatorSettings.method = CalculatorSettings::methodREC;
	}
}

void ProgramSettings::readSkewTripleSettings(const libconfig::Setting& stg)
{
	m_engineSettings.qxRange = readRange(stg["qxrange"]);
	m_engineSettings.qzRange = readRange(stg["qzrange"]);
	m_calculatorSettings.nbMCCalls = stg["nbMCCalls"];
}

void ProgramSettings::readMisfitDislocations(const libconfig::Setting& stg)
{
	size_t nbDislocations;

	if(stg.isList())
	{
		nbDislocations = stg.getLength();
		m_sampleSettings.misfit.resize(nbDislocations);
		for(size_t i = 0; i < nbDislocations; ++i)
		{
			m_sampleSettings.misfit[i].rho = stg[i]["rho"];
			m_sampleSettings.misfit[i].b_x = stg[i]["b_x"];
			m_sampleSettings.misfit[i].b_z = stg[i]["b_z"];
		}
	}
}

void ProgramSettings::readThreadingDislocations(const libconfig::Setting& stg)
{
	size_t nbDislocations;

	if(stg.isList())
	{
		nbDislocations = stg.getLength();
		m_sampleSettings.threading.resize(nbDislocations);
		for(size_t i = 0; i < nbDislocations; ++i)
		{
			m_sampleSettings.threading[i].rho = stg[i]["rho"];
			m_sampleSettings.threading[i].rc = stg[i]["rc"];
			m_sampleSettings.threading[i].b_edge = stg[i]["b_edge"];
			m_sampleSettings.threading[i].b_screw = stg[i]["b_screw"];
		}
	}
}

void ProgramSettings::printMisfitDislocations() const
{
	for (size_t i = 0; i < m_sampleSettings.misfit.size(); ++i)
	{
		std::cout << "\tSet:\t" << i << std::endl;
		std::cout << "\t\tBurgers vector (bx):\t"
				<< m_sampleSettings.misfit[i].b_x << std::endl;
		std::cout << "\t\tBurgers vector (bz):\t"
				<< m_sampleSettings.misfit[i].b_z << std::endl;
		std::cout << "\t\tDensity :\t" << m_sampleSettings.misfit[i].rho
				<< std::endl;
	}
}

void ProgramSettings::printThreadingDislocations() const
{
	for (size_t i = 0; i < m_sampleSettings.threading.size(); ++i)
	{
		std::cout << "\tSet:\t" << i << std::endl;
		std::cout << "\t\tBurgers vector (edge):\t"
				<< m_sampleSettings.threading[i].b_edge << std::endl;
		std::cout << "\t\tBurgers vector (screw):\t"
				<< m_sampleSettings.threading[i].b_screw << std::endl;
		std::cout << "\t\tDensity:\t" << m_sampleSettings.threading[i].rho
				<< std::endl;
		std::cout << "\t\tCorrelation radius:\t" << m_sampleSettings.threading[i].rc
				<< std::endl;
	}
}
