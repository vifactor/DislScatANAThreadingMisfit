//============================================================================
// Name        : DislScatANAMisfitHex.cpp
// Author      : Viktor Kopp
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include "ProgramSettings.h"
#include "MCCalculatorCoplanarTriple.h"
#include "MCCalculatorSkewDouble.h"
#include "RecCalculatorSkewDouble.h"
#include "MCCalculatorSkewTriple.h"

using namespace std;

int main()
{
	const int * Q;
    double Qperp, Qpar, Qx, Qy, Qz, Qnorm;
	double cosAlpha, sinAlpha, sinPhi, cosPhi, sinPsi, cosPsi, sinThetaB,
			cosThetaB, lambda;
	double omega;

    std::vector<double> qx, qz;
    ofstream fout;
    gsl_rng * rng;

    ProgramSettings programSettings;

	MCSampleHex  * sample;
	MCCalculatorCoplanarTriple * calculator_coplanar_triple;
	MCCalculatorSkewDouble * calculator_skew_double;
	RecCalculatorSkewDouble * calculator_rec_skew_double;
	MCCalculatorSkewTriple * calculator_skew_triple;

	try
	{
		programSettings.read("default.cfg");
		programSettings.print();
	}catch(const ProgramSettings::Exception& ex)
	{
		std::cout << ex.what() << std::endl;
		return -1;
	}
    
    //save settings file with new name
    boost::filesystem::path backupcfgfile;
    backupcfgfile = programSettings.getEngineSettings().outfile;
    backupcfgfile.replace_extension(".~cfg");
    boost::filesystem::copy_file(programSettings.getConfigfile(),
            backupcfgfile,
            boost::filesystem::copy_option::overwrite_if_exists);

	Q = programSettings.getCalculatorSettings().Q;
	Qperp = 2 * M_PI * sqrt(2.0 / 3 * (Q[0] * Q[0] + Q[1] * Q[1] + Q[2] * Q[2]))
	  / programSettings.getSampleSettings().a0;
	Qpar = 2 * M_PI * Q[3]
	                  / programSettings.getSampleSettings().c0;

	rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, 836);

	sample = new MCSampleHex(programSettings.getSampleSettings().thickness,
			programSettings.getSampleSettings().width);

	/*add threading layers*/
	for (size_t i = 0; i < programSettings.getSampleSettings().threading.size();
			++i)
	{
		sample->addThreadingLayer(
				programSettings.getSampleSettings().threading[i].rho * 1e-14,
				programSettings.getSampleSettings().threading[i].b_edge,
				programSettings.getSampleSettings().threading[i].b_screw,
				programSettings.getSampleSettings().threading[i].rc, Qperp, Qpar,
				programSettings.getSampleSettings().nu);
	}
	//std::cout << sample->T_threading(1, 0.926) << std::endl;
	//std::cout << sample->T_threading(10, 0.926) << std::endl;
	//std::cout << sample->T_threading(100, 0.926) << std::endl;

	if(programSettings.getEngineSettings().m_geometry == ProgramSettings::EngineSettings::geomCOPLANAR)
	{
		Qx = Qperp;
		Qy = 0.0;
		Qz = Qpar;
	    /*one misfit interfaces*/
	    for (size_t i = 0; i < programSettings.getSampleSettings().misfit.size();
			    ++i)
	    {
		    sample->addMisfitInterface(
				    programSettings.getSampleSettings().misfit[i].rho * 1e-7,
				    programSettings.getSampleSettings().misfit[i].b_x,
				    0.0,//by
				    programSettings.getSampleSettings().misfit[i].b_z,
				    Qx, Qy, Qz,
				    programSettings.getSampleSettings().nu,
				    programSettings.getSampleSettings().thickness);
	    }
		std::cout << "Q_lab:\t" << Qx << "\t" << Qy << "\t" << Qz << std::endl;
		calculator_coplanar_triple = new MCCalculatorCoplanarTriple(sample, rng,
				programSettings.getCalculatorSettings().nbMCCalls);
		calculator_coplanar_triple->setResolution(programSettings.getCalculatorSettings().qresolX,
				programSettings.getCalculatorSettings().qresolZ);

		programSettings.getEngineSettings().qxRange.toVector(qx);
		programSettings.getEngineSettings().qzRange.toVector(qz);

		fout.open(programSettings.getEngineSettings().outfile.c_str());
		fout << "#qx\tqz\tintens" << std::endl;
		for (size_t i = 0; i < qx.size(); ++i)
		{
			for (size_t j = 0; j < qz.size(); ++j)
			{
				cout << i << "\t" << j << endl;
				fout << qx[i] << "\t" << qz[j] << "\t"
						<< calculator_coplanar_triple->I(qx[i], qz[j]) << endl;
			}
		}
		fout.close();
		delete calculator_coplanar_triple;
	}
	if(programSettings.getEngineSettings().m_geometry == ProgramSettings::EngineSettings::geomSKEW)
	{
		lambda = programSettings.getCalculatorSettings().lambda;
		Qnorm = sqrt(Qperp * Qperp + Qpar * Qpar);
		sinThetaB = Qnorm * lambda / (4 * M_PI);
		cosThetaB = sqrt(1.0 - gsl_pow_2(sinThetaB));
		sinPsi = Qpar / Qnorm;
		sinPhi = sinPsi * sinThetaB;
		cosPsi = sqrt(1.0 - gsl_pow_2(sinPsi));
		cosPhi = sqrt(1.0 - gsl_pow_2(sinPhi));
		cosAlpha = sinThetaB * cosPsi / cosPhi;
		sinAlpha = sqrt(1.0 - gsl_pow_2(cosAlpha));

		if(programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffDOUBLE)
		{
			Qx = Qperp * cosAlpha;
			Qy = Qperp * sinAlpha;
			Qz = Qpar;

			std::cout << "Q_lab:\t" << Qx << "\t" << Qy << "\t" << Qz << std::endl;
            
     	    for (size_t i = 0; i < programSettings.getSampleSettings().misfit.size();
			    ++i)
	        {
		        sample->addMisfitInterface(
				        programSettings.getSampleSettings().misfit[i].rho * 1e-7,
				        programSettings.getSampleSettings().misfit[i].b_x,
				        0.0,//by
				        programSettings.getSampleSettings().misfit[i].b_z,
				        Qx, Qy, Qz,
				        programSettings.getSampleSettings().nu,
				        programSettings.getSampleSettings().thickness);
	        }
			programSettings.getEngineSettings().qxRange.toVector(qx);

			fout.open(programSettings.getEngineSettings().outfile.c_str());
			fout << "#q\tomega[deg]\tintens" << std::endl;
			if(programSettings.getCalculatorSettings().method == ProgramSettings::CalculatorSettings::methodMC)
			{

				calculator_skew_double = new MCCalculatorSkewDouble(sample,
						asin(sinAlpha), asin(sinPhi), rng,
						programSettings.getCalculatorSettings().nbMCCalls);

				calculator_skew_double->setResolution(
						programSettings.getCalculatorSettings().qresolX,
						programSettings.getCalculatorSettings().qresolZ);
				for (size_t i = 0; i < qx.size(); ++i)
				{
					fout << qx[i] << "\t"
							<< qx[i] / (Qnorm * cosThetaB) * 180.0 / M_PI
							<< "\t" << calculator_skew_double->I(qx[i]) << endl;
				}
				/*delete calculator_skew_double;*/
				delete calculator_skew_double;
			}
			if(programSettings.getCalculatorSettings().method == ProgramSettings::CalculatorSettings::methodREC)
			{
				calculator_rec_skew_double = new RecCalculatorSkewDouble(sample,
						asin(sinAlpha), asin(sinPhi),
						programSettings.getCalculatorSettings().sampling);

				calculator_rec_skew_double->setResolution(
						programSettings.getCalculatorSettings().qresolX,
						programSettings.getCalculatorSettings().qresolZ);
				for (size_t i = 0; i < qx.size(); ++i)
				{
					fout << qx[i] << "\t"
							<< qx[i] / (Qnorm * cosThetaB) * 180.0 / M_PI
							<< "\t" << calculator_rec_skew_double->I(qx[i])
							<< endl;
				}
				/*delete calculator_skew_double;*/
				delete calculator_rec_skew_double;
			}
			fout.close();

		}
		if(programSettings.getEngineSettings().m_diffractometry == ProgramSettings::EngineSettings::diffTRIPLE)
		{
			Qx = Qperp;
			Qy = 0.0;
			Qz = Qpar;
                        
            for (size_t i = 0; i < programSettings.getSampleSettings().misfit.size();
			    ++i)
	        {
		        sample->addMisfitInterface(
				        programSettings.getSampleSettings().misfit[i].rho * 1e-7,
				        programSettings.getSampleSettings().misfit[i].b_x,
				        0.0,//by
				        programSettings.getSampleSettings().misfit[i].b_z,
				        Qx, Qy, Qz,
				        programSettings.getSampleSettings().nu,
				        programSettings.getSampleSettings().thickness);
	        }

			calculator_skew_triple = new MCCalculatorSkewTriple(sample,
					asin(sinPsi), rng,
					programSettings.getCalculatorSettings().nbMCCalls);
			calculator_skew_triple->setResolution(
					programSettings.getCalculatorSettings().qresolX,
					programSettings.getCalculatorSettings().qresolZ);

			/*TEMPORARY MODIFIED*/
			/*programSettings.getEngineSettings().qxRange.toVector(qx);
			programSettings.getEngineSettings().qzRange.toVector(qz);*/

			fout.open(programSettings.getEngineSettings().outfile.c_str());
			fout << "#om\tqx\tqz\tintens" << std::endl;

			/*for (size_t i = 0; i < qx.size(); ++i)
			{
				omega = 180.0 / M_PI * qx[i] / Qnorm;
				for (size_t j = 0; j < qz.size(); ++j)
				{
					cout << i << "\t" << j << endl;
					fout << omega << "\t" << qx[i] << "\t" << qz.at(j) << "\t" << calculator_skew_triple->I(qx[i], qz[j])
							<< endl;
				}
			}*/
			double x, z, omega_rad;
			for(size_t i = 0; i < 100; ++i)
			{
				//omega in degrees
				omega = -1.0 + i * 2.0 /100;
				omega_rad = M_PI * omega / 180;
				x = Qnorm * sin(omega_rad);
				z = Qnorm * (cos(omega_rad) - 1);
				cout << i << "\t" << x << "\t" << z << endl;
				fout << omega << "\t" << x << "\t" << z << "\t" << calculator_skew_triple->I(x, z)
						<< endl;
			}

			fout.close();
			delete calculator_skew_triple;
		}
	}

	gsl_rng_free(rng);
	delete sample;

	std::cout << "Done." << std::endl;

	return EXIT_SUCCESS;
}
