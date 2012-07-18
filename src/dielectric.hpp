#ifndef DIELECTRIC_HPP_
#define DIELECTRIC_HPP_
#include "stixVars.hpp"
#include <complex>
#include <armadillo>

class dielectric
{
		public:
				dielectric(StixVars);
				dielectric(std::vector<HotPlasmaSpecies>, double,int,arma::cx_colvec,std::vector<float>);

		private:
				std::complex<double> I;
		public:
				arma::cx_mat stix;
				arma::cx_mat stixRotated;
				arma::mat rotQ;
				std::vector< std::complex<double> > roots;

				void rotateEpsilon(std::vector<float> b);
				void coldRoots(double, double, double);

};

#endif
