#ifndef DIELECTRIC_HPP_
#define DIELECTRIC_HPP_
#include "stixVars.hpp"
#include <complex>
#include <armadillo>

class dielectric
{
		public:
				dielectric(StixVars);
				dielectric(std::vector<HotPlasmaSpecies>, double,int,arma::cx_colvec);

		private:
				std::complex<double> I;
		public:
				arma::cx_mat epsilon;

				arma::mat rotQ;
				std::vector< std::complex<double> > roots;

				void rotate(arma::mat _rot);
				void rotate(arma::cx_mat _rot);

				void coldRoots(double, double, double);
};

#endif
