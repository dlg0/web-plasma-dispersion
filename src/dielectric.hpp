#ifndef DIELECTRIC_HPP_
#define DIELECTRIC_HPP_
#include "stixVars.hpp"
#include <complex>
#include <armadillo>
#include "rotation.hpp"

class dielectric
{
		public:
				dielectric(StixVars _s);
				dielectric(std::vector<HotPlasmaSpecies> species, std::complex<double> _omega_c, int _l, 
								std::complex<double> _ky, std::complex<double> _kz, RotationMatrix _R);

		private:
				std::complex<double> I;
		public:
				arma::cx_mat epsilon;

				RotationMatrix R;
				std::vector< std::complex<double> > roots;

				std::complex<double> K0,K1,K2,K3,K4,K5;
				std::complex<double> Ka0,Ka1,Ka2,Ka3,Ka4,Ka5;

				arma::cx_colvec k_xyz, k_abp;
				double omega;
				std::complex<double> omega_c;
				int l;
				std::vector<HotPlasmaSpecies> species;

				void rotate(arma::mat _rot);
				void rotate(arma::cx_mat _rot);
				//void determinant(double kxIn_re, double kxIn_im, double *kxOu_re, double *kxOu_im);
				std::complex<double> determinant(std::complex<double> _kx);
				std::complex<double> determinant(std::complex<double> _kx, std::complex<double> _ky, std::complex<double> _kz, RotationMatrix _R, double _omega, int _coldVersion);
				int populateSwansonKs(std::complex<double> _kx);
				//void coldRoots(double, double, double);
};

#endif
