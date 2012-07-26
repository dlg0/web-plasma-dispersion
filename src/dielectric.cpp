#include "dielectric.hpp"
#include "constants.hpp"
#include <math.h>
#include "arrays.hpp"
#include "stixVars.hpp"
#include <complex>
#include <vector>
#include <numeric>
#include <armadillo>
#include <o2scl/poly.h>
#include "matpack.h"
#include "rotation.hpp"
#include <iostream>

extern "C" { void __bessel_mod_MOD_besiexp( double*, double*, int*,
			   	double[], double[], double[],
			   	double[], double[], double[] ); }

#define MAXL_ 10

dielectric::dielectric ( StixVars _s )
{
		I = std::complex<double>(0.0,1.0);
		epsilon.zeros(3,3);

		epsilon(0,0) = _s.S;
		epsilon(0,1) = -I*_s.D ;
		epsilon(0,2) = 0;

		epsilon(1,0) = I*_s.D;
		epsilon(1,1) = _s.S ;
		epsilon(1,2) = 0;

		epsilon(2,0) = 0;
		epsilon(2,1) = 0 ;
		epsilon(2,2) = _s.P;
}

dielectric::dielectric ( std::vector<HotPlasmaSpecies> _s, double _omega, int _l, 
				std::complex<double> _ky, std::complex<double> _kz, RotationMatrix _R )
{
		I = std::complex<double>(0.0,1.0);
		epsilon.zeros(3,3);

		omega = _omega;
		omega_c = std::complex<double>(omega,0.0);
		l = _l;
		species = _s;

		R = _R;

		k_xyz = arma::zeros<arma::cx_colvec>(3);
		k_abp = arma::zeros<arma::cx_colvec>(3);

		k_xyz(1) = _ky;
		k_xyz(2) = _kz;

}

void dielectric::populateSwansonKs(std::complex<double> _kx)
{

		k_xyz(0) = _kx;
		k_abp = R.xyz2abp * k_xyz;

		//k_xyz.print("k_xyz: ");
		//k_abp.print("k_abp: ");

		K0 = std::complex<double>(0.0,0.0);
		K1 = std::complex<double>(0.0,0.0);
		K2 = std::complex<double>(0.0,0.0);
		K3 = std::complex<double>(0.0,0.0);
		K4 = std::complex<double>(0.0,0.0);
		K5 = std::complex<double>(0.0,0.0);
	
		double expBes_r[MAXL_];
		double expBesPrime_r[MAXL_];
		double expBesOverGam_r[MAXL_];	
		double expBes_i[MAXL_];
		double expBesPrime_i[MAXL_];
		double expBesOverGam_i[MAXL_];	

		std::complex<double> kPer = sqrt(pow(k_abp(0),2)+pow(k_abp(1),2));
		std::complex<double> cosPsi = k_abp(0)/kPer;
		std::complex<double> sinPsi = k_abp(1)/kPer;

		epsilon.zeros(3,3); // initialize the dielectric to zero

		for(int s=0;s<species.size();s++) // species loop
		{
			double rho = species[s].vTh / omega;
			std::complex<double> lambda = pow(kPer,2)*pow(rho,2) / 2.0;

			std::complex<double> factor1 = pow(species[s].wp,2)/(omega_c*k_abp(2)*species[s].vTh);
#ifdef _DEBUG
			std::cout<<"factor1 components: "<<std::endl;
			std::cout<<"wp: "<<species[s].wp<<std::endl;
			std::cout<<"omega_c: "<<omega_c<<std::endl;
			std::cout<<"kp: "<<k_abp(2)<<std::endl;
			std::cout<<"vTh: "<<species[s].vTh<<std::endl;
			std::cout<<"density: "<<species[s].n<<std::endl;
#endif
			std::complex<double> factor2 = kPer*pow(species[s].wp,2)/(2.0*k_abp(2)*omega_c*species[s].wc);
			double e_swan = species[s].z/abs(species[s].z);

			double lambda_r = lambda.real();
			double lambda_i = lambda.imag();

			__bessel_mod_MOD_besiexp(&lambda_r,&lambda_i,&l,
					expBes_r,expBes_i,
					expBesPrime_r,expBesPrime_i,
					expBesOverGam_r,expBesOverGam_i);					

			arma::cx_colvec expBes(l+1);
			arma::cx_colvec expBesPrime(l+1);
			arma::cx_colvec expBesOverGam(l+1);

			for(int n=0;n<=l;n++)
			{
					expBes(n) = std::complex<double>(expBes_r[n],expBes_i[n]);
					expBesPrime(n) = std::complex<double>(expBesPrime_r[n],expBesPrime_i[n]);
					expBesOverGam(n) = std::complex<double>(expBesOverGam_r[n],expBesOverGam_i[n]);
#ifdef _DEBUG
					std::cout << "expBes("<<n<<"): " << expBes(abs(n)) <<std::endl;
					std::cout << "expBesPrime("<<n<<"): " << expBesPrime(abs(n)) <<std::endl;
					std::cout << "expBesOverGam("<<n<<"): " << expBesOverGam(abs(n));
					std::cout<<", lambda: " << lambda <<", rho: " << rho << ", vTh: " << species[s].vTh << std::endl;
#endif
			}

			std::complex<double> K0_s,K1_s,K2_s,K3_s,K4_s,K5_s;

			K0_s = std::complex<double>(0.0,0.0);
			K1_s = std::complex<double>(0.0,0.0);
			K2_s = std::complex<double>(0.0,0.0);
			K3_s = std::complex<double>(0.0,0.0);
			K4_s = std::complex<double>(0.0,0.0);
			K5_s = std::complex<double>(0.0,0.0);

			for(int n=-l;n<=l;n++) // harmonic number loop
			{

					int nabs = abs(n);
					std::complex<double> n_c(n,0.0);

					std::complex<double> zeta_n = ( omega_c + n*species[s].wc ) /  (k_abp(2)*species[s].vTh);
					std::complex<double> w = MATPACK::Faddeeva_2(zeta_n);
					//std::complex<double> w = MATPACK::Faddeeva(zeta_n);
					std::complex<double> Z = sqrt(_pi)*I*w; // Z Function
					std::complex<double> Zprime = -2.0*(1.0+zeta_n*Z);
#ifdef _DEBUG
					std::cout << "Z: " << Z << " , zeta_n: " << zeta_n << std::endl;
#endif
					std::complex<double> InExp = expBes(nabs);
					std::complex<double> InPrimeExp = expBesPrime(nabs);
					std::complex<double> InOverLambdaExp = expBesOverGam(nabs);

					K0_s += lambda*(InExp-InPrimeExp)*Z;
					K1_s += pow(n_c,2)*InOverLambdaExp*Z;
					K2_s += n_c*(InExp-InPrimeExp)*Z;
					K3_s += InExp*zeta_n*Zprime;
					K4_s += n_c*InOverLambdaExp*Zprime;
					K5_s += (InExp-InPrimeExp)*Zprime;
#ifdef _DEBUG
					std::cout<<std::endl<<std::endl;
					std::cout<<"n: "<<n<<std::endl;
					std::cout<<"s: "<<s<<std::endl;
					std::cout<<"lambda: "<<lambda<<std::endl;
					std::cout<<"vTh: "<<species[s].vTh<<std::endl;
					std::cout<<"k_abp(2): "<<k_abp(2)<<std::endl;
					std::cout << "K0: " << K0_s << std::endl;
					std::cout << "K1: " << K1_s << std::endl;
					std::cout << "K2: " << K2_s << std::endl;
					std::cout << "K3: " << K3_s << std::endl;
					std::cout << "K4: " << K4_s << std::endl;
					std::cout << "K5: " << K5_s << std::endl;
					std::cout << "Zprime: "<<Zprime<<std::endl;
					std::cout << "1/omega: "<<1.0/omega_c<<std::endl;
					std::cout << "1/omega cold: "<<InExp*zeta_n*Zprime/(k_abp(2)*species[s].vTh)<<std::endl;
#endif
			}

			K0_s = factor1*K0_s;
			K1_s = factor1*K1_s;
			K2_s = e_swan*factor1*K2_s;
			K3_s = factor1*K3_s;
			K4_s = factor2*K4_s;
			K5_s = e_swan*factor2*K5_s;
#ifdef _DEBUG
			std::cout << "factor1: " << factor1 << std::endl;
			std::cout << "factor2: " << factor2 << std::endl;
			std::cout << "e_swan: " << e_swan << std::endl;
			std::cout <<"sinPsi: "<<sinPsi<<std::endl;
			std::cout <<"cosPsi: "<<cosPsi<<std::endl;
#endif
			// rememeber the above is the dielectric 
			// and below we turn it into a conductivity
			
			K0 += K0_s;
			K1 += K1_s;
			K2 += K2_s;
			K3 += K3_s;
			K4 += K4_s;
			K5 += K5_s;

		}

		K0 = 2.0*K0;
		K1 = 1.0+K1;
		K2 = I*K2;
		K3 = 1.0-K3;
		K4 = K4;
		K5 = I*K5;

		Ka0 = K0*pow(omega_c,2)/pow(_c,2);
		Ka1 = K1*pow(omega_c,2)/pow(_c,2);
		Ka2 = K2*pow(omega_c,2)/pow(_c,2);
		Ka3 = K3*pow(omega_c,2)/pow(_c,2);
		Ka4 = K4*pow(omega_c,2)/pow(_c,2);
		Ka5 = K5*pow(omega_c,2)/pow(_c,2);

		epsilon(0,0) += K1+pow(sinPsi,2)*K0;
		epsilon(0,1) += K2-cosPsi*sinPsi*K0;
		epsilon(0,2) += cosPsi*K4+sinPsi*K5;
	
		epsilon(1,0) += -K2-cosPsi*sinPsi*K0;
		epsilon(1,1) += K1+pow(cosPsi,2)*K0;
		epsilon(1,2) += sinPsi*K4-cosPsi*K5;

		epsilon(2,0) += cosPsi*K4-sinPsi*K5;
		epsilon(2,1) += sinPsi*K4+cosPsi*K5;
		epsilon(2,2) += K3;
	
		//arma::mat ident;
		//ident.zeros(3,3);
		//ident(0,0) = 1;
		//ident(1,1) = 1;
		//ident(2,2) = 1;

		//stix = -(stix-ident)*omega_c*_e0*I;
		//stix.print("Hot conductivity stix:");
}

//void dielectric::determinant (double kxIn_re, double kxIn_im, double *kxOu_re, double *kxOu_im)
std::complex<double> dielectric::determinant (std::complex<double> _kx)
{
	//std::complex<double> kx = std::complex<double>(kxIn_re,kxIn_im);
	std::complex<double> kx = _kx;//std::complex<double>(kxIn_re,kxIn_im);

	//k_xyz(0) = kx;
	//k_abp = R.xyz2abp * k_xyz;

	populateSwansonKs(kx);

	std::complex<double> kPer = sqrt(pow(k_abp(0),2)+pow(k_abp(1),2));
	std::complex<double> gamma = pow(k_abp(2),2)-Ka1;

	// pg 177 of Swanson
	
	//std::complex<double> determinant =
	//		(gamma*(gamma-Ka0+pow(kPer,2))+pow(Ka2,2))*Ka3
	//		+pow(kPer,2)*((gamma-Ka0+pow(kPer,2))*Ka1-pow(Ka2,2))
	//		+Ka4*(gamma-Ka0+pow(kPer,2))*(2.0*kPer*k_abp(2)+Ka4)
	//		-Ka5*(gamma*Ka5+2.0*Ka2*(kPer*k_abp(2)+Ka4));

	arma::cx_mat WaveEqn = arma::zeros<arma::cx_mat>(3,3);
	std::complex<double> w2_c2 = pow(omega_c,2)/pow(_c,2);

	WaveEqn(0,0) = epsilon(0,0)*w2_c2-pow(k_abp(2),2)-pow(k_abp(1),2);
	WaveEqn(0,1) = epsilon(0,1)*w2_c2+k_abp(0)*k_abp(1);
	WaveEqn(0,2) = epsilon(0,2)*w2_c2+k_abp(0)*k_abp(2);

	WaveEqn(1,0) = epsilon(1,0)*w2_c2+k_abp(1)*k_abp(0);
	WaveEqn(1,1) = epsilon(1,1)*w2_c2-pow(k_abp(2),2)-pow(k_abp(0),2);
	WaveEqn(1,2) = epsilon(1,2)*w2_c2+k_abp(1)*k_abp(2);

	WaveEqn(2,0) = epsilon(2,0)*w2_c2+k_abp(2)*k_abp(0);
	WaveEqn(2,1) = epsilon(2,1)*w2_c2+k_abp(2)*k_abp(1);
	WaveEqn(2,2) = epsilon(2,2)*w2_c2-pow(kPer,2);

	return(arma::det(WaveEqn));
	//std::complex<double> determinantA2 = arma::det(WaveEqn2);

	//double determinantA_ = arma::det(arma::real(WaveEqn));

	//*kxOu_re = std::real(determinantA)/_c;	
	//*kxOu_im = std::imag(determinantA)/_c;
}

void dielectric::rotate ( arma::mat _rot )
{
		epsilon = _rot * epsilon * arma::inv(_rot);
}

void dielectric::rotate ( arma::cx_mat _rot )
{
		epsilon = _rot * epsilon * arma::inv(_rot);
}



