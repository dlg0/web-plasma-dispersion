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

dielectric::dielectric ( StixVars s )
{
		I = std::complex<double>(0.0,1.0);
		epsilon.zeros(3,3);

		epsilon(0,0) = s.S;
		epsilon(0,1) = -I*s.D ;
		epsilon(0,2) = 0;

		epsilon(1,0) = I*s.D;
		epsilon(1,1) = s.S ;
		epsilon(1,2) = 0;

		epsilon(2,0) = 0;
		epsilon(2,1) = 0 ;
		epsilon(2,2) = s.P;
}

dielectric::dielectric ( std::vector<HotPlasmaSpecies> _s, double _omega, int _l, arma::cx_colvec _k_abp )
{
		std::complex<double> K0,K1,K2,K3,K4,K5;
		std::complex<double> I(0.0,1.0);

		double expBes_r[MAXL_];
		double expBesPrime_r[MAXL_];
		double expBesOverGam_r[MAXL_];	
		double expBes_i[MAXL_];
		double expBesPrime_i[MAXL_];
		double expBesOverGam_i[MAXL_];	

		std::complex<double> omega_c(_omega,_omega*0.0);
		std::complex<double> kPer = sqrt(pow(_k_abp(0),2)+pow(_k_abp(1),2));
		std::complex<double> cosPsi = _k_abp(0)/kPer;
		std::complex<double> sinPsi = _k_abp(1)/kPer;

		epsilon.zeros(3,3); // initialize the dielectric to zero

		for(int s=0;s<_s.size();s++) // species loop
		{
			double rho = _s[s].vTh / _omega;
			std::complex<double> lambda = pow(kPer,2)*pow(rho,2) / 2.0;

			std::complex<double> factor1 = pow(_s[s].wp,2)/(omega_c*_k_abp(2)*_s[s].vTh);
#ifdef _DEBUG
			std::cout<<"factor1 components: "<<std::endl;
			std::cout<<"wp: "<<_s[s].wp<<std::endl;
			std::cout<<"omega_c: "<<omega_c<<std::endl;
			std::cout<<"kp: "<<_k_abp(2)<<std::endl;
			std::cout<<"vTh: "<<_s[s].vTh<<std::endl;
			std::cout<<"density: "<<_s[s].n<<std::endl;
#endif
			std::complex<double> factor2 = kPer*pow(_s[s].wp,2)/(2.0*_k_abp(2)*omega_c*_s[s].wc);
			double e_swan = _s[s].z/abs(_s[s].z);

			double lambda_r = lambda.real();
			double lambda_i = lambda.imag();

			__bessel_mod_MOD_besiexp(&lambda_r,&lambda_i,&_l,
					expBes_r,expBes_i,
					expBesPrime_r,expBesPrime_i,
					expBesOverGam_r,expBesOverGam_i);					

			arma::cx_colvec expBes(_l+1);
			arma::cx_colvec expBesPrime(_l+1);
			arma::cx_colvec expBesOverGam(_l+1);

			for(int n=0;n<=_l;n++)
			{
					expBes(n) = std::complex<double>(expBes_r[n],expBes_i[n]);
					expBesPrime(n) = std::complex<double>(expBesPrime_r[n],expBesPrime_i[n]);
					expBesOverGam(n) = std::complex<double>(expBesOverGam_r[n],expBesOverGam_i[n]);
#ifdef _DEBUG
					std::cout << "expBes("<<n<<"): " << expBes(abs(n)) <<std::endl;
					std::cout << "expBesPrime("<<n<<"): " << expBesPrime(abs(n)) <<std::endl;
					std::cout << "expBesOverGam("<<n<<"): " << expBesOverGam(abs(n));
					std::cout<<", lambda: " << lambda <<", rho: " << rho << ", vTh: " << _s[s].vTh << std::endl;
#endif
			}

			K0 = std::complex<double>(0.0,0.0);
			K1 = std::complex<double>(0.0,0.0);
			K2 = std::complex<double>(0.0,0.0);
			K3 = std::complex<double>(0.0,0.0);
			K4 = std::complex<double>(0.0,0.0);

			for(int n=-_l;n<=_l;n++) // harmonic number loop
			{

					int nabs = abs(n);
					std::complex<double> n_c(n,0.0);

					std::complex<double> zeta_n = ( omega_c + n*_s[s].wc ) /  (_k_abp(2)*_s[s].vTh);
					std::complex<double> w = MATPACK::Faddeeva_2(zeta_n);
					std::complex<double> Z = sqrt(_pi)*I*w; // Z Function
					std::complex<double> Zprime = -2.0*(1.0+zeta_n*Z);
#ifdef _DEBUG
					std::cout << "Z: " << Z << " , zeta_n: " << zeta_n << std::endl;
#endif
					std::complex<double> InExp = expBes(nabs);
					std::complex<double> InPrimeExp = expBesPrime(nabs);
					std::complex<double> InOverLambdaExp = expBesOverGam(nabs);

					K0 += lambda*(InExp-InPrimeExp)*Z;
					K1 += pow(n_c,2)*InOverLambdaExp*Z;
					K2 += n_c*(InExp-InPrimeExp)*Z;
					K3 += InExp*zeta_n*Zprime;
					K4 += n_c*InOverLambdaExp*Zprime;
					K5 += (InExp-InPrimeExp)*Zprime;
#ifdef _DEBUG
					std::cout << "K0: " << K0 << std::endl;
					std::cout << "K1: " << K1 << std::endl;
					std::cout << "K2: " << K2 << std::endl;
					std::cout << "K3: " << K3 << std::endl;
					std::cout << "K4: " << K4 << std::endl;
					std::cout << "K5: " << K5 << std::endl;
					std::cout << "Zprime: "<<Zprime<<std::endl;
#endif
			}

			K0 = 2.0*factor1*K0;
			K1 = 1.0+factor1*K1;
			K2 = I*e_swan*factor1*K2;
			K3 = 1.0-factor1*K3;
			K4 = factor2*K4;
			K5 = I*e_swan*factor2*K5;
#ifdef _DEBUG
			std::cout << "factor1: " << factor1 << std::endl;
			std::cout << "factor2: " << factor2 << std::endl;
			std::cout << "e_swan: " << e_swan << std::endl;
			std::cout <<"sinPsi: "<<sinPsi<<std::endl;
			std::cout <<"cosPsi: "<<cosPsi<<std::endl;
#endif
			// rememeber the above is the dielectric 
			// and below we turn it into a conductivity
			
			epsilon(0,0) += K1+pow(sinPsi,2)*K0;
			epsilon(0,1) += K2-cosPsi*sinPsi*K0;
			epsilon(0,2) += cosPsi*K4+sinPsi*K5;
	
			epsilon(1,0) += -K2-cosPsi*sinPsi*K0;
			epsilon(1,1) += K1+pow(cosPsi,2)*K0;
			epsilon(1,2) += sinPsi*K4-cosPsi*K5;

			epsilon(2,0) += cosPsi*K4-sinPsi*K5;
			epsilon(2,1) += sinPsi*K4+cosPsi*K5;
			epsilon(2,2) += K3;
	
		}

		arma::mat ident;
		ident.zeros(3,3);
		ident(0,0) = 1;
		ident(1,1) = 1;
		ident(2,2) = 1;

		//stix = -(stix-ident)*omega_c*_e0*I;
		//stix.print("Hot conductivity stix:");

}

void dielectric::rotate ( arma::mat _rot )
{
		epsilon = _rot * epsilon * arma::inv(_rot);
}

void dielectric::rotate ( arma::cx_mat _rot )
{
		epsilon = _rot * epsilon * arma::inv(_rot);
}

void dielectric::coldRoots( double w, double kp, double kz)
{

		// see rsfxc_1D.nb for the derivation of these polynomial coeffs.

		std::complex<double> k4 = pow(w,2)/pow(_c,2) * epsilon(0,0);

		std::complex<double> k3 = pow(w,2)/pow(_c,2) 
				* ( kp * ( epsilon(0,1) + epsilon(1,0) ) 
						+ kz * ( epsilon(0,2) + epsilon(2,0) ) );

		std::complex<double> k2 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(kp,2) * ( epsilon(0,0) + epsilon(1,1) ) 
						+ pow(_c,2) * kz * kp * ( epsilon(1,2) + epsilon(2,1) ) 
						+ pow(_c,2) * pow(kz,2) * ( epsilon(0,0) + epsilon(2,2) ) 
						+ pow(w,2) * ( epsilon(0,1) * epsilon(1,0) 
										+ epsilon(0,2) * epsilon(2,0)
										- epsilon(0,0) * ( epsilon(1,1) + epsilon(2,2) ) ) );

		std::complex<double> k1 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(kz,2) * kp * ( epsilon(0,1) + epsilon(1,0) )
						+ pow(_c,2) * pow(kz,3) * ( epsilon(0,2) + epsilon(2,0) )
						+ kz * ( pow(_c,2) * pow(kp,2) * ( epsilon(0,2) + epsilon(2,0) )
							+ pow(w,2) * ( epsilon(0,1) * epsilon(1,2) 
										- epsilon(1,1) * ( epsilon(0,2) + epsilon(2,0) )
										+ epsilon(1,0) * epsilon(2,1) ) ) 
						+ kp * ( pow(w,2) * ( epsilon(1,2) * epsilon(2,0)
												+ epsilon(0,2) * epsilon(2,1) )
									+ ( epsilon(0,1) + epsilon(1,0) ) 
										* ( pow(_c,2) * pow(kp,2) - pow(w,2) * epsilon(2,2) ) ) );
		
		std::complex<double> k0 = pow(w,2)/pow(_c,6) * ( pow(_c,4) * pow(kp,4) * epsilon(1,1) 
							+ pow(_c,4) * kz * pow(kp,3) * ( epsilon(1,2) + epsilon(2,1) ) 
							+ pow(_c,2) * kz * kp 
                                * ( pow(_c,2) * pow(kz,2) * ( epsilon(1,2) + epsilon(2,1) ) 
									+ pow(w,2) * ( epsilon(0,2) * epsilon(1,0) 
													+ epsilon(0,1) * epsilon(2,0) 
													- epsilon(0,0) * 
														( epsilon(1,2) + epsilon(2,1) ) ) ) 
							+ pow(_c,4) * pow(kz,4) * epsilon(2,2) 
							+ pow(_c,2) * pow(w,2) * pow(kz,2) 
								* ( epsilon(0,2) * epsilon(2,0) 
									+ epsilon(1,2) * epsilon(2,1) 
									- ( epsilon(0,0) + epsilon(1,1) ) * epsilon(2,2) ) 
							+pow(w,4) * ( epsilon(0,2) * ( -epsilon(1,1) * epsilon(2,0) 
															+ epsilon(1,0) * epsilon(2,1) ) 
								+ epsilon(0,1) * ( epsilon(1,2) * epsilon(2,0) 
													- epsilon(1,0) * epsilon(2,2) ) 
								+ epsilon(0,0) * ( -epsilon(1,2) * epsilon(2,1) 
													+ epsilon(1,1) * epsilon(2,2) ) ) 
							+ pow(_c,2) * pow(kp,2) * ( pow(_c,2) * pow(kz,2) * 
								( epsilon(1,1) + epsilon(2,2) ) 
								+ pow(w,2) * ( epsilon(0,1) * epsilon(1,0) 
											+ epsilon(1,2) * epsilon(2,1) 
											- epsilon(1,1) * 
												( epsilon(0,0) + epsilon(2,2) ) ) ) );
#ifdef _DEBUG
		std::cout << "Coeffs: " << std::endl;
		std::cout << "k0: " << k0 << std::endl;
		std::cout << "k1: " << k1 << std::endl;
		std::cout << "k2: " << k2 << std::endl;
		std::cout << "k3: " << k3 << std::endl;
		std::cout << "k4: " << k4 << std::endl;
#endif
		o2scl::simple_quartic_complex quart;
		roots.clear();
		roots.resize(4);
		std::complex<double> rootsTmp[4];
		quart.solve_c(k4,k3,k2,k1,k0,rootsTmp[0],rootsTmp[1],rootsTmp[2],rootsTmp[3]);

		for(int i=0;i<4;i++)
		{
				roots[i] = rootsTmp[i];
				std::cout << "Root: " << roots[i] << std::endl;
		}

		//std::sort(rootsVec.begin(),rootsVec.end());
}

