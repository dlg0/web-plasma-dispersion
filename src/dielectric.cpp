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

dielectric::dielectric ( StixVars s )
{
		I = std::complex<double>(0.0,1.0);
		stix.zeros(3,3);

		stix(0,0) = s.S;
		stix(0,1) = I*s.D ;
		stix(0,2) = 0;

		stix(1,0) = -I*s.D;
		stix(1,1) = s.S ;
		stix(1,2) = 0;

		stix(2,0) = 0;
		stix(2,1) = 0 ;
		stix(2,2) = s.P;

		stix.print();
}

dielectric::dielectric ( std::vector<HotPlasmaSpecies> _s, double _omega, int _l, 
				arma::cx_colvec _k_cyl, std::vector<float> bUnit_cyl )
{
		std::complex<double> K0,K1,K2,K3,K4,K5;
		std::complex<double> I(0.0,1.0);

		RotationMatrix _rotQ(bUnit_cyl);
		rotQ = _rotQ.rotQ;
		rotQ.print("RotQ:");

		arma::cx_colvec k_abp = rotQ * _k_cyl; // get k in alp,bet,prl (_abp)

		for(int s=0;s<_s.size();s++) // species loop
		{
			for(int n=-_l;n<=_l;n++) // harmonic number loop
			{
					std::complex<double> zeta_n = ( _omega - n*_s[s].wc ) /  (k_abp(2)*_s[s].vTh);
					std::complex<double> w = MATPACK::Faddeeva_2(zeta_n);
					std::complex<double> Z = sqrt(_pi)*I*w; // Z Function
					std::complex<double> Zprime = -2.0*(1.0+zeta_n*Z);
					std::cout << "Z: " << Z << std::endl;
			}

		}
}

void dielectric::rotateEpsilon ( std::vector<float> bUnit_car )
{
	// stix rotated

		RotationMatrix _rotQ(bUnit_car);
		rotQ = _rotQ.rotQ;
		stixRotated = rotQ * stix * arma::inv(rotQ);
		stixRotated.print();

}

void dielectric::coldRoots( double w, double kp, double kz)
{

		// see rsfxc_1D.nb for the derivation of these polynomial coeffs.

		std::complex<double> k4 = pow(w,2)/pow(_c,2) * stixRotated(0,0);

		std::complex<double> k3 = pow(w,2)/pow(_c,2) 
				* ( kp * ( stixRotated(0,1) + stixRotated(1,0) ) 
						+ kz * ( stixRotated(0,2) + stixRotated(2,0) ) );

		std::complex<double> k2 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(kp,2) * ( stixRotated(0,0) + stixRotated(1,1) ) 
						+ pow(_c,2) * kz * kp * ( stixRotated(1,2) + stixRotated(2,1) ) 
						+ pow(_c,2) * pow(kz,2) * ( stixRotated(0,0) + stixRotated(2,2) ) 
						+ pow(w,2) * ( stixRotated(0,1) * stixRotated(1,0) 
										+ stixRotated(0,2) * stixRotated(2,0)
										- stixRotated(0,0) * ( stixRotated(1,1) + stixRotated(2,2) ) ) );

		std::complex<double> k1 = pow(w,2)/pow(_c,4) * ( pow(_c,2) * pow(kz,2) * kp * ( stixRotated(0,1) + stixRotated(1,0) )
						+ pow(_c,2) * pow(kz,3) * ( stixRotated(0,2) + stixRotated(2,0) )
						+ kz * ( pow(_c,2) * pow(kp,2) * ( stixRotated(0,2) + stixRotated(2,0) )
							+ pow(w,2) * ( stixRotated(0,1) * stixRotated(1,2) 
										- stixRotated(1,1) * ( stixRotated(0,2) + stixRotated(2,0) )
										+ stixRotated(1,0) * stixRotated(2,1) ) ) 
						+ kp * ( pow(w,2) * ( stixRotated(1,2) * stixRotated(2,0)
												+ stixRotated(0,2) * stixRotated(2,1) )
									+ ( stixRotated(0,1) + stixRotated(1,0) ) 
										* ( pow(_c,2) * pow(kp,2) - pow(w,2) * stixRotated(2,2) ) ) );
		
		std::complex<double> k0 = pow(w,2)/pow(_c,6) * ( pow(_c,4) * pow(kp,4) * stixRotated(1,1) 
							+ pow(_c,4) * kz * pow(kp,3) * ( stixRotated(1,2) + stixRotated(2,1) ) 
							+ pow(_c,2) * kz * kp 
                                * ( pow(_c,2) * pow(kz,2) * ( stixRotated(1,2) + stixRotated(2,1) ) 
									+ pow(w,2) * ( stixRotated(0,2) * stixRotated(1,0) 
													+ stixRotated(0,1) * stixRotated(2,0) 
													- stixRotated(0,0) * 
														( stixRotated(1,2) + stixRotated(2,1) ) ) ) 
							+ pow(_c,4) * pow(kz,4) * stixRotated(2,2) 
							+ pow(_c,2) * pow(w,2) * pow(kz,2) 
								* ( stixRotated(0,2) * stixRotated(2,0) 
									+ stixRotated(1,2) * stixRotated(2,1) 
									- ( stixRotated(0,0) + stixRotated(1,1) ) * stixRotated(2,2) ) 
							+pow(w,4) * ( stixRotated(0,2) * ( -stixRotated(1,1) * stixRotated(2,0) 
															+ stixRotated(1,0) * stixRotated(2,1) ) 
								+ stixRotated(0,1) * ( stixRotated(1,2) * stixRotated(2,0) 
													- stixRotated(1,0) * stixRotated(2,2) ) 
								+ stixRotated(0,0) * ( -stixRotated(1,2) * stixRotated(2,1) 
													+ stixRotated(1,1) * stixRotated(2,2) ) ) 
							+ pow(_c,2) * pow(kp,2) * ( pow(_c,2) * pow(kz,2) * 
								( stixRotated(1,1) + stixRotated(2,2) ) 
								+ pow(w,2) * ( stixRotated(0,1) * stixRotated(1,0) 
											+ stixRotated(1,2) * stixRotated(2,1) 
											- stixRotated(1,1) * 
												( stixRotated(0,0) + stixRotated(2,2) ) ) ) );

		std::cout << "Coeffs: " << std::endl;
		std::cout << "k0: " << k0 << std::endl;
		std::cout << "k1: " << k1 << std::endl;
		std::cout << "k2: " << k2 << std::endl;
		std::cout << "k3: " << k3 << std::endl;
		std::cout << "k4: " << k4 << std::endl;

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

