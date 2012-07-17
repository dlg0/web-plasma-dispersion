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
				std::vector<complex<double> > _k_cyl )
{
		complex<double> K0,K1,K2,K3,K4,K5;

		for(int s=0;s<_s.size();s++) // species loop
		{
			for(int n=-_l;n<=_l;n++) // harmonic number loop
			{
					complex<double> zeta_n = ( _omega - n*_s[s].wc ) /  (_k_cyl[1]*_s[s].vTh);
			}

		}
}

void dielectric::rotateEpsilon ( std::vector<float> bUnit_car )
{
    // get vector perp to both z axis and b

		arma::vec zaxis_ = arma::zeros<arma::vec>(3);
		zaxis_(2) = 1;
		arma::vec bUnit_ = arma::zeros<arma::vec>(3);
		bUnit_(0) = bUnit_car[0];
		bUnit_(1) = bUnit_car[1];
		bUnit_(2) = bUnit_car[2];

		arma::vec perp_ = arma::cross(zaxis_,bUnit_);

    // check angle between perp and b

		double perpDotB_ = arma::dot(perp_,bUnit_);
		double perpMag_ = arma::norm(perp_,2);
		double bUnitMag_ = arma::norm(bUnit_,2);

		arma::vec perpUnit_(perp_ / perpMag_);

    // get angle between z axis and b

    	float theta = acos ( bUnit_(2) );

    // calculate the quaternions

    	float q0  = cos ( theta / 2.0 );
    	float q1  = sin ( theta / 2.0 ) * perpUnit_(0);
    	float q2  = sin ( theta / 2.0 ) * perpUnit_(1);
    	float q3  = sin ( theta / 2.0 ) * perpUnit_(2);

    // construct the rotation matrix

		arma::mat rotQ = arma::zeros<arma::mat>(3,3);

    	rotQ(0,0) = pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2); 
		rotQ(0,1) = 2*(q1*q2-q0*q3);
		rotQ(0,2) = 2*(q1*q3+q0*q2);

		rotQ(1,0) = 2*(q2*q1+q0*q3);
		rotQ(1,1) = pow(q0,2)-pow(q1,2)+pow(q2,2)-pow(q3,2);
		rotQ(1,2) = 2*(q2*q3-q0*q1);

		rotQ(2,0) = 2*(q3*q1-q0*q2);
		rotQ(2,1) = 2*(q3*q2+q0*q1);
		rotQ(2,2) = pow(q0,2)-pow(q1,2)-pow(q2,2)+pow(q3,2);

	// stix rotated

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

