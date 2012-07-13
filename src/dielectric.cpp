#include "dielectric.hpp"
#include "constants.hpp"
#include <math.h>
#include "arrays.hpp"
#include "stixVars.hpp"
#include <complex>
#include <vector>
#include <numeric>
#include <armadillo>

dielectric::dielectric ( StixVars s )
{
		I = (0,1);
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
