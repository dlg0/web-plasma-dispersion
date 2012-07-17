#include "rotation.hpp"
#include <armadillo>
#include <vector>
#include <math.h>

RotationMatrix::RotationMatrix(std::vector<float> bUnit_cyl )
{

    // get vector perp to both z axis and b

		arma::vec zaxis_ = arma::zeros<arma::vec>(3);
		zaxis_(2) = 1;
		arma::vec bUnit_ = arma::zeros<arma::vec>(3);
		bUnit_(0) = bUnit_cyl[0];
		bUnit_(1) = bUnit_cyl[1];
		bUnit_(2) = bUnit_cyl[2];

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

		rotQ = arma::zeros<arma::mat>(3,3);

    	rotQ(0,0) = pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2); 
		rotQ(0,1) = 2*(q1*q2-q0*q3);
		rotQ(0,2) = 2*(q1*q3+q0*q2);

		rotQ(1,0) = 2*(q2*q1+q0*q3);
		rotQ(1,1) = pow(q0,2)-pow(q1,2)+pow(q2,2)-pow(q3,2);
		rotQ(1,2) = 2*(q2*q3-q0*q1);

		rotQ(2,0) = 2*(q3*q1-q0*q2);
		rotQ(2,1) = 2*(q3*q2+q0*q1);
		rotQ(2,2) = pow(q0,2)-pow(q1,2)-pow(q2,2)+pow(q3,2);
}

