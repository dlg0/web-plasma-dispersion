#include "rotation.hpp"
#include <armadillo>
#include <vector>
#include <math.h>

RotationMatrix::RotationMatrix(arma::colvec _b_rtz, arma::cx_colvec _k_rtz )
{

	// Rotation from rtz to abp
	// ------------------------
	
    // get vector perp to both z axis and b

		arma::vec zaxis_ = arma::zeros<arma::vec>(3);
		zaxis_(2) = 1;
		arma::vec bUnit_rtz = arma::zeros<arma::vec>(3);
		bUnit_rtz(0) = _b_rtz(0)/arma::norm(_b_rtz,2);
		bUnit_rtz(1) = _b_rtz(1)/arma::norm(_b_rtz,2);
		bUnit_rtz(2) = _b_rtz(2)/arma::norm(_b_rtz,2);

		arma::vec perp_ = arma::cross(zaxis_,bUnit_rtz);

    // check angle between perp and b

		double perpDotB_ = arma::dot(perp_,bUnit_rtz);
		double perpMag_ = arma::norm(perp_,2);
		double bUnitMag_ = arma::norm(bUnit_rtz,2);

		arma::vec perpUnit_(perp_ / perpMag_);

    // get angle between z axis and b

    	float theta = acos ( bUnit_rtz(2) );

    // calculate the quaternions

    	float q0  = cos ( theta / 2.0 );
    	float q1  = sin ( theta / 2.0 ) * perpUnit_(0);
    	float q2  = sin ( theta / 2.0 ) * perpUnit_(1);
    	float q3  = sin ( theta / 2.0 ) * perpUnit_(2);

    // construct the rotation matrix

		rtz2abp = arma::zeros<arma::mat>(3,3);

    	rtz2abp(0,0) = pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2); 
		rtz2abp(0,1) = 2*(q1*q2-q0*q3);
		rtz2abp(0,2) = 2*(q1*q3+q0*q2);

		rtz2abp(1,0) = 2*(q2*q1+q0*q3);
		rtz2abp(1,1) = pow(q0,2)-pow(q1,2)+pow(q2,2)-pow(q3,2);
		rtz2abp(1,2) = 2*(q2*q3-q0*q1);

		rtz2abp(2,0) = 2*(q3*q1-q0*q2);
		rtz2abp(2,1) = 2*(q3*q2+q0*q1);
		rtz2abp(2,2) = pow(q0,2)-pow(q1,2)-pow(q2,2)+pow(q3,2);

		abp2rtz = arma::zeros<arma::mat>(3,3);
		abp2rtz = arma::inv(rtz2abp);


	// Now initialize the rotation matrix from stx to abp
	// --------------------------------------------------
	
		arma::cx_colvec k_abp = rtz2abp * _k_rtz; // get k in alp,bet,prl (_abp)

		stx2abp = arma::zeros<arma::cx_mat>(3,3);

		std::complex<double> kPerp = sqrt(pow(k_abp(0),2)+pow(k_abp(1),2));
		std::complex<double> sinPsi = k_abp(0)/kPerp; 
		std::complex<double> cosPsi = k_abp(1)/kPerp; 

		stx2abp(0,0) = cosPsi;
		stx2abp(0,1) = -sinPsi;
		stx2abp(0,2) = 0;

		stx2abp(1,0) = sinPsi;
		stx2abp(1,1) = cosPsi;
		stx2abp(1,2) = 0;

		stx2abp(2,0) = 0;
		stx2abp(2,1) = 0;
		stx2abp(2,2) = 1;

		abp2stx = arma::zeros<arma::cx_mat>(3,3);
		abp2stx = arma::inv(stx2abp);

		stx2abp.print("Stx2ABP: ");
}

