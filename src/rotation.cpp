#include "rotation.hpp"
#include <armadillo>
#include <vector>
#include <math.h>
#include "constants.hpp"

RotationMatrix::RotationMatrix()
{
}

RotationMatrix::RotationMatrix(arma::colvec _b_xyz)
{

	
	// vectors are:
	//
	// 	alpha,beta,parallel -> abp
	// 	x,y,z
	//
	//	y cross p = a
	//
	//	u are unit vectors

    // get vector perp to both y axis and b

		arma::vec xu_xyz = arma::zeros<arma::vec>(3);
		xu_xyz(0) = 1;
		arma::vec yu_xyz = arma::zeros<arma::vec>(3);
		yu_xyz(1) = 1;
		arma::vec zu_xyz = arma::zeros<arma::vec>(3);
		zu_xyz(2) = 1;


		arma::vec pu_xyz = arma::zeros<arma::vec>(3);
		pu_xyz(0) = _b_xyz(0)/arma::norm(_b_xyz,2);
		pu_xyz(1) = _b_xyz(1)/arma::norm(_b_xyz,2);
		pu_xyz(2) = _b_xyz(2)/arma::norm(_b_xyz,2);

		arma::vec au_xyz = arma::cross(yu_xyz,pu_xyz);
		arma::vec bu_xyz = arma::cross(pu_xyz,au_xyz);

		au_xyz.print("au_xyz: ");
		bu_xyz.print("bu_xyz: ");
		pu_xyz.print("pu_xyz: ");

    // get angle between z axis and p

    	float theta = acos ( arma::dot(zu_xyz,pu_xyz) );

		std::cout<<"theta: "<<theta*180.0/_pi<<std::endl;

	// rotation is equivilant to that rotation that will
	// make the xyz axes lie on to of abp axes. Here we
	// rotate around b by angle theta.

    // calculate the quaternions

    	float q0  = cos ( theta / 2.0 );
    	float q1  = sin ( theta / 2.0 ) * bu_xyz(0);
    	float q2  = sin ( theta / 2.0 ) * bu_xyz(1);
    	float q3  = sin ( theta / 2.0 ) * bu_xyz(2);

    // construct the rotation matrix

		xyz2abp = arma::zeros<arma::mat>(3,3);

    	xyz2abp(0,0) = pow(q0,2)+pow(q1,2)-pow(q2,2)-pow(q3,2); 
		xyz2abp(0,1) = 2*(q1*q2-q0*q3);
		xyz2abp(0,2) = 2*(q1*q3+q0*q2);

		xyz2abp(1,0) = 2*(q2*q1+q0*q3);
		xyz2abp(1,1) = pow(q0,2)-pow(q1,2)+pow(q2,2)-pow(q3,2);
		xyz2abp(1,2) = 2*(q2*q3-q0*q1);

		xyz2abp(2,0) = 2*(q3*q1-q0*q2);
		xyz2abp(2,1) = 2*(q3*q2+q0*q1);
		xyz2abp(2,2) = pow(q0,2)-pow(q1,2)-pow(q2,2)+pow(q3,2);

		abp2xyz = arma::zeros<arma::mat>(3,3);
		abp2xyz = arma::inv(xyz2abp);

		xyz2abp.print("xyz2abp: ");

	//// Now initialize the rotation matrix from stx to abp
	//// --------------------------------------------------
	//
	//	arma::cx_colvec k_abp = xyz2abp * _k_xyz; // get k in alp,bet,prl (_abp)

	//	stx2abp = arma::zeros<arma::cx_mat>(3,3);

	//	std::complex<double> kPerp = sqrt(pow(k_abp(0),2)+pow(k_abp(1),2));
	//	std::complex<double> sinPsi = k_abp(0)/kPerp; 
	//	std::complex<double> cosPsi = k_abp(1)/kPerp; 

	//	stx2abp(0,0) = cosPsi;
	//	stx2abp(0,1) = -sinPsi;
	//	stx2abp(0,2) = 0;

	//	stx2abp(1,0) = sinPsi;
	//	stx2abp(1,1) = cosPsi;
	//	stx2abp(1,2) = 0;

	//	stx2abp(2,0) = 0;
	//	stx2abp(2,1) = 0;
	//	stx2abp(2,2) = 1;

	//	abp2stx = arma::zeros<arma::cx_mat>(3,3);
	//	abp2stx = arma::inv(stx2abp);

	//	stx2abp.print("Stx2ABP: ");
}

