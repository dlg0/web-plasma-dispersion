#ifndef ROTATION_HPP_
#define ROTATION_HPP_

#include <armadillo>

class RotationMatrix { 
		public:
				RotationMatrix();
				RotationMatrix(arma::colvec _b_xyz);
		private:
		public:
			arma::mat xyz2abp, abp2xyz;
			//arma::cx_mat stx2abp, abp2stx;
};

#endif
