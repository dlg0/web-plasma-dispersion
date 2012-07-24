#ifndef ROTATION_HPP_
#define ROTATION_HPP_

#include <armadillo>

class RotationMatrix
{
		public:
				RotationMatrix(arma::colvec _b_rtz);//, arma::cx_colvec _k_rtz);
		private:
		public:
			arma::mat rtz2abp, abp2rtz;
			arma::cx_mat stx2abp, abp2stx;
};

#endif
