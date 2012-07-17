#ifndef ROTATION_HPP_
#define ROTATION_HPP_

#include <armadillo>

class RotationMatrix
{
		public:
				RotationMatrix(std::vector<float>);
		private:
		public:
			arma::mat rotQ;

};

#endif
