#include "stixVars.hpp"
#include "constants.hpp"
#include <math.h>
#include <vector>

PlasmaSpecies::PlasmaSpecies(double _z, double _amu, double _n, float _bMag)
{
	z = _z;
	amu = _amu;
	n = _n;
	bMag = _bMag;
	q = z * _e;	
	m = amu * _mi;
	wp = sqrt( n*pow(q,2) / (m*_e0) );
	wc = q * bMag / m;
}

StixVars::StixVars(double _omega, std::vector<PlasmaSpecies> _s)
{
	double _R = 0.0;
	double _L = 0.0;
	double _P = 0.0;

	for(int s=0;s<_s.size();s++)
	{
		_R += pow(_s[s].wp,2) / ( _omega * ( _omega + _s[s].wc ) );
		_L += pow(_s[s].wp,2) / ( _omega * ( _omega - _s[s].wc ) );
		_P += pow(_s[s].wp,2) / pow(_omega,2);
	}

	R = 1.0 - _R;
	L = 1.0 - _L;
	P = 1.0 - _P;
	S = 0.5 * ( R + L );
	D = 0.5 * ( R - L );
}