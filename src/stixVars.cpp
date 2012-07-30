#include "stixVars.hpp"
#include "constants.hpp"
#include <math.h>
#include <vector>
#include <cmath>

PlasmaSpecies::PlasmaSpecies(double _z, double _amu, double _n, float _bMag)
{
	z = _z;
	amu = _amu;
	n = _n;
	bMag = _bMag;
	q = z * _e;	
	m = amu * _mi;
	wp = sqrt( n*pow(q,2) / (m*_e0) );
	wc = std::abs(q) * bMag / m; // Swanson definition, i.e., always positive.
	e_swan = z/std::abs(z);
}

HotPlasmaSpecies::HotPlasmaSpecies(double _z, double _amu, double _n, float _bMag, float _TeV, int _maxHarmN)
{
	z = _z;
	amu = _amu;
	n = _n;
	bMag = _bMag;
	q = z * _e;	
	m = amu * _mi;
	wp = sqrt( n*pow(q,2) / (m*_e0) );
	wc = std::abs(q) * bMag / m; // Swanson definition, i.e., always positive.
	e_swan = z/std::abs(z);
	T_eV = _TeV;
	vTh = sqrt ( 2.0 * T_eV*_e / m);
	maxHarmN = _maxHarmN;
}

StixVars::StixVars(double _omega, std::vector<PlasmaSpecies> _s)
{
	double _R = 0.0;
	double _L = 0.0;
	double _P = 0.0;

	// These are Swanson variables. It DOES matter!

	for(int s=0;s<_s.size();s++)
	{
		_R += pow(_s[s].wp,2) / ( _omega * ( _omega + _s[s].e_swan*_s[s].wc ) );
		_L += pow(_s[s].wp,2) / ( _omega * ( _omega - _s[s].e_swan*_s[s].wc ) );
		_P += pow(_s[s].wp,2) / pow(_omega,2);
	}

	R = 1.0 - _R;
	L = 1.0 - _L;
	P = 1.0 - _P;
	S = 0.5 * ( R + L );
	D = 0.5 * ( R - L );
}
