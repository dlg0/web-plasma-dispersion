#ifndef STIXVARS_HPP_
#define STIXVARS_HPP_
#include <vector>

class PlasmaSpecies
{
		public:
				PlasmaSpecies(double, double, double, float);
		private:
		public:
				double z, amu; // Atomic charge and mass
				double wp, wc; // Plasma and Cyclotron angular frequencies 
				double q, m, n, bMag; // Charge, mass, density, bMag


};

class StixVars
{
		public:
				StixVars(double _omega, std::vector<PlasmaSpecies> _s);
		private:
		public:
				double R, L, S, D, P;
};


#endif
