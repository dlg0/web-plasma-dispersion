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
				double e_swan; // q/|q| and wc is ALWAYS postive
};

class HotPlasmaSpecies
{
		public:
				HotPlasmaSpecies(double, double, double, float, float);
		private:
		public:
				double z, amu; // Atomic charge and mass
				double wp, wc; // Plasma and Cyclotron angular frequencies 
				double q, m, n, bMag; // Charge, mass, density, bMag
				double T_eV, vTh; // Temperature, Thermal velocity
				double e_swan; // q/|q| and wc is ALWAYS postive
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
