#ifndef STIXVARS_HPP_
#define STIXVARS_HPP_
#include <vector>
#include <complex>

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
				HotPlasmaSpecies(double, double, double, float, float,int);
		private:
		public:
				double z, amu; // Atomic charge and mass
				double wp, wc; // Plasma and Cyclotron angular frequencies 
				double q, m, n, bMag; // Charge, mass, density, bMag
				double T_eV, vTh; // Temperature, Thermal velocity
				double e_swan; // q/|q| and wc is ALWAYS postive
				int maxHarmN;
};


class StixVars
{
		public:
				StixVars(std::complex<double> _omega_c, std::vector<PlasmaSpecies> _s);
		private:
		public:
				std::complex<double> R, L, S, D, P;
};


#endif
