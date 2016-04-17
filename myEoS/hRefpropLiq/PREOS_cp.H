/*===============================================================*/
/*return heat capacity at T and rho*/
/*===============================================================*/
{
	scalar delta_Cp, Z, dadt, ddadt;
	scalar Ap, Bp, M, N;

	scalar compress;

	// specific universal gas constant
	scalar Rs = this->R();
	
	// define ideal gas heat capacity coefficient*/
	scalar CC1 = 453.577;
	scalar CC2 = 1.65014;
	scalar CC3 = -1.24814e-3;
	scalar CC4 = 3.78201e-7;
	scalar CC5 = 0.00;

	// calculate the redueced temperature 

	scalar kvar = 0.37464 + 1.54226*omega_- 0.26992*sqr(omega_); // dimensionless 

	scalar Tr = T / Tc_;
	scalar a = 0.45724*sqr(Rs)*sqr(Tc_)/Pc_;   // universal gas constant is in molar base 
    	scalar b = 0.07780*Rs*Tc_/Pc_;             // same as above 

	scalar alpha =
        sqr
        (
          1.0
          + kvar
          * (1.0 - sqrt(Tr))
        );
	
	scalar afun = alpha*a; // mass base 

	scalar cp_ideal_gas = (CC1 + T*(CC2 + T*(CC3 + T*(CC4 + T*CC5)))); // mass base 	

	ddadt = (kvar*a*(1. + kvar)) / (2.*sqrt(T*T*T*Tc_));   // mass base
	dadt = kvar*a*((kvar / Tc_) - (1. + kvar) / sqrt(T*Tc_)); // mass base
	
	compress = this->Z(p, T);
	
	Ap = (afun*p) / (sqr(Rs)*sqr(T)); 
	Bp = (b*p) / (Rs*T);

	M = (compress*compress + 2.*Bp*compress - Bp*Bp) / (compress - Bp);
	N = (dadt)*(Bp / (b*Rs));

	delta_Cp = ddadt*(T / (2.828427*b))*log((compress + 2.414214*Bp) / (compress - 0.414214*Bp)) + (Rs*(M - N)*(M - N)) / (M*M - 2.*Ap*(compress + Bp)) - Rs; // mass base

	heatCapacity = cp_ideal_gas*this->W()/1000. + delta_Cp*this->W()/1000.;                 /*PREOS_Cp_ideal_gas(temp) + delta_Cp;*/ /*j/kg-k*/

}
