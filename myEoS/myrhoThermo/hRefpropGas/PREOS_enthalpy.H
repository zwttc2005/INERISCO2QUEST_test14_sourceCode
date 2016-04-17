/*=====================================================*/
/*return specific enthalpy at T and rho*/
/*=====================================================*/
{
	scalar  Bp;
	scalar  dadt;
	scalar  delta_h;
	scalar  compress;	

	// define ideal gas heat capacity coefficient*/
	scalar  CC1 = 453.577;
	scalar  CC2 = 1.65014;
	scalar  CC3 = -1.24814e-3;
	scalar  CC4 = 3.78201e-7;
	scalar  CC5 = 0.00;

	// define specific universal gas constant
	scalar Rs = this->R();

	// define reference enthalpy 
	scalar  H_REF = 13192.75;

	scalar  kvar = 0.37464+1.54226*omega_-0.26992*sqr(omega_);

	scalar  a = 0.45724*sqr(Rs)*sqr(Tc_)/Pc_;   // universal gas constant is in molar base 
    	scalar  b = 0.07780*Rs*Tc_/Pc_;             // sane as above 
    	scalar  Tr = T/Tc_;
	
	scalar  alpha =
        sqr
        (
          1.0
          + kvar
          * (1.0 - sqrt(Tr))
        );
	
	scalar  afun = alpha*a; // mass base 
	
	scalar  H_ideal_gas = T*(CC1 + T*(0.5*CC2 + T*(0.333333*CC3 + T*(0.25*CC4 + T*0.2*CC5))));	

	dadt = kvar*a*((kvar / Tc_) - (1. + kvar) / sqrt(T*Tc_)); // mass base
	compress = this->Z(p, T);	

	Bp = (b*p) / (Rs*T);

	delta_h = Rs*T*(compress - 1.) + ((T*dadt - afun) / (2.828427*b))*log((compress + 2.414214*Bp) / (compress - 0.414214*Bp)); // pre kg 

	enthalpy =  (delta_h*this->W()/1000.) + (H_ideal_gas*this->W()/1000.) + H_REF;
}
