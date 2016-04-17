/* this is to calculate the viscousity and thermo-conductivity of the vapour phase */

// for simplicity, we may assume during the depressurisation, those properties does not change much


scalar totalArea;
scalar totalVolume;

scalar TvapInit = thermo1.T()[celln/2];
scalar TliqInit = thermo2.T()[celln/2]; 
scalar rhovapInit = rho1[celln/2];
scalar rholiqInit = rho2[celln/2];

// pipeline characteristics 
scalar length = 40.0;
scalar thickness = 0.011;
scalar id = 0.0508;
scalar od = id + 2.*thickness; 
scalar rhowall = 7950.0;
scalar cpwall = 60.0;
scalar kwall = 40.0;
scalar htcWall = 1./(0.5*id/kwall*Foam :: log((od)/id));
totalArea = 3.1415926*id*length;
totalVolume = 3.1415926*Foam :: pow (id,2)/4.*length;

// bulk properties
scalar muliq, muvap, kliq, kvap;
scalar Revap, Prvap, Reliq, Prliq;

// vapour phase heat capacity
volScalarField cpVap("cpVap", thermo1.Cpv());

// calculate the viscousity and thermoconductivity
{ 
	long ierr = 0; 
	char herr[255];
	double eta=NULL, tcx=NULL;   	
	double xvap[]={1.0}, xliq[]={1.0}, T = NULL, rho = NULL;

	T = TvapInit;
	rho = rhovapInit/44.01;
	TRNPRPdll (T,rho,xvap,eta,tcx,ierr,herr,255); // 255 is compeletly random 

	muvap = eta/1000000.; // upa.s to pa.s
	kvap = tcx; // W/m-K

	//TRNPRPdll (TliqInit,rholiqInit,xliq,eta,tcx,ierr,herr,255);

	//muliq = eta/1000000.	
	//kliq = tcx;
}


for(i = 0; i < celln+1; i++)
{
	scalar TwallS = Twall[i];
	scalar qV_wallS;
	scalar htcVapS;	

	htcVapS = 0.023*Foam :: pow((alpha1[i]*(rho1[i]*mag(U1[celln/2])*id/muvap)),0.8) 
	* Foam :: pow((cpVap[i]*muvap/kvap),0.4)*kvap/(1.0*id);

	htcVapS = 1./(1./htcVapS+1./htcWall);
	
	// heat transfer to the vapour phase from the wall
	if (TwallS > thermo1.T()[i])
	{	
		qV_wallS = 5.0*alpha1[i]*htcVapS*totalArea/totalVolume*(TwallS - thermo1.T()[i]);
	}
	else
	{
		qV_wallS = 0.0;
	}

	htcVap[i] = htcVapS;
	qV_wall[i] = qV_wallS;

	// updating pipeline wall temperature 
	TwallS = TwallS - (qV_wallS*dt)/(rhowall*cpwall*((Foam :: pow((0.5*id+thickness),2)-Foam :: pow((0.5*id),2))/(Foam :: pow (0.5*id,2))));
	Twall[i] = TwallS; 
}

//Info << "wall temperature" << Twall << nl << endl;
//Info << "heat transfer"  << qV_wall << nl << endl;
