#include <stdio.h>
#include <math.h>
#include <iostream>

struct density_rho0{
  double rho0, c_den[100], a_den[100]; //the overall normalization factor for nuclear charge densities
};

struct ZBQL0001{
  double ZBQLIX[43], B, C;
};

struct Brem_spect{
  double Eg[500], Br[500];
};


struct density_rho0 density_rho0;
struct ZBQL0001 zbql0001;
struct Brem_spect brem_spect;

double Brem(bool brem_init, bool cobrems, double E0, double Egamma);
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF, double phi_JT, double k1[3], double k2[3]);//units of nb/sr^2
double FF2(double q2, int ztgt, bool nuc_FF);
void analysis(double E0, double mtgt, double k1[3], double k2[3], double ktgt[3], double w_mumu, double t, double missing_mass, double m_part, bool pion_hypothesis);
void density_init(int ztgt);
double FF(double Q2, int ztgt);;
void ZBQLINI(int SEED);
//void ZBQLU01(double DUMMY);
//void ZBQLUAB(double A, double B);
//void ZBQLexp(double MU);
//void ZBQLNOR(double MU, double SIGMA);
//void ZBQLBIN(int N, double P);
//void ZBQLGEO(double P);
//void ZBQLPOI(double MU);
//void ZBQLGAM(double G, double H);
//void ZBQLBET1(double NU1, double NU2);
//void ZBQLWEI(double A, double B);
//void ZBQLNB(double R, double P);
//void ZBQLPAR(double A, double B);
//void ZBQLLG(double X);


//
//	To compile: "gfortran -ffixed-line-length-132 lepton_event_v17_4.f -o lepton_event_v17_4.exe”
//	To run: “./lepton_event_v17_4.exe”
//
//	version 6 of code weights the distribution of events by xs= acos( costheta ).  This is the most efficient version
//	of the code to run.  The nuclear and atomic form factors are included. 
//	version 8 allows for normal event weighting, or non-standard event weighting, theta=x**2, switchable with ‘standard_phase_space’. 
//		Non-standard is the most efficient. 
//	version 9: histograms J_T phi angle, uses simplified version of W_pol. 
//	version 10: histograms cross section arrays
//	version 11: histogram integrated cross section as a function of x and phi1. 
//	version 12: introduced logicals to do electron, muon, and muon with pion mass hypothesis histogramming. 
//	version 13: added some histogram options and controls
//	version 14: added some more controls, got rid of things in the code that weren’t useful 
//	version 15: brought over some features from the pi+pi- code, fixed some bugs.   Can control shape of the proton form factor. 
//	version 16: consistently use transverse momentum transfer squared everywhere, changed name to lepton_event
//	Uses rejection sampling: 
//	version 17: added feature to run one kinematic point
//	version 17_1: use weighting dcos theta/dx = (1-cos theta)**N , got rid of the option for bending the FF
//	version 17_2: use weighting theta = x**phase_space, with phase_space an integer >1 .   Set phase_space = 0 for dcos theta/dx =1. 
//		      makes a guess for the largest cross section, which seems to be correct
//	version 17_3: add option for log t plot, control proton rms radius independent of dipole FF parameter
//	version 17_4: read brem. file
//
//
double E0, theta_min, theta_max;
double pi, cos_max, cos_min, cross_sum, cross_max, cross_test, x_min, x_max, theta1, phi1, costheta1, theta2, phi2, costheta2, cross_section, total_xscn, k1[3], k2[3], m_e, m_part, m_muon, pol, x, q2, failure, w_mumu, xs1, xs2, xs_max, xs_min, Atgt[100], tgtlen[100], radlen[100], hbarc, W, jacobian, m_pi, phi_JT, delta_w, delta_x, delta_phi, mtgt, ktgt[3], Elab[2], missing_mass, t, delta_t, x_value, q2_T, W_pol, W_unpol, JS, JT[2], Rexp, xmax, theta1_max, theta2_max, phi1_max, phi2_max, temp, delta_log_t, frac_delta_t, delta_Egamma, data_array[200], error_array[200], Egamma_max, E_hi, E_lo, E_coherent, Egamma;
int iseed;
int i, itest[4], nevent, j, nfail, bad_max, i_array, j_array, ztgt, phase_space; 
bool hist_w, hist_x, hist_t, hist_phi_JT, hist_Egamma, output_event, hist_log_t, nuc_FF, muon, electron, pion_hypothesis, brem_init, cobrems;
//
//	Standard CPP configuration 
//	double ztgt = 82, E0 = 5.5, pol = 1.0, theta_min = 0.80,theta_max = 5.3; //Standard CPP configuration, min angle to TOF and max angle to MWPC
//
//	Standard GlueX configuration
//	data ztgt,E0,E_coherent,pol,theta_min,theta_max /1,11.0,8.7,1.0,0.90,13.12/	//standard GlueX config., min and max angles in deg. to TOF
int main(){
    ztgt = 1, E0 = 11.0, E_coherent = 8.7, pol = 1.0, theta_min = 0.090, theta_max = 13.12;///standard GlueX config., min and max angles in deg. to TOF;
//
//	Set tagging interval
    E_hi = 8.8, E_lo = 8.6;
//
    itest[0] = 100000, itest[1] = 1000000, itest[2] = 10000000, itest[3] = 100000000, nevent = 10;
    m_e = 0.000511, m_muon = 0.105658, m_pi = 0.139570, hbarc = 0.197;
//	
//  Histogram parameters
    delta_w = 0.02, delta_t = 0.0002, delta_x = 0.02, delta_phi = 5.0, frac_delta_t = 0.2, delta_Egamma = 0.05; //units of GeV, GeV^2, //DIMENSIONless, degrees,
//											fractional bin width in t, GeV
//
//  Target information
    Atgt[0] = 1.0, tgtlen[0] = 0.0338, radlen[0] = 63.04;//tgtlen = # Rad lengths, radlen = rad length of material in g/cm^2;
//								this target is based on 30 cm LH2
    Atgt[5] = 12.0, tgtlen[5] = 0.05, radlen[5] = 42.70, density_rho0.c_den[5] = 2.45, density_rho0.a_den[5] = 0.524;//RL is in g/cm^2;
    Atgt[13] = 28.0, tgtlen[13] = 0.05, radlen[13] = 21.82, density_rho0.c_den[13] = 3.14, density_rho0.a_den[13] = 0.537;//units of c_den and a_den are fm;
    Atgt[19] = 40.0, tgtlen[19] = 0.05, radlen[19] = 16.12, density_rho0.c_den[19] = 3.51, density_rho0.a_den[19] = 0.563;//target length is in units of RL;
    Atgt[25] = 56.0, tgtlen[25] = 0.05, radlen[25] = 13.84, density_rho0.c_den[25] = 3.971,density_rho0.a_den[25] = 0.5935;//RL is g/cm^2;
    Atgt[49] = 116.0, tgtlen[49] = 0.05, radlen[49] = 8.82, density_rho0.c_den[49] = 5.416, density_rho0.a_den[49] = 0.552;
    Atgt[81] = 208.0, tgtlen[81] = 0.05, radlen[81] = 6.37;
//
//COMMON/density_rho0/rho0, c_den, a_den//the overall normalization factor for nuclear charge densities
//
    double ZBQLUAB, zlo, zhi;
    double zlo = 0.0, zhi = 1.0;
//
//      BLOCK DATA ZBQLBD01
//
//       Initializes seed array etc. for random number generator.
//       The values below have themselves been generated using the
//       NAG generator.
//
    double ZBQLIX[43], B, C;
//COMMON/ZBQL0001/ZBQLIX, B, C
//      INTEGER I
    ZBQLIX[0] = 8.001441;
    ZBQLIX[1] = 5.5321801;
    ZBQLIX[2] = 1.69570999;
    ZBQLIX[3] = 2.88589930;
    ZBQLIX[4] = 2.91581871;
    ZBQLIX[5] = 1.03842493;
    ZBQLIX[6] = 7.9952507;
    ZBQLIX[7] = 3.81202335;
    ZBQLIX[8] = 3.11575334;
    ZBQLIX[9] = 4.02878631;
    ZBQLIX[10] = 2.49757109;
    ZBQLIX[11] = 1.15192595;
    ZBQLIX[12] = 2.10629619;
    ZBQLIX[13] = 3.99952890;
    ZBQLIX[14] = 4.12280521;
    ZBQLIX[15] = 1.33873288;
    ZBQLIX[16] = 7.1345525;
    ZBQLIX[17] = 2.23467704;
    ZBQLIX[18] = 2.82934796;
    ZBQLIX[19] = 9.9756750;
    ZBQLIX[20] = 1.68564303;
    ZBQLIX[21] = 2.86817366;
    ZBQLIX[22] = 1.14310713;
    ZBQLIX[23] = 3.47045253;
    ZBQLIX[24] = 9.3762426;
    ZBQLIX[25] = 1.09670477;
    ZBQLIX[26] = 3.20029657;
    ZBQLIX[27] = 3.26369301;
    ZBQLIX[28] = 9.441177;
    ZBQLIX[29] = 3.53244738;
    ZBQLIX[30] = 2.44771580;
    ZBQLIX[31] = 1.59804337;
    ZBQLIX[32] = 2.07319904;
    ZBQLIX[33] = 3.37342907;
    ZBQLIX[34] = 3.75423178;
    ZBQLIX[35] = 7.0893571;
    ZBQLIX[36] = 4.26059785;
    ZBQLIX[37] = 3.95854390;
    ZBQLIX[38] = 2.0081010;
    ZBQLIX[39] = 5.9250059;
    ZBQLIX[40] = 1.62176640;
    ZBQLIX[41] = 3.20429173;
    ZBQLIX[42] = 2.63576576;

    double B = 4.294967291;
    double C = 0.0;
//
    iseed = 0;//set equal to 0 for the // program to do a call to the system clock to initialize random number generator;
//			!for random number initialization
    ZBQLINI(iseed);
//
//  Initializations
//
    temp = -1;
    pi = acos(temp);
//
//Find the delta log t step
    delta_log_t = log10(1. + frac_delta_t);

//Set target mass in GeV
    mtgt = Atgt[ztgt]*0.931494;
    if (ztgt == 1) mtgt = 0.93828;
//
    do{
        data_array[i] = 0;
        error_array[i] = 0;
        i++;
    }while(i < 200);

//Initialize Brem. distribution: select 1/Egamma or coherent Brems. file
    cobrems = true;//set true for scanfing coherent Brems. file, false for using a 1/Egamma distribution;
    if(cobrems == true) 
{ // scanf coherent Brem. file
    brem_init = true;
    temp = Brem(brem_init, cobrems, E0, Egamma);//scanf coherent brems file, then set brem_init = false;
}

//Start logical assignments

//Only one of the histogram logicals can be true, the rest must be false.   The event output can be turned on independent of histograming. 
    hist_w = false;
    hist_x = false;
    hist_t = false;
    hist_phi_JT = false;
    hist_log_t = false;
    hist_Egamma = false;
    output_event = true;
//Logical assignments
    phase_space = 4;//theta = x**phase_space, with int phase_space >1 . Note: phase_space = 4 seems to be fastest.;
//Set phase_space=0 for standard dcos theta/dx =1
    Rexp = float(phase_space);
    nuc_FF = true;
    muon = false;//this is how you change the particle type
//electron = .not.muon;
    if (muon == true) 
    {
        m_part = m_muon;
        pion_hypothesis = false; //set true if you want calculated invariant masses with pion assumption;
    }else{
        m_part = m_e;
        pion_hypothesis = false;//should always have this set false since we can distinguish e from mu;
    }
//
    if (nuc_FF) 
    {// setup the nuclear form factors
        density_rho0.c_den[ztgt] = density_rho0.c_den[ztgt]/hbarc;
        density_rho0.a_den[ztgt] = density_rho0.a_den[ztgt]/hbarc;
        density_init(ztgt); //for initializing the nuclear form factor(density_rho0.rho0 is density const.);
    }
//
    if((hist_w) || (hist_x) || (hist_t) || (hist_phi_JT) || (hist_log_t) || (hist_Egamma)) fopen("lepton_v17_4_hist.txt","w");
    if (output_event) fopen("lepton_v17_4_event.txt", "w");
//
//	End logical assignments
//
//***********************************
//	Evaluate the cross section at one kinematic point at coherent peak
    x = 0.5;
    theta1 = theta_min*(pi/180);
    theta2 = theta_min*(pi/180);
    phi1 = 90*(pi/180);
    phi2 = 270*(pi/180);
    cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, phi_JT, k1, k2);
    std::cout << " cross section nb/sr^2 = " << cross_section << "\n";
//*************************************
    theta_min = theta_min*pi/180;//switch to radians;
    theta_max = theta_max*pi/180;
    cos_max = cos(theta_min);
    cos_min = cos(theta_max);
// Limits on xs
    if (phase_space == 0){//dcos theta/dx = 1
        xs_max = cos_max;
        xs_min = cos_min;
    }else{
        xs_max = pow(theta_max, 1/Rexp);
        xs_min = pow(theta_min, 1/Rexp);
    }
//
//***************************************************************************
//
    j = 0;
    i = 0;
    do{//loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the maximum
//cross section*Brem converges
    cross_max = 0;
    do{//find maximum cross section in allowed phase space
        Egamma = ZBQLUAB(E_lo, E_hi);//get tagged photon energy
        x = 0.5;//x_min + (x_max - x_min)*ZBQLUAB(zlo, zhi)//make a guess for the energy fraction
        phi1 = 90*pi/180;//2.*pi*ZBQLUAB(zlo, zhi)//make a guess for phi1
        phi2 = 270*pi/180;//2.*pi*ZBQLUAB(zlo, zhi)//make a guess for phi2
        xs2 = ZBQLUAB(xs_min, xs_max);
        if (phase_space == 0) {// dcos theta/dx = 1
            theta1 = theta_min;//make a guess for theta1
            theta2 = acos(xs2);
            jacobian = 1;
        }else{
            theta1 = theta_min;//make a guess for theta1
            xs1 = pow(theta_min, (1/Rexp));
            theta2 = pow(xs2, phase_space);
            jacobian = (Rexp * pow(xs1, (phase_space - 1)) * sin(pow(xs1, phase_space))) * (Rexp * pow(xs2,(phase_space - 1)) * sin(pow(xs2,phase_space)));
        }
        cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, phi_JT, k1, k2)*jacobian*Brem(brem_init, cobrems, E0, Egamma);
        if(cross_section > cross_max){
            cross_max = cross_section;
            Egamma_max = Egamma;
            xmax = x;
            theta1_max = theta1*180/pi;
            theta2_max = theta2*180/pi;
            phi1_max = phi1*180/pi;
            phi2_max = phi2*180/pi;
            i++;
        }
    }while(i < itest[j]);
    std::cout << "test events " << itest[j] <<  " maximum xsctn*Brem " <<  cross_max << "\n";
    std::cout << "Egamma max " << Egamma_max << " x max " << xmax << " theta1 max " << theta1_max << " theta2 max " << theta2_max << " phi1 max " << phi1_max << " phi2 max " <<  phi2_max + "\n";
    j++;
    }while(j < 4);
//
//**********************************************************************************************************
// Loop over 4 samplings of the phase space at coherent peak, each a factor of x10 larger, to see if the integrated cross section converges
// Set limits on energy fraction
    j = 0;
    i = 0;
    x_max = (E_coherent - m_part)/E_coherent;
    x_min = m_part/E_coherent;
//
    do{//loop over 4 samplings of the phase space, each a factor of x10 larger, to see if the integrated
//cross section at the coherent peak converges
    cross_sum = 0;
    do{
        x = x_min + (x_max - x_min)*ZBQLUAB(zlo, zhi);//energy fraction
        phi1 = 2*pi*ZBQLUAB(zlo, zhi);
        phi2 = 2*pi*ZBQLUAB(zlo, zhi);
        xs1 = ZBQLUAB(xs_min, xs_max);
        xs2 = ZBQLUAB(xs_min, xs_max);
        if (phase_space == 0){// dcos theta/dx = 1
            theta1 = acos(xs1);
            theta2 = acos(xs2);
            jacobian = 1;
        }else{
            theta1 = pow(xs1,phase_space);
            theta2 = pow(xs2,phase_space);
            jacobian = (Rexp * pow(xs1,(phase_space - 1)) * sin(pow(xs1,phase_space))) * (Rexp*pow(xs2,(phase_space-1))*sin(pow(xs2,phase_space)));
        }
        cross_section = xsctn(E_coherent, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, phi_JT, k1, k2)*jacobian;
        cross_sum = cross_sum + cross_section;
        i++;
    }while(0 < itest(j));
    total_xscn = cross_sum/float(itest[j])*pow((xs_max - xs_min), 2)*pow((2*pi), 2) * (x_max - x_min);
    std::cout << "test events " << itest[j] << " Egamma " <<  E_coherent << " total cross section nb " <<  total_xscn << "\n";
    j++;
    }while(j < 4);


//*************************************************************************
//
//	Start event generation
//
//Use the widest possible range in x by using the maximum accepted tagged photon energy, then test it
    x_max = (E_hi - m_part)/E_hi;//largest possible x
    x_min = m_part/E_hi;//smallest possible x
//
    nfail = 0;
    bad_max = 0;
//
    do{
    g100:

        Eg2.63576576D8amma = ZBQLUAB(E_lo, E_hi);//get tagged photon energy
        x = ZBQLUAB(x_min, x_max);//energy fraction
  //	Test x to make sure it's within the allowed range for the photon energy Egamma
        if (x >= ((Egamma - m_part)/Egamma)) || (x <= (m_part/Egamma)){goto g100} // x is out of range, try again
        phi1 = 2.*pi*ZBQLUAB(zlo, zhi);
        phi2 = 2.*pi*ZBQLUAB(zlo, zhi);
        xs1 = ZBQLUAB(xs_min, xs_max);
        xs2 = ZBQLUAB(xs_min, xs_max);
//
        if(phase_space == 0){// dcos theta/dx = 1
            theta1 = acos(xs1);
            theta2 = acos(xs2);
            jacobian = 1;
        }else{
            theta1 = pow(xs1, phase_space);
            theta2 = pow(xs2, phase_space);
            jacobian = (Rexp * pow(xs1,(phase_space_1)) * sin(pow(xs1,phase_space))) * (Rexp*pow(xs2,(phase_space-1)) * sin(pow(xs2,phase_space)));
        }

        cross_section = xsctn(Egamma, ztgt, x, theta1, phi1, theta2, phi2, pol, m_part, nuc_FF, phi_JT, k1, k2)*jacobian*Brem(brem_init, cobrems, E0, Egamma);
        if (cross_section > cross_max){
            bad_max = bad_max + 1;//an occurrence of cross section larger than cross_max, not supposed to happen
            std::cout <<  "bad max cross section= " <<  cross_section << "\n";
        }
        cross_test = cross_max*ZBQLUAB(zlo, zhi);
        if (cross_test > cross_section){//selection fails
            nfail = nfail + 1;
            goto g100;
        }
//	Event selection succeeds: 
//
        analysis(Egamma, mtgt, k1, k2, ktgt, w_mumu, t, missing_mass, m_part, pion_hypothesis);//analyze the event;
//							
//		Do the histogramming
//
        if(hist_w) i_array = nint((w_mumu - 0.200)/delta_w);//w distribution;
        if(hist_x) i_array = nint(x/delta_x);//x distribution
        if(hist_t) i_array = nint(t/delta_t);//t distribution;
        if(hist_phi_JT) i_array = nint(phi_JT*180/pi/delta_phi); //JT phi distribution in degrees;
        if(hist_log_t) i_array = nint((log10(t) + 6)/delta_log_t);// 10^ - 6 GeV^2 is bin 0;
        if(hist_Egamma) i_array = nint(Egamma/delta_Egamma);//photon energy distribution;
        if(i_array < 0) i_array = 0;
        if (i_array > 200) i_array = 200;
        data_array[i_array] = data_array[i_array] + 1;
//
//	3-momentum event output
        if(output_event) {printf(4, 200) Egamma, k1[1], k1[2], k1[3], k2[1], k2[2], k2[3], ktgt[1], ktgt[2], ktgt[3]};
g200:
// format(2x, f6.3, 1x, 9(f10.6, 1x))
//
        if (mod(i, 100) == 0) {std::cout << ' event # ' << i}
    }while(i < nevent);


//
//   Event generation ends
//	
    float failure = float(nfail)/float(nevent);
    std::cout <<  "Failures per event = " << failure << " Events with cross section exceeding max xsctn = " << bad_max << "\n";
//
    int i = 0;
    do{//printf out the arrays
        if(hist_w){
        x_value = float(i)*delta_w + 0.200;
        error_array[i] = sqrt(data_array[i]);
        printf(2, 20) x_value, data_array[i], error_array[i];  // format(2x, e10.4, 2x, e10.4, 2x, e10.4)
        }else if (hist_x){
            x_value = float(i)*delta_x;
            error_array[i] = sqrt(data_array[i]);
            printf(2, 20) x_value, data_array[i], error_array[i];  // format(2x, e10.4, 2x, e10.4, 2x, e10.4)
        }else if(hist_t){
            x_value = float(i)*delta_t;
            error_array[i] = sqrt(data_array[i]);
            printf(2, 20) x_value, data_array[i], error_array[i];  // format(2x, e10.4, 2x, e10.4, 2x, e10.4)
        }else if(hist_phi_JT){
            x_value = float(i)*delta_phi;
            if((x_value == 0.) || (x_value == 360.)) data_array[i] = 2.*data_array[i]; // this is a binning problem
            error_array[i] = sqrt(data_array[i]);
            printf(2, 20) x_value, data_array[i], error_array[i];  // format(2x, e10.4, 2x, e10.4, 2x, e10.4)
        }else if(hist_log_t){
            x_value = pow(10,(float[i]*delta_log_t-6));
            error_array[i] = sqrt(data_array[i])/(pow(10,(float[i+1]*delta_log_t-6)) - pow(10,(float[i-1]*delta_log_t - 6))) * 2;
            data_array[i] = data_array(i)/(pow(10,(float[i+1]*delta_log_t-6)) - pow(10,(float[i-1]*delta_log_t - 6))) * 2;
            printf(2, 20) x_value, data_array[i], error_array[i];  // format(2x, e10.4, 2x, e10.4, 2x, e10.4)
        }else if(hist_Egamma){
            x_value = float(i)*delta_Egamma;
            error_array[i] = sqrt(data_array[i]);
            printf(2, 20) x_value, data_array[i], error_array[i];  // format(2x, e10.4, 2x, e10.4, 2x, e10.4)
        }
    }while(i < 200);


//

//
exit(0);
};


//
//*******************************************************************
//
// --------------------------------------------
double Brem(bool brem_init, bool cobrems, double E0, double Egamma){
    double Eg[500], Br[500];
    int i, imax, ipoint;
    bool brem_init, cobrems;
    //COMMON/Brem_spect/Eg, Br
    //	
    //Brem = 0;
    //
    if (brem_init){//fopen and scanf coherent Brem file
      fopen(unit = 2, file = "CobremsDistribution.dat");
      i = 1;
g10:
        
      scanf(2, *, end = 20) brem_spect.Eg[i], brem_spect.Br[i];
        //		print *, Eg(i), Br(i)
      i = i + 1;
      goto g10;
g20:
        
      imax = i - 1;
      close(2);
      brem_init = false//done with initialization;
      return;
    }
    //
    if (.not.cobrems) 
    {
        Brem = E0/Egamma;
        return;
    }
    //
    if (cobrems) 
    { //return coherent brems distribution
        ipoint = nint((Egamma + .02)/.04);
        Brem = brem_spect.Br[ipoint];
        return;
    }
}


//	
//******************************************************************
// The reference for this is Bakmaev et al., Physics Letters B 660 (2008) 494-500, eqn. 23
// Had to correct a mistake in Eqn. 23 of their paper.   The cos(phi1 + phi2) term should be 
// multiplied by 2.  You can see this by comparing Wp in Eqn. 22 with the vector current part of Eqn. 23
//
// --------------------------------------------
double xsctn(double E0, int ztgt, double x, double theta1, double phi1, double theta2, double phi2, double pol, double m_part, bool nuc_FF, double phi_JT, double k1, double k2)//units of nb/sr^2;
{
    //implicit none
    double Z, k1[3], k2[3], W_unpol, W_pol, q2_T;
    double  alpha, hbarc;
    double pi, E1, E2, k1_mag, k2_mag, p1, p2, q2, c1, c2, JS, JT[2], FF2, FF_nuc, FF_TFM, phi_JT;
    int i;
    //
    double alpha = 7.297352e-3, hbarc = 0.197;
    //
    pi = acos(-1);
    //
    Z = float(ztgt);
    E1 = E0 * x;
    E2 = E0 * (1 - x);
    k1_mag = sqrt(pow(E1, 2) - pow(m_part, 2));
    k2_mag = sqrt(pow(E2, 2) - pow(m_part, 2));
    k1[1] = k1_mag * sin(theta1) * cos(phi1);
    k1[2] = k1_mag * sin(theta1) * sin(phi1);
    k1[3] = k1_mag * cos(theta1);
    k2[1] = k2_mag * sin(theta2) * cos(phi2);
    k2[2] = k2_mag * sin(theta2) * sin(phi2);
    k2[3] = k2_mag * cos(theta2);
    //
    p1 = sqrt(pow(k1[1], 2) + pow(k1[2],2));//transverse momenta of muon #1, GeV
    p2 = sqrt(pow(k2[1],2) + pow(k2[2],2));//transverse momenta of muon #2, GeV
    q2_T = pow((k1[1] + k2[1]),2) + pow((k1[2] + k2[2]),2); // this is transverse momentum transfer squared
    q2 = q2_T + pow((E0 - k1[3] - k2[3]), 2); // this is 3 - momentum transfer squared
    c1 = pow(p1, 2) + pow(m_part, 2);
    c2 = pow(p2, 2) + pow(m_part, 2);
    JS = 1/c1 - 1/c2;//units of 1/GeV^2//scalar current, units of GeV^ - 2;
    i = 0;
    do{
      JT[i] = k1[i]/c1 + k2[i]/c2;//vector current, units of GeV^ - 1;
      i++;
    }while(i < 2);


    phi_JT = acos(JT[1]/sqrt(pow(JT[1],2) + pow(JT[2], 2)));//phi angle of JT wrt to x axis, radians
    if (JT[2] < 0){phi_JT = 2.*pi - phi_JT};
//
    W_unpol = pow(m_part,2) * pow(JS,2) + (pow(x,2) + pow((1 -x),2)) * (pow(JT[1],2) + pow(JT[2], 2));
//     	W_pol = -2.*x*(1.-x)*((p2/c2)**2*cos(2.*phi2)+(p1/c1)**2*cos(2.*phi1)+2.*(p2/c2)*(p1/c1)*cos(phi1+phi2)) !the Bakmaev expression
//	xsctn=2.*alpha**3*Z**2*E0**4*x**2*(1.-x)**2/(pi**2*q2_T**2)*(W_unpol+pol*W_pol) ! note the absence of cos(2phi_JT) in this expression
//     &	*hbarc**2/100.*1.e9*FF2(q2_T,ztgt,nuc_FF) !units of nb/sr^2  The denominator uses the transverse 3-momentum transfer^2, 
    W_pol = -2 * x * (1 - x) * (pow(JT[1], 2) + pow(JT[2], 2));//this is my reduction of the Bakmaev equations;
    xsctn = 2 * pow(alpha,3) * pow(Z, 2) * pow(E0, 4) * pow(x, 2) * pow((1 - x), 2)/(pow(pi,2) * pow(q2_T, 2) * (W_unpol + pol * cos(2 * phi_JT) * W_pol));
//this contains the cos(2phi_JT) term*hbarc**2/100.*1.e9*FF2(q2_T, ztgt, nuc_FF) //units of nb/sr^2 The denominator uses the transverse 3 - momentum transfer^2
//
    return;
};


//
//****************************************************************************************
//
// --------------------------------------------
double FF2(double q2, int ztgt, bool nuc_FF);
{
    double hbarc, z, FF_nuc, FF_TFM, alpha[3], b[3], b0, m_e, c, FF;
    int i;
    //bool nuc_FF;
    double m_e = 0.511e-3, hbarc = 0.197;
    double alpha[1] = 0.1, alpha[2] = 0.55, alpha[3] = 0.35, b[1] = 6.0, b[2] = 1.2, b[3] = 0.3;
    z = float(ztgt);
    c = m_e*pow(z,0.333);
    FF_nuc = FF(q2, ztgt);
    FF_TFM = 1;
    int i = 0;
    do{
      FF_TFM = FF_TFM - alpha[i]*q2/(q2 + pow((b[i]*c),2));
      i++;
    }while{i < 3);
    if(nuc_FF) 
      {
	FF2 = pow((FF_nuc - FF_TFM),2);
      }else
      {
	FF2 = 1;
      }
    return;
};
//
//****************************************************************************************
//
// --------------------------------------------
void analysis(double E0, double mtgt, double k1, double k2, double ktgt, double w_mumu, double t, doublle missing_mass, double m_part, bool pion_hypothesis)
{
    // implicit none
    bool pion_hypothesis;
    double  k1[3], k2[3], w_mumu, t, m_part, ktgt[3];
    double m_pi = 0.139570;
    pi = acos(-1);
    E1 = sqrt(pow(k1[1], 2) + pow(k1[2], 2) + pow(k1[3], 2) + pow(m_part, 2));//lepton energies
    E2 = sqrt(k2[1]**2 + k2[2]**2 + k2[3]**2 + m_part**2);
    ks[1] = k1[1] + k2[1]//lepton summed momentum;
    ks[2] = k1[2] + k2[2];
    ks[3] = k1[3] + k2[3];
    ktgt[1] = -ks[1]//target momentum;
    ktgt[2] = -ks[2];
    ktgt[3] = E0 - ks[3];
    missing_mass = sqrt(pow((E0 + mtgt - E1 - E2), 2) - pow(ktgt[1], 2) - pow(ktgt[2],2) - pow(ktgt[3], 2));
    t = pow(ks[1], 2) + pow(ks[2], 2) + pow((E0 - ks[3]), 2) - pow((E0 - E1 - E2), 2);//4 - momentum transfer squared to nucleus, this is positive;
    //
    //		mu mu invariant mass, possibly with pion hypothesis
    m_x = m_part;
    if (pion_hypothesis){m_x = m_pi}
    E1 = sqrt(pow(k1[1], 2) + pow(k1[2],2) + pow(k1[3], 2) + pow(m_x, 2));//need to put in the mass hypothesis
    E2 = sqrt(pow(k2[1], 2) + pow(k2[2], 2) + pow(k2[3], 2) + pow(m_x, 2));
    w_mumu = sqrt(pow((E1 + E2), 2) - pow(ks[1], 2) - pow(ks[2], 2) - pow(ks[3],2));
    //
    return;
};





//c******************************************************************
//c
// --------------------------------------------
void density_init(int ztgt)
{
    // implicit none
    double pi, c_den[100], a_den[100],rho0, w;
    int i;
    //COMMON/density_rho0/rho0, c_den, a_den
    //c
    pi = acos(-1);
    //c
    rho0 = 0;
    if((ztgt == 82) || (ztgt == 1)) return;
    //c	These equations have to do with Fermi distribution, reference? 
    w = 4 * pi * c_den[ztgt]/3 * (pow((pi * a_den[ztgt]), 2) + pow(c_den[ztgt], 2));
    do{
      w = w + 8 * pi * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * exp(-float(i)*c_den[ztgt]/a_den[ztgt])/pow(float(i), 3);
    i++;
    }while(i < 10);
    density_rho0.rho0 = 1/w;
    return;
};


//
//******************************************************************

// --------------------------------------------
double FF(double Q2, int ztgt);
{
    // implicit none
    double Q2, q02, hbarc, Q, gamma, r[12], A[12], rho0, c_den[100], a_den[100], pi, norm, proton_rms;
    int i, ztgt;
    //c
    double q02 = 0.71, proton_rms = 0.879;//proton dipole form factor parameter GeV^2, proton rms radius fm
    double hbarc = 0.197;
    //COMMON/density_rho0/rho0, c_den, a_den
    //c
    double r[0] = 0.1, A[0] = 0.003845;
    double r[1] = 0.7, A[1] = 0.009724;
    double r[2] = 1.6, A[2] = 0.033093;
    double r[3] = 2.1, A[3] = 0.000120;
    double r[4] = 2.7, A[4] = 0.083107;
    double r[5] = 3.5, A[5] = 0.080869;
    double r[6] = 4.2, A[6] = 0.139957;
    double r[7] = 5.1, A[7] = 0.260892;
    double r[8] = 6.0, A[8] = 0.336013;
    double r[9] = 6.6, A[9] = 0.033637;
    double r[10] = 7.6, A[10] = 0.018729;
    double r[11] = 8.7, A[12] = 0.000020;
    double gamma = 1.388;
    //
    //  Select the FF
    //
    pi = acos(-1);
    q = sqrt(Q2);
    //
    if (ztgt == 1) 
    {//proton
        FF = 1./(1. + Q2/q02)**2 + 2.*Q2/q02 - 1./6.*Q2*proton_rms**2/hbarc**2;
    }else if(ztgt == 82){//lead
        FF = 0;
	do{
	  FF = FF + A[i] * (pow(gamma,2) * cos(Q * r[i]/hbarc) + 2 * r[i] * hbarc/Q * sin(Q * r[i]/hbarc))/(pow(gamma, 2) + 2 * pow(r[i], 2)) * exp(-Q2/4 * pow(gamma, 2)/pow(hbarc, 2));
	  i++;
	}(while i < 12);
    }else{ //for everything else use 2 - parameter fermi model, reference ?
      FF = 4 * pow(pi,2) * rho0 * pow(a_den[ztgt], 3)/(pow((q * a_den[ztgt]), 2) * pow((sinh(pi * q * a_den[ztgt])),2)) * (pi * q * a_den[ztgt] * cosh(pi * q * a_den[ztgt]) * sin(q * c_den[ztgt]) - q * c_den[ztgt] * cos(q * c_den[ztgt]) * sinh(pi * q * a_den[ztgt]));
      i = 0;
      do{
	FF = FF + 8 * pi * rho0 * pow(a_den[ztgt], 3) * pow((-1), (i - 1)) * float(i) * exp(-float(i) * c_den[ztgt]/a_den[ztgt])/pow((pow(float(i), 2) + pow((q * a_den[ztgt]),2)),2);
	i++;
      }while(i < 10);
    }
return;
};


//
//*******************************************************************
//********	FILE: randgen.f				***********
//********	AUTHORS: Richard Chandler		***********
//********		 (richard@stats.ucl.ac.uk)	***********
//********		 Paul Northrop 			***********
//********		 (northrop@stats.ox.ac.uk)	***********
//********	LAST MODIFIED: 26/8/03			***********
//********	See file randgen.txt for details	***********
//*******************************************************************


//******************************************************************
//******************************************************************
// --------------------------------------------
void ZBQLINI(int SEED)
{
    //******************************************************************
    //*       To initialize the random number generator - either
    //*       repeatably or nonrepeatably. Need double precision
    //*       variables because integer storage can't handle the
    //*       numbers involved
    //******************************************************************
    //*	ARGUMENTS
    //*	=========
    //*	SEED	(integer, input). User-input number which generates
    //*		elements of the array ZBQLIX, which is subsequently used 
    //*		in the random number generation algorithm. If SEED=0,
    //*		the array is seeded using the system clock if the 
    //*		FORTRAN implementation allows it.
    //******************************************************************
    //*	PARAMETERS
    //*	==========
    //*	LFLNO	(integer). Number of lowest file handle to try when
    //*		opening a temporary file to copy the system clock into.
    //*		Default is 80 to keep out of the way of any existing
    //*		open files (although the program keeps searching till
    //*		it finds an available handle). If this causes problems,
    //*               (which will only happen if handles 80 through 99 are 
    //*               already in use), decrease the default value.
    //******************************************************************
    const int LFLNO = 80;
    
    //******************************************************************
    //*	VARIABLES
    //*	=========
    //*	SEED	See above
    //*	ZBQLIX	Seed array for the random number generator. Defined
    //*		in ZBQLBD01
    //*	B,C	Used in congruential initialisation of ZBQLIX
    //*	SS,MM,}	System clock secs, mins, hours and days
    //*	HH,DD }
    //*	FILNO	File handle used for temporary file
    //*	INIT	Indicates whether generator has already been initialised
    //*
    int SS, MM, HH, DD, FILNO, I;
    int INIT;
    double zbql0001.ZBQLIX(43), zbql0001.B, zbql0001.C;
    double TMPVAR1, DSS, DMM, DHH, DDD;
    
    //COMMON/ZBQL0001/ZBQLIX, B, C
    SAVE INIT;
    
    //*
    //*	Ensure we don't call this more than once in a program
    //*
    if (INIT >= 1) 
    {
        if (INIT == 1) 
        {
            printf(*, 1)  // format(//5X, "****WARNING**** You have called routine ZBQLINI ", "more than", /5X, "once. I""m ignoring any subsequent calls.", //);
            INIT = 2;
        }
        return;
    } else
    {
        INIT = 1;
    }
    //*
    //*       If SEED = 0, cat the contents of the clock into a file
    //*       and transform to obtain ZQBLIX(1), then use a congr.
    //*       algorithm to set remaining elements. Otherwise take
    //*       specified value of SEED.
    //*
    //*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //*>>>>>>>	NB FOR SYSTEMS WHICH DO NOT SUPPORT THE  >>>>>>>
    //*>>>>>>>	(NON-STANDARD) 'CALL SYSTEM' COMMAND,    >>>>>>>
    //*>>>>>>>	THIS WILL NOT WORK, AND THE FIRST CLAUSE >>>>>>>
    //*>>>>>>>	OF THE FOLLOWING IF BLOCK SHOULD BE	 >>>>>>>
    //*>>>>>>>	COMMENTED OUT.				 >>>>>>>
    //*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (SEED == 0) 
    {
        //*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        //*>>>>>>>	COMMENT OUT FROM HERE IF YOU DON'T HAVE  >>>>>>>
        //*>>>>>>>	'CALL SYSTEM' CAPABILITY ...		 >>>>>>>
        //*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        SYSTEM(" date +%S%M%H%j > zbql1234.tmp");
        //*
        //*       Try all file numbers for LFLNO to 999 
        //*
        FILNO = LFLNO;
g10:
        fopen(FILNO, FILE = "zbql1234.tmp", ERR = 11);
        goto g12;
g11:
        FILNO = FILNO + 1;
        if (FILNO > 999) 
        {
            printf(*, 2)  // format(//5X, "**** ERROR **** In routine ZBQLINI, I couldn""t", " find an", /5X, "available file number. To rectify the problem, decrease the ", "value of", /5X, "the parameter LFLNO at the start of this routine (in file ", "randgen.f)", /5X, "and recompile. Any number less than 100 should work.");
            return;
        }
        goto g10;
g12:
        scanf(FILNO, "(3(I2),I3)") SS, MM, HH, DD;
        CLOSE(FILNO);
        SYSTEM("rm zbql1234.tmp");
        DSS = DINT((DBLE(SS)/6.0D1)*zbql0001.B);
        DMM = DINT((DBLE(MM)/6.0D1)*zbql0001.B);
        DHH = DINT((DBLE(HH)/2.4D1)*zbql0001.B);
        DDD = DINT((DBLE(DD)/3.65D2)*zbql0001.B);
        TMPVAR1 = DMOD(DSS + DMM + DHH + DDD, zbql0001.B);
        //*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        //*<<<<<<<<	... TO HERE (END OF COMMENTING OUT FOR 	  <<<<<<<
        //*<<<<<<<<	USERS WITHOUT 'CALL SYSTEM' CAPABILITY	  <<<<<<<
        //*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    } else
    {
        TMPVAR1 = DMOD(DBLE(SEED), zbql0001.B);
    }
    zbql0001.ZBQLIX(1) = TMPVAR1;
    for(I=2; I<=43; I++)
    {
        TMPVAR1 = zbql0001.ZBQLIX(I - 1)*3.0269D4;
        TMPVAR1 = DMOD(TMPVAR1, zbql0001.B);
        zbql0001.ZBQLIX(I) = TMPVAR1;
        
    }
    
    
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLU01(DUMMY)
{
    //*
    //*       Returns a uniform random number between 0 & 1, using
    //*       a Marsaglia-Zaman type subtract-with-borrow generator.
    //*       Uses double precision, rather than integer, arithmetic 
    //*       throughout because MZ's integer constants overflow
    //*       32-bit integer storage (which goes from -2^31 to 2^31).
    //*       Ideally, we would explicitly truncate all integer 
    //*       quantities at each stage to ensure that the double
    //*       precision representations do not accumulate approximation
    //*       error; however, on some machines the use of DNINT to
    //*       accomplish this is *seriously* slow (run-time increased
    //*       by a factor of about 3). This double precision version 
    //*       has been tested against an integer implementation that
    //*       uses long integers (non-standard and, again, slow) -
    //*       the output was identical up to the 16th decimal place
    //*       after 10^10 calls, so we're probably OK ...
    //*
    double ZBQLU01, DUMMY, zbql0001.B, zbql0001.C, zbql0001.ZBQLIX(43), X, B2, BINV;
    int CURPOS, ID22, ID43;
    
    //COMMON/ZBQL0001/ZBQLIX, B, C
    SAVE/ZBQL0001/;
    SAVE CURPOS, ID22, ID43;
    DATA CURPOS, ID22, ID43/1, 22, 43/;
    
    B2 = zbql0001.B;
    BINV = 1.0D0/zbql0001.B;
g5:
    X = zbql0001.ZBQLIX(ID22) - zbql0001.ZBQLIX(ID43) - zbql0001.C;
    if (X < 0.0D0) 
    {
        X = X + zbql0001.B;
        zbql0001.C = 1.0D0;
    } else
    {
        zbql0001.C = 0.0D0;
    }
    zbql0001.ZBQLIX(ID43) = X;
    //*
    //*     Update array pointers. Do explicit check for bounds of each to
    //*     avoid expense of modular arithmetic. If one of them is 0 the others
    //*     won't be
    //*
    CURPOS = CURPOS - 1;
    ID22 = ID22 - 1;
    ID43 = ID43 - 1;
    if (CURPOS == 0) 
    {
        CURPOS = 43;
    } else
    {if (ID22 == 0) THEN
        ID22 = 43;
    } else
    {if (ID43 == 0) THEN
        ID43 = 43;
    }
    //*
    //*     The integer arithmetic there can yield X=0, which can cause 
    //*     problems in subsequent routines (e.g. ZBQLEXP). The problem
    //*     is simply that X is discrete whereas U is supposed to 
    //*     be continuous - hence if X is 0, go back and generate another
    //*     X and return X/B^2 (etc.), which will be uniform on (0,1/B). 
    //*
    if (X < BINV) 
    {
        B2 = B2*zbql0001.B;
        goto g5;
    }
    
    ZBQLU01 = X/B2;
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLUAB(A, B)
{
    //*
    //*       Returns a random number uniformly distributed on (A,B)
    //*
    double A, B, ZBQLU01, ZBQLUAB;
    
    //*
    //*       Even if A > B, this will work as B-A will then be -ve
    //*
    if (A != B) 
    {
        ZBQLUAB = A + ((B - A)*ZBQLU01(0.0D0));
    } else
    {
        ZBQLUAB = A;
        printf(*, 1)  // format(/5X, "****WARNING**** (function ZBQLUAB) Upper and ", "lower limits on uniform", /5X, "distribution are identical", /);
    }
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLexp(MU)
{
    //*
    //*       Returns a random number exponentially distributed with
    //*       mean MU
    //*
    double MU, ZBQLEXP, ZBQLU01;
    
    ZBQLEXP = 0.0D0;
    
    if (MU < 0.0D0) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLEXP", /);
        return;
    }
    
    ZBQLEXP = -DLOG(ZBQLU01(0.0D0))*MU;
    
    
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLNOR(MU, SIGMA)
{
    //*
    //*       Returns a random number Normally distributed with mean
    //*       MU and standard deviation |SIGMA|, using the Box-Muller
    //*       algorithm
    //*
    double THETA, R, ZBQLNOR, ZBQLU01, PI, MU, SIGMA;
    double SPARE;
    int STATUS;
    SAVE STATUS, SPARE, PI;
    DATA STATUS/ - 1/;
    
    if (STATUS == -1) PI = 4.0D0*DATAN(1.0D0);
    
    if (STATUS <= 0) 
    {
        THETA = 2.0D0*PI*ZBQLU01(0.0D0);
        R = Dsqrt(-2.0D0*DLOG(ZBQLU01(0.0D0)));
        ZBQLNOR = (R*Dcos(THETA));
        SPARE = (R*Dsin(THETA));
        STATUS = 1;
    } else
    {
        ZBQLNOR = SPARE;
        STATUS = 0;
    }
    
    ZBQLNOR = MU + (SIGMA*ZBQLNOR);
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLBIN(N, P)
{
    //*
    //*       Returns a random number binomially distributed (N,P)
    //*
    double P, ZBQLBET1;
    double PP, PPP, G, Y, TINY;
    int N, ZBQLBIN, ZBQLGEO, IZ, NN;
    
    TINY = 1.0D - 8;
    ZBQLBIN = 0;
    
    if (.NOT.((P >= 0.0D0)) && ((P <= 1.0D0))) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLBIN", /);
        return;
    } else
    {if (N <= 0) THEN
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLBIN", /);
        return;
    }
    //*
    //*	First step: if NP > 10, say, things will be expensive, and 
    //*	we can get into the right ballpark by guessing a value for
    //*	ZBQLBIN (IZ, say), and simulating Y from a Beta distribution 
    //*	with parameters IZ and NN-IZ+1 (NN starts off equal to N).
    //*	If Y is less than PP (which starts off as P) then the IZth order 
    //*	statistic from NN U(0,1) variates is less than P, and we know 
    //*	that there are at least IZ successes. In this case we focus on 
    //*	the remaining (NN-IZ) order statistics and count how many are
    //*	less than PP, which is binomial (NN-IZ,(PP-Y)/(1-Y)). 
    //*	Otherwise, if Y is greater than PP there must be less 
    //*	than IZ successes, so we can count the number of order statistics
    //*	under PP, which is binomial (IZ-1,P/Y). When we've got NN*PP
    //*	small enough, we go to the next stage of the algorithm and 
    //*	generate the final bits directly.
    //*
    NN = N;
    PP = P;
g10:
    IZ = INT(DBLE(NN)*PP) + 1;
    if ((IZ > 10)) && ((IZ < NN - 10)) 
    {
        Y = ZBQLBET1(DBLE(IZ), DBLE(NN - IZ + 1));
        if (Y < PP) 
        {
            ZBQLBIN = ZBQLBIN + IZ;
            NN = NN - IZ;
            PP = (PP - Y)/(1.0D0 - Y);
        } else
        {
            NN = IZ - 1;
            PP = PP/Y;
        }
        goto g10;
    }
    //*
    //*	PP is the probability of the binomial we're currently
    //*	simulating from. For the final part, we simulate either number 
    //*	of failures or number of success, depending which is cheaper.
    //*      
g20:
    if (PP > 0.5) 
    {
        PPP = 1.0D0 - PP;
    } else
    {
        PPP = PP;
    }
    
    G = 0;
    IZ = 0;
    //*
    //*     ZBQLGEO falls over for miniscule values of PPP, so ignore these
    //*     (tiny probability of any successes in this case, anyway)
    //* 
    if (PPP > TINY) 
    {
g30:
        G = G + ZBQLGEO(PPP);
        if (G <= NN) 
        {
            IZ = IZ + 1;
            goto g30;
        }
    }
    
    if (PP > 0.5) IZ = NN - IZ;
    ZBQLBIN = ZBQLBIN + IZ;
    
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLGEO(P)
{
    //*
    //*       Returns a random number geometrically distributed with 
    //*       parameter P ie. mean 1/P
    //* 
    
    double P, ZBQLU01, U, TINY;
    int ZBQLGEO;
    
    TINY = 1.0D - 12;
    ZBQLGEO = 0;
    
    if (.NOT.((P >= 0.0D0)) && ((P <= 1.0D0))) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLGEO", /);
        return;
    }
    
    if (P > 0.9D0) 
    {
g10:
        ZBQLGEO = ZBQLGEO + 1;
        U = ZBQLU01(0.0D0);
        if (U > P) goto g10;
    } else
    {
        U = ZBQLU01(0.0D0);
        //*
        //*	For tiny P, 1-p will be stored inaccurately and log(1-p) may
        //*	be zero. In this case approximate log(1-p) by -p
        //*
        if (P > TINY) 
        {
            ZBQLGEO = 1 + INT(DLOG(U)/DLOG(1.0D0 - P));
        } else
        {
            ZBQLGEO = 1 + INT(-DLOG(U)/P);
        }
    }
    
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLPOI(MU)
{
    //*
    //*       Returns a random number Poisson distributed with mean MU
    //*
    
    double ZBQLU01, X, Y, MU, PI;
    double ZBQLLG, ZBQLGAM, MU1, TMP1, TMP2, T;
    int ZBQLPOI, ZBQLBIN, K, INIT;
    SAVE INIT, PI;
    DATA INIT/0/;
    
    if (INIT == 0) 
    {
        PI = 4.0D0*DATAN(1.0D0);
        INIT = 1;
    }
    
    ZBQLPOI = 0;
    
    if (MU < 0.0D0) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLPOI", /);
        return;
    }
    //*
    //*      For small MU, generate exponentials till their sum exceeds 1
    //*      (equivalently, uniforms till their product falls below e^-MU)
    //*
    if (MU <= 1.0D3) 
    {
        MU1 = MU;
        //*
        //*     For values of MU less than 1000, use order statistics - the Kth
        //*     event in a Poisson process of rate MU has a Gamma distribution
        //*     with parameters (MU,K); if it's greater than 1 we know that there 
        //*     are less than K events in (0,1) (and the exact number is binomial)
        //*     and otherwise the remaining number is another Poisson. Choose K so
        //*     that we'll get pretty close to 1 in the first go but are unlikely
        //*     to overshoot it.
        //*
g19:
        if (MU1 > 1.0D1) 
        {
            K = INT(MU1 - Dsqrt(MU1));
            Y = ZBQLGAM(DBLE(K), MU1);
            if (Y > 1.0D0) 
            {
                ZBQLPOI = ZBQLPOI + ZBQLBIN(K - 1, (1.0D0/Y));
                return;
            }
            ZBQLPOI = ZBQLPOI + K;
            MU1 = MU1*(1.0D0 - Y);
            goto g19;
        }
        Y = Dexp(-MU1);
        X = 1.0D0;
g20:
        X = X*ZBQLU01(0.0D0);
        if (X > Y) 
        {
            ZBQLPOI = ZBQLPOI + 1;
            goto g20;
        }
        //*
        //*     For really huge values of MU, use rejection sampling as in 
        //*     Press et al (1992) - large numbers mean some accuracy may be
        //*     lost, but it caps the execution time.
        //*
    } else
    {
        TMP1 = Dsqrt(2.0D0*MU);
        TMP2 = ZBQLLG(MU + 1.0D0) - (MU*DLOG(MU));
g30:
        Y = DTAN(PI*ZBQLU01(0.0D0));
        ZBQLPOI = INT(MU + (TMP1*Y));
        if (ZBQLPOI < 0) goto g30;
        X = DBLE(ZBQLPOI);
        T = (X*DLOG(MU) - ZBQLLG(X + 1.0D0)) + TMP2;
        if (DABS(T) < 1.0D2) 
        {
            T = 0.9D0*(1.0D0 + (Y*Y))*Dexp(T);
            if (ZBQLU01(0.0D0) > T) goto g30;
        } else
        {
            T = DLOG(0.9D0*(1.0D0 + (Y*Y))) + T;
            if (DLOG(ZBQLU01(0.0D0)) > T) goto g30;
        }
    }
    
    
};


//******************************************************************
// --------------------------------------------
unknown_retval ZBQLGAM(G, H)
{
    //*
    //*       Returns a random number with a gamma distribution with mean
    //*       G/H and variance G/(H^2). (ie. shape parameter G & scale
    //*       parameter H)
    //*
    double C, D, R, ZBQLGAM, ZBQLU01, G, H, A, z1, z2, B1, B2, M;
    double U1, U2, U, V, TEST, X;
    double c1, c2, c3, c4, c5, w;
    
    ZBQLGAM = 0.0D0;
    
    if ((G <= 0.0D0)) || ((H < 0.0D0)) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLGAM", /5X, "(both parameters must be positive)", /);
        return;
    }
    
    if (G < 1.0D0) 
    {
g889:
        u = ZBQLU01(0.0d0);
        v = ZBQLU01(0.0d0);
        if (u > exp(1.0d0)/(g + exp(1.0d0))) goto g891;
        ZBQLGAM = ((g + exp(1.0d0))*u/exp(1.0d0))**(1.0d0/g);
        if (v > exp(-ZBQLGAM)) 
        {
            goto g889;
        } else
        {
            goto g892;
        }
g891:
        ZBQLGAM = -log((g + exp(1.0d0))*(1.0d0 - u)/(g*exp(1.0d0)));
        if (v > ZBQLGAM**(g - 1.0)) goto g889;
g892:
        ZBQLGAM = ZBQLGAM/h;
        return;
    } else
    {if (G < 2.0D0) THEN
        M = 0.0D0;
    } else
    {if (g > 10.0d0) then
        c1 = g - 1.0d0;
        c2 = (g - 1.0d0/(6.0d0*g))/c1;
        c3 = 2.0d0/c1;
        c4 = c3 + 2.0d0;
        c5 = 1.0d0/sqrt(g);
g777:
        u = ZBQLU01(0.0d0);
        v = ZBQLU01(0.0d0);
        if (g > 2.50d0) 
        {
            u = v + c5*(1.0d0 - 1.860d0*u);
        }
        if (u <= 0.0d0) || (u >= 1.0d0) goto g777;
        w = c2*v/u;
        if (c3*u + w + 1.0d0/w <= c4) goto g778;
        if (c3*log(u) - log(w) + w >= 1.0d0) goto g777;
g778:
        ZBQLGAM = c1*w/h;
        return;
    } else
    {
        M = -(G - 2.0D0);
    }
    R = 0.50D0;
    a = ((g - 1.0d0)/exp(1.0d0))**((g - 1.0d0)/(r + 1.0d0));
    C = (R*(M + G) + 1.0D0)/(2.0D0*R);
    D = M*(R + 1.0D0)/R;
    z1 = C - Dsqrt(C*C - D);
    //*
    //*     On some systems (e.g. g77 0.5.24 on Linux 2.4.24), C-DSQRT(C*C)
    //*     is not exactly zero - this needs trapping if negative.
    //*
    if ((Z1 - M < 0.0D0)) && ((Z1 - M > -1.0D - 12)) Z1 = M;
    z2 = C + Dsqrt(C*C - D);
    B1 = (z1*(z1 - M)**(R*(G - 1.0D0)/(R + 1.0D0)))*Dexp(-R*(z1 - M)/(R + 1.0D0));
    B2 = (z2*(z2 - M)**(R*(G - 1.0D0)/(R + 1.0D0)))*Dexp(-R*(z2 - M)/(R + 1.0D0));
g50:
    U1 = ZBQLU01(0.0D0);
    U2 = ZBQLU01(0.0D0);
    U = A*U1;
    V = B1 + (B2 - B1)*U2;
    X = V/(U**R);
    if (X <= M) goto g50;
    TEST = ((X - M)**((G - 1)/(R + 1)))*exp(-(X - M)/(R + 1.0D0));
    if (U <= TEST) 
    {
        ZBQLGAM = (X - M)/H;
    } else
    {
        goto g50;
    }
    
    
};


//***************************************************************
// --------------------------------------------
unknown_retval ZBQLBET1(NU1, NU2)
{
    //*
    //*       Returns a random number, beta distributed with degrees
    //*       of freedom NU1 and NU2.
    //*
    double NU1, NU2, ZBQLGAM, ZBQLBET1, ZBQLU01, X1, X2;
    
    ZBQLBET1 = 0.0D0;
    
    if ((NU1 <= 0.0)) || ((NU2 <= 0.0)) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLBET1", /5X, "(both degrees of freedom must be positive)", /);
        return;
    }
    //*      
    //*       If parameters are too small, gamma subroutine tends to return zero
    //*       as all the probability goes to the origin and we get rounding
    //*       errors, even with double precision. In this case, we use Johnk's
    //*       method, suitably scaled to avoid rounding errors as much as possible.
    //*
    
    if ((NU1 < 0.9D0)) && ((NU2 < 0.9D0)) 
    {
g10:
        X1 = ZBQLU01(0.0D0);
        X2 = ZBQLU01(0.0D0);
        if ((X1**(1.0D0/NU1)) + (X2**(1.0D0/NU2)) > 1.0D0) goto g10;
        X1 = (DLOG(X2)/NU2) - (DLOG(X1)/NU1);
        ZBQLBET1 = (1.0D0 + Dexp(X1))**(-1);
        if (ZBQLBET1 > 1.0D0) goto g10;
    } else
    {
        X1 = ZBQLGAM(NU1, 1.0D0);
        X2 = ZBQLGAM(NU2, 1.0D0);
        ZBQLBET1 = X1/(X1 + X2);
    }
     ;
    return;
    
    
    
};


//***************************************************************
// --------------------------------------------
unknown_retval ZBQLWEI(A, B)
{
    //*
    //*       Returns a random number, Weibull distributed with shape parameter
    //*       A and location parameter B, i.e. density is
    //*	f(x) = ( A/(B**A) ) * x**(A-1) * EXP( -(x/B)**A )
    //*
    double A, B, ZBQLU01, ZBQLWEI, U;
    
    ZBQLWEI = 0.0D0;
    
    if ((A <= 0.0)) || ((B <= 0.0)) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLWEI", /5X, "(both parameters must be positive)", /);
        return;
    }
    
    U = ZBQLU01(0.0D0);
    ZBQLWEI = B*((-DLOG(U))**(1.0D0/A));
    
    
};


//***************************************************************
// --------------------------------------------
unknown_retval ZBQLNB(R, P)
{
    //*
    //*       Returns a pseudo-random number according to a Negative
    //*	Binomial distribution with parameters (R,P). NB these are
    //*	both DOUBLE - it copes with non-integer R as well. The
    //*       form of the distribution is *not* the no. of trials to 
    //*       the Rth success - see documentation for full spec.
    //*
    double R, P, ZBQLGAM, Y;
    int ZBQLNB, ZBQLPOI;
    
    ZBQLNB = 0;
    
    if ((R <= 0.0D0)) || ((P <= 0.0D0)) || ((P >= 1.0D0)) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLNB");
        return;
    }
    
    Y = ZBQLGAM(R, 1.0D0);
    Y = Y*P/(1.0D0 - P);
    ZBQLNB = ZBQLPOI(Y);
    
    
};


//***************************************************************
// --------------------------------------------
unknown_retval ZBQLPAR(A, B)
{
    //*
    //*     Returns a random number, Pareto distributed with parameters
    //*     A and B. The density is A*(B**A) / (B+X)**(A+1) for X > 0.
    //*     (this is slightly nonstandard - see documentation in 
    //*     randgen.txt). The algorithm is straightforward - it uses the
    //*     inverse CDF method.
    //*
    double A, B, ZBQLPAR, ZBQLU01, U;
    
    ZBQLPAR = 0.0D0;
    
    if ((A <= 0.0D0)) || ((B <= 0.0D0)) 
    {
        printf(*, 1)  // format(/5X, "****ERROR**** Illegal parameter value in ", " ZBQLPAR", /5X, "(both parameters must be positive)", /);
        return;
    }
    
    U = ZBQLU01(0.0D0);
    ZBQLPAR = B*(U**(-1.0D0/A) - 1.0D0);
    
    
};


//***************************************************************
// --------------------------------------------
unknown_retval ZBQLLG(X)
{
    //*
    //*     Returns log(G(X)) where G is the Gamma function. The algorithm is
    //*     that given in Press et al (1992), Section 6.1, although this
    //*     version also allows for arguments less than 1.
    //*
    double X, Z, Z2, ZBQLLG, PI, RLN2P, C(0:6), TMP, SUM;
    int INIT, I;
    SAVE INIT, C, RLN2P, PI;
    DATA INIT/0/;
    DATA(C(I), I = 0, 6)/1.000000000190015D0, 76.18009172947146D0,       -86.50532032941677D0, 24.01409824083091D0,             -1.231739572450155D0, 0.1208650973866179D - 2,             -0.5395239384953D - 5/;
    
    if (INIT == 0) 
    {
        PI = 4.0D0*DATAN(1.0D0);
        RLN2P = 0.5D0*DLOG(2.0D0*PI);
        INIT = 1;
    }
    //*
    //*     Compute for x > 1, then use transformation if necessary. Z is
    //*     our working argument.
    //*
    if (X >= 1.0D0) 
    {
        Z = X;
    } else
    {
        Z = 2.0D0 - X;
        Z2 = 1.0D0 - X;
    }
    
    if (DABS(Z - 1.0D0) < 1.0D - 12) 
    {
        ZBQLLG = 0.0D0;
        return;
    }
    
    TMP = Z + 4.5D0;
    TMP = ((Z - 0.5D0)*DLOG(TMP)) - TMP + RLN2P;
    
    SUM = C(0);
    for(I=1; I<=6; I++)
    {
        SUM = SUM + (C(I)/(Z + DBLE(I - 1)));
        
    }
    ZBQLLG = TMP + DLOG(SUM);
    //*
    //*     Transformation required if X<1
    //*
    if (X < 1.0D0) 
    {
        TMP = PI*Z2;
        ZBQLLG = DLOG(TMP/Dsin(TMP)) - ZBQLLG;
    }
    
};
