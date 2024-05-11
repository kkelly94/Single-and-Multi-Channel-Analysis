#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>

using namespace std;

int main()
{
	fstream textfile;

	double input[13];
	int xx = 1;
	textfile.open("input.txt");
	while (!textfile.eof())
	{
		textfile >> input[xx];
		//cout << input[x] << endl;
		if (xx < 13)
		{
			xx++;
		}
		else
		{
			break;
		}
	}
	textfile.close();

	const double PI = 3.1415926535897;

	double L = input[1];
	double Le = input[2];
	double D = input[3];
	double Pitch = input[4];
	double Tin = input[5];
	double mdot = input[6];
	double Pnominal = input[7];
	double qprimemax = input[8];
	double PWR = input[9]; //0 pwr 1 bwr
	double Dci = input[10];
	double Dfo = input[11];
	double SC = input[12]; //1 on 0 off
	//Added properties for htc calc.
	double pd=Pitch/D;
	double D_h=D*((4/PI)*pow(pd, 2)-1);
	//cout<<"DH "<<D_h<<endl;
	double A=pow(Pitch,2)-(PI/4*pow(D,2));
	double G=mdot/A;
	//cout<<G<<endl;
	/*
	cout << L << endl;
	cout << Le << endl;
	cout << D << endl;
	cout << Pitch << endl;
	cout << Tin << endl;
	cout << mdot << endl;
	cout << Pnominal << endl;
	cout << qprimemax << endl;
	cout << PWR << endl;
	cout << Dci << endl;
	cout << Dfo << endl;
	cout << SC << endl;
	cout << L << endl;
	cout << Pitch << endl;
	*/

	//Creates array for property table ****************************************************************************************
	double prop[11][373];
	int i = 0;
	int j = 0;
	textfile.open("proptable.txt");
	while (!textfile.eof() && j < 373)
	{
		if (i < 10)
		{
			textfile >> prop[i][j];
			//cout << prop[i][j] << endl;
			i++;
		}
		else if ((i = 10) && (j < 373))
		{
			i = 0;
			j++;
		}
		else
		{
			cout << "how?" << endl;
			break;
		}
	}
	textfile.close();

    ofstream out;
    out.open("output.txt");
	//Linear interpolation to set properties for a given temperature **********************************************************************
	double T_mn=Tin;
	double T_mn2=Tin;
	int m = 0;
	int n=0;
    int m2 = 0;
	int n2=0;
    int m3 = 0;
	int n3=0;
	int cell=1;
	int xtrack=0;
	double x_e=0;

	//Tsat properties

		int tempfinder3=0;
	double tempfinder4=0;
	while (Pnominal > tempfinder4)
    {
        tempfinder3++;
        tempfinder4= prop[1][tempfinder3];
    }
    //cout<<"Psat " << prop[1][tempfinder3];

	double T_sat = (Pnominal - prop[1][tempfinder3-1])/(prop[1][tempfinder3]-prop[1][tempfinder3-1])+prop[0][tempfinder3-1];//Interpolates temperature value
	//cout<<"Tsat " << T_sat<<endl;
	while (T_sat > m3)
	{
		m3++;
	}

	n3 = m3 - 2;//This is the order! n first, m second. Otherwise, m will be redifined and give a false n value
	m3 = m3 - 1;

		double PX2 = prop[2][m3];
	double QX2 = prop[2][n3];
	double vol_fsat = ((T_sat - m3)*(PX2 - QX2) + QX2);

	double PX3 = prop[3][m3];
	double QX3 = prop[3][n3];
	double vol_gsat = ((T_sat - m3)*(PX3 - QX3) + QX3);
	//cout << "vol_g " << vol_g << endl;

	double PX4 = prop[4][m3];
	double QX4 = prop[4][n3];
	double h_fsat = ((T_sat - m3)*(PX4 - QX4)+ QX4);
	//cout << setprecision (9) <<"h_fsat "<< h_fsat << endl;

	double PX5 = prop[5][m3];
	double QX5 = prop[5][n3];
	double h_gsat = ((T_sat - m3)*(PX5 - QX5) + QX5);
	//cout << "h_gsat " <<h_gsat << endl;

	double PX6 = prop[6][m3];
	double QX6 = prop[6][n3];
	double mu_fsat = ((T_sat - m3)*(PX6 - QX6)+ QX6);
	//cout << "mu_f " << mu_f << endl;

	double PX7 = prop[7][m3];
	double QX7 = prop[7][n3];
	double k_fsat = ((T_sat - m3)*(PX7 - QX7) + QX7);
	//cout << "k_f " << k_f << endl;

	double PX8 = prop[8][m3];
	double QX8 = prop[8][n3];
	double Pr_fsat = ((T_sat - m3)*(PX8 - QX8)+ QX8);
	//cout <<"Pr_f " << Pr_f << endl;

	double PX9 = prop[9][m3];
	double QX9 = prop[9][n3];
	double mu_gsat = ((T_sat - m3)*(PX9 - QX9) + QX9);

	double z = -L/2+L/400;
	while (z< (-L/2+400*L/400))
   {
    out<< setprecision (9) << "Cell# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<cell<<endl;
    cell++;
   //Enthalpy at z+.5 for calculating T at z+1
	while (T_mn > m)
	{
		m++;
	}

	n = m - 2;//This is the order! n first, m second. Otherwise, m will be redifined and give a false n value
	m = m - 1;//This is the order! n first, m second. Otherwise, m will be redifined and give a false n value

	out<<"z " << z << endl;

	double P4 = prop[4][m];
	double Q4 = prop[4][n];
	double h_f = ((T_mn - m)*(P4 - Q4)+ Q4);



    //double bear = cos(PI*(z+(L/800)));
	double qprimez=cos(PI*(z-(L/800))/Le)*qprimemax; //The idea is to have z be channel location that we are calculating properties for, so z-L/800 gives the .5 fraction before that
	double qprimez2=cos(PI*(z)/Le)*qprimemax;
	//cout << "COSINE VALUE IS: " <<bear << endl;
	//cout<<"q'(z) " <<qprimez << endl;
	double h_fn = h_f+(qprimez*(L/400))/mdot;
	//double h_fn2 = h_f+(qprimez2*L/400)/mdot;
	//double cat=qprimez*(L/400)/mdot;

			int tempfinder=0;
	double tempfinder2=0;
	while (h_fn > tempfinder2)
    {
        tempfinder++;
        tempfinder2= prop[4][tempfinder];
    }
    //cout<<prop[0][tempfinder]<< endl;
    //cout << prop[0][tempfinder-1]<<endl;
	T_mn = (h_fn-prop[4][tempfinder-1])/(prop[4][tempfinder]-prop[4][tempfinder-1])*(prop[0][tempfinder]-prop[0][tempfinder-1])+prop[0][tempfinder-1];  //Interpolates temperature value

	//Outlet Properties
    T_mn2=T_mn;

	while (T_mn2 > m2)
	{
		m2++;
	}

	n2 = m2 - 2;//This is the order! n first, m second. Otherwise, m will be redifined and give a false n value
	m2 = m2 - 1;//This is the order! n first, m second. Otherwise, m will be redifined and give a false n value

    //cout<<"T_mn2 " <<T_mn2 <<endl;

	double PP2 = prop[2][m2];
	double QQ2 = prop[2][n2];
	double vol_f2 = ((T_mn2 - m2)*(PP2 - QQ2) + QQ2);
	//cout << "vol_f2 " << vol_f2 << endl;

	double PP3 = prop[3][m2];
	double QQ3 = prop[3][n2];
	double vol_g2 = ((T_mn2 - m2)*(PP3 - QQ3) + QQ3);
	//cout << "vol_g " << vol_g << endl;

    double PP6 = prop[6][m2];
	double QQ6 = prop[6][n2];
	double mu_f2 = ((T_mn2 - m2)*(PP6 - QQ6)+ QQ6);
	//cout << "mu_f2 " << mu_f2 << endl;
    double PP4 = prop[4][m2];
	double QQ4 = prop[4][n2];
	double h_f2 = ((T_mn2 - m2)*(PP4 - QQ4)+ QQ4);
	//cout <<"h_f2 "<< h_f2 << endl;
    double PP5 = prop[5][m2];
	double QQ5 = prop[5][n2];
	double h_g2 = ((T_mn2 - m2)*(PP5 - QQ5) + QQ5);

    double PP7 = prop[7][m2];
	double QQ7 = prop[7][n2];
	double k_f2 = ((T_mn2 - m2)*(PP7 - QQ7) + QQ7);

    double PP8 = prop[8][m2];
	double QQ8 = prop[8][n2];
	double Pr_f2 = ((T_mn2 - m2)*(PP8 - QQ8)+ QQ8);

	double PP9 = prop[9][m2];
	double QQ9 = prop[9][n2];
	double mu_g2 = ((T_mn2 - m2)*(PP9 - QQ9) + QQ9);


	//cout<< cat<< endl;
	//out<<"h_f2 "<<h_f2<< endl;
	double Re=G*D_h/mu_f2;
	double Nu=(.023*(pow(Re, .8))*(pow (Pr_f2, .333))*(1.826*Pitch/D-1.0430));   //Im not sure about this equation, dunno where the "(1.826*Pitch/D-1.0430)" comes from
    double htc=Nu*k_f2/D_h;
    //cout<<"Re " << Re << endl;
    //cout << "htc " << htc << endl;





    //cout<<prop[0][tempfinder]<< endl;
    //cout << prop[0][tempfinder-1]<<endl;



    //cout<<prop[0][tempfinder]<< endl;
    //cout << prop[0][tempfinder-1]<<endl;


    out<< setprecision (9) <<"T_mn "<<T_mn<< endl;

    double k_c=15;
	double R_co=D/2;
	double R_ci=Dci/2;
	double R_fo=Dfo/2;
	double T_co=qprimez2/(2*PI*R_co*htc)+T_mn2; //subbing in schubrings values
	double T_ci=qprimez2/(2*PI*k_c)*log(R_co/R_ci)+T_co;
    out<<"T_co "<<T_co<< endl;
	out<<"T_ci "<<T_ci<< endl;
	//double moose =log(R_co/R_ci);
//COPIED CODE STARTS HERE
	double d_eff = R_ci - R_fo;
	double k_gas = 0;
	double htc_g = 0;
	double T_fo = 0;
	double T_fo_k = T_fo + 273.15;
	double T_fo_guess = 300;
	double T_fo_guess_k = T_fo_guess+ 273.15;
	double T_ci_k = T_ci + 273.15;
	while (abs(T_fo_k - T_fo_guess_k) > 0.000001)
	{
		//T_fo_guess_k = T_fo_guess_k + 0.01;
		k_gas = (15.8*0.0001)*pow((((T_fo_guess_k) + (T_ci_k)) / 2), 0.79);
		//cout << k_gas << endl;
		htc_g = (k_gas / d_eff) + (5.67E-8)*((pow(T_fo_guess_k, 4) - pow(T_ci_k, 4)) / ((T_fo_guess_k) - (T_ci_k)));
		//cout << htc_g << endl;
		T_fo_k = qprimez2 / (PI*(R_ci + R_fo)*htc_g) + T_ci_k;
		//cout << T_fo_k << endl;

		if (T_fo_k > T_fo_guess_k)
		{
			T_fo_guess_k = T_fo_guess_k + (T_fo_k - T_fo_guess_k) / 2;
			//cout << T_fo_guess_k << endl;
		}
		else if (T_fo_k < T_fo_guess_k)
		{
			T_fo_guess_k = T_fo_guess_k - (T_fo_k - T_fo_guess_k) / 2;
			//cout << T_fo_guess_k << endl;
		}
	}
	//cout << T_fo_guess_k << endl;
	T_fo = T_fo_guess_k - 273.15;
	out <<"T_fo is: "<< T_fo << endl;

	//8.Solving for T_max Good
	double T_max = 0; // needs to be solved iteratively
	double M1 = 3824 * log(402.4 + T_fo) + (6.1256E-11) / 4 * pow(T_fo + 273, 4) ;
	double T_max_guess = 300;

	double qprimez_guess = 0;
	while (abs(qprimez2 - qprimez_guess)>0.01)
	{
		double M2 = 3824 * log(402.4 + T_max_guess) + (6.1256E-11) / 4 * pow(T_max_guess + 273, 4);
		qprimez_guess = 4 * PI*(M2 - M1);
		//cout << qprimez_guess << endl;
		T_max_guess = T_max_guess + 0.0001;

	}
	out <<"T_max is: " <<T_max_guess << endl;




//COPIED CODE ENDS*/

	//cout << moose <<endl;
	//double T_fo2=qprimez/(PI*(R_ci+R_fo)*htc_g)+289.73536;//+T_ci;
	//cout<<"T_fo "<<T_fo2<< endl;
	double dPfric=pow((D_h*G/mu_f2),-.18)*(.1339+.09059*((Pitch/D)-1)-.09926*pow(((Pitch/D)-1),2))*pow(G,2)*vol_f2/(D_h*2)*(L/400);
	//cout<<"dPfric " << dPfric << endl;
/*		double dPaccel=pow((G,2)*(vol_f-vol_f2);
	cout<<"dPaccel " << dPaccel << endl; NYI*/
		double dPgrav=9.81/vol_f2*(L/400);
	//cout<<"dPgrav " << dPgrav << endl;
	double dP=dPgrav+dPfric;
	out << "dP " << dP << endl;

    double qprimed = qprimez2/(2*PI*R_co);
		x_e = (h_f2-h_fsat)/(h_gsat-h_fsat);
		out<<"xe right after calc " <<x_e<<endl;

    if ((x_e<0 && SC==0)|| x_e>0)
    {
        T_mn2=T_sat;
        vol_f2=vol_fsat;
        vol_g2=vol_gsat;
        mu_g2=mu_gsat;
        mu_f2=mu_fsat;
        h_f2=h_fsat;
        h_g2=h_gsat;
        Pr_f2=Pr_fsat;
        k_f2=k_fsat;
    }

		double x;
		double x_in=x;
		x=0;
		double X_tt;
		//double T_mnboil;d
		bool flag = true;
		bool boil =false;
if (x_e < 0 && (SC==0))
{
    x_e = 0;
    x=x_e;
}
else if (x_e>0 && (SC==0))
{
    x=x_e;
}
//assigns first equillibrium quality
    double x_e0;
    if (xtrack <1)
    {
        x_e0=0;
    }
	if (T_co>T_sat)
    {
        xtrack++;
    }

    /*out<<"xtraxk " <<xtrack<<endl;
    out<<T_sat<<endl;
    out<<T_co<<endl;
    out<<SC<<endl;*/
    if (xtrack == 1 && flag == true)
    {
      x_e0 = x_e;
      out<<"x_e0 " <<x_e0<<endl;
      flag = false;
      boil=true;
    }

    if (xtrack<2)
    {
        x_e=0;
    }
    // for subcooled
    if ((SC == 1) && xtrack>1)
    {
    x=x_e-x_e0*exp(x_e/x_e0-1);
    }
    else
    {
    x=0;
    }
//Recalculation of T_co to T_max with 2 phase


if (x>0)
{
    X_tt=pow(((1-x)/x),.9)*pow((vol_f2/vol_g2),.5)*pow((mu_f2/mu_g2),.1);
    //cout<<"X_tt "<<X_tt<<endl;
    double Re4=G*D_h/mu_f2;
	double Nu4=(.023*(pow(Re4, .8))*(pow (Pr_f2, .333))*(1.826*Pitch/D-1.0430));
    double htc_lo=Nu4*k_f2/D_h;
    //cout<<htc<<endl;
    double htc_2phase=htc_lo*(7400*qprimed/(G*(h_g2-h_f2))+1.11*pow(X_tt, -.66));
    //cout<<htc_2phase<<endl;

    T_co=qprimez2/(2*PI*R_co*htc_2phase)+T_mn2;
	T_ci=qprimez2/(2*PI*k_c)*log(R_co/R_ci)+T_co;

//COPIED CODE STARTS HERE

	k_gas = 0;
	htc_g = 0;
	T_fo = 0;
	T_fo_k = T_fo + 273.15;
	T_fo_guess = 300;
	T_fo_guess_k = T_fo_guess+ 273.15;
	T_ci_k = T_ci + 273.15;
	while (abs(T_fo_k - T_fo_guess_k) > 0.0000001)
	{
		//T_fo_guess_k = T_fo_guess_k + 0.01;
		k_gas = (15.8*0.0001)*pow((((T_fo_guess_k) + (T_ci_k)) / 2), 0.79);
		//cout << k_gas << endl;
		htc_g = (k_gas / d_eff) + (5.67E-8)*((pow(T_fo_guess_k, 4) - pow(T_ci_k, 4)) / ((T_fo_guess_k) - (T_ci_k)));
		//cout << htc_g << endl;
		T_fo_k = qprimez2 / (PI*(R_ci + R_fo)*htc_g) + T_ci_k;
		//cout << T_fo_k << endl;

		if (T_fo_k > T_fo_guess_k)
		{
			T_fo_guess_k = T_fo_guess_k + (T_fo_k - T_fo_guess_k) / 2;
			//cout << T_fo_guess_k << endl;
		}
		else if (T_fo_k < T_fo_guess_k)
		{
			T_fo_guess_k = T_fo_guess_k - (T_fo_k - T_fo_guess_k) / 2;
			//cout << T_fo_guess_k << endl;
		}
	}
	//cout << T_fo_guess_k << endl;
	T_fo = T_fo_guess_k - 273.15;


	//8.Solving for T_max Good
	T_max = 0; // needs to be solved iteratively
	M1 = 3824 * log(402.4 + T_fo) + (6.1256E-11) / 4 * pow(T_fo + 273, 4) ;
	T_max_guess = 300;

	qprimez_guess = 0;
	while (abs(qprimez2 - qprimez_guess)>0.01)
	{
		double M2 = 3824 * log(402.4 + T_max_guess) + (6.1256E-11) / 4 * pow(T_max_guess + 273, 4);
		qprimez_guess = 4 * PI*(M2 - M1);
		//cout << qprimez_guess << endl;
		T_max_guess = T_max_guess + 0.0001;

	}

    }

	//dp/dz HEM for 2phase
	double Re3=G*D_h/mu_f2;
	double rho_m = 1/((x*vol_g2)+(1-x)*(vol_f2));
	double f=pow(Re3,(-0.18))*(0.1339+0.09059*((Pitch/D)-1)-0.09926*(pow(((Pitch/D)-1),2)));
	/*cout<<"ff "<<f<<endl;
	cout<<"G "<<G<<endl;
    cout<<"Relo "<<Re_lo<<endl;
    cout<<"vol_m "<<vol_m<<endl;
	cout<<"D_h "<<D_h<<endl;*/
	double boilfric=(f/D_h*pow(G,2)/(2*rho_m))*L/400;
	double boilaccel=pow(G,2)*(vol_g2-vol_f2)*(x-x_in);
	double boilgrav=9.81*rho_m*L/400;
	/*out<<"f " <<f<<endl;
	out<<"rho_m "<<rho_m<<endl;
	out<<"dPfric " << boilfric<<endl;
    out<<"dPaccel " << boilaccel<<endl;
    out<<"dPgrav " << boilgrav<<endl;*/
    dP=(boilfric+boilaccel+boilgrav);

    //CHF calculations (bowring with weisman Psi), notify user when CHFR drops below 1.3 PWR 1.9 BWR
    double P_r=.145*Pnominal/1e6;
    double F_1 = pow(P_r, -.368)*exp(.648*(1-P_r));
    double F_2 = F_1/(pow(P_r, -.448)*exp(.245*(1-P_r)));
    double F_3 =pow(P_r, .219);

if (P_r<1) //PWR=1 means BWR
{
    F_1 = (1/1.917)*(pow(P_r,18.942)*(exp(20.89*(1-P_r))) + 0.917);
    F_2 = 1.309*F_1/((pow(P_r,1.316)*(exp(2.444*(1-P_r))) + 0.309));
    F_3 = (1/1.667)*(pow(P_r,17.023)*exp(16.658*(1-P_r))+0.667);
}
    double F_4 =F_3*pow(P_r, 1.649);

double n=2-.5*(P_r);
double A=(2.317*((h_g2-h_f2)*D_h*G/4)*F_1)/(1+.0143*F_2*G*pow(D_h,.5));
double B = G*D_h/4;
double C = (.077*F_3*D_h*G)/(1+.347*F_4*pow((G/1356), n));
double qbowring = (A-B*(h_g2-h_f2)*x)/C*(1.826*Pitch/D-1.0430);
double CHFR = qbowring/(qprimed);
/*cout<<" F1 "<<F_1 <<endl;
cout<<"F2  "<<F_2 <<endl;
cout<<"F3  "<< F_3<<endl;
cout<<"F4  "<<F_4 <<endl;*/
/*out<<"A  "<<A <<endl;
out<<"B  "<< B<<endl;
out<<"C  "<<C <<endl;
//cout<<"qbowring  "<<qbowring <<endl;
//cout<<"q''  "<<qprimed <<endl;*/
if(x==0)
{
    CHFR=0;
}
out<<"CHFR is : "<<CHFR<<endl;

if(PWR==0 && CHFR<1.3 && x!=0)
{
    out<<"warning, below CHFR limit"<<endl;
}
if(PWR==1 && CHFR<1.9 && x!=0)
{
    out<<"warning, below CHFR limit"<<endl;
}


if(x<0){
    out<< setprecision (9) <<"T_mn "<<T_mn<< endl;
}
if(x>0)
{
 out<<"T_mn "<<T_mn2<<endl;
}
    out<<"dP boil:" << dP << endl;
    out<<"T_co boil:"<<T_co<< endl;
	out<<"T_ci boil:"<<T_ci<< endl;
    out <<"T_fo boil: "<< T_fo << endl;
    out <<"T_max boil: " <<T_max_guess << endl;
    out<< " x " <<x<<endl;
    out << "x_e "<<x_e<<endl;

z+=L/400;
}

	return 0;
}
