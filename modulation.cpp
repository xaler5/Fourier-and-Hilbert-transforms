#include<iostream>
#include<complex>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#include<hilbert.h>
#define PI 3.14159265358979323846
using namespace std;

//AM-MODULATION (doubled bandwith)
void am_modulation(int N, double *O, double *s, double f, double t, double T)
{
for(int i=0; i<N; i++){O[i]=(s[i]*cos(2*PI*f*t)); t+=T;}
t=0;
}

//SINGLE SIDE-BAND MODULATION
void single_sideband_modulation(int N, double *p, double *o, double *s, double f,double t, double T,double *&a, double *&b, long double *hil_out)
{
 //CALCULATING FIRST DFT IN ORDER TO BE ABLE TO CALCULATE ALSO HILBERT TRANSFORM AND USE IT IN MODULATION
 int const c=N;
 double *m = new double [c]; //Auxiliar variable used in hilbert transform
 double *l = new double [c]; //Auxiliar variable used in hilbert transform
 forwardDFT(s,N,a,b);
 hilbert(m,l,c,hil_out,a,b);
 for(int i=0; i<N; i++){a[i]=0; b[i]=0;} //Reset of a,b vectors
//SSB-MODULATION
 for(int i=0; i<N; i++){o[i]=((s[i]*cos(2*PI*f*t))-(hil_out[i]*sin(2*PI*f*t))); p[i]=((s[i]*sin(2*PI*f*t))+(hil_out[i]*cos(2*PI*f*t))); t+=T;}	
t=0;
}


double reorder(const int c,int N,double *o, double *O,double *k)
{
double *g= new double [c];
for(int i=0; i<c; i++){g[i]=sqrt((o[i]*o[i]) + (O[i]*O[i]))/N;}
for(int i=0; i<5000; i++){k[i]=g[4999-i];}
for(int i=5000; i<c; i++){k[i]=g[14999-i];}
}

void AMdemodulation(double *m,int N,double *o,double F, double t,double T)
{
for(int i=0; i<N; i++){m[i]=((2*(o[i]*cos(2*PI*F*t)))/(1+cos(4*PI*F*t))); t+=T;}
}

void SSBdemodulation(double *m,int N,double *o,double F,double t, double T,long double *hil_out)
{
for(int i=0; i<N; i++){m[i]=(((2*(o[i]*cos(2*PI*F*t)))+(hil_out[i]*sin(4*PI*F*t)))/(1+cos(4*PI*F*t))); t+=T;}
}

int main()
{
double ff=3; //Max frequency present in the signal
double f=200*ff; //Set carrier frequency
int N=10000;    //Set number of samples
double t=0; //Set starting time
double T=0.0001; // Set sampling rate
int const c=N;
double *reOUTPUT1 = new double [c]; //SSB-MODULATED OUTPUT
double *imOUTPUT1 = new double [c]; //SSB-MODULATED OUTPUT
double *OUTPUT2 = new double [c];   //AM-MODULATED OUTPUT
double *signal = new double [c];
//reading the input datas and filling the signal variable
std::ifstream input("PATH TO BE READ");
for(std::string line; getline(input,line);)
{for(int i=0; i<c; i++){input >> signal[i];}}
//for(int i=0; i<c; i++){signal[i]=sin(2*PI*ff*t); t+=T;} //UNCOMMENT TO BUILD YOUR OWN SIGNAL
t=0;
double *a = new double [c];
double *b = new double [c]; 
double *a2 = new double [c];
double *b2 = new double [c];
double *ddsb = new double [c]; //DemodulatedDSB
double *dssb = new double [c]; //DemodulatedDSB
long double *hil_out = new long double [c]; //Output of hilbert transform

am_modulation(N,OUTPUT2,signal,f,t,T);
forwardDFT(OUTPUT2,N,a,b); 

single_sideband_modulation(N,imOUTPUT1,reOUTPUT1,signal,f,t,T,a2,b2,hil_out);
forwardDFTcomplex(imOUTPUT1,reOUTPUT1,N,a2,b2); //FT of the SSB-modulated signal

AMdemodulation(ddsb,N,OUTPUT2,f,t,T);
SSBdemodulation(dssb,N,reOUTPUT1,f,t,T,hil_out);

int *r= new int [c];
for(int i=0; i<c; i++){r[i]=-5000+i;}
double *k= new double [c];
double *K= new double [c];
reorder(c,N,a2,b2,k);
reorder(c,N,a,b,K);

//WRITING THE OUTPUTS
ofstream output;
output.open("PATH TO BE WRITTEN");
output << "SSB-modulated signal"<<";"<<"DFT of SSB signal" << ";" << "AM-modulated signal" << ";" << "DFT of AM signal" << ";" << "Frequencies" << ";" << "Signal" << ";" << "DSBdemodulated"<<";"<<"SSBdemodulated"<<endl;
for(int i=0; i<c; i++)
{
output << reOUTPUT1[i] << ";" << k[i]/*sqrt((a2[i]*a2[i]) + (b2[i]*b2[i]))/N*/ << ";" << OUTPUT2[i] << ";" << K[i]/*sqrt((a[i]*a[i]) + (b[i]*b[i]))/N*/ << ";" << r[i] << ";" << signal[i] << ";" << ddsb[i] << ";" << dssb[i] << endl;
}
output.close();
system("pause");
return 0;
};
