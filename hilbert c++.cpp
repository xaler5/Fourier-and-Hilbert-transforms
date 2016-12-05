#include<iostream>
#include<complex>
#include<fstream>
#include<sstream>
#include<vector>
#include<stdlib.h>
#include<math.h>
#define PI 3.14159265358979323846
using namespace std;

void forwardDFT(const double *s, const int &N, double *&a, double *&b)
{
    for (int k = 0; k < N; ++k) {
        a[k] = b[k] = 0;
        for (int x = 0; x < N; ++x) {
            a[k] += s[x] * cos(2 * PI * k * x / N);
            b[k] -= s[x] * sin(2 * PI * k * x / N);
        }
    }
}

void inverseDFT( const double *a, const double *b, const int &N, long double *&s)
{
    for (int x = 0; x < N; ++x) {
        s[x] = a[0];
        for (int k = 1; k < N ; ++k) {
            s[x] += (a[k] * cos(2* PI * k * x / N) ) - (b[k] * sin(2* PI * k * x / N));
			
        }
    }
}

void hilbert( double *m, double *l, const int &c, long double *&hil_out,const double *a, const double *b)
{
for(int i=1; i<(c/2); i++){m[i]=b[i]*(2);l[i]=-a[i]*2;}
for(int i=(c/2)+1; i<c; i++){m[i]=0;l[i]=0;}
inverseDFT(m,l,c,hil_out);
}

void filter(double *M, double *n, double *D, double *d,double *s,double *F,int N)
{
for(int i=0; i<N; i++)
{
//FIND ALL MAX's and put them into M and n vectors
if((s[i]>s[i+1])&&(s[i]>s[i-1])){if(s[i]>0){M[i]=s[i];} else {n[i]=s[i];}}
if((s[i]==s[i+1])&&(s[i+1]>s[i+2])){if(s[i]>0){M[i]=s[i];} else {n[i]=s[i];}}
//Cut out all values below a certain limit(=what we consider noise by looking at the signal) and Open a window for each peak above the limit and construct the filtered signal F
if(M[i]>0.3){for(int j=0; j<20; j++){F[i-j+10]=s[i-j+10];F[i+j]=s[i+j];}}
}
}

int main () 
{
int s=5000; //NUMBERS OF SAMPLES DATA FROM ECG
const int c=s;
double *sig_in = new double [c]; //Input Signal
long double *sig_out = new long double [c]; //Inverse Fourier transform variable
long double *hil_out = new long double [c]; //Hilber transform variable
double *a = new double [c]; //real part of FT variable
double *b = new double [c]; //imag part of FT variable
double *F = new double [c]; //filtered signal
double *m = new double [c]; //Auxiliar variable used in hilbert transform
double *l = new double [c]; //Auxiliar variable used in hilbert transform
double *M = new double [c]; //max
double *n = new double [c]; //min
double *D = new double [c]; //Derivatives
double *d = new double [c]; //Differentiations
//READING THE INPUT AND PUTTING THE VALUES INSIDE sig_in VECTOR
std::ifstream input("C:\\Users\\Dario\\Desktop\\Tesi\\cartel1.csv");
for(std::string line; getline(input,line);)
{for(int i=0; i<c; i++){input >> sig_in[i];}}
//EVALUATING DFT, IDFT, FILTERED SIGNAL, HILBERT TRANSFORM
forwardDFT(sig_in,c,a,b);
inverseDFT(a,b,c,sig_out);
//FILTERING THE SIGNAL
filter(M,n,D,d,sig_in,F,s);
forwardDFT(F,c,a,b); //REDEFINING A and B with the new filtered signal
hilbert(m,l,c,hil_out,a,b); //filtered hilbert

//WRITING THE OUTPUTS (real part of ft, imag part of ft, abs value of ft, inverse Ft normalized, Hilbert transform (evaluated with the filtered signal)
ofstream output;
output.open("C:\\Users\\Dario\\Desktop\\Tesi\\excel\\file.csv");
output <<"Real" << ";" << "Imag" <<";" << "FT" << ";" << "Inverse FT"<< ";" <<"Filtered signal" << ";" << "Hilbert transform" << endl;
for(int i=0; i<s; i++)
{
output <<a[i]<<";"<<b[i]<<";"<< sqrt((a[i]*a[i]) + (b[i]*b[i]))<< ";" <<sig_out[i]/s << ";" << F[i] << ";" << hil_out[i]/s << endl;
}
output.close();
system("pause");
return 0;
}