#include<iostream>
#include<complex>
using namespace std;

int log2(int N)  //logarithm base 2 of an inteer
{
int k = N, i = 0;
while(k) {k >>= 1; i++;}
return i-1;
}

int power2(int n)// function check if the argument is a power of 2
{                  // if not it generates h which is the next power of 2 after the argument
 int h=n;          
 if(!(n > 0 && (n & (n - 1)))){ return n;}
 else                                   
 {                
 while (h & (h-1)) 
 {
 h = h & (h-1);
 }
 h = h << 1;
 }
 return h;
};

int inversebit(int N, int n)//N is the samples number and n an integer index 
{                           //(in bit) that represents each sample
  int j, r = 0;             //this function assume as value the reverse bit order
  for(j = 1; j <= log2(N); j++) 
  { 
  if(n & (1 << (log2(N) - j)))
  {r |= 1 << (j - 1);}
  }
return r; 
};

void order(complex<double>* vec, int n)//reorder the elements in reverse bit order 
{                                      
const int c=n;                                      
complex<double> vec2[c];
 for(int i = 0; i < n; i++)
 {vec2[i] = vec[inversebit(n, i)];}
 for(int j = 0; j < n; j++)
 {vec[j] = vec2[j];}
}

void FFT(complex<double>* vec1, int N, double t)
{
 complex<double> W[N / 2];        //vector of the weights
 W[1] = polar(1., -2. * M_PI / N); //polar gives us (1*cost(theta),1*sin(theta) 
 W[0] = 1;                       
 for(int i = 2; i < N / 2; i++)
 {
 W[i] = pow(W[1], i);             //calculating each weight and assigning it to the vector
 }
 int n = 1;
 int a = N / 2;
 for(int j = 0; j < log2(N); j++)//It starts the iteration of the Danielson-Lanczos algorithm
 {                               //log2N iterations
 for(int i = 0; i < N; i++)      //each of them N times ---> Nlog2(N)
 {
  if(!(i & n)) 
  {
  complex<double> temp = vec1[i];                            
  complex<double> Temp = W[(i * a) % (n * a)] * vec1[i + n]; 
  vec1[i] = temp + Temp;      // NB having a period of N/2 of the weights with opposite signs make it possible to have
  vec1[i + n] = temp - Temp;  // N/2 additions and N/2 substractions 
  }
 }
 n *= 2;   
 a = a / 2; 
 };
for(int o = 0; o<N; o++)
{vec1[o] *= t;} //Normalization by the sampling time
}

int main () 
{
double t; //sampling time
int s; //number of initial samples
cout << "Insert sample time interval (s): " << endl;
cin >> t;
cout << "Insert number of samples: " << endl; cin >> s;
int r=power2(s); 
int d=r-s;        
const int c=r;
complex<double> vec1[c];
cout << "Insert data in the form 'real-part  imaginary-part': " << endl;
 for(int i = 0; i < s; i++) //filling the data vector 
 {
 cout << "Insert component number:  " << i+1 << endl;
 double real, imag;
 cin >> real; //>> imag;
 complex<double> p(real, 0);
 vec1[i]=p;
 p=0; real=0; imag=0;
 }
 for(int j=s; j<r-1; j++) //Adding zeros
 {
 vec1[j]=0;
 }  
order(vec1,r);
FFT(vec1,r,t);
for(int j=0; j<r; j++) {cout << vec1[j] << " ";} //transformed vector
cout << " " << endl;
system("pause");
return 0;
};
