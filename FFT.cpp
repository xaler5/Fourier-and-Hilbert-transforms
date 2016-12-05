#include<iostream>
#include<complex>
using namespace std;

int log2(int N)  //logaritmo in base 2 di un intero
{
int k = N, i = 0;
while(k) {k >>= 1; i++;}
return i-1;
}

int potenza2(int n)// la funzione controlla se il numero in argomento è una potenza di due
{                  // e in caso non lo sia genera h che è la potenza di due successiva a n dato
 int h=n;          // infine la funzione assume come valore o n o la potenza più vicina h
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

int invertibit(int N, int n)//N è numero di samples e n è l'indice intero 
{                           //(in bit) che rappresenta i singoli N
  int j, r = 0;             //la funzione assume come valore il reverse bit di n
  for(j = 1; j <= log2(N); j++) 
  { 
  if(n & (1 << (log2(N) - j)))
  {r |= 1 << (j - 1);}
  }
return r; 
};

void ordina(complex<double>* vec, int n)//dispone gli elementi del vettore ordinandoli 
{                                      //in ordine dato dal reverse bit precedente
const int c=n;                                      
complex<double> vec2[c];
 for(int i = 0; i < n; i++)
 {vec2[i] = vec[invertibit(n, i)];}
 for(int j = 0; j < n; j++)
 {vec[j] = vec2[j];}
}

void FFT(complex<double>* vec1, int N, double t)
{
 complex<double> W[N / 2];        //vettore dei pesi NB sono sempre in numero N/2 (di cui uno è sempre 1).
 W[1] = polar(1., -2. * M_PI / N); //polar restituisce (1*cost(theta),1*sin(theta) 
 W[0] = 1;                        //ovvero la forma cartesiana
 for(int i = 2; i < N / 2; i++)
 {
 W[i] = pow(W[1], i);             //calcolo i vari W e li assegno al vettore dei pesi
 }
 int n = 1;
 int a = N / 2;
 for(int j = 0; j < log2(N); j++)//Comincia l'iterazione dell'algoritmo di Danielson-Lanczos
 {                               //log2N iterazioni
 for(int i = 0; i < N; i++)   //ognuna per N volte ---> Nlog2(N)
 {
  if(!(i & n)) 
  {
  complex<double> temp = vec1[i];                            
  complex<double> Temp = W[(i * a) % (n * a)] * vec1[i + n]; 
  vec1[i] = temp + Temp;      // NB il ripetersi con periodo N/2 dei pesi ma con segno opposto fa si che ci
  vec1[i + n] = temp - Temp;  // siano N/2 addizioni e N/2 sottrazioni
  }
 }
 n *= 2;    //Ad ogni passo n aumenta.. ad esempio con N=8 n=1,2,4
 a = a / 2; //e lo stesso per a.. con N=8 a=4,2,1
 };
for(int o = 0; o<N; o++)
{vec1[o] *= t;} //moltiplico il vettore trasformato (DFT) per il tempo di campionamento
}                                 //e ottengo la FFT definitiva

int main () 
{
double t; //tempo di campionamento
int s; //numero samples iniziali
cout << "Inserire il tempo di campionamento (s): " << endl;
cin >> t;
cout << "Inserire il numero di samples: " << endl; cin >> s;
int r=potenza2(s); //r=s se s  è potenza di 2 altrimenti è la potenza di due successiva
int d=r-s;         //d è 0 se r=s altrimenti è il numero di 0 da aggiungere al set di dati
const int c=r;
complex<double> vec1[c];
cout << "Inserire i dati nella forma 'reale immaginaria': " << endl;
 for(int i = 0; i < s; i++) //RIEMPIO IL VETTORE CON I DATI 
 {
 cout << "inserisci la componente " << i+1 << endl;
 double real, imag;
 cin >> real; //>> imag;
 complex<double> p(real, 0);
 vec1[i]=p;
 p=0; real=0; imag=0;
 }
 for(int j=s; j<r-1; j++) //RIEMPIO DI 0 i rimanenti termini
 {
 vec1[j]=0;
 }  
//for(int j=0; j<r; j++) {cout << vec1[j] << " ";} //vettore non invertito
//cout << " " << endl;
ordina(vec1,r);
//for(int j=0; j<r; j++) {cout << vec1[j] << " ";} //vettore invertito
//cout << " " << endl;
FFT(vec1,r,t);
for(int j=0; j<r; j++) {cout << vec1[j] << " ";} //vettore trasformato
cout << " " << endl;
system("pause");
return 0;
};
