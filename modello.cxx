using namespace std;

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <time.h>

double gen_gauss (double mean, double sigma, double xi, double xf);
double cov(int N, double* x, double* y);

//eq modello: dx_{i}=-x_{i}/tau+sum_{k=!i}C_{i,k}x_{k}+dB_{i}

int main(){

//Parametri necessari: pesi Cij (matrice nxn), vettore errore (n elementi)
double tau; 
int n; //dimensione matrice Cij
double T; // durata simulazione
double dt;

n=100; 
T=1000; 
int m=10; //prova: uso 10 valori per la cov 
int l0=0; //Valore di partenza per il campione di tempi per la cov
int l=l0;

double *Cij_mat, *xt, *xt1, *xt1_fin, *er, *Xij_mat;
double **Cij, **Xij; //Xij Matrice che contiene elementi dell'evoluzione temporale fino a t=j*tau per ogni regione (100 righe x m colonne) per fare cov
Cij_mat = (double*) new double[n*n];
Cij=(double**) new double*[n];
Xij_mat = (double*) new double[m*n];
Xij=(double**) new double*[n];
xt=(double*) new double[n];
xt1=(double*) new double[n];
xt1_fin=(double*) new double[n];
er=(double*) new double[n];
for (int i=0; i<n; i++){
    Cij[i]=&Cij_mat[n*i];
    Xij[i]=&Xij_mat[m*i];
}
ifstream f ("example_real_Sigma.txt"); // file txt con parametri per simulazione e sigma (tutto in colonna)
if (f.good()){
    for (int i=0; i<n; i++){
       f >> er[i];
    }
}
f.close();

//Leggere la matrice dei pesi dimensione n
ifstream fin ("example_real.txt");
if (fin.good()){
    int j=0; int i=0;
    while (!fin.eof()){
        for(int i=0; i<n; i++){
            fin >> Cij[j][i];
        }
        j++;
    }
}
fin.close();
tau=abs(1./Cij[0][0]);
dt=0.005*tau;
int N=round(T/dt); //step per la simulazione
cout << N << endl;
//dati iniziali creati da programma 
for (int i=0; i<n; i++){
    xt[i]=1.;
}

ofstream fout ("ris_evol.txt");
/*ris_sim.txt file con risultati simulazione con file matriciale ij
elemento di matrice x_{i, j}= evoluzione temporale al tempo t=j*dt della regione i-esima*/
ofstream outFC ("FC1.txt");  //Matrice di correlazione

fout << "#t\t";
for(int i=0; i<n; i++){
    fout << "reg " << i+1 << '\t'; 
}
fout << endl;
double sum=0.;

//metodo Eulero
for (int j=0; j<N; j++){ // ciclo temporale
    if(j%100==1) fout << j*dt << '\t';
    for (int i=0; i<n; i++){    //ciclo spaziale
        if (j%100==1) fout << xt[i] << '\t';
        for (int k=0; k<n; k++){
            if(k!=i) sum+=Cij[i][k]*xt[k];
        }
        xt1[i]= xt[i]+(-xt[i]/tau+sum)*dt;
        sum=0;
        xt1_fin[i]=xt1[i]+gen_gauss(0, sqrt(er[i]*dt), -3*sqrt(er[i]*dt), 3*sqrt(er[i]*dt));
    }
    if (floor((j-1)*dt)==l && floor(j*dt)==l+1 && l<l0+m){ //floor approssima al più piccolo intero, round approssima all'intero più vicino
        for(int i=0; i<n; i++){
            Xij[i][l-l0]=xt[i];
        }
        l++;
    }

    //Aggiornamento
    for(int i=0; i<n; i++){
        xt[i]=xt1_fin[i];
    }
    if(j%100==1) fout << endl;
    if (j%1000==1) cout << "fatto " << j << endl;
}

//Calcolo correlazione per tempi da 1s a 10s campionati ogni secondo 
for (int s=0; s<n; s++){
    for (int i=0; i<n; i++){
        outFC << cov(m, Xij[i], Xij[s])/sqrt(cov(m, Xij[i], Xij[i])*cov(m, Xij[s], Xij[s])) << '\t'; 
    }
    outFC << endl;
}

outFC.close();
fout.close();

delete [] xt;
delete [] xt1;
delete [] xt1_fin;
delete [] er;
delete [] Cij;
delete [] Cij_mat;
delete [] Xij;
delete [] Xij_mat;

    return 0;
}


double gen_gauss (double mean, double sigma, double xi, double xf){
    double g, pmax, xgauss, x;
    double r[2];
    pmax= 1/sqrt(2*M_PI*sigma*sigma);
    srand(time(NULL));
    do{
        for(int j=0; j<2; j++){
            r[j]=(double) rand()/RAND_MAX; //genera numero casuale tra 0 e 1 
        }
        xgauss=0;
        x= (xf-xi)*r[0]+xi;
        g=exp(-pow(x-mean, 2)/2*pow(sigma, 2))/sqrt(2*M_PI*pow(sigma, 2));
    } while(r[1]>=g/pmax);
    xgauss=x;
    return xgauss;
}

//covarianza campionaria
double cov(int N, double* x, double* y){
    double covxy=0;
    double sumx=0;
    double sumy=0;
    for(int i=0; i<N; i++){
        sumx+=x[i];
        sumy+=y[i];
    }
    double xm=sumx/N;
    double ym=sumy/N;
    double sumxy=0;
    for(int i=0; i<N; i++){
        sumxy+=(x[i]-xm)*(y[i]-ym);
    }
    covxy=sumxy/(N-1);
    return covxy;
}