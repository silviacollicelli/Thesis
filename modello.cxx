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
tau=1/0.681; 
dt=0.005*tau; 
T=1000; 

double *Cij_mat, *xt, *xt1, *xt1_fin, *er;
double **Cij; 
Cij_mat = (double*) new double[n*n];
Cij=(double**) new double*[n];
xt=(double*) new double[n];
xt1=(double*) new double[n];
xt1_fin=(double*) new double[n];
er=(double*) new double[n];
for (int i=0; i<n; i++){
    Cij[i]=&Cij_mat[i*n];
}

int m=round(tau/dt); 

double *Xij_mat;
double **Xij; //Matrice che contiene elementi dell'evoluzione temporale fino a t=j*tau per ogni regione (100 righe x m colonne) per fare cov
Xij_mat = (double*) new double[m*n];
Xij=(double**) new double*[m];
for (int i=0; i<m; i++){
    Xij[i]=&Xij_mat[n*i];
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

//dati iniziali creati da programma     Possono andare bene?
for (int i=0; i<n; i++){
    xt[i]=1.;
}

ofstream fout ("ris_sim.txt");
/*ris_sim.txt file con risultati simulazione con file matriciale ij
elemento di matrice x_{i, j}= evoluzione temporale al tempo t=j*dt della regione i-esima*/
ofstream out200 ("FC200.txt");  //Matrice di correlazione

int N=T/dt; //step per la simulazione

//matrice RUMORI GAUSSIANI n righe x N colonne
double *G_mat; double **G;
G_mat = (double*) new double[N*n];
G=(double**) new double*[N];
for (int i=0; i<N; i++){
    G[i]=&G_mat[n*i];
}
for(int i=0; i<n; i++){
    for (int j=0; j<N; j++){
        G[i][j]=gen_gauss(0, er[i]*dt, -3*er[i]*dt, 3*er[i]*dt);
    }
    cout << "fatto " << i << endl; 
}

fout << "#t\t";
for(int i=0; i<n; i++){
    fout << "reg " << i+1 << '\t'; 
}
fout << endl;
double sum=0.;

//metodo Eulero
for (int j=0; j<N; j++){ // ciclo temporale
    fout << j*dt << '\t';
    for (int i=0; i<n; i++){    //ciclo spaziale
        fout << xt[i] << '\t';
        for (int k=0; k<n; k++){
            if(k!=i) sum+=Cij[i][k]*xt[k];
        }
        xt1[i]= xt[i]+(-xt[i]/tau+sum)*dt;
        sum=0;
        xt1_fin[i]=xt1[i]+G[i][j];
    }
    if (j>=20000 && j<=20200){ 
        for(int i=0; i<n; i++){
            Xij[i][j-20000]=xt[i];
        }
    }

    //Aggiornamento
    for(int i=0; i<n; i++){
        xt[i]=xt1[i];
    }
    fout << endl;
}

//per cov
//Calcolo cov per t=200*tau

for (int l=0; l<n; l++){
    for (int i=0; i<n; i++){
        out200 << cov(m, Xij[i], Xij[l]) << '\t'; 
    }
    out200 << endl;
}

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