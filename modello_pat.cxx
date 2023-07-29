//PROGRAMMA PER SIMULAZIONE SERIE TEMPORALE E CONSEGUENTE ANALISI DI FC, DFC PER PAZIENTI MALATI

using namespace std;

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <string>

double gen_gauss (double mean, double sigma, double xi, double xf);
double cov(int N, double* x, double* y);

//eq modello: dx_{i}=-x_{i}/tau+sum_{k=!i}C_{i,k}x_{k}+dB_{i}

int main(){

double tau; 
int n; //dimensione matrice Cij
double T; // durata simulazione
double dt;

n=119; 
T=1000; 

int l0=0; //Valore di partenza per il campione di tempi per la cov
int l=l0;

double *Cij_mat, *xt, *xt1, *xt1_fin, *er, *erij_mat;
double **Cij, **erij;

Cij_mat = (double*) new double[n*n];
Cij=(double**) new double*[n];
erij_mat = (double*) new double[n*n];
erij=(double**) new double*[n];
double **X= new double*[(int)T]; //Xij Matrice che contiene elementi dell'evoluzione temporale ad ogni t intero per ogni regione T*n 
xt=(double*) new double[n];
xt1=(double*) new double[n];
xt1_fin=(double*) new double[n];
er=(double*) new double[n];
for (int i=0; i<n; i++){
    Cij[i]=&Cij_mat[n*i];
    erij[i]=&erij_mat[n*i];
}

for (int i=0; i<(int)T; i++){
    X[i]=new double[n];
}

double a=0; 
string name1_0= ".\\data\\patients\\AMCJ\\FCS_0";
string name2_0= ".\\data\\patients\\AMCS\\FCS_0";
string name1= ".\\data\\patients\\AMCJ\\FCS_"; 
string name2= ".\\data\\patients\\AMCS\\FCS_"; 

string namef1= "_A_nov_J.txt";
string namef2= "_A_nov_Sigma.txt";

int np=127; int npmax=196;
int np_start=24;
ifstream AMCp[2][np];
string names[np];

int r=0;
int y=0;

//APERTURA FILE INPUT CON MATRICI J E SIGMA

for (int s=0; s<2; s++){
    if(s==0){   //S=0 indica matrici J
        for (int i=np_start; i<np+np_start; i++){
            if (i+r<100) {
                AMCp[s][i-np_start].open(name1_0 + to_string(i+r) + namef1);
                names[i-np_start]=name1_0 + to_string(i+r) + namef1;
                //cout << i-np_start << "\t" << name1_0 + to_string(i+r) + namef1;
            }
            else {
                AMCp[s][i-np_start].open(name1 + to_string(i+r) + namef1);
                names[i-np_start]=name1 + to_string(i+r) + namef1;
                //cout << i-np_start << "\t" << name1 + to_string(i+r) + namef1;

            }
            if (!AMCp[s][i-np_start]){
                i=i-1;
                r=r+1;
                //cout << "\tNO";
            }
            //cout << endl;
        }
    }
    //return 0;
    if(s==1){   //s=1 indica matrice sigma
        for (int i=np_start; i<np+np_start; i++){
            if (i+r<100) AMCp[s][i-np_start].open(name2_0 + to_string(i+r) + namef2);
            else AMCp[s][i-np_start].open(name2 + to_string(i+r) + namef2);
            if (!AMCp[s][i-np_start]){
                i=i-1;
                r=r+1;
            }
        }
    }
    r=0;
}

ofstream outV[3];
outV[0].open("Vp_short.txt");
outV[1].open("Vp_medium.txt");
outV[2].open("Vp_long.txt");


for (int my=0; my<np; my++){    //Ciclo sui file

cout << "SIM paziente MALATO\tN " << my << "\t" << names[my] << endl;
if(my==8 || my==56 || my==81 || my==82 || my==87 || my==98 || my==111 || my ==121 || my==124 || my==125) {
    cout << "File danneggiato"<< endl;  //File con valori nan o serie temporale nulla
}
else{
if (AMCp[1][my].good()){
   int j=0; int i=0;
    while (!AMCp[1][my].eof()){
        for(int i=0; i<n; i++){
            AMCp[1][my] >> erij[j][i];
            if (i==j) er[i]= erij[j][i];
        }
        j++;
    }
    AMCp[1][my].close();
}

if (AMCp[0][my].good()){
    int j=0; int i=0;
    while (!AMCp[0][my].eof()){
        for(int i=0; i<n; i++){
            AMCp[0][my] >> Cij[j][i];
        }
        j++;
    }
}

tau=abs(1./Cij[0][0]);
dt=0.005*tau;
int N=round(T/dt); //step per la simulazione

//dati iniziali creati da programma 
for (int i=0; i<n; i++){
    xt[i]=0.;
}

double sum=0.;
int count=0;

//metodo Eulero
for (int j=0; j<N; j++){ // ciclo temporale
    for (int i=0; i<n; i++){    //ciclo spaziale     
        sum=0;
        for (int k=0; k<n; k++){
            if(k!=i) sum+=Cij[i][k]*xt[k];
        }
        xt1[i]= xt[i]+(-xt[i]/tau+sum)*dt;
        double noise=gen_gauss(0, sqrt(er[i]*dt), -3*sqrt(er[i]*dt), 3*sqrt(er[i]*dt));
        xt1_fin[i]= xt1[i]+ noise; 
    }
    
    if (floor((j-1)*dt)==l && floor(j*dt)==l+1){ //floor approssima al più piccolo intero, round approssima all'intero più vicino
        for(int i=0; i<n; i++){
            X[l][i]=xt1_fin[i];
        }
        l++;
    }

    //Aggiornamento
    for(int i=0; i<n; i++){
        xt[i]=xt1_fin[i];
    }
    //if (j%10000==1) cout << "fatto " << j << endl;
}

/*ofstream TS("ts.txt");  //evoluzione temporale per tempi interi

for(int i=0; i<(int)T; i++) {
    for(int j=0; j<n;j++) {
        TS << X[i][j] << "\t";
    }
    TS << endl;
}*/

//TS.close();

//cout << "\tended simulations" << endl;

//Calcolo di FC e dFC per tre gruppi di finestre temporali: short, medium e long
int m=0;
int n_durations = 3; 
double fc;

for(int q=0; q<3; q++) {
  
	int min_l, max_l;   //lunghezze max e min delle finestre temporali

	if(q==0)  {
		min_l=3;
		max_l=8;
	}
		
	else if(q==1)  {
		min_l=9;
		max_l=32;
	}

	else if(q==2)  {
		min_l=33;
		max_l=105;
	}
  
    for(m=min_l; m<=max_l; m++) {
 
        int w=T/m; //numero di matrici FC calcolabili
	    int T1 = w*m;
 
        double **F_mat = new double *[w];   //w matrici n*n

    	for (int r=0; r<w; r++){
          F_mat[r]=new double [n*n];
        }      
       
        //cout <<"m = " << m  << " w " << w << endl;
        
        for(int r=0; r<w; r++){

    	    double **xij;
		    xij = new double *[n];
            int tmin=r*m;

            for (int i=0; i<n; i++){
                xij[i] = new double[m];
      
               for (int s=0; s<m; s++){
                   xij[i][s] = X[tmin+s][i];
	    	   }
	        }

	        int cc=0;
            for (int s=0; s<n; s++){
                for (int i=0; i<n; i++){
                    fc = cov(m, xij[i], xij[s])/sqrt(cov(m, xij[i], xij[i])*cov(m, xij[s], xij[s]));
			        F_mat[r][cc] = fc; //ogni matrice sta in un'unica riga
	                //cout << F_mat[r][cc] << "\t";
	                cc++; //outFC << "\t";s
	                
                }
            }

            for (int i=0; i<n; i++){
                delete [] xij[i];
		    }
            delete [] xij;

        } // end for r   
         
        double *v = new double [w-1];   //dfc speed
        for(int r=1; r<w; r++){
	        v[r-1] = 1-cov(n*n, F_mat[r], F_mat[r-1])/sqrt(cov(n*n, F_mat[r], F_mat[r])*cov(n*n, F_mat[r-1], F_mat[r-1]));
	        outV[q] << v[r-1] << endl; 
        }

     
        for (int r=0; r<w; r++){
            delete [] F_mat[r];
        }      
		delete [] F_mat;
    
    } // end for m
    
} // end for q

} //end for else
} // end for my
    return 0;
}


double gen_gauss (double mean, double sigma, double xi, double xf){
    double g, pmax, xgauss, x;
    double r[2];
    pmax= 1/sqrt(2*M_PI*sigma*sigma);
    
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
