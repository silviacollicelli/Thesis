//PROGRAMMA PER SIMULAZIONE SERIE TEMPORALE E CONSEGUENTE ANALISI DI FC, DFC PER PAZIENTI SANI

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

//Parametri necessari: pesi Cij (matrice nxn), vettore errore (n elementi)
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
string name1a_0= ".\\data\\controls_session1\\AMCJ\\FCS_00";
string name2a_0= ".\\data\\controls_session1\\AMCS\\FCS_00";
string name1a= ".\\data\\controls_session1\\AMCJ\\FCS_0"; 
string name2a= ".\\data\\controls_session1\\AMCS\\FCS_0"; 

string name1b_0= ".\\data\\controls_session2\\AMCJ\\FCS_00";
string name2b_0= ".\\data\\controls_session2\\AMCS\\FCS_00";
string name1b= ".\\data\\controls_session2\\AMCJ\\FCS_0"; 
string name2b= ".\\data\\controls_session2\\AMCS\\FCS_0";

string namef1a= "_AMC_nov_J.txt";
string namef2a= "_AMC_nov_Sigma.txt";
string namef1b= "_AMC2_nov_J.txt";
string namef2b= "_AMC2_nov_Sigma.txt";

int n1=31;  int n1max=36;
int n2=26;  int n2max=34;
int n2_start=2;
string names[n1+n2];
ifstream AMCh[2][n1+n2];

int r=0;
int y=0;
//APERTURA FILE

for (int s=0; s<2; s++){
    if(s==0){
        for (int i=0; i<n1; i++){
            if (i<9) {
                AMCh[s][i].open(name1a_0 + to_string(i+r+1) + namef1a);
                names[i]=name1a_0 + to_string(i+r+1) + namef1a;
            }
            else {
                AMCh[s][i].open(name1a + to_string(i+r+1) + namef1a);
                names[i]=name1a + to_string(i+r+1) + namef1a;
            }
            if (!AMCh[s][i]){
                i=i-1;
                r=r+1;
            }
        }
        for (int i=n2_start; i<n2+n2_start; i++){
            if (i+y<10) {
                AMCh[s][i+n1-n2_start].open(name1b_0 + to_string(i+y) + namef1b);
                names[i+n1-n2_start]=name1b_0 + to_string(i+y) + namef1b;
            }
            else {
                AMCh[s][i+n1-n2_start].open(name1b + to_string(i+y) + namef1b);
                names[i+n1-n2_start]=name1b + to_string(i+y) + namef1b;
            }
            if (!AMCh[s][i+n1-n2_start]){
                i=i-1;
                y=y+1;
            }
        }
    }

    r=0;
    y=0;
    if(s==1){
        for (int i=0; i<n1; i++){
            if (i<9) AMCh[s][i].open(name2a_0 + to_string(i+r+1) + namef2a);
            else AMCh[s][i].open(name2a + to_string(i+r+1) + namef2a);
            if (!AMCh[s][i]){
                i=i-1;
                r=r+1;
            }
        }
        for (int i=n2_start; i<n2+n2_start; i++){
            if (i+y<10) AMCh[s][i+n1-n2_start].open(name2b_0 + to_string(i+y) + namef2b);
            else AMCh[s][i+n1-n2_start].open(name2b + to_string(i+y) + namef2b);
            if (!AMCh[s][i+n1-n2_start]){
                i=i-1;
                y=y+1;
            }
        }
    }
    r=0;
    y=0;
}

ofstream outV[3];
outV[0].open(".\\V\\Vh_short.txt");
outV[1].open(".\\V\\Vh_medium.txt");
outV[2].open(".\\V\\Vh_long.txt");
ofstream outFC(".\\FC\\FC_hs.txt");
ofstream outFC2(".\\FC\\FC_hs105.txt");
ofstream ts ("time.txt");

for (int my=0; my<n1+n2; my++){    //Ciclo sui file

cout << "SIM paziente SANO\tN " << my << "\t" << names[my] << endl;

if (AMCh[1][my].good()){
   int j=0; int i=0;
    while (!AMCh[1][my].eof()){
        for(int i=0; i<n; i++){
            AMCh[1][my] >> erij[j][i];
            if (i==j) er[i]= erij[j][i];
        }
        j++;
    }
    AMCh[1][my].close();
}

//Leggere la matrice dei pesi dimensione n
if (AMCh[0][my].good()){
    int j=0; int i=0;
    while (!AMCh[0][my].eof()){
        for(int i=0; i<n; i++){
            AMCh[0][my] >> Cij[j][i];
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
/*for (int j=0; j<n; j++){
    for (int i=0; i<n; i++){
        if (i==j) Cij[j][i]=-2.6;
        else Cij[j][i]=-my*0.01-0.1;
    }
}
*/

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
            if (my==0) ts << X[l][i] << '\t';
        }
        if (my==0) ts << endl;
        l++;
    }

    for(int i=0; i<n; i++){
        xt[i]=xt1_fin[i];
    }
}
//if (my==0) cout << "fatto" << endl;

/*ofstream TS("ts.txt");  //evoluzione temporale per tempi interi

for(int i=0; i<(int)T; i++) {
    for(int j=0; j<n;j++) {
        TS << X[i][j] << "\t";
    }
    TS << endl;
}*/

//TS.close();


int m=0;
int n_durations = 3; //tre gruppi di finestre temporali: short, medium e long
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
	                cc++;
	                if (my==0 && q==1 && m==min_l+11 && r==1) outFC << fc << "\t";
	                if (my==0 && q==2 && m==max_l && r==1) outFC2 << fc << "\t";
                }
                if (my==0 && q==1 && m==min_l+11 && r==1) outFC << endl;
                if (my==0 && q==2 && m==max_l && r==1) outFC2 << endl;
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
