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
/*int m=30; //uso 10 valori per la cov 
int w=T/m; //numero di matrici FC calcolabili
T = w*m;*/

//m=10s->short window size range
int l0=0; //Valore di partenza per il campione di tempi per la cov
int l=l0;

double *Cij_mat, *xt, *xt1, *xt1_fin, *er;

//double **F_mat = new double *[w];

double **Cij;
Cij_mat = (double*) new double[n*n];
Cij=(double**) new double*[n];
double **X= new double*[(int)T]; //Xij Matrice che contiene elementi dell'evoluzione temporale fino a t=j*tau per ogni regione wxnxm 
xt=(double*) new double[n];
xt1=(double*) new double[n];
xt1_fin=(double*) new double[n];
er=(double*) new double[n];
for (int i=0; i<n; i++){
    Cij[i]=&Cij_mat[n*i];
}

for (int i=0; i<(int)T; i++){
    X[i]=new double[n];
    //F_mat[i]=new double [n*n];
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
        //count=floor(l/m);
        
        for(int i=0; i<n; i++)
        {
            X[l][i]=xt1_fin[i];
        }
        l++;
    }

    //Aggiornamento
    for(int i=0; i<n; i++){
        xt[i]=xt1_fin[i];
    }
    if (j%1000==1) cout << "fatto " << j << endl;
}

ofstream TS("ts.txt");

for(int i=0; i<(int)T; i++) {
    for(int j=0; j<n;j++) {
        TS << X[i][j] << "\t";
    }
    TS << endl;
}

TS.close();


cout << "ended simulations" << endl;

int m=0;


int n_durations = 3;
double fc;


for(int q=0; q<3; q++) {
  
	int min_l, max_l;
	ofstream outV;

	if(q==0)  {
		outV.open("V_short.txt");  //Matrice di correlazione
		min_l=3;
		max_l=8;
	}
		
	else if(q==1)  {
		outV.open("V_medium.txt");  //Matrice di correlazione
		min_l=9;
		max_l=32;
	}

	else if(q==2)  {
		outV.open("V_long.txt");  //Matrice di correlazione
		min_l=33;
		max_l=105;
	}

  
    for(m=min_l; m<=max_l; m++) {
 
        int w=T/m; //numero di matrici FC calcolabili
	    int T1 = w*m;
 
        double **F_mat = new double *[w];

    	for (int r=0; r<w; r++){
          F_mat[r]=new double [n*n];
        }      
       
        cout <<"m = " << m  << " w " << w << endl;
        
     
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
			        F_mat[r][cc] = fc;
	                //cout << F_mat[r][cc] << "\t";
	                cc++; //outFC << "\t";s
	                
                }
                //cout << endl;
            }
          
           //cout << endl;
           
           

            for (int i=0; i<n; i++){
                delete [] xij[i];
		    }
            delete [] xij;

        } // end for r   

         
        double *v = new double [w-1];

        for(int r=1; r<w; r++){
	        v[r-1] = 1-cov(n*n, F_mat[r], F_mat[r-1])/sqrt(cov(n*n, F_mat[r], F_mat[r])*cov(n*n, F_mat[r-1], F_mat[r-1]));
	        outV << v[r-1] << endl; 
        }

     
        for (int r=0; r<w; r++){
            delete [] F_mat[r];
        }      
		delete [] F_mat;
    
    } // end for m
    

} // end for q


/*

double fc;

//Calcolo correlazione per tempi da 1s a 10s campionati ogni secondo 
for(int r=0; r<w; r++){
	int cc=0;
    for (int s=0; s<n; s++){
        for (int i=0; i<n; i++){
            fc = cov(m, Xij[r][i], Xij[r][s])/sqrt(cov(m, Xij[r][i], Xij[r][i])*cov(m, Xij[r][s], Xij[r][s]));
			F_mat[r][cc] = fc;
	        cc++; //outFC << "\t";s 
        }
    }
   //outFC<< endl;
}

//outFC.close();
fout.close();

double *v = new double [w-1];

ofstream outV ("V.txt");  //Matrice di correlazione

for(int r=1; r<w; r++){
	v[r-1] = 1-cov(n*n, F_mat[r], F_mat[r-1])/sqrt(cov(n*n, F_mat[r], F_mat[r])*cov(n*n, F_mat[r-1], F_mat[r-1]));
	outV << v[r-1] << endl;
}


outV.close();

delete [] xt;
delete [] xt1;
delete [] xt1_fin;
delete [] er;
delete [] Cij;
delete [] Cij_mat;
for (int i=0; i<w; i++){
    for (int j=0; j<n; j++) {
        delete [] Xij[i][j];
    }
    delete [] Xij[i];
    delete [] F_mat[i];
}
delete [] Xij;
delete [] F_mat;
*/


    return 0;
}


double gen_gauss (double mean, double sigma, double xi, double xf){
    double g, pmax, xgauss, x;
    double r[2];
    pmax= 1/sqrt(2*M_PI*sigma*sigma);
    
    //cout << "time " << time(NULL) << endl; 
    //srand(time(NULL));
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
