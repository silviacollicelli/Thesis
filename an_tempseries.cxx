using namespace std;

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <string>

//ANALISI DATI PER PAZIENTI SANI
//Ho serie temporale di dati per ogni regione e devo calcolare dfc speed 
double gen_gauss (double mean, double sigma, double xi, double xf);
double cov(int N, double* x, double* y);

int main(){

//Parametri necessari: pesi Cij (matrice nxn), vettore errore (n elementi)
int n=119; //dimensione matrice Cij

int l0=0; //Valore di partenza per il campione di tempi per la cov
int l=l0;

string name1= ".\\data\\controls_session1\\AMC\\FCS_0";
string name2= ".\\data\\controls_session2\\AMC\\FCS_0";
string name1_0= ".\\data\\controls_session1\\AMC\\FCS_00"; 
string name2_0= ".\\data\\controls_session2\\AMC\\FCS_00"; 

string namef1= "_AMC_nov.txt";
string namef2= "_AMC2_nov.txt";

int n1=31;  int n1max=36;
int n2=26;  int n2max=34;
ifstream AMCh[n1+n2];

int r=0;
int y=0;
for (int i=0; i<n1; i++){
    if (i<9) AMCh[i].open(name1_0 + to_string(i+r+1) + namef1);
    else AMCh[i].open(name1 + to_string(i+r+1) + namef1);
    if (!AMCh[i]){
        i=i-1;
        r=r+1;
    }
}
r=0;
for (int i=n1; i<n1+n2; i++){
    if (i-n1<8) AMCh[i].open(name2_0 + to_string(i+y+1-n1) + namef2);
    else AMCh[i].open(name2 + to_string(i+y+1-n1) + namef2);
    if (!AMCh[i]){
        i=i-1;
        y=y+1;
    }
}


ofstream outVh[3];
outVh[0].open("Vh_data_short.txt");
outVh[1].open("Vh_data_medium.txt");
outVh[2].open("Vh_data_long.txt");

string line;
int N=0;

//Ho raggruppato insieme sessione 1 e sessione 2
for (int my=0; my<n1+n2; my++){
    if(my<n1) cout << "Sessione 1\nPAZIENTE " << my << endl;
    else cout << "Sessione 2\nPAZIENTE " << my-n1 << endl; 
double a;
    if (AMCh[my].good()){
        N=0;
        int j=0; int i=0;
        while (getline(AMCh[my], line)){
            N++;
        }
    }
    double **X= new double*[N];
    for (int i=0; i<N; i++){
        X[i]=new double[n];
    }
    int j=0; int i=0;    
    AMCh[my].clear();
    AMCh[my].seekg(0, AMCh[my].beg);
    while (!AMCh[my].eof()){
        for(int i=0; i<n; i++){
            AMCh[my] >> X[j][i];
        }
        j++;
    }

double sum=0.;
int count=0;
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
 
        int w=N/m; //numero di matrici FC calcolabili
	    int T1 = w*m;
 
        double **F_mat = new double *[w];   //w matrici n*n

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
      
               for (int s=0; s<m; s++){     //Considero i punti della serie temporali a distanza di 1s uno dall'altro
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
	        outVh[q] << v[r-1] << endl; 
        }

     
        for (int r=0; r<w; r++){
            delete [] F_mat[r];
        }      
		delete [] F_mat;
    
    } // end for m
    
} // end for q

cout << "done " << my << endl;
} // end for my

    return 0;
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
