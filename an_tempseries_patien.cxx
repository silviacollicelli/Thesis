using namespace std;

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <string>

//ANALISI DATI PER PAZIENTI MALATI
//Ho serie temporale di dati per ogni regione e devo calcolare dfc speed 
double gen_gauss (double mean, double sigma, double xi, double xf);
double cov(int N, double* x, double* y);

int main(){

//Parametri necessari: pesi Cij (matrice nxn), vettore errore (n elementi)
int n=119; //dimensione matrice Cij

int l0=0; //Valore di partenza per il campione di tempi per la cov
int l=l0;

string name3= ".\\data\\patients\\AMC\\FCS_";
string name3_0= ".\\data\\patients\\AMC\\FCS_0"; 
string namef3= "_A_nov.txt";

int np=127; int npmax=196;
ifstream AMCp[np];

int r=0;
int y=0;
int np_start=24;
string names[np];
for (int i=np_start; i<np+np_start; i++){
    if (i+r<100) {
        AMCp[i-np_start].open(name3_0 + to_string(i+r) + namef3);
        names[i-np_start]=name3_0 + to_string(i+r) + namef3;
        //cout << i-np_start << "\t" << name3_0 + to_string(i+r) + namef3;
    }
    else {
        AMCp[i-np_start].open(name3 + to_string(i+r) + namef3);
        names[i-np_start]=name3 + to_string(i+r) + namef3;
        // cout << i-np_start << "\t" << name3 + to_string(i+r) + namef3;

    }
    if (!AMCp[i-np_start]){
        i=i-1;
        r=r+1;
        //cout << "\tNO";
    }
    //cout << endl;
}
//return 0;
ofstream outVp[3];
outVp[0].open(".\\V\\Vp_data_short.txt");
outVp[1].open(".\\V\\Vp_data_medium.txt");
outVp[2].open(".\\V\\Vp_data_long.txt");
ofstream outFC(".\\FC\\FC_pd.txt");
ofstream outFC2(".\\FC\\FC_pd105.txt");

string line;
int N=0;
double q;
for (int my=0; my<np; my++){
    cout << "PAZIENTE malato N\t" << my << "\tFile " << names[my] << endl;
    if(my==8 || my==56 || my==81 || my==82 || my==87 || my==98 || my==111 || my ==121 || my==124 || my==125) cout << "File danneggiato" << endl;
    else{
        if (AMCp[my].good()){
            N=0;
            while (getline(AMCp[my], line)){
                N++;
            }
        }
        double **X= new double*[N];
        for (int i=0; i<N; i++){
            X[i]=new double[n];
        }
        int j=0;    
        AMCp[my].clear();
        AMCp[my].seekg(0, AMCp[my].beg);
        //if (my==56) cout << N << endl;
        while (!AMCp[my].eof()){
            for(int i=0; i<n; i++){
                AMCp[my] >> X[j][i];
                //if (my==56 && j==0) cout << "\t" << X[j][i] << endl; 
            }
            j++;
            //if (my==56) cout << j << endl;
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
       
        //cout <<"m = " << m  << " w " << w << endl;
        
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
                    if (my==0 && q==1 && m==min_l+11 && r==1) outFC << fc << "\t";
                    if (my==0 && q==2 && m==max_l && r==1) outFC2 << fc << "\t";
	                //cout << F_mat[r][cc] << "\t";
	                cc++; //outFC << "\t";s
	                
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
	        if(isnan(v[r-1])) {
                cout << my << "\tErrore NaN" << endl; 
                break;
            }
            outVp[q] << v[r-1] << endl; 

        }

     
        for (int r=0; r<w; r++){
            delete [] F_mat[r];
        }      
		delete [] F_mat;
    
    } // end for m
    
} // end for q

    for (int r=0; r<n; r++){
        delete [] X[r];
    }      
    delete [] X;

    }//End for else
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
