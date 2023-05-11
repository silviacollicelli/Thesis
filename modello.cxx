using namespace std;

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include<time.h>

double gen_gauss (double mean, double sigma, double xi, double xf);

//eq modello: dx_{i}=-x_{i}/tau+sum_{k=!i}C_{i,k}x_{k}+dB_{i}

int main(){

//Parametri necessari: pesi Cij (matrice nxn), tau_{x}, dati iniziali vettore x (n elementi), vettore errore (n elementi)
double tau=1; //parametro del modello
int n; //dimensione matrice Cij
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


//Leggere la matrice dei pesi dimensione n
ifstream fin ("pesi_Cij.txt");
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

ifstream f ("x0_iniziale.txt"); //file txt con elementi in colonna
if (f.good()){
    for (int i=0; i<n; i++){
       f >> xt[i];
    }
}

ifstream s ("sigma.txt"); // file txt con errore gaussiano per ogni regione iesima colonna
if (s.good()){
    for (int i=0; i<n; i++){
       f >> er[i];
    }
}

ofstream fout ("ris_sim.txt");
/*ris_sim.txt file con risultati simulazione con 
file matriciale ij
riga i= evoluzione temporale di tale regione
colonna xj= regione i-esima cervello
elemento di matrice x_{i, j}= evoluzione temporale al tempo t=j*dt della regione i-esima
*/

//PARAMETRI SIMULAZIONE
double T=10; // durata simulazione
double dt=0.05*tau;
int N=T/dt; //step per la simulazione
double sum=0.;

//metodo Eulero
//sigma vettore? allora sigma è uguale per la stessa regione e quindi per ogni step dell'evoluzione temporale
for (int j=0; j<N; j++){ // ciclo temporale
    if (j==0) {
        fout << "t\t";
        for(int i=0; i<n; i++){
            fout << "reg " << i+1 << '\t'; 
        }
        fout << endl; 
    }
    fout << j*dt << '\t';
    for (int i=0; i<n; i++){    //ciclo spaziale
        fout << xt[i] << '\t';
        for (int k=0; k<n; k++){
            if(k!=i) sum+=Cij[i][k]*xt[k];
        }
        xt1[i]= xt[i]+(-xt[i]/tau+sum)*dt;
        sum=0;
        //Generazione numero casuale che segue distr gaussiana
        double a, b; //estremi distribuzione (+-3 sigma?)
        a=xt1[i]-3*er[i];
        b=xt1[i]+3*er[i];
        xt1_fin[i]=xt1[i]+gen_gauss(0, er[i]*dt, a, b);
    }
    for(int i=0; i<n; i++){
        xt[i]=xt1[i];
    }
    fout << endl;
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