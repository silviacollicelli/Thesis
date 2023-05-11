using namespace std;

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include<time.h>

int main(){
    ofstream fout ("prova_gauss.txt");
    double xgauss[100]; 
    double r[2];
    double g, pmax, xf, xi, x;
    double sigma=1;
    double mean=10;
    xi=0; xf=20; 
    pmax= 1/sqrt(2*M_PI*sigma*sigma);
    srand(time(NULL));
    for (int i=0; i<100; i++){
        do{
            for(int j=0; j<2; j++){
                r[j]=(double) rand()/RAND_MAX; //genera numero casuale tra 0 e 1 
                if (i==0) cout << r[j] << endl;
            }
            xgauss[i]=0;
            x= (xf-xi)*r[0]+xi;
            g=exp(-pow(x-mean, 2)/2*pow(sigma, 2))/sqrt(2*M_PI*pow(sigma, 2));
        } while(r[1]>=g/pmax);
    xgauss[i]=x;
    fout << xgauss[i] << endl;
    }
   
}
