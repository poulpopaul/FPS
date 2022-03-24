#include<iostream>
#include"time.h"
using namespace std;
#include"projet.h"
#include<fstream>
int main(){
int Nx = 20;
double densite = 0.3;
double rc = 2.5;
double h = 0.01;
double delta_r = 1; ////IL va jouer sur la fréquence de calcul des voisins donc ne doit pas etre trop petit ni trop grand
system1 sys(Nx,densite,rc,delta_r);
sys.init_system(1.0);
sys.integration_neighbour(h,100);
//sys.integration(h,100);
  return 0;
}
