#include<iostream>
#include"time.h"
using namespace std;
#include"projet.h"
#include<fstream>
int main(){
int Nx = 20;
double densite = 0.3;
double rc = 5;
double h = 0.01;
system1 sys(Nx,densite,rc,1);
sys.init_system(1.0);
cout << "end init";
sys.integration_neighbour(h,1000);
  return 0;
}
