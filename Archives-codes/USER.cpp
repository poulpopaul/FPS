#include<iostream>
#include"time.h"
using namespace std;
#include"projet.h"
#include<fstream>

int main(){

/* System initialization
##########################################
*/

int Nx = 20;
double h = 1E-4;
double rc = 2.5;
double delta_r = 0.5;
double v0 = 1.;
double size = 0;

system1 sys(Nx,0.3,rc,delta_r);
sys.mag = false;
sys.cutoff = 1E3;
sys.init_system(v0,size);


/* Creation of the files
##########################################
*/

int nb_part = sys.list_particle.size();
fstream fichx, fichy, ficht, fiche, fichep, fichec;
fichx.open("x.txt", ios::out);
fichy.open("y.txt", ios::out);
ficht.open("t.txt", ios::out);
//fiche.open("e.txt", ios::out);
//fichec.open("ec.txt", ios::out);
//fichep.open("ep.txt", ios::out);
fichx << nb_part << endl;
fichy << nb_part << endl;


/* Simulation and file filling
##########################################
*/


int N_it1 = 100; //main loop (file filling). Please N_it1 >=20 for % display
int N_it2 = 100;
int N_it_tot = N_it1 * N_it2; 

clock_t begin = clock();

int p = 0;
for (int c = 0; c<N_it1; c++){
  for (int k = 0; k<nb_part; k++){
    fichx << sys.list_particle[k].X.x << endl;
    fichy << sys.list_particle[k].X.y << endl;
    ficht << sys.list_particle[k].type << endl;
  }

sys.compute_E_kin();
//fiche<< sys.energy << endl;
//fichec<< sys.E_kin << endl;
//fichep<< sys.E_pot << endl;


sys.integration(h,N_it2);


/* Graphic interface
##########################################
*/

if(c%int((N_it1/20.)) == 0 || c == N_it1 -1){
  system("clear");
  cout << "Simulation in progress..."<<endl<<endl;
  for (int i = 0; i<p*0.2; i++){
    cout << "#";
  }
  for (int i = 0; i<(100-p)*0.2; i++){
    cout << ".";
  }
  cout <<" ||  "<<p<<"%"<<endl<<endl;
  p+= 5;
  }
}

clock_t end = clock();

unsigned long sec = (end -  begin) / CLOCKS_PER_SEC;

cout << "*********************************" <<endl;
cout<<"   ***************************"<<endl<<endl;
cout << "Total number of iterations: " << N_it_tot << endl<<endl;
if (sec <1){unsigned long millis = 1000*(end -  begin) / CLOCKS_PER_SEC;
cout << "Computation time: "<<millis<<" ms" << endl<<endl;
}
else{cout << "Computation time: "<<sec<<" s" << endl<<endl;}
cout<<"   ***************************"<<endl;
cout << "*********************************" <<endl<<endl;



system("python3 anim.py"); 

/* Creates and save a GIF animation: By default the number of frames is N_it1 
(be careful, you have to lower it if N_it1 becomes big, otherwise  the creation
will be very long and the GIF very heavy). A maximum of 1000 frames is ok (45 MB). */

return 0;
}
