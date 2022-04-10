#include<iostream>
#include"time.h"
using namespace std;
#include"projet.h"
#include<fstream>
int main(){

int Nx = 40;
double density = 0.3;
double rc = 2.5;
double h = 0.00001;
double delta_r = 0.5;
double size = 0;

system1 sys(Nx,density,rc,delta_r);
sys.init_system(5.,size);
int taille = sys.list_particle.size();

fstream fichx, fichy, ficht, fiche, fichep, fichec;
fichx.open("x.txt", ios::out);
fichy.open("y.txt", ios::out);
ficht.open("t.txt", ios::out);
fiche.open("e.txt", ios::out);
fichec.open("ec.txt", ios::out);
fichep.open("ep.txt", ios::out);
fichx << taille << endl;
fichy << taille << endl;

clock_t begin = clock();

//cout << "Atteinte de l'équilibre...";
int N_it = 10000;
/*int p = 0;
for (int c = 0; c<N_it; c++){
for (int k = 0; k<taille; k++){
    fichx << sys.list_particle[k].X.x << endl;
    fichy << sys.list_particle[k].X.y << endl;
    ficht << sys.list_particle[k].type << endl;
  }
sys.compute_E_kin();
fichep<< sys.E_pot << endl;
fichec<< sys.E_kin << endl;
fiche<< sys.energy << endl;
sys.integration(h,1,false);
if(c%int((N_it/10.)) == 0){
system("clear");

cout << "Atteinte de l'équilibre..."<<endl;
for (int i = 0; i<p*0.1; i++){
  cout << "##";
}
for (int i = 0; i<(100-p)*0.1; i++){
  cout << "  ";
}
cout <<" ||  "<<p<<"%"<<endl;
p+= 10;
}
}
*/
  int p = 0;
for (int c = 0; c<N_it; c++){
for (int k = 0; k<taille; k++){
    fichx << sys.list_particle[k].X.x << endl;
    fichy << sys.list_particle[k].X.y << endl;
    ficht << sys.list_particle[k].type << endl;
  }
sys.compute_E_kin();
fiche<< sys.energy << endl;
fichec<< sys.E_kin << endl;
fichep<< sys.E_pot << endl;
sys.integration(h,10,true);

if(c%int((N_it/10.)) == 0 || c == N_it -1){
system("clear");

cout << "Simulation en cours"<<endl;
for (int i = 0; i<p*0.1; i++){
  cout << "##";
}
for (int i = 0; i<(100-p)*0.1; i++){
  cout << "  ";
}
cout <<" ||  "<<p<<"%"<<endl;
p+= 10;
}



}

clock_t end = clock();
unsigned long sec = (end -  begin) / CLOCKS_PER_SEC;
cout << "Temps de simulation: "<<sec<<"s" << endl;
//cout << "Time with neighbours: " << millis<<" ms"<<endl;
//for (int k = 0; k<sys.list_particle.size();k++){
  //cout << sys.list_particle[k].X.x << "/" << sys.list_particle[k].X.y << endl;
//}
//cout << sys.L << endl;
//system("python3 anim.py");
  return 0;
}
