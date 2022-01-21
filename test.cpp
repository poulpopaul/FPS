#include<iostream>
#include"time.h"
using namespace std;
#include"projet.h"
#include<fstream>
int main(){
//srand( (unsigned)time( NULL ) );
int Nx = 20;
double densite = 0.3;
double rc = 2.5;
double h = 0.01;
system1 sys(Nx,densite,rc,5);
sys.init_system(1);
sys.integration_neighbour(h,1000);

///int n = sys.list_particle.size();
///fstream fichx, fichy;
///fichx.open("x.txt", ios::out);
///fichy.open("y.txt", ios::out);
 // system1 S(10,1,1,1);
 // //S.init_system(1.0);
 // S.construct_neighbour_list();
 // S.compute_force_with_neighbour();
  ///sys.init_particle(0,X,V);
  ///sys.list_particle[0].A = A;
  ///sys.verlet_neighbour(h,h*0.5);
  //sys.init_particle(1,X1,V1);
  //sys.move_particle(&sys.list_particle[0],0,vect(2,2));
///for (int k = 0; k<n; k++){
///  sys.list_particle[k].display();
///  fichx << sys.list_particle[k].X.x << endl;
///  fichy << sys.list_particle[k].X.y << endl;
///}
///sys.list_particle[0].display();
///cout << "ecriture dans les fichiers terminee" <<endl;
  //cout << "list particle : ";
  //S.list_particle[2].display();
  //cout <<endl<< "----------"<<endl;
  //int i = int(X.x/S.lc);
  //int j = int(X.y/S.lc);
  //cell C = S.Grid.get_cell(i,j);
  //cout << "particles' index in cell : ";
  //C.display();
  //particle* p = &(S.list_particle[2]); // reference !
  //int index = 2;
  //vect X1(5,6);
  //S.move_particle(p, index, X1);
  //(*p).display();
  //S.list_particle[2].display();
  //cout << "ok !" <<endl;
  return 0;
}
