#include<iostream>
using namespace std;
#include"projet.h"
int main(){
int Nx = 20;
double densite = 0.3;
double rc = 2.5;
double h = 0.01;
system1 sys(Nx,densite,rc,0);
sys.init_system(1.0);
sys.integration_neighbour(h,10);
 //vect V(1,1);
 //vect X(3,4);
 // system1 S(10,1,1,1);
 // //S.init_system(1.0);
 // S.construct_neighbour_list();
 // S.compute_force_with_neighbour();
  //S.init_particle(2,X,V);
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
