//#include<iostream>
//using namespace std;
#include"projet.h"
int main(){
 vect V(1,1);
 vect X(3,4);
  system1 S(10,1,1,1);
  S.init_particle(2,X,V);
  cout << "list particle : ";
  S.list_particle[2].display();
  cout <<endl<< "----------"<<endl;
  int i = int(X.x/S.lc);
  int j = int(X.y/S.lc);
  cell C = S.Grid.get_cell(i,j);
  cout << "particles' index in cell : ";
  C.display();
  return 0;

}
