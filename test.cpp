#include<iostream>
using namespace std;
#include"projet.h"

int main(){
  grid G(3,3);
  //cell *cell_ = G.tab[0];
  //cell C = cell(0,0);
  //C.add_particle(8);
  //cell_[0] = C;

  //cell_[0].display();
  //G.tab[0][0].display();
  cell C = G.get_cell(0,0);
  //C.add_particle(5);
 //C.display();
  cout << "test" << endl;
  return 0;

}
