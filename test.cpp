#include<iostream>
using namespace std;
#include"projet.h"

int main(){
  grid G(3,3);
  cell C = G.get_cell(0,0);
  C.add_particle(5);
  C.display();
  return 0;

}
