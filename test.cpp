#include<iostream>
using namespace std;
#include"projet.h"

int main(){
  grid G = grid(100,100);
  cell C = G.get_cell(50,50);  
  C.display();
  cout << "test" << endl;
  return 0;

}
