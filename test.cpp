#include<iostream>
using namespace std;
#include"projet.h"

int main(){
  vector B(2,3);
  vector A;//bite
  vector C(1,10);
  B.display();
  A.display();
  vector D = A+B+C;
  D.display();
  vector E = B*C;
  E.display();
  vector F = A-B+C;
  F.display();
  vector G = 5*B;
  G.display();
  return 0;

}
