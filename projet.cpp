#include<iostream>
using namespace std;
#include"projet.h"

vector::vector():
  x(0),y(0)
{
}

vector::vector(double x_,double y_):
  x(x_),y(y_)
{
}

void vector::display(){
  cout << "(" << this->x << "," << this->y << ")" << endl;
}

//vector::~vector()
//{
//}

vector operator+(vector v1, vector v2){
  return vector(v1.x + v2.x, v1.y + v2.y);
}

vector operator-(vector v1, vector v2){
  return vector(v1.x - v2.x, v1.y - v2.y);
}

vector operator*(vector v1, vector v2){
  return vector(v1.x * v2.x, v1.y * v2.y);
}

vector operator*(vector v1, double l){
  return vector(v1.x * l, v1.y * l);
}

vector operator*(double l, vector v1){
  return vector(v1.x * l, v1.y * l);
}

particule::particule()
{
  X = vector();
  V = vector();
  A = vector();
}

particule::particule(vector X_, vector V_, vector A_):
  X(X_),V(V_),A(A_)
{
}

//particule::~particule(){
//}

