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
  cout << "(" << this->x << "," << this->y << ")";
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

particle::particle()
{
  X = vector();
  V = vector();
  A = vector();
}

particle::particle(vector X_, vector V_, vector A_):
  X(X_),V(V_),A(A_)
{
}

void particle::display(){
  cout << "[";
  X.display();
  cout << ",";
  V.display();
  cout << ",";
  A.display();
  cout << "]";
}
//particule::~particule(){
//}

