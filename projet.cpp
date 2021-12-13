#include<iostream>
#include <vector>
using namespace std;
#include"projet.h"

vect::vect():
  x(0),y(0)
{
}

vect::vect(double x_,double y_):
  x(x_),y(y_)
{
}

void vect::display(){
  cout << "(" << this->x << "," << this->y << ")";
}

//vect::~vect()
//{
//}

vect operator+(vect v1, vect v2){
  return vect(v1.x + v2.x, v1.y + v2.y);
}

vect operator-(vect v1, vect v2){
  return vect(v1.x - v2.x, v1.y - v2.y);
}

vect operator*(vect v1, vect v2){
  return vect(v1.x * v2.x, v1.y * v2.y);
}

vect operator*(vect v1, double l){
  return vect(v1.x * l, v1.y * l);
}

vect operator*(double l, vect v1){
  return vect(v1.x * l, v1.y * l);
}

particle::particle()
{
  X = vect();
  V = vect();
  A = vect();
}

particle::particle(vect X_, vect V_, vect A_):
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

cell::cell():
x(0), y(0), n(0)
{
}

cell::cell(unsigned int x_, unsigned int y_):
x(x_), y(y_), n(0)
{
}

void cell::add_particle(unsigned int index){
  list_particle.push_back(index);
  n++;
}

void cell::remove_particle(unsigned int index){
  unsigned int i = 0;
  bool b = true;
  while(i<this->n && b){
    if (this->list_particle[i] == index)
    {
      if (this->n == 1 || i == this->n - 1)
      {
        list_particle.pop_back();
        this->n--;
      }
      else {
        list_particle[i] = list_particle[this->n - 1];
        list_particle.pop_back();
        this->n--;
      }
    }
    else {i++;}
  }
}

void cell::display(){
  for (unsigned int i = 0; i < n;i++){
    cout << list_particle[i] << endl;
  }
}

grid::grid():
  L(0),N(0)
{
}

grid::grid(int L_, int N_):
  L(L_), N(N_)
{
 cell **tab = new cell*[N_];
 for(int i = 0; i<N_; i++){
   tab[i] = new cell[N_];
   for(int j = 0; j<N_; j++){
     tab[i][j] = cell();
   }
 }
}

cell grid::get_cell(int i, int j){
int x = 0;
int y = 0;
if (i < 0){
  i += this->N;
  x = - this->L;
}
if (i >= this->N){
  i -= this->N;
  x = this->L;
}
if (j<0){
  j+= this->N;
  y = -this->N;
}
if (j>= this->N){
  j -= this->N;
  y = this->L;
}

cell cell_ = this->tab[i][j];
cell_.x = x;
cell_.y = y;
return cell_;
}
