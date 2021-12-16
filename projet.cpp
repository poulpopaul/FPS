#include<iostream>
#include <vector>
#include<cstring>
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

vect::vect(const vect & v){
  x = v.x;
  y = v.y;
 // vect_name = new char[strlen(v.vect_name)];
  strcpy(vect_name, v.vect_name);
}

void vect::display(){
  cout << "(" << this->x << "," << this->y << ")";
}

vect::~vect()
{
  delete vect_name;
}

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

particle::particle():
X(vect()),V(vect()),A(vect())
{
}

particle::particle(vect X_, vect V_, vect A_):
  X(X_),V(V_),A(A_)
{
}

particle::particle(const particle &p): vect(p)
{
X = p.X;
V = p.V;
A = p.A;
//particle_name = new char[strlen(p.particle_name)];
strcpy(particle_name, p.particle_name);
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
particle::~particle(){
  delete particle_name;
}

cell::cell():
x(0), y(0), n(0)
{
}

cell::cell(unsigned int x_, unsigned int y_):
x(x_), y(y_), n(0)
{
}


cell::cell(const cell &c): particle(c)
{
//cell_name = new char[strlen(c.cell_name)];
strcpy(cell_name, c.cell_name);
x = c.x;
y = c.y;
n = c.n;
list_particle = c.list_particle;
}

cell::~cell(){
  delete cell_name;
 list_particle.clear();
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
  L(1),N(1)
{
  tab = new cell*[1];
  tab[1] = new cell[1];
}

grid::grid(int L_, int N_):
  L(L_), N(N_)
{
 tab = new cell*[N_];
 for (int i = 0; i<N_;i++){
   tab[i] = new cell[N_];
 }
}

grid::grid(const grid &g): cell(g)
{
  L = g.L;
  N = g.N;
  tab = g.tab;
}

cell grid::get_cell(int i,int j){
int x_ = 0;
int y_ = 0;
if (i < 0){
  i += this->N;
  x_ = - this->L;
}
if (i >= this->N){
  i -= this->N;
  x_ = this->L;
}
if (j<0){
  j+= this->N;
  y_ = -this->N;
}
if (j>= this->N){
  j -= this->N;
  y_ = this->L;
}
(this->tab[i][j]).x = x_;
(this->tab[i][j]).y = y_;
return this->tab[i][j];
}