#include<iostream>
#include <vector>
#include<cstring>
#include<cmath>
#include<time.h>
#include <stdlib.h>
#include<fstream>
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
  //strcpy(vect_name, v.vect_name);
}

void vect::display(){
  cout << "(" << this->x << "," << this->y << ")";
}

vect::~vect()
{
  //delete vect_name;
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
particle_name = new char[strlen(p.particle_name)];
strcpy(particle_name, p.particle_name);
}

void particle::set_x(double x_)
{
  X.x = x_;
}

void particle::set_y(double y_)
{
  X.y = y_;
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
  //delete particle_name;
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
 // delete cell_name;
 //list_particle.clear();
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

grid::~grid(){}

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

void grid::set_cell(int i,int j,int index){ //allow to update list_particle easily
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
tab[i][j].add_particle(index);
}

system1::system1():
Nx(1),density(1),rc(1),deltaR(0.1),rc2(rc*rc),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(PI/double(density))*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/double(rv))),lc(L/(double)Nc)
{
 Grid = grid(Nc,L);
}

system1::system1(int Nx_, double density_, double rc_, double deltaR_):
Nx(Nx_),N(Nx_*Nx_),density(density_),rc(rc_),rc2(rc_*rc_),deltaR(deltaR_),rv(rc_+deltaR),rv2((rc_+deltaR)*(rc_+deltaR)),L(sqrt(PI/double(density))*0.5*Nx_),area(L*L),half_L(0.5*L),Nc(int(L/double(rv))),lc(L/(double)Nc)
{
  Grid = grid(Nc,L);
  for (int i = 0; i < N; i++){
    list_particle.push_back(particle());
  }
}

system1::system1(const system1& s) : Nx(s.Nx),N(Nx*Nx),density(s.density),rc(s.rc),rc2(rc*rc),deltaR(s.deltaR),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(PI/(double)density)*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/rv)),lc(L/(double)Nc)
{
  strcpy(system1_name, s.system1_name);
  Grid = s.Grid;
  list_particle = s.list_particle;
}

system1::~system1(){
  delete system1_name;
}

void system1::init_particle(int index,vect X,vect V){
  this->list_particle[index].X = X;
  this->list_particle[index].V = V;
  double i = int(X.x/double(this->lc));
  double j = int(X.y/double(this->lc));
  this -> Grid.set_cell(int(i),int(j),index); //remplace get_cell + add_particle !
}

void system1::init_system(double velocity){
  double dx = this->L/double(this->Nx);
  cout << "initial distance : " << dx << endl;
   if (dx <= this->diameter){     throw "Density: too high!";//cout << "too high density" << endl;
   }
   else{
  double dy = dx;
  double x = dx*0.5;
  double y = x;
  //cout << "y avant : " << y << endl;
  double px = 0.0;
   double py = 0.0;
   srand( (unsigned)time( NULL ) );
   for (int k = 0; k<this->N; k++){
  //srand( (unsigned)time( NULL ) );
   //  cout << k << "/" << this->N << endl;
     
     double a = (rand()/double(RAND_MAX))*PI*2.0;  
    double vx = velocity*cos(a);
    double vy = velocity*sin(a);
    px += vx;
    py += vy;
    //cout << "a : " << a<< endl;
    //cout << "vy :" << vy << endl;
    vect X1(x,y);
    //X1.display();
    vect V1(vx,vy);
    init_particle(k,X1,V1);
    x += dx;
   // cout << "y apres :" << y<<endl;
     if (x>this->L){
      x = 0;//x*0.5;
       y += dy;
      // cout << "y 2apres :" << y << endl;
    }
  }
    for (int k = 0; k < this->L; k++){
      //cout << k << "/" << this->L << endl;
      particle p = this->list_particle[k];
      vect S(px/double(this->N),py/double(this->N));
      p.V = p.V - S;  
    }
  compute_force();
  construct_neighbour_list();
   }
}

void system1::move_particle(particle* p, int index, vect X1){ // reference in order to really change p ! 
  if (X1.x < 0){X1 = X1 + vect(this->L,0);}
  if (X1.x > this->L){X1 = X1 - vect(this->L, 0);}
  if (X1.y < 0){X1 = X1 + vect(0,this->L);}
  if (X1.y > this->L){X1 = X1 - vect(0,this->L);}
  int i = int(p->X.x/double(this->lc));
  int j = int(p->X.y/double(this->lc));
  int i1 = int(X1.x/double(this->lc));
 int j1 = int(X1.y/double(this->lc));
  if (i != i1 || j!=j1){
    cell C = this->Grid.get_cell(i,j);
    C.remove_particle(index);
    cell C1 = this->Grid.get_cell(i1,j1);
    C1.add_particle(index);
  }
  //p->display();

  //X.display();
  p->X = X1;

  //X1.display();
  //p->display();
  //cout << endl;
  //p->display();
  //cout << endl;
}

void system1::construct_neighbour_list(){
  vector<vect> list_neighbour = {};
  double move_max = 0.0;
  int Nc = this->Nc;
  for (int i = 0; i< Nc; i++){
    for (int j = 0; j< Nc; j++){
      cell C = this->Grid.get_cell(i,j);
      for (int k = -1; k<2;k++){
        for (int l = -1; l<2;l++){
          cell C1 = this->Grid.get_cell(i+k,j+l);
          for (int index = 0; index<C.n; index++){
            for (int index1 = 0; index1<C1.n; index1++){
              if (C1.list_particle[index1] < C.list_particle[index]){
                particle p = this->list_particle[index];
                particle p1 = this->list_particle[index1];
                double dx = p1.X.x+ C1.x -p.X.x;
                double dy = p1.X.y+ C1.y - p.X.y;
                double r2 = dx*dx + dy*dy;
                if (r2 <= this->rv2){
                  list_neighbour.push_back(vect(index1,index));
                }
              }
            }
          }
        }

      }
    }
  }
}

void system1::compute_force(){
  for (int k = 0; k<this->N; k++){
    particle p = this->list_particle[k];
    p.A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
  for (int i = 0; i<this->Nc;i++){
    //cout << i+1 << "/" << Nc <<endl;
    for (int j = 0; j<this->Nc;j++){
      //cout << "|||" << j+1 <<"/" << Nc << endl;
      cell c = this->Grid.get_cell(i,j);
      for (int k = -1;k<2;k++){
        for (int l = -1;l<2;l++){
          cell c1 = this->Grid.get_cell(i+k,j+l);
          for (int index=0;index < c.n;index++){
            for (int index1=0;index1 < c1.n;index1++){
              if (c1.list_particle[index1] < c.list_particle[index]){
                particle p = this->list_particle[index];
                particle p1 = this->list_particle[index1];
                double dx = p1.X.x+ c.x -p.X.x;
                double dy = p1.X.y+ c.y - p.X.y;
                double r2 = dx*dx + dy*dy;
                if (r2 <= this->rv2){
                    double ir2 = 1.0/double(r2);
                    double ir6 = ir2*ir2*ir2;
                    double v = 24.0*ir6*(ir6-0.5);
                    double f = 2.0*v*ir2;
                    double fx = f*dx;
                    double fy = f*dy;
                    p1.A = p1.A + vect(fx,fy);
                    p.A = p.A - vect(fx,fy);
                    this->E_pot += 4.0*ir6*(ir6-1.0);
                    this->viriel +=v;
                }
              }
            }
          } 
        }
      }
    }
  }
}

void system1::compute_force_with_neighbour(){
  for (int k = 0; k<this->N; k++){
    particle p = this->list_particle[k];
    p.A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
  for (int i = 0; i<this->list_neighbour.size();i++){
    particle p = this->list_particle[list_neighbour[i].x];
    particle p1 = this->list_particle[list_neighbour[i].y];
    vect dX(p1.X.x-p.X.x,p1.X.y-p.X.y);
    if (dX.x >= this->half_L){
      dX = dX - vect(this->L,0);
    }
    if (dX.x < -this->half_L){
      dX = dX + vect(this->L,0);
    }
    if (dX.y >= this->half_L){
      dX = dX - vect(0,this->L);
    }
    if (dX.y < -this->half_L){
      dX = dX + vect(0,this->L);
    }
    double r2 = dX.x*dX.x + dX.y*dX.y;
    if (r2 <= this->rv2){
      double ir2 = 1.0/double(r2);
      double ir6 = ir2*ir2*ir2;
      double v = 24.0*ir6*(ir6-0.5);
      double f = 2.0*v*ir2;
      double fx = f*dX.x;
      double fy = f*dX.y;
      p1.A = p1.A + vect(fx,fy);
      p.A = p.A - vect(fx,fy);
      this->E_pot += 4.0*ir6*(ir6-1.0);
      this->viriel +=v;
    }
  }
}

void system1::compute_E_kin(){
  this->E_kin = 0.0;
  for (int i = 0; i<this->N;i++){
    particle p = this->list_particle[i];
    this->E_kin += 0.5*(p.V.x*p.V.x + p.V.y*p.V.y);
  }
  this->pressure = (this->viriel + this->E_kin)/double(this->area);
  this->E_kin /= double(this->N);
  this->energy = this->E_kin+this->E_pot/double(this->N);
  this->sum_temp += this->E_kin;
  this->sum_temp2 += this->E_kin * this->E_kin;
  this->sum_pressure += this->pressure;
  this->sum_pressure2 += this->pressure * this->pressure;
  this->counter += 1;
}

void system1::init_mean(){
  this->counter = 0.0;
  this->sum_temp = 0.0;
  this->sum_pressure = 0.0;
  this->sum_temp2 = 0.0;
  this->sum_pressure2 = 0.0;
}


void system1::mean_temp(vect* v){
  double Tm = this->sum_temp/double(this->counter);
  v->x = Tm;
  v->y = sqrt(this->sum_temp2/double(this->counter) - Tm*Tm);
}

void system1::mean_pressure(vect* v){
  double Pm= this->sum_pressure/double(this->counter);
  v->x = Pm;
  v->y = sqrt(this->sum_pressure2/double(this->counter) - Pm*Pm);
}

void system1::adjust_v(double T){
  double Tm = this->sum_temp/double(this->counter);
  double f = sqrt(T/double(Tm));
  for (int k = 0; k<this->N;k++){
    particle p = this->list_particle[k];
    p.V = p.V*f;
  }
}

void system1::verlet(double h, double hd2){
  for (int k =0; k<this->N;k++){
    particle p =this->list_particle[k];
    p.V = p.V + hd2*p.A;
    vect X1(p.X.x + h*p.V.x,p.X.y+h*p.V.y);
    this->move_particle(&p,k,X1);
  }
  this->compute_force();
  for (int k =0; k<this->N;k++){
    particle p =this->list_particle[k];
    p.V = p.V + hd2*p.A;
  }
}


void system1::verlet_neighbour(double h, double hd2){
  for (int k =0; k<this->N;k++){
    particle p =this->list_particle[k];
    p.V = p.V + hd2*p.A;
    vect X1(p.X.x + h*p.V.x,p.X.y+h*p.V.y);
    this->move_particle(&list_particle[k],k,vect(p.X.x + h*p.V.x,p.X.y+h*p.V.y));
  }
  this->compute_force_with_neighbour();
  double v2max = 0.0;
  for (int k =0; k<this->N;k++){
    particle p =this->list_particle[k];
    p.V = p.V + hd2*p.A;
    double v2 = p.V.x * p.V.x + p.V.y*p.V.y;
    if (v2 < v2max){
      v2max = v2;
    }
    this->move_max += sqrt(v2max)*h;
    if (this->move_max*2.0 > this->deltaR){
      this->construct_neighbour_list();
      this->move_max = 0.0;
    }
  }
}


void system1::integration(double h,unsigned int n){
  double hd2 = h/2.0;
  for (int i = 0; i<n;i++){
    this->verlet(h,hd2);
  }
}

void system1::integration_neighbour(double h,unsigned int n){
  int taille = this->list_particle.size();
  fstream fichx, fichy;
  fichx.open("x.txt", ios::out);
  fichy.open("y.txt", ios::out);
  fichx << taille << endl;
  fichy << taille << endl;
  double hd2 = h*0.5;
  for (int i = 0; i<n;i++){
    for (int k = 0; k<taille; k++){
  //sys.list_particle[k].display();
      fichx << this->list_particle[k].X.x << endl;
      fichy << this->list_particle[k].X.y << endl;
    }
    //cout << i << "/" << n<< endl;
    this->verlet_neighbour(h,hd2);
  }
}
