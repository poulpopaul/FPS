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
  //vect_name = new char[strlen(v.vect_name)];
  //strcpy(vect_name, v.vect_name);
}

vect& vect::operator=(const vect&v){
  x = v.x;
  y = v.y;
  return *this;
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
X(vect(0,0)),V(vect(0,0)),A(vect(0,0))
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
//strcpy(particle_name, p.particle_name);
}

particle& particle::operator=(const particle&p){
  X = p.X;
  V = p.V;
  A = p.A;
  return *this;
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
 delete particle_name;
}

cell::cell():
x(0.0), y(0.0), n(0)
{
}

cell::cell(double x_, double y_):
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
list_index_particle = c.list_index_particle;
}

cell& cell::operator=(const cell&c){
  x = c.x;
  
  y = c.y;
  n = c.n;
  list_index_particle = c.list_index_particle;
  return *this;
}

cell::~cell(){
  delete cell_name;
}

void cell::add_particle(int index){
  this->list_index_particle.push_back(index);
  this->n++;
}

void cell::remove_particle(int index){
  int i = 0;
  bool b = true;
  while(i<this->n && b){
    if (this->list_index_particle[i] == index)
    {
      if (this->n == 1)// || i == (this->n-1))
      {
        this->list_index_particle.pop_back();
        this->n--;
      }
      else {
        // cout << list_index_particle[i];
         // cout <<"   ";
        int a  = list_index_particle[this->n-1];
        this->list_index_particle[i] = a;
        this->list_index_particle.pop_back();
        this->list_index_particle = {};
        // cout << list_index_particle[i];
        //cout << endl;
      }
      this->n--;
      b = false;
    }
    else {i++;}
  }
}

void cell::display(){
  for (int i = 0; i < this->n;i++){
    cout << this->list_index_particle[i] << endl;
  }
}

grid::grid():
  L(1),Nc(1)
{
  tab = new cell*[1];
  tab[1] = new cell[1];
}

grid::grid(int L_, int N_):
  L(L_), Nc(N_)
{
 tab = new cell*[N_];
 for (int i = 0; i<N_;i++){
   tab[i] = new cell[N_];
 }
}

grid::grid(const grid &g): cell(g)
{
  // tab = new cell*[Nc];
  // for (int i = 0; i<Nc;i++){
  // tab[i] = new cell[Nc];
  //}
  L = g.L;
  Nc = g.Nc;
  tab = g.tab;
}

grid& grid::operator=(const grid&g){
  L = g.L;
  Nc = g.Nc;
  tab = g.tab;
  return *this;
}

grid::~grid(){
  //delete tab;
  delete grid_name;
}

void grid::remove(int i,int j,int index){
  int c = 0;
  bool b = true;
  while(i<this->n && b){
    if (this->tab[i][j].list_index_particle[c] == index)
    {
      if (this->n == 1)// || i == this->n - 1)
      {
        this->tab[i][j].list_index_particle.pop_back();
        this->tab[i][j].n--;
        this->tab[i][j].list_index_particle = {};
      }
      else {
        // cout << list_index_particle[i];
         // cout <<"   ";
        int a  = this->tab[i][j].list_index_particle[this->tab[i][j].n-1];
        this->tab[i][j].list_index_particle[c] = a;
        this->tab[i][j].list_index_particle.pop_back();
        this->tab[i][j].list_index_particle = {};
        //this->list_index_particle = {};
        // cout << list_index_particle[i];
        //cout << endl;
      }
      this->n--;
      b = false;
    }
    else {c++;}
  }
 // this->tab[i][j].remove_particle(index);
}

cell grid::get_cell(int i,int j){
int x_ = 0;
int y_ = 0;
if (i < 0){
  i += this->Nc;
  x_ = - this->L;
}
if (i >= this->Nc){
  i -= this->Nc;
  x_ = this->L;
}
if (j<0){
  j+= this->Nc;
  y_ = -this->L;
}
if (j>= this->Nc){
  j -= this->Nc;
  y_ = this->L;
}
this->tab[i][j].x = x_;
this->tab[i][j].y = y_;
return this->tab[i][j];
}

void grid::set_cell(int i,int j,int index){ //allow to update list_particle easily
int x_ = 0;
int y_ = 0;
while (i < 0){
  i += this->Nc;
  x_ -= this->L;
}
while (i >= this->Nc){
  i -= this->Nc;
  x_ += this->L;
}
while (j<0){
  j+= this->Nc;
  y_ -= this->L;
}
while (j>= this->Nc){
  j -= this->Nc;
  y_ += this->L;
}
this->tab[i][j].x = x_;
this->tab[i][j].y = y_;
tab[i][j].add_particle(index);
}

system1::system1():
Nx(1),density(1),rc(1),rc2(rc*rc),deltaR(0.1),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(PI/double(density))*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/double(rv))),lc(L/(double)Nc)
{
  // Grid = grid(Nc,L);
}

system1::system1(int Nx_, double density_, double rc_, double deltaR_):
Nx(Nx_),N(Nx_*Nx_),density(density_),rc(rc_),rc2(rc_*rc_),deltaR(deltaR_),rv(rc_+deltaR),rv2((rc_+deltaR)*(rc_+deltaR)),L(sqrt(PI/double(density))*0.5*Nx_),area(L*L),half_L(0.5*L),Nc(int(L/double(rv))),lc(L/(double)Nc)
{
  // Grid = grid(Nc,L);
  for (int i = 0; i < N; i++){
    list_particle.push_back(particle());
  }
}

system1::system1(const system1& s) : grid(s),Nx(s.Nx),N(Nx*Nx),density(s.density),rc(s.rc),rc2(rc*rc),deltaR(s.deltaR),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(PI/(double)density)*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/rv)),lc(L/(double)Nc),list_particle(s.list_particle)
{
  //strcpy(system1_name, s.system1_name);
  //  tab = new cell*[Nc];
  //for (int i = 0; i<Nc;i++){
  // tab[i] = new cell[Nc];
  //}
  //  Grid = s.Grid;
  //list_particle = s.list_particle;
  // Grid = grid(Nc,L);
  //for (int i = 0; i < N; i++){
  // list_particle.push_back(particle());
  //}
}

system1& system1::operator=(const system1&s){
 Grid = s.Grid;
 list_particle = s.list_particle;
  return *this;
}



system1::~system1(){
  delete system1_name;
}

void system1::init_particle(int index,vect X,vect V){
  this->list_particle[index].X = X;
  this->list_particle[index].V = V;
 // this->list_particle[index].A = vect(0,0);
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
  double px = 0.0;
  double py = 0.0;
   srand( (unsigned)time( NULL ) );
   for (int k = 0; k<this->N; k++){
      double a = (rand()/double(RAND_MAX))*PI*2.0;  
      double vx = velocity*cos(a);
      double vy = velocity*sin(a);
      px += vx;
      py += vy;
      init_particle(k,vect(x,y),vect(vx,vy));
      x += dx;
      if (x > this->L){
        x = dx*0.5;
        y += dy;
      }
    }
    for (int k = 0; k < this->N; k++){
      particle p = this->list_particle[k];
      vect S(px/double(this->N),py/double(this->N));
      this->list_particle[k].V = p.V - S;  
    }
   }
compute_force();
construct_neighbour_list();
}

void system1::move_particle(particle* p, int index, vect* X1){ // reference in order to really change p ! 
  //if (X1->x < 0){*X1 = *X1 + vect(this->L,0);}
  //if (X1->x > this->L){*X1 = *X1 - vect(this->L, 0);}
  //if (X1->y < 0){*X1 = *X1 + vect(0,this->L);}
  //if (X1->y > this->L){*X1 = *X1 - vect(0,this->L);}
  int i = int(p->X.x/double(this->lc));
  int j = int(p->X.y/double(this->lc));
  int i1 = int(X1->x/double(this->lc));
  int j1 = int(X1->y/double(this->lc));
  //cout << "(i,j) : ("<<i<<","<<j<<"),("<<i1<<","<<j1<<")"<<endl;
  if (i != i1 || j!=j1){
    int x_1 = 0;
    int y_1 = 0;
    if (i1 < 0){
      i1 += this->Nc;
      x_1 = - this->L;
    }
    if (i1 >= this->Nc){
       i1 -= this->Nc;
       x_1 = this->L;
     }
    if (j1<0){
      j1+= this->Nc;
      y_1 = -this->L;
    }
    if (j1>= this->Nc){
      j1 -= this->Nc;
      y_1 = this->L;
    }
    this->Grid.tab[i1][j1].x = x_1;
    this->Grid.tab[i1][j1].y = y_1;
    
   // cout << "size : " << Grid.get_cell(i,j).n;
   // cout << " "<<endl;;
   // this->Grid.tab[i][j].display();
   // this->Grid.tab[i][j].list_index_particle.pop_back();
   // cout <<"ok"<<endl;
    //this->Grid.set_cell(i1,j1,index);
  /////  AJOUTER LES CONDITIONS LIMITES ///////
  this->Grid.tab[i][j].remove_particle(index);
 ///this->Grid.remove(i,j,index);
 // this->Grid.set_cell(i1,j1,index);
  this->Grid.tab[i1][j1].list_index_particle.push_back(index);
 // this->Grid.tab[i1][j1].n += 1;

//remove(i,j,index);
 // cout << "move" << endl;
   //this->Grid.get_cell(i,j).remove_particle(index);
   //this->Grid.tab[i][j].display();
   // cout << " size 1: " << Grid.get_cell(i,j).n;
    //cout << endl;
  }
  p->X = *X1;
}

void system1::construct_neighbour_list(){
  this->list_neighbour = {};
  this->move_max = 0.0;
  for (int i = 0; i< this->Nc; i++){
    for (int j = 0; j< this->Nc; j++){
      //cell C = this->Grid.tab[i][j];
      for (int k = -1; k<2;k++){
        for (int l = -1; l<2;l++){
          int i1 = i+k;
          int j1 = j+l;
          int x_1 = 0;
          int y_1 = 0;
          if (i1 < 0){
            i1 += this->Nc;
            x_1 = - this->L;
          }
          if (i1 >= this->Nc){
            i1 -= this->Nc;
            x_1 = this->L;
          }
          if (j1<0){
            j1+= this->Nc;
            y_1 = -this->L;
          }
          if (j1>= this->Nc){
            j1 -= this->Nc;
            y_1 = this->L;
          }
          this->Grid.tab[i1][j1].x = x_1;
          this->Grid.tab[i1][j1].y = y_1;
          //cell C1 = this->Grid.tab[i+k][j+l];
          for (int index = 0; index<this->Grid.tab[i][j].list_index_particle.size(); index++){
            for (int index1 = 0; index1<this->Grid.tab[i1][j1].list_index_particle.size(); index1++){
              if (this->Grid.tab[i1][j1].list_index_particle[index1] < this->Grid.tab[i][j].list_index_particle[index]){
                //particle p = this->list_particle[index];
                //particle p1 = this->list_particle[index1];
                double dx = this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].X.x+ this->Grid.tab[i1][j1].x -this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].X.x;
                double dy = this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].X.y+ this->Grid.tab[i1][j1].y - this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].X.y;
                double r2 = dx*dx + dy*dy;
               if (r2 <= this->rv2){
                  this->list_neighbour.push_back(vect(this->Grid.tab[i][j].list_index_particle[index],this->Grid.tab[i1][j1].list_index_particle[index1]));
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
    this->list_particle[k].A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
  for (int i = 0; i<this->Nc;i++){
    for (int j = 0; j<this->Nc;j++){
      //cell c = this->Grid.tab[i][j];
      for (int k = -1;k<2;k++){
        for (int l = -1;l<2;l++){
          int i1 = i+k;
          int j1 = j+l;
          int x_ = 0;
          int y_ = 0;
          if (i1 < 0){
            i1 += this->Nc;
            x_ = - this->L;
          }
          if (i1 >= this->Nc){
            i1 -= this->Nc;
            x_ = this->L;
          }
          if (j1<0){
            j1+= this->Nc;
            y_ = -this->L;
          }
          if (j1>= this->Nc){
            j1 -= this->Nc;
            y_ = this->L;
          }
          this->Grid.tab[i1][j1].x = x_;
          this->Grid.tab[i1][j1].y = y_;
          //cell c1 = this->Grid.tab[i+k][j+l];
          for (int index=0;index < this->Grid.tab[i][j].list_index_particle.size();index++){
            for (int index1=0;index1 < this->Grid.tab[i1][j1].list_index_particle.size();index1++){
              if (this->Grid.tab[i1][j1].list_index_particle[index1] < this->Grid.tab[i][j].list_index_particle[index]){
               // cout << c1.list_index_particle[index1];
               // cout << c.list_index_particle[index];
                //particle p = this->list_particle[index];
                //particle p1 = this->list_particle[index1];
                double dx = this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].X.x+ this->Grid.tab[i1][j1].x -this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].X.x;
                double dy = this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].X.y+ this->Grid.tab[i1][j1].y - this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].X.y;
                double r2 = dx*dx + dy*dy;
                if (r2 < this->rc2){
                    double ir2 = 1.0/double(r2);
                    double ir6 = ir2*ir2*ir2;
                    double v = 24.0*ir6*(ir6-0.5);
                    double f = 2.0*v*ir2;
                    double fx = f*dx;
                    double fy = f*dy;
                    //cout << "fx :" << fx << endl;
                   // vect(fx,fy).display();
                   //this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].display();
                    this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].A.x = this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].A.x + fx;
                    this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].A.y = this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].A.y + fy;
                    this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].A.x = this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].A.x - fx;
                    this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].A.y = this->list_particle[this->Grid.tab[i][j].list_index_particle[index]].A.y - fy;
                   // this->list_particle[index].A = this->list_particle[index1].A - vect(fx,fy);
                  // this->list_particle[this->Grid.tab[i1][j1].list_index_particle[index1]].display();
                  // cout << endl;
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
    this->list_particle[k].A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
  int size = this->list_neighbour.size();
  for (int i = 0; i<size;i++){
    particle p = this->list_particle[list_neighbour[i].x];
    particle p1 = this->list_particle[list_neighbour[i].y];
    vect dX(p1.X.x-p.X.x,p1.X.y-p.X.y);
    while (dX.x >= this->half_L){
      dX.x = dX.x - this->L;
    }
    while (dX.x < -this->half_L){
      dX.x = dX.x + this->L;
    }
    while (dX.y >= this->half_L){
      dX.y = dX.y - this->L;
    }
    while (dX.y < -this->half_L){
      dX.y = dX.y + this->L;
    }
    double r2 = dX.x*dX.x + dX.y*dX.y;
    if (r2 <= this->rc2){
      double ir2 = 1.0/double(r2);
      double ir6 = ir2*ir2*ir2;
      double v = 24.0*ir6*(ir6-0.5);
      double f = 2.0*v*ir2;
      double fx = f*dX.x;
      double fy = f*dX.y;
      this->list_particle[list_neighbour[i].y].A.x = p1.A.x + fx;
      this->list_particle[list_neighbour[i].y].A.y = p1.A.y + fy;
      this->list_particle[list_neighbour[i].x].A.x = p1.A.x - fx;
      this->list_particle[list_neighbour[i].x].A.y = p1.A.y - fy;
     // this->list_particle[list_neighbour[i].x].A = p.A - vect(fx,fy);
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
    this->move_particle(&list_particle[k],k,&X1);
  }
  this->compute_force();
  for (int k =0; k<this->N;k++){
    particle p =this->list_particle[k];
    p.V = p.V + hd2*p.A;
  }
}


void system1::verlet_neighbour(double h, double hd2){
  for (int k =0; k<this->N;k++){
    this->list_particle[k].V = this->list_particle[k].V + hd2*this->list_particle[k].A;
    vect X1(this->list_particle[k].X.x + h*this->list_particle[k].V.x,this->list_particle[k].X.y+h*this->list_particle[k].V.y);
    this->move_particle(&list_particle[k],k,&X1);
  }
  this->compute_force_with_neighbour();
  double v2max = 0.0;
  for (int k =0; k<this->N;k++){
    this->list_particle[k].V = this->list_particle[k].V + hd2*this->list_particle[k].A;
    double v2 = this->list_particle[k].V.x * this->list_particle[k].V.x + this->list_particle[k].V.y*this->list_particle[k].V.y;
    if (v2 > v2max){
      v2max = v2;
    }
    this->move_max += sqrt(v2max)*h;
    if (this->move_max*2.0 > this->deltaR){
      this->construct_neighbour_list();
      this->move_max = 0.0;
    }
  }
}


void system1::integration(double h,int n){
  double hd2 = h/2.0;
  for (int i = 0; i<n;i++){
    this->verlet(h,hd2);
  }
}

void system1::integration_neighbour(double h,int n){
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
    this->verlet(h,hd2);
   // this->list_particle[0].display();
    //cout << endl;
    //cout << "&&&&&&&&&&&&&&&&&&&&&&" << endl;
  }
}
