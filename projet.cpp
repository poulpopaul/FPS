#include <vector>
#include<cstring>
#include<math.h>
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
}

vect& vect::operator=(const vect&v){
  x = v.x;
  y = v.y;
  return *this;
}

vect::~vect()
{
  delete vect_name;
}



particle::particle():
X(vect(0,0)),V(vect(0,0)),A(vect(0,0)),type(1)
{
}

particle::particle(vect X_, vect V_, vect A_,double type_):
  X(X_),V(V_),A(A_),type(type_)
{
}

particle::particle(const particle &p): vect(p)
{
X = p.X;
V = p.V;
A = p.A;
type = p.type;
}

particle& particle::operator=(const particle&p){
  X = p.X;
  V = p.V;
  A = p.A;
  type = p.type;
  return *this;
}

particle::~particle(){
 delete particle_name;
}

cell::cell():
move_x(0.0), move_y(0.0), n(0)
{
}

cell::cell(double x_, double y_):
move_x(x_), move_y(y_), n(0)
{
}


cell::cell(const cell &c): particle(c)
{
move_x = c.move_x;
move_y = c.move_y;
n = c.n;
list_index_particle = c.list_index_particle;
}

cell& cell::operator=(const cell&c){
  move_x = c.move_x;
  move_y = c.move_y;
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
  // remove the index in the list_index_particle: the particle index is no longer in the cell
  int i = 0;
  bool b = true;
  while(i<this->n && b){
    if (this->list_index_particle[i] == index)
    {
      if (this->n == 1)
      {
        this->list_index_particle.pop_back();
        this->n--;
      }
      else {
        int a  = list_index_particle[this->n-1];
        this->list_index_particle[i] = a;
        this->list_index_particle.pop_back();
      }
      this->n--;
      b = false;
    }
    else {i++;}
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
  delete tab;
  delete grid_name;
}

void grid::get_cell(int* i,int* j){
int x_ = 0;
int y_ = 0;
if (*i < 0){
  *i += this->Nc;
  x_ = - this->L;
}
if (*i >= this->Nc){
  *i -= this->Nc;
  x_ = this->L;
}
if (*j<0){
  *j+= this->Nc;
  y_ = -this->L;
}
if (*j>= this->Nc){
  *j -= this->Nc;
  y_ = this->L;
}
this->tab[*i][*j].move_x = x_;
this->tab[*i][*j].move_y = y_;
}

void grid::set_cell(int i,int j,int index){ //allow to update list_particle
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
this->tab[i][j].move_x = x_;
this->tab[i][j].move_y = y_;
tab[i][j].add_particle(index);
}

system1::system1():
Nx(1),density(1),rc(1),rc2(rc*rc),deltaR(0.1),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(M_PI/double(density))*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/double(rv))),lc(L/(double)Nc)
{
}

system1::system1(int Nx_, double density_, double rc_, double deltaR_):
Nx(Nx_),N(Nx_*Nx_),density(density_),rc(rc_),rc2(rc_*rc_),deltaR(deltaR_),rv(rc_+deltaR),rv2((rc_+deltaR)*(rc_+deltaR)),L(sqrt(M_PI/double(density))*0.5*Nx_),area(L*L),half_L(0.5*L),Nc(int(L/double(rv))),lc(L/(double)Nc)
{
  for (int i = 0; i < N; i++){
    list_particle.push_back(particle());
  }
}

system1::system1(const system1& s) : grid(s),Nx(s.Nx),N(Nx*Nx),density(s.density),rc(s.rc),rc2(rc*rc),deltaR(s.deltaR),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(M_PI/(double)density)*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/rv)),lc(L/(double)Nc),list_particle(s.list_particle)
{
  Grid = s.Grid;
  list_particle = s.list_particle;
}

system1& system1::operator=(const system1&s){
 Grid = s.Grid;
 list_particle = s.list_particle;
  return *this;
}


system1::~system1(){
  delete system1_name;
}

void system1::init_particle(int index,vect X,vect V,double type_){
  this->list_particle[index].type = type_;
  this->list_particle[index].X = X;
  this->list_particle[index].V = V;
  double i = int(X.x/double(this->lc));
  double j = int(X.y/double(this->lc));
  this -> Grid.set_cell(int(i),int(j),index);
}

void system1::init_system(double velocity){
  double dx = this->L/double(this->Nx);
  cout << "initial distance : " << dx << endl;
   //if (dx <= diameter){   cout << "Density too high !"<<endl;
   //throw "Density: too high!";
   //}
  // else{
  double dy = dx;
  double x = dx*0.5;
  double y = x;
  double px = 0.0;
  double py = 0.0;
   srand( (unsigned)time( NULL ) );
   for (int k = 0; k<this->N; k++){
     double t = 5; 
      double a = (rand()/double(RAND_MAX))*M_PI*2.0;  
      double vx = velocity*cos(a);
      double vy = velocity*sin(a);
      px += vx;
      py += vy;
      if ((x - 0.5*L)*(x - 0.5*L) + (y - 0.5*L)*(y - 0.5*L) <= L*2){t = 0.1;}
      init_particle(k,vect(x,y),vect(vx,vy),t);
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
  // }
compute_force();
construct_neighbour_list();
}

void system1::move_particle(particle* p, int index, vect* X1){
  if (X1->x < 0){X1->x += L;}
  if (X1->x > this->L){X1->x -= L;}
  if (X1->y < 0){X1->y += L;}
  if (X1->y > this->L){X1->y -= L;}
  int i = int(p->X.x/double(this->lc));
  int j = int(p->X.y/double(this->lc));
  int i1 = int(X1->x/double(this->lc));
  int j1 = int(X1->y/double(this->lc));
  if ((i != i1 || j!=j1) && i1>=0 && j1 >= 0 && i1 < this->Nc && j1< this->Nc){
    this->Grid.tab[i1][j1].move_x = X1->x;
    this->Grid.tab[i1][j1].move_y = X1->y;
    this->Grid.tab[i][j].remove_particle(index);
    this->Grid.tab[i1][j1].list_index_particle.push_back(index);
  }
  p->X = *X1;
}

void system1::construct_neighbour_list(){
 // clock_t begin = clock();
  this->list_neighbour = {};
  int* i1 = new int;
  int* j1 = new int;
  for (int i = 0; i< this->Nc; i++){
    for (int j = 0; j< this->Nc; j++){
      cell *c;
      c = &this->Grid.tab[i][j];
      for (int k = -1; k<2;k++){
        for (int l = -1; l<2;l++){
          *i1 = i+k;
          *j1 = j+l;
          this->Grid.get_cell(i1,j1);
          cell *c1;
          c1 = &this->Grid.tab[*i1][*j1];
          for (int id0 : c->list_index_particle){
            for (int id1 : c1->list_index_particle){
              if (id1 < id0){
                double dx = this->list_particle[id1].X.x+ c1->move_x -this->list_particle[id0].X.x;
                double dy = this->list_particle[id1].X.y+ c1->move_y - this->list_particle[id0].X.y;
                double r2 = dx*dx + dy*dy;
               if (r2 <= this->rv2){
                  this->list_neighbour.push_back(vect(id0,id1));
                }
              }
            }
          }
        }
      }
    }
  }
delete i1;
delete j1;
//clock_t end = clock();
//unsigned long millis = (end -  begin) * 1000 / CLOCKS_PER_SEC;
//cout << "cn: " << millis;
} 

void system1::compute_force(){
  for (int k = 0; k<this->N; k++){
    this->list_particle[k].A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
  int* i1 = new int;
  int* j1 = new int;
  for (int i = 0; i<this->Nc;i++){
    for (int j = 0; j<this->Nc;j++){
      for (int k = -1;k<2;k++){
        for (int l = -1;l<2;l++){
          *i1 = i+k;
          *j1 = j+l;
          this->Grid.get_cell(i1,j1);
          cell *c1;
          c1 = &this->Grid.tab[*i1][*j1];
          for (int id0 : this->Grid.tab[i][j].list_index_particle){
            for (int id1 : this->Grid.tab[*i1][*j1].list_index_particle){
              if (id1 < id0){
                double dx = this->list_particle[id1].X.x+ c1->move_x -this->list_particle[id0].X.x;
                double dy = this->list_particle[id1].X.y+ c1->move_y-this->list_particle[id0].X.y;
                double r2 = dx*dx + dy*dy;
                double ir2 = 1.0/double(r2);
                    double ir6 = ir2*ir2*ir2;
                    double v = 24.0*ir6*(ir6-0.5);
                if (r2 < this->rc2 && r2 > this->diameter){
                     if (abs(v)>10){
                        v=1;
                      }
                      double f = 2.0*v*ir2;
                      double fx = f*dx;
                      double fy = f*dy;
                      this->list_particle[id1].A.x = this->list_particle[id1].A.x + fx*this->list_particle[id1].type;
                      this->list_particle[id1].A.y = this->list_particle[id1].A.y + fy*this->list_particle[id1].type;
                      this->list_particle[id0].A.x = this->list_particle[id0].A.x - fx*this->list_particle[id0].type;
                      this->list_particle[id0].A.y = this->list_particle[id0].A.y - fy*this->list_particle[id0].type;
                      this->E_pot += 4.0*ir6*(ir6-1.0);
                      this->viriel +=v;
                }
                if (r2 <this->diameter/10){
                  if (abs(v)>100){
                        v=1;
                  }
                    if (list_particle[id1].type == list_particle[id0].type){
                      this->list_particle[id1].V.x = -this->list_particle[id1].V.x;//this->list_particle[id0].type;
                      this->list_particle[id1].V.y = -this->list_particle[id1].V.y;//this->list_particle[id0].type;
                      this->list_particle[id0].V.x = -this->list_particle[id0].V.x; //* this->list_particle[id1].type;
                      this->list_particle[id0].V.y = -this->list_particle[id0].V.y; //* this->list_particle[id1].type;
                    }
                    if (list_particle[id1].type < list_particle[id0].type){
                      this->list_particle[id1].V.x = -this->list_particle[id1].V.x;//this->list_particle[id0].type;
                      this->list_particle[id1].V.y = -this->list_particle[id1].V.y;//this->list_particle[id0].type;
                      this->list_particle[id0].V.x = -this->list_particle[id0].V.x; //* this->list_particle[id1].type;
                      this->list_particle[id0].V.y = -this->list_particle[id0].V.y; //* this->list_particle[id1].type;
                    }
                    if (list_particle[id1].type > list_particle[id0].type){
                      this->list_particle[id1].V.x = -this->list_particle[id1].V.x;//this->list_particle[id0].type;
                      this->list_particle[id1].V.y = -this->list_particle[id1].V.y;//this->list_particle[id0].type;
                      this->list_particle[id0].V.x = -this->list_particle[id0].V.x; //* this->list_particle[id1].type;
                      this->list_particle[id0].V.y = -this->list_particle[id0].V.y; //* this->list_particle[id1].type;
                    }
                }
              } 
            }
          }
        }
      }
    }
  }
delete i1;
delete j1;  
}


void system1::compute_force_with_neighbour(){
  //clock_t begin = clock();
  for (int k = 0; k<this->N; k++){
    this->list_particle[k].A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
  //cout <<"taille : "<< list_neighbour.size() << endl;
  for (vect V : this->list_neighbour){
    particle *p,*p1;
    p = &this->list_particle[V.x];
    p1 = &this->list_particle[V.y];
    vect dX = p1->X - p->X;
    if (dX.x >= this->half_L){
      dX.x = dX.x - this->L;
    }
    if (dX.x < -this->half_L){
      dX.x = dX.x + this->L;
    }
    if (dX.y >= this->half_L){
      dX.y = dX.y - this->L;
    }
    if (dX.y < -this->half_L){
      dX.y = dX.y + this->L;
    }
    double r2 = dX.x*dX.x + dX.y*dX.y;
    if (r2 < this->rc2){
      double ir2 = 1.0/double(r2);
      double ir6 = ir2*ir2*ir2;
      double v = 24.0*ir6*(ir6-0.5);
      if (abs(v) > 5){v=1;}
      double f = 2.0*v*ir2;
      double fx = f*dX.x;
      double fy = f*dX.y;
      p->A = p->A - vect(fx,fy);
      p1->A = p1->A + vect(fx,fy);
      this->E_pot += 4.0*ir6*(ir6-1.0);
      this->viriel +=v;
    }
  }
//clock_t end = clock();
//unsigned long millis = (end -  begin) * 1000 / CLOCKS_PER_SEC;
//cout << "f: " << millis;
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
    particle *p;
    p = &list_particle[k];
    if ((p->A.x)*(p->A.x) + (p->A.y)*(p->A.y) > 10){p->A.x = 1;
    p->A.y = 1;}
    p->V = p->V + hd2*p->A;
    vect X1(p->X.x + h*p->V.x,p->X.y+h*p->V.y);
    this->move_particle(p,k,&X1);
  }
  
  this->compute_force();

  for (int k =0; k<this->N;k++){
    particle *p;
    p = &list_particle[k];
    p->V = p->V + hd2*p->A;
  }
}


void system1::verlet_neighbour(double h, double hd2){
  int nb = 0;
  for (int k =0; k<this->N;k++){
    particle *p;
    p = &list_particle[k];
    p->V = p->V + hd2* p->A;
    vect X1(p->X.x + h* p->V.x , p->X.y+h* p->V.y);
    this->move_particle(p,k,&X1);
  }
  this->compute_force_with_neighbour();
  double v2max = 0.0;
  for (int i =0; i<this->N;i++){
    particle *p;
    p = &list_particle[i];
    p->V = p->V + hd2*p->A;
    double v2 = p->V.x * p->V.x + p->V.y*p->V.y;
    if (v2 > v2max){
      v2max = v2;
    }
    this->move_max += sqrt(v2max)*h;
    
    if (this->move_max*2.0 > this->deltaR){
      this->construct_neighbour_list();
      nb++;
      this->move_max = 0.0;
    }
  }
cout <<"nb: "<< nb << endl;
}


void system1::integration(double h,int n){
  clock_t begin = clock();
  int taille = this->list_particle.size();
  fstream fichx, fichy, ficht;
  fichx.open("x.txt", ios::out);
  fichy.open("y.txt", ios::out);
  ficht.open("t.txt", ios::out);
  fichx << taille << endl;
  fichy << taille << endl;
  double hd2 = h*0.5;
  for (int i = 0; i<n;i++){
    system("clear");
    cout << i+1 << " / "<< n << endl;
    for (int k = 0; k<taille; k++){
      fichx << this->list_particle[k].X.x << endl;
      fichy << this->list_particle[k].X.y << endl;
      ficht << this->list_particle[k].type << endl;
    }
    this->verlet(h,hd2);
  }
clock_t end = clock();
unsigned long millis = (end -  begin) * 1000 / CLOCKS_PER_SEC;
cout << "Time without neighbours: " << millis<<" ms"<<endl;
}

void system1::integration_neighbour(double h,int n){
  clock_t begin = clock();
  double hd2 = h*0.5;
  int taille = this->list_particle.size();
  fstream fichx, fichy,ficht;
  fichx.open("x.txt", ios::out);
  fichy.open("y.txt", ios::out);
  ficht.open("t.txt", ios::out);
  fichx << taille << endl;
  fichy << taille << endl;
  for (int i = 0; i<n;i++){
    cout << i << endl;
    for (int k = 0; k<taille; k++){
      fichx << this->list_particle[k].X.x << endl;
      fichy << this->list_particle[k].X.y << endl;
      ficht << this->list_particle[k].type << endl;
    }
    this->verlet_neighbour(h,hd2);
  }
clock_t end = clock();
unsigned long millis = (end -  begin) * 1000 / CLOCKS_PER_SEC;
cout << "Time with neighbours: " << millis<<" ms"<<endl;
}
