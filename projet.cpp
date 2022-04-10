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
X(vect(0,0)),V(vect(0,0)),A(vect(0,0)),type(0)
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


grid::grid(int N_, int L_):
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
// transforms the coordinates respecting the periodic boundary conditions and adds the offset to the "move_x/y" variables
  int x_ = 0;
  int y_ = 0;
  if (*i < 0){
    *i += this->Nc;
    x_ = -this->L;
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

void grid::set_cell(int i,int j,int index){ 
// transforms the coordinates respecting the periodic boundary conditions and adds a new particle in the cell 
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


system1::system1(const system1& s) : grid(s),Nx(s.Nx),N(Nx*Nx),density(s.density),rc(s.rc),rc2(rc*rc),deltaR(s.deltaR),rv(rc+deltaR),rv2((rc+deltaR)*(rc+deltaR)),L(sqrt(M_PI/(double)density)*0.5*Nx),area(L*L),half_L(0.5*L),Nc(int(L/rv)),lc(L/(double)Nc),list_particle(s.list_particle),cutoff(s.cutoff)
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
// initializes a new particle of index "index", position X, velocity V and type "type_" 
  this->list_particle[index].type = type_;
  this->list_particle[index].X = X;
  this->list_particle[index].V = V;
  double i = int(X.x/double(this->lc));
  double j = int(X.y/double(this->lc));
  this -> Grid.set_cell(int(i),int(j),index);
}


void system1::init_system(double velocity,double size){
/* Initializes the system: we place the spheres regularly spaced with a random direction velocity, whose norm is given. 
A correction on the velocities is made so that the total momentum is zero. At the center of the square, we place the 
particles of type "1" (the magnetic ones) in a disk of radius sqrt(size*L) */
  double dx = this->L/double(this->Nx);
  if (dx <= diameter){throw "Density: too high!";}
  else{
    double dy = dx;
    double x = dx*0.5;
    double y = x;
    double px = 0.0;
    double py = 0.0;
    srand( (unsigned)time( NULL ) );
    for (int k = 0; k<this->N ; k++){
      double t = 0; 
      double a = (rand()/double(RAND_MAX))*M_PI*2.0;  
      double vx = velocity*cos(a);
      double vy = velocity*sin(a);
      px += vx;
      py += vy;
      if ((x - 0.5*L)*(x - 0.5*L) + (y - 0.5*L)*(y - 0.5*L) < size*L){t = 1;}
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
  }
compute_force();
construct_neighbour_list();
}


void system1::move_particle(particle* p, int index, vect* X1){
  /* Changes the position of the particle in X1. If this position corresponds 
  to a distance that exceeds 2*L (the particle is too fast) then it will leave
  the square, which we will avoid */
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
// Builds the list of nearest neighbor pairs of particles
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
} 


vect system1::compute_force_mag(particle p){
/* Calculates the sum of magnetic interactions with all particles of type "1" (not only the neighbors). 
Returns the force vector */ 
  double fx = 0.0;
  double fy = 0.0;
  for (int i = 0; i<this->Nc;i++){
    for (int j = 0; j<this->Nc;j++){
      this->Grid.get_cell(&i,&j);
      cell *c;
      c = &this->Grid.tab[i][j];
      for (int index : c->list_index_particle){
        if (this->list_particle[index].type == 1){
          double dx = this->list_particle[index].X.x+ c->move_x - p.X.x;
          double dy = this->list_particle[index].X.y+ c->move_y- p.X.y;
          double r2 = dx*dx + dy*dy;
          if (r2 > this->diameter/100.){
            double ir4 = (1/r2) * (1/r2);
            fx += +0.01 * ir4 * dx/sqrt(dx*dx + dy*dy);
            fy += +0.01 * ir4 * dy/sqrt(dx*dx + dy*dy);
          }
        }
      }
    }
  }
return vect(fx,fy);
}


void system1::compute_force(){
/* Calculates the force associated with the Lennard Jones potential 
and modifies the acceleration for particles that are less than "rc" away. 
Also takes into account the magnetic force for particles of type "1".
"mag" is a boolean deciding if we take into account the magnetic force or not */
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
            if (this->mag){
              if (this->list_particle[id0].type == 1){
                vect F = compute_force_mag(this->list_particle[id0]);
                if (sqrt(F.x * F.x + F.y * F.y) <= this->cutoff){
                this->list_particle[id0].A.x = this->list_particle[id0].A.x + F.x;
                this->list_particle[id0].A.y = this->list_particle[id0].A.y + F.y;
                }
              }
            }
            for (int id1 : this->Grid.tab[*i1][*j1].list_index_particle){
              if (id1 < id0){
                vect dX = this->list_particle[id1].X - this->list_particle[id0].X;
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
                double dx = dX.x;
                double dy = dX.y;
                //double dx = this->list_particle[id1].X.x+ c1->move_x -this->list_particle[id0].X.x;
                //double dy = this->list_particle[id1].X.y+ c1->move_y-this->list_particle[id0].X.y;
                double r2 = dx*dx + dy*dy;
                double ir2 = 1.0/double(r2);
                double ir6 = ir2*ir2*ir2;
                double v=0;
                double f=0;
                if (r2 < this->rc2){
                      if (this->list_particle[id0].type == 1 && this->list_particle[id1].type == 1){
                        v = (24.0/3.)*ir6*((ir6/0.0134)- 0.5/(0.3659));
                        f = 2.0*v*sqrt(ir2);
                        this->E_pot += 4.0*ir6*((ir6/0.0134)-1.0/0.3659);
                      }
                      if (this->list_particle[id0].type == 0 && this->list_particle[id1].type == 0){
                        v = (24.0)*ir6*(ir6-0.5);
                        f = 2.0*v*sqrt(ir2);
                        this->E_pot += 4.0*ir6*(ir6-1.0);
                      }
                      if (this->list_particle[id0].type != this->list_particle[id1].type){
                        v = (24.0/1.73)*ir6*(ir6/0.0233 -0.5/0.1527);
                        f = 2.0*v*sqrt(ir2);
                        this->E_pot += 4.0*ir6*((ir6/0.0233)-1.0/0.1527);
                      }
                      if (f> this->cutoff){f = 1;}
                      double fx = f*dx;
                      double fy = f*dy;
                      this->list_particle[id1].A.x = this->list_particle[id1].A.x + fx;
                      this->list_particle[id1].A.y = this->list_particle[id1].A.y + fy;
                      this->list_particle[id0].A.x = this->list_particle[id0].A.x - fx;
                      this->list_particle[id0].A.y = this->list_particle[id0].A.y - fy;
                      this->viriel +=v;
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
/* Calculates the force associated with the Lennard Jones potential 
and modifies the acceleration for the particles that are nearest neighbors. 
It is a test function to decrease the computation time. 
Thus, the magnetic force has not been implemented here. */
  for (int k = 0; k<this->N; k++){
    this->list_particle[k].A = vect(0,0);
    this->E_pot = 0.0;
    this->viriel = 0.0;
  }
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
    double ir2 = 1.0/double(r2);
    double ir6 = ir2*ir2*ir2;
    double v = 0;
    double f= 0;
    if (r2 < this->rc2){
      if (p->type == 1 && p1->type == 1){
        v = (24.0/3.)*ir6*((ir6/0.0134)- 0.5/(0.3659));
        f = 2.0*v*sqrt(ir2);
        this->E_pot += 4.0*ir6*((ir6/0.0134)-1.0/0.3659);
      }
      if (p->type == 0 && p1->type == 0){
        v = (24.0)*ir6*(ir6-0.5);
        f = 2.0*v*sqrt(ir2);
        this->E_pot += 4.0*ir6*(ir6-1.0);
      }
      if (p->type != p1->type){
        v = (24.0/1.73)*ir6*(ir6/0.0233 -0.5/0.1527);
        f = 2.0*v*sqrt(ir2);
        this->E_pot += 4.0*ir6*((ir6/0.0233)-1.0/0.1527);
      }
      if (f> this->cutoff){f = 1;}
      double fx = f*dX.x;
      double fy = f*dX.y;
      p1->A.x = p1->A.x + fx;
      p1->A.y = p1->A.y + fy;
      p->A.x = p->A.x - fx;
      p->A.y = p->A.y - fy;
      this->viriel +=v;
    }
  }
}


void system1::compute_E_kin(){
/*Calculate the kinetic and total energy, the instantaneous pressure 
as well as the sums allowing the calculation of the average values then */
  this->E_kin = 0.0;
  for (int i = 0; i<this->N;i++){
    particle p = this->list_particle[i];
    this->E_kin += 0.5*(p.V.x*p.V.x + p.V.y*p.V.y);
  }
  this->pressure = (this->viriel + this->E_kin)/double(this->area);
  this->E_kin /= double(this->N);
  this->energy = this->E_kin+(this->E_pot/double(this->N));
  this->E_pot /= double(this->N);
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
// Calculates the mean and standard deviation of temperature
  double Tm = this->sum_temp/double(this->counter);
  v->x = Tm;
  v->y = sqrt(this->sum_temp2/double(this->counter) - Tm*Tm);
}


void system1::mean_pressure(vect* v){
//Calculates the mean and standard deviation of pressure
  double Pm= this->sum_pressure/double(this->counter);
  v->x = Pm;
  v->y = sqrt(this->sum_pressure2/double(this->counter) - Pm*Pm);
}


void system1::adjust_v(double T){
// Performs a speed adjustment to obtain a given temperature
  double Tm = this->sum_temp/double(this->counter);
  double f = sqrt(T/double(Tm));
  for (int k = 0; k<this->N;k++){
    this->list_particle[k].V = this->list_particle[k].V*f;
  }
}


void system1::verlet(double h, double hd2){
// Performs the elementary step of the Verlet method, using the cells to calculate the forces
  for (int k =0; k<this->N;k++){
    particle *p;
    p = &list_particle[k];
    p->V = p->V + hd2*p->A;
    vect X1(p->X.x + h*p->V.x, p->X.y + h*p->V.y);
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
/* Performs the elementary step of the Verlet method, using the list of neighbors. 
In the second evaluation of the velocities, the maximum velocity is also calculated, 
in order to increment the maximum displacement of the spheres. If this exceeds delta_r, 
the neighbor list is reconstructed. */
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
    double v2 = p->V.x * p->V.x + p->V.y * p->V.y;
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
  this->adjust_v(0.5);
}


void system1::integration(double h,int n){
  double hd2 = h*0.5;
  for (int i = 0; i<n;i++){
    this->verlet(h,hd2);
    }
}

void system1::integration_neighbour(double h,int n){
  double hd2 = h*0.5;
  for (int i = 0; i<n;i++){
    this->verlet_neighbour(h,hd2);
  }
}
