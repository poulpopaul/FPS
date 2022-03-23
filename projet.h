#ifndef projet_h
#define projet_h

#include<iostream>
#include<vector>
#include<cstring>
using namespace std;

class vect{
  // vect creation for position, velocity and acceleration.
 public:
  char * vect_name = new char[15];
  double x;
  double y;
  vect();
  vect(const vect &);
  vect& operator=(const vect&);
  vect(double x_, double y_);
  void display();
  virtual ~vect();
  friend vect operator+(vect v1, vect v2);
  friend vect operator-(vect v1, vect v2);
  friend vect operator*(vect v1, vect v2);
  friend vect operator*(vect v1, double l);
  friend vect operator*(double l, vect v1);
};



class particle: public vect {
  // A particule is descibred by its coordinates, velocity and acceleration.
 public:
  vect X;
  vect V;
  vect A;
  char * particle_name = new char[15];
  particle();
  particle(const particle &);
  particle(vect X_, vect V_, vect A_);
  void set_x(double x_);
  void set_y(double y_);
  void display();
  virtual ~particle();
  particle& operator=(const particle&);
};

class cell: public particle{
  // A cell is an elementary square containing particles.
  // Thus, we must be able to add or remove  particles from the cell.
  // A cell will be described by its (x,y) coordinates in a grid.
public:
vector<int> list_index_particle;
double x;
double y;
int n = 0;
char * cell_name = new char[15];
cell();
//~cell();
cell(double x_, double y_);
cell(const cell &);
virtual ~cell();
void add_particle(int index);
void remove_particle(int index);
void display();
cell& operator=(const cell&);
};

class grid: public cell{
  public:
  int L;
  int Nc;
  char * grid_name = new char[15];
  cell **tab = new cell*[Nc];
  grid();
  grid( const grid&g);
  grid(int L_, int N_);
  virtual ~grid();
  cell get_cell(int i, int j);
  void set_cell(int i,int j, int index); //allow to update list_particle easily
  void remove(int i,int j,int index);
  grid& operator=(const grid&);



};
class system1 : public grid{
public:
    char * system1_name = new char[15];
    const double PI  =3.141592653589793238463;
    int Nx;
    int N;
    double density;
    double rc;
    double rc2;
    double deltaR;
    double rv;
    double rv2;
    double radius = 0.5;
    double diameter = 1;
    double L;
    double area;
    double half_L;
    double Nc;
    double lc;
    grid Grid = grid(Nc,L);;
    vector<particle> list_particle; //= {}
    double energy = 0;
    double viriel = 0;
    double E_kin = 0;
    double E_pot = 0;
    double pressure = 0;
    double counter = 0;
    double sum_temp = 0;
    double sum_pressure = 0;
    double sum_temp2 = 0;
    double sum_pressure2 = 0;
    vector<vect> list_neighbour = {};
    double move_max = 0.0;
    system1();
    system1(int Nx_, double density_, double rc_, double deltaR_);
    system1(const system1&s);
    virtual ~system1();
    void init_particle(int i,vect X,vect V);
    void init_system(double velocity);
    void move_particle(particle* p, int index, vect* X1);
    void construct_neighbour_list();
    void compute_force();
    void compute_force_with_neighbour();
    void compute_E_kin();
    void init_mean();
    void mean_temp(vect* v);
    void mean_pressure(vect* v);
    void adjust_v(double T);
    void verlet(double h, double hd2);
    void verlet_neighbour(double h, double hd2);
    void integration(double h,int n);
    void integration_neighbour(double h,int n);
    system1& operator=(const system1&);
    };

#endif
