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
  vect(double x_, double y_);
  void display();
  ~vect();
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
  void display();
  ~particle();
};

class cell: public particle{
  // A cell is an elementary square containing particles.
  // Thus, we must be able to add or remove  particles from the cell.
  // A cell will be described by its (x,y) coordinates in a grid.
public:
vector<unsigned int> list_particle = {};
unsigned int x;
unsigned int y;
unsigned int n = 0;
char * cell_name = new char[15];
cell();
cell(unsigned int x_, unsigned int y_);
cell(const cell &);
~cell();
void add_particle(unsigned int index);
void remove_particle(unsigned int index);
void display();
};

class grid: public cell{
  public:
  int L;
  int N;
  cell **tab = new cell*[N];
  grid();
  grid( const grid&);
  grid(int L_, int N_);
  cell get_cell(int i, int j);


};


#endif
