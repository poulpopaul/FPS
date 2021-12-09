#ifndef projet_h
#define projet_h

class vector{
  // Vector creation for position, velocity and acceleration.
 public:
  double x;
  double y;
  vector();
  vector(double x_, double y_);
  void display();
  //  ~vector();
  friend vector operator+(vector v1, vector v2);
  friend vector operator-(vector v1, vector v2);
  friend vector operator*(vector v1, vector v2);
  friend vector operator*(vector v1, double l);
  friend vector operator*(double l, vector v1);
  
};



class particle: public vector {
  // A particule is descibred by its coordinates, velocity and acceleration.
 public:
  vector X;
  vector V;
  vector A;
  particle();
  particle(vector X_, vector V_, vector A_);
  void display();
  //  ~particule();
};

//class cell: public particule
//{
  // A cell is an elementary square containing particles. Thus, we must be able to add or remove  particles from the cell.
  
  
//};

#endif
