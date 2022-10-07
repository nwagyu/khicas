// -*- mode:C++ ; compile-command: "g++ -I.. -g -c Equation.cc" -*-
#ifndef _MP_H
#define _MP_H
#include "config.h"
#include "giacPCH.h"
extern "C" {
  void c_draw_rectangle(int x,int y,int w,int h,int c);
  void c_fill_rect(int x,int y,int w,int h,int c);
  void c_draw_line(int x0,int y0,int x1,int y1,int c);
  void c_draw_circle(int xc,int yc,int r,int color,bool q1,bool q2,bool q3,bool q4);
  void c_draw_filled_circle(int xc,int yc,int r,int color,bool left,bool right);
  void c_draw_polygon(int * x,int *y ,int n,int color);
  void c_draw_filled_polygon(int * x,int *y, int n,int xmin,int xmax,int ymin,int ymax,int color);
  void c_draw_arc(int xc,int yc,int rx,int ry,int color,double theta1, double theta2);
  void c_draw_filled_arc(int x,int y,int rx,int ry,int theta1_deg,int theta2_deg,int color,int xmin,int xmax,int ymin,int ymax,bool segment);
  void c_set_pixel(int x,int y,int c);
  int c_draw_string(int x,int y,int c,int bg,const char * s,bool fake);
  int c_draw_string_small(int x,int y,int c,int bg,const char * s,bool fake);
  int c_draw_string_medium(int x,int y,int c,int bg,const char * s,bool fake);
  int select_item(const char ** ptr,const char * title,bool askfor1=true);
  // C conversion to gen from atomic data type 
  unsigned long long c_double2gen(double); 
  unsigned long long c_int2gen(int);
  // linalg on double matrices
  void doubleptr2matrice(double * x,int n,giac::matrice & m);
  bool matrice2doubleptr(const giac::matrice &M,double *x); // x must have enough space
  bool r_inv(double *,int n);
  bool r_rref(double *,int n,int m);
  double r_det(double *,int);
  struct double_pair {
    double r,i;
  } ;
  typedef struct double_pair c_complex;
  bool matrice2c_complexptr(const giac::matrice &M,c_complex *x);
  void c_complexptr2matrice(c_complex * x,int n,int m,giac::matrice & M);  
  bool c_inv(c_complex *,int n);
  bool c_rref(c_complex *,int n,int m);
  c_complex c_det(c_complex *,int);
  bool c_egv(c_complex * x,int n); // eigenvectors
  bool c_eig(c_complex * x,c_complex * d,int n); // x eigenvect, d reduced mat
  bool c_proot(c_complex * x,int n); // poly root
  bool c_pcoeff(c_complex * x,int n); // root->coeffs
  bool c_fft(c_complex * x,int n,bool inverse); // FFT
  void c_sprint_double(char * s,double d);
  int os_get_pixel(int x,int y);
  void sync_screen();
  void turtle_freeze();
  void c_turtle_forward(double d);
  void c_turtle_left(double d);
  void c_turtle_up(int i);
  void c_turtle_goto(double x,double y);
  void c_turtle_cap(double x);
  void c_turtle_crayon(int i);
  void c_turtle_rond(int x,int y,int z);
  void c_turtle_disque(int x,int y,int z,int centered);
  void c_turtle_fill(int i);
  void c_turtle_fillcolor(double r,double g,double b,int entier);
  void c_turtle_getposition(double * x,double * y);
}
#endif // MP_H
