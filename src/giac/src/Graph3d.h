// -*- mode:C++ ; compile-command: "g++ -I.. -g -c Graph3d.cc" -*-
#ifndef _GRAPH3D_H
#define _GRAPH3D_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "giacPCH.h"
#include "giac.h"
#include <fstream>
#include <string>
#include <stdio.h>
#ifdef HAVE_LIBFLTK_GL
#include <FL/Fl_Menu.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Tile.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Text_Buffer.H>
#include <FL/Fl_Text_Editor.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#endif

#ifdef HAVE_LC_MESSAGES
#include <locale.h>
#endif
#include "giacintl.h"
#include "Xcas1.h"
#include "Input.h"
#include "Editeur.h"

#ifndef NO_NAMESPACE_XCAS
namespace xcas {
#endif // ndef NO_NAMESPACE_XCAS
#ifdef HAVE_LIBFLTK

#ifdef HAVE_LIBFLTK_GL

  // translate giac GL constant to open GL constant
  unsigned gl_translate(unsigned i);
  // utilities for matrix 4x4 represented as a double[16] 
  // written in columns
  void mult4(double * colmat,double * vect,double * res);
  void mult4(double * colmat,float * vect,double * res);
  void mult4(double * c,double k,double * res);
  double det4(double * c);
  void inv4(double * c,double * res);
  // return in i and j the distance to the BOTTOM LEFT of the window
  // use window()->h()-j for the FLTK coordinates in this window
  void dim32dim2(double * view,double * proj,double * model,double x0,double y0,double z0,double & i,double & j,double & depth);
  // quaternion for the rotation of axis (x,y,z) angle theta
  quaternion_double rotation_2_quaternion_double(double x, double y, double z,double theta);

  // native 3d rendering types/defs
  typedef double float3d;
  // typedef float float3d;
  struct double3 {
    float3d x,y,z;
    double3(double x_,double y_,double z_):x(x_),y(y_),z(z_){};
    double3():x(0),y(0),z(0){};
  };

  struct int4 {
    int u,d,du,dd;
    int4(int u_,int d_,int du_,int dd_):u(u_),d(d_),du(du_),dd(dd_) {}
  };
  
  struct int2 {
    int i,j;
    int2(int i_,int j_):i(i_),j(j_) {}
    int2(): i(0),j(0) {}
  };
  inline bool operator < (const int2 & a,const int2 & b){ if (a.i!=b.i) return a.i<b.i; return a.j<b.j;}
  inline bool operator == (const int2 & a,const int2 & b){ return a.i==b.i && a.j==b.j;}

  struct int2_double2 {
    int i,j;
    double arg,norm;
  };
  inline bool operator < (const int2_double2 & a,const int2_double2 & b){ if (a.arg!=b.arg) return a.arg<b.arg; return a.norm<b.norm;}

#define giac3d_default_upcolor 65535
#define giac3d_default_downcolor 12345
#define giac3d_default_downupcolor 18432 // 12297
#define giac3d_default_downdowncolor 22539
#define COLOR_GREEN _GREEN
#define COLOR_RED _RED
#define COLOR_BLUE _BLUE
#define COLOR_CYAN _CYAN
#define COLOR_BLACK _BLACK
#define COLOR_YELLOW _YELLOW
#define COLOR_MAGENTA _MAGENTA
#define COLOR_WHITE _WHITE

  // end native 3d rendering types/defs

  class Graph3d : public Graph2d3d {
  public:
    Graph3d(int x,int y,int width, int height, const char* title,History_Pack * hp_);
    virtual ~Graph3d();
    double theta_z,theta_x,theta_y; // rotations
    double delta_theta; // rotation increment
    int draw_mode; // for sphere drawing
    FILE * printing; // non 0 if printing with gl2ps
    void * glcontext;
    // save values of the projection and modelview matrices
    double proj[16],model[16],proj_inv[16],model_inv[16]; 
    double view[4];
    int dragi,dragj;
    double depth;
    bool push_in_area;
    bool below_depth_hidden;
    unsigned char * screenbuf;
    virtual void resize(int X,int Y,int W,int H);
    virtual void draw();
    virtual void orthonormalize();
    void display(); 
    void displaypolyg(const std::vector<int2> & polyg,const int2 & IJmin,int color,int & Px,int & Py) const ;
    // internally callled by draw, maybe multiple times when printing
    void print(); // assumes that printing is assigned to a FILE *
    virtual int in_handle(int event);
    void indraw(const giac::vecteur & v);
    void indraw(const giac::gen & g);
    void legende_draw(const giac::gen & g,const std::string & s,int mode);
    void draw_string(const std::string & s);
    virtual const char * latex(const char * filename);
    // i,j,z -> x,y
    virtual void find_xyz(double i,double j,double depth,double & x,double & y,double & z)  ;
    // x,y,z -> FLTK coordinates i,j
    void find_ij(double x,double y,double z,double & i,double & j,double & depth) ;
    void current_normal(double & a,double &b,double &c) ;
    void normal2plan(double & a,double &b,double &c);
    virtual void geometry_round(double x,double y,double z,double eps,giac::gen & tmp,const giac::context *) ;
    int opengl2png(const std::string & filename);
    /* 
       build-in 3d rendering
    */
    void switch_renderer();
    void reset_view();
    virtual void update_rotation();
    void glinter1(double z,double dz,
		  double *zmin,double *zmax,double ZMIN,double ZMAX,
		  int ih,int lcdz,
		  int upcolor,int downcolor,int diffusionz,int diffusionz_limit,bool interval
			 ) const;
    void glinter(double a,double b,double c,
		 double xscale,double xc,double yscale,double yc,
		 double *zmin,double *zmax,double ZMIN,double ZMAX,
		 int i,int horiz,int j,int w,int h,int lcdz,
		 int upcolor,int downcolor,int diffusionz,int diffusionz_limit,bool interval
		 ) const;
  void glinter(double x1,double x2,double x3,
	       double y1,double y2,double y3,
	       double z1,double z2,double z3,
	       double xscale,double xc,double yscale,double yc,
	       double *zmin,double *zmax,double ZMIN,double ZMAX,
	       int i,int horiz,int j,int w,int h,int lcdz,
	       int upcolor,int downcolor,int diffusionz,int diffusionz_limit,bool interval
	       ) const;

    bool indraw3d(int w,int h,int lcdz,const giac::context*,int upcolor,int downcolor,int downupcolor,int downdowncolor) ;
    void draw_decorations(const giac::gen & title_tmp);
    void xyz2ij(const double3 & d,int &i,int &j) const; // d not transformed
    void xyz2ij(const double3 & d,double &i,double &j) const; // d not transformed
    void xyz2ij(const double3 & d,double &i,double &j,double3 & d3) const; // d not transformed, d3 is
    void XYZ2ij(const double3 & d,int &i,int &j) const; // d is transformed
    void draw3d(const giac::context *); // 3d rendering engine if opengl is false (e.g. for mac os)
    void addpolyg(vector<int2> & polyg,double x,double y,double z,int2 & IJmin) const ;
    void adddepth(vector<int2> & polyg,const double3 &A,const double3 &B,int2 & IJmin) const;
    double transform3d[16],invtransform3d[16];
    int show_edges,lcdz,default_upcolor,default_downcolor,default_downupcolor,default_downdowncolor,LCD_WIDTH_PX,LCD_HEIGHT_PX;
    short int precision,diffusionz,diffusionz_limit;
    bool opengl,doprecise,hide2nd,interval,solid3d,displayhint;
    double Ai,Aj,Bi,Bj,Ci,Cj,Di,Dj,Ei,Ej,Fi,Fj,Gi,Gj,Hi,Hj; // visualization cube coordinates
    std::vector< std::vector< std::vector<float3d> > > surfacev;
    std::vector<double3> plan_pointv; // point in plan 
    std::vector<double3> plan_abcv; // plan equation z=a*x+b*y+c
    std::vector<bool> plan_filled;
    std::vector<double3> sphere_centerv;
    std::vector<double> sphere_radiusv;
    giac::vecteur sphere_quadraticv; // matrix of the transformed quad. form
    std::vector<bool> sphere_isclipped;
    std::vector< std::vector<double3> > polyedrev;
    std::vector<double3> polyedre_abcv;
    std::vector<double> polyedre_xyminmax;
    std::vector<bool> polyedre_faceisclipped,polyedre_filled;
    std::vector<double3> linev; // 2 double3 per object
    std::vector<short> linetypev;
    std::vector<const char *> lines; // legende
    std::vector< std::vector<double3> > curvev;
    std::vector<double3> pointv; 
    std::vector<const char *> points; // legende
    std::vector<int4> hyp_color,plan_color,sphere_color,polyedre_color,line_color,curve_color,point_color;

  };

  class Geo3d : public Graph3d {
  public:
    virtual FL_EXPORT void draw();
    virtual int in_handle(int event);
    Geo3d(int x,int y,int width, int height, History_Pack * _hp);
  };
#else
  class Graph3d : public Graph2d3d {
  public:
    void print(){}; // assumes that printing is assigned to a FILE *
    Graph3d(int x,int y,int width, int height, const char* title,History_Pack * hp_):Graph2d3d(x,y,width,height,title,hp_),printing(0) {};
    double theta_z,theta_x,theta_y; // rotations
    double delta_theta; // rotation increment
    int draw_mode; // for sphere drawing
    FILE * printing; // non 0 if printing with gl2ps
    void * glcontext;
    // save values of the projection and modelview matrices
    double proj[16],model[16],proj_inv[16],model_inv[16]; 
    double view[4];
    int dragi,dragj;
    bool push_in_area;
    double depth;
    bool below_depth_hidden;
    unsigned char * screenbuf;
    // virtual void resize(int X,int Y,int W,int H);
    // virtual void draw();
    // virtual void orthonormalize();
    // void display(); 
    // void indraw(const giac::vecteur & v);
    // void indraw(const giac::gen & g);
    // void legende_draw(const giac::gen & g,const std::string & s,int mode);
    // void draw_string(const std::string & s);
    // virtual const char * latex(const char * filename);
    // i,j,z -> x,y
    // virtual void find_xyz(double i,double j,double depth,double & x,double & y,double & z)  ;
    // x,y,z -> FLTK coordinates i,j
    // void find_ij(double x,double y,double z,double & i,double & j,double & depth) ;
    void current_normal(double & a,double &b,double &c) {};
    void normal2plan(double & a,double &b,double &c){};
    // virtual void geometry_round(double x,double y,double z,double eps,giac::gen & tmp,const giac::context *) ;
    int opengl2png(const std::string & filename){ return 0;};
  };

  class Geo3d : public Graph3d {
  public:
    Geo3d(int x,int y,int width, int height, History_Pack * _hp):Graph3d(x,y,width,height,0,_hp){};
    virtual FL_EXPORT void draw(){};
  };

  inline void mult4(double * colmat,double * vect,double * res){};
  inline void mult4(double * colmat,float * vect,double * res){};
  inline void mult4(double * c,double k,double * res){};
  inline quaternion_double rotation_2_quaternion_double(double x, double y, double z,double theta){};

#endif // HAVE_LIBFLTK_GL
#endif // HAVE_LIBFLTK

#ifndef NO_NAMESPACE_XCAS
} // namespace xcas
#endif // ndef NO_NAMESPACE_XCAS

#endif // _GRAPH3D_H
