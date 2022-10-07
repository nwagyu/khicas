/* -*- mode:C++ ; compile-command: "g++ -I.. -I../include -DHAVE_CONFIG_H -DIN_GIAC -DGIAC_GENERIC_CONSTANTS -fno-strict-aliasing -g -c mp.cc -Wall" -*- */
#include "giacPCH.h"
#include "mp.h"
/*
 *  Copyright (C) 2020 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

using namespace std;
using namespace giac;

context * caseval_context(){
  return (context *) caseval("caseval contextptr");
}
void c_draw_rectangle(int x,int y,int w,int h,int c){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_line(x,y,x+w,y,c,contextptr);
  draw_line(x+w,y,x+w,y+h,c,contextptr);
  draw_line(x,y+h,x+w,y+h,c,contextptr);
  draw_line(x,y,x,y+h,c,contextptr);
}
void c_draw_line(int x0,int y0,int x1,int y1,int c){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_line(x0,y0,x1,y1,c,contextptr);
}
void c_draw_circle(int xc,int yc,int r,int color,bool q1,bool q2,bool q3,bool q4){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_circle(xc,yc,r,color,q1,q2,q3,q4,contextptr);
}
void c_draw_filled_circle(int xc,int yc,int r,int color,bool left,bool right){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_filled_circle(xc,yc,r,color,left,right,contextptr);
}
void c_convert(int *x,int*y,vector< vector<int> > & v){
  for (int i=0;i<v.size();++i,++x,++y){
    v[i].push_back(*x);
    v[i].push_back(*y);
  }
}
void c_draw_polygon(int * x,int *y ,int n,int color){
  giac::freeze=true;
  vector< vector<int> > v(n);
  c_convert(x,y,v);
  context * contextptr=caseval_context();
  draw_polygon(v,color,contextptr);
}
void c_draw_filled_polygon(int * x,int *y, int n,int xmin,int xmax,int ymin,int ymax,int color){
  giac::freeze=true;
  vector< vector<int> > v(n);
  c_convert(x,y,v);
  context * contextptr=caseval_context();
  draw_filled_polygon(v,xmin,xmax,ymin,ymax,color,contextptr);
}
void c_draw_arc(int xc,int yc,int rx,int ry,int color,double theta1, double theta2){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_arc(xc,yc,rx,ry,color,theta1,theta2,contextptr);
}
void c_draw_filled_arc(int x,int y,int rx,int ry,int theta1_deg,int theta2_deg,int color,int xmin,int xmax,int ymin,int ymax,bool segment){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_filled_arc(x,y,rx,ry,theta1_deg,theta2_deg,color,xmin,xmax,ymin,ymax,segment,contextptr);
}
void c_set_pixel(int x,int y,int c){
  giac::freeze=true;
  context * contextptr=caseval_context();
  set_pixel(x,y,c,contextptr);
}
void c_fill_rect(int x,int y,int w,int h,int c){
  giac::freeze=true;
  context * contextptr=caseval_context();
  draw_rectangle(x,y,w,h,c,contextptr);
}

int c_draw_string(int x,int y,int c,int bg,const char * s,bool fake){
  giac::freeze=true;
  context * contextptr=caseval_context();
  _draw_string(makevecteur(string2gen(s,false),x,y,c),contextptr);
  return 0;//os_draw_string(x,y,c,bg,s,fake);
}
int c_draw_string_small(int x,int y,int c,int bg,const char * s,bool fake){
  giac::freeze=true;
  context * contextptr=caseval_context();
  _draw_string(makevecteur(string2gen(s,false),x,y,c),contextptr);
  return 0;//os_draw_string_small(x,y,c,bg,s,fake);
}
int c_draw_string_medium(int x,int y,int c,int bg,const char * s,bool fake){
  giac::freeze=true;
  context * contextptr=caseval_context();
  _draw_string(makevecteur(string2gen(s,false),x,y,c),contextptr);
  return 0;//os_draw_string_medium(x,y,c,bg,s,fake);
}

int os_get_pixel(int x,int y){
  context * contextptr=caseval_context();
  return _get_pixel(makevecteur(x,y),contextptr).val;
}

ulonglong double2gen(double d){
  giac::gen g(d);
  return *(ulonglong *) &g;
}

ulonglong int2gen(int d){
  giac::gen g(d);
  return *(ulonglong *) &g;
}

//bool freezeturtle=false;
void turtle_freeze(){
  freezeturtle=true;
}

void doubleptr2matrice(double * x,int n,int m,giac::matrice & M){
  M.resize(n);
  for (int i=0;i<n;++i){
    M[i]=giac::vecteur(m);
    giac::vecteur & w=*M[i]._VECTptr;
    for (int j=0;j<m;++j){
      w[j]=*x;
      ++x;
    }
  }
}

// x must have enough space!
bool matrice2doubleptr(const giac::matrice &M,double *x){
  int n=M.size();
  if (n==0 || M.front().type!=giac::_VECT)
    return false;
  int m=M.front()._VECTptr->size();
  for (int i=0;i<n;++i){
    if (M[i].type!=giac::_VECT || M[i]._VECTptr->size()!=m)
      return false;
    giac::vecteur & w=*M[i]._VECTptr;
    for (int j=0;j<m;++j){
      giac::gen g =giac::evalf_double(w[j],1,giac::context0);
      if (g.type!=giac::_DOUBLE_)
	return false;
      *x=g._DOUBLE_val;
      ++x;
    }
  }
  return true;
}

bool r_inv(double * x,int n){
  giac::matrice M(n);
  doubleptr2matrice(x,n,n,M);
  M=giac::minv(M,giac::context0);
  return matrice2doubleptr(M,x);
}


bool r_rref(double * x,int n,int m){
  giac::matrice M(n);
  doubleptr2matrice(x,n,m,M);
  giac::gen g=giac::_rref(M,giac::context0);
  if (g.type!=giac::_VECT)
    return false;
  return matrice2doubleptr(*g._VECTptr,x);
}

double r_det(double *x,int n){
  giac::matrice M(n);
  doubleptr2matrice(x,n,n,M);
  giac::gen g=giac::mdet(M,giac::context0);
  g=giac::evalf_double(g,1,giac::context0);
  double d=1.0,e=1.0;
  if (g.type!=_DOUBLE_)
    return 0.0/(d-e);
  return g._DOUBLE_val;
}

void c_complexptr2matrice(c_complex * x,int n,int m,giac::matrice & M){
  M.resize(n);
  for (int i=0;i<n;++i){
    if (m==0){
      M[i]=gen(x->r,x->i);
      ++x;
      continue;
    }
    M[i]=giac::vecteur(m);
    giac::vecteur & w=*M[i]._VECTptr;
    for (int j=0;j<m;++j){
      w[j]=gen(x->r,x->i);
      ++x;
    }
  }
}

c_complex gen2c_complex(giac::gen & g){
  double d=1.0,e=1.0;
  c_complex c={0,0};
  if (g.type!=giac::_DOUBLE_ && g.type!=giac::_CPLX)
    c.r=c.i=0.0/(d-e);
  else {
    if (g.type==giac::_DOUBLE_)
      c.r=g._DOUBLE_val;
    else {
      if (g.subtype!=3)
	c.r=c.i=0.0/(d-e);
      c.r=g._CPLXptr->_DOUBLE_val;
      c.i=(g._CPLXptr+1)->_DOUBLE_val;
    }
  }
  return c;
}

// x must have enough space!
bool matrice2c_complexptr(const giac::matrice &M,c_complex *x){
  int n=M.size();
  if (n==0)
    return false;
  if (M.front().type!=giac::_VECT){
    for (int i=0;i<n;++i){
      giac::gen g =giac::evalf_double(M[i],1,giac::context0);
      if (g.type!=giac::_DOUBLE_ && g.type!=giac::_CPLX)
	return false;
      *x=gen2c_complex(g);
      ++x;
    }
    return true;
  }
  int m=M.front()._VECTptr->size();
  for (int i=0;i<n;++i){
    if (M[i].type!=giac::_VECT || M[i]._VECTptr->size()!=m)
      return false;
    giac::vecteur & w=*M[i]._VECTptr;
    for (int j=0;j<m;++j){
      giac::gen g =giac::evalf_double(w[j],1,giac::context0);
      if (g.type!=giac::_DOUBLE_ && g.type!=giac::_CPLX)
	return false;
      *x=gen2c_complex(g);
      ++x;
    }
  }
  return true;
}

bool c_inv(c_complex * x,int n){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,n,M);
  M=giac::minv(M,giac::context0);
  return matrice2c_complexptr(M,x);
}

bool c_proot(c_complex * x,int n){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,0,M);
  M=giac::proot(M);
  return matrice2c_complexptr(M,x);
}

bool c_pcoeff(c_complex * x,int n){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,0,M);
  M=giac::pcoeff(M);
  return matrice2c_complexptr(M,x);
}

bool c_fft(c_complex * x,int n,bool inverse){
#if 1
  complex<double> * X=(complex<double> *) x;
  double theta=2*M_PI/n;
  if (!inverse)
    theta=-theta;
  fft2(X,n,theta);
  if (inverse){
    for (int i=0;i<n;++i)
      X[i]=X[i]/double(n);
  }
  return true;
#else
  giac::matrice M(n);
  c_complexptr2matrice(x,n,0,M);
  gen g=inverse?giac::_ifft(M,giac::context0):giac::_fft(M,giac::context0);
  if (g.type!=_VECT)
    return false;
  return matrice2c_complexptr(*g._VECTptr,x);
#endif
}

bool c_egv(c_complex * x,int n){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,n,M);
  gen g=giac::_egv(M,giac::context0);
  if (!ckmatrix(g))
    return false;
  return matrice2c_complexptr(*g._VECTptr,x);
}

bool c_eig(c_complex * x,c_complex * d,int n){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,n,M);
  gen g=giac::_jordan(M,giac::context0);
  if (g.type!=_VECT || g._VECTptr->size()!=2 || !ckmatrix(g[0]) || !ckmatrix(g[1]))
    return false;
  return matrice2c_complexptr(*g[0]._VECTptr,x) && matrice2c_complexptr(*g[1]._VECTptr,d);
}

bool c_rref(c_complex * x,int n,int m){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,m,M);
  giac::gen g=giac::_rref(M,giac::context0);
  if (g.type!=giac::_VECT)
    return false;
  return matrice2c_complexptr(*g._VECTptr,x);
}

c_complex c_det(c_complex *x,int n){
  giac::matrice M(n);
  c_complexptr2matrice(x,n,n,M);
  giac::gen g=giac::mdet(M,giac::context0);
  g=giac::evalf_double(g,1,giac::context0);
  return gen2c_complex(g);
}

void c_sprint_double(char * s,double d){
  giac::sprint_double(s,d);
}

void c_turtle_forward(double d){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  //const context * contextptr=caseval_context();
  giac::_avance(d,cascontextptr);
}

void c_turtle_left(double d){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  giac::_tourne_gauche(d,cascontextptr);
}

void c_turtle_up(int i){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  if (i)
    giac::_leve_crayon(0,cascontextptr);
  else
    giac::_baisse_crayon(0,cascontextptr);
}

void c_turtle_goto(double x,double y){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  giac::_position(makesequence(x,y),cascontextptr);
}

void c_turtle_cap(double x){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  giac::_cap(x,cascontextptr);
}

void c_turtle_crayon(int i){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  giac::_crayon(i,cascontextptr);
}

void c_turtle_rond(int x,int y,int z){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  giac::_rond(makesequence(x,y,z),cascontextptr);
}

void c_turtle_disque(int x,int y,int z,int centre){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  if (centre)
    giac::_disque_centre(makesequence(x,y,z),cascontextptr);
  else
    giac::_disque(makesequence(x,y,z),cascontextptr);
}

void c_turtle_fill(int i){
  gen arg(vecteur(0));
  if (i==0) 
    arg.subtype=_SEQ__VECT;
  context * cascontextptr=(context *)caseval("caseval contextptr");
  giac::_polygone_rempli(arg,cascontextptr);
}

void c_turtle_fillcolor(double r,double g,double b,int entier){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  if (entier)
    giac::_polygone_rempli(makesequence(int(r),int(g),int(b)),cascontextptr);
  else
    giac::_polygone_rempli(makesequence(r,g,b),cascontextptr);
}

void c_turtle_getposition(double * x,double * y){
  context * cascontextptr=(context *)caseval("caseval contextptr");
  gen arg(vecteur(0)); arg.subtype=_SEQ__VECT;
  giac::gen g=giac::_position(arg,cascontextptr);
  if (g.type==_VECT && g._VECTptr->size()==2){
    gen a=g._VECTptr->front(),b=g._VECTptr->back();
    a=evalf_double(a,1,cascontextptr);
    b=evalf_double(b,1,cascontextptr);
    *x=a._DOUBLE_val;
    *y=b._DOUBLE_val;
  }
}
