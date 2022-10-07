// -*- mode:C++ ; compile-command: "g++-3.4 -I. -I.. -I../include -I../../giac/include -g -c Graph3d.cc -DIN_GIAC -DHAVE_CONFIG_H" -*-
#include "Graph3d.h"
/*
 *  Copyright (C) 2006,2014 B. Parisse, Institut Fourier, 38402 St Martin d'Heres
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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define __CARBONSOUND__
#ifdef HAVE_LIBFLTK_GL
#include <FL/fl_ask.H>
#include <FL/fl_ask.H>
#include <FL/Fl.H>
#include <FL/Fl_File_Chooser.H>
#include <FL/Fl_Value_Input.H>
#include <FL/Fl_Return_Button.H>
#include <FL/Fl_Check_Button.H>
#include <fstream>
#include "vector.h"
#include <algorithm>
#include <fcntl.h>
#include <cmath>
#include <time.h> // for nanosleep
#include <stdio.h>
#include <dirent.h>
#include <sys/stat.h> // auto-recovery function
#include "path.h"
#include "misc.h"
#include "Equation.h"
#include "Editeur.h"
#include "Xcas1.h"
#include "Cfg.h"
#include "Print.h"
#include "gl2ps.h"

#ifdef __APPLE__
//#include <OpenGL/gl.h>
//#include <AGL/agl.h>
#include <FL/x.H>
#include <FL/fl_draw.H>
#define __APPLE_QUARTZ__ 1
#include "Fl_Gl_Choice.H"
extern Fl_Gl_Choice * gl_choice;
// be sure to remove static declaration in gl_start.cxx in fltk src directory
#endif

#ifndef HAVE_PNG_H
#undef HAVE_LIBPNG
#endif
#ifdef HAVE_LIBPNG
#include <png.h>
#endif

using namespace std;
using namespace giac;

#ifndef NO_NAMESPACE_XCAS
namespace xcas {
#endif // ndef NO_NAMESPACE_XCAS

  std::map<std::string,Fl_Image *> texture_cache;

  int Graph3d::opengl2png(const std::string & filename){
    if (!opengl){
      fl_alert(gettext("Supports only OpenGL rendering"));
      return -1;
    }
#ifdef HAVE_LIBPNG
    if (!screenbuf)
      return -1;
    int i;
    // unsigned rowbytes = w()*4;
    unsigned char *rows[h()];
    for (i = 0; i < h(); i++) {
      rows[i] = &screenbuf[(h() - i - 1)*4*w()];
    }
    int res= write_png(filename.c_str(), rows, w(), h(), PNG_COLOR_TYPE_RGBA, 8);
    if (res!=-1){
      string command="pngtopnm "+filename+" | pnmtops > "+remove_extension(filename)+".ps &";
      cerr << command << '\n';
      system_no_deprecation(command.c_str()); 
      command="pngtopnm "+filename+" | pnmtojpeg > "+remove_extension(filename)+".jpg &";
      cerr << command << '\n';
      system_no_deprecation(command.c_str()); 
    }
#endif
    return -1;
  }


  double giac_max(double i,double j){
    return i>j?i:j;
  }

  quaternion_double rotation_2_quaternion_double(double x, double y, double z,double theta){
    double t=theta*M_PI/180;
    double qx,qy,qz,qw,s=std::sin(t/2),c=std::cos(t/2);
    qx=x*s;
    qy=y*s;
    qz=z*s;
    qw=c;
    double n=std::sqrt(qx*qx+qy*qy+qz*qz+qw*qw);
    return quaternion_double(qw/n,qx/n,qy/n,qz/n);
  }

  void Graph3d::resize(int X,int Y,int W,int H){
    int oldh=h(),oldw=w();
    Graph2d3d::resize(X,Y,W,H);
    if (screenbuf){
      if (oldh==H && oldw==W)
	return;
      delete screenbuf;
    }
    screenbuf = new unsigned char[H*W*4];
  }

  Graph3d::Graph3d(int x,int y,int width, int height,const char* title,History_Pack * hp_): 
    Graph2d3d(x,y,width, height, title,hp_),
    theta_z(-110),theta_x(-13),theta_y(-95),
    delta_theta(5),draw_mode(GL_QUADS),printing(0),glcontext(0),dragi(0),dragj(0),push_in_area(false),depth(0),below_depth_hidden(false),screenbuf(0) {
    current_depth=0;
    displayhint=true;
#if 1 // ndef WIN32 // def HYPERSURFACE3
    if (gnuplot_opengl)
      opengl=true;
    else {
      opengl=false;
      theta_z=theta_x=theta_y=0;
    }
#else
    opengl=true;
#endif
    q=euler_deg_to_quaternion_double(theta_z,theta_x,theta_y);
    // end();
    // mode=0;
    display_mode |= 0x80;
    display_mode |= 0x200;
    box(FL_FLAT_BOX);
    couleur=_POINT_WIDTH_5;
    legende_size=max(min(legende_size,width/4),width/6);
    resize(x,y,width-legende_size,height);
    add_mouse_param_group(x,y,width,height);
    // 8 light initialization
    for (int i=0;i<8;++i)
      reset_light(i);
    if (!opengl){
      diffusionz=5; diffusionz_limit=5; hide2nd=false; interval=false;
      default_upcolor=giac3d_default_upcolor;
      default_downcolor=giac3d_default_downcolor;
      default_downupcolor=giac3d_default_downupcolor;
      default_downdowncolor=giac3d_default_downdowncolor;
      doprecise=false;
      lcdz= h()/4;
      autoscale(false);
      update_rotation();
      precision=1;
    }

  }

  // t angle in radians -> r,g,b
  void arc_en_ciel2(double t,int & r,int & g,int &b){
    int k=int(t/2/M_PI*126);
    arc_en_ciel(k,r,g,b);
  }

  bool get_glvertex(const vecteur & v,double & d1,double & d2,double & d3,double realiscmplx,double zmin,GIAC_CONTEXT){
    if (v.size()==3){
      gen tmp;
      tmp=evalf_double(v[0],2,contextptr);
      if (tmp.type!=_DOUBLE_) return false;
      d1=tmp._DOUBLE_val;
      tmp=evalf_double(v[1],2,contextptr);
      if (tmp.type!=_DOUBLE_) return false;
      d2=tmp._DOUBLE_val;
      tmp=evalf_double(v[2],2,contextptr);
      if (realiscmplx){
	double arg=0;
	if (realiscmplx<0){
	  d3=evalf_double(im(tmp,contextptr),2,contextptr)._DOUBLE_val-zmin;
	  arg=-d3*realiscmplx;
	}
	else {
	  if (tmp.type==_DOUBLE_){
	    d3=tmp._DOUBLE_val;
	    if (d3<0){
	      arg=M_PI;
	      d3=-d3;
	    }
	  }
	  else {
	    if (tmp.type==_CPLX && tmp._CPLXptr->type==_DOUBLE_ && (tmp._CPLXptr+1)->type==_DOUBLE_){
	      double r=tmp._CPLXptr->_DOUBLE_val;
	      double i=(tmp._CPLXptr+1)->_DOUBLE_val;
	      arg=std::atan2(i,r);
	      d3=std::sqrt(r*r+i*i);
	    }
	    else 
	      return false;
	  }
	} // end realiscmplx>0
	// set color corresponding to argument
	int r,g,b;
	arc_en_ciel2(arg,r,g,b);
	glColor3f(r/255.,g/255.,b/255.);
	// glColor4i(r,g,b,int(std::log(d3+1)));
      } // end if (realiscmplx)
      else {
	if (tmp.type!=_DOUBLE_) return false;
	d3=tmp._DOUBLE_val;
      }
      return true;
    }
    return false;
  }

  bool get_glvertex(const gen & g,double & d1,double & d2,double & d3,double realiscmplx,double zmin,GIAC_CONTEXT){
    if (g.type!=_VECT)
      return false;
    return get_glvertex(*g._VECTptr,d1,d2,d3,realiscmplx,zmin,contextptr);
  }

  bool glvertex(const vecteur & v,double realiscmplx,double zmin,GIAC_CONTEXT){
    double d1,d2,d3;
    if (get_glvertex(v,d1,d2,d3,realiscmplx,zmin,contextptr)){
      glVertex3d(d1,d2,d3);
      return true;
    }
    return false;
  }

  void glraster(const vecteur & v){
    if (v.size()==3){
      double d1=evalf_double(v[0],2,context0)._DOUBLE_val;
      double d2=evalf_double(v[1],2,context0)._DOUBLE_val;
      double d3=evalf_double(v[2],2,context0)._DOUBLE_val;
      glRasterPos3d(d1,d2,d3);
    }
  }

  void glraster(const gen & g){
    if (g.type==_VECT)
      glraster(*g._VECTptr);
  }

  // draw s at g with mode= 0 (upper right), 1, 2 or 3
  void Graph3d::legende_draw(const gen & g,const string & s,int mode){
    context * contextptr=hp?hp->contextptr:get_context(this);
    gen gf=evalf_double(g,1,contextptr);
    if (gf.type==_VECT && gf._VECTptr->size()==3){
      double Ax=gf[0]._DOUBLE_val;
      double Ay=gf[1]._DOUBLE_val;
      double Az=gf[2]._DOUBLE_val;
      double Ai,Aj,Ad;
      find_ij(Ax,Ay,Az,Ai,Aj,Ad);
      int di=3,dj=1;
      find_dxdy(s,mode,labelsize(),di,dj);
      find_xyz(Ai+di,Aj+dj,Ad,Ax,Ay,Az);
      glRasterPos3d(Ax,Ay,Az);
      // string s1(s);
      // cst_greek_translate(s1);
      // draw_string(s1);
      draw_string(s);
    }
  }

  void glnormal(const vecteur & v){
    if (v.size()==3){
      double d1=evalf_double(v[0],2,context0)._DOUBLE_val;
      double d2=evalf_double(v[1],2,context0)._DOUBLE_val;
      double d3=evalf_double(v[2],2,context0)._DOUBLE_val;
      glNormal3d(d1,d2,d3);
    }
  }

  void glnormal(double d1,double d2,double d3,double e1,double e2,double e3,double f1,double f2,double f3){
    double de1(e1-d1),de2(e2-d2),de3(e3-d3),df1(f1-d1),df2(f2-d2),df3(f3-d3);
    glNormal3d(de2*df3-de3*df2,de3*df1-de1*df3,de1*df2-de2*df1);
  }

  void gltranslate(const vecteur & v){
    if (v.size()==3){
      double d1=evalf_double(v[0],2,context0)._DOUBLE_val;
      double d2=evalf_double(v[1],2,context0)._DOUBLE_val;
      double d3=evalf_double(v[2],2,context0)._DOUBLE_val;
      glTranslated(d1,d2,d3);
    }
  }


  void iso3d(const double &i1,const double &i2,const double &i3,const double &j1,const double &j2,const double &j3,const double &k1,const double &k2,const double &k3,double & x,double & y,double & z){
    double X=x,Y=y,Z=z;
    x=i1*X+j1*Y+k1*Z;
    y=i2*X+j2*Y+k2*Z;
    z=i3*X+j3*Y+k3*Z;
  }


#define LI 64
#define LH 64
  GLubyte image[LI][LH][3];

  void makeImage(void) {
    int i,j,c;
    for( i = 0 ; i < LI ; i++ ) {
      for( j = 0 ; j < LH ; j++ ) {
	c = (((i&0x8)==0)^
	     ((j&0x8)==0))*255;
	image[i][j][0] =(GLubyte) c;
	image[i][j][1] =(GLubyte) c;
	image[i][j][2] =(GLubyte) c; } }
  }

  bool test_enable_texture(Fl_Image * texture){
    if (!texture)
      return false;
    int depth=-1;
    if (texture->count()==1)
      depth=texture->d();
    if (depth==3 || depth==4){
      // define texture
      // makeImage();
      char * ptr=(char *)texture->data()[0];
      /*
      int W=texture->w(),H=texture->h();
      for (int y=0;y<H;++y){
	for (int x=0;x<W;x++)
	  cerr << unsigned(*(ptr+(x+y*W)*depth)) << " " << unsigned(*(ptr+(x+y*W)*depth+1)) << " " << unsigned(*(ptr+(x+y*W)*depth+2)) << ", ";
	cerr << '\n';
      }
      */
      // texture->w() and texture->h() must be a power of 2!!
      glTexImage2D(GL_TEXTURE_2D,0,depth,
		   // LI,LH,
		   texture->w(),texture->h(),
		   0,
		   // GL_RGB,
		   (depth==3?GL_RGB:GL_RGBA),
		   GL_UNSIGNED_BYTE,
		   // &image[0][0][0]);
		   ptr);
      // not periodically
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP); 
      // adjust to nearest 
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,
		      // GL_NEAREST);
		      GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,
		      // GL_NEAREST);
		      GL_LINEAR);
      //
      glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,
		GL_MODULATE);
      //GL_REPLACE);
      glEnable(GL_TEXTURE_2D);
      return true;
    }
    return false;
  }

  // Sphere centered at center, radius radius, i,j,k orthonormal, ntheta/nphi
  // number of subdivisions
  // mode=GL_QUADS for example
  void glsphere(const vecteur & center,const gen & radius,const vecteur & i0,const vecteur & j0,const vecteur & k0,int ntheta,int nphi,int mode,Fl_Image * texture,GIAC_CONTEXT){
    test_enable_texture(texture);
    double c1=evalf_double(center[0],1,contextptr)._DOUBLE_val; // center
    double c2=evalf_double(center[1],1,contextptr)._DOUBLE_val;
    double c3=evalf_double(center[2],1,contextptr)._DOUBLE_val;
    double r=evalf_double(radius,1,contextptr)._DOUBLE_val;
    double i1=evalf_double(i0[0],1,contextptr)._DOUBLE_val;
    double i2=evalf_double(i0[1],1,contextptr)._DOUBLE_val;
    double i3=evalf_double(i0[2],1,contextptr)._DOUBLE_val;
    double j1=evalf_double(j0[0],1,contextptr)._DOUBLE_val;
    double j2=evalf_double(j0[1],1,contextptr)._DOUBLE_val;
    double j3=evalf_double(j0[2],1,contextptr)._DOUBLE_val;
    double k1=evalf_double(k0[0],1,contextptr)._DOUBLE_val;
    double k2=evalf_double(k0[1],1,contextptr)._DOUBLE_val;
    double k3=evalf_double(k0[2],1,contextptr)._DOUBLE_val;
    double dtheta=2*M_PI/ntheta; // longitude
    double dphi=M_PI/nphi; // latitude
    double jsurtheta=0,isurphi=0,djsurtheta=1.0/ntheta,disurphi=1.0/nphi;
    double parallele1[ntheta+1],parallele2[ntheta+1],parallele3[ntheta+1];
    double x,y,z,oldx,oldy,oldz,X,Z;
    // Set initial parallel to the North pole
    for (int j=0;j<=ntheta;++j){
      parallele1[j]=k1; parallele2[j]=k2; parallele3[j]=k3;
    }
    // parallele1/2/3 contains the coordinate of the previous parallel
    for (int i=1;i<=nphi;i++,isurphi+=disurphi){
      // longitude theta=0, latitude phi=i*dphi
      X=std::sin(i*dphi); Z=std::cos(i*dphi); 
      oldx=X; oldy=0; oldz=Z;
      iso3d(i1,i2,i3,j1,j2,j3,k1,k2,k3,oldx,oldy,oldz);
      double * par1j=parallele1,* par2j=parallele2,*par3j=parallele3;
      jsurtheta=0;
      if (mode==GL_QUADS){
	glBegin(GL_QUAD_STRIP);
	for (int j=0;j<=ntheta;++par1j,++par2j,++par3j,jsurtheta+=djsurtheta){
	  glNormal3d(*par1j,*par2j,*par3j);
	  if (texture)
	    glTexCoord2f(jsurtheta,isurphi);
	  glVertex3d(c1+r*(*par1j),c2+r*(*par2j),c3+r*(*par3j));
	  if (texture)
	    glTexCoord2f(jsurtheta,isurphi+disurphi);
	  glVertex3d(c1+r*oldx,c2+r*oldy,c3+r*oldz);
	  *par1j=oldx; *par2j=oldy; *par3j=oldz;
	  // theta=j*dtheta, phi=i*dphi
	  ++j;
	  x=X*std::cos(j*dtheta); y=X*std::sin(j*dtheta); z=Z;
	  iso3d(i1,i2,i3,j1,j2,j3,k1,k2,k3,x,y,z);
	  oldx=x; oldy=y; oldz=z;
	}
	glEnd();
      }
      else 
      {
	for (int j=0;j<ntheta;){
	  glBegin(mode);
	  glNormal3d(parallele1[j],parallele2[j],parallele3[j]);
	  if (texture)
	    glTexCoord2f(double(j)/(ntheta),double(i-1)/(nphi));
	  glVertex3d(c1+r*parallele1[j],c2+r*parallele2[j],c3+r*parallele3[j]);
	  if (texture)
	    glTexCoord2f(double(j)/(ntheta),double(i)/(nphi));
	  glVertex3d(c1+r*oldx,c2+r*oldy,c3+r*oldz);
	  parallele1[j]=oldx; parallele2[j]=oldy; parallele3[j]=oldz;
	  ++j;
	  // theta=j*dtheta, phi=i*dphi
	  x=X*std::cos(j*dtheta); y=X*std::sin(j*dtheta); z=Z;
	  iso3d(i1,i2,i3,j1,j2,j3,k1,k2,k3,x,y,z);
	  if (texture)
	    glTexCoord2f(double(j+1)/(ntheta),double(i-1)/(nphi));
	  glVertex3d(c1+r*x,c2+r*y,c3+r*z);
	  if (texture)
	    glTexCoord2f(double(j+1)/(ntheta),double(i)/(nphi));
	  glVertex3d(c1+r*parallele1[j],c2+r*parallele2[j],c3+r*parallele3[j]);
	  glEnd();
	  oldx=x; oldy=y; oldz=z;
	}
      }
      parallele1[ntheta]=oldx; parallele2[ntheta]=oldy; parallele3[ntheta]=oldz;
    }
    if (texture)
      glDisable(GL_TEXTURE_2D);
  }

  // if all z values of surfaceg are pure imaginary return a negative int
  // else return 1
  bool find_zscale(const gen & surface,double & zmin,double & zmax){
    if (surface.type!=_VECT)
      return true;
    if (surface.subtype!=_POLYEDRE__VECT && surface.subtype!=_GROUP__VECT && surface._VECTptr->size()==3){
      if (!is_zero(re(surface._VECTptr->back(),context0)))
	return false;
      gen s3=evalf_double(surface._VECTptr->back()/cst_i,2,context0);
      if (s3.type==_DOUBLE_){
	double s3d = s3._DOUBLE_val;
	if (s3d<zmin) zmin=s3d;
	if (s3d>zmax) zmax=s3d;
      }
      return true;
    }
    const_iterateur it=surface._VECTptr->begin(),itend=surface._VECTptr->end();
    for (;it!=itend;++it){
      if (!find_zscale(*it,zmin,zmax))
	return false;
    }
    return true;
  }
  
  struct float2 {
    float f,a;
  } ;
  double absarg(const gen & g,double & argcolor){
    if (g.type==_DOUBLE_){
      double d=g._DOUBLE_val;
      if (d>=0){ argcolor=0;  return d; }
      argcolor=-M_PI; return -d;
    }
    double x=g._CPLXptr->_DOUBLE_val,y=(g._CPLXptr+1)->_DOUBLE_val;
    argcolor=std::atan2(x,y);
    double n=std::sqrt(x*x+y*y); // will be encoded in a float, no overflow care
    return n;
  }

  void glsurface(const gen & surfaceg,int draw_mode,Fl_Image * texture,GIAC_CONTEXT){
    if (!ckmatrix(surfaceg,true))
      return;
    test_enable_texture(texture);
    double realiscmplx=has_i(surfaceg),zmin=1e300,zmax=-1e300;
    matrice & surface = *surfaceg._VECTptr;
    if (realiscmplx && !texture){
      if (find_zscale(surface,zmin,zmax))
	realiscmplx=2*M_PI/(zmin-zmax);
    }
    int n=surface.size();
    if (surfaceg.subtype==_POLYEDRE__VECT){
      // implicit surface drawing with the given triangulation
      if (draw_mode==GL_QUADS || draw_mode==GL_TRIANGLES){
	glBegin(GL_TRIANGLES);
	for (int i=0;i<n;++i){
	  vecteur & v=*surface[i]._VECTptr;
	  const_iterateur it=v.begin();
	  double d1,d2,d3,e1,e2,e3,f1,f2,f3;
	  if (get_glvertex(*it,d1,d2,d3,realiscmplx,zmin,contextptr) &&
	      get_glvertex(*(it+1),e1,e2,e3,realiscmplx,zmin,contextptr) &&
	      get_glvertex(*(it+2),f1,f2,f3,realiscmplx,zmin,contextptr)){
	    glnormal(d1,d2,d3,e1,e2,e3,f1,f2,f3);
	    glVertex3d(d1,d2,d3);
	    glVertex3d(e1,e2,e3);
	    glVertex3d(f1,f2,f3);
	  }
	}
	glEnd();
      }
      else {
	glBegin(draw_mode);
	for (int i=0;i<n;++i){
	  vecteur & v=*surface[i]._VECTptr;
	  const_iterateur it=v.begin();
	  double d1,d2,d3,e1,e2,e3,f1,f2,f3;
	  if (get_glvertex(*it,d1,d2,d3,realiscmplx,zmin,contextptr) &&
	      get_glvertex(*(it+1),e1,e2,e3,realiscmplx,zmin,contextptr) &&
	      get_glvertex(*(it+2),f1,f2,f3,realiscmplx,zmin,contextptr)){
	    glnormal(d1,d2,d3,e1,e2,e3,f1,f2,f3);
	    glVertex3d(d1,d2,d3);
	    glVertex3d(e1,e2,e3);
	    glVertex3d(f1,f2,f3);
	    glVertex3d(d1,d2,d3);	    
	  }
	}
	glEnd();
      }
      if (texture)
	glDisable(GL_TEXTURE_2D);
      return;
    }
    if (surface.front()._VECTptr->size()<2){
      if (texture)
	glDisable(GL_TEXTURE_2D);
      return;
    }
    gen a,b,c,d;
    double xt=0,yt,dxt=double(1)/n,dyt;
    for (int i=1;i<n;++i,xt+=dxt){
      const vecteur & vprec=*surface[i-1]._VECTptr;
      const vecteur & v=*surface[i]._VECTptr;
      const_iterateur itprec=vprec.begin(),it=v.begin(),itend=v.end();
      yt=0; dyt=double(1)/(itend-it);
      a=*itprec; b=*it;
      double a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3;
      /*
      if (draw_mode==GL_QUADS){
	get_glvertex(b,b1,b2,b3,realiscmplx,zmin,contextptr);
	glBegin(GL_QUAD_STRIP);
	for (;it!=itend;yt+=dyt){
	  a=*itprec;
	  get_glvertex(a,a1,a2,a3,realiscmplx,zmin,contextptr);
	  ++it; ++itprec;
	  if (it==itend){
	    if (texture)
	      glTexCoord2f(xt,yt);
	    glVertex3f(a1,a2,a3);
	    glVertex3f(b1,b2,b3);
	    break;
	  }
	  b=*it;
	  get_glvertex(b,c1,c2,c3,realiscmplx,zmin,contextptr);
	  if (texture)
	    glTexCoord2f(xt,yt);
	  glnormal(a1,a2,a3,b1,b2,b3,c1,c2,c3);
	  glVertex3f(a1,a2,a3);
	  glVertex3f(b1,b2,b3);
	  b1=c1; b2=c2; b3=c3;
	}
	glEnd();
      }
      else 
      */
      {
	get_glvertex(a,a1,a2,a3,realiscmplx,zmin,contextptr);
	get_glvertex(b,b1,b2,b3,realiscmplx,zmin,contextptr);
	++it; ++itprec;
	for (;it!=itend;++it,++itprec,yt+=dyt){
	  c = *itprec;
	  d = *it;
	  get_glvertex(c,c1,c2,c3,realiscmplx,zmin,contextptr);
	  get_glvertex(d,d1,d2,d3,realiscmplx,zmin,contextptr);
	  glBegin(draw_mode);
	  if (texture)
	    glTexCoord2f(xt,yt);
	  glnormal(a1,a2,a3,b1,b2,b3,c1,c2,c3);
	  glVertex3f(a1,a2,a3); // itprec
	  glVertex3f(b1,b2,b3); // it
	  glVertex3f(d1,d2,d3); // it
	  glVertex3f(c1,c2,c3); // itprec
	  glEnd();
	  if (draw_mode==GL_QUADS){
	    glnormal(a1,a2,a3,c1,c2,c3,b1,b2,b3);
	    glBegin(draw_mode);
	    glVertex3f(a1,a2,a3);
	    glVertex3f(c1,c2,c3);
	    glVertex3f(d1,d2,d3);
	    glVertex3f(b1,b2,b3);
	    glEnd();
	  }
	  a=c; a1=c1; a2=c2; a3=c3;
	  b=d; b1=d1; b2=d2; b3=d3;
	}
      }
    }
    if (texture)
      glDisable(GL_TEXTURE_2D);
  }

  // surface without grid evaluation, should not happen!
  void glsurface(const vecteur & point,const gen & uv,double umin,double umax,double vmin,double vmax,int nu,int nv,int draw_mode,GIAC_CONTEXT){
    double u=umin,v=vmin,deltau=(umax-umin)/nu,deltav=(vmax-vmin)/nv;
    vecteur prevline(nv+1); //gen prevline[nv+1];//,line[nv+1];
    vecteur curuv(2);
    gen old,current;
    curuv[0]=u;
    for (int j=0;j<=nv;j++,v += deltav){
      curuv[1] = v;
      prevline[j]=subst(point,uv,curuv,false,contextptr);
    }
    u += deltau;
    for (int i=1;i<=nu;i++,u+=deltau){
      v=vmin;
      curuv = makevecteur(u,v);
      v += deltav;
      old = subst(point,uv,curuv,false,contextptr);
      for (int j=0;j<nv;v +=deltav){
	glBegin(draw_mode);
	curuv[1] = v;
	current = subst(point,uv,curuv,false,contextptr);
	glvertex(*prevline[j]._VECTptr,0,0,contextptr);
	prevline[j]=old;
	glvertex(*old._VECTptr,0,0,contextptr);
	++j;
	glvertex(*current._VECTptr,0,0,contextptr);
	glvertex(*prevline[j]._VECTptr,0,0,contextptr);
	old=current;
	glEnd();
      }
      prevline[nv]=old;
      /* for (int j=0;j<=nv;++j)
	 cerr << prevline[j] << " " ;
	cerr << '\n';
      */
    }
  }

  unsigned int line_stipple(unsigned int i){
    switch (i){
    case 1: case 4:
      return 0xf0f0;
    case 2: case 5:
      return 0xcccc;
    case 3: case 6:
      return 0xaaaa;
    default:
      return 0xffff;
    }
  }

  bool is_approx_zero(const gen & dP,double window_xmin,double window_xmax,double window_ymin,double window_ymax,double window_zmin,double window_zmax){
    bool closed=false;
    if (dP.type==_VECT && dP._VECTptr->size()==3){
      closed=true;
      double dPx=evalf_double(dP[0],2,context0)._DOUBLE_val;
      if (fabs(dPx)>(window_xmax-window_xmin)*1e-6)
	closed=false;
      double dPy=evalf_double(dP[1],2,context0)._DOUBLE_val;
      if (fabs(dPy)>(window_ymax-window_ymin)*1e-6)
	closed=false;
      double dPz=evalf_double(dP[2],2,context0)._DOUBLE_val;
      if (fabs(dPz)>(window_zmax-window_zmin)*1e-6)
	closed=false;
    }
    return closed;
  }

  // translate giac GL constant to open GL constant
  unsigned gl_translate(unsigned i){
    switch (i){
    case _GL_LIGHT0:
      return GL_LIGHT0;
    case _GL_LIGHT1:
      return GL_LIGHT1;
    case _GL_LIGHT2:
      return GL_LIGHT2;
    case _GL_LIGHT3:
      return GL_LIGHT3;
    case _GL_LIGHT4:
      return GL_LIGHT4;
    case _GL_LIGHT5:
      return GL_LIGHT5;
    case _GL_AMBIENT:
      return GL_AMBIENT;
    case _GL_SPECULAR:
      return GL_SPECULAR;
    case _GL_DIFFUSE:
      return GL_DIFFUSE;
    case _GL_POSITION:
      return GL_POSITION;
    case _GL_SPOT_DIRECTION:
      return GL_SPOT_DIRECTION;
    case _GL_SPOT_EXPONENT:
      return GL_SPOT_EXPONENT;
    case _GL_SPOT_CUTOFF:
      return GL_SPOT_CUTOFF;
    case _GL_CONSTANT_ATTENUATION:
      return GL_CONSTANT_ATTENUATION;
    case _GL_LINEAR_ATTENUATION:
      return GL_LINEAR_ATTENUATION;
    case _GL_QUADRATIC_ATTENUATION:
      return GL_QUADRATIC_ATTENUATION;
    case _GL_LIGHT_MODEL_AMBIENT:
      return GL_LIGHT_MODEL_AMBIENT;
    case _GL_LIGHT_MODEL_LOCAL_VIEWER:
      return GL_LIGHT_MODEL_LOCAL_VIEWER;
    case _GL_LIGHT_MODEL_TWO_SIDE:
      return GL_LIGHT_MODEL_TWO_SIDE;
#ifndef WIN32
    case _GL_LIGHT_MODEL_COLOR_CONTROL:
      return GL_LIGHT_MODEL_COLOR_CONTROL;
#endif
    case _GL_SMOOTH:
      return GL_SMOOTH;
    case _GL_FLAT:
      return GL_FLAT;
    case _GL_SHININESS:
      return GL_SHININESS;
    case _GL_FRONT:
      return GL_FRONT;
    case _GL_BACK:
      return GL_BACK;
    case _GL_FRONT_AND_BACK:
      return GL_FRONT_AND_BACK;
    case _GL_AMBIENT_AND_DIFFUSE:
      return GL_AMBIENT_AND_DIFFUSE;
    case _GL_EMISSION:
      return GL_EMISSION;
#ifndef WIN32
    case _GL_SEPARATE_SPECULAR_COLOR:
      return GL_SEPARATE_SPECULAR_COLOR;
    case _GL_SINGLE_COLOR:
      return GL_SINGLE_COLOR;
#endif
    case _GL_BLEND:
      return GL_BLEND;
    case _GL_SRC_ALPHA:
      return GL_SRC_ALPHA;
    case _GL_ONE_MINUS_SRC_ALPHA:
      return GL_ONE_MINUS_SRC_ALPHA;
    case _GL_COLOR_INDEXES:
      return GL_COLOR_INDEXES;
    }
    cerr << "No GL equivalent for " << i << '\n';
    return i;
  }

  void tran4(double * colmat){
    giac::swapdouble(colmat[1],colmat[4]);
    giac::swapdouble(colmat[2],colmat[8]);
    giac::swapdouble(colmat[3],colmat[12]);
    giac::swapdouble(colmat[6],colmat[9]);
    giac::swapdouble(colmat[7],colmat[13]);
    giac::swapdouble(colmat[11],colmat[14]);    
  }

  void get_texture(const gen & attrv1,Fl_Image * & texture){
    // set texture
    if (attrv1.type==_STRNG){
      std::map<std::string,Fl_Image *>::const_iterator it,itend=texture_cache.end();
      it=texture_cache.find(*attrv1._STRNGptr);
      if (it!=itend){
	texture=it->second;
	// texture->uncache();
      }
      else {
	texture=Fl_Shared_Image::get(attrv1._STRNGptr->c_str());
	if (texture){
	  int W=texture->w(),H=texture->h();
	  // take a power of 2 near w/h
	  W=1 << min(int(std::log(double(W))/std::log(2.)+.5),8);
	  H=1 << min(int(std::log(double(H))/std::log(2.)+.5),8);
	  texture=texture->copy(W,H);
	  texture_cache[*attrv1._STRNGptr]=texture;
	}
      }
    }
  }

  int gen2int(const gen & g){
    if (g.type==_INT_)
      return g.val;
    if (g.type==_DOUBLE_)
      return int(g._DOUBLE_val);
    setsizeerr(gettext("Unable to convert to int")+g.print());
    return -1;
  }

  void Graph3d::indraw(const giac::gen & g){
    context * contextptr=hp?hp->contextptr:get_context(this);
    if (g.type==_VECT)
      indraw(*g._VECTptr);
    if (g.is_symb_of_sommet(at_animation)){
      indraw(get_animation_pnt(g,animation_instructions_pos));
      return;
    }
    if (!g.is_symb_of_sommet(at_pnt)){
      update_infos(g,contextptr);
      return;
    }
    gen & f=g._SYMBptr->feuille;
    if (f.type!=_VECT)
      return;
    vecteur & v = *f._VECTptr;
    gen v0=v[0];
    bool est_hyperplan=v0.is_symb_of_sommet(at_hyperplan);
    string legende;
    vecteur style(get_style(v,legende));
    int styles=style.size();
    // color
    bool hidden_name = false;
    int ensemble_attributs=style.front().val;
    if (style.front().type==_ZINT){
      ensemble_attributs = mpz_get_si(*style.front()._ZINTptr);
      hidden_name=true;
    }
    else
      hidden_name=ensemble_attributs<0;
    int couleur=ensemble_attributs & 0x0000ff;
    int width           =(ensemble_attributs & 0x00070000) >> 16; // 3 bits
    int epaisseur_point =((ensemble_attributs & 0x00380000) >> 19)+1; // 3 bits
    int type_line       =(ensemble_attributs & 0x01c00000) >> 22; // 3 bits
    int type_point      =(ensemble_attributs & 0x0e000000) >> 25; // 3 bits
    int labelpos        =(ensemble_attributs & 0x30000000) >> 28; // 2 bits
    bool fill_polygon   =(ensemble_attributs & 0x40000000) >> 30;
    hidden_name = hidden_name || legende.empty();
    Fl_Image * texture=0;
    glLineWidth(width+1);
    gl2psLineWidth(width);
    // FIXME line_stipple disabled because of printing
    if (!printing)
      glLineStipple(1,line_stipple(type_line));
    else
      glLineStipple(1,0xffff);
    glPointSize(epaisseur_point);
    gl2psPointSize(epaisseur_point);
    if (styles<=2){
      glEnable(GL_COLOR_MATERIAL);
      glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);
    }
    else 
      glDisable(GL_COLOR_MATERIAL);
    GLfloat tab[4]={0,0,0,1};
    glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,tab);
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,50);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,tab);
    GLfloat tab1[4]={0.2,0.2,0.2,1};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT,tab1);
    GLfloat tab2[4]={0.8,0.8,0.8,1};
    glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE,tab2);
    if (est_hyperplan){
      glPolygonMode(GL_FRONT_AND_BACK,fill_polygon?GL_FILL:GL_LINE);
    }
    else {
      glPolygonMode(GL_FRONT_AND_BACK,fill_polygon?GL_FILL:GL_LINE);
      // glPolygonMode(GL_BACK,GL_POINT);
      // glMaterialfv(GL_BACK,GL_EMISSION,tab);
      // glMaterialfv(GL_BACK,GL_SPECULAR,tab);
      // glMaterialfv(GL_BACK,GL_AMBIENT_AND_DIFFUSE,tab);
    }
    // glGetFloatv(GL_CURRENT_COLOR,tab);

    // set other style attributs : 
    /* material = [gl_front|gl_back|gl_front_and_back,gl_shininess,valeur] or
       material = [gl_front|gl_back|gl_front_and_back, GL_AMBIENT| GL_DIFFUSE |
       GL_SPECULAR | GL_EMISSION | GL_AMBIENT_AND_DIFFUSE | GL_COLOR_INDEXES,
      [r,g,b,a] ]
      last arg is [ambient,diffuse,specular] for GL_COLOR_INDEXES
      or
      material=[gl_texture,"image_filename",...]
    */
    for (int i=2;i<styles;++i){
      gen & attr=style[i];
      if (attr.type==_VECT){
	vecteur & attrv=*attr._VECTptr;
	if (attrv.size()>=2 && attrv.front().type==_INT_ && attrv.front().val==_GL_TEXTURE){
	  get_texture(attrv[1],texture);
	  continue;
	}
	if (attrv.size()==2 && attrv.front().type==_INT_ && attrv.front().val==_GL_MATERIAL){
	  gen attrm =evalf_double(attrv.back(),1,contextptr);
	  if (debug_infolevel)
	    cerr << "Setting material " << attrm << '\n';
	  if (attrm.type==_VECT && attrm._VECTptr->size()<=3 ){
	    gen attrv0=attrv.back()._VECTptr->front();
	    if (attrv0.type==_INT_ && attrv0.val==_GL_TEXTURE){
	      gen attrv1=(*attrv.back()._VECTptr)[1];
	      get_texture(attrv1,texture);
	      continue;
	    }
	    vecteur & attrmv = *attrm._VECTptr;
	    if (attrmv.back().type==_VECT && attrmv.back()._VECTptr->size()==4){
	      vecteur & w=*attrmv.back()._VECTptr;
	      GLfloat tab[4]={w[0]._DOUBLE_val,w[1]._DOUBLE_val,w[2]._DOUBLE_val,w[3]._DOUBLE_val};
	      glMaterialfv(gl_translate(gen2int(attrmv[0])),gl_translate(gen2int(attrmv[1])),tab);
	      continue;
	    }
	    if (attrmv.back().type==_VECT && attrmv.back()._VECTptr->size()==3){
	      vecteur & w=*attrmv.back()._VECTptr;
	      GLfloat tab[3]={w[0]._DOUBLE_val,w[1]._DOUBLE_val,w[2]._DOUBLE_val};
	      glMaterialfv(gl_translate(gen2int(attrmv[0])),GL_COLOR_INDEXES,tab);
	      continue;
	    }
	    if (attrmv.back().type==_DOUBLE_ && attrmv[1]._DOUBLE_val==_GL_SHININESS){
	      // glMaterialf(gl_translate(gen2int(attrmv[0])),GL_SHININESS,float(attrmv[2]._DOUBLE_val));
	      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,float(attrmv[2]._DOUBLE_val));
	      continue;
	    }
	  }
	}
      }
    }
    if (texture){
      fill_polygon=true;
      couleur=FL_WHITE;
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    }
    bool hidden_line = fill_polygon && (width==7 || (display_mode & 0x8) || texture );
    xcas_color(couleur,true);
    if (debug_infolevel){
      cerr << "opengl displaying " << g << '\n';
      GLint b;
      GLfloat posf[4],direcf[4],ambient[4],diffuse[4],specular[4],emission[4],shini[1];
      GLfloat expo,cutoff;
      double pos[4],direc[4];
      glGetIntegerv(GL_BLEND,&b);
      cerr << "blend " << b << '\n';
      for (int i=0;i<8;++i){
	glGetIntegerv(GL_LIGHT0+i,&b);
	if (b){
	  glGetLightfv(GL_LIGHT0+i,GL_SPOT_EXPONENT,&expo);
	  glGetLightfv(GL_LIGHT0+i,GL_SPOT_CUTOFF,&cutoff);
	  glGetLightfv(GL_LIGHT0+i,GL_POSITION,posf);
	  glGetLightfv(GL_LIGHT0+i,GL_SPOT_DIRECTION,direcf);
	  direcf[3]=0; // 4 was out of bound
	  mult4(model_inv,posf,pos);
	  tran4(model);
	  mult4(model,direcf,direc);
	  tran4(model);
	  glGetLightfv(GL_LIGHT0+i,GL_AMBIENT,ambient);
	  glGetLightfv(GL_LIGHT0+i,GL_DIFFUSE,diffuse);
	  glGetLightfv(GL_LIGHT0+i,GL_SPECULAR,specular);
	  cerr << "light " << i << ": " <<
	    " pos " << pos[0] << "," << pos[1] << "," << pos[2] << "," << pos[3] << 
	    " dir " << direc[0] << "," << direc[1] << "," << direc[2] << "," << direc[3] << 
	    " ambient " << ambient[0] << "," << ambient[1] << "," << ambient[2] << "," << ambient[3] << 
	    " diffuse " << diffuse[0] << "," << diffuse[1] << "," << diffuse[2] << "," << diffuse[3] << 
	    " specular " << specular[0] << "," << specular[1] << "," << specular[2] << "," << specular[3] << 
	    " exponent " << expo << " cutoff " << cutoff << '\n';
	}
      }
      // material colors
      glGetMaterialfv(GL_FRONT,GL_AMBIENT,ambient);
      glGetMaterialfv(GL_FRONT,GL_DIFFUSE,diffuse);
      glGetMaterialfv(GL_FRONT,GL_SPECULAR,specular);
      glGetMaterialfv(GL_FRONT,GL_EMISSION,emission);
      glGetMaterialfv(GL_FRONT,GL_SHININESS,shini);
      cerr << "front " << ": " <<
	" ambient " << ambient[0] << "," << ambient[1] << "," << ambient[2] << "," << ambient[3] << 
	" diffuse " << diffuse[0] << "," << diffuse[1] << "," << diffuse[2] << "," << diffuse[3] << 
	" specular " << specular[0] << "," << specular[1] << "," << specular[2] << "," << specular[3] << 
	" emission " << emission[0] << "," << emission[1] << "," << emission[2] << "," << emission[3] << 
	" shininess " << shini[0] << '\n';
      glGetMaterialfv(GL_BACK,GL_AMBIENT,ambient);
      glGetMaterialfv(GL_BACK,GL_DIFFUSE,diffuse);
      glGetMaterialfv(GL_BACK,GL_SPECULAR,specular);
      glGetMaterialfv(GL_BACK,GL_EMISSION,emission);
      glGetMaterialfv(GL_BACK,GL_SHININESS,shini);
      cerr << "back " << ": " <<
	" ambient " << ambient[0] << "," << ambient[1] << "," << ambient[2] << "," << ambient[3] << 
	" diffuse " << diffuse[0] << "," << diffuse[1] << "," << diffuse[2] << "," << diffuse[3] << 
	" specular " << specular[0] << "," << specular[1] << "," << specular[2] << "," << specular[3] << 
	" emission " << emission[0] << "," << emission[1] << "," << emission[2] << "," << emission[3] << 
	" shininess " << shini[0] << '\n';	
    }
    if (est_hyperplan){
      vecteur P,n;
      if (!hyperplan_normal_point(v0,n,P))
	return;
      P=*evalf_double(P,1,contextptr)._VECTptr;
      vecteur Porig(P);
      double n1=evalf_double(n[0],1,contextptr)._DOUBLE_val;
      double n2=evalf_double(n[1],1,contextptr)._DOUBLE_val;
      double n3=evalf_double(n[2],1,contextptr)._DOUBLE_val;
      if (fill_polygon){
	vecteur v1,v2;
	if (!normal3d(n,v1,v2))
	  return;
	double v11=evalf_double(v1[0],1,contextptr)._DOUBLE_val;
	double v12=evalf_double(v1[1],1,contextptr)._DOUBLE_val;
	double v13=evalf_double(v1[2],1,contextptr)._DOUBLE_val;
	double v21=evalf_double(v2[0],1,contextptr)._DOUBLE_val;
	double v22=evalf_double(v2[1],1,contextptr)._DOUBLE_val;
	double v23=evalf_double(v2[2],1,contextptr)._DOUBLE_val;
	// moves P along v1 so that one coordinate is in clip
	if (std::abs(v11)>=std::abs(v12) && std::abs(v11)>=std::abs(v13)){
	  P=addvecteur(P,multvecteur((window_xmin-P[0])/v11,v1));
	}
	if (std::abs(v12)>std::abs(v11) && std::abs(v12)>=std::abs(v13)){
	  P=addvecteur(P,multvecteur((window_ymin-P[1])/v12,v1));
	}
	if (std::abs(v13)>std::abs(v11) && std::abs(v13)>std::abs(v12)){
	  P=addvecteur(P,multvecteur((window_zmin-P[2])/v13,v1));
	}
	// find a large constant so that the points are outside clipping
	double v31=v11+v21,v32=v12+v22,v33=v13+v23,v41=v11-v21,v42=v12-v22,v43=v13-v23;
	double nv3=std::sqrt(v31*v31+v32*v32+v33*v33);
	double nv4=std::sqrt(v41*v41+v42*v42+v43*v43);
	double lambda=1;
	if (std::abs(v31)>0.1*nv3)
	  lambda=giac_max(lambda,(window_xmax-window_xmin)/std::abs(v31));
	if (std::abs(v32)>0.1*nv3)
	  lambda=giac_max(lambda,(window_ymax-window_ymin)/std::abs(v32));
	if (std::abs(v33)>0.1*nv3)
	  lambda=giac_max(lambda,(window_zmax-window_zmin)/std::abs(v33));
	if (std::abs(v41)>0.1*nv4)
	  lambda=giac_max(lambda,(window_xmax-window_xmin)/std::abs(v41));
	if (std::abs(v42)>0.1*nv4)
	  lambda=giac_max(lambda,(window_ymax-window_ymin)/std::abs(v42));
	if (std::abs(v43)>0.1*nv4)
	  lambda=giac_max(lambda,(window_zmax-window_zmin)/std::abs(v43));
	lambda *= 2;
	if (display_mode & 0x8){
	  // if lighting enabled, we must draw small quads otherwise
	  // light will be incorrectly rendered (vertices too far)
	  // divide v1 and v2 by 10 (100 quads)
	  P=subvecteur(P,multvecteur(lambda,addvecteur(v1,v2))); // base point
	  unsigned hyperplan_light_rep=10;
	  v1=multvecteur(2*double(lambda)/hyperplan_light_rep,v1);
	  v2=multvecteur(2*double(lambda)/hyperplan_light_rep,v2);
	  for (unsigned j=1;j<hyperplan_light_rep;++j){
	    vecteur P1=addvecteur(P,multvecteur(int(j-1),v2));
	    for (unsigned i=1;i<hyperplan_light_rep;++i){
	      vecteur P2(addvecteur(P1,v1));
	      vecteur P3(addvecteur(P2,v2));
	      vecteur P4(addvecteur(P1,v2));
	      glBegin(GL_QUADS);
	      glNormal3d(n1,n2,n3);
	      glvertex(P1,0,0,contextptr);
	      glvertex(P2,0,0,contextptr);
	      glvertex(P3,0,0,contextptr);
	      glvertex(P4,0,0,contextptr);
	      glNormal3d(-n1,-n2,-n3);
	      glvertex(P1,0,0,contextptr);
	      glvertex(P4,0,0,contextptr);
	      glvertex(P3,0,0,contextptr);
	      glvertex(P2,0,0,contextptr);
	      glEnd();
	      // cerr << P1 << "," << P2 << "," << P3 << "," << P4 << '\n';
	      P1=P2;
	    }
	  }
	} // end if (display_mode & 0x8)
	else {
	  // Plan representation: face P+lambda*(-v1-v2), P+lambda*(-v1+v2), ...
	  vecteur P1(subvecteur(P,multvecteur(lambda,addvecteur(v1,v2))));
	  vecteur P2(subvecteur(P,multvecteur(lambda,addvecteur(v1,-v2))));
	  vecteur P4(addvecteur(P,multvecteur(lambda,addvecteur(v1,-v2))));
	  vecteur P3(addvecteur(P,multvecteur(lambda,addvecteur(v1,v2))));
	  glBegin(GL_QUADS);
	  // normal not needed (light not enabled)
	  glNormal3d(n1,n2,n3);
	  glvertex(P1,0,0,contextptr);
	  glvertex(P2,0,0,contextptr);
	  glvertex(P3,0,0,contextptr);
	  glvertex(P4,0,0,contextptr);
	  glNormal3d(-n1,-n2,-n3);
	  glvertex(P1,0,0,contextptr);
	  glvertex(P4,0,0,contextptr);
	  glvertex(P3,0,0,contextptr);
	  glvertex(P2,0,0,contextptr);
	  glEnd();
	}
      } // end fill_polygon
      else { // use equation
	glNormal3d(n1,n2,n3);
	if (std::abs(n1)>=std::abs(n2) && std::abs(n1)>=std::abs(n3)){
	  // x=a*y+b*z+c
	  double a=-n2/n1, b=-n3/n1;
	  double c=evalf_double(P[0]-a*P[1]-b*P[2],1,contextptr)._DOUBLE_val;
	  double dy=(window_ymax-window_ymin)/10;
	  for (int j=0;j<=10;++j){
	    double y=window_ymin+j*dy;
	    glBegin(GL_LINES);
	    glVertex3d(a*y+b*window_zmin+c,y,window_zmin);
	    glVertex3d(a*y+b*window_zmax+c,y,window_zmax);
	    glEnd();
	  }
	  double dz=(window_zmax-window_zmin)/10;
	  for (int j=0;j<=10;++j){
	    double z=window_zmin+j*dz;
	    glBegin(GL_LINES);
	    glVertex3d(a*window_ymin+b*z+c,window_ymin,z);
	    glVertex3d(a*window_ymax+b*z+c,window_ymax,z);
	    glEnd();
	  }
	}
	if (std::abs(n2)>std::abs(n1) && std::abs(n2)>=std::abs(n3)){
	  // y=a*x+b*z+c
	  double a=-n1/n2, b=-n3/n2;
	  double c=evalf_double(P[1]-a*P[0]-b*P[2],1,contextptr)._DOUBLE_val;
	  double dx=(window_xmax-window_xmin)/10;
	  for (int j=0;j<=10;++j){
	    double x=window_xmin+j*dx;
	    glBegin(GL_LINES);
	    glVertex3d(x,a*x+b*window_zmin+c,window_zmin);
	    glVertex3d(x,a*x+b*window_zmax+c,window_zmax);
	    glEnd();
	  }
	  double dz=(window_zmax-window_zmin)/10;
	  for (int j=0;j<=10;++j){
	    double z=window_zmin+j*dz;
	    glBegin(GL_LINES);
	    glVertex3d(window_xmin,a*window_xmin+b*z+c,z);
	    glVertex3d(window_xmax,a*window_xmax+b*z+c,z);
	    glEnd();
	  }
	}
	if (std::abs(n3)>std::abs(n1) && std::abs(n3)>std::abs(n2)){
	  // z=a*x+b*y+c
	  double a=-n1/n3, b=-n2/n3;
	  double c=evalf_double(P[2]-a*P[0]-b*P[1],1,contextptr)._DOUBLE_val;
	  double dx=(window_xmax-window_xmin)/10;
	  for (int j=0;j<=10;++j){
	    double x=window_xmin+j*dx;
	    glBegin(GL_LINES);
	    glVertex3d(x,window_ymin,a*x+b*window_ymin+c);
	    glVertex3d(x,window_ymax,a*x+b*window_ymax+c);
	    glEnd();
	  }
	  double dy=(window_ymax-window_ymin)/10;
	  for (int j=0;j<=10;++j){
	    double y=window_ymin+j*dy;
	    glBegin(GL_LINES);
	    glVertex3d(window_xmin,y,a*window_xmin+b*y+c);
	    glVertex3d(window_xmax,y,a*window_xmax+b*y+c);
	    glEnd();
	  }
	}
      }
      if (!hidden_name && show_names) 
	legende_draw(Porig,legende,labelpos);
      return;
    }
    if (v0.is_symb_of_sommet(at_hypersphere)){
      gen & f=v0._SYMBptr->feuille;
      if (f.type==_VECT && f._VECTptr->size()>=2){
	vecteur & v=*f._VECTptr;
	// Check that center is a 3-d point
	if (v.front().type==_VECT && v.front()._VECTptr->size()==3){
	  // Radius
	  gen r=v[1];
	  if (r.type==_VECT && r._VECTptr->size()==3)
	    r=l2norm(*r._VECTptr,contextptr); 
	  // Direction of axis, parallels and meridiens
	  vecteur dir1(makevecteur(1,0,0)),dir2(makevecteur(0,1,0)),dir3(makevecteur(0,0,1));
	  if (v.size()>=3 && v[2].type==_VECT && v[2]._VECTptr->size()==3){
	    dir3=*v[2]._VECTptr;
	    dir3=divvecteur(dir3,sqrt(dotvecteur(dir3,dir3),contextptr));
	    if (v.size()>=4 && v[3].type==_VECT && v[3]._VECTptr->size()==3){
	      dir1=*v[3]._VECTptr;
	      dir1=divvecteur(dir1,sqrt(dotvecteur(dir1,dir1),contextptr));
	    }
	    else {
	      if (!is_zero(dir3[0]) || !is_zero(dir3[1]) ){
		dir1=makevecteur(-dir3[1],dir3[0],0);
		dir1=divvecteur(dir1,sqrt(dotvecteur(dir1,dir1),contextptr));
	      }
	    }
	    dir2=cross(dir3,dir1,contextptr);
	  }
	  // optional discretisation info for drawing
	  if (v.size()>=5 && v[4].type==_INT_){
	    nphi=max(absint(v[4].val),3);
	    ntheta=ntheta;
	  }
	  if (v.size()>=6 && v[5].type==_INT_){
	    ntheta=max(absint(v[5].val),3);
	  }
	  // Now make the sphere
	  if (fill_polygon){
	    glsphere(*v.front()._VECTptr,r,dir1,dir2,dir3,ntheta,nphi,GL_QUADS,texture,contextptr);
	  }
	  if (!hidden_line){
	    if (fill_polygon)
	      glColor3d(0,0,0);
	    glsphere(*v.front()._VECTptr,r,dir1,dir2,dir3,min(ntheta,36),min(nphi,36),GL_LINE_LOOP,texture,contextptr);
	    xcas_color(couleur,true);
	  }
	  if (!hidden_name && show_names) legende_draw(v.front(),legende,labelpos);
	}
      }
      return;
    }
    if (v0.is_symb_of_sommet(at_curve) && v0._SYMBptr->feuille.type==_VECT && !v0._SYMBptr->feuille._VECTptr->empty()){
      gen & f = v0._SYMBptr->feuille._VECTptr->front();
      // f = vect[ pnt,var,xmin,xmax ]
      if (f.type==_VECT && f._VECTptr->size()>=4){
	vecteur & vf = *f._VECTptr;
	if (vf.size()>4){
	  gen poly=vf[4];
	  if (ckmatrix(poly)){
	    const_iterateur it=poly._VECTptr->begin(),itend=poly._VECTptr->end();
	    bool closed=fill_polygon && !(display_mode & 0x8); // ?is_approx_zero(poly._VECTptr->back()-poly._VECTptr->front(),window_xmin,window_xmax,window_ymin,window_ymax,window_zmin,window_zmax):false;
	    if (it->_VECTptr->size()==3){
	      glBegin(closed?GL_POLYGON:GL_LINE_STRIP);
	      // if (closed) ++it;
	      for (;it!=itend;++it){
		if (!glvertex(*it->_VECTptr,0,0,contextptr)){
		  glEnd();
		  glBegin(closed?GL_POLYGON:GL_LINE_STRIP);
		}
	      }
	      glEnd();
	      if (!hidden_name && show_names && !poly._VECTptr->empty()) legende_draw(poly._VECTptr->front(),legende,labelpos);
	      return ;
	    }
	  }
	}
	gen point=vf[0];
	gen var=vf[1];
	gen mini=vf[2];
	gen maxi=vf[3];
	bool closed=false;
	if (fill_polygon && !(display_mode & 0x8) )
	  closed=is_approx_zero(subst(point,var,mini,false,contextptr)-subst(point,var,maxi,false,contextptr),window_xmin,window_xmax,window_ymin,window_ymax,window_zmin,window_zmax);
	int n=nphi*10;
	gen delta=(maxi-mini)/n;
	glBegin(closed?GL_POLYGON:GL_LINE_STRIP);
	for (int i=closed?1:0;i<=n;++i){
	  if (!glvertex(*subst(point,var,mini,false,contextptr)._VECTptr,0,0,contextptr)){
	    glEnd();
	    glBegin(closed?GL_POLYGON:GL_LINE_STRIP);
	  }
	  mini = mini + delta;
	}
	glEnd();
	if (!hidden_name && show_names) legende_draw(mini,legende,labelpos);
      }
      return;
    }
    if (v0.is_symb_of_sommet(at_hypersurface)){
      gen & f = v0._SYMBptr->feuille;
      if (f.type!=_VECT || f._VECTptr->size()<3)
	return;
      gen & tmp = f._VECTptr->front();
      if (tmp.type!=_VECT || tmp._VECTptr->size()<4)
	return;
      if (tmp._VECTptr->size()>4){
	if (!fill_polygon)
	  glsurface((*tmp._VECTptr)[4],GL_LINE_LOOP,texture,contextptr);
	else {
	  glsurface((*tmp._VECTptr)[4],GL_QUADS,texture,contextptr);
	  // glLineWidth(width+2);
	  if (!hidden_line){
	    glColor3d(0,0,0);
	    glsurface((*tmp._VECTptr)[4],GL_LINE_LOOP,texture,contextptr);
	    xcas_color(couleur,true);
	  }
	}
	// if (!hidden_name && show_names) legende_draw(legende.c_str());
	return;
      }
      gen point = tmp._VECTptr->front(); // [x(u,v),y(u,v),z(u,v)]
      gen vars = (*tmp._VECTptr)[1]; // [u,v]
      gen mini = (*tmp._VECTptr)[2]; // [umin,vmin]
      gen maxi = (*tmp._VECTptr)[3]; // [umax,vmax]
      if (!check3dpoint(point) || 
	  vars.type!=_VECT || vars._VECTptr->size()!=2 
	  || mini.type!=_VECT || mini._VECTptr->size()!=2 
	  || maxi.type!=_VECT || maxi._VECTptr->size()!=2 )
	return;
      double umin=evalf_double(mini._VECTptr->front(),1,contextptr)._DOUBLE_val;
      double vmin=evalf_double(mini._VECTptr->back(),1,contextptr)._DOUBLE_val;
      double umax=evalf_double(maxi._VECTptr->front(),1,contextptr)._DOUBLE_val;
      double vmax=evalf_double(maxi._VECTptr->back(),1,contextptr)._DOUBLE_val;
      if (fill_polygon){
	glsurface(*point._VECTptr,vars,umin,umax,vmin,vmax,ntheta,nphi,GL_QUADS,contextptr);
	if (!hidden_line){
	  glColor3d(0,0,0);
	  glsurface(*point._VECTptr,vars,umin,umax,vmin,vmax,ntheta,nphi,GL_LINE_LOOP,contextptr);
	  xcas_color(couleur,true);
	}
      }
      else
	glsurface(*point._VECTptr,vars,umin,umax,vmin,vmax,ntheta,nphi,GL_LINE_LOOP,contextptr);
      // if (!hidden_name && show_names) legende_draw(?,legende,labelpos);
      return;
    }
    if (v0.type==_VECT && v0.subtype==_POINT__VECT && v0._VECTptr->size()==3 ){
      gen A(evalf_double(v0,1,contextptr));
      if (A.type==_VECT && A._VECTptr->size()==3 && type_point!=4){
	double xA=A._VECTptr->front()._DOUBLE_val;
	double yA=(*A._VECTptr)[1]._DOUBLE_val;
	double zA=A._VECTptr->back()._DOUBLE_val;
	double iA,jA,depthA;
	find_ij(xA,yA,zA,iA,jA,depthA);
	glLineWidth(1+epaisseur_point/2);
	switch(type_point){ 
	case 0:
	  glBegin(GL_LINES);
	  find_xyz(iA-epaisseur_point,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA-epaisseur_point,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  glEnd();
	  break;
	case 1:
	  // 1 losange, 
	  glBegin(GL_LINE_LOOP);
	  find_xyz(iA-epaisseur_point,jA,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  glEnd();
	  break;
	case 2:	  // 2 croix verticale, 
	  glBegin(GL_LINES);
	  find_xyz(iA,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA-epaisseur_point,jA,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  glEnd();
	  break;
	case 3: // 3 carre. 
	  glBegin(GL_LINE_LOOP);
	  find_xyz(iA-epaisseur_point,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA-epaisseur_point,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  glEnd();
	  break;
	case 5: // 5 triangle, 
	  glBegin(GL_LINE_LOOP);
	  find_xyz(iA+epaisseur_point,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA-epaisseur_point,jA,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  glEnd();
	  break;
	case 6:  // 6 etoile, 
	  glBegin(GL_LINES);
	  find_xyz(iA,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA-epaisseur_point/2,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point/2,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA-epaisseur_point/2,jA-epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  find_xyz(iA+epaisseur_point/2,jA+epaisseur_point,depthA,xA,yA,zA);
	  glVertex3d(xA,yA,zA);
	  glEnd();
	  break;
	default: 	  // 7 point
	  glBegin(GL_POINTS);
	  glvertex(*v0._VECTptr,0,0,contextptr);
	  glEnd();
	}
	glLineWidth(width+1);
      }
      if (!hidden_name && show_names) 
	legende_draw(v0,legende,labelpos);
      return;
    }
    if (v0.type==_VECT && v0.subtype!=_POINT__VECT){
      vecteur & vv0=*v0._VECTptr;
      vecteur w,lastpnt;
      if (v0.subtype==_POLYEDRE__VECT){
	// each element of v is a face 
	const_iterateur it=vv0.begin(),itend=vv0.end();
	for (;it!=itend;++it){
	  if (it->type==_VECT){ 
	    w = *evalf_double(*it,1,contextptr)._VECTptr;
	    int s=w.size();
	    if (s<3)
	      continue;
	    double d1,d2,d3,e1,e2,e3,f1,f2,f3;	    
	    if (get_glvertex(w[0],d1,d2,d3,0,0,contextptr) &&
		get_glvertex(w[1],e1,e2,e3,0,0,contextptr) &&
		get_glvertex(w[2],f1,f2,f3,0,0,contextptr) ){
	      glnormal(d1,d2,d3,e1,e2,e3,f1,f2,f3);
	      if (fill_polygon){
		glBegin(GL_POLYGON);
		glVertex3d(d1,d2,d3);
		glVertex3d(e1,e2,e3);
		glVertex3d(f1,f2,f3);
		for (int j=3;j<s;++j){
		  if (w[j].type!=_VECT || w[j]._VECTptr->size()!=3)
		    glvertex(vecteur(3,0),0,0,contextptr);
		  else 
		    glvertex( lastpnt=*w[j]._VECTptr,0,0,contextptr );
		}
		glEnd();
		xcas_color(couleur,true);
		// FIXME?? back face
		/*
		glnormal(f1,f2,f3,e1,e2,e3,d1,d2,d3);
		glBegin(GL_POLYGON);
		for (int j=s-1;j>=3;--j){
		  if (w[j].type!=_VECT || w[j]._VECTptr->size()!=3)
		    glvertex(vecteur(3,0),0,0,contextptr);
		  else 
		    glvertex( *w[j]._VECTptr ,0,0,contextptr);
		}
		glVertex3d(f1,f2,f3);
		glVertex3d(e1,e2,e3);
		glVertex3d(d1,d2,d3);
		glEnd();
		*/
	      }
	    }
	  }
	}
	for (it=vv0.begin();it!=itend;++it){
	  if (it->type==_VECT){ 
	    w = *evalf_double(*it,1,contextptr)._VECTptr;
	    int s=w.size();
	    if (s<3)
	      continue;
	    double d1,d2,d3,e1,e2,e3,f1,f2,f3;	    
	    if (get_glvertex(w[0],d1,d2,d3,0,0,contextptr) &&
		get_glvertex(w[1],e1,e2,e3,0,0,contextptr) &&
		get_glvertex(w[2],f1,f2,f3,0,0,contextptr) ){
	      glnormal(d1,d2,d3,e1,e2,e3,f1,f2,f3);
	      if (!hidden_line){
		glColor3d(0,0,0);
		glBegin(GL_LINE_LOOP);
		glVertex3d(d1,d2,d3);
		glVertex3d(e1,e2,e3);
		glVertex3d(f1,f2,f3);
		for (int j=3;j<s;++j){
		  if (w[j].type!=_VECT || w[j]._VECTptr->size()!=3)
		    glvertex(vecteur(3,0),0,0,contextptr);
		  else 
		    glvertex( lastpnt=*w[j]._VECTptr,0,0,contextptr );
		}
		glEnd();
		xcas_color(couleur,true);
	      }
	    }
	  } // end it->type==_VECT
	} // end for (;it!=itend;)
	if (!hidden_name && show_names) legende_draw(lastpnt,legende,labelpos);
	return;
      } // end polyedre
      int s=vv0.size();
      if (s==2 && check3dpoint(vv0.front()) && check3dpoint(vv0.back()) ){
	// segment, half-line, vector or line
	vecteur A(*evalf_double(vv0.front(),1,contextptr)._VECTptr),B(*evalf_double(vv0.back(),1,contextptr)._VECTptr);
	vecteur dir(subvecteur(B,A));
	double lambda=1,nu;
	double d1=evalf_double(dir[0],1,contextptr)._DOUBLE_val;
	double d2=evalf_double(dir[1],1,contextptr)._DOUBLE_val;
	double d3=evalf_double(dir[2],1,contextptr)._DOUBLE_val;
	double nd=std::sqrt(d1*d1+d2*d2+d3*d3);
	if (std::abs(d1)>=std::abs(d2) && std::abs(d1)>=std::abs(d3)){
	  nu=(window_xmin-evalf_double(A[0],1,contextptr)._DOUBLE_val)/d1;
	}
	if (std::abs(d2)>std::abs(d1) && std::abs(d2)>=std::abs(d3)){
	  nu=(window_ymin-evalf_double(A[1],1,contextptr)._DOUBLE_val)/d2;
	}
	if (std::abs(d3)>std::abs(d1) && std::abs(d3)>std::abs(d2)){
	  nu=(window_zmin-evalf_double(A[2],1,contextptr)._DOUBLE_val)/d3;
	}
	if (std::abs(d1)>=0.1*nd)
	  lambda=giac_max(lambda,(window_xmax-window_xmin)/std::abs(d1));
	if (std::abs(d2)>=0.1*nd)
	  lambda=giac_max(lambda,(window_ymax-window_ymin)/std::abs(d2));
	if (std::abs(d3)>=0.1*nd)
	  lambda=giac_max(lambda,(window_ymax-window_ymin)/std::abs(d3));
	lambda *= 2;
	glBegin(GL_LINE_STRIP);
	if (v0.subtype==_LINE__VECT){
	  // move A in clip
	  A=addvecteur(A,multvecteur(nu,dir));
	  B=subvecteur(B,multvecteur(nu,dir));
	  glvertex(subvecteur(A,multvecteur(lambda,dir)),0,0,contextptr);
	  glvertex(addvecteur(A,multvecteur(lambda,dir)),0,0,contextptr);
	}
	else {
	  glvertex(A,0,0,contextptr);
	  if (v0.subtype==_HALFLINE__VECT)
	    glvertex(addvecteur(A,multvecteur(lambda,dir)),0,0,contextptr);
	  else
	    glvertex(B,0,0,contextptr);
	}
	glEnd();
	if (v0.subtype==_VECTOR__VECT){
	  double xB=evalf_double(B[0],1,contextptr)._DOUBLE_val;
	  double yB=evalf_double(B[1],1,contextptr)._DOUBLE_val;
	  double zB=evalf_double(B[2],1,contextptr)._DOUBLE_val;
	  double xA=evalf_double(A[0],1,contextptr)._DOUBLE_val;
	  double yA=evalf_double(A[1],1,contextptr)._DOUBLE_val;
	  double zA=evalf_double(A[2],1,contextptr)._DOUBLE_val;
	  /* 2-d code */
	  double iA,jA,depthA,iB,jB,depthB,di,dj,dij;
	  find_ij(xB,yB,zB,iB,jB,depthB);
	  find_ij(xA,yA,zA,iA,jA,depthA);
	  di=iA-iB; dj=jA-jB;
	  dij=std::sqrt(di*di+dj*dj);
	  if (dij){
	    dij /= min(5,int(dij/10))+width;
	    di/=dij;
	    dj/=dij;
	    double dip=-dj,djp=di;
	    di*=std::sqrt(3.0);
	    dj*=std::sqrt(3.0);
	    double iC=iB+di+dip,jC=jB+dj+djp;
	    double iD=iB+di-dip,jD=jB+dj-djp;
	    double xC,yC,zC,xD,yD,zD;
	    find_xyz(iC,jC,depthB,xC,yC,zC);
	    find_xyz(iD,jD,depthB,xD,yD,zD);
	    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	    glBegin(GL_POLYGON);
	    glVertex3d(xB,yB,zB);
	    glVertex3d(xC,yC,zC);
	    glVertex3d(xD,yD,zD);
	    glEnd();
	  }
	}
	if (!hidden_name && show_names && v0.subtype!=_GROUP__VECT) 
	  legende_draw(multvecteur(0.5,addvecteur(A,B)),legende,labelpos);
	return;
      }
      if (s>2){ // polygon
	bool closed=vv0.front()==vv0.back();
	double d1,d2,d3,e1,e2,e3,f1,f2,f3;	    
	if (get_glvertex(vv0[0],d1,d2,d3,0,0,contextptr) &&
	    get_glvertex(vv0[1],e1,e2,e3,0,0,contextptr) &&
	    get_glvertex(vv0[2],f1,f2,f3,0,0,contextptr) ){
	  if (test_enable_texture(texture) && s<=5 && closed){
	    glBegin(GL_QUADS);
	    glnormal(d1,d2,d3,e1,e2,e3,f1,f2,f3);
	    glTexCoord2f(0,0);
	    glVertex3d(d1,d2,d3);
	    glTexCoord2f(0,1);
	    glVertex3d(e1,e2,e3);
	    glTexCoord2f(1,1);
	    glVertex3d(f1,f2,f3);
	    if (s==5){
	      glTexCoord2f(1,0);
	      glvertex(*vv0[3]._VECTptr,0,0,contextptr);
	    }
	    else
	      glVertex3d(d1,d2,d3);
	    glEnd();
	    glDisable(GL_TEXTURE_2D);
	  }
	  else {
	    glBegin(closed?GL_POLYGON:GL_LINE_STRIP);
	    glnormal(d1,d2,d3,e1,e2,e3,f1,f2,f3);
	    const_iterateur it=vv0.begin(),itend=vv0.end();
	    if (closed)
	      ++it;
	    for (;it!=itend;++it){
	      if (check3dpoint(*it))
		glvertex(*it->_VECTptr,0,0,contextptr);
	    }
	    glEnd();
	  }
	  if (!hidden_name && show_names && !vv0.empty()) 
	    legende_draw(vv0.front(),legende,labelpos);
	}
	return;
      }
    }
  }

  void Graph3d::indraw(const vecteur & v){
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
      indraw(*it);
  }

  void mult4(double * colmat,double * vect,double * res){
    res[0]=colmat[0]*vect[0]+colmat[4]*vect[1]+colmat[8]*vect[2]+colmat[12]*vect[3];
    res[1]=colmat[1]*vect[0]+colmat[5]*vect[1]+colmat[9]*vect[2]+colmat[13]*vect[3];
    res[2]=colmat[2]*vect[0]+colmat[6]*vect[1]+colmat[10]*vect[2]+colmat[14]*vect[3];
    res[3]=colmat[3]*vect[0]+colmat[7]*vect[1]+colmat[11]*vect[2]+colmat[15]*vect[3];
  }

  void mult4(double * colmat,float * vect,double * res){
    res[0]=colmat[0]*vect[0]+colmat[4]*vect[1]+colmat[8]*vect[2]+colmat[12]*vect[3];
    res[1]=colmat[1]*vect[0]+colmat[5]*vect[1]+colmat[9]*vect[2]+colmat[13]*vect[3];
    res[2]=colmat[2]*vect[0]+colmat[6]*vect[1]+colmat[10]*vect[2]+colmat[14]*vect[3];
    res[3]=colmat[3]*vect[0]+colmat[7]*vect[1]+colmat[11]*vect[2]+colmat[15]*vect[3];
  }

  void mult4(double * c,double k,double * res){
    for (int i=0;i<16;i++)
      res[i]=k*c[i];
  }
  
  double det4(double * c){
    return c[0]*c[5]*c[10]*c[15]-c[0]*c[5]*c[14]*c[11]-c[0]*c[9]*c[6]*c[15]+c[0]*c[9]*c[14]*c[7]+c[0]*c[13]*c[6]*c[11]-c[0]*c[13]*c[10]*c[7]-c[4]*c[1]*c[10]*c[15]+c[4]*c[1]*c[14]*c[11]+c[4]*c[9]*c[2]*c[15]-c[4]*c[9]*c[14]*c[3]-c[4]*c[13]*c[2]*c[11]+c[4]*c[13]*c[10]*c[3]+c[8]*c[1]*c[6]*c[15]-c[8]*c[1]*c[14]*c[7]-c[8]*c[5]*c[2]*c[15]+c[8]*c[5]*c[14]*c[3]+c[8]*c[13]*c[2]*c[7]-c[8]*c[13]*c[6]*c[3]-c[12]*c[1]*c[6]*c[11]+c[12]*c[1]*c[10]*c[7]+c[12]*c[5]*c[2]*c[11]-c[12]*c[5]*c[10]*c[3]-c[12]*c[9]*c[2]*c[7]+c[12]*c[9]*c[6]*c[3];
  }

  void inv4(double * c,double * res){
    res[0]=c[5]*c[10]*c[15]-c[5]*c[14]*c[11]-c[10]*c[7]*c[13]-c[15]*c[9]*c[6]+c[14]*c[9]*c[7]+c[11]*c[6]*c[13];
    res[1]=-c[1]*c[10]*c[15]+c[1]*c[14]*c[11]+c[10]*c[3]*c[13]+c[15]*c[9]*c[2]-c[14]*c[9]*c[3]-c[11]*c[2]*c[13];
    res[2]=c[1]*c[6]*c[15]-c[1]*c[14]*c[7]-c[6]*c[3]*c[13]-c[15]*c[5]*c[2]+c[14]*c[5]*c[3]+c[7]*c[2]*c[13];
    res[3]=-c[1]*c[6]*c[11]+c[1]*c[10]*c[7]+c[6]*c[3]*c[9]+c[11]*c[5]*c[2]-c[10]*c[5]*c[3]-c[7]*c[2]*c[9];
    res[4]=-c[4]*c[10]*c[15]+c[4]*c[14]*c[11]+c[10]*c[7]*c[12]+c[15]*c[8]*c[6]-c[14]*c[8]*c[7]-c[11]*c[6]*c[12];
    res[5]=c[0]*c[10]*c[15]-c[0]*c[14]*c[11]-c[10]*c[3]*c[12]-c[15]*c[8]*c[2]+c[14]*c[8]*c[3]+c[11]*c[2]*c[12];
    res[6]=-c[0]*c[6]*c[15]+c[0]*c[14]*c[7]+c[6]*c[3]*c[12]+c[15]*c[4]*c[2]-c[14]*c[4]*c[3]-c[7]*c[2]*c[12];
    res[7]=c[0]*c[6]*c[11]-c[0]*c[10]*c[7]-c[6]*c[3]*c[8]-c[11]*c[4]*c[2]+c[10]*c[4]*c[3]+c[7]*c[2]*c[8];
    res[8]=c[4]*c[9]*c[15]-c[4]*c[13]*c[11]-c[9]*c[7]*c[12]-c[15]*c[8]*c[5]+c[13]*c[8]*c[7]+c[11]*c[5]*c[12];
    res[9]=-c[0]*c[9]*c[15]+c[0]*c[13]*c[11]+c[9]*c[3]*c[12]+c[15]*c[8]*c[1]-c[13]*c[8]*c[3]-c[11]*c[1]*c[12];
    res[10]=c[0]*c[5]*c[15]-c[0]*c[13]*c[7]-c[5]*c[3]*c[12]-c[15]*c[4]*c[1]+c[13]*c[4]*c[3]+c[7]*c[1]*c[12];
    res[11]=-c[0]*c[5]*c[11]+c[0]*c[9]*c[7]+c[5]*c[3]*c[8]+c[11]*c[4]*c[1]-c[9]*c[4]*c[3]-c[7]*c[1]*c[8];
    res[12]=-c[4]*c[9]*c[14]+c[4]*c[13]*c[10]+c[9]*c[6]*c[12]+c[14]*c[8]*c[5]-c[13]*c[8]*c[6]-c[10]*c[5]*c[12];
    res[13]=c[0]*c[9]*c[14]-c[0]*c[13]*c[10]-c[9]*c[2]*c[12]-c[14]*c[8]*c[1]+c[13]*c[8]*c[2]+c[10]*c[1]*c[12];
    res[14]=-c[0]*c[5]*c[14]+c[0]*c[13]*c[6]+c[5]*c[2]*c[12]+c[14]*c[4]*c[1]-c[13]*c[4]*c[2]-c[6]*c[1]*c[12];
    res[15]=c[0]*c[5]*c[10]-c[0]*c[9]*c[6]-c[5]*c[2]*c[8]-c[10]*c[4]*c[1]+c[9]*c[4]*c[2]+c[6]*c[1]*c[8];
    double det=det4(c);
    mult4(res,1/det,res);
  }

  void dim32dim2(double * view,double * proj,double * model,double x0,double y0,double z0,double & i,double & j,double & dept){
    double vect[4]={x0,y0,z0,1},res1[4],res2[4];
    mult4(model,vect,res1);
    mult4(proj,res1,res2);
    i=res2[0]/res2[3]; // x and y are in [-1..1]
    j=res2[1]/res2[3];
    dept=res2[2]/res2[3];
    i=view[0]+(i+1)*view[2]/2;
    j=view[1]+(j+1)*view[3]/2;
    // x and y are the distance to the BOTTOM LEFT of the window
  }

  void dim22dim3(double * view,double * proj_inv,double * model_inv,double i,double j,double depth_,double & x,double & y,double & z){
    i=(i-view[0])*2/view[2]-1;
    j=(j-view[1])*2/view[3]-1;
    double res2[4]={i,j,depth_,1},res1[4],vect[4];
    mult4(proj_inv,res2,res1);
    mult4(model_inv,res1,vect);
    x=vect[0]/vect[3];
    y=vect[1]/vect[3];
    z=vect[2]/vect[3];
  }

  void Graph3d::find_ij(double x,double y,double z,double & i,double & j,double & depth_) {
    dim32dim2(view,proj,model,x,y,z,i,j,depth_);
#ifdef __APPLE__
    j=this->y()+h()-j;
    i=i+this->x();
#else
    j=window()->h()-j;
#endif
    // cout << i << " " << j <<  '\n';
  }

  double sqrt3over2=std::sqrt(double(3.0))/2;

  bool find_xmin_dx(double x0,double x1,double & xmin,double & dx){
    if (x0>=x1)
      return false;
    double x0x1=x1-x0;
    dx=std::pow(10,std::floor(std::log10(x0x1)));
    if (x0x1/dx>6)
      dx *= 2;    
    if (x0x1/dx<1.3)
      dx /= 5;
    if (x0x1/dx<3)
      dx /=2;
    if (!dx)
      return false;
    xmin=std::ceil(x0/dx)*dx;
    return true;
  }

  void Graph3d::draw_string(const string & s){
    if (printing){
      gl2psText(s.c_str(),"Helvetica",labelsize());
    }
    else {
#if !defined(__APPLE__) || !defined(INT128) // (does not work with textures on mac 64 bits)
      gl_font(FL_HELVETICA,labelsize());
      gl_draw(s.c_str());
#endif
    }
  }

  void normalize(double & a,double &b,double &c){
    double n=std::sqrt(a*a+b*b+c*c);
    a /= n;
    b /= n;
    c /= n;
  }

  void Graph3d::normal2plan(double & a,double &b,double &c){
    context * contextptr=hp?hp->contextptr:get_context(this);
    a /= std::pow(window_xmax-window_xmin,2);
    b /= std::pow(window_ymax-window_ymin,2);
    c /= std::pow(window_zmax-window_zmin,2);
    normalize(a,b,c);
    if (std::abs(a)<=1e-3){
      a=0;
      if (std::abs(b)<=1e-3){
	b=0;
	c=1;
      }
      else {
	gen gcb(float2rational(c/b,1e-3,contextptr));
	gen cbn,cbd;
	fxnd(gcb,cbn,cbd);
	if (cbn.type==_INT_ && cbd.type==_INT_ &&  std::abs(double(cbn.val)/cbd.val-c/a)<1e-2){
	  b=cbn.val;
	  c=cbd.val;
	}
      }
    }
    else {
      gen gba(float2rational(b/a,1e-3,contextptr));
      gen ban,bad,can,cad;
      fxnd(gba,ban,bad);
      if (ban.type==_INT_ && bad.type==_INT_  && std::abs(double(ban.val)/bad.val-b/a)<1e-2){
	gen gca(float2rational(c/a,1e-3,contextptr));
	fxnd(gca,can,cad);
	if (can.type==_INT_ && cad.type==_INT_ &&  std::abs(double(can.val)/cad.val-c/a)<1e-2){
	  int g=gcd(cad.val,bad.val);
	  // b/a=ban/bad=ban*(cad/g)/ppcm
	  int ai,bi,ci;
	  ai=(cad.val*bad.val/g);
	  bi=(ban.val*cad.val/g);
	  ci=(can.val*bad.val/g);
	  if (std::abs(ai)<=13 && std::abs(bi)<=13 && std::abs(ci)<13){
	    a=ai; b=bi; c=ci;
	  }
	}
      }
    }
  }

  void Graph3d::current_normal(double & a,double & b,double & c) {
    double res1[4]={0,0,1,1},vect[4];
    mult4(opengl?model_inv:invtransform3d,res1,vect);
    a=vect[0]/vect[3]-(window_xmax+window_xmin)/2;
    b=vect[1]/vect[3]-(window_ymax+window_ymin)/2;
    c=vect[2]/vect[3]-(window_zmax+window_zmin)/2;
    if (std::abs(a)<1e-3*(window_xmax-window_xmin))
      a=0;
    if (std::abs(b)<1e-3*(window_ymax-window_ymin))
      b=0;
    if (std::abs(c)<1e-3*(window_zmax-window_zmin))
      c=0;
  }

  void round0(double & x,double xmin,double xmax){
    if (std::abs(x)<1e-3*(xmax-xmin))
      x=0;
  }

  void Graph3d::display(){
    glEnable(GL_NORMALIZE);
    glEnable(GL_LINE_STIPPLE);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0,1.0);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
    glEnable(GL_CLIP_PLANE4);
    glEnable(GL_CLIP_PLANE5);
    // cout << glIsEnabled(GL_CLIP_PLANE0) << '\n';
    glDepthFunc(GL_LESS);
    glEnable(GL_DEPTH_TEST);
    glShadeModel((display_mode & 0x10)?GL_FLAT:GL_SMOOTH);
    bool lighting=display_mode & 0x8;
    if (lighting)
      glClearColor(0, 0, 0, 0);
    else
      glClearColor(1, 1, 1, 1);
    // clear the color and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);      
    // view transformations
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();   
    glLoadIdentity();
    bool notperspective=display_mode & 0x4;
    if (notperspective)
      // glFrustum(-sqrt3over2,sqrt3over2,-sqrt3over2,sqrt3over2,sqrt3over2,3*sqrt3over2);
      glOrtho(-sqrt3over2,sqrt3over2,-sqrt3over2,sqrt3over2,-sqrt3over2,sqrt3over2);
    else
      glFrustum(-0.5,0.5,-0.5,0.5,sqrt3over2,3*sqrt3over2);
    // put the visualisation cube inside above visualisation
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();   
    glLoadIdentity();
    double dx=(window_xmax-window_xmin),dy=(window_ymax-window_ymin),dz=(window_zmax-window_zmin);
    if (dx==0) { dx=1; ++window_xmax; }
    if (dy==0) { dy=1; ++window_ymax; }
    if (dz==0) { dz=1; ++window_zmax; }
    double x,y,z,theta;
    get_axis_angle_deg( (dragi || dragj)?(quaternion_double(dragi*180/h(),0,0)*rotation_2_quaternion_double(1,0,0,dragj*180/w())*q):q,x,y,z,theta);
    // cerr << theta << " " << x << "," << y << "," << z << '\n';
    if (!notperspective)
      glTranslated(0,0,-2*sqrt3over2);
    glRotated(theta,x,y,z);
    glScaled(1/dx,1/dy,1/dz);
    glTranslated(-(window_xmin+window_xmax)/2,-(window_ymin+window_ymax)/2,-(window_zmin+window_zmax)/2);
    // glRotated(theta_y,0,0,1);
    double plan0[]={1,0,0,0.5};
    double plan1[]={-1,0,0,0.5};
    double plan2[]={0,1,0,0.5};
    double plan3[]={0,-1,0,0.5};
    double plan4[]={0,0,1,0.5};
    double plan5[]={0,0,-1,0.5};
    plan0[3]=-window_xmin+dx/256;
    plan1[3]=window_xmax+dx/256;
    plan2[3]=-window_ymin+dy/256;
    plan3[3]=window_ymax+dy/256;
    plan4[3]=-window_zmin+dz/256;
    plan5[3]=window_zmax+dz/256;
    glGetDoublev(GL_PROJECTION_MATRIX,proj); // projection matrix in columns
    /*
    for (int i=0;i<15;++i)
      proj[i]=0;
    proj[0]=1/sqrt3over2;
    proj[5]=proj[0];
    proj[10]=-proj[0];
    proj[15]=1;
    */
    inv4(proj,proj_inv);
    glGetDoublev(GL_MODELVIEW_MATRIX,model); // modelview matrix in columns
    inv4(model,model_inv);
    glGetDoublev(GL_VIEWPORT,view);
    if (debug_infolevel>=2){
      double check[16];
      mult4(model,&model_inv[0],&check[0]);
      mult4(model,&model_inv[4],&check[4]);
      mult4(model,&model_inv[8],&check[8]);
      mult4(model,&model_inv[12],&check[12]);
      for (int i=0;i<16;++i){
	cout << model[i] << ",";
	if (i%4==3) cout << '\n';
      }
      cout << '\n';
      for (int i=0;i<16;++i){
	cout << model_inv[i] << ",";
	if (i%4==3) cout << '\n';
      }
      cout << '\n';
      for (int i=0;i<16;++i){
	cout << check[i] << ",";
	if (i%4==3) cout << '\n';
      }
      cout << '\n';
    }
    // drax
    bool fbox=(display_mode & 0x100);
    bool triedre=(display_mode & 0x200);
    glLineStipple(1,0xffff);
    glLineWidth(1);
    gl2psLineWidth(1);
    gl_color(FL_RED);
    // glColor3f(1,0,0);
    if (show_axes || triedre){
      glLineWidth(3);
      gl2psLineWidth(3);
      glBegin(GL_LINES);
      glVertex3d(0,0,0);
      glVertex3d(1,0,0);
      glEnd();
    }
    if (show_axes){
      glLineWidth(1);
      gl2psLineWidth(1);
      glBegin(GL_LINES);
      glVertex3d(0,0,0);
      glVertex3d(window_xmax,0,0);
      glEnd();
    }
    if (show_axes && !printing){
      glBegin(GL_LINES);
      glVertex3d(window_xmin,window_ymin,window_zmin);
      glVertex3d(window_xmax,window_ymin,window_zmin);
      glVertex3d(window_xmin,window_ymax,window_zmin);
      glVertex3d(window_xmax,window_ymax,window_zmin);
      glVertex3d(window_xmin,window_ymin,window_zmax);
      glVertex3d(window_xmax,window_ymin,window_zmax);
      glVertex3d(window_xmin,window_ymax,window_zmax);
      glVertex3d(window_xmax,window_ymax,window_zmax);
      glEnd();
    }
    gl_color(FL_GREEN);
    // glColor3f(0,1,0);
    if (show_axes || triedre){
      glLineWidth(3);
      gl2psLineWidth(3);
      glBegin(GL_LINES);
      glVertex3d(0,0,0);
      glVertex3d(0,1,0);
      glEnd();
    }
    if (show_axes){
      glLineWidth(1);
      gl2psLineWidth(1);
      glBegin(GL_LINES);
      glVertex3d(0,0,0);
      glVertex3d(0,window_ymax,0);
      glEnd();
    }
    if (show_axes && !printing){
      glBegin(GL_LINES);
      glVertex3d(window_xmin,window_ymin,window_zmin);
      glVertex3d(window_xmin,window_ymax,window_zmin);
      glVertex3d(window_xmax,window_ymin,window_zmin);
      glVertex3d(window_xmax,window_ymax,window_zmin);
      glVertex3d(window_xmin,window_ymin,window_zmax);
      glVertex3d(window_xmin,window_ymax,window_zmax);
      glVertex3d(window_xmax,window_ymin,window_zmax);
      glVertex3d(window_xmax,window_ymax,window_zmax);
      glEnd();
    }
    gl_color(FL_BLUE);
    // glColor3f(0,0,1);
    if (show_axes || triedre){
      glLineWidth(3);
      gl2psLineWidth(3);
      glBegin(GL_LINES);
      glVertex3d(0,0,0);
      glVertex3d(0,0,1);
      glEnd();
    }
    if (show_axes){
      glLineWidth(1);
      gl2psLineWidth(1);
      glBegin(GL_LINES);
      glVertex3d(0,0,0);
      glVertex3d(0,0,window_zmax);
      glEnd();
    }
    if (show_axes && !printing){
      glBegin(GL_LINES);
      glVertex3d(window_xmin,window_ymin,window_zmin);
      glVertex3d(window_xmin,window_ymin,window_zmax);
      glVertex3d(window_xmax,window_ymax,window_zmin);
      glVertex3d(window_xmax,window_ymax,window_zmax);
      glVertex3d(window_xmin,window_ymax,window_zmin);
      glVertex3d(window_xmin,window_ymax,window_zmax);
      glVertex3d(window_xmax,window_ymin,window_zmin);
      glVertex3d(window_xmax,window_ymin,window_zmax);
      glEnd();
    }
    if(show_axes){ // maillage
      glColor3f(1,0,0);
      glRasterPos3d(1,0,0);
      draw_string(x_axis_name.empty()?"x":x_axis_name);
      glColor3f(0,1,0);
      glRasterPos3d(0,1,0);
      draw_string(y_axis_name.empty()?"y":y_axis_name);
      glColor3f(0,0,1);
      glRasterPos3d(0,0,1);
      draw_string(z_axis_name.empty()?"z":z_axis_name);
      if (fbox){
	double xmin,dx,x,ymin,dy,y,zmin,dz,z;
	find_xmin_dx(window_xmin,window_xmax,xmin,dx);
	find_xmin_dx(window_ymin,window_ymax,ymin,dy);
	find_xmin_dx(window_zmin,window_zmax,zmin,dz);
	glLineStipple(1,0x3333);
	glBegin(GL_LINES);
	gl_color(FL_CYAN);
	for (z=zmin;z<=window_zmax;z+=2*dz){
	  for (y=ymin;y<=window_ymax;y+=2*dy){
	    glVertex3d(window_xmin,y,z);
	    glVertex3d(window_xmax,y,z);
	  }
	}
	for (z=zmin;z<=window_zmax;z+=2*dz){
	  for (x=xmin;x<=window_xmax;x+=2*dx){
	    glVertex3d(x,window_ymin,z);
	    glVertex3d(x,window_ymax,z);
	  }
	}
	for (x=xmin;x<=window_xmax;x+=2*dx){
	  for (y=ymin;y<=window_ymax;y+=2*dy){
	    glVertex3d(x,y,window_zmin);
	    glVertex3d(x,y,window_zmax);
	  }
	}
	glEnd();
	gl_color((display_mode & 0x8)?FL_WHITE:FL_BLACK);
	glPointSize(3);
	for (x=xmin;x<=window_xmax;x+=dx){
	  round0(x,window_xmin,window_xmax);
	  glBegin(GL_POINTS);
	  glVertex3d(x,window_ymin,window_zmin);
	  glEnd();
	  glRasterPos3d(x,window_ymin,window_zmin);
	  string tmps=giac::print_DOUBLE_(x,2)+x_axis_unit;
	  draw_string(tmps);
	}
	for (y=ymin;y<=window_ymax;y+=dy){
	  round0(y,window_ymin,window_ymax);
	  glBegin(GL_POINTS);
	  glVertex3d(window_xmin,y,window_zmin);
	  glEnd();
	  glRasterPos3d(window_xmin,y,window_zmin);
	  string tmps=giac::print_DOUBLE_(y,2)+y_axis_unit;
	  draw_string(tmps);
	}
	for (z=zmin;z<=window_zmax;z+=dz){
	  round0(z,window_zmin,window_zmax);
	  glBegin(GL_POINTS);
	  glVertex3d(window_xmin,window_ymin,z);
	  glEnd();
	  glRasterPos3d(window_xmin,window_ymin,z);
	  string tmps=giac::print_DOUBLE_(z,2)+z_axis_unit;
	  draw_string(tmps);
	}
      }
    }
    glClipPlane(GL_CLIP_PLANE0,plan0);
    glClipPlane(GL_CLIP_PLANE1,plan1);
    glClipPlane(GL_CLIP_PLANE2,plan2);
    glClipPlane(GL_CLIP_PLANE3,plan3);
    glClipPlane(GL_CLIP_PLANE4,plan4);
    glClipPlane(GL_CLIP_PLANE5,plan5);
    // mouse plan
    double normal_a,normal_b,normal_c;
    current_normal(normal_a,normal_b,normal_c);
    normal2plan(normal_a,normal_b,normal_c);
    double plan_x0,plan_y0,plan_z0,plan_t0;
    find_xyz(0,0,depth,plan_x0,plan_y0,plan_z0);
    plan_t0=normal_a*plan_x0+normal_b*plan_y0+normal_c*plan_z0;
    if (std::abs(plan_t0)<std::abs(window_zmax-window_zmin)/1000)
      plan_t0=0;
    if (show_axes && !printing && parent() && dynamic_cast<Figure *>(parent())){
      gl_color((display_mode & 0x8)?FL_WHITE:FL_BLACK);
      glLineWidth(2);
      gl2psLineWidth(2);
      double xa,ya,za,xb,yb,zb;
      // show mouse plan intersections with the faces of the clip planes
      // example with the face z=Z:
      // a*x+b*y=plan_t0-normal_c*Z
      // get the 4 intersections, keep only the valid ones
      glLineStipple(1,0x3333);
      if (is_clipped(normal_a,window_xmin,window_xmax,normal_b,window_ymin,window_ymax,plan_t0-normal_c*window_zmin,xa,ya,xb,yb)){
	glBegin(GL_LINE_STRIP);
	glVertex3d(xa,ya,window_zmin+dz/256);
	glVertex3d(xb,yb,window_zmin+dz/256);
	glEnd();
      }
      if (is_clipped(normal_a,window_xmin,window_xmax,normal_b,window_ymin,window_ymax,plan_t0-normal_c*window_zmax,xa,ya,xb,yb)){
	glBegin(GL_LINE_STRIP);
	glVertex3d(xa,ya,window_zmax-dz/256);
	glVertex3d(xb,yb,window_zmax-dz/256);
	glEnd();
      }
      if (is_clipped(normal_a,window_xmin,window_xmax,normal_c,window_zmin,window_zmax,plan_t0-normal_b*window_ymin,xa,za,xb,zb)){
	glBegin(GL_LINE_STRIP);
	glVertex3d(xa,window_ymin+dy/256,za);
	glVertex3d(xb,window_ymin+dy/256,zb);
	glEnd();
      }
      if (is_clipped(normal_a,window_xmin,window_xmax,normal_c,window_zmin,window_zmax,plan_t0-normal_b*window_ymax,xa,za,xb,zb)){
	glBegin(GL_LINE_STRIP);
	glVertex3d(xa,window_ymax-dy/256,za);
	glVertex3d(xb,window_ymax-dy/256,zb);
	glEnd();
      }
      if (is_clipped(normal_b,window_ymin,window_ymax,normal_c,window_zmin,window_zmax,plan_t0-normal_a*window_xmin,ya,za,yb,zb)){
	glBegin(GL_LINE_STRIP);
	glVertex3d(window_xmin+dx/256,ya,za);
	glVertex3d(window_xmin+dx/256,yb,zb);
	glEnd();
      }
      if (is_clipped(normal_b,window_ymin,window_ymax,normal_c,window_zmin,window_zmax,plan_t0-normal_a*window_xmax,ya,za,yb,zb)){
	glBegin(GL_LINE_STRIP);
	glVertex3d(window_xmax-dx/256,ya,za);
	glVertex3d(window_xmax-dx/256,yb,zb);
	glEnd();
      }
      // same for window_zmax, etc.
      glLineWidth(1);
      gl2psLineWidth(1);
      glLineStipple(1,0xffff);
    }
    if (lighting){
      /*
      GLfloat mat_specular[] = { 1.0,1.0,1.0,1.0 };
      GLfloat mat_shininess[] = { 50.0 };
      glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
      glMaterialfv(GL_FRONT,GL_SHININESS,mat_shininess);
      */
      glEnable(GL_LIGHTING);
      glEnable(GL_CULL_FACE);
      glEnable(GL_COLOR_MATERIAL);
      static GLfloat l_pos[8][4],l_dir[8][3],ambient[8][4],diffuse[8][4],specular[8][4];
      for (int i=0;i<8;++i){
	if (!light_on[i]){
	  glDisable(GL_LIGHT0+i);
	  continue;
	}
	glEnable(GL_LIGHT0+i);
	l_pos[i][0]=light_x[i];
	l_pos[i][1]=light_y[i];
	l_pos[i][2]=light_z[i];
	l_pos[i][3]=light_w[i];
	glLightfv(GL_LIGHT0+i,GL_POSITION,l_pos[i]);
	l_dir[i][0]=light_spot_x[i];
	l_dir[i][1]=light_spot_y[i];
	l_dir[i][2]=light_spot_z[i];
	glLightfv(GL_LIGHT0+i,GL_SPOT_DIRECTION,l_dir[i]);
	glLightf(GL_LIGHT0+i,GL_SPOT_EXPONENT,light_spot_exponent[i]);
	glLightf(GL_LIGHT0+i,GL_SPOT_CUTOFF,light_spot_cutoff[i]);
	glLightf(GL_LIGHT0+i,GL_CONSTANT_ATTENUATION,light_0[i]);
	glLightf(GL_LIGHT0+i,GL_LINEAR_ATTENUATION,light_1[i]);
	glLightf(GL_LIGHT0+i,GL_QUADRATIC_ATTENUATION,light_2[i]);
	ambient[i][0]=light_ambient_r[i];
	ambient[i][1]=light_ambient_g[i];
	ambient[i][2]=light_ambient_b[i];
	ambient[i][3]=light_ambient_a[i];
	glLightfv(GL_LIGHT0+i,GL_AMBIENT,ambient[i]);
	diffuse[i][0]=light_diffuse_r[i];
	diffuse[i][1]=light_diffuse_g[i];
	diffuse[i][2]=light_diffuse_b[i];
	diffuse[i][3]=light_diffuse_a[i];
	glLightfv(GL_LIGHT0+i,GL_DIFFUSE,diffuse[i]);
	specular[i][0]=light_specular_r[i];
	specular[i][1]=light_specular_g[i];
	specular[i][2]=light_specular_b[i];
	specular[i][3]=light_specular_a[i];
	glLightfv(GL_LIGHT0+i,GL_SPECULAR,specular[i]);
	if (debug_infolevel>=2){
	  GLfloat posf[4],direcf[4],ambient[4],diffuse[4],specular[4];
	  GLfloat expo,cutoff;
	  double pos[4],direc[4];
	  glGetLightfv(GL_LIGHT0+i,GL_SPOT_EXPONENT,&expo);
	  glGetLightfv(GL_LIGHT0+i,GL_SPOT_CUTOFF,&cutoff);
	  glGetLightfv(GL_LIGHT0+i,GL_POSITION,posf);
	  glGetLightfv(GL_LIGHT0+i,GL_SPOT_DIRECTION,direcf);
	  direcf[3]=0; // 4 was out of bound
	  mult4(model_inv,posf,pos);
	  tran4(model);
	  mult4(model,direcf,direc);
	  tran4(model);
	  glGetLightfv(GL_LIGHT0+i,GL_AMBIENT,ambient);
	  glGetLightfv(GL_LIGHT0+i,GL_DIFFUSE,diffuse);
	  glGetLightfv(GL_LIGHT0+i,GL_SPECULAR,specular);
	  cerr << "light " << i << ": " <<
	    " pos " << pos[0] << "," << pos[1] << "," << pos[2] << "," << pos[3] << 
	    " dir " << direc[0] << "," << direc[1] << "," << direc[2] << "," << direc[3] << 
	    " ambient " << ambient[0] << "," << ambient[1] << "," << ambient[2] << "," << ambient[3] << 
	    " diffuse " << diffuse[0] << "," << diffuse[1] << "," << diffuse[2] << "," << diffuse[3] << 
	    " specular " << specular[0] << "," << specular[1] << "," << specular[2] << "," << specular[3] << 
	    " exponent " << expo << " cutoff " << cutoff << '\n';
	}
      }
    }
    if (display_mode & 0x20){
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    }
    // now draw each object
    if ( (display_mode & 2) && !animation_instructions.empty())
      indraw(animation_instructions[animation_instructions_pos % animation_instructions.size()]);
    if ( display_mode & 0x40 )
      indraw(trace_instructions);
    if (display_mode & 1)
      indraw(plot_instructions);
    if (display_mode & 0x8){
      glDisable(GL_LIGHTING);
      for (int i=0;i<8;++i)
	glDisable(GL_LIGHT0+i);
    }
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    gen plot_tmp,title_tmp;
    History_Pack * hp =get_history_pack(this);
    context * contextptr=hp?hp->contextptr:get_context(this);
    find_title_plot(title_tmp,plot_tmp,contextptr);
    indraw(plot_tmp);
    if (mode==1 && pushed && push_in_area && in_area){
      // draw segment between push and current
      gl_color(FL_RED);
      double x,y,z;
      glBegin(GL_LINES);
      find_xyz(push_i,push_j,push_depth,x,y,z);
      glVertex3d(x,y,z);
      find_xyz(current_i,current_j,current_depth,x,y,z);
      glVertex3d(x,y,z);
      glEnd();
    }
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
    /*
    if( show_axes){
      gl_color((display_mode & 0x8)?FL_WHITE:FL_BLACK);
      glRasterPos3d(window_xmin,window_ymin,window_zmin);
      string tmps=giac::print_DOUBLE_(window_xmin,2)+","+giac::print_DOUBLE_(window_ymin,2)+","+giac::print_DOUBLE_(window_zmin,2);
      draw_string(tmps);
      glRasterPos3d(window_xmax,window_ymax,window_zmax);
      tmps=giac::print_DOUBLE_(window_xmax,2)+","+giac::print_DOUBLE_(window_ymax,2)+","+giac::print_DOUBLE_(window_zmax,2);
      draw_string(tmps);
    }
    */
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gl_color((display_mode & 0x8)?FL_WHITE:FL_BLACK);
    double td=fl_width(title.c_str());
    glRasterPos3d(-0.4*td/w(),-1,depth-0.001);
    string mytitle(title);
    if (!is_zero(title_tmp) && function_final.type==_FUNC)
      mytitle=gen(symbolic(*function_final._FUNCptr,title_tmp)).print(contextptr);
    if (!mytitle.empty())
      draw_string(mytitle);
    glRasterPos3d(-1,-1,depth-0.001);
    if (!args_help.empty() && args_tmp.size()<= args_help.size()){
      draw_string(gettext("Click ")+args_help[giacmax(1,args_tmp.size())-1]);
    }
    glRasterPos3d(-0.98,0.87,depth-0.001);
    if (show_axes && !printing){
      string tmps=gettext("mouse plan ")+giac::print_DOUBLE_(normal_a,3)+"x+"+giac::print_DOUBLE_(normal_b,3)+"y+"+giac::print_DOUBLE_(normal_c,3)+"z="+ giac::print_DOUBLE_(plan_t0,3);
      // cerr << tmps << '\n';
      draw_string(tmps); // +" Z="+giac::print_DOUBLE_(-depth,3));
    }
    if (displayhint){
      glRasterPos3d(-0.98,-0.99,depth-0.001);
      draw_string(gettext("click k: kill 3d view, o: native rendering"));
    }
    if (below_depth_hidden){
      /* current mouse position
	double i=Fl::event_x();
	double j=Fl::event_y();
	j=window()->h()-j;
	i=(i-view[0])*2/view[2]-1;
	j=(j-view[1])*2/view[3]-1;
      */
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      xcas_color((display_mode & 0x8)?FL_WHITE:FL_BLACK,true);
      glBegin(GL_POLYGON);
      glVertex3d(-1,-1,depth);
      glVertex3d(-1,1,depth);
      glVertex3d(1,1,depth);
      glVertex3d(1,-1,depth);
      glEnd();
    }
    glDisable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glFlush();
    glFinish();
    if (screenbuf && !printing){
      glReadBuffer(GL_FRONT);
      glReadPixels(this->x(), window()->h()-this->y()-this->h(), w(), h(), GL_RGBA, GL_UNSIGNED_BYTE, screenbuf);
    }
  }

  void Graph3d::print(){
    if (!printing)
      return;
    // EPS output
    glViewport(0,0, w(), h()); // lower left
    int buf=1024*1024;
    for (;;buf+=1024*1024){
      gl2psBeginPage(label(),"xcasnew",0,GL2PS_EPS,GL2PS_BSP_SORT,GL2PS_SIMPLE_LINE_OFFSET | GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT,GL_RGBA,0,0,0,0,0,buf,printing,"output");
      gl2psEnable(GL2PS_LINE_STIPPLE);
      display();
      if (gl2psEndPage()!=GL2PS_OVERFLOW){
	break;
      }
    }    
  }

  bool find_clip_box(Fl_Widget * wid,int & x,int & y,int & w,int & h){
    if (!wid || !wid->window())
      return false;
    fl_push_no_clip();
    vector<Fl_Widget *> p;
    for (;wid;wid=wid->parent())
      p.push_back(wid);
    for (int i=p.size()-2;i>=0;--i){
      fl_clip_box(p[i]->x(),p[i]->y(),p[i]->w(),p[i]->h(),x,y,w,h);
      fl_push_clip(x,y,w,h);
    }
    for (int i=p.size()-1;i>=0;--i)
      fl_pop_clip();
    return true;
  }

  /*  
      native 3d rendering
   */
  void do_transform(const double mat[16],double x,double y,double z,double & X,double & Y,double &Z){
    X=mat[0]*x+mat[1]*y+mat[2]*z+mat[3];
    Y=mat[4]*x+mat[5]*y+mat[6]*z+mat[7];
    Z=mat[8]*x+mat[9]*y+mat[10]*z+mat[11];
    // double t=mat[12]*x+mat[13]*y+mat[14]*z+mat[15];
    // X/=t; Y/=t; Z/=t;
  }
  void Graph3d::find_xyz(double i,double j,double depth_,double & x,double & y,double & z) {
    if (!opengl){ // FIXME
      update_rotation();
      // cout <<  i << " " << j << "  " <<  this->x() << " " <<  this->y() << " " << this->w() << " " << this->h() << " \n" ;
      i -= this->x();
      j -= this->y();
      int horiz=LCD_WIDTH_PX/2,vert=horiz/2;//LCD_HEIGHT_PX/2;
      double lcdz= LCD_HEIGHT_PX/4;
      double xmin=-1,ymin=-1,xmax=1,ymax=1,xscale=0.6*(xmax-xmin)/horiz,yscale=0.6*(ymax-ymin)/vert;
      double Z=current_depth; // -1..1
      double I=i-horiz;
      double J=j-LCD_HEIGHT_PX/2+lcdz*Z;
      double X=yscale*J-xscale*I;
      double Y=yscale*J+xscale*I;
      do_transform(invtransform3d,X,Y,Z,x,y,z);
      return;
    }
#ifdef __APPLE__
    j=this->y()+h()-j;
    i=i-this->x();
#else
    j=window()->h()-j;
#endif
    dim22dim3(view,proj_inv,model_inv,i,j,depth_,x,y,z);
  }

  void Graph3d::xyz2ij(const double3 &d,int &i,int &j) const {
    double X,Y,Z;
    do_transform(transform3d,d.x,d.y,d.z,X,Y,Z);
    i=LCD_WIDTH_PX/2+(Y-X)/4.8*LCD_WIDTH_PX;
    j=LCD_HEIGHT_PX/2-Z*lcdz+(Y+X)/9.6*LCD_WIDTH_PX;    
  }
  
  void Graph3d::xyz2ij(const double3 &d,double &i,double &j) const {
    double X,Y,Z;
    do_transform(transform3d,d.x,d.y,d.z,X,Y,Z);
    i=LCD_WIDTH_PX/2+(Y-X)/4.8*LCD_WIDTH_PX;
    j=LCD_HEIGHT_PX/2-Z*lcdz+(Y+X)/9.6*LCD_WIDTH_PX;    
  }
  
  void Graph3d::xyz2ij(const double3 &d,double &i,double &j,double3 & d3) const {
    do_transform(transform3d,d.x,d.y,d.z,d3.x,d3.y,d3.z);
    i=LCD_WIDTH_PX/2+(d3.y-d3.x)/4.8*LCD_WIDTH_PX;
    j=LCD_HEIGHT_PX/2-d3.z*lcdz+(d3.y+d3.x)/9.6*LCD_WIDTH_PX;    
  }
  
  void Graph3d::XYZ2ij(const double3 &d,int &i,int &j) const {
    double X=d.x,Y=d.y,Z=d.z;
    i=LCD_WIDTH_PX/2+(Y-X)/4.8*LCD_WIDTH_PX;
    j=LCD_HEIGHT_PX/2-Z*lcdz+(Y+X)/9.6*LCD_WIDTH_PX;    
  }

  // returns true if filled, false otherwise
  bool get_colors(gen attr,int & upcolor,int & downcolor,int & downupcolor,int & downdowncolor){
    if (attr.is_symb_of_sommet(at_pnt)){
      attr=attr[1];
    }
    if (attr.type==_INT_ && (attr.val & 0xffff)!=56){
      upcolor=attr.val &0xffff;
      int color=rgb565to888(upcolor);
      int r=(color&0xff0000)>>16,g=(color & 0xff00)>>8,b=192;
      r >>= 2;
      g >>= 2;
      downcolor=rgb888to565((r<<16)|(g<<8)|b);
      r >>= 1;
      g >>= 1;
      downupcolor=rgb888to565((r<<16)|(g<<8)|b);
      r >>= 2;
      g >>= 2;
      downdowncolor=rgb888to565((r<<16)|(g<<8)|b);
    }
    if (attr.type==_INT_)
      return attr.val & 0x40000000;
    return false;
  }

  void Graph3d::update_rotation(){
    context * contextptr=hp?hp->contextptr:get_context(this);
    solid3d=false;
    double rx,ry,rz,theta;
    get_axis_angle_deg( (dragi || dragj)?(rotation_2_quaternion_double(0,0,1,dragi*180/h())*rotation_2_quaternion_double(0,1,0,dragj*180/w())*q):q,rx,ry,rz,theta);
    // get_axis_angle_deg( (dragi || dragj)?(quaternion_double(dragi*180/h(),0,0)*rotation_2_quaternion_double(1,0,0,dragj*180/w())*q):q,rx,ry,rz,theta);
    // rx=-0.51; ry=-.197; rz=-.835; theta=327.88;
    // rx=0;ry=1;rz=0; //theta=60;
    // cout << rx << " " << ry << " " << rz << " " << theta << "\n";
    // rx=0;ry=0;rz=1;
    // theta=10;
    theta *= M_PI/180;
    double r2=rx*rx+ry*ry+rz*rz,invr2=1/r2;
    double r=std::sqrt(r2);
    double c=std::cos(theta),s=std::sin(theta);
    // mkisom([[rx,ry,rz],theta],1)
    // 1/r2*[[rx*rx+ry*ry*c+rz*rz*c,rx*ry-rx*ry*c-rz*s*r,rx*rz-rx*rz*c+ry*s*r],[rx*ry-rx*ry*c+rz*s*r,ry*ry+rx*rx*c+rz*rz*c,ry*rz-rx*s*r-ry*rz*c],[rx*rz-rx*rz*c-ry*s*r,ry*rz+rx*s*r-ry*rz*c,rz*rz+rx*rx*c+ry*ry*c]]
    double a11=invr2*(rx*rx+ry*ry*c+rz*rz*c),a12=invr2*(rx*ry-rx*ry*c-rz*s*r),a13=invr2*(rx*rz-rx*rz*c+ry*s*r);
    double a21=invr2*(rx*ry-rx*ry*c+rz*s*r),a22=invr2*(ry*ry+rx*rx*c+rz*rz*c),a23=invr2*(ry*rz-rx*s*r-ry*rz*c);
    double a31=invr2*(rx*rz-rx*rz*c-ry*s*r),a32=invr2*(ry*rz+rx*s*r-ry*rz*c),a33=invr2*(rz*rz+rx*rx*c+ry*ry*c);
    double xt=(window_xmin+window_xmax)/2,xs=2.0/(window_xmax-window_xmin),yt=(window_ymin+window_ymax)/2,ys=2.0/(window_ymax-window_ymin),zt=(window_zmin+window_zmax)/2,zs=2.0/(window_zmax-window_zmin);
    double mat[16]={a11*xs,a12*ys,a13*zs,-a11*xt*xs-a12*yt*ys-a13*zt*zs,
		    a21*xs,a22*ys,a23*zs,-a21*xt*xs-a22*yt*ys-a23*zt*zs,
		    a31*xs,a32*ys,a33*zs,-a31*xt*xs-a32*yt*ys-a33*zt*zs,
		    0,0,0,1
    };
    for (int i=0;i<sizeof(mat)/sizeof(double);++i)
      transform3d[i]=mat[i];
    inv4(transform3d,invtransform3d);
    surfacev.clear();
    polyedrev.clear(); polyedre_xyminmax.clear(); polyedre_abcv.clear(); polyedre_faceisclipped.clear(); polyedre_filled.clear();
    plan_pointv.clear(); plan_abcv.clear(); plan_filled.clear();
    sphere_centerv.clear(); sphere_radiusv.clear(); sphere_quadraticv.clear();
    linev.clear(); linetypev.clear(); curvev.clear();
    pointv.clear(); points.clear();
    plan_color.clear();sphere_color.clear();polyedre_color.clear();polyedre_faceisclipped.clear();line_color.clear();curve_color.clear(); hyp_color.clear(); point_color.clear();
    // rotate+translate+scale g
    vecteur v;
    aplatir(native_renderer_instructions,v);
    for (int i=0;i<v.size();++i){
      int u=default_upcolor,d=default_downcolor,du=default_downupcolor,dd=default_downdowncolor;
      const char * ptr=0;
      bool fill_polyedre=false;
      if (v[i].is_symb_of_sommet(at_pnt)){
	vecteur & attrv=*v[i]._SYMBptr->feuille._VECTptr;
	if (attrv.size()>1){
	  gen attr=attrv[1];
	  fill_polyedre=get_colors(attr,u,d,du,dd);
	  if (fill_polyedre) solid3d=true;
	  if (attrv.size()>2){
	    attr=attrv[2];
	    if (attr.type==_STRNG)
	      ptr=attr._STRNGptr->c_str();
	    if (attr.type==_IDNT)
	      ptr=attr._IDNTptr->id_name;
	  }
	}
      }
      gen G=remove_at_pnt(v[i]);
      if (G.is_symb_of_sommet(at_curve)){
	gen f=G._SYMBptr->feuille;
	f=f[0];
	f=f[4];
	if (ckmatrix(f)){
	  vecteur v=*f._VECTptr;
	  int n=v.size();
	  curvev.push_back(vector<double3>(n));
	  curve_color.push_back(int4(u,d,du,dd));
	  vector<double3> & cur=curvev.back();
	  for (int k=0;k<n;++k){
	    gen P=v[k];
	    if (P.type==_VECT && P._VECTptr->size()==3){
	      double X,Y,Z;
	      do_transform(transform3d,P[0]._DOUBLE_val,P[1]._DOUBLE_val,P[2]._DOUBLE_val,X,Y,Z);
	      cur[k]=double3(X,Y,Z);
	    }
	  }
	}
	continue;
      }
      if (G.type==_VECT && G.subtype==_POINT__VECT){
	gen m=evalf_double(G,1,contextptr);
	if (m.type==_VECT && m._VECTptr->size()==3){
	  vecteur & A=*m._VECTptr;
	  if (A[0].type==_DOUBLE_ && A[1].type==_DOUBLE_ && A[2].type==_DOUBLE_ ){
	    double X,Y,Z; int I,J;
	    do_transform(transform3d,A[0]._DOUBLE_val,A[1]._DOUBLE_val,A[2]._DOUBLE_val,X,Y,Z);
	    xyz2ij(double3(A[0]._DOUBLE_val,A[1]._DOUBLE_val,A[2]._DOUBLE_val),I,J);
	    pointv.push_back(double3(I,J,Z));
	    points.push_back(ptr);
	    point_color.push_back(int4(u,d,du,dd));
	  }
	}
	continue;
      }
      bool line=G.subtype==_LINE__VECT,halfline=G.subtype==_HALFLINE__VECT,segment= G.subtype==_GROUP__VECT;
      if (G.type==_VECT && G._VECTptr->size()>=2 && (line || halfline || segment)){
	for (int n=1;n<G._VECTptr->size();++n){
	  gen a=evalf_double((*G._VECTptr)[n-1],1,contextptr),b=evalf_double((*G._VECTptr)[n],1,contextptr);
	  if (a.type==_VECT && b.type==_VECT && a._VECTptr->size()==3 && b._VECTptr->size()==3){
	    vecteur & A=*a._VECTptr;
	    vecteur & B=*b._VECTptr;
	    if (A[0].type==_DOUBLE_ && A[1].type==_DOUBLE_ && A[2].type==_DOUBLE_ && B[0].type==_DOUBLE_ && B[1].type==_DOUBLE_ && B[2].type==_DOUBLE_ ){
	      lines.push_back(ptr);
	      double x=A[0]._DOUBLE_val,y=A[1]._DOUBLE_val,z=A[2]._DOUBLE_val;
#if 0 // ndef OLD_LINE_RENDERING
	      double3 prev(x,y,z);
	      linev.push_back(prev);
#endif
	      double X,Y,Z;
	      do_transform(transform3d,x,y,z,X,Y,Z);
	      double3 M(X,Y,Z);
	      x=B[0]._DOUBLE_val;y=B[1]._DOUBLE_val;z=B[2]._DOUBLE_val;
#if 0 // ndef OLD_LINE_RENDERING
	      linev.push_back(double3(x-prev.x,y-prev.y,z-prev.z));
#endif
	      do_transform(transform3d,x,y,z,X,Y,Z);
	      double3 N(X,Y,Z);
	      double3 v(N.x-M.x,N.y-M.y,N.z-M.z);
	      linev.push_back(M); linev.push_back(v);
	      linetypev.push_back(G.subtype);
	      line_color.push_back(int4(u,d,du,dd));
	    }
	  }
	}
	continue;
      }
      if (G.is_symb_of_sommet(at_hypersphere)){
	solid3d=true;
	vecteur hyp=*G._SYMBptr->feuille._VECTptr;
	gen c=evalf_double(hyp[0],1,contextptr);
	double x=c[0]._DOUBLE_val,y=c[1]._DOUBLE_val,z=c[2]._DOUBLE_val,X,Y,Z;
	do_transform(transform3d,x,y,z,X,Y,Z);
	sphere_centerv.push_back(double3(X,Y,Z));
	gen R=evalf(hyp[1],1,contextptr);
	double r=R._DOUBLE_val;
	sphere_radiusv.push_back(r);
	double * mat=invtransform3d;
	matrice qmat(makevecteur(
				 makevecteur(mat[0],mat[4],mat[8]),
				 makevecteur(mat[1],mat[5],mat[9]),
				 makevecteur(mat[2],mat[6],mat[10])
				 ));
	qmat=mmult(mtran(qmat),qmat);
	sphere_quadraticv.push_back(qmat);
	sphere_color.push_back(int4(u,d,du,dd));
	bool isclipped=x>=window_xmin+r && x<=window_xmax-r && y>=window_ymin+r && y<=window_ymax-r && z>=window_zmin+r && z<=window_zmax-r;
	// check if distance of center to window_x/y/xmin/max is <=R
	sphere_isclipped.push_back(isclipped);
	continue;
      }
      if (G.is_symb_of_sommet(at_hyperplan)){
	plan_filled.push_back(fill_polyedre);
	vecteur hyp=*G._SYMBptr->feuille._VECTptr;
	gen hyp1=evalf_double(hyp[1],1,contextptr);
	vecteur & hyp1v=*hyp1._VECTptr;
	double A,B,C,x,y,z,X,Y,Z;
	A=hyp1v[0]._DOUBLE_val;B=hyp1v[1]._DOUBLE_val;C=hyp1v[2]._DOUBLE_val;
	do_transform(transform3d,A,B,C,X,Y,Z);
	gen hyp0=evalf_double(hyp[0],1,contextptr);
	vecteur & hyp0v=*hyp0._VECTptr;
	x=hyp0v[0]._DOUBLE_val;y=hyp0v[1]._DOUBLE_val;z=hyp0v[2]._DOUBLE_val;
	double * mat=invtransform3d;
	A=mat[0]*x+mat[4]*y+mat[8]*z;
	B=mat[1]*x+mat[5]*y+mat[9]*z;
	C=mat[2]*x+mat[6]*y+mat[10]*z;
	double AB=std::abs(A)+std::abs(B);
	if (std::abs(C)<1e-2*AB){
	  // almost vertical plane, equation A*(x-X)+B*(y-Y)=0
	  continue;
	}
	A/=C; B/=C;
	plan_pointv.push_back(double3(X,Y,Z));	
	plan_abcv.push_back(double3(-A,-B,Z+A*X+B*Y));
	plan_color.push_back(int4(u,d,du,dd));
	continue;
      }
      if (G.is_symb_of_sommet(at_hypersurface)){
	solid3d=true;
	const vecteur & hyp=*G._SYMBptr->feuille._VECTptr;
	gen hyp0=hyp[0];
	const vecteur & hyp0v=*hyp0._VECTptr;
	gen h=hyp0v[4];
	if (h.type==_VECT && h.subtype==_POLYEDRE__VECT)
	  G=h;
	else if (ckmatrix(h,true)){
	  bool cplx=has_i(h); // 4d hypersurface, encode color in a float+int
	  double argcplx;
	  surfacev.push_back(vector< vector<float3d> >(0));
	  vector< vector<float3d> > & S=surfacev.back();
	  const vecteur & V=*h._VECTptr;
	  S.reserve(V.size());
	  for (int j=0;j<V.size();++j){
	    gen Vj=V[j];
	    const vecteur & vj=*Vj._VECTptr;
	    S.push_back(vector<float3d>(0));
	    if (vj.empty())
	      continue;
	    vector<float3d> &S_=S.back();
	    if (vj[0].type==_VECT && vj[0]._VECTptr->size()==3){
	      S_.reserve(3*vj.size());
	      for (int k=0;k<vj.size();++k){
		if (vj[k].type!=_VECT || vj[k]._VECTptr->size()!=3)
		  continue;
		const vecteur & cur=*vj[k]._VECTptr;
		double X,Y,Z;
		if (cplx)
		  do_transform(mat,cur[0]._DOUBLE_val,cur[1]._DOUBLE_val,absarg(cur[2],argcplx),X,Y,Z);
		else
		  do_transform(mat,cur[0]._DOUBLE_val,cur[1]._DOUBLE_val,cur[2]._DOUBLE_val,X,Y,Z);		
		S_.push_back(X); S_.push_back(Y);
		if (cplx){
		  float2 * fptr=(float2 *) &Z;
		  fptr->f = Z;
		  fptr->a = argcplx;
		  S_.push_back(Z);
		}
		else
		  S_.push_back(Z);
	      }
	    }
	    else {
	      S_.reserve(vj.size());
	      for (int k=0;k<vj.size();k+=3){
		double X,Y,Z;
		do_transform(mat,vj[k]._DOUBLE_val,vj[k+1]._DOUBLE_val,vj[k+2]._DOUBLE_val,X,Y,Z);
		//vj[k]=X; vj[k+1]=Y; vj[k+2]=Z;
		S_.push_back(X); S_.push_back(Y); S_.push_back(Z);
	      }
	    }
	  }
	  hyp_color.push_back(cplx?int4(0,0,0,0):int4(u,d,du,dd));
	  continue;
	} // end quad hypersurface
      } // end hypersurface
      if (G.type==_VECT && G.subtype==_GROUP__VECT){
	G=gen(makevecteur(G),_POLYEDRE__VECT);
      }
      if (G.type==_VECT && G.subtype==_POLYEDRE__VECT){
	vecteur p=*G._VECTptr;
	polyedrev.reserve(polyedrev.size()+p.size());
	polyedre_color.reserve(polyedre_color.size()+p.size());
	polyedre_xyminmax.reserve(polyedre_xyminmax.size()+2*p.size());
	for (int j=0;j<p.size();++j){
	  bool is_clipped=false;
	  gen g=p[j];
	  if (g.type==_VECT){
	    vector<double3> cur;
	    vecteur w=*g._VECTptr;
	    cur.reserve(w.size()+1);
	    for (int k=0;k<w.size();++k){
	      gen P=evalf_double(w[k],1,contextptr);
	      if (P.type==_VECT && P._VECTptr->size()==3){
		double x=P[0]._DOUBLE_val,y=P[1]._DOUBLE_val,z=P[2]._DOUBLE_val;
		if (is_clipped && (x<window_xmin || x>window_xmax || y<window_ymin || y>window_ymax || z<window_zmin || z>window_zmax) )
		  is_clipped=false;
		double X,Y,Z;
		do_transform(transform3d,x,y,z,X,Y,Z);
		cur.push_back(double3(X,Y,Z));
	      }
	    }
	    if (cur.size()>=3){
	      double Z3;int l=0;
	      for (;l<cur.size()-2;++l){
		double x0=cur[l].x,y0=cur[l].y,z0=cur[l].z;
		double x1=cur[l+1].x,y1=cur[l+1].y,z1=cur[l+1].z;
		double x2=cur[l+2].x,y2=cur[l+2].y,z2=cur[l+2].z;
		double X1=x1-x0,Y1=y1-y0,Z1=z1-z0;
		double X2=x2-x0,Y2=y2-y0,Z2=z2-z0;
		double X3=Y1*Z2-Y2*Z1,Y3=X2*Z1-X1*Z2;Z3=X1*Y2-X2*Y1;
		double prec=1e-3;
#if 1
		if ( (X3!=0 || Y3!=0) && std::abs(Z3)<prec*(std::abs(X3)+std::abs(Y3))){
		  Z3=1.0001*prec*(std::abs(X3)+std::abs(Y3));
		}
#endif
		if (std::abs(Z3)>prec*(std::abs(X3)+std::abs(Y3))){
		  // X3*(x-x0)+Y3*(y-y0)+Z3*(z-z0)=0
		  X3/=Z3; Y3/=Z3;
		  polyedre_abcv.push_back(double3(-X3,-Y3,z0+X3*x0+Y3*y0));
		  break;
		}
	      }
	      if (l==cur.size()-2)
		continue;
	      cur.push_back(cur.front());
	      double facemin=1e306,facemax=-1e306;
	      if (fill_polyedre){
		for (int l=1;l<cur.size();++l){
		  double xy=cur[l].y-cur[l].x;
		  if (xy<facemin)
		    facemin=xy;
		  if (xy>facemax)
		    facemax=xy;
		  // replace unused z coordinate by slope
		  // cur[l].z=(cur[l].y-cur[l-1].y)/(cur[l].x-cur[l-1].x);
		}
	      }
	      polyedrev.push_back(vector<double3>(0)); polyedrev.back().swap(cur); // polyedrev.push_back(cur);
	      polyedre_color.push_back(int4(u,d,du,dd));
	      polyedre_xyminmax.push_back(facemin);
	      polyedre_xyminmax.push_back(facemax);
	      polyedre_faceisclipped.push_back(is_clipped);
	      polyedre_filled.push_back(fill_polyedre);
	    } // end cur.size()>=3
	  } // end g.type==_VECT
	}
	continue;
      }      
    }
  }

  void Graph3d::displaypolyg(const std::vector<int2> & polyg,const int2 & IJmin,int color,int & Px,int & Py) const {
    if (polyg.empty())
      return;
    // sort list of arguments
    vector<int2_double2> p;
    for (int k=0;k<polyg.size();++k){
      const int2 & cur=polyg[k];
      if (cur==IJmin){
	int2_double2 id={cur.i,cur.j,0,0};
	p.push_back(id);
      } else {
	double di=cur.i-IJmin.i,dj=cur.j-IJmin.j;
	int2_double2 id={cur.i,cur.j,std::atan2(di,dj),di*di+dj*dj};
	p.push_back(id);
      }
    }
    sort(p.begin(),p.end());
    // draw polygon
    vector< vector<int> > P;
    for (int k=0;k<p.size();++k){
      vector<int> vi(2);
      vi[0]=p[k].i;
      vi[1]=p[k].j;
      P.push_back(vi);
    }
    xcas_color(color & 0xffff);
    int Ps=P.size();
    for (int i=0;i<Ps;++i){
      const vector<int> & A=P[(i?i:Ps)-1];
      const vector<int> & B=P[i];
      fl_line(x()+A[0],y()+A[1],x()+B[0],y()+B[1]);
    }
    Px=P[0][0];
    Py=P[0][1];
  }

  void Graph3d::adddepth(vector<int2> & polyg,const double3 &A,const double3 &B,int2 & IJmin) const {
    if ((A.z-current_depth)*(B.z-current_depth)>0)
      return;
    double t=(current_depth-A.z)/(B.z-A.z);
    double x=A.x+t*(B.x-A.x);
    double y=A.y+t*(B.y-A.y);
    int I,J;
    XYZ2ij(double3(x,y,current_depth),I,J);
    int2 IJ(I,J);
    polyg.push_back(IJ);
    if (IJ<IJmin)
      IJmin=IJ;
  }

  void Graph3d::addpolyg(vector<int2> & polyg,double x,double y,double z,int2 & IJmin) const {
    int I,J;
    xyz2ij(double3(x,y,z),I,J);
    int2 IJ(I,J);
    polyg.push_back(IJ);
    if (IJ<IJmin)
      IJmin=IJ;
  }

  void normalize(double & x, double & y){
    double fx=fabs(x),fy=fabs(y);
    if (fx>fy){
      y/=fx;
      x/=fx;
    }
    else {
      x/=fy;
      y/=fy;
    }
    double n=std::sqrt(x*x+y*y);
    x/=n; y/=n;
  }

  int diffuse(int color_orig,double diffusionz){
    if (diffusionz<1.1)
      return color_orig;
    int color=rgb565to888(color_orig);
    int r=(color&0xff0000)>>16,g=(color & 0xff00)>>8,b=192;
    double attenuate=(.3*(diffusionz-1));
    attenuate=1.0/(1+attenuate*attenuate);
    r*=attenuate; g*=attenuate; b*=attenuate;
    return rgb888to565((r<<16)|(g<<8)|b);
  }
  
  void Graph3d::glinter1(double z,double dz,
		double *zmin,double *zmax,double ZMIN,double ZMAX,
		int ih,int lcdz,
		int upcolor,int downcolor,int diffusionz,int diffusionz_limit,bool interval
		) const {
    if (ZMIN<z && z<ZMAX)
      return;
    // lcdz tests below: avoid marking too large regions
    if (*zmax<*zmin || z<*zmin-lcdz || z>*zmax+lcdz)
      *zmax=*zmin=z;
    bool diffus=diffusionz<diffusionz_limit;
    double deltaz;
    if (interval){
      bool intervalonly=false;
      if (z<0) {
	// return;
	z=0; intervalonly=true;
      }
      if (z>=LCD_HEIGHT_PX) {
	// return;
	z=LCD_HEIGHT_PX-1; intervalonly=true;
      }
      deltaz=diffus?1:diffusionz;
      if (z>*zmax+deltaz){
	if (diffus){
	  drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),diffuse(downcolor,diffusionz));
	  if (!intervalonly)
	    os_set_pixel(ih,z,downcolor);
	}
	else {
	  drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),_BLACK);
	  // draw interval
	  int nstep=int(z-*zmax)/diffusionz;
	  double zstep=(z-*zmax)/nstep;
	  for (double zz=*zmax+zstep;zz<=z;zz+=zstep)
	    os_set_pixel(ih,zz,downcolor);
	}
	*zmax=z;
	return;
      }
      else if (z<*zmin-deltaz){
	if (diffus){
	  drawRectangle(ih,z,1,std::ceil(*zmin-z),diffuse(upcolor,diffusionz));
	  if (!intervalonly)
	    os_set_pixel(ih,z,upcolor);
	}
	else {
	  drawRectangle(ih,z,1,std::ceil(*zmin-z),_BLACK);
	  // draw interval
	  int nstep=int(*zmin-z)/diffusionz;
	  double zstep=(z-*zmin)/nstep; // zstep<0
	  for (double zz=*zmin+zstep;zz>=z;zz+=zstep)
	    os_set_pixel(ih,zz,upcolor);
	}
	*zmin=z;
	return;
      }
    } // end if interval
    if (z>=0 &&  z<=LCD_HEIGHT_PX){
      int color=-1;
      if (diffus){
	if (z<=*zmin){
	  // mark all points with diffuse color from upcolor
	  drawRectangle(ih,z,1,std::ceil(*zmin-z),diffuse(upcolor,std::min(double(diffusionz),std::max(-dz,1.0))));
	  color=upcolor;
	  *zmin=z;
	}
	if (z>=*zmax){
	  // mark all points with diffuse color from downcolor
	  drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),diffuse(downcolor,std::min(double(diffusionz),std::max(dz,1.0))));
	  *zmax=z;
	}
	return;
      }
      if (z>*zmax){ // mark only 1 point
	color=downcolor;
	drawRectangle(ih,*zmax+1,1,z-*zmax-1,_BLACK);
	*zmax=z;
      }
      if (z<*zmin){ // mark 1 point
	color=upcolor;
	// drawRectangle(ih,z+1,1,*zmin-z-1,_BLACK);
	*zmin=z;
      }
      if (color>=0) os_set_pixel(ih,z,color);
    }
  }
  
  void Graph3d::glinter(double a,double b,double c,
	       double xscale,double xc,double yscale,double yc,
	       double *zmin,double *zmax,double ZMIN,double ZMAX,
	       int i,int horiz,int j,int w,int h,int lcdz,
	       int upcolor,int downcolor,int diffusionz,int diffusionz_limit,bool interval
	       ) const {
    double dz=lcdz*(a+b)*yscale-1;
    //if (dz<-10 || dz>10) cout << "dz=" << dz << "\n";
    // plane equation solved
    if (//0
	h==1 && w==1
	){
      int ih=i+horiz;
      double x = yscale*j-xscale*i + xc;
      // if (x<xmin) continue;
      double y = yscale*j+xscale*i + yc;
      // if (y<ymin) continue;
      double z = (a*x+b*y+c);
      z=LCD_HEIGHT_PX/2+j-lcdz*z;
      if (ZMIN<z && z<ZMAX)
	return;
      bool intervalonly=false;
      // lcdz tests below: avoid marking too large regions
      if (*zmax<*zmin || z<*zmin-lcdz || z>*zmax+lcdz)
	*zmax=*zmin=z;
      if (0 && (*zmax<50 || *zmin<50 || z<50))
	cout << *zmax << " "; // debug
      bool diffus=diffusionz<diffusionz_limit;
      double deltaz;
      if (interval){
	if (z<0) {
	  // return;
	  z=0; intervalonly=true;
	}
	if (z>=LCD_HEIGHT_PX) {
	  // return;
	  z=LCD_HEIGHT_PX-1; intervalonly=true;
	}
	deltaz=diffus?1:diffusionz;
	if (z>*zmax+deltaz){
	  if (diffus){
	    drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),diffuse(downcolor,diffusionz));
	    if (!intervalonly)
	      os_set_pixel(ih,z,downcolor);
	  }
	  else {
	    drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),_BLACK);
	    // draw interval
	    int nstep=int(z-*zmax)/diffusionz;
	    double zstep=(z-*zmax)/nstep;
	    for (double zz=*zmax+zstep;zz<=z;zz+=zstep)
	      os_set_pixel(ih,zz,downcolor);
	  }
	  *zmax=z;
	  return;
	}
	else if (z<*zmin-deltaz){
	  if (diffus){
	    drawRectangle(ih,z,1,std::ceil(*zmin-z),diffuse(upcolor,diffusionz));
	    if (!intervalonly)
	      os_set_pixel(ih,z,upcolor);
	  }
	  else {
	    drawRectangle(ih,z,1,std::ceil(*zmin-z),_BLACK);
	    // draw interval
	    int nstep=int(*zmin-z)/diffusionz;
	    double zstep=(z-*zmin)/nstep; // zstep<0
	    for (double zz=*zmin+zstep;zz>=z;zz+=zstep)
	      os_set_pixel(ih,zz,upcolor);
	  }
	  *zmin=z;
	  return;
	}
      } // end if interval
      if (z>=0 &&  z<=LCD_HEIGHT_PX){
	int color=-1;
	if (diffus){
	  if (z<=*zmin){
	    // mark all points with diffuse color from upcolor
	    drawRectangle(ih,z,1,std::ceil(*zmin-z),diffuse(upcolor,std::min(double(diffusionz),std::max(-dz,1.0))));
	    color=upcolor;
	    *zmin=z;
	  }
	  if (z>=*zmax){
	    // mark all points with diffuse color from downcolor
	    drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),diffuse(downcolor,std::min(double(diffusionz),std::max(dz,1.0))));
	    *zmax=z;
	  }
	  return;
	}
	if (z>=*zmax){ // mark only 1 point
	  color=downcolor;
	  drawRectangle(ih,*zmax+1,1,z-*zmax-1,_BLACK);
	  *zmax=z;
	}
	if (z<=*zmin){ // mark 1 point
	  color=upcolor;
	// drawRectangle(ih,z+1,1,*zmin-z-1,_BLACK);
	  *zmin=z;
	}
	if (color>=0) os_set_pixel(ih,z,color);
      }
      return; // end h==1 and w==1
    }
    for (int I=i;I<i+w;++I,++zmax,++zmin){
      int ih=I+horiz;
      double x = yscale*j-xscale*I + xc;
      // if (x<xmin) continue;
      double y = yscale*j+xscale*I + yc;
      // if (y<ymin) continue;
      double z = (a*x+b*y+c);
      bool intervalonly=false;
      z=LCD_HEIGHT_PX/2+j-lcdz*z;
      if (ZMIN<z && z<ZMAX)
	return;
      if (interval){
	if (z<0) {
	  z=0; intervalonly=true;
	}
	if (z>=LCD_HEIGHT_PX) {
	  z=LCD_HEIGHT_PX-1; intervalonly=true;
	}
      }
      if (i==0)
	; // cout << "i=" << i << " j=" << j << ", zmin=" << *zmin << " z=" << z << " zmax=" << *zmax << " dz=" << dz << ", a=" << a << " b=" << b << " c=" << c <<"\n";
      if (*zmax<*zmin || z<*zmin-lcdz || z>*zmax+lcdz)
	*zmax=*zmin=z;
      int deltaz=(diffusionz<diffusionz_limit?1:diffusionz);
      if (interval && z>*zmax+deltaz){
	if (//0
	    diffusionz<diffusionz_limit
	    )	  
	  drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),diffuse(downcolor,diffusionz));
	else {
	  drawRectangle(ih,*zmax,1,std::ceil(z-*zmax),_BLACK);
	  // draw interval
	  int nstep=int(z-*zmax)/diffusionz;
	  double zstep=(z-*zmax)/nstep;
	  for (double zz=*zmax+zstep;zz<=z;zz+=zstep)
	    os_set_pixel(ih,zz,downcolor);
	}
	if (intervalonly){
	  *zmax=z;
	  continue;
	}
      }
      if (interval && z<*zmin-deltaz){
	if (//0
	    diffusionz<diffusionz_limit
	    )
	  drawRectangle(ih,z,1,std::ceil(*zmin-z),diffuse(upcolor,diffusionz));
	else {
	  drawRectangle(ih,z,1,std::ceil(*zmin-z),_BLACK);
	  // draw interval
	  int nstep=int(*zmin-z)/diffusionz;
	  double zstep=(z-*zmin)/nstep; // zstep<0
	  for (double zz=*zmin+zstep;zz>=z;zz+=zstep)
	    os_set_pixel(ih,zz,upcolor);
	}
	if (intervalonly){
	  *zmin=z;
	  continue;
	}
      }
      if (//x-(h-1)*yscale>xmin && y-(h-1)*yscale>ymin &&
	  z>=0 && z+(h-1)*dz>=0 && z<=LCD_HEIGHT_PX && z+(h-1)*dz<=LCD_HEIGHT_PX
	  ){
	int color=-1;
	if ( (h>1 || diffusionz<diffusionz_limit) && dz>0 && z>=*zmax){
	  // mark all points with downcolor
	  *zmax=z+(h-1)*dz;
	  color=downcolor;
	  if (diffusionz>=diffusionz_limit && dz>diffusionz){
	    drawRectangle(ih,z,1,std::ceil(*zmax-z),_BLACK);
	    // draw interval
	    int nstep=int(std::ceil((*zmax-z)/diffusionz));
	    double zstep=(*zmax-z)/nstep;
	    for (int i=0;i<=nstep;++i)
	      os_set_pixel(ih,z+i*zstep,color);		    
	    continue;
	  }
	}
	if ( (h>1 || diffusionz<diffusionz_limit) && dz<0 && z<=*zmin){
	  // mark all points with upcolor
	  *zmin=z+(h-1)*dz;
	  color=upcolor;
	  if (diffusionz>=diffusionz_limit && dz<-diffusionz){
	    drawRectangle(ih,z,1,std::ceil(z-*zmin),_BLACK);
	    // draw interval
	    int nstep=int(std::ceil((z-*zmin)/diffusionz));
	    double zstep=(*zmin-z)/nstep;
	    for (int i=0;i<=nstep;++i)
	      os_set_pixel(ih,z+i*zstep,color);
	    continue;
	  }
	}
	if (color>=0){
	  if (diffusionz<diffusionz_limit){
	    if (dz>0)
	      drawRectangle(ih,z,1,std::ceil(h*dz),diffuse(color,std::min(double(diffusionz),std::max(dz,1.0))));
	    else
	      drawRectangle(ih,z-std::ceil(-h*dz),1,std::ceil(-h*dz),diffuse(color,std::min(double(diffusionz),std::max(-dz,1.0))));
	    continue;
	  }
	  if (dz>1)
	    drawRectangle(ih,z,1,std::ceil(h*dz),_BLACK);
	  os_set_pixel(ih,z,color);
	  if (h==1) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==2) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==3) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==4) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==5) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==6) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==7) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  if (h==8) continue;
	  z += dz;
	  os_set_pixel(ih,z,color);
	  continue;
	}
	if (dz<=0 && z>=*zmax && z+(h-1)*dz>=*zmin){ // mark only 1 point
	  *zmax=z;
	  color=downcolor;
	}
	if (dz>=0 && z<=*zmin && z+(h-1)*dz<=*zmax){ // mark 1 point
	  *zmin=z;
	  color=upcolor;
	}
	if (color>=0){
	  os_set_pixel(ih,z,color);
	  continue;
	}
      }
      for (int J=0;J<h
	     // && x>=xmin && y>=ymin
	     ;++J,z+=dz,x-=yscale,y-=yscale){
	int color=-1;
	if (z>*zmax){
	  drawRectangle(i,*zmax,1,z-*zmax,_BLACK);
	  *zmax=z;
	  color=downcolor;
	}
	if (z<*zmin){
	  drawRectangle(i,z,1,*zmin-z,_BLACK);
	  *zmin=z;
	  color=upcolor;
	}
	if (z<=-0.5 || z>=LCD_HEIGHT_PX)
	  continue;
	if (color>=0)
	  os_set_pixel(ih,z,color); // drawRectangle(i,z,w,h,color);
      }
    }
  }

  void find_abc(double x1,double x2,double x3,
		double y1,double y2,double y3,
		double z1,double z2,double z3,
		double &a,double &b,double &c){
    // solve([a*x1+b*y1+c=z1,a*x2+b*y2+c=z2,a*x3+b*y3+c=z3],[a,b,c])
    // double d=(x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2);
    double d=(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
    if (d==0) return;
    d=1/d;
    double z12=z2-z1,z23=z3-z2,z31=z1-z3;
    //double a=(-y1*z2+y1*z3+y2*z1-y2*z3-y3*z1+y3*z2)/d;
    a=d*(y1*z23+y2*z31+y3*z12);
    // double b=(x1*z2-x1*z3-x2*z1+x2*z3+x3*z1-x3*z2)/d;
    b=-d*(x1*z23+x2*z31+x3*z12);
    //double c=(x1*y2*z3-x1*y3*z2-x2*y1*z3+x2*y3*z1+x3*y1*z2-x3*y2*z1)/d;
    c=d*(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1));
  }

  void Graph3d::glinter(double x1,double x2,double x3,
	       double y1,double y2,double y3,
	       double z1,double z2,double z3,
	       double xscale,double xc,double yscale,double yc,
	       double *zmin,double *zmax,double ZMIN,double ZMAX,
	       int i,int horiz,int j,int w,int h,int lcdz,
	       int upcolor,int downcolor,int diffusionz,int diffusionz_limit,bool interval
	       ) const{
    double a,b,c;
    find_abc(x1,x2,x3,y1,y2,y3,z1,z2,z3,a,b,c);
    glinter(a,b,c,xscale,xc,yscale,yc,zmin,zmax,ZMIN,ZMAX,i,horiz,j,w,h,lcdz,upcolor,downcolor,diffusionz,diffusionz_limit,interval);
  }
  
  void update12(bool & found,bool &found2,
		double x1,double x2,double x3,double y1,double y2,double y3,double z1,double z2,double z3,
		int upcolor,int downcolor,int downupcolor,int downdowncolor,
		double & curx1, double &curx2, double &curx3, double &cury1, double &cury2, double &cury3, double &curz1, double &curz2, double &curz3, 
		double &cur2x1, double &cur2x2, double &cur2x3, double &cur2y1, double &cur2y2, double &cur2y3, double &cur2z1, double &cur2z2, double &cur2z3,
		int & u,int & d,int & du,int & dd){
    if (found){
      if (z1+z2+z3<curz1+curz2+curz3){
	// no need to update cur, perhaps cur2?
	if (found2 && cur2z1+cur2z2+cur2z3<z1+z2+z3)
	  return;
	found2=true;
	cur2x1=x1; cur2x2=x2; cur2x3=x3;
	cur2y1=y1; cur2y2=y2; cur2y3=y3;
	cur2z1=z1; cur2z2=z2; cur2z3=z3;
	du=downupcolor; dd=downdowncolor;
	return;
      }
      else { // need to update cur, and perhaps cur2
	if (!found2 || curz1+curz2+curz3<cur2z1+cur2z2+cur2z3){
	  found2=true;
	  cur2x1=curx1; cur2x2=curx2; cur2x3=curx3;
	  cur2y1=cury1; cur2y2=cury2; cur2y3=cury3;
	  cur2z1=curz1; cur2z2=curz2; cur2z3=curz3;
	  du=downupcolor; dd=downdowncolor;
	}
	curx1=x1; curx2=x2; curx3=x3;
	cury1=y1; cury2=y2; cury3=y3;
	curz1=z1; curz2=z2; curz3=z3;
	u=upcolor; d=downcolor;
	return;
      }
    }
    found=true;
    curx1=x1; curx2=x2; curx3=x3;
    cury1=y1; cury2=y2; cury3=y3;
    curz1=z1; curz2=z2; curz3=z3;
    u=upcolor; d=downcolor;
  }
	      
  void update12(bool & found,bool &found2,
		double a,double b,double c,double z,
		int upcolor,int downcolor,int downupcolor,int downdowncolor,
		double & cura1,double & curb1,double & curc1,double & curz1,
		double & cura2,double & curb2,double & curc2,double & curz2,
		int & u,int & d,int & du,int & dd){
    if (!found || z>curz1){
      if (found){ // update cur2
	found2=true;
	cura2=cura1; curb2=curb1; curc2=curc1; curz2=curz1;
      }
      found=true;
      cura1=a; curb1=b; curc1=c; curz1=z;
      u=upcolor; d=downcolor;
      return;
    }
    if (z>curz2){
      found2=true;
      cura2=a; curb2=b; curc2=c; curz2=z;
      du=downupcolor; dd=downdowncolor;
    }
  }
	        
  // intersect plane x-y=xy with line m+t*v
  // m.x+t*v.x-m.y-t*v.y=-xy
  double intersect(const double3 & m,const double3 & v,double xy){
    return (-xy+m.y-m.x)/(v.x-v.y);
  }

  // 2d coordinates of m+t*v
  void grmtv2ij(const Graph3d & gr,const double3 & m,const double3 & v,double t,int & i,int & j){
    double x=m.x+t*v.x;
    double y=m.y+t*v.y;
    double z=m.z+t*v.z;
    gr.XYZ2ij(double3(x,y,z),i,j);
  }    
  
  /* generate colors
#include <iostream>

using namespace std;

// t angle in radians -> r,g,b
int arc_en_ciel(int k,double intensite){
  int r,g,b;
  k = (k+126)%126;
  if (k<21){
    r=251; g=0; b=12*k;
  }
  if (k>=21 && k<42){
    r=251-(12*(k-21)); g=0; b=251;
  } 
  if (k>=42 && k<63){
    r=0; g=(k-42)*12; b=251;
  } 
  if (k>=63 && k<84){
    r=0; g=251; b=251-(k-63)*12;
  } 
  if (k>=84 && k<105){
    r=(k-84)*12; g=251; b=0;
  } 
  if (k>=105 && k<126){
    r=251; g=251-(k-105)*12; b=0;
  }
  r*=intensite; g*=intensite; b*=intensite;
  if (r==0) r=8;
  return (((r*32)/256)<<11) | (((g*64)/256)<<5) | (b*32/256);
}


int main(){
  for (int i=0;i<126;++i){
    int u=arc_en_ciel(i,1),d=arc_en_ciel(i,.75),du=arc_en_ciel(i,.5),dd=arc_en_ciel(i,.25);
    cout << "{" << u << "," << d << "," << du << "," << dd << "},\n";
  }
}
  */

  const int4 tabcolorcplx[]={
{63488,47104,30720,14336},
{63489,47105,30720,14336},
{63491,47106,30721,14336},
{63492,47107,30722,14337},
{63494,47108,30723,14337},
{63495,47109,30723,14337},
{63497,47110,30724,14338},
{63498,47111,30725,14338},
{63500,47113,30726,14339},
{63501,47114,30726,14339},
{63503,47115,30727,14339},
{63504,47116,30728,14340},
{63506,47117,30729,14340},
{63507,47118,30729,14340},
{63509,47119,30730,14341},
{63510,47120,30731,14341},
{63512,47122,30732,14342},
{63513,47123,30732,14342},
{63515,47124,30733,14342},
{63516,47125,30734,14343},
{63518,47126,30735,14343},
{63519,47127,30735,14343},
{59423,45079,28687,14343},
{57375,43031,28687,14343},
{53279,40983,26639,12295},
{51231,38935,24591,12295},
{47135,34839,22543,10247},
{45087,32791,22543,10247},
{40991,30743,20495,10247},
{38943,28695,18447,8199},
{34847,26647,16399,8199},
{32799,24599,16399,8199},
{28703,22551,14351,6151},
{26655,20503,12303,6151},
{22559,16407,10255,4103},
{20511,14359,10255,4103},
{16415,12311,8207,4103},
{14367,10263,6159,2055},
{10271,8215,4111,2055},
{8223,6167,4111,2055},
{4127,4119,2063,7},
{2079,2071,15,7},
{2079,2071,2063,2055},
{2175,2135,2095,2055},
{2271,2199,2159,2087},
{2367,2263,2191,2119},
{2463,2359,2255,2151},
{2559,2423,2287,2151},
{2655,2487,2351,2183},
{2751,2551,2383,2215},
{2847,2647,2447,2247},
{2943,2711,2479,2247},
{3039,2775,2543,2279},
{3135,2839,2575,2311},
{3231,2935,2639,2343},
{3327,2999,2671,2343},
{3423,3063,2735,2375},
{3519,3127,2767,2407},
{3615,3223,2831,2439},
{3711,3287,2863,2439},
{3807,3351,2927,2471},
{3903,3415,2959,2503},
{3999,3511,3023,2535},
{4063,3575,3055,2535},
{4061,3574,3054,2535},
{4060,3573,3054,2535},
{4058,3572,3053,2534},
{4057,3571,3052,2534},
{4055,3569,3051,2533},
{4054,3568,3051,2533},
{4052,3567,3050,2533},
{4051,3566,3049,2532},
{4049,3565,3048,2532},
{4048,3564,3048,2532},
{4046,3563,3047,2531},
{4045,3562,3046,2531},
{4043,3560,3045,2530},
{4042,3559,3045,2530},
{4040,3558,3044,2530},
{4039,3557,3043,2529},
{4037,3556,3042,2529},
{4036,3555,3042,2529},
{4034,3554,3041,2528},
{4033,3553,3040,2528},
{4032,3552,3040,2528},
{4032,3552,992,480},
{8128,5600,3040,480},
{10176,7648,5088,2528},
{14272,9696,7136,2528},
{16320,11744,7136,2528},
{20416,13792,9184,4576},
{22464,15840,11232,4576},
{26560,19936,13280,6624},
{28608,21984,13280,6624},
{32704,24032,15328,6624},
{34752,26080,17376,8672},
{38848,28128,19424,8672},
{40896,30176,19424,8672},
{44992,32224,21472,10720},
{47040,34272,23520,10720},
{51136,38368,25568,12768},
{53184,40416,25568,12768},
{57280,42464,27616,12768},
{59328,44512,29664,14816},
{63424,46560,31712,14816},
{65472,48608,31712,14816},
{65376,48512,31648,14784},
{65280,48448,31616,14784},
{65184,48384,31552,14752},
{65088,48320,31520,14720},
{64992,48224,31456,14688},
{64896,48160,31424,14688},
{64800,48096,31360,14656},
{64704,48032,31328,14624},
{64608,47936,31264,14592},
{64512,47872,31232,14592},
{64416,47808,31168,14560},
{64320,47744,31136,14528},
{64224,47648,31072,14496},
{64128,47584,31040,14496},
{64032,47520,30976,14464},
{63936,47456,30944,14432},
{63840,47360,30880,14400},
{63744,47296,30848,14400},
{63648,47232,30784,14368},
{63552,47168,30752,14336},
  };
  
  struct hypertriangle_t {
    const int4 * colorptr; // hypersurface color 
    double xmin,xmax,ymin,ymax; // minmax values intersection with plane y-x=Cte
    double a,b,c; // plane equation of triangle
    double zG; // altitude for gravity center 
  }  ; // data struct for hypesurface triangulation cache

  void compute(double yx,double3 * cur,hypertriangle_t & res){
    double xmin=1e307,xmax=-1e307,ymin=1e307,ymax=-1e307;
    for (int l=0;l<4;++l){
      int prev=l==0?3:l-1;
      double3 & d3=cur[prev];
      double x0=d3.x,y0=d3.y,x1=cur[l].x,y1=cur[l].y;
      double yx0=y0-x0,yx1=y1-x1,m=yx1-yx0;
      if (m==0){
	if (yx==yx1){
	  if (x0>xmax) xmax=x0; if (x0<xmin) xmin=x0;
	  if (x1>xmax) xmax=x1; if (x1<xmin) xmin=x1;
	  if (y0>ymax) ymax=y0; if (y0<ymin) ymin=y0;
	  if (y1>ymax) ymax=y1; if (y1<ymin) ymin=y1;
	}
	continue;
      }
      double t=(yx-yx0)/m;
      if (t>=0 && t<=1){
	double X=x0+t*(x1-x0),Y=y0+t*(y1-y0);
	if (X>xmax) xmax=X; if (X<xmin) xmin=X;
	if (Y>ymax) ymax=Y; if (Y<ymin) ymin=Y;
      }
    }
    res.zG=(cur[0].z+cur[1].z+cur[2].z+cur[3].z)/4;
    res.xmin=xmin; res.xmax=xmax; res.ymin=ymin; res.ymax=ymax;
    find_abc(cur[0].x,cur[1].x,cur[2].x,
	     cur[0].y,cur[1].y,cur[2].y,
	     cur[0].z,cur[1].z,cur[2].z,
	     res.a,res.b,res.c);
  }

  
  
  void update_hypertri(const vector<hypertriangle_t> & hypertriangles,double x,double y,
		       bool & found,bool &found2,
		       double3 & curabc1,double & curz1,
		       double3 & curabc2,double & curz2,
		       int & upcolor,int & downcolor,int & downupcolor,int & downdowncolor){
    vector<hypertriangle_t>::const_iterator it=hypertriangles.begin(),itend=hypertriangles.end();
    for (;it!=itend;++it){
      if (x<it->xmin){
	++it;
	if (it==itend) break;
	if (x<it->xmin){
	  ++it;
	  if (it==itend) break;
	  if (x<it->xmin){
	    ++it;
	    if (it==itend) break;
	  }
	}
      }
      else if (x>it->xmax){
	++it;
	if (it==itend) break;
	if (x>it->xmax){
	  ++it;
	  if (it==itend) break;
	  if (x>it->xmax){
	    ++it;
	    if (it==itend) break;
	  }
	}
      }
      const hypertriangle_t & cur=*it;
      if (x<cur.xmin || x>cur.xmax || y<cur.ymin || y>cur.ymax)
	continue;
      if (!found || cur.zG>curz1){
	if (found){
	  found2=true;
	  curabc2=curabc1;
	  curz2=curz1;
	}
	found=true;
	curabc1.x=cur.a; curabc1.y=cur.b; curabc1.z=cur.c;
	curz1=cur.zG;
	upcolor=cur.colorptr->u; downcolor=cur.colorptr->d;
	continue;
      }
      if (cur.zG>curz2){
	found2=true;
	curabc2.x=cur.a; curabc2.y=cur.b; curabc2.z=cur.c;
	curz2=cur.zG;
	downupcolor=cur.colorptr->du; downdowncolor=cur.colorptr->dd;
	continue;
      }
    } // end loop on k
  }
  
  bool is_inside(const vector<double3> & v,double x,double y){
    int n=0;
    for (int i=1;i<v.size();++i){
      const double3 & prev=v[i-1];
      const double3 & cur=v[i];
      double prevx=prev.x,prevy=prev.y,curx=cur.x,cury=cur.y,m=cur.z;
      if (prevx==curx){
	if (x==curx && (y-prevy)*(cury-y)>0) // on vertical edge
	  return false;
	continue;
      }
      if (x==prevx) 
	continue;
      if ((x-prevx)*(curx-x)<0)
	continue;
      //double Y=cury+m*(x-curx);
      double Y=cury+(cur.y-prev.y)/(cur.x-prev.x)*(x-curx);
      if (Y>=y)
	++n; 
    }
    if (n%2)
      return true;
    return false;
  }

  void Graph3d::draw_decorations(const giac::gen & title_tmp){
    if (args_tmp.empty()){ // add selected names
      char s[256]; strcpy(s,modestr.c_str());
      int pos=0,modestrsize=modestr.size();
      pos += modestrsize;
      s[pos++]=' ';
      if (mode!=0 && drag_name.type==_IDNT){
	strcpy(s+pos,drag_name._IDNTptr->id_name);
	pos += strlen(drag_name._IDNTptr->id_name);
      }
      else {
	if (1 || mode==0 || mode==255){ // print selected names 
	  vecteur v;
	  if (mode!=255) v=selected_names(true,false);
	  int vs=v.size();
	  if (!vs){ // print current coordinates
	    double i=Fl::event_x()-x(),j=Fl::event_y()-y(),x,y,z;
	    find_xyz(i,j,current_depth,x,y,z);
	    round3(z,window_zmin,window_zmax);
	    sprintf(s+pos," %.3g,%.3g,%.3g",x,y,z);
	    pos=strlen(s);
	  }
	  for (int i=0;i<vs && pos<100;++i){
	    if (v[i].type==_IDNT){
	      strcpy(s+pos,v[i]._IDNTptr->id_name);
	      pos += strlen(v[i]._IDNTptr->id_name);
	      if (i<vs-1){
		s[pos++]=',';
	      }
	    }
	  }
	}
      }
      s[pos]=0;
      os_draw_string_small(LCD_WIDTH_PX-fl_width(s),LCD_HEIGHT_PX-14,FL_WHITE,FL_BLACK,s,false);
    }
    if (mode && title.size()<100 && (!title.empty() || !is_zero(title_tmp))){
      std::string mytitle;
      if (!is_zero(title_tmp) && function_final.type==_FUNC){
	if (function_final==at_point && title_tmp.is_symb_of_sommet(at_point))
	  mytitle=title_tmp.print(context0);
	else
	  mytitle=function_final._FUNCptr->ptr()->s+('('+title_tmp.print(context0)+')'); // gen(symbolic(*function_final._FUNCptr,title_tmp)).print(contextptr);
      }
      else
	mytitle=title;
      if (!mytitle.empty()){
	int dt=int(fl_width(mytitle.c_str()));
	if (dt>LCD_WIDTH_PX)
	  dt=LCD_WIDTH_PX;
	os_draw_string_small((LCD_WIDTH_PX-dt)/2,LCD_HEIGHT_PX-14,FL_WHITE,FL_BLACK,mytitle.c_str(),false);
      }
    }
  }

  
  // hpersurface encoded as a matrix
  // with lines containing 3 coordinates per point
  bool Graph3d::indraw3d(int w,int h,int lcdz,GIAC_CONTEXT,
			  int upcolor_,int downcolor_,int downupcolor_,int downdowncolor_)  {
    if (w>9) w=9; if (w<1) w=1;
    if (h>9) h=9; if (h<1) h=1;
    // save zmin/zmax on the stack (4K required)
    const int jmintabsize=2048;
    short int *jmintab=(short int *)alloca(jmintabsize*sizeof(short int)), * jmaxtab=(short int *)alloca(jmintabsize*sizeof(short int)); // assumes LCD_WIDTH_PX<=jmintabsize
    for (int i=0;i<jmintabsize;++i){
      jmintab[i]=LCD_HEIGHT_PX;
      jmaxtab[i]=0;
    }
    vecteur attrv(native_renderer_instructions);
    std::vector< std::vector< vector<float3d> >::const_iterator > hypv; // 3 iterateurs per hypersurface
    int upcolor,downcolor,downupcolor,downdowncolor;
    for (int i=0;i<int(attrv.size());++i){
      gen attr=attrv[i];
      upcolor=upcolor_;downcolor=downcolor_;downupcolor=downupcolor_;downdowncolor=downdowncolor_;
      get_colors(attrv[i],upcolor,downcolor,downupcolor,downdowncolor);
    }
    for (int i=0;i<int(surfacev.size());++i){
      hypv.push_back(surfacev[i].begin());
      hypv.push_back(surfacev[i].end());
    }
    int horiz=LCD_WIDTH_PX/2,vert=horiz/2;//LCD_HEIGHT_PX/2;
    int imin=Ai,imax=Ai,itmp;
    // 12 segments from cube visualization
    double segments_x1[12]={Ai,Bi,Ei,Fi,Ai,Bi,Ci,Di,Ai,Ci,Ei,Gi};
    double segments_x2[12]={Ci,Di,Gi,Hi,Ei,Fi,Gi,Hi,Bi,Di,Fi,Hi};
    double segments_y1[12]={Aj,Bj,Ej,Fj,Aj,Bj,Cj,Dj,Aj,Cj,Ej,Gj};
    double segments_y2[12]={Cj,Dj,Gj,Hj,Ej,Fj,Gj,Hj,Bj,Dj,Fj,Hj};
    double segments_m[12];
    for (int i=0;i<12;++i){
      segments_m[i]=(segments_y2[i]-segments_y1[i])/(segments_x2[i]-segments_x1[i]);
      itmp=segments_x1[i];
      if (itmp<imin) imin=itmp; if (itmp>imax) imax=itmp;
      itmp=segments_x2[i];
      if (itmp<imin) imin=itmp; if (itmp>imax) imax=itmp;      
    }
    double xmin=-1,ymin=-1,xmax=1,ymax=1,xscale=0.6*(xmax-xmin)/horiz,yscale=0.6*(ymax-ymin)/vert,x,y,z,xc=(xmin+xmax)/2,yc=(ymin+ymax)/2;
    drawRectangle(0,0,imin,LCD_HEIGHT_PX,COLOR_BLACK); // clear    
    drawRectangle(imax,0,LCD_WIDTH_PX-imax,LCD_HEIGHT_PX,COLOR_BLACK); // clear
    sync_screen();
    int count=0;
    vector<int> polyedrei; polyedrei.reserve(polyedrev.size()); // cache for polyedres polygons edges
    vector<double> polyedrexmin,polyedrexmax,polyedreymin,polyedreymax;
    polyedrexmin.reserve(polyedrev.size());polyedrexmax.reserve(polyedrev.size());
    polyedreymin.reserve(polyedrev.size());polyedreymax.reserve(polyedrev.size());
    vector<hypertriangle_t> hypertriangles;
    for (int i=imin-horiz;i<imax-horiz;i+=w,++count){
    //for (int i=-horiz;i<horiz;i+=w){
      drawRectangle(i+horiz,0,w,LCD_HEIGHT_PX,COLOR_BLACK); // clear
      // find min and max values for j using vertical intersections of line x=i
      // with segments
      int ih=i+horiz+w/2,jmin=RAND_MAX,jmax=-RAND_MAX;
      for (int k=0;k<12;++k){
	if ( !is_inf(segments_m[k]) && (ih-segments_x1[k])*(segments_x2[k]-ih)>=0 ){
	  double y=segments_y1[k]+segments_m[k]*(ih-segments_x1[k]);
	  if (y<=jmin)
	    jmin=std::floor(y);
	  if (y>=jmax)
	    jmax=std::ceil(y);
	}
      }
      if (jmin>jmax) continue;
      if (jmin<0) jmin=0;
      if (jmax>LCD_HEIGHT_PX) jmax=LCD_HEIGHT_PX;
      jmin -= LCD_HEIGHT_PX/2;
      jmax -= LCD_HEIGHT_PX/2;      
      double yx=2*xscale*(i+(w-1)/2.0)+yc-xc;
      // poledrev indices for yx, and xmin/xmax/ymin/ymax values
      // (xmin/xmax should be enough, except limit cases)
      polyedrei.clear(); polyedrexmin.clear(); polyedrexmax.clear(); polyedreymin.clear(); polyedreymax.clear();
      for (int k=0;k<int(polyedrev.size());++k){
	double facemin=polyedre_xyminmax[2*k],facemax=polyedre_xyminmax[2*k+1];
	if (yx<facemin || yx>facemax)
	  continue;
	polyedrei.push_back(k);
	vector<double3> & cur=polyedrev[k];
	double xmin=1e307,xmax=-1e307,ymin=1e307,ymax=-1e307;
	for (int l=0;l<int(cur.size());++l){
	  double3 & d3=cur[l?l-1:cur.size()-1];
	  double x0=d3.x,y0=d3.y,x1=cur[l].x,y1=cur[l].y;
	  double yx0=y0-x0,yx1=y1-x1,m=yx1-yx0;
	  if (m==0){
	    if (yx==yx1){
	      if (x0>xmax) xmax=x0; if (x0<xmin) xmin=x0;
	      if (x1>xmax) xmax=x1; if (x1<xmin) xmin=x1;
	      if (y0>ymax) ymax=y0; if (y0<ymin) ymin=y0;
	      if (y1>ymax) ymax=y1; if (y1<ymin) ymin=y1;
	    }
	    continue;
	  }
	  double t=(yx-yx0)/m;
	  if (t>=0 && t<=1){
	    double X=x0+t*(x1-x0),Y=y0+t*(y1-y0);
	    if (X>xmax) xmax=X; if (X<xmin) xmin=X;
	    if (Y>ymax) ymax=Y; if (Y<ymin) ymin=Y;
	  }
	}
	polyedrexmin.push_back(xmin);
	polyedrexmax.push_back(xmax);
	polyedreymin.push_back(ymin);
	polyedreymax.push_back(ymax);
      }
      // hypersurfaces: find triangles
      hypertriangles.clear();
      double hyperxymax=-1e307,hyperxymin=1e307;
      double3 tri[4]; 
      for (int k=0;k<int(hypv.size());k+=2){
	bool cplx=hyp_color[k].u==0 && hyp_color[k].d==0 && hyp_color[k].du==0 && hyp_color[k].dd==0;
	vector< vector<float3d> >::const_iterator sbeg=hypv[k],send=hypv[k+1],sprec,scur;
	vector<float3d>::const_iterator itprec,itcur,itprecend;
	for (sprec=sbeg,scur=sprec+1;scur<send;++sprec,++scur){
	  itprec=sprec->begin(); 
	  itprecend=sprec->end();
	  itcur=scur->begin();
	  double yx1,yx2=*(itprec+1)-*itprec,yx3,yx4=*(itcur+1)-*itcur;
	  for (itprec+=3,itcur+=3;itprec<itprecend;itprec+=3,itcur+=3){
	    yx1=yx2;
	    yx2=*(itprec+1)-*itprec;
	    yx3=yx4;
	    yx4=*(itcur+1)-*itcur;
	    if (yx<yx1 && yx<yx2 && yx<yx3 && yx<yx4){
	      for (;;){
		// per iteration: 2 incr, 1 test, 2 read, 2 comp, && , test
		itprec+=3;itcur+=3;
		if (itprec<itprecend && yx<(yx2=*(itprec+1)-*itprec) && yx<(yx4=*(itcur+1)-*itcur)){
		  itprec+=3;itcur+=3;
		  if (itprec<itprecend && yx<(yx2=*(itprec+1)-*itprec) && yx<(yx4=*(itcur+1)-*itcur)){
		    itprec+=3;itcur+=3;
		    if (itprec<itprecend && yx<(yx2=*(itprec+1)-*itprec) && yx<(yx4=*(itcur+1)-*itcur)){
		      itprec+=3;itcur+=3;
		      if (itprec<itprecend && yx<(yx2=*(itprec+1)-*itprec) && yx<(yx4=*(itcur+1)-*itcur)){
			continue;
		      }
		    }
		  }
		}
		break;
	      }
	      if (yx<yx2 && yx<yx4) continue;
	    } // end yx<yxk
	    else if (yx>yx1 && yx>yx2 && yx>yx3 && yx>yx4){
	      for (;;){
		// per iteration: 2 incr, 1 test, 2 read, 2 comp, && , test
		itprec+=3;itcur+=3;
		if (itprec<itprecend && yx>(yx2=*(itprec+1)-*itprec) && yx>(yx4=*(itcur+1)-*itcur)){
		  itprec+=3;itcur+=3;
		  if (itprec<itprecend && yx>(yx2=*(itprec+1)-*itprec) && yx>(yx4=*(itcur+1)-*itcur)){
		    itprec+=3;itcur+=3;
		    if (itprec<itprecend && yx>(yx2=*(itprec+1)-*itprec) && yx>(yx4=*(itcur+1)-*itcur)){
		      itprec+=3;itcur+=3;
		      if (itprec<itprecend && yx>(yx2=*(itprec+1)-*itprec) && yx>(yx4=*(itcur+1)-*itcur)){
			continue;
		      }
		    }
		  }
		}
		break;
	      }
	      if (yx>yx2 && yx>yx4) continue;
	    }
	    // found one quad intersecting plane
	    double x1=*(itprec-3),x2=*(itprec),x3=*(itcur-3),x4=*(itcur);
	    double y1=*(itprec-2),y2=*(itprec+1),y3=*(itcur-2),y4=*(itcur+1);
	    double z1=*(itprec-1),z2=*(itprec+2),z3=*(itcur-1),z4=*(itcur+2);
	    double a1,a2,a3,a4;
	    if (cplx){
	      a1 = ((float2 *)&z1)->a;
	      z1 = ((float2 *)&z1)->f;
	      a2 = ((float2 *)&z2)->a;
	      z2 = ((float2 *)&z2)->f;
	      a3 = ((float2 *)&z3)->a;
	      z3 = ((float2 *)&z3)->f;
	      a4 = ((float2 *)&z4)->a;
	      z4 = ((float2 *)&z4)->f;
	    }
	    yx1=y1-x1; yx2=y2-x2; yx3=y3-x3; yx4=y4-x4;
	    tri[0]=double3(x1,y1,z1);
	    tri[1]=double3(x2,y2,z2);
	    tri[2]=double3(x4,y4,z4);
	    tri[3]=double3(x3,y3,z3);
	    double x123=(x1+x2+x3+x4)/4,y123=(y1+y2+y3+y4)/4,z123=(z1+z2+z3+z4)/4,X,Y,Z;
	    double xy123=x123+y123;
	    if (xy123<hyperxymin) hyperxymin=xy123;
	    if (xy123>hyperxymax) hyperxymax=xy123;
	    do_transform(invtransform3d,x123,y123,z123,X,Y,Z);
	    if ( Z>=window_zmin && Z<=window_zmax &&
		X>=window_xmin && X<=window_xmax && Y>=window_ymin && Y<=window_ymax ){
	      hypertriangle_t res; 
	      if (cplx){
		int idx=(a1+M_PI)*sizeof(tabcolorcplx)/(sizeof(int4)*2*M_PI);
		if (idx<0 || idx >=sizeof(tabcolorcplx)/(sizeof(int4))){
		  idx = 0;
		}
		// arg of z1*exp(i*a1)+z2*exp(i*a2)+z3*exp(i*a3)+z4*exp(i*a4)
		// does not work because z1,... are transformed
		//CERR << idx << " ";
		res.colorptr=&tabcolorcplx[idx];
	      }
	      else
		res.colorptr=&hyp_color[k];
	      compute(yx,tri,res);
	      hypertriangles.push_back(res);
	    }
	  }
	}
      }
      vector<int> spheres(sphere_centerv.size()); // is plane y-x=yx intersecting sphere, vector<bool> does not work with Keil
      for (int k=0;k<int(sphere_centerv.size());++k){
	const double3 & c=sphere_centerv[k];
	double xc=c.x,yc=c.y;
	double r=sphere_radiusv[k];
	const matrice & m=*sphere_quadraticv[k]._VECTptr;
	const vecteur & m0=*m[0]._VECTptr;
	const vecteur & m1=*m[1]._VECTptr;
	const vecteur & m2=*m[2]._VECTptr;
	double m00=m0[0]._DOUBLE_val,m01=m0[1]._DOUBLE_val,m02=m0[2]._DOUBLE_val,m11=m1[1]._DOUBLE_val,m12=m1[2]._DOUBLE_val,m22=m2[2]._DOUBLE_val;
	/* q0:=m00*x^2+2*m01*x*y+2*m02*x*z+m11*y^2+2*m12*y*z+m22*z^2; q:=subst(q0,[x,y],[x-xc,y-yc]);
	   a,b,c:=coeffs(q(y=yx+x)-r^2,x);
	   delta:=b^2-4*a*c; 
	   // if delta<0 for all z, there is no intersection
	   // delta is a second order polynomial in z, check discriminant
	   A,B,C:=coeffs(delta,z);
	   D:=B^2-4*A*C;  // if D<0 no intersection
	*/
	double A=4*m02*m02+4*m12*m12-4*m00*m22-8*m01*m22+8*m02*m12-4*m11*m22,
	  B=-8*m00*m12*xc+8*m00*m12*yc-8*m00*m12*yx+8*m01*m02*xc-8*m01*m02*yc+8*m01*m02*yx-8*m01*m12*xc+8*m01*m12*yc-8*m01*m12*yx+8*m02*m11*xc-8*m02*m11*yc+8*m02*m11*yx,
	  C=4*m01*m01*xc*xc+4*m01*m01*yc*yc+4*m01*m01*yx*yx-4*m00*m11*xc*xc-4*m00*m11*yc*yc-4*m00*m11*yx*yx-8*m01*m01*xc*yc+8*m01*m01*xc*yx-8*m01*m01*yc*yx+8*m00*m11*xc*yc-8*m00*m11*xc*yx+8*m00*m11*yc*yx+4*m00*r*r+8*m01*r*r+4*m11*r*r,
	  D=B*B-4*A*C;
	spheres[k]=D>=0;
      }
      double zmin[10]={220.220,220,220,220,220,220,220,220,220},
	zmax[10]={0,0,0,0,0,0,0,0,0,0},
	zmin2[10]={220.220,220,220,220,220,220,220,220,220},
	zmax2[10]={0,0,0,0,0,0,0,0,0,0}	; // initialize for these vertical lines
      double3 curabc1,curabc2; 
      double curz1=-1e306,curz2=1e306;
      int u,d,du,dd;
      // loop earlier if there are only hypersurfaces
      bool only_hypertri=true;
      for (int ki=0;ki<int(polyedrei.size());++ki){
	if (polyedrexmin[ki]<=polyedrexmax[ki]){ only_hypertri=false; break; }
      }
      for (int k=0;k<int(sphere_centerv.size());++k){
	if (spheres[k]){ only_hypertri=false; break; }
      }
      for (int k=0;k<int(plan_abcv.size());++k){
	if (plan_filled[k]){ only_hypertri=false; break; }
      }
      if (only_hypertri){
	if (hypertriangles.empty()) 
	  goto suite3d;
	int effjmax=(hyperxymax-xc-yc)/yscale/2.0,effjmin=(hyperxymin-xc-yc)/yscale/2.0;
	if (effjmax+1<jmax)
	  jmax=effjmax+1;
	if (effjmin-1>jmin)
	  jmin=effjmin-1;
	x = yscale*(jmax-(h-1)/2.0)-xscale*(i+(w-1)/2.0) + xc;
	y = yscale*(jmax-(h-1)/2.0)+xscale*(i+(w-1)/2.0) + yc;
	for (int j=jmax;j>=jmin;j-=h,x-=yscale*h,y-=yscale*h){
	  bool found=false,found2=false;
	  update_hypertri(hypertriangles,x,y,found,found2,curabc1,curz1,curabc2,curz2,upcolor,downcolor,downupcolor,downdowncolor);
	  if (!found) continue;
	  if (h==1 && w==1){
	    if (found2 && !hide2nd){
	      double dz=lcdz*(curabc2.x+curabc2.y)*yscale-1;
	      // if (y<ymin) continue;
	      double z = (curabc2.x*x+curabc2.y*y+curabc2.z);
	      z=LCD_HEIGHT_PX/2+j-lcdz*z;
	      glinter1(z,dz,
		       zmin2,zmax2,zmin[0],zmax[0],
		       ih,lcdz,
		       downupcolor,downdowncolor,diffusionz,diffusionz_limit,interval);
	    }
	    double dz=lcdz*(curabc1.x+curabc1.y)*yscale-1;
	    // if (y<ymin) continue;
	    double z = (curabc1.x*x+curabc1.y*y+curabc1.z);
	    z=LCD_HEIGHT_PX/2+j-lcdz*z;

	    glinter1(z,dz,
		     zmin,zmax,1e307,-1e307,
		     ih,lcdz,
		     upcolor,downcolor,diffusionz,diffusionz_limit,interval);
	  }
	  else {
	    if (found2 && !hide2nd)
	    glinter(curabc2.x,curabc2.y,curabc2.z,xscale,xc,yscale,yc,zmin2,zmax2,zmin[0],zmax[0],i,horiz,j,w,h,lcdz,downupcolor,downdowncolor,diffusionz,diffusionz_limit,interval);
	    glinter(curabc1.x,curabc1.y,curabc1.z,xscale,xc,yscale,yc,zmin,zmax,1e307,-1e307,i,horiz,j,w,h,lcdz,upcolor,downcolor,diffusionz,diffusionz_limit,interval);
	  }
	}
      }
      else {
	x = yscale*(jmax-(h-1)/2.0)-xscale*(i+(w-1)/2.0) + xc;
	y = yscale*(jmax-(h-1)/2.0)+xscale*(i+(w-1)/2.0) + yc;
	for (int j=jmax;j>=jmin;j-=h,x-=yscale*h,y-=yscale*h){
	  if (0 && i==-35 && j==-44)
	    u=0; // debug
	  // x = yscale*(j-(h-1)/2.0)-xscale*(i+(w-1)/2.0) + xc;
	  // y = yscale*(j-(h-1)/2.0)+xscale*(i+(w-1)/2.0) + yc;
	  bool found=false,found2=false;
	  if (x+y>=hyperxymin && x+y<=hyperxymax)
	    update_hypertri(hypertriangles,x,y,found,found2,curabc1,curz1,curabc2,curz2,upcolor,downcolor,downupcolor,downdowncolor);
	  for (int ki=0;ki<int(polyedrei.size());++ki){
	    int k=polyedrei[ki];
	    vector<double3> & cur=polyedrev[k];
	    if (
#if 1
		x>=polyedrexmin[ki] && x<=polyedrexmax[ki] && y>=polyedreymin[ki] && y<=polyedreymax[ki]
#else
		inside(cur,x,y)
#endif
		){
	      const double3 & abc=polyedre_abcv[k];
	      const int4 & color=polyedre_color[k];
	      // std::cout << k << " " << x << " " << y << " " << color.u << "\n";
	      double a=abc.x,b=abc.y,c=abc.z;
	      z=a*x+b*y+c;
	      bool is_clipped=polyedre_faceisclipped[k];
	      if (!is_clipped){
		double X,Y,Z;
		do_transform(invtransform3d,x,y,z,X,Y,Z);
		is_clipped=X>=window_xmin && X<=window_xmax && Y>=window_ymin && Y<=window_ymax && Z>=window_zmin && Z<=window_zmax;
	      }
	      if (is_clipped){
		update12(found,found2,
			 a,b,c,z,
			 color.u,color.d,color.du,color.dd,
			 curabc1.x,curabc1.y,curabc1.z,curz1,
			 curabc2.x,curabc2.y,curabc2.z,curz2,
			 upcolor,downcolor,downupcolor,downdowncolor);
	      }
	    } // end if inside(cur,x,y)
	  }
	  for (int k=0;k<int(sphere_centerv.size());++k){
	    if (!spheres[k]) continue;
	    const double3 & c=sphere_centerv[k];
	    double R=sphere_radiusv[k];
	    const matrice & m=*sphere_quadraticv[k]._VECTptr;
	    const vecteur & m0=*m[0]._VECTptr;
	    const vecteur & m1=*m[1]._VECTptr;
	    const vecteur & m2=*m[2]._VECTptr;
	    double v0=x-c.x,v1=y-c.y;
	    double a=m2[2]._DOUBLE_val,b=2*(m0[2]._DOUBLE_val*v0+m1[2]._DOUBLE_val*v1),C=(m0[0]._DOUBLE_val*v0+2*m0[1]._DOUBLE_val*v1)*v0+m1[1]._DOUBLE_val*v1*v1-R*R;
	    double delta=b*b-4*a*C;
	    if (delta<0)
	      continue;
	    const int4 & color=sphere_color[k];
	    delta=std::sqrt(delta);
	    double sol1,sol2;
	    if (b>0){
	      sol1=(-b-delta)/2/a;
	      sol2=2*C/(-b-delta); // (-b+delta)/2/a;
	    }
	    else {
	      sol1=2*C/(-b+delta);//(-b-delta)/2/a;
	      sol2=(-b+delta)/2/a;
	    }
	    double v2=sol1;
	    z=v2+c.z;
	    bool is_clipped=sphere_isclipped[k];
	    if (!is_clipped){
	      double X,Y,Z;
	      do_transform(invtransform3d,x,y,z,X,Y,Z);
	      is_clipped=X>=window_xmin && X<=window_xmax && Y>=window_ymin && Y<=window_ymax && Z>=window_zmin && Z<=window_zmax;
	    }
	    if (is_clipped){
	      double w0=v0*m0[0]._DOUBLE_val+v1*m1[0]._DOUBLE_val+v2*m2[0]._DOUBLE_val;
	      double w1=v0*m0[1]._DOUBLE_val+v1*m1[1]._DOUBLE_val+v2*m2[1]._DOUBLE_val;
	      double w2=v0*m0[2]._DOUBLE_val+v1*m1[2]._DOUBLE_val+v2*m2[2]._DOUBLE_val;
	      double a=-w0/w2,b=-w1/w2,c=z-(a*x+b*y);
	      update12(found,found2,
		       a,b,c,z,
		       color.u,color.d,color.du,color.dd,
		       curabc1.x,curabc1.y,curabc1.z,curz1,
		       curabc2.x,curabc2.y,curabc2.z,curz2,
		       upcolor,downcolor,downupcolor,downdowncolor);
	    }
	    if (delta<=0) continue; // delta==0, twice the same point
	    v2=sol2;
	    z=v2+c.z;
	    is_clipped=sphere_isclipped[k];
	    if (!is_clipped){
	      double X,Y,Z;
	      do_transform(invtransform3d,x,y,z,X,Y,Z);
	      is_clipped=X>=window_xmin && X<=window_xmax && Y>=window_ymin && Y<=window_ymax && Z>=window_zmin && Z<=window_zmax;
	    }
	    if (is_clipped){
	      double w0=v0*m0[0]._DOUBLE_val+v1*m1[0]._DOUBLE_val+v2*m2[0]._DOUBLE_val;
	      double w1=v0*m0[1]._DOUBLE_val+v1*m1[1]._DOUBLE_val+v2*m2[1]._DOUBLE_val;
	      double w2=v0*m0[2]._DOUBLE_val+v1*m1[2]._DOUBLE_val+v2*m2[2]._DOUBLE_val;
	      double a=-w0/w2,b=-w1/w2,c=z-(a*x+b*y);
	      update12(found,found2,
		       a,b,c,z,
		       color.u,color.d,color.du,color.dd,
		       curabc1.x,curabc1.y,curabc1.z,curz1,
		       curabc2.x,curabc2.y,curabc2.z,curz2,
		       upcolor,downcolor,downupcolor,downdowncolor);
	    }
	  } // end hypersphere loop
	  for (int k=0;k<int(plan_abcv.size());++k){
	    if (!plan_filled[k])
	      continue;
	    double3 abc=plan_abcv[k];
	    int4 color=plan_color[k];
	    // z=a*x+b*y+c
	    double z=abc.x*x+abc.y*y+abc.z,X,Y,Z;
	    do_transform(invtransform3d,x,y,z,X,Y,Z);
	    if (X>=window_xmin && X<=window_xmax && Y>=window_ymin && Y<=window_ymax && Z>=window_zmin && Z<=window_zmax)
	      update12(found,found2,
		       abc.x,abc.y,abc.z,z,
		       color.u,color.d,color.du,color.dd,
		       curabc1.x,curabc1.y,curabc1.z,curz1,
		       curabc2.x,curabc2.y,curabc2.z,curz2,
		       upcolor,downcolor,downupcolor,downdowncolor);
	  } // end hyperplan loop
	  if (found){
	    if (found2){
	      if (!hide2nd)
		glinter(curabc2.x,curabc2.y,curabc2.z,xscale,xc,yscale,yc,zmin2,zmax2,zmin[0],zmax[0],i,horiz,j,w,h,lcdz,downupcolor,downdowncolor,diffusionz,diffusionz_limit,interval);
	      glinter(curabc1.x,curabc1.y,curabc1.z,xscale,xc,yscale,yc,zmin,zmax,1e307,-1e307,i,horiz,j,w,h,lcdz,upcolor,downcolor,diffusionz,diffusionz_limit,interval);
	    }
	    else
	      glinter(curabc1.x,curabc1.y,curabc1.z,xscale,xc,yscale,yc,zmin,zmax,1e307,-1e307,i,horiz,j,w,h,lcdz,upcolor,downcolor,diffusionz,diffusionz_limit,interval);
	  }
	  else {
	    //std::cout << "not inside " << i << " " << j << " " << x << " " << y << "\n";	      
	  }
	} // end pixel vertical loop on j
      } // end else only_hypertri
    suite3d:
      // update jmintab/jmaxtab
      if (i+horiz+w<jmintabsize){
	for (int I=0;I<w;++I){
	  jmintab[i+horiz+I]=zmin[I];
	  jmaxtab[i+horiz+I]=zmax[I];
	}
      }
      // now render line/segments/curves: find intersection with plane
      // y-x=yc-xc+xscale*(2*i+I), 0<=I<w
      for (int j=0;j<int(curvev.size());++j){
	vector<double3> & cur=curvev[j];
	int s=cur.size();
	if (s<2) continue;
	int4 color=curve_color[j];
	double xy=yc-xc+xscale*2*i;
	for (int l=0;l<s-1;++l){
	  double3 m=cur[l];
	  double3 n=cur[l+1];
	  double3 v(n.x-m.x,n.y-m.y,n.z-m.z);
	  if (std::abs(v.y-v.x)<1e-4*(std::abs(v.y)+std::abs(v.x))){
	    v.y=v.x+1e-4*(std::abs(v.y)+std::abs(v.x));
	  }
	  double t1=intersect(m,v,xy);
	  if (t1<0 || t1>1)
	    continue;
	  double dt=2*xscale/(v.y-v.x); // di==1
	  double x1=m.x+t1*v.x;
	  double y1=m.y+t1*v.y;
	  double z1=m.z+t1*v.z;
	  double X1,Y1,Z1;
	  do_transform(invtransform3d,x1,y1,z1,X1,Y1,Z1);
	  if (X1<window_xmin || X1>window_xmax || Y1<window_ymin || Y1>window_ymax || Z1<window_zmin || Z1>window_zmax)
	    continue;
	  z1=LCD_HEIGHT_PX/2-lcdz*z1+(x1+y1)/2/yscale;
	  double dz=-dt*v.z*lcdz+dt*(v.x+v.y)/2/yscale;
	  int horiz=LCD_WIDTH_PX/2;
	  for (int k=0;k<w;++k){
	    double Z1=z1,Z2=z1+dz;
	    z1=Z2;
	    if (Z1>Z2) std::swap(Z1,Z2);
	    // line [ (i+k,Z1), (i+k,Z2) ]
	    if (Z2<zmin[k]){
	      drawRectangle(i+horiz+k,Z1,1,std::ceil(Z2-Z1),color.u);
	      continue;
	    }
	    if (Z1>zmax[k]){
	      drawRectangle(i+horiz+k,Z1,1,std::ceil(Z2-Z1),color.d);
	      continue;
	    }
	    drawRectangle(i+horiz+k,Z1,1,std::ceil(Z2-Z1),color.du);
	    if (Z1<zmin[k])
	      drawRectangle(i+horiz+k,Z1,1,std::ceil(zmin[k]-Z1),color.u);
	    if (Z2>zmax[k])
	      drawRectangle(i+horiz+k,zmax[k],1,std::ceil(Z2-zmax[k]),color.d);
	  }
	} // end l loop on curve discretization
      } // end loop on curves
#ifdef OLD_LINE_RENDERING
      for (int j=0;j<linetypev.size();j++){
	double3 m=linev[2*j];
	double3 v=linev[2*j+1];
	int4 color=line_color[j];
	if (std::abs(v.y-v.x)<1e-4*(std::abs(v.y)+std::abs(v.x))){
	  v.y=v.x+1e-6*(std::abs(v.y)+std::abs(v.x));
	}
	double xy=yc-xc+xscale*2*i;
	double t1=intersect(m,v,xy);
	// for halfline, t1 must be >=0, for segments between 0 and 1
	if (linetypev[j]==_HALFLINE__VECT && t1<0)
	  continue;
	if (linetypev[j]==_GROUP__VECT && (t1<0 || t1>1))
	  continue;
	double dt=2*xscale/(v.y-v.x); // di==1
	double x1=m.x+t1*v.x;
	double y1=m.y+t1*v.y;
	double z1=m.z+t1*v.z;
	double X1,Y1,Z1;
	do_transform(invtransform3d,x1,y1,z1,X1,Y1,Z1);
	// int dbgi,dbgj; xyz2ij(double3(x1,y1,z1),dbgi,dbgj);
	/// double x2=x1+dt*v.x,y2=y1+dt*v.y,z2=z1+dt*v.z;
	if (X1<window_xmin || X1>window_xmax || Y1<window_ymin || Y1>window_ymax || Z1<window_zmin || Z1>window_zmax)
	  continue;
	z1=LCD_HEIGHT_PX/2-lcdz*z1+(x1+y1)/2/yscale;
	double dz=-dt*v.z*lcdz+dt*(v.x+v.y)/2/yscale;
	int horiz=LCD_WIDTH_PX/2;
	for (int k=0;k<w;++k){
	  double Z1=z1,Z2=z1+dz;
	  z1=Z2;
	  if (Z1>Z2) std::swap(Z1,Z2);
	  // line [ (i+k,Z1), (i+k,Z2) ]
	  if (Z2<zmin[k]){
	    drawRectangle(i+horiz+k,Z1,1,std::ceil(Z2-Z1),color.u);
	    continue;
	  }
	  if (Z1>zmax[k]){
	    drawRectangle(i+horiz+k,Z1,1,std::ceil(Z2-Z1),color.d);
	    continue;
	  }
	  drawRectangle(i+horiz+k,Z1,1,std::ceil(Z2-Z1),color.du);
	  if (Z1<zmin[k])
	    drawRectangle(i+horiz+k,Z1,1,std::ceil(zmin[k]-Z1),color.u);
	  if (Z2>zmax[k])
	    drawRectangle(i+horiz+k,zmax[k],1,std::ceil(Z2-zmax[k]),color.d);	    
	}
      } // end lines rendering
#endif
      // points rendering
      for (int j=0;j<int(pointv.size());++j){
	const double3 & m=pointv[j];
	if (m.x<i+horiz || m.x>=i+horiz+w)
	  continue;
	const int4 & c=point_color[j];
	int k=m.x-i-horiz,color=-1;
	double mz=LCD_HEIGHT_PX/2-lcdz*m.z;
	double dz=(zmax[k]-zmin[k])*1e-3;
	if (mz>=zmax[k]-dz)
	  color=c.u;
	else if (mz<=zmin[k]+dz)
	  color=c.u; // c.d?
	else color=c.du;
	drawRectangle(m.x,m.y,3,3,color);
	if (points[j]){
	  // int dx=RAND_MAX+os_draw_string(-RAND_MAX,0,color,0,points[j],false); // fake print
	  int dx=os_draw_string_small(0,0,color,0,points[j],true); // fake print
	  os_draw_string_small(m.x-dx,m.y,color,0,points[j],false); 
	}
      } // end points rendering
    } // end pixel horizontal loop on i
#ifndef OLD_LINE_RENDERING
    // new line rendering
    for (int j=0;j<int(linetypev.size());j++){
      double mx,my,mz,vx,vy,vz;
      double3 m=linev[2*j];
      do_transform(invtransform3d,m.x,m.y,m.z,mx,my,mz);
      double3 v=linev[2*j+1];
      do_transform(invtransform3d,m.x+v.x,m.y+v.y,m.z+v.z,vx,vy,vz);
      vx -= mx; vy -= my; vz -= mz;
      int4 color=line_color[j];
      // find min/max parameter value for hidden/visible
      double tmax=1e306,tmin=-1e306;
      if (linetypev[j]==_HALFLINE__VECT){
	tmin=0;
      }
      if (linetypev[j]==_GROUP__VECT){
	tmin=0;
	tmax=1;
      }
      if (vx!=0){ // intersect with x=window_xmin and x=window_xmax
	// m.x+v.x*t=window_xmin/max
	double t1=(window_xmin-mx)/vx;
	double t2=(window_xmax-mx)/vx;
	if (t1>t2)
	  std::swap<double>(t1,t2);
	if (t1>tmin)
	  tmin=t1;
	if (t2<tmax)
	  tmax=t2;
      }
      if (vy!=0){ // intersect with y=window_ymin and y=window_ymax
	double t1=(window_ymin-my)/vy;
	double t2=(window_ymax-my)/vy;
	if (t1>t2)
	  std::swap<double>(t1,t2);
	if (t1>tmin)
	  tmin=t1;
	if (t2<tmax)
	  tmax=t2;
      }
      if (vz!=0){ // intersect with z=window_zmin and z=window_zmax
	double t1=(window_zmin-mz)/vz;
	double t2=(window_zmax-mz)/vz;
	if (t1>t2)
	  std::swap<double>(t1,t2);
	if (t1>tmin)
	  tmin=t1;
	if (t2<tmax)
	  tmax=t2;
      }
      // find which intersection we want
      bool usetmax=v.x+v.y==0?v.z>=0:v.x+v.y>0; 
      double tmin_=usetmax?tmin:tmax,
	tmax_=usetmax?tmin:tmax;
      for (int k=0;k<int(plan_abcv.size());++k){
	if (!plan_filled[k]) continue;
	// z >= z_plan=a*x+b*y+c where (x,y,z)=m+t*v
	// m.z+t*v.z >= a*m.x+t*a*v.x+b*m.y+t*b*v.y+c
	// t*(v.z-a*v.x-b.v.y) >= a*m.x+b*m.y+c-m.z
	double3 abc=plan_abcv[k];
	double A=v.z-abc.x*v.x-abc.y*v.y,B=abc.x*m.x+abc.y*m.y+abc.z-m.z;
	if (A==0) continue;
	double t=B/A;
	if (t<=tmin || t>=tmax)
	  continue;
	if (tmax_<t)
	  tmax_=t;
	if (tmin_>t)
	  tmin_=t;
      }
      for (int k=0;k<int(sphere_centerv.size());++k){
	// sphere interesect line 
	double3 c=sphere_centerv[k];
	double R=sphere_radiusv[k];
	int4 color=sphere_color[k];
	double cx,cy,cz;
	do_transform(invtransform3d,c.x,c.y,c.z,cx,cy,cz);
	// (mx+t*vx-cx)^2+(my+t*vy-cy)^2+(mz+t*vz-cz)^2=R^2
	double a=vx*vx+vy*vy+vz*vz,b=(vx*(mx-cx)+vy*(my-cy)+vz*(mz-cz)),C=(mx-cx)*(mx-cx)+(my-cy)*(my-cy)+(mz-cz)*(mz-cz)-R*R;
	double deltaprime=b*b-a*C;
	if (deltaprime>0){
	  deltaprime=std::sqrt(deltaprime);
	  double t1=(-b-deltaprime)/a,t2=(-b+deltaprime)/a;
	  if (t1>t2)
	    std::swap(t1,t2);
	  if (t1<tmin || t1>tmax)
	    t1=t2;
	  if (t2<tmin || t2>tmax)
	    t2=t1;
	  double t=usetmax?t2:t1;
	  if (t<=tmin || t>=tmax)
	    continue;
	  if (tmin_>t)
	    tmin_=t;
	  if (tmax_<t)
	    tmax_=t;	  
	}
      }
      vector<double> interpoly;
      for (int k=0;k<int(polyedre_abcv.size());++k){
	if (!polyedre_filled[k]) continue;
	double3 abc=polyedre_abcv[k];
	double a=abc.x,b=abc.y,c=abc.z;
	// intersect z=a*x+b*y+c with line m+t*v
	// m.z+t*v.z=a*(m.x+t*v.x)+b*(m.y+t*v.y)+c
	double A=v.z-a*v.x-b*v.y,B=a*m.x+b*m.y+c-m.z;
	if (A!=0){
	  double t=B/A;
	  double x=m.x+t*v.x;
	  double y=m.y+t*v.y;
	  if (t>tmin && t<tmax && is_inside(polyedrev[k],x,y)){
	    interpoly.push_back(t);
	  }
	}
      }
      if (interpoly.size()>=1){
	sort(interpoly.begin(),interpoly.end());
	double t=usetmax?interpoly.back():interpoly.front();
	if (tmin_>t)
	  tmin_=t;
	if (tmax_<t)
	  tmax_=t;
      }
      vecteur sv(native_renderer_instructions);
      for (int k=0;k<int(sv.size());++k){
	gen surf=remove_at_pnt(sv[k]);
	if (surf.is_symb_of_sommet(at_hypersurface)){
	  const vecteur & hyp=*surf._SYMBptr->feuille._VECTptr;
	  if (hyp.size()>2 && !is_undef(hyp[1])){
	    gen eq=hyp[1],vars=hyp[2];
	    if (_is_polynomial(makesequence(eq,vars[0]),contextptr)==1 && _is_polynomial(makesequence(eq,vars[1]),contextptr)==1 && _is_polynomial(makesequence(eq,vars[2]),contextptr)==1){
	      vecteur V(makevecteur(mx+vx*t__IDNT_e,my+vy*t__IDNT_e,mz+vz*t__IDNT_e));
	      gen eqt=subst(eq,vars,V,false,contextptr);
	      gen gradeq=subst(derive(eq,vars,contextptr),vars,V,false,contextptr); // not used
	      vecteur tval_=gen2vecteur(_sort(solve(eqt,t__IDNT_e,0,contextptr),contextptr));
	      int s=tval_.size();
	      vector<double> tval;
	      for (int l=0;l<s;++l){
		gen tg=evalf_double(tval_[l],1,contextptr);
		if (tg.type==_DOUBLE_){
		  double t=tg._DOUBLE_val;
		  if (t>=tmin && t<=tmax){
		    if (tval.size()>1)
		      tval.back()=t;
		    else
		      tval.push_back(t);
		  }
		}
	      }
#if 1
	      if (!tval.empty()){
		double t=usetmax?tval.back():tval.front();
		if (tmin_>t)
		  tmin_=t;
		if (tmax_<t)
		  tmax_=t;
	      }
#else
	      if (tval.size()==1){
		// grad of eq dot [vx,vy.vz] is positive if entering surface
		gen gradt=subst(gradeq,t__IDNT_e,tval[0],false,contextptr);
		gen dotp=dotvecteur(gen2vecteur(gradt),makevecteur(vx,vy,vz));
		if (is_positive(dotp,contextptr)){
		  if (tmin_>tval[0])
		    tmin_=tval[0];
		  tmax_=tmax;
		}
		else {
		  if (tmax_<tval[0])
		    tmax_=tval[0];
		  tmin_=tmin;
		}
	      }
	      else if (tval.size()==2){
		if (tmin_>tval[0])
		  tmin_=tval[0];
		if (tmax_<tval[1])
		  tmax_=tval[1];
	      }
#endif
	    } // end polynomial hypersurface
	  }
	} // end hypersurface	
      }
      int i1,j1,i2,j2;
      grmtv2ij(*this,m,v,tmin,i1,j1);
      grmtv2ij(*this,m,v,tmax,i2,j2);
      if (tmin==tmin_ && tmax==tmax_){
	int curcolor=color.u;
	if (0 && i1<jmintabsize){
	  if (jmaxtab[i1]>=jmintab[i1] && j1>=jmaxtab[i1] )
	    curcolor=color.d;
	}
	drawLine(i1,j1,i2,j2,curcolor);
      }
      else {
	int curcolor=color.d;
	if (0 && i1<jmintabsize){
	  if (jmaxtab[i1]>=jmintab[i1] && j1<jmintab[i1])
	    curcolor=color.u;
	}
	drawLine(i1,j1,i2,j2,curcolor | 0x400000);
      }
      bool name=show_names && lines[j];
      if (tmin<tmin_){
	grmtv2ij(*this,m,v,tmin,i1,j1);
	grmtv2ij(*this,m,v,tmin_,i2,j2);
	int curcolor=color.u;
	if (0 && i1<jmintabsize){
	  if (jmaxtab[i1]>=jmintab[i1] && j1>=jmaxtab[i1] )
	    curcolor=color.d;
	}
	drawLine(i1,j1,i2,j2,curcolor);
	if (name){
	  os_draw_string(i1,j1,color.u,0,lines[j]);
	  name=false;
	}
	if (tmin_!=tmax) os_fill_rect(i2-2,j2-2,4,4,curcolor);
      }
      if (tmax>tmax_){
	grmtv2ij(*this,m,v,tmax_,i1,j1);
	grmtv2ij(*this,m,v,tmax,i2,j2);
	int curcolor=color.u;
	if (0 && i1<jmintabsize){
	  if (jmaxtab[i1]>=jmintab[i1] && j1>=jmaxtab[i1] )
	    curcolor=color.d;
	}
	drawLine(i1,j1,i2,j2,curcolor);
	if (name){
	  os_draw_string(i2,j2,color.u,0,lines[j]);
	  name=false;
	}
	if (tmax_!=tmin) os_fill_rect(i1-2,j1-2,4,4,curcolor);
      }
    }
#endif // OLD_LINE_RENDERING
    sync_screen();
    return true;
  }
  
  void vpush(vecteur & w,const gen & g){
    vecteur v(gen2vecteur(g));
    w.reserve(w.size()+v.size());
    for (int i=0;i<v.size();++i)
      w.push_back(v[i]);
  }

  void Graph3d::draw3d(GIAC_CONTEXT){
    native_renderer_instructions.clear();
    if ( (display_mode & 2) && !animation_instructions.empty()){
      vpush(native_renderer_instructions,animation_instructions[animation_instructions_pos % animation_instructions.size()]);
    }
    if ( display_mode & 0x40 ){
      vpush(native_renderer_instructions,trace_instructions);
    }
    if (display_mode & 1)
      vpush(native_renderer_instructions,plot_instructions);
    gen plot_tmp,title_tmp;
    find_title_plot(title_tmp,plot_tmp,contextptr);
    vpush(native_renderer_instructions,plot_tmp);
    update_rotation();
    LCD_WIDTH_PX=w();
    LCD_HEIGHT_PX=h();
    double3 A(window_xmin,window_ymin,window_zmin),
      B(window_xmin,window_ymin,window_zmax),
      C(window_xmax,window_ymin,window_zmin),
      D(window_xmax,window_ymin,window_zmax),
      E(window_xmin,window_ymax,window_zmin),
      F(window_xmin,window_ymax,window_zmax),
      G(window_xmax,window_ymax,window_zmin),
      H(window_xmax,window_ymax,window_zmax);
    double3 A3,B3,C3,D3,E3,F3,G3,H3;
    xyz2ij(A,Ai,Aj,A3);
    xyz2ij(B,Bi,Bj,B3);
    xyz2ij(C,Ci,Cj,C3);
    xyz2ij(D,Di,Dj,D3);
    xyz2ij(E,Ei,Ej,E3);
    xyz2ij(F,Fi,Fj,F3);
    xyz2ij(G,Gi,Gj,G3);
    xyz2ij(H,Hi,Hj,H3);
    indraw3d(precision,precision,lcdz,contextptr,default_upcolor,default_downcolor,default_downupcolor,default_downdowncolor);
    if (hp)
      draw_decorations(title_tmp);
    if (show_edges){
      // polyhedrons
      for (int k=0;k<int(polyedrev.size());++k){
	const vector<double3> & cur=polyedrev[k]; // current face
	const int4 & col=polyedre_color[k];
	for (int l=1;l<int(cur.size());++l){
	  const double3 & p=cur[l?l-1:cur.size()-1];
	  const double3 & c=cur[l];
	  // is edge visible?
	  double3 m(p.x/2+c.x/2,p.y/2+c.y/2,p.z/2+c.z/2);
	  double xy=m.x+m.y;
	  int mi,mj;
	  XYZ2ij(m,mi,mj);
	  int kk,jmin=RAND_MAX,jmax=-RAND_MAX;
	  for (kk=0;kk<int(polyedrev.size());++kk){
	    if (k==kk)
	      continue;
	    const vector<double3> & Cur=polyedrev[kk];
	    int ll;
	    // first check if point is in face
	    for (ll=1;ll<int(Cur.size());++ll){
	      const double3 & P=Cur[ll?ll-1:Cur.size()-1];
	      const double3 & C=Cur[ll];
	      double3 M(P.x/2+C.x/2,P.y/2+C.y/2,P.z/2+C.z/2);
	      if (M.x==m.x && M.y==m.y && M.z==m.z){
		break; // edge PC has same midpoint, will ignore face
	      }
	    }
	    if (ll<int(Cur.size())) // point is in face, ignore face
	      continue;
	    double3 M0; bool found1st=false;
	    for (ll=1;ll<int(Cur.size());++ll){
	      const double3 & P=Cur[ll?ll-1:Cur.size()-1];
	      const double3 & C=Cur[ll];
	      // intersect plane y-x=m.y-m.x with PC edge P+t*PC
	      double PCx=C.x-P.x,PCy=C.y-P.y,dPC=PCy-PCx;
	      // P.y-P.x + t*dPC=m.y-m.x
	      if (dPC==0) // edge is parallel
		continue;
	      double t=((m.y-m.x)+(P.x-P.y))/dPC;
	      if (t<0 || t>1)
		continue;
	      double x=P.x+t*PCx;
	      double y=P.y+t*PCy;
	      double z=P.z+t*(C.z-P.z);
	      if (!found1st){
		M0=double3(x,y,z);
		found1st=true;
		continue;
	      }
	      if (x==M0.x && y==M0.y && z==M0.z)
		continue;
	      // segment([x,y,z],M0) has same y-x as m,
	      // find segment position for same y+x as m [x,y,z]+t*(M0-[x,y,z])
	      // yx=x+y+t*(M0.x-x+M0.y-y)
	      double M0xy=M0.x-x+M0.y-y;
	      int i1,j1,i2,j2; // N.B. i1,i2 should be the same as mi
	      if (std::abs(M0xy)<1e-14)
		t=-1;
	      else
		t=(xy-x-y)/M0xy;
	      if (t<=0 || t>=1){
		if (x+y<=xy) // segment is behind midpoint m
		  continue;
		XYZ2ij(M0,i1,j1); 
		XYZ2ij(double3(x,y,z),i2,j2);
		if (j1>j2) swapint(j1,j2);
		if (jmin>j1) jmin=j1;
		if (jmax<j2) jmax=j2;
		if (jmin<mj && mj<jmax){
		  break;
		}
		continue;
	      }
	      // find segment part that might mask midpoint m
	      double X = x+t*(M0.x-x);
	      double Y = y+t*(M0.y-y);
	      double Z = z+t*(M0.z-z);
	      XYZ2ij(double3(X,Y,Z),i1,j1);
	      if (x+y<=xy)
		XYZ2ij(M0,i2,j2);
	      else
		XYZ2ij(double3(x,y,z),i2,j2);
	      if (j1>j2) swapint(j1,j2);
	      if (jmin>j1) jmin=j1;
	      if (jmax<j2) jmax=j2;
	      if (jmin<mj && mj<jmax){
		break;
	      }
	    } // end for
	    if (ll<int(Cur.size())){
	      // means edge is not visible
	      break;
	    }
	  }
	  // polyedre attribute: filled/not filled
	  bool filled=polyedre_filled[k];
	  bool hidden=kk<int(polyedrev.size());
	  if (filled && hidden)
	    continue;
	  int i1,j1,i2,j2;
	  XYZ2ij(p,i1,j1);
	  XYZ2ij(c,i2,j2);
	  if (i1>i2 || (i1==i2 && j1>j2)){
	    swapint(i1,i2); swapint(j1,j2);
	  }
	  drawLine(i1,j1,i2,j2,
		   // col.d | 0x400000
		   col.u | ((hidden || filled)?0x400000:0x20000)
		   );
	}
      }
    }
    if (show_axes){
      // cube A,B,C,D,E,F,G,H
      // X
      drawLine(Ai,Aj,Ci,Cj,COLOR_RED | 0x800000);
      drawLine(Bi,Bj,Di,Dj,COLOR_RED | 0x800000);
      drawLine(Ei,Ej,Gi,Gj,COLOR_RED | 0x800000);
      drawLine(Fi,Fj,Hi,Hj,COLOR_RED | 0x800000);
      // Y
      drawLine(Ai,Aj,Ei,Ej,COLOR_GREEN | 0x800000);
      drawLine(Bi,Bj,Fi,Fj,COLOR_GREEN | 0x800000);
      drawLine(Ci,Cj,Gi,Gj,COLOR_GREEN | 0x800000);
      drawLine(Di,Dj,Hi,Hj,COLOR_GREEN | 0x800000);
      // Z
      drawLine(Ai,Aj,Bi,Bj,COLOR_BLUE | 0x800000);
      drawLine(Ci,Cj,Di,Dj,COLOR_BLUE | 0x800000);
      drawLine(Ei,Ej,Fi,Fj,COLOR_BLUE | 0x800000);
      drawLine(Gi,Gj,Hi,Hj,COLOR_BLUE | 0x800000);
      // planes
      // current_depth
      if (dynamic_cast<Geo3d *>(this)){
	vector<int2> polyg; int2 IJmin(RAND_MAX,RAND_MAX);
	// x: A3-C3, B3-D3; E3-G3,F3-H3
	adddepth(polyg,A3,C3,IJmin);
	adddepth(polyg,B3,D3,IJmin);
	adddepth(polyg,E3,G3,IJmin);
	adddepth(polyg,F3,H3,IJmin);
	// y: A3-E3; B3-F3; C3-G3, D3-H3
	adddepth(polyg,A3,E3,IJmin);
	adddepth(polyg,B3,F3,IJmin);
	adddepth(polyg,C3,G3,IJmin);
	adddepth(polyg,D3,H3,IJmin);
	// z: A3-B3, C3-D3, E3-F3, G3-H3
	adddepth(polyg,A3,B3,IJmin);
	adddepth(polyg,C3,D3,IJmin);
	adddepth(polyg,E3,F3,IJmin);
	adddepth(polyg,G3,H3,IJmin);
	int Px,Py;
	displaypolyg(polyg,IJmin,
		     // FL_CYAN,
		     FL_YELLOW | 0x400000,
		     Px,Py);
      }
      vecteur attrv(native_renderer_instructions);
      for (int i=0;i<attrv.size();++i){
	gen attr=attrv[i];
	gen cur=remove_at_pnt(attr);
	int upcolor=44444;
	const char * nameptr=0;
	if (attr.is_symb_of_sommet(at_pnt)){
	  if (show_names && attr._SYMBptr->feuille.type==_VECT && attr._SYMBptr->feuille._VECTptr->size()==3){
	    gen name=attr._SYMBptr->feuille._VECTptr->back();
	    if (name.type==_IDNT)
	      nameptr=name._IDNTptr->id_name;
	    if (name.type==_STRNG)
	      nameptr=name._STRNGptr->c_str();
	  }
	  attr=attr._SYMBptr->feuille[1];
	  if (attr.type==_INT_ && (attr.val & 0xffff)!=56){
	    upcolor=attr.val &0xffff;
	  }
	}
	if (cur.is_symb_of_sommet(at_hyperplan)){
	  vecteur & w=*cur._SYMBptr->feuille._VECTptr;
	  gen m=evalf_double(w[1],1,contextptr),n=evalf_double(w[0],1,contextptr);
	  double a=n[0]._DOUBLE_val,b=n[1]._DOUBLE_val,c=n[2]._DOUBLE_val;
	  double x0=m[0]._DOUBLE_val,y0=m[1]._DOUBLE_val,z0=m[2]._DOUBLE_val;
	  // a*(x-x0)+b*(y-y0)+c*(z-z0)=0
	  // replace 2 coordinates of M with window_xyzminmax and find last coord
	  vector<int2> polyg; int2 IJmin(RAND_MAX,RAND_MAX);
	  // x
	  if (a!=0){
	    double x=x0-1/a*(b*(window_ymin-y0)+c*(window_zmin-z0));
	    if (x>=window_xmin && x<=window_xmax)
	      addpolyg(polyg,x,window_ymin,window_zmin,IJmin);
	    x=x0-1/a*(b*(window_ymin-y0)+c*(window_zmax-z0));
	    if (x>=window_xmin && x<=window_xmax)
	      addpolyg(polyg,x,window_ymin,window_zmax,IJmin);
	    x=x0-1/a*(b*(window_ymax-y0)+c*(window_zmin-z0));
	    if (x>=window_xmin && x<=window_xmax)
	      addpolyg(polyg,x,window_ymax,window_zmin,IJmin);
	    x=x0-1/a*(b*(window_ymax-y0)+c*(window_zmax-z0));
	    if (x>=window_xmin && x<=window_xmax)
	      addpolyg(polyg,x,window_ymax,window_zmax,IJmin);
	  }
	  // y
	  if (b!=0){
	    double y=y0-1/b*(a*(window_xmin-x0)+c*(window_zmin-z0));
	    if (y>=window_ymin && y<=window_ymax)
	      addpolyg(polyg,window_xmin,y,window_zmin,IJmin);
	    y=y0-1/b*(a*(window_xmin-x0)+c*(window_zmax-z0));
	    if (y>=window_ymin && y<=window_ymax)
	      addpolyg(polyg,window_xmin,y,window_zmax,IJmin);
	    y=y0-1/b*(a*(window_xmax-x0)+c*(window_zmin-z0));
	    if (y>=window_ymin && y<=window_ymax)
	      addpolyg(polyg,window_xmax,y,window_zmin,IJmin);
	    y=y0-1/b*(a*(window_xmax-x0)+c*(window_zmax-z0));
	    if (y>=window_ymin && y<=window_ymax)
	      addpolyg(polyg,window_xmax,y,window_zmax,IJmin);
	  }
	  // z
	  if (c!=0){
	    double z=z0-1/c*(a*(window_xmin-x0)+b*(window_ymin-y0));
	    if (z>=window_zmin && z<=window_zmax)
	      addpolyg(polyg,window_xmin,window_ymin,z,IJmin);
	    z=z0-1/c*(a*(window_xmin-x0)+b*(window_ymax-y0));
	    if (z>=window_zmin && z<=window_zmax)
	      addpolyg(polyg,window_xmin,window_ymax,z,IJmin);
	    z=z0-1/c*(a*(window_xmax-x0)+b*(window_ymin-y0));
	    if (z>=window_zmin && z<=window_zmax)
	      addpolyg(polyg,window_xmax,window_ymin,z,IJmin);
	    z=z0-1/c*(a*(window_xmax-x0)+b*(window_ymax-y0));
	    if (z>=window_zmin && z<=window_zmax)
	      addpolyg(polyg,window_xmax,window_ymax,z,IJmin);
	  }
	  // sort list of arguments
	  vector<int2_double2> p;
	  for (int k=0;k<polyg.size();++k){
	    int2 & cur=polyg[k];
	    if (cur==IJmin){
	      int2_double2 id={cur.i,cur.j,0,0};
	      p.push_back(id);
	    } else {
	      double di=cur.i-IJmin.i,dj=cur.j-IJmin.j;
	      int2_double2 id={cur.i,cur.j,atan2(di,dj),di*di+dj*dj};
	      p.push_back(id);
	    }
	  }
	  sort(p.begin(),p.end());
	  // draw polygon
	  vector< vector<int> > P;
	  int ps=p.size();
	  drawLine(p[ps-1].i,p[ps-1].j,p[0].i,p[0].j,upcolor);	  
	  for (int k=1;k<ps;++k){
	    drawLine(p[k-1].i,p[k-1].j,p[k].i,p[k].j,upcolor);
	  }
	  if (nameptr){
	    int x=os_draw_string_small(0,0,0,upcolor,nameptr,true);
	    os_draw_string_small(p[0].i-x,p[0].j,upcolor,0,nameptr);
	  }
	}
      }
      // frame
      double xi=Ci-Ai,xj=Cj-Aj;
      normalize(xi,xj);
      drawLine(20,20,20+20*xi,20+20*xj,COLOR_RED);
      os_draw_string_small(20+20*xi,20+20*xj,COLOR_RED,COLOR_BLACK,"x");
      double yi=Ei-Ai,yj=Ej-Aj;
      normalize(yi,yj);
      drawLine(20,20,20+20*yi,20+20*yj,COLOR_GREEN);
      os_draw_string_small(20+20*yi,20+20*yj,COLOR_GREEN,COLOR_BLACK,"y");
      double zi=Bi-Ai,zj=Bj-Aj;
      normalize(zi,zj);
      drawLine(20,20,20+20*zi,20+20*zj,COLOR_BLUE);
      os_draw_string_small(20+20*zi,20+20*zj,COLOR_BLUE,COLOR_BLACK,"z");
    } // end show_axes
      // now handle legend([x,y],string)
    const vecteur & V=native_renderer_instructions;
    for (int i=0;i<V.size();++i){
      gen attr=V[i];
      if (attr.is_symb_of_sommet(at_pnt)){
	attr=attr._SYMBptr->feuille;
	if (attr.type==_VECT && attr._VECTptr->size()>1){
	  int color=65535;
	  gen attr0=attr._VECTptr->front();
	  attr=attr[1];
	  if (attr.type==_INT_ && (attr.val & 0xffff)!=0){
	    color=attr.val &0xffff;
	  }
	  if (attr0.is_symb_of_sommet(at_legende)){
	    gen leg=attr0._SYMBptr->feuille;
	    if (leg.type==_VECT && leg._VECTptr->size()>=2){
	      gen pos=leg._VECTptr->front();
	      leg=leg[1];
	      if (pos.type==_VECT && pos._VECTptr->size()==2 && leg.type==_STRNG){
		gen x=pos._VECTptr->front(),y=pos._VECTptr->back();
		if (x.type==_INT_ && y.type==_INT_)
		  os_draw_string(x.val,y.val,color,0,leg._STRNGptr->c_str());
	      }
	    }
	  }
	}
      }
    }
    if (mode==1 && pushed && push_in_area && in_area){
      // draw segment between push and current
      double x,y,z;
      find_xyz(push_i,push_j,push_depth,x,y,z);
#if 0 // FIXME 3d geometry
      gl_color(FL_RED);
      glBegin(GL_LINES);
      glVertex3d(x,y,z);
      find_xyz(current_i,current_j,current_depth,x,y,z);
      glVertex3d(x,y,z);
      glEnd();
#endif
    }
    if (displayhint)
      os_draw_string_small(0,h()-16,COLOR_WHITE,COLOR_BLACK,gettext("click o: OpenGL"));
  }  

  // if printing is true, we call gl2ps to make an eps file
  void Graph3d::draw(){
    // cerr << "graph3d" << '\n';
    int clip_x,clip_y,clip_w,clip_h;
#ifdef __APPLE__
    if (!find_clip_box(this,clip_x,clip_y,clip_w,clip_h))
#endif
      fl_clip_box(x(),y(),w(),h(),clip_x,clip_y,clip_w,clip_h);
    // cerr << clip_x << " " << clip_y<< " " << clip_w<< " " << clip_h << '\n';
    if (opengl && printing){
      fprintf(printing,"%s","\nGS\n\nCR\nCS\n% avant eps -> BeginDocument\n/b4_Inc_state save def\n%save state for cleanup\n/dict_count countdictstack def\n/op_count count 1 sub def\n%count objects on op stack\nuserdict begin\n%make userdict current dict\n/showpage { } def\n%redefine showpage to be null\n0 setgray 0 setlinecap\n1 setlinewidth 0 setlinejoin\n10 setmiterlimit [] 0 setdash newpath\n/languagelevel where\n%if not equal to 1 then\n{pop languagelevel                %set strokeadjust and\n1 ne\n%overprint to their defaults\n{false setstrokeadjust false setoverprint\n} if\n} if\n");
      // Translate by previous widget h()
      fprintf(printing,"[1 0 0 -1 0 0] concat\n%d %d translate\n",x(),-h()-y());
      fprintf(printing,"%s","\n%%BeginDocument: out.eps\n");
      print();
      fprintf(printing,"%s","\n%%EndDocument\ncount op_count sub {pop} repeat\ncountdictstack dict_count sub {end} repeat  %clean up dict stack\nb4_Inc_state restore\n% apres eps\n\nGR\n");
      return;
    }
    Fl_Window * win = window();
    if (!win)
      return;
    if (!hp)
      hp=geo_find_history_pack(this);
    context * contextptr = hp?hp->contextptr:0;
    int locked=0;
#ifdef HAVE_LIBPTHREAD
    locked=pthread_mutex_trylock(&interactive_mutex);
#endif
    bool b,block;
    if (!locked){
      b=io_graph(contextptr);
      io_graph(contextptr)=false;
      block=block_signal;
      block_signal=true;
    }
    if (opengl){
#if defined __APPLE__ && !defined GRAPH_WINDOW
      GLContext context;
      if (!glcontext){ // create context
	GLContext shared_ctx = 0;
	context = fl_create_gl_context(Fl_Window::current(), gl_choice);//aglCreateContext( gl_choice->pixelformat, shared_ctx);
	if (!context){ 
	  if (!locked){
	    block_signal=block;
	    io_graph(contextptr)=b;
#ifdef HAVE_LIBPTHREAD
	    pthread_mutex_unlock(&interactive_mutex);
#endif
	  }
	  return ;
	}
	glcontext = (void *) context;
      }
      else
	context = (GLContext) glcontext;
      fl_set_gl_context(Fl_Window::current(), context);
      
      //aglSetCurrentContext(context);
      //GLint rect[] = { clip_x, win->h()-clip_y-clip_h, clip_w, clip_h};
      //aglSetInteger( (GLContext)context, AGL_BUFFER_RECT, rect );
      //aglEnable( (GLContext)context, AGL_BUFFER_RECT );
      //OpaqueWindowPtr * winid=(OpaqueWindowPtr*) win->window_ref();
      //aglSetWindowRef( context, winid );
      glEnable(GL_DEPTH_TEST);
      glLoadIdentity();
      glEnable(GL_SCISSOR_TEST);
      glScissor(clip_x, win->h()-clip_y-clip_h, clip_w, clip_h); // lower left
      glViewport(clip_x, win->h()-clip_y-clip_h, clip_w, clip_h); // lower left
      //glViewport(0,0,w(),h());
      // glDrawBuffer(GL_FRONT);
      gl_font(FL_HELVETICA,labelsize());
      // GLint viewport[4];
      // glGetIntergerv(GL_VIEWPORT,viewport);
      //fl_push_clip(clip_x,clip_y,clip_w,clip_h);
      display();
#else
      gl_start();
      glEnable(GL_SCISSOR_TEST);
      // glViewport(clip_x, win->h()-clip_y-clip_h, clip_w, clip_h); // lower left
#ifdef GRAPH_WINDOW
      glViewport(0,0, w(), h()); // lower left
#else
      glScissor(clip_x, win->h()-clip_y-clip_h, clip_w, clip_h); // lower left
      glViewport(x(),win->h()-y()-h(), w(), h()); // lower left
#endif
      gl_font(FL_HELVETICA,labelsize());
      // GLint viewport[4];
      // glGetIntergerv(GL_VIEWPORT,viewport);
      //fl_push_clip(clip_x,clip_y,clip_w,clip_h);
      display();
      gl_finish();
#endif
    }
    else {
      fl_push_clip(clip_x,clip_y,clip_w,clip_h);
      draw3d(contextptr);
      fl_pop_clip();
    }
    struct timezone tz;
    gettimeofday(&animation_last,&tz);
    if (!paused)
      ++animation_instructions_pos;
    //fl_pop_clip();
    if (!locked){
      block_signal=block;
      io_graph(contextptr)=b;
#ifdef HAVE_LIBPTHREAD
      pthread_mutex_unlock(&interactive_mutex);
#endif
    }
  }

  const char * Graph3d::latex(const char * filename_){
    const char * filename=0;
    if ( (filename=Graph2d3d::latex(filename_)) ){
      ofstream of(filename);
      if (of){
	string epsfile(remove_extension(filename)+".eps");
	FILE * f=fopen(epsfile.c_str(),"w");
	if (f){
	  of << "% Generated by xcas\n" << '\n';
	  of << giac::tex_preamble << '\n';
	  of << "\\includegraphics[bb=0 0 400 "<< h() <<"]{" << remove_path(epsfile) << "}" << '\n' << '\n' ;
	  of << giac::tex_end << '\n';
	  of.close();
	  int gw=w();
	  resize(x(),y(),400,h());
	  printing=f;
	  print();
	  printing=0;
	  resize(x(),y(),gw,h());
	  fclose(f);	
	}
      }
      
      /*
      FILE * f=fopen(filename,"w");
      if (f){
	int gw=w();
	resize(x(),y(),400,h());
	fprintf(f,"%s","% Generated by xcas\n");
	fprintf(f,"%s",giac::tex_preamble.c_str()) ;
	fprintf(f,"\n\n{\n\\setlength{\\unitlength}{1pt}\n\\begin{picture}(%d,%d)(%d,%d)\n{\\special{\"\n",w(),h(),0,0);
	printing=f;
	print();
	printing=0;
	resize(x(),y(),gw,h());
	fprintf(f,"%s","}\n}\\end{picture}\n}\n\n");
	fprintf(f,"%s",giac::tex_end.c_str());
	fclose(f);	
      }
      */
    }
    return filename;
  }

  Graph3d::~Graph3d(){ 
#if defined __APPLE__ && !defined GRAPH_WINDOW
    if (glcontext){
      fl_delete_gl_context((GLContext) glcontext);
      //aglSetCurrentContext( NULL );
      //aglSetWindowRef((GLContext) glcontext, NULL );    
      //aglDestroyContext((GLContext)glcontext);
    }
#endif
  }

  void Graph3d::reset_view() {
    if (opengl){
      theta_x=-13;
      theta_y=-95;
      theta_z=-110; 
    }
    else 
      theta_x=theta_y=theta_z=0;
    q=euler_deg_to_quaternion_double(theta_z,theta_x,theta_y);
  }

  void Graph3d::switch_renderer() {
#if 0 // def WIN32
    fl_alert(gettext("Not supported under Windows. Sticking to OpenGL."));
    return;
#endif
    opengl=!opengl;
    reset_view();
  }

  int Graph3d::in_handle(int event){
    int res=Graph2d3d::in_handle(event);
    if (event==FL_FOCUS){
      if (!paused)
	--animation_instructions_pos;
      redraw();
    }
    if (event==FL_FOCUS || event==FL_UNFOCUS)
      return 1;
    if (event==FL_KEYBOARD){
      theta_z -= int(theta_z/360)*360;
      theta_x -= int(theta_x/360)*360;
      theta_y -= int(theta_y/360)*360;
      switch (Fl::event_text()?Fl::event_text()[0]:0){
      case 'y':
	theta_z += delta_theta;
	q=q*euler_deg_to_quaternion_double(delta_theta,0,0);
	redraw();
	return 1;
      case 'z':
	theta_x += delta_theta;
	q=q*euler_deg_to_quaternion_double(0,delta_theta,0);
	redraw();
	return 1;
      case 'x':
	theta_y += delta_theta;
	q=q*euler_deg_to_quaternion_double(0,0,delta_theta);
	redraw();
	return 1;
      case 'Y':
	theta_z -= delta_theta;
	q=q*euler_deg_to_quaternion_double(-delta_theta,0,0);
	redraw();
	return 1;
      case 'Z':
	theta_x -= delta_theta;
	q=q*euler_deg_to_quaternion_double(0,-delta_theta,0);
	redraw();
	return 1;
      case 'X':
	theta_y -= delta_theta;
	q=q*euler_deg_to_quaternion_double(0,0,-delta_theta);
	redraw();
	return 1;
      case 'u': case 'f':
	depth += 0.01;
	if (!opengl){
	  if (depth>0.99) depth=0.99;
	  current_depth=depth;
	}
	mouse_position->redraw();
	redraw();
	return 1;
      case 'd': case 'n':
	depth -= 0.01;
	if (!opengl){
	  if (depth<-0.99) depth=-0.99;
	  current_depth=depth;
	}
	mouse_position->redraw();
	redraw();
	return 1;
      case 'h': case 'H': // hide below depth
	below_depth_hidden=true;
	redraw();
	return 1;
      case 's': case 'S': // show below depth
	below_depth_hidden=false;
	redraw();
	return 1;
      case 'k': case 'K': 
	{// kill 3d window (for Mac)
	  int hp_pos;
	  History_Pack * hp=get_history_pack(this,hp_pos);
	  if (hp){
	    Fl_Widget * g=hp->child(hp_pos);
	    if (Fl_Group * gr=dynamic_cast<Fl_Group *>(g)){
	      g=gr->child(0);
	      hp->add_entry(hp_pos+1,g);
	    }
	    hp->remove_entry(hp_pos);
	  }
	  //delete this;
	  return 1;	
	}
      case 'o': case 'O':
	switch_renderer();
	redraw();
	return 1;
      case 'p': case 'P': // switch display hint, useful for printing
	displayhint=!displayhint;
	redraw();
	return 1;
      }      
    }
    current_i=Fl::event_x();
    current_j=Fl::event_y();
    current_depth=depth;
    double x,y,z;
    find_xyz(current_i,current_j,current_depth,x,y,z);
    //std::cout << current_i << " " << current_j << " " << current_depth << " " << x << " " << y << " " << z << '\n';
    in_area=(x>=window_xmin) && (x<=window_xmax) &&
      (y>=window_ymin) && (y<=window_ymax) && 
      (z>=window_zmin) && (z<=window_zmax);
    if (event==FL_MOUSEWHEEL){
      if (!Fl::event_inside(this))
	return 0;
      if (!opengl){
	depth -= Fl::e_dy*0.01;
	if (depth<-0.99) depth=-0.99;
	if (depth>0.99) depth=0.99;
	current_depth=depth;
	redraw();
	return 1;
      }
      if (show_axes){ // checking in_area seems to difficult, especially if mouse plan is not viewed
	depth = int(1000*(depth-Fl::e_dy *0.01)+.5)/1000.0;
	mouse_position->redraw();
	redraw();
	return 1;
      }
      else {
	if (Fl::e_dy<0)
	  zoom(0.8);
	else
	  zoom(1.25);
	return 1;
      }
    }
    if ( (event==FL_PUSH || push_in_area))
      ;
    else
      in_area=false;
    if (event==FL_PUSH){
      if (this!=Fl::focus()){
	Fl::focus(this);
	handle(FL_FOCUS);
      }
      push_i=current_i;
      push_j=current_j;
      push_depth=current_depth;
      push_in_area = in_area;
      pushed = true;
      return 1;
    }
    if (!(display_mode & 0x80) && push_in_area && (event==FL_DRAG || event==FL_RELEASE)){
      double x1,y1,z1,x2,y2,z2;
      find_xyz(current_i,current_j,current_depth,x1,y1,z1);
      find_xyz(push_i,push_j,push_depth,x2,y2,z2);
      double newx=x1-x2, newy=y1-y2, newz=z1-z2;
      round3(newx,window_xmin,window_xmax);
      round3(newy,window_ymin,window_ymax);      
      round3(newz,window_zmin,window_zmax);
      window_xmin -= newx;
      window_xmax -= newx;
      window_ymin -= newy;
      window_ymax -= newy;
      window_zmin -= newz;
      window_zmax -= newz;
      push_i = current_i;
      push_j = current_j;
      redraw();
      if (event==FL_RELEASE)
	pushed=false;
      return 1;
    }
    if (event==FL_DRAG){
      if (push_in_area)
	return 0;
      dragi=current_i-push_i;
      dragj=current_j-push_j;
      redraw();
      return 1;
    }
    if (event==FL_RELEASE){
      pushed = false;
      if (push_in_area)
	return 0;
      dragi=current_i-push_i;
      dragj=current_j-push_j;
      if (paused && absint(dragi)<4 && absint(dragj)<4)
	++animation_instructions_pos;
      else {
	if (opengl)
	  q=quaternion_double(dragi*180/h(),0,0)*rotation_2_quaternion_double(1,0,0,dragj*180/w())*q;
	else
	  q=rotation_2_quaternion_double(0,0,1,dragi*180/h())*rotation_2_quaternion_double(0,1,0,dragj*180/w())*q;
      }
      dragi=dragj=0;
      redraw();
      return 1;
    }
    return res;
  }

  Geo3d::Geo3d(int x,int y,int width, int height, History_Pack * _hp): Graph3d(x,y,width,height,0,_hp) {
    displayhint=false;
    hp=_hp?_hp:geo_find_history_pack(this);
    if (hp){
      hp->eval_below=true;
      hp->eval_below_once=false;
      show_mouse_on_object=true;
      update_rotation();
      //default_color(FL_WHITE,hp->contextptr);
    }
  }

  int Geo3d::in_handle(int event){
    if (!event)
      return 1;
    if (!hp)
      hp=geo_find_history_pack(this);
    // cerr << event << " " << mode << '\n';
    int res=Graph3d::in_handle(event);
    if (event==FL_UNFOCUS){
      return 1;
    }
    if (event==FL_FOCUS){
      return 1;
    }
    if (int gres=geo_handle(event))
      return gres;
    // event not handeld by geometry
    return res;
  }

  void Geo3d::draw(){
    Graph3d::draw();
  }

  // set evryone to x
  void Graph3d::orthonormalize(){ 
    double dx=window_xmax-window_xmin,dy=window_ymax-window_ymin,dz=window_zmax-window_zmin;
    double d=giac_max(dx,giac_max(dy,dz));
    zoomx(d/dx);
    zoomy(d/dy);
    zoomz(d/dz);
    redraw();
  }

  void Graph3d::geometry_round(double x,double y,double z,double eps,gen & tmp,GIAC_CONTEXT) {
    tmp= geometry_round_numeric(x,y,z,eps,approx);
    if (tmp.type==_VECT)
      tmp.subtype = _POINT__VECT;
    selected=nearest_point(plot_instructions,tmp,eps,contextptr);
    // if there is a point inside selected, stop there, 
    // otherwise find a point that is near the line
    int pos=findfirstpoint(selection2vecteur(selected));
    if (pos<0){
      // line passing through tmp      
      double a,b,c;
      double i,j,k;
      find_ij(x,y,z,i,j,k);
      k--;
      find_xyz(i,j,k,a,b,c);
      gen line(makevecteur(tmp,makevecteur(a,b,c)),_LINE__VECT);
      /*
	current_normal(a,b,c);
	vecteur v(makevecteur(a,b,c));
	gen line(makevecteur(tmp,tmp+v),_LINE__VECT);
      */
      vector<int> sel2=nearest_point(plot_instructions,line,eps,contextptr);
      pos=findfirstpoint(selection2vecteur(sel2));
      if (pos>=0){
	selected.insert(selected.begin(),sel2[pos]);
      }
      else {
	vector<int>::const_iterator it=sel2.begin(),itend=sel2.end();
	for (;it!=itend;++it)
	  selected.push_back(*it);
      }
      if (selected.empty()){
	// add hyperplans
	const_iterateur it=plot_instructions.begin(),itend=plot_instructions.end();
	for (int pos=0;it!=itend;++it,++pos){
	  gen tmp=remove_at_pnt(*it);
	  if (tmp.is_symb_of_sommet(at_hyperplan)){
	    vecteur v=interdroitehyperplan(line,tmp,contextptr);
	    if (!v.empty() && !is_undef(v.front())){
	      gen inters=evalf_double(remove_at_pnt(v.front()),1,contextptr);
	      if (inters.type==_VECT && inters._VECTptr->size()==3){
		vecteur & xyz=*inters._VECTptr;
		if (xyz[0].type==_DOUBLE_ && xyz[1].type==_DOUBLE_ && xyz[2].type==_DOUBLE_){
		  double x=xyz[0]._DOUBLE_val;
		  double y=xyz[1]._DOUBLE_val;
		  double z=xyz[2]._DOUBLE_val;
		  if (x>=window_xmin && x<=window_xmax &&
		      y>=window_ymin && y<=window_ymax &&
		      z>=window_zmin && z<=window_zmax)
		    selected.push_back(pos);
		}
	      }
	    }
	  }
	}
      }
    }
  }

#ifndef NO_NAMESPACE_XCAS
} // namespace giac
#endif // ndef NO_NAMESPACE_XCAS

#endif // HAVE_LIBFLTK_GL
