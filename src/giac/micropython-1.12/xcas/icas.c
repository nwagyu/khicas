#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
/* micropython alone test file */

void console_output(const char * s){
  printf("%s",s);
}

const char * console_input(){
  printf("%s",">>> ");
  static char buf[1024];
  int i=0;
  for (;i<sizeof(buf)-1;++i){
    char c=getc(stdin);
    if (c==10 || c==13)
      break;
    buf[i]=c;
  }
  buf[i]=0;
  return buf;
}

const char * read_file(const char * filename){
  FILE * f =fopen(filename,"r");
  static char buf[256*1024];
  fscanf(f,"%s",buf);
  return buf;
}

int getkey(int allow_suspend){
  return getc(stdin);
}

int file_exists(const char * filename){
  if (access(filename,R_OK))
    return 0;
  return 1;
}

int turtle_freeze=0;

void sync_screen(){
}

const char * caseval(const char * s){
  return s;
}

void c_sprint_double(char * s,double d){
  sprintf("%.14g",s,d);
}

void c_draw_rectangle(int x,int y,int w,int h,int c){
}
void c_draw_line(int x0,int y0,int x1,int y1,int c){
}
void c_fill_rect(int x,int y,int w,int h,int c){
}
int c_draw_string(int x,int y,int c,int bg,const char * s,int fake){
  return 0;
}
int c_draw_string_small(int x,int y,int c,int bg,const char * s,int fake){
  return 0;
}
int c_draw_string_medium(int x,int y,int c,int bg,const char * s,int fake){
  return 0; 
}
void c_draw_circle(int xc,int yc,int r,int color,int q1,int q2,int q3,int q4){
}
void c_draw_filled_circle(int xc,int yc,int r,int color,int left,int right){
}
void c_draw_polygon(int * x,int *y ,int n,int color){
}
void c_draw_filled_polygon(int * x,int *y, int n,int xmin,int xmax,int ymin,int ymax,int color){
}
int os_get_pixel(int x,int y){
  return 0; 
}
void c_set_pixel(int x,int y,int c){
}
void c_draw_arc(int xc,int yc,int rx,int ry,int color,double theta1, double theta2){
}
void c_draw_filled_arc(int x,int y,int rx,int ry,int theta1_deg,int theta2_deg,int color,int xmin,int xmax,int ymin,int ymax,int segment){
}
struct double_pair {
  double r,i;
} ;
typedef struct double_pair c_complex;

int c_inv(c_complex * x,int n){
  return 0;
}

c_complex c_det(c_complex *x,int n){
  return x[0];
}

int c_rref(c_complex * x,int n,int m){
  return 0;
}

int c_egv(c_complex * x,int n){
  return 0;
}

int c_eig(c_complex * x,c_complex * d,int n){
  return 0;
}

int c_proot(c_complex * x,int n){
  return 0;
}

int c_pcoeff(c_complex * x,int n){
  return 0;
}

int c_fft(c_complex * x,int n,int inverse){
  return 0;
}
int ctrl_c_interrupted(){
  return 0;
}


#include "py/mpstate.h"
#include "py/nlr.h"
#include "py/compile.h"
#include "py/runtime.h"
#include "py/builtin.h"
#include "py/repl.h"
#include "py/gc.h"
#include "py/stackctrl.h"
#include "py/mphal.h"
#include "py/mpthread.h"
#include "py/mperrno.h"
#include "extmod/misc.h"
#include "genhdr/mpversion.h"
#include "input.h"


MP_NOINLINE int main_(int argc, char **argv);

int main(int argc, char **argv) {
    #if MICROPY_PY_THREAD
    mp_thread_init();
    #endif
    // We should capture stack top ASAP after start, and it should be
    // captured guaranteedly before any other stack variables are allocated.
    // For this, actual main (renamed main_) should not be inlined into
    // this function. main_() itself may have other functions inlined (with
    // their own stack variables), that's why we need this main/main_ split.
    mp_stack_ctrl_init();
    return main_(argc, argv);
}
