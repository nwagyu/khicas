#include "os.h"
#include "giac.h"
#include "kdisplay.h"

using namespace giac;
using namespace std;

__attribute__((noinline)) int main_(int argc,char ** argv){
  cfg_register_fileext("xw", "khicas");
  cfg_register_fileext("py", "khicas");
  lcd_init(lcd_type()); // clrscr();
#ifdef MICROPY_LIB
  mp_stack_ctrl_init();
  // python_heap=micropy_init(128*1024,1024*1024);
  // if (!python_heap) return 1;
#endif
  vx_var=identificateur("x");
  giac::context * cptr=(giac::context *)caseval("caseval contextptr");
#if 1
  xcas::console_main(cptr,argc>1?argv[1]:"session");
#ifdef MICROPY_LIB
  mp_deinit();
  if (python_heap)
    free(python_heap);
#endif
  return 0;
#else
  // identificateur x("x");
  // COUT << _factor(pow(x,4,cptr)-1,cptr) << endl;
  COUT << "Enter expression to eval, 0 to quit, ?command for help " << endl;
  for (unsigned i=0;;++i){
    //console_cin.foreground_color(nio::COLOR_RED);
    COUT << i << ">> ";
    //console_cin.foreground_color(nio::COLOR_BLACK);
    string s;
    CIN >> s;
#ifdef MICROPY_LIB
    micropy_eval(s.c_str()); 
#else
    //console_cin.foreground_color(nio::COLOR_GREEN);
    ctrl_c=interrupted=false;
    giac::gen g(s,cptr);
    if (g==0)
      return 0;
    // COUT << "type " << g.type << endl;
    // if (g.type==_SYMB) COUT << g._SYMBptr->sommet << endl;
    // COUT << "before eval " << g << endl; wait_key_pressed();
    g=eval(g,cptr);
    //console_cin.foreground_color(nio::COLOR_BLUE);
    COUT << g << endl;
#endif
  }
#endif
#ifdef MICROPY_LIB
  mp_deinit();
  if (python_heap)
    free(python_heap);
#endif
  return 0;
}

int main(int argc,char ** argv){
#ifdef MICROPY_LIB
  mp_stack_ctrl_init();
#endif
  main_(argc,argv);
}

