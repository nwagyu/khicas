#include "giacPCH.h"
#include <os.h>
#include <lauxlib.h>
#include "luabridge.h"
#include "kdisplay.h"

using namespace giac;
//const char std::endl='\n';
//nio::console * std::console_cin_ptr=0;
using namespace std;

static int initialized=0;

void luagiac_free(){
#if 0
  if (nspire_exam_mode==1){
    set_exam_mode(3,contextptr); exam_mode=0;
  }
#endif
#ifdef MICROPY_LIB
  python_free();
#endif
  xcas::Console_Free();
  giac::release_globals();
}

struct bridge_bidon_t{
  bridge_bidon_t(){
  }
  ~bridge_bidon_t() {
    if (initialized)
      luagiac_free();
  }
};

bridge_bidon_t bridge_bidon;

void luagiac_init(){
  unsigned green=*(unsigned *) 0x90110b04;
  unsigned red=*(unsigned *) 0x90110b0c;
  if (green || red)
    nspire_exam_mode=1; 
  nspirelua=1;
  giac::context * contextptr=(giac::context *)giac::caseval("caseval contextptr");
#ifdef MICROPY_LIB
  giac::micropy_ptr=micropy_ck_eval;
#endif
  freeze=true;
  python_heap=0;
  xcas::sheetptr=0;
  shutdown=do_shutdown;
  // SetQuitHandler(save_session); // automatically save session when exiting
  if (!turtleptr){
    turtle();
    giac::_efface_logo(giac::vecteur(0),contextptr);
  }
  giac::caseval("floor"); // init xcas parser for Python syntax coloration (!)
  int key;
  xcas::Console_Init(contextptr);
  rand_seed(millis(),contextptr);
  giac::angle_radian(os_get_angle_unit()==0,contextptr);
#if 0
  if (nspire_exam_mode){ // must save LED state for restoration at end
    set_exam_mode(2,contextptr); exam_mode=0;
  }
#endif
  // GetKey(&key);
  // Console_Disp(1,contextptr);
}

const char * giac_caseval(const char * s){
  if (!initialized){
    luagiac_init();
    initialized=1;
  }
  vx_var=identificateur("x");
  //static nio::console console_cin;
  //console_cin_ptr=&console_cin;
  //giac::debug_infolevel=2;
  int intmask = TCT_Local_Control_Interrupts(-1); // disable
  const char * res=giac::caseval(s);
  reset_gc();
  TCT_Local_Control_Interrupts(intmask); // restore (0 to enable)
  return res;
}

