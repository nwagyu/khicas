// -*- compile-command Debug : "cl /Zi /MDd /I. /EHsc giac_oo.cpp giac_oo.lib mpfr.lib mpir.lib" -*-
// -*- compile-command Release : "cl /MD /I. /EHsc giac_oo.cpp giac_oo.lib mpfr.lib mpir.lib" -*-

#include <giac/config.h>
#include <giac/gen.h>
#include <giac/prog.h>

int main(){
  std::string s;
  giac::context ctx;
  for (;;){
    std::cout << "Expression : " ;
    std::getline(std::cin,s);
    giac::gen g(s,&ctx);
    if (is_zero(g))
      return 0;
    g=giac::protecteval(g,1,&ctx);
    std::cout << "Sortie giac : " << g.print(&ctx) << std::endl;
  }
}
