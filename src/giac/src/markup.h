/*  markup.h
 *
 *  copyright (c) 2019 Luka Marohnić
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

#ifndef MARKUP_H
#define MARKUP_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "first.h"
#include "gen.h"

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

void enable_texmacs_compatible_latex_export(bool yes);
string export_latex(const gen &g,GIAC_CONTEXT);
bool has_improved_latex_export(const gen &g,string &s,bool override_texmacs,GIAC_CONTEXT);
string export_mathml(const gen &g,GIAC_CONTEXT);
string export_mathml_presentation(const gen &g,GIAC_CONTEXT);
string export_mathml_content(const gen &g,GIAC_CONTEXT);
gen _export_mathml(const gen &g,GIAC_CONTEXT);
gen _xml_print(const gen &g,GIAC_CONTEXT);

#if 0
enum EqElementType {
  EQ_ELEMENT_TEXT=1,
  EQ_ELEMENT_APPLY=2,
  EQ_ELEMENT_POSTFIX=3,
  EQ_ELEMENT_OP=4,
  EQ_ELEMENT_INDEX=5,
  EQ_ELEMENT_POWER=6,
  EQ_ELEMENT_ROOT=7,
  EQ_ELEMENT_OVER=8,
  EQ_ELEMENT_TABLE=9,
  EQ_ELEMENT_LIST=10,
  EQ_ELEMENT_BIG=11
};

enum EqTextType {
  EQ_TEXT_NORMAL=1,
  EQ_TEXT_ITALIC=2,
  EQ_TEXT_BOLD=3,
  EQ_TEXT_SCRIPT=4,
  EQ_TEXT_FRAKTUR=5,
  EQ_TEXT_BLACKBOARD=6,
  EQ_TEXT_GREEK=7
};

enum EqBracketType {
  EQ_BRACKET_ROUND=1,
  EQ_BRACKET_SQUARE=2,
  EQ_BRACKET_CURLY=3,
  EQ_BRACKET_STRAIGHT=4
};

enum EqOperatorType {
  EQ_OPERATOR_PLUS=1,
  EQ_OPERATOR_MINUS=2,
  EQ_OPERATOR_CDOT=3,
  EQ_OPERATOR_INVISIBLE_TIMES=4,
  EQ_OPERATOR_TIMES=5,
  EQ_OPERATOR_CIRC=6,
  EQ_OPERATOR_POINTPLUS=7,
  EQ_OPERATOR_POINTMINUS=8,
  EQ_OPERATOR_POINTPROD=9,
  EQ_OPERATOR_POINTDIV=10,
  EQ_OPERATOR_AND=11,
  EQ_OPERATOR_OR=12,
  EQ_OPERATOR_XOR=13,
  EQ_OPERATOR_UNION=14,
  EQ_OPERATOR_INTERSECT=15,
  EQ_OPERATOR_SETMINUS=16,
  EQ_OPERATOR_USER=17
};

enum EqBigOperatorType {
    EQ_BIG_SUM=1,
    EQ_BIG_PROD=2
};

enum EqSelectionType {
  EQ_SELECTION_ATOM=1,
  EQ_SELECTION_SINGLE=2,
  EQ_SELECTION_MULTIPLE=3
};

class EqElement {
protected:
  EqElement *_parent;
  int _left_bracket_type;
  int _right_bracket_type;
  int _selection_type;
  gen _expr;
  double _rel_x;
  double _rel_y;
  double _width;
  double _height;
  double _depth;
  bool _is_overlined;
public:
  EqElement(const gen &e,EqElement *p) {
    _expr=e;
    _parent=p;
    _left_bracket_type=0;
    _right_bracket_type=0;
    _selection_type=0;
    _width=_height=_depth=0;
    _rel_x=_rel_y=0;
    _is_overlined=false;
  }
  EqElement* const parent() { return _parent; }
  virtual int type() const { return 0; }
  void parenthesize(int btype=EQ_BRACKET_ROUND) {
    _left_bracket_type=btype;
    _right_bracket_type=btype;
  }
  void set_left_bracket(int btype) { _left_bracket_type=btype; }
  void set_right_bracket(int btype) { _right_bracket_type=btype; }
  int left_bracket() const { return _left_bracket_type; }
  int right_bracket() const { return _right_bracket_type; }
  int selection_type() const { return _selection_type; }
  void set_selection_type(int t) { _selection_type=t; }
  double rel_x() const { return _rel_x; }
  double rel_y() const { return _rel_y; }
  double width() const { return _width; }
  double height() const { return _height; }
  double depth() const { return _depth; }
  void set_size(double w,double h,double d) {
    _width=w;
    _height=h;
    _depth=d;
  }
  void set_position(double x,double y) {
    _rel_x=x;
    _rel_y=y;
  }
  void set_overlined(bool enable) { _is_overlined=enable; }
  bool is_overlined() const { return _is_overlined; }
};

class EqElementText : public EqElement {
private:
  string _text;
  int _text_type;
public:
  EqElementText(const gen &e,EqElement *p,const string &s,int t) : EqElement(e,p) {
    _text=s;
    _text_type=t;
  }
  int type() const { return EQ_ELEMENT_TEXT; }
  const string &text() const { return _text; }
  int text_type() const { return _text_type; }
};

class EqElementApply : public EqElement {
private:
  EqElement *_func;
  EqElement *_arg;
  bool _insert_space;
public:
  EqElementApply(const gen &e,EqElement *p,EqElement *f,EqElement *a,bool isp=true) : EqElement(e,p) {
    _func=f;
    _arg=a;
    _insert_space=isp;
    _selection_type=EQ_SELECTION_ATOM;
  }
  ~EqElementApply() {
    delete _func;
    delete _arg;
  }
  int type() const { return EQ_ELEMENT_APPLY; }
  EqElement* function() const { return _func; }
  EqElement* argument() const { return _arg; }
  bool insert_space() const { return _insert_space; }
};

class EqElementOp: public EqElement {
private:
  int _nops;
  EqElement **_ops;
  int *_opr;
  string _opr_user;
public:
  EqElementOp(const gen &e,EqElement *p,int n=2) : EqElement(e,p) {
    assert(n>1);
    _nops=n;
    _ops=(EqElement**)malloc(_nops*sizeof(EqElement*));
    _opr=new int[_nops-1];
  }
  ~EqElementOp() {
    for (int i=0;i<_nops;++i) {
      delete _ops[i];
    }
    free(_ops);
    delete[] _opr;
  }
  int type() const { return EQ_ELEMENT_OP; }
  int nops() const { return _nops; }
  void set_nth_operand(EqElement *e,int n) { assert(n<_nops && n>=0); _ops[n]=e; }
  void set_nth_operator(int opr,int n=1) { assert(n<_nops && n>=1); _opr[n-1]=opr; }
  void set_user_operator(const string &s) {
    assert(_nops==2);
    delete[] _opr;
    _opr=NULL;
    _opr_user=s;
  }
  EqElement* nth_operand(int n) const { assert(n<_nops && n>=0); return _ops[n]; }
  int nth_operator(int n=1) const {
    assert(n<_nops && n>=1);
    if (_opr==NULL) return 0;
    else return _opr[n-1];
  }
};

class EqElementPostfix : public EqElement {
private:
  EqElement *_post;
  EqElement *_arg;
public:
  EqElementPostfix(const gen &e,EqElement *p,EqElement *f,EqElement *a) : EqElement(e,p) {
    _post=f;
    _arg=a;
    _selection_type=EQ_SELECTION_ATOM;
  }
  ~EqElementPostfix() {
    delete _post;
    delete _arg;
  }
  int type() const { return EQ_ELEMENT_POSTFIX; }
  EqElement* postfix() const { return _post; }
  EqElement* argument() const { return _arg; }
};

class EqElementIndex : public EqElement {
private:
  EqElement *_base;
  EqElement *_sub;
public:
  EqElementIndex(const gen &e,EqElement *p,EqElement *b,EqElement *s) : EqElement(e,p) {
    _base=b;
    _sub=s;
  }
  ~EqElementIndex() {
    delete _base;
    delete _sub;
  }
  int type() const { return EQ_ELEMENT_INDEX; }
  EqElement* base() const { return _base; }
  EqElement* index() const { return _sub; }
};

class EqElementPower : public EqElement {
private:
  EqElement *_base;
  EqElement *_sup;
public:
  EqElementPower(const gen &e,EqElement *p,EqElement *b,EqElement *s) : EqElement(e,p) {
    _base=b;
    _sup=s;
  }
  ~EqElementPower() {
    delete _base;
    delete _sup;
  }
  int type() const { return EQ_ELEMENT_POWER; }
  EqElement* base() const { return _base; }
  EqElement* exponent() const { return _sup; }
};

class EqElementRoot : public EqElement {
private:
  EqElement *_radic;
  EqElement *_deg;
public:
  EqElementRoot(const gen &e,EqElement *p,EqElement *r,EqElement *d=NULL) : EqElement(e,p) {
    _radic=r;
    _radic->set_overlined(true);
    _deg=d;
  }
  ~EqElementRoot() {
    delete _radic;
    if (_deg!=NULL) delete _deg;
  }
  int type() const { return EQ_ELEMENT_ROOT; }
  EqElement* radicand() const { return _radic; }
  EqElement* degree() const { return _deg; }
  bool is_square_root() const { return _deg==NULL; }
};

class EqElementOver : public EqElement {
private:
  EqElement *_numer;
  EqElement *_denom;
public:
  EqElementOver(const gen &e,EqElement *p,EqElement *num,EqElement *den) : EqElement(e,p) {
    _numer=num;
    _denom=den;
    _selection_type=EQ_SELECTION_SINGLE;
  }
  ~EqElementOver() {
    delete _numer;
    delete _denom;
  }
  int type() const { return EQ_ELEMENT_OVER; }
  EqElement* numerator() const { return _numer; }
  EqElement* denominator() const { return _denom; }
};

class EqElementBig : public EqElement {
private:
  EqElement *_lbound;
  EqElement *_ubound;
  EqElement *_formula;
  int _op_type;
public:
  EqElementBig(const gen &e,EqElement *p,EqElement *f,EqElement *lb,EqElement *ub,int t) : EqElement(e,p) {
    assert(f!=NULL);
    _formula=f;
    _lbound=lb;
    _ubound=ub;
    _op_type=t;
    _selection_type=EQ_SELECTION_SINGLE;
  }
  ~EqElementBig() {
    delete _formula;
    if (_lbound!=NULL) delete _lbound;
    if (_ubound!=NULL) delete _ubound;
  }
  int type() const { return EQ_ELEMENT_BIG; }
  EqElement *lower_bound() const { return _lbound; }
  EqElement *upper_bound() const { return _ubound; }
  EqElement *formula() const { return _formula; }
  int operator_type() const { return _op_type; }
};
#endif

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC

#endif // MARKUP_H
