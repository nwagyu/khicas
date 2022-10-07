/*
 * optimization.cc
 *
 * Copyright 2021-2022 Luka Marohnić
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "giacPCH.h"
#include "giac.h"
#include "optimization.h"
#include "lpsolve.h"
#include "graphe.h"
#include "signalprocessing.h"
#include <bitset>
#include <string>
#ifdef HAVE_PARI_PARI_H
#include <pari/pari.h>
#endif

using namespace std;

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

/* Create a temporary symbol "name<count>". */
gen temp_symb(const string &name,int count,GIAC_CONTEXT) {
    gen ret=identificateur(" "+name+(count>=0?print_INT_(count):""));
    _purge(ret,contextptr);
    return ret;
}

/* Make a matrix of zeros. */
matrice zero_mat(int rows,int cols,GIAC_CONTEXT) {
    return *_matrix(makesequence(rows,cols,0),contextptr)._VECTptr;
}

/* Copy matrix SRC to DEST, resizing the latter if necessary. */
void copy_matrice(const matrice &src,matrice &dest) {
    int m=mrows(src),n=mcols(src);
    if (int(dest.size())!=m) dest.resize(m);
    const_iterateur it=src.begin(),itend=src.end();
    iterateur jt=dest.begin();
    for (;it!=itend;++it,++jt) {
        if (jt->type!=_VECT)
            *jt=vecteur(it->_VECTptr->begin(),it->_VECTptr->end());
        else {
            if (int(jt->_VECTptr->size())!=n) jt->_VECTptr->resize(n);
            std::copy(it->_VECTptr->begin(),it->_VECTptr->end(),jt->_VECTptr->begin());
        }
    }
}

/* Return A+B (add two matrices). */
matrice madd(const matrice &a,const matrice &b) {
    matrice ret;
    ret.reserve(a.size());
    const_iterateur at=a.begin(),bt=b.begin();
    for (;at!=a.end()&&bt!=b.end();++at,++bt) {
        ret.push_back(addvecteur(*at->_VECTptr,*bt->_VECTptr));
    }
    return ret;
}
void madd_inplace(matrice &a,const matrice &b,bool subt) {
    iterateur at=a.begin(),ajt;
    const_iterateur bt=b.begin(),bjt;
    for (;at!=a.end()&&bt!=b.end();++at,++bt) {
        ajt=at->_VECTptr->begin();
        bjt=bt->_VECTptr->begin();
        for (;ajt!=at->_VECTptr->end()&&bjt!=bt->_VECTptr->end();++ajt,++bjt) {
            if (subt) *ajt-=*bjt; else *ajt+=*bjt;
        }
    }
}

/* Multiply each element of matrix B by A. */
matrice mscale(const gen &a,const matrice &b) {
    matrice ret;
    ret.reserve(b.size());
    for (const_iterateur it=b.begin();it!=b.end();++it) {
        ret.push_back(multvecteur(a,*it->_VECTptr));
    }
    return ret;
}
void mscale_inplace(const gen &a,matrice &b) {
    for (iterateur it=b.begin();it!=b.end();++it) {
        for (iterateur jt=it->_VECTptr->begin();jt!=it->_VECTptr->end();++jt) {
            *jt=*jt*a;
        }
    }
}

/* Return A-B (subtract two matrices). */
matrice msub(const matrice &a,const matrice &b) {
    return madd(a,mscale(-1,b));
}
void msub_inplace(matrice &a,const matrice &b) {
    madd_inplace(a,b,true);
}

/* Return the index of the first minimal/maximal element in V */
int argminmax(const vecteur &v,bool min,GIAC_CONTEXT) {
    gen val=min?plus_inf:minus_inf;
    int ret=-1;
    const_iterateur it,itstart=v.begin(),itend=v.end();
    for (it=itstart;it!=itend;++it) {
        if (is_strictly_greater(min?val:*it,min?*it:val,contextptr)) {
            ret=it-itstart;
            val=*it;
        }
    }
    return ret;
}
int argmin(const vecteur &v,GIAC_CONTEXT) {
    return argminmax(v,true,contextptr);
}
int argmax(const vecteur &v,GIAC_CONTEXT) {
    return argminmax(v,false,contextptr);
}

gen to_algebraic(const gen &g,vecteur &s,GIAC_CONTEXT) {
    if (g.type==_IDNT || is_fully_numeric(g,num_mask_withfrac|num_mask_withint))
        return g;
    const_iterateur it,itend;
    if (g.type==_VECT) {
        vecteur res;
        res.reserve(g._VECTptr->size());
        for (it=g._VECTptr->begin(),itend=g._VECTptr->end();it!=itend;++it)
            res.push_back(to_algebraic(*it,s,contextptr));
        return change_subtype(res,g.subtype);
    }
    if (g.is_symb_of_sommet(at_plus)) {
        gen ret(0);
        for (it=g._SYMBptr->feuille._VECTptr->begin(),itend=g._SYMBptr->feuille._VECTptr->end();it!=itend;++it)
            ret+=to_algebraic(*it,s,contextptr);
        return ret;
    }
    if (g.is_symb_of_sommet(at_prod)) {
        gen ret(1);
        for (it=g._SYMBptr->feuille._VECTptr->begin(),itend=g._SYMBptr->feuille._VECTptr->end();it!=itend;++it)
            ret=ret*to_algebraic(*it,s,contextptr);
        return ret;
    }
    if (g.is_symb_of_sommet(at_inv))
        return symb_inv(to_algebraic(g._SYMBptr->feuille,s,contextptr));
    if (g.is_symb_of_sommet(at_div))
        return to_algebraic(g._SYMBptr->feuille._VECTptr->front(),s,contextptr)/
               to_algebraic(g._SYMBptr->feuille._VECTptr->back(),s,contextptr);
    if (g.is_symb_of_sommet(at_neg))
        return -to_algebraic(g._SYMBptr->feuille,s,contextptr);
    if (g.is_symb_of_sommet(at_pow) && is_fully_numeric(g._SYMBptr->feuille._VECTptr->back(),num_mask_withfrac|num_mask_withint))
        return symb_pow(to_algebraic(g._SYMBptr->feuille._VECTptr->front(),s,contextptr),g._SYMBptr->feuille._VECTptr->back());
    gen u=temp_symb("toalg_tmp",s.size(),contextptr);
    vecteur cond;
    if (g.is_symb_of_sommet(at_exp))
        cond.push_back(symb_superieur_strict(u,0));
    else if (g.is_symb_of_sommet(at_cosh))
        cond.push_back(symb_superieur_egal(u,1));
    else if (g.is_symb_of_sommet(at_cos) || g.is_symb_of_sommet(at_sin) || g.is_symb_of_sommet(at_tanh))
        cond.push_back(makevecteur(-1,1));
    else if (g.is_symb_of_sommet(at_atan) || g.is_symb_of_sommet(at_asin))
        cond.push_back(makevecteur(-cst_pi/2,cst_pi/2));
    else if (g.is_symb_of_sommet(at_acos))
        cond.push_back(makevecteur(0,cst_pi));
    set_assumptions(u,cond,vecteur(0),false,contextptr);
    s.push_back(u);
    return u;
}
bool is_algebraic(const gen &g,GIAC_CONTEXT) {
    vecteur s;
    to_algebraic(g,s,contextptr);
    if (s.empty())
        return true;
    _purge(s,contextptr);
    return false;
}

/* Remove elements from V with indices in IND. */
void remove_elements_with_indices(vecteur &v,const set<int> &ind) {
    set<int>::const_reverse_iterator it=ind.rbegin(),itend=ind.rend();
    for (;it!=itend;++it)
        v.erase(v.begin()+*it);
}
void remove_elements_with_indices(vecteur &v,const set<pair<int,gen> > &ind) {
    set<pair<int,gen> >::const_reverse_iterator it=ind.rbegin(),itend=ind.rend();
    for (;it!=itend;++it)
        v.erase(v.begin()+it->first);
}
/* Return true iff g is an mapleconversion integer (corresponding to v) */
bool is_mcint(const gen &g,int v) {
    return g.is_integer() && g.subtype==_INT_MAPLECONVERSION && (v<0 || g.val==v);
}

/* Return the index of g in v, or -1 if g is not contained in v. */
int indexof(const gen &g,const vecteur &v) {
    const_iterateur it=find(v.begin(),v.end(),g);
    if (it==v.end())
        return -1;
    return it-v.begin();
}

/* Return the list of indices of linearly dependent rows in matrix m.
 * The returned list is sorted in ascending order. */
vector<int> linearly_dependent_rows(const matrice &m,GIAC_CONTEXT) {
    int nc=mcols(m);
    std::set<int> found;
    vector<int> res;
    res.reserve(m.size());
    matrice mf=mtran(*_rref(mtran(*_evalf(m,contextptr)._VECTptr),contextptr)._VECTptr);
    for (const_iterateur it=mf.begin();it!=mf.end();++it) {
        int i=0;
        for (;i<nc && is_zero(it->_VECTptr->at(i));++i);
        if (i<nc && is_one(it->_VECTptr->at(i)) && found.find(i)==found.end())
            found.insert(i);
        else res.push_back(it-mf.begin());
    }
    std::sort(res.begin(),res.end());
    return res;
}

/* Return true if variables in V are all identifiers. */
bool ckvars(const vecteur &v,GIAC_CONTEXT) {
    vecteur vp;
    for (const_iterateur it=v.begin();it!=v.end();++it) {
        if (it->type!=_IDNT || _eval(*it,contextptr)!=*it)
            return false;
    }
    return true;
}

/* Set assumptions on VARS w.r.t. BNDS. If OPEN=true, exclude bounds from range. */
void bound_variables(const vecteur &vars,const vecteur &bnds,bool open,GIAC_CONTEXT) {
    assert(vars.size()==bnds.size() && (vars.empty() || (ckmatrix(bnds) && bnds.front()._VECTptr->size()==2)));
    for (int i=vars.size();i-->0;) {
        const gen &v=vars[i],&lb=bnds[i]._VECTptr->front(),&ub=bnds[i]._VECTptr->back();
        if (!is_inf(lb) && !is_inf(ub))
            assume_t_in_ab(v,lb,ub,open,open,contextptr);
        else if (!is_inf(lb))
            giac_assume(open?symb_superieur_strict(v,lb):symb_superieur_egal(v,lb),contextptr);
        else if (!is_inf(ub))
            giac_assume(open?symb_inferieur_strict(v,ub):symb_inferieur_egal(v,ub),contextptr);
    }
}

/* Get the assumptions on identifier G: obtain the domain, feasible interval(s)
 * and excluded points. */
bool get_assumptions(const gen &g,int &dom,matrice &intervals,vecteur &excluded,GIAC_CONTEXT) {
    dom=-1;
    intervals.clear();
    excluded.clear();
    gen res=_about(g,contextptr);
    if (res.type==_IDNT && res==g)
        return true; // no assumptions
    if (res.type!=_VECT || res.subtype!=_ASSUME__VECT)
        return false;
    if (res._VECTptr->size()==1) {
        if (!res._VECTptr->front().is_integer())
            return false;
        dom=res._VECTptr->front().val;
        return true;
    }
    if (res._VECTptr->size()!=3)
        return false;
    if (res._VECTptr->front().is_integer())
        dom=res._VECTptr->front().val; // 1 - real, 2 - integer, 4 - complex, 10 - rational
    intervals=*res._VECTptr->at(1)._VECTptr;
    excluded=*res._VECTptr->at(2)._VECTptr;
    return true;
}

/* Set (additional) assumptions on identifier G: must be either inequality
 * constraints or bounds in form [L,U]. */
void set_assumptions(const gen &g,const vecteur &cond,const vecteur &excluded,bool additionally,GIAC_CONTEXT) {
    if (cond.empty())
        return;
    vecteur args;
    for (const_iterateur it=cond.begin();it!=cond.end();++it) {
        if (it->type==_VECT) {
            assert(it->_VECTptr->size()==2);
            const gen &lb=it->_VECTptr->front(),&ub=it->_VECTptr->back();
            gen l=is_inf(lb)?undef:(contains(excluded,lb)?symb_superieur_strict(g,lb):symb_superieur_egal(g,lb));
            gen r=is_inf(ub)?undef:(contains(excluded,ub)?symb_inferieur_strict(g,ub):symb_inferieur_egal(g,ub));
            if (!is_undef(l) && !is_undef(r))
                args.push_back(symb_and(l,r));
            else if (!is_undef(l))
                args.push_back(l);
            else if (!is_undef(r))
                args.push_back(r);
        } else args.push_back(*it);
    }
    gen a=symbolic(at_ou,change_subtype(args,_SEQ__VECT));
    if (additionally)
        giac_additionally(a,contextptr);
    else giac_assume(a,contextptr);
}

void purge_variables(const vecteur &vars,GIAC_CONTEXT) {
/* Purge identifiers in VARS. */
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        assert(it->type==_IDNT);
        _purge(*it,contextptr);
    }
}

/* Sort a list of elements such that letters and numbers are treated separately.
 * Useful for sorting lists of numbered identifiers, such as [x1,x2,...,x10,x11,...]. */
vecteur sort_identifiers(const vecteur &v,GIAC_CONTEXT) {
    if (v.empty())
        return v;
    vecteur strv=*_apply(makesequence(at_string,v),contextptr)._VECTptr,snv;
    for (const_iterateur it=strv.begin();it!=strv.end();++it) {
        if (it->type!=_STRNG)
            return v;
        const string &s=*it->_STRNGptr;
        size_t i=s.find_last_not_of("0123456789");
        string str_id=i==string::npos?"":s.substr(0,i+1),str_num=i==string::npos?s:s.substr(i+1);
        gen id=str_id.size()==0?undef:string2gen(str_id,false);
        gen num=_expr(string2gen(str_num,false),contextptr);
        vecteur sn;
        if (!is_undef(id))
            sn.push_back(id);
        if (!is_undef(num))
            sn.push_back(num);
        else if (!is_undef(id))
            sn.push_back(0);
        snv.push_back(sn);
    }
    return *mtran(*_sort(mtran(makevecteur(snv,v)),contextptr)._VECTptr)[1]._VECTptr;
}

/* Convert bounds given in matrix form [[L1,U1],[L2,U2],..] to the sequence
 * x1=L1..U1,x2=L2..U2,.. */
gen _box_constraints(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()!=2 ||
            g._VECTptr->front().type!=_VECT || g._VECTptr->back().type!=_VECT ||
            g._VECTptr->front()._VECTptr->empty() || g._VECTptr->back()._VECTptr->empty())
        return generrtype("Expected a sequence of two nonempty lists");
    const vecteur &x=*g._VECTptr->front()._VECTptr;
    const vecteur &b=*g._VECTptr->back()._VECTptr;
    if (!ckmatrix(b,false) || b.size()!=x.size() || b.front()._VECTptr->size()!=2)
        return generrdim("Invalid list of bounds");
    matrice B=mtran(b);
    gen intrv=_zip(makesequence(at_interval,B.front(),B.back()),contextptr);
    return change_subtype(_zip(makesequence(at_equal,x,intrv),contextptr),_SEQ__VECT);
}
static const char _box_constraints_s []="box_constraints";
static define_unary_function_eval (__box_constraints,&_box_constraints,_box_constraints_s);
define_unary_function_ptr5(at_box_constraints,alias_at_box_constraints,&__box_constraints,0,true)

/* Return the sum of elements in V (exclude the last one if DROP_LAST=true). */
int ipdiff::sum_ivector(const ivector &v,bool drop_last) {
    int res=0;
    for (ivector_iter it=v.begin();it!=v.end()-drop_last?1:0;++it) {
        res+=*it;
    }
    return res;
}

/* Return 1 if G is a <(=) contraint, 2 if G is a >(=) constraint, and 0 otherwise. */
int which_ineq(const gen &g) {
    if (g.is_symb_of_sommet(at_inferieur_egal) || g.is_symb_of_sommet(at_inferieur_strict))
        return 1;
    if (g.is_symb_of_sommet(at_superieur_egal) || g.is_symb_of_sommet(at_superieur_strict))
        return 2;
    return 0;
}

/* Compute the intersection RES of two sorted lists A and B, return the size of RES. */
int intersect(const vector<int> &a,const vector<int> &b,vector<int> &res) {
    res.resize(a.size()>b.size()?a.size():b.size());
    vector<int>::iterator it=set_intersection(a.begin(),a.end(),b.begin(),b.end(),res.begin());
    res.resize(it-res.begin());
    return res.size();
}

/* Compute the union RES of two sorted lists A and B, return the size of RES. */
int unite(const vector<int> &a,const vector<int> &b,vector<int> &res) {
    res.resize(a.size()+b.size());
    vector<int>::iterator it=set_union(a.begin(),a.end(),b.begin(),b.end(),res.begin());
    res.resize(it-res.begin());
    return res.size();
}

/* Simplify expression G. */
gen simp(const gen &g,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        vecteur res;
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            res.push_back(simp(*it,contextptr));
        }
        return change_subtype(res,g.subtype);
    }
    log_output_redirect lor(contextptr);
    gen ret=recursive_normal(g,contextptr);
#if 0
    if (_evalf(g,contextptr).is_approx()) {
        gen ret_simp=_simplify(g,contextptr);
        if (!has_rootof(ret_simp))
            return ret_simp;
        if (has_rootof(ret)) {
            ret=ratnormal(g,contextptr);
            return has_rootof(ret)?ret_simp:ret;
        }
    }
#endif
    return ret;
}

/* Return true if the expression E contains a non-elementary derivative w.r.t. X,
 * which must be an identifier or a list of identifiers. */
bool has_diff(const gen &e,const gen &x) {
    if (x.type==_VECT) {
        for (const_iterateur it=x._VECTptr->begin();it!=x._VECTptr->end();++it) {
            if (has_diff(e,*it))
                return true;
        }
        return false;
    }
    assert(x.type==_IDNT);
    if (e.type==_VECT) {
        for (const_iterateur it=e._VECTptr->begin();it!=e._VECTptr->end();++it) {
            if (has_diff(*it,x))
                return true;
        }
        return false;
    }
    if (e.type!=_SYMB)
        return false;
    if ((e.is_symb_of_sommet(at_derive) && e._SYMBptr->feuille.type==_VECT &&
            e._SYMBptr->feuille._VECTptr->size()>=2 && e._SYMBptr->feuille._VECTptr->at(1)==x) ||
            (e.is_symb_of_sommet(at_of) && e._SYMBptr->feuille._VECTptr->front().is_symb_of_sommet(at_derive) &&
                e._SYMBptr->feuille._VECTptr->back().type==_VECT && contains(*e._SYMBptr->feuille._VECTptr->back()._VECTptr,x)))
        return true;
    return has_diff(e._SYMBptr->feuille,x);
}

/* Return true iff G contains any of the following:
 * - absolute values
 * - piecewise, when and ifte functions
 * - Heaviside and Dirac functions
 * - max, min, floor, ceil, round functions
 * - sign function
 */
bool has_breaks(const gen &g) {
    if (g.type==_VECT) {
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (has_breaks(*it))
                return true;
        }
        return false;
    }
    if (g.type!=_SYMB)
        return false;
    const gen &s=g._SYMBptr->sommet;
    const gen &f=g._SYMBptr->feuille;
    if (s==at_abs || s==at_piecewise || s==at_when || s==at_ifte || s==at_Heaviside || s==at_Dirac || s==at_min ||
            s==at_max || s==at_floor || s==at_ceil || s==at_round || s==at_sign)
        return true;
    return has_breaks(f);
}

/* Brent's algorithm (implemented by John Burkardt), code below */
double local_min_rc(double &a,double &b,int &status,double value);

/* Return true if G = A*X^3+B*X^2+C*X+D. */
bool is_cubic_wrt(const gen &g,const gen &x,gen &a,gen &b,gen &c,gen &d,GIAC_CONTEXT) {
    if (is_quadratic_wrt(derive(g,x,contextptr),x,a,b,c,contextptr)) {
        a=a/3;
        b=b/2;
        d=g-a*x*x*x-b*x*x-c*x;
        return true;
    }
    return false;
}

/* Return true if G = A*X^4+B*X^3+C*X^2+D*x+E. */
bool is_quartic_wrt(const gen &g,const gen &x,gen &a,gen &b,gen &c,gen &d,gen &e,GIAC_CONTEXT) {
    if (is_cubic_wrt(derive(g,x,contextptr),x,a,b,c,d,contextptr)) {
        a=a/4;
        b=b/3;
        c=c/2;
        e=g-a*x*x*x*x-b*x*x*x-c*x*x-d*x;
        return true;
    }
    return false;
}

/* Return true if G is a quadratic with positive leading coefficient and
 * negative discriminant. */
bool is_positive_quadratic(const gen &g,bool strict,GIAC_CONTEXT) {
    vecteur p=*_lname(g,contextptr)._VECTptr;
    gen a,b,c,d;
    unsigned mt=100;
    for (const_iterateur it=p.begin();it!=p.end();++it) {
        if (!is_quadratic_wrt(g,*it,a,b,c,contextptr))
            continue;
        d=4*a*c-b*b;
        if (!is_exactly_zero(a) && is_positive_safe(a,false,mt,contextptr) && is_positive_safe(d,strict,mt,contextptr))
            return true;
    }
    return false;
}

/* Return true if G is a (strictly) positive quartic. */
bool is_positive_quartic(const gen &g,bool strict,GIAC_CONTEXT) {
    vecteur p=*_lname(g,contextptr)._VECTptr;
    gen a,b,c,d,e,dsc,dsc0,D,P,R;
    unsigned mt=100;
    for (const_iterateur it=p.begin();it!=p.end();++it) {
        if (!is_quartic_wrt(g,*it,a,b,c,d,e,contextptr) || is_positive_safe(-a,false,mt,contextptr))
            continue;
        dsc=((256*a*e-192*b*d-128*c*c)*e*e+(144*c*e-27*d*d)*d*d)*a*a+(144*c*e*e-6*d*d*e)*a*b*b
            +(-80*c*d*e+18*d*d*d)*a*b*c+(16*c*e-4*d*d)*a*c*c*c+(-27*b*e*e+18*c*d*e-4*d*d*d)*b*b*b+b*b*c*c*d*d;
        dsc0=c*c-3*b*d+12*a*e;
        D=64*a*a*a*e-16*a*a*c*c+16*a*b*b*c-16*a*a*b*d-3*b*b*b*b;
        P=8*a*c-3*b*b;
        R=b*b*b+8*d*a*a-4*a*b*c;
        bool pos_dsc=is_positive_safe(dsc,true,mt,contextptr),zero_dsc=is_exactly_zero(dsc);
        bool pos_P=is_positive_safe(P,true,mt,contextptr),neg_P=is_positive_safe(-P,true,mt,contextptr),pos_D=is_positive_safe(D,true,mt,contextptr);
        bool zero_R=is_exactly_zero(R),zero_D=is_exactly_zero(D),zero_dsc0=is_exactly_zero(dsc0);
        if ((pos_dsc && pos_D && pos_P) || (zero_dsc && zero_D && zero_R && pos_P) ||
                (!strict && (zero_dsc && ((pos_D || (pos_P && (!zero_D || !zero_R))) || (zero_D && (neg_P || zero_dsc0))))))
            return true;
    }
    return false;
}

/* Return true if G does not change sign (but may be zero). */
bool is_const_sign(const gen &g,bool pos,GIAC_CONTEXT) {
    return is_positive_safe(pos?g:-g,false,100,contextptr) ||
           is_positive_quadratic(pos?g:-g,false,contextptr) ||
           is_positive_quartic(pos?g:-g,false,contextptr);
}

/* Return true if g is *strictly* positive/negative given pos=true/false.
 * LIN should be set to true. */
bool is_const_sign_strict(const gen &g,bool pos,bool lin,GIAC_CONTEXT) {
    if (is_real_number(g,contextptr) && ((pos && is_strictly_positive(to_real_number(g,contextptr),contextptr)) ||
                                         (!pos && is_strictly_positive(-to_real_number(g,contextptr),contextptr))))
        return true;
    if (g.is_symb_of_sommet(at_neg))
        return is_const_sign_strict(g._SYMBptr->feuille,!pos,lin,contextptr);
    if (g.is_symb_of_sommet(at_abs))
        return is_const_sign_strict(g._SYMBptr->feuille,true,lin,contextptr) || is_const_sign_strict(g._SYMBptr->feuille,false,lin,contextptr);
    gen eps=exact(epsilon(contextptr),contextptr);
    if ((pos?is_greater(g,eps,contextptr):is_greater(-eps,g,contextptr)) ||
            is_positive_quadratic(pos?g:-g,true,contextptr) || is_positive_quartic(pos?g:-g,true,contextptr))
        return true;
    if (!lin && g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT) {
        const vecteur &f=*g._SYMBptr->feuille._VECTptr;
        vecteur s(f.size(),0);
        for (const_iterateur it=f.begin();it!=f.end();++it) {
            gen fac=*it,exponent;
            bool power=false;
            if (it->is_symb_of_sommet(at_pow) && (exponent=_abs(it->_SYMBptr->feuille._VECTptr->back(),contextptr)).is_integer()) {
                power=true;
                fac=it->_SYMBptr->feuille._VECTptr->front();
            }
            if (is_const_sign_strict(fac,true,false,contextptr))
                s[it-f.begin()]=1;
            else if (is_const_sign_strict(fac,false,false,contextptr))
                s[it-f.begin()]=power && exponent.val%2==0?1:-1;
        }
        gen ps=_product(s,contextptr);
        if (is_zero(ps))
            return false;
        return pos?is_one(ps):is_minus_one(ps);
    }
    vecteur terms(1,g);
    if (g.is_symb_of_sommet(at_plus) && g._SYMBptr->feuille.type==_VECT)
        terms=*g._SYMBptr->feuille._VECTptr;
    bool surplus=false;
    gen rest(0);
    for (const_iterateur it=terms.begin();it!=terms.end();++it) {
        gen l=lin?_lin(*it,contextptr):*it;
        if ((pos && (l.is_symb_of_sommet(at_exp) || l.is_symb_of_sommet(at_cosh))) ||
                ((l.is_symb_of_sommet(at_atan) || l.is_symb_of_sommet(at_asin) || l.is_symb_of_sommet(at_atanh) || l.is_symb_of_sommet(at_asinh))
                    && is_const_sign_strict(l._SYMBptr->feuille,pos,true,contextptr)) ||
                (pos && l.is_symb_of_sommet(at_acosh) && is_const_sign_strict(1-l._SYMBptr->feuille,true,true,contextptr)) ||
                (pos && (l.is_symb_of_sommet(at_ln) || l.is_symb_of_sommet(at_acosh))
                    && is_const_sign_strict(l._SYMBptr->feuille-1,true,true,contextptr)) ||
                (is_real_number(*it,contextptr) && is_strictly_positive((pos?1:-1)*to_real_number(*it,contextptr),contextptr)) ||
                (lin && is_const_sign_strict(l,pos,false,contextptr)))
            surplus=true;
        else if (surplus && is_const_sign(*it,pos,contextptr))
            ;
        else rest+=*it;
    }
    if (!surplus)
        return false;
    return is_const_sign(rest,pos,contextptr);
}
bool is_definitely_positive(const gen &g,GIAC_CONTEXT) {
    return is_const_sign_strict(g,true,true,contextptr);
}

/* Return true if im(g)=0. */
bool has_imag(const gen &g,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (has_imag(*it,contextptr))
                return true;
        }
        return false;
    }
    if (is_inf(g) || is_real_number(g,contextptr))
        return false;
    gen gr,gi;
    reim(g,gr,gi,contextptr);
    return !is_zero(simp(gi,contextptr),contextptr);
}

/* Make piecewise sub-expressions nested, e.g. replace
 * piecewise(cond1,g1,cond2,g2,...,otherwise) by
 * piecewise(cond1,g1,piecewise(cond2,g2,piecewise(...,piecewise(condn,gn,otherwise))...)). */
gen make_piecewise_nested(const gen &g) {
    if (g.type==_VECT) {
        vecteur res;
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            res.push_back(make_piecewise_nested(*it));
        }
        return change_subtype(res,g.subtype);
    } else if (g.type==_SYMB) {
        const gen &f=g._SYMBptr->feuille;
        if (g.is_symb_of_sommet(at_piecewise) && f.type==_VECT && f._VECTptr->size()>3) {
            const gen &cond=f._VECTptr->front(),&h=f._VECTptr->at(1);
            vecteur rest(f._VECTptr->begin()+2,f._VECTptr->end());
            return symbolic(at_piecewise,makesequence(make_piecewise_nested(cond),make_piecewise_nested(h),
                            make_piecewise_nested(symbolic(at_piecewise,change_subtype(rest,_SEQ__VECT)))));
        }
        return symbolic(g._SYMBptr->sommet,make_piecewise_nested(f));
    } else return g;
}

/* Replace the first found piecewise expression in G with PCW and store its
 * arguments in FEU. */
bool find_piecewise(gen &g,const gen &pcw,vecteur &feu) {
    if (g.type==_VECT) {
        for (iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (find_piecewise(*it,pcw,feu))
                return true;
        }
    } else if (g.type==_SYMB) {
        gen &f=g._SYMBptr->feuille;
        if (g._SYMBptr->sommet==at_piecewise) {
            feu=*f._VECTptr;
            g=pcw;
            return true;
        }
        return find_piecewise(f,pcw,feu);
    }
    return false;
}

/* Replace the first found absolute value in G with ASYMB and store its
 * argument in VAL. */
bool find_abs(gen &g,const gen &asymb,gen &val) {
    if (g.type==_VECT) {
        for (iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (find_abs(*it,asymb,val))
                return true;
        }
    } else if (g.type==_SYMB) {
        gen &f=g._SYMBptr->feuille;
        if (g._SYMBptr->sommet==at_abs) {
            val=f;
            g=asymb;
            return true;
        }
        return find_abs(f,asymb,val);
    }
    return false;
}

/* Return true if E is rational with respect to variables in VARS. */
bool is_rational_wrt_vars(const gen &e,const vecteur &vars) {
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        vecteur l(rlvarx(e,*it));
        if (l.size()>1)
            return false;
    }
    return true;
}

/* Return true if E is linear with respect to variables in VARS. */
bool is_linear_wrt_vars(const gen &e,const vecteur &vars,GIAC_CONTEXT) {
    gen a,b;
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        if (!is_linear_wrt(e,*it,a,b,contextptr) || !is_constant_wrt_vars(a,vars,contextptr))
            return false;
    }
    return true;
}

/* Return true if E is quadratic with respect to variables in VARS. */
bool is_quadratic_wrt_vars(const gen &e,const vecteur &vars,GIAC_CONTEXT) {
    gen a,b,c;
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        if (!is_quadratic_wrt(e,*it,a,b,c,contextptr) || !is_constant_wrt_vars(a,vars,contextptr) || !is_linear_wrt_vars(b,vars,contextptr))
            return false;
    }
    return true;
}

/* Return an uniform random real between A and B. */
gen rand_uniform(const gen &a,const gen &b,GIAC_CONTEXT) {
    gen fa=to_real_number(a,contextptr),fb=to_real_number(b,contextptr);
    if (is_zero(b-a,contextptr))
        return fa;
    return _rand(is_greater(b,a,contextptr)?makesequence(fa,fb):makesequence(fb,fa),contextptr);
}

/* Generate a normal variable with parameters MU=0 and SIGMA.
 * If ABSOLUT=true, return the absolute value. */
gen rand_normal(const gen &sigma,bool absolut,GIAC_CONTEXT) {
    assert(is_positive(sigma,contextptr));
    gen s=to_real_number(sigma,contextptr);
    if (is_zero(sigma,contextptr))
        return 0;
    gen r=_randNorm(makesequence(0,s),contextptr);
    return absolut?_abs(r,contextptr):r;
}

/* Generate a multi-dimensional normal random variable.
 * SIGMA is either a single number or a vector of the same length as MU. */
void rand_multinormal(const vecteur &mu,const gen &sigma,vecteur &res,GIAC_CONTEXT) {
    int n=mu.size(),i;
    res.resize(n);
    iterateur it=res.begin(),itend=res.end();
    const_iterateur mt=mu.begin();
    if (sigma.type==_VECT) {
        const_iterateur st=sigma._VECTptr->begin();
        for (;it!=itend;++it,++st,++mt) {
            *it=_randNorm(makesequence(*mt,*st),contextptr);
        }
    } else for (;it!=itend;++it,++mt) {
        *it=_randNorm(makesequence(*mt,sigma),contextptr);
    }
}

/* Generate categorical random integer variable. */
int rand_categorical(const vecteur &weights,GIAC_CONTEXT) {
    gen c=rand_uniform(0,_sum(weights,contextptr),contextptr),s=0;
    int ret=0;
    for (const_iterateur it=weights.begin();it!=weights.end();++it,++ret) {
        s+=*it;
        if (is_greater(s,c,contextptr))
            return ret;
    }
    return ret-1;
}

/* TODO */
gen eval_continuous(const gen &g,const vecteur &x,const vecteur &bnds,const vecteur &a,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        if (ckmatrix(g)) {
            const matrice &mat=*g._VECTptr;
            int nr=mat.size(),nc=mat.front()._VECTptr->size();
            matrice ret=zero_mat(nr,nc,contextptr);
            for (int i=0;i<nr;++i) {
                for (int j=0;j<nc;++j) {
                    ret[i]._VECTptr->at(j)=eval_continuous(mat[i][j],x,bnds,a,contextptr);
                }
            }
            return ret;
        } else {
            vecteur ret;
            for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
                ret.push_back(eval_continuous(*it,x,bnds,a,contextptr));
            }
            return change_subtype(ret,g.subtype);
        }
    }
    assert(x.size()==a.size());
    if (is_constant_wrt_vars(g,x,contextptr))
        return g;
    return simp(subst(g,x,a,false,contextptr),contextptr);
    /* TODO: implement multivariable limits */
}

/* Filter critical points in CV. Discard those being complex, those not satisfying
 * the inequality constraints and those for which f (and its extension) is not defined. */
void filter_cpts(vecteur &cv,const vecteur &vars,const vecteur &bnds,const vecteur &ineq,bool open,const gen &f,GIAC_CONTEXT) {
    bound_variables(vars,bnds,open,contextptr);
    for (int j=cv.size();j-->0;) {
        gen val,&cpt=cv[j].type==_VECT && !cv[j]._VECTptr->empty() && cv[j]._VECTptr->front().type==_VECT?cv[j]._VECTptr->front():cv[j];
        vecteur cp=gen2vecteur(cpt);
        cp.resize(vars.size());
        if (has_inf_or_undef(cp) || has_imag(cp,contextptr)) {
            cv.erase(cv.begin()+j);
            continue;
        }
        if (is_constant_wrt_vars(cp,vars,contextptr)) {
            int i=ineq.size();
            for (;i-->0;) {
                gen g=simp(subst(ineq[i],vars,cp,false,contextptr),contextptr);
                if (!(open?is_strictly_positive(-g,contextptr)
                          :((is_approx(g) && is_zero(_pow(makesequence(_abs(g,contextptr),1.2),contextptr),contextptr))
                            || is_positive(-g,contextptr))))
                    break;
            }
            if (i>=0 || is_undef(val=simp(eval_continuous(f,vars,bnds,cp,contextptr),contextptr)) || has_imag(val,contextptr)) {
                cv.erase(cv.begin()+j);
                continue;
            }
        } else if (is_undef(val=simp(subst(f,vars,cp,false,contextptr),contextptr)) ||
                !is_constant_wrt_vars(val,vars,contextptr) || has_imag(val,contextptr)) {
            cv.erase(cv.begin()+j);
            continue;
        }
        cpt=simp(cpt,contextptr);
    }
    purge_variables(vars,contextptr);
}

/* Call _solve without printing any messages.
 * If an error occurs, an empty list is returned. */
vecteur solve_quiet(const gen &e,const gen &x,GIAC_CONTEXT) {
    gen sol;
    try {
        log_output_redirect lor(contextptr);
        sol=_solve(makesequence(e,x),contextptr);
        if (sol.type!=_VECT)
            return vecteur(0);
    } catch (const std::runtime_error &e) {
        return vecteur(0);
    }
    return *sol._VECTptr;
}

int var_index=0;

/* Remove strictly positive or strictly negative factors from G and return the result. */
gen remove_nonzero_factors(const gen &g,GIAC_CONTEXT) {
    gen f(g);
    bool is_neg=false;
    if (f.is_symb_of_sommet(at_neg)) {
        f=f._SYMBptr->feuille;
        is_neg=true;
    }
    if (f.is_symb_of_sommet(at_prod) && f._SYMBptr->feuille.type==_VECT) {
        const vecteur &fv=*f._SYMBptr->feuille._VECTptr;
        gen p(1);
        for (const_iterateur jt=fv.begin();jt!=fv.end();++jt) {
            if (is_const_sign_strict(*jt,true,true,contextptr) || is_const_sign_strict(*jt,false,true,contextptr))
              continue;
            else p=*jt*p;
        }
        f=p;
    } else if (is_const_sign_strict(f,true,true,contextptr) || is_const_sign_strict(f,false,true,contextptr))
        f=1;
    return (is_neg?-1:1)*f;
}

/* Solve a system of equations.
 * This function is based on _solve but handles cases where a variable
 * is found inside trigonometric, hyperbolic or exponential functions. */
vecteur solve2(const vecteur &e_orig,const vecteur &vars_orig,GIAC_CONTEXT) {
    int m=e_orig.size(),n=vars_orig.size(),i=0;
    vecteur e_orig_simp=*expexpand(expand(_pow2exp(e_orig,contextptr),contextptr),contextptr)._VECTptr;
    for (iterateur it=e_orig_simp.begin();it!=e_orig_simp.end();++it) {
        if (is_equal(*it))
            *it=equal2diff(*it);
        gen f=(it->type==_SYMB?_factor(*it,contextptr):*it);
        gen num=remove_nonzero_factors(_numer(f,contextptr),contextptr);
        gen den=remove_nonzero_factors(_denom(f,contextptr),contextptr);
        *it=num/den;
    }
    for (;i<m;++i) {
        if (!is_rational_wrt_vars(e_orig_simp[i],vars_orig))
            break;
    }
    if (n==1 || i==m)
        return solve_quiet(e_orig_simp,vars_orig,contextptr);
    vecteur e(*halftan(_texpand(hyp2exp(e_orig_simp,contextptr),contextptr),contextptr)._VECTptr);
    vecteur lv(*exact(lvar(_evalf(lvar(e),contextptr)),contextptr)._VECTptr);
    vecteur deps(n),depvars(n,gen(0));
    vecteur vars(vars_orig);
    const_iterateur it=lv.begin();
    for (;it!=lv.end();++it) {
        i=0;
        for (;i<n;++i) {
            if (is_undef(vars[i]))
                continue;
            if (*it==(deps[i]=vars[i]) ||
                    *it==(deps[i]=exp(vars[i],contextptr)) ||
                    is_exactly_zero(simp(*it-(deps[i]=tan(vars[i]/2,contextptr)),contextptr))) {
                vars[i]=undef;
                depvars[i]=temp_symb("depvar",i,contextptr);
                break;
            }
        }
        if (i==n)
            break;
    }
    if (it!=lv.end() || contains(depvars,gen(0)))
        return solve_quiet(e_orig_simp,vars_orig,contextptr);
    vecteur e_subs=subst(e,deps,depvars,false,contextptr);
    vecteur sol=solve_quiet(e_subs,depvars,contextptr);
    vecteur ret;
    for (const_iterateur it=sol.begin();it!=sol.end();++it) {
        vecteur r(n);
        i=0;
        for (;i<n;++i) {
            gen c(it->_VECTptr->at(i));
            if (deps[i].type==_IDNT)
                r[i]=c;
            else if (deps[i].is_symb_of_sommet(at_exp) && is_strictly_positive(c,contextptr))
                r[i]=simp(ln(c,contextptr),contextptr);
            else if (deps[i].is_symb_of_sommet(at_tan))
                r[i]=simp(2*atan(c,contextptr),contextptr);
            else
                break;
        }
        if (i==n)
            ret.push_back(r);
    }
    return ret;
}

/* Solve the system of equations e=0 w.r.t. variables in v, including the solutions of a
 * continuous extension too. It is assumed that e is exact. */
vecteur zeros_ext(const vecteur &e,const vecteur &v,const vecteur &bnds,GIAC_CONTEXT) {
    vecteur e_numer=*_apply(makesequence(at_numer,e),contextptr)._VECTptr;
    vecteur e_denom=*_apply(makesequence(at_denom,e),contextptr)._VECTptr;
    vecteur sol=solve2(e_numer,v,contextptr),res;
    for (const_iterateur it=sol.begin();it!=sol.end();++it) {
        bool ok=true;
        if (!is_approx(*it)) {
            for (const_iterateur jt=e.begin();ok && jt!=e.end();++jt) {
                if (!is_constant_wrt_vars(e_denom[jt-e.begin()],v,contextptr) &&
                        is_zero(simp(subst(e_denom[jt-e.begin()],v,*(it->_VECTptr),false,contextptr),contextptr))) {
                    ok=is_zero(simp(eval_continuous(*jt,v,bnds,*(it->_VECTptr),contextptr),contextptr));
                }
            }
        }
        if (ok)
            res.push_back(*it);
    }
    return res;
}

/* Return true if G = (A and B) for some A, B. */
bool is_conjunction(const gen &g) {
    return g.is_symb_of_sommet(at_and) || g.is_symb_of_sommet(at_et);
}

/* Get bounds on variables in VARS using the inequalities INEQ. */
bool get_variable_bounds(const vecteur &vars,vecteur &ineq,vecteur &bnds,GIAC_CONTEXT) {
    gen a,b,c,vmin,vmax;
    bnds=zero_mat(vars.size(),2,contextptr);
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        vecteur &bnd=*bnds[it-vars.begin()]._VECTptr;
        vmin=bnd.front()=minus_inf;
        vmax=bnd.back()=plus_inf;
        vecteur as;
        vector<int> ind;
        for (const_iterateur jt=ineq.begin();jt!=ineq.end();++jt) {
            if ((is_linear_wrt(*jt,*it,a,b,contextptr) || (is_quadratic_wrt(*jt,*it,a,b,c,contextptr) && is_constant_wrt_vars(c,vars,contextptr))) &&
                    is_constant_wrt_vars(a,vars,contextptr) && is_constant_wrt_vars(b,vars,contextptr) && _lname(*jt,contextptr)._VECTptr->size()==1) {
                as.push_back(symb_inferieur_egal(*jt,0));
                ind.push_back(jt-ineq.begin());
            }
        }
        if (!as.empty()) {
            as=solve_quiet(as,*it,contextptr); // solve a system of inequalities in one variable
            if (as.size()==1 && as.front().type!=_IDNT) {
                std::sort(ind.begin(),ind.end());
                for (vector<int>::const_reverse_iterator jt=ind.rbegin();jt!=ind.rend();++jt) {
                    ineq.erase(ineq.begin()+*jt);
                }
                const gen &s=as.front();
                if (s.is_symb_of_sommet(at_inferieur_egal) && s._SYMBptr->feuille._VECTptr->front()==*it)
                    vmax=s._SYMBptr->feuille._VECTptr->back();
                else if (s.is_symb_of_sommet(at_superieur_egal) && s._SYMBptr->feuille._VECTptr->front()==*it)
                    vmin=s._SYMBptr->feuille._VECTptr->back();
                else if (is_conjunction(s) && s._SYMBptr->feuille._VECTptr->size()==2 &&
                        s._SYMBptr->feuille._VECTptr->front().is_symb_of_sommet(at_superieur_egal) &&
                        s._SYMBptr->feuille._VECTptr->front()._SYMBptr->feuille._VECTptr->front()==*it &&
                        s._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_inferieur_egal) &&
                        s._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->front()==*it) {
                    vmin=s._SYMBptr->feuille._VECTptr->front()._SYMBptr->feuille._VECTptr->back();
                    vmax=s._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->back();
                } else if (is_conjunction(s) && s._SYMBptr->feuille._VECTptr->size()==2 &&
                        s._SYMBptr->feuille._VECTptr->front().is_symb_of_sommet(at_inferieur_egal) &&
                        s._SYMBptr->feuille._VECTptr->front()._SYMBptr->feuille._VECTptr->front()==*it &&
                        s._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_superieur_egal) &&
                        s._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->front()==*it) {
                    vmax=s._SYMBptr->feuille._VECTptr->front()._SYMBptr->feuille._VECTptr->back();
                    vmin=s._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->back();
                }
                if (!is_inf(vmin))
                    ineq.push_back(vmin-*it);
                if (!is_inf(vmax))
                    ineq.push_back(*it-vmax);
            } else if (as.empty())
                return false; // infeasible
        } else continue;
        bnd.front()=vmin;
        bnd.back()=vmax;
    }
    return true;
}

/*
 *
 ******* IPDIFF CLASS IMPLEMENTATION *******
 *
 */

ipdiff::ipdiff(const gen &f_orig,const vecteur &g_orig,const vecteur &vars_orig,GIAC_CONTEXT) {
    ctx=contextptr;
    f=f_orig;
    g=g_orig;
    vars=vars_orig;
    ord=0;
    nconstr=g.size();
    nvars=vars.size()-nconstr;
    assert(nvars>0);
    pdv[ivector(nvars,0)]=f; // make the zeroth order derivative initially available
}
void ipdiff::ipartition(int m,int n,ivectors &c,const ivector &p) {
    for (int i=0;i<n;++i) {
        if (!p.empty() && p[i]!=0)
            continue;
        ivector r;
        if (p.empty())
            r.resize(n,0);
        else r=p;
        for (int j=0;j<m;++j) {
            ++r[i];
            int s=sum_ivector(r);
            if (s==m && find(c.begin(),c.end(),r)==c.end())
                c.push_back(r);
            else if (s<m)
                ipartition(m,n,c,r);
            else break;
        }
    }
}
ipdiff::diffterms ipdiff::derive_diffterms(const diffterms &terms,ivector &sig) {
    while (!sig.empty() && sig.back()==0) {
        sig.pop_back();
    }
    if (sig.empty())
        return terms;
    int k=sig.size()-1,p;
    diffterms tv;
    ivector u(nvars+1,0);
    for (diffterms::const_iterator it=terms.begin();it!=terms.end();++it) {
        int c=it->second;
        diffterm t(it->first);
        const ivector_map &h_orig=it->first.second;
        ++t.first.at(k);
        tv[t]+=c;
        --t.first.at(k);
        ivector_map h(h_orig);
        for (ivector_map::const_iterator jt=h_orig.begin();jt!=h_orig.end();++jt) {
            ivector v=jt->first;
            if ((p=jt->second)==0)
                continue;
            if (p==1)
                h.erase(h.find(v));
            else
                --h[v];
            ++v[k];
            ++h[v];
            t.second=h;
            tv[t]+=c*p;
            --h[v];
            --v[k];
            ++h[v];
        }
        t.second=h_orig;
        for (int i=0;i<nconstr;++i) {
            ++t.first.at(nvars+i);
            u[k]=1;
            u.back()=i;
            ++t.second[u];
            tv[t]+=c;
            --t.first.at(nvars+i);
            --t.second[u];
            u[k]=0;
        }
    }
    --sig.back();
    return derive_diffterms(tv,sig);
}
const gen &ipdiff::get_pd(const pd_map &pds,const ivector &sig) const {
    try {
        return pds.at(sig);
    }
    catch (out_of_range &e) {
        return undef;
    }
}
const gen &ipdiff::differentiate(const gen &e,pd_map &pds,const ivector &sig) {
    const gen &pd=get_pd(pds,sig);
    if (!is_undef(pd))
        return pd;
    vecteur v(1,e);
    bool do_derive=false;
    assert(vars.size()<=sig.size());
    for (int i=0;i<int(vars.size());++i) {
        if (sig[i]>0) {
            v=mergevecteur(v,vecteur(sig[i],vars[i]));
            do_derive=true;
        }
    }
    if (do_derive)
        return pds[sig]=_derive(_feuille(v,ctx),ctx);
    return e;
}
void ipdiff::compute_h(const vector<diffterms> &grv,int order) {
    if (g.empty())
        return;
    ivectors hsigv;
    matrice A;
    vecteur b(g.size()*grv.size(),gen(0));
    gen t;
    int grv_sz=grv.size();
    for (int i=0;i<nconstr;++i) {
        for (int j=0;j<grv_sz;++j) {
            vecteur eq(g.size()*grv_sz,gen(0));
            const diffterms &grvj=grv[j];
            for (diffterms::const_iterator it=grvj.begin();it!=grvj.end();++it) {
                ivector sig(it->first.first),hsig;
                sig.push_back(i);
                t=gen(it->second)*differentiate(g[i],pdg,sig);
                for (ivector_map::const_iterator ht=it->first.second.begin();ht!=it->first.second.end();++ht) {
                    if (ht->second==0)
                        continue;
                    const ivector &sigh=ht->first;
                    if (sum_ivector(sigh,true)<order) {
                        gen h(get_pd(pdh,sigh));
                        assert(!is_undef(h));
                        t=t*pow(h,ht->second);
                    }
                    else {
                        assert(ht->second==1);
                        hsig=sigh;
                    }
                }
                if (hsig.empty())
                    b[grv_sz*i+j]-=t;
                else {
                    int k=0,hsigv_sz=hsigv.size();
                    for (;k<hsigv_sz;++k) {
                        if (hsigv[k]==hsig)
                            break;
                    }
                    eq[k]+=t;
                    if (k==hsigv_sz)
                        hsigv.push_back(hsig);
                }
            }
            A.push_back(*simp(eq,ctx)._VECTptr);
        }
    }
    matrice B;
    B.push_back(*simp(b,ctx)._VECTptr);
    matrice invA=*_inv(A,ctx)._VECTptr;
    vecteur sol(*mtran(mmult(invA,mtran(B))).front()._VECTptr);
    for (int i=0;i<int(sol.size());++i) {
        pdh[hsigv[i]]=simp(sol[i],ctx);
    }
}
void ipdiff::find_nearest_terms(const ivector &sig,diffterms &match,ivector &excess) {
    excess=sig;
    int i;
    for (map<ivector,diffterms>::const_iterator it=cterms.begin();it!=cterms.end();++it) {
        ivector ex(nvars,0);
        for (i=0;i<nvars;++i) {
            if ((ex[i]=sig[i]-it->first.at(i))<0)
                break;
        }
        if (i<nvars)
            continue;
        if (sum_ivector(ex)<sum_ivector(excess)) {
            excess=ex;
            match=it->second;
        }
    }
}
void ipdiff::raise_order(int order) {
    if (g.empty())
        return;
    ivectors c;
    ivector excess,init_f(nvars+nconstr,0);
    diffterm init_term;
    init_term.first=init_f;
    diffterms init_terms;
    init_terms[init_term]=1;
    vector<diffterms> grv;
    for (int k=ord+1;k<=order;++k) {
        grv.clear();
        c.clear();
        ipartition(k,nvars,c);
        for (ivectors::const_iterator it=c.begin();it!=c.end();++it) {
            diffterms terms=init_terms;
            find_nearest_terms(*it,terms,excess);
            if (sum_ivector(excess)>0) {
                terms=derive_diffterms(terms,excess);
                cterms[*it]=terms;
            }
            grv.push_back(terms);
        }
        compute_h(grv,k);
    }
    ord=order;
}
void ipdiff::compute_pd(int order,const ivector &sig) {
    gen pd;
    ivectors c;
    ipartition(order,nvars,c);
    for (ivectors::const_iterator ct=c.begin();ct!=c.end();++ct) {
        if (!sig.empty() && sig!=*ct)
            continue;
        if (g.empty()) {
            differentiate(f,pdv,sig);
            continue;
        }
        diffterms &terms=cterms[*ct];
        pd=gen(0);
        for (diffterms::const_iterator it=terms.begin();it!=terms.end();++it) {
            ivector sig(it->first.first);
            gen t(gen(it->second)*differentiate(f,pdf,sig));
            if (!is_exactly_zero(t)) {
                for (ivector_map::const_iterator jt=it->first.second.begin();jt!=it->first.second.end();++jt) {
                    if (jt->second==0)
                        continue;
                    gen h(get_pd(pdh,jt->first));
                    assert(!is_undef(h));
                    t=t*pow(h,jt->second);
                }
                pd+=t;
            }
        }
        pdv[*ct]=simp(pd,ctx);
    }
}
void ipdiff::gradient(vecteur &res) {
    if (nconstr==0)
        res=*_grad(makesequence(f,vars),ctx)._VECTptr;
    else {
        res.resize(nvars);
        ivector sig(nvars,0);
        if (ord<1) {
            raise_order(1);
            compute_pd(1);
        }
        for (int i=0;i<nvars;++i) {
            sig[i]=1;
            res[i]=derivative(sig);
            sig[i]=0;
        }
    }
}
void ipdiff::hessian(matrice &res) {
    if (nconstr==0)
        res=*_hessian(makesequence(f,vars),ctx)._VECTptr;
    else {
        res.clear();
        ivector sig(nvars,0);
        if (ord<2) {
            raise_order(2);
            compute_pd(2);
        }
        for (int i=0;i<nvars;++i) {
            vecteur r(nvars);
            ++sig[i];
            for (int j=0;j<nvars;++j) {
                ++sig[j];
                r[j]=derivative(sig);
                --sig[j];
            }
            res.push_back(r);
            --sig[i];
        }
    }
}
const gen &ipdiff::derivative(const ivector &sig) {
    if (nconstr==0)
        return differentiate(f,pdf,sig);
    int k=sum_ivector(sig); // the order of the derivative
    if (k>ord) {
        raise_order(k);
        compute_pd(k,sig);
    }
    return get_pd(pdv,sig);
}
const gen &ipdiff::derivative(const vecteur &dvars) {
    ivector sig(nvars,0);
    int j;
    for (const_iterateur it=dvars.begin();it!=dvars.end();++it) {
        if ((j=indexof(*it,vars))<0)
            return undef;
        ++sig[j];
    }
    return derivative(sig);
}
void ipdiff::partial_derivatives(int order,pd_map &pdmap) {
    if (nconstr>0 && ord<order) {
        raise_order(order);
        compute_pd(order);
    }
    ivectors c;
    ipartition(order,nvars,c);
    for (ivectors::const_iterator it=c.begin();it!=c.end();++it) {
        pdmap[*it]=derivative(*it);
    }
}
gen ipdiff::taylor_term(const vecteur &a,int k,bool shift) {
    assert(k>=0);
    if (k==0)
        return subst(f,vars,a,false,ctx);
    ivectors sigv;
    ipartition(k,nvars,sigv);
    gen term(0);
    if (nconstr>0) while (k>ord) {
        raise_order(ord+1);
        compute_pd(ord);
    }
    for (ivectors::const_iterator it=sigv.begin();it!=sigv.end();++it) {
        gen pd;
        if (g.empty()) {
            vecteur args(1,f);
            for (int i=0;i<nvars;++i) {
                for (int j=0;j<it->at(i);++j) {
                    args.push_back(vars[i]);
                }
            }
            pd=_derive(_feuille(args,ctx),ctx);
        }
        else
            pd=derivative(*it);
        pd=subst(pd,vars,a,false,ctx);
        for (int i=0;i<nvars;++i) {
            int ki=it->at(i);
            if (ki==0)
                continue;
            pd=pd*(shift?pow(vars[i]-a[i],ki):pow(vars[i],ki))/factorial(ki);
        }
        term+=pd;
    }
    return term;
}
gen ipdiff::taylor(const vecteur &a,int order) {
    assert(order>=0);
    gen T(0);
    for (int k=0;k<=order;++k) {
        T+=taylor_term(a,k);
    }
    return T;
}
/* END OF IPDIFF CLASS */

void vars_arrangements(matrice J,ipdiff::ivectors &arrs,GIAC_CONTEXT) {
    int m=J.size(),n=J.front()._VECTptr->size();
    assert(n<=32 && m<n);
    matrice tJ(mtran(J));
    ulong N=std::pow(2,n);
    vector<ulong> sets(comb(n,m).val);
    int i=0;
    for (ulong k=1;k<N;++k) {
        bitset<32> b(k);
        if (b.count()==(size_t)m)
            sets[i++]=k;
    }
    matrice S;
    ipdiff::ivector arr(n);
    for (vector<ulong>::const_iterator it=sets.begin();it!=sets.end();++it) {
        for (i=0;i<n;++i) arr[i]=i;
        N=std::pow(2,n);
        for (i=n;i-->0;) {
            N/=2;
            if ((*it & N)!=0) {
                arr.erase(arr.begin()+i);
                arr.push_back(i);
            }
        }
        S.clear();
        for (ipdiff::ivector::const_iterator it=arr.end()-m;it!=arr.end();++it) {
            S.push_back(tJ[*it]);
        }
        if (!is_zero(mdet(S,contextptr)))
            arrs.push_back(arr);
    }
}

/* Computes the Jacobian of G w.r.t. VARS, returns false on failure. */
bool jacobian(const vecteur &g,const vecteur &vars,matrice &J,GIAC_CONTEXT) {
    J.resize(g.size());
    gen gr;
    for (const_iterateur it=g.begin();it!=g.end();++it) {
        try {
            gr=_grad(makesequence(*it,vars),contextptr);
        } catch (const std::runtime_error &e) {
            return false;
        }
        if (gr.type==_VECT && gr._VECTptr->size()==vars.size() && !has_inf_or_undef(gr) && !has_diff(gr,vars))
            J[it-g.begin()]=gr;
        else return false;
    }
    return true;
}

bool ck_jacobian(vecteur &g,vecteur &vars,GIAC_CONTEXT) {
    assert(!g.empty());
    matrice J;
    if (!jacobian(g,vars,J,contextptr))
        return false;
    int m=g.size();
    int n=vars.size()-m;
    if (_rank(J,contextptr).val<m)
        return false;
    J=mtran(J);
    J.erase(J.begin(),J.begin()+n);
    return !is_zero(mdet(J,contextptr),contextptr);
}

/*
 * 'implicitdiff' differentiates function(s) defined by equation(s) or a
 * function f(x1,x2,...,xn,y1,y2,...,ym) where y1,...,ym are functions of
 * x1,x2,...,xn defined by m equality constraints.
 *
 * Usage
 * ^^^^^
 *      implicitdiff(f,constr,depvars,diffvars)
 *      implicitdiff(f,constr,vars,order=<posint>,[P])
 *      implicitdiff(constr,[depvars],y,diffvars)
 *
 * Parameters
 * ^^^^^^^^^^
 *      - f         : expression
 *      - constr    : (list of) equation(s)
 *      - depvars   : (list of) dependent variable(s), each of them given
 *                    either as a symbol, e.g. y, or a function, e.g. y(x,z)
 *      - diffvars  : sequence of variables w.r.t. which the differentiation
 *                    will be carried out
 *      - vars      : list all variables on which f depends such that
 *                  : dependent variables come after independent ones
 *      - P         : (list of) coordinate(s) to compute derivatives at
 *      - y         : (list of) dependent variable(s) to differentiate w.r.t.
 *                    diffvars, each of them given as a symbol
 *
 * The return value is partial derivative specified by diffvars. If
 * 'order=m' is given as the fourth argument, all partial derivatives of
 * order m will be computed and returned as vector for m=1, matrix for m=2 or
 * table for m>2. The first two cases produce gradient and hessian of f,
 * respectively. For m>2, the partial derivative
 * pd=d^m(f)/(d^k1(x1)*d^k2(x2)*...*d^kn(xn)) is saved under key [k1,k2,...kn].
 * If P is specified, pd(P) is saved.
 *
 * Examples
 * ^^^^^^^^
 * implicitdiff(x^2*y+y^2=1,y,x)
 *      >> -2*x*y/(x^2+2*y)
 * implicitdiff(R=P*V/T,P,T)
 *      >> P/T
 * implicitdiff([x^2+y=z,x+y*z=1],[y(x),z(x)],y,x)
 *      >> (-2*x*y-1)/(y+z)
 * implicitdiff([x^2+y=z,x+y*z=1],[y(x),z(x)],[y,z],x)
 *      >> [(-2*x*y-1)/(y+z),(2*x*z-1)/(y+z)]
 * implicitdiff(y=x^2/z,y,x)
 *      >> 2x/z
 * implicitdiff(y=x^2/z,y,z)
 *      >> -x^2/z^2
 * implicitdiff(y^3+x^2=1,y,x)
 *      >> -2*x/(3*y^2)
 * implicitdiff(y^3+x^2=1,y,x,x)
 *      >> (-8*x^2-6*y^3)/(9*y^5)x+3y-z,2x^2+y^2=z,[x,y,z]
 * implicitdiff(a*x^3*y-2y/z=z^2,y(x,z),x)
 *      >> -3*a*x^2*y*z/(a*x^3*z-2)
 * implicitdiff(a*x^3*y-2y/z=z^2,y(x,z),x,z)
 *      >> (12*a*x^2*y-6*a*x^2*z^3)/(a^2*x^6*z^2-4*a*x^3*z+4)
 * implicitdiff([-2x*z+y^2=1,x^2-exp(x*z)=y],[y(x),z(x)],y,x)
 *      >> 2*x/(y*exp(x*z)+1)
 * implicitdiff([-2x*z+y^2=1,x^2-exp(x*z)=y],[y(x),z(x)],[y,z],x)
 *      >> [2*x/(y*exp(x*z)+1),(2*x*y-y*z*exp(x*z)-z)/(x*y*exp(x*z)+x)]
 * implicitdiff([a*sin(u*v)+b*cos(w*x)=c,u+v+w+x=z,u*v+w*x=z],[u(x,z),v(x,z),w(x,z)],u,z)
 *      >> (a*u*x*cos(u*v)-a*u*cos(u*v)+b*u*x*sin(w*x)-b*x*sin(w*x))/
 *         (a*u*x*cos(u*v)-a*v*x*cos(u*v)+b*u*x*sin(w*x)-b*v*x*sin(w*x))
 * implicitdiff(x*y,-2x^3+15x^2*y+11y^3-24y=0,y(x),x$2)
 *      >> (162*x^5*y+1320*x^4*y^2-320*x^4-3300*x^3*y^3+800*x^3*y+968*x^2*y^4-1408*x^2*y^2+
 *          512*x^2-3630*x*y^5+5280*x*y^3-1920*x*y)/(125*x^6+825*x^4*y^2-600*x^4+
 *          1815*x^2*y^4-2640*x^2*y^2+960*x^2+1331*y^6-2904*y^4+2112*y^2-512)
 * implicitdiff((x-u)^2+(y-v)^2,[x^2/4+y^2/9=1,(u-3)^2+(v+5)^2=1],[v(u,x),y(u,x)],u,x)
 *      >> (-9*u*x-4*v*y+27*x-20*y)/(2*v*y+10*y)
 * implicitdiff(x*y*z,-2x^3+15x^2*y+11y^3-24y=0,[x,z,y],order=1)
 *      >> [(2*x^3*z-5*x^2*y*z+11*y^3*z-8*y*z)/(5*x^2+11*y^2-8),x*y]
 * implicitdiff(x*y*z,-2x^3+15x^2*y+11y^3-24y=0,[x,z,y],order=2,[1,-1,0])
 *      >> [[64/9,-2/3],[-2/3,0]]
 * pd:=implicitdiff(x*y*z,-2x^3+15x^2*y+11y^3-24y=0,[x,z,y],order=4,[0,z,0]);pd[4,0,0]
 *      >> -2*z
 */
gen _implicitdiff(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()<2)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    const gen &f=gv[0];
    if (int(gv.size())<3)
        return generr("Too few arguments");
    int ci=gv[0].type!=_VECT && !is_equal(gv[0])?1:0;
    vecteur freevars,depvars,diffdepvars;
    gen_map diffvars;
    // get the constraints as a list of vanishing expressions
    vecteur constr=gen2vecteur(gv[ci]);
    for (int i=0;i<int(constr.size());++i) {
        if (is_equal(constr[i]))
            constr[i]=equal2diff(constr[i]);
    }
    int m=constr.size();
    int dvi=3;
    if (ci==0) {
        if (gv[ci+1].type==_VECT)
            diffdepvars=gen2vecteur(gv[ci+2]);
        else dvi=2;
    }
    bool compute_all=false;
    int order=0;
    if (ci==1 && is_equal(gv[dvi])) {
        vecteur &v=*gv[dvi]._SYMBptr->feuille._VECTptr;
        if (v.front()!=at_order || !v.back().is_integer() || v.back().val<=0)
            return generrtype("Expected order=<posint>");
        order=v.back().val;
        compute_all=true;
    }
    // get dependency specification
    vecteur deplist=gen2vecteur(gv[ci+1]);
    if (compute_all) {
        // vars must be specified as x1,x2,...,xn,y1,y2,...,ym
        int nd=deplist.size();
        if (nd<=m)
            return generrdim("Too few variables");
        for (int i=0;i<nd;++i) {
            if (i<nd-m)
                freevars.push_back(deplist[i]);
            else depvars.push_back(deplist[i]);
        }
    } else {
        // get (in)dependent variables
        for (const_iterateur it=deplist.begin();it!=deplist.end();++it) {
            if (it->type==_IDNT)
                depvars.push_back(*it);
            else if (it->is_symb_of_sommet(at_of)) {
                vecteur fe(*it->_SYMBptr->feuille._VECTptr);
                depvars.push_back(fe.front());
                if (fe.back().type==_VECT) {
                    for (int i=0;i<int(fe.back()._VECTptr->size());++i) {
                        gen &x=fe.back()._VECTptr->at(i);
                        if (!contains(freevars,x))
                            freevars.push_back(x);
                    }
                } else freevars.push_back(fe.back());
            } else return generrtype("Invalid variable specification");
        }
        // get diffvars
        for (const_iterateur it=gv.begin()+dvi;it!=gv.end();++it) {
            gen v(eval(*it,contextptr));
            gen x;
            if (v.type==_IDNT)
                diffvars[(x=v)]+=1;
            else if (v.type==_VECT && v.subtype==_SEQ__VECT)
                diffvars[(x=v._VECTptr->front())]+=v._VECTptr->size();
            else return generrtype("Invalid variable specification");
            if (!contains(freevars,x))
                freevars.push_back(x);
        }
    }
    int n=freevars.size();  // number of independent variables
    if (m!=int(depvars.size()))
        return generrdim("Invalid number of dependent variables");
    vecteur vars(mergevecteur(freevars,depvars));  // list of all variables
    // check whether the conditions of implicit function theorem hold
    if (!ck_jacobian(constr,vars,contextptr))
        return generr("Conditions of implicit function theorem are violated");
    // build partial derivative specification 'sig'
    ipdiff::ivector sig(n,0); // sig[i]=k means: derive k times with respect to ith independent variable
    ipdiff ipd(f,constr,vars,contextptr);
    if (compute_all) {
        vecteur pt(0);
        if (int(gv.size())>4) {
            pt=gen2vecteur(gv[4]);
            if (int(pt.size())!=n+m)
                return generrdim("The given point does not match the number of variables");
        }
        ipdiff::pd_map pdv;
        ipd.partial_derivatives(order,pdv);
        if (order==1) {
            vecteur gr;
            ipd.gradient(gr);
            return pt.empty()?gr:simp(subst(gr,vars,pt,false,contextptr),contextptr);
        }
        else if (order==2) {
            matrice hess;
            ipd.hessian(hess);
            return pt.empty()?hess:simp(subst(hess,vars,pt,false,contextptr),contextptr);
        }
        else {
            ipdiff::ivectors c;
            ipdiff::ipartition(order,n,c);
            gen_map ret_pdv;
            for (ipdiff::ivectors::const_iterator it=c.begin();it!=c.end();++it) {
                vecteur v;
                for (int i=0;i<n;++i) {
                    v.push_back(gen(it->at(i)));
                }
                ret_pdv[v]=pt.empty()?pdv[*it]:simp(subst(pdv[*it],vars,pt,false,contextptr),contextptr);
            }
            return ret_pdv;
        }
    }
    for (gen_map::const_iterator it=diffvars.begin();it!=diffvars.end();++it) {
        int i=0;
        for (;i<n;++i) {
            if (it->first==freevars[i]) {
                sig[i]=it->second.val;
                break;
            }
        }
        assert(i<n);
    }
    // compute the partial derivative specified by 'sig'
    order=ipdiff::sum_ivector(sig);
    if (ci==1)
        return simp(ipd.derivative(sig),contextptr);
    vecteur ret;
    if (diffdepvars.empty()) {
        assert(m==1);
        diffdepvars=vecteur(1,depvars.front());
    }
    for (const_iterateur it=diffdepvars.begin();it!=diffdepvars.end();++it) {
        if (!contains(depvars,*it)) {
            // variable *it is not in depvars, so it's treated as a constant
            ret.push_back(0);
            continue;
        }
        ipdiff tmp(*it,constr,vars,contextptr);
        ret.push_back(simp(tmp.derivative(sig),contextptr));
    }
    return ret.size()==1?ret.front():ret;
}
static const char _implicitdiff_s []="implicitdiff";
static define_unary_function_eval (__implicitdiff,&_implicitdiff,_implicitdiff_s);
define_unary_function_ptr5(at_implicitdiff,alias_at_implicitdiff,&__implicitdiff,0,true)

iterateur find_cpt(vecteur &cpts,const vecteur &cand,GIAC_CONTEXT) {
    iterateur ret;
    for (ret=cpts.begin();ret!=cpts.end();++ret) {
        if (ret->type!=_VECT || ret->_VECTptr->empty())
            continue;
        if (ret->_VECTptr->front().type==_VECT) {
            const vecteur &c=*ret->_VECTptr->front()._VECTptr;
            if (c.size()!=cand.size())
                continue;
            const_iterateur it;
            for (it=c.begin();it!=c.end();++it) {
                if (!is_zero(simp(*it-cand[it-c.begin()],contextptr),contextptr))
                    break;
            }
            if (it==c.end())
                break;
        }
    }
    return ret;
}

/* classification using the bordered Hessian and Theorem 1 of David Spring (1985) */
int critical_point_class(const matrice &hess,int n,int m,GIAC_CONTEXT) {
    vecteur s;
    int i,j,k;
    for (k=1;k<=n;++k) {
        matrice M;
        for (i=0;i<2*m+k;++i) {
            const vecteur &row=*hess[i]._VECTptr;
            M.push_back(vecteur(row.begin(),row.begin()+2*m+k));
        }
        s.push_back(simp(pow(gen(-1),m)*mdet(M,contextptr),contextptr));
    }
    if (contains(s,gen(0)))
        return _CPCLASS_UNDECIDED; // paranoid check, a counterexample exists for Theorem 1 of D.S.
    for (i=n-1;i>=0 && is_zero(s[i],contextptr);--i);
    if (i<0) return _CPCLASS_UNDECIDED; // the sequence is trivial
    for (j=0;j<=i && is_strictly_positive(s[j],contextptr);++j);
    if (j>i) return i==n-1?_CPCLASS_MIN:_CPCLASS_UNDECIDED;
    for (j=0;j<=i && is_strictly_positive(pow(gen(-1),j+1)*s[j],contextptr);++j);
    if (j>i) return i==n-1?_CPCLASS_MAX:_CPCLASS_UNDECIDED;
    return _CPCLASS_SADDLE;
}

/* Find those variables in VARS which appear in E and store them in X. */
void find_vars(const gen &e,const vecteur &vars,vecteur &x,GIAC_CONTEXT) {
    vecteur s=*_lname(e,contextptr)._VECTptr;
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        if (contains(s,*it))
            x.push_back(*it);
    }
}

bool get_parameter_assumptions(const vecteur &parm,gen_map &parm_asmp,GIAC_CONTEXT) {
    for (const_iterateur it=parm.begin();it!=parm.end();++it) {
        int dom;
        matrice intervals;
        vecteur excluded;
        get_assumptions(*it,dom,intervals,excluded,contextptr);
        if (dom==_ZINT || dom==_FRAC || dom==_CPLX)
            return false; // domain not supported, only real parameters are valid
        parm_asmp[*it]=makevecteur(intervals,excluded);
    }
    return true;
}

bool get_additional_parameter_assumptions(const vecteur &cpt,const vecteur &vars,const vecteur &ineq,bool open,gen_map &conds,GIAC_CONTEXT) {
    vecteur rest_conds,parm=*_lname(cpt,contextptr)._VECTptr;
    conds.clear();
    for (int i=parm.size();i-->0;) {
        if (contains(vars,parm[i]))
            parm.erase(parm.begin()+i);
    }
    if (parm.empty())
        return true;
    for (const_iterateur it=ineq.begin();it!=ineq.end();++it) {
        gen g=simp(subst(*it,vars,cpt,false,contextptr),contextptr);
        vecteur s;
        find_vars(g,parm,s,contextptr);
        if (s.empty())
            continue;
        gen h=open?symb_inferieur_strict(g,0):symb_inferieur_egal(g,0);
        if (s.size()==1) {
            gen &c=conds[s.front()];
            if (c.type==_VECT)
                c._VECTptr->push_back(h);
            else c=vecteur(1,h);
        } else rest_conds.push_back(h);
    }
    for (gen_map::iterator it=conds.begin();it!=conds.end();++it) {
        assert(it->second.type==_VECT);
        vecteur sol=solve_quiet(it->second,it->first,contextptr);
        it->second._VECTptr->clear();
        if (sol.empty()) {
            *logptr(contextptr) << gettext("Warning") << ": " << gettext("failed to solve ") << it->second << gettext(" for ")
                                << it->first << ", " << gettext("parameter seems to have no value") << "\n";
            return false; // infeasible
        } else if (sol.size()==1 && sol.front()==it->first)
            ; // no assumption needed
        else {
            for (const_iterateur jt=sol.begin();jt!=sol.end();++jt) {
                if (jt->type==_SYMB && (which_ineq(*jt)!=0 || is_conjunction(*jt))) {
                    gen e=_evalb(*jt,contextptr);
                    if (e.is_integer() && e.subtype==_INT_BOOLEAN) {
                        if (e.val==0)
                            return false; // infeasible
                        continue; // no assumption needed
                    }
                    it->second._VECTptr->push_back(*jt);
                }
            }
        }
    }
    if (!rest_conds.empty())
        conds[_lname(rest_conds,contextptr)]=rest_conds;
    return true;
}

void restore_parameter_assumptions(const gen_map &parm_asmp,GIAC_CONTEXT) {
    for (gen_map::const_iterator jt=parm_asmp.begin();jt!=parm_asmp.end();++jt) {
        set_assumptions(jt->first,*jt->second._VECTptr->front()._VECTptr,*jt->second._VECTptr->back()._VECTptr,false,contextptr);
    }
}

void find_local_extrema(vecteur &cpts,const gen &f,const vecteur &g,const vecteur &vars,const ipdiff::ivector &arr,const vecteur &ineq,
                        const vecteur &bnds,bool open,const vecteur &initial,const gen_map &parm_asmp,int ord,bool approx_hompol,bool cont_ext,GIAC_CONTEXT) {
    int nv=vars.size(),m=g.size(),n=nv-m,cls;
    if (ord==0 && m>0) { // apply the method of Lagrange
        gen L(f);
        vecteur multipliers(m),allinitial;
        if (!initial.empty())
            allinitial=mergevecteur(vecteur(m,0),initial);
        for (int i=m;i-->0;) {
            multipliers[i]=temp_symb("lambda",++var_index,contextptr);
            _purge(multipliers[i],contextptr);
            L+=-multipliers[i]*g[i];
        }
        vecteur allvars=mergevecteur(vars,multipliers),allbnds=bnds;
        for (int i=multipliers.size();i-->0;) allbnds.push_back(makevecteur(minus_inf,plus_inf));
        matrice cv,bhess;
        gen tmpgr=_grad(makesequence(L,allvars),contextptr);
        if (tmpgr.type==_VECT) {
            vecteur &gr=*tmpgr._VECTptr;
            if (allinitial.empty()) {
                bound_variables(vars,bnds,open,contextptr);
                cv=cont_ext?zeros_ext(gr,allvars,allbnds,contextptr):solve2(gr,allvars,contextptr);
                purge_variables(vars,contextptr);
            } else {
                gen tmpfsol=_fsolve(makesequence(gr,allvars,allinitial),contextptr);
                if (tmpfsol.type==_VECT) {
                    vecteur &fsol=*tmpfsol._VECTptr;
                    if (!fsol.empty())
                        cv.push_back(fsol);
                }
            }
            filter_cpts(cv,vars,bnds,ineq,open,f,contextptr);
            for (iterateur it=cv.begin();it!=cv.end();++it) {
                *it=mergevecteur(vecteur(it->_VECTptr->begin()+nv,it->_VECTptr->end()),
                                 vecteur(it->_VECTptr->begin(),it->_VECTptr->begin()+nv));
            }
            allvars=mergevecteur(vecteur(allvars.begin()+nv,allvars.end()),vecteur(allvars.begin(),allvars.begin()+nv));
            allbnds=mergevecteur(vecteur(allbnds.begin()+nv,allbnds.end()),vecteur(allbnds.begin(),allbnds.begin()+nv));
            if (!cv.empty()) {
                gen tmphess=_hessian(makesequence(L,allvars),contextptr);
                if (ckmatrix(tmphess))
                    bhess=*tmphess._VECTptr; // bordered Hessian
            }
            gen s;
            for (const_iterateur it=cv.begin();it!=cv.end();++it) {
                vecteur cpt=vecteur(it->_VECTptr->begin()+m,it->_VECTptr->end());
                gen_map asmp;
                if (!get_additional_parameter_assumptions(cpt,vars,ineq,true,asmp,contextptr))
                    continue;
                for (gen_map::const_iterator jt=asmp.begin();jt!=asmp.end();++jt) {
                    set_assumptions(jt->first,*jt->second._VECTptr,vecteur(0),true,contextptr);
                }
                cls=_CPCLASS_UNDECIDED;
                if (!bhess.empty())
                    cls=critical_point_class(*eval_continuous(bhess,allvars,allbnds,*it->_VECTptr,contextptr)._VECTptr,n,m,contextptr);
                iterateur cit=find_cpt(cpts,cpt,contextptr);
                if (cit!=cpts.end()) {
                    if (cls!=_CPCLASS_UNDECIDED)
                        cit->_VECTptr->at(1)=cls;
                } else cpts.push_back(makevecteur(cpt,cls,asmp));
                restore_parameter_assumptions(parm_asmp,contextptr);
            }
        }
    } else if (ord>0) { // use implicit differentiation instead of Lagrange multipliers
        vecteur gr,taylor_terms,a(nv),cpt_arr(nv);
        ipdiff ipd(f,g,vars,contextptr);
        ipd.gradient(gr);
        matrice cv,hess;
        vecteur eqv=mergevecteur(gr,g);
        if (initial.empty()) {
            bound_variables(vars,bnds,open,contextptr);
            cv=cont_ext?zeros_ext(eqv,vars,bnds,contextptr):solve2(eqv,vars,contextptr);
            purge_variables(vars,contextptr);
        } else {
            gen tmpfsol=_fsolve(makesequence(eqv,vars,initial),contextptr);
            if (tmpfsol.type==_VECT) {
                vecteur &fsol=*tmpfsol._VECTptr;
                if (!fsol.empty())
                    cv.push_back(fsol);
            }
        }
        filter_cpts(cv,vars,bnds,ineq,open,f,contextptr);
        if (!cv.empty()) {
            if (nv==1) {
                gen d;
                for (const_iterateur it=cv.begin();it!=cv.end();++it) {
                    assert(it->type==_VECT);
                    gen_map asmp;
                    if (!get_additional_parameter_assumptions(*it->_VECTptr,vars,ineq,true,asmp,contextptr))
                        continue;
                    cls=_CPCLASS_UNDECIDED;
                    for (int k=2;k<=ord;++k) {
                        d=simp(eval_continuous(derive(f,vars.front(),k,contextptr),vars,bnds,*(it->_VECTptr),contextptr),contextptr);
                        if (is_zero(d,contextptr))
                            continue;
                        if ((k%2)!=0)
                            cls=_CPCLASS_SADDLE;
                        else cls=is_strictly_positive(d,contextptr)?_CPCLASS_MIN:_CPCLASS_MAX;
                        break;
                    }
                    cpts.push_back(makevecteur(it->_VECTptr->front(),cls,asmp));
                }
            } else {
                vecteur fvars(vars),ip;
                fvars.resize(n);
                gen l=temp_symb("lambda",++var_index,contextptr);
                bool approx_hp;
                gen pmin,pmax,sp(-1);
                for (int j=0;j<n;++j) {
                    sp+=pow(vars[j],2);
                    ip.push_back(sqrt(fraction(1,n),contextptr));
                }
                if (ord>1)
                    ipd.hessian(hess);
                for (int i=0;i<nv;++i) {
                    a[i]=temp_symb("a",i,contextptr);
                }
                for (const_iterateur it=cv.begin();it!=cv.end();++it) {
                    const vecteur &cpt=*it->_VECTptr;
                    for (int j=0;j<nv;++j) {
                        cpt_arr[arr[j]]=it->_VECTptr->at(j);
                    }
                    gen_map asmp;
                    if (!get_additional_parameter_assumptions(cpt_arr,vars,ineq,true,asmp,contextptr))
                        continue;
                    iterateur ct=find_cpt(cpts,cpt_arr,contextptr);
                    if (ct==cpts.end()) {
                        cpts.push_back(makevecteur(cpt_arr,0,asmp));
                        ct=cpts.begin()+cpts.size()-1;
                    }
                    gen &cpt_class=ct->_VECTptr->at(1);
                    if (ord==1 || !is_exactly_zero(cpt_class))
                        continue;
                    cls=_CPCLASS_UNDECIDED;
                    if (ord>=2) {
                        for (gen_map::const_iterator jt=asmp.begin();jt!=asmp.end();++jt) {
                            set_assumptions(jt->first,*jt->second._VECTptr,vecteur(0),true,contextptr);
                        }
                        for (int k=2;k<=ord;++k) {
                            if (int(taylor_terms.size())<k-1)
                                taylor_terms.push_back(ipd.taylor_term(a,k,false));
                            if (is_zero(expand(taylor_terms[k-2],contextptr),contextptr))
                                break;
                            vecteur csts=*_lname(taylor_terms[k-2],contextptr)._VECTptr;
                            for (int j=csts.size();j-->0;) {
                                if (contains(vars,csts[j]) || contains(a,csts[j]))
                                    csts.erase(csts.begin()+j);
                            }
                            if (!csts.empty()) {
                                vecteur csol=solve_quiet(symb_equal(taylor_terms[k-2],0),csts,contextptr);
                                vecteur sub=subst(csol,a,*it,false,contextptr);
                                if (!is_constant_wrt_vars(sub,csts,contextptr))
                                    sub=solve_quiet(sub,csts,contextptr);
                                if (!sub.empty()) {
                                    *logptr(contextptr)
                                        << gettext("Warning") << ": " << gettext("assuming ") << csts << gettext(" not in ")
                                        << simp(sub,contextptr) << "\n";
                                }
                            }
                            gen p=simp(subst(taylor_terms[k-2],a,*it,false,contextptr),contextptr);
                            if (is_zero(p,contextptr))
                                continue;
                            else if (k%2)
                                cls=_CPCLASS_SADDLE;
                            else {
                                gen gL,tmp(undef);
                                if ((approx_hp=approx_hompol)) {
                                    vecteur hpvars=*_lname(p,contextptr)._VECTptr;
                                    const_iterateur hit=hpvars.begin();
                                    for (;hit!=hpvars.end() && contains(fvars,*hit);++hit);
                                    if (hit==hpvars.end()) {
                                        gL=_fMin(makesequence(p,vecteur(1,symb_equal(sp,0)),fvars,ip),contextptr);
                                        if (gL.type==_VECT && int(gL._VECTptr->size())==n) {
                                            tmp=vecteur(1,_epsilon2zero(subst(p,fvars,*gL._VECTptr,false,contextptr),contextptr));
                                            gL=_fMax(makesequence(p,vecteur(1,symb_equal(sp,0)),fvars,ip),contextptr);
                                            if (gL.type==_VECT && int(gL._VECTptr->size())==n)
                                                tmp._VECTptr->push_back(_epsilon2zero(subst(p,fvars,*gL._VECTptr,false,contextptr),contextptr));
                                            else tmp=undef;
                                        }
                                    }
                                }
                                if (is_undef(tmp)) {
                                    approx_hp=false;
                                    fvars.push_back(l);
                                    gL=_grad(makesequence(p-l*sp,fvars),contextptr);
                                    if (gL.type==_VECT)
                                        tmp=solve_quiet(gL,fvars,contextptr);
                                }
                                if (tmp.type!=_VECT || tmp._VECTptr->empty() || (!approx_hp && !ckmatrix(tmp))) {
                                    if (k>2) break;
                                    // apply the second order derivative test
                                    matrice chess=*simp(eval_continuous(hess,vars,bnds,*it->_VECTptr,contextptr),contextptr)._VECTptr;
                                    if (is_undef(chess))
                                        break;
                                    if (_lname(chess,contextptr)._VECTptr->empty()) // no symbols here
                                        chess=*_evalf(chess,contextptr)._VECTptr;
                                    gen res=_eigenvals(chess,contextptr);
                                    if (res.type==_VECT) {
                                        vecteur eigvals=*res._VECTptr;
                                        gen e=undef;
                                        for (const_iterateur et=eigvals.begin();et!=eigvals.end();++et) {
                                            if (is_zero(*et,contextptr)) {
                                                cls=_CPCLASS_UNDECIDED;
                                                break;
                                            } else if (is_undef(e)) {
                                                e=*et;
                                                cls=is_positive(e,contextptr)?_CPCLASS_MIN:_CPCLASS_MAX;
                                            } else if (is_strictly_positive(-e*(*et),contextptr))
                                                cls=_CPCLASS_SADDLE;
                                        }
                                    }
                                    break;
                                }
                                if (approx_hp) {
                                    pmin=tmp._VECTptr->front();
                                    pmax=tmp._VECTptr->back();
                                } else {
                                    vecteur lst;
                                    const_iterateur mt=tmp._VECTptr->begin();
                                    for (;mt!=tmp._VECTptr->end();++mt) {
                                        lst.push_back(simp(subst(p,fvars,*mt->_VECTptr,false,contextptr),contextptr));
                                    }
                                    pmin=_min(lst,contextptr);
                                    pmax=_max(lst,contextptr);
                                    fvars.pop_back();
                                }
                                if (is_zero(pmin,contextptr))
                                    cls=_CPCLASS_POSSIBLE_MIN;
                                else if (is_zero(pmax,contextptr))
                                    cls=_CPCLASS_POSSIBLE_MAX;
                                else if (is_strictly_positive(pmin,contextptr))
                                    cls=_CPCLASS_MIN;
                                else if (is_strictly_positive(pmax,contextptr))
                                    cls=_CPCLASS_SADDLE;
                                else
                                    cls=_CPCLASS_MAX;
                            }
                            break;
                        }
                        restore_parameter_assumptions(parm_asmp,contextptr);
                    }
                    cpt_class=gen(cls);
                }
            }
        }
    }
}

/* Find critical points of f subject to constr=0 and ineq<=0 w.r.t. vars in bnds.
 * If an initial point is provided, find a single critical point using fsolve.
 * Classify critical points with respect to ord and approx_hp.
 * If cont_ext is TRUE, then ord must be 1 and critical points will also contain those
 * of an extension of f by continuity.
 * This function returns FALSE if the process is interrupted by user, otherwise it returns TRUE. */
bool find_critical(const gen &f,const vecteur &constr_orig,const vecteur &vars,const vecteur &bnds,const vecteur &ineq_orig,bool open,
                   const vecteur &initial,const gen_map &parm_asmp,int ord,bool approx_hp,bool cont_ext,vecteur &cpts,GIAC_CONTEXT) {
    ipdiff::ivectors arrs;
    vecteur constr=constr_orig,ineq=ineq_orig;
    if (ord>0 && !constr.empty()) {
        log_output_redirect lor(contextptr);
        matrice J;
        if (!jacobian(constr,vars,J,contextptr))
            return true;
        vector<int> ind=linearly_dependent_rows(J,contextptr);
        for (vector<int>::const_reverse_iterator it=ind.rbegin();it!=ind.rend();++it) {
            ineq.push_back(constr[*it]);
            ineq.push_back(-constr[*it]);
            constr.erase(constr.begin()+*it);
            J.erase(J.begin()+*it);
        }
        if (constr.size()>=vars.size())
            return true; // too many constraints
        vars_arrangements(J,arrs,contextptr);
    } else {
        ipdiff::ivector arr(vars.size());
        for (int i=vars.size();i-->0;) {
            arr[i]=i;
        }
        arrs.push_back(arr);
    }
    vecteur tmp_vars(vars.size()),tmp_bnds(vars.size()),tmp_initial(initial.size());
    /* iterate through all possible variable arrangements */
    bool undetected=false;
    for (ipdiff::ivectors::const_iterator ait=arrs.begin();ait!=arrs.end();++ait) {
        const ipdiff::ivector &arr=*ait;
        for (ipdiff::ivector::const_iterator it=arr.begin();it!=arr.end();++it) {
            tmp_vars[it-arr.begin()]=vars[*it];
            tmp_bnds[it-arr.begin()]=bnds[*it];
            if (!initial.empty())
                tmp_initial[it-arr.begin()]=initial[*it];
        }
        try {
            find_local_extrema(cpts,f,constr,tmp_vars,arr,ineq,tmp_bnds,open,tmp_initial,parm_asmp,ord,approx_hp,cont_ext,contextptr);
        } catch (const std::runtime_error &e) {
            purge_variables(tmp_vars,contextptr);
            undetected=true;
            if (strstr(e.what(),"user interruption")!=NULL || interrupted || ctrl_c)
                throw e;
        }
    }
    if (undetected)
        print_warning("some critical points may have been undetected",contextptr);
    return true;
}

/* Return the list of implied inequalities considering the domain of G.
 * Inequalities are computed in form h(<=0) and stored in INEQ. */
void implied_inequalities(const gen &g,const vecteur &vars,vecteur &ineq,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            implied_inequalities(*it,vars,ineq,contextptr);
        }
    } else if (g.type==_SYMB) {
        const gen &f=g._SYMBptr->feuille,&s=g._SYMBptr->sommet;
        if (!is_constant_wrt_vars(f,vars,contextptr)) {
            if (s==at_sqrt || s==at_ln)
                ineq.push_back(symb_superieur_egal(f,0));
            else if (s==at_pow && is_real_number(f._VECTptr->back(),contextptr) && !f._VECTptr->back().is_integer())
                ineq.push_back(symb_superieur_egal(f._VECTptr->front(),0));
            else if (s==at_asin || s==at_acos || s==at_atanh) {
                ineq.push_back(symb_inferieur_egal(f,1));
                ineq.push_back(symb_superieur_egal(f,-1));
            } else if (s==at_acosh)
                ineq.push_back(symb_superieur_egal(f,1));
        }
        implied_inequalities(f,vars,ineq,contextptr);
    }
}

/* Find global minimum and maximum recursively, return false if interrupted. */
bool find_global_extrema_recursively(const gen &f,const vecteur &g,const vecteur &h,const vecteur &vars,const vecteur &bnds,
                                     const gen_map &parm_asmp,const map<int,bool> &sel_orig,const map<int,gen> &fx_orig,gen &mn,gen &mx,vecteur &min_loc,vecteur &ass,GIAC_CONTEXT) {
    vecteur constr_orig=h;
    for (map<int,bool>::const_iterator it=sel_orig.begin();it!=sel_orig.end();++it) {
        if (it->second)
            constr_orig.push_back(g[it->first]);
    }
    int m=g.size(),i=-1;
    if (!sel_orig.empty())
        i=sel_orig.rbegin()->first+1;
    for (;i<m;++i) {
        map<int,gen> fx=fx_orig;
        map<int,bool> sel=sel_orig;
        vecteur constr=constr_orig,vars_rd,bnds_rd,ineq=g,cpts,cv;
        gen a,b,f_rd=f;
        if (i>=0) {
            const_iterateur jt=vars.begin();
            for (;jt!=vars.end();++jt) {
                if (is_linear_wrt(g[i],*jt,a,b,contextptr) &&
                        is_constant_wrt_vars(a,vars,contextptr) && !is_zero(a,contextptr) && is_constant_wrt_vars(b,vars,contextptr))
                    break;
            }
            if (jt!=vars.end()) { // there is a variable with constant value
                int j=jt-vars.begin();
                gen jval=-b/a;
                map<int,gen>::const_iterator kt=fx.find(j);
                if (kt!=fx.end() && !is_zero(simp(kt->second-jval,contextptr)))
                    continue; // infeasible
                if (kt==fx.end())
                    fx[j]=jval;
                sel[i]=false;
            } else {
                constr.push_back(g[i]);
                sel[i]=true;
            }
        }
        vector<int> vars_ind;
        for (const_iterateur it=vars.begin();it!=vars.end();++it) {
            int j=it-vars.begin();
            if (fx.find(j)==fx.end()) {
                vars_rd.push_back(*it);
                bnds_rd.push_back(bnds[j]);
                vars_ind.push_back(j);
            }
        }
        for (map<int,gen>::const_iterator it=fx.begin();it!=fx.end();++it) {
            f_rd=subst(f_rd,vars[it->first],it->second,false,contextptr);
            constr=subst(constr,vars[it->first],it->second,false,contextptr);
            ineq=subst(ineq,vars[it->first],it->second,false,contextptr);
        }
        constr=*simp(constr,contextptr)._VECTptr;
        f_rd=simp(f_rd,contextptr);
        int ci=constr.size();
        for (;ci-->0;) {
            const gen &c=constr[ci];
            if (is_real_number(c,contextptr) && !is_zero(to_real_number(c,contextptr),contextptr))
                break; // infeasible
            if (is_constant_wrt_vars(c,vars,contextptr))
                constr.erase(constr.begin()+ci);
        }
        if (ci>=0)
            continue;
        bool stopped=false;
        try {
            if (is_constant_wrt_vars(f_rd,vars_rd,contextptr)) {
                bound_variables(vars_rd,bnds_rd,false,contextptr);
                if (!vars_rd.empty()) {
                    vecteur sol=solve2(constr,vars_rd,contextptr);
                    for (const_iterateur it=sol.begin();it!=sol.end();++it) {
                        cpts.push_back(makevecteur(*it,0,gen_map()));
                    }
                } else cpts=vecteur(1,makevecteur(vecteur(0),0,gen_map()));
            } else {
                if (!find_critical(f_rd,constr,vars_rd,bnds_rd,ineq,false,vecteur(0),parm_asmp,1,false,true,cpts,contextptr))
                    return false;
                if (cpts.empty() && !constr.empty()) {
                    bound_variables(vars_rd,bnds_rd,false,contextptr);
                    vecteur sol=solve2(constr,vars_rd,contextptr);
                    for (const_iterateur it=sol.begin();it!=sol.end();++it) {
                        cpts.push_back(makevecteur(*it,0,gen_map()));
                    }
                }
            }
        } catch (const std::runtime_error &e) {
            if (strstr(e.what(),"user interruption")!=NULL || interrupted || ctrl_c) {
                interrupted=ctrl_c=false;
                print_error("stopped by user interruption",contextptr);
                stopped=true;
            } else print_warning("some critical points may have been undetected",contextptr);
        }
        purge_variables(vars_rd,contextptr);
        if (stopped)
            return false;
        for (iterateur it=cpts.begin();it!=cpts.end();++it) {
            gen &cpt=it->_VECTptr->front();
            if (cpt.type!=_VECT)
                cpt=vecteur(1,cpt);
            vecteur fullcpt(vars.size(),undef);
            for (map<int,gen>::const_iterator jt=fx.begin();jt!=fx.end();++jt) {
                fullcpt[jt->first]=jt->second;
            }
            int vi=0;
            for (iterateur jt=fullcpt.begin();jt!=fullcpt.end();++jt) {
                if (is_undef(*jt))
                    *jt=cpt._VECTptr->at(vi++);
            }
            cpt=fullcpt;
        }
        filter_cpts(cpts,vars,bnds,g,false,f,contextptr);
        // update global min and max
        for (const_iterateur it=cpts.begin();it!=cpts.end();++it) {
            vecteur cpt=*it->_VECTptr->front()._VECTptr;
            gen asmp=it->_VECTptr->size()==3?it->_VECTptr->back():gen_map();
            gen val=simp(eval_continuous(f,vars,bnds,cpt,contextptr),contextptr);
            vecteur gg=*simp(subst(g,vars,cpt,false,contextptr),contextptr)._VECTptr,vars_gg;
            bool is_feas=true,append_domain=false;
            for (int j=gg.size();j-->0;) {
                if (gg[j].is_integer() && gg[j].subtype==_INT_BOOLEAN) {
                    if (gg[j].val==0) {
                        is_feas=false;
                        break;
                    }
                    gg.erase(gg.begin()+j);
                }
                if (is_positive(-gg[j],contextptr))
                    gg.erase(gg.begin()+j);
            }
            if (!is_feas)
                continue;
            find_vars(gg,vars,vars_gg,contextptr);
            for (const_iterateur jt=cpt.begin();jt!=cpt.end();++jt) {
                const gen &var=vars[jt-cpt.begin()];
                if (*jt==var && contains(vars_gg,var))
                    append_domain=true;
            }
            if (append_domain) {
                if (vars_gg.size()==1) {
                    vecteur sol=solve_quiet(gg,vars_gg.front(),contextptr);
                    if (sol.empty())
                        continue; // infeasible
                    if (sol.front()==vars_gg.front())
                        sol.clear();
                    cpt=mergevecteur(cpt,sol.empty()?*_zip(makesequence(at_inferieur_egal,gg,vecteur(gg.size(),0)),contextptr)._VECTptr:sol);
                } else cpt=mergevecteur(cpt,*_zip(makesequence(at_inferieur_egal,gg,vecteur(gg.size(),0)),contextptr)._VECTptr);
            }
            if (!is_undef(mn) && is_zero(simp(val-mn,contextptr),contextptr)) {
                if (!contains(min_loc,cpt)) {
                    min_loc.push_back(cpt);
                    ass.push_back(asmp);
                }
            } else if (is_undef(mn) || is_strictly_greater(mn,val,contextptr)) {
                mn=val;
                min_loc=vecteur(1,cpt);
                ass=vecteur(1,asmp);
            }
            if (is_undef(mx) || is_strictly_greater(val,mx,contextptr))
                mx=val;
            if (!is_undef(mn) && !is_strictly_greater(mn,val,contextptr) && !is_greater(val,mn,contextptr)) { // unable to sort
                print_error("failed to minimize",contextptr);
                vecteur prm=*_lname(makevecteur(mn,val),contextptr)._VECTptr;
                if (!prm.empty()) {
                    *logptr(contextptr) << gettext("The minimum depends on: ");
                    for (const_iterateur it=prm.begin();it!=prm.end();++it) {
                        if (it!=prm.begin())
                            *logptr(contextptr) << ",";
                        *logptr(contextptr) << *it;
                    }
                }
                return false;
            }
        }
        // recurse
        if (i>=0 && !find_global_extrema_recursively(f,g,h,vars,bnds,parm_asmp,sel,fx,mn,mx,min_loc,ass,contextptr))
            return false;
    }
    return true;
}

bool parse_variables(const gen &g,vecteur &vars,vecteur &ineq,vecteur &initial,GIAC_CONTEXT) {
    vecteur varlist=gen2vecteur(g);
    for (const_iterateur it=varlist.begin();it!=varlist.end();++it) {
        if (is_equal(*it)) {
            vecteur &ops=*it->_SYMBptr->feuille._VECTptr;
            gen &v=ops.front(), &rh=ops.back();
            vars.push_back(v);
            if (rh.is_symb_of_sommet(at_interval)) {
                vecteur &range=*rh._SYMBptr->feuille._VECTptr;
                if (!is_inf(range.front())) {
                    if (!is_real_number(range.front(),contextptr))
                        return false;
                    ineq.push_back(symb_inferieur_egal(range.front()-v,0));
                }
                if (!is_inf(range.back())) {
                    if (!is_real_number(range.back(),contextptr))
                        return false;
                    ineq.push_back(symb_inferieur_egal(v-range.back(),0));
                }
            } else initial.push_back(rh);
        } else vars.push_back(*it);
    }
    return ckvars(vars,contextptr) &&
            (initial.empty() || (initial.size()==vars.size() && _lname(initial,contextptr)._VECTptr->empty()));
}

bool minimize(const gen &f_orig,const vecteur &constr,const vecteur &vars,const gen_map &parm_asmp,const vecteur &initial,
              gen &mn,gen &mx,vecteur &loc,vecteur &ass,GIAC_CONTEXT) {
    // handle conjunctions
    vecteur g;
    gen f=f_orig;
    for (const_iterateur it=constr.begin();it!=constr.end();++it) {
        if (is_conjunction(*it) && it->_SYMBptr->feuille.type==_VECT)
            g=mergevecteur(g,*it->_SYMBptr->feuille._VECTptr);
        else g.push_back(*it);
    }
    // handle disjunctions
    for (const_iterateur it=g.begin();it!=g.end();++it) {
        if (it->is_symb_of_sommet(at_ou) && it->_SYMBptr->feuille.type==_VECT) {
            for (const_iterateur jt=it->_SYMBptr->feuille._VECTptr->begin();jt!=it->_SYMBptr->feuille._VECTptr->end();++jt) {
                vecteur gj=g;
                gj[it-g.begin()]=*jt;
                if (!minimize(f,gj,vars,parm_asmp,initial,mn,mx,loc,ass,contextptr))
                    return false;
            }
            return true;
        }
    }
    vecteur feu,g1,g2;
    gen ts=temp_symb("tpw",-1,contextptr),val,f1,f2;
    if (find_piecewise(f,ts,feu)) {
        if (feu.size()==2) {
            g.push_back(feu.front());
            f=subst(f,ts,feu.back(),false,contextptr);
            return minimize(f,g,vars,parm_asmp,initial,mn,mx,loc,ass,contextptr);
        } else if (feu.size()==3) {
            g1=g2=g;
            g1.push_back(feu.front());
            g2.push_back(_eval(symb_not(feu.front()),contextptr));
            f1=subst(f,ts,feu.at(1),false,contextptr);
            f2=subst(f,ts,feu.at(2),false,contextptr);
            return minimize(f1,g1,vars,parm_asmp,initial,mn,mx,loc,ass,contextptr) &&
                   minimize(f2,g2,vars,parm_asmp,initial,mn,mx,loc,ass,contextptr);
        } else {
            print_error("invalid piecewise expression",contextptr);
            return false;
        }
    }
    if (find_abs(f,ts,val)) {
        g1=g2=g;
        g1.push_back(symb_superieur_egal(val,0));
        g2.push_back(symb_inferieur_egal(val,0));
        f1=subst(f,ts,val,false,contextptr);
        f2=subst(f,ts,-val,false,contextptr);
        return minimize(f1,g1,vars,parm_asmp,initial,mn,mx,loc,ass,contextptr) &&
                minimize(f2,g2,vars,parm_asmp,initial,mn,mx,loc,ass,contextptr);
    }
    implied_inequalities(f,vars,g,contextptr);
    vecteur h;
    int wi;
    for (int i=g.size();i-->0;) {
        g[i]=_eval(g[i],contextptr);
        if (is_equal(g[i])) {
            h.push_back(equal2diff(g[i]));
            g.erase(g.begin()+i);
        } else if ((wi=which_ineq(g[i]))!=0) {
            const vecteur &s=*g[i]._SYMBptr->feuille._VECTptr;
            g[i]=wi==1?s[0]-s[1]:s[1]-s[0];
        } else if (g[i].is_integer() && g[i].subtype==_INT_BOOLEAN) {
            if (g[i].val==0)
                return true; // infeasible
            g.erase(g.begin()+i);
        } else {
            h.push_back(g[i]);
            g.erase(g.begin()+i);
        }
    }
    vecteur bnds;
    if (!get_variable_bounds(vars,g,bnds,contextptr))
        return true; // infeasible
    if (!initial.empty()) { // approx local minimization
        gen args=makesequence(f,mergevecteur(
                                *_zip(makesequence(at_inferieur_egal,g,vecteur(g.size(),0)),contextptr)._VECTptr,
                                *_zip(makesequence(at_equal,h,vecteur(h.size(),0)),contextptr)._VECTptr),
                                symb_equal(change_subtype(_NLP_INITIALPOINT,_INT_MAPLECONVERSION),
                                            _zip(makesequence(at_equal,vars,initial),contextptr)));
        gen asol=_nlpsolve(args,contextptr);
        if (asol.type==_VECT && asol._VECTptr->size()==2) {
            mn=asol._VECTptr->front();
            if (mn.type==_STRNG) // failed to find minimum
                return false;
            vecteur lc(vars.size());
            gen pos=asol._VECTptr->at(1),v;
            if (pos.type!=_VECT || pos._VECTptr->size()!=vars.size())
                return false;
            for (const_iterateur it=pos._VECTptr->begin();it!=pos._VECTptr->end();++it) {
                if (!is_equal(*it))
                    return false;
                int j=indexof(it->_SYMBptr->feuille._VECTptr->front(),vars);
                if (j<0)
                    return false;
                lc[j]=it->_SYMBptr->feuille._VECTptr->back();
            }
            loc.push_back(lc);
            ass.push_back(gen_map());
        } else return false;
    } else {
        if (!find_global_extrema_recursively(f,g,h,vars,bnds,parm_asmp,map<int,bool>(),map<int,gen>(),mn,mx,loc,ass,contextptr))
            return false; // interrupted
    }
    return true;
}

bool extrema(const gen &f,const vecteur &h,const vecteur &constr,const vecteur &vars,const vecteur &initial,
             const gen_map &parm_asmp,int ord,bool approx_hp,vecteur &cpts,GIAC_CONTEXT) {
    // handle conjunctions
    vecteur g;
    for (const_iterateur it=constr.begin();it!=constr.end();++it) {
        if (is_conjunction(*it) && it->_SYMBptr->feuille.type==_VECT)
            g=mergevecteur(g,*it->_SYMBptr->feuille._VECTptr);
        else g.push_back(*it);
    }
    // handle disjunctions
    for (const_iterateur it=g.begin();it!=g.end();++it) {
        if (it->is_symb_of_sommet(at_ou) && it->_SYMBptr->feuille.type==_VECT) {
            for (const_iterateur jt=it->_SYMBptr->feuille._VECTptr->begin();jt!=it->_SYMBptr->feuille._VECTptr->end();++jt) {
                vecteur gj=g;
                gj[it-g.begin()]=*jt;
                if (!extrema(f,h,gj,vars,initial,parm_asmp,ord,approx_hp,cpts,contextptr))
                    return false;
            }
            return true;
        }
    }
    for (int i=g.size();i-->0;) {
        if (g[i].is_integer() && g[i].subtype==_INT_BOOLEAN) {
            if (g[i].val==0)
                return true; // infeasible
            g.erase(g.begin()+i);
            continue;
        }
        int wi=which_ineq(g[i]);
        assert(wi!=0);
        const vecteur &s=*g[i]._SYMBptr->feuille._VECTptr;
        g[i]=wi==1?s[0]-s[1]:s[1]-s[0];
    }
    vecteur bnds;
    if (!get_variable_bounds(vars,g,bnds,contextptr))
        return true; // infeasible
    return find_critical(f,h,vars,bnds,g,true,initial,parm_asmp,ord,approx_hp,true,cpts,contextptr);
}

vecteur make_temp_vars(const vecteur &vars,vecteur &constr,GIAC_CONTEXT) {
    vecteur tmpvars(vars.size());
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        int i=it-vars.begin(),dom;
        gen &tv=tmpvars[i];
        tv=temp_symb("var",i,contextptr);
        matrice intrv;
        vecteur exc,adc;
        get_assumptions(*it,dom,intrv,exc,contextptr);
        if (dom>0 && dom!=1)
            *logptr(contextptr) << gettext("Warning") << ": " << gettext("assuming that ") << *it << gettext(" is real");
        for (const_iterateur jt=intrv.begin();jt!=intrv.end();++jt) {
            const vecteur &b=*jt->_VECTptr;
            assert(b.size()==2);
            adc.push_back(symb_and(symb_superieur_egal(*it,b.front()),symb_inferieur_egal(*it,b.back())));
        }
        if (adc.size()>1)
            constr.push_back(symbolic(at_ou,adc));
        else if (adc.size()==1)
            constr.push_back(adc.front());
    }
    return tmpvars;
}

/*
 * Function 'minimize' minimizes a multivariate continuous real function on a
 * closed and bounded region using the method of Lagrange multipliers. The
 * feasible region is specified by bounding variables or by adding one or more
 * (in)equality constraints.
 *
 * Usage
 * ^^^^^
 *     minimize(obj,[constr],vars,[opt])
 *
 * Parameters
 * ^^^^^^^^^^
 *   - obj                 : objective function to minimize
 *   - constr (optional)   : single equality or inequality constraint or
 *                           a list of constraints, if constraint is given as
 *                           expression it is assumed that it is equal to zero
 *   - vars                : single variable or a list of problem variables, where
 *                           optional bounds of a variable may be set by appending '=a..b'
 *   - location (optional) : the option keyword 'locus' or 'coordinates' or 'point'
 *
 * Objective function must be continuous in all points of the feasible region,
 * which is assumed to be closed and bounded. If one of these condinitions is
 * not met, the final result may be incorrect.
 *
 * When the fourth argument is specified, point(s) at which the objective
 * function attains its minimum value are also returned as a list of vector(s).
 * The keywords 'locus', 'coordinates' and 'point' all have the same effect.
 * For univariate problems, a vector of numbers (x values) is returned, while
 * for multivariate problems it is a vector of vectors, i.e. a matrix.
 *
 * The function attempts to obtain the critical points in exact form, if the
 * parameters of the problem are all exact. It works best for problems in which
 * the lagrangian function gradient and the constraints are rational expressions.
 * Points at which the function is not differentiable are also considered critical.
 * This function handles univariate piecewise functions.
 *
 * If no critical points were obtained, the return value is undefined.
 *
 * Examples
 * ^^^^^^^^
 * minimize(sin(x),[x=0..4])
 *    >> sin(4)
 * minimize(asin(x),x=-1..1)
 *    >> -pi/2
 * minimize(x^2+cos(x),x=0..3)
 *    >> 1
 * minimize(x^4-x^2,x=-3..3,locus)
 *    >> [-1/4,[(sqrt(2))/2,-(sqrt(2))/2]]
 * minimize(abs(x),x=-1..1)
 *    >> 0
 * minimize(x-abs(x),x=-1..1)
 *    >> -2
 * minimize(abs(exp(-x^2)-1/2),x=-4..4)
 *    >> 0
 * minimize(piecewise(x<=-2,x+6,x<=1,x^2,3/2-x/2),x=-3..2)
 *    >> 0
 * minimize(x^2-3x+y^2+3y+3,[x=2..4,y=-4..-2],point)
 *    >> -1,[[2,-2]]
 * minimize(2x^2+y^2,x+y=1,[x,y])
 *    >> 2/3
 * minimize(2x^2-y^2+6y,x^2+y^2<=16,[x,y])
 *    >> -40
 * minimize(x*y+9-x^2-y^2,x^2+y^2<=9,[x,y])
 *    >> -9/2
 * minimize(sqrt(x^2+y^2)-z,[x^2+y^2<=16,x+y+z=10],[x,y,z])
 *    >> -4*sqrt(2)-6
 * minimize(x*y*z,x^2+y^2+z^2=1,[x,y,z])
 *    >> -sqrt(3)/9
 * minimize(sin(x)+cos(x),x=0..20,coordinates)
 *    >> -sqrt(2),[5*pi/4,13*pi/4,21*pi/4]
 * minimize((1+x^2+3y+5x-4*x*y)/(1+x^2+y^2),x^2/4+y^2/3=9,[x,y])
 *    >> -2.44662490691
 * minimize(x^2-2x+y^2+1,[x+y<=0,x^2<=4],[x,y],locus)
 *    >> 1/2,[[1/2,-1/2]]
 * minimize(x^2*(y+1)-2y,[y<=2,sqrt(1+x^2)<=y],[x,y])
 *    >> -4
 * minimize(4x^2+y^2-2x-4y+1,4x^2+y^2=1,[x,y])
 *    >> -sqrt(17)+2
 * minimize(cos(x)^2+cos(y)^2,x+y=pi/4,[x,y],locus)
 *    >> (-sqrt(2)+2)/2,[[5*pi/8,-3*pi/8]]
 * minimize(x^2+y^2,x^4+y^4=2,[x,y])
 *    >> 1.41421356237
 * minimize(z*x*exp(y),z^2+x^2+exp(2y)=1,[x,y,z])
 *    >> -sqrt(3)/9
 * minimize(y,x^2-x^4-y^5=0,[x=-1..1,y=0..1])
 *    >> 0
 * minimize(x*y,x^2-x^4-y^5=0,[x=-1..1,y=0..1])
 *    >> -sqrt(7)*42^(1/5)/9
 */
gen _minimize(const gen &args,GIAC_CONTEXT) {
    if (args.type==_STRNG && args.subtype==-1) return args;
    if (args.type!=_VECT || args.subtype!=_SEQ__VECT || args._VECTptr->size()>4)
        return gentypeerr(contextptr);
    vecteur &argv=*args._VECTptr,constr,vars,ineq,initial;
    bool location=false;
    int nargs=argv.size();
    if (argv.back()==at_coordonnees || argv.back()==at_lieu || argv.back()==at_point) {
        location=true;
        --nargs;
    }
    if (nargs==3)
        constr=gen2vecteur(argv[1]);
    if (!parse_variables(argv[nargs-1],vars,constr,initial,contextptr))
        return generr("Invalid specification of variables");
    gen f=make_piecewise_nested(argv[0]);
    if (is_constant_wrt_vars(f,vars,contextptr))
        return location?makevecteur(f,vars):f; // function is constant
    vecteur parm=*_lname(makevecteur(f,constr),contextptr)._VECTptr;
    for (int i=parm.size();i-->0;) {
        if (contains(vars,parm[i]))
            parm.erase(parm.begin()+i);
    }
    if (!initial.empty() && !parm.empty())
        return generr("Cannot perform local search with symbolic constants");
    gen_map parm_asmp;
    if (!get_parameter_assumptions(parm,parm_asmp,contextptr))
        return generr("Symbolic constants must be real");
    vecteur tmpvars=make_temp_vars(vars,constr,contextptr);
    f=subst(f,vars,tmpvars,false,contextptr);
    constr=subst(constr,vars,tmpvars,false,contextptr);
    gen mn(undef),mx(undef);
    vecteur loc,ass;
    if (!minimize(f,constr,tmpvars,parm_asmp,initial,mn,mx,loc,ass,contextptr)) { // interrupted or failed
        return undef;
    }
    if (is_undef(mn)) {
        print_error("failed to minimize, no critical points found",contextptr);
        return undef;
    }
    if (is_inf(mn))
        return generr("Objective function is unbounded");
    if (location) {
        loc=subst(loc,tmpvars,vars,false,contextptr);
#if 1
        for (int i=loc.size();i-->0;) {
            const gen_map &asmp=*ass[i]._MAPptr;
            if (!asmp.empty()) {
                *logptr(contextptr) << gettext("Critical point") << " " << loc << " " << gettext("exists under assumptions") << " ";
                for (gen_map::const_iterator it=asmp.begin();it!=asmp.end();++it) {
                    for (const_iterateur jt=it->second._VECTptr->begin();jt!=it->second._VECTptr->end();++jt) {
                        if (jt!=it->second._VECTptr->begin())
                            *logptr(contextptr) << ",";
                        *logptr(contextptr) << subst(*jt,tmpvars,vars,false,contextptr);
                    }
                }
                *logptr(contextptr) << "\n";
            }
        }
#endif
        return makevecteur(mn,loc);
    }
    return mn;
}
static const char _minimize_s []="minimize";
static define_unary_function_eval (__minimize,&_minimize,_minimize_s);
define_unary_function_ptr5(at_minimize,alias_at_minimize,&__minimize,0,true)

/*
 * 'maximize' takes the same arguments as the function 'minimize', but
 * maximizes the objective function. See 'minimize' for details.
 *
 * Examples
 * ^^^^^^^^
 * maximize(cos(x),x=1..3)
 *    >> cos(1)
 * maximize(piecewise(x<=-2,x+6,x<=1,x^2,3/2-x/2),x=-3..2)
 *    >> 4
 * maximize(x^2-3x+y^2+3y+3,[x=2..4,y=-4..-2])
 *    >> 11
 * maximize(x*y*z,x^2+2*y^2+3*z^2<=1,[x,y,z],point)
 *    >> sqrt(2)/18,[[-sqrt(3)/3,sqrt(6)/6,-1/3],[sqrt(3)/3,-sqrt(6)/6,-1/3],
 *                   [-sqrt(3)/3,-sqrt(6)/6,1/3],[sqrt(3)/3,sqrt(6)/6,1/3]]
 * maximize(x^2-x*y+2*y^2,[x=-1..0,y=-1/2..1/2],coordinates)
 *    >> 2,[[-1,1/2]]
 * maximize(x*y,[x+y^2<=2,x>=0,y>=0],[x,y],locus)
 *    >> 4*sqrt(6)/9,[[4/3,sqrt(6)/3]]
 * maximize(y^2-x^2*y,y<=x,[x=0..2,y=0..2])
 *    >> 4
 * maximize(2x+y,4x^2+y^2=8,[x,y])
 *    >> 4
 * maximize(x^2*(y+1)-2y,[y<=2,sqrt(1+x^2)<=y],[x,y])
 *    >> 5
 * maximize(4x^2+y^2-2x-4y+1,4x^2+y^2=1,[x,y])
 *    >> sqrt(17)+2
 * maximize(3x+2y,2x^2+3y^2<=3,[x,y])
 *    >> sqrt(70)/2
 * maximize(x*y,[2x+3y<=10,x>=0,y>=0],[x,y])
 *    >> 25/6
 * maximize(x^2+y^2+z^2,[x^2/16+y^2+z^2=1,x+y+z=0],[x,y,z])
 *    >> 8/3
 * assume(a>0);maximize(x^2*y^2*z^2,x^2+y^2+z^2=a^2,[x,y,z])
 *    >> a^6/27
 */
gen _maximize(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()<2)
        return gentypeerr(contextptr);
    vecteur gv(*g._VECTptr);
    gv[0]=-gv[0];
    gen res=_minimize(_feuille(gv,contextptr),contextptr);
    if (res.type==_VECT && res._VECTptr->size()>0) {
        res._VECTptr->front()=ratnormal(-res._VECTptr->front(),contextptr);
    } else if (res.type!=_VECT)
        res=ratnormal(-res,contextptr);
    return res;
}
static const char _maximize_s []="maximize";
static define_unary_function_eval (__maximize,&_maximize,_maximize_s);
define_unary_function_ptr5(at_maximize,alias_at_maximize,&__maximize,0,true)

/*
 * 'extrema' attempts to find all points of strict local minima/maxima of a
 * smooth (uni/multi)variate function (subject to equality constraints).
 *
 * Usage
 * ^^^^^
 *     extrema(expr,[constr],vars,[order=n||lagrange])
 *
 * Parameters
 * ^^^^^^^^^^
 *   - expr               : differentiable expression
 *   - constr (optional)  : (list of) equality constraint(s)
 *   - vars               : (list of) problem variable(s)
 *   - order (optional)   : specify 'order=<nonnegative integer>' to
 *                          bound the order of the derivatives being
 *                          inspected when classifying the critical points
 *
 * The number of constraints must be less than the number of variables. When
 * there are more than one constraint/variable, they must be specified in
 * form of list.
 *
 * Variables may be specified with symbol, e.g. 'var', or by using syntax
 * 'var=a..b', which restricts the variable 'var' to the open interval (a,b),
 * where a and b are real numbers or +/-infinity. If variable list includes a
 * specification of initial point, such as, for example, [x=1,y=0,z=2], then
 * numeric solver is activated to find critical point in the vicinity of the
 * given point. In this case, the single critical point, if found, is examined.
 *
 * The function attempts to find the critical points in exact form, if the
 * parameters of the problem are all exact. It works best for problems in which
 * the gradient of lagrangian function consists of rational expressions. The
 * result may be inexact, however, if exact solutions could not be obtained.
 *
 * For classifying critical points, the bordered hessian method is used first.
 * It is only a second order test, so it may be inconclusive in some cases. In
 * these cases function looks at higher-order derivatives, up to order
 * specified by 'order' option (the extremum test). Set 'order' to 1
 * to use only the bordered hessian test or 0 to output critical points without
 * attempting to classify them. Setting 'order' to 2 or higher activates
 * checking for saddle points and inspecting higher derivatives (up to 'order')
 * to determine the nature of some or all unclassified critical points.
 * By default, 'order' equals to 5.
 *
 * The return value is a sequence with two elements: list of strict local
 * minima and list of strict local maxima. If only critical points are
 * requested (by setting 'order' to 0), the output consists of a single
 * list, as no classification was attempted. For univariate problems the points
 * are real numbers, while for multivariate problems they are specified as
 * lists of coordinates, so "lists of points" are in fact matrices with rows
 * corresponding to points in multivariate cases, i.e. vectors in univariate
 * cases.
 *
 * The function prints out information about saddle/inflection points and also
 * about critical points for which no decision could be made, so that the user
 * can inspect candidates for local extrema by plotting the graph, for example.
 *
 * Examples
 * ^^^^^^^^
 * extrema(-2*cos(x)-cos(x)^2,x)
 *    >> [0],[pi]
 * extrema((x^3-1)^4/(2x^3+1)^4,x=0..inf)
 *    >> [1],[]
 * extrema(x/2-2*sin(x/2),x=-12..12)
 *    >> [2*pi/3,-10*pi/3],[10*pi/3,-2*pi/3]
 * extrema(x-ln(abs(x)),x)
 *    >> [1],[]
 * assume(a>=0);extrema(x^2+a*x,x)
 *    >> [-a/2],[]
 * extrema(x^7+3x^6+3x^5+x^4+2x^2-x,x)
 *    >> [0.225847362349],[-1.53862319761]
 * extrema((x^2+x+1)/(x^4+1),x)
 *    >> [],[0.697247087784]
 * extrema(x^2+exp(-x),x)
 *    >> [LambertW(1/2)],[]
 * extrema(exp(-x)*ln(x),x)
 *    >> [],[exp(LambertW(1))]
 * extrema(tan(x)*(x^3-5x^2+1),x=-0.5)
 *    >> [-0.253519032024],[]
 * extrema(tan(x)*(x^3-5x^2+1),x=0.5)
 *    >> [],[0.272551772027]
 * extrema(exp(x^2-2x)*ln(x)*ln(1-x),x=0.5)
 *    >> [],[0.277769149124]
 * extrema(ln(2+x-sin(x)^2),x=0..2*pi)
 *    >> [],[] (needed to compute third derivative to drop critical points pi/4 and 5pi/4)
 * extrema(x^3-2x*y+3y^4,[x,y])
 *    >> [[12^(1/5)/3,(12^(1/5))^2/6]],[]
 * extrema((2x^2-y)*(y-x^2),[x,y])  //Peano surface
 *    >> [],[] (saddle point at origin)
 * extrema(5x^2+3y^2+x*z^2-z*y^2,[x,y,z])
 *    >> [],[] (possible local minimum at origin, in fact saddle)
 * extrema(3*atan(x)-2*ln(x^2+y^2+1),[x,y])
 *    >> [],[[3/4,0]]
 * extrema(x*y,x+y=1,[x,y])
 *    >> [],[[1/2,1/2]]
 * extrema(sqrt(x*y),x+y=2,[x,y])
 *    >> [],[[1,1]]
 * extrema(x*y,x^3+y^3=16,[x,y])
 *    >> [],[[2,2]]
 * extrema(x^2+y^2,x*y=1,[x=0..inf,y=0..inf])
 *    >> [[1,1]],[]
 * extrema(ln(x*y^2),2x^2+3y^2=8,[x,y])
 *    >> [],[[2*sqrt(3)/3,-4/3],[2*sqrt(3)/3,4/3]]
 * extrema(y^2+4y+2x-x^2,x+2y=2,[x,y])
 *    >> [],[[-2/3,4/3]]
 * assume(a>0);extrema(x/a^2+a*y^2,x+y=a,[x,y])
 *    >> [[(2*a^4-1)/(2*a^3),1/(2*a^3)]],[]
 * extrema(6x+3y+2z,4x^2+2y^2+z^2=70,[x,y,z])
 *    >> [[-3,-3,-4]],[[3,3,4]]
 * extrema(x*y*z,x+y+z=1,[x,y,z])
 *    >> [],[[1/3,1/3,1/3]]
 * extrema(x*y^2*z^2,x+y+z=5,[x,y,z])
 *    >> [],[[1,2,2]]
 * extrema(4y-2z,[2x-y-z=2,x^2+y^2=1],[x,y,z])
 *    >> [[2*sqrt(13)/13,-3*sqrt(13)/13,(7*sqrt(13)-26)/13]],
 *       [[-2*sqrt(13)/13,3*sqrt(13)/13,(-7*sqrt(13)-26)/13]]
 * extrema((x-3)^2+(y-1)^2+(z-1)^2,x^2+y^2+z^2=4,[x,y,z])
 *    >> [[6*sqrt(11)/11,2*sqrt(11)/11,2*sqrt(11)/11]],
 *       [[-6*sqrt(11)/11,-2*sqrt(11)/11,-2*sqrt(11)/11]]
 * extrema(x+3y-z,2x^2+y^2=z,[x,y,z])
 *    >> [],[[1/4,3/2,19/8]]
 * extrema(2x*y+2y*z+x*z,x*y*z=4,[x,y,z])
 *    >> [[2,1,2]],[]
 * extrema(x+y+z,[x^2+y^2=1,2x+z=1],[x,y,z])
 *    >> [[sqrt(2)/2,-sqrt(2)/2,-sqrt(2)+1]],[[-sqrt(2)/2,sqrt(2)/2,sqrt(2)+1]]
 * assume(a>0);extrema(x+y+z,[y^2-x^2=a,x+2z=1],[x,y,z])
 *    >> [[-sqrt(3)*sqrt(a)/3,2*sqrt(3)*sqrt(a)/3,(sqrt(3)*sqrt(a)+3)/6]],
 *       [[sqrt(3)*sqrt(a)/3,-2*sqrt(3)*sqrt(a)/3,(-sqrt(3)*sqrt(a)+3)/6]]
 * extrema((x-u)^2+(y-v)^2,[x^2/4+y^2/9=1,(u-3)^2+(v+5)^2=1],[u,v,x,y])
 *    >> [[2.35433932354,-4.23637555425,0.982084902545,-2.61340692712]],
 *       [[3.41406613851,-5.91024679428,-0.580422508346,2.87088778158]]
 * extrema(x2^6+x1^3+4x1+4x2,x1^5+x2^4+x1+x2=0,[x1,x2])
 *    >> [[-0.787457596325,0.758772338924],[-0.784754836317,-1.23363062357]],
 *       [[0.154340184382,-0.155005038065]]
 * extrema(x*y,-2x^3+15x^2*y+11y^3-24y=0,[x,y])
 *    >> [[sqrt(17)*sqrt(3*sqrt(33)+29)/17,sqrt(-15*sqrt(33)+127)*sqrt(187)/187],
 *        [-sqrt(17)*sqrt(3*sqrt(33)+29)/17,-sqrt(-15*sqrt(33)+127)*sqrt(187)/187],
 *        [sqrt(-3*sqrt(33)+29)*sqrt(17)/17,-sqrt(15*sqrt(33)+127)*sqrt(187)/187],
 *        [-sqrt(-3*sqrt(33)+29)*sqrt(17)/17,sqrt(15*sqrt(33)+127)*sqrt(187)/187]],
 *       [[1,1],[-1,-1],[0,0]]
 * extrema(x2^4-x1^4-x2^8+x1^10,[x1,x2],order=1)
 *    >> [[0,0],[0,-(1/2)^(1/4)],[0,(1/2)^(1/4)],[-(2/5)^(1/6),0],[-(2/5)^(1/6),-(1/2)^(1/4)],
 *        [-(2/5)^(1/6),(1/2)^(1/4)],[(2/5)^(1/6),0],[(2/5)^(1/6),-(1/2)^(1/4)],[(2/5)^(1/6),(1/2)^(1/4)]]
 * extrema(x2^4-x1^4-x2^8+x1^10,[x1,x2])
 *    >> [[],[]]
 * extrema(x2^6+x1^3+4x1+4x2,x1^5+x2^4+x1+x2=0,[x1,x2])
 *    >> [[-0.787457596325,0.758772338924],[-0.784754836317,-1.23363062357]],
 *       [[0.154340184382,-0.155005038065]]
 * extrema(x2^6+x1^3+2x1^2-x2^2+4x1+4x2,x1^5+x2^4+x1+x2=0,[x1,x2])
 *    >> [[-0.662879934158,-1.18571027742],[0,0]],[[0.301887394815,-0.314132376868]]
 * extrema(3x^2-2x*y+y^2-8y,[x,y])
 *    >> [[2,6]],[]
 * extrema(x^3+3x*y^2-15x-12y,[x,y])
 *    >> [[2,1]],[[-2,-1]]
 * extrema(4x*y-x^4-y^4,[x,y])
 *    >> [],[[1,1],[-1,-1]]
 * extrema(x*sin(y),[x,y])
 *    >> [],[]
 * extrema(x^4+y^4,[x,y])
 *    >> [[0,0]],[]
 * extrema(x^3*y-x*y^3,[x,y])  //dog saddle at origin
 *    >> [],[]
 * extrema(x^2+y^2+z^2,x^4+y^4+z^4=1,[x,y,z])
 *    >> [[0,0,1],[0,0,-1]],[]
 * extrema(3x+3y+8z,[x^2+z^2=1,y^2+z^2=1],[x,y,z])
 *    >> [[-3/5,-3/5,-4/5]],[[3/5,3/5,4/5]]
 * extrema(2x^2+y^2,x^4-x^2+y^2=5,[x,y])
 *    >> [[0,-sqrt(5)],[0,sqrt(5)]],
 *       [[-1/2*sqrt(6),-1/2*sqrt(17)],[-1/2*sqrt(6),1/2*sqrt(17)],
 *        [1/2*sqrt(6),-1/2*sqrt(17)],[1/2*sqrt(6),1/2*sqrt(17)]]
 * extrema((3x^4-4x^3-12x^2+18)/(12*(1+4y^2)),[x,y])
 *    >> [[2,0]],[[0,0]]
 * extrema(x-y+z,[x^2+y^2+z^2=1,x+y+2z=1],[x,y,z])
 *    >> [[(-2*sqrt(70)+7)/42,(4*sqrt(70)+7)/42,(-sqrt(70)+14)/42]],
 *       [[(2*sqrt(70)+7)/42,(-4*sqrt(70)+7)/42,(sqrt(70)+14)/42]]
 * extrema(ln(x)+2*ln(y)+3*ln(z)+4*ln(u)+5*ln(v),x+y+z+u+v=1,[x,y,z,u,v])
 *    >> [],[[1/15,2/15,1/5,4/15,1/3]]
 * extrema(x*y*z,-2x^3+15x^2*y+11y^3-24y=0,[x,y,z])
 *    >> [],[]
 * extrema(x+y-exp(x)-exp(y)-exp(x+y),[x,y])
 *    >> [],[[ln(-1/2*(-sqrt(5)+1)),ln(-1/2*(-sqrt(5)+1))]]
 * extrema(x^2*sin(y)-4*x,[x,y])    // has two saddle points
 *    >> [],[]
 * extrema((1+y*sinh(x))/(1+y^2+tanh(x)^2),[x,y])
 *    >> [],[[0,0]]
 * extrema((1+y*sinh(x))/(1+y^2+tanh(x)^2),y=x^2,[x,y])
 *    >> [[1.42217627369,2.02258535346]],[[8.69443642205e-16,7.55932246971e-31]]
 * extrema(x^2*y^2,[x^3*y-2*x^2+3x-2y^2=0],[x=0..inf,y])
 *    >> [[3/2,0]],[[0.893768095046,-0.5789326693514]]
 * extrema(alpha*x^2+beta*y^2+gamma*z^2,[a1*x+a2*y+a3*z=c,b1*x+b2*y+b3*z=d],[x,y,z])
 * assume(a>0):;extrema(sqrt(a+x+y)/(1+y+sqrt(a+x)),[x=0..inf,y=0..inf])
 *    >> [[],[]]    // saddle point: [x=(-a+1/4),y=1]
 * extrema((1+(x+y+1)^2*(19-14x+3x^2-14y+6x*y+3y^2))*(30+(2x-3y)^2*(18-32x+12x^2+48y-36x*y+27y^2)),[x,y])
 *    >> [[[0,-1],[6/5,4/5],[9/5,1/5],[-3/5,-2/5]],[[4/5,1/5]]]   // has four saddle points
 */
gen _extrema(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    gen f=gv[0];
    vecteur constr;
    int ord=5; // will not compute the derivatives of order higher than ORD
    int ngv=gv.size();
    bool approx_hp=false; // use nlpsolve to determine images of homogeneous polynomials on spheres
    if (gv.back()==at_lagrange) {
        ord=0; // use Lagrange method
        --ngv;
    } else if (is_equal(gv.back())) {
        vecteur &v=*gv.back()._SYMBptr->feuille._VECTptr;
        if (v[0]==at_order && is_integer(v[1])) {
            if ((ord=v[1].val)<1)
                return generr("Expected a positive integer");
            --ngv;
        }
    }
    if (ngv<2 || ngv>3)
        return generr("Wrong number of input arguments");
    // get the variables
    vecteur vars,initial,ineq;
    // parse variables and their ranges, if given
    if (!parse_variables(gv[ngv-1],vars,ineq,initial,contextptr))
        return generrtype("Invalid specification of variables");
#if 0
    if (vars.size()>20)
        return generrdim("Too many variables");
#endif
    if (!initial.empty() && initial.size()!=vars.size())
        return generrdim("Invalid initial point specification");
    if (ngv==3) {
        // get the constraint(s)
        if (gv[1].type==_VECT)
            constr=*gv[1]._VECTptr;
        else
            constr=vecteur(1,gv[1]);
    }
    if (ord==0 && constr.empty())
        return generr("At least one constraint is required for Lagrange method");
    for (iterateur it=constr.begin();it!=constr.end();++it) {
        if (is_equal(*it))
            *it=equal2diff(*it);
        else if (which_ineq(*it)!=0)
            return generr("Inequality constraints are not supported");
    }
    vecteur parm=*_lname(makevecteur(f,constr),contextptr)._VECTptr;
    for (int i=parm.size();i-->0;) {
        if (contains(vars,parm[i]))
            parm.erase(parm.begin()+i);
    }
    if (!initial.empty() && !parm.empty())
        return generr("Cannot perform local search with symbolic constants");
    gen_map parm_asmp;
    if (!get_parameter_assumptions(parm,parm_asmp,contextptr))
        return generr("Symbolic constants must be real");
    vecteur tmpvars=make_temp_vars(vars,ineq,contextptr);
    f=subst(f,vars,tmpvars,false,contextptr);
    constr=subst(constr,vars,tmpvars,false,contextptr);
    ineq=subst(ineq,vars,tmpvars,false,contextptr);
    vecteur cpts;
    if (!extrema(f,constr,ineq,tmpvars,initial,parm_asmp,ord,approx_hp,cpts,contextptr)) {
        ;
    }
    for (int i=cpts.size();i-->0;) {
        if (!is_constant_wrt_vars(cpts[i]._VECTptr->front(),tmpvars,contextptr))
            cpts.erase(cpts.begin()+i); // not a strict local extremum, discard it
    }
    if (ord==1) { // return the list of critical points
        vecteur cv;
        for (const_iterateur it=cpts.begin();it!=cpts.end();++it) {
            cv.push_back(it->_VECTptr->front());
        }
        return cv;
    }
    // return sequence of minima and maxima in separate lists and report non- or possible extrema
    vecteur minv(0),maxv(0);
    for (const_iterateur it=cpts.begin();it!=cpts.end();++it) {
        gen dispt(vars.size()==1?symb_equal(vars[0],it->_VECTptr->front()):_zip(makesequence(at_equal,vars,it->_VECTptr->front()),contextptr));
        int cls=it->_VECTptr->at(1).val;
        const gen_map &asmp=*it->_VECTptr->back()._MAPptr;
        switch(cls) {
        case _CPCLASS_MIN:
            minv.push_back(it->_VECTptr->front());
            if (!asmp.empty())
                *logptr(contextptr) << dispt << ": " << gettext("local minimum");
            break;
        case _CPCLASS_MAX:
            maxv.push_back(it->_VECTptr->front());
            if (!asmp.empty())
                *logptr(contextptr) << dispt << ": " << gettext("local maximum");
            break;
        case _CPCLASS_SADDLE:
            *logptr(contextptr) << dispt << ": " << gettext("saddle point");
            break;
        case _CPCLASS_POSSIBLE_MIN:
        case _CPCLASS_POSSIBLE_MAX:
            *logptr(contextptr) << dispt << ": " << gettext("indeterminate critical point");
            break;
        case _CPCLASS_UNDECIDED:
            *logptr(contextptr) << dispt << ": " << gettext("unclassified critical point");
            break;
        }
#if 1
        if (!asmp.empty()) {
            bool msg=false;
            for (gen_map::const_iterator it=asmp.begin();it!=asmp.end();++it) {
                for (const_iterateur jt=it->second._VECTptr->begin();jt!=it->second._VECTptr->end();++jt) {
                    if (!msg) {
                        *logptr(contextptr) << gettext(" under condition(s) ");
                        msg=true;
                    }
                    if (jt!=it->second._VECTptr->begin())
                        *logptr(contextptr) << ",";
                    *logptr(contextptr) << subst(*jt,tmpvars,vars,false,contextptr);
                }
            }
        }
#endif
        if ((cls!=_CPCLASS_MIN && cls!=_CPCLASS_MAX) || !asmp.empty())
            *logptr(contextptr) << "\n";
    }
    return change_subtype(makevecteur(minv,maxv),_LIST__VECT);
}
static const char _extrema_s []="extrema";
static define_unary_function_eval (__extrema,&_extrema,_extrema_s);
define_unary_function_ptr5(at_extrema,alias_at_extrema,&__extrema,0,true)

/* Compute the value of expression F(X) (or |F(X)| if ABSOLUTE is true) for X=A. */
gen compf(const gen &f,identificateur &x,gen &a,bool absolute,GIAC_CONTEXT) {
    gen val(subst(f,x,a,false,contextptr));
    return _evalf(absolute?_abs(val,contextptr):val,contextptr);
}

/* find zero of expression F(X) for X in [A,B] using Brent solver. */
gen find_zero(const gen &f,identificateur &x,gen &a,gen &b,GIAC_CONTEXT) {
    gen var(symb_equal(x,symb_interval(a,b)));
    gen tmpsol=_fsolve(makesequence(f,var,_BRENT_SOLVER),contextptr);
    if (tmpsol.type==_VECT) {
        vecteur &sol=*tmpsol._VECTptr;
        return sol.empty()?(a+b)/2:sol[0];
    } else return (a+b)/2;
}

/* Find the point of maximum/minimum of unimodal expression f(x) in [a,b] using the
 * golden-section search. */
gen find_peak(const gen &f,identificateur &x,gen &a_orig,gen &b_orig,GIAC_CONTEXT) {
    gen a(a_orig),b(b_orig);
    gen c(b-(b-a)/GOLDEN_RATIO),d(a+(b-a)/GOLDEN_RATIO);
    while (is_strictly_greater(_abs(c-d,contextptr),epsilon(contextptr),contextptr)) {
        gen fc(compf(f,x,c,true,contextptr)),fd(compf(f,x,d,true,contextptr));
        if (is_strictly_greater(fc,fd,contextptr))
            b=d;
        else
            a=c;
        c=b-(b-a)/GOLDEN_RATIO;
        d=a+(b-a)/GOLDEN_RATIO;
    }
    return (a+b)/2;
}

/* Compute n Chebyshev nodes in [a,b]. */
vecteur chebyshev_nodes(gen &a,gen &b,int n,GIAC_CONTEXT) {
    vecteur nodes(1,a);
    for (int i=1;i<=n;++i) {
        nodes.push_back(_evalf((a+b)/2+(b-a)*symbolic(at_cos,((2*i-1)*cst_pi/(2*n)))/2,contextptr));
    }
    nodes.push_back(b);
    return *_sort(nodes,contextptr)._VECTptr;
}

/*
 * Implementation of Remez method for minimax polynomial approximation of a
 * continuous bounded function, which is not necessary differentiable in all
 * points of (a,b).
 *
 * Source: http://homepages.rpi.edu/~tasisa/remez.pdf
 *
 * Usage
 * ^^^^^
 *      minimax(expr,var=a..b,n,[opts])
 *
 * Parameters
 * ^^^^^^^^^^
 *      - expr            : expression to be approximated on [a,b]
 *      - var             : variable
 *      - a,b             : bounds of the function domain
 *      - n               : degree of the minimax approximation polynomial
 *      - opts (optional) : sequence of options
 *
 * This function uses 'compf', 'find_zero' and 'find_peak'. It does not use
 * derivatives to determine points of local extrema of error function, but
 * instead implements the golden search algorithm to find these points in the
 * exchange phase of Remez method.
 *
 * The returned polynomial may have degree lower than n, because the latter is
 * decremented during execution of the algorithm if there is no unique solution
 * for coefficients of a nth degree polynomial. After each decrement, the
 * algorithm is restarted. If the degree of resulting polynomial is m<n, it
 * means that polynomial of degree between m and n cannot be obtained by using
 * this implementation.
 *
 * In 'opts' one may specify 'limit=<posint>' which limits the number of
 * iterations. By default, it is unlimited.
 *
 * Be aware that, in some cases, the result with high n may be unsatisfying,
 * producing larger error than the polynomials for smaller n. This happens
 * because of the rounding errors. Nevertheless, a good approximation of an
 * "almost" smooth function can usually be obtained with n less than 30. Highly
 * oscillating functions containing sharp cusps and spikes will probably be
 * approximated poorly.
 *
 * Examples
 * ^^^^^^^^
 * minimax(x*exp(-x),x=0..10,24)
 * minimax(x*sin(x),x=0..10,25)
 * minimax(ln(2+x-sin(x)^2),x=0..2*pi,20)
 * minimax(cos(x^2-x+1),x=-2..2,40)
 * minimax(atan(x),x=-5..5,25)
 * minimax(tanh(sin(9x)),x=-1/2..1/2,40)
 * minimax(abs(x),x=-1..1,20)
 * minimax(abs(x)*sqrt(abs(x)),x=-2..2,15)
 * minimax(min(1/cosh(3*sin(10x)),sin(9x)),x=-0.3..0.4,25)
 * minimax(when(x==0,0,exp(-1/x^2)),x=-1..1,30)
 */
gen _minimax(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    if (gv.size()<3)
        return generr("Too few arguments");
    if (!is_equal(gv[1]) || !is_integer(gv[2]))
        return generrtype("Expected an equality and an integer");
    // detect parameters
    vecteur s(*gv[1]._SYMBptr->feuille._VECTptr);
    if (s[0].type!=_IDNT || !s[1].is_symb_of_sommet(at_interval))
        return generrtype("Expected x=a..b as the second argument");
    identificateur x(*s[0]._IDNTptr);
    s=*s[1]._SYMBptr->feuille._VECTptr;
    gen a(_evalf(s[0],contextptr)),b(_evalf(s[1],contextptr));
    if (!is_strictly_greater(b,a,contextptr))
        return generrtype("Invalid bounds in the second argument");
    const gen &f=gv[0];
    int n=gv[2].val;
    gen threshold(1.02);  // threshold for stopping criterion
    // detect options
    int limit=0;
    //bool poly=true;
    for (const_iterateur it=gv.begin()+3;it!=gv.end();++it) {
        if (is_equal(*it)) {
            vecteur &p=*it->_SYMBptr->feuille._VECTptr;
            if (p[0]==at_limit) {
                if (!p[1].is_integer() || p[1].val<=0)
                    return generrtype("Expected a positive integer");
                limit=p[1].val;
                continue;
            }
        }
        return generr("Invalid option");
    }
    // create Chebyshev nodes to start with
    vecteur nodes(chebyshev_nodes(a,b,n,contextptr)),sol;
    gen p,best_p,best_emax,emax,emin,tmpsol;
    int iteration_count=0;
    while (true) { // iterate the algorithm
        iteration_count++;
        if (n<1 || (limit>0 && iteration_count>limit))
            break;
        // compute polynomial p
        matrice m;
        vecteur fv;
        for (int i=0;i<n+2;++i) {
            fv.push_back(_evalf(subst(f,x,nodes[i],false,contextptr),contextptr));
            vecteur r;
            for (int j=0;j<n+1;++j) {
                r.push_back(j==0?gen(1):pow(nodes[i],j));
            }
            r.push_back(pow(gen(-1),i));
            m.push_back(r);
        }
        sol.clear();
        if (!is_zero(mdet(m,contextptr))) {
            try {
                tmpsol=_LSQ(makesequence(m,fv),contextptr);
                if (ckmatrix(tmpsol))
                    sol=*mtran(*tmpsol._VECTptr).front()._VECTptr;
            } catch (const std::runtime_error &e) { ; }
        }
        if (sol.empty()) {
            // decrease n and start over
            nodes=chebyshev_nodes(a,b,--n,contextptr);
            continue;
        }
        p=gen(0);
        for (int i=0;i<n+1;++i) {
            p+=sol[i]*pow(x,i);
        }
        // compute the error function and its zeros
        gen e(f-p);
        vecteur zv(1,a);
        for (int i=0;i<n+1;++i) {
            zv.push_back(find_zero(e,x,nodes[i],nodes[i+1],contextptr));
        }
        zv.push_back(b);
        // remez exchange:
        // determine points of local extrema of error function e
        vecteur ev(n+2,0);
        for (int i=0;i<n+2;++i) {
            if (i>0 && i<n+1) {
                nodes[i]=find_peak(e,x,zv[i],zv[i+1],contextptr);
                ev[i]=compf(e,x,nodes[i],true,contextptr);
                continue;
            }
            gen e1(compf(e,x,zv[i],true,contextptr)),e2(compf(e,x,zv[i+1],true,contextptr));
            if (is_greater(e1,e2,contextptr)) {
                nodes[i]=zv[i];
                ev[i]=e1;
            } else {
                nodes[i]=zv[i+1];
                ev[i]=e2;
            }
        }
        // compute minimal and maximal absolute error
        emin=_min(ev,contextptr);
        emax=_max(ev,contextptr);
        if (is_exactly_zero(best_emax) || is_strictly_greater(best_emax,emax,contextptr)) {
            best_p=p;
            best_emax=emax;
        }
        // emin >= E is required to continue, also check
        // if the threshold is reached, i.e. the difference
        // between emax and emin is at least fifty times
        // smaller than emax
        if (is_strictly_greater(sol.back(),emin,contextptr) ||
                is_greater(threshold*emin,emax,contextptr)) {
            break;
        }
    }
    *logptr(contextptr) << gettext("max. absolute error: ") << best_emax << "\n";
    return best_p;
}
static const char _minimax_s []="minimax";
static define_unary_function_eval (__minimax,&_minimax,_minimax_s);
define_unary_function_ptr5(at_minimax,alias_at_minimax,&__minimax,0,true)

/*
 *
 ******* TPROB CLASS IMPLEMENTATION *******
 * 
 */

/* The constructor.
 * - s: supply vector
 * - d: demand vector
 * - m: a symbol representing infinite cost */
tprob::tprob(const vecteur &s,const vecteur &d,const gen &m,GIAC_CONTEXT) {
    eps=epsilon(contextptr)/(s.size()*d.size()+1);
    ctx=contextptr;
    supply=s;
    demand=d;
    M=m;
}
/* North-West-Corner method giving the initial feasible solution to the
 * transportatiom problem with given supply and demand vectors. It handles
 * degeneracy cases (assignment problems, for example, always have degenerate
 * solutions). */
void tprob::north_west_corner(matrice &feas) {
    int m=supply.size(),n=demand.size();
    feas=zero_mat(m,n,ctx);
    int i=0,j=0;
    while (i<m && j<n) {
        const gen &s=supply[i],&d=demand[j];
        gen u,v;
        for (int k=0;k<i;++k) {
            v+=_epsilon2zero(feas[k][j],ctx);
        }
        for (int k=0;k<j;++k) {
            u+=_epsilon2zero(feas[i][k],ctx);
        }
        gen a=min(s-u,d-v,ctx);
        feas[i]._VECTptr->at(j)=a;
        int k=i+j;
        if (u+a==s)
            ++i;
        if (v+a==d)
            ++j;
        if (i<m && j<n && i+j-k==2) // avoid degeneracy
            feas[i-1]._VECTptr->at(j)=eps;
    }
}
/* Stepping stone path method for determining a closed path "jumping" from one
 * positive element of X to another in the same row or column. */
tprob::ipairs tprob::stepping_stone_path(ipairs &path_orig,const matrice &X) {
    ipairs path(path_orig);
    int I=path.back().first,J=path.back().second;
    int m=X.size(),n=X.front()._VECTptr->size();
    if (path.size()>1 && path.front().second==J)
        return path;
    bool hrz=path.size()%2==1;
    for (int i=0;i<(hrz?n:m);++i) {
        int cnt=0;
        for (ipairs::const_iterator it=path.begin();it!=path.end();++it) {
            if ((hrz && it->second==i) || (!hrz && it->first==i))
                ++cnt;
        }
        if (cnt<2 && !is_exactly_zero(X[hrz?I:i][hrz?i:J])) {
            path.push_back(make_pair(hrz?I:i,hrz?i:J));
            ipairs fullpath(stepping_stone_path(path,X));
            if (!fullpath.empty())
                return fullpath;
            path.pop_back();
        }
    }
    return ipairs(0);
}
/* Implementation of the MODI (modified ditribution) method.
 * It handles forbidden routes and degenerate solutions.
 * - P_orig: the cost matrix
 * - X: the solution matrix (must be initialized to a feasible solution)
 * After the function returns, X will contain an optimal solution. */
void tprob::modi(const matrice &P_orig,matrice &X) {
    matrice P;
    copy_matrice(P_orig,P);
    int m=X.size(),n=X.front()._VECTptr->size();
    vecteur u(m),v(n);
    if (M.type==_IDNT)
        P=subst(P,M,_inv(eps,ctx),false,ctx);
    for (int i=0;i<m;++i)
        u[i]=i==0?gen(0):temp_symb("u",++var_index,ctx);
    for (int j=0;j<n;++j)
        v[j]=temp_symb("v",++var_index,ctx);
    vecteur vars(mergevecteur(vecteur(u.begin()+1,u.end()),v));
    while (true) {
        vecteur eqv;
        for (int i=0;i<m;++i) {
            for (int j=0;j<n;++j) {
                if (!is_exactly_zero(X[i][j]))
                    eqv.push_back(u[i]+v[j]-P[i][j]);
            }
        }
        vecteur sol=*_linsolve(makesequence(eqv,vars),ctx)._VECTptr;
        vecteur U(1,0),V(sol.begin()+m-1,sol.end());
        U=mergevecteur(U,vecteur(sol.begin(),sol.begin()+m-1));
        gen cmin(0);
        bool optimal=true;
        int I,J;
        for (int i=0;i<m;++i) {
            for (int j=0;j<n;++j) {
                if (is_exactly_zero(X[i][j])) {
                    gen c(P[i][j]-U[i]-V[j]);
                    if (is_strictly_greater(cmin,c,ctx)) {
                        cmin=c;
                        optimal=false;
                        I=i;
                        J=j;
                    }
                }
            }
        }
        if (optimal)
            break;
        ipairs path;
        path.push_back(make_pair(I,J));
        path=stepping_stone_path(path,X);
        gen d(X[path.at(1).first][path.at(1).second]);
        for (ipairs::const_iterator it=path.begin()+3;it<path.end();it+=2) {
            d=min(d,X[it->first][it->second],ctx);
        }
        for (int i=0;i<int(path.size());++i) {
            gen &Xij=X[path.at(i).first]._VECTptr->at(path.at(i).second);
            gen x(Xij+(i%2?-d:d));
            bool has_zero=false;
            for (ipairs::const_iterator it=path.begin();it!=path.end();++it) {
                if (is_exactly_zero(X[it->first][it->second])) {
                    has_zero=true;
                    break;
                }
            }
            if ((!is_exactly_zero(x) && is_strictly_greater(fraction(1,2),x,ctx)) ||
                    (is_exactly_zero(x) && has_zero))
                x=eps;
            Xij=x;
        }
    }
    X=*exact(_epsilon2zero(X,ctx),ctx)._VECTptr;
}
void tprob::solve(const matrice &cost_matrix,matrice &sol) {
    north_west_corner(sol);
    modi(cost_matrix,sol);
}
/* END OF TPROB CLASS */

/*
 * Function 'tpsolve' solves a transportation problem using MODI method.
 *
 * Usage
 * ^^^^^
 *      tpsolve(supply,demand,cost_matrix)
 *
 * Parameters
 * ^^^^^^^^^^
 *      - supply      : source capacity (vector of m positive integers)
 *      - demand      : destination demand (vector of n positive integers)
 *      - cost_matrix : real matrix C=[c_ij] of type mXn where c_ij is cost of
 *                      transporting an unit from ith source to jth destination
 *                      (a nonnegative number)
 *
 * Supply and demand vectors should contain only positive integers. Cost matrix
 * must be consisted of nonnegative real numbers, which do not have to be
 * integers. There is a possibility of adding a certain symbol to cost matrix,
 * usually M, to indicate the "infinite cost", effectively forbidding the
 * transportation on a certain route. The notation of the symbol may be chosen
 * arbitrarily, but must be used consistently within a single problem.
 *
 * The return value is a sequence of total (minimal) cost and matrix X=[x_ij]
 * of type mXn where x_ij is equal to number of units which have to be shipped
 * from ith source to jth destination, for all i=1,2,..,m and j=1,2,..,n.
 *
 * This function uses 'north_west_corner' to determine initial feasible
 * solution and then applies MODI method to optimize it (function 'modi', which
 * uses 'stepping_stone_path'). Also, it is capable of handling degeneracy of
 * the initial solution and during iterations of MODI method.
 *
 * If the given problem is not balanced, i.e. if supply exceeds demand or vice
 * versa, dummy supply/demand points will be automatically added to the
 * problem, augmenting the cost matrix with zeros. Resulting matrix will not
 * contain dummy point.
 *
 * Examples
 * ^^^^^^^^
 * Balanced transportation problem:
 *  tpsolve([12,17,11],[10,10,10,10],[[500,750,300,450],[650,800,400,600],[400,700,500,550]])
 *      >> 20200,[[0,0,2,10],[0,9,8,0],[10,1,0,0]]
 * Non-balanced transportation problem:
 *  tpsolve([7,10,8,8,9,6],[9,6,12,8,10],[[36,40,32,43,29],[28,27,29,40,38],[34,35,41,29,31],[41,42,35,27,36],[25,28,40,34,38],[31,30,43,38,40]])
 *      >> [[0,0,2,0,5],[0,0,10,0,0],[0,0,0,0,5],[0,0,0,8,0],[9,0,0,0,0],[0,6,0,0,0]]
 * Transportation problem with forbidden routes:
 *  tpsolve([95,70,165,165],[195,150,30,45,75],[[15,M,45,M,0],[12,40,M,M,0],[0,15,25,25,0],[M,0,M,12,0]])
 *      >> [[20,0,0,0,75],[70,0,0,0,0],[105,0,30,30,0],[0,150,0,15,0]]
 * Assignment problem:
 *  tpsolve([1,1,1,1],[1,1,1,1],[[10,12,9,11],[5,10,7,8],[12,14,13,11],[8,15,11,9]])
 *      >> [[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]]
 */
gen _tpsolve(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    if (gv.size()<3)
        return generr("Too few arguments");
    if (gv[0].type!=_VECT || gv[1].type!=_VECT ||
            gv[2].type!=_VECT || !ckmatrix(*gv[2]._VECTptr))
        return generrtype("Expected two vectors and a matrix");
    vecteur supply(*gv[0]._VECTptr),demand(*gv[1]._VECTptr);
    if (!is_integer_vecteur(supply,true) || !is_integer_vecteur(demand,true))
        return generrtype("Supply and demand quantites must be integers");
    matrice P(*gv[2]._VECTptr);
    if (is_strictly_greater(0,_min(_min(P,contextptr),contextptr),contextptr))
        return generr("All costs must be non-negative");
    vecteur sy(*_lname(P,contextptr)._VECTptr);
    if (sy.size()>1)
        return generr("At most one symbol is allowed in the cost matrix");
    int m=supply.size(),n=demand.size();
    if (m!=int(P.size()) || n!=int(P.front()._VECTptr->size()))
        return generrdim("Cost matrix dimensions do not match supply and demand");
    gen M(sy.size()==1 && sy[0].type==_IDNT?sy[0]:0);
    gen ts(_sum(supply,contextptr)),td(_sum(demand,contextptr));
    if (ts!=td) {
        if (is_greater(ts,td,contextptr)) {
            demand.push_back(ts-td);
            P=mtran(P);
            P.push_back(vecteur(m,0));
            P=mtran(P);
        } else {
            supply.push_back(td-ts);
            P.push_back(vecteur(n,0));
        }
    }
    matrice X;
    tprob tp(supply,demand,M,contextptr);
    tp.solve(P,X);
    if (is_strictly_greater(ts,td,contextptr)) {
        X=mtran(X);
        X.pop_back();
        X=mtran(X);
    } else if (is_strictly_greater(td,ts,contextptr))
        X.pop_back();
    gen cost(0);
    for (int i=0;i<m;++i) {
        for (int j=0;j<n;++j) {
            cost+=P[i][j]*X[i][j];
        }
    }
    if (contains(*_lname(cost,contextptr)._VECTptr,M))
        return generr("The problem has no feasible solution");
    return makesequence(cost,X);
}
static const char _tpsolve_s []="tpsolve";
static define_unary_function_eval (__tpsolve,&_tpsolve,_tpsolve_s);
define_unary_function_ptr5(at_tpsolve,alias_at_tpsolve,&__tpsolve,0,true)

void add_identifiers(const gen &source,vecteur &dest,GIAC_CONTEXT) {
    vecteur v(*_lname(source,contextptr)._VECTptr);
    for (const_iterateur it=v.begin();it!=v.end();++it) {
        if (_eval(*it,contextptr).type==_IDNT && !contains(dest,*it))
            dest.push_back(*it);
    }
    dest=sort_identifiers(dest,contextptr);
}

/*
 *
 ******* NLP PROBLEM CLASS IMPLEMENTATION ********
 *
 */
nlp_problem::optimum::optimum() {
    res=_NLP_PENDING;
    f=undef;
}
/* Initialize a NLP problem: MIN f(x) SUBJECT TO list of (in)equality constraints g(x), x in R^n. */
nlp_problem::nlp_problem(const gen &f,const vecteur &g,const vecteur &x,GIAC_CONTEXT) : ctx(contextptr)
{
    _msg_level=1;
    _initialized=_presolved=_infeas=_convex_obj=_convex_constr=_debug=_iter_data=false;
    _linobj=_linconstr=_compact=_smooth_obj=_smooth_constr=false;
    feasible_region=NULL;
    vars=obj_vars=x;
    obj=f;
    best_obj_val=undef;
    obj_shift=0;
    _penalty_scheme=_fevalc=_iter_count=0;
    _is_intvar.resize(vars.size(),false);
    _sol_stat.resize(6);
    debug("Constructing NLP problem...");
    reset_timer();
    /* Store equality and inequality constraints in G separately, ignore other elements of G.
       Inequalities are stored as nonnegativity constraints. */
    for (const_iterateur it=g.begin();it!=g.end();++it) {
        if (is_equal(*it)) {
            eq_cons.push_back(equal2diff(*it));
        } else switch (which_ineq(*it)) {
            case 1:
                ineq_cons.push_back(it->_SYMBptr->feuille._VECTptr->back()-it->_SYMBptr->feuille._VECTptr->front());
                break;
            case 2:
                ineq_cons.push_back(it->_SYMBptr->feuille._VECTptr->front()-it->_SYMBptr->feuille._VECTptr->back());
                break;
            default:
                break;
        }
    }
}
nlp_problem::~nlp_problem() {
    debug("Deleting NLP problem...");
    if (feasible_region!=NULL)
        delete feasible_region;
    purge_bounded_vars();
}
/* Presolve (optional) - Compute gradients - Determine bounding rect - Check convexity */
void nlp_problem::initialize(int method,const meth_parm &parm) {
    if (parm.presolve && !(_presolved=preprocess())) {
        _infeas=true;
        return;
    }
    int n=var_count();
    if (n>0) {
        make_bounded_vars(false);
        remove_redundant_constraints();
        if (_infeas) return;
        if (method<_NLP_NELDER_MEAD) {
            obj_gradient=*_grad(makesequence(obj,vars),ctx)._VECTptr;
            if (_smooth_obj=!has_diff(obj_gradient,vars) && !has_breaks(obj_gradient))
                debug("Objective function is differentiable");
            else debug ("Objective function is not differentiable");
            eq_jacobian.reserve(eq_cons.size());
            ineq_jacobian.reserve(ineq_cons.size());
            _linconstr=_smooth_constr=_convex_constr=true;
        } else _linconstr=_smooth_obj=_convex_obj=_smooth_constr=_convex_constr=false;
        for (const_iterateur it=eq_cons.begin();_smooth_constr && it!=eq_cons.end();++it) {
            eq_jacobian.push_back(_grad(makesequence(*it,vars),ctx));
            _smooth_constr=!has_diff(eq_jacobian.back(),vars) && !has_breaks(eq_jacobian.back());
            if (_smooth_constr && _linconstr && !is_linear_wrt_vars(*it,vars,ctx)) {
                _linconstr=_convex_constr=false;
                if (parm.convex==1)
                    warn("nonlinear equality constraint found, dropping convexity assumption");
            }
        }
        for (const_iterateur it=ineq_cons.begin();_smooth_constr && it!=ineq_cons.end();++it) {
            ineq_jacobian.push_back(_grad(makesequence(*it,vars),ctx));
            _smooth_constr=!has_diff(ineq_jacobian.back(),vars) && !has_breaks(ineq_jacobian.back());
            if (_smooth_constr && _linconstr && !is_linear_wrt_vars(*it,vars,ctx))
                _linconstr=false;
        }
        if (!_smooth_constr) {
            eq_jacobian.clear();
            ineq_jacobian.clear();
            _linconstr=_convex_constr=false;
        } else {
            if (feasible_region!=NULL)
                delete feasible_region;
            feasible_region=new region(*this,parm.tol,100);
        }
        // detect feasible region
        if (bd_rect.empty())
            bd_rect=_smooth_constr && parm.border?feasible_region->bounding_rect():mtran(makevecteur(lb,ub));
        if (_infeas)
            return;
        _compact=!has_inf_or_undef(bd_rect);
        // Is the objective function convex?
        if (_smooth_obj) {
            debug("Checking convexity of the objective...");
            _linobj=is_linear_wrt_vars(obj,vars,ctx);
            if (parm.convex>=0) {
                if (jacobian(obj_gradient,vars,obj_hessian,ctx)) {
                    _convex_obj=_have_hessian=!has_diff(obj_hessian,vars) && !has_breaks(obj_hessian);
                } else {
                    _convex_obj=_have_hessian=false;
                }
                if (parm.convex==1 && !_convex_obj)
                    warn("the objective is not twice differentiable, dropping convexity assumption");
            } else {
                make_bounded_vars(true);
                obj_hessian=zero_mat(n,n,ctx);
                if (_linobj) {
                    _convex_obj=_have_hessian=true;
                } else _convex_obj=is_expression_convex(obj_gradient,obj_hessian,_have_hessian,false);
            }
            if (!_have_hessian)
                obj_hessian.clear();
            if (_convex_obj)
                debug("Objective function is convex under the bound constraints");
        }
        // Are inequality constraints convex?
        if (_convex_constr && parm.convex==-1) {
            debug("Checking convexity of the constraints...");
            matrice constr_hessian=zero_mat(n,n,ctx);
            make_bounded_vars(true);
            bool hess_ok;
            for (const_iterateur it=ineq_jacobian.begin();_convex_constr && it!=ineq_jacobian.end();++it) {
                _convex_constr=is_numericv(*it->_VECTptr,num_mask_withfrac|num_mask_withint) ||
                                is_expression_convex(multvecteur(-1,*it->_VECTptr),constr_hessian,hess_ok,true);
            }
        } else _convex_constr=parm.convex==1;
        // Is the problem convex?
        if (_convex_obj && _convex_constr)
            msg(parm.convex==1?"Assuming that the problem is convex":"Problem is convex");
    }
    _initialized=true;
}
/* This routine requires BND_VARS to be initialized. */
bool nlp_problem::is_expression_convex(const vecteur &gradient,matrice &H,bool &hess_ok,bool lt) {
    hess_ok=true;
    int n=var_count();
    for (int i=0;hess_ok&&i<n;++i) {
        const gen &dxi=gradient[i];
        for (int j=0;j<=i;++j) {
            gen dxixj=derive(dxi,vars[j],ctx);
            if (has_inf_or_undef(dxixj) || has_diff(dxixj,vars[j]) || has_breaks(dxixj)) {
                hess_ok=false;
                break;
            }
            H[i]._VECTptr->at(j)=dxixj;
            if (!lt && i!=j)
                H[j]._VECTptr->at(i)=dxixj;
        }
    }
    return hess_ok && is_positive_semidefinite(subst(H,vars,bnd_vars,false,ctx),0.5,100,ctx);
}
/* Make bounded copies of variables using assumptions */
void nlp_problem::make_bounded_vars(bool use_bd_rect) {
    int nv=var_count();
    purge_bounded_vars();
    bnd_vars.resize(nv);
    for (int i=0;i<nv;++i) {
        bnd_vars[i]=temp_symb("xbnd",i,ctx);
        const gen &l=use_bd_rect?bd_rect[i][0]:lb[i];
        const gen &u=use_bd_rect?bd_rect[i][1]:ub[i];
        if (!is_inf(l) || !is_inf(u))
            set_assumptions(bnd_vars[i],vecteur(1,makevecteur(l,u)),vecteur(0),false,ctx);
    }
}
void nlp_problem::remove_redundant_constraints() {
    assert(bnd_vars.size()==vars.size());
    int rm=0,ch=0;
    debug("Finding redundant constraints...");
    for (int i=ineq_count();i-->0;) {
        const gen &c=ineq_cons[i];
        gen p=subst(c,vars,bnd_vars,false,ctx);
        //if (!is_real_number(p,ctx)) continue;
        //p=to_real_number(p,ctx);
        if (is_positive_safe(p,false,100,ctx)) {
            _ss << "Removing redundant constraint " << c << ""; debug();
            rm++;
            ineq_cons.erase(ineq_cons.begin()+i);
        }
        else if (is_definitely_positive(-p,ctx)) {
            _ss << "Constraint " << c << " is not feasible"; debug();
            _infeas=true;
            break;
        }
#if 0
        else if (is_positive(-p,ctx)) {
            _ss << "Replacing " << c << ">=0 by equality constraint"; debug();
            ch++;
            eq_cons.push_back(p);
            ineq_cons.erase(ineq_cons.begin()+i);
        }
#endif
    }
    if (rm>0) { _ss << gettext("Removed") << " " << rm << " " << gettext("redundant constraints"); msg(); }
    else debug("No redundant constraints found");
}
pair<int,int> nlp_problem::make_leqv_sys(const vector<int> &ind,const vector<lineq> &leqv) {
    vector<int> isect,iunit,tmp;
    for (graphe::ivector_iter jt=ind.begin();jt!=ind.end();++jt) {
        const lineq &lineq=leqv[*jt];
        if (jt==ind.begin()) {
            isect=lineq.lvars;
            iunit=lineq.lvars;
        } else {
            intersect(isect,lineq.lvars,tmp);
            isect=tmp;
            unite(iunit,lineq.lvars,tmp);
            iunit=tmp;
        }
    }
    return make_pair(isect.size(),-int(iunit.size()));
}
bool nlp_problem::eliminate_equalities(const vector<lineq> &leqv,vecteur &lsol,vecteur &leq_vars) {
    int n=leqv.size();
    graphe G(ctx,false);
    G.add_nodes(n);
    vector<int> tmp;
    for (int i=0;i<n;++i) {
        const lineq &nl1=leqv[i];
        for (int j=i+1;j<n;++j) {
            const lineq &nl2=leqv[j];
            if (intersect(nl1.lvars,nl2.nlvars,tmp)==0 && intersect(nl1.nlvars,nl2.lvars,tmp)==0)
                G.add_edge(i,j);
        }
    }
    G.find_maximal_cliques();
    graphe::ivectors mc=G.maximal_cliques();
    if (mc.empty())
        return false;
    vector<pair<pair<int,int>,int> > mca,rem;
    for (int i=0;i<int(mc.size());++i) {
        vector<int> &clq=mc[i];
        int sz=clq.size();
        pair<int,int> szp=make_leqv_sys(clq,leqv);
        while (sz>-szp.second) {
            rem.clear();
            for (int j=clq.size();j-->0;) {
                vector<int> clqj=clq;
                clqj.erase(clqj.begin()+j);
                rem.push_back(make_pair(make_leqv_sys(clqj,leqv),j));
            }
            std::sort(rem.begin(),rem.end());
            clq.erase(clq.begin()+rem.rbegin()->second);
            szp=rem.rbegin()->first;
            --sz;
        }
        if (sz==0) {
            mc.erase(mc.begin()+i);
            --i;
        } else mca.push_back(make_pair(szp,i));
    }
    if (mca.empty())
        return false;
    std::sort(mca.begin(),mca.end());
    const graphe::ivector &c=mc[mca.rbegin()->second];
    leq_vars.clear();
    vecteur eqv;
    vector<int> rc;
    for (graphe::ivector_iter it=c.begin();it!=c.end();++it) {
        const lineq &lineq=leqv[*it];
        eqv.push_back(lineq.e);
        rc.push_back(lineq.index);
        for (vector<int>::const_iterator jt=lineq.lvars.begin();jt!=lineq.lvars.end();++jt) {
            const gen &v=vars[*jt];
            if (!contains(leq_vars,v))
                leq_vars.push_back(v);
        }
    }
    lsol.clear();
    vecteur sol=*_solve(makesequence(eqv,leq_vars),ctx)._VECTptr;
    if (!sol.empty()) {
        assert(sol.front().type==_VECT);
        lsol=*sol.front()._VECTptr;
        std::sort(rc.begin(),rc.end());
        for (vector<int>::const_reverse_iterator it=rc.rbegin();it!=rc.rend();++it) {
            eq_cons.erase(eq_cons.begin()+*it);
        }
    }
    return true;
}
bool nlp_problem::subs_fxvars(const fxvars &fv) {
    obj=subst(obj,fv.names,fv.values,false,ctx);
    for (int i=eq_cons.size();i-->0;) {
        gen &c=eq_cons[i];
        c=simp(subst(c,fv.names,fv.values,false,ctx),ctx);
        if (is_zero(c,ctx))
            eq_cons.erase(eq_cons.begin()+i);
        else if (is_real_number(c,ctx))
            return false;
    }
    for (int i=ineq_cons.size();i-->0;) {
        gen &c=ineq_cons[i];
        c=simp(subst(c,fv.names,fv.values,false,ctx),ctx);
        if (is_positive_safe(c,false,100,ctx))
            ineq_cons.erase(ineq_cons.begin()+i);
        else if (is_real_number(c,ctx))
            return false;
    }
    return true;
}
/* Return the list of constraints (optionally including variable bounds).
 * If all_geq is true, then the list contains lhs of constraints
 * written in form f_i(x) >= 0 (used by COBYLA). */
vecteur nlp_problem::constraints(bool incl_bnds,bool all_geq) const {
    vecteur res;
    int n=var_count();
    res.reserve(2*eq_cons.size()+ineq_cons.size()+2*n);
    if (!all_geq) {
        vecteur eqs=eq_cons.empty()?vecteur(0):*_zip(makesequence(at_equal,eq_cons,vecteur(eq_cons.size(),0)),ctx)._VECTptr;
        vecteur ineqs=ineq_cons.empty()?vecteur(0):*_zip(makesequence(at_superieur_egal,ineq_cons,vecteur(ineq_cons.size(),0)),ctx)._VECTptr;
        res=mergevecteur(eqs,ineqs);
    } else {
        if (!eq_cons.empty())
            res=mergevecteur(eq_cons,multvecteur(-1,eq_cons));
        if (!ineq_cons.empty())
            res=mergevecteur(res,ineq_cons);
    }
    if (!incl_bnds)
        return res;
    for (int i=0;i<n;++i) {
        if (!is_inf(lb[i]))
            res.push_back(all_geq?vars[i]-lb[i]:symb_superieur_egal(vars[i],lb[i]));
        if (!is_inf(ub[i]))
            res.push_back(all_geq?ub[i]-vars[i]:symb_inferieur_egal(vars[i],ub[i]));
    }
    return res;
}
/* Return lists of variables and the respective values. */
void nlp_problem::get_vars_and_values(const set<pair<int,gen> > &lst,vecteur &v,vecteur &w) const {
    set<pair<int,gen> >::const_iterator it=_fixed_vars.begin(),itend=_fixed_vars.end();
    v.clear(); w.clear();
    v.reserve(_fixed_vars.size());
    w.reserve(_fixed_vars.size());
    for (;it!=itend;++it) {
        v.push_back(vars[it->first]);
        w.push_back(it->second);
    }
}
/* Return G where fixed variables are substituted with the respective values. */
gen nlp_problem::subs_fixed_vars(const gen &g) const {
    if (_fixed_vars.empty()) return g;
    vecteur fv,fval;
    get_vars_and_values(_fixed_vars,fv,fval);
    return subst(g,fv,fval,false,ctx);
}
vecteur nlp_problem::subs_fixed_vars(const vecteur &v) const {
    if (_fixed_vars.empty()) return v;
    vecteur fv,fval;
    get_vars_and_values(_fixed_vars,fv,fval);
    return subst(v,fv,fval,false,ctx);
}
/* Find inactive equality constraints among those in EQ.
 * Return false if there are violated constraints, else return true. */
bool nlp_problem::find_inactive_eq(const vecteur &eq) {
    const_iterateur it=eq.begin(),itend=eq.end();
    _inactive_eq_indices.clear();
    int i=0;
    for (;it!=itend;++it,++i) {
        if (!is_real_number(*it,ctx)) continue;
        if (is_zero(to_real_number(*it,ctx),ctx))
            _inactive_eq_indices.insert(i);
        else {
            _inactive_eq_indices.clear();
            return false;
        }
    }
    return true;
}
/* Find inactive inequality constraints among those in INEQ.
 * Return false if there are violated constraints, else return true. */
bool nlp_problem::find_inactive_ineq(const vecteur &ineq) {
    const_iterateur it=ineq.begin(),itend=ineq.end();
    _inactive_ineq_indices.clear();
    int i=0;
    for (;it!=itend;++it,++i) {
        if (is_positive_safe(*it,false,100,ctx))
            _inactive_ineq_indices.insert(i);
        else if (is_definitely_positive(-*it,ctx)) {
            _inactive_ineq_indices.clear();
            return false;
        }
    }
    return true;
}
/* Preprocess the problem. */
bool nlp_problem::preprocess() {
    int pass=0,oldvarcount=var_count(),rmvars=0;
    fxvars fv(vars.size());
    for (int j=vars.size();j-->0;) { // remove fixed variables
        if (is_zero(ub[j]-lb[j],ctx)) {
            fv.names.push_back(vars[j]);
            fv.values.push_back((lb[j]+ub[j])/2);
            fv.indices.push_back(j);
            vars.erase(vars.begin()+j);
            lb.erase(lb.begin()+j);
            ub.erase(ub.begin()+j);
            if (!initp.empty())
                initp=*_delcols(makesequence(initp,j),ctx)._VECTptr;
            if (!bd_rect.empty())
                bd_rect.erase(bd_rect.begin()+j);
            ++rmvars;
        } else if (is_strictly_greater(lb[j],ub[j],ctx))
            return false; // infeasible
    }
    if (rmvars>0) {
        _ss << rmvars << " " << gettext("variable(s) fixed by preprocessor");
        msg();
    }
    if (!fv.names.empty()) {
        if (!subs_fxvars(fv))
            return false;
        fx_subs.push(fv);            
    }
    bool changed;
    do { // substitute variables
        changed=false;
        fxvars fxv(vars.size());
        vector<lineq> leqv;
        for (const_iterateur ct=eq_cons.begin();ct!=eq_cons.end();++ct) {
            gen d=exact(*ct,ctx),a,b;
            vecteur lin;
            for (const_iterateur it=vars.begin();it!=vars.end();++it) {
                if (!is_intvar(it-vars.begin()) && is_linear_wrt(d,*it,a,b,ctx) &&
                        !is_exactly_zero(a) && is_constant_wrt_vars(a,vars,ctx))
                    lin.push_back(*it);
            }
            if (lin.empty())
                continue;
            vecteur v=*_lname(d,ctx)._VECTptr;
            for (const_iterateur it=lin.begin();it!=lin.end();++it) {
                v.erase(find(v.begin(),v.end(),*it));
            }
            lineq lineq(d,ct-eq_cons.begin());
            for (const_iterateur it=lin.begin();it!=lin.end();++it) {
                lineq.lvars.push_back(indexof(*it,vars));
            }
            for (const_iterateur it=v.begin();it!=v.end();++it) {
                lineq.nlvars.push_back(indexof(*it,vars));
            }
            std::sort(lineq.lvars.begin(),lineq.lvars.end());
            std::sort(lineq.nlvars.begin(),lineq.nlvars.end());
            leqv.push_back(lineq);
        }
        if (leqv.empty())
            continue;
        vecteur lsol,leq_vars;
        if (!eliminate_equalities(leqv,lsol,leq_vars))
            continue;
        if (!lsol.empty()) {
            for (const_iterateur it=lsol.begin();it!=lsol.end();++it) {
                int i,j=it-lsol.begin();
                if (is_constant_wrt(*it,leq_vars[j],ctx)) {
                    fxv.names.push_back(leq_vars[j]);
                    fxv.values.push_back(*it);
                    fxv.indices.push_back(i=indexof(leq_vars[j],vars));
                    if (!is_inf(lb[i]))
                        ineq_cons.push_back(*it-lb[i]);
                    if (!is_inf(ub[i]))
                        ineq_cons.push_back(ub[i]-*it);
                }
            }
            for (int j=vars.size();j-->0;) {
                if (find(fxv.indices.begin(),fxv.indices.end(),j)!=fxv.indices.end()) {
                    vars.erase(vars.begin()+j);
                    lb.erase(lb.begin()+j);
                    ub.erase(ub.begin()+j);
                    _is_intvar.erase(_is_intvar.begin()+j);
                    if (!initp.empty())
                        initp=*_delcols(makesequence(initp,j),ctx)._VECTptr;
                    if (!bd_rect.empty())
                        bd_rect.erase(bd_rect.begin()+j);
                }
            }
        } else return false;
        if (changed=!fxv.names.empty()) {
            if (!subs_fxvars(fxv))
                return false;
            fx_subs.push(fxv);
        }
    } while (changed && (++pass>0));
    if (var_count()<oldvarcount-rmvars) {
        _ss << oldvarcount-var_count()-rmvars << " variable(s) substituted by preprocessor";
        msg();
    }
    return true;
}
/* Postprocess the (optimal) point x0 (may be vector of points and/or
 * vectors of points) for returning to user.
 * This works only once! (subs info is destroyed) */
void nlp_problem::postprocess_point(gen &x0,const vecteur &oldvars,const fxvars &fx,vecteur &csol) {
    if (x0.type!=_VECT)
        return;
    vecteur &pt=*x0._VECTptr;
    if (is_numericv(pt)) {
        vecteur fxvals=oldvars.empty()?fx.values:subst(fx.values,oldvars,pt,false,ctx);
        for (int i=fx.indices.size();i-->0;) csol[fx.indices[i]]=fxvals[i];
        int i=0;
        for (const_iterateur it=pt.begin();it!=pt.end();++it) {
            while (!is_undef(csol[i])) ++i;
            csol[i]=*it;
        }
        x0=csol;
    } else for (iterateur it=pt.begin();it!=pt.end();++it) {
        postprocess_point(*it,oldvars,fx,csol);
    }
}
void nlp_problem::postprocess(gen &x0) {
    int i;
    while (!fx_subs.empty()) {
        fxvars fx=fx_subs.top();
        vecteur csol(fx.n_oldvars,undef),oldvars=vars,fxvals;
        vars=vecteur(fx.n_oldvars,undef);
        for (i=fx.indices.size();i-->0;) vars[fx.indices[i]]=fx.names[i];
        i=0;
        for (const_iterateur it=oldvars.begin();it!=oldvars.end();++it) {
            while (!is_undef(vars[i])) ++i;
            vars[i]=*it;
        }
        postprocess_point(x0,oldvars,fx,csol);
        fx_subs.pop();
    }
}
/* NLP_REGION CLASS IMPLEMENTATION */
nlp_problem::region::region(nlp_problem &p,double tol,int itlim) : prob(p), ctx(p.ctx) {
    ftol=tol;
    maxiter=itlim;
    ipt=NULL;
    n=prob.var_count();
    x_lb=prob.lb;
    x_ub=prob.ub;
    c=mergevecteur(prob.eq_cons,prob.ineq_cons);
    m=c.size();
    dc=mergevecteur(prob.eq_jacobian,prob.ineq_jacobian);
    cx.resize(m);
    dir.resize(m);
    c.resize(m); cx.resize(m); dc.resize(m);
    kl.resize(n); ku.resize(n);
    sl.resize(n); su.resize(n);
    int neq=prob.eq_cons.size();
    for (int i=0;i<m;++i) {
        cx[i]=_lname(c[i],ctx);
        dir[i]=i<neq?0:1;
    }
}
/* Constraint consensus algorithm: moves x0 to the feasible region, returns true iff successful.
 * SOURCE: Ibrahim and Chinneck 2005 (in "Feasibility and Infeasibility in Optimization", p. 70) */
bool nlp_problem::region::mfeas(vecteur &x0,bool store) {
    int NINF,iter=1;
    vecteur p(x0),t(n);
    gen evb;
    for (;iter<=maxiter;++iter) {
        NINF=0;
        std::fill(kl.begin(),kl.end(),0);
        std::fill(sl.begin(),sl.end(),0);
        std::fill(ku.begin(),ku.end(),0);
        std::fill(su.begin(),su.end(),0);
        for (int i=0;i<m;++i) { // for each constraint c
            vecteur nabla_c=subst(*dc[i]._VECTptr,prob.vars,p,false,ctx);
            gen viol,ec=subst(c[i],prob.vars,p,false,ctx);
            viol=dir[i]==0?_abs(ec,ctx):(dir[i]<0?max(0,ec,ctx):max(0,-ec,ctx));
            if (!is_zero(viol,ctx)) {
                vecteur fv=multvecteur((dir[i]==0?_sign(-ec,ctx):dir[i])*viol/pow(l2norm(nabla_c,ctx),2),nabla_c);
                gen d=l2norm(fv,ctx);
                if (is_strictly_greater(d,ftol,ctx)) {
                    NINF++;
                    for (int j=0;j<n;++j) { // for each variable in c
                        if (!contains(*cx[i]._VECTptr,prob.vars[j]))
                            continue;
                        if (is_strictly_positive(fv[j],ctx)) {
                            ku[j]++;
                            if (is_strictly_greater(fv[j],su[j],ctx))
                                su[j]=fv[j];
                        } else if (!is_zero(fv[j],ctx)) {
                            kl[j]++;
                            if (is_strictly_greater(sl[j],fv[j],ctx))
                                sl[j]=fv[j];
                        }
                    }
                }
            }
        }
        if (NINF==0 || (ipt!=NULL && ipt->filter_accepts(ipt->vars2x(p)))) {
            if (store)
                x0=p;
            return true;
        }
        std::fill(t.begin(),t.end(),0);
        for (int j=0;j<n;++j) {
            if (kl[j]==ku[j])
                t[j]=(sl[j]+su[j])/2;
            else if (ku[j]>kl[j])
                t[j]=su[j];
            else t[j]=sl[j];
        }
        if (is_zero(l2norm(t,ctx),ctx))
            return false; // movement tolerance exceeded
        p=addvecteur(p,t);
        for (int j=0;j<n;++j) {
            if (is_strictly_greater(p[j],x_ub[j],ctx))
                p[j]=x_ub[j];
            else if (is_strictly_greater(x_lb[j],p[j],ctx))
                p[j]=x_lb[j];
        }
    }
    return false; // max number of iterations exceeded
}
/* Randomized standard heuristic */
void nlp_problem::region::rs_heur(vecteur &x0) const {
    x0.resize(n);
    gen D=rand_uniform(0.5,1.0,ctx);
    for (int j=0;j<n;++j) {
        if (!is_inf(x_lb[j]) && !is_inf(x_ub[j])) { // doubly bounded
            gen r=(x_ub[j]-x_lb[j])/2;
            if (is_greater(D,r,ctx))
                D=rand_uniform(r/10,(9*r)/10,ctx);
            x0[j]=(x_lb[j]+x_ub[j])/2+D;
        } else if (!is_inf(x_lb[j])) // single lower bound
            x0[j]=x_lb[j]+D;
        else if (!is_inf(x_ub[j])) // single upper bound
            x0[j]=x_ub[j]-D;
        else x0[j]=D; // unbounded in both directions
    }
}
gen nlp_problem::region::bisection(const vecteur &x0,int i,const gen &a,const gen &b,int d) {
    assert(is_strictly_greater(b,a,ctx) && d!=0 && i>=0 && i<int(x0.size()));
    if (is_greater(ftol,b-a,ctx))
        return d>0?b:a;
    vecteur p(x0);
    gen q=(a+b)/2;
    p[i]=q;
    gen old_b=d>0?x_lb[i]:x_ub[i];
    if (d>0) x_lb[i]=q; else x_ub[i]=q;
    bool succ=mfeas(p,false);
    if (d>0) x_lb[i]=old_b; else x_ub[i]=old_b;
    if (succ)
        return bisection(x0,i,d>0?q:a,d>0?b:q,d);
    return bisection(x0,i,d>0?a:q,d>0?q:b,d);
}
bool nlp_problem::region::strides(const vecteur &x0_orig,int i,int d,int lim,gen &a,gen &b) {
    assert(d!=0 && lim>0 && !x0_orig.empty() && i>=0 && i<int(x0_orig.size()));
    if (d>0) b=plus_inf; else a=minus_inf;
    vecteur x0(x0_orig);
    gen old_x=x0[i],old_b=d>0?x_lb[i]:x_ub[i];
    double step=0.618033;
    bool ret=false;
    for (int j=0;j<lim;++j) {
        gen dx(step*(1<<j));
        if (d>0) x_lb[i]=old_x+dx; else x_ub[i]=old_x-dx;
        x0[i]=d>0?x_lb[i]:x_ub[i];
        if (!mfeas(x0,false)) {
            a=d>0?old_x:x_ub[i];
            b=d>0?x_lb[i]:old_x;
            ret=true;
            break;
        }
        old_x=d>0?x_lb[i]:x_ub[i];
    }
    if (d>0) x_lb[i]=old_b; else x_ub[i]=old_b;
    return ret;
}
/* Find bound of the variable x[i] in direction d (d<0 is lower bound, d>0 is upper bound) */
bool nlp_problem::region::find_bound(const vecteur &x0,int i,int d,gen &bnd) {
    if (prob._linconstr) {
        vecteur sol;
        gen optval;
        int res=prob.minimize_linear((d<0?1:-1)*prob.vars[i],prob.constraints(true),sol,optval);
        if (res==_NLP_OPTIMAL) {
            bnd=sol[i];
            return true;
        }
        return false;
    }
    gen a,b;
    if (!strides(x0,i,d,25,a,b))
        return false;
    bnd=bisection(x0,i,a,b,d);
    return true;
}
/* Compute bounding rectangle of the feasible area */
matrice nlp_problem::region::bounding_rect() {
    matrice bd=mtran(makevecteur(x_lb,x_ub));
    if (has_inf_or_undef(bd)) {
        prob.msg("Estimating variable bounds...");
        vecteur x0;
        rs_heur(x0);
        gen bnd;
        int i=0;
        for (;i<n;++i) {
            gen &lb=bd[i]._VECTptr->front(),&ub=bd[i]._VECTptr->back();
            if (is_inf(lb) && find_bound(x0,i,-1,bnd))
                lb=bnd-ftol;
            if (is_inf(ub) && find_bound(x0,i,1,bnd))
                ub=bnd+ftol;
            if (is_strictly_positive(lb-ub,ctx)) {
                prob._infeas=true;
                prob.warn("problem appears to be infeasible, try again with bounded variables or initial rectangle");
                break;
            }
        }
        if (has_inf_or_undef(bd))
            prob.msg("Feasible region is unbounded");
        else
            prob.msg("Feasible region is bounded");
    }
    return bd;
}
/* END NLP_REGION CLASS */

/* COBYLA optimizer. */
bool nlp_problem::cobyla(const vecteur &cons,const meth_parm &parm,optima_t &res) {
    vecteur fcv=makevecteur(obj,cons,vars),x0;
    cobyla_gc gc={fcv,ctx};
    int maxiter=parm.max_iter,cres,iter,tries,max_tries=5;
    while (maxiter>0 && restart(x0)) {
        optimum_t opt;
        try { // run COBYLA
            tries=0;
            while (maxiter>0 && tries++<max_tries) {
                iter=maxiter;
                cres=giac_cobyla(&gc,x0,iter,parm.eps,0);
                maxiter-=iter;
                if (cres!=COBYLA_ROUNDING)
                    break;
            }
            if (cres==COBYLA_MAXFUN)
                maxiter=0;
            if (cres<0) // error
                opt.res=_NLP_ERROR;
            else {
                opt.x=x0;
                save_point(x0);
                opt.res=cres>0?_NLP_FAILED:_NLP_OPTIMAL;
            }
        } catch (const std::runtime_error &err) { // an error raised in COBYLA or user interruption
            msg(err.what());
            return false;
        }
        res.push_back(opt);
    }
    return true;
}
/* Evaluate objective Hessian at point x0 (if full=false, evaluate only the lower triangle). */
bool nlp_problem::subst_hessian(const vecteur &x0,matrice &H,bool full) const {
    if (!_have_hessian || mrows(H)!=var_count() || mcols(H)!=var_count())
        return false;
    const_iterateur it=obj_hessian.begin(),itend=obj_hessian.end(),jt,jtend;
    int i,j;
    for (i=0;it!=itend;++it,++i) {
        jt=it->_VECTptr->begin();
        jtend=it->_VECTptr->end();
        for (j=0;jt!=jtend;++jt,++j) {
            gen hij=subst(*jt,vars,x0,false,ctx);
            if (!is_real_number(hij,ctx))
                return false;
            hij=to_real_number(hij,ctx);
            H[i]._VECTptr->at(j)=hij;
            if (full && i!=j)
                H[j]._VECTptr->at(i)=hij;
        }
    }
    return true;
}
/* Minimize O subject to C and problem variables by using LPSOLVE. */
int nlp_problem::minimize_linear(const gen &o,const vecteur &c,vecteur &sol,gen &optval) const {
    gen res;
    vecteur iv;
    int i=0,n=var_count();
    for (;i<n;++i) if (is_intvar(i)) iv.push_back(vars[i]);
    try {
        log_output_redirect lor(ctx);
        res=_lpsolve(makesequence(o,c,symb_equal(change_subtype(_LP_INTEGERVARIABLES,_INT_MAPLECONVERSION),iv)),ctx);
    } catch (const std::runtime_error &e) {
        return _NLP_ERROR;
    }
    if (res.type!=_VECT)
        return _NLP_ERROR;
    if (res._VECTptr->empty())
        return _NLP_INFEAS;
    if (res._VECTptr->size()!=2 || res._VECTptr->back().type!=_VECT)
        return _NLP_ERROR;
    optval=res._VECTptr->front();
    if (is_inf(optval))
        return _NLP_UNBOUNDED;
    const vecteur &tmpsol=*res._VECTptr->back()._VECTptr;
    sol.resize(vars.size(),undef);
    int j;
    for (const_iterateur it=tmpsol.begin();it!=tmpsol.end();++it) {
        if (!is_equal(*it) || (j=indexof(it->_SYMBptr->feuille._VECTptr->front(),vars))<0)
            continue;
        sol[j]=it->_SYMBptr->feuille._VECTptr->back();
    }
    return has_inf_or_undef(sol)?_NLP_ERROR:_NLP_OPTIMAL;
}
/* Return true if G is a vector of intervals. */
bool is_vect_of_intervals(const gen &g,vecteur &lb,vecteur &ub,GIAC_CONTEXT) {
    if (g.type!=_VECT || g._VECTptr->empty())
        return false;
    lb.resize(g._VECTptr->size());
    ub.resize(g._VECTptr->size());
    for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
        int i=it-g._VECTptr->begin();
        if (!it->is_symb_of_sommet(at_interval))
            return false;
        const vecteur &lr=*it->_SYMBptr->feuille._VECTptr;
        if (!is_real_number(lr.front(),contextptr) || !is_real_number(lr.back(),contextptr))
            return false;
        lb[i]=to_real_number(lr.front(),contextptr);
        ub[i]=to_real_number(lr.back(),contextptr);
    }
    return true;
}
/* Compute initial size from parm.init_size.
 * If parm.init_size is zero, return uniform random number between LO*D and HI*D,
 * where D is the range of the smallest dimension. Otherwise, return parm.init_size.
 * If there is not a variable with lower and upper bound, and parm.init_size is zero,
 * then return 1 (each dimension has an inifinite range).
 */
double nlp_problem::initial_size(const meth_parm &parm,double lo,double hi) const {
    if (parm.size>0)
        return parm.size;
    int n=var_count();
    double ret;
    ret=rand_uniform(lo,hi,ctx).to_double(ctx);
    bool has_bounded_dim=false;
    gen min_dim_range=plus_inf;
    for (int i=0;i<n;++i) {
        gen rng=ub[i]-lb[i];
        if (is_inf(rng))
            continue;
        has_bounded_dim=true;
        if (is_strictly_greater(min_dim_range,rng,ctx))
            min_dim_range=rng;
    }
    return has_bounded_dim?ret*to_real_number(min_dim_range,ctx).to_double(ctx):1.0;
}
void nlp_problem::randomize_ith_component(vecteur &x,int i,bool use_varbnds) const {
    const gen &lo=use_varbnds?lb[i]:bd_rect[i][0];
    gen hi=(use_varbnds?ub[i]:bd_rect[i][1])+(is_intvar(i)?1:0);
    if (is_inf(lo) && is_inf(hi))
        x[i]=rand_normal(1.0,false,ctx);
    else if (is_inf(hi))
        x[i]=lo+rand_normal(1.0,true,ctx);
    else if (is_inf(lo))
        x[i]=hi-rand_normal(1.0,true,ctx);
    else x[i]=rand_uniform(lo,hi,ctx);
}
/* Generate a random initial point X which conforms to the bounds of variables.
 * This procedure avoids previously stored points. If convexity clusters are available,
 * it attempts to generate points near their centers.
 * NSAMP - maximum number of samples that will be generated in case there are stored points.
 * SIGMA-SCALE - standard-deviation scaling along the dimensions if there are stored points.
 * Return true on success, else return false. */
bool nlp_problem::make_random_initial_point(vecteur &x,gen *fval,int nsamp,double sigma_scale) {
    int n=vars.size(),iter=0,ccsize=cclust_centers.size(),k;
    vecteur samples,x0(n);
    while (nsamp>(int)samples.size() && ++iter<nsamp*10) {
        if (ccsize==0) {
            if (saved_points.empty()) {
                /* No previously saved points, generate an uniformly random point in bd_rect,
                * for unbounded dimensions resp. those bounded only from one side,
                * use normal distribution with mean zero resp. the bound. */
                for (int i=0;i<n;++i)
                    randomize_ith_component(x0,i,false);
            } else {
                /* stick to stored points */
                if (initp_sigma.empty()) { // compute sigma for each dimension
                    initp_sigma.resize(n);
                    for (k=0;k<n;++k) {
                        gen &sk=initp_sigma[k];
                        if (is_inf(sk=std::sqrt(sigma_scale)*(bd_rect[k][1]-bd_rect[k][0])))
                            sk=1;
                    }
                }
                k=_rand(saved_points.size(),ctx).val;
                rand_multinormal(*saved_points[k]._VECTptr,initp_sigma,x0,ctx);
            }
        } else { // have clusters, try generating a point within a randomly chosen cluster
            if (cclust_freqs.empty())
                cclust_freqs.resize(cclust_centers.size(),0);
            vecteur z=*_find(makesequence(0,cclust_freqs),ctx)._VECTptr;
            /* Choose uniformly between cluster centers that haven't been selected before,
            * if there are not such clusters select k-th cluster with probability
            * proportional to 1/freq[k], where freq[k] is the choice frequency of k-th cluster. */
            if (z.empty())
                k=rand_categorical(*_inv(z,ctx)._VECTptr,ctx);
            else k=z[_rand(z.size(),ctx).val].val;
            // random normal variable with cluster center as mean and cluster size as sigma
            rand_multinormal(*cclust_centers[k]._VECTptr,cclust_sizes[k],x0,ctx);
            // assure the point resides within the variable bounds
        }
        if (clamp(x0) && ccsize==0 && !saved_points.empty())
            continue;
        if (fval==NULL || compute_obj_val(x0,*fval)) {
            samples.push_back(makevecteur(x0,fval==NULL?undef:*fval));
            if (ccsize>0) {
                cclust_freqs[k]+=1;
                break;
            }
            if (saved_points.empty())
                break;
        }
    }
    if (samples.empty())
        return false;
    if (samples.size()==1) {
        x=*samples.front()._VECTptr->front()._VECTptr;
        if (fval!=NULL)
            *fval=samples.front()._VECTptr->back();
    } else { // multiple samples
        assert(!saved_points.empty());
        gen minp=plus_inf,pval,e;
        const_iterateur it=samples.begin(),itend=samples.end(),jt,jtend;
        for (;it!=itend;++it) {
            pval=0;
            const vecteur &pt=*it->_VECTptr->front()._VECTptr;
            for (jt=saved_points.begin(),jtend=saved_points.end();jt!=jtend;++jt) {
                vecteur d=subvecteur(pt,*jt->_VECTptr);
                e=0;
                for (k=0;k<n;++k) e+=pow(d[k]/initp_sigma[k],2);
                pval+=exp(-e/2,ctx);
            }
            if (is_strictly_greater(minp,pval,ctx)) {
                x=pt;
                minp=pval;
                if (fval!=NULL)
                    *fval=it->_VECTptr->back();
            }
        }
    }
    return true;
}
/* Make next initial point.
 * If there are initial points provided by the user, they will be iterated over.
 * Otherwise, make_random_initial_point is used to generate one.
 * If this method returns false, the respective optimization process should terminate. */
bool nlp_problem::restart(vecteur &x,gen *fval) {
    int max_tries=10,tries=0;
    bool ret=true;
    if (initp_iter!=initp.end()) {
        // fetch the next initial point among those provided in INITP
        x=*(initp_iter++)->_VECTptr;
        if (fval!=NULL)
            ret=compute_obj_val(x,*fval);
    } else {
        if (!initp.empty())
            return false; // done iterating over the user-provided initial points
        // attempt to generate a random initial point
        while (tries++<max_tries) {
            ret=make_random_initial_point(x,fval);
            if (ret) {
                save_point(x);
                break;
            }
        }
    }
    if (ret) { _ss << "Starting from point " << x; debug(); }
    return ret;
}
/* Compute objective function value at X and store it to VAL.
 * Penalty function is integrated according to the penalty scheme and weights.
 * Return true if the computed value is a real number. */
bool nlp_problem::compute_obj_val(const vecteur &x_orig,gen &val,bool print_err) {
    vecteur x(x_orig);
    int i=0,n=vars.size();
    if (n!=(int)x.size())
        return false;
    for (;i<n;++i) if (is_intvar(i)) x[i]=INT(x[i]);
    try {
        val=_evalf(subst(obj,vars,x,false,ctx),ctx);
        const_iterateur it,pit,itend;
        gen s,b;
        switch (_penalty_scheme) {
        case 0: // no penalty function
            break;
        case 1: // Lagrangian
            for (it=eq_cons.begin(),itend=eq_cons.end(),pit=penalty_eq.begin();it!=itend;++it,++pit)
                if (!is_zero(*pit,ctx)) val+=*pit*subst(_abs(*it,ctx),vars,x,false,ctx);
            for (it=ineq_cons.begin(),itend=ineq_cons.end(),pit=penalty_ineq.begin();it!=itend;++it,++pit)
                if (!is_zero(*pit,ctx)) val+=*pit*max(0,-subst(*it,vars,x,false,ctx));
            break;
        case 2: // multiplicative penalty function
            val+=obj_shift;
            for (it=eq_cons.begin(),itend=eq_cons.end(),pit=penalty_eq.begin();it!=itend;++it,++pit) {
                reim(*pit,s,b,ctx);
                val=val*pow(1+s*subst(_abs(*it,ctx),vars,x,false,ctx),b,ctx);
            }
            for (it=ineq_cons.begin(),itend=ineq_cons.end(),pit=penalty_ineq.begin();it!=itend;++it,++pit) {
                reim(*pit,s,b,ctx);
                val=val*pow(1+s*max(0,-subst(*it,vars,x,false,ctx),ctx),b,ctx);
            }
            break;
        default:
            assert(false); // not implemented
        }
        if (is_real_number(val,ctx)) {
            val=to_real_number(val,ctx);
            _fevalc++;
            return true;
        }
    } catch (const std::runtime_error &e) {
        if (print_err) { _ss << e.what(); err(); }
    }
    val=undef;
    if (print_err)
        err("failed to compute objective value");
    return false;
}
/* Clamp X0 to fit the bounds of variables. */
bool nlp_problem::clamp(vecteur &x0,bool rnd) const {
    int n=vars.size();
    gen nb,xval;
    bool clamped=false,lviol,uviol;
    for (int i=0;i<n;++i) {
        xval=is_intvar(i)?INT(x0[i]):x0[i];
        lviol=is_strictly_greater(lb[i],xval,ctx);
        uviol=is_strictly_greater(xval,ub[i],ctx);
        if (!clamped)
            clamped=lviol||uviol;
        if (rnd && (lviol||uviol))
            randomize_ith_component(x0,i,true);
        if (!rnd) {
            if (lviol) x0[i]=lb[i]; else if (uviol) x0[i]=ub[i];
        }
    }
    return clamped;
}
/* Return true if at least FR*100% points in X are within the specified
 * relative tolerance TOL from X[BEST]. */
bool is_within_tolerance(const vecteur &x,int best,double tol,double fr,GIAC_CONTEXT) {
    const vecteur &v=*x[best]._VECTptr;
    gen vnorm=l2norm(v,contextptr);
    int cnt=0;
    for (const_iterateur it=x.begin();it!=x.end();++it) {
        if (is_strictly_greater(tol,l2norm(subvecteur(*it->_VECTptr,v),contextptr)/(1+vnorm),contextptr))
            cnt++;
    }
    return cnt>=fr*x.size();
}
/* Generate initial simplex at X0 with the given SIZE (Haftka and Guerdal 1993).
 * Note that simplex will be clamped to fit variable bounds.
 * Return TRUE on success, else return false (when objective is undefined at some vertex). */
bool nlp_problem::initialize_simplex(const vecteur &x0,double size) {
    nm_simplex[0]=x0;
    int n=var_count();
    double q=size*(std::sqrt(n+1)-1)/(n*M_SQRT2),p=size*(std::sqrt(n+1)+n-1)/(n*M_SQRT2);
    for (int i=1;i<=n;++i) {
        if (nm_simplex[i].type==_VECT)
            std::fill(nm_simplex[i]._VECTptr->begin(),nm_simplex[i]._VECTptr->end(),q);
        else nm_simplex[i]=vecteur(n,q);
        vecteur &Si=*nm_simplex[i]._VECTptr;
        Si[i-1]=p;
        addvecteur(Si,x0,Si);
        clamp(Si); if (!compute_obj_val(Si,nm_yarr[i])) return false;
    }
    return true;
}
/* Initialize all penalties to P. */
void nlp_problem::initialize_penalties(const gen &p) {
    if (penalty_eq.size()!=eq_cons.size())
        penalty_eq.resize(eq_cons.size());
    std::fill(penalty_eq.begin(),penalty_eq.end(),p);
    if (penalty_ineq.size()!=ineq_cons.size())
        penalty_ineq.resize(ineq_cons.size());
    std::fill(penalty_ineq.begin(),penalty_ineq.end(),p);
}
/* Nelder-Mead algorithm iteration, which saves the current best point.
 * if NO_GSL is true, then GSL implementation will not be used (e.g. when some variables are bounded).
 * Return values: negative (error), 0 (continue), 1 (converged/small), 2 (flat/stalled), 3 (degenerated). */
int nlp_problem::nm_iteration(vecteur &best_point,gen &fval,double &ssize,const meth_parm &parm,bool no_gsl) {
    _iter_count++;
#ifdef HAVE_LIBGSL
    /* GSL NM implementation 'nmsimplex2' (linear complexity) */
    if (!no_gsl) {
        if (nm_s==NULL)
            return -1; // error
        int status=gsl_multimin_fminimizer_iterate(nm_s); // make a move
        if (status)
            return 2; // cannot continue
        best_point=gsl_vector2vecteur(gsl_multimin_fminimizer_x(nm_s));
        fval=gsl_multimin_fminimizer_minimum(nm_s);
        double size=gsl_multimin_fminimizer_size(nm_s); // simplex size
        if (gsl_multimin_test_size(size,parm.eps)==GSL_SUCCESS)
            return 1; // converged
        return 0; // continue iterating
    }
#endif
    /* Giac NM implementation, detects small/flat/degenerated simplices */
    matrice E; // simplex edges
    vector<int> p;
    int n=var_count(),i,opt_ind;
    sortperm(nm_yarr,p,ctx); // sort from lowest to highest
    apply_permutation(nm_simplex,p,true);
    apply_permutation(nm_yarr,p,true);
    _base_vertex_index=p[_base_vertex_index]; // update index of the base vertex
    vecteur xl=*nm_simplex.front()._VECTptr,xh=*nm_simplex.back()._VECTptr,xs=*nm_simplex[n-1]._VECTptr;
    gen yl=nm_yarr.front(),yh=nm_yarr.back(),ys=nm_yarr[n-1],yr,ye,yc,dc,rng;
    vecteur xm=*_mean(vecteur(nm_simplex.begin(),nm_simplex.begin()+n),ctx)._VECTptr;
    vecteur xr=addvecteur(xm,multvecteur(nm_alpha,subvecteur(xm,xh))); // reflection point
    clamp(xr); if (!compute_obj_val(xr,yr)) return -1;
    ssize=0;
    if (is_strictly_greater(yl,yr,ctx)) {
        vecteur xe=addvecteur(xm,multvecteur(nm_beta,subvecteur(xr,xm))); // expansion point
        clamp(xe); if (!compute_obj_val(xe,ye)) return -1;
        if (is_strictly_greater(yr,ye,ctx)) { // expansion
            nm_simplex.back()=xe; nm_yarr.back()=ye;
        } else { // reflection
            nm_simplex.back()=xr; nm_yarr.back()=yr;
        }
    } else if (is_strictly_greater(yr,ys,ctx)) {
        if (is_greater(yh,yr,ctx)) {
            xh=xr; yh=yr; nm_simplex.back()=xr; nm_yarr.back()=yr;
        }
        vecteur xc=addvecteur(xm,multvecteur(nm_gamma,subvecteur(xh,xm))); // contraction point
        if (!compute_obj_val(xc,yc)) return -1;
        if (is_strictly_greater(yc,yh,ctx)) { // shrinking
            gen maxd=minus_inf,d;
            for (i=1;i<n+1;++i) {
                vecteur &Si=*nm_simplex[i]._VECTptr;
                Si=addvecteur(multvecteur(1-nm_delta,Si),multvecteur(nm_delta,xl));
                if (!compute_obj_val(Si,nm_yarr[i])) return -1;
            }
        } else { // contraction
            nm_simplex.back()=xc; nm_yarr.back()=yc;
        }
    } else { // reflection
        nm_simplex.back()=xr; nm_yarr.back()=yr;
    }
    // store best point and the corresponding objective value
    opt_ind=argmin(nm_yarr,ctx);
    best_point=*nm_simplex[opt_ind]._VECTptr;
    fval=nm_yarr[opt_ind];
    // get simplex edges
    E.clear();
    const vecteur &bv=*nm_simplex[_base_vertex_index]._VECTptr;
    for (i=0;i<=n;++i) {
        if (i!=_base_vertex_index)
            E.push_back(subvecteur(*nm_simplex[i]._VECTptr,bv));
    }
    assert(is_squarematrix(E));
    for (const_iterateur it=E.begin();it!=E.end();++it) {
        const_iterateur jt=it->_VECTptr->begin(),jtend=it->_VECTptr->end();
        dc=0;
        for (i=0;jt!=jtend;++jt,++i) {
            rng=(!bd_rect.empty() && !has_inf_or_undef(bd_rect[i]))?bd_rect[i][1]-bd_rect[i][0]:1;
            dc+=_abs(*jt/rng,ctx);
        }
        if (is_strictly_greater(dc,ssize,ctx))
            ssize=_evalf(dc,ctx).to_double(ctx);
    }
    if (is_strictly_greater(parm.eps,ssize,ctx))
        return 1;
    if (is_zero(_max(nm_yarr,ctx)-_min(nm_yarr,ctx),ctx))
        return 2;
    vecteur Enorm=*_apply(makesequence(at_l2norm,E),ctx)._VECTptr;
    if (is_strictly_greater(parm.tol,_min(Enorm,ctx)/_max(Enorm,ctx),ctx) ||
            is_strictly_greater(parm.tol,_abs(_det(E,ctx),ctx)/_product(Enorm,ctx),ctx))
        return 3;
#if 0
    if (is_strictly_greater(parm.eps,_stddev(nm_yarr,ctx),ctx)) // original Nelder-Mead termination criterion
        return 1;
#endif
    return 0; // continue iterating
}
/* Nelder-Mead minimizer. */
bool nlp_problem::nelder_mead(const meth_parm &parm,optima_t &res) {
    int maxiter=parm.max_iter,n=var_count(),status,degeneracy_count=0;
    if (n<=0 || maxiter<=0)
        return false;
    bool no_gsl=has_bounded_variables(),has_ips=!initp.empty(),fixed_penalty=false;
    double nm_initsize,ssize=-1;
    vecteur ip,bp;
    gen fval;
#ifdef HAVE_LIBGSL
    gsl_vector *nm_ss=NULL,*nm_x0=NULL;
    gsl_multimin_function nm_func;
    if (!no_gsl) {
        msg("Using GSL minimizer...");
        nm_s=gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2,n);
        nm_func.n=n;
        nm_func.f=&nlp_problem::gsl_my_f;
        nm_func.params=static_cast<void*>(this);
        nm_ss=gsl_vector_alloc(n);
        gsl_vector_set_all(nm_ss,initial_size(parm,0.02,0.10));
        initialize_penalties(1.0/parm.eps);
        fixed_penalty=true;
    } else {
#endif
        nm_alpha=parm.reflect_ratio;
        nm_beta=parm.expand_ratio;
        nm_gamma=parm.contract_ratio;
        nm_delta=parm.shrink_ratio;
        if (nm_alpha<=0 || nm_beta<=1 || nm_gamma<=0 || nm_gamma>=1 || nm_delta<=0 || nm_delta>=1) {
            err("invalid parameter(s) for Nelder & Mead algorithm");
            return false;
        }
        nm_initsize=initial_size(parm,0.02,0.1);
        nm_simplex.resize(n+1);
        nm_yarr.resize(n+1);
        initialize_penalties(is_real_number(parm.penalty,ctx)&&is_positive(parm.penalty,ctx)?to_real_number(parm.penalty,ctx):0);
#ifdef HAVE_LIBGSL
    }
#endif
    bool start_from_last_point=false;
    gen best_L,L;
    const_iterateur cit,citend;
    iterateur pit;
    gen pstep=parm.step<=0?1e-2:parm.step;
    bool has_constraints=eq_count()+ineq_count()>0;
    while (_iter_count<maxiter) {
        optimum_t opt;
        /* 1) provide initial point */
        if (start_from_last_point) {
            ip=bp;
            start_from_last_point=false;
        } else { // start from a new initial point
            if (!restart(ip))
                break;
        }
        /* 2) initialize simplex */
#ifdef HAVE_LIBGSL
        if (!no_gsl) {
            nm_x0=vecteur2gsl_vector(ip,ctx);
            gsl_multimin_fminimizer_set(nm_s,&nm_func,nm_x0,nm_ss);
        } else {
#endif
            if (ssize<0) ssize=nm_initsize;
            else { // simplex has to be inflated, need to move away from the bounds
                for (int i=0;i<n;++i) {
                    if (is_greater(2*ssize,ub[i]-lb[i],ctx))
                        ip[i]=(lb[i]+ub[i])/2.0;
                    else if (is_greater(ssize,ip[i]-lb[i],ctx))
                        ip[i]=lb[i]+ssize;
                    else if (is_greater(ssize,ub[i]-ip[i],ctx))
                        ip[i]=ub[i]-ssize;
                }
                _ss << "Inflating simplex at point " << ip << "..."; debug();
            }
            if (!initialize_simplex(ip,ssize)) {
                err("simplex initialization failed, aborting...");
                opt.res=_NLP_ERROR;
                break;
            } else _base_vertex_index=0;
#ifdef HAVE_LIBGSL
        }
#endif
        /* 3) iterate */
        status=2;
        while (_iter_count<maxiter && (status=nm_iteration(bp,fval,ssize,parm,no_gsl))==0) {
            if (fixed_penalty || !has_constraints) continue;
            /* handle constraints using adaptive linear penalty (Luersen, Riche & Guyon, 2003) */
            // update penalty parameters
            for (cit=eq_cons.begin(),pit=penalty_eq.begin(),citend=eq_cons.end();cit!=citend;++cit,++pit)
                *pit+=pstep*_abs(subst(*cit,vars,bp,false,ctx),ctx);
            for (cit=ineq_cons.begin(),pit=penalty_ineq.begin(),citend=ineq_cons.end();cit!=citend;++cit,++pit)
                *pit+=pstep*max(0,-subst(*cit,vars,bp,false,ctx),ctx);
            // recompute objective values at the simplex vertices
            for (cit=nm_simplex.begin(),citend=nm_simplex.end(),pit=nm_yarr.begin();cit!=citend;++cit,++pit)
                compute_obj_val(*cit->_VECTptr,*pit);
        }
        switch (status) {
        case 1: // converged
            debug("Converged");
            opt.res=_NLP_OPTIMAL;
            break;
        case 2: // flat simplex/stalled
            debug("Stalled or simplex is flat");
            opt.res=_NLP_FAILED;
            break;
        case 3: // degenerated simplex
            _ss << "Simplex is degenerated (size: " << ssize << ")"; debug();
            if (++degeneracy_count==3) {
                opt.res=_NLP_FAILED;
                pstep=10*pstep; // increase step-size for constraints
                _ss << "Increasing step-size to " << pstep << "..."; debug();
                break;
            }
            start_from_last_point=true;
            continue; // restart from the best point
        default: // error
            debug("Aborted");
            opt.res=_NLP_ERROR;
            break;
        }
        degeneracy_count=0;
        ssize=-1; // restart with simplex of initial size
        if (status>0) {
            _ss << "Best point: " << bp; debug();
            opt.x=bp;
            save_point(bp);
        }
        /* 4) store the best point as a solution candidate */
        res.push_back(opt);
    }
#ifdef HAVE_LIBGSL
    if (!no_gsl) {
        gsl_multimin_fminimizer_free(nm_s);
        gsl_vector_free(nm_ss);
        gsl_vector_free(nm_x0);
    }
#endif
    _ss << "Penalties: " << mergevecteur(penalty_eq,penalty_ineq); debug();
    return true;
}
/* Initialize the population for differential evolution. */
bool nlp_problem::de_initialize(iterateur &it,const iterateur &itend,iterateur &kt,iterateur &ft,int &best_k) {
    bool ok;
    int i=0,ips=initp.size();
    for (;it!=itend;++it,++kt,++ft,++i) {
        *it=i;
        if (i<ips) {
            *kt=initp[i];
            ok=compute_obj_val(*kt->_VECTptr,*ft);
        } else *kt=vecteur(0);
        if ((i>=ips || !ok) && !make_random_initial_point(*kt->_VECTptr,&*ft))
            return false;
        if (is_strictly_greater(best_obj_val,*ft,ctx)) {
            best_obj_val=*ft;
            best_k=i;
        }
    }
    return true;
}
/* Differential evolution minimizer. */
bool nlp_problem::differential_evolution(const meth_parm &parm,optimum_t &opt) {
    int conv_iter=10,m=parm.search_points,maxiter=parm.max_iter+conv_iter;
    double crossp=parm.cross_prob,tol=parm.tol,eps=parm.eps,s=parm.scale;
    int n=var_count(),k,i,ia,ib,ic;
    if (m<=0)
        m=std::max((int)initp.size(),10*n);
    if (s<=0) { // compute an "optimal" scale (Zaharie, 2002)
        s=std::sqrt((2-crossp)/(2*n)+5e-2/crossp);
        _ss << "Optimal scale factor: " << s; debug();
    }
    if (m<=0 || maxiter<=0 || crossp<0 || crossp>1 || s<0 || s>2 || tol<=0 || eps<=0)
        return false;
    /* create initial population and constraint penalties */
    best_obj_val=plus_inf;
    int best_k;
    initialize_penalties(parm.penalty.type==_CPLX?parm.penalty:makecomplex(10,1));
    vecteur x(m),fval(m),ind(m);
    iterateur it=ind.begin(),itend=ind.end(),kt=x.begin(),ktend,ft=fval.begin();
    if (!de_initialize(it,itend,kt,ft,best_k)) {
        err("failed to initialize population");
        return false;
    }
    if (eq_count()>0 || ineq_count()>0) // compute the objective shift
        obj_shift=max(-_min(fval,ctx),0,ctx)+1000;
    /* iterate DE */
    const_iterateur jt;
    bool converged=false,changed;
    vecteur best_obj_vals,pts,y,abc;
    best_obj_vals.reserve(maxiter);
    y.reserve(n);
    abc.reserve(3);
    set<int> viol_eq,viol_ineq;
    set<int>::const_iterator sit,sitend;
    while (_iter_count++<maxiter) {
        if (_iter_data) pts.push_back(x);
        best_obj_vals.push_back(best_obj_val);
        /* stopping criteria */
        if (_iter_count>=conv_iter &&
                is_greater(eps,(best_obj_vals.back()-best_obj_vals[_iter_count-conv_iter])/(1+best_obj_vals.back()),ctx) &&
                is_within_tolerance(x,best_k,tol,0.67,ctx)) {
            debug("Converged");
            converged=true;
            if (!is_feasible(*x[best_k]._VECTptr,tol,false,&viol_eq,&viol_ineq)) { // restart with larger penalty coefficients
                _ss << int(viol_eq.size()+viol_ineq.size()) << " "
                    << gettext("constraints violated, restarting with larger penalties..."); msg();
                for (sit=viol_eq.begin(),sitend=viol_eq.end();sit!=sitend;++sit) {
                    gen &p=penalty_eq[*sit];
                    p=makecomplex(10*re(p,ctx),im(p,ctx)+1);
                }
                for (sit=viol_ineq.begin(),sitend=viol_ineq.end();sit!=sitend;++sit) {
                    gen &p=penalty_ineq[*sit];
                    p=makecomplex(10*re(p,ctx),im(p,ctx)+1);
                }
                it=ind.begin(); itend=ind.end(); kt=x.begin(); ft=fval.begin();
                if (de_initialize(it,itend,kt,ft,best_k))
                    continue;
                err("failed to re-initialize the population");
                opt.res=_NLP_INFEAS;
            }
            break;
        }
        /* change the population */
        for (kt=x.begin(),ktend=x.end(),ft=fval.begin(),k=0;kt!=ktend;++kt,++ft,++k) {
            vecteur &xk=*kt->_VECTptr;
            /* pick three random agents different than xk */
            do abc=*_rand(makesequence(3,ind),ctx)._VECTptr; while (contains(abc,k));
            ia=abc[0].val; ib=abc[1].val; ic=abc[2].val;
            /* mutation */
            y=addvecteur(*x[ia]._VECTptr,multvecteur(s,subvecteur(*x[ib]._VECTptr,*x[ic]._VECTptr)));
            changed=false;
            for (it=y.begin(),itend=y.end(),jt=xk.begin();it!=itend;++it,++jt) {
                if (is_greater(crossp,rand_uniform(0,1,ctx),ctx)) {
                    changed=true;
                    *it=*jt;
                }
            }
            if (!changed) { // force a change anyway
                i=_rand(n,ctx).val;
                gen ival=y[i];
                y=xk;
                y[i]=ival;
            }
            clamp(y,true);
            /* fitness */
            gen val;
            if (compute_obj_val(y,val)) {
                if (is_greater(*ft,val,ctx)) {
                    xk=y;
                    *ft=val;
                    if (is_strictly_greater(best_obj_val,val,ctx)) {
                        best_obj_val=val;
                        best_k=k;
                    }
                }
            }
        }
    }
    obj_shift=0;
    if (_iter_data)
        iteration_data.push_back(pts);
    opt.x=*x[best_k]._VECTptr;
    for (i=0;i<n;++i) if (is_intvar(i)) opt.x[i]=INT(opt.x[i]);        
    //opt.f=best_obj_val;
    if (opt.res==_NLP_PENDING)
        opt.res=(converged?_NLP_OPTIMAL:_NLP_FAILED);
    return true;
}
#ifdef HAVE_LIBGSL
double nlp_problem::gsl_my_f(const gsl_vector *x,void *params) {
    nlp_problem *prob=static_cast<nlp_problem*>(params);
    gen fval;
    if (!prob->compute_obj_val(gsl_vector2vecteur(x),fval))
        return GSL_NAN;
    return fval.to_double(prob->ctx);
}
void nlp_problem::gsl_my_df(const gsl_vector *x,void *params,gsl_vector *df) {
    nlp_problem *prob=static_cast<nlp_problem*>(params);
    vecteur Dval=subst(prob->obj_gradient,prob->vars,gsl_vector2vecteur(x),false,prob->ctx);
    for (iterateur it=Dval.begin();it!=Dval.end();++it)
        *it=evalf_double(*it,1,prob->ctx);
    const iterateur it=Dval.begin(),itend=Dval.end();
    if (vecteur2gsl_vector(it,itend,df,prob->ctx)!=GSL_SUCCESS)
        prob->debug("Failed to compute gradient for GSL minimizer");
}
void nlp_problem::gsl_my_fdf(const gsl_vector *x,void *params,double *f,gsl_vector *df) {
    *f=gsl_my_f(x,params);
    gsl_my_df(x,params,df);
}
double nlp_problem::gsl_my_uf(double x,void *params) {
    nlp_problem *prob=static_cast<nlp_problem*>(params);
    gen fval;
    if (!prob->compute_obj_val(vecteur(1,x),fval))
        return GSL_NAN;
    return fval.to_double(prob->ctx);
}
/* This method uses the vector Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
 * This is a quasi-Newton method which builds up an approximation to the second
 * derivatives of the function f using the difference between successive gradient
 * vectors. By combining the first and second derivatives the algorithm is able
 * to take Newton-type steps towards the function minimum, assuming quadratic
 * behavior in that region. */
bool nlp_problem::gsl_bfgs(const meth_parm &parm,optima_t &res) {
    int N=vars.size(),maxiter=parm.max_iter,status,restarts,max_restarts=0;
    gsl_vector *xstart=gsl_vector_alloc(N);
    gsl_multimin_fdfminimizer *s=gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2,N);
    gsl_multimin_function_fdf my_func;
    my_func.n=N;
    my_func.f=&nlp_problem::gsl_my_f;
    my_func.df=&nlp_problem::gsl_my_df;
    my_func.fdf=&nlp_problem::gsl_my_fdf;
    my_func.params=static_cast<void*>(this);
    vecteur x0;
    while (_iter_count<maxiter && restart(x0)) {
        optimum_t opt;
        if (vecteur2gsl_vector(x0.begin(),x0.end(),xstart,ctx)!=GSL_SUCCESS)
            opt.res=_NLP_ERROR;
        else {
            opt.res=_NLP_FAILED;
            gsl_multimin_fdfminimizer_set(s,&my_func,xstart,parm.step<=0?1.0:parm.step,0.1);
            restarts=0;
            do {
                status=gsl_multimin_fdfminimizer_iterate(s);
                ++_iter_count;
                opt.x=gsl_vector2vecteur(gsl_multimin_fdfminimizer_x(s));
                opt.f=gsl_multimin_fdfminimizer_minimum(s);
                if (status) {
                    if (status==GSL_ENOPROG) { // GSL was unable to make improvement
                        gsl_multimin_fdfminimizer_restart(s);
                        if (restarts++<max_restarts) {
                            _ss << "Unable to make improvement, restarting... (attempt " << restarts << ")";
                            debug();
                            continue;
                        }
                        //warn("unable to make improvement, aborting search");
                        break; // cannot continue;                        
                    }
                    _ss << gettext("GSL minimizer") << " " << gsl_multimin_fdfminimizer_name(s)
                        << " " << gettext("failed with error code") << " " << status;
                    err();
                    opt.res=_NLP_ERROR;
                    break;
                }
                restarts=0;
                status=gsl_multimin_test_gradient(s->gradient,parm.eps);
                if (status==GSL_SUCCESS) {
                    _ss << "Converged to point " << opt.x; debug();
                    opt.res=_NLP_OPTIMAL;
                    break;
                }
            } while (_iter_count<maxiter);
        }
        res.push_back(opt);
    }
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(xstart);
    return true;
}
#endif // HAVE_LIBGSL

/* 
 * INTERIOR POINT IMPLEMENTATION (IPT_SOLVER CLASS)
 */
nlp_problem::ipt_solver::ipt_solver(nlp_problem &p,const meth_parm &parm) : prob(p), ctx(p.ctx) {
    eps_tol=parm.eps;
    feas_tol=parm.tol;
    maxiter=parm.max_iter;
    initialized=false;
    mu=temp_symb("mu",-1,ctx);
#ifdef HAVE_LIBLAPACK
    factorization=NULL;
    lapack_work=NULL;
    lapack_ipiv=NULL;
    lapack_rhs=NULL;
#endif
    init_parameters();
}
nlp_problem::ipt_solver::~ipt_solver() {
#ifdef HAVE_LIBLAPACK
    if (factorization!=NULL) delete[] factorization;
    if (lapack_work!=NULL) delete[] lapack_work;
    if (lapack_ipiv!=NULL) delete[] lapack_ipiv;
    if (lapack_rhs!=NULL) delete[] lapack_rhs;
#endif
}
/* initialize internal parameters with values from the Waechter & Biegler paper */
void nlp_problem::ipt_solver::init_parameters() {
    kappa_eps=10;
    kappa_mu=0.2;
    kappa_Sigma=1e10;
    kappa_soc=0.1; // not specified in the paper, guessing...
    kappa_1=1e-2;
    tau_min=0.99;
    theta_mu=1.5;
    delta=1;
    s_max=100;
    s_theta=1.1;
    s_phi=2.3;
    eta_phi=1e-4;
    gamma_alpha=0.05;
    gamma_theta=1e-5;
    gamma_phi=1e-5;
    p_max=4;
    bardelta_w_min=1e-20;
    bardelta_w_0=1e-4;
    bardelta_w_max=1e40;
    bardelta_c=1e-8;
    kappa_w_minus=1.0/3.0;
    kappa_w_plus=0;
    barkappa_w_plus=100;
    kappa_c=0.25;
    kappa_resto=0.9;
    lambda_max=1e3;
    theta_min=theta_max=undef;
}
/* substitute variables in EXPR */
vecteur nlp_problem::ipt_solver::subs_vars(const vecteur &expr,const vecteur &x0,const vecteur &lambda0,const vecteur &z0,double muj) {
    vecteur res=expr;
    if (!x0.empty()) res=subst(res,x,x0,false,ctx);
    if (!lambda0.empty()) res=subst(res,lambda,lambda0,false,ctx);
    if (!z0.empty()) res=subst(res,z,z0,false,ctx);
    if (muj>=0) res=subst(res,mu,muj,false,ctx);
    return res;
}
/* compute the barrier function value */
gen nlp_problem::ipt_solver::barrier_fval(const vecteur &x0) const {
    return _evalf(subst(f,x,x0,false,ctx)-mu_j*ln(_product(x0,ctx),ctx),ctx);
}
/* compute constraint violation */
gen nlp_problem::ipt_solver::constraint_violation(const vecteur &x0) const {
    return _evalf(_l1norm(subst(c,x,x0,false,ctx),ctx),ctx);
}
/* initialize filter */
void nlp_problem::ipt_solver::init_filter() {
    F.clear();
    F.resize(1,makevecteur(theta_max,minus_inf));
}
/* update filter with a new pair (THETA,PHI) */
void nlp_problem::ipt_solver::update_filter(const gen &b_theta,const gen &b_phi) {
    F.push_back(makevecteur(b_theta,b_phi));
}
/* Optimality error for the barrier problem with mu = mu0 */
gen nlp_problem::ipt_solver::E(const gen &mu0,const vecteur &x0,const vecteur &lambda0,const vecteur &z0) {
    gen s_d=max(s_max,(_l1norm(lambda0,ctx)+_l1norm(z0,ctx))/(nvars+ncons),ctx)/s_max;
    gen s_c=max(s_max,_l1norm(z0,ctx)/nvars,ctx)/s_max;
    vecteur e1v=subvecteur(grad_f,ncons>0?subvecteur(z0,multmatvecteur(A,lambda0)):z0),cx;
    e1v.resize(nvars-nslacks);
    e1v=subst(e1v,x,x0,false,ctx);
    cx=subst(c,x,x0,false,ctx);
    gen ret=_max(makevecteur(linfn(e1v)/s_d,linfn(subvecteur(pprod(x0,z0),vecteur(nvars,mu0)))/s_c,linfn(cx)),ctx);
    if (!is_real_number(ret,ctx)) {
        prob.debug("Failed to compute merit function");
        return undef;
    }
    return to_real_number(ret,ctx);
}
/* Inertia correction of the KKT matrix */
bool nlp_problem::ipt_solver::IC(const matrice &mat) {
    int n_pos,n_neg,n_zero,i,j;
    delta_w=delta_c=0;
    bool is_first=true;
    while (true) {
#ifdef HAVE_LIBLAPACK
        const_iterateur it=mat.begin(),itend=mat.end(),jt,jtend;
        gen e;
        for (i=0;it!=itend;++it,++i) {
            jt=it->_VECTptr->begin();
            jtend=it->_VECTptr->end();
            for(j=0;jt!=jtend;++jt,++j) {
                e=evalf_double(*jt,1,ctx);
                if (i==j) e+=i<nvars?delta_w:-delta_c;
                if (e.type!=_DOUBLE_)
                    return false;
                factorization[i+j*(nvars+ncons)]=e._DOUBLE_val;
            }
        }
        if (!solve_indef(factorization,&lapack_work,lapack_ipiv,NULL,nvars+ncons,0,&n_pos,&n_neg,&n_zero))
            return false;
#else
        copy_matrice(mat,factorization);
        for (i=0;i<nvars;++i) factorization[i]._VECTptr->at(i)+=delta_w;
        for (i=0;i<ncons;++i) factorization[i+nvars]._VECTptr->at(i+nvars)-=delta_c;
        if (!solve_indef(factorization,NULL,ldl_perm,&n_pos,&n_neg,&n_zero,ctx))
            return false;
#endif
        if (n_zero==0 && n_pos==nvars)
            break;
        if (is_first) {
            delta_w=delta_w_last>0?std::max(bardelta_w_min,kappa_w_minus*delta_w_last):bardelta_w_0;
            delta_c=n_zero>0?bardelta_c*std::pow(mu_j,kappa_c):0;
            is_first=false;
            continue;
        }
        delta_w=delta_w_last>0?kappa_w_plus*delta_w:barkappa_w_plus*delta_w;
        if (delta_w>bardelta_w_max)
            return false;
    }
    return true;
}
/* return true iff (THETA,PHI) is accepted by the current filter */
bool nlp_problem::ipt_solver::filter_accepts(const gen &theta,const gen &phi) const {
    const_iterateur it=F.begin(),itend=F.end();
    for (;it!=itend;++it) {
        if (is_greater(theta,it->_VECTptr->front(),ctx) && is_greater(phi,it->_VECTptr->back(),ctx))
            return false;
    }
    return true;
}
/* initialize the symbolic KKT system for the main iteration */
bool nlp_problem::ipt_solver::kkt_init() {
    matrice W(hess_f);
    const_iterateur it=A.begin(),itend=A.end();
    for (;it!=itend;++it) {
        matrice J;
        if (!jacobian(*it->_VECTptr,x,J,ctx)) {
            prob.err("failed to compute second-order constraint derivatives");
            return false;
        }
        W=madd(W,mscale(lambda[it-A.begin()],J));
    }
    sigma=pprod(*inv(x,ctx)._VECTptr,z);
    if (ncons>0) {
        matrice tmp=mtran(A);
        kkt_mat=*_blockmatrix(makesequence(2,2,makevecteur(W,tmp,A,zero_mat(ncons,ncons,ctx))),ctx)._VECTptr;
        A=tmp;
    } else kkt_mat=W;
    for (int i=0;i<nvars;++i) kkt_mat[i]._VECTptr->at(i)+=sigma[i];
    kkt_rhs=multvecteur(-1,mergevecteur(subvecteur(ncons>0?addvecteur(grad_f,multmatvecteur(A,lambda)):grad_f,z),c));
    return true;
}
/* initialize lambda values for the main iteration */
vecteur nlp_problem::ipt_solver::init_lambda(const vecteur &x0,const vecteur &z0) {
    if (ncons==0)
        return vecteur(0);
    matrice mat=subst(*_blockmatrix(makesequence(2,2,makevecteur(midn(nvars),A,mtran(A),zero_mat(ncons,ncons,ctx))),ctx)._VECTptr,x,x0,false,ctx);
    vecteur rhs=multvecteur(-1,mergevecteur(subvecteur(subst(grad_f,x,x0,false,ctx),z0),vecteur(ncons,0)));
    vecteur sol;
    if (!solve_lsq(mat,rhs,sol,ctx))
        return vecteur(ncons,0);
    vecteur lambda0(sol.begin()+nvars,sol.end());
    if (is_strictly_greater(linfn(lambda0),lambda_max,ctx))
        return vecteur(ncons,0);
    return lambda0;
}
/* linear system solver using LSQ */
bool nlp_problem::solve_lsq(const matrice &lhs,const vecteur &rhs,vecteur &sol,GIAC_CONTEXT) {
    gen res//=_linsolve(makesequence(lhs,rhs),contextptr)
    ;
    //return res.type==_VECT && (sol=*res._VECTptr).size()==rhs.size();
    try {
        log_output_redirect lor(contextptr);
        res=_LSQ(makesequence(lhs,rhs),contextptr);
    } catch (const std::runtime_error &e) {
        return false;
    }
    return ckmatrix(res) && res._VECTptr->front()._VECTptr->size()==1 && (sol=*mtran(*res._VECTptr).front()._VECTptr).size()==rhs.size();
}
/* Feasibility restoration phase algorithm.
 * We simply solve the relaxed Lagrangian sub-problem without
 * nonnegativity constraints, starting from the current point X0.
 * If the solution is not componentwise nonnegative, or the solution
 * does not significantly reduce the infeasibility, the routine
 * returns FAILED, otherwise it returns OPTIMAL. */
int nlp_problem::ipt_solver::feas_restoration(vecteur &x0,vecteur &lambda0,vecteur &z0,const gen &viol_old,gen &viol) {
    //prob.debug("Attempting feasibility restoration...");
    vecteur d_R=*_apply(makesequence(at_min,mtran(makevecteur(vecteur(nvars,1.0),inv(_abs(x0,ctx),ctx)))),ctx)._VECTptr;
    vecteur d_R2=pprod(d_R,d_R);
    vecteur L_eq=multvecteur(2,pprod(d_R2,subvecteur(x,x0)));
    if (ncons>0) L_eq=addvecteur(L_eq,multmatvecteur(A,lambda));
    gen sol;
    try {
        sol=_fsolve(makesequence(mergevecteur(L_eq,c),mergevecteur(x,lambda),mergevecteur(x0,lambda0)),ctx);
    } catch (const std::runtime_error &e) {
        prob.err(e.what());
        return _NLP_ERROR;
    }
    if (sol.type!=_VECT || int(sol._VECTptr->size())<nvars)
        return _NLP_ERROR;
    vecteur x_k(sol._VECTptr->begin(),sol._VECTptr->begin()+nvars);
    if (!is_strictly_positive(_min(x_k,ctx),ctx)) {
        viol=viol_old;
        return _NLP_FAILED;
    }
    /* compute new values of z and lambda */
    vecteur d_x=subvecteur(x_k,x0),d_z;
    compute_dz(x_k,d_x,z0,d_z);
    gen alpha_z=fraction_to_boundary(z0,d_z);
    z0=addvecteur(z0,multvecteur(alpha_z,d_z));
    lambda0=init_lambda(x0,z0);
    x0=x_k;
    viol=constraint_violation(x0);
    bool yes=filter_accepts(viol,barrier_fval(x0)) && is_greater(kappa_resto*viol_old,viol,ctx);
    return yes?_NLP_OPTIMAL:_NLP_FAILED;
}
/* computing the movement in the direction of the z-variables */
void nlp_problem::ipt_solver::compute_dz(const vecteur &x0,const vecteur &dx0,const vecteur &z0,vecteur &dz) {
    dz=subs_vars(subvecteur(multvecteur(mu_j,*inv(x0,ctx)._VECTptr),addvecteur(z0,pprod(sigma,dx0))),x0,vecteur(0),z0);
}
/* reducing the step size */
gen nlp_problem::ipt_solver::fraction_to_boundary(const vecteur &y,const vecteur &d_y) {
    gen alpha(1);
    const_iterateur it=d_y.begin(),itend=d_y.end();
    for (;it!=itend;++it) {
        if (is_strictly_positive(-(*it),ctx)) {
            gen cand=-tau_j*y[it-d_y.begin()]/(*it);
            if (is_strictly_greater(alpha,cand,ctx))
                alpha=cand;
        }
    }
    return alpha;
}
bool nlp_problem::ipt_solver::solve_kkt(const matrice &mat,const vecteur &rh,vecteur &sol) {
#ifdef HAVE_LIBLAPACK
    const_iterateur it=rh.begin(),itend=rh.end();
    for (int i=0;it!=itend;++it,++i) {
        gen e=evalf_double(*it,1,ctx);
        if (e.type!=_DOUBLE_) return _NLP_ERROR;
        lapack_rhs[i]=e._DOUBLE_val;
    }
    if (!solve_indef(factorization,&lapack_work,lapack_ipiv,lapack_rhs,nvars+ncons,1,NULL,NULL,NULL))
        return false;
    iterateur jt=sol.begin(),jtend=sol.end();
    for (int i=0;jt!=jtend;++jt,++i) *jt=lapack_rhs[i];
#else
    std::copy(ldl_perm.begin(),ldl_perm.end(),sol.begin());
    if (!solve_indef(factorization,&rh,sol,NULL,NULL,NULL,ctx))
        return false;
#endif
    return true;
}
/* Interior-point iterative solver */
int nlp_problem::ipt_solver::ls_filter_barrier_method(vecteur &x_k,vecteur &lambda_k,vecteur &z_k,double mu0) {
    /* A-1. Initialize */
    int j=0,k=-1,heur_iter=0,heur_threshold=4,skip_iter=0,skip_threshold=4;
    delta_w_last=0;
    init_filter();
    mu_j=mu0;
    tau_j=std::max(tau_min,1-mu_j);
    bool prev_iter_full_step=false,converged=false,feas_restored=false,first_trial_stepsize_rejected,accepted;
    vecteur pts,p_k,n_k,z_pk,z_nk,d_p,d_n,d_zp,d_zn;
    while (prob._iter_count++<maxiter) {
        ++k;
        if (prob._iter_data) pts.push_back(x_k);
        /* A-2. Check convergence for the overall problem */
        gen merit=E(0,x_k,lambda_k,z_k);
        if (is_undef(merit)) return _NLP_ERROR;
        //sprintf(buf,"Iteration %d: %g",k,merit.to_double(ctx));
        //prob.debug(buf);
        if (is_greater(eps_tol,merit,ctx)) {
            converged=true;
            break;
        }
        /* A-3. Check convergence for the barrier problem */
        while (true) {
            merit=E(mu_j,x_k,lambda_k,z_k);
            if (is_undef(merit)) return _NLP_ERROR;
            if (is_strictly_greater(merit,kappa_eps*mu_j,ctx))
                break;
            mu_j=std::max(eps_tol/10,std::min(kappa_mu*mu_j,std::pow(mu_j,theta_mu)));
            tau_j=std::max(tau_min,1-mu_j);
            ++j;
            init_filter();
            if (k>0) break;
        }
        /* A-4. Update KKT matrix and compute the search direction */
        matrice kktm=subs_vars(kkt_mat,x_k,lambda_k,z_k,mu_j);
        theta_k=constraint_violation(x_k);
        if (is_undef(theta_min)) theta_min=1e-4*max(1,theta_k,ctx);
        if (is_undef(theta_max)) theta_max=1e4*max(1,theta_k,ctx);
        phi_k=barrier_fval(x_k);
        if (!IC(kktm)) {
            prob.warn("KKT matrix is ill-conditioned");
            if (ncons==0)
                return _NLP_FAILED;
            update_filter((1-gamma_theta)*theta_k,subst(f,x,x_k,false,ctx)-gamma_phi*theta_k);
            switch (feas_restoration(x_k,lambda_k,z_k,theta_k,theta_resto)) {
            case _NLP_ERROR:
                //prob.debug("Feasibility restoration phase failed");
                return _NLP_ERROR;
            case _NLP_FAILED:
                return is_greater(theta_resto,feas_tol,ctx)?_NLP_INFEAS:_NLP_FAILED;
            default:
                //prob.debug("Feasibility restored successfully");
                continue;
            }
        }
        vecteur kkt_sol(nvars+ncons),rhs=subs_vars(kkt_rhs,x_k,lambda_k,z_k,mu_j);
        if (!solve_kkt(kktm,rhs,kkt_sol)) {
            //prob.debug("Failed to obtain search directions");
            return _NLP_ERROR;
        }
        vecteur d_x(kkt_sol.begin(),kkt_sol.begin()+nvars),x_kn;
        vecteur d_lambda(kkt_sol.begin()+nvars,kkt_sol.end()),d_z;
        compute_dz(x_k,d_x,z_k,d_z);
        /* A-5. Backtracking line search */
        gen alpha_max=fraction_to_boundary(x_k,d_x),alpha_z=fraction_to_boundary(z_k,d_z);
        /* A-5.1. Initialize the line search */
        gen alpha_kl=alpha_max,alpha_k;
        int l=0;
        bool skip=false;
        vecteur grad_phi=subvecteur(subst(subst(grad_f,x,x_k,false,ctx),mu,mu_j,false,ctx),z_k);
        gen sp=-scalarproduct(grad_phi,d_x,ctx);
        gen alpha_min=ncons==0?gamma_alpha:gamma_alpha*(is_strictly_positive(sp,ctx)?
            _min(makesequence(gamma_theta,gamma_phi*theta_k/sp,delta*pow(theta_k,s_theta,ctx)/pow(sp,s_phi,ctx)),ctx):gen(gamma_theta));
        alpha_min=max(min(alpha_min,alpha_max/100,ctx),eps_tol,ctx);
        first_trial_stepsize_rejected=false;
        while (true) {
            if (ctrl_c || interrupted) {
                ctrl_c=interrupted=false;
                throw std::runtime_error("Stopped by user interruption");
            }
            /* A-5.2. Compute the new trial point */
            vecteur x_kl=addvecteur(x_k,multvecteur(alpha_kl,d_x));
            gen theta=constraint_violation(x_kl),phi=barrier_fval(x_kl);
            bool switching_condition=is_strictly_positive(sp,ctx) && is_strictly_greater(alpha_kl*pow(sp,s_phi,ctx),delta*pow(theta_k,s_theta,ctx),ctx);
            /* A-5.3. Check acceptability to the filter */
            accepted=false;
            while (!accepted) {
                if (filter_accepts(theta,phi)) {
                    /* A-5.4. Check sufficient decrease with respect to the current iterate */
                    if ((is_greater(theta_min,theta_k,ctx) && switching_condition && is_greater(phi_k,phi+eta_phi*alpha_kl*sp,ctx)) ||
                            ((is_strictly_greater(theta_k,theta_min,ctx) || !switching_condition) &&
                                (is_greater((1-gamma_theta)*theta_k,theta,ctx) || is_greater(phi_k-gamma_theta*theta_k,phi,ctx)))) {
                        alpha_k=alpha_kl;
                        x_kn=x_kl;
                        accepted=true;
                    } else break;
                } else if (heur_iter==heur_threshold) { // apply accelerating heuristic
                    //prob.debug("Applying accelerating heuristic...");
                    theta_max=0.1*theta_max;
                    init_filter();
                    heur_iter=0;
                } else break;
            }
            if (accepted)
                break;
            if (l==0 && ncons>0 && is_greater(theta,theta_k,ctx)) {
                //prob.debug("Applying second-order correction...");
                /* A-5.5. Initialize the second-order correction */
                int p=1; // SOC counter
                vecteur c_soc=addvecteur(multvecteur(alpha_kl,subst(c,x,x_k,false,ctx)),subst(c,x,x_kl,false,ctx));
                gen theta_old_soc=theta_k;
                while (true) {
                    /* A-5.6. Compute the second-order correction */
                    rhs.resize(nvars);
                    rhs=mergevecteur(rhs,c_soc);
                    if (!solve_kkt(kktm,rhs,kkt_sol))
                        return _NLP_ERROR;
                    vecteur d_x_cor(kkt_sol.begin(),kkt_sol.begin()+nvars);
                    gen alpha_soc=fraction_to_boundary(x_k,d_x_cor);
                    vecteur x_soc=addvecteur(x_k,multvecteur(alpha_soc,d_x_cor));
                    /* A-5.7. Check acceptability to the filter (in SOC) */
                    theta=constraint_violation(x_soc);
                    phi=barrier_fval(x_soc);
                    if (filter_accepts(theta,phi)) {
                        /* A-5.8. Check sufficient decrease with respect to the current iterate (in SOC) */
                        if ((is_greater(theta_min,theta_k,ctx) && switching_condition && is_greater(phi_k,phi+eta_phi*alpha_kl*sp,ctx)) ||
                                ((is_strictly_greater(theta_k,theta_min,ctx) || !switching_condition) &&
                                    (is_greater((1-gamma_theta)*theta_k,theta,ctx) || is_greater(phi_k-gamma_theta*theta_k,phi,ctx)))) {
                            alpha_k=alpha_soc;
                            x_kn=x_soc;
                            accepted=true;
                            break;
                        }
                    }
                    /* A-5.9. Next second-order correction */
                    if (p==p_max || is_strictly_greater(theta,kappa_soc*theta_old_soc,ctx))
                        break;
                    p++;
                    c_soc=addvecteur(multvecteur(alpha_soc,c_soc),subst(c,x,x_soc,false,ctx));
                    theta_old_soc=theta;
                }
                if (accepted)
                    break;
            }
            /* A-5.10. Choose the new trial step size */
            alpha_kl=0.5*alpha_kl;
            first_trial_stepsize_rejected=true;
            l++;
            if (is_strictly_greater(alpha_min,alpha_kl,ctx)) {
                update_filter((1-gamma_theta)*theta_k,subst(subst(f,x,x_k,false,ctx),mu,mu_j,false,ctx)-gamma_phi*theta_k);
                switch (feas_restoration(x_k,lambda_k,z_k,theta_k,theta_resto)) {
                case _NLP_ERROR:
                    //prob.debug("Feasibility restoration phase failed");
                    return _NLP_ERROR;
                case _NLP_FAILED:
                    if (is_greater(theta_resto,feas_tol,ctx))
                        return _NLP_INFEAS;
                    break;
                default:
                    //prob.debug("Feasibility restored successfully");
                    break;
                }
                skip=true;
                break;
            }
        }
        if (first_trial_stepsize_rejected)
            heur_iter++;
        else heur_iter=0;
        if (skip) {
            if (skip_iter++==skip_threshold)
                return _NLP_FAILED;
            continue;
        } else skip_iter=0;
        /* A-6. Accept the trial point */
        int xsize=x_k.size();
        vecteur d_rel(xsize);
        for (int i=0;i<xsize;++i) d_rel[i]=(x_k[i]-x_kn[i])/(1+_abs(x_k[i],ctx));
        bool fullstep=is_strictly_greater(10*epsilon(ctx),linfn(d_rel),ctx);
        if (fullstep) { // handle small search directions
            alpha_k=alpha_max;
            x_kn=addvecteur(x_k,multvecteur(alpha_k,d_x));
        }
        vecteur lambda_kn=ncons>0?addvecteur(lambda_k,multvecteur(alpha_k,d_lambda)):vecteur(0);
        vecteur z_kn=addvecteur(z_k,multvecteur(alpha_z,d_z));
        /* A-7. Augment the filter if necessary */
        if (is_positive(-sp,ctx) || is_greater(delta*pow(theta_k,s_theta,ctx),alpha_k*pow(sp,s_phi,ctx),ctx) ||
                is_strictly_greater(barrier_fval(x_kn),phi_k-eta_phi*alpha_k*sp,ctx))
            update_filter((1-gamma_theta)*theta_k,subst(subst(f,x,x_k,false,ctx),mu,mu_j,false,ctx)-gamma_phi*theta_k);
        /* A-8. Update search directions and continue with the next iteration */
        if (fullstep) {
            if (prev_iter_full_step) {
                if (mu_j>eps_tol/10) {
                    //prob.debug("Search directions are too small, terminating...");
                    converged=true;
                    break;
                }
                else return _NLP_FAILED;
            } else prev_iter_full_step=true;
        } else prev_iter_full_step=false;
        x_k=x_kn;
        lambda_k=lambda_kn;
        z_k=z_kn;
    }
    if (prob._iter_data)
        iteration_data.push_back(pts);
    if (!converged) {
        //prob.debug("Iteration limit exceeded");
        return _NLP_FAILED;
    }
    return _NLP_OPTIMAL;
}
/* convert nlp_problem variables to ipt_solver variables */
vecteur nlp_problem::ipt_solver::vars2x(const vecteur &x0,const vecteur &x0_orig) const {
    int pos=0;
    vecteur ret(nvars);
    const_iterateur jt=x0.begin(),jtend=x0.end();
    for (;jt!=jtend;++jt) {
        if (dup_ind.find(pos)!=dup_ind.end()) {
            if (x0_orig.empty()) {
                if (is_positive(*jt,ctx)) {
                    ret[pos]=*jt+kappa_1;
                    ret[pos+1]=kappa_1;
                } else {
                    ret[pos]=kappa_1;
                    ret[pos+1]=-*jt+kappa_1;
                }
            } else {
                gen a=x0_orig[pos],b=x0_orig[pos+1],val1=(a+b+*jt)/2,val2=(a+b-*jt)/2,mv=min(val1,val2,ctx);
                if (is_positive(-mv,ctx)) {
                    val1+=kappa_1-mv;
                    val2+=kappa_1-mv;
                }
                ret[pos]=max(kappa_1,val1,ctx);
                ret[pos+1]=max(kappa_1,val2,ctx);
            }
            pos+=2;
        } else {
            gen val=_solve(makesequence(symb_equal(x_subs[jt-x0.begin()],*jt),x[pos]),ctx)._VECTptr->front();
            ret[pos]=max(kappa_1,val,ctx);
            pos++;
        }
    }
    int pos0=pos;
    vecteur ns(x.begin(),x.begin()+pos0),nsval(ret.begin(),ret.begin()+pos0);
    for (;pos<nvars;++pos) {
        gen val=_solve(makesequence(symb_equal(subst(c[pos-pos0],ns,nsval,false,ctx),0),x[pos]),ctx)._VECTptr->front();
        ret[pos]=max(kappa_1,val,ctx);
    }
    return ret;
}
/* convert ipt_solver variables to nlp_problem variables */
vecteur nlp_problem::ipt_solver::x2vars(const vecteur &x0) const {
    return subst(x_subs,x,x0,false,ctx);
}
/* the main optimization routine */
bool nlp_problem::ipt_solver::optimize(optima_t &res) {
    if (!initialized) {
        /* This code transforms the problem into min f(x) s.t. h(x)=0, x>=0.
         * Inequalities come first in the list C of constraints, sclacks are used for
         * conversion to equalities.
         * Unbounded variables are replaced by x_pos-x_neg, where x_pos and x_neq
         * are both nonnegative. Additionally, the constraints Jacobian J_C is computed. */
        int nv_orig=prob.var_count(),nv=nv_orig;
        if (nv==0) {
            optimum_t trivial;
            trivial.res=_NLP_OPTIMAL;
            prob.insert_fixed_vars(trivial.x);
            res.push_back(trivial);
            return true;            
        }
        x.reserve(2*nv);
        x_subs.resize(nv);
        vecteur c_eq=prob.subs_fixed_vars(prob.eq_cons),c_ineq=prob.subs_fixed_vars(prob.ineq_cons),lb=prob.lb,ub=prob.ub;
        matrice J_eq=prob.subs_fixed_vars(prob.eq_jacobian),J_ineq=prob.subs_fixed_vars(prob.ineq_jacobian);
        prob._inactive_eq_indices.clear();
        prob._inactive_ineq_indices.clear();
        if (!prob._fixed_vars.empty() && (!prob.find_inactive_eq(c_eq) || !prob.find_inactive_ineq(c_ineq))) {
            // the problem is infeasible
            optimum_t opt;
            opt.res=_NLP_INFEAS;
            res.push_back(opt);
            return true;
        }
        remove_elements_with_indices(lb,prob._fixed_vars);
        remove_elements_with_indices(ub,prob._fixed_vars);
        remove_elements_with_indices(c_eq,prob._inactive_eq_indices);
        remove_elements_with_indices(c_ineq,prob._inactive_ineq_indices);
        remove_elements_with_indices(J_eq,prob._inactive_eq_indices);
        remove_elements_with_indices(J_ineq,prob._inactive_ineq_indices);
        int ineq_size=c_ineq.size();
        ncons=c_eq.size()+c_ineq.size();
        c=mergevecteur(c_ineq,c_eq);
        A=mergevecteur(J_ineq,J_eq);
        grad_f=prob.subs_fixed_vars(prob.obj_gradient);
        remove_elements_with_indices(grad_f,prob._fixed_vars);
        hess_f=prob.subs_fixed_vars(prob.obj_hessian);
        if (!prob._fixed_vars.empty()) {
            if (!A.empty()) {
                A=mtran(A);
                remove_elements_with_indices(A,prob._fixed_vars);
                A=mtran(A);
            }
            remove_elements_with_indices(hess_f,prob._fixed_vars);
            hess_f=mtran(hess_f);
            remove_elements_with_indices(hess_f,prob._fixed_vars);
            hess_f=mtran(hess_f);
        }
        gen L,U;
        vecteur vars(prob.vars);
        remove_elements_with_indices(vars,prob._fixed_vars);
        const_iterateur it=vars.begin(),itend=vars.end();
        for (;it!=itend;++it) {
            int i=it-vars.begin();
            if (!is_inf(L=lb[i])) {
                x.push_back(temp_symb("x",i,ctx));
                x_subs[i]=L+x.back();
                if (!is_inf(U=ub[i])) {
                    c.insert(c.begin()+ineq_size,U-*it);
                    A.insert(A.begin()+ineq_size,vecteur(nv_orig,0));
                    A[ineq_size]._VECTptr->at(i)=-1;
                    ++ineq_size;
                    ++ncons;
                }
            } else if (!is_inf(U=ub[i])) {
                x.push_back(temp_symb("x",i,ctx));
                x_subs[i]=U-x.back();
                grad_f[i]=-grad_f[i];
                for (int j=0;j<ncons;++j) A[j]._VECTptr->at(i)=-A[j][i];
                for (int j=0;j<nv_orig;++j) {
                    hess_f[j]._VECTptr->at(i)=-hess_f[j][i];
                    if (i!=j) hess_f[i]._VECTptr->at(j)=hess_f[j][i];
                }
            } else { // unbounded variable
                x.push_back(temp_symb("x_pos",i,ctx));
                x.push_back(temp_symb("x_neg",i,ctx));
                x_subs[i]=x[x.size()-2]-x.back();
                dup_ind.insert(x.size()-2);
                ++nv;
            }
        }
        if (ncons>0 && !A.empty()) A=mtran(A);
        for (set<int>::const_iterator it=dup_ind.begin();it!=dup_ind.end();++it) {
            grad_f.insert(grad_f.begin()+*it+1,-grad_f[*it]);
            hess_f.insert(hess_f.begin()+*it+1,multvecteur(-1,*hess_f[*it]._VECTptr));
            if (ncons>0) A.insert(A.begin()+*it+1,multvecteur(-1,*A[*it]._VECTptr));
        }
        grad_f.resize(nv+ineq_size,0);
        hess_f=mtran(ineq_size>0?mergevecteur(hess_f,zero_mat(ineq_size,mcols(hess_f),ctx)):hess_f);
        for (set<int>::const_iterator it=dup_ind.begin();it!=dup_ind.end();++it) {
            hess_f.insert(hess_f.begin()+*it+1,multvecteur(-1,*hess_f[*it]._VECTptr));
        }
        if (ineq_size>0) {
            hess_f=mergevecteur(hess_f,zero_mat(ineq_size,mcols(hess_f),ctx));
            A=mergevecteur(A,zero_mat(ineq_size,mcols(A),ctx));
        }
        if (ncons>0 && !A.empty()) A=mtran(subst(A,vars,x_subs,false,ctx));
        grad_f=subst(grad_f,vars,x_subs,false,ctx);
        hess_f=subst(hess_f,vars,x_subs,false,ctx);
        if (ncons>0) c=subst(c,vars,x_subs,false,ctx);
        nslacks=ineq_size;
        for (int i=0;i<nslacks;++i) {
            x.push_back(temp_symb("slack",i,ctx));
            c[i]-=x.back();
            if (ncons>0) A[i]._VECTptr->at(nv+i)=-1;
        }
        nvars=nv+nslacks;
        z.resize(nvars);
        for (int i=0;i<nvars;++i) z[i]=temp_symb("z",i,ctx);
        lambda.resize(ncons);
        for (int i=0;i<ncons;++i) lambda[i]=temp_symb("lambda",i,ctx);
        f=prob.subs_fixed_vars(prob.obj);
        f=subst(f,vars,x_subs,false,ctx);
        initialized=true;
#ifdef HAVE_LIBLAPACK
        int fsize=nvars+ncons;
        factorization=new double[fsize*fsize];
        lapack_ipiv=new int[fsize];
        lapack_rhs=new double[fsize];
#endif
    }
    if (!kkt_init()) {
        prob.err("failed to initialize KKT system");
        return false;
    }
    vecteur ip;
    double mu0=0.1;
    while (prob._iter_count<maxiter && prob.restart(ip)) {
        optimum_t opt;
        remove_elements_with_indices(ip,prob._fixed_vars);
        vecteur x_k=vars2x(ip),z_k(nvars,1),lambda_k=init_lambda(x_k,z_k);
        opt.res=ls_filter_barrier_method(x_k,lambda_k,z_k,mu0);
        opt.x=x2vars(x_k);
        prob.insert_fixed_vars(opt.x);
        if (opt.res<_NLP_INFEAS)
            prob.save_point(opt.x);
        res.push_back(opt);
    }
    return true;
}
/* END OF IPT_SOLVER CLASS IMPLEMENTATION */

/* Return true if there are no constraints/bounds. */
bool nlp_problem::is_unconstrained() const {
    if (eq_count()+ineq_count()>0)
        return false;
    for (int i=0;i<var_count();++i) {
        if (!is_inf(lb[i]) || !is_inf(ub[i]))
            return false;
    }
    return true;
}
/* return true if some variables are bounded from either side, else return false. */
bool nlp_problem::has_bounded_variables() const {
    int n=var_count();
    for (int i=0;i<n;++i) {
        if (!is_inf(lb[i]) || !is_inf(ub[i]))
            return true;
    }
    return false;
}
/* Compute the index of the branching variable for B&B algorithm.
 * Possible strategies are: 0 (random), 1 (first), 2 (last), 3 (pseudocost).
 * If there is no variable to branch upon, return false, else return true. */
bool nlp_problem::branch_var(const meth_parm &parm,const vecteur &sol,const gen &cost,int strategy,int &var) {
    vector<int> iv;
    iv.reserve(sol.size());
    const_iterateur it=sol.begin(),itend=sol.end();
    vector<bool>::const_iterator jt=_is_intvar.begin();
    vector<int>::const_iterator kt,ktend;
    int i=0;
    for (;it!=itend;++it,++jt,++i) {
        if (*jt && is_strictly_greater(_abs(*it-_round(*it,ctx),ctx),parm.tol,ctx))
            iv.push_back(i);
    }
    if (iv.empty())
        return false; // no fractional variables
    gen pc,maxpc=minus_inf;
    switch (strategy) {
    case 0:
        var=iv[_rand(iv.size(),ctx).val];
        break;
    case 1:
        var=iv.front();
        break;
    case 2:
        var=iv.back();
        break;
    case 3:
        if (_min(pc_lo_count,ctx).val==0 || _min(pc_hi_count,ctx).val==0) // there is a variable not being branched upon yet
            var=iv.front(); // branch upon the first fractional variable
        else for (kt=iv.begin(),ktend=iv.end();kt!=ktend;++kt) { // use pseudocost branching
            i=var2intvar(*kt);
            pc=cost+max(pc_lo[i]*(sol[*kt]-_floor(sol[*kt],ctx)),pc_hi[i]*(_ceil(sol[*kt],ctx)-sol[*kt]),ctx);
            if (is_strictly_greater(pc,maxpc,ctx)) {
                maxpc=pc;
                var=*kt;
            }
        }
        break;
    default:
        assert(false); // should be unreachable
    }
    return true;
}
/* Return true if X is a feasible point with respect to the constraints. */
bool nlp_problem::is_feasible(const vecteur &x,double feas_tol,bool check_bounds,set<int> *viol_eq,set<int> *viol_ineq,set<int> *viol_bnd) const {
    gen max_infeas(0),infeas;
    const_iterateur it,itend;
    int i,n=var_count();
    if (viol_eq!=NULL) viol_eq->clear();
    if (viol_ineq!=NULL) viol_ineq->clear();
    for (i=0,it=eq_cons.begin(),itend=eq_cons.end();it!=itend;++it,++i) {
        infeas=_abs(subst(*it,vars,x,false,ctx),ctx);
        if (is_strictly_greater(infeas,feas_tol,ctx)) {
            max_infeas=max(max_infeas,infeas);
            if (viol_eq!=NULL)
                viol_eq->insert(i);
        }
    }
    for (i=0,it=ineq_cons.begin(),itend=ineq_cons.end();it!=itend;++it,++i) {
        infeas=max(0,-subst(*it,vars,x,false,ctx),ctx);
        if (is_strictly_greater(infeas,feas_tol,ctx)) {
            max_infeas=max(max_infeas,infeas);
            if (viol_ineq!=NULL)
                viol_ineq->insert(i);
        }
    }
    if (check_bounds) {
        if (viol_bnd!=NULL) viol_bnd->clear();
        for (i=0;i<n;++i) {
            infeas=max(0,max(lb[i]-x[i],x[i]-ub[i],ctx),ctx);
            if (is_strictly_greater(infeas,max_infeas,ctx)) {
                max_infeas=max(max_infeas,infeas);
                if (viol_bnd!=NULL)
                    viol_bnd->insert(i);
            }
        }
    }
    return is_strictly_greater(feas_tol,max_infeas,ctx);
}
/* Return the number of integer variables. */
int nlp_problem::intvar_count() const {
    int ret=0;
    for(vector<bool>::const_iterator it=_is_intvar.begin();it!=_is_intvar.end();++it)
        if (*it) ++ret;
    return ret;
}
/* Return true iff all integer variables are binary. */
bool nlp_problem::is_intvars_all_binary() const {
    int i=0;
    for (vector<bool>::const_iterator it=_is_intvar.begin();it!=_is_intvar.end();++it,++i) {
        if (*it && (!is_zero(lb[i]) || !is_one(ub[i])))
            return false;
    }
    return true;
}
/* Return true iff SOL satisfies integrality constraints. */
bool nlp_problem::is_intsol(const meth_parm &parm,const vecteur &sol) const {
    int i=0;
    const_iterateur it=sol.begin(),itend=sol.end();
    for (;it!=itend;++it,++i) {
        if (is_intvar(i) && is_strictly_greater(_abs(*it-_round(*it,ctx),ctx),parm.tol,ctx))
            return false;
    }
    return true;
}
/* Find fixed variables (with mutually equal lower and upper bound). */
void nlp_problem::find_fixed_vars() {
    _fixed_vars.clear();
    int n=var_count(),i=0;
    for (;i<n;++i) {
        if (is_zero(ub[i]-lb[i],ctx))
            _fixed_vars.insert(make_pair(i,(ub[i]+lb[i])/2.0));
    }
}
/* Insert fixed variables in X. */
void nlp_problem::insert_fixed_vars(vecteur &x) const {
    set<pair<int,gen> >::const_iterator it=_fixed_vars.begin(),itend=_fixed_vars.end();
    for (;it!=itend;++it)
        x.insert(x.begin()+it->first,it->second);
}
/* Branch & bound subproblem solving routine. */
bool nlp_problem::bb_solve_subproblem(const meth_parm &parm,vecteur &active_nodes,gen &fval,optima_t &sub_res,optimum_t &incumbent) {
    find_fixed_vars();
    initp_iter=initp.begin();
    _iter_count=0;
    ipt_solver *ipt=new ipt_solver(*this,parm);
    sub_res.clear();
    bool ret=false;
    int i;
    optimum_t r;
    ++_spc;
    clear_saved_points();
    if (ipt->optimize(sub_res) && !sub_res.empty()) {
        find_best_solution(parm,sub_res,r.x,r.f,r.res);
        if ((r.res==_NLP_OPTIMAL || r.res==_NLP_FAILED) && compute_obj_val(r.x,fval)) {
            ret=true; // pseudocosts will be updated
            r.f=fval;
            bool crit=incumbent.res==_NLP_PENDING || is_strictly_greater(incumbent.f,fval,ctx);
            if (is_intsol(parm,r.x)) { // a feasible solution was found
                for (i=active_nodes.size();i-->0;) { // delete nodes with too large lower bound
                    if (is_greater(active_nodes[i]._VECTptr->back(),fval,ctx))
                        active_nodes.erase(active_nodes.begin()+i);
                }
                if (crit) {
                    _ss << gettext("New incumbent solution found") << ": " << fval; msg();
                    incumbent=r;
                }
            } else if (crit)
                active_nodes.push_back(makevecteur(lb,ub,r.x,fval)); // append new active node
        }
    }
    delete ipt;
    return ret;
}
/* Return the index of VAR as an integer variable. */
int nlp_problem::var2intvar(int var) const {
    vector<bool>::const_iterator it=_is_intvar.begin();
    int i=0,j=0;
    for (;j!=var;++it,++j) if (*it) ++i;
    return i;
}
/* Update pseudocost for integer variable VAR. */
void nlp_problem::update_pseudocost(int var,vecteur &pc,vecteur &pc_count,const gen &df,const gen &dx) {
    int i=var2intvar(var);
    pc[i]=(pc_count[i]*pc[i]+df/dx)*_inv(pc_count[i]+=1,ctx);
}
/* Branch & bound algorithm for MINLP using the interior-point solver.
 * CONT_SOL is the solution if the continuous relaxation.
 * INCUMBENT status must be set to PENDING prior the call. */
bool nlp_problem::branch_and_bound(const meth_parm &parm,const vecteur &cont_sol,optimum_t &incumbent) {
    incumbent.res=_NLP_PENDING;
    incumbent.x.clear();
    incumbent.f=plus_inf;
    int niv=intvar_count(),bv;
    _spc=0;
    pc_lo.resize(niv,0); pc_lo_count.resize(niv,0);
    pc_hi.resize(niv,0); pc_hi_count.resize(niv,0);
    vecteur active_nodes;
    gen fval,minbnd,oldb,s_fval,s_xf,s_xl,s_xu;
    if (!compute_obj_val(cont_sol,fval))
        return false;
    vecteur root_node=makevecteur(lb,ub,cont_sol,fval);
    active_nodes.push_back(root_node);
    iterateur nt,ntend,sel;
    optima_t sub_res;
    while (!active_nodes.empty()) {
        if (interrupted || ctrl_c) {
            interrupted=ctrl_c=false;
            throw std::runtime_error(gettext("Stopped by user interruption"));
        }
        // node selection (strategy: best local bound)
        minbnd=plus_inf;
        sel=active_nodes.end();
        for (nt=active_nodes.begin(),ntend=active_nodes.end();nt!=ntend;++nt) {
            if (is_strictly_greater(minbnd,nt->_VECTptr->back(),ctx)) {
                minbnd=nt->_VECTptr->back();
                sel=nt;
            }
        }
        if (sel==active_nodes.end()) {
            err("node selection failed");
            break;
        }
        const vecteur &s_x=*sel->_VECTptr->at(2)._VECTptr;
        s_fval=sel->_VECTptr->back();
        // branching variable selection (default strategy: pseudocost)
        assert(branch_var(parm,s_x,s_fval,3,bv));
        s_xf=s_x[bv];
        s_xl=_floor(s_xf,ctx);
        s_xu=_ceil(s_xf,ctx);
        lb=*sel->_VECTptr->at(0)._VECTptr;
        ub=*sel->_VECTptr->at(1)._VECTptr;
        oldb=ub[bv];
        active_nodes.erase(sel); // remove parent subproblem
        // branch down
        ub[bv]=s_xl;
        if (is_greater(ub[bv],lb[bv],ctx) && bb_solve_subproblem(parm,active_nodes,fval,sub_res,incumbent))
            update_pseudocost(bv,pc_lo,pc_lo_count,fval-s_fval,s_xf-s_xl);
        // branch up
        ub[bv]=oldb;
        lb[bv]=s_xu;
        if (is_greater(ub[bv],lb[bv],ctx) && bb_solve_subproblem(parm,active_nodes,fval,sub_res,incumbent))
            update_pseudocost(bv,pc_hi,pc_hi_count,fval-s_fval,s_xu-s_xf);
    }
    return true;
}
/* Outer Approximation algorithm implementation (Duran & Grossman, Fletcher & Leyffer).
 * CONT_SOL is the solution if the continuous relaxation.
 * INCUMBENT status must be set to PENDING prior the call. */
bool nlp_problem::outer_approximation(const meth_parm &parm,const vecteur &cont_sol,optimum_t &incumbent) {
    gen UB=plus_inf,LB=minus_inf,LB_prev,gap=plus_inf;
    vecteur moa_sol,constr(eq_cons),ipt_opt(cont_sol),s,cobyla_cons,ccons,cobyla_vars;
    gen alpha=temp_symb("alpha",-1,ctx),beta=temp_symb("beta",-1,ctx),fval,ipt_optval,c;
    int res,i,n=vars.size(),status,ones;
    ipt_solver *ipt=NULL;
    optima_t ipt_res;
    const_iterateur it,itend,jt;
    bool append_f_lin=true,isbin=is_intvars_all_binary(),ret=true;
    if (!compute_obj_val(cont_sol,ipt_optval))
        return false;
    for (i=0;i<n;++i) {
        if (!is_inf(lb[i]))
            constr.push_back(symb_superieur_egal(vars[i],lb[i]));
        if (!is_inf(ub[i]))
            constr.push_back(symb_inferieur_egal(vars[i],ub[i]));
    }
    if (isbin)
        debug("All integer variables are binary, replacing infeasibility problems with integer cuts");
    // make cobyla constraints
    cobyla_cons=mergevecteur(mergevecteur(ineq_cons,eq_cons),multvecteur(-1,eq_cons));
    cobyla_cons=mergevecteur(cobyla_cons,subvecteur(vars,lb));
    cobyla_cons=mergevecteur(cobyla_cons,subvecteur(ub,vars));
    // iterate OA
    _spc=0;
    initp.clear();
    while(true) {
        if (interrupted || ctrl_c) {
            interrupted=ctrl_c=false;
            throw std::runtime_error(gettext("Stopped by user interruption"));
        }
        // solve master problem
        if (!ipt_opt.empty()) {
            s=subvecteur(vars,ipt_opt);
            if (append_f_lin) {
                c=ipt_optval;
                c+=scalarproduct(subst(obj_gradient,vars,ipt_opt,false,ctx),s,ctx);
                constr.push_back(symb_inferieur_egal(c,alpha));
            }
            for (it=ineq_cons.begin(),itend=ineq_cons.end(),jt=ineq_jacobian.begin();it!=itend;++it,++jt) {
                c=subst(*it,vars,ipt_opt,false,ctx);
                c+=scalarproduct(subst(*jt->_VECTptr,vars,ipt_opt,false,ctx),s,ctx);
                constr.push_back(symb_superieur_egal(c,0));
            }
        }
        vars.push_back(alpha);
        LB_prev=LB;
        res=minimize_linear(alpha,constr,moa_sol,LB);
        vars.pop_back();
        moa_sol.pop_back();
        ++_spc;
        _ss << "Lower bound: " << LB; debug();
        gap=(UB-LB)/(1+_abs(LB,ctx));
        if (res!=_NLP_OPTIMAL || is_strictly_greater(parm.eps,gap,ctx))
            break;
        // fix integer variables and solve the reduced NLP
        for (i=0;i<n;++i) {
            if (is_intvar(i))
                _fixed_vars.insert(make_pair(i,moa_sol[i]));
        }
        ipt=new ipt_solver(*this,parm);
        ipt_res.clear();
        _iter_count=0;
        if (!ipt->optimize(ipt_res)) {
            ret=false;
            break;
        }
        if (ipt_res.empty())
            break;
        ipt_optval=undef;
        status=_NLP_PENDING;
        find_best_solution(parm,ipt_res,ipt_opt,ipt_optval,status);
        if ((status==_NLP_OPTIMAL || status==_NLP_FAILED) && compute_obj_val(ipt_opt,ipt_optval)) {
            if (is_strictly_greater(UB,ipt_optval,ctx)) {
                UB=ipt_optval;
                gap=(UB-LB)/(1+_abs(LB,ctx)); // UB has changed, update the gap value
                incumbent.f=ipt_optval;
                incumbent.x=ipt_opt;
                incumbent.res=status;
                _ss << gettext("New incumbent solution found") << ": " << UB << ", "
                    << gettext("gap") << ": " << 100.0*gap << "%";
                msg();
            }
            append_f_lin=true;
        } else { // reduced NLP is infeasible
            if (isbin) { // append integer cut to MILP constraints
                debug("NLP is infeasible, adding integer cut...");
                ones=0;
                c=0;
                for (it=moa_sol.begin(),itend=moa_sol.end(),i=0;it!=itend;++it,++i) {
                    if (!is_intvar(i)) continue;
                    if (is_one(_round(*it,ctx))) {
                        ++ones;
                        c+=vars[i];
                    } else c-=vars[i];
                }
                constr.push_back(symb_inferieur_egal(c,ones-1));
                ipt_opt.clear();
            } else {
                // solve the corresponding NLPF using COBYLA
                debug("NLP is infeasible, solving feasibility problem...");
                if (_fixed_vars.size()==vars.size()) {
                    ipt_opt=moa_sol;
                } else {
                    ccons=subs_fixed_vars(cobyla_cons);
                    for (i=ccons.size();i-->0;) {
                        if (is_positive_safe(ccons[i],false,100,ctx))
                            ccons.erase(ccons.begin()+i);
                    }
                    addvecteur(ccons,vecteur(ccons.size(),beta),ccons);
                    cobyla_vars=vars;
                    remove_elements_with_indices(cobyla_vars,_fixed_vars);
                    cobyla_vars.push_back(beta);
                    vecteur fcv=makevecteur(beta,ccons,cobyla_vars),x0=moa_sol;
                    remove_elements_with_indices(x0,_fixed_vars);
                    x0.push_back(0);
                    cobyla_gc gc={fcv,ctx};
                    int maxiter=parm.max_iter,cres,iter,tries,max_tries=5;
                    try {
                        tries=0;
                        while (maxiter>0 && tries++<max_tries) {
                            iter=maxiter;
                            cres=giac_cobyla(&gc,x0,iter,parm.eps,0);
                            maxiter-=iter;
                            if (cres!=COBYLA_ROUNDING)
                                break;
                        }
                        if (cres<0) { // error
                            ret=false;
                            break;
                        }
                        x0.pop_back();
                        insert_fixed_vars(x0);
                        ipt_opt=x0;
                    } catch (const std::runtime_error &e) { // an error raised in COBYLA or user interruption
                        if (strstr(e.what(),"user interruption")!=NULL || interrupted || ctrl_c) {
                            throw e;
                        } else {
                            msg(e.what());
                            ret=false;
                        }
                        break;
                    }
                }
                append_f_lin=false;
            }
        }
        _fixed_vars.clear();
        delete ipt;
        ipt=NULL;
    }
    if (ipt!=NULL) delete ipt;
    _fixed_vars.clear();
    return ret;
}
/* Global minimization using METHOD with parameters PARM.
 * IP must be a list of initial point(s) for the specified method (possibly empty)
 * or [interval,L,U] where L and U are lists of lower and upper bounds, respectively,
 * which is used as a starting (hyper)rectangle.
 * The location of the global minimum is stored to OPT before returning. */
int nlp_problem::optimize(int method,const meth_parm &parm,const vecteur &ip,vecteur &opt,gen &optval) {
    int ret=_NLP_PENDING,i;
    if (var_count()==0) { // no variables, objective is constant
        opt.clear();
        optval=obj;
        return _NLP_OPTIMAL;
    }
    // there is at least one independent variable
    _msg_level=parm.msg_level;
    optval=undef;
    if (method==0) {
        err("no method specified");
        return _NLP_ERROR;
    }
    if (!_initialized) {
        _has_initial_rect=false;
        reset_timer();
        if (!ip.empty() && ip.front()==at_interval) {
            assert(ip.size()==3);
            bd_rect=mtran(makevecteur(ip[1],ip[2]));
            initp.clear();
            _has_initial_rect=true;
        } else initp=ip;
        initialize(method,parm);
        _ss << "-- Initialization time: " << elapsed_secs()*1000.0 << " ms"; debug();
        if (_infeas)
            return _NLP_INFEAS;
        if (var_count()>0) {
            _ss << gettext("Tolerance") << ": " << gettext("optimality") << " = " << parm.eps
                << ", " << gettext("feasibility") << " = " << parm.tol; msg();
        }
    } else if (!ip.empty())
        initp=mergevecteur(initp,ip);
    if (var_count()>0) {
        int max_tries=5,tries=0;
        bool direct_meth=method>_NLP_INTERIOR_POINT;
        bool is_unconstr=is_unconstrained();
        vecteur cons;
        initp_iter=initp.begin();   // initialize the initial point iterator
        _penalty_scheme=0;          // no inclusion of constraints into the objective
        _iter_count=0;              // restart iteration counter
        optval=undef;               // no solution (yet)
        if (method==_NLP_AUTOMATIC) { // guess the method
            if (var_count()==1) {
                /* Univariate minimizer */
                if (eq_count()+ineq_count()>0 || is_inf(lb[0]) || is_inf(ub[0]))
                    ret=optimize(is_intvar(0)?_NLP_DIFFERENTIAL_EVOLUTION:_NLP_NELDER_MEAD,parm,vecteur(0),opt,optval);
                else { // in-segment optimization
                    msg("Applying univariate minimization...");
                    if (is_strictly_greater(lb.front(),ub.front(),ctx))
                        return _NLP_INFEAS;
                    opt.resize(1);
                    if (is_intvar(0)) { // integer variable
                        int min_ival=_ceil(lb.front(),ctx).val,max_ival=_floor(ub.front(),ctx).val;
                        gen fval;
                        for (i=min_ival;i<=max_ival;++i) {
                            if (!compute_obj_val(vecteur(1,i),fval))
                                continue;
                            if (is_undef(optval) || is_strictly_greater(optval,fval,ctx)) {
                                optval=fval;
                                opt.front()=i;
                            }
                        }
                        if (is_undef(optval)) {
                            opt.clear();
                            return _NLP_INFEAS;
                        }
                        return _NLP_OPTIMAL;
                    }
    #ifdef HAVE_LIBGSL
                    // use GSL univariate minimizer (Brent minimization algorithm)
                    debug("Using GSL minimizer");
                    gsl_min_fminimizer *brent=gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
                    gsl_function F;
                    F.function=&gsl_my_uf;
                    F.params=static_cast<void*>(this);
                    vecteur x0;
                    double left=_evalf(lb.front(),ctx).to_double(ctx),right=_evalf(ub.front(),ctx).to_double(ctx);
                    double fval,lfv,rfv,x0d,prev_fval,best_fval=DBL_MAX,best_x0d;
                    int status;
                    bool improv_fail;
                    while (_iter_count<parm.max_iter && restart(x0)) {
                        x0d=_evalf(x0.front(),ctx).to_double(ctx);
                        gsl_min_fminimizer_set(brent,&F,x0d,left,right);
                        improv_fail=false;
                        prev_fval=DBL_MAX;
                        do {
                            status=gsl_min_fminimizer_iterate(brent);
                            ++_iter_count;
                            x0d=gsl_min_fminimizer_x_minimum(brent);
                            left=gsl_min_fminimizer_x_lower(brent);
                            right=gsl_min_fminimizer_x_upper(brent);
                            fval=gsl_min_fminimizer_f_minimum(brent);
                            lfv=gsl_min_fminimizer_f_lower(brent);
                            rfv=gsl_min_fminimizer_f_upper(brent);
                            if (lfv<fval) {
                                x0d=left;
                                fval=lfv;
                            } else if (rfv<fval) {
                                x0d=right;
                                fval=rfv;
                            }
                            if (status==GSL_FAILURE) {
                                _ss << gettext("failed to improve the current interval") << " [" << left << "," << right << "]";
                                warn();
                                improv_fail=true; // update incumbent if necessary and terminate
                            }
                            if (status==GSL_EBADFUNC) {
                                warn("objective function has singular point(s), the problem may be unbounded");
                                break; // start over
                            }
                            if (!improv_fail)
                                status=gsl_min_test_interval(left,right,parm.tol,0.0);
                            if (improv_fail || status==GSL_SUCCESS) {
                                if (!improv_fail) {
                                    save_point(vecteur(1,x0d));
                                    if (prev_fval-fval>parm.eps)
                                        status=GSL_CONTINUE;
                                }
                                if (fval<best_fval) {
                                    best_fval=fval;
                                    best_x0d=x0d;
                                }
                            }
                            prev_fval=fval;
                        } while (!improv_fail && status==GSL_CONTINUE && _iter_count<parm.max_iter);
                    }
                    gsl_min_fminimizer_free(brent);
                    if (best_fval<DBL_MAX) {
                        optval=best_fval;
                        opt.front()=best_x0d;
                    } else opt.front()=x0d; // last point
                    ret=is_undef(optval)?_NLP_FAILED:_NLP_OPTIMAL;
    #else
                    // no GSL, use global segment optimization with Brent's method (initial points are ignored!)
                    try {
                        opt[0]=find_minimum(obj,lb[0],ub[0],true,parm.eps,parm.max_iter,ctx);
                    } catch (const std::runtime_error &e) {
                        _ss << "Univariate minimizer error: " << e.what();
                        debug();
                        err("Univariate minimizer error");
                        return _NLP_ERROR;
                    }
                    if (is_undef(opt[0]) || !compute_obj_val(opt,optval)) {
                        err("Brent's method failed");
                        return _NLP_ERROR;
                    }
                    ret=_NLP_OPTIMAL;
    #endif
                }
            } else if (_linconstr && _linobj) { // use linear programming solver
                msg("Applying the method of linear programming...");
                reset_timer();
                ret=minimize_linear(obj,constraints(true),opt,optval);
                if (ret==_NLP_ERROR)
                    err("lpsolve failed");
            } else if (_smooth_obj && _smooth_constr) {
                if (eq_count()==0 && ineq_count()==0 && _compact && parm.cluster<0 && initp.empty() && !has_intvars())
                    ret=optimize(_NLP_DIFFERENTIAL_EVOLUTION,parm,vecteur(0),opt,optval);
                else ret=optimize(_NLP_INTERIOR_POINT,parm,vecteur(0),opt,optval);
            } else if (has_intvars() || (_compact && eq_count()==0 && ineq_count()==0 && initp.empty())) {
                // use differential evolution for MINLP or box-constrained NLP
                ret=optimize(_NLP_DIFFERENTIAL_EVOLUTION,parm,vecteur(0),opt,optval);
            } else {
                // try Nelder-Mead, if it fails use differential evolution
                ret=optimize(_NLP_NELDER_MEAD,parm,vecteur(0),opt,optval);
                if (ret==_NLP_FAILED || ret==_NLP_ERROR) {
                    warn("Nelder-Mead algorithm failed, switching to differential evolution");
                    if (_iter_data) iteration_data.clear();
                    ret=optimize(_NLP_DIFFERENTIAL_EVOLUTION,parm,vecteur(0),opt,optval);
                }
            }
            return ret;
        }
        // solve by the specified method
        reset_timer();
        gen val;
        optima_t res,res_improv;
        optimum_t o;
        res.reserve(initp.size());
        optima_t::iterator opt_iter;
        vecteur failed_initp;
        failed_initp.reserve(initp.size());
        vector<int> failed_ind;
        failed_ind.reserve(initp.size());
        bool has_failed;
        // optimize by using the specified method
        switch (method) {
        case _NLP_AUTOMATIC:
            assert(false); // this case is handled above!
        case _NLP_NELDER_MEAD:
            msg("Applying the method of Nelder & Mead...");
            if (has_intvars()) {
                err("Nelder-Mead method does not support integrality constraints");
                return _NLP_ERROR;
            }
            _penalty_scheme=1;
            if (!nelder_mead(parm,res))
                return _NLP_ERROR;
            break;
        case _NLP_DIFFERENTIAL_EVOLUTION:
            msg("Applying the method of differential evolution...");
            _penalty_scheme=2;
            if (!differential_evolution(parm,o))
                return _NLP_ERROR;
            res.push_back(o);
            break;
        case _NLP_COBYLA:
            msg("Applying COBYLA algorithm...");
            if (!cobyla(constraints(true,true),parm,res))
                return _NLP_ERROR;
            break;
        case _NLP_INTERIOR_POINT:
            if (parm.cluster>=0 && _have_hessian) {
                vecteur ip_tmp;
                int n=var_count(),nip=parm.search_points;
                if (nip<=0) nip=_round(1000*std::sqrt(n),ctx).val;
                if (!initp.empty())
                    warn("ignoring provided initial points");
                initp.resize(nip);
                for (iterateur ip_it=initp.begin();ip_it!=initp.end();++ip_it) {
                    make_random_initial_point(ip_tmp);
                    *ip_it=ip_tmp;
                }
    #ifdef HAVE_LIBGSL
                gsl_matrix *gslm=gsl_matrix_alloc(n,n);
                size_t j,k;
    #endif
                for (i=nip;i-->0;) {
                    matrice H=subst(obj_hessian,vars,*initp[i]._VECTptr,false,ctx);
    #ifdef HAVE_LIGBSL
                    j=0;
                    for (const_iterateur jt=H.begin();jt!=H.end();++jt,++j) {
                        k=0;
                        for (const_iterateur kt=jt->_VECTptr->begin();kt!=jt->_VECTptr->end();++kt,++k) {
                            gsl_matrix_set(gslm,j,k,kt->to_double(ctx));
                        }
                    }
                    if (gsl_linalg_cholesky_decomp1(gslm)==GSL_EDOM)
                        initp.erase(initp.begin()+i);
                }
                gsl_matrix_free(gslm);
    #else
                    if (!is_positive_definite_inplace(H,true,ctx))
                        initp.erase(initp.begin()+i);
                }
    #endif
                if (!initp.empty()) {
                    _ss << "Clustering " << (int)initp.size() << " random PDH points...";
                    msg();
                    if ((parm.cluster>0 && (int)initp.size()>parm.cluster) ||
                            (parm.cluster==0 && (int)initp.size()>10*n)) { // cluster
                        vector<int> cluster_indices;
                        int nc=parm.cluster,max_nc=initp.size();
                        vecteur cc;
                        if (kmeans(initp,nc,2,max_nc,100,cluster_indices,&cc,NULL,ctx)!=0)
                            warn("clustering failed");
                        else {
                            cclust_centers.resize(cc.size());
                            cclust_sizes.resize(cc.size());
                            const_iterateur cct=cc.begin(),cctend=cc.end();
                            for (i=0;cct!=cctend;++cct,++i) {
                                cclust_centers[i]=cct->_VECTptr->front();
                                cclust_sizes[i]=cct->_VECTptr->back();
                            }
                        }
                    }
                    _ss << "Found " << (int)cclust_centers.size() << " PDH clusters";
                    msg();
                } else warn("objective function appears to be concave");
                initp=cclust_centers;
                initp_iter=initp.begin();
            }
    #ifdef HAVE_LIBGSL
            if (is_unconstr) {
                msg("Applying BFGS method...");
                if (!gsl_bfgs(parm,res))
                    return _NLP_ERROR;
            } else {
    #endif
                msg("Applying interior-point method...");
                if (!_have_hessian || !_smooth_constr) {
                    err("Could not compute objective hessian and/or constraint gradients");
                    return _NLP_ERROR;
                }
                { ipt_solver *ipt=new ipt_solver(*this,parm);
                if (!ipt->optimize(res))
                    ret=_NLP_ERROR;
                delete ipt; }
    #ifdef HAVE_LIBGSL
            }
    #endif
            break;
        default:
            err("method not implemented");
            return _NLP_ERROR;
        }
        _penalty_scheme=0;
        res_improv.clear();
        failed_initp.clear();
        failed_ind.clear();
        for (opt_iter=res.begin();opt_iter!=res.end();++opt_iter) {
            if (opt_iter->res==_NLP_UNBOUNDED)
                return opt_iter->res;
            if (opt_iter->res==_NLP_ERROR || opt_iter->res==_NLP_PENDING)
                continue;
            if (method==_NLP_NELDER_MEAD || method==_NLP_DIFFERENTIAL_EVOLUTION) {
                failed_initp.push_back(opt_iter->x);
                failed_ind.push_back(opt_iter-res.begin());
            }
        }
        if (parm.postprocess && !failed_initp.empty()) {
            // attempt to improve solutions
            initp=failed_initp;
            initp_iter=initp.begin();
            if (method==_NLP_DIFFERENTIAL_EVOLUTION && has_intvars()) {
                // append integrality constraints
                for (i=0;i<var_count();++i) {
                    if (is_intvar(i)) {
                        ineq_cons.push_back(_neg(pow(sin(cst_pi*vars[i],ctx),2),ctx));
                        _linconstr=false;
                        if (_smooth_constr)
                            ineq_jacobian.push_back(_grad(makesequence(ineq_cons.back(),vars),ctx));
                    }
                }
            }
            if (_have_hessian && _smooth_constr && !has_intvars()) { // use interior-point method
                debug("Using interior-point method to improve solutions...");
                ipt_solver *ipt=new ipt_solver(*this,parm);
                _iter_count=0;
                if (!ipt->optimize(res_improv)) {
                    debug("Interior-point algorithm failed");
                    res_improv.clear();
                }
                delete ipt;
            } else { // use COBYLA
                debug("Using COBYLA to improve solutions...");
                if (!cobyla(constraints(true,true),parm,res_improv)) {
                    debug("COBYLA failed");
                    res_improv.clear();
                }
            }
            vector<int>::const_iterator it=failed_ind.begin();
            int improv_count=0;
            for (opt_iter=res_improv.begin();opt_iter!=res_improv.end();++opt_iter,++it) {
                if (opt_iter->res==_NLP_OPTIMAL) {
                    res[*it]=*opt_iter;
                    ++improv_count;
                }
            }
            if (improv_count>0) { _ss << gettext("Improved") << " " << improv_count << " "
                                    << gettext("solution(s)"); msg(); }
        }
        find_best_solution(parm,res,opt,optval,ret);
        switch (ret) {
        case _NLP_INFEAS:
        case _NLP_ERROR:
            return ret;
        case _NLP_PENDING:
            return _NLP_ERROR;
        default:
            break;
        }
        // print solver summary
        _ss << gettext("Total number of restarts:") << " " << res.size() << " ("
            << _sol_stat[0] << " " << gettext("optimal") << ", "
            << _sol_stat[1] << " " << gettext("failed") << ", "
            << _sol_stat[2] << " " << gettext("infeasible") << ", "
            << _sol_stat[4] << " " << gettext("error(s)") << ", "
            << _sol_stat[5] << " " << gettext("pending") << ")";
        msg();
        if (method==_NLP_INTERIOR_POINT && has_intvars() && !is_intsol(parm,opt)) {
            // do branch & bound for MINLP
            initp.clear();
            optimum_t intsol;
            if (is_convex() && _compact) {
                msg("Applying outer approximation method to find integer solutions...");
                if (!outer_approximation(parm,opt,intsol) || intsol.res==_NLP_ERROR)
                    return _NLP_ERROR;
                _ss << gettext("Finished with") << " " << _spc << " " << gettext("major iterations"); msg();
            } else {
                msg("Applying branch & bound algorithm to find integer solutions...");
                if (!branch_and_bound(parm,opt,intsol) || intsol.res==_NLP_ERROR)
                    return _NLP_ERROR;
                _ss << gettext("Solved") << " " << _spc << " " << gettext("subproblems"); msg();
            }
            _fixed_vars.clear();
            optval=undef;
            opt.clear();
            ret=intsol.res;
            if (ret==_NLP_PENDING || ret==_NLP_INFEAS)
                return _NLP_INFEAS;
            opt=intsol.x;
            optval=intsol.f;
            if (is_undef(optval) && !compute_obj_val(opt,optval))
                ret=_NLP_ERROR;
        }
        _ss << gettext("Optimization time") << ": " << elapsed_secs()*1000.0 << " ms"; msg();
        if (_fevalc>0) { _ss << gettext("Objective function evaluations") << ": " << _fevalc; msg(); }
    } else {
        msg("Solution found by preprocessor");
        opt.clear();
        optval=obj;
        ret=_NLP_OPTIMAL;
    }
    // postsolve
    if (_presolved) {
        gen pt=_iter_data?makevecteur(opt,iteration_data):opt;
        debug("Post-processing...");
        postprocess(pt);
        opt=_iter_data?*pt._VECTptr->front()._VECTptr:*pt._VECTptr;
        if (_iter_data)
            iteration_data=*pt._VECTptr->back()._VECTptr;
    }
    return ret;
}
/* Filter out the best solution. */
void nlp_problem::find_best_solution(const meth_parm &parm,optima_t &res,vecteur &opt,gen &optval,int &ret) {
    std::fill(_sol_stat.begin(),_sol_stat.end(),0);
    optima_t::iterator opt_iter=res.begin(),opt_itend=res.end();
    for (;opt_iter!=opt_itend;++opt_iter) {
        ++_sol_stat[opt_iter->res];
        switch (opt_iter->res) {
        case _NLP_PENDING:
        case _NLP_ERROR:
            continue;
        default:
            break;
        }
        if (is_undef(opt_iter->f) && !compute_obj_val(opt_iter->x,opt_iter->f)) {
            opt_iter->res=_NLP_ERROR;
            continue;
        }
        if (opt_iter->res==_NLP_INFEAS || !is_feasible(opt_iter->x,parm.tol,true)) {
            opt_iter->res=_NLP_INFEAS;
            if (ret==_NLP_PENDING)
                ret=opt_iter->res;
            continue;
        }
        if (!is_undef(optval)) {
            if (is_greater(parm.eps,_abs(optval-opt_iter->f,ctx),ctx)) {
                if (opt_iter->res==_NLP_FAILED && ret==_NLP_OPTIMAL)
                    continue;
            } else if (is_greater(opt_iter->f,optval,ctx))
                continue;
        }
        opt=opt_iter->x;
        optval=opt_iter->f;
        ret=opt_iter->res;
    }
}
/* Rudimentary AMPL modeling language parser, which is able to load problems from minlplib.org.
 * Variables are written to V, objective(s) to O and constraints to C.
 * Returns true iff successful, otherwise it returns false with message MSG. */
bool nlp_problem::ampl_load(const std::string &filename,vecteur &v,vecteur &o,vecteur &c,string &msg,GIAC_CONTEXT) {
    ifstream ifs(filename.c_str());
    if (!ifs.is_open()) {
        msg=gettext("AMPL parser: failed to open file for reading");
        return false; // failed to open file
    }
    string line,token;
    size_t pos,p;
    int status=0,st,block=0; // 1 - variable, 2 - objective(s), 3 - constraints
    bool ismin=true;
    string e;
    vecteur x;
    gen g;
    while (getline(ifs,line)) {
        pos=line.find('#');
        if (pos!=string::npos)
            line=line.substr(0,pos);
        pos=line.find_first_not_of(" \t");
        if (pos==string::npos)
            line="";
        else line=line.substr(pos);
        pos=line.find_last_not_of(" \t");
        if (pos!=string::npos)
            line=line.substr(0,pos+1);
        if (line.find("var ")==0) { block=1; line=line.substr(4); }
        else if (line.find("minimize ")==0) { block=2; line=line.substr(9); }
        else if (line.find("maximize ")==0) { block=2; line=line.substr(9); ismin=false; }
        else if (line.find("subject to ")==0) { block=3; line=line.substr(11); }
        else {
            if (line=="var") { block=1; continue; }
            if (line=="minimize") { block=2; continue; }
            if (line=="maximize") { block=2; ismin=false; continue; }
            if (line=="subject to") { block=3; continue; }
        }
        if (block==0)
            continue;
        // tokenize line
        do {
            pos=line.find_first_of(" \t");
            token=line.substr(0,pos);
            if (pos!=string::npos)
                line.erase(0,pos+1);
            if (token.empty())
                continue;
            switch (block) {
            case 1:
                switch (status) {
                case 0: // initialize variable
                    if (token.at(token.size()-1)==';') {
                        token.erase(token.size()-1,1);
                        v.push_back(makevecteur(identificateur(token),0,minus_inf,plus_inf,undef));
                        status=0;
                    } else {
                        x=makevecteur(identificateur(token),0,minus_inf,plus_inf,undef);
                        status=1;
                    }
                    break;
                case 1:
                    if (token=="integer")
                        x[1]=1;
                    else if (token=="binary")
                        x[1]=2;
                    else if (token==">=") status=2;
                    else if (token=="<=") status=3;
                    else if (token==":=") status=4;
                    else {
                        msg=gettext("AMPL parser: invalid variable specification");
                        ifs.close();
                        return false;
                    }
                    break;
                case 2: // read lower bound
                case 3: // read upper bound
                case 4: // read initial value
                    st=status;
                    if (token.at(token.size()-1)==',') {
                        token.erase(token.size()-1,1);
                        status=1;
                    } else if (token.at(token.size()-1)==';') {
                        token.erase(token.size()-1,1);
                        status=0;
                    } else { ifs.close(); return false; }
                    try {
                        x[st]=_expr(string2gen(token,false),contextptr);
                    } catch (const std::runtime_error &e) {
                        msg=string(gettext("AMPL parser"))+string(": ")+e.what();
                        ifs.close();
                        return false;
                    }
                    if (status==0)
                        v.push_back(x);
                    break;
                default:
                    assert(false);
                }
                break;
            case 2: // read objective(s)
            case 3: // read constraints
                switch (status) {
                case 0:
                    if (token.at(token.size()-1)==':')
                        status=1;
                    break;
                case 1:
                    if (token.at(token.size()-1)==';') {
                        token.erase(token.size()-1,1);
                        e+=token;
                        try {
                            g=_expr(string2gen(e,false),contextptr);
                        } catch (const std::runtime_error &e) {
                            msg=string(gettext("AMPL parser"))+string(": ")+e.what();
                            ifs.close();
                            return false;
                        }
                        if (block==2) o.push_back(g);
                        else c.push_back(g);
                        status=0;
                        e.clear();
                    } else {
                        e+=token+" ";
                    }
                    break;
                default:
                    assert(false);
                }
                break;
            }
        } while (pos!=string::npos);
    }
    ifs.close();
    return true;
}
/* initialize method parameters with default values */
void nlp_problem::meth_parm::init_defaults(double eps) {
    cross_prob=0.5;
    scale=0;            // the optimal scale will be computed (DE)
    reflect_ratio=1.0;
    expand_ratio=2.0;
    contract_ratio=shrink_ratio=0.5;
    search_points=0;    // use 10*n points, where n is the number of variables
    penalty=undef;      // penalty factor for constraint violation in direct methods (undef means automatic)
    max_iter=500;       // maximum number of iterations allowed
    tol=1e-5;           // feasibility tolerance
    this->eps=eps;      // optimality tolerance
    size=0;             // initial simplex size for NM (automatic)
    step=0;             // initial step-size for BGFS or adaptive penalty step for NM (automatic)
    postprocess=true;   // attempt to improve failed direct-method solutions using COBYLA
    presolve=true;      // the problem will be presolved by default
    border=true;        // attempt to estimate bounds on variables
    convex=-1;          // check whether the problem is convex automatically
    cluster=-1;         // no clustering (0 - automatic number of clusters, integer k - k clusters)
    msg_level=1;        // 0 - no output
                        // 1 - show only errors (default)
                        // 2 - show warnings and errors
                        // 3 - show all messages
                        // 4 - debug mode
}
/* END OF NLP_PROBLEM CLASS */

/* Nonlinear programming problem solver.
 * Usage: nlpsolve(Expr(obj)||Func(obj),Lst(constr),[Seq(bounds)],[opts]) */
gen _nlpsolve(const gen &g,GIAC_CONTEXT) {
    gen fnameg(undef);
    vecteur constr,vars,initp,lbv,ubv,intvars,binvars;
    gen obj;
    if (g.type==_STRNG) {
        if (g.subtype==-1)
            return g;
        fnameg=g;
    }
    if (g.type==_VECT && g.subtype==_SEQ__VECT && !g._VECTptr->empty() && g._VECTptr->front().type==_STRNG) {
        fnameg=g._VECTptr->front();
    }
    if (!is_undef(fnameg)) {
        string fname=gen2string(fnameg),ampl_msg;
        vecteur rest_args=g.type==_VECT?vecteur(g._VECTptr->begin()+1,g._VECTptr->end()):vecteur(0),pv,po,pc;
        if (!nlp_problem::ampl_load(fname,pv,po,pc,ampl_msg,contextptr))
            return generr(ampl_msg.c_str(),false);
        if (po.empty())
            obj=0;
        else if (po.size()>1)
            return generr("Multi-objective problems are not supported");
        else obj=po.front();
        *logptr(contextptr) << gettext("Loaded a problem with") << " " << pv.size() << " " << gettext("variables and")
                            << " " << pc.size() << " " << gettext("constraints") << "\n";
        lbv.resize(pv.size(),minus_inf);
        ubv.resize(pv.size(),plus_inf);
        initp.resize(pv.size(),undef);
        int i=0;
        const_iterateur it,itend;
        for (it=pv.begin(),itend=pv.end();it!=itend;++it,++i) {
            assert(it->type==_VECT && it->_VECTptr->size()==5 && it->_VECTptr->at(1).is_integer());
            const vecteur &vspec=*it->_VECTptr;
            vars.push_back(vspec[0]);
            switch (vspec[1].val) {
                case 0: break;
                case 1: intvars.push_back(vspec[0]); break;
                case 2: binvars.push_back(vspec[0]); break;
                default: assert(false);
            }
            lbv[i]=vspec[2];
            ubv[i]=vspec[3];
            initp[i]=vspec[4];
        }
        if (has_inf_or_undef(initp))
            initp.clear();
        else initp=vecteur(1,initp);
        vecteur arguments=makevecteur(obj,pc);
        arguments=mergevecteur(arguments,*_box_constraints(makesequence(vars,mtran(makevecteur(lbv,ubv))),contextptr)._VECTptr);
        if (!initp.empty())
            arguments.push_back(symb_equal(change_subtype(_NLP_INITIALPOINT,_INT_MAPLECONVERSION),_zip(makesequence(at_equal,vars,initp.front()),contextptr)));
        if (!intvars.empty())
            arguments.push_back(symb_equal(change_subtype(_NLP_INTEGERVARIABLES,_INT_MAPLECONVERSION),intvars));
        if (!binvars.empty())
            arguments.push_back(symb_equal(change_subtype(_NLP_BINARYVARIABLES,_INT_MAPLECONVERSION),binvars));
        arguments=mergevecteur(arguments,rest_args);
        return _nlpsolve(change_subtype(arguments,_SEQ__VECT),contextptr);
    }
    if (g.type==_SYMB) // unconstrained minimization without any options and bounds
        return _nlpsolve(makesequence(g,vecteur(0)),contextptr);
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()<2)
        return generrtype("Expected a sequence of two or more arguments");
    const vecteur &gv=*g._VECTptr;
    obj=gv.front();
    add_identifiers(obj,vars,contextptr);
    bool maximize=false,minlp=false,nonneg=false;
    int num_mask=num_mask_withint | num_mask_withfrac;
    int machine_accuracy=-_round(std::log10(epsilon(contextptr)),contextptr).val;
    int wc=0,method=_NLP_AUTOMATIC;
    double eps=std::pow(10,-(machine_accuracy*2)/3);
    nlp_problem::meth_parm parm;
    parm.init_defaults(eps);
    const_iterateur it=gv.begin();
    /* get the constraints, if any */
    if (gv[1].type==_VECT) {
        constr=*gv[1]._VECTptr;
        // parse constraints, detect infeasible/redundant ones, convert expr=a..b to expr>=a and expr<=b
        for (int i=constr.size();i-->0;) {
            gen &c=constr[i];
            if (c.is_integer() && c.subtype==_INT_BOOLEAN) {
                if (c.val==0)
                    return vecteur(0); // infeasible
                constr.erase(constr.begin()+i);
                continue;
            }
            if (!is_equal(c) && which_ineq(c)==0)
                return gensizeerr(gettext("Invalid constraint ")+print_INT_(1+i));
            if (is_equal(c) && c._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_interval)) {
                gen a=c._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->front();
                gen b=c._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->back();
                gen lh=c._SYMBptr->feuille._VECTptr->front();
                c=symb_inferieur_egal(lh,b);
                constr.insert(constr.begin()+i,symb_superieur_egal(lh,a));
                continue;
            }
        }
        add_identifiers(constr,vars,contextptr);
        ++it;
    }
    if (vars.empty()) {
        print_error("no variables found, returning the objective value",contextptr);
        return change_subtype(makevecteur(obj,vecteur(0)),_LIST__VECT);
    }
    lbv.resize(vars.size(),minus_inf);
    ubv.resize(vars.size(),plus_inf);
    /* parse options and bounds */
    while (++it!=gv.end()) {
        bool invalid=false;
        if (*it==at_maximize || is_mcint(*it,_NLP_MAXIMIZE)) {
            maximize=true;
        } else if (*it==at_cluster) {
            parm.cluster=0;
        } else if (*it==at_convex) {
            parm.convex=1;
        } else if (is_equal(*it)) {
            const gen &lh=it->_SYMBptr->feuille._VECTptr->front();
            const gen &rh=it->_SYMBptr->feuille._VECTptr->back();
            if (lh==at_assume && rh.is_integer() && (rh.val==_ZINT || rh.subtype==_INT_MAPLECONVERSION)) {
                // parse assumptions
                switch (rh.val) {
                case _NONNEGINT: case _NLP_NONNEGINT:
                    intvars=vars;
                case _NLP_NONNEGATIVE:
                    nonneg=true;
                    break;
                case _NLP_INTEGER: case _ZINT:
                    intvars=vars;
                    break;
                case _NLP_BINARY:
                    binvars=vars;
                    break;
                default:
                    invalid=true;
                    break;
                }
            } else if (lh==at_maximize && rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                maximize=(bool)rh.val;
            } else if (lh==at_border && rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                parm.border=(bool)rh.val;
            } else if (lh==at_convex && rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                parm.convex=rh.val;
            } else if (lh==at_cluster && rh.is_integer()) {
                if (rh.subtype==_INT_BOOLEAN)
                    parm.cluster=rh.val-1;
                else if (rh.val<1)
                    return generr("Expected a positive integer");
                else parm.cluster=rh.val;
            } else if (lh.is_integer() && lh.val==_NLP_INITIALPOINT && rh.type==_VECT) {
                initp.resize(vars.size(),undef);
                vecteur &pnt=*rh._VECTptr,ipl,ipu;
                const_iterateur jt=pnt.begin();
                int iv;
                for (;jt!=pnt.end();++jt) {
                    if (is_equal(*jt) && (iv=indexof(jt->_SYMBptr->feuille._VECTptr->front(),vars))>=0)
                        initp.at(iv)=jt->_SYMBptr->feuille._VECTptr->back();
                    else break;
                }
                if (jt!=pnt.end() || has_inf_or_undef(initp))
                    return generr("Invalid initial point specification");
                if (is_vect_of_intervals(initp,ipl,ipu,contextptr)) {
                    if (!is_numericv(ipl,num_mask) || !is_numericv(ipu,num_mask))
                        return generr("Invalid initial area specification");
                    initp=makevecteur(at_interval,ipl,ipu);
                } else {
                    if (ckmatrix(initp)) {
                        if (!is_numericm(initp,num_mask))
                            return generr("Initial points must be numeric");
                        initp=mtran(initp);
                    } else {
                        if (!is_numericv(initp,num_mask))
                            return generr("Initial point must be numeric");
                        initp=vecteur(1,initp);
                    }
                }
            } else if (is_mcint(lh,_NLP_ITERATIONLIMIT) && rh.is_integer() && rh.val>0) {
                parm.max_iter=rh.val;
            } else if (is_mcint(lh,_NLP_METHOD) && (rh.type==_STRNG || rh.type==_VECT)) {
                string method_name;
                if (rh.type==_STRNG) {
                    method_name=*rh._STRNGptr;
                } else {
                    const vecteur &meth=*rh._VECTptr;
                    if (meth.empty() || meth.front().type!=_STRNG)
                        return generr("Invalid method specification");
                    method_name=*meth.front()._STRNGptr;
                    /* parse method parameters */
                    for (const_iterateur it=meth.begin()+1;it!=meth.end();++it) {
                        if (!is_equal(*it) || it->_SYMBptr->feuille._VECTptr->front().type!=_STRNG)
                            return generr("Invalid method parameter specification");
                        string param_name=*it->_SYMBptr->feuille._VECTptr->front()._STRNGptr;
                        const gen &pval=it->_SYMBptr->feuille._VECTptr->back();
                        if (param_name=="penalty") {
                            if (pval.type==_CPLX) {
                                if (!is_strictly_positive(re(pval,contextptr),contextptr) || !is_strictly_positive(im(pval,contextptr),contextptr))
                                    return generr("Penalty coefficients must be positive");
                                parm.penalty=_evalf(pval,contextptr);
                            } else {
                                if (!is_real_number(pval,contextptr) || is_positive(-pval,contextptr))
                                    return generr("Expected a positive real number");
                                parm.penalty=to_real_number(pval,contextptr).to_double(contextptr);
                            }
                        } else if (param_name=="cross-probability") {
                            if (!is_real_number(pval,contextptr) || !is_positive(pval,contextptr) || !is_greater(1,pval,contextptr))
                                return generr("Expected a real number in [0,1]");
                            parm.cross_prob=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="scale") {
                            if (!is_real_number(pval,contextptr) || !is_positive(pval,contextptr) || !is_greater(2,pval,contextptr))
                                return generr("Expected a real number in [0,2]");
                            parm.scale=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="contract-ratio") {
                            if (!is_real_number(pval,contextptr) || is_positive(-pval,contextptr) || is_positive(pval-1,contextptr))
                                return generr("Expected a real number in (0,1)");
                            parm.contract_ratio=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="shrink-ratio") {
                            if (!is_real_number(pval,contextptr) || is_positive(-pval,contextptr) || is_positive(pval-1,contextptr))
                                return generr("Expected a real number in (0,1)");
                            parm.shrink_ratio=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="reflect-ratio") {
                            if (!is_real_number(pval,contextptr) || is_positive(-pval,contextptr))
                                return generr("Expected a positive real number");
                            parm.reflect_ratio=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="expand-ratio") {
                            if (!is_real_number(pval,contextptr) || is_positive(1-pval,contextptr))
                                return generr("Expected a real number greater than 1");
                            parm.expand_ratio=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="search-points") {
                            if (!pval.is_integer() || pval.val<1)
                                return generr("Expected a positive integer");
                            parm.search_points=pval.val;
                        } else if (param_name=="postprocess") {
                            if (!pval.is_integer() || pval.subtype!=_INT_BOOLEAN)
                                return generr("Expected a boolean value");
                            parm.postprocess=(bool)pval.val;
                        } else if (param_name=="size") {
                            if (!is_real_number(pval,contextptr) || is_positive(-pval,contextptr))
                                return generr("Expected a positive real number");
                            parm.size=to_real_number(pval,contextptr).to_double(contextptr);
                        } else if (param_name=="step") {
                            if (!is_real_number(pval,contextptr) || is_positive(-pval,contextptr))
                                return generr("Expected a positive real number");
                            parm.step=to_real_number(pval,contextptr).to_double(contextptr);
                        } else return generr("Unrecognized method parameter");
                    }
                }
                if (method_name=="differential-evolution" || method_name=="de")
                    method=_NLP_DIFFERENTIAL_EVOLUTION;
                else if (method_name=="automatic" || method_name=="auto")
                    method=_NLP_AUTOMATIC;
                else if (method_name=="nelder-mead" || method_name=="nm")
                    method=_NLP_NELDER_MEAD;
                else if (method_name=="interior-point" || method_name=="ipt")
                    method=_NLP_INTERIOR_POINT;
                else if (method_name=="cobyla" || method_name=="cbl")
                    method=_NLP_COBYLA;
                else return generr("Method not supported");
            } else if (is_mcint(lh,_NLP_PRESOLVE) && rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                parm.presolve=(bool)rh.val;
            } else if (is_mcint(lh,_NLP_MAXIMIZE) && rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                maximize=(bool)rh.val;
            } else if (is_mcint(lh,_NLP_PRECISION) && is_real_number(rh,contextptr) && is_strictly_positive(rh,contextptr)) {
                parm.eps=to_real_number(rh,contextptr).to_double(contextptr);
            } else if (is_mcint(lh,_NLP_TOLERANCE) && is_real_number(rh,contextptr) && is_strictly_positive(rh,contextptr)) {
                parm.tol=to_real_number(rh,contextptr).to_double(contextptr);
            } else if (is_mcint(lh,_NLP_VERBOSE) && rh.is_integer()) {
                parm.msg_level=rh.subtype==_INT_BOOLEAN?((bool)rh.val?3:0):rh.val;
                if (parm.msg_level<0 || parm.msg_level>4)
                    invalid=true;
            } else if (is_mcint(lh,_NLP_INTEGERVARIABLES) && rh.type==_VECT) {
                intvars=*rh._VECTptr;
            } else if (is_mcint(lh,_NLP_BINARYVARIABLES) && rh.type==_VECT) {
                binvars=*rh._VECTptr;
            } else if (rh.is_symb_of_sommet(at_interval)) {
                vecteur vs=gen2vecteur(lh);
                const gen &lb=rh._SYMBptr->feuille._VECTptr->front();
                const gen &ub=rh._SYMBptr->feuille._VECTptr->back();
                if (is_strictly_greater(lb,ub,contextptr))
                    return gensizeerr(gettext("Infeasible bounds on variable(s) ")+_feuille(vs,contextptr).print(contextptr)+".");
                for (const_iterateur vt=vs.begin();vt!=vs.end();++vt) {
                    int j=indexof(*vt,vars);
                    if (j<0)
                        return generr("Expected a variable");
                    lbv[j]=max(lb,lbv[j],contextptr);
                    ubv[j]=min(ub,ubv[j],contextptr);
                }
            } else invalid=true;
        } else invalid=true;
        if (invalid)
            return gensizeerr(gettext("Invalid argument ")+print_INT_(it-gv.begin())+".");
    }
    // convert trivial inequalities to variable bounds
    for (int i=constr.size();i-->0;) {
        gen &c=constr[i],a,b;
        int dir=which_ineq(c);
        if (dir>0) {
            dir=2*dir-3;
            gen d=c._SYMBptr->feuille._VECTptr->front()-c._SYMBptr->feuille._VECTptr->back();
            for (const_iterateur it=vars.begin();it!=vars.end();++it) {
                int j=it-vars.begin();
                gen fa;
                if (is_linear_wrt(d,*it,a,b,contextptr) && is_real_number(a,contextptr) &&
                        !is_zero((fa=to_real_number(a,contextptr)),contextptr) && is_real_number(b,contextptr)) {
                    if (!is_positive(fa,contextptr))
                        dir=-dir;
                    if (dir>0)
                        lbv[j]=max(lbv[j],-b/a,contextptr);
                    else ubv[j]=min(ubv[j],-b/a,contextptr);
                    constr.erase(constr.begin()+i);
                    break;
                }
            }
        }
    }
    if (maximize)
        obj=-obj;
    if (nonneg) for (iterateur it=lbv.begin();it!=lbv.end();++it) {
        *it=max(0,*it,contextptr);
    }
    // handle integral variables
    intvars=mergevecteur(intvars,binvars);
    minlp=!intvars.empty();
    set<int> intvars_indices;
    for (const_iterateur it=intvars.begin();it!=intvars.end();++it) {
        int j=indexof(*it,vars);
        if (j<0) // specified variable not found among problem variables
            return generr("Invalid specification of integer variables");
        intvars_indices.insert(j);
        if (contains(binvars,*it)) { // tighten bounds on binary variable
            lbv[j]=max(lbv[j],0,contextptr);
            ubv[j]=min(ubv[j],1,contextptr);
        }
        // make bounds integral
        lbv[j]=_ceil(lbv[j],contextptr);
        ubv[j]=_floor(ubv[j],contextptr);
    }
    //if (!minlp) parm.convex=0; // check convexity only for MINLP
    // solve the problem
    nlp_problem prob(obj,constr,vars,contextptr);
    prob.set_variable_bounds(lbv,ubv);
    for (set<int>::const_iterator iit=intvars_indices.begin();iit!=intvars_indices.end();++iit) {
        prob.set_integer_variable(*iit);
    }
    prob.enable_debug_messages(parm.msg_level>=4);
    vecteur optsol;
    gen optval;
    int res;
    try {
        res=prob.optimize(method,parm,initp,optsol,optval);
    } catch (const std::runtime_error &e) {
        return generr(e.what(),false);
    }
    switch (res) {
    case _NLP_OPTIMAL:
        return change_subtype(makevecteur(maximize?-optval:optval,_zip(makesequence(at_equal,vars,optsol),contextptr)),_LIST__VECT);
    case _NLP_FAILED:
        return change_subtype(makevecteur(string2gen(gettext("Failed to optimize at given precision"),false),
                                          _zip(makesequence(at_equal,vars,optsol),contextptr)),_LIST__VECT);
    case _NLP_ERROR:
        return generr("Failed to find a solution");
    case _NLP_INFEAS:
        if (optsol.size()!=vars.size())
            return generr("The problem is infeasible");
        return change_subtype(makevecteur(string2gen(gettext("Solution is infeasible"),false),
                                          _zip(makesequence(at_equal,vars,optsol),contextptr)),_LIST__VECT);
    case _NLP_UNBOUNDED:
        return maximize?plus_inf:minus_inf;
    default: // not implemented
        break;
    }
    return undef; // should not be reachable
}
static const char _nlpsolve_s []="nlpsolve";
static define_unary_function_eval (__nlpsolve,&_nlpsolve,_nlpsolve_s);
define_unary_function_ptr5(at_nlpsolve,alias_at_nlpsolve,&__nlpsolve,0,true)

gen compute_invdiff(int n,int k,vecteur &xv,vecteur &yv,map<tprob::ipair,gen> &invdiff,GIAC_CONTEXT) {
    tprob::ipair I=make_pair(n,k);
    assert(n<=k);
    gen res(invdiff[I]);
    if (!is_zero(res,contextptr))
        return res;
    if (n==0)
        return invdiff[I]=yv[k];
    if (n==1)
        return invdiff[I]=(xv[k]-xv[0])/(yv[k]-yv[0]);
    gen d1(compute_invdiff(n-1,n-1,xv,yv,invdiff,contextptr));
    gen d2(compute_invdiff(n-1,k,xv,yv,invdiff,contextptr));
    return invdiff[I]=(xv[k]-xv[n-1])/(d2-d1);
}

gen thiele(int k,vecteur &xv,vecteur &yv,identificateur &var,map<tprob::ipair,gen> &invdiff,GIAC_CONTEXT) {
    if (k==int(xv.size()))
        return gen(0);
    gen phi(compute_invdiff(k,k,xv,yv,invdiff,contextptr));
    return (var-xv[k-1])/(phi+thiele(k+1,xv,yv,var,invdiff,contextptr));
}

/*
 * 'thiele' computes rational interpolation for the given list of points using
 * Thiele's method with continued fractions.
 *
 * Source: http://www.astro.ufl.edu/~kallrath/files/pade.pdf (see page 19)
 *
 * Usage
 * ^^^^^
 *      thiele(data,[v])
 * or   thiele(data_x,data_y,[v])
 *
 * Parameters
 * ^^^^^^^^^^
 *      - data      : list of points [[x1,y1],[x2,y2],...,[xn,yn]]
 *      - v         : identifier or a symbolic expression
 *      - data_x    : list of x coordinates [x1,x2,...,xn]
 *      - data_y    : list of y coordinates [y1,y2,...,yn]
 *
 * The return value is an expression R(v), where R is rational interpolant of
 * the given set of points.
 *
 * Note that the interpolant may have singularities in
 * [min(data_x),max(data_x)].
 *
 * Example
 * ^^^^^^^
 * Function f(x)=(1-x^4)*exp(1-x^3) is sampled on interval [-1,2] in 13
 * equidistant points:
 *
 * data_x:=[-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2],
 * data_y:=[0.0,2.83341735599,2.88770329586,2.75030303645,2.71828182846,
 *          2.66568510781,2.24894558809,1.21863761951,0.0,-0.555711613283,
 *         -0.377871362418,-0.107135851128,-0.0136782294833]
 *
 * To obtain rational function passing through these points, input:
 *      thiele(data_x,data_y,x)
 * Output:
 *      (-1.55286115659*x^6+5.87298387514*x^5-5.4439152812*x^4+1.68655817708*x^3
 *       -2.40784868317*x^2-7.55954205222*x+9.40462512097)/(x^6-1.24295718965*x^5
 *       -1.33526268624*x^4+4.03629272425*x^3-0.885419321*x^2-2.77913222418*x+3.45976823393)
 */
gen _thiele(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return generrtype("Expected a matrix or a sequence of arguments");
    vecteur xv,yv;
    gen x0=identificateur("x");
    int opt_pos=1;
    if (g.subtype==_SEQ__VECT) {
        const vecteur &gv=*g._VECTptr;
        if (gv[0].type!=_VECT)
            return generrtype("Expected a list or a matrix");
        if (ckmatrix(gv[0],false)) {
            matrice m(mtran(*gv[0]._VECTptr));
            if (m.size()!=2)
                return generr("Expected a two-column matrix");
            xv=*m[0]._VECTptr;
            yv=*m[1]._VECTptr;
        } else {
            if (gv.size()<2)
                return generr("Expected a pair of lists");
            if (gv[1].type!=_VECT)
                return generrtype("Expected a list");
            if (gv[0]._VECTptr->size()!=gv[1]._VECTptr->size())
                return generrdim("Expected two lists of equal size");
            xv=*gv[0]._VECTptr;
            yv=*gv[1]._VECTptr;
            opt_pos=2;
        }
        if (gv.size()-opt_pos>=1)
            x0=gv[opt_pos];
    } else {
        if (!ckmatrix(*g._VECTptr,false))
            return generrtype("Expected a matrix");
        matrice m(mtran(*g._VECTptr));
        if (m.size()!=2)
            return generrdim("Expected a two-column matrix");
        xv=*m[0]._VECTptr;
        yv=*m[1]._VECTptr;
    }
    if (xv.size()<2)
        return generr("At least two sample points are required");
    gen x=temp_symb("x",-1,contextptr);
    map<tprob::ipair,gen> invdiff;
    gen rat(yv[0]+thiele(1,xv,yv,*x._IDNTptr,invdiff,contextptr));
    // detect singularities
    gen den(_denom(rat,contextptr));
    matrice sing;
    if (*_lname(den,contextptr)._VECTptr==vecteur(1,x)) {
        for (int i=0;i<int(xv.size())-1;++i) {
            gen y1(_evalf(subst(den,x,xv[i],false,contextptr),contextptr));
            gen y2(_evalf(subst(den,x,xv[i+1],false,contextptr),contextptr));
            if (is_positive(-y1*y2,contextptr))
                sing.push_back(makevecteur(xv[i],xv[i+1]));
        }
    }
    if (!sing.empty()) {
        *logptr(contextptr) << gettext("Warning") << ": " << gettext("interpolant has singularities in ");
        for (int i=0;i<int(sing.size());++i) {
            *logptr(contextptr) << "(" << sing[i][0] << "," << sing[i][1] << ")";
            if (i<int(sing.size())-1)
                *logptr(contextptr) << (i<int(sing.size())-2?", ":" and ");
        }
        *logptr(contextptr) << "\n";
    }
    if (x0.type==_VECT) {
        vecteur res;
        res.reserve(x0._VECTptr->size());
        for (const_iterateur it=x0._VECTptr->begin();it!=x0._VECTptr->end();++it) {
            res.push_back(simp(subst(rat,x,*it,false,contextptr),contextptr));
        }
        return res;
    }
    return simp(subst(rat,x,x0,false,contextptr),contextptr);
}
static const char _thiele_s []="thiele";
static define_unary_function_eval (__thiele,&_thiele,_thiele_s);
define_unary_function_ptr5(at_thiele,alias_at_thiele,&__thiele,0,true)

/* Returns the trigonometric polynomial in variable x passing through points
 * with ordinate componets in 'data' and the abscissa components equally spaced between
 * a and b (the first being equal a and the last being equal to b).
 */
gen triginterp(const vecteur &data,const gen &a,const gen &b,const identificateur &x,GIAC_CONTEXT) {
    int n=data.size();
    if (n<2)
        return generr("At least two data points are required");
    int N=(n%2)==0?n/2:(n-1)/2;
    gen T=(b-a)*fraction(n,n-1),twopi=2*_IDNT_pi(),X;
    matrice cos_coeff=zero_mat(N,n,contextptr);
    matrice sin_coeff=zero_mat(N,n,contextptr);
    for (int k=0;k<n;++k) {
        X=twopi*(a/T+fraction(k,n));
        for (int j=1;j<=N;++j) {
            cos_coeff[j-1]._VECTptr->at(k)=cos(j*X,contextptr);
            sin_coeff[j-1]._VECTptr->at(k)=sin(j*X,contextptr);
        }
    }
    gen tp=_mean(data,contextptr);
    for (int j=0;j<N;++j) {
        gen c=fraction(((n%2)==0 && j==N-1)?1:2,n);
        gen ak=_evalc(trig2exp(scalarproduct(data,*cos_coeff[j]._VECTptr,contextptr),contextptr),contextptr);
        gen bk=_evalc(trig2exp(scalarproduct(data,*sin_coeff[j]._VECTptr,contextptr),contextptr),contextptr);
        tp+=simp(c*ak,contextptr)*cos(simp((j+1)*twopi/T,contextptr)*x,contextptr);
        tp+=simp(c*bk,contextptr)*sin(simp((j+1)*twopi/T,contextptr)*x,contextptr);
    }
    return tp;
}

gen _triginterp(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    vecteur &args=*g._VECTptr;
    if (args.size()<2)
        return generr("Too few arguments");
    if (args.front().type!=_VECT)
        return generrtype("Expected a list");
    vecteur &data=*args.front()._VECTptr;
    gen x,ab,a,b,&vararg=args.at(1);
    if (is_equal(vararg) && (x=_lhs(vararg,contextptr)).type==_IDNT && (ab=_rhs(vararg,contextptr)).is_symb_of_sommet(at_interval)) {
        a=_lhs(ab,contextptr);
        b=_rhs(ab,contextptr);
    } else if (args.size()==4 && (x=args.back()).type==_IDNT) {
        a=args.at(1);
        b=args.at(2);
    } else return generr("Expected x=a..b or x,a,b");
    gen tp=triginterp(data,a,b,*x._IDNTptr,contextptr);
    if (is_approx(data) || is_approx(a) || is_approx(b))
        tp=_evalf(tp,contextptr);
    return tp;
}
static const char _triginterp_s []="triginterp";
static define_unary_function_eval (__triginterp,&_triginterp,_triginterp_s);
define_unary_function_ptr5(at_triginterp,alias_at_triginterp,&__triginterp,0,true)

/* Compute a rational interpolation of the given points using
 * the method of Floater and Hormann.
 * Usage: ratinterp(Mtrx(M) || Vect(X),[Vect(Y)],[Var(x) || Vect(P)],[d])
 * - If P is given, then the interpolant is evaluated at these points
 * - d is a nonnegative integer not larger than |X| (default: d=|X|)
 */
gen _ratinterp(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return generrtype("Expected a matrix or a sequence of arguments");
    vecteur X,Y;
    gen x=temp_symb("x",-1,contextptr),x0=identificateur("x");
    int d=-1,opts_at=2;
    /* parse input arguments */
    if (g.subtype==_SEQ__VECT) {
        const vecteur &gv=*g._VECTptr;
        if (gv.front().type!=_VECT)
            return generrtype("Expected a list or a matrix");
        if (ckmatrix(*gv.front()._VECTptr,false)) {
            const matrice &M=*gv.front()._VECTptr;
            if (M.empty() || M.front()._VECTptr->size()!=2)
                return generr("Expected a two-column matrix");
            matrice tM=mtran(M);
            X=*tM.front()._VECTptr;
            Y=*tM.back()._VECTptr;
            opts_at=1;
        } else {
            if (gv.size()<2)
                return generr("Expected a pair of lists");
            if (gv[1].type!=_VECT)
                return generrtype("Expected a list");
            X=*gv[0]._VECTptr;
            Y=*gv[1]._VECTptr;
            if (X.empty() || X.size()!=Y.size())
                return generrdim("Expected two lists of equal size");
        }
        if (gv.size()-opts_at>=1)
            x0=gv[opts_at];
        if (gv.size()-opts_at>=2) {
            const gen &arg=gv[opts_at+1];
            if (!arg.is_integer() || arg.val<0 || arg.val>=int(X.size()))
                return generr("Expected an integer 0<=d<|X|");
            d=arg.val;
        }
    } else {
        if (!ckmatrix(*g._VECTptr,false))
            return generrtype("Expected a matrix");
        const matrice &M=*g._VECTptr;
        if (M.empty() || M.front()._VECTptr->size()!=2)
            return generr("Expected a two-column matrix");
        matrice tM=mtran(M);
        X=*tM.front()._VECTptr;
        Y=*tM.back()._VECTptr;
    }
    int n=X.size()-1;
    if (n<1)
        return generr("At least two sample points are required");
    if (d<0)
        d=0;
    /* compute the interpolant corresponding to d in barycentric form */
    vecteur w(n+1,0);
    vecteur D(n+1);
    for (int k=0;k<=n;++k) {
        for (int i=0;i<=n-d;++i) {
            if (k-d<=i && i<=k) {
                gen p(i%2?-1:1);
                for (int j=i;j<=i+d;++j) {
                    if (j!=k)
                        p=p*_inv(X[k]-X[j],contextptr);
                }
                w[k]+=p;
            }
        }
        D[k]=w[k]/(x-X[k]);
    }
    gen N=scalarproduct(D,Y,contextptr),r=N/_sum(D,contextptr);
    if (x0.type==_VECT) {
        vecteur res;
        res.reserve(x0._VECTptr->size());
        for (const_iterateur it=x0._VECTptr->begin();it!=x0._VECTptr->end();++it) {
            res.push_back(simp(subst(r,x,*it,false,contextptr),contextptr));
        }
        return res;
    }
    return simp(subst(r,x,x0,false,contextptr),contextptr);
}
static const char _ratinterp_s []="ratinterp";
static define_unary_function_eval (__ratinterp,&_ratinterp,_ratinterp_s);
define_unary_function_ptr5(at_ratinterp,alias_at_ratinterp,&__ratinterp,0,true)

/* select a good bandwidth for kernel density estimation using a direct plug-in method (DPI),
 * Gaussian kernel is assumed */
double select_bandwidth_dpi(const vector<double> &data,double sd) {
    int n=data.size();
    double g6=1.23044723*sd,s=0,t,t2;
    for (vector<double>::const_iterator it=data.begin();it!=data.end();++it) {
        for (vector<double>::const_iterator jt=it+1;jt!=data.end();++jt) {
            t=(*it-*jt)/g6;
            t2=t*t;
            s+=(2*t2*(t2*(t2-15)+45)-30)*std::exp(-t2/2);
        }
    }
    s-=15*n;
    double g4=g6*std::pow(-(6.0*n)/s,1/7.0);
    s=0;
    for (vector<double>::const_iterator it=data.begin();it!=data.end();++it) {
        for (vector<double>::const_iterator jt=it+1;jt!=data.end();++jt) {
            t=(*it-*jt)/g4;
            t2=t*t;
            s+=(2*t2*(t2-6)+6)*std::exp(-t2/2);
        }
    }
    s+=3*n;
    return std::pow(double(n)/(M_SQRT2*s),0.2)*g4;
}

gen fft_sum(const vecteur &c,const vecteur &k,int M,GIAC_CONTEXT) {
    return _scalar_product(makesequence(c,_mid(makesequence(_convolution(makesequence(c,k),contextptr),M,M),contextptr)),contextptr);
}

/* faster bandwidth DPI selector using binned data and FFT */
double select_bandwidth_dpi_bins(int n,const vecteur &c,double d,double sd,GIAC_CONTEXT) {
    int M=c.size();
    vecteur k(2*M+1);
    double g6=1.23044723*sd,s=0,t,t2;
    for (int i=0;i<=2*M;++i) {
        t=d*double(i-M)/g6;
        t2=t*t;
        k[i]=gen((2*t2*(t2*(t2-15)+45)-30)*std::exp(-t2/2));
    }
    s=to_real_number(fft_sum(c,k,M,contextptr),contextptr).to_double(contextptr);
    double g4=g6*std::pow(-(6.0*n)/s,1/7.0);
    for (int i=0;i<=2*M;++i) {
        t=d*double(i-M)/g4;
        t2=t*t;
        k[i]=gen((2*t2*(t2-6)+6)*std::exp(-t2/2));
    }
    s=to_real_number(fft_sum(c,k,M,contextptr),contextptr).to_double(contextptr);
    return std::pow(double(n)/(M_SQRT2*s),0.2)*g4;
}

/* kernel density estimation with Gaussian kernel */
gen kernel_density(const vector<double> &data,double bw,double sd,int bins,double a,double b,int interp,const gen &x,GIAC_CONTEXT) {
    int n=data.size();
    double SQRT_2PI=std::sqrt(2.0*M_PI);
    if (bins<=0) { // return density as a sum of exponential functions, usable for up to few hundred samples
        if (bw<=0)
            bw=select_bandwidth_dpi(data,sd);
        double fac=bw*n*SQRT_2PI;
        gen res(0),h(2.0*bw*bw);
        for (vector<double>::const_iterator it=data.begin();it!=data.end();++it) {
            res+=exp(-pow(x-gen(*it),2)/h,contextptr);
        }
        return res/gen(fac);
    }
    /* FFT method, constructs an approximation on [a,b] with the specified number of bins.
     * If interp>0, interpolation of order interp is performed and the density is returned piecewise. */
    assert(b>a && bins>0);
    double d=(b-a)/(bins-1);
    vecteur c(bins,0);
    int index;
    for (vector<double>::const_iterator it=data.begin();it!=data.end();++it) {
        index=(int)((*it-a)/d+0.5);
        if (index>=0 && index<bins) c[index]+=1;
    }
    if (bw<=0) { // select bandwidth
        if (n<=1000)
            bw=select_bandwidth_dpi(data,sd);
        else bw=select_bandwidth_dpi_bins(n,c,d,sd,contextptr);
        *logptr(contextptr) << gettext("selected bandwidth: ") << bw << "\n";
    }
    int L=std::min(bins-1,(int)std::floor(1+4*bw/d));
    vecteur k(2*L+1);
    for (int i=0;i<=2*L;++i) {
        k[i]=gen(1.0/(n*bw*SQRT_2PI)*std::exp(-std::pow(d*double(i-L)/bw,2)/2.0));
    }
    gen res=_mid(makesequence(_convolution(makesequence(c,k),contextptr),L,bins),contextptr);
    if (interp>0) { // interpolate the obtained points
        int pos0=0;
        if (x.type!=_IDNT) {
            double xd=to_real_number(x,contextptr).to_double(contextptr);
            if (xd<a || xd>=b || (pos0=std::floor((xd-a)/d))>bins-2)
                return 0;
            if (interp==1) {
                const gen &y1=res._VECTptr->at(pos0),&y2=res._VECTptr->at(pos0+1),x1=a+pos0*d;
                return y1+(x-x1)*(y2-y1)/gen(d);
            }
        }
        vecteur pos(bins);
        for (int i=0;i<bins;++i) pos[i]=a+d*i;
        gen X=x.type==_IDNT?x:temp_symb("X",-1,contextptr);
        vecteur p=*_spline(makesequence(pos,res,X,interp),contextptr)._VECTptr;
        vecteur args(0);
        if (x.type==_IDNT)
            args.reserve(2*bins+1);
        for (int i=0;i<bins;++i) {
            if (x.type==_IDNT) {
                args.push_back(i+1<bins?symb_inferieur_strict(X,pos[i]):symb_inferieur_egal(X,pos[i]));
                args.push_back(i==0?0:p[i-1]);
            } else if (i==pos0) res=simp(subst(p[i],X,x,false,contextptr),contextptr);
            if (i+1<bins && interp>1 && !solve_quiet(p[i],symb_equal(X,symb_interval(pos[i],pos[i+1])),contextptr).empty())
                *logptr(contextptr) << gettext("Warning") << ": " << gettext("interpolated density has negative values in") << " ["
                                    << pos[i] << "," << pos[i+1] << "]\n";
        }
        if (x.type!=_IDNT) return res;
        args.push_back(0);
        res=symbolic(at_piecewise,change_subtype(args,_SEQ__VECT));
        return res;
    }
    return res;
}

bool parse_interval(const gen &feu,double &a,double &b,GIAC_CONTEXT) {
    const vecteur &v=*feu._VECTptr;
    gen l=v.front(),r=v.back(),fl,fr;
    if (!is_real_number(l,contextptr) || !is_real_number(r,contextptr) ||
            !is_strictly_greater(fr=to_real_number(r,contextptr),fl=to_real_number(l,contextptr),contextptr))
        return false;
    a=fl.to_double(contextptr); b=fr.to_double(contextptr);
    return true;
}

gen _kernel_density(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return generrtype("Expected a list or a sequence of arguments");
    gen x=identificateur("x");
    double a=0,b=0,bw=0,sd,d,sx=0,sxsq=0;
    int bins=100,interp=1,method=_KDE_METHOD_LIST,bw_method=_KDE_BW_METHOD_DPI;
    if (g.subtype==_SEQ__VECT) {
        // parse options
        for (const_iterateur it=g._VECTptr->begin()+1;it!=g._VECTptr->end();++it) {
            if (is_equal(*it)) {
                const gen &opt=it->_SYMBptr->feuille._VECTptr->front();
                const gen &v=it->_SYMBptr->feuille._VECTptr->back();
                if (opt==_KDE_BANDWIDTH) {
                    if (v==at_select)
                        bw_method=_KDE_BW_METHOD_DPI;
                    else if (v==at_gauss || v==at_normal || v==at_normald)
                        bw_method=_KDE_BW_METHOD_ROT;
                    else {
                        gen fv;
                        if (!is_real_number(v,contextptr) || !is_strictly_positive(fv=to_real_number(v,contextptr),contextptr))
                            return generr("Invalid bandwidth specification");
                        bw=fv.to_double(contextptr);
                    }
                } else if (opt==_KDE_BINS) {
                    if (!v.is_integer() || v.val<=0)
                        return generr("Number of bins must be positive");
                    bins=v.val;
                } else if (opt==at_range) {
                    if (v.type==_VECT) {
                        if (v._VECTptr->size()!=2 || !parse_interval(v,a,b,contextptr))
                            return generr("Invalid range specification");
                    } else if (!v.is_symb_of_sommet(at_interval) ||
                               !parse_interval(v._SYMBptr->feuille,a,b,contextptr))
                        return generr("Invalid range specification");
                } else if (opt==at_output || opt==at_Output) {
                    if (v==at_exact)
                        method=_KDE_METHOD_EXACT;
                    else if (v==at_piecewise)
                        method=_KDE_METHOD_PIECEWISE;
                    else if (v==_MAPLE_LIST)
                        method=_KDE_METHOD_LIST;
                    else return generr("Invalid method specification");
                } else if (opt==at_interp) {
                    if (!v.is_integer() || (interp=v.val)<1)
                        return generr("Invalid interpolation method specification");
                } else if (opt==at_spline) {
                    if (!v.is_integer() || (interp=v.val)<1)
                        return generr("Invalid interpolation method specification");
                    method=_KDE_METHOD_PIECEWISE;
                } else if (opt.type==_IDNT) {
                    x=opt;
                    if (!v.is_symb_of_sommet(at_interval) || !parse_interval(v._SYMBptr->feuille,a,b,contextptr))
                        return generr("Invalid range specification");
                } else if (opt==at_eval) x=v;
                else return generr("Unrecognized option");
            } else if (it->type==_IDNT) x=*it;
            else if (it->is_symb_of_sommet(at_interval)) {
                if (!parse_interval(it->_SYMBptr->feuille,a,b,contextptr))
                    return generr("Invalid range specification");
            } else if (*it==at_exact)
                method=_KDE_METHOD_EXACT;
            else if (*it==at_piecewise)
                method=_KDE_METHOD_PIECEWISE;
            else return generr("Unrecognized option");
        }
    }
    if (x.type!=_IDNT) {
        if (!is_real_number(x,contextptr))
            return generr("Expected a real constant");
        if (method==_KDE_METHOD_LIST)
            return generr("Invalid method");
    }
    vecteur &data=g.subtype==_SEQ__VECT?*g._VECTptr->front()._VECTptr:*g._VECTptr;
    int n=data.size();
    if (n<2)
        return generr("Data too small");
    vector<double> ddata(n);
    gen e;
    for (const_iterateur it=data.begin();it!=data.end();++it) {
        if (!is_real_number(*it,contextptr))
            return generr("Expected a real constant");
        d=ddata[it-data.begin()]=to_real_number(*it,contextptr).to_double(contextptr);
        sx+=d;
        sxsq+=d*d;
    }
    sd=std::sqrt(1/double(n-1)*(sxsq-1/double(n)*sx*sx));
    if (bw_method==_KDE_BW_METHOD_ROT) { // Silverman's rule of thumb
        double iqr=to_real_number(_quartile3(data,contextptr)-_quartile1(data,contextptr),contextptr).to_double(contextptr);
        bw=1.06*std::min(sd,iqr/1.34)*std::pow(double(data.size()),-0.2);
        *logptr(contextptr) << gettext("Selected bandwidth: ") << bw << "\n";
    }
    if (bins>0 && a==0 && b==0) {
        a=to_real_number(_min(data,contextptr),contextptr).to_double(contextptr)-3*bw;
        b=to_real_number(_max(data,contextptr),contextptr).to_double(contextptr)+3*bw;
    }
    if (method==_KDE_METHOD_EXACT)
        bins=0;
    else if (method==_KDE_METHOD_LIST) {
        if (bins<1)
            return generr("Number of bins must be positive");
        interp=0;
    } else if (method==_KDE_METHOD_PIECEWISE) {
        if (bins<1 || interp<1)
            return generr("Number of bins and interpolation specification must be positive");
    }
    return kernel_density(ddata,bw,sd,bins,a,b,interp,x,contextptr);
}
static const char _kernel_density_s []="kernel_density";
static define_unary_function_eval (__kernel_density,&_kernel_density,_kernel_density_s);
define_unary_function_ptr5(at_kernel_density,alias_at_kernel_density,&__kernel_density,0,true)

static const char _kde_s []="kde";
static define_unary_function_eval (__kde,&_kernel_density,_kde_s);
define_unary_function_ptr5(at_kde,alias_at_kde,&__kde,0,true)

/* maximum likelihood estimation for Weibull distribution */
gen weibull_mle(const vecteur &S,const gen &k0,const gen &eps,GIAC_CONTEXT) {
    gen s0(0),s1(0),s2(0),L(0),t,l,n(S.size());
    for (const_iterateur it=S.begin();it!=S.end();++it) {
        L+=(l=ln(*it,contextptr));
        t=exp(k0*l,contextptr);
        s0+=t; s1+=t*l; s2+=t*sq(l);
    }
    L=L/n;
    gen ik0=_inv(k0,contextptr),k=k0-(s1-s0*(ik0+L))/(s2+s0*sq(ik0)-s1*(ik0+L));
    if (is_positive(-k,contextptr))
        return undef;
    if (is_greater(_abs(k-k0,contextptr),eps,contextptr))
        return weibull_mle(S,k,eps,contextptr);
    gen lambda=exp(_inv(k,contextptr)*ln(s0/n,contextptr),contextptr);
    return symbolic(at_weibulld,makesequence(k,lambda));
}

/* maximum likelihood estimation for Cauchy distribution */
gen cauchy_mle(const vecteur &S,const gen &x0_init,const gen &gama_init,const gen &eps,GIAC_CONTEXT) {
    matrice m=zero_mat(2,2,contextptr);
    gen n(S.size()),d,dsq,den,densq,gsq=sq(gama_init),F(0),G(-n/2);
    gen &Fx0=m.front()._VECTptr->front(),&Fgama=m.front()._VECTptr->back(),
            &Gx0=m.back()._VECTptr->front(),&Ggama=m.back()._VECTptr->back();
    for (const_iterateur it=S.begin();it!=S.end();++it) {
        d=*it-x0_init; dsq=sq(d); den=gsq+dsq; densq=sq(den);
        F+=d/den; G+=gsq/den;
        Fx0+=(dsq-gsq)/densq; Fgama+=d/densq; Ggama+=dsq/densq;
    }
    Gx0=2*gsq*Fgama;
    Fgama=-2*gama_init*Fgama;
    Ggama=2*gama_init*Ggama;
    gen tmpdelta=_linsolve(makesequence(m,makevecteur(-F,-G)),contextptr);
    if (tmpdelta.type==_VECT) {
        vecteur &delta=*tmpdelta._VECTptr;
        if (delta.empty())
            return undef;
        gen x0=x0_init+delta.front(),gama=gama_init+delta.back();
        if (is_greater(l2norm(delta,contextptr),eps,contextptr))
            return cauchy_mle(S,x0,gama,eps,contextptr);
        return symbolic(at_cauchyd,makesequence(x0,gama));
    }
    return undef;
}

/* fit distribution of the given type to the given data using the method of maximum likelihood */
gen _fitdistr(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    if (g._VECTptr->size()!=2)
        return generrdim("Expected a sequence of two arguments");
    if (g._VECTptr->front().type!=_VECT)
        return generr("Expected a list");
    vecteur &S=*g._VECTptr->front()._VECTptr;
    int n=S.size(); // number of samples
    if (n<2)
        return generr("At least two samples are required");
    const gen &dist=g._VECTptr->back();
    gen N(n);
    /* compute sample mean and variance */
    gen mean(0),var(0),ev;
    for (const_iterateur it=S.begin();it!=S.end();++it) {
        if (!is_real_number(*it,contextptr))
            return generr("Expected a real constant");
        mean+=to_real_number(*it,contextptr);
    }
    mean=mean/(N-1);
    for (const_iterateur it=S.begin();it!=S.end();++it) {
        var+=sq(*it-mean);
    }
    var=var/(N-1);
    gen sdev=sqrt(var,contextptr);
    /* fit the parameters of the specified distribution by the method of maximum likelihood */
    if (dist==at_normal || dist==at_normald || dist==at_NORMALD) {
        return symbolic(at_normald,makesequence(mean,sdev));
    } else if (dist==at_poisson || dist==at_POISSON) {
        for (const_iterateur it=S.begin();it!=S.end();++it) {
            if (!it->is_integer() || it->val<0) return generr("Expected a nonnegative integer");
        }
        return symbolic(at_poisson,mean);
    } else if (dist==at_exp || dist==at_EXP || dist==at_exponential || dist==at_exponentiald) {
        for (const_iterateur it=S.begin();it!=S.end();++it) {
            if (is_positive(-*it,contextptr)) return generr("Expected a positive number");
        }
        return symbolic(at_exponentiald,_inv(mean,contextptr));
    } else if (dist==at_geometric) {
        for (const_iterateur it=S.begin();it!=S.end();++it) {
            if (!it->is_integer() || it->val<1) return generr("Expected a positive integer");
        }
        return symbolic(at_geometric,_inv(mean,contextptr));
    } else if (dist==at_gammad || dist==at_Gamma) {
        gen slog(0);
        for (const_iterateur it=S.begin();it!=S.end();++it) {
            if (is_positive(-*it,contextptr)) return generr("Expected a positive number");
            slog+=ln(*it,contextptr);
        }
        gen a_init=sq(mean)/var,aidn=temp_symb("a",-1,contextptr);
        gen e=ln(aidn,contextptr)-Psi(aidn,contextptr)-ln(mean,contextptr)+slog/N;
        gen a=_solve(makesequence(e,symb_equal(aidn,a_init),_NEWTON_SOLVER),contextptr);
        return symbolic(at_gammad,makesequence(a,a/mean));
    } else if (dist==at_betad || dist==at_Beta) {
        gen slog(0),s1log(0),aidn=temp_symb("a",-1,contextptr),bidn=temp_symb("b",-1,contextptr);
        for (const_iterateur it=S.begin();it!=S.end();++it) {
            if (!is_greater(*it,0,contextptr) || is_strictly_greater(*it,1,contextptr))
                return generr("Expected a number in [0,1]");
            slog+=ln(*it,contextptr); s1log+=ln(1-*it,contextptr);
        }
        gen fac=(var+sq(mean)-mean)/var,a_init=-mean*fac,b_init=(mean-1)*fac;
        gen e1=Psi(aidn,contextptr)-Psi(aidn+bidn,contextptr)-slog/N;
        gen e2=Psi(bidn,contextptr)-Psi(aidn+bidn,contextptr)-s1log/N;
        gen tmpsol=_fsolve(makesequence(makevecteur(e1,e2),makevecteur(aidn,bidn),makevecteur(a_init,b_init),_NEWTONJ_SOLVER),contextptr);
        if (tmpsol.type!=_VECT || tmpsol._VECTptr->size()!=2)
            return generr("Unable to compute parameters for Beta distribution");
        vecteur &sol=*tmpsol._VECTptr;
        return symbolic(at_betad,change_subtype(sol,_SEQ__VECT));
    } else if (dist==at_cauchy || dist==at_cauchyd) {
        gen x0_init=_median(S,contextptr);
        gen gama_init=(_quartile3(S,contextptr)-_quartile1(S,contextptr))/2;
        return cauchy_mle(S,x0_init,gama_init,1e-5,contextptr);
    } else if (dist==at_weibull || dist==at_weibulld) {
        if (is_zero(var)) return generr("Expected a nonzero variance");
        gen kidn=temp_symb("k",-1,contextptr),slog(0);
        for (const_iterateur it=S.begin();it!=S.end();++it) {
            if (is_positive(-*it,contextptr)) return generr("Expected a positive number");
            slog+=ln(*it,contextptr);
        }
        gen e=_Gamma(1+gen(2)/kidn,contextptr)/_Gamma(1+gen(1)/kidn,contextptr)-var/sq(mean)-1;
        gen k_init=_fsolve(makesequence(e,symb_equal(kidn,max(1,_inv(var,contextptr),contextptr)),_NEWTON_SOLVER),contextptr);
        return weibull_mle(S,k_init,1e-5,contextptr);
    }
    return generr("Unrecognized probability distribution");
}
static const char _fitdistr_s []="fitdistr";
static define_unary_function_eval (__fitdistr,&_fitdistr,_fitdistr_s);
define_unary_function_ptr5(at_fitdistr,alias_at_fitdistr,&__fitdistr,0,true)

/* evaluates the function f(x,y,y') in the specified point and returns the result as a floating point number. */
double eval_func(const gen &f,const vecteur &vars,const gen &x,const gen &y,const gen &dy,bool &errflag,GIAC_CONTEXT) {
    gen e=subst(f,vars,makevecteur(x,y,dy),false,contextptr);
    if (!is_real_number(e,contextptr)) {
        errflag=false;
        return 0;
    }
    return to_real_number(e,contextptr).to_double(contextptr);
}

/* approximate the solution of the following boundary-value problem:
 * y''=f(x,y,y'), a<=x<=b, y(a)=alpha, y(b)=beta.
 * The solution is stored in the lists X, w1 and w2 such that
 * X=[x0=a,x1,x2,..,xN=b] and w1[k]=y(xk), w2[k]=y'(xk) for k=0,1,...,N.
 * Return value: 0 on success, 1 if the maximum number of iterations M is exceeded and 2 on computation failure.
 */
int shooting(const gen &f,const gen &x_idn,const gen &y_idn,const gen &dy_idn,const gen &TK_orig,
             const gen &x1,const gen &x2,const gen &y1,const gen &y2,
             int N,double tol,int M,vecteur &X,vecteur &Y,vecteur &dY,GIAC_CONTEXT) {
    gen dfy=derive(f,y_idn,contextptr),dfdy=derive(f,dy_idn,contextptr);
    vecteur vars=makevecteur(x_idn,y_idn,dy_idn);
    double a=x1.to_double(contextptr),b=x2.to_double(contextptr),alpha=y1.to_double(contextptr),beta=y2.to_double(contextptr);
    double h=(b-a)/N,x,k11,k12,k21,k22,k31,k32,k41,k42,dk11,dk12,dk21,dk22,dk31,dk32,dk41,dk42,u1,u2,fv1,fv2,w1i,w2i;
    double TK=is_undef(TK_orig)?(beta-alpha)/(b-a):TK_orig.to_double(contextptr);
    vector<double> w1(N+1),w2(N+1);
    int k=1;
    dfy=simp(dfy,contextptr);
    dfdy=simp(dfdy,contextptr);
    bool ef=true;
    while (k<=M) {
        w1[0]=alpha;
        w2[0]=TK;
        u1=0;
        u2=1;
        for (int i=0;i<N;++i) {
            x=a+i*h;
            w1i=w1[i];
            w2i=w2[i];
            k11=h*w2i;
            k12=h*eval_func(f,vars,x,w1i,w2i,ef,contextptr);
            k21=h*(w2i+k12/2);
            k22=h*eval_func(f,vars,x+h/2,w1i+k11/2,w2i+k12/2,ef,contextptr);
            k31=h*(w2i+k22/2);
            k32=h*eval_func(f,vars,x+h/2,w1i+k21/2,w2i+k22/2,ef,contextptr);
            k41=h*(w2i+k32);
            k42=h*eval_func(f,vars,x+h,w1i+k31,w2i+k32,ef,contextptr);
            w1[i+1]=w1i+(k11+2*k21+2*k31+k41)/6;
            w2[i+1]=w2i+(k12+2*k22+2*k32+k42)/6;
            dk11=h*u2;
            dk12=h*(eval_func(dfy,vars,x,w1i,w2i,ef,contextptr)*u1+
                    eval_func(dfdy,vars,x,w1i,w2i,ef,contextptr)*u2);
            dk21=h*(u2+dk12/2);
            fv1=eval_func(dfy,vars,x+h/2,w1i,w2i,ef,contextptr);
            fv2=eval_func(dfdy,vars,x+h/2,w1i,w2i,ef,contextptr);
            dk22=h*(fv1*(u1+dk11/2)+fv2*(u2+dk12/2));
            dk31=h*(u2+dk22/2);
            dk32=h*(fv1*(u1+dk21/2)+fv2*(u2+dk22/2));
            dk41=h*(u2+dk32);
            dk42=h*(eval_func(dfy,vars,x+h,w1i,w2i,ef,contextptr)*(u1+dk31)+
                    eval_func(dfdy,vars,x+h,w1i,w2i,ef,contextptr)*(u2+dk32));
            u1+=(dk11+2*dk21+2*dk31+dk41)/6;
            u2+=(dk12+2*dk22+2*dk32+dk42)/6;
            if (!ef) return 2;
        }
        if (std::abs(w1[N]-beta)<=tol) {
            X.resize(N+1);
            Y.resize(N+1);
            dY.resize(N+1);
            for (int i=0;i<=N;++i) {
                X[i]=a+i*h;
                Y[i]=w1[i];
                dY[i]=w2[i];
            }
            return 0; // success
        }
        TK-=(w1[N]-beta)/u1;
        ++k;
    }
    return 1; // max number of iterations is exceeded
}

/* the finite-difference method as an slower but more stable alternative to shooting method */
int finitediff(const gen &f,const gen &x_idn,const gen &y_idn,const gen &dy_idn,const gen &x1,const gen &x2,
               const gen &y1,const gen &y2,int N,double tol,int M,vecteur &X,vecteur &Y,GIAC_CONTEXT) {
    gen dfy=derive(f,y_idn,contextptr),dfdy=derive(f,dy_idn,contextptr);
    vecteur vars=makevecteur(x_idn,y_idn,dy_idn);
    double a=x1.to_double(contextptr),b=x2.to_double(contextptr),alpha=y1.to_double(contextptr),beta=y2.to_double(contextptr);
    double h=(b-a)/(N+1),fac=(beta-alpha)/(b-a)*h,x,t;
    vector<double> W(N+2,alpha),A(N+1),B(N+1),C(N+1),D(N+1),U(N+1),L(N+1),Z(N+1);
    vecteur V(N+1);
    W[N+1]=beta;
    for (int i=1;i<=N;++i) W[i]+=i*fac;
    int k=1;
    bool ef=true;
    dfy=simp(dfy,contextptr);
    dfdy=simp(dfdy,contextptr);
    while (k<=M) {
        x=a+h;
        t=(W[2]-alpha)/(2*h);
        A[1]=2+h*h*eval_func(dfy,vars,x,W[1],t,ef,contextptr);
        B[1]=h*eval_func(dfdy,vars,x,W[1],t,ef,contextptr)/2-1;
        D[1]=W[2]+alpha-2*W[1]-h*h*eval_func(f,vars,x,W[1],t,ef,contextptr);
        for (int i=2;i<N;++i) {
            x=a+i*h;
            t=(W[i+1]-W[i-1])/(2*h);
            A[i]=2+h*h*eval_func(dfy,vars,x,W[i],t,ef,contextptr);
            B[i]=h*eval_func(dfdy,vars,x,W[i],t,ef,contextptr)/2-1;
            C[i]=-h*eval_func(dfdy,vars,x,W[i],t,ef,contextptr)/2-1;
            D[i]=W[i+1]+W[i-1]-2*W[i]-h*h*eval_func(f,vars,x,W[i],t,ef,contextptr);
        }
        x=b-h;
        t=(beta-W[N-1])/(2*h);
        A[N]=2+h*h*eval_func(dfy,vars,x,W[N],t,ef,contextptr);
        C[N]=-h*eval_func(dfdy,vars,x,W[N],t,ef,contextptr)/2-1;
        D[N]=W[N-1]+beta-2*W[N]-h*h*eval_func(f,vars,x,W[N],t,ef,contextptr);
        if (!ef) return 2;
        L[1]=A[1];
        U[1]=B[1]/A[1];
        Z[1]=D[1]/L[1];
        for (int i=2;i<N;++i) {
            L[i]=A[i]-C[i]*U[i-1];
            U[i]=B[i]/L[i];
            Z[i]=(D[i]-C[i]*Z[i-1])/L[i];
        }
        L[N]=A[N]-C[N]*U[N-1];
        Z[N]=(D[N]-C[N]*Z[N-1])/L[N];
        V[N]=Z[N];
        W[N]+=Z[N];
        for (int i=N;i-->1;) {
            V[i]=gen(Z[i])-gen(U[i])*V[i+1];
            W[i]+=Z[i]-U[i]*V[i+1].to_double(contextptr);
        }
        if (is_greater(tol,l2norm(V,contextptr),contextptr)) {
            X.resize(N+2);
            Y.resize(N+2);
            for (int i=0;i<=N+1;++i) {
                X[i]=a+i*h;
                Y[i]=W[i];
            }
            return 0; // success
        }
        ++k;
    }
    return 1; // maximum number of iterations is exceeded
}

gen _bvpsolve(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    gen tk(undef);
    int maxiter=RAND_MAX;
    if (gv.size()<3)
        return generr("Too few input arguments");
    gen arg2=_eval(gv[1],contextptr),arg3=_eval(gv[2],contextptr);
    if (arg2.type!=_VECT || arg3.type!=_VECT || arg2._VECTptr->size()!=2 || arg3._VECTptr->size()<2 || arg3._VECTptr->size()>3)
        return generr("Invalid input arguments");
    gen f=idnteval(gv.front(),contextptr);
    const gen &t=arg2._VECTptr->front(),&y=arg2._VECTptr->back();
    gen y1=arg3._VECTptr->at(0),y2=arg3._VECTptr->at(1);
    if (arg3._VECTptr->size()==3 && !is_real_number(arg3._VECTptr->at(2),contextptr))
        return generr("Expected a real constant");
    tk=to_real_number(arg3._VECTptr->at(2),contextptr);
    if (y.type!=_IDNT)
        return generr("Expected an identifier");
    if (!is_equal(t) || t._SYMBptr->feuille._VECTptr->front().type!=_IDNT ||
            !t._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_interval))
        return generrtype("Expected a variable with range");
    const gen &x=t._SYMBptr->feuille._VECTptr->front();
    const vecteur &rng=*t._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr;
    gen x1=rng.front(),x2=rng.back();
    if (!is_real_number(x1,contextptr) || !is_real_number(x2,contextptr) ||
            !is_real_number(y1,contextptr) || !is_real_number(y2,contextptr))
        return generr("Expected a real constant");
    x1=to_real_number(x1,contextptr); x2=to_real_number(x2,contextptr);
    if (!is_strictly_greater(x2,x1,contextptr))
        return gensizeerr("It should be x1<x2");
    y1=to_real_number(y1,contextptr); y2=to_real_number(y2,contextptr);
    int N=100;
    /* parse options */
    int output_type=_BVP_LIST;
    for (const_iterateur it=gv.begin()+3;it!=gv.end();++it) {
        if (is_equal(*it)) {
            const gen &lh=it->_SYMBptr->feuille._VECTptr->front();
            const gen &rh=it->_SYMBptr->feuille._VECTptr->back();
            if (lh==at_output || lh==at_Output) {
                if (rh==_MAPLE_LIST)
                    output_type=_BVP_LIST;
                else if (rh==at_derive)
                    output_type=_BVP_DIFF;
                else if (rh==at_piecewise)
                    output_type=_BVP_PIECEWISE;
                else if (rh==at_spline)
                    output_type=_BVP_SPLINE;
                else return generr("Invalid output specification");
            } else if (lh==at_limit) {
                if (!rh.is_integer() || (maxiter=rh.val)<1)
                    return generr("Maximum number of iterations must be a positive integer");
            } else return generr("Unrecognized option");
        } else if (it->is_integer()) {
            if ((N=it->val)<2)
                return generr("Too small subdivision size, N>=2 is required");
        } else return generr("Unrecognized option");
    }
    gen dy=temp_symb("dy",-1,contextptr);
    gen F=subst(f,derive(symb_of(y,x),x,contextptr),dy,false,contextptr);
    F=subst(F,symb_of(y,x),y,false,contextptr);
    F=subst(F,symbolic(at_derive,y),dy,false,contextptr);
    vecteur X,Y,dY;
    double tol=_epsilon(change_subtype(vecteur(0),_SEQ__VECT),contextptr).to_double(contextptr);
    int ec=shooting(F,x,y,dy,tk,x1,x2,y1,y2,N,tol,maxiter,X,Y,dY,contextptr);
    if (ec==1) {
        print_error("maximum number of iterations exceeded",contextptr);
        return undef;
    }
    if (ec==2) {
        print_error("shooting method failed to converge",contextptr);
        if (is_undef(tk))
            *logptr(contextptr) << ", " << gettext("try to set an initial guess for ") << y << "'(" << x1 << ")";
        *logptr(contextptr) << "\n";
        if (N>=3 && (output_type==_BVP_LIST || output_type==_BVP_PIECEWISE)) {
            *logptr(contextptr) << gettext("Trying the finite-difference method instead") << "\n";
            ec=finitediff(F,x,y,dy,x1,x2,y1,y2,N-1,tol,maxiter,X,Y,contextptr);
            if (ec==2) {
                print_error("failed to converge",contextptr);
                return undef;
            }
            if (ec==1) {
                print_error("maximum number of iterations exceeded",contextptr);
                return undef;
            }
        } else return undef;
    }
    vecteur res,coeff;
    matrice m=zero_mat(4,4,contextptr);
    m[0]._VECTptr->at(3)=m[1]._VECTptr->at(3)=m[2]._VECTptr->at(2)=m[3]._VECTptr->at(2)=1;
    switch (output_type) {
    case _BVP_LIST:
        res.reserve(N+1);
        for (int i=0;i<=N;++i)
            res.push_back(makevecteur(X[i],Y[i]));
        break;
    case _BVP_DIFF:
        res.reserve(N+1);
        for (int i=0;i<=N;++i)
            res.push_back(makevecteur(X[i],Y[i],dY[i]));
        break;
    case _BVP_PIECEWISE:
    case _BVP_SPLINE:
        res.reserve(2*(N+1)+1);
        for (int i=0;i<=N;++i) {
            res.push_back(i<N?symb_inferieur_strict(x,X[i]):symb_inferieur_egal(x,X[i]));
            if (i==0)
                res.push_back(0);
            else if (output_type==_BVP_PIECEWISE)
                res.push_back(Y[i-1]+(x-X[i-1])*(Y[i]-Y[i-1])/(X[i]-X[i-1]));
            else {
                m[0]._VECTptr->at(0)=pow(X[i-1],3);
                m[0]._VECTptr->at(1)=pow(X[i-1],2);
                m[0]._VECTptr->at(2)=X[i-1];
                m[1]._VECTptr->at(0)=pow(X[i],3);
                m[1]._VECTptr->at(1)=pow(X[i],2);
                m[1]._VECTptr->at(2)=X[i];
                m[2]._VECTptr->at(0)=3*m[0][1];
                m[2]._VECTptr->at(1)=2*X[i-1];
                m[3]._VECTptr->at(0)=3*m[1][1];
                m[3]._VECTptr->at(1)=2*X[i];
                coeff=*_linsolve(makesequence(m,makevecteur(Y[i-1],Y[i],dY[i-1],dY[i])),contextptr)._VECTptr;
                res.push_back(pow(x,3)*coeff[0]+pow(x,2)*coeff[1]+x*coeff[2]+coeff[3]);
            }
        }
        res.push_back(0);
        return symbolic(at_piecewise,change_subtype(res,_SEQ__VECT));
    default:
        break;
    }
    return res;
}
static const char _bvpsolve_s []="bvpsolve";
static define_unary_function_eval (__bvpsolve,&_bvpsolve,_bvpsolve_s);
define_unary_function_ptr5(at_bvpsolve,alias_at_bvpsolve,&__bvpsolve,_QUOTE_ARGUMENTS,true)

gen strip_constants(const gen &g,GIAC_CONTEXT) {
    if (g.is_symb_of_sommet(at_neg))
        return g._SYMBptr->feuille;
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT) {
        const vecteur &feu=*g._SYMBptr->feuille._VECTptr;
        gen ret(1);
        for (const_iterateur it=feu.begin();it!=feu.end();++it) {
            if (is_real_number(*it,contextptr))
                continue;
            ret=ret*strip_constants(*it,contextptr);
        }
        return ret;
    }
    if (g.is_symb_of_sommet(at_inv))
        return _inv(strip_constants(g._SYMBptr->feuille,contextptr),contextptr);
    return g;
}

/* return the expression for conjugate points */
gen _conjugate_equation(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    if (gv.size()!=5)
        return generr("Invalid number of input arguments");
    const gen &y0=gv[0],&parm=gv[1],&pvals=gv[2],&t=gv[3],&a=gv[4];
    if (y0.type!=_SYMB || t.type!=_IDNT || parm.type!=_VECT || parm._VECTptr->size()!=2 ||
            pvals.type!=_VECT || pvals._VECTptr->size()!=2 || !is_real_number(a,contextptr))
        return generr("Invalid input arguments");
    const gen &alpha=parm._VECTptr->front(),&beta=parm._VECTptr->back();
    if (alpha.type!=_IDNT || beta.type!=_IDNT)
        return generr("Expected pair of identifiers");
    gen y1=derive(y0,alpha,contextptr),y2=derive(y0,beta,contextptr);
    gen ret=_collect(simp(subst(y1*subst(y2,t,a,false,contextptr)-y2*subst(y1,t,a,false,contextptr),
                                     parm,pvals,false,contextptr),contextptr),contextptr);
    return strip_constants(ret,contextptr);
}
static const char _conjugate_equation_s []="conjugate_equation";
static define_unary_function_eval (__conjugate_equation,&_conjugate_equation,_conjugate_equation_s);
define_unary_function_ptr5(at_conjugate_equation,alias_at_conjugate_equation,&__conjugate_equation,0,true)

static int cnst_count=0;

/* return the (list of) Euler-Lagrange equation(s) for functional L(u,du,t) */
gen _euler_lagrange(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    gen L,t=identificateur("x");
    vecteur u=makevecteur(identificateur("y"));
    if (g.type!=_VECT) {
        L=idnteval(g,contextptr);
        if (!contains(lidnt(L),t))
            t=t__IDNT_e;
    } else {
        if (g.subtype!=_SEQ__VECT)
            return gentypeerr(contextptr);
        const vecteur &gv=*g._VECTptr;
        L=idnteval(gv.front(),contextptr);
        if (gv.size()>1) {
            gen arg=_eval(gv[1],contextptr);
            if (gv.size()==2 && arg.is_symb_of_sommet(at_of)) {
                u.front()=arg._SYMBptr->feuille._VECTptr->front();
                t=arg._SYMBptr->feuille._VECTptr->back();
            } else if (gv.size()==2 && arg.type==_VECT && !arg._VECTptr->empty()) {
                u.clear();
                t=undef;
                for (const_iterateur it=arg._VECTptr->begin();it!=arg._VECTptr->end();++it) {
                    if (!it->is_symb_of_sommet(at_of))
                        return generr("Expected a function");
                    u.push_back(it->_SYMBptr->feuille._VECTptr->front());
                    if (is_undef(t))
                        t=it->_SYMBptr->feuille._VECTptr->back();
                    else if (t!=it->_SYMBptr->feuille._VECTptr->back())
                        return generr("Expected a single parameter variable");
                }
            } else t=arg;
            if (t.type!=_IDNT)
                return generr("Expected a variable");
        }
        if (gv.size()>2) {
            gen arg=_eval(gv[2],contextptr);
            if (arg.type==_IDNT)
                u.front()=arg;
            else if (arg.type==_VECT)
                u=*arg._VECTptr;
            else return generr("Expected a variable or a list of variables");
        }
    }
    L=idnteval(L,contextptr);
    int n=u.size();
    vecteur du(n),Du(n),Dut(n),DU(n),D2U(n),d2u(n),ut(n);
    for (int i=0;i<n;++i) {
        if (u[i].type!=_IDNT)
            return generr("Expected an identifier");
        ut[i]=symb_of(u[i],t);
        du[i]=temp_symb("du",i,contextptr);
        d2u[i]=temp_symb("d2u",i,contextptr);
        Du[i]=symbolic(at_derive,u[i]);
        Dut[i]=symb_of(Du[i],t);
        DU[i]=derive(ut[i],t,contextptr);
        D2U[i]=derive(ut[i],t,2,contextptr);
    }
    L=subst(L,Dut,du,false,contextptr);
    L=subst(L,Du,du,false,contextptr);
    L=subst(L,DU,du,false,contextptr);
    L=subst(L,ut,u,false,contextptr);
    vecteur ret;
    if (n==1 && !depend(L,*t._IDNTptr)) {
        ret.push_back(symb_equal(simp(du[0]*derive(L,du[0],contextptr)-L,contextptr),
                identificateur("K_"+print_INT_(cnst_count++))));
    }
    for (int i=0;i<n;++i) {
        if (!depend(L,*u[i]._IDNTptr)) {
            ret.push_back(symb_equal(simp(derive(L,du[i],contextptr),contextptr),
                                     identificateur("K_"+print_INT_(cnst_count++))));
        } else {
            gen eq=derive(L,u[i],contextptr);
            eq-=derive(subst(derive(L,du[i],contextptr),makevecteur(du[i],u[i]),
                             makevecteur(DU[i],ut[i]),false,contextptr),t,contextptr);
            eq=subst(eq,makevecteur(DU[i],D2U[i]),makevecteur(du[i],d2u[i]),false,contextptr);
            eq=symb_equal(simp(subst(eq,ut[i],u[i],false,contextptr),contextptr),0);
            vecteur sol;
            if (depend(eq,*d2u[i]._IDNTptr) && (sol=solve_quiet(eq,d2u[i],contextptr)).empty())
                eq=symb_equal(d2u[i],simp(sol.front(),contextptr));
            ret.push_back(eq);
        }
    }
    ret=subst(ret,ut,u,false,contextptr);
    gen tmprs=radsimp(ret,contextptr);
    if (tmprs.type==_VECT && tmprs._VECTptr->size()==ret.size())
        ret=*tmprs._VECTptr;
    ret=subst(subst(ret,u,ut,false,contextptr),mergevecteur(du,d2u),
                    mergevecteur(DU,D2U),false,contextptr);
    return ret.size()==1?ret.front():ret;
}
static const char _euler_lagrange_s []="euler_lagrange";
static define_unary_function_eval_quoted (__euler_lagrange,&_euler_lagrange,_euler_lagrange_s);
define_unary_function_ptr5(at_euler_lagrange,alias_at_euler_lagrange,&__euler_lagrange,_QUOTE_ARGUMENTS,true)

gen parse_functional(const gen &L,const gen &t,const gen &y,const gen &dy,GIAC_CONTEXT) {
    assert(t.type==_IDNT && y.type==_IDNT && dy.type==_IDNT);
    gen ret=subst(L,symb_of(symbolic(at_derive,y),t),dy,false,contextptr);
    ret=subst(ret,symbolic(at_derive,y),dy,false,contextptr);
    ret=subst(ret,derive(symb_of(y,t),t,contextptr),dy,false,contextptr);
    ret=subst(ret,symb_of(y,t),y,false,contextptr);
    return ret;
}

/* return the Jacobi equation(s) for the functional L(y,y',a<=t<=b) and stationary function Y,
 * use h as unknown with h(a)=0 */
gen _jacobi_equation(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    gen L,t,y,Y,h,a;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    L=idnteval(gv.front(),contextptr);
    vecteur args=*_eval(vecteur(gv.begin()+1,gv.end()),contextptr)._VECTptr;
    if (args.size()==4) {
        if (!args[0].is_symb_of_sommet(at_of))
            return generr("Expected a function");
        y=args[0]._SYMBptr->feuille._VECTptr->front();
        t=args[0]._SYMBptr->feuille._VECTptr->back();
        Y=args[1]; h=args[2]; a=args[3];
    } else if (args.size()==5) {
        t=args[0]; y=args[1]; Y=args[2]; h=args[3]; a=args[4];
    } else return generr("Invalid number of input arguments");
    if (t.type!=_IDNT || h.type!=_IDNT || y.type!=_IDNT || !is_real_number(a,contextptr))
        return generr("Invalid input arguments");
    L=idnteval(L,contextptr);
    gen dy=temp_symb("dy",-1,contextptr);
    L=parse_functional(L,t,y,dy,contextptr);
    gen dY=derive(Y,t,contextptr);
    gen Ldydy=subst(simp(derive(L,dy,2,contextptr),contextptr),
                    makevecteur(y,dy),makevecteur(Y,dY),false,contextptr);
    gen Lyy=subst(simp(derive(L,y,2,contextptr),contextptr),
                  makevecteur(y,dy),makevecteur(Y,dY),false,contextptr);
    gen Lydy=subst(simp(derive(derive(L,y,contextptr),dy,contextptr),contextptr),
                   makevecteur(y,dy),makevecteur(Y,dY),false,contextptr);
    gen d=simp(Lyy-derive(Lydy,t,contextptr),contextptr);
    if (is_zero(Ldydy))
        return change_subtype(is_zero(d),_INT_BOOLEAN);
    gen jeq=simp(derive(Ldydy,t,contextptr),contextptr)*derive(symb_of(h,t),t,contextptr);
    jeq+=Ldydy*derive(symb_of(h,t),t,2,contextptr)-d*symb_of(h,t);
    jeq=symb_equal(jeq,0);
    gen sol=_deSolve(makesequence(makevecteur(jeq,symb_equal(symb_of(h,a),0)),makevecteur(t,h)),contextptr);
    if (sol.type==_VECT && sol._VECTptr->size()==1 && sol._VECTptr->front().type==_STRNG)
        return jeq;
    return makesequence(jeq,sol.is_symb_of_sommet(at_prod) || sol.is_symb_of_sommet(at_neg)?simp(sol,contextptr):sol);
}
static const char _jacobi_equation_s []="jacobi_equation";
static define_unary_function_eval_quoted (__jacobi_equation,&_jacobi_equation,_jacobi_equation_s);
define_unary_function_ptr5(at_jacobi_equation,alias_at_jacobi_equation,&__jacobi_equation,_QUOTE_ARGUMENTS,true)

gen makevars(const gen &e,const gen &t,const vecteur &depvars,const vecteur &diffvars,GIAC_CONTEXT) {
    if (e.is_symb_of_sommet(at_of) && e._SYMBptr->feuille._VECTptr->back()==t) {
        const gen &u=e._SYMBptr->feuille._VECTptr->front();
        for (const_iterateur it=depvars.begin();it!=depvars.end();++it) {
            if (*it==u) return u;
            if (u==symbolic(at_derive,*it)) return diffvars[it-depvars.begin()];
        }
    } else if (e.is_symb_of_sommet(at_derive)) {
        const gen &feu=e._SYMBptr->feuille;
        if (feu.type!=_VECT || (feu._VECTptr->size()==2 && feu._VECTptr->at(1)==t)) {
            gen f=makevars(feu.type==_VECT?feu._VECTptr->front():feu,t,depvars,diffvars,contextptr);
            for (const_iterateur it=depvars.begin();it!=depvars.end();++it) {
                if (*it==f) return diffvars[it-depvars.begin()];
            }
        }
    } else if (e.type==_SYMB) {
        const gen &feu=e._SYMBptr->feuille;
        if (feu.type==_VECT) {
            vecteur nf;
            nf.reserve(feu._VECTptr->size());
            for (const_iterateur it=feu._VECTptr->begin();it!=feu._VECTptr->end();++it) {
                nf.push_back(makevars(*it,t,depvars,diffvars,contextptr));
            }
            return symbolic(e._SYMBptr->sommet,change_subtype(nf,_SEQ__VECT));
        }
        return symbolic(e._SYMBptr->sommet,makevars(feu,t,depvars,diffvars,contextptr));
    }
    return _eval(e,contextptr);
}

gen apply_sign(const gen &g,const gen &simp,GIAC_CONTEXT) {
    gen e=_sign(g,contextptr);
    if (!e.is_symb_of_sommet(at_sign))
        return e;
    gen arg=_apply(makesequence(simp,vecteur(1,e._SYMBptr->feuille)),contextptr)._VECTptr->front();
    arg=_factor(arg,contextptr);
    if (is_zero(simp(arg-g,contextptr)))
        return e;
    return apply_sign(arg,simp,contextptr);
}

/* return the sign of the expression, simplified */
gen determine_sign(const gen &e_orig,const gen &simp,GIAC_CONTEXT) {
    gen e=_apply(makesequence(simp,vecteur(1,e_orig)),contextptr)._VECTptr->front();
    if (e.type==_SYMB)
        return apply_sign(_factor(e,contextptr),simp,contextptr);
    return _sign(e,contextptr);
}

gen strip_sign(const gen &g) {
    if (g.is_symb_of_sommet(at_neg))
        return -strip_sign(g._SYMBptr->feuille);
    if (g.is_symb_of_sommet(at_sign))
        return g._SYMBptr->feuille;
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT) {
        gen ret(1);
        const vecteur &v=*g._SYMBptr->feuille._VECTptr;
        for (const_iterateur it=v.begin();it!=v.end();++it)
            ret=ret*strip_sign(*it);
        return ret;
    }
    return g;
}

/* return the condition(s) under which the given function is convex */
gen _convex(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    gen f=gv.front(),t(undef),TRUE=change_subtype(1,_INT_BOOLEAN),FALSE=change_subtype(0,_INT_BOOLEAN);
    if (gv.size()<2)
        return generr("Too few input arguments");
    f=idnteval(f,contextptr);
    vecteur vars;
    gen vars_arg=idnteval(gv[1],contextptr);
    if (vars_arg.type==_VECT)
        vars=*vars_arg._VECTptr;
    else vars.push_back(vars_arg);
    if (vars.empty())
        return generr("List of variables is empty");
    bool dosimp=false;
    if (gv.size()>2) {
        if (gv[2]!=at_simplify)
            return generr("Invalid optional argument");
        dosimp=true;
    }
    vecteur fvars,depvars;
    for (iterateur it=vars.begin();it!=vars.end();++it) {
        if (it->is_symb_of_sommet(at_of)) {
            const gen &u=it->_SYMBptr->feuille._VECTptr->front();
            const gen &v=it->_SYMBptr->feuille._VECTptr->back();
            if (v.type!=_IDNT || u.type!=_IDNT)
                return generr("Expected an identifier");
            if (is_undef(t))
                t=v;
            else if (t!=v)
                return generr("Expected a single parameter variable");
            if (!contains(depvars,u))
                depvars.push_back(u);
        } else if (it->type==_IDNT) {
            if (!contains(fvars,*it))
                    fvars.push_back(*it);
        } else return generr("Invalid variable specification");
    }
    vecteur diffvars,diffs;
    int cnt=0;
    for (const_iterateur it=depvars.begin();it!=depvars.end();++it) {
        diffvars.push_back(temp_symb("tmp",++cnt,contextptr));
        diffs.push_back(derive(symb_of(*it,t),t,contextptr));
    }
    gen F=is_undef(t)?_eval(f,contextptr):makevars(f,t,depvars,diffvars,contextptr);
    vecteur allvars=mergevecteur(fvars,mergevecteur(depvars,diffvars));
    gen simp_func=dosimp?at_simplify:at_ratnormal;
    if (allvars.size()==1) { // univariate function
        gen e=determine_sign(derive(F,allvars.front(),2,contextptr),simp_func,contextptr);
        if (is_one(e))
            return TRUE;
        if (is_minus_one(e))
            return FALSE;
        e=strip_sign(e);
        e=dosimp?_simplify(e,contextptr):_ratnormal(e,contextptr);
        return eval(e.is_symb_of_sommet(at_neg)?symb_inferieur_egal(e._SYMBptr->feuille,0):symb_superieur_egal(e,0),contextptr);
    }
    // multivariate case
    log_output_redirect lor(contextptr);
    int n=allvars.size();
    matrice H=*_hessian(makesequence(F,allvars),contextptr)._VECTptr;
    if (has_diff(H,allvars))
        return generr("Function is not twice differentiable");
    vector<int> perm;
    bool sing;
    int mt=is_numericm(H,num_mask_withfrac|num_mask_withint)?2:0;
    if (!ldl(H,perm,mt,sing,0,contextptr))
        return generr("Factorization failed");
    gen a11,a12,a22,simb;
    vecteur cond,res;
    cond.reserve(n);
    res.reserve(n);
    for (int i=0;i<n;++i) {
        if (i<n-1 && (!(a12=H[0][i+1]).is_integer() || a12.val!=0)) {
            // 2x2 block (all principal minors must be nonnegative!)
            cond.push_back(determine_sign(a11=H[i][i],simp_func,contextptr));
            cond.push_back(determine_sign(a22=H[i+1][i+1],simp_func,contextptr));
            cond.push_back(determine_sign(a11*a22-pow(a12,2),simp_func,contextptr));
            ++i;
        } else {
            // 1x1 block
            cond.push_back(determine_sign(H[i][i],simp_func,contextptr));
        }
    }
    for (const_iterateur it=cond.begin();it!=cond.end();++it) {
        if (is_minus_one(*it))
            return FALSE;
        if (is_one(*it) || is_zero(*it))
            continue;
        res.push_back(strip_sign(*it));
    }
    if (res.empty())
        return TRUE;
    for (int i=0;i<int(res.size());++i) {
        for (int j=res.size();j-->i+1;) {
            if (is_zero(simp(res[j]-res[i],contextptr)))
                res.erase(res.begin()+j);
        }
    }
    for (int i=res.size();i-->0;) {
        gen &r=res[i];
        r=dosimp?_simplify(r,contextptr):_ratnormal(r,contextptr);
        if (is_positive_safe(r,false,100,contextptr))
            res.erase(res.begin()+i);
        else if (is_definitely_positive(-r,contextptr))
            return FALSE;
        r=symb_superieur_egal(subst(r,diffvars,diffs,false,contextptr),0);
    }
    if (res.empty())
        return TRUE;
    if (dosimp && (simb=_lname(res,contextptr)).type==_VECT && simb._VECTptr->size()==1) {
        /* Attempt to simplify the result by solving the
         * system of inequalities depending on a single parameter. */
        vecteur s=solve_quiet(res,simb._VECTptr->front(),contextptr);
        if (!s.empty())
            res=s;
    }
    return eval(res,contextptr);
}
static const char _convex_s []="convex";
static define_unary_function_eval_quoted (__convex,&_convex,_convex_s);
define_unary_function_ptr5(at_convex,alias_at_convex,&__convex,_QUOTE_ARGUMENTS,true)

gen _numdiff(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    if (gv.size()<3)
        return generr("Too few input arguments");
    if (gv[0].type!=_VECT || gv[1].type!=_VECT)
        return generrtype("First two arguments must be lists");
    if (gv[0]._VECTptr->size()!=gv[1]._VECTptr->size() || gv[0]._VECTptr->empty())
        return generr("Lists X and Y must have equal sizes");
    vecteur X=*gv[0]._VECTptr;
    vecteur Y=*gv[1]._VECTptr;
    gen x0=gv[2];
    int M=1,N=X.size()-1;
    vecteur rest(1,M);
    if (gv.size()>=4) {
        rest=vecteur(gv.begin()+3,gv.end());
        for (const_iterateur it=rest.begin();it!=rest.end();++it) {
            if (!it->is_integer() || it->val<0)
                return generrtype("Expected a nonnegative integer");
        }
        M=_max(rest,contextptr).val;
    }
    vecteur res(rest.size());
    if (is_fully_numeric(_evalf(makevecteur(X,Y,x0),contextptr)) &&
            (is_approx(X) || is_approx(Y) || is_approx(x0))) { // the input is numeric and not fully exact
        X=*_evalf(X,contextptr)._VECTptr;
        Y=*_evalf(Y,contextptr)._VECTptr;
        x0=_evalf(x0,contextptr);
    } else { // the input is exact or symbolic
        X=*exact(X,contextptr)._VECTptr;
        Y=*exact(Y,contextptr)._VECTptr;
        x0=exact(x0,contextptr);
    }
#if 1 // Sadiq & Viswanath (2014)
    FDWeights fdw(X,M,contextptr);
    fdw.setz0(x0);
    iterateur rt=res.begin();
    for (const_iterateur it=rest.begin();it!=rest.end();++it,++rt) {
        *rt=0;
        int k=0;
        for (const_iterateur jt=Y.begin();jt!=Y.end();++jt,++k) {
            *rt+=fdw(k,it->val)*(*jt);
        }
    }
#else // Fornberg algorithm
    vecteur del((N+1)*(N+1)*(M+1),0);
    del[0]=1;
    gen c1(1),c2,c3;
    for (int n=1;n<=N;++n) {
        c2=1;
        int m_max=n<M?n:M;
        for (int ni=0;ni<n;++ni) {
            c3=X[n]-X[ni];
            c2=c2*c3;
            for (int m=0;m<=m_max;++m) {
                del[(n*(N+1)+ni)*(M+1)+m]=((X[n]-x0)*del[((n-1)*(N+1)+ni)*(M+1)+m]-
                m*del[((n-1)*(N+1)+ni)*(M+1)+(m==0?m:m-1)])/c3;
            }
        }
        for (int m=0;m<=m_max;++m) {
            del[(n*(N+1)+n)*(M+1)+m]=c1/c2*(m*del[((n-1)*(N+1)+n-1)*(M+1)+(m==0?m:m-1)]-
            (X[n-1]-x0)*del[((n-1)*(N+1)+n-1)*(M+1)+m]);
        }
        c1=c2;
    }
    // compute approximations of the derivatives
    iterateur rt=res.begin();
    for (const_iterateur it=rest.begin();it!=rest.end();++it,++rt) {
        *rt=0;
        for (int ni=0;ni<=N;++ni) {
            *rt+=del[(N*(N+1)+ni)*(M+1)+it->val]*Y[ni];
        }
    }
#endif
    gen result=res.size()==1?res.front():res;
    return result;
}
static const char _numdiff_s []="numdiff";
static define_unary_function_eval (__numdiff,&_numdiff,_numdiff_s);
define_unary_function_ptr5(at_numdiff,alias_at_numdiff,&__numdiff,0,true)

bool parse_equations(const gen &eq,vecteur &eqv) {
    if (eq.type==_VECT) {
        for (const_iterateur it=eq._VECTptr->begin();it!=eq._VECTptr->end();++it) {
            if (!parse_equations(*it,eqv))
                return false;
        }
    } else {
        gen e;
        if (is_equal(eq)) {
            const gen &lhs=eq._SYMBptr->feuille._VECTptr->front();
            const gen &rhs=eq._SYMBptr->feuille._VECTptr->back();
            e=lhs-rhs;
        } else e=eq;
        if (e.type!=_SYMB)
            return false;
        eqv.push_back(e);
    }
    return true;
}

static int integer_placeholder_count=-1;
static int positive_integer_placeholder_count=-1;

gen make_integer_placeholder(const string &base,GIAC_CONTEXT) {
    if (integer_placeholder_count<0)
        integer_placeholder_count=array_start(contextptr);
    gen ph=identificateur(string("_")+base+print_INT_(integer_placeholder_count++));
    return ph;
}

gen make_positive_integer_placeholder(const string &base,GIAC_CONTEXT) {
    if (positive_integer_placeholder_count<0)
        positive_integer_placeholder_count=array_start(contextptr);
    gen ph=identificateur(string("_")+base+print_INT_(positive_integer_placeholder_count++));
    return ph;
}

/* return true iff a is divisible by b */
bool is_divisible(const gen &a,const gen &b,GIAC_CONTEXT) {
    assert(a.is_integer() && b.is_integer());
    return is_zero(_irem(makesequence(a,b),contextptr));
}

/* if g is a homogeneous polynomial in vars of degree d>0, return d, else return 0 */
int is_homogeneous_poly(const gen &g_orig,const vecteur &vars,gen_map *coef,GIAC_CONTEXT) {
    gen g=expand(g_orig,contextptr);
    if (is_exactly_zero(g))
        return 0;
    vecteur cf(1,g),deg(1,vecteur(0));
    for (const_iterateur it=vars.begin();it!=vars.end();++it) {
        vecteur cf_next,deg_next;
        for (const_iterateur jt=cf.begin();jt!=cf.end();++jt) {
            vecteur d=*deg[jt-cf.begin()]._VECTptr;
            vecteur cfx=*_coeff(makesequence(*jt,*it),contextptr)._VECTptr;
            for (const_iterateur kt=cfx.begin();kt!=cfx.end();++kt) {
                if (is_exactly_zero(*kt))
                    continue;
                cf_next.push_back(*kt);
                d.push_back(cfx.end()-kt-1);
                deg_next.push_back(d);
                d.pop_back();
            }
        }
        cf=cf_next;
        deg=deg_next;
    }
    if (!is_constant_wrt_vars(cf,vars,contextptr))
        return 0;
    assert(!deg.empty() && ckmatrix(deg,false));
    vecteur sum_deg=*_sum(mtran(deg),contextptr)._VECTptr;
    assert(is_integer_vecteur(sum_deg,true));
    gen d=sum_deg.front();
    for (const_iterateur it=sum_deg.begin()+1;it!=sum_deg.end();++it) {
        if (!is_exactly_zero(*it-d))
            return 0;
    }
    if (coef!=NULL) {
        int n=cf.size();
        for (int i=0;i<n;++i) {
            (*coef)[deg[i]]=cf[i];
        }
    }
    return d.val;
}

/* return the root of the square part of g */
gen square_part(const gen &g,GIAC_CONTEXT) {
    assert(g.is_integer());
    vecteur f=*_ifactors(g,contextptr)._VECTptr;
    int n=f.size()/2;
    gen ret(1);
    for (int i=0;i<n;++i) {
        gen p=f[2*i],r=f[2*i+1];
        assert(r.is_integer());
        for (int j=0;j<r.val/2;++j)
            ret=ret*p;
    }
    return ret;
}

/* return true iff g = a*x^2 + b*x*y + c*y^2 + d*x + e*y + f */
bool is_binary_quadratic_wrt(const gen &g,const gen &x,const gen &y,
                             gen &a,gen &b,gen &c,gen &d,gen &e,gen &f,GIAC_CONTEXT) {
    gen A,B,C,D,E,F;
    if (!is_quadratic_wrt(g,x,A,B,C,contextptr)) {
        A=0;
        if (!is_linear_wrt(g,x,B,C,contextptr) || is_exactly_zero(B))
            return false;
    } else if (!is_constant_wrt(A,y,contextptr))
        return false;
    a=A;
    if (!is_exactly_zero(B)) {
        if (!is_linear_wrt(B,y,D,E,contextptr))
            return false;
        b=D;
        d=E;
    } else b=d=0;
    if (!is_exactly_zero(C)) {
        if (is_constant_wrt(C,y,contextptr)) {
            c=e=0;
            f=C;
        } else if (is_linear_wrt(C,y,D,E,contextptr)) {
            c=0;
            e=D;
            f=E;
        } else if (is_quadratic_wrt(C,y,D,E,F,contextptr)) {
            c=D;
            e=E;
            f=F;
        } else return false;
    } else c=e=f=0;
    return true;
}

/* convert rationals in coef to integers by multiplying them with
 * lcm of the denominators and subsequently dividing by gcd */
bool integralize(const vecteur &coef,vecteur &coef_int,GIAC_CONTEXT) {
    if (coef.empty())
        return true;
    for (const_iterateur it=coef.begin()+1;it!=coef.end();++it) {
        if (!is_rational(*it))
            return false;
    }
    coef_int=multvecteur(_lcm(_apply(makesequence(at_denom,coef),contextptr),contextptr),coef);
    assert(is_integer_vecteur(coef_int,true));
    coef_int=divvecteur(coef_int,_gcd(coef_int,contextptr));
    return is_integer_vecteur(coef_int,true);
}

/* Rosser's algorithm for an integer solution of a1 x1 + a2 x2 + ... + an xn = b,
 * return true iff it exists and in that case write it to sol
 * (source: Esmaeili et al. 2001) */
bool rosser(const vecteur &a_orig,const gen &b,vecteur &sol,GIAC_CONTEXT) {
    int n=a_orig.size();
    vecteur a=a_orig,ind(n);
    for (int i=0;i<n;++i) {
        if (is_strictly_positive(-a[i],contextptr)) {
            a[i]=-a[i];
            ind[i]=-(i+1);
        } else ind[i]=i+1;
    }
    matrice z=mtran(*_revlist(_sort(_zip(makesequence(a,ind),contextptr),contextptr),contextptr)._VECTptr);
    a=*z.front()._VECTptr;
    vector<int> iv=vecteur_2_vector_int(*z.back()._VECTptr);
    matrice C=midn(n);
    C=mtran(mergevecteur(matrice(1,a),C));
    while (!is_zero(C[1][0])) {
        gen f=_floor(fraction(C[0][0],C[1][0]),contextptr);
        C[0]=subvecteur(*C[0]._VECTptr,multvecteur(f,*C[1]._VECTptr));
        for (int i=0;i<n-1;++i) {
            if (!is_greater(C[i][0],C[i+1][0],contextptr)) {
                gen tmp=C[i+1];
                C[i+1]=C[i];
                C[i]=tmp;
            } else break;
        }
    }
    vecteur u(C[0]._VECTptr->begin()+1,C[0]._VECTptr->end());
    gen q=b/C[0][0];
    if (!q.is_integer())
        return false;
    vecteur x=multvecteur(q,u);
    sol.resize(n);
    for (int i=0;i<n;++i) {
        int j=iv[i];
        if (j<0)
            sol[-j-1]=-x[i];
        else sol[j-1]=x[i];
    }
    return true;
}

/* solve system of linear diophantine equations AX=b with integer coefficients
 * using the method of Esmaeili et al. (Numer. Math 90: 101-115, 2001),
 * return true iff there exists a solution (q is arbitrary integer vector) */
bool ABS_diophantine(const matrice &A,const vecteur &b,const vecteur &q,vecteur &sol,GIAC_CONTEXT) {
    assert(!A.empty());
    int m=A.size(),n=A.front()._VECTptr->size();
    assert(m<=n);
    vecteur x(n,0),z(n);
    matrice H=midn(n);
    for (int i=0;i<m;++i) {
        const vecteur &a=*A[i]._VECTptr;
        gen tau=scalarproduct(a,x,contextptr)-b[i];
        vecteur s=multmatvecteur(H,a);
        if (is_zero(s,contextptr)) {
            if (is_zero(tau))
                continue;
            return false;
        }
        gen delta=_gcd(s,contextptr),alpha;
        if (!(alpha=tau/delta).is_integer())
            return false;
        assert(rosser(s,delta,z,contextptr));
        vecteur p=multmatvecteur(mtran(H),z);
        x=subvecteur(x,multvecteur(alpha,p));
        H=subvecteur(H,divvecteur(mmult(mtran(matrice(1,s)),matrice(1,multvecteurmat(z,H))),delta));
    }
    H=mtran(H);
    sol=addvecteur(x,multmatvecteur(H,q));
    return true;
}

/* return true iff list of equations eqv is equivalent to linear system A*x=b */
bool is_linear_system(const vecteur &eqv,const vecteur &x,matrice &A,vecteur &b,GIAC_CONTEXT) {
    A.clear();
    b.clear();
    int n=x.size(),m=eqv.size(),i,j;
    A.reserve(m);
    b.resize(m);
    gen p,q;
    vecteur c;
    for (const_iterateur it=eqv.begin();it!=eqv.end();++it) {
        i=it-eqv.begin();
        vecteur eq(n);
        gen e=*it;
        for (const_iterateur jt=x.begin();jt!=x.end();++jt) {
            j=jt-x.begin();
            if (is_constant_wrt(e,*jt,contextptr))
                eq[j]=0;
            else if (is_linear_wrt(e,*jt,p,q,contextptr) && is_rational(p)) {
                eq[j]=p;
                e=q;
            } else return false;
        }
        if (!is_rational(e))
            return false;
        eq.push_back(-e);
        assert(integralize(eq,c,contextptr));
        b[i]=c.back();
        c.pop_back();
        A.push_back(c);
    }
    return true;
}

/* return +-1 if g is Pell equation x^2-d*y^2=+-1
 * for positive integer d, otherwise return 0 */
int is_pell_equation(const gen &g,const gen &x,const gen &y,gen &d,GIAC_CONTEXT) {
    gen a,b,c,f;
    if (!is_quadratic_wrt(g,x,a,b,c,contextptr) || !is_one(a) || !is_exactly_zero(b))
        return 0;
    if (!is_quadratic_wrt(c,y,a,b,f,contextptr) || !is_exactly_zero(b) || !is_one(_abs(f,contextptr)))
        return 0;
    d=-a;
    if (!d.is_integer() || is_positive(-d,contextptr) || is_perfect_square(d))
        return 0;
    return is_positive(f,contextptr)?-1:1;
}

/* solve Pell equation x^2-d*y^2=+-1, using the algorithm of Alexandre Junod (2015) */
vecteur solve_pell_equation(const gen &d,bool rhs_positive,int nsols,GIAC_CONTEXT) {
    if (is_perfect_square(d))
        return vecteur(0);
    int maxiter=32767,i=-1;
    gen a0,b0,c0,a1,b1,c1,a,b,c,q;
    /* initialize */
    a0=0; b0=1; c0=d;
    a=1; b=0; c=1;
    vecteur sols;
    while (++i<maxiter) {
        if (i>0) {
            /* found a solution */
            if (rhs_positive) {
                if (is_divisible(2*a,c,contextptr))
                    sols.push_back(makevecteur(2*a*a/c+gen((i+1)%2!=0?-1:1),2*a*b/c));
            } else if (i%2!=0 && is_one(c))
                sols.push_back(makevecteur(a,b));
            if (nsols==int(sols.size()))
                break;
        }
        a1=a; b1=b; c1=c;
        q=_floor((sqrt(d-c0*c,contextptr)+sqrt(d,contextptr))/c,contextptr);
        a=q*a+a0;
        b=q*b+b0;
        c=2*q*sqrt(d-c0*c,contextptr)+c0-q*q*c;
        a0=a1; b0=b1; c0=c1;
    }
    if (i>=maxiter)
        print_warning("maximum number of iterations exceeded",contextptr);
    return sols;
}

bool is_negative_pell_insoluble(const gen &d,const gen &u1,GIAC_CONTEXT) {
    if (is_one(_is_prime(d,contextptr)))
        return !is_divisible(d-1,4,contextptr);
    if (is_divisible(d-1,4,contextptr) || is_divisible(d-2,4,contextptr))
        return !is_divisible(u1+1,2*d,contextptr);
    return false;
}

gen pell_fundamental(const gen &d,bool rhs_positive,GIAC_CONTEXT) {
    if (!rhs_positive) {
        vecteur fs=solve_pell_equation(d,true,1,contextptr);
        if (fs.empty())
            return undef;
        if (is_negative_pell_insoluble(d,fs.front(),contextptr))
            return vecteur(0);
    }
    vecteur sols=solve_pell_equation(d,rhs_positive,1,contextptr);
    if (sols.empty())
        return undef;
    return sols.front();
}

gen pell_general(const gen &d,bool rhs_positive,const gen &ph,bool output_components,GIAC_CONTEXT) {
    gen fund=pell_fundamental(d,rhs_positive,contextptr);
    if (is_undef(fund) || (fund.type==_VECT && fund._VECTptr->empty()))
        return fund;
    gen n=is_undef(ph)?make_positive_integer_placeholder("N",contextptr):ph;
    const gen &u1=fund._VECTptr->front(),&v1=fund._VECTptr->back();
    gen sqr_d=sqrt(d,contextptr),f1,f2;
    if (rhs_positive) {
        f1=pow(u1+v1*sqr_d,n,contextptr); f2=pow(u1-v1*sqr_d,n,contextptr);
    } else {
        f1=pow(u1+v1*sqr_d,2*n+1,contextptr); f2=pow(u1-v1*sqr_d,2*n+1,contextptr);
    }
    if (output_components)
        return makevecteur(f1,f2,sqr_d);
    return makevecteur((f1+f2)/2,(f1-f2)/(2*sqr_d));
}

vecteur tengely_zeros(const gen &F,const gen &u,const gen &t,const gen &x,const gen &p,GIAC_CONTEXT) {
    gen s1=_zeros(makesequence(F-_pow(makesequence(u-t,p),contextptr),x),contextptr);
    gen s2=_zeros(makesequence(F-_pow(makesequence(u+t,p),contextptr),x),contextptr);
    vecteur s=mergevecteur(s1.type==_VECT?*s1._VECTptr:vecteur(0),s2.type==_VECT?*s2._VECTptr:vecteur(0));
    s.push_back(0);
    return s;
}

vecteur tengely_T_set(const gen &t,const gen &D,GIAC_CONTEXT) {
    vecteur dv=*_divis(D,contextptr)._VECTptr,res;
    for (const_iterateur it=dv.begin();it!=dv.end();++it) {
        const gen &d=*it;
        gen lb=-2*d*t,ub=2*d*t;
        gen n=lb.is_integer()?lb+1:_ceil(lb,contextptr);
        while (is_strictly_greater(ub,n,contextptr)) {
            if (!contains(res,n/d))
                res.push_back(n/d);
            n=n+1;
        }
    }
    return res;
}

/* factorize g over Q and return the list of factors with multiplicities */
vecteur factorize_over_Q(const gen &g,GIAC_CONTEXT) {
    bool old=contextptr->globalptr->_withsqrt_;
    contextptr->globalptr->_withsqrt_=false;
    vecteur fv=*_factors(g,contextptr)._VECTptr;
    assert(!fv.empty() && fv.size()%2==0);
    contextptr->globalptr->_withsqrt_=old;
    return fv;
}

/* return true if (multivariate) polynomial g with
 * rational coefficients is irreducible in Q */
bool is_irreducible_in_Q(const gen &g,GIAC_CONTEXT) {
    vecteur fv=factorize_over_Q(g,contextptr);
    int cnt=0;
    for (const_iterateur it=fv.begin();it!=fv.end();it=it+2) {
        if (!it->is_integer())
            ++cnt;
    }
    return cnt==1;
}

bool is_tengely(const gen &g,const gen &x,const gen &y,gen &F,gen &G,GIAC_CONTEXT) {
    vecteur xcf=*_coeff(makesequence(expand(g,contextptr),x),contextptr)._VECTptr;
    if (xcf.size()<2)
        return false;
    gen rest=xcf.back();
    xcf.pop_back();
    vecteur ycf=*_coeff(makesequence(expand(rest,contextptr),y),contextptr)._VECTptr;
    if (ycf.size()<2)
        return false;
    vecteur v=mergevecteur(xcf,ycf),vi;
    if (!is_constant_wrt_vars(v,makevecteur(x,y),contextptr) || !integralize(v,vi,contextptr))
        return false;
    xcf=vecteur(vi.begin(),vi.begin()+xcf.size());
    ycf=vecteur(vi.begin()+xcf.size(),vi.end());
    xcf.push_back(0);
    if (!is_one(-xcf.front()*ycf.front()))
        return false;
    F=_r2e(makesequence(xcf,x),contextptr);
    G=_r2e(makesequence(ycf,y),contextptr);
    if (is_one(_lcoeff(makesequence(F,x),contextptr)))
        G=-G;
    else F=-F;
    if (!is_irreducible_in_Q(F-G,contextptr))
        return false;
    gen n=xcf.size()-1,m=ycf.size()-1;
    return is_greater(m,n,contextptr) && !is_one(_gcd(makesequence(n,m),contextptr));
}

/* solve Diophantine equation in form F(x)=G(y) with gcd(deg(G),deg(F))>1,
 * where F and G are monic polynomials such that F(X)-G(Y) is irreducible
 * in Q[X,Y] and also deg(G)>=deg(F), using Tengely's algorithm */
vecteur solve_tengely(const gen &F,const gen &G,const gen &x,const gen &y,GIAC_CONTEXT) {
    gen n=_degree(makesequence(F,x),contextptr);
    gen m=_degree(makesequence(G,y),contextptr);
    gen d=_gcd(makesequence(n,m),contextptr);
    assert(d.is_integer() && d.val>1);
    gen p=_ifactors(d,contextptr)._VECTptr->front();
    gen pol=change_subtype(gen(_POLY1__VECT),_INT_MAPLECONVERSION);
    gen Fp=_pow(makesequence(F,_inv(p,contextptr)),contextptr);
    gen Gp=_pow(makesequence(G,_inv(p,contextptr)),contextptr);
    gen u=_series(makesequence(Fp,symb_equal(x,plus_inf),n/p,pol),contextptr);
    gen v=_series(makesequence(Gp,symb_equal(y,plus_inf),m/p,pol),contextptr);
    u=_ratnormal(u,contextptr);
    v=_ratnormal(v,contextptr);
    vecteur cu=*_coeff(makesequence(u,x),contextptr)._VECTptr;
    vecteur cv=*_coeff(makesequence(v,y),contextptr)._VECTptr;
    gen c0=cu.back()-cv.back();
    cu.pop_back();
    cv.pop_back();
    vecteur cf=mergevecteur(cu,cv);
    cf.push_back(c0);
    gen D=_lcm(_apply(makesequence(at_denom,cf),contextptr),contextptr);
    gen xlo,xhi,ylo,yhi,xl,xh,yl,yh,t=_inv(2*D,contextptr),numeq,last_numeq(undef);
    while (true) {
        vecteur xs=tengely_zeros(F,u,t,x,p,contextptr);
        vecteur ys=tengely_zeros(G,v,t,y,p,contextptr);
        xlo=_min(xs,contextptr);
        xhi=_max(xs,contextptr);
        ylo=_min(ys,contextptr);
        yhi=_max(ys,contextptr);
        numeq=_floor(xhi-xlo,contextptr)+_floor(yhi-ylo,contextptr)+(p==2?8*D*t:4*D*t+1);
        if (!is_undef(last_numeq) && !is_strictly_greater(last_numeq,numeq,contextptr))
            break;
        last_numeq=numeq;
        xl=_ceil(xlo,contextptr); if (!xl.is_integer()) xl=_ceil(_evalf(xlo,contextptr),contextptr);
        xh=_floor(xhi,contextptr); if (!xh.is_integer()) xh=_ceil(_evalf(xhi,contextptr),contextptr);
        yl=_ceil(ylo,contextptr); if (!yl.is_integer()) yl=_ceil(_evalf(ylo,contextptr),contextptr);
        yh=_floor(yhi,contextptr); if (!yh.is_integer()) yh=_ceil(_evalf(yhi,contextptr),contextptr);
        t=2*t;
    }
    vecteur sols;
    gen k,sol,soly,s,res;
    for (k=yl;is_greater(yh,k,contextptr);k=k+1) {
        sol=_rationalroot(subst(G,y,k,false,contextptr)-F,contextptr);
        if (sol.type==_VECT) for (const_iterateur it=sol._VECTptr->begin();it!=sol._VECTptr->end();++it) {
            if (it->is_integer() && !contains(sols,s=makevecteur(*it,k)))
                sols.push_back(s);
        }
    }
    for (k=xl;is_greater(xh,k,contextptr);k=k+1) {
        sol=_rationalroot(subst(F,x,k,false,contextptr)-G,contextptr);
        if (sol.type==_VECT) for (const_iterateur it=sol._VECTptr->begin();it!=sol._VECTptr->end();++it) {
            if (it->is_integer() && !contains(sols,s=makevecteur(k,*it)))
                sols.push_back(s);
        }
    }
    vecteur T=tengely_T_set(t,D,contextptr);
    for (int i=0;i<2;++i) {
        for (const_iterateur it=T.begin();it!=T.end();++it) {
            res=_resultant(makesequence(F-G,u+(i==0?-v:v)-*it,y),contextptr);
            sol=_rationalroot(res,contextptr);
            if (sol.type==_VECT) for (const_iterateur jt=sol._VECTptr->begin();jt!=sol._VECTptr->end();++jt) {
                if (jt->is_integer()) {
                    soly=_rationalroot(subst(F,x,*jt,false,contextptr)-G,contextptr);
                    if (soly.type==_VECT) for (const_iterateur kt=soly._VECTptr->begin();kt!=soly._VECTptr->end();++kt) {
                        if (kt->is_integer() && !contains(sols,s=makevecteur(*jt,*kt)))
                            sols.push_back(s);
                    }
                }
            }
        }
        if (is_strictly_greater(p,2,contextptr))
            break;
    }
    return sols;
}

/* return true iff g is a polynomial with integer coefficients w.r.t. vars */
bool is_int_poly(const gen &g,const vecteur &vars,GIAC_CONTEXT) {
    if (vars.empty() || is_constant_wrt_vars(g,vars,contextptr))
        return g.is_integer();
    if (!is_rational_wrt_vars(g,vars) || !is_one(_denom(g,contextptr)))
        return false;
    vecteur cf=*_coeff(makesequence(expand(g,contextptr),vars.front()),contextptr)._VECTptr;
    vecteur rest_vars(vars.begin()+1,vars.end());
    for (const_iterateur it=cf.begin();it!=cf.end();++it) {
        if (it->is_integer())
            continue;
        if (!is_int_poly(*it,rest_vars,contextptr))
            return false;
    }
    return true;
}

/* return true iff g is a polynomial in Q[x] with either:
 * - at least two distinct irreducible factors (over Q),
 * - exactly one irreducible factor with degree > 2,
 * - two (conjugate) complex roots. */
bool is_thue_polynomial(const gen &g,const gen &x,GIAC_CONTEXT) {
    vecteur cf=*_coeff(makesequence(expand(g,contextptr),x),contextptr)._VECTptr;
    for (const_iterateur it=cf.begin();it!=cf.end();++it) {
        if (!is_rational(*it))
            return false;
    }
    vecteur fv=factorize_over_Q(g,contextptr);
    int cnt=0;
    for (const_iterateur it=fv.begin();it!=fv.end();it+=2) {
        if (!it->is_integer())
            ++cnt;
    }
    if (cnt>=2 || (cnt==1 && is_greater(_degree(makesequence(g,x),contextptr),3,contextptr)))
        return true;
    gen z=_cZeros(makesequence(g,x),contextptr);
    if (z.type==_VECT) for (const_iterateur it=z._VECTptr->begin();it!=z._VECTptr->end();++it) {
        if (!is_exactly_zero(im(*it,contextptr)))
            return true;
    }
    return false;
}

/* return true iff g=P(X,Y)-a=0 is a Thue equation in x and y with finitely many solutions,
 * provides univariate polynomial P(x,1) and the free term a (nonzero integer) */
bool is_thue_equation(const gen &g,const gen &x,const gen &y,gen &P,gen &a,GIAC_CONTEXT) {
    vecteur vars=makevecteur(x,y);
    a=-_ratnormal(subst(g,vars,vecteur(vars.size(),0),false,contextptr),contextptr);
    if (!is_rational(a) || is_homogeneous_poly(g+a,vars,NULL,contextptr)==0)
        return false;
    P=subst(g+a,y,1,false,contextptr);
    return is_thue_polynomial(P,x,contextptr);
}

/* return true iff x^2+A=0 (mod B) is soluble */
bool is_soluble(const gen &A,const gen &B,GIAC_CONTEXT) {
    gen x(0);
    for (;is_strictly_greater(B,x,contextptr);x=x+1) {
        if (is_divisible(x*x+A,B,contextptr))
            return true;
    }
    return false;
}

/* return true iff a, b, and c do not have the same sign,
 * ax^2+by^2+cz^2=0 (mod 8) has a non-trivial solution (modulo 8) and
 * x^2+bc=0 (mod a), x^2+ac=0 (mod b), and x^2+ab=0 (mod c) are also soluble,
 * where a, b, and c are square-free and pairwise coprime
 * (this is equivalent to ax^2+by^2+cz^2=0 being locally soluble everywhere) */
bool is_locally_soluble_everywhere(const gen &a,const gen &b,const gen &c,GIAC_CONTEXT) {
    if (is_strictly_positive(a*b, contextptr) &&
        is_strictly_positive(a*c, contextptr) &&
        is_strictly_positive(b*c, contextptr))
        return false;
    bool is_soluble_in_Q2=false;
    for (int x=0;!is_soluble_in_Q2 && x<8;++x) {
        for (int y=0;!is_soluble_in_Q2 && y<8;++y) {
            for (int z=0;z<8;++z) {
                gen X(x),Y(y),Z(z);
                if (is_divisible(a*X*X+b*Y*Y+c*Z*Z,8,contextptr)) {
                    is_soluble_in_Q2=true;
                    break;
                }
            }
        }
    }
    if (!is_soluble_in_Q2)
        return false;
    return is_soluble(b*c,abs(a,contextptr),contextptr) &&
            is_soluble(a*c,abs(b,contextptr),contextptr) &&
            is_soluble(a*b,abs(c,contextptr),contextptr);
}

/* use Lagrange's method of descent to determine a solution of w^2=Ax^2+By^2 which is
 * locally soluble everywhere, return [x0,y0,w0] */
vecteur descent(const gen &A,const gen &B,GIAC_CONTEXT) {
    if (is_strictly_greater(abs(A,contextptr),abs(B,contextptr),contextptr)) {
        vecteur sol=descent(B,A,contextptr);
        return makevecteur(sol[1],sol[0],sol[2]);
    }
    if (is_one(A))
        return makevecteur(1,0,1);
    if (is_one(B))
        return makevecteur(0,1,1);
    gen ub=abs(B,contextptr)/gen(2),r(0);
    for (;is_greater(ub,r,contextptr);r=r+1) {
        if (is_divisible(r*r-A,B,contextptr))
            break;
    }
    assert(is_greater(ub,r,contextptr));
    gen Q=_ratnormal((r*r-A)/B,contextptr);
    gen d=square_part(Q,contextptr),d2=d*d;
    vecteur sol=descent(A,_ratnormal(Q/d2,contextptr),contextptr);
    return makevecteur((r*sol[0]-sol[2])/(r*r-A),sol[1]/(B*d),(-A*sol[0]+r*sol[2])/(r*r-A));
}

/* convert a ternary quadratic form into a*x^2+b*y^2+c*z^2 where a, b and c are square-free and pairwise coprime */
bool ternary_quadratic_normal_form(const gen &g_orig,const vecteur &vars,vecteur &P,vecteur &q,GIAC_CONTEXT) {
    /* diagonalize the form */
    gen g=_gauss(makesequence(g_orig,vars,0),contextptr);
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()!=3 ||
        g._VECTptr->at(1).type!=_VECT || g._VECTptr->at(2).type!=_VECT)
        return false;
    vecteur cf=*g._VECTptr->at(1)._VECTptr;
    P=*g._VECTptr->at(2)._VECTptr;
    if (!integralize(cf,q,contextptr))
        return false;
    /* make the coefficients square-free */
    for (iterateur it=q.begin();it!=q.end();++it) {
        gen d=square_part(*it,contextptr),&p=P[it-q.begin()];
        p=multvecteur(d,*p._VECTptr);
        *it=*it/(d*d);
    }
    /* make the coefficients pairwise coprime */
    for (int c=0;c<3;++c) {
        int i=c,j=(c+1)%3,k=(c+2)%3;
        gen d=_gcd(makesequence(q[i],q[j]),contextptr);
        q[i]=q[i]/d;
        q[j]=q[j]/d;
        q[k]=q[k]*d;
        P[k]=divvecteur(*P[k]._VECTptr,d);
    }
    return true;
}

/* return the general solution of Q(x,y,z)=0, where Q is a ternary quadratic form */
gen solve_ternary_quadratic_form(const gen &Q,const vecteur &vars,const vecteur &fvars,bool all_sols,GIAC_CONTEXT) {
    assert(fvars.size()>=2);
    const gen &p=fvars[0],&q=fvars[1],&Z3=fvars[2];
    /* normalize the quadratic form Q */
    vecteur R,P;
    if (!ternary_quadratic_normal_form(Q,vars,P,R,contextptr)) {
        print_error("failed to normalize quadratic form",contextptr);
        return undef;
    }
    const gen &a=R[0],&b=R[1],&c=R[2];
    /* check local solubility */
    if (!is_locally_soluble_everywhere(a,b,c,contextptr))
        return vecteur(0);
    /* find a particular integral solution */
    vecteur s0=descent(-a*c,-b*c,contextptr),sol;
    s0[2]=s0[2]/c; //get z from w
    s0=*_linsolve(makesequence(P,s0),contextptr)._VECTptr;
    if (!integralize(s0,sol,contextptr)) {
        print_error("particular solution is not rational",contextptr);
        return undef;
    }
    if (all_sols) {
        /* find the general solution */
        gen r=temp_symb("r",-1,contextptr),A,B;
        gen pol=_ratnormal(subst(Q,vars,makevecteur(r*sol[0],r*sol[1]+p,r*sol[2]+q),false,contextptr),contextptr);
        assert(is_linear_wrt(pol,r,A,B,contextptr));
        vecteur cr(5),ci;
        cr[0]=-_coeff(makesequence(A,p,1),contextptr);
        cr[1]=-_coeff(makesequence(A,q,1),contextptr);
        cr[2]=_coeff(makesequence(B,p,2),contextptr);
        cr[3]=_coeff(makesequence(B,makevecteur(p,q),makevecteur(1,1)),contextptr);
        cr[4]=_coeff(makesequence(B,q,2),contextptr);
        assert(integralize(cr,ci,contextptr));
        vecteur gsol(3);
        gen p2=pow(p,2,contextptr),q2=pow(q,2,contextptr);
        gsol[0]=sol[0]*ci[2]*p2+sol[0]*ci[3]*p*q+sol[0]*ci[4]*q2;
        gsol[1]=(ci[2]*sol[1]+ci[0])*p2+(ci[3]*sol[1]+ci[1])*p*q+ci[4]*sol[1]*q2;
        gsol[2]=ci[2]*sol[2]*p2+(ci[3]*sol[2]+ci[0])*p*q+(ci[4]*sol[2]+ci[1])*q2;
        return divvecteur(multvecteur(Z3,gsol),_Gcd(change_subtype(gsol,_SEQ__VECT),contextptr));
    }
    return sol;
}

/* remove solutions of x^2+bxy+cy^2+f=0 associated with the other ones
 * by using the criterion of Andreescu (2015, p. 56) */
void purge_nonfundamental_solutions(vecteur &sols,const gen &a,const gen &b,const gen &c,const gen &f,GIAC_CONTEXT) {
    gen D=b*b-4*a*c,N=-4*a*f;
    for (int i=0;i<int(sols.size());++i) {
        const gen &x1=sols[i]._VECTptr->front(),&y1=sols[i]._VECTptr->back();
        gen u1=2*a*x1+b*y1,v1=y1;
        for (int j=sols.size();j-->i+1;) {
            const gen &x2=sols[j]._VECTptr->front(),&y2=sols[j]._VECTptr->back();
            gen u2=2*a*x2+b*y2,v2=y2;
            if (is_divisible(u1*u2-D*v1*v2,N,contextptr) && is_divisible(v1*u2-u1*v2,N,contextptr)) {
                if (is_positive(v1,contextptr) && (is_strictly_positive(-v2,contextptr) || is_greater(v2,v1,contextptr)))
                    sols.erase(sols.begin()+j);
                else {
                    sols.erase(sols.begin()+i);
                    --i;
                    break;
                }
            }
        }
    }
}

/* find continued fraction expansion of a quadratic irrationality g, double the period if odd-length */
vecteur contfrac_expansion(const gen &g,GIAC_CONTEXT) {
    assert(!is_rational(g));
    int k=0;
    vecteur ex,period;
    while (true) {
        ex=*_dfc(makesequence(g,pow(10,++k,contextptr)),contextptr)._VECTptr;
        if (ex.back().type==_VECT) {
            period=*ex.back()._VECTptr;
            ex.pop_back();
            break;
        } // note that the loop will go forever if g is not a quadratic irrationality!
    }
    if (period.size()%2==0)
        return ex;
    return mergevecteur(ex,period);
}

/* find solutions to ax^2+bxy+cy^2+f=0 where 4f^2<b^2-4ac and b^2-4ac>0 is not a perfect square */
vecteur solve_contfrac(const gen &a,const gen &b,const gen &c,const gen &f,bool only_first,GIAC_CONTEXT) {
    gen d=b*b-4*a*c,sqr_d=sqrt(d,contextptr);
    assert(is_strictly_positive(d,contextptr) && !sqr_d.is_integer() && is_strictly_greater(d,4*f*f,contextptr));
    vecteur t=makevecteur((-b+sqr_d)/(2*a),(-b-sqr_d)/(2*a)),sols;
    for (const_iterateur it=t.begin();it!=t.end();++it) {
        vecteur ex=contfrac_expansion(*it,contextptr),p,q;
        int n=ex.size();
        p.resize(n+2,0); q.resize(n+2,0);
        p[1]=q[0]=1;
        for (int k=0;k<n;++k) {
            gen x=ex[k]*p[k+1]+p[k],y=ex[k]*q[k+1]+q[k],g=gcd(x,y);
            x=x/g; y=y/g;
            if (is_zero(a*x*x+b*x*y+c*y*y+f)) {
                if (only_first)
                    return makevecteur(x,y);
                sols.push_back(makevecteur(x,y));
            }
            p[k+2]=x;
            q[k+2]=y;
        }
    }
    return sols;
}

/* find solutions to ax^2+bxy+cy^2+f=0 where b^2-4ac>0 is not a perfect square, 4f^2>=b^2-4ac and gcd(a,b,f)=1 */
vecteur solve_abf_coprime(const gen &a,const gen &b,const gen &c,const gen &f,GIAC_CONTEXT) {
    vecteur sols;
    for (gen s=0;is_strictly_greater(_abs(f,contextptr),s,contextptr);s=s+1) {
        if (is_divisible(a*s*s+b*s+c,f,contextptr)) {
            vecteur sol=solve_contfrac((a*s*s+b*s+c)/f,-2*a*s-b,a*f,1,false,contextptr);
            for (const_iterateur it=sol.begin();it!=sol.end();++it) {
                const vecteur yz=*(it->_VECTptr);
                gen x=s*yz[0]-f*yz[1],y=yz[0];
                sols.push_back(makevecteur(x,y));
            }
        }
    }
    return mergevecteur(sols,multvecteur(-1,sols));
}

/* find solutions to ax^2+bxy+cy^2+f=0 where b^2-4ac>0 is not a perfect square */
vecteur solve_homogeneous_binary_quadratic(const vecteur &cf,GIAC_CONTEXT) {
    gen a=cf[0],b=cf[1],c=cf[2],f=cf[3];
    gen d=b*b-4*a*c,g=_gcd(makesequence(a,b,c),contextptr);
    if (!is_divisible(f,g,contextptr))
        return vecteur(0); // no solutions
    a=a/g; b=b/g; c=c/g; f=f/g;
    if (is_strictly_greater(d,4*f*f,contextptr))
        return solve_contfrac(a,b,c,f,false,contextptr);
    vecteur dv=*_idivis(f,contextptr)._VECTptr,sols,sol;
    for (const_iterateur it=dv.begin();it!=dv.end();++it) {
        if (!is_perfect_square(*it))
            continue;
        gen F=f/(*it),i,m;
        if (is_strictly_greater(d,4*F*F,contextptr)) {
            sol=solve_contfrac(a,b,c,F,false,contextptr);
        } else if (is_one(_gcd(makesequence(a,b,F),contextptr))) {
            sol=solve_abf_coprime(a,b,c,F,contextptr);
        } else if (is_one(_gcd(makesequence(b,c,F),contextptr))) {
            sol=solve_abf_coprime(c,b,a,F,contextptr);
            sol=*_apply(makesequence(at_revlist,sol),contextptr)._VECTptr;
        } else {
            bool found=false;
            for (i=0;!found && is_strictly_greater(_abs(F,contextptr),i,contextptr);i=i+1) {
                for (m=0;is_greater(i+1,m,contextptr);m=m+1) {
                    if (is_one(gcd(i,m)) && is_one(gcd(a*i*i+b*i*m+c*m*m,F))) {
                        found=true;
                        break;
                    }
                }
            }
            if (found) {
                vecteur nj=*_iabcuv(makesequence(i,-m,1),contextptr)._VECTptr;
                const gen &n=nj.front(),&j=nj.back();
                sol=solve_abf_coprime(a*i*i+b*i*m+c*m*m,2*a*i*j+b*i*n+b*j*m+2*c*m*n,a*j*j+b*j*n+c*n*n,F,contextptr);
                for (iterateur jt=sol.begin();jt!=sol.end();++jt) {
                    gen X=jt->_VECTptr->front(),Y=jt->_VECTptr->back();
                    *jt=makevecteur(i*X+j*Y,m*X+n*Y);
                }
            }
        }
        if (!sol.empty())
            sols=mergevecteur(sols,multvecteur(sqrt(*it,contextptr),sol));
    }
    return sols;
}

vecteur solve_linrec(const vecteur &r1,const vecteur &r2,const gen &x1,const gen &y1,const gen &n,GIAC_CONTEXT) {
    gen u=temp_symb("u",-1,contextptr),v=temp_symb("v",-1,contextptr);
    gen un=symb_of(u,n),vn=symb_of(v,n),unp=symb_of(u,n+1),vnp=symb_of(v,n+1),u1=symb_of(u,0),v1=symb_of(v,0);
    gen rec1=symb_equal(unp,r1[0]*un+r1[1]*vn+r1[2]),rec2=symb_equal(vnp,r2[0]*un+r2[1]*vn+r2[2]);
    gen icond1=symb_equal(u1,x1),icond2=symb_equal(v1,y1);
    gen args=makesequence(makevecteur(rec1,rec2),makevecteur(un,vn),makevecteur(icond1,icond2));
    return *_rsolve(args,contextptr)._VECTptr->front()._VECTptr;
}

/* solve general binary quadratic equation using the methods of Dario Alejandro Alpern at
 * https://www.alpertron.com.ar/METHODS.HTM */
gen solve_binary_quadratic(const vecteur &cf,gen &ph,bool set_ph,bool &alt_sols,bool fund_only,GIAC_CONTEXT) {
    gen a=cf[0],b=cf[1],c=cf[2],d=cf[3],e=cf[4],f=cf[5];
    if (is_zero(a) && is_zero(b) && is_zero(c)) {
        /* linear case */
        gen g=_gcd(makesequence(d,e),contextptr);
        if (!is_divisible(f,g,contextptr))
            return vecteur(0); // no solutions
        d=d/g; e=e/g; f=f/g;
        gen uv=_iabcuv(makesequence(d,e,-f),contextptr);
        assert(uv.type==_VECT && uv._VECTptr->size()==2);
        gen u=uv._VECTptr->front(),v=uv._VECTptr->back();
        gen t=is_undef(ph)?make_integer_placeholder("Z",contextptr):ph;
        if (is_undef(ph) && set_ph)
            ph=t;
        return makevecteur(u+e*t,v-d*t);
    }
    vecteur sol;
    if (is_zero(a) && is_zero(c) && !is_zero(b)) {
        /* simple hyperbolic case */
        if (is_zero(d*e-b*f)) {
            bool c1=is_divisible(e,b,contextptr),c2=is_divisible(d,b,contextptr);
            if (!c1 && !c2)
                return vecteur(0); // no solutions
            gen t=is_undef(ph)?make_integer_placeholder("Z",contextptr):ph;
            if (c1) sol.push_back(makevecteur(-e/b,t));
            if (c2) sol.push_back(makevecteur(t,-d/b));
        } else {
            vecteur dv=*_idivis(d*e-b*f,contextptr)._VECTptr;
            dv=mergevecteur(dv,multvecteur(-1,dv)); // also include negative divisors
            for (const_iterateur it=dv.begin();it!=dv.end();++it) {
                const gen &D=*it;
                gen x=(D-e)/b,y=((d*e-b*f)/D-d)/b;
                if (x.is_integer() && y.is_integer())
                    sol.push_back(makevecteur(x,y));
            }
        }
        return sol;
    }
    gen dsc=b*b-4*a*c;
    if (is_positive(dsc,contextptr)) {
        /* elliptical case */
        gen A=dsc,B=2*(b*e-2*c*d),C=e*e-4*c*f;
        dsc=B*B-4*A*C;
        if (!is_positive(dsc,contextptr))
            return vecteur(0); // no solutions
        if (!is_positive(A,contextptr)) {
            A=-A; B=-B; C=-C;
        }
        gen r1=(-B-sqrt(dsc,contextptr))/(2*A),r2=(-B+sqrt(dsc,contextptr))/(2*A);
        assert(is_greater(r2,r1,contextptr));
        A=c;
        for (gen x=_ceil(r1,contextptr);is_greater(r2,x,contextptr);x=x+1) {
            B=e+b*x; C=a*x*x+d*x+f;
            if (!is_positive(dsc=B*B-4*A*C,contextptr))
                continue;
            gen y1=(-B-sqrt(dsc,contextptr))/(2*A),y2=(-B+sqrt(dsc,contextptr))/(2*A);
            if (y1.is_integer())
                sol.push_back(makevecteur(x,y1));
            if (!is_zero(dsc) && y2.is_integer())
                sol.push_back(makevecteur(x,y2));
        }
        return sol;
    }
    if (is_zero(dsc)) {
        /* parabolic case */
        gen sgn(is_positive(a*b,contextptr)?1:-1);
        gen g=(is_positive(a,contextptr)?1:-1)*_gcd(makesequence(a,c),contextptr);
        a=a/g; b=b/g; c=c/g;
        gen sqr_a=sqrt(a,contextptr),sqr_c=sgn*sqrt(c,contextptr);
        assert(sqr_a.is_integer() && sqr_c.is_integer());
        gen r=sqr_c*d-sqr_a*e;
        if (is_zero(r)) {
            dsc=d*d-4*a*f*g;
            if (!is_positive(dsc,contextptr))
                return vecteur(0); // no solutions
            gen u1=(-d+sqrt(dsc,contextptr))/(2*sqr_a*g),u2=(-d-sqrt(dsc,contextptr))/(2*sqr_a*g);
            if (u1.is_integer()) {
                gen s=solve_binary_quadratic(makevecteur(0,0,0,sqr_a,sqr_c,-u1),ph,true,alt_sols,false,contextptr);
                if (s.type==_VECT && !s._VECTptr->empty())
                    sol.push_back(s);
            }
            if (!is_zero(dsc) && u2.is_integer()) {
                gen s=solve_binary_quadratic(makevecteur(0,0,0,sqr_a,sqr_c,-u2),ph,false,alt_sols,false,contextptr);
                if (s.type==_VECT && !s._VECTptr->empty())
                    sol.push_back(s);
            }
        } else {
            gen t(undef);
            for (gen u(0);is_greater(_abs(r,contextptr),u,contextptr);u=u+1) {
                gen numx=sqr_c*g*u*u+e*u+sqr_c*f,numy=sqr_a*g*u*u+d*u+sqr_a*f;
                if (is_divisible(numx,r,contextptr) && is_divisible(numy,r,contextptr)) {
                    if (is_undef(t))
                        t=is_undef(ph)?make_integer_placeholder("Z",contextptr):ph;
                    gen x=-sqr_c*g*r*pow(t,2,contextptr)-t*(e+2*sqr_c*g*u)-numx/r;
                    gen y=sqr_a*g*r*pow(t,2,contextptr)+t*(d+2*sqr_a*g*u)+numy/r,contextptr;
                    sol.push_back(makevecteur(x,y));
                }
            }
        }
        return sol;
    }
    if (is_strictly_positive(dsc,contextptr)) {
        /* tedious hyperbolic case */
        if (is_zero(d) && is_zero(e)) { // binary quadratic form
            if (!is_perfect_square(dsc)) {
                vecteur psols=solve_homogeneous_binary_quadratic(makevecteur(a,b,c,f),contextptr);
                if (!psols.empty()) {
                    purge_nonfundamental_solutions(psols,a,b,c,f,contextptr);
                    gen n;
                    if (!fund_only)
                        n=is_undef(ph)?make_integer_placeholder("Z",contextptr):ph;
                    for (const_iterateur it=psols.begin();it!=psols.end();++it) {
                        gen x=it->_VECTptr->front(),y=it->_VECTptr->back();
                        if (!fund_only) {
                            /* find other solutions by solving a system of linear recurrences */
                            vecteur rs=solve_contfrac(1,b,a*c,-1,true,contextptr);
                            const gen &r=rs.front(),&s=rs.back();
                            vecteur rsol=solve_linrec(makevecteur(r,-c*s,0),makevecteur(a*s,r+b*s,0),x,y,n,contextptr);
                            x=rsol.front(); y=rsol.back();
                        }
                        sol.push_back(makevecteur(x,y));
                    }
                }
            } else { // simple factorization
                if (is_zero(a)) {
                    gen s=solve_binary_quadratic(makevecteur(c,b,a,f),ph,false,alt_sols,false,contextptr);
                    if (s.type==_VECT)
                        sol=*_apply(makesequence(at_revlist,s),contextptr)._VECTptr;
                } else {
                    gen k=sqrt(dsc,contextptr);
                    if (is_zero(f)) {
                        gen s=solve_binary_quadratic(makevecteur(0,0,0,2*a,b+k,0),ph,true,alt_sols,false,contextptr);
                        if (s.type==_VECT && !s._VECTptr->empty())
                            sol.push_back(s);
                        s=solve_binary_quadratic(makevecteur(0,0,0,2*a,b-k,0),ph,false,alt_sols,false,contextptr);
                        if (s.type==_VECT && !s._VECTptr->empty())
                            sol.push_back(s);
                    } else {
                        vecteur u=*_idivis(-4*a*f,contextptr)._VECTptr;
                        u=mergevecteur(u,multvecteur(-1,u)); // also include negative divisors
                        for (const_iterateur it=u.begin();it!=u.end();++it) {
                            gen y=(*it+4*a*f/(*it))/(2*k);
                            gen x=(*it-(b+k)*y)/(2*a);
                            if (x.is_integer() && y.is_integer())
                                sol.push_back(makevecteur(x,y));
                        }
                    }
                }
            }
        } else { // the non-homogeneous case
            gen G=2*a*e-b*d,s,g;
            g=_gcd(makesequence(dsc,G),contextptr);
            gen F=4*a*(-f*dsc-a*e*e+b*d*e-c*d*d)/g;
            s=solve_binary_quadratic(makevecteur(-dsc/g,0,g,0,0,F),ph,false,alt_sols,true,contextptr);
            if (s.type==_VECT) {
                for (const_iterateur it=s._VECTptr->begin();it!=s._VECTptr->end();++it) {
                    const gen &x1=it->_VECTptr->front(),&y1=it->_VECTptr->back();
                    gen y=-(g*y1-G)/dsc;
                    gen x=(x1-b*y-d)/(2*a);
                    vecteur Z=*_lname(makevecteur(x,y),contextptr)._VECTptr;
                    if (Z.empty()) {
                        if (x.is_integer() && y.is_integer())
                            sol.push_back(makevecteur(x,y));
                    } else {
                        assert(Z.size()==1); // there should be exactly one global variable
                        const gen &t=Z.front();
                        gen A,B,C,D;
                        assert(is_linear_wrt(x,t,A,B,contextptr) && is_linear_wrt(y,t,C,D,contextptr));
                        if (!B.is_integer() || !D.is_integer())
                            continue;
                        gen pA=_numer(A,contextptr),pC=_numer(C,contextptr);
                        gen qA=_denom(A,contextptr),qC=_denom(C,contextptr);
                        A=pA/_gcd(makesequence(pA,qA),contextptr);
                        C=pC/_gcd(makesequence(pC,qC),contextptr);
                        sol.push_back(makevecteur(A*t+B,C*t+D));
                    }
                }
                if (!fund_only && !is_perfect_square(dsc) && !sol.empty()) {
                    /* find other solutions by solving a system of linear recurrences */
                    gen n=is_undef(ph)?make_integer_placeholder("Z",contextptr):ph;
                    for (iterateur it=sol.begin();it!=sol.end();++it) {
                        gen x=it->_VECTptr->front(),y=it->_VECTptr->back();
                        vecteur rsv=solve_contfrac(1,b,a*c,-1,false,contextptr);
                        gen P,Q,R,S,K,L,r,s;
                        for (const_iterateur jt=rsv.begin();jt!=rsv.end();++jt) {
                            r=jt->_VECTptr->front(); s=jt->_VECTptr->back();
                            P=r,Q=-c*s; R=a*s; S=r+b*s;
                            K=-(c*d*(P+S-2)+e*(b-b*r-2*a*c*s))/dsc;
                            L=d*s-(d*(b-b*r-2*a*c*s)+a*e*(P+S-2))/dsc;
                            if (K.is_integer() && L.is_integer())
                                break;
                        }
                        if (!K.is_integer() || !L.is_integer()) {
                            r=rsv.front()._VECTptr->front(); s=rsv.front()._VECTptr->back();
                            P=r*r-a*c*s*s; Q=-c*s*(2*r+b*s); K=-c*d*s*s-e*r*s;
                            R=a*s*(2*r+b*s); S=r*r+2*b*r*s+(b*b-a*c)*s*s; L=d*s*(r+b*s)-a*e*s*s;
                        }
                        *it=solve_linrec(makevecteur(P,Q,K),makevecteur(R,S,L),x,y,n,contextptr);
                    }
                }
            }
        }
        return sol;
    }
    return undef; // unable to solve the equation
}

/* Diophantine equation solver, which can solve:
 *  - (systems of) linear equation(s) (Esmaeili et al., 2001)
 *  - Pell equations (Junod, 2015)TRUE
 *  - Thue equations (using PARI) - only in MinGW!
 *  - homogeneous equations Q(x,y,z)=0, where Q is ternary quadratic form (Smart, 1998)
 *  - general quadratic equations in two variables (Alpern)
 *  - certain polynomial equations of the type f(x)=g(y) (Tengely, 2003)
 */
gen _isolve(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    gen eq;
    vecteur vars,fvars,opts,eqv;
    if (g.type==_VECT && g.subtype==_SEQ__VECT) {
        const vecteur &gv=*g._VECTptr;
        eq=gv.front();
        opts=vecteur(gv.begin()+1,gv.end());
    } else eq=g;
    if (!parse_equations(eq,eqv))
        return generrtype("Failed to parse input equation(s)");
    vars=sort_identifiers(*_lname(eqv,contextptr)._VECTptr,contextptr);
    gen sol(undef);
    matrice A;
    vecteur Fr;
    bool find_alternate_sols=false,output_sequences=true;
    int n=vars.size(),m=eqv.size();
    /* parse options */
    for (const_iterateur it=opts.begin();it!=opts.end();++it) {
        if (is_equal(*it)) {
            const gen &ls=it->_SYMBptr->feuille._VECTptr->front();
            const gen &rs=it->_SYMBptr->feuille._VECTptr->back();
            if (ls==at_seq && rs.is_integer() && rs.subtype==_INT_BOOLEAN)
                output_sequences=(bool)rs.val;
        } else if (it->type==_IDNT) {
            fvars.push_back(*it);
        } else if (it->type==_VECT) {
            for (const_iterateur jt=it->_VECTptr->begin();jt!=it->_VECTptr->end();++jt) {
                if (jt->type==_IDNT)
                    fvars.push_back(*jt);
                else return generrtype("Expected a variable");
            }
        }
    }
    /* solve equation(s) */
    if (n==1) {
        /* univariate case */
        sol=solve_quiet(eqv,vars.front(),contextptr);
        for (int i=sol._VECTptr->size();i-->0;) {
            if (!sol._VECTptr->at(i).is_integer())
                sol._VECTptr->erase(sol._VECTptr->begin()+i);
        }
    } else if (is_linear_system(eqv,vars,A,Fr,contextptr) && (n>2 || m>1)) {
        /* linear equation with at least three indeterminates or a system of linear equations */
        gen rk=_rank(A,contextptr);
        if (!rk.is_integer() || rk.val<m)
            return generr("System matrix is not full rank");
        vecteur tmp(n);
        for (int i=0;i<n;++i) {
            tmp[i]=temp_symb("tmp_",i,contextptr);
        }
        sol=vecteur(0);
        if (ABS_diophantine(A,Fr,tmp,*sol._VECTptr,contextptr)) {
            vecteur tmp_vars=*_lname(sol,contextptr)._VECTptr;
            for (const_iterateur it=tmp_vars.begin();it!=tmp_vars.end();++it) {
                int i=it-tmp_vars.begin();
                sol=subst(sol,*it,i<int(fvars.size())?fvars[i]:make_integer_placeholder("Z",contextptr),
                          false,contextptr);
            }
        }
    } else if (n==2 && m==1) {
        /* single equation with two indeterminates */
        gen &x=vars.front(),&y=vars.back(),&eq=eqv.front(),a,b,c,d,e,f;
        vecteur t(6);
        int pd1=0,pd2=0;
        bool y_first=false;
        if ((pd1=is_pell_equation(eq,x,y,d,contextptr))!=0 ||
                (pd2=is_pell_equation(eq,y,x,d,contextptr))!=0) {
            /* Pell equation */
            if (output_sequences) {
                gen ph=fvars.empty()?undef:fvars.front();
                sol=pell_general(d,pd1+pd2>0,ph,false,contextptr);
            } else sol=pell_fundamental(d,pd1+pd2>0,contextptr);
            if (pd1==0)
                sol=_apply(makesequence(at_revlist,sol),contextptr);
        } else if (is_binary_quadratic_wrt(eq,x,y,a,b,c,d,e,f,contextptr) &&
                   integralize(makevecteur(a,b,c,d,e,f),t,contextptr)) {
            gen ph(fvars.empty()?undef:fvars.front());
            sol=solve_binary_quadratic(t,ph,false,find_alternate_sols,!output_sequences,contextptr);
        } else if (is_thue_equation(eq,x,y,t[0],t[1],contextptr) ||
                   (y_first=is_thue_equation(eq,y,x,t[0],t[1],contextptr))) {
            /* Thue equation */
            vecteur Pcf=*_coeff(makesequence(t[0],y_first?y:x),contextptr)._VECTptr,cf_int;
            Pcf.push_back(t[1]);
            assert(integralize(Pcf,cf_int,contextptr));
            if (is_zero(cf_int.back())) {
                cf_int.pop_back();
                gen P=_r2e(makesequence(cf_int,x),contextptr);
                gen rsol=_rationalroot(P,contextptr);
                gen ph=fvars.empty()?undef:fvars.front();
                sol=vecteur(0);
                if (rsol.type==_VECT) for (const_iterateur it=rsol._VECTptr->begin();it!=rsol._VECTptr->end();++it) {
                    if (is_zero(*it))
                        continue;
                    if (it->is_integer()) {
                        a=*it; b=1;
                    } else {
                        a=_numer(*it,contextptr);
                        b=_denom(*it,contextptr);
                    }
                    if (is_undef(ph))
                        ph=make_integer_placeholder("Z",contextptr);
                    sol._VECTptr->push_back(makevecteur(a*ph,b*ph));
                }
            } else {
#ifdef __MINGW_H
#ifdef HAVE_LIBPARI
                GEN alpha=gen2GEN(cf_int.back(),vecteur(0),contextptr);
                cf_int.pop_back();
                gen R=expand(_r2e(makesequence(cf_int,x),contextptr),contextptr);
                GEN P=gen2GEN(R,vecteur(1,x),contextptr);
                GEN tnf=thueinit(P,1,DEFAULTPREC);
                sol=GEN2gen(thue(tnf,alpha,NULL),vecteur(0));
#else
                print_warning("PARI library is required for solving Thue equations",contextptr);
#endif
#else
                print_error("failed to solve Thue equation",contextptr);
#endif
            }
            if (sol.type==_VECT && y_first)
                sol=_apply(makesequence(at_revlist,sol),contextptr);
            find_alternate_sols=true;
        } else if (is_tengely(eq,x,y,t[0],t[1],contextptr) || (y_first=is_tengely(eq,y,x,t[0],t[1],contextptr))) {
            /* polynomial equation of type f(x)=g(y) */
            sol=solve_tengely(t[0],t[1],y_first?y:x,y_first?x:y,contextptr);
            if (y_first)
                sol=_apply(makesequence(at_revlist,sol),contextptr);
        }
        if (find_alternate_sols && sol.type==_VECT) for (int i=sol._VECTptr->size();i-->0;) {
            const vecteur &s=*sol._VECTptr->at(i)._VECTptr;
            matrice alt_s;
            alt_s.push_back(makevecteur(s[1],s[0]));
            alt_s.push_back(makevecteur(-s[1],s[0]));
            alt_s.push_back(makevecteur(s[1],-s[0]));
            alt_s.push_back(makevecteur(-s[1],-s[0]));
            alt_s.push_back(makevecteur(-s[0],s[1]));
            alt_s.push_back(makevecteur(s[0],-s[1]));
            alt_s.push_back(makevecteur(-s[0],-s[1]));
            for (const_iterateur it=alt_s.begin();it!=alt_s.end();++it) {
                if (is_zero(_ratnormal(subst(eq,vars,*it,false,contextptr),contextptr)) && !contains(sol,*it))
                    sol._VECTptr->push_back(*it);
            }
        }
    } else if (n==3 && m==1) {
        /* three variables, single equation */
        const gen &eq=eqv.front();
        if (is_homogeneous_poly(eq,vars,NULL,contextptr)==2) {
            /* ternary quadratic form */
            while (fvars.size()<3) {
                fvars.push_back(make_integer_placeholder("Z",contextptr));
            }
            sol=solve_ternary_quadratic_form(eq,vars,fvars,output_sequences,contextptr);
        }
    }
    if (is_undef(sol)) {
        /* check whether a system of polynomial equations with integer coefficients is given, try solving */
        bool is_int_sys=true;
        for (const_iterateur it=eqv.begin();it!=eqv.end();++it) {
            if (!is_int_poly(*it,vars,contextptr)) {
                is_int_sys=false;
                break;
            }
        }
        if (is_int_sys) {
            sol=solve_quiet(eqv,vars,contextptr);
            if (sol.type==_VECT) {
                for (int i=sol._VECTptr->size();i-->0;) {
                    const vecteur &isol=*(sol._VECTptr->at(i))._VECTptr;
                    for (const_iterateur it=isol.begin();it!=isol.end();++it) {
                        if (it->is_integer())
                            continue;
                        if (!is_int_poly(*it,vars,contextptr)) {
                            sol._VECTptr->erase(sol._VECTptr->begin()+i);
                            break;
                        }
                    }
                }
            }
            vecteur sph=*_lname(sol,contextptr)._VECTptr;
            if (!sph.empty()) {
                for (const_iterateur it=sph.begin();it!=sph.end();++it) {
                    int i=it-sph.begin();
                    gen ph=i<int(fvars.size())?fvars[i]:make_integer_placeholder("Z",contextptr);
                    sol=subst(sol,*it,ph,false,contextptr);
                }
            }
        } else return generrtype("Failed to solve equation over integers");
    }
    if (n>1 && sol.type==_VECT) {
        /* format the output */
        if (ckmatrix(*sol._VECTptr,false) && int(sol._VECTptr->front()._VECTptr->size())==n) {
            for (iterateur it=sol._VECTptr->begin();it!=sol._VECTptr->end();++it) {
                *it=_zip(makesequence(at_equal,vars,*it),contextptr);
            }
        } else if (!sol._VECTptr->empty())
            sol=_zip(makesequence(at_equal,vars,sol),contextptr);
        sol=change_subtype(sol,_LIST__VECT);
    }
    return sol;
}
static const char _isolve_s []="isolve";
static define_unary_function_eval (__isolve,&_isolve,_isolve_s);
define_unary_function_ptr5(at_isolve,alias_at_isolve,&__isolve,0,true)

vecteur yvalues(const gen &f,const gen &x,const vecteur &xi,GIAC_CONTEXT) {
    vecteur ret(xi.size());
    for (const_iterateur it=xi.begin();it!=xi.end();++it) {
        ret[it-xi.begin()]=simp(subst(f,x,*it,false,contextptr),contextptr);
    }
    return ret;
}

bool inner_product(const gen &f,const gen &g,const gen &x,const gen &a,const gen &b,const gen &weight,gen &res,GIAC_CONTEXT) {
    bool ret=true;
    try {
        log_output_redirect lor(contextptr);
        res=_gaussquad(makesequence(f*g*weight,x,a,b),contextptr);
    } catch (const std::runtime_error &e) {
        print_error(e.what(),contextptr);
        ret=false;
    }
    if (res.type==_VECT && res._VECTptr->size()==2) {
        res=_mean(res,contextptr);
    }
    if (!is_real_number(res,contextptr)) {
        print_error("inner-product computation failed",contextptr);
        ret=false;
    } else res=to_real_number(res,contextptr);
    return ret;
}

gen fitpoly(const gen &f,const vecteur &vardata,const gen &tol_orig,int maxdeg,int d,GIAC_CONTEXT) {
    vecteur xdata,ydata;
    gen t,tol(tol_orig);
    assert(vardata.size()==3);
    const gen &x=vardata[0],&a=vardata[1],&b=vardata[2];
    if (ckmatrix(f)) {
        int nc=f._VECTptr->front()._VECTptr->size();
        assert(f._VECTptr->size()==2 && nc>=2);
        xdata.reserve(nc);
        ydata.reserve(nc);
        for (int i=0;i<nc;++i) {
            const gen &xi=f._VECTptr->front()._VECTptr->at(i),&yi=f._VECTptr->back()._VECTptr->at(i);
            if (!is_real_number(xi,contextptr) || !is_real_number(yi,contextptr)) {
                print_error("data must be numeric",contextptr);
                return undef;
            }
            if ((!is_undef(a) && is_strictly_greater(a,xi,contextptr)) || (!is_undef(b) && is_strictly_greater(xi,b,contextptr)))
                continue;
            xdata.push_back(xi);
            ydata.push_back(yi);
        }
    } else {
        assert(!is_undef(a) && !is_undef(b));
        t=temp_symb("t",-1,contextptr);
    }
    int m=ydata.size(),n,N=giacmax(63,maxdeg+1);
    vecteur q,Fq,qq,rho,qy,qy0,c,lp,cheby_zeros,cheby_zeros_fvals;
    gen cf,fp(0),f2(0),f2_actual(1),mapx=2*(x-a)/(b-a)-1,corr=gen(2)/(M_PI_2*(b-a));
    if (m==0) {
        if (d<0 && (!inner_product(f,f,x,a,b,_inv(sqrt(1-pow(mapx,2),contextptr),contextptr),f2,contextptr) ||
                    !inner_product(f,f,x,a,b,1,f2_actual,contextptr)))
            return undef;
        tol=tol*corr*f2/f2_actual;
        f2=f2*corr;
        c.reserve(maxdeg+1);
        cheby_zeros.reserve(N+1);
        cheby_zeros_fvals.reserve(N+1);
        for (int j=0;j<=N;++j) {
            cheby_zeros.push_back((b-a)*(1+cos((2.0*j+1)*M_PI_2/(N+1),contextptr))/2+a);
            cheby_zeros_fvals.push_back(subst(f,x,cheby_zeros.back(),false,contextptr));
        }
    }
    // find the optimal degree n
    const_iterateur czf_iter,cz_iter,cz_itend=cheby_zeros.end();
    for (n=0;(m==0 || n<m) && n<=maxdeg+(d<0?0:d);++n) {
        if (m==0) {
            gen T=_e2r(makesequence(subst(_tchebyshev1(makesequence(n,t),contextptr),t,mapx,false,contextptr),x),contextptr);
            cf=0;
            for (cz_iter=cheby_zeros.begin(),czf_iter=cheby_zeros_fvals.begin();cz_iter!=cz_itend;++cz_iter,++czf_iter) {
                cf+=*czf_iter*_peval(makesequence(T,*cz_iter),contextptr);
            }
            cf=(2.0*cf)/(N+1);
            fp+=pow(cf,2)/(n==0?2:1);
            c.push_back(cf);
            lp.push_back(T);
            if ((d>=0 && n==d) || (d<0 && is_greater(tol,_abs(f2-fp,contextptr),contextptr)))
                break; // the desired precision or degree is reached
        } else {
            if (n==0)
                q.push_back(1);
            else if (n==1)
                q.push_back(x-_sum(xdata,contextptr)/qq[0]);
            else {
                vecteur yval=yvalues(x*q[n-1],x,xdata,contextptr);
                gen alpha=scalarproduct(yval,qy,contextptr)/qq[n-1];
                gen beta=scalarproduct(yval,qy0,contextptr)/qq[n-2];
                q.push_back(expand((x-alpha)*q[n-1]-beta*q[n-2],contextptr));
            }
            qy0=n==1?vecteur(xdata.size(),1):qy;
            qy=yvalues(q[n],x,xdata,contextptr);
            Fq.push_back(scalarproduct(ydata,qy,contextptr));
            qq.push_back(scalarproduct(qy,qy,contextptr));
            rho.push_back(max(0,n==0?scalarproduct(ydata,ydata,contextptr)-pow(Fq[0],2)/qq[0]:rho[n-1]-pow(Fq[n],2)/qq[n],contextptr));
            if (is_zero(rho[n]/m,contextptr) || (d<=0 && n==-d))
                break; // near-perfect fit or the desired degree is reached
            if (d>0 && n>=d-1 && is_greater(tol,_stddev(vecteur(rho.begin()+n-d+1,rho.end()),contextptr)/_mean(rho,contextptr),contextptr)) {
                n-=d-1;
                break;
            }
        }
    }
    if (m>0 && n>=m)
        n=m-1;
    if (n>maxdeg) {
        n=maxdeg;
        print_warning("degree limit reached",contextptr);
    }
    // build the polynomial
    gen p(0);
    for (int i=0;i<=n;++i) p+=(m==0?c[i]*lp[i]/(i==0?2:1):q[i]*Fq[i]/qq[i]);
    if (m==0) p=_r2e(makesequence(p,x),contextptr);
    return _evalf(expand(symb_normal(p),contextptr),contextptr);
}

gen _fitpoly(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur x,y;
    gen var=temp_symb("x",-1,contextptr),f(undef),a(undef),b(undef),p,tol(0);
    bool has_var=false;
    int maxdeg=15,d=5,deg=-1;
    if (g.type!=_VECT)
        return generrtype("Expected a matrix or a sequence of arguments");
    if (g.subtype==_SEQ__VECT) {
        const vecteur &args=*g._VECTptr;
        size_t argn=2;
        if (args.size()<2)
            return generrdim("Too few input arguments");
        if (args[0].type!=_VECT && args[0].type!=_SYMB)
            return generrtype("Expected a list of x-values, expression, or matrix");
        if (args[0].type==_SYMB) {
            f=args[0];
            vecteur parm;
            if ((parm=*_lname(f,contextptr)._VECTptr).size()!=1)
                return generr("Expected an univariate expression");
            var=parm.front();
            --argn;
        } else if (ckmatrix(args[0])) {
            if (args[0]._VECTptr->front()._VECTptr->size()!=2)
                return generrdim("Expected a two-column matrix");
            matrice data=mtran(*args[0]._VECTptr);
            x=*data.front()._VECTptr;
            y=*data.back()._VECTptr;
            --argn;
        } else {
            if (args[1].type!=_VECT)
                return generrtype("Expected a list of y-values");
            x=*args[0]._VECTptr;
            y=*args[1]._VECTptr;
        }
        if (argn<args.size()) {
            const gen &v=args[argn];
            bool iseq;
            if (v.type==_IDNT) {
                var=v;
                has_var=true;
                argn++;
            } else if ((iseq=is_equal(v) && (v._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_interval) ||
                                             v._SYMBptr->feuille._VECTptr->front().type==_IDNT)) ||
                                             v.is_symb_of_sommet(at_interval)) {
                if (iseq) {
                    if ((var=v._SYMBptr->feuille._VECTptr->front()).type!=_IDNT)
                        return generr("Expected a variable");
                    if (!v._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_interval))
                        return generr("Expected a range");
                    has_var=true;
                }
                const gen &intrv=iseq?v._SYMBptr->feuille._VECTptr->back():v;
                a=intrv._SYMBptr->feuille._VECTptr->front();
                b=intrv._SYMBptr->feuille._VECTptr->back();
                if (!is_real_number(a,contextptr) || !is_real_number(b,contextptr) || is_greater(a,b,contextptr))
                    return generr("Invalid range");
                a=to_real_number(a,contextptr);
                b=to_real_number(b,contextptr);
                argn++;
            }
        }
        // parse options
        for (const_iterateur it=args.begin()+argn;it!=args.end();++it) {
            if (is_equal(*it)) {
                const gen &lh=it->_SYMBptr->feuille._VECTptr->front();
                const gen &rh=it->_SYMBptr->feuille._VECTptr->back();
                if (lh==at_limit) {
                    if (!rh.is_integer() || rh.val<1)
                        return generr("Expected a positive integer");
                    maxdeg=rh.val;
                } else if (lh==at_threshold) {
                    if (!is_real_number(rh,contextptr) || is_positive(-rh,contextptr))
                        return generr("Expected a positive real number");
                    tol=to_real_number(rh,contextptr);
                } else if (lh==at_len || lh==at_length) {
                    if (!rh.is_integer() || rh.val<1)
                        return generr("Expected a positive integer");
                    d=rh.val;
                } else if (lh==at_degree) {
                    if (!rh.is_integer() || rh.val<0)
                        return generr("Expected a nonnegative integer");
                    deg=rh.val;
                } else return generr("Option not supported");
            } else if (is_real_number(*it,contextptr) && is_strictly_positive(*it,contextptr)) {
                tol=to_real_number(*it,contextptr);
            } else return generr("Option not supported");
        }
    } else {
        if (!ckmatrix(g) || g._VECTptr->front()._VECTptr->size()!=2)
            return generrdim("Expected a two-column matrix");
        matrice data=mtran(*g._VECTptr);
        x=*data.front()._VECTptr;
        y=*data.back()._VECTptr;
    }
    vecteur vardata=makevecteur(var,a,b);
    if (is_undef(f)) {
        if (x.size()!=y.size())
            return generr("Sizes of x-values and y-values do not match");
        if (tol==0)
            tol=0.01;
        int m=x.size();
        if (m<2)
            return generrdim("Too few input data points");
        if (d>m)
            d=m;
        p=fitpoly(exact(makevecteur(x,y),contextptr),vardata,tol,maxdeg,deg>=0?-deg:d,contextptr);
        if (!has_var)
            p=_e2r(makesequence(p,var),contextptr);
    } else {
        if (is_undef(a) || is_undef(b))
            return generr("Variable range is required");
        if (tol==0)
            tol=1e-5;
        p=fitpoly(f,vardata,tol,maxdeg,deg,contextptr);
    }
    return p;
}
static const char _fitpoly_s []="fitpoly";
static define_unary_function_eval (__fitpoly,&_fitpoly,_fitpoly_s);
define_unary_function_ptr5(at_fitpoly,alias_at_fitpoly,&__fitpoly,0,true)

// Routine for solving K-means problem
// Source: https://people.math.sc.edu/Burkardt/cpp_src/asa136/asa136.html
//
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

void optra ( double a[], int m, int n, double c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int live[], int *indx )

//****************************************************************************80
//
//  Purpose:
//
//    OPTRA carries out the optimal transfer stage.
//
//  Discussion:
//
//    This is the optimal transfer stage.
//
//    Each point is re-allocated, if necessary, to the cluster that
//    will induce a maximum reduction in the within-cluster sum of
//    squares.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by John Hartigan, Manchek Wong.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hartigan, Manchek Wong,
//    Algorithm AS 136:
//    A K-Means Clustering Algorithm,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, pages 100-108.
//
//  Parameters:
//
//    Input, double A(M,N), the points.
//
//    Input, int M, the number of points.
//
//    Input, int N, the number of spatial dimensions.
//
//    Input/output, double C(K,N), the cluster centers.
//
//    Input, int K, the number of clusters.
//
//    Input/output, int IC1(M), the cluster to which each 
//    point is assigned.
//
//    Input/output, int IC2(M), used to store the cluster 
//    which each point is most likely to be transferred to at each step.
//
//    Input/output, int NC(K), the number of points in 
//    each cluster.
//
//    Input/output, double AN1(K).
//
//    Input/output, double AN2(K).
//
//    Input/output, int NCP(K).
//
//    Input/output, double D(M).
//
//    Input/output, int ITRAN(K).
//
//    Input/output, int LIVE(K).
//
//    Input/output, int *INDX, the number of steps since a 
//    transfer took place.
//
{
  double al1;
  double al2;
  double alt;
  double alw;
  double da;
  double db;
  double dc;
  double dd;
  double de;
  double df;
  int i;
  int j;
  int l;
  int l1;
  int l2;
  int ll;
  double r2;
  double rr;
//
//  If cluster L is updated in the last quick-transfer stage, it
//  belongs to the live set throughout this stage.   Otherwise, at
//  each step, it is not in the live set if it has not been updated
//  in the last M optimal transfer steps.
//
  for ( l = 1; l <= k; l++ )
  {
    if ( itran[l-1] == 1)
    {
      live[l-1] = m + 1;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    *indx = *indx + 1;
    l1 = ic1[i-1];
    l2 = ic2[i-1];
    ll = l2;
//
//  If point I is the only member of cluster L1, no transfer.
//
    if ( 1 < nc[l1-1]  )
    {
//
//  If L1 has not yet been updated in this stage, no need to
//  re-compute D(I).
//
      if ( ncp[l1-1] != 0 )
      {
        de = 0.0;
        for ( j = 1; j <= n; j++ )
        {
          df = a[i-1+(j-1)*m] - c[l1-1+(j-1)*k];
          de = de + df * df;
        }
        d[i-1] = de * an1[l1-1];
      }
//
//  Find the cluster with minimum R2.
//
      da = 0.0;
      for ( j = 1; j <= n; j++ )
      {
        db = a[i-1+(j-1)*m] - c[l2-1+(j-1)*k];
        da = da + db * db;
      }
      r2 = da * an2[l2-1];

      for ( l = 1; l <= k; l++ )
      {
//
//  If LIVE(L1) <= I, then L1 is not in the live set.   If this is
//  true, we only need to consider clusters that are in the live set
//  for possible transfer of point I.   Otherwise, we need to consider
//  all possible clusters.
//
        if ( ( i < live[l1-1] || i < live[l2-1] ) && l != l1 && l != ll )
        {
          rr = r2 / an2[l-1];

          dc = 0.0;
          for ( j = 1; j <= n; j++ )
          {
            dd = a[i-1+(j-1)*m] - c[l-1+(j-1)*k];
            dc = dc + dd * dd;
          }

          if ( dc < rr )
          {
            r2 = dc * an2[l-1];
            l2 = l;
          }
        }
      }
//
//  If no transfer is necessary, L2 is the new IC2(I).
//
      if ( d[i-1] <= r2 )
      {
        ic2[i-1] = l2;
      }
//
//  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
//  L2, and update IC1(I) and IC2(I).
//
      else
      {
        *indx = 0;
        live[l1-1] = m + i;
        live[l2-1] = m + i;
        ncp[l1-1] = i;
        ncp[l2-1] = i;
        al1 = ( double ) ( nc[l1-1] );
        alw = al1 - 1.0;
        al2 = ( double ) ( nc[l2-1] );
        alt = al2 + 1.0;
        for ( j = 1; j <= n; j++ )
        {
          c[l1-1+(j-1)*k] = ( c[l1-1+(j-1)*k] * al1 - a[i-1+(j-1)*m] ) / alw;
          c[l2-1+(j-1)*k] = ( c[l2-1+(j-1)*k] * al2 + a[i-1+(j-1)*m] ) / alt;
        }
        nc[l1-1] = nc[l1-1] - 1;
        nc[l2-1] = nc[l2-1] + 1;
        an2[l1-1] = alw / al1;
        if ( 1.0 < alw )
        {
          an1[l1-1] = alw / ( alw - 1.0 );
        }
        else
        {
          an1[l1-1] = r8_huge ( );
        }
        an1[l2-1] = alt / al2;
        an2[l2-1] = alt / ( alt + 1.0 );
        ic1[i-1] = l2;
        ic2[i-1] = l1;
      }
    }

    if ( *indx == m )
    {
      return;
    }
  }
//
//  ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
//  decreased by M before re-entering OPTRA.
//
  for ( l = 1; l <= k; l++ )
  {
    itran[l-1] = 0;
    live[l-1] = live[l-1] - m;
  }

  return;
}
//****************************************************************************80

void qtran ( double a[], int m, int n, double c[], int k, int ic1[], 
  int ic2[], int nc[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int *indx )

//****************************************************************************80
//
//  Purpose:
//
//    QTRAN carries out the quick transfer stage.
//
//  Discussion:
//
//    This is the quick transfer stage.
//
//    IC1(I) is the cluster which point I belongs to.
//    IC2(I) is the cluster which point I is most likely to be
//    transferred to.
//
//    For each point I, IC1(I) and IC2(I) are switched, if necessary, to
//    reduce within-cluster sum of squares.  The cluster centers are
//    updated after each step.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by John Hartigan, Manchek Wong.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Hartigan, Manchek Wong,
//    Algorithm AS 136:
//    A K-Means Clustering Algorithm,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, pages 100-108.
//
//  Parameters:
//
//    Input, double A(M,N), the points.
//
//    Input, int M, the number of points.
//
//    Input, int N, the number of spatial dimensions.
//
//    Input/output, double C(K,N), the cluster centers.
//
//    Input, int K, the number of clusters.
//
//    Input/output, int IC1(M), the cluster to which each 
//    point is assigned.
//
//    Input/output, int IC2(M), used to store the cluster 
//    which each point is most likely to be transferred to at each step.
//
//    Input/output, int NC(K), the number of points in 
//    each cluster.
//
//    Input/output, double AN1(K).
//
//    Input/output, double AN2(K).
//
//    Input/output, int NCP(K).
//
//    Input/output, double D(M).
//
//    Input/output, int ITRAN(K).
//
//    Input/output, int INDX, counts the number of steps 
//    since the last transfer.
//
{
  double al1;
  double al2;
  double alt;
  double alw;
  double da;
  double db;
  double dd;
  double de;
  int i;
  int icoun;
  int istep;
  int j;
  int l1;
  int l2;
  double r2;
//
//  In the optimal transfer stage, NCP(L) indicates the step at which
//  cluster L is last updated.   In the quick transfer stage, NCP(L)
//  is equal to the step at which cluster L is last updated plus M.
//
  icoun = 0;
  istep = 0;

  for ( ; ; )
  {
    for ( i = 1; i <= m; i++ )
    {
      icoun = icoun + 1;
      istep = istep + 1;
      l1 = ic1[i-1];
      l2 = ic2[i-1];
//
//  If point I is the only member of cluster L1, no transfer.
//
      if ( 1 < nc[l1-1] )
      {
//
//  If NCP(L1) < ISTEP, no need to re-compute distance from point I to
//  cluster L1.   Note that if cluster L1 is last updated exactly M
//  steps ago, we still need to compute the distance from point I to
//  cluster L1.
//
        if ( istep <= ncp[l1-1] )
        {
          da = 0.0;
          for ( j = 1; j <= n; j++ )
          {
            db = a[i-1+(j-1)*m] - c[l1-1+(j-1)*k];
            da = da + db * db;
          }
          d[i-1] = da * an1[l1-1];
        }
//
//  If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will be no transfer of
//  point I at this step.
//
        if ( istep < ncp[l1-1] || istep < ncp[l2-1] )
        {
          r2 = d[i-1] / an2[l2-1];

          dd = 0.0;
          for ( j = 1; j <= n; j++ )
          {
            de = a[i-1+(j-1)*m] - c[l2-1+(j-1)*k];
            dd = dd + de * de;
          }
//
//  Update cluster centers, NCP, NC, ITRAN, AN1 and AN2 for clusters
//  L1 and L2.   Also update IC1(I) and IC2(I).   Note that if any
//  updating occurs in this stage, INDX is set back to 0.
//
          if ( dd < r2 )
          {
            icoun = 0;
            *indx = 0;
            itran[l1-1] = 1;
            itran[l2-1] = 1;
            ncp[l1-1] = istep + m;
            ncp[l2-1] = istep + m;
            al1 = ( double ) ( nc[l1-1] );
            alw = al1 - 1.0;
            al2 = ( double ) ( nc[l2-1] );
            alt = al2 + 1.0;
            for ( j = 1; j <= n; j++ )
            {
              c[l1-1+(j-1)*k] = ( c[l1-1+(j-1)*k] * al1 - a[i-1+(j-1)*m] ) / alw;
              c[l2-1+(j-1)*k] = ( c[l2-1+(j-1)*k] * al2 + a[i-1+(j-1)*m] ) / alt;
            }
            nc[l1-1] = nc[l1-1] - 1;
            nc[l2-1] = nc[l2-1] + 1;
            an2[l1-1] = alw / al1;
            if ( 1.0 < alw )
            {
              an1[l1-1] = alw / ( alw - 1.0 );
            }
            else
            {
              an1[l1-1] = r8_huge ( );
            }
            an1[l2-1] = alt / al2;
            an2[l2-1] = alt / ( alt + 1.0 );
            ic1[i-1] = l2;
            ic2[i-1] = l1;
          }
        }
      }
//
//  If no re-allocation took place in the last M steps, return.
//
      if ( icoun == m )
      {
        return;
      }
    }
  }
}
//****************************************************************************80

void kmns ( double a[], int m, int n, double c[], int k, int ic1[], int nc[], 
  int iter, double wss[], int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    KMNS carries out the K-means algorithm.
//
//  Discussion:
//
//    This routine attempts to divide M points in N-dimensional space into 
//    K clusters so that the within cluster sum of squares is minimized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by John Hartigan, Manchek Wong.
//    C++ version by John Burkardt
//
//  Reference:
//
//    John Hartigan, Manchek Wong,
//    Algorithm AS 136:
//    A K-Means Clustering Algorithm,
//    Applied Statistics,
//    Volume 28, Number 1, 1979, pages 100-108.
//
//  Parameters:
//
//    Input, double A(M,N), the points.
//
//    Input, int M, the number of points.
//
//    Input, int N, the number of spatial dimensions.
//
//    Input/output, double C(K,N), the cluster centers.
//
//    Input, int K, the number of clusters.
//
//    Output, int IC1(M), the cluster to which each point 
//    is assigned.
//
//    Output, int NC(K), the number of points in each cluster.
//
//    Input, int ITER, the maximum number of iterations allowed.
//
//    Output, double WSS(K), the within-cluster sum of squares
//    of each cluster.
//
//    Output, int *IFAULT, error indicator.
//    0, no error was detected.
//    1, at least one cluster is empty after the initial assignment.  A better
//       set of initial cluster centers is needed.
//    2, the allowed maximum number of iterations was exceeded.
//    3, K is less than or equal to 1, or greater than or equal to M.
//
{
  double aa;
  double *an1;
  double *an2;
  double *d;
  double da;
  double db;
  double dc;
  double dt[2];
  int i;
  int *ic2;
  int ii;
  int ij;
  int il;
  int indx;
  int *itran;
  int j;
  int l;
  int *live;
  int *ncp;
  double temp;

  *ifault = 0;

  if ( k <= 1 || m <= k )
  {
    *ifault = 3;
    return;
  }
  ic2 = new int[m];
  an1 = new double[k];
  an2 = new double[k];
  ncp = new int[k];
  d = new double[m];
  itran = new int[k];
  live = new int[k];
//
//  For each point I, find its two closest centers, IC1(I) and
//  IC2(I).  Assign the point to IC1(I).
//
  for ( i = 1; i <= m; i++ )
  {
    ic1[i-1] = 1;
    ic2[i-1] = 2;

    for ( il = 1; il <= 2; il++ )
    {
      dt[il-1] = 0.0;
      for ( j = 1; j <= n; j++ )
      {
        da = a[i-1+(j-1)*m] - c[il-1+(j-1)*k];
        dt[il-1] = dt[il-1] + da * da;
      }
    }

    if ( dt[1] < dt[0] )
    {
      ic1[i-1] = 2;
      ic2[i-1] = 1;
      temp = dt[0];
      dt[0] = dt[1];
      dt[1] = temp;
    }

    for ( l = 3; l <= k; l++ )
    {
      db = 0.0;
      for ( j = 1; j <= n; j++ )
      {
        dc = a[i-1+(j-1)*m] - c[l-1+(j-1)*k];
        db = db + dc * dc;
      }

      if ( db < dt[1] )
      {
        if ( dt[0] <= db )
        {
          dt[1] = db;
          ic2[i-1] = l;
        }
        else
        {
          dt[1] = dt[0];
          ic2[i-1] = ic1[i-1];
          dt[0] = db;
          ic1[i-1] = l;
        }
      }
    }
  }
//
//  Update cluster centers to be the average of points contained within them.
//
  for ( l = 1; l <= k; l++ )
  {
    nc[l-1] = 0;
    for ( j = 1; j <= n; j++ )
    {
      c[l-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    l = ic1[i-1];
    nc[l-1] = nc[l-1] + 1;
    for ( j = 1; j <= n; j++ )
    {
      c[l-1+(j-1)*k] = c[l-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }
//
//  Check to see if there is any empty cluster at this stage.
//
  *ifault = 1;

  for ( l = 1; l <= k; l++ )
  {
    if ( nc[l-1] == 0 )
    {
      *ifault = 1;
      return;
    }

  }

  *ifault = 0;

  for ( l = 1; l <= k; l++ )
  {
    aa = ( double ) ( nc[l-1] );

    for ( j = 1; j <= n; j++ )
    {
      c[l-1+(j-1)*k] = c[l-1+(j-1)*k] / aa;
    }
//
//  Initialize AN1, AN2, ITRAN and NCP.
//
//  AN1(L) = NC(L) / (NC(L) - 1)
//  AN2(L) = NC(L) / (NC(L) + 1)
//  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
//           = 0 otherwise
//
//  In the optimal-transfer stage, NCP(L) stores the step at which
//  cluster L is last updated.
//
//  In the quick-transfer stage, NCP(L) stores the step at which
//  cluster L is last updated plus M.
//
    an2[l-1] = aa / ( aa + 1.0 );

    if ( 1.0 < aa )
    {
      an1[l-1] = aa / ( aa - 1.0 );
    }
    else
    {
      an1[l-1] = r8_huge ( );
    }
    itran[l-1] = 1;
    ncp[l-1] = -1;
  }

  indx = 0;
  *ifault = 2;

  for ( ij = 1; ij <= iter; ij++ )
  {
//
//  In this stage, there is only one pass through the data.   Each
//  point is re-allocated, if necessary, to the cluster that will
//  induce the maximum reduction in within-cluster sum of squares.
//
    optra ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &indx );
//
//  Stop if no transfer took place in the last M optimal transfer steps.
//
    if ( indx == m )
    {
      *ifault = 0;
      break;
    }
//
//  Each point is tested in turn to see if it should be re-allocated
//  to the cluster to which it is most likely to be transferred,
//  IC2(I), from its present cluster, IC1(I).   Loop through the
//  data until no further change is to take place.
//
    qtran ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &indx );
//
//  If there are only two clusters, there is no need to re-enter the
//  optimal transfer stage.
//
    if ( k == 2 )
    {
      *ifault = 0;
      break;
    }
//
//  NCP has to be set to 0 before entering OPTRA.
//
    for ( l = 1; l <= k; l++ )
    {
      ncp[l-1] = 0;
    }

  }
//
//  Compute the within-cluster sum of squares for each cluster.
//
  for ( l = 1; l <= k; l++ )
  {
    wss[l-1] = 0.0;
    for ( j = 1; j <= n; j++ )
    {
      c[l-1+(j-1)*k] = 0.0;
    }
  }

  for ( i = 1; i <= m; i++ )
  {
    ii = ic1[i-1];
    for ( j = 1; j <= n; j++ )
    {
      c[ii-1+(j-1)*k] = c[ii-1+(j-1)*k] + a[i-1+(j-1)*m];
    }
  }

  for ( j = 1; j <= n; j++ )
  {
    for ( l = 1; l <= k; l++ )
    {
      c[l-1+(j-1)*k] = c[l-1+(j-1)*k] / ( double ) ( nc[l-1] );
    }
    for ( i = 1; i <= m; i++ )
    {
      ii = ic1[i-1];
      da = a[i-1+(j-1)*m] - c[ii-1+(j-1)*k];
      wss[ii-1] = wss[ii-1] + da * da;
    }
  }

  delete [] ic2;
  delete [] an1;
  delete [] an2;
  delete [] ncp;
  delete [] d;
  delete [] itran;
  delete [] live;

  return;
}

void transpose_array(double *ar,int k,int n) {
    double tmp;
    for (int i=0;i<k;++i) {
        for (int j=0;j<n;++j) {
            tmp=ar[i+j*k];
            ar[i+j*k]=ar[j+i*n];
            ar[j+i*n]=tmp;
        }
    }
}

/* Initialize cluster centers k data points randomly */
void init_cluster_centers(const double *x,const int *ind,int m,int n,int k,double *d,GIAC_CONTEXT) {
    map<int,vecteur> partition;
    for (int i=0;i<m;++i) partition[ind==NULL?1:ind[i]].push_back(i);
    vecteur ch;
    ch.reserve(k);
    int lc,maxsize=0;
    for (map<int,vecteur>::iterator it=partition.begin();it!=partition.end();++it) {
        if (int(it->second.size())>maxsize) {
            maxsize=it->second.size();
            lc=it->first;
        }
        int pos=_rand(it->second.size(),contextptr).val;
        ch.push_back(it->second[pos]);
        it->second.erase(it->second.begin()+pos);
    }
    if (int(ch.size())<k) {
        ch=mergevecteur(ch,*_rand(makesequence(k-int(ch.size()),partition[lc]),contextptr)._VECTptr);
    }
    assert(int(ch.size())==k);
    for (int i=0;i<k;++i) {
        for (int j=0;j<n;++j) {
            d[i+j*k]=x[ch[i].val+j*m];
        }
    }
}

/* Compute list RES of cluster indices for data points using an efficient k-means
 * implementation AS 136 (Hartigan & Wong). Set K to the desired number of clusters
 * or to a negative value, in which case the cluster index -K is computed.
 * MIN_K and MAX_K are bounds for K. The total number of clusters is written to K.
 * Cluster centers and average sizes [center,size] are stored to CC unless the latter is NULL.
 * The list of values of used clustering indices is stored in IND unless the latter is NULL. */
int kmeans(const matrice &data,int &k,int min_k,int &max_k,int maxiter,vector<int> &res,vecteur *cc,matrice *ind,GIAC_CONTEXT) {
    int m=mrows(data),n=mcols(data); // m data points as vectors of size n
    double *x=new double[m*n];       // data points
    int *b=new int[m];               // cluster indices for data points
    double *d=NULL;                  // cluster centers
    double *dev=NULL;                // sums of squared deviations from centers
    int *e=NULL;                     // number of points in each cluster
    int ifault;
    /* input data points */
    for (int i=0;i<m;++i) {
        for (int j=0;j<n;++j) {
            x[i+j*m]=to_real_number(data[i][j],contextptr).to_double(contextptr);
        }
    }
    res.resize(m);
    int tries=0;
    if (k>0) { // the number of clusters is fixed (k)
        d=new double[k*n];
        dev=new double[k];
        e=new int[k];
        do { // try to generate initial centers at most 10 times
            init_cluster_centers(x,NULL,m,n,k,d,contextptr);
            kmns(x,m,n,d,k,b,e,maxiter,dev,&ifault);
        } while (ifault==1 && tries++<10);
        if (ifault!=0)
            return ifault;
        for (int i=0;i<m;++i) res[i]=b[i]-1;
    } else { // guess the optimal number of clusters
        double W,W_prev=0,score;
        int best_k=0,meth=-k;
        int *b_best=meth==_HARTIGAN_CRITERION?new int[m]:NULL;
        dev=new double[max_k];
        e=new int[max_k];
        if (b_best!=NULL) for (int i=0;i<m;++i) b_best[i]=i+1;
        bool rough_optimum_found=false;
        cluster_crit ccr(data,at_longueur2,contextptr);
        map<int,vector<int> > b_map;
        vecteur kvals;
        kvals.reserve(max_k-min_k+1);
        for (k=min_k;k<=max_k;++k) {
            if (d==NULL)
                d=new double[max_k*n];
            tries=0;
            do { // try to generate initial centers at most 10 times
                init_cluster_centers(x,k==min_k||tries>0?NULL:b,m,n,k,d,contextptr);
                kmns(x,m,n,d,k,b,e,maxiter,dev,&ifault);
            } while (ifault==1 && tries++<10);
            if (ifault!=0) return ifault;
            if (!rough_optimum_found) {
                W=0;
                for (int i=0;i<k;++i) W+=dev[i];
                if (W>0 && W_prev>0) {
                    score=(W_prev/W-1)*(m-k);
                    if (score<=10) {
                        rough_optimum_found=true;
                        if (meth==_HARTIGAN_CRITERION) {
                            k--;
                            break;
                        }
                        max_k=std::min(max_k,2*k);
                    }
                }
                W_prev=W;
                if (meth==_HARTIGAN_CRITERION) {
                    for (int i=0;i<m;++i) b_best[i]=b[i];
                    if (W==0) break;
                }
            }
            if (meth>0) {
                ccr.init(k,b,d,dev);
                ccr.compute_indices(meth);
                kvals.push_back(k);
                vector<int> lst(m);
                for (int i=0;i<m;++i) lst[i]=b[i]-1;
                b_map.insert(make_pair(k,lst));
            }
        }
        if (meth==_HARTIGAN_CRITERION) {
            if (k>max_k) k=max_k;
            for (int i=0;i<m;++i) res[i]=(rough_optimum_found?b_best[i]:b[i])-1;
            delete[] b_best;
        } else {
            k=ccr.ncopt();
            res=b_map[k];
        }
        if (ind!=NULL) { // store clustering indices
            vecteur indval=ccr.get_index_values();
            gen ndef=string2gen(gettext("not defined"),false);
            for (iterateur it=indval.begin();it!=indval.end();++it) {
                if (is_undef(*it))
                    *it=ndef;
            }
            *ind=mtran(makevecteur(kvals,indval));
        }
    }
    if (cc!=NULL && d!=NULL) { // store cluster centers & sizes
        cc->resize(k);
        vecteur xc(n);
        for (int i=0;i<k;++i) {
            for (int j=0;j<n;++j) xc[j]=d[i+j*k];
            cc->at(i)=makevecteur(xc,std::sqrt(dev[i]/e[i]));
        }
    }
    if (d!=NULL) delete[] d;
    delete[] dev;
    delete[] e;
    delete[] x;
    delete[] b;
    return 0;
}

gen comp_dist(const gen &a,const gen &b,const gen &dist_func,GIAC_CONTEXT) {
    gen val;
    try {
        val=_evalf(dist_func(makesequence(a,b),contextptr),contextptr);
    } catch (const std::runtime_error &e) {
        return undef;
    }
    if (!is_real_number(val,contextptr))
        return undef;
    return to_real_number(val,contextptr);
}

/*
 ****** IMPLEMENTATION OF CLUSTER_CRIT CLASS *******
 *
 * Various cluster indices for obtaining optimal clustering
 * SOURCE: Bernard Desgraupes, "Clustering Indices" (2017) -- package clusterCrit for R
 * https://cran.r-project.org/web/packages/clusterCrit/vignettes/clusterCrit.pdf
 */
cluster_crit::cluster_crit(const matrice &data,const gen &dist_func,GIAC_CONTEXT)
    : ctx(contextptr),pts(data),distf(dist_func)
{
    N=mrows(pts);
    p=mcols(pts);
    K=0;
    _mu.resize(1,_mean(pts,ctx)); // global barycenter
    _ci.resize(N,-1);
    _mu.reserve(N+1);
    _n.reserve(N);
    _WGSSk.reserve(N);
    _WG.reserve(N);
    _T=midn(p);
    _BG=midn(p);
    _TSS.resize(p);
    _BGSS.resize(p);
    ET=0;
    for (int i=0;i<N;++i) ET+=comp_dist(pts[i],mu(),distf,ctx);
}
gen cluster_crit::Delta(int k,int l) {
    pair<int,int> kl=make_pair(k<l?k:l,k<l?l:k);
    map<pair<int,int>,gen>::const_iterator it=Delta_map.find(kl);
    if (it==Delta_map.end())
        it=Delta_map.insert(make_pair(kl,comp_dist(mu(k),mu(l),(distf==at_longueur2?at_longueur:distf),ctx))).first;
    return it->second;
}
void cluster_crit::init(int k,const int *indices,const double *bc,const double *wgss) {
    assert(k>0 && K<=N);
    K=k;
    _n.resize(K);
    std::fill(_n.begin(),_n.end(),0);
    for (int i=0;i<N;++i) _n[_ci[i]=indices[i]-1]++;
    _mu.resize(K+1);
    _WG.resize(K);
    _WGSSk.resize(K);
    _WGSS=0;
    _BGSS[0]=_TSS[0]=undef;
    for (int i=0;i<K;++i) {
        _WGSS+=(_WGSSk[i]=wgss[i]);
        if (_mu[i+1].type!=_VECT)
            _mu[i+1]=vecteur(p);
        for (int j=0;j<p;++j) {
            _mu[i+1]._VECTptr->at(j)=bc[i+j*K];
        }
    }
    _DB_min=_DB_max=undef;
    have_T=have_WG=have_BG=false;
    Delta_map.clear();
}
const vecteur &cluster_crit::T() {
    if (!have_T) {
        for (int i=0;i<p;++i) {
            _TSS[i]=0;
            for (int j=i;j<p;++j) {
                _T[i]._VECTptr->at(j)=0;
                for (int l=0;l<N;++l) {
                    _T[i]._VECTptr->at(j)+=(pts[l][i]-mu()[i])*(pts[l][j]-mu()[j]);
                }
                if (i!=j)
                    _T[j]._VECTptr->at(i)=_T[i][j];
                else _TSS[i]=_T[i][i];
            }
        }
        have_T=true;
    }
    return _T;
}
gen cluster_crit::TSS(int k) {
    if (is_undef(_TSS.front())) {
        for (int j=0;j<p;++j) {
            _TSS[j]=0;
            for (int l=0;l<N;++l) {
                _TSS[j]+=pow(pts[l][j]-mu()[j],2);
            }
        }
    }
    return k<0?_sum(_TSS,ctx):_TSS[k];
}
const vecteur &cluster_crit::WG(int l) {
    if (!have_WG) {
        _WG.resize(K);
        _WGsum.clear();
        for (int k=0;k<K;++k) {
            if (_WG[k].type!=_VECT)
                _WG[k]=zero_mat(p,p,ctx);
            matrice &WGk=*_WG[k]._VECTptr;
            for (int i=0;i<p;++i) {
                for (int j=i;j<p;++j) {
                    WGk[i]._VECTptr->at(j)=0;
                    for (int l=0;l<N;++l) {
                        if (ci(l)==k)
                            WGk[i]._VECTptr->at(j)+=(pts[l][i]-mu(k)[i])*(pts[l][j]-mu(k)[j]);
                    }
                    if (i!=j)
                        WGk[j]._VECTptr->at(i)=_T[i][j];
                }
            }
            _WGsum=_WGsum.empty()?WGk:madd(_WGsum,WGk);
        }
        have_WG=true;
    }
    return l<0?_WGsum:*_WG[l]._VECTptr;
}
const vecteur &cluster_crit::BG() {
    if (!have_BG) {
        for (int i=0;i<p;++i) {
            _BGSS[i]=0;
            for (int j=i;j<p;++j) {
                _BG[i]._VECTptr->at(j)=0;
                for (int k=0;k<K;++k) {
                    _BG[i]._VECTptr->at(j)+=n(k)*(mu(k)[i]-mu()[i])*(mu(k)[j]-mu()[j]);
                }
                if (i!=j)
                    _BG[j]._VECTptr->at(i)=_BG[i][j];
                else _BGSS[i]=_BG[i][i];
            }
        }
        have_BG=true;
    }
    return _BG;
}
gen cluster_crit::BGSS(int k) {
    if (is_undef(_BGSS.front())) {
        for (int i=0;i<p;++i) {
            _BGSS[i]=0;
            for (int k=0;k<K;++k) {
                _BGSS[i]+=n(k)*pow(mu(k)[i]-mu()[i],2);
            }
        }
    }
    return k<0?_sum(_BGSS,ctx):_BGSS[k];
}
gen cluster_crit::DB(int dir) {
    if (is_undef(_DB_min) || is_undef(_DB_max)) {
        _DB_min=plus_inf;
        _DB_max=minus_inf;
        for (int k=0;k<K;++k) {
            for (int l=k+1;l<K;++l) {
                gen D=Delta(k,l);
                _DB_min=min(_DB_min,D,ctx);
                _DB_max=max(_DB_max,D,ctx);
            }
        }
    }
    return dir<0?_DB_min:_DB_max;
}
gen cluster_crit::banfeld_raftery() {
    gen ret=0;
    for (int k=0;k<K;++k) {
        if (is_zero(WGSS(k),ctx))
            return undef;
        ret+=n(k)*ln(WGSS(k)/n(k),ctx);
    }
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::calinski_harabasz() {
    if (K<=1)
        return undef;
    gen ret=((N-K)*BGSS())/((K-1)*WGSS());
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::davies_bouldin() {
    vecteur M(K,minus_inf);
    gen m;
    for (int k=0;k<K;++k) {
        for (int l=0;l<K;++l) {
            if (k!=l) M[k]=max((WGSS(k)/n(k)+WGSS(l)/n(l))/Delta(k,l),M[k],ctx);
        }
    }
    gen ret=_mean(M,ctx);
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::pbm() {
    gen ret=pow(ET*DB(1)/(K*WGSS()),2);
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::ratkowsky_lance() {
    gen Rbar=0;
    for (int j=0;j<p;++j) Rbar+=BGSS(j)/TSS(j);
    gen ret=sqrt(Rbar/p,ctx);
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::ray_turi() {
    gen ret=WGSS()/(N*pow(DB(-1),2));
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::scott_symons() {
    gen ret=0;
    for (int k=0;k<K;++k) {
        gen d=_det(WG(k),ctx);
        if (is_zero(d,ctx))
            return undef;
        ret+=n(k)*log(d/n(k),ctx);
    }
    return is_real_number(ret,ctx)?to_real_number(ret,ctx):undef;
}
gen cluster_crit::det_ratio() {
    gen den=_det(WG(),ctx);
    if (is_zero(den,ctx))
        return undef;
    return _det(T(),ctx)/den;
}
gen cluster_crit::log_det_ratio() {
    gen num=_det(T(),ctx),den=_det(WG(),ctx);
    if (is_zero(num,ctx) || is_zero(den,ctx))
        return undef;
    return N*ln(num/den,ctx);
}
gen cluster_crit::ksq_detW() {
    return K*K*_det(WG(),ctx);
}
gen cluster_crit::log_ss_ratio() {
    gen num=BGSS(),den=WGSS();
    if (is_zero(num,ctx) || is_zero(den,ctx))
        return undef;
    return ln(num/den,ctx);
}
gen cluster_crit::trace_W() {
    return WGSS();
}
gen cluster_crit::trace_WiB() {
    gen C;
    try {
        C=_evalf(symbolic(at_trace,mmult(*_inv(WG(),ctx)._VECTptr,BG())),ctx);
    } catch (const std::runtime_error &e) {
        return undef;
    }
    return is_real_number(C,ctx)?to_real_number(C,ctx):undef;
}
gen cluster_crit::ball_hall() {
    gen ret=0;
    for (int k=0;k<K;++k) ret+=WGSS(k)/n(k);
    return ret/K;
}
void cluster_crit::update(int pos,const gen &val,int dir) {
    assert(dir!=0);
    _vmap[pos].push_back(val);
    if (is_undef(val))
        return;
    bool elbow=(dir==-2||dir==2);
    gen V=val;
    if (elbow) {
        const vecteur &v=_vmap[pos];
        int vs=v.size();
        if (vs<3)
            return;
        V=v[vs-1]-2*v[vs-2]+v[vs-3];
    }
    map<int,pair<int,double> >::const_iterator it=_optvals.find(pos);
    if (!is_real_number(V,ctx))
        return;
    double Vd=to_real_number(V,ctx).to_double(ctx);
    if (it==_optvals.end() || (dir<0 && it->second.second>Vd) || (dir>0 && it->second.second<Vd)) {
        _optvals[pos]=make_pair(elbow?K-1:K,Vd);
        _updated=true;
    }
}
bool cluster_crit::compute_indices(int crit) {
    int pos=0;
    _updated=false;
    if (crit & _BANFELD_RAFTERY_INDEX)     update(pos++,banfeld_raftery(),-1);
    if (crit & _DAVIES_BOULDIN_INDEX)      update(pos++,davies_bouldin(),-1);
    if (crit & _RAY_TURI_INDEX)            update(pos++,ray_turi(),-1);
    if (crit & _SCOTT_SYMONS_INDEX)        update(pos++,scott_symons(),-1);
    if (crit & _CALINSKI_HARABASZ_INDEX)   update(pos++,calinski_harabasz(),1);
    if (crit & _PBM_INDEX)                 update(pos++,pbm(),1);
    if (crit & _RATKOWSKY_LANCE_INDEX)     update(pos++,ratkowsky_lance(),1);
    if (crit & _DET_RATIO_INDEX)           update(pos++,det_ratio(),-2);
    if (crit & _LOG_DET_RATIO_INDEX)       update(pos++,log_det_ratio(),-2);
    if (crit & _LOG_SS_RATIO_INDEX)        update(pos++,log_ss_ratio(),-2);
    if (crit & _BALL_HALL_INDEX)           update(pos++,ball_hall(),2);
    if (crit & _KSQ_DETW_INDEX)            update(pos++,ksq_detW(),2);
    if (crit & _TRACE_W_INDEX)             update(pos++,trace_W(),2);
    if (crit & _TRACE_WIB_INDEX)           update(pos++,trace_WiB(),2);
    return _updated;
}
int cluster_crit::name2index(const string &name) {
    if (name=="banfeld-raftery")                return _BANFELD_RAFTERY_INDEX;
    if (name=="davies-bouldin")                 return _DAVIES_BOULDIN_INDEX;
    if (name=="ray-turi")                       return _RAY_TURI_INDEX;
    if (name=="scott-symons")                   return _SCOTT_SYMONS_INDEX;
    if (name=="calinski-harabasz")              return _CALINSKI_HARABASZ_INDEX;
    if (name=="pbm")                            return _PBM_INDEX;
    if (name=="ratkowsky-lance")                return _RATKOWSKY_LANCE_INDEX;
    if (name=="det")                            return _DET_RATIO_INDEX;
    if (name=="log-det")                        return _DET_RATIO_INDEX;
    if (name=="ksq-detW")                       return _KSQ_DETW_INDEX;
    if (name=="log-ss")                         return _LOG_SS_RATIO_INDEX;
    if (name=="ball-hall")                      return _BALL_HALL_INDEX;
    if (name=="trace-W")                        return _TRACE_W_INDEX;
    if (name=="trace-WiB")                      return _TRACE_WIB_INDEX;
    if (name=="all")                            return _ALL_CLUSTER_INDICES;
    return -1;
}
int cluster_crit::optimal_number_of_clusters(const map<int,pair<int,double> > &optvals) {
    map<int,int> sugg;
    for (map<int,pair<int,double> >::const_iterator it=optvals.begin();it!=optvals.end();++it) {
        sugg[it->second.first]++;
    }
    int best_k=0,num_k;
    for (map<int,int>::const_iterator it=sugg.begin();it!=sugg.end();++it) {
        if (best_k==0 || it->second>num_k) {
            best_k=it->first;
            num_k=it->second;
        }
    }
    return best_k;
}
/* END OF CLUSTER_CRIT CLASS */

/*
 ****** IMPLEMENTATION OF HCLUST CLASS *******
 *
 * Hierarchical agglomerative clustering algorithm
 * 
 * SOURCE: Daniel Muellner, "Modern hierarchical, agglomerative clustering algorithms" (2011)
 * https://arxiv.org/pdf/1109.2378.pdf
 */
hclust::hclust(const matrice &data,const gen &dist_func,GIAC_CONTEXT)
    : ctx(contextptr),points(data),distf(dist_func)
{
    N=data.size();
    assert(N>0);
    S.resize(N);
    vector<int>::iterator it=S.begin(),itend=S.end();
    for (int i=0;it!=itend;++it,++i) *it=i;
    _dist_cache=NULL;
}
hclust::~hclust() {
    if (_dist_cache!=NULL)
        delete[] _dist_cache;
}
double hclust::dist(int i,int j,bool cache) {
    if (i==j)
        return 0;
    gen d;
    if (cache) {
        int I=i<j?i:j,J=i<j?j:i;
        double &dc=dist_cache(I,J);
        if (dc<0) {
            d=comp_dist(points[i],points[j],distf,ctx);
            if (is_undef(d))
                return -1;
            dc=d.to_double(ctx);
        }
        return dc;
    }
    d=comp_dist(points[i],points[j],distf,ctx);
    return is_undef(d)?-1:d.to_double(ctx);
}
double hclust::formula(int a,int b,int x,int na,int nb,int nx,double dab,int meth) {
    double dax=dist(a,x,true),dbx=dist(b,x,true);
    if (dax<0 || dbx<0 || dab<0 || dax+dbx<dab)
        return -1; // bad distance function
    switch (meth) {
    case _SINGLE_LINKAGE:
        return std::min(dax,dbx);
    case _COMPLETE_LINKAGE:
        return std::max(dax,dbx);
    case _AVERAGE_LINKAGE:
        return (na*dax+nb*dbx)/double(na+nb);
    case _WEIGHTED_LINKAGE:
        return (dax+dbx)/2.0;
    case _WARD_LINKAGE:
        return ((na+nx)*dax+(nb+nx)*dbx-nx*dab)/(na+nb+nx);
    case _CENTROID_LINKAGE:
        return (na*dax+nb*dbx)/(na+nb)-na*nb*dab/std::pow(na+nb,2);
    case _MEDIAN_LINKAGE:
        return (dax+dbx)/2.0-dab/4.0;
    default:
        break;
    }
    return giac::nan(); // should not be reachable
}
int hclust::k_low(const dendrogram &dg) {
    double dprev=0,max_step=-1,d;
    int k=0,i;
    for (dendrogram::const_iterator it=dg.begin();it!=dg.end();++it) {
        i=it-dg.begin();
        d=it->first;
        if (d-dprev>max_step) {
            max_step=d-dprev;
            k=N-i;
        }
        dprev=d;
    }
    return k;
}
bool hclust::mst_linkage(dendrogram &dg,bool cache_distances) {
    int c=S.front();
    dg.clear();
    dg.reserve(N-1);
    double *D;
    try { D=new double[2*N]; } catch (const std::bad_alloc &ba) { return false; }
    bool prev_first=true;
    for (vector<int>::const_iterator it=S.begin();it!=S.end();++it) {
        if (*it!=c)
            D[*it]=DBL_MAX;
    }
    for (int i=1;i<N;++i) {
        double minD=DBL_MAX,Dval;
        int n=-1,offset=prev_first?0:N;
        S.erase(std::find(S.begin(),S.end(),c));
        for (vector<int>::const_iterator it=S.begin();it!=S.end();++it) {
            Dval=std::min(D[*it+offset],dist(*it,c,cache_distances));
            if (Dval<0) {
                delete[] D;
                return false;
            }
            if (Dval<minD) {
                n=*it;
                minD=Dval;
            }
            D[*it+N-offset]=Dval;
        }
        assert(n>=0);
        drow r;
        r.first=D[n+N-offset];
        r.second=make_pair(c,n);
        dg.push_back(r);
        c=n;
        prev_first=!prev_first;
    }
    delete[] D;
    return true;
}
bool hclust::nn_chain_linkage(dendrogram &dg,int meth) {
    dg.clear();
    dg.reserve(N-1);
    vector<int> chain;
    chain.reserve(N);
    map<int,double> size;
    vector<int>::const_iterator it;
    vector<int>::iterator et;
    for (it=S.begin();it!=S.end();++it) size[*it]=1;
    int a,b,c,x;
    while (S.size()>1) {
        if (chain.size()<=3) {
            a=S[0];
            chain.resize(1);
            chain[0]=a;
            b=S[1];
        } else {
            a=chain[chain.size()-4];
            b=chain[chain.size()-3];
            chain.resize(chain.size()-3);
        }
        do {
            double mind=DBL_MAX,d;
            for (it=S.begin();it!=S.end();++it) {
                if ((x=*it)!=a && ((d=dist(x,a,true))<mind || (x==b && d==mind))) {
                    c=x;
                    mind=d;
                }
            }
            b=a;
            a=c;
            chain.push_back(a);
        } while (chain.size()<3 || a!=chain[chain.size()-3]);
        double dab=dist(a,b,true);
        if (dab<0)
            return false;
        dg.push_back(make_pair(dab,make_pair(a,b))); // append (d(a,b),a,b) to dendrogram
        et=std::lower_bound(S.begin(),S.end(),b);
        assert(et!=S.end());
        S.erase(et);
        for (it=S.begin();it!=S.end();++it) { // update distances
            if ((x=*it)==a)
                continue;
            double f=formula(a,b,x,size[a],size[b],size[x],dab,meth);
            if (f<0)
                return false;
            dist_cache(x<a?x:a,x<a?a:x)=f;
        }
        size[a]+=size[b];
    }
    return true;
}
int hclust::linkage(dendrogram &dg,vector<int> &ind,int meth,int K,int min_k,int max_k) {
    int dcs=(N*(N-1))/2;
    if (meth!=_SINGLE_LINKAGE || K<0) {
        if (_dist_cache!=NULL)
            delete[] _dist_cache;
        try { _dist_cache=new double[dcs]; } catch (const std::bad_alloc &ba) { return -1; }
        for (int i=0;i<dcs;++i) _dist_cache[i]=-1;
    }
    switch (meth) {
    case _SINGLE_LINKAGE:
        if (!mst_linkage(dg,K<0))
            return -1;
        break;
    case _COMPLETE_LINKAGE:
    case _AVERAGE_LINKAGE:
    case _WEIGHTED_LINKAGE:
    case _WARD_LINKAGE:
        if (!nn_chain_linkage(dg,meth))
            return -1;
        break;
    default:
        // TODO: implement CENTROID and MEDIAN linkage
        assert(false);
    }
    // agglomerate
    std::stable_sort(dg.begin(),dg.end());
    graphe::unionfind UF(N);
    for (int i=0;i<N;++i) UF.make_set(i);
    int k=N,best_k=0,k_lb=std::max(min_k,K<=0?k_low(dg):0),a,b,c,indf=-K;
    double indval;
    ind.resize(N,-1);
    pair<int,double> ep=make_pair(0,0.0);
    set<int> cind;
    if (indf>0) {
        if (meth!=_SINGLE_LINKAGE) { // force recomputation of distances
            for (int i=0;i<dcs;++i) _dist_cache[i]=-1;
        }
        if (indf & _SILHOUETTE_INDEX) {
            slh_a.resize(N);
            slh_b.resize(N);
            slh_s.resize(N);
        }
    }
    if (indf>0) {
        _kvals.clear();
        _kvals.reserve(max_k-min_k);
    }
    for (dendrogram::iterator it=dg.begin();it!=dg.end();++it) {
        if (k==(K>0?K:k_lb))
            break;
        a=it->second.first;
        b=it->second.second;
        --k;
        if (indf>0 && k<=max_k) { // use index function
            int ci,cj;
            if (k<max_k) {
                int ca=UF.find(a),cb=UF.find(b);
                if (indf & _SILHOUETTE_INDEX) {
                    for (int i=0;i<N;++i) {
                        ci=UF.find(i);
                        map<int,pair<int,double> > &bm=slh_b[i];
                        if (ci==ca || ci==cb) {
                            map<int,pair<int,double> >::iterator jt=bm.find(ci==ca?cb:ca);
                            assert(jt!=bm.end());
                            pair<int,double> &p=slh_a[i];
                            p.first+=jt->second.first;
                            p.second+=jt->second.second;
                            bm.erase(jt);
                        }
                    }
                }
                UF.unite(a,b);
                c=UF.find(a);
                assert(c==ca||c==cb);
                if (c==cb) { // swap ca and cb
                    cb=ca;
                    ca=c;
                }
                if (indf & _SILHOUETTE_INDEX) {
                    for (int i=0;i<N;++i) {
                        ci=UF.find(i);
                        map<int,pair<int,double> > &bm=slh_b[i];
                        if (ci!=c) {
                            map<int,pair<int,double> >::iterator jta=bm.find(ca),jtb=bm.find(cb);
                            assert(jta!=bm.end() && jtb!=bm.end());
                            jta->second.first+=jtb->second.first;
                            jta->second.second+=jtb->second.second;
                            bm.erase(jtb);
                        }
                    }
                }
            } else {
                UF.unite(a,b);
                if (indf & _SILHOUETTE_INDEX)
                    std::fill(slh_a.begin(),slh_a.end(),ep);
                for (int i=0;i<N;++i) {
                    ci=UF.find(i);
                    cind.insert(ci);
                    for (int j=i+1;j<N;++j) {
                        cj=UF.find(j);
                        double d=dist(i,j,true);
                        if (ci==cj) {
                            if (indf & _SILHOUETTE_INDEX) {
                                pair<int,double> &p=slh_a[i],&q=slh_a[j];
                                ++p.first;
                                ++q.first;
                                p.second+=d;
                                q.second+=d;
                            }
                        } else {
                            if (indf & _SILHOUETTE_INDEX) {
                                map<int,pair<int,double> > &bv=slh_b[i],&bw=slh_b[j];
                                pair<int,double> &p=bv[cj],&q=bw[ci];
                                ++p.first;
                                ++q.first;
                                p.second+=d;
                                q.second+=d;
                            }
                        }
                    }
                }
            }
            int pos=0;
            bool updated=false;
            if (indf & _SILHOUETTE_INDEX) {
                std::fill(slh_s.begin(),slh_s.end(),ep);
                double ai,bi,minb,si;
                for (int i=0;i<N;++i) {
                    ci=UF.find(i);
                    minb=DBL_MAX;
                    map<int,pair<int,double> > &bv=slh_b[i];
                    pair<int,double> &ap=slh_a[i];
                    ai=ap.first==0?0:ap.second/ap.first;
                    for (map<int,pair<int,double> >::const_iterator jt=bv.begin();jt!=bv.end();++jt) {
                        bi=jt->second.second/jt->second.first;
                        if (bi<minb) minb=bi;
                    }
                    si=ap.first==0?0:(minb-ai)/std::max(minb,ai);
                    pair<int,double> &s=slh_s[ci];
                    ++s.first;
                    s.second+=si;
                }
                indval=DBL_MIN;
                for (vector<pair<int,double> >::const_iterator it=slh_s.begin();it!=slh_s.end();++it) {
                    if (it->first==0)
                        continue;
                    double S=it->second/it->first;
                    if (S>indval) indval=S;
                }
                _vmap[pos].push_back(indval);
                map<int,pair<int,double> >::const_iterator mt=_optvals.find(pos);
                if (mt==_optvals.end() || mt->second.second<indval) {
                    _optvals[pos]=make_pair(k,indval);
                    updated=true;
                }
                pos++;
            }
            _kvals.push_back(k);
            if (updated) {
                int optk=cluster_crit::optimal_number_of_clusters(_optvals);
                if (optk!=best_k) {
                    for (vector<int>::iterator jt=ind.begin();jt!=ind.end();++jt) {
                        *jt=UF.find(jt-ind.begin());
                    }
                    best_k=optk;
                }
            }
        } else UF.unite(a,b);
    }
    map<int,int> ci_map;
    map<int,int>::iterator jt;
    for (vector<int>::iterator it=ind.begin();it!=ind.end();++it) {
        if (*it<0) *it=UF.find(it-ind.begin());
        jt=ci_map.find(*it);
        if (jt==ci_map.end())
            jt=ci_map.insert(make_pair(*it,ci_map.size())).first;
        *it=jt->second;
    }
    return std::max(k,best_k);
}
/* the following methods are adapted from:
 * https://www.centerspace.net/drawing-dendrograms */
const vector<int> &hclust::order(dendrogram &dg) {
    vector<int> parent(2*N-1,-1);
    int pos=0,lab=N;
    for (dendrogram::iterator it=dg.begin();it!=dg.end();++it) {
        int a=it->second.first,b=it->second.second,pa=a,pb=b;
        while (parent[pa]>=0) pa=parent[pa];
        while (parent[pb]>=0) pb=parent[pb];
        it->second=make_pair(pa,pb);
        parent[pa]=parent[pb]=lab;
        ++lab;
    }
    ord.resize(N);
    reorder_children(dg);
    walk_children(dg,pos,N-2);
    return ord;
}
void hclust::walk_children(const dendrogram &dg,int &pos,int row) {
    int node=dg[row].second.first;
    if (node>=N) walk_children(dg,pos,node-N);
    if (node<N) ord[pos++]=node;
    node=dg[row].second.second;
    if (node>=N) walk_children(dg,pos,node-N);
    if (node<N) ord[pos++]=node;
}
void hclust::reorder_children(dendrogram &dg) {
    for (int i=0;i<N-1;++i) {
        int l1=dg[i].second.first;
        int l2=dg[i].second.second;
        if ((l1<N && l2<N) || (l1>=N && l2>=N)) {
            if (l1>l2) swap(dg,i);
        } else if (l1<N && l2>=N) swap(dg,i);
    }
}
void hclust::swap(dendrogram &dg,int row) {
    int tmp=dg[row].second.first;
    dg[row].second.first=dg[row].second.second;
    dg[row].second.second=tmp;
}
int hclust::name2index(const string &name) {
    if (name=="silhouette")     return _SILHOUETTE_INDEX;
    if (name=="mcclain-rao")    return _MCCLAIN_RAO_INDEX;
    if (name=="dunn")           return _DUNN_INDEX;
    if (name=="all")            return _ALL_HCLUST_INDICES;
    return -1;
}
/* end of the centerspace code slice */
/* END OF HCLUST CLASS */

vecteur data2list(const gen &g,GIAC_CONTEXT) {
    vecteur v=gen2vecteur(g);
    for (unsigned i=0;i<v.size();++i) {
        gen &tmp=v[i];
        if (tmp.is_symb_of_sommet(at_pnt))
            tmp=complex2vecteur(remove_at_pnt(tmp),contextptr);
    }
    return v;
}

int lev_size(const gen &g) {
    if (g.type==_VECT)
        return g._VECTptr->size();
    if (g.type==_STRNG)
        return g._STRNGptr->length();
    assert(false);
}

bool lev_is_equal(const gen &g1,const gen &g2,int i,int j) {
    if (g1.type==_VECT && g2.type==_VECT)
        return g1._VECTptr->at(i)==g2._VECTptr->at(j);
    if (g1.type==_STRNG && g2.type==_STRNG)
        return g1._STRNGptr->at(i)==g2._STRNGptr->at(j);
    assert(false);
}

/* Modified from https://en.wikipedia.org/wiki/Levenshtein_distance */
int LevenshteinDistance(const gen &source, const gen &target) {
    int m=lev_size(source),n=lev_size(target),del_cost,ins_cost,sub_cost;
    if (m>n) return LevenshteinDistance(target,source);
    vector<int> v0(n+1),v1(n+1);
    for (int i=0;i<=n;++i) v0[i]=i;
    for (int i=0;i<m;++i) {
        v1[0]=i+1;
        for (int j=0;j<n;++j) {
            del_cost=v0[j+1]+1;
            ins_cost=v1[j]+1;
            sub_cost=v0[j]+(lev_is_equal(source,target,i,j)?0:1);
            v1[j+1]=std::min(std::min(del_cost,ins_cost),sub_cost);
        }
        v0.swap(v1);
    }
    return v0[n];
}

/* Return the Levenshtein distance between two vectors/strings */
gen _levenshtein(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()!=2)
        return gentypeerr(contextptr);
    const gen &arg1=g._VECTptr->front(),&arg2=g._VECTptr->back();
    if ((arg1.type==_VECT && arg2.type==_VECT) || (arg1.type==_STRNG && arg2.type==_STRNG))
        return LevenshteinDistance(arg1,arg2);
    return generrtype("Both arguments must be either lists or strings");
}
static const char _levenshtein_s []="levenshtein";
static define_unary_function_eval (__levenshtein,&_levenshtein,_levenshtein_s);
define_unary_function_ptr5(at_levenshtein,alias_at_levenshtein,&__levenshtein,0,true)

int cyclic_colormap(int c) {
    int col=c%16; // 16 colors are available for distinguishing between clusters (cyclic colormap)
    switch (col) { // make some colors better for displaying
    case 3: col=94; break;
    case 6: col=222; break;
    case 7: col=92; break;
    }
    return col;
}

gen clustering_result(const vecteur &data,int NP,int p,int out_what,int k,int disp,const vector<int> &res,GIAC_CONTEXT) {
    vecteur ret;
    int as=array_start(contextptr);
    char buf[256];
    switch (out_what) {
    case 1: // return list of clusters
        ret.resize(k,undef);
        for (vector<int>::const_iterator it=res.begin();it!=res.end();++it) {
            int i=it-res.begin();
            gen &c=ret[*it];
            if (is_undef(c))
                c=vecteur(1,data[i]);
            else c._VECTptr->push_back(data[i]);
        }
        return change_subtype(ret,_LIST__VECT);    
    case 2: // return cluster indices for data points
        for (vector<int>::const_iterator it=res.begin();it!=res.end();++it) {
            ret.push_back(*it+as);
        }
        return ret;
    case 3: // plot clusters in different colors
        assert(p==2 || p==3);
        sprintf(buf,"%s: %d",gettext("Number of clusters"),k);
        ret.push_back(symb_equal(change_subtype(_TITLE,_INT_PLOT),string2gen(buf,false))); // show title
        ret.push_back(symb_equal(change_subtype(gen(_AXES),_INT_PLOT),0)); // hide axes
        for (int i=0;i<NP;++i) {
            int pt_attrib=disp!=0?disp:(_POINT_POINT | (p==2?_POINT_WIDTH_2:_POINT_WIDTH_4));
            gen attr=change_subtype(pt_attrib | cyclic_colormap(res[i]),_INT_COLOR);
            ret.push_back(symbolic(at_point,makesequence(data[i],symb_equal(at_display,attr))));
        }
        return ret;
    case 4:
        return k;
    default:
        break;
    }
    return undef; // should not be reachable
}

gen _kmeans(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data=data2list(g.type==_VECT && g.subtype==_SEQ__VECT?g._VECTptr->front():g,contextptr);
    if (data.empty())
        return generr("No data found");
    if (!is_fully_numeric(data,num_mask_withint | num_mask_withfrac))
        return generr("Data must be numeric");
    int NP=data.size(),p=data.front().type==_VECT?data.front()._VECTptr->size():1;
    int out_what=2,min_k=2,max_k=NP-1,k=0,maxiter=100,disp=0;
    if (g.subtype==_SEQ__VECT) for (const_iterateur it=g._VECTptr->begin()+1;it!=g._VECTptr->end();++it) { // parse options
        if (is_equal(*it)) {
            const gen &lh=it->_SYMBptr->feuille._VECTptr->front(),&rh=it->_SYMBptr->feuille._VECTptr->back();
            if (lh==at_count) {
                if (rh.is_integer()) {
                    if ((k=rh.val)<=0) {
                        print_error("invalid value for count option",contextptr);
                        return generr("Expected a positive integer");
                    }
                } else if (rh.is_symb_of_sommet(at_interval)) {
                    const gen &lb=rh._SYMBptr->feuille._VECTptr->front();
                    const gen &ub=rh._SYMBptr->feuille._VECTptr->back();
                    if (!lb.is_integer() || !ub.is_integer() || (min_k=lb.val)<=0 || (max_k=ub.val)<min_k) {
                        print_error("Invalid value for count option",contextptr);
                        return generr("Expected a range of positive integers");
                    }
                } else {
                    print_error("invalid value for count option",contextptr);
                    return generrtype("Expected an integer or range");
                }
            } else if (lh==at_count_sup) {
                if (!rh.is_integer()) {
                    print_error("invalid value for count_sup option",contextptr);
                    return generr("Expected a positive integer");
                }
                max_k=std::min(rh.val,NP-1);
            } else if (lh==at_count_inf) {
                if (!rh.is_integer()) {
                    print_error("invalid value for count_inf option",contextptr);
                    return generr("Expected a positive integer");
                }
                min_k=std::max(rh.val,2);
            } else if (lh==at_limit) {
                if (!rh.is_integer() || (maxiter=rh.val)<1) {
                    print_error("invalid value for limit option",contextptr);
                    return generr("Expected a positive integer");
                }
            } else if (lh==at_display) {
                if (!rh.is_integer() && rh.subtype!=_INT_COLOR) {
                    print_error("invalid value for display option",contextptr);
                    return generr("Expected an integer");
                }
                disp=rh.val;
            } else if (lh==at_index) {
                if (rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                    k=(bool)rh.val?-_CALINSKI_HARABASZ_INDEX:0;
                } else if (rh.type==_STRNG) {
                    int cr=cluster_crit::name2index(*rh._STRNGptr);
                    if (cr<0) return generr("Unknown index name");
                    k=-cr;
                } else if (is_string_list(rh)) {
                    int cr_aggreg=0,cr;
                    for (const_iterateur jt=rh._VECTptr->begin();jt!=rh._VECTptr->end();++jt) {
                        cr=cluster_crit::name2index(*jt->_STRNGptr);
                        if (cr<0) return generr("Unknown index name");
                        cr_aggreg|=cr;
                    }
                    k=-cr_aggreg;
                } else {
                    print_error("invalid value for index option",contextptr);
                    return generrtype("Expected a string or list of strings");
                }
            } else if (lh==at_output) {
                if (rh==at_part) out_what=1;                    // partition
                else if (is_mcint(rh,_MAPLE_LIST)) out_what=2;  // list of cluster indices
                else if (rh==at_plot) out_what=3;               // plot clusters (only for 2D/3D data)
                else if (rh==at_count) out_what=4;              // return the number of clusters
                else if (rh==at_index) out_what=5;              // return values of the first index used
                else return generr("Invalid output specification");
            }
        } else if (it->is_integer()) {
            if ((k=it->val)<=0) {
                print_error("invalid number of clusters",contextptr);
                return generr("Expected a positive integer");
            }
        } else return generr("Invalid option");
    }
    if (out_what==3 && p!=2 && p!=3)
        return generrdim("Plot output works only with 2D or 3D data");
    if (out_what==5 && k>=0)
        return generr("No index specified for output");
    max_k=std::min(max_k,NP-1);
    if (min_k>max_k)
        return generr("Invalid k-range");
    if (min_k==max_k)
        k=max_k;
    vector<int> res;
    matrice ind;
    int code=kmeans(data,k,min_k,max_k,maxiter,res,NULL,out_what==5?&ind:NULL,contextptr);
    if (code==1)
        print_error("failed to generate good initial cluster centers",contextptr);
    else if (code==2)
        print_error("maximum number of iterations exceeded",contextptr);
    if (code!=0)
        return vecteur(0);
    if (out_what==5)
        return ind;
    return clustering_result(data,NP,p,out_what,k,disp,res,contextptr);
}
static const char _kmeans_s []="kmeans";
static define_unary_function_eval (__kmeans,&_kmeans,_kmeans_s);
define_unary_function_ptr5(at_kmeans,alias_at_kmeans,&__kmeans,0,true)

gen make_2D_line_segment(const vecteur &P,const vecteur &Q,int style,GIAC_CONTEXT) {
    vecteur attr(1,style),seg=makevecteur(makecomplex(P[0],P[1]),makecomplex(Q[0],Q[1]));
    return pnt_attrib(gen(seg,_GROUP__VECT),attr,contextptr);
}

gen _cluster(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data=data2list(g.type==_VECT && g.subtype==_SEQ__VECT?g._VECTptr->front():g,contextptr);
    if (data.empty())
        return generr("No data found");
    gen dist_func=is_string_list(data)?at_levenshtein:at_longueur2;
    int NP=data.size(),p=data.front().type==_VECT?data.front()._VECTptr->size():1;
    int out_what=2,k=0,min_k=0,max_k=16,lmeth=_SINGLE_LINKAGE,lab=1,indf=0,disp=0;
    if (g.subtype==_SEQ__VECT) for (const_iterateur it=g._VECTptr->begin()+1;it!=g._VECTptr->end();++it) { // parse options
        if (is_equal(*it)) {
            const gen &lh=it->_SYMBptr->feuille._VECTptr->front(),&rh=it->_SYMBptr->feuille._VECTptr->back();
            if (lh==at_count) {
                if (rh.is_integer()) {
                    if ((k=rh.val)<=0) {
                        print_error("invalid value for count option",contextptr);
                        return generr("Expected a positive integer");
                    }
                } else if (rh.is_symb_of_sommet(at_interval)) {
                    const gen &lb=rh._SYMBptr->feuille._VECTptr->front();
                    const gen &ub=rh._SYMBptr->feuille._VECTptr->back();
                    if (!lb.is_integer() || !ub.is_integer() || (min_k=lb.val)<=0 || (max_k=ub.val)<min_k) {
                        print_error("Invalid value for count option",contextptr);
                        return generr("Expected a range of positive integers");
                    }
                } else {
                    print_error("invalid value for count option",contextptr);
                    return generrtype("Expected an integer or range");
                }
            } else if (lh==at_count_sup) {
                if (!rh.is_integer()) {
                    print_error("invalid value for count_sup option",contextptr);
                    return generr("Expected a positive integer");
                }
                max_k=std::min(rh.val,NP-1);
            } else if (lh==at_display) {
                if (!rh.is_integer() && rh.subtype!=_INT_COLOR) {
                    print_error("invalid value for display option",contextptr);
                    return generr("Expected an integer");
                }
                disp=rh.val;
            } else if (lh==at_count_inf) {
                if (!rh.is_integer()) {
                    print_error("invalid value for count_inf option",contextptr);
                    return generr("Expected a positive integer");
                }
                min_k=std::max(rh.val,2);
            } else if (lh==at_type) {
                if (rh.type!=_STRNG) {
                    print_error("invalid value for type option (linkage)",contextptr);
                    return generrtype("Expected a string");
                }
                const string &lt=*rh._STRNGptr;
                if (lt=="single")           lmeth=_SINGLE_LINKAGE;
                else if (lt=="complete")    lmeth=_COMPLETE_LINKAGE;
                else if (lt=="average")     lmeth=_AVERAGE_LINKAGE;
                else if (lt=="weighted")    lmeth=_WEIGHTED_LINKAGE;
                else if (lt=="ward")        lmeth=_WARD_LINKAGE;
#if 0
                else if (lt=="centroid")    lmeth=_CENTROID_LINKAGE;
                else if (lt=="median")      lmeth=_MEDIAN_LINKAGE;
#endif
                else return generr("Unknown linkage type");
            } else if (lh==at_longueur) {
                if (rh.type!=_FUNC && (rh.type!=_SYMB || rh._SYMBptr->sommet!=at_program))
                    return generrtype("Invalid distance function specification");
                dist_func=rh;
            } else if (lh==at_output) { // set desired output
                if (rh==at_part) out_what=1;                    // partition of the set of points
                else if (is_mcint(rh,_MAPLE_LIST)) out_what=2;  // list of cluster indices for each point
                else if (rh==at_plot) out_what=3;               // plot clusters (works only with 2D or 3D euclidean points)
                else if (rh==at_count) out_what=4;              // return the number of clusters
                else if (is_mcint(rh,_GT_TREE)) out_what=5;     // dendrogram drawing
                else if (rh==at_index) out_what=6;              // return values of the last index used
                else return generr("Invalid output specification");
            } else if (lh.is_integer() && lh.subtype==_INT_PLOT && lh.val==_LABELS) {
                if (rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                    lab=(bool)rh.val?1:0;
                else if (rh==at_index)
                    lab=2;
                else return generr("Invalid labels specification");
            } else if (lh==at_index) {
                if (rh.is_integer() && rh.subtype==_INT_BOOLEAN) {
                    indf=(bool)rh.val?_SILHOUETTE_INDEX:0;
                } else if (rh.type==_STRNG) {
                    indf=hclust::name2index(*rh._STRNGptr);
                    if (indf<0) return generr("Unknown index name");
                } else if (is_string_list(rh)) {
                    int cr;
                    indf=0;
                    for (const_iterateur jt=rh._VECTptr->begin();jt!=rh._VECTptr->end();++jt) {
                        cr=hclust::name2index(*jt->_STRNGptr);
                        if (cr<0) return generr("Unknown index name");
                        indf|=cr;
                    }
                } else {
                    print_error("invalid value for index option",contextptr);
                    return generrtype("Expected a string or list of strings");
                }
            } else return generr("Invalid option");
        } else if (it->is_integer()) {
            if ((k=it->val)<=0) {
                print_error("invalid number of clusters",contextptr);
                return generr("Expected a positive integer");
            }
        } else if (it->type==_FUNC) {
            dist_func=*it;
        } else return generr("Invalid option");
    }
    if (out_what==3 && p!=2 && p!=3)
        return generrdim("Plot output works only with 2D or 3D data");
    if (out_what==6 && indf<=0)
        return generr("No index specified for output");
    if (max_k<min_k)
        return generr("Invalid k-range");
    if (max_k==min_k)
        k=min_k;
    hclust hc(data,dist_func,contextptr);
    hclust::dendrogram dg;
    vector<int> indices;
    if ((k=hc.linkage(dg,indices,lmeth,k>0?k:-indf,min_k,max_k))<=0) {
        print_error("linkage method failed",contextptr);
        return generr("Bad distance function");
    }
    if (out_what==6)
        return hc.get_index_values();
    if (out_what==5) { // draw the dendrogram with horizontal orientation (better for labels placement)
        vector<int> ord=hc.order(dg);
        vecteur drawing;
        vecteur pos(2*NP-1,undef);
        int as=array_start(contextptr),ulab=NP;
        for (int i=0;i<NP;++i) {
            int j=ord[i];
            pos[j]=makevecteur(makevecteur(0,i+1),indices[j]);
            if (lab>0) {
                gen label=_string(lab==1?data[j]:gen(j+as),contextptr);
                gen attr=symb_equal(at_display,change_subtype(_QUADRANT2,_INT_COLOR));
                const vecteur &p=*pos[j]._VECTptr->front()._VECTptr;
                drawing.push_back(_legende(makesequence(makecomplex(p[0],p[1]),label,attr),contextptr));
            }
        }
        for (hclust::dendrogram::const_iterator it=dg.begin();it!=dg.end();++it) {
            int a=it->second.first,b=it->second.second;
            double h=it->first;
            const vecteur &A=*pos[a]._VECTptr->front()._VECTptr,&B=*pos[b]._VECTptr->front()._VECTptr;
            int col=pos[a]._VECTptr->back().val,next_ind;
            if (pos[b]._VECTptr->back().val!=col) {
                next_ind=ulab;
                col=50; // light gray
            } else {
                next_ind=col;
                col=cyclic_colormap(col);
            }
            pos[ulab]=makevecteur(makevecteur(h,(A[1]+B[1])/2),next_ind);
            drawing.push_back(make_2D_line_segment(A,makevecteur(h,A[1]),col,contextptr));
            drawing.push_back(make_2D_line_segment(B,makevecteur(h,B[1]),col,contextptr));
            drawing.push_back(make_2D_line_segment(makevecteur(h,A[1]),makevecteur(h,B[1]),col,contextptr));
            ulab++;
        }
        return drawing;
    }
    return clustering_result(data,NP,p,out_what,k,disp,indices,contextptr);
}
static const char _cluster_s []="cluster";
static define_unary_function_eval (__cluster,&_cluster,_cluster_s);
define_unary_function_ptr5(at_cluster,alias_at_cluster,&__cluster,0,true)

/* Return TRUE iff g is a list of strings */
bool is_string_list(const gen &g) {
    if (g.type!=_VECT || g._VECTptr->empty())
        return false;
    for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
        if (it->type!=_STRNG)
            return false;
    }
    return true;
}

// Brent's method for finding local minimizer of a univariate function
// SOURCE: https://people.math.sc.edu/Burkardt/cpp_src/local_min_rc/local_min_rc.html
//
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double local_min_rc ( double &a, double &b, int &status, double value )

//****************************************************************************80
//
//  Purpose:
//
//    LOCAL_MIN_RC seeks a minimizer of a scalar function of a scalar variable.
//
//  Discussion:
//
//    This routine seeks an approximation to the point where a function
//    F attains a minimum on the interval (A,B).
//
//    The method used is a combination of golden section search and
//    successive parabolic interpolation.  Convergence is never much
//    slower than that for a Fibonacci search.  If F has a continuous
//    second derivative which is positive at the minimum (which is not
//    at A or B), then convergence is superlinear, and usually of the
//    order of about 1.324...
//
//    The routine is a revised version of the Brent local minimization
//    algorithm, using reverse communication.
//
//    It is worth stating explicitly that this routine will NOT be
//    able to detect a minimizer that occurs at either initial endpoint
//    A or B.  If this is a concern to the user, then the user must
//    either ensure that the initial interval is larger, or to check
//    the function value at the returned minimizer against the values
//    at either endpoint.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
//    David Kahaner, Cleve Moler, Steven Nash,
//    Numerical Methods and Software,
//    Prentice Hall, 1989,
//    ISBN: 0-13-627258-4,
//    LC: TA345.K34.
//
//  Parameters
//
//    Input/output, double &A, &B.  On input, the left and right
//    endpoints of the initial interval.  On output, the lower and upper
//    bounds for an interval containing the minimizer.  It is required
//    that A < B.
//
//    Input/output, int &STATUS, used to communicate between
//    the user and the routine.  The user only sets STATUS to zero on the first
//    call, to indicate that this is a startup call.  The routine returns STATUS
//    positive to request that the function be evaluated at ARG, or returns
//    STATUS as 0, to indicate that the iteration is complete and that
//    ARG is the estimated minimizer.
//
//    Input, double VALUE, the function value at ARG, as requested
//    by the routine on the previous call.
//
//    Output, double LOCAL_MIN_RC, the currently considered point.
//    On return with STATUS positive, the user is requested to evaluate the
//    function at this point, and return the value in VALUE.  On return with
//    STATUS zero, this is the routine's estimate for the function minimizer.
//
//  Local parameters:
//
//    C is the squared inverse of the golden ratio.
//
//    EPS is the square root of the relative machine precision.
//
{
  static double arg;
  static double c;
  static double d;
  static double e;
  static double eps;
  static double fu;
  static double fv;
  static double fw;
  static double fx;
  static double midpoint;
  static double p;
  static double q;
  static double r;
  static double tol;
  static double tol1;
  static double tol2;
  static double u;
  static double v;
  static double w;
  static double x;
//
//  STATUS (INPUT) = 0, startup.
//
  if ( status == 0 )
  {
    if ( b <= a )
    {
      status = -1;
      return 0;
    }
    c = 0.5 * ( 3.0 - std::sqrt ( 5.0 ) );

    eps = std::sqrt ( r8_epsilon ( ) );
    tol = r8_epsilon ( );

    v = a + c * ( b - a );
    w = v;
    x = v;
    e = 0.0;

    status = 1;
    arg = x;

    return arg;
  }
//
//  STATUS (INPUT) = 1, return with initial function value of FX.
//
  else if ( status == 1 )
  {
    fx = value;
    fv = fx;
    fw = fx;
  }
//
//  STATUS (INPUT) = 2 or more, update the data.
//
  else if ( 2 <= status )
  {
    fu = value;

    if ( fu <= fx )
    {
      if ( x <= u )
      {
        a = x;
      }
      else
      {
        b = x;
      }
      v = w;
      fv = fw;
      w = x;
      fw = fx;
      x = u;
      fx = fu;
    }
    else
    {
      if ( u < x )
      {
        a = u;
      }
      else
      {
        b = u;
      }

      if ( fu <= fw || w == x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
//
//  Take the next step.
//
  midpoint = 0.5 * ( a + b );
  tol1 = eps * fabs ( x ) + tol / 3.0;
  tol2 = 2.0 * tol1;
//
//  If the stopping criterion is satisfied, we can exit.
//
  if ( fabs ( x - midpoint ) <= ( tol2 - 0.5 * ( b - a ) ) )
  {
    status = 0;
    return arg;
  }
//
//  Is golden-section necessary?
//
  if ( fabs ( e ) <= tol1 )
  {
    if ( midpoint <= x )
    {
      e = a - x;
    }
    else
    {
      e = b - x;
    }
    d = c * e;
  }
//
//  Consider fitting a parabola.
//
  else
  {
    r = ( x - w ) * ( fx - fv );
    q = ( x - v ) * ( fx - fw );
    p = ( x - v ) * q - ( x - w ) * r;
    q = 2.0 * ( q - r );
    if ( 0.0 < q )
    {
      p = - p;
    }
    q = fabs ( q );
    r = e;
    e = d;
//
//  Choose a golden-section step if the parabola is not advised.
//
    if (
      ( fabs ( 0.5 * q * r ) <= fabs ( p ) ) ||
      ( p <= q * ( a - x ) ) ||
      ( q * ( b - x ) <= p ) )
    {
      if ( midpoint <= x )
      {
        e = a - x;
      }
      else
      {
        e = b - x;
      }
      d = c * e;
    }
//
//  Choose a parabolic interpolation step.
//
    else
    {
      d = p / q;
      u = x + d;

      if ( ( u - a ) < tol2 )
      {
        d = tol1 * r8_sign ( midpoint - x );
      }

      if ( ( b - u ) < tol2 )
      {
        d = tol1 * r8_sign ( midpoint - x );
      }
    }
  }
//
//  F must not be evaluated too close to X.
//
  if ( tol1 <= fabs ( d ) )
  {
    u = x + d;
  }
  if ( fabs ( d ) < tol1 )
  {
    u = x + tol1 * r8_sign ( d );
  }
//
//  Request value of F(U).
//
  arg = u;
  status = status + 1;

  return arg;
}

/* Find local/global minimum of f in (a,b) using (a modification of) Brent's method */
gen find_minimum(const gen &f,const gen &a,const gen &b,bool global,double eps,int maxiter,GIAC_CONTEXT) {
    if (!is_real_number(a,contextptr) || !is_real_number(b,contextptr) || maxiter<1 || eps<=0)
        return undef;
    bool isfunc=f.type==_FUNC||(f.type==_SYMB && f._SYMBptr->sommet==at_program);
    gen val,x=isfunc?undef:_lname(f,contextptr);
    if (!is_undef(x) && (x.type!=_VECT || x._VECTptr->size()!=1))
        return undef;
    if (x.type==_VECT) x=x._VECTptr->front();
    double xL=a.to_double(contextptr);
    double xU=b.to_double(contextptr);
    if (xL>=xU)
        return undef;
    double arg,min_arg,min_fval=DBL_MAX,value;
    std::queue<pair<double,double> > intq;
    intq.push(make_pair(xL,xU));
    while (!intq.empty()) {
        int status=0,iter=0;
        pair<double,double> bnds=intq.front();
        intq.pop();
        double x1=bnds.first,x2=bnds.second;
        while (iter++<maxiter) {
            arg=local_min_rc(x1,x2,status,value);
            if (status<0)
                return undef;
            if (status==0 || x2-x1<eps)
                break;
            try {
                val=isfunc?_evalf(f(arg,contextptr),contextptr):subst(f,x,arg,false,contextptr);
            } catch (const std::runtime_error &e) {
                *logptr(contextptr) << e.what() << endl;
                return generr("Failed to compute function value");
            }
            if (!is_real_number(val,contextptr))
                return generr("Function values must be real");
            value=to_real_number(val,contextptr).to_double(contextptr);
        }
        if (!global) {
            min_arg=arg;
            continue;
        }
        if (arg<bnds.first+eps||arg>bnds.second-eps) {
            if (value<min_fval) {
                min_arg=arg;
                min_fval=value;
            }
            continue;
        }
        if (value>=min_fval)
            continue;
        min_arg=arg;
        min_fval=value;
        intq.push(make_pair(bnds.first,arg));
        intq.push(make_pair(arg,bnds.second));
    }
    return min_arg;
}

/* Minimize a scalar univariate function on the interval (a,b) using Brent algorithm.
 * Usage: find_minimum(f,a,b,[optargs])
 * Sequence OPTARGS may contain real number EPS in (0,1) and/or a positive integer MAXITER.
 */
gen _find_minimum(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return generrtype("Expected a sequence of arguments");
    const vecteur &args=*g._VECTptr;
    int argc=args.size();
    if (argc<3 || argc>5)
        return generrdim("Wrong number of input arguments");
    const gen &f=args[0];
    gen a=args[1],b=args[2];
    if (!is_real_number(a,contextptr) || !is_real_number(b,contextptr))
        return generr("Range bounds must be real");
    a=to_real_number(a,contextptr);
    b=to_real_number(b,contextptr);
    if (is_greater(a,b,contextptr))
        return generr("Invalid range");
    double eps=std::pow(epsilon(contextptr),2.0/3.0);
    int maxiter=100;
    for (const_iterateur it=args.begin()+3;it!=args.end();++it) {
        if (!is_real_number(*it,contextptr) || is_positive(-*it,contextptr))
            return generr("Expected a positive number");
        if (it->is_integer())
            maxiter=it->val;
        else {
            eps=to_real_number(*it,contextptr).to_double(contextptr);
            if (is_greater(eps,plus_one,contextptr))
                return generr("Expected a real number in (0,1)");
        }
    }
    return find_minimum(f,a,b,true,eps,maxiter,contextptr);
}
static const char _find_minimum_s []="find_minimum";
static define_unary_function_eval (__find_minimum,&_find_minimum,_find_minimum_s);
define_unary_function_ptr5(at_find_minimum,alias_at_find_minimum,&__find_minimum,0,true)

/* Return true if MAT is positive definite,
 * set isnum to true if MAT is fully numeric */
bool is_positive_definite_inplace(vecteur &m,bool isnum,GIAC_CONTEXT) {
    int n=int(m.size()),j,k,l;
    for (j=0;j<n;j++) {
        gen s,d;
        for (l=j;l<n;l++) {
	        s=0;
	        for (k=0;k<j;k++) {
                const gen &kk=m[k][k];
	            if ((isnum && is_positive(-kk,contextptr)) ||
                        (!isnum && !is_definitely_positive(kk,contextptr)))
	                return false;
            	s+=m[l][k]*m[j][k]/kk;
	        }
            d=m[l][j]-s;
            m[l]._VECTptr->at(j)=isnum?d:ratnormal(d,contextptr);
        }
    }
    return true;
}
bool is_positive_definite(const vecteur &mat,int isnum,GIAC_CONTEXT) {
    bool isn=isnum<0?is_numericm(mat):(isnum==1);
#ifdef HAVE_LIBGSL
    if (isn && is_zero(im(mat,contextptr),contextptr)) {
        gsl_matrix *m=matrice2gsl_matrix(*re(mat,contextptr)._VECTptr,contextptr);
        int res=gsl_linalg_cholesky_decomp(m); // should be cholesky_decomp1 in GSL 2.7
        gsl_matrix_free(m);
        return res!=GSL_EDOM;
    }
#endif // HAVE_LIBGSL
    int n=int(mat.size()),i=0;
    vecteur m;
    m.reserve(n);
    for (const_iterateur it=mat.begin();it!=mat.end();++it) {
        const vecteur &row=*it->_VECTptr;
        m.push_back(vecteur(row.begin(),row.begin()+(++i)));
    }
    return is_positive_definite_inplace(m,isn,contextptr);
}
/* Return true if MAT is positive semidefinite.
 * TIME-LIMIT, if positive, is applied when determining LDL decomposition.
 * MAX-TAILLE, if positive, limits the complexity of the elements of block-diagonal matrix D. */
bool is_positive_semidefinite(const matrice &mat,double time_limit,unsigned max_taille,GIAC_CONTEXT) {
    matrice M;
    copy_matrice(mat,M);
    vector<int> perm;
    bool sing;
    log_output_redirect lor(contextptr);
    int mt=is_numericm(mat,num_mask_withfrac|num_mask_withint)?2:0,i,n=M.size();
    bool res=ldl(M,perm,mt,sing,time_limit,contextptr);
    if (!res)
        return false;
    gen din,dli;
    for (i=0;i<n;++i) {
        const gen &di=M[i][i];
        if (i<n-1 && (!(dli=M[0][i+1]).is_integer() || dli.val!=0)) { // 2x2 block
            if (!is_positive_safe(di,false,max_taille,contextptr) ||
                    !is_positive_safe(din=M[i+1][i+1],false,max_taille,contextptr) ||
                    !is_positive_safe(di*din-pow(dli,2),false,max_taille,contextptr))
                return false;
            ++i;
        } else { // 1x1 block
            if (!is_positive_safe(di,false,max_taille,contextptr))
                return false;
        }
    }
    return true;
}

/* Return 1 if g is a positive definite matrix, else return 0 */
gen _isposdef(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (!is_squarematrix(g))
        return generr("Expected a square matrix");
    return is_positive_definite(*g._VECTptr,-1,contextptr)?1:0;
}
static const char _isposdef_s []="isposdef";
static define_unary_function_eval (__isposdef,&_isposdef,_isposdef_s);
define_unary_function_ptr5(at_isposdef,alias_at_isposdef,&__isposdef,0,true)

vecteur symbol_array(const vector<int> &sz,vector<int> &ind,int pos,const vector<int> &hashpos,const string &tmpl,int as,GIAC_CONTEXT) {
    int nd=sz.size();
    vecteur ret(sz[pos]);
    iterateur rt=ret.begin(),rtend=ret.end();
    for (int i=0;rt!=rtend;++rt,++i) {
        ind[pos]=i;
        if (pos==nd-1) {
            string name(tmpl);
            vector<int>::const_reverse_iterator it=hashpos.rbegin(),itend=hashpos.rend(),jt=ind.rbegin();
            for (;it!=itend;++it,++jt) {
                name.replace(*it,1,print_INT_(*jt+as));
            }
            *rt=identificateur(name);
            if (_eval(*rt,contextptr).type!=_IDNT)
                *logptr(contextptr) << gettext("Warning") << ": " << gettext("variable") << " "
                                    << name << " " << gettext("should be purged") << "\n";
        } else *rt=symbol_array(sz,ind,pos+1,hashpos,tmpl,as,contextptr);
    }
    return ret;
}

gen _symbol_array(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return generrtype("Expected a sequence of arguments");
    const vecteur &args=*g._VECTptr;
    if (args.size()<2)
        return generrdim("Too few input arguments");
    if (args.front().type!=_STRNG)
        return generrtype("First argument must be a string");
    string tmpl=*args.front()._STRNGptr;
    vecteur d(args.begin()+1,args.end());
    if (!is_integer_vecteur(d,true) || is_strictly_greater(1,_min(d,contextptr),contextptr))
        return generr("Expected a sequence of positive integers");
    int nd=d.size(),cnt=0,as=array_start(contextptr);
    vector<int> sz=vecteur_2_vector_int(d),ind(nd,0),hashpos(tmpl.size(),-1);
    string::const_iterator st,ststart=tmpl.begin(),stend=tmpl.end();
    char hash_ch='%';
    for (st=ststart;st!=stend;++st) {
        if (*st==hash_ch)
            hashpos[cnt++]=st-ststart;
        else if ((st==ststart && isdigit(*st)) || (!isalnum(*st) && *st!='_'))
            break;
    }
    hashpos.resize(cnt);
    if (st!=stend || (cnt>0 && cnt!=nd))
        return generr("Invalid name template");
    if (cnt==0) {
        string suffix(nd,hash_ch);
        hashpos.resize(nd);
        for (int i=nd-1;i>0;i--) {
            suffix.insert(suffix.begin()+i,'_');
        }
        if (nd>1) tmpl+="_";
        for (int i=0;i<nd;++i) {
            hashpos[i]=tmpl.size()+2*i;
        }
        tmpl+=suffix;
    }
    return symbol_array(sz,ind,0,hashpos,tmpl,as,contextptr);
}
static const char _symbol_array_s []="symbol_array";
static define_unary_function_eval (__symbol_array,&_symbol_array,_symbol_array_s);
define_unary_function_ptr5(at_symbol_array,alias_at_symbol_array,&__symbol_array,0,true)

void PermuSort::sort(vector<int> &p) {
    comparator comp(*this);
    std::stable_sort(p.begin(),p.end(),comp);
}
void sortperm(const vecteur &v,vector<int> &p,GIAC_CONTEXT) {
    int n=v.size(),i=array_start(contextptr);
    if ((int)p.size()!=n) p.resize(n);
    vector<int>::iterator it=p.begin(),itend=p.end();
    for (;it!=itend;++it,++i) *it=i;
    PermuSort ps(v,contextptr);
    ps.sort(p);
}

/* With 1 argument: return the permutation that sorts the given vector.
 * With 2 arguments: first arg must be a vector, the second a permutation,
 *                   then return the first sorted according to the second. */
gen _sortperm(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return gentypeerr(contextptr);
    if (g.subtype==_SEQ__VECT) {
        const vecteur &args=*g._VECTptr;
        if (args.size()!=2)
            return gendimerr(contextptr);
        if (args.front().type!=_VECT || args.back().type!=_VECT)
            return gentypeerr(contextptr);
        vector<int> p;
        if (!is_permu(*args.back()._VECTptr,p,contextptr))
            return gentypeerr(contextptr);
        vecteur v=*args.front()._VECTptr;
        apply_permutation(v,p,false);
        return gen(v,args.front().subtype);    
    }
    const vecteur &v=*g._VECTptr;
    if (v.empty()) return vecteur(0);
    if (v.size()==1) return vecteur(1,array_start(contextptr));
    vector<int> p;
    sortperm(v,p,contextptr);
    return vector_int_2_vecteur(p);
}
static const char _sortperm_s []="sortperm";
static define_unary_function_eval (__sortperm,&_sortperm,_sortperm_s);
define_unary_function_ptr5(at_sortperm,alias_at_sortperm,&__sortperm,0,true)

/*
 * FDWEIGHTS CLASS IMPLEMENTATION
 * Fast and accurate numerical differentiation
 */
FDWeights::FDWeights(const vecteur &grid_points,int diff_order,GIAC_CONTEXT) : ctx(contextptr),z(grid_points) {
    N=grid_points.size();
    M=diff_order;
    w_lagrange.resize(N);
    L=zero_mat(N+2,M+1,ctx);
    R=zero_mat(N+2,M+1,ctx);
    w=zero_mat(N,M+1,ctx);
    lagrange_weights();
}
void FDWeights::lagrange_weights() {
    iterateur it=w_lagrange.begin(),itend=w_lagrange.end();
    const_iterateur zit=z.begin(),zjt,zjtend=z.end();
    for (;it!=itend;++it,++zit) {
        *it=1;
        for (zjt=z.begin();zjt!=zjtend;++zjt) if (zit!=zjt) *it=*it*(*zit-*zjt);
        *it=_inv(*it,ctx);
    }
}
void FDWeights::find_weights(const vecteur &c,int k) {
    gen f=w_lagrange[k];
    vecteur &wk=*w[k]._VECTptr;
    iterateur it=wk.begin(),itend=wk.end();
    const_iterateur ct=c.begin();
    int m=0;
    for (;it!=itend;++it,++ct,++m) {
        *it=*ct*f;
        f=(m+1)*f;
    }
}
inline void FDWeights::multbinom(const vecteur &a,vecteur &b,const gen &zeta) {
    b[0]=-zeta*a[0];
    const_iterateur it=a.begin();
    iterateur jt=b.begin()+1,jtend=b.end();
    for (;jt!=jtend;++it,++jt) *jt=-*(it+1)*zeta+*it;
}
inline void FDWeights::convolve(const vecteur &a,const vecteur &b,vecteur &c) {
    int m=0,i;
    iterateur it=c.begin(),itend=c.end();
    for (;it!=itend;++it,++m) {
        *it=0;
        for (i=0;i<=m;++i) *it+=a[m-i]*b[i];
    }
}
void FDWeights::setz0(const gen &zeta) {
    vecteur zeta_v(N);
    iterateur it=zeta_v.begin(),itend=zeta_v.end();
    const_iterateur jt=z.begin(),jtend;
    for (;it!=itend;++it,++jt) *it=*jt-zeta;
    vecteur &L0=*L[0]._VECTptr,&Rlast=*R[N+1]._VECTptr;
    std::fill(L0.begin(),L0.end(),0);
    std::fill(Rlast.begin(),Rlast.end(),0);
    L0.front()=Rlast.front()=1;
    int k;
    for (k=1;k<N;++k) multbinom(*L[k-1]._VECTptr,*L[k]._VECTptr,zeta_v[k-1]);
    for (k=N+1;k-->2;) multbinom(*R[k+1]._VECTptr,*R[k]._VECTptr,zeta_v[k-1]);
    vecteur c(M+1);
    for (k=1;k<=N;++k) {
        convolve(*L[k-1]._VECTptr,*R[k+1]._VECTptr,c);
        find_weights(c,k-1);
    }
}
/* END OF FDWEIGHTS CLASS */

double siman_energy_func(void *xp) {
    siman_config *x=static_cast<siman_config*>(xp);
    assert(x!=NULL);
    if (!x->is_valid())
        throw std::runtime_error(gettext("Invalid configuration"));
    gen e=x->energy();
    if (e.type!=_DOUBLE_)
        throw std::runtime_error(gettext("energy function should return a real number"));
    return e.DOUBLE_val();
}

double siman_distance_func(void *xp,void *yp) {
    siman_config *x=static_cast<siman_config*>(xp);
    siman_config *y=static_cast<siman_config*>(yp);
    assert(x!=NULL && y!=NULL);
    if (!x->is_valid() || !y->is_valid())
        throw std::runtime_error(gettext("Invalid configuration"));
    gen d=x->distance(*y);
    if (d.type!=_DOUBLE_)
        throw std::runtime_error(gettext("distance function should return a real number"));
    return d.DOUBLE_val();
}

void siman_modify_func(const gsl_rng *r,void *xp,double step_size) {
    siman_config *x=static_cast<siman_config*>(xp);
    assert(x!=NULL);
    if (!x->is_valid())
        throw std::runtime_error(gettext("Invalid configuration"));
    x->modify(step_size);
    if (!x->is_valid())
        throw std::runtime_error(gettext("failed to modify configuration"));
}

void siman_copy_func(void *source,void *dest) {
    siman_config *src=static_cast<siman_config*>(source);
    siman_config *dst=static_cast<siman_config*>(dest);
    assert(src!=NULL && dst!=NULL);
    if (!src->is_valid() || !dst->is_valid())
        throw std::runtime_error(gettext("Invalid configuration"));
    dst->assign(*src);
}

void *siman_copy_constructor(void *xp) {
    siman_config *x=static_cast<siman_config*>(xp),*copy;
    assert(x!=NULL);
    copy=new siman_config(*x);
    return static_cast<void*>(copy);
}

void siman_destroy_func(void *xp) {
    siman_config *x=static_cast<siman_config*>(xp);
    assert(x!=NULL);
    delete x;
}

gen _siman(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
#ifdef HAVE_LIBGSL
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &args=*g._VECTptr;
    if (args.size()<5 || args.size()>7)
        return generrdim("Wrong number of input arguments");
    const gen &x0=args[0],&efunc=args[1],&distfunc=args[2],&modfunc=args[3],&tp=args[4];
    if (!efunc.is_symb_of_sommet(at_program) || efunc._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=1)
        return generrtype("Energy function should accept one argument");
    if (!distfunc.is_symb_of_sommet(at_program) || distfunc._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=2)
        return generrtype("Distance function should accept two arguments");
    if (!modfunc.is_symb_of_sommet(at_program))
        return generrtype("Invalid step function specification");
    int mod_nargs=modfunc._SYMBptr->feuille._VECTptr->front()._VECTptr->size(),mask=num_mask_withfrac|num_mask_withint;
    int n_tries=10,iters_fixed_T=100;
    double step_size=0;
    if (tp.type!=_VECT || tp._VECTptr->size()!=4 || !is_numericv(*tp._VECTptr,mask))
        return generrtype("Invalid specification of cooling schedule");
    const vecteur &cooling_params=*tp._VECTptr;
    if (!is_zero__VECT(*im(cooling_params,contextptr)._VECTptr,contextptr) || !is_strictly_positive(_min(cooling_params,contextptr),contextptr))
        return generr("Cooling schedule parameters must be positive real numbers");
    if (args.size()>5) { // set iteration parameters
        if (args[5].type!=_VECT || args[5]._VECTptr->size()!=2 || !is_integer_vecteur(*args[5]._VECTptr))
            return generrtype("Invalid specification of iteration limits");
        n_tries=args[5]._VECTptr->front().val;
        iters_fixed_T=args[5]._VECTptr->back().val;
        if (n_tries<1)
            return generr("Invalid number of tries per step");
        if (iters_fixed_T<1)
            return generr("Invalid number of iterations per temperature value");
    }
    if (args.size()>6) { // set max step size
        if (!is_real_number(args[6],contextptr) || (step_size=to_real_number(args[6],contextptr).to_double(contextptr))<=0)
            return generr("Invalid maximum step-size specification");
        if (modfunc._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=2)
            return generrdim("Step function should accept two arguments");
    } else if (modfunc._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=1)
        return generrdim("Step function should accept one argument");
    gsl_siman_params_t params;
    params.n_tries=n_tries;
    params.iters_fixed_T=iters_fixed_T;
    params.step_size=step_size;
    params.k=to_real_number(cooling_params[0],contextptr).to_double(contextptr);
    params.t_initial=to_real_number(cooling_params[1],contextptr).to_double(contextptr);
    params.mu_t=to_real_number(cooling_params[2],contextptr).to_double(contextptr);
    params.t_min=to_real_number(cooling_params[3],contextptr).to_double(contextptr);
    siman_config *cfg=new siman_config(x0,efunc,distfunc,modfunc,contextptr);
    gsl_rng *r=NULL;
    try {
        gsl_rng_env_setup();
        r=gsl_rng_alloc(gsl_rng_default);
        gsl_siman_solve(r,static_cast<void*>(cfg),siman_energy_func,siman_modify_func,siman_distance_func,
                        NULL,siman_copy_func,siman_copy_constructor,siman_destroy_func,0,params);
    } catch (const std::runtime_error &e) {
        delete cfg;
        if (r!=NULL) gsl_rng_free(r);
        return generr(e.what(),false);
    }
    gen ret=cfg->obj();
    if (r!=NULL) gsl_rng_free(r);
    delete cfg;
    return ret;
#else
    return generr("GSL is required for simulated annealing");
#endif
}
static const char _siman_s []="simulated_annealing";
static define_unary_function_eval (__siman,&_siman,_siman_s);
define_unary_function_ptr5(at_siman,alias_at_siman,&__siman,0,true)

#ifndef NO_NAMESPACE_GIAC
}
#endif // ndef NO_NAMESPACE_GIAC
