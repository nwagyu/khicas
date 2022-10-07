/*
 * lpsolve.cc
 *
 * Copyright 2021 Luka Marohnić
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
 */

#include "giacPCH.h"
#include "giac.h"
#include "lpsolve.h"
#include "optimization.h"
#include <ctime>
#include <set>

#ifndef DBL_MAX
#define DBL_MAX 1.79769313486e+308
#endif

using namespace std;

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif //ndef NO_NAMESPACE_GIAC

/* 
 * A variant of "is_positive" which hopefully does not hang.
 * Implemented here because optimization.cc includes pari.h which redefines taille!
 */
bool is_positive_safe(const gen &g,bool strict,unsigned max_taille,GIAC_CONTEXT) {
    vecteur s;
    gen alg=g;//to_algebraic(g,s,contextptr); // FIXME
    if (max_taille>0 && taille(alg,max_taille+1)>=max_taille)
        return false;
    bool ret=strict?is_strictly_positive(alg,contextptr):is_positive(alg,contextptr);
    if (!s.empty())
        _purge(s,contextptr);
    return ret;
}

/*
 * Return the number of digits of an unsigned integer.
 */
int numdigits(unsigned i) {
    return 1+(i>0?(int)std::log10((double)i):0);
}

/*
 * Return true iff g is a (vector of) real number(s), +inf, or -inf.
 */
bool is_realcons(const gen &g,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        const vecteur &v = *g._VECTptr;
        const_iterateur it=v.begin(),itend=v.end();
        for (;it!=itend;++it) {
            if (!is_realcons(*it,contextptr))
                return false;
        }
        return true;
    }
    return (is_inf(g) || g.type==_REAL || _evalf(g,contextptr).type==_DOUBLE_);
}

/*
 * If g is an interval, store its bounds to pair p and return true.
 * If g is not an interval, return false.
 */
bool interval2pair(const gen &g,pair<gen,gen> &p,GIAC_CONTEXT) {
    if (g.type!=_SYMB || !g.is_symb_of_sommet(at_interval))
        return false;  //g is not an interval
    const vecteur &v=*g._SYMBptr->feuille._VECTptr;
    p=make_pair(v[0],v[1]);
    return is_realcons(v,contextptr);
}

/*
 * Return true iff g is a vector of identifiers.
 */
bool is_idnt_vecteur(const gen &g,GIAC_CONTEXT) {
    if (g.type!=_VECT)
        return false;
    const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
    for (;it!=itend;++it) {
        if (it->type!=_IDNT)
            return false;
    }
    return true;
}

/*
 * Make a singleton vector with 1 at position j and 0 at other positions.
 */
vecteur singleton(int n,int j,bool negative=false) {
    vecteur v(n,0);
    v[j]=negative?-1:1;
    return v;
}

/*
 * Insert column c in matrix at position j.
 */
void insert_column(matrice &m,const vecteur &c,int j) {
    assert(m.size()==c.size());
    for (int i=m.size();i-->0;) {
        vecteur &row=*m[i]._VECTptr;
        row.insert(j>=0?row.begin()+j:row.end()+j,c[i]);
    }
}

/*
 * Append a column c to matrix m.
 */
void append_column(matrice &m,const vecteur &c) {
    assert(m.size()==c.size());
    for (iterateur it=m.begin();it!=m.end();++it) {
        it->_VECTptr->push_back(c[it-m.begin()]);
    }
}

/*
 * Remove the jth column from the matrix m.
 * If j<0, count from the last column towards the first.
 */
void remove_column(matrice &m,int j) {
    for (iterateur it=m.begin();it!=m.end();++it) {
        it->_VECTptr->erase(j+(j>=0?it->_VECTptr->begin():it->_VECTptr->end()));
    }
}

/*
 * Get the jth column from matrix m.
 * If j<0, count from the last column towards the the first.
 */
vecteur jth_column(const matrice &m,int j) {
    int n=m.front()._VECTptr->size();
    vecteur col(m.size());
    for (const_iterateur it=m.begin();it!=m.end();++it) {
        col[it-m.begin()]=it->_VECTptr->at(j>=0?j:n+j);
    }
    return col;
}

/*
 * Multiply the coefficients in v_orig by LCM of
 * denominators and then divide by their GCD.
 */
void integralize(vecteur &v,GIAC_CONTEXT) {
    vecteur vd;
    for (const_iterateur it=v.begin();it!=v.end();++it) {
        if (!is_zero(*it,contextptr))
            vd.push_back(_denom(*it,contextptr));
    }
    if (vd.empty())
        return;
    v=multvecteur(abs(_lcm(vd,contextptr),contextptr),v);
    v=divvecteur(v,abs(_gcd(v,contextptr),contextptr));
}

void lp_variable::assign(const lp_variable &other) {
    _is_integral=other._is_integral;
    _sign_type=other._sign_type;
    _range=other._range;
    _name=other._name;
    _subs_coef=other._subs_coef;
    pseudocost[0]=other.pseudocost[0];
    pseudocost[1]=other.pseudocost[1];
    nbranch[0]=other.nbranch[0];
    nbranch[1]=other.nbranch[1];
}

/*
 * LP variable constructor. By default, it is a nonnegative variable
 * unrestricted from above.
 */
lp_variable::lp_variable() {
    _is_integral=false;
    _sign_type=_LP_VARSIGN_POS;
    _range=lp_range();
    _range.set_lb(0);
    fill_n(nbranch,2,0);
}

/*
 * Update lower (dir=0) or upper (dir=1) pseudocost.
 */
void lp_variable::update_pseudocost(double delta,double fr,int dir) {
    if (fr==0) return;
    double sigma=pseudocost[dir]*nbranch[dir];
    sigma+=delta/(dir==0?fr:(1-fr));
    pseudocost[dir]=sigma/(++nbranch[dir]);
}

/*
 * Return score, a positive value based on pseudocost values. Variable with the
 * best (highest) score is selected for branching.
 */
double lp_variable::score(double fr) const {
    if (nbranch[0]==0 || nbranch[1]==0)
        return 0;
    double qlo=fr*pseudocost[0],qhi=(1-fr)*pseudocost[1];
    return (1-LP_SCORE_FACTOR)*std::min(qlo,qhi)+LP_SCORE_FACTOR*std::max(qlo,qhi);
}

/* Set the type of this variable. */
void lp_variable::set_type(int t,GIAC_CONTEXT) {
    switch (t) {
    case _LP_BINARYVARIABLES:
        tighten_lbound(0,contextptr);
        tighten_ubound(1,contextptr);
    case _LP_INTEGERVARIABLES:
        _is_integral=true;
        break;
    }
}

/*
 * Range constructor: by default, it contains no restriction.
 */
lp_range::lp_range () {
    lbound=minus_inf;
    ubound=plus_inf;
}

/*
 * Settings constructor loads some sensible defaults.
 */
lp_settings::lp_settings() {
    verbose=false;
    status_report_freq=0.2;
    solver=_LP_SIMPLEX;
    precision=_LP_PROB_DEPENDENT;
    presolve=1;
    maximize=false;
    acyclic=true;
    relative_gap_tolerance=0.0;
    has_binary_vars=false;
    varselect=-1;
    nodeselect=-1;
    depth_limit=0;
    node_limit=0;
    iteration_limit=0;
    time_limit=0;
    max_cuts=5;
    use_heuristic=true;
}

/*
 * Stats constructor initializes the status container for the problem being
 * solved. It is used to monitor the progress and to summarize when done.
 */
lp_stats::lp_stats() {
    subproblems_examined=0;
    cuts_applied=0;
    cut_improvement=0;
    max_active_nodes=0;
    mip_gap=-1; //negative means undefined
}

/*
 * Pivot on element with coordinates I,J.
 */
void lp_node::pivot_ij(matrice &m,int I,int J,bool negate) {
    pivot_elm=m[I][J];
    vecteur &mI=*m[I]._VECTptr;
    iterateur it=mI.begin(),itend=mI.end(),jt;
    vector<intgen>::iterator rcbeg=pivot_row.begin(),rc;
    int i,j,jprev,jmax;
    for (j=jprev=0,rc=rcbeg;it!=itend;++it,++j) {
        if (j!=J && !is_zero(*it,prob->ctx)) {
            *(rc++)=make_pair(j-jprev,*it=*it/pivot_elm);
            jprev=j;
        }
    }
    jmax=rc-rcbeg;
    mI[J]=gen(negate?-1:1)/pivot_elm;
    for (it=m.begin(),itend=m.end(),i=0;it!=itend;++it,++i) {
        vecteur &mi=*it->_VECTptr;
        gen &miJ=mi[J];
        if (i!=I && !is_zero(miJ,prob->ctx)) {
            j=jmax;
            for (rc=rcbeg,jt=mi.begin();j-->0;++rc)
                *(jt+=rc->first)-=miJ*rc->second;
            miJ=gen(negate?1:-1)*miJ/pivot_elm;
        }
    }
}

static clock_t srbt;

/*
 * Change basis by choosing entering and leaving variables and swapping them.
 * Return true iff there is no need to change basis any further.
 */
bool lp_node::change_basis(matrice &m,const vecteur &u,vector<bool> &is_slack,ints &basis,ints &cols) {
    // ev, lv: indices of entering and leaving variables
    // ec, lr: 'entering' column and 'leaving' row in matrix m, respectively
    int ec,ev,lr,lv,nc=cols.size(),nr=basis.size(),i,j;
    gen ratio,mincoeff=0;
    // choose a variable to enter basis
    ev=-1;
    const vecteur &last=*m.back()._VECTptr;
    vector<int>::const_iterator it;
    const_iterateur jt,jtend;
    for (j=0,it=cols.begin(),jt=last.begin();j<nc;++j,++it,++jt) {
        if ((use_bland && !is_positive(*jt,prob->ctx) &&
                (ev<0 || *it+(is_slack[*it]?nv:0)<ev+(is_slack[ev]?nv:0))) ||
                (!use_bland && is_strictly_greater(mincoeff,*jt,prob->ctx))) {
            ec=j;
            ev=*it;
            mincoeff=*jt;
        }
    }
    if (ev<0) // current solution is optimal
        return true;
    // choose a variable to leave basis
    mincoeff=plus_inf;
    lv=-1;
    bool hits_ub,ub_subs;
    for (i=0,jt=m.begin();i<nr;++i,++jt) {
        const gen &a=jt->_VECTptr->at(ec),&b=jt->_VECTptr->back();
        j=basis[i];
        if (is_strictly_positive(a,prob->ctx) && is_greater(mincoeff,ratio=b/a,prob->ctx))
            hits_ub=false;
        else if (!is_positive(a,prob->ctx) && !is_inf(u[j]) &&
                    is_greater(mincoeff,ratio=(b-u[j])/a,prob->ctx))
            hits_ub=true;
        else continue;
        if (is_strictly_greater(mincoeff,ratio,prob->ctx)) {
            lv=-1;
            mincoeff=ratio;
        }
        if (lv<0 || (use_bland && j+(is_slack[j]?nv:0)<lv+(is_slack[lv]?nv:0))) {
            lv=j;
            lr=i;
            ub_subs=hits_ub;
        }
    }
    if (lv<0 && is_inf(u[ev])) { // solution is unbounded
        optimum=minus_inf;
        return true;
    }
    if (prob->settings.acyclic)
        use_bland=is_zero(mincoeff,prob->ctx); // Bland's rule
    if (lv<0 || is_greater(mincoeff,u[ev],prob->ctx)) {
        for (jt=m.begin(),jtend=m.end();jt!=jtend;++jt) {
            gen &a=jt->_VECTptr->at(ec);
            a=-a;
            jt->_VECTptr->back()+=u[ev]*a;
        }
        if (ev<nv)
            is_slack[ev]=!is_slack[ev];
        return false;
    }
    // swap variables: basic leaves, nonbasic enters
    if (ub_subs) {
        m[lr]._VECTptr->back()-=u[lv];
        if (lv<nv)
            is_slack[lv]=!is_slack[lv];
    }
    pivot_ij(m,lr,ec,ub_subs);
    basis[lr]=ev;
    cols[ec]=lv;
    return false;
}

/*
 * Simplex algorithm that handles upper bounds of the variables. The solution x
 * satisfies 0<=x<=u. An initial basis must be provided.
 *
 * Basis is a vector of integers B and B[i]=j means that j-th variable is basic
 * and appears in i-th row (constraint). Basic columns are not kept in matrix,
 * which contains only the columns of nonbasic variables. A nonbasic variable
 * is assigned to the respective column with integer vector 'cols': cols[i]=j
 * means that i-th column of the matrix is associated with the j-th variable. The
 * algorithm uses the upper-bounding technique when pivoting and an adaptation
 * of Bland's rule to prevent cycling. i-th element of 'is_slack' is true iff
 * the i-th (nonbasic) variable is at its upper bound.
 *
 * If limit>0, the simplex algorithm will terminate after that many iterations.
 */
void lp_node::simplex_reduce_bounded(matrice &m,const vecteur &u,vector<bool> &is_slack,
                                     ints &basis,ints &cols,int phase,const gen &obj_ct) {
    int nr=basis.size(),nc=cols.size();
    int limit=prob->settings.iteration_limit;
    int &icount=prob->iteration_count;
    // iterate the simplex method
    use_bland=false;
    optimum=undef;
    char buffer[256],numbuf[8];
    double obj0=-1;
    while (true) {
        ++icount;
        if (limit>0 && icount>limit)
            break;
        if (phase>0) {
            clock_t now=clock();
            double obj=_evalf(phase==1?-m[nr][nc]:(prob->settings.maximize?-1:1)*(obj_ct-m[nr][nc]),prob->ctx).DOUBLE_val();
            if (phase==1 && obj0<obj) {
                if (is_zero(m[nr][nc],prob->ctx))
                    break;
                obj0=obj;
            }
            if (CLOCKS_PER_SEC/double(now-srbt)<=0.5) { // display progress info
                sprintf(numbuf,"%d",icount);
                string ns(numbuf);
                while (ns.length()<7) ns.insert(0,1,' ');
                sprintf(buffer," %s %s  obj: %g%s",phase==1?" ":"*",ns.c_str(),phase==1?100*obj/obj0:obj,phase==1?"%":"");
                prob->message(buffer);
                srbt=now;
            }
        }
        if (change_basis(m,u,is_slack,basis,cols) || (std::abs(phase)==1 && is_zero(m[nr][nc],prob->ctx)))
            break;
    }
    if (is_undef(optimum))
        optimum=m[nr][nc];
}

/*
 * Remove a free (empty-column) variable corresponding to the j-th column.
 */
bool lp_node::remove_free_variable(int j,const vecteur &l,const vecteur &u,ints &cols,const vecteur &obj,gen &obj_ct) {
    gen val;
    int v=cols[j];
    if (is_positive(obj[v],prob->ctx))
        val=l[v];
    else if (is_inf(u[v]))
        return false; // unbounded
    else val=u[v];
    vector<intgen> rcol(1,make_pair(v,val));
    removed_cols.push(rcol);
    obj_ct+=obj[v]*val;
    cols.erase(cols.begin()+j);
    return true;
}

/*
 * Preprocess a B&B-tree node.
 */
int lp_node::preprocess(matrice &m,vecteur &bv,vecteur &l,vecteur &u,ints &cols,vecteur &obj,gen &obj_ct) {
    bool changed;
    do {
        changed=false;
        for (int i=m.size();i-->0;) { // detect empty and singleton rows
            const vecteur &a=*m[i]._VECTptr;
            const gen &b=bv[i];
            int nz=0,j0=-1;
            for (int j=cols.size();j-->0;) {
                if (!is_zero(a[j],prob->ctx)) {
                    ++nz;
                    j0=j;
                }
            }
            if (nz==0) { // empty row detected
                if (!is_zero(b,prob->ctx))
                    return _LP_INFEASIBLE;
            } else if (nz==1) { // singleton row detected
                assert(j0>=0);
                int v=cols[j0];
                if (is_strictly_greater(b/a[j0],u[v],prob->ctx) || is_strictly_greater(l[v],b/a[j0],prob->ctx))
                    return _LP_INFEASIBLE;
                u[v]=l[v]=b/a[j0];
            } else continue;
            m.erase(m.begin()+i);
            bv.erase(bv.begin()+i);
            if (m.empty()) return 0;
            changed=true;
        }
        for (int i=m.size();i-->0;) { // detect forcing constraints
            const vecteur &a=*m[i]._VECTptr;
            const gen &b=bv[i];
            gen g(0),h(0);
            for (int j=cols.size();j-->0;) {
                if (is_zero(a[j],prob->ctx))
                    continue;
                int v=cols[j];
                g+=a[j]*(is_positive(a[j],prob->ctx)?l[v]:u[v]);
                h+=a[j]*(is_positive(a[j],prob->ctx)?u[v]:l[v]);
            }
            if (is_undef(g) || is_undef(h))
                continue;
            if (is_strictly_greater(b,h,prob->ctx) || is_strictly_greater(g,b,prob->ctx))
                return _LP_INFEASIBLE;
            if (is_inf(g) || is_inf(h))
                continue;
            if (is_zero(g-b,prob->ctx) || is_zero(h-b,prob->ctx)) { // forcing constraint detected
                for (int j=cols.size();j-->0;) {
                    if (is_zero(a[j],prob->ctx))
                        continue;
                    int v=cols[j];
                    if (is_positive(a[j]*gen(is_zero(g-b,prob->ctx)?1:-1),prob->ctx))
                        u[v]=l[v];
                    else l[v]=u[v];
                }
                m.erase(m.begin()+i);
                bv.erase(bv.begin()+i);
                if (m.empty()) return 0;
                changed=true;
                continue;
            }
            // detect free singleton columns, tighten variable bounds
            for (int j=cols.size();j-->0;) {
                if (is_zero(a[j],prob->ctx))
                    continue;
                int v=cols[j];
                gen vl=(is_positive(a[j],prob->ctx)?b-h:b-g)/a[j]+u[v];
                gen vu=(is_positive(a[j],prob->ctx)?b-g:b-h)/a[j]+l[v];
                int k=0;
                for (;k<int(m.size()) && (k==i || is_zero(m[k][j],prob->ctx));++k);
                if (k==int(m.size()) && is_greater(vl,l[v],prob->ctx) && is_greater(u[v],vu,prob->ctx)) {
                    // free column singleton detected
                    vector<intgen> rcol(1,make_pair(v,b/a[j]));
                    obj_ct+=obj[v]*b/a[j];
                    for (k=cols.size();k-->0;) {
                        if (k!=j && !is_zero(a[k],prob->ctx)) {
                            rcol.push_back(make_pair(cols[k],-a[k]/a[j]));
                            obj[cols[k]]-=obj[v]*a[k]/a[j];
                        }
                    }
                    removed_cols.push(rcol);
                    cols.erase(cols.begin()+j);
                    m.erase(m.begin()+i);
                    bv.erase(bv.begin()+i);
                    if (m.empty()) return 0;
                    remove_column(m,j);
                    if (m.front()._VECTptr->empty()) return 0;
                    changed=true;
                    break;
                } else { // tighten bounds on var
                    l[v]=max(l[v],vl,prob->ctx);
                    u[v]=min(u[v],vu,prob->ctx);
                }
            }
        }
        // integralize bounds, check for impossible bounds
        for (ints::const_iterator it=cols.begin();it!=cols.end();++it) {
            if (prob->variables[*it].is_integral()) {
                l[*it]=_ceil(l[*it],prob->ctx);
                u[*it]=_floor(u[*it],prob->ctx);
            }
            if (is_strictly_greater(l[*it],u[*it],prob->ctx))
                return _LP_INFEASIBLE;
        }
        // remove empty columns and fixed variables
        for (int j=cols.size();j-->0;) {
            int v=cols[j];
            if (is_zero__VECT(jth_column(m,j),prob->ctx)) { // j-th column is empty
                if (!remove_free_variable(j,l,u,cols,obj,obj_ct))
                    return _LP_UNBOUNDED;
                remove_column(m,j);
                if (m.front()._VECTptr->empty()) return 0;
                changed=true;
            } else if (is_zero(u[v]-l[v],prob->ctx)) { // variable v is fixed
                if (!is_zero(l[v]),prob->ctx) for (int i=bv.size();i-->0;) {
                    bv[i]-=m[i][j]*l[v];
                }
                vector<intgen> rcol(1,make_pair(v,l[v]));
                removed_cols.push(rcol);
                obj_ct+=obj[v]*l[v];
                cols.erase(cols.begin()+j);
                remove_column(m,j);
                if (m.front()._VECTptr->empty()) return 0;
                changed=true;
            }
        }
    } while (changed);
    return 0;
}

/*
 * Solve the relaxed subproblem corresponding to this node.
 *
 * This function uses two-phase simplex method and applies suitable Gomory
 * mixed integer cuts generated after (each re)optimization. Weak GMI cuts are
 * discarded either because of small away or because not being parallel enough
 * to the objective. Cuts with too large coefficients (when integralized) are
 * discarded too because they slow down the computational process. Generated
 * cuts are kept in the problem structure to be used by child suboproblems
 * during the branch&bound algorithm.
 */
int lp_node::solve_relaxation() {
    int nrows=prob->constr.nrows(),ncols=prob->constr.ncols(),bs,nc;
    matrice m,cuts;
    vecteur obj=prob->objective.first,l(ncols),u(ncols),br,row,b,lh;
    gen rh,obj_ct(prob->objective.second),mgn;
    ints cols(ncols),basis;
    vector<bool> is_slack(ncols,false);
    bool is_mip=prob->has_integral_variables();
    // determine the upper and the lower bound
    for (int j=0;j<ncols;++j) {
        cols[j]=j;
        const lp_variable &var=prob->variables[j];
        const lp_range &rng=ranges[j];
        l[j]=max(var.lb(),rng.lb(),prob->ctx);
        u[j]=min(var.ub(),rng.ub(),prob->ctx);
        if (is_strictly_greater(l[j],u[j],prob->ctx))
            return _LP_INFEASIBLE;
    }
    // populate matrix with constraint coefficients
    m=*_matrix(makesequence(nrows,ncols,0),prob->ctx)._VECTptr;
    for (int i=0;i<nrows;++i) for (int j=0;j<ncols;++j)
        m[i]._VECTptr->at(j)=prob->constr.lhs[i][j];
    b=prob->constr.rhs;
    // preprocess
    if (prob->settings.presolve==1 && is_mip) {
        int res=preprocess(m,b,l,u,cols,obj,obj_ct);
        if (res!=0)
            return res;
        if (m.empty()) for (int j=cols.size();j-->0;) {
            if (!remove_free_variable(j,l,u,cols,obj,obj_ct))
                return _LP_UNBOUNDED;
        }
        nrows=m.size();
    }
    nc=cols.size();
    // shift tableau variables so that l<=x<=u becomes 0<=x'<=u'
    ints shifted;
    for (ints::const_iterator it=cols.begin();it!=cols.end();++it) {
        int i=it-cols.begin(),j=*it;
        if (!is_zero(l[j],prob->ctx)) {
            b=subvecteur(b,multvecteur(l[j],jth_column(m,i)));
            u[j]-=l[j];
            obj_ct+=obj[j]*l[j];
            shifted.push_back(j);
        }
    }
    // append b to the tableau
    if (!m.empty())
        append_column(m,b);
    // assure that the right-hand side column is nonnegative
    iterateur it=m.begin(),itend=m.end();
    for (;it!=itend;++it) {
        if (!is_positive(it->_VECTptr->back(),prob->ctx)) {
            *it=multvecteur(-1,*(it->_VECTptr));
        }
    }
    // optimize and cut, repeat
    srbt=clock();
    int pass=0,cuts0;
    gen opt1;
    double imp=0;
    if (nc>0) while (true) { // apply two-phase simplex algorithm
        ++pass;
        nv=ncols+nrows;
        nc=cols.size();
        if (pivot_row.size()<=cols.size())
            pivot_row.resize(nc+1);
        br=vecteur(nc+1,0);
        bs=basis.size();
        basis.resize(nrows);
        u.resize(ncols+nrows-bs);
        for (int i=bs;i<nrows;++i) {
            br=subvecteur(br,*m[i]._VECTptr);
            basis[i]=ncols+i-bs;
            u[ncols+i-bs]=plus_inf;
        }
        m.push_back(br);
        assert(ckmatrix(m));
        // phase 1: minimize the sum of artificial variables
        simplex_reduce_bounded(m,u,is_slack,basis,cols,is_mip || !prob->settings.verbose?-1:1,obj_ct);
        if (!is_zero(optimum,prob->ctx))
            return _LP_INFEASIBLE; // at least one artificial variable is basic and positive
        m.pop_back();
        // push artificial variables out of the basis (requires that m is full-rank)
        for (int i=0;i<nrows;++i) {
            int j=basis[i];
            if (j<ncols)
                continue;
            int k=0;
            for (;k<nc && (is_zero(m[i][k],prob->ctx) || cols[k]>=ncols);++k);
            if (k==nc)
                return _LP_ERROR;
            pivot_ij(m,i,k);
            basis[i]=cols[k];
            cols[k]=j;
        }
        // remove artificial columns from m
        for (int i=nc;i-->0;) {
            if (cols[i]>=ncols) {
                remove_column(m,i);
                cols.erase(cols.begin()+i);
                --nc;
            }
        }
        // append the bottom row to maximize -obj
        br=vecteur(nc+1,0);
        for (int i=0;i<nc;++i) {
            int j=cols[i];
            br[i]=obj[j];
            if (is_slack[j]) {
                br.back()-=br[i]*u[j];
                br[i]=-br[i];
            }
        }
        for (int i=0;i<nrows;++i) {
            int j=basis[i];
            if (is_slack[j])
                br.back()-=obj[j]*u[j];
            br=subvecteur(br,multvecteur(is_slack[j]?-obj[j]:obj[j],*m[i]._VECTptr));
        }
        m.push_back(br);
        assert(ckmatrix(m));
        u.resize(ncols);
        // phase 2: optimize the objective
        simplex_reduce_bounded(m,u,is_slack,basis,cols,is_mip || !prob->settings.verbose?-2:2,obj_ct);
        if (is_inf(optimum))
            return _LP_UNBOUNDED; // the solution is unbounded
        if (pass==1) {
            opt1=optimum;
            cuts0=prob->stats.cuts_applied;
        }
        else if (!is_zero(opt1,prob->ctx))
            imp=((opt1-optimum)/_abs(opt1,prob->ctx)).to_double(prob->ctx)*100;
        // get the solution
        solution=vecteur(ncols+nrows,0);
        for (int i=0;i<nrows;++i) {
            solution[basis[i]]=m[i]._VECTptr->back();
        }
        for (int j=solution.size();j-->0;) {
            if (is_slack[j])
                solution[j]=u[j]-solution[j];
        }
        if (prob->stats.cuts_applied-cuts0>=prob->settings.max_cuts)
            break;
        m.pop_back(); // remove the bottom row
        // attempt to generate GMI cuts
        cuts.clear();
        for (int i=0;i<nrows;++i) {
            int j0=basis[i];
            vecteur eq(*m[i]._VECTptr);
            gen f0=fracpart(eq.back()),fj,sp(0),eqnorm(0);
            if (j0>=prob->nv() || !prob->variables[j0].is_integral() || is_zero(f0,prob->ctx) ||
                    min(f0,1-f0,prob->ctx).to_double(prob->ctx)<LP_MIN_AWAY)
                continue;
            eq.pop_back();
            for (int k=0;k<nc;++k) {
                int j=cols[k];
                if (j<prob->nv() && prob->variables[j].is_integral()) {
                    fj=fracpart(eq[k]);
                    eq[k]=is_strictly_greater(fj,f0,prob->ctx)?(fj-1)/(f0-1):fj/f0;
                } else eq[k]=eq[k]/(is_positive(eq[k],prob->ctx)?f0:f0-1);
                if (j<prob->nv()) {
                    sp+=eq[k]*obj[j];
                    eqnorm+=eq[k]*eq[k];
                }
            }
            if (is_zero__VECT(eq,prob->ctx))
                continue; // not integer feasible anyway
            if (std::abs(sp.to_double(prob->ctx)/(prob->objective_norm*std::sqrt(eqnorm.to_double(prob->ctx))))<LP_MIN_PARALLELISM)
                continue; // this cut is not parallel enough to the objective
            integralize(eq,prob->ctx);
            mgn=_max(_abs(eq,prob->ctx),prob->ctx);
            if (mgn.to_double(prob->ctx)>LP_MAX_MAGNITUDE)
                continue; // this cut has too large coefficients
            cuts.push_back(eq);
        }
        if (cuts.empty())
            break;
        // append GMI cuts to the simplex tableau and reoptimize
        for (int cc=cuts.size();cc-->0;) {
            insert_column(m,vecteur(nrows,0),-1);
            l.push_back(0);
            u.push_back(plus_inf);
            is_slack.push_back(false);
            obj.push_back(0);
            cols.push_back(ncols++);
        }
        const_iterateur ct=cuts.begin(),itend=cuts.end();
        for (;ct!=itend;++ct) {
            int cc=ct-cuts.begin();
            vecteur cut=mergevecteur(*(ct->_VECTptr),singleton(cuts.size(),cc,true));
            cut.push_back(1);
            m.push_back(cut);
            ++nrows;
            if (++prob->stats.cuts_applied==prob->settings.max_cuts+cuts0)
                break;
        }
        assert(ckmatrix(m));
    }
    pivot_row.clear();
    prob->stats.cut_improvement+=imp;
    // undo shifting variables
    for (ints::const_iterator it=shifted.begin();it!=shifted.end();++it) {
        solution[*it]+=l[*it];
    }
    // restore removed variables
    while (!removed_cols.empty()) {
        const vector<intgen> &rcol=removed_cols.top();
        int v=rcol.front().first;
        for (vector<intgen>::const_iterator it=rcol.begin();it!=rcol.end();++it) {
            solution[v]+=(it==rcol.begin()?1:solution[it->first])*it->second;
        }
        removed_cols.pop();
        if (is_strictly_greater(l[v],solution[v],prob->ctx) || is_strictly_greater(solution[v],u[v],prob->ctx))
            return _LP_INFEASIBLE;
    }
    solution.resize(prob->nv());
    // compute objective value
    if (nc>0)
        optimum=obj_ct-optimum;
    else {
        obj.resize(prob->nv());
        optimum=prob->objective.second+scalarproduct(prob->objective.first,solution,prob->ctx);
    }
    // compute some useful data for branch&bound
    opt_approx=optimum.to_double(prob->ctx);
    infeas=0;
    most_fractional=-1;
    gen p,ifs,max_ifs(-1);
    ints mf_cand;
    for (int i=0;i<prob->nv();++i) {
        if (!prob->variables[i].is_integral())
            continue;
        if (is_zero(p=fracpart(solution[i]),prob->ctx))
            continue;
        ifs=min(p,1-p,prob->ctx);
        fractional_vars[i]=p.to_double(prob->ctx);
        infeas+=ifs;
        if (is_greater(ifs,max_ifs,prob->ctx)) {
            if (is_strictly_greater(ifs,max_ifs,prob->ctx))
                mf_cand.clear();
            mf_cand.push_back(i);
            max_ifs=ifs;
        }
    }
    most_fractional=_rand(vector_int_2_vecteur(mf_cand),prob->ctx).val;
    return _LP_SOLVED;
}

/*
 * Return true iff a feasible solution sol can be obtained by rounding heuristic
 */
bool lp_node::rounding_heuristic(vecteur &sol,gen &cost) const {
    map<int,bool> fv;
    sol=solution;
    for (map<int,double>::const_iterator it=fractional_vars.begin();it!=fractional_vars.end();++it) {
        int i=it->first;
        gen &s=sol[i];
        s=_round(s,prob->ctx);
        fv[i]=is_greater(s,solution[i],prob->ctx);
    }
    double viol;
    int vc,j,last_j=-1,pass=0;
    vecteur solc;
    pair<int,double> vp,min_vp,cur_vp;
    cur_vp=prob->constr.violated_constraints(sol,prob->ctx);
    if (cur_vp.first>0) do {
        solc=sol;
        min_vp=cur_vp;
        j=-1;
        for (map<int,bool>::const_iterator it=fv.begin();it!=fv.end();++it) {
            if (last_j==it->first)
                continue;
            solc[it->first]+=it->second?-1:1;
            vp=prob->constr.violated_constraints(solc,prob->ctx);
            if (vp<min_vp) {
                j=it->first;
                min_vp=vp;
            }
            solc[it->first]+=it->second?1:-1;
        }
        if (j<0)
            break;
        sol[j]+=fv[j]?-1:1;
        fv[j]=!fv[j];
        last_j=j;
    } while ((cur_vp=min_vp).first>0);
    if (cur_vp.first>0)
        return false;
    cost=scalarproduct(sol,prob->objective.first,prob->ctx)+prob->objective.second;
    return true;
}

double lp_node::get_fractional_var(int index) const {
    if (fractional_vars.find(index)!=fractional_vars.end())
        return fractional_vars.at(index);
    return 0;
}

void lp_node::assign(const lp_node &other) {
    prob=other.prob;
    depth=other.depth;
    ranges=other.ranges;
    optimum=other.optimum;
    solution=other.solution;
    opt_approx=other.opt_approx;
    infeas=other.infeas;
    most_fractional=other.most_fractional;
    fractional_vars=other.fractional_vars;
    removed_cols=other.removed_cols;
    pivot_row.resize(other.pivot_row.size());
}

/*
 * Return the fractional part of g, i.e. [g]=g-floor(g). It is always 0<=[g]<1.
 */
gen lp_node::fracpart(const gen &g) const {
    return g-_floor(g,prob->ctx);
}

/*
 * Initialize child node with copy of ranges and depth increased by one.
 */
void lp_node::init_child(lp_node &child) {
    child.depth=depth+1;
    child.ranges=this->ranges;
}

/*
 * Return true iff this is a (mixed) integer problem.
 */
bool lp_problem::has_integral_variables() {
    for (vector<lp_variable>::const_iterator it=variables.begin();it!=variables.end();++it) {
        if (it->is_integral())
            return true;
    }
    return false;
}

/*
 * Return true iff the problem has floating-point coefficients.
 */
bool lp_problem::has_approx_coefficients() {
    if (is_approx(objective.first) ||
            objective.second.is_approx() ||
            is_approx(constr.lhs) ||
            is_approx(constr.rhs))
        return true;
    for (vector<lp_variable>::const_iterator it=variables.begin();it!=variables.end();++it) {
        if (it->lb().is_approx() || it->ub().is_approx())
            return true;
    }
    return false;
}

/*
 * Set the objective function parameters.
 */
void lp_problem::set_objective(const vecteur &v,const gen &ft) {
    objective.first=v;
    objective.second=ft;
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it) {
        obj_approx.push_back(abs(*it,ctx).to_double(ctx));
    }
}

/*
 * Display a message.
 */
void lp_problem::message(const char *msg,int type) {
    if (type==1 || settings.verbose) {
        switch (type) {
        case 0:
            break;
        case 1:
            *logptr(ctx) << gettext("Error") << ": ";
            break;
        case 2:
            *logptr(ctx) << gettext("Warning") << ": ";
        }
        *logptr(ctx) << gettext(msg) << "\n";
    }
}

/*
 * Duplicate the jth column by inserting a copy to the position j.
 */
void lp_constraints::duplicate_column(int index) {
    assert(index<ncols());
    vecteur col=jth_column(lhs,index);
    insert_column(lhs,col,index);
}

/*
 * Change signs of the coefficients in the jth column.
 */
void lp_constraints::negate_column(int index) {
    for (int i=0;i<nrows();++i) {
        vecteur &lh=*lhs[i]._VECTptr;
        lh[index]=-lh[index];
    }
}

/*
 * Subtract v from rhs column of constraints.
 */
void lp_constraints::subtract_from_rhs_column(const vecteur &v) {
    assert(int(v.size())==nrows());
    for (int i=0;i<nrows();++i) {
        rhs[i]-=v[i];
    }
}

/*
 * Append the constraint lh<rel>rh.
 */
void lp_constraints::append(const vecteur &lh,const gen &rh,int relation_type) {
    assert(nrows()==0 || int(lh.size())==ncols());
    lhs.push_back(lh);
    rhs.push_back(rh);
    rv.push_back(relation_type);
}

/*
 * Set the constraint with specified index.
 */
void lp_constraints::set(int index,const vecteur &lh,const gen &rh,int relation_type) {
    assert(index<nrows());
    lhs[index]=lh;
    rhs[index]=rh;
    rv[index]=relation_type;
}

/*
 * Get left and right side of the constraint with specified index.
 */
void lp_constraints::get_lr(int index,vecteur &lh,gen &rh) {
    assert(index<nrows());
    lh=*lhs[index]._VECTptr;
    rh=rhs[index];
}

/*
 * Get the constraint with specified index.
 */
void lp_constraints::get(int index,vecteur &lh,gen &rh,int &relation_type) {
    get_lr(index,lh,rh);
    relation_type=rv[index];
}

/*
 * Divide the constraint by g.
 */
void lp_constraints::div(int index,const gen &g,GIAC_CONTEXT) {
    assert(index<nrows());
    lhs[index]=divvecteur(*lhs[index]._VECTptr,g);
    rhs[index]=rhs[index]/g;
    if (!is_positive(g,contextptr))
        rv[index]*=-1;
}

/*
 * Subtract the vector [v,g] from the constraint with specified index.
 */
void lp_constraints::subtract(int index,const vecteur &v,const gen &g) {
    assert(index<nrows());
    lhs[index]=subvecteur(*lhs[index]._VECTptr,v);
    rhs[index]-=g;
}

/*
 * Remove the constraint with specified index.
 */
void lp_constraints::remove(int index) {
    lhs.erase(lhs.begin()+index);
    rhs.erase(rhs.begin()+index);
    rv.erase(rv.begin()+index);
}

/*
 * Delete linearly independent equality constraints, if any.
 */
int lp_constraints::remove_linearly_dependent(GIAC_CONTEXT) {
    matrice m;
    ints pos;
    std::set<int> ri,found;
    for (int i=0;i<nrows();++i) {
        if (rv[i]==_LP_EQ) {
            pos.push_back(i);
            m.push_back(mergevecteur(*lhs[i]._VECTptr,vecteur(1,rhs[i])));
        }
    }
    if (!pos.empty()) {
        matrice mf=mtran(*_rref(mtran(*_evalf(m,contextptr)._VECTptr),contextptr)._VECTptr);
        const_iterateur it=mf.begin(),itend=mf.end();
        for (;it!=itend;++it) {
            int i=0;
            for (;i<=ncols() && is_zero(it->_VECTptr->at(i));++i);
            if (i<=ncols() && is_one(it->_VECTptr->at(i)) && found.find(i)==found.end())
                found.insert(i);
            else ri.insert(it-mf.begin());
        }
        for (std::set<int>::const_reverse_iterator it=ri.rbegin();it!=ri.rend();++it) {
            remove(pos[*it]);
        }
    }
    return ri.size();
}

/*
 * Return the number of violated constraints and compute the infeasibility |b-Ax|.
 */
pair<int,double> lp_constraints::violated_constraints(const vecteur &x,GIAC_CONTEXT) {
    assert(int(x.size())>=ncols_orig);
    vecteur xv(x.begin(),x.begin()+ncols_orig);
    vecteur d=*_epsilon2zero(subvecteur(multmatvecteur(lhs_f,xv),rhs_f),contextptr)._VECTptr;
    for (int i=rv_orig.size();i-->0;) {
        if (rv_orig[i]!=_LP_EQ && is_positive(rv_orig[i]*d[i],contextptr))
            d[i]=0;
    }
    int vc=int(d.size())-_count_eq(makesequence(0,d),contextptr).val;
    double viol=_l2norm(d,contextptr).to_double(contextptr);
    if (viol<LP_FEAS_TOL) vc=0;
    return make_pair(vc,viol);
}

/*
 * Return the total numer of nonzeros in the constraint matrix.
 */
int lp_constraints::nonzeros() {
    int nz=0;
    for (int i=0;i<nrows();++i) {
        for (int j=0;j<ncols();++j) {
            if (!is_zero(lhs[i][j]))
                ++nz;
        }
        if (!is_zero(rhs[i]))
            ++nz;
    }
    return nz;
}

/*
 * Add identifiers from g to variable_identifiers.
 */
void lp_problem::add_identifiers_from(const gen &g) {
    vecteur vars(*_lname(g,ctx)._VECTptr);
    const_iterateur it=vars.begin(),itend=vars.end();
    for (;it!=itend;++it) {
        if (!contains(variable_identifiers,*it))
            variable_identifiers.push_back(*it);
    }
}

/*
 * Return variable index corresponding to the given identifier.
 */
int lp_problem::get_variable_index(const identificateur &idnt) {
    int n=variable_identifiers.size();
    for (int i=0;i<n;++i) {
        if (*variable_identifiers[i]._IDNTptr==idnt)
            return i;
    }
    return -1;
}

/*
 * Create the problem variables, continuous and unrestricted by default.
 */
void lp_problem::create_variables(int n) {
    int m=variables.size();
    if (m>=n)
        return;
    variables.resize(n);
    for (int i=m;i<n;++i) {
        lp_variable var;
        var.set_lb(minus_inf);
        variables[i]=var;
    }
}

/*
 * Tighten both upper and lower bound of the i-th variable.
 */
void lp_problem::tighten_variable_bounds(int i,const gen &l,const gen &u) {
    lp_variable &var=variables[i];
    var.tighten_lbound(l,ctx);
    var.tighten_ubound(u,ctx);
}

/*
 * Output solution in form [x1=v1,x2=v2,...,xn=vn] where xk are variable
 * identifiers and vk are solution values.
 */
vecteur lp_problem::output_solution(bool sort_vars) {
    if (variable_identifiers.empty())
        return solution;
    vecteur v=variable_identifiers;
    if (sort_vars) {
        vecteur sv=sort_identifiers(v,ctx),ret;
        ret.reserve(solution.size());
        const_iterateur it=sv.begin(),itend=sv.end();
        for (;it!=itend;++it) {
            const_iterateur jt=std::find(v.begin(),v.end(),*it);
            assert(jt!=v.end());
            ret.push_back(symb_equal(*it,solution[jt-v.begin()]));
        }
        return ret;
    }
    return *_zip(makesequence(at_equal,v,solution),ctx)._VECTptr;
}

/*
 * Determine coeffcients of linear combination g of variables x in
 * variable_identifiers. varcoeffs C and freecoeff c are filled such that
 * g=C*x+c.
 */
bool lp_problem::lincomb_coeff(const gen &g,vecteur &varcoeffs,gen &freecoeff) {
    gen e(g),a;
    varcoeffs.clear();
    const_iterateur it=variable_identifiers.begin(),itend=variable_identifiers.end();
    for (;it!=itend;++it) {
        a=0;
        if (is_constant_wrt(e,*it,ctx) || (is_linear_wrt(e,*it,a,e,ctx) && is_realcons(a,ctx)))
            varcoeffs.push_back(a);
        else return false;
    }
    return is_realcons(freecoeff=e,ctx);
}

/*
 * Add slack variables to the problem if necessary (i.e. convert all
 * inequalities to equalities).
 */
void lp_problem::add_slack_variables() {
    ints posv;
    constr.ncols_orig=constr.ncols();
    constr.rv_orig=constr.rv;
    constr.lhs_f=*_evalf(constr.lhs,ctx)._VECTptr;
    constr.rhs_f=*_evalf(constr.rhs,ctx)._VECTptr;
    for (int i=0;i<nc();++i) {
        if (constr.rv[i]==_LP_EQ)
            continue;
        append_column(constr.lhs,singleton(nc(),i,constr.rv[i]>0));
        constr.rv[i]=_LP_EQ;
        variables.push_back(lp_variable()); //add slack variable
        posv.push_back(i);
    }
    objective.first.resize(nv(),0);
#if 0
    //determine types of slack variables
    vecteur lh;
    gen rh;
    for (int k=0;k<int(posv.size());++k) {
        int i=posv[k],j=nv0+k;
        lp_variable &var=variables[j];
        constr.get_lr(i,lh,rh);
        if (is_exact(lh) && is_exact(rh)) {
            gen den(_denom(rh,ctx));
            for (int l=0;l<nv0;++l) {
                if (is_zero(lh[l]))
                    continue;
                if (!variables[l].is_integral()) {
                    den=undef;
                    break;
                }
                if (!lh[l].is_integer())
                    den=_lcm(makesequence(den,_denom(lh[l],ctx)),ctx);
            }
            if (!is_undef(den)) {
                for (iterateur it=lh.begin();it!=lh.end();++it) {
                    *it=(*it)*den;
                }
                rh=den*rh;
                var.set_integral(true);
            }
        }
    }
#endif
}

/*
 * Make all decision variables bounded below by negating variables unrestricted
 * below and replacing variables unrestricted from both above and below with
 * the difference of two nonnegative variables. This process is reversed after
 * an optimal solution is found.
 */
void lp_problem::make_all_vars_bounded_below() {
    for (int i=nv();i-->0;) {
        lp_variable &var=variables[i];
        if (var.range().is_unrestricted_below()) {
            if (var.range().is_unrestricted_above()) {
                var.set_lb(0);
                lp_variable negvar(var);
                var.set_sign_type(_LP_VARSIGN_POS_PART);
                negvar.set_sign_type(_LP_VARSIGN_NEG_PART);
                variables.insert(variables.begin()+i,negvar);
                objective.first.insert(objective.first.begin()+i,-objective.first[i]);
                constr.duplicate_column(i);
            } else {
                var.set_lb(-var.ub());
                var.set_ub(plus_inf);
                var.set_sign_type(_LP_VARSIGN_NEG);
                objective.first[i]=-objective.first[i];
            }
            constr.negate_column(i);
        }
    }
}

/*
 * Force all problem parameters to be exact.
 */
void lp_problem::make_problem_exact() {
    objective.first=*exact(objective.first,ctx)._VECTptr;
    objective.second=exact(objective.second,ctx);
    constr.lhs=*exact(constr.lhs,ctx)._VECTptr;
    constr.rhs=*exact(constr.rhs,ctx)._VECTptr;
    for (vector<lp_variable>::iterator it=variables.begin();it!=variables.end();++it) {
        it->set_lb(exact(it->lb(),ctx));
        it->set_ub(exact(it->ub(),ctx));
    }
}

/*
 * Report status.
 */
void lp_problem::report_status(const char *msg) {
    char buf[16];
    sprintf(buf,"%d: ",stats.subproblems_examined);
    int nd=numdigits((unsigned)stats.subproblems_examined);
    string str(msg);
    str.insert(0,buf);
    while (nd<8) {
        str.insert(str.begin(),(char)32);
        ++nd;
    }
    message(str.c_str());
}

/*
 * Return true iff there is a variable with lb>ub.
 */
bool lp_problem::has_infeasible_var() const {
    for (vector<lp_variable>::const_iterator it=variables.begin();it!=variables.end();++it) {
        if (is_strictly_greater(it->lb(),it->ub(),ctx))
            return true;
    }
    return false;
}

/*
 * Remove the j-th variable from the problem.
 * Use subs_coef for substitution.
 */
void lp_problem::remove_variable(int j) {
    const lp_variable &var=variables[j];
    vecteur &obj=objective.first;
    vector<intgen>::const_iterator it=var.subs_coef().begin(),itend=var.subs_coef().end();
    if (it!=itend && it->first<0) {
        objective.second+=obj[j]*it->second;
        ++it;
    }
    for (;it!=itend;++it) {
        if (it->first!=j)
            obj[it->first]+=obj[j]*it->second;
    }
    removed_vars.push(var);
    variables.erase(variables.begin()+j);
    remove_column(constr.lhs,j);
    obj.erase(obj.begin()+j);
    removed_cols.push(j);
}

/*
 * Deduce the optimal value val of a (implied) free variable.
 * c is the respective coefficient in the objective.
 */
int lp_variable::find_opt_free(const gen &c,gen &val,GIAC_CONTEXT) const {
    if (is_zero(c)) {
        if (!is_inf(lb())) {
            val=is_integral()?_ceil(lb(),contextptr):lb();
        } else if (!is_inf(ub())) {
            val=is_integral()?_floor(ub(),contextptr):ub();
        } else val=0;
    } else if (is_positive(c,contextptr)) {
        if (!is_inf(lb())) {
            val=is_integral()?_ceil(lb(),contextptr):lb();
        } else return _LP_UNBOUNDED;
    } else {
        if (!is_inf(ub())) {
            val=is_integral()?_floor(ub(),contextptr):ub();
        } else return _LP_UNBOUNDED;
    }
    if (!is_greater(val,lb(),contextptr) || !is_greater(ub(),val,contextptr))
        return _LP_INFEASIBLE;
    return 0;
}

/*
 * Find implied integral variables and do basic preprocessing afterwards.
 */
void lp_problem::find_implied_integers() {
    vecteur cns;
    vector<bool> is_eq;
    int nc=constr.ncols(),nr=constr.nrows();
    for (int i=0;i<nr;++i) {
        vecteur row=*constr.lhs[i]._VECTptr;
        row.push_back(constr.rhs[i]);
        integralize(row,ctx);
        cns.push_back(row);
        is_eq.push_back(constr.rv[i]==_LP_EQ);
    }
    bool changed;
    do {
        changed=false;
        for (int j=0;j<nc;++j) {
            lp_variable &var=variables[j];
            if (var.is_integral())
                continue;
            for (int i=0;i<nr;++i) {
                if (is_one(cns[i]._VECTptr->at(j))) {
                    int k=0;
                    for (;k<nc && (k==j || is_zero(cns[i]._VECTptr->at(k),ctx) || variables[k].is_integral());++k);
                    if (k==nc) {
                        if (is_eq[i]) {
                            var.set_integral(true);
                            var.tighten_lbound(_ceil(var.lb(),ctx),ctx);
                            var.tighten_ubound(_floor(var.ub(),ctx),ctx);
                            imp_int_count++;
                            changed=true;
                            break;
                        }
                    }
                }
            }
        }
    } while (changed);
    if (imp_int_count>0)
    preprocess(false);
}

/*
 * Preprocess the problem.
 * Handles empty, singleton, forcing rows, and singleton columns.
 */
int lp_problem::preprocess(bool find_imp_int) {
    // basic preprocessing
    bool changed;
    do {
        changed=false;
        for (int i=constr.nrows();i-->0;) { // detect empty and singleton rows
            const vecteur &a=*constr.lhs[i]._VECTptr;
            const gen &b=constr.rhs[i];
            int nz=0,j0=-1,r=constr.rv[i];
            for (int j=0;j<constr.ncols();++j) {
                if (!is_zero(a[j],ctx)) {
                    ++nz;
                    j0=j;
                }
            }
            if (nz==0) { // empty row detected
                if ((r==_LP_EQ && !is_zero(b,ctx)) ||
                        (r==_LP_LEQ && !is_positive(b,ctx)) ||
                        (r==_LP_GEQ && is_strictly_positive(b,ctx)))
                    return _LP_INFEASIBLE;
            } else if (nz==1) { // singleton row detected
                assert(j0>=0);
                lp_variable &var=variables[j0];
                gen val=b/a[j0];
                if (!is_positive(a[j0],ctx))
                    r*=-1;
                if (r<=0) var.tighten_ubound(val,ctx);
                if (r>=0) var.tighten_lbound(val,ctx);
            } else continue;
            constr.remove(i);
            changed=true;
        }
        if (has_infeasible_var())
            return _LP_INFEASIBLE;
        for (int i=constr.nrows();i-->0;) { // detect forcing constraints
            const vecteur &a=*constr.lhs[i]._VECTptr;
            const gen &b=constr.rhs[i];
            int r=constr.rv[i];
            gen g(0),h(0);
            for (int j=0;j<constr.ncols();++j) {
                if (is_zero(a[j],ctx))
                    continue;
                const lp_variable &var=variables[j];
                g+=a[j]*(is_positive(a[j],ctx)?var.lb():var.ub());
                h+=a[j]*(is_positive(a[j],ctx)?var.ub():var.lb());
            }
            if (is_undef(g) || is_undef(h))
                continue;
            if ((r>=0 && is_strictly_greater(b,h,ctx)) || (r<=0 && is_strictly_greater(g,b,ctx)))
                return _LP_INFEASIBLE;
            if (r!=_LP_EQ || is_inf(g) || is_inf(h))
                continue;
            if (is_zero(g-b,ctx) || is_zero(h-b,ctx)) { // forcing constraint detected
                for (int j=0;j<constr.ncols();++j) {
                    if (is_zero(a[j],ctx))
                        continue;
                    lp_variable &var=variables[j];
                    if (is_positive(a[j]*gen(is_zero(g-b,ctx)?1:-1),ctx))
                        var.tighten_ubound(var.lb(),ctx);
                    else var.tighten_lbound(var.ub(),ctx);
                }
                constr.remove(i);
                changed=true;
                continue;
            }
            // detect free singleton columns, tighten variable bounds
            for (int j=constr.ncols();j-->0;) {
                if (is_zero(a[j],ctx))
                    continue;
                lp_variable &var=variables[j];
                gen l=(is_positive(a[j],ctx)?b-h:b-g)/a[j]+var.ub();
                gen u=(is_positive(a[j],ctx)?b-g:b-h)/a[j]+var.lb();
                int k=0;
                for (;k<constr.nrows() && (k==i || is_zero(constr.lhs[k][j],ctx));++k);
                if (k==constr.nrows() && is_greater(l,var.lb(),ctx) && is_greater(var.ub(),u,ctx)) {
                    // free column singleton detected
                    bool intg=true;
                    if (!is_zero(b,ctx)) var.push_subs_coef(-1,b/a[j]);
                    intg=intg && (b/a[j]).is_integer();
                    for (k=0;k<constr.ncols();++k) {
                        if (k!=j && !is_zero(a[k],ctx)) {
                            var.push_subs_coef(k,-a[k]/a[j]);
                            intg=intg && variables[k].is_integral() && (a[k]/a[j]).is_integer();
                        }
                    }
                    if (var.is_integral() && !intg) {
                        var.clear_subs_coef();
                        continue;
                    }
                    remove_variable(j);
                    constr.remove(i);
                    changed=true;
                    break;
                } else { // tighten bounds on var
                    if (var.tighten_lbound(l,ctx) || var.tighten_ubound(u,ctx))
                        changed=true;
                }
            }
        }
        if (has_infeasible_var())
            return _LP_INFEASIBLE;
        // remove fixed variables and empty columns
        for (int j=constr.ncols();j-->0;) {
            lp_variable &var=variables[j];
            gen val;
            if (is_zero__VECT(jth_column(constr.lhs,j),ctx)) { // empty column
                int res=var.find_opt_free(objective.first[j],val,ctx);
                if (res!=0)
                    return res;
                if (!is_zero(val,ctx)) var.push_subs_coef(-1,val);
                remove_variable(j);
                changed=true;
            } else if (var.is_fixed()) { // fixed variable
                val=var.lb();
                if (var.is_integral() && !val.is_integer())
                    return _LP_INFEASIBLE;
                for (int i=0;i<constr.nrows();++i) {
                    constr.rhs[i]-=constr.lhs[i][j]*val;
                }
                if (!is_zero(val,ctx)) var.push_subs_coef(-1,val);
                remove_variable(j);
                changed=true;
            }
        }
    } while (changed);
#if 0
    if (find_imp_int && constr.nrows()>0 && constr.ncols()>0)
        find_implied_integers();
#endif
    return 0;
}

/*
 * Postprocess, affecting the solution.
 */
void lp_problem::postprocess() {
    while (!removed_cols.empty()) {
        int j=removed_cols.top();
        lp_variable &var=removed_vars.top();
        solution.insert(solution.begin()+j,0);
        gen &s=solution[j];
        vector<intgen>::const_iterator it=var.subs_coef().begin(),itend=var.subs_coef().end();
        for (;it!=itend;++it)
            s+=(it->first<0?1:solution[it->first])*it->second;
        variables.insert(variables.begin()+j,var);
        removed_cols.pop();
        removed_vars.pop();
    }
}

/*
 * Print the constraint matrix dimensions and sparsity.
 */
void lp_problem::print_constraint_matrix_dim(char *buffer) {
    sprintf(buffer,gettext("Constraint matrix has %d rows, %d columns, and %d nonzeros"),constr.nrows(),constr.ncols()+1,constr.nonzeros());
}

/*
 * Solve the problem using the specified settings.
 */
int lp_problem::solve(bool make_exact) {
    stats=lp_stats();
    char buffer[256];
    print_constraint_matrix_dim(buffer);
    message(buffer);
    if (make_exact)
        make_problem_exact();
    if (objective.first.size()!=variables.size())
        objective.first.resize(variables.size(),0);
    imp_int_count=0;
    if (settings.presolve>0) {
        message(gettext("Preprocessing..."));
        int res=preprocess();
        if (res!=0)
            return res; // the problem is either infeasible or unbounded
    }
    int rr=constr.remove_linearly_dependent(ctx);
    if (rr>0) {
        sprintf(buffer,gettext("Removed %d redundant equality constraint(s)"),rr);
        message(buffer);
    }
    if (settings.presolve>0 || rr>0) {
        print_constraint_matrix_dim(buffer);
        message(buffer);
    }
    // print information about variables and constraints
    int cvc=0,ivc=0,bvc=0,ecc=0,inecc=0;
    for (vector<lp_variable>::const_iterator it=variables.begin();it!=variables.end();++it) {
        if (it->is_integral())
            ++ivc;
        else ++cvc;
        if (it->is_binary())
            ++bvc;
    }
    sprintf(buffer,gettext("Variables: %d continuous, %d integer (%d binary)"),cvc,ivc,bvc);
    message(buffer);
    if (imp_int_count>0) {
        sprintf(buffer,gettext("Found %d implied integer variables"),imp_int_count);
        message(buffer);
    }
    for (ints::const_iterator it=constr.rv.begin();it!=constr.rv.end();++it) {
        if (*it==0) ecc++;
        else inecc++;
    }
    sprintf(buffer,gettext("Constraints: %d equalities, %d inequalities"),ecc,inecc);
    message(buffer);
    if (settings.maximize) { // convert to minimization problem
        objective.first=multvecteur(-1,objective.first);
        objective.second=-objective.second;
    }
    int nvars=variables.size();
    iteration_count=0;
    if (nvars>0) {
        if (constr.nrows()==0) { // the constraint matrix is empty
            solution.resize(nvars);
            gen val;
            for (int i=0;i<nvars;++i) {
                int res=variables[i].find_opt_free(objective.first[i],val,ctx);
                if (res!=0)
                    return res;
                solution[i]=val;
                objective.second+=val*objective.first[i];
            }
            optimum=objective.second;
        } else { // use simplex method
            make_all_vars_bounded_below();
            add_slack_variables();
            objective_norm=_l2norm(objective.first,ctx).to_double(ctx);
            int result;
            optimum=undef;
            double opt_approx;
            lp_node root(this,true);
            root.resize_ranges(nv());
            root.set_depth(0);
            message(gettext("Optimizing..."));
            if ((result=root.solve_relaxation())!=_LP_SOLVED)
                return result;
            if (root.is_integer_feasible()) {
                solution=root.get_solution();
                optimum=root.get_optimum();
            } else {
                message(gettext("Starting branch & bound..."));
                if (settings.nodeselect<0)
                    settings.nodeselect=_LP_BEST_LOCAL_BOUND;
                double root_optimum=root.get_optimum().to_double(ctx);
                double root_infeas=root.get_infeas().to_double(ctx);
                root.resize_ranges(nv());
                vector<lp_node> active_nodes(1,root);
                clock_t t=clock(),t0=t,now;
                int n,j,k,nsel=settings.nodeselect,vsel=settings.varselect;
                double opt_lbound,fr,max_score;
                bool depth_exceeded=false,incumbent_updated,is_use_pseudocost=false;
                map<double,int> candidates;
                set<pair<pair<double,double>,int> > cand;
                ints best_cand;
                pair<double,double> qual;
                while (!active_nodes.empty()) {
                    if (settings.node_limit>0 && stats.subproblems_examined>=settings.node_limit) {
                        message(gettext("node limit exceeded"),2);
                        break;
                    }
                    if (settings.iteration_limit>0 && iteration_count>settings.iteration_limit) {
                        message(gettext("simplex iteration limit exceeded"),2);
                        break;
                    }
                    n=active_nodes.size();
                    if (n>stats.max_active_nodes)
                        stats.max_active_nodes=n;
                    opt_lbound=DBL_MAX;
                    bool dive=nsel==_LP_DEPTHFIRST || (nsel==_LP_HYBRID && is_undef(optimum));
                    k=nsel==_LP_BREADTHFIRST?0:(dive?int(active_nodes.size())-1:-1);
                    double proj,iopt=is_undef(optimum)?0:optimum.to_double(ctx);
                    for (vector<lp_node>::const_iterator it=active_nodes.begin();it!=active_nodes.end();++it) {
                        double bnd=it->get_opt_approx();
                        if (opt_lbound>bnd)
                            opt_lbound=bnd;
                        if (nsel==_LP_BEST_LOCAL_BOUND || (nsel==_LP_HYBRID && !is_undef(optimum)))
                            qual=make_pair(bnd,-it->get_depth());
                        else if (nsel==_LP_BEST_PROJECTION) {
                            proj=bnd+(iopt-root_optimum)*it->get_infeas().to_double(ctx)/root_infeas;
                            qual=make_pair(proj,bnd);
                        } else continue;
                        cand.insert(make_pair(qual,it-active_nodes.begin()));
                    }
                    if (!cand.empty()) {
                        best_cand.clear();
                        const pair<double,double> &best_qual=cand.begin()->first;
                        for (set<pair<pair<double,double>,int> >::const_iterator it=cand.begin();it!=cand.end();++it) {
                            if (it->first==best_qual)
                                best_cand.push_back(it->second);
                            else break;
                        }
                        cand.clear();
#if 0
                        for (int i=best_cand.size();i-->0;) {
                            if (best_cand.size()==1)
                                break;
                            if (!active_nodes[best_cand[i]].is_up_branch())
                                best_cand.erase(best_cand.begin()+i);
                        }
#endif
                    }
                    if (!best_cand.empty() && k<0) {
                        int idx;
                        do idx=_rand(best_cand.size(),ctx).val; while (idx<0 || idx>=int(best_cand.size()));
                        k=best_cand[idx];
                    }
                    if (k<0) {
                        message(gettext("node selection strategy failed"),1);
                        break;
                    }
                    j=-1;
                    if (vsel==_LP_PSEUDOCOST || vsel<0) {
                        max_score=0;
                        for (int i=0;i<nv();++i) {
                            if (!active_nodes[k].is_var_fractional(i))
                                continue;
                            lp_variable &var=variables[i];
                            double score=var.score(active_nodes[k].get_fractional_var(i));
                            if (score==0) {
                                j=-1;
                                break;
                            }
                            if (score>max_score) {
                                j=i;
                                max_score=score;
                            }
                        }
                        if (j>=0 && !is_use_pseudocost) {
                            report_status("Switched to pseudocost-based branching");
                            is_use_pseudocost=true;
                        }
                    }
                    if (j<0) switch (vsel) {
                    case -1:
                    case _LP_MOSTFRACTIONAL:
                    case _LP_PSEUDOCOST:
                        j=active_nodes[k].get_most_fractional();
                        break;
                    case _LP_FIRSTFRACTIONAL:
                        j=active_nodes[k].get_first_fractional_var();
                        break;
                    case _LP_LASTFRACTIONAL:
                        j=active_nodes[k].get_last_fractional_var();
                        break;
                    default: assert(false);
                    }
                    if (j<0) {
                        message(gettext("branching variable selection strategy failed"),1);
                        break;
                    }
                    if (!is_undef(optimum)) {
                        stats.mip_gap=is_zero(optimum,ctx)?-opt_lbound:(opt_approx-opt_lbound)/std::abs(opt_approx);
                        if (stats.mip_gap<=settings.relative_gap_tolerance) {
                            if (settings.relative_gap_tolerance>0)
                                message(gettext("integrality gap threshold reached"),2);
                            break;
                        }
                    }
                    incumbent_updated=false;
                    double p=active_nodes[k].get_fractional_var(j);
                    for (int i=0;i<2;++i) {
                        lp_node child_node(this);
                        active_nodes[k].init_child(child_node);
                        if (settings.depth_limit>0 && child_node.get_depth()>settings.depth_limit) {
                            if (!depth_exceeded) {
                                message (gettext("depth limit exceeded"),2);
                                depth_exceeded=true;
                            }
                            break;
                        }
                        if (i==(dive?0:1)) {
                            child_node.set_up_branch(true);
                            child_node.get_range(j).tighten_lbound(_ceil(active_nodes[k].get_solution()[j],ctx),ctx);
                            if (is_strictly_positive(child_node.get_range(j).lb(),ctx)) {
                                switch (variables[j].sign_type()) {
                                case _LP_VARSIGN_POS_PART:
                                    child_node.get_range(j-1).set_ub(0);
                                    break;
                                case _LP_VARSIGN_NEG_PART:
                                    child_node.get_range(j+1).set_ub(0);
                                }
                            }
                        } else {
                            child_node.get_range(j).tighten_ubound(_floor(active_nodes[k].get_solution()[j],ctx),ctx);
                            child_node.set_up_branch(false);
                        }
                        ++stats.subproblems_examined;
                        if (child_node.solve_relaxation()==_LP_SOLVED) {
#if 1
                            pair<int,double> viol=constr.violated_constraints(child_node.get_solution(),ctx);
                            if (viol.first!=0) {
                                sprintf(buffer,gettext("relaxation solution is infeasible with violation %g, please report"),viol.second);
                                message(buffer,1);
                                return _LP_ERROR;
                            }
#endif
                            variables[j].update_pseudocost(std::abs(child_node.get_opt_approx()-active_nodes[k].get_opt_approx()),p,dive?1-i:1);
                            if (is_undef(optimum) || is_strictly_greater(optimum,child_node.get_optimum(),ctx)) {
                                if (child_node.is_integer_feasible()) { // incumbent found
                                    sprintf(buffer,"Incumbent solution found: %g",(settings.maximize?-1:1)*child_node.get_opt_approx());
                                    if (!is_undef(optimum))
                                        sprintf(buffer+strlen(buffer)," (improvement: %g%%)",(opt_approx-child_node.get_opt_approx())/std::abs(opt_approx)*100.0);
                                    report_status(buffer);
                                    solution=child_node.get_solution();
                                    optimum=child_node.get_optimum();
                                    opt_approx=child_node.get_opt_approx();
                                    incumbent_updated=true;
                                } else {
                                    active_nodes.push_back(child_node);
                                    if (settings.use_heuristic) {
                                        vecteur heur_sol;
                                        gen heur_cost;
                                        if (child_node.rounding_heuristic(heur_sol,heur_cost)) {
                                            if (is_undef(optimum) || is_strictly_greater(optimum,heur_cost,ctx)) {
                                                solution=heur_sol;
                                                optimum=heur_cost;
                                                opt_approx=heur_cost.to_double(ctx);
                                                incumbent_updated=true;
                                                sprintf(buffer,"Solution found by heuristic: %g",(settings.maximize?-1:1)*opt_approx);
                                                report_status(buffer);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (depth_exceeded)
                        break;
                    // pruning
                    active_nodes.erase(active_nodes.begin()+k);
                    if (incumbent_updated) {
                        for (int i=active_nodes.size();i-->0;) {
                            if (is_greater(active_nodes[i].get_optimum(),optimum,ctx))
                                active_nodes.erase(active_nodes.begin()+i);
                        }
                    }
                    now=clock();
                    if (settings.time_limit>0 && 1e3*double(now-t0)/CLOCKS_PER_SEC>settings.time_limit) {
                        message(gettext("time limit exceeded"),2);
                        break;
                    }
                    if (CLOCKS_PER_SEC/double(now-t)<=settings.status_report_freq) { //report status
                        sprintf(buffer,"%d nodes active, bound: %g",
                                     (int)active_nodes.size(),opt_lbound*(settings.maximize?-1:1));
                        if (stats.mip_gap>=0)
                            sprintf(buffer+strlen(buffer),", gap: %g%%",stats.mip_gap*100);
                        report_status(buffer);
                        t=clock();
                    }
                }
                if (active_nodes.empty())
                    report_status("Tree is empty");
                if (!is_undef(optimum)) {
                    //show branch&bound summary
                    sprintf(buffer,gettext("Summary:\n * %d subproblem(s) examined\n * max. tree size: %d nodes\n * %d GMI cut(s) applied"),
                            stats.subproblems_examined,stats.max_active_nodes,stats.cuts_applied);
                    if (stats.cuts_applied>0)
                        sprintf(buffer+strlen(buffer),gettext(" (average improvement: %g%%)"),stats.cut_improvement/stats.cuts_applied);
                    message(buffer);
                }
            }
            if (is_undef(optimum))
                return _LP_INFEASIBLE;
            for (int i=nv();i-->0;) {
                lp_variable &var=variables[i];
                switch (var.sign_type()) {
                case _LP_VARSIGN_NEG:
                    solution[i]=-solution[i];
                    break;
                case _LP_VARSIGN_NEG_PART:
                    solution[i+1]-=solution[i];
                    solution.erase(solution.begin()+i);
                    break;
                }
            }
            solution.resize(nvars);
        }
    } else {
        if (settings.presolve>0) {
            message(gettext("Solution found by preprocessor"));
            optimum=objective.second;
        } else {
            optimum=undef;
            solution.clear();
        }
    }
    if (settings.presolve>0)
        postprocess();
    if (settings.maximize)
        optimum=-optimum;
    return _LP_SOLVED;
}

#ifdef HAVE_LIBGLPK

/*
 * Create GLPK problem from constr.
 */
glp_prob *lp_problem::glpk_initialize() {
    glp_prob *glp=glp_create_prob();
    glp_add_rows(glp,nc());
    glp_add_cols(glp,nv());
    glp_set_obj_dir(glp,settings.maximize?GLP_MAX:GLP_MIN);
    const vecteur &obj=objective.first;
    int obj_sz=obj.size();
    for (int i=0;i<=nv();++i) {
        glp_set_obj_coef(glp,i,(i==0?objective.second:(i<=obj_sz?obj[i-1]:gen(0))).to_double(ctx));
        if (i>0) {
            lp_variable &var=variables[i-1];
            glp_set_col_kind(glp,i,var.is_integral()?GLP_IV:GLP_CV);
            int bound_type=GLP_FR;
            if (!var.range().is_unrestricted_below() && var.range().is_unrestricted_above())
                bound_type=GLP_LO;
            else if (var.range().is_unrestricted_below() && !var.range().is_unrestricted_above())
                bound_type=GLP_UP;
            else if (!var.range().is_unrestricted_below() && !var.range().is_unrestricted_above())
                bound_type=is_zero(var.ub()-var.lb(),ctx)?GLP_FX:GLP_DB;
            double lo=var.range().is_unrestricted_below()?0:var.lb().to_double(ctx);
            double hi=var.range().is_unrestricted_above()?0:var.ub().to_double(ctx);
            if (var.is_integral()) {
                lo=ceil(lo);
                hi=floor(hi);
                if (bound_type==GLP_DB && lo==hi)
                    bound_type=GLP_FX;
            }
            if (bound_type==GLP_DB && lo==0 && hi==1 && var.is_integral())
                glp_set_col_kind(glp,i,GLP_BV);
            else glp_set_col_bnds(glp,i,bound_type,lo,hi);
        }
    }
    int n=constr.nrows()*constr.ncols();
    int *ia=new int[n+1],*ja=new int[n+1]; int k=0;
    double *ar=new double[n+1];
    gen a;
    for (int i=0;i<constr.nrows();++i) {
        for (int j=0;j<constr.ncols();++j) {
            if (is_zero(a=constr.lhs[i][j],ctx))
                continue;
            ++k;
            ia[k]=i+1;
            ja[k]=j+1;
            ar[k]=a.to_double(ctx);
        }
        double rh=constr.rhs[i].to_double(ctx);
        switch (constr.rv[i]) {
        case _LP_EQ:
            glp_set_row_bnds(glp,i+1,GLP_FX,rh,0.0);
            break;
        case _LP_LEQ:
            glp_set_row_bnds(glp,i+1,GLP_UP,0.0,rh);
            break;
        case _LP_GEQ:
            glp_set_row_bnds(glp,i+1,GLP_LO,rh,0.0);
            break;
        }
    }
    glp_load_matrix(glp,k,ia,ja,ar);
    delete[] ia;
    delete[] ja;
    delete[] ar;
    return glp;
}

/*
 * Solve LP problem using GLPK implementation of simplex method.
 */
int lp_problem::glpk_simplex(glp_prob *prob) {
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev=settings.verbose?GLP_MSG_ON:GLP_MSG_ERR;
    if (settings.iteration_limit>0)
        parm.it_lim=settings.iteration_limit;
    if (settings.time_limit>0)
        parm.tm_lim=settings.time_limit;
    parm.presolve=settings.presolve>0?GLP_ON:GLP_OFF;
    return glp_simplex(prob,&parm);
}

/*
 * Solve LP problem using GLPK implementation of interior point method.
 */
int lp_problem::glpk_interior_point(glp_prob *prob) {
    glp_iptcp parm;
    glp_init_iptcp(&parm);
    parm.msg_lev=settings.verbose?GLP_MSG_ON:GLP_MSG_ERR;
    return glp_interior(prob,&parm);
}

void lp_problem::glpk_callback(glp_tree *tree,void *info) {
    if (interrupted || ctrl_c) {
        interrupted=ctrl_c=false;
        glp_ios_terminate(tree);
    }
    lp_problem *prob=static_cast<lp_problem*>(info);
    switch (glp_ios_reason(tree)) {
    case GLP_ISELECT:
        if (prob->use_blb)
            glp_ios_select_node(tree,glp_ios_best_node(tree));
        break;
    case GLP_IBINGO:
        if (prob->settings.nodeselect==_LP_HYBRID && !prob->use_blb) {
            prob->message(gettext("Switching to best-local-bound node selection"));
            prob->use_blb=true;
        }
        break;
    default:
        break;
    }
}

/*
 * Solve MIP problem using GLPK implementation of branch and cut method.
 */
int lp_problem::glpk_branchcut(glp_prob *prob) {
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.msg_lev=settings.verbose?GLP_MSG_ON:GLP_MSG_ERR;
    if (settings.time_limit>0)
        parm.tm_lim=settings.time_limit;
    parm.out_frq=(int)(1000.0/settings.status_report_freq);
    parm.mip_gap=settings.relative_gap_tolerance;
    parm.gmi_cuts=settings.max_cuts>0?GLP_ON:GLP_OFF;
    parm.mir_cuts=settings.max_cuts>0?GLP_ON:GLP_OFF;
    parm.clq_cuts=settings.has_binary_vars?GLP_ON:GLP_OFF;
    parm.cov_cuts=settings.has_binary_vars?GLP_ON:GLP_OFF;
    parm.presolve=settings.presolve>0?GLP_ON:GLP_OFF;
    parm.pp_tech=settings.presolve==1?GLP_PP_ALL:(settings.presolve==2?GLP_PP_ROOT:GLP_PP_NONE);
    parm.sr_heur=settings.use_heuristic?GLP_ON:GLP_OFF;
    parm.cb_func=&glpk_callback;
    parm.cb_info=static_cast<void*>(this);
    switch (settings.varselect) {
    case _LP_FIRSTFRACTIONAL:
        parm.br_tech=GLP_BR_FFV;
        break;
    case _LP_LASTFRACTIONAL:
        parm.br_tech=GLP_BR_LFV;
        break;
    case _LP_MOSTFRACTIONAL:
        parm.br_tech=GLP_BR_MFV;
        break;
    case _LP_PSEUDOCOST:
        parm.br_tech=GLP_BR_PCH;
        break;
    default:
        parm.br_tech=GLP_BR_DTH;
    }
    switch (settings.nodeselect) {
    case _LP_HYBRID:
    case _LP_DEPTHFIRST:
        parm.bt_tech=GLP_BT_DFS;
        break;
    case _LP_BREADTHFIRST:
        parm.bt_tech=GLP_BT_BFS;
        break;
    case _LP_BEST_PROJECTION:
        parm.bt_tech=GLP_BT_BPH;
        break;
    case _LP_BEST_LOCAL_BOUND:
    default:
        parm.bt_tech=GLP_BT_BLB;
    }
    use_blb=false;
    return glp_intopt(prob,&parm);
}

int lp_problem::term_hook(void *info,const char *s) {
    *static_cast<std::ostream*>(info) << s;
    return true;
}

#endif

/*
 * Solve the problem using the GLPK library.
 */
int lp_problem::glpk_solve() {
#ifndef HAVE_LIBGLPK
    message(gettext("GLPK library not found, using the native solver"),2);
    return solve(true);
#else
    /* redirect GLPK output to logptr */
    int term_old;
    if (settings.verbose)
        glp_term_hook(term_hook,static_cast<void*>(logptr(ctx)));
    else term_old=glp_term_out(GLP_OFF);
    glp_prob *prob=glpk_initialize();
    int result=0,solution_status;
    bool is_mip=has_integral_variables();
    switch (settings.solver) {
    case _LP_SIMPLEX:
        if (is_mip) {
            result=glpk_branchcut(prob);
            if (result==GLP_EMIPGAP || result==GLP_ETMLIM || result==GLP_ESTOP)
                result=0;
            if (result==0)
                solution_status=glp_mip_status(prob);
        } else {
            result=glpk_simplex(prob);
            if (result==GLP_ETMLIM || result==GLP_EITLIM)
                result=0;
            if (result==0)
                solution_status=glp_get_status(prob);
        }
        break;
    case _LP_INTERIOR_POINT:
        if ((result=glpk_interior_point(prob))==0)
            solution_status=glp_ipt_status(prob);
        break;
    }
    if (result==0) { //solving process finished successfully
        switch (solution_status) {
        case GLP_OPT:
            solution_status=_LP_SOLVED;
            break;
        case GLP_FEAS:
            message(gettext("the solution is not necessarily optimal"),2);
            solution_status=_LP_SOLVED;
            break;
        case GLP_INFEAS:
            message(gettext("the solution is infeasible"),2);
            solution_status=_LP_SOLVED;
            break;
        case GLP_NOFEAS:
            solution_status=_LP_INFEASIBLE;
            break;
        case GLP_UNBND:
            solution_status=_LP_UNBOUNDED;
            break;
        case GLP_UNDEF:
            solution_status=_LP_ERROR;
            break;
        }
        if (solution_status==_LP_SOLVED) { //get solution and optimum
            solution.resize(nv());
            switch (settings.solver) {
            case _LP_SIMPLEX:
                if (is_mip) {
                    optimum=glp_mip_obj_val(prob);
                    for (int i=0;i<nv();++i) {
                        solution[i]=glp_mip_col_val(prob,i+1);
                        if (variables[i].is_integral())
                            solution[i]=_round(solution[i],ctx);
                    }
                }
                else {
                    optimum=glp_get_obj_val(prob);
                    for (int i=0;i<nv();++i) {
                        solution[i]=glp_get_col_prim(prob,i+1);
                    }
                }
                break;
            case _LP_INTERIOR_POINT:
                optimum=glp_ipt_obj_val(prob);
                for (int i=0;i<nv();++i) {
                    solution[i]=glp_ipt_col_prim(prob,i+1);
                }
                break;
            }
        }
    } else switch (result) { // GLPK returned an error
    case GLP_EBADB:   message(gettext("invalid basis"),1); break;
    case GLP_ESING:   message(gettext("singular matrix"),1); break;
    case GLP_ECOND:   message(gettext("ill-conditioned matrix"),1); break;
    case GLP_EBOUND:  message(gettext("invalid bounds"),1); break;
    case GLP_EFAIL:   message(gettext("solver failed"),1); break;
    case GLP_EOBJLL:  message(gettext("objective lower limit reached"),1); break;
    case GLP_EOBJUL:  message(gettext("objective upper limit reached"),1); break;
    case GLP_ENOPFS:  message(gettext("no primal feasible solution"),1); break;
    case GLP_ENODFS:  message(gettext("no dual feasible solution"),1); break;
    case GLP_EROOT:   message(gettext("root LP optimum not provided"),1); break;
    case GLP_ENOFEAS: message(gettext("no primal/dual feasible solution"),1); break;
    case GLP_ENOCVG:  message(gettext("no convergence"),1); break;
    case GLP_EINSTAB: message(gettext("numerical instability"),1); break;
    case GLP_EDATA:   message(gettext("invalid data"),1); break;
    case GLP_ERANGE:  message(gettext("result out of range"),1); break;
    default:          message(gettext("unknown GLPK error"),1); break;
    }
    glp_delete_prob(prob);
    if (settings.verbose)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
    if (result==0)
        return solution_status;
    return _LP_ERROR;
#endif
}

/*
 * Load problem from file in LP or (possibly gzipped) MPS format. File types
 * are distinguished by examining the filename extension. If it is .gz or .mps,
 * it is assumed that the file is in (gzipped) MPS format. If it is .lp, it is
 * assumed that file is in LP format.
 */
bool lp_problem::glpk_load_from_file(const char *fname) {
#ifndef HAVE_LIBGLPK
    message(gettext("loading CPLEX and MPS files requires GLPK library"),1);
    return false;
#endif
#ifdef HAVE_LIBGLPK
    glp_prob *prob=glp_create_prob();
    bool print_messages=true;
    int result,term_old;
    string str(fname);
    string ext(str.substr(str.find_last_of(".") + 1));
    if (print_messages)
        glp_term_hook(term_hook,static_cast<void*>(logptr(ctx)));
    else term_old=glp_term_out(GLP_OFF);
    if (ext=="mps" || ext=="gz") // MPS format (modern)
        result=glp_read_mps(prob,GLP_MPS_FILE,NULL,fname);
    else if (ext=="lp") // CPLEX LP format
        result=glp_read_lp(prob,NULL,fname);
    else {
        message(gettext("file format not supported"),1);
        result=1;
    }
    if (result==0) { // successfully loaded
        char buffer[256];
        int obj_dir=glp_get_obj_dir(prob),len,t,k;
        int nr=glp_get_num_rows(prob),n=glp_get_num_cols(prob);
        if (nr*n>LP_CONSTR_MAXSIZE) {
            message(gettext("constraint matrix too large"),1);
            return false;
        }
        int *ind=new int[1+n];
        double *val=new double[1+n];
        create_variables(n);
        constr.lhs=*_matrix(makesequence(nr,nv()),ctx)._VECTptr;
        constr.rhs=vecteur(nr);
        constr.rv=ints(nr);
        variable_identifiers=vecteur(nv());
        settings.maximize=(obj_dir==GLP_MAX);
        objective.second=glp_get_obj_coef(prob,0);
        objective.first=vecteur(nv());
        for (int j=0;j<nv();++j) {
            objective.first[j]=glp_get_obj_coef(prob,j+1);
            lp_variable &var=variables[j];
            t=glp_get_col_type(prob,j+1);
            k=glp_get_col_kind(prob,j+1);
            var.set_name(string(glp_get_col_name(prob,j+1)));
            variable_identifiers[j]=identificateur(var.name());
            if (t!=GLP_UP && t!=GLP_FR)
                var.set_lb(glp_get_col_lb(prob,j+1));
            if (t!=GLP_LO && t!=GLP_FR)
                var.set_ub(glp_get_col_ub(prob,j+1));
            if (k!=GLP_CV)
                var.set_integral(true);
            if (k==GLP_BV) {
                var.tighten_lbound(0,ctx);
                var.tighten_ubound(1,ctx);
            }
        }
        for (int i=nr;i>0;--i) {
            if ((len=glp_get_mat_row(prob,i,ind,val))>0) {
                vecteur &row=*constr.lhs[i-1]._VECTptr;
                for (int j=1;j<=len;++j) {
                    row[ind[j]-1]=val[j];
                }
            }
            t=glp_get_row_type(prob,i);
            if (t==GLP_FR) {
                sprintf(buffer,gettext("row bounds not set, discarding constraint %d"),i);
                message(buffer,2);
                constr.remove(i-1);
                continue;
            }
            if (t==GLP_FX) {
                constr.rhs[i-1]=glp_get_row_lb(prob,i);
                constr.rv[i-1]=_LP_EQ;
            }
            else {
                if (t==GLP_LO || t==GLP_DB) {
                    constr.rhs[i-1]=glp_get_row_lb(prob,i);
                    constr.rv[i-1]=_LP_GEQ;
                }
                if (t==GLP_DB)
                    constr.append(*constr.lhs[i-1]._VECTptr,glp_get_row_ub(prob,i),_LP_LEQ);
                if (t==GLP_UP || t==GLP_DB) {
                    constr.rhs[i-1]=glp_get_row_ub(prob,i);
                    constr.rv[i-1]=_LP_LEQ;
                }
            }
        }
        delete[] ind;
        delete[] val;
    }
    glp_delete_prob(prob);
    if (print_messages)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
    return result==0;
#endif
}

bool lp_problem::assign_variable_types(const gen &g,int t) {
    pair<gen,gen> range;
    int i,i0=array_start(ctx);
    if (g.type==_VECT) {
        const_iterateur it=g._VECTptr->begin(),itend=g._VECTptr->end();
        for (;it!=itend;++it) {
            if (!assign_variable_types(*it,t))
                return false;
        }
    } else if ((g.type==_IDNT && (i=get_variable_index(*g._IDNTptr))>=0) || (g.is_integer() && (i=g.val-i0)>=0))
        variables[i].set_type(t,ctx);
    else if (interval2pair(g,range,ctx) && range.first.is_integer() && range.second.is_integer()) {
        for (i=range.first.val;i<=range.second.val;++i) {
            variables[i-i0].set_type(t,ctx);
        }
    } else return false;
    return true;
}

bool parse_limit(const gen &g,int &lim,GIAC_CONTEXT) {
    if (is_positive(g,contextptr)) {
        if (g.is_integer())
            lim=g.val;
        else if (is_inf(g))
            lim=0;
        else return false;
    } else return false;
    return true;
}

bool parse_options_and_bounds(const_iterateur &it,const_iterateur &itend,lp_problem &prob) {
    for (;it!=itend;++it) {
        if (*it==at_maximize)
            prob.settings.maximize=true;
        else if (*it==at_minimize)
            prob.settings.maximize=false;
        else if (is_mcint(*it)) {
            switch(it->val) {
            case _LP_MAXIMIZE:
                prob.settings.maximize=true;
                break;
            case _LP_VERBOSE:
                prob.settings.verbose=true;
                break;
            case _LP_PRESOLVE:
                prob.settings.presolve=true;
                break;
            }
        } else if (it->is_symb_of_sommet(at_equal)) {
            //parse the argument in form "option=value"
            pair<gen,gen> bounds;
            gen a,lh(_lhs(*it,prob.ctx)),rh(_rhs(*it,prob.ctx));
            if (lh.type==_IDNT && interval2pair(rh,bounds,prob.ctx)) {
                int vi=prob.get_variable_index(*lh._IDNTptr);
                if (vi>=0)
                    prob.tighten_variable_bounds(vi,bounds.first,bounds.second);
                else return false;
            } else if (is_idnt_vecteur(lh,prob.ctx) && interval2pair(rh,bounds,prob.ctx)) {
                const_iterateur it=lh._VECTptr->begin(),itend=lh._VECTptr->end();
                for (;it!=itend;++it) {
                    int vi=prob.get_variable_index(*(it->_IDNTptr));
                    if (vi>=0)
                        prob.tighten_variable_bounds(vi,bounds.first,bounds.second);
                    else return false;
                }
            } else if (lh==at_maximize && rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                prob.settings.maximize=(bool)rh.val;
            else if (lh==at_assume && rh.is_integer())
                prob.settings.assumption=rh.val;
            else if (is_mcint(lh)) {
                switch(lh.val) {
                case _LP_INTEGERVARIABLES:
                case _LP_BINARYVARIABLES:
                    if (!prob.assign_variable_types(rh,lh.val))
                        return false;
                    break;
                case _LP_DEPTHLIMIT:
                    if (!parse_limit(rh,prob.settings.depth_limit,prob.ctx))
                        return false;
                    break;
                case _LP_NODE_LIMIT:
                    if (!parse_limit(rh,prob.settings.node_limit,prob.ctx))
                        return false;
                    break;
                case _LP_ITERATION_LIMIT:
                    if (!parse_limit(rh,prob.settings.iteration_limit,prob.ctx))
                        return false;
                    break;
                case _LP_MAX_CUTS:
                    if (!parse_limit(rh,prob.settings.max_cuts,prob.ctx))
                        return false;
                    break;
                case _LP_TIME_LIMIT:
                    if (!parse_limit(rh,prob.settings.time_limit,prob.ctx))
                        return false;
                    break;
                case _LP_NODESELECT:
                    if (is_mcint(rh))
                        prob.settings.nodeselect=rh.val;
                    else return false;
                    if (prob.settings.nodeselect<_LP_DEPTHFIRST || prob.settings.nodeselect>_LP_HYBRID)
                        return false;
                    break;
                case _LP_VARSELECT:
                    if (is_mcint(rh))
                        prob.settings.varselect=rh.val;
                    else return false;
                    if (prob.settings.varselect<_LP_FIRSTFRACTIONAL || prob.settings.varselect>_LP_PSEUDOCOST)
                        return false;
                    break;
                case _LP_METHOD:
                    if (rh==at_exact)
                        prob.settings.precision=_LP_EXACT;
                    else if (rh==at_float)
                        prob.settings.precision=_LP_INEXACT;
                    else if (is_mcint(rh)) {
                        if (rh.val==_LP_INTERIOR_POINT) {
                            prob.settings.precision=_LP_INEXACT;
                            prob.settings.solver=_LP_INTERIOR_POINT;
                        } else if (rh.val==_LP_SIMPLEX) {
                            prob.settings.precision=_LP_PROB_DEPENDENT;
                            prob.settings.solver=_LP_SIMPLEX;
                        } else return false;
                    } else return false;
                    break;
                case _LP_GAP_TOLERANCE:
                    a=_evalf(rh,prob.ctx);
                    if (a.type==_DOUBLE_ && is_positive(a,prob.ctx))
                        prob.settings.relative_gap_tolerance=a.DOUBLE_val();
                    else return false;
                    break;
                case _LP_MAXIMIZE:
                    if (rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                        prob.settings.maximize=(bool)rh.val;
                    else return false;
                    break;
                case _LP_VERBOSE:
                    if (rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                        prob.settings.verbose=(bool)rh.val;
                    else return false;
                    break;
                case _LP_ASSUME:
                    if (rh.is_integer())
                        prob.settings.assumption=rh.val;
                    else return false;
                    break;
                case _LP_PRESOLVE:
                    if (rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                        prob.settings.presolve=rh.val;
                    else if (rh==at_maple_root)
                        prob.settings.presolve=2;
                    else return false;
                    break;
                case _LP_HEURISTIC:
                    if (rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                        prob.settings.use_heuristic=(bool)rh.val;
                    else return false;
                    break;
                case _GT_ACYCLIC:
                    if (rh.is_integer() && rh.subtype==_INT_BOOLEAN)
                        prob.settings.acyclic=(bool)rh.val;
                    else return false;
                    break;
                default:
                    break;
                }
            } else return false;
        } else return false;
    }
    return true;
}

int parse_constraints(bool is_matrix_form,const_iterateur &it,const_iterateur &itend,lp_problem &prob) {
    if (is_matrix_form) {
        if (it->type!=_VECT)
            return _LP_ERR_TYPE;
        const vecteur arg=*it->_VECTptr;
        int n=arg.size(),mask=num_mask_withfrac|num_mask_withint;
        if ((n%2)!=0)
            return _LP_ERR_SIZE;
        if (n>0) {
            //constraints are given in form [A,b,Aeq,beq] such that
            //A.x<=b and/or Aeq.x=beq (Aeq and beq are optional)
            if (!((arg[0].type==_VECT && arg[0]._VECTptr->empty()) || ckmatrix(arg[0])) || arg[1].type!=_VECT)
                return _LP_ERR_TYPE;
            const matrice &A=*arg[0]._VECTptr;
            const vecteur &b=*arg[1]._VECTptr;
            if (!is_numericm(A,mask) || !is_numericv(b,mask))
                return _LP_ERR_TYPE;
            int len=A.size();
            if (len!=int(b.size()) || (!A.empty() && mcols(A)!=prob.nv()))
                return _LP_ERR_DIM;
            for (int i=0;i<len;++i) {
                prob.constr.append(*A[i]._VECTptr,b[i],_LP_LEQ);
            }
            if (n>2) { //there are equality constraints
                if (!ckmatrix(arg[2]) || arg[3].type!=_VECT)
                    return _LP_ERR_TYPE;
                const matrice &Aeq=*arg[2]._VECTptr;
                const vecteur &beq=*arg[3]._VECTptr;
                if (!is_numericm(Aeq,mask) || !is_numericv(beq,mask))
                    return _LP_ERR_TYPE;
                if ((len=Aeq.size())!=int(beq.size()) || mcols(Aeq)!=prob.nv())
                    return _LP_ERR_DIM;
                for (int i=0;i<len;++i) {
                    prob.constr.append(*Aeq[i]._VECTptr,beq[i],_LP_EQ);
                }
            }
            if (!ckmatrix(prob.constr.lhs))
                return _LP_ERR_DIM;
            if (prob.nv()==0)
                prob.create_variables(prob.constr.ncols());
            else if (prob.nv()!=prob.constr.ncols())
                return _LP_ERR_DIM;
        }
        //if the third arg is [bl,bu] then assume bl<=x<=bu
        ++it;
        if (it!=itend && it->type==_VECT) {
            matrice &bounds=*it->_VECTptr;
            if (!ckmatrix(bounds) || bounds.size()!=2 || int(bounds.front()._VECTptr->size())<prob.nv())
                return _LP_ERR_DIM;
            if (prob.nv()==0)
                prob.create_variables(bounds.front()._VECTptr->size());
            for (int i=0;i<prob.nv();++i) {
                prob.tighten_variable_bounds(i,bounds[0][i],bounds[1][i]);
            }
            ++it;
        }
        if (prob.nv()==0)
            return _LP_ERR_SIZE;
        if (prob.objective.first.empty())
            prob.set_objective(vecteur(prob.nv(),0),0);
    } else { //the problem is given in symbolic form
        if (it->type==_VECT) {
            const vecteur &arg=*it->_VECTptr;
            prob.add_identifiers_from(arg);
            ++it;
            vecteur sides,c;
            pair<gen,gen> bounds;
            gen r;
            int rel;
            const_iterateur jt=arg.begin(),itend=arg.end();
            for (;jt!=itend;++jt) {
                if (jt->is_symb_of_sommet(at_equal))
                    rel=_LP_EQ;
                else if (jt->is_symb_of_sommet(at_inferieur_egal))
                    rel=_LP_LEQ;
                else if (jt->is_symb_of_sommet(at_superieur_egal))
                    rel=_LP_GEQ;
                else
                    return _LP_ERR_TYPE;
                sides=*jt->_SYMBptr->feuille._VECTptr;
                if (jt->is_symb_of_sommet(at_equal) && interval2pair(sides.back(),bounds,prob.ctx)) {
                    if (!prob.lincomb_coeff(sides.front(),c,r))
                        return _LP_ERR_TYPE;
                    prob.constr.append(c,bounds.first-r,_LP_GEQ);
                    prob.constr.append(c,bounds.second-r,_LP_LEQ);
                    continue;
                }
                if (!prob.lincomb_coeff(sides.front()-sides.back(),c,r))
                    return _LP_ERR_TYPE;
                prob.constr.append(c,-r,rel);
            }
        }
        prob.create_variables(prob.variable_identifiers.size());
    }
    if (!parse_options_and_bounds(it,itend,prob))
        return _LP_ERR_TYPE;
    return 0;
}

gen _lpsolve(const gen &args,GIAC_CONTEXT) {
    gen g(args);
    if (g.type==_STRNG)
        g=makesequence(g);
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gensizeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    vecteur obj;
    const_iterateur it=gv.begin(),itend=gv.end();
    lp_problem prob(contextptr); //create LP problem with default settings
    bool is_matrix_form=false,sort_vars=true;
    if (it->type==_STRNG) { //problem is given in file
        int len=_size(*it,contextptr).val;
        string fname(it->_STRNGptr->begin(),it->_STRNGptr->begin()+len);
        if (!prob.glpk_load_from_file(fname.c_str()))
            return undef;
        ++it;
        if (it->type==_VECT)
            return gentypeerr(contextptr);
        sort_vars=false;
    } else { //problem is entered manually
        is_matrix_form=it->type==_VECT;
        if (is_matrix_form) {
            obj=*it->_VECTptr;
            if (!obj.empty()) {
                prob.set_objective(obj,0);
                prob.create_variables(obj.size());
            }
        } else {
            prob.add_identifiers_from(*it);
            gen obj_ct;
            if (!prob.lincomb_coeff(*it,obj,obj_ct))
                return gentypeerr(contextptr);
            prob.set_objective(obj,obj_ct);
        }
        ++it;
    }
    switch (parse_constraints(is_matrix_form,it,itend,prob)) {
    case _LP_ERR_SIZE:
        return gensizeerr(contextptr);
    case _LP_ERR_TYPE:
        return gentypeerr(contextptr);
    case _LP_ERR_DIM:
        return gendimerr(contextptr);
    default:
        if (prob.nc()==0) {
            for (int i=0;i<prob.nv();++i) {
                lp_variable &var=prob.variables[i];
                vecteur row(prob.nv(),0);
                if (!is_inf(var.lb())) {
                    row[i]=1;
                    prob.constr.append(row,var.lb(),_LP_GEQ);
                    break;
                }
                else if (!is_inf(var.ub())) {
                    row[i]=1;
                    prob.constr.append(row,var.ub(),_LP_LEQ);
                    break;
                }
            }
            if (prob.nc()==0) {
                prob.message(gettext("no constraints detected"),1);
                return gensizeerr(contextptr);
            }
        }
        break;
    }
    for (int i=0;i<prob.nv();++i) {
        lp_variable &var=prob.variables[i];
        switch (prob.settings.assumption) {
        case _ZINT:
        case _LP_INTEGER:
            var.set_type(_LP_INTEGERVARIABLES,contextptr);
            break;
        case _LP_BINARY:
            var.set_type(_LP_BINARYVARIABLES,contextptr);
            break;
        case _NONNEGINT:
        case _LP_NONNEGINT:
            var.set_type(_LP_INTEGERVARIABLES,contextptr);
        case _LP_NONNEGATIVE:
            var.tighten_lbound(0,contextptr);
            break;
        }
    }
    bool is_solver_exact;
    switch (prob.settings.precision) {
    case _LP_EXACT:
        is_solver_exact=true;
        break;
    case _LP_INEXACT:
        is_solver_exact=false;
        break;
    case _LP_PROB_DEPENDENT:
        is_solver_exact=!prob.has_approx_coefficients();
        break;
    }
    vector<lp_variable>::const_iterator vt=prob.variables.begin();
    for (;vt!=prob.variables.end();++vt) {
        if (vt->is_integral() && is_zero(vt->lb(),contextptr) && is_one(vt->ub())) {
            prob.settings.has_binary_vars=true;
            break;
        }
    }
    //solve the problem
    switch (is_solver_exact?prob.solve(true):prob.glpk_solve()) {
    case _LP_INFEASIBLE:
        prob.message(gettext("problem has no feasible solutions"),1);
        return vecteur(0);
    case _LP_UNBOUNDED:
        prob.message(gettext("problem has unbounded solution"),1);
        return makevecteur(prob.settings.maximize?plus_inf:minus_inf,vecteur(0));
    case _LP_ERROR:
        return undef;
    default: //_LP_SOLVED
        return gen(makevecteur(_eval(prob.optimum,contextptr),_eval(prob.output_solution(sort_vars),contextptr)),_LIST__VECT);
    }
}
static const char _lpsolve_s[]="lpsolve";
static define_unary_function_eval (__lpsolve,&_lpsolve,_lpsolve_s);
define_unary_function_ptr5(at_lpsolve,alias_at_lpsolve,&__lpsolve,0,true)

#ifndef NO_NAMESPACE_GIAC
}
#endif //ndef NO_NAMESPACE_GIAC
