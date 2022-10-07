/*
 * lpsolve.h
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

#ifndef __LPSOLVE_H
#define __LPSOLVE_H

#include "config.h"
#include "gen.h"
#include "unary.h"
#include <stack>
#ifdef HAVE_LIBGLPK
#include <glpk.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif //ndef NO_NAMESPACE_GIAC

#define LP_SCORE_FACTOR 0.1667
#define LP_MIN_AWAY 0.08
#define LP_MIN_PARALLELISM 0.86
#define LP_MAX_MAGNITUDE 1e6
#define LP_CONSTR_MAXSIZE 1e8
#define LP_FEAS_TOL 1e-5

typedef vector<int> ints;
typedef pair<int,gen> intgen;

enum lp_results {
    _LP_SOLVED     = 0,
    _LP_INFEASIBLE = 1,
    _LP_UNBOUNDED  = 2,
    _LP_ERROR      = 3
};

enum lp_parsing_errors {
    _LP_ERR_SIZE = 1,
    _LP_ERR_TYPE = 2,
    _LP_ERR_DIM  = 3
};

enum lp_precision_types {
    _LP_EXACT,
    _LP_INEXACT,
    _LP_PROB_DEPENDENT
};

enum lp_relation_types {
    _LP_LEQ = -1,
    _LP_EQ  =  0,
    _LP_GEQ =  1
};

enum lp_variable_sign_types {
    _LP_VARSIGN_POS,
    _LP_VARSIGN_NEG,
    _LP_VARSIGN_POS_PART,
    _LP_VARSIGN_NEG_PART
};

struct lp_settings {
    //solver parameters
    int solver;
    int precision;
    int presolve;
    bool maximize;
    int assumption;
    int iteration_limit;
    bool has_binary_vars;
    bool acyclic;
    //branch&bound parameters
    int varselect;
    int nodeselect;
    double relative_gap_tolerance;
    int depth_limit;
    int node_limit;
    int time_limit; //in miliseconds
    int max_cuts;
    bool use_heuristic;
    //message report parameters
    bool verbose;
    double status_report_freq;
    lp_settings();
};

struct lp_stats {
    int subproblems_examined;
    int cuts_applied;
    double cut_improvement;
    int max_active_nodes;
    double mip_gap;
    lp_stats();
};

class lp_range {
    gen lbound;
    gen ubound;
    void assign(const lp_range &other) { lbound=other.lb(); ubound=other.ub(); }
public:
    lp_range();
    lp_range(const lp_range &other) { assign(other); }
    ~lp_range() { }
    lp_range &operator =(const lp_range &other) { assign(other); return *this; }
    const gen &lb() const { return lbound; }
    const gen &ub() const { return ubound; }
    void set_lb(const gen &L) { lbound=L; }
    void set_ub(const gen &U) { ubound=U; }
    bool tighten_lbound(const gen &l,GIAC_CONTEXT)
        { if (is_strictly_greater(l,lbound,contextptr)) { lbound=l; return true; } return false; }
    bool tighten_ubound(const gen &u,GIAC_CONTEXT)
        { if (is_strictly_greater(ubound,u,contextptr)) { ubound=u; return true; } return false; }
    bool is_unrestricted_below() const { return is_inf(lbound); }
    bool is_unrestricted_above() const { return is_inf(ubound); }
};

class lp_variable {
    bool _is_integral;
    int _sign_type;
    lp_range _range;
    string _name;
    vector<intgen> _subs_coef;
    double pseudocost[2];
    int nbranch[2];
    void assign(const lp_variable &other);
public:
    lp_variable();
    lp_variable(const lp_variable &other) { assign(other); }
    ~lp_variable() { }
    lp_variable &operator =(const lp_variable &other) { assign(other); return *this; }
    bool is_integral() const { return _is_integral; }
    void set_integral(bool yes) { _is_integral = yes; }
    bool is_binary() const { return _is_integral && is_zero(_range.lb()) && is_one(_range.ub()); }
    int sign_type() const { return _sign_type; }
    void set_sign_type(int type) { _sign_type=type; }
    const lp_range &range() const { return _range; }
    bool is_fixed() const { return is_zero(range().ub()-range().lb()); }
    const string &name() const { return _name; }
    void set_name(const string &s) { _name=s; }
    void set_type(int t,GIAC_CONTEXT);
    void set_lb(const gen &L) { _range.set_lb(L); }
    void set_ub(const gen &U) { _range.set_ub(U); }
    const gen &lb() const { return _range.lb(); }
    const gen &ub() const { return _range.ub(); }
    bool tighten_lbound(const gen &L,GIAC_CONTEXT)
        { return _range.tighten_lbound(_is_integral?_ceil(L,contextptr):L,contextptr); }
    bool tighten_ubound(const gen &U,GIAC_CONTEXT)
        { return _range.tighten_ubound(_is_integral?_floor(U,contextptr):U,contextptr); }
    void update_pseudocost(double delta,double fr,int dir);
    double score(double fr) const;
    void push_subs_coef(int j,const gen &g) { _subs_coef.push_back(make_pair(j,g)); }
    void clear_subs_coef() { _subs_coef.clear(); }
    const vector<intgen> &subs_coef() const { return _subs_coef; }
    int find_opt_free(const gen &c,gen &val,GIAC_CONTEXT) const;
};

struct lp_constraints {
    matrice lhs,lhs_f;
    vecteur rhs,rhs_f;
    ints rv,rv_orig;
    int ncols_orig;
    int nrows() { return lhs.size(); }
    int ncols() { return lhs.empty()?0:lhs.front()._VECTptr->size(); }
    int nonzeros();
    void append(const vecteur &lh,const gen &rh,int relation_type);
    void duplicate_column(int index);
    void negate_column(int index);
    void subtract_from_rhs_column(const vecteur &v);
    void set(int index,const vecteur &lh,const gen &rh,int relation_type);
    void get(int index,vecteur &lh,gen &rh,int &relation_type);
    void get_lr(int index,vecteur &lh,gen &rh);
    void div(int index,const gen &g,GIAC_CONTEXT);
    void subtract(int index,const vecteur &v,const gen &g);
    void remove(int index);
    int remove_linearly_dependent(GIAC_CONTEXT);
    pair<int,double> violated_constraints(const vecteur &x,GIAC_CONTEXT);
};

class lp_node;

struct lp_problem {
    const context *ctx;
    pair<vecteur,gen> objective;
    double objective_norm;
    vector<double> obj_approx;
    vector<lp_variable> variables;
    vecteur variable_identifiers;
    lp_constraints constr;
    lp_settings settings;
    lp_stats stats;
    vecteur solution;
    gen optimum;
    stack<lp_variable> removed_vars;
    stack<int> removed_cols;
    int iteration_count;
    bool use_blb;
    int imp_int_count;
    lp_problem(GIAC_CONTEXT) {
        ctx=contextptr;
        settings=lp_settings();
    }
    int nc() { return constr.lhs.size(); }
    int nv() { return variables.size(); }
    void message(const char* msg,int type=0);
    void report_status(const char* msg);
    void add_identifiers_from(const gen &g);
    bool assign_variable_types(const gen &g,int t);
    int get_variable_index(const identificateur &idnt);
    void set_objective(const vecteur &v,const gen &ft);
    void create_variables(int n);
    void make_problem_exact();
    void add_slack_variables();
    void tighten_variable_bounds(int i,const gen &l,const gen &u);
    void make_all_vars_bounded_below();
    bool has_integral_variables();
    bool has_approx_coefficients();
    bool lincomb_coeff(const gen &g,vecteur &varcoeffs,gen &freecoeff);
    int preprocess(bool find_imp_int=true);
    bool has_infeasible_var() const;
    void remove_variable(int j);
    void find_implied_integers();
    void postprocess();
    void print_constraint_matrix_dim(char *buffer);
    int solve(bool make_exact);
    vecteur output_solution(bool sort_vars);
    //GLPK routines
#ifdef HAVE_LIBGLPK
    glp_prob *glpk_initialize();
    int glpk_simplex(glp_prob *prob);
    int glpk_interior_point(glp_prob *prob);
    int glpk_branchcut(glp_prob *prob);
    static int term_hook(void *info,const char *s);
    static void glpk_callback(glp_tree *tree,void *info);
#endif
    int glpk_solve();
    bool glpk_load_from_file(const char *fname);
};

class lp_node {
    lp_problem *prob;
    int depth;
    vector<lp_range> ranges;
    gen optimum;
    vecteur solution;
    double opt_approx;
    gen infeas;
    int most_fractional;
    map<int,double> fractional_vars;
    bool use_bland,up_branch,is_root;
    stack<vector<intgen> > removed_cols;
    vector<intgen> pivot_row;
    gen pivot_elm;
    int nv;
    void assign(const lp_node &other);
    gen fracpart(const gen &g) const;
    int preprocess(matrice &m,vecteur &bv,vecteur &l,vecteur &u,ints &cols,vecteur &obj,gen &obj_ct);
    bool remove_free_variable(int j,const vecteur &l,const vecteur &u,ints &cols,const vecteur &obj,gen &obj_ct);
    bool change_basis(matrice &m,const vecteur &u,vector<bool> &is_slack,ints &basis,ints &cols);
    void pivot_ij(matrice &m,int I,int J,bool negate=false);
    void simplex_reduce_bounded(matrice &m,const vecteur &u,vector<bool> &is_slack,
                                ints &basis,ints &cols,int phase,const gen &obj_ct);
public:
    lp_node(lp_problem *p,bool isroot=false) { prob=p; is_root=isroot; }
    lp_node(const lp_node &other) { assign(other); }
    ~lp_node() { }
    lp_node &operator =(const lp_node &other) { assign(other); return *this; }
    lp_problem *get_prob() const { return prob; }
    void resize_ranges(size_t s) { ranges.resize(s); }
    const vector<lp_range> &get_ranges() const { return ranges; }
    lp_range &get_range(int i) { return ranges[i]; }
    int get_depth() const { return depth; }
    void set_depth(int d) { depth=d; }
    const vecteur &get_solution() const { return solution; }
    const gen &get_optimum() const { return optimum; }
    double get_opt_approx() const { return opt_approx; }
    const gen &get_infeas() const { return infeas; }
    bool is_integer_feasible() const { return is_zero(infeas); }
    bool is_var_fractional(int index) const { return fractional_vars.find(index)!=fractional_vars.end(); }
    double get_fractional_var(int index) const;
    int get_most_fractional() const { return most_fractional; }
    int get_first_fractional_var() const { return fractional_vars.begin()->first; }
    int get_last_fractional_var() const { return fractional_vars.rbegin()->first; }
    void init_child(lp_node &child);
    void set_up_branch(bool yes) { up_branch=yes; }
    bool is_up_branch() const { return up_branch; }
    int solve_relaxation();
    bool rounding_heuristic(vecteur &sol,gen &cost) const;
};

gen _lpsolve(const gen &args,GIAC_CONTEXT);
extern const unary_function_ptr * const  at_lpsolve;

#ifndef NO_NAMESPACE_GIAC
} //namespace giac
#endif //ndef NO_NAMESPACE_GIAC
#endif //__LPSOLVE_H
