/*
 * optimization.h
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

#ifndef __OPTIMIZATION_H
#define __OPTIMIZATION_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "first.h"
#include "gen.h"
#include "unary.h"
#include <set>
#include <stack>
#include <sstream>
#ifdef HAVE_LIBGSL
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_rng.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

#define GOLDEN_RATIO 1.61803398875   // used in golden-ratio search method

enum critical_point_classification {
    _CPCLASS_UNDECIDED      = 0,
    _CPCLASS_MIN            = 1,
    _CPCLASS_MAX            = 2,
    _CPCLASS_POSSIBLE_MIN   = 3,
    _CPCLASS_POSSIBLE_MAX   = 4,
    _CPCLASS_SADDLE         = 5
};

enum kernel_density_estimation_method {
    _KDE_METHOD_EXACT,
    _KDE_METHOD_PIECEWISE,
    _KDE_METHOD_LIST
};

enum bandwidth_selection_method {
    _KDE_BW_METHOD_DPI,
    _KDE_BW_METHOD_ROT
};

enum bvp_output_type {
    _BVP_LIST,
    _BVP_DIFF,
    _BVP_PIECEWISE,
    _BVP_SPLINE
};

enum hclust_linkage_type {
    _SINGLE_LINKAGE   = 1,
    _AVERAGE_LINKAGE  = 2,
    _COMPLETE_LINKAGE = 3,
    _WEIGHTED_LINKAGE = 4,
    _WARD_LINKAGE     = 5,
    _MEDIAN_LINKAGE   = 6,
    _CENTROID_LINKAGE = 7
};
enum hclust_index_function {
    _SILHOUETTE_INDEX     = 1,
    _MCCLAIN_RAO_INDEX    = 2,
    _DUNN_INDEX           = 4,
    _ALL_HCLUST_INDICES   = 7
};

enum clustering_indices {
    _HARTIGAN_CRITERION         = 0,
    _BANFELD_RAFTERY_INDEX      = 1,      // min
    _DAVIES_BOULDIN_INDEX       = 2,      // min
    _RAY_TURI_INDEX             = 4,      // min
    _SCOTT_SYMONS_INDEX         = 8,      // min
    _CALINSKI_HARABASZ_INDEX    = 16,     // max
    _PBM_INDEX                  = 32,     // max
    _RATKOWSKY_LANCE_INDEX      = 64,     // max
    _DET_RATIO_INDEX            = 128,    // min diff
    _LOG_DET_RATIO_INDEX        = 256,    // min diff
    _LOG_SS_RATIO_INDEX         = 512,    // min diff
    _KSQ_DETW_INDEX             = 1024,   // max diff
    _BALL_HALL_INDEX            = 2048,   // max diff
    _TRACE_W_INDEX              = 4096,   // max diff
    _TRACE_WIB_INDEX            = 8192,   // max diff
    _ALL_CLUSTER_INDICES        = 16383
};

enum nlp_method {
    _NLP_AUTOMATIC              = 1,
    _NLP_INTERIOR_POINT         = 2,
    _NLP_NELDER_MEAD            = 3,
    _NLP_DIFFERENTIAL_EVOLUTION = 4,
    _NLP_COBYLA                 = 5
};

enum nlp_solution_status {
    _NLP_OPTIMAL,
    _NLP_FAILED,
    _NLP_INFEAS,
    _NLP_UNBOUNDED,
    _NLP_ERROR,
    _NLP_PENDING
};

class ipdiff {
    /* IPDIFF CLASS (Implicit Partial DIFFerentiation)
     * This class is used for implicit differentiation of f with respect to g=0.
     * Here, f:R^(n+m)->R, g:R^(n+m)->R^m with rank(g')=m and J(x)=f(x,h(x)) where x in R^n and h is
     * the implicit function guaranteed to exist (in a neighborhood of some point x0) by the
     * Implicit Function Theorem if rank(g')=m. Function J may now be differentiated wrt x1,x2,...,xn. */
public:
    typedef std::vector<int> ivector;
    typedef ivector::const_iterator ivector_iter;
    typedef std::vector<ivector> ivectors;
    typedef std::map<ivector,int> ivector_map;
    typedef std::pair<ivector,ivector_map> diffterm;
    typedef std::map<diffterm,int> diffterms;
    typedef std::map<ivector,gen> pd_map;
private:
    vecteur vars; // list of variables [x1,x2,..,xn,y1,y2,..,ym] where x are free and y are dependent
    gen f; // the expression f(x1,x2,..,xn,y1,y2,..,ym)
    vecteur g; // list [g1,g2,..,gm], it is assumed that gj(x1,x2,..,xn,y1,y2,..,ym)=0 for all j=1,2,..,m
    const context *ctx; // giac context
    pd_map pdv; // partial derivatives of J
    pd_map pdf; // partial derivatives of f
    pd_map pdg; // partial derivatives of g
    pd_map pdh; // partial derivatives of h
    std::map<ivector,diffterms> cterms;
    int ord; // the maximal order of the currently computed partial derivatives of h
    int nvars; // the number of free variables
    int nconstr; // the number of constraints = the number of dependent variables
    diffterms derive_diffterms(const diffterms &terms,ivector &sig);
    void find_nearest_terms(const ivector &sig,diffterms &match,ivector &excess);
    const gen &differentiate(const gen &e,pd_map &pds,const ivector &sig);
    void compute_h(const std::vector<diffterms> &grv,int order);
    void compute_pd(int order,const ivector &sig=ivector(0));
    void raise_order(int order); // compute the missing partial derivatives of h
    const gen &get_pd(const pd_map &pds,const ivector &sig) const;
public:
    /* construct the ipdiff instance for function f_orig(vars_orig) subject to
     * the constraints in g_orig, the Jacobian of g_orig must have maximal rank */
    ipdiff(const gen &f_orig,const vecteur &g_orig,const vecteur &vars_orig,GIAC_CONTEXT);
    /* return the partial derivative of J represented by sig */
    const gen &derivative(const ivector &sig);
    /* return the partial derivative of J with respect to the variables in dvars */
    const gen &derivative(const vecteur &dvars);
    /* store all partial derivatives of the specified order to pd_map */
    void partial_derivatives(int order,pd_map &pdmap);
    /* return 1/k!*(sum(i=1)^n (xi-ai)*d/dxi)^k f(x)|a */
    gen taylor_term(const vecteur &a,int k,bool shift=true);
    /* return the Taylor polynomial, with the specified order, of J in the vicinity of a */
    gen taylor(const vecteur &a,int order);
    /* store grad(J) as vector res */
    void gradient(vecteur &res);
    /* store Hessian(J) as matrix res */
    void hessian(matrice &res);
    /* generate all partitions of m into n terms and store them to vector c */
    static void ipartition(int m,int n,ivectors &c,const ivector &p=ivector(0));
    /* return the sum the elements of a vector v of integers, ignore the last element if drop_last=true */
    static int sum_ivector(const ivector &v,bool drop_last=false);
};

class tprob {
    /* TPROB CLASS (Transportation PROBlem)
     * The class implementing the MODI method for balanced TP with degeneracy handling */
public:
    typedef std::pair<int,int> ipair;
    typedef std::vector<ipair> ipairs;
private:
    const context *ctx;
    vecteur demand;
    vecteur supply;
    gen eps; // epsilon
    gen M; // symbol for marking forbidden routes in cost matrix
    void north_west_corner(matrice &feas);
    ipairs stepping_stone_path(ipairs &path_orig,const matrice &X);
    void modi(const matrice &P_orig,matrice &X);
public:
    /* construct the TP with supply s and demand d, m marks forbidden routes */
    tprob(const vecteur &s,const vecteur &d,const gen &m,GIAC_CONTEXT);
    /* solve the transportation problem with the given cost matrix, output in sol */
    void solve(const matrice &cost_matrix,matrice &sol);
};

bool is_mcint(const gen &g,int v=-1);
vecteur sort_identifiers(const vecteur &v,GIAC_CONTEXT);
std::vector<int> linearly_dependent_rows(const matrice &m,GIAC_CONTEXT);
bool get_assumptions(const gen &g,int &dom,matrice &intervals,vecteur &excluded,GIAC_CONTEXT);
void set_assumptions(const gen &g,const vecteur &cond,const vecteur &excluded,bool additionally,GIAC_CONTEXT);
bool is_cubic_wrt(const gen &g,const gen &x,gen &a,gen &b,gen &c,gen &d,GIAC_CONTEXT);
bool is_quartic_wrt(const gen &g,const gen &x,gen &a,gen &b,gen &c,gen &d,gen &e,GIAC_CONTEXT);
bool is_definitely_positive(const gen &g,GIAC_CONTEXT);
bool is_positive_definite(const vecteur &mat,int isnum,GIAC_CONTEXT);
bool is_positive_definite_inplace(vecteur &m,bool isnum,GIAC_CONTEXT);
bool is_positive_semidefinite(const matrice &mat,double time_limit,unsigned max_taille,GIAC_CONTEXT);
matrice madd(const matrice &a,const matrice &b);
void madd_inplace(matrice &a,const matrice &b,bool subt=false);
matrice msub(const matrice &a,const matrice &b);
void msub_inplace(matrice &a,const matrice &b);
matrice mscale(const gen &a,const matrice &b);
void mscale_inplace(const gen &a,matrice &b);
void copy_matrice(const matrice &src,matrice &dest);
bool is_string_list(const gen &g);
vecteur data2list(const gen &g,GIAC_CONTEXT);
gen find_minimum(const gen &f,const gen &a,const gen &b,bool global,double eps,int maxiter,GIAC_CONTEXT);
void sortperm(const vecteur &v,std::vector<int> &p,GIAC_CONTEXT);
int argmin(const vecteur &v,GIAC_CONTEXT);
int argmax(const vecteur &v,GIAC_CONTEXT);
int kmeans(const matrice &data,int &k,int min_k,int &max_k,int maxiter,std::vector<int> &res,vecteur *cc,matrice *ind,GIAC_CONTEXT);
void remove_elements_with_indices(vecteur &v,const set<int> &ind);
void remove_elements_with_indices(vecteur &v,const set<pair<int,gen> > &ind);
gen to_algebraic(const gen &g,vecteur &S,GIAC_CONTEXT);
bool is_algebraic(const gen &g,GIAC_CONTEXT);
extern bool is_positive_safe(const gen &g,bool strict,unsigned max_taille,GIAC_CONTEXT);

class nlp_problem {
    gen obj,best_obj_val,obj_shift;
    int _fevalc; // objective function evaluation counter
    matrice bd_rect;
    vecteur vars,obj_vars,bnd_vars;
    // problem data
    vecteur lb,ub,pc_lo,pc_lo_count,pc_hi,pc_hi_count;
    vecteur eq_cons,ineq_cons;
    vecteur obj_gradient,obj_hessian;
    matrice eq_jacobian,ineq_jacobian;
    // internal data
    matrice saved_points;
    matrice nm_simplex; // simplex vertices in Nelder-Mead
    matrice nm_yarr;    // objective function vaues at simplex vertices
    double  nm_alpha,nm_beta,nm_gamma,nm_delta;
    vecteur initp_sigma;
    vecteur cclust_centers,cclust_sizes,cclust_freqs;
    vecteur iteration_data;
    vecteur penalty_eq,penalty_ineq; // adaptive penalty factors
    vecteur initp;
    const_iterateur initp_iter; // initial-point iterator
    const context *ctx;
    int _msg_level; // 0 - quiet, 1 - errors, 2 - errors/warnings, 3 - all output
    bool _smooth_obj,_convex_obj,_smooth_constr,_convex_constr,_have_hessian;
    bool _linconstr,_linobj;
    bool _presolved,_initialized;
    bool _has_initial_rect;
    bool _compact;
    bool _infeas;
    bool _iter_data;
    bool _debug;
    int _penalty_scheme; // 0 - no penalty, 1 - Lagrangian, 2 - multiplicative (see Lampinen & Zelinka, 2000)
    int _iter_count;
    int _base_vertex_index;
    int _spc;
    vector<bool> _is_intvar;
    set<pair<int,gen> > _fixed_vars;
    set<int> _inactive_eq_indices,_inactive_ineq_indices;
    stringstream _ss;
    vector<int> _sol_stat;
    class fxvars {
        void assign(const fxvars &other) { n_oldvars=other.n_oldvars; indices=other.indices; names=other.names; values=other.values; }
    public:
        int n_oldvars;
        std::vector<int> indices;
        vecteur names;
        vecteur values;
        fxvars(int nov) { n_oldvars=nov; }
        fxvars(const fxvars &other) { assign(other); }
        fxvars &operator=(const fxvars &other) { assign(other); return *this; }
    };
    std::stack<fxvars> fx_subs;
    class lineq {
        void assign(const lineq &other) { lvars=other.lvars; nlvars=other.nlvars; index=other.index; e=other.e; }
    public:
        std::vector<int> lvars;
        std::vector<int> nlvars;
        int index;
        gen e;
        lineq(const gen &g,int i) { e=g; index=i; }
        lineq(const lineq &other) { assign(other); }
        lineq &operator=(const lineq &other) { assign(other); return *this; }
    };
    static std::pair<int,int> make_leqv_sys(const std::vector<int> &ind,const std::vector<lineq> &leqv);
    bool eliminate_equalities(const std::vector<lineq> &leqv,vecteur &lsol,vecteur &leq_vars);
    bool subs_fxvars(const fxvars &fv);
public:
    typedef struct optimum {
        int res;
        vecteur x;
        gen f;
        optimum();
    } optimum_t;
    typedef std::vector<optimum_t> optima_t;
    struct meth_parm {
        double eps;             // optimality tolerance
        double tol;             // feasibility tolerance
        gen penalty;            // penalty factors/exponents for constraints (PENALTY[i] corresponds to the i-th constraint)
        double scale;           // scaling factor for differential evolution (DE)
        double cross_prob;      // cross-probability for mutation in DE
        double reflect_ratio;   // ratio used for point reflection in Nelder-Mead (NM)
        double expand_ratio;    // ratio used for polytope expansion in NM
        double contract_ratio;  // ratio used for polytope contraction in NM
        double shrink_ratio;    // ratio used for polytope shrinking in NM
        double size;            // initial simplex size for NM
        double step;            // step-size for NM and BFGS
        int search_points;      // number of starting points, used in multistart and DE
        int max_iter;           // maximum number of iterations
        int msg_level;          // filtering messages
        bool presolve;          // whether to presolve the problem first
        bool postprocess;       // whether to refine global solution using local methods
        int convex;             // is this a convex problem? (saves the time required for checking)
        int cluster;            // number of clusters of initial points (set to 1 for automatic clustering)
        bool border;            // setting this to false disables automatic bound finding (default: true)
        void init_defaults(double eps);
    };
    class ipt_solver {
        /* this class implements interior-point method for solving: min  f(x) s.t. c(x)=0, x>=0
         * SOURCE: Waechter & Biegler, 2004 */
        // fixed internal parameters
        double kappa_eps;       // positive
        double kappa_mu;        // in (0,1)
        double kappa_Sigma;     // >= 1
        double kappa_soc;       // in (0,1)
        double kappa_1;         // positive
        double tau_min;         // positive
        double theta_mu;        // in (1,2)
        gen    theta_min;       // positive
        gen    theta_max;       // positive
        double delta;           // positive
        double s_max;           // positive
        double s_theta;         // > 1
        double s_phi;           // > 1
        double eta_phi;         // in (0,1/2)
        double gamma_alpha;     // in (0,1]
        double gamma_theta;     // in (0,1)
        double gamma_phi;       // in (0,1)
        int p_max;              // positive
        double bardelta_w_min;  // positive
        double bardelta_w_0;    // > bardelta_w_min
        double bardelta_w_max;  // > bardelta_w_0
        double bardelta_c;      // positive
        double kappa_w_minus;   // in (0,1)
        double kappa_w_plus;    // > 1
        double barkappa_w_plus; // > kappa_w_plus
        double kappa_c;         // positive
        double kappa_resto;     // positive
        double lambda_max;      // positive
        // other parameters
        double eps_tol,feas_tol;
        int maxiter;
        double delta_w,delta_w_last,delta_c;
        double mu_j,tau_j;
#ifdef HAVE_LIBLAPACK
        double *factorization,*lapack_work,*lapack_rhs;
        int *lapack_ipiv;
#else
        matrice factorization;
        vecteur ldl_perm;
#endif
        gen f,mu,theta_k,theta_resto,phi_k;
        std::set<int> dup_ind;
        const context *ctx;
        nlp_problem &prob;
        int nvars,ncons,nslacks;
        vecteur c,grad_f,hess_f;
        matrice A;  // transposed Jacobian of constraints
        vecteur x,x_subs,lambda,z;
        matrice kkt_mat;
        vecteur kkt_rhs,sigma;
        matrice F; // filter
        bool initialized;
        char buf[256];
        vecteur iteration_data;
        vecteur pprod(const vecteur &vec1,const vecteur &vec2) { return *_pointprod(makesequence(vec1,vec2),ctx)._VECTptr; }
        gen linfn(const vecteur &v) { return _max(_abs(v,ctx),ctx); }
        vecteur subs_vars(const vecteur &vec,const vecteur &x0,const vecteur &lambda0,const vecteur &z0,double muj=-1);
        void compute_dz(const vecteur &x0,const vecteur &dx0,const vecteur &z0,vecteur &dz);
        gen E(const gen &mu0,const vecteur &x0,const vecteur &lambda0,const vecteur &z0);
        bool IC(const matrice &mat);
        void init_filter();
        void update_filter(const gen &b_theta,const gen &b_phi);
        vecteur init_lambda(const vecteur &x0,const vecteur &z0);
        void compute_np(const gen &mu0,const vecteur &x0,vecteur &n0,vecteur &p0) const;
        gen constraint_violation(const vecteur &x0) const;
        gen barrier_fval(const vecteur &x0) const;
        bool filter_accepts(const gen &theta,const gen &phi) const;
        bool kkt_init();
        int feas_restoration(vecteur &x0,vecteur &lambda0,vecteur &z0,const gen &viol_old,gen &viol);
        gen fraction_to_boundary(const vecteur &y,const vecteur &d_y);
        int ls_filter_barrier_method(vecteur &x_k,vecteur &lambda_k,vecteur &z_k,double mu0);
        void init_parameters();
        bool solve_kkt(const matrice &mat,const vecteur &rh,vecteur &sol);
    public:
        ipt_solver(nlp_problem &p,const meth_parm &parm);
        ~ipt_solver();
        /* Return TRUE if x0 is accepted by the currently active filter */
        bool filter_accepts(const vecteur &x0) const { return filter_accepts(constraint_violation(x0),barrier_fval(x0)); }
        /* run interior point method for each point in x0 (or in x0) and return results as list of pairs (retval,x) */
        bool optimize(optima_t &res);
        const vecteur &get_iteration_data() const { return iteration_data; }
        vecteur x2vars(const vecteur &x0) const;
        vecteur vars2x(const vecteur &x0,const vecteur &x0_orig=vecteur(0)) const;
    };
private:
    class region {
        nlp_problem &prob;
        vecteur x_lb,x_ub,c,dc,cx,sl,su;
        std::vector<int> dir,kl,ku;
        double ftol;
        int n,m,maxiter;
        const context *ctx;
        const ipt_solver *ipt;
        char buf[256];
        void rs_heur(vecteur &x0) const;
        bool mfeas(vecteur &x0,bool store);
        gen bisection(const vecteur &x0,int i,const gen &a,const gen &b,int d);
        bool strides(const vecteur &x0_orig,int i,int d,int lim,gen &a,gen &b);
        bool find_bound(const vecteur &x0,int i,int d,gen &bnd);
    public:
        /* tol: feasibility tolerance, itlim: iteration limit for strides/bisection */
        region(nlp_problem &p,double tol,int itlim);
        /* attach/detach ipt solver for premature termination on filter acceptance */
        void attach_ipt(const ipt_solver *ipts) { ipt=ipts; }
        void detach_ipt() { ipt=NULL; }
        /* make the point x0 feasible, returns TRUE if successful */
        bool make_feasible(vecteur &x0) { return mfeas(x0,true); }
        /* generate a feasible point, returns TRUE if successful */
        bool generate_feasible_point(vecteur &x0) { rs_heur(x0); return make_feasible(x0); }
        /* approximate the bounding hyperrectangle for this region, may return infinity bounds */
        matrice bounding_rect();
    };
    region *feasible_region;
    clock_t _start_time;
    void err(const std::string &s,bool translate=true) const
        { if (_msg_level>0) *logptr(ctx) << gettext("Error") << ": " << (translate?gettext(s.c_str()):s.c_str()) << "\n"; }
    void err() { err(_ss.str(),false); _ss.str(std::string()); }
    void warn(const std::string &s,bool translate=true) const
        { if (_msg_level>1) *logptr(ctx) << gettext("Warning") << ": " << (translate?gettext(s.c_str()):s.c_str()) << "\n"; }
    void warn() { warn(_ss.str(),false); _ss.str(std::string()); }
    void msg(const std::string &s,bool translate=true) const
        { if (_msg_level>2) *logptr(ctx) << (translate?gettext(s.c_str()):s.c_str()) << "\n"; }
    void msg() { msg(_ss.str(),false); _ss.str(std::string()); }
    void debug(const std::string &s) const
        { if (!_debug) return; if (_iter_count>0) *logptr(ctx) << _iter_count << ": "; *logptr(ctx) << s << "\n"; }
    void debug() { debug(_ss.str()); _ss.str(std::string()); }
    gen make_label(const vecteur &pt,const std::string &str,int disp=0) {
        return _legende(makesequence(_point(pt,ctx),string2gen(str,false),symb_equal(at_display,change_subtype(disp,_INT_COLOR))),ctx);
    }
    bool compute_obj_val(const vecteur &x_orig,gen &val,bool print_err=false);
    bool subst_hessian(const vecteur &x0,matrice &H,bool full=true) const;
    bool is_expression_convex(const vecteur &gradient,matrice &H,bool &hess_ok,bool lt);
    bool clamp(vecteur &x0,bool rnd=false) const;
    void make_bounded_vars(bool use_bd_rect);
    void purge_bounded_vars() { if (!bnd_vars.empty()) _purge(change_subtype(bnd_vars,_SEQ__VECT),ctx); }
    void remove_redundant_constraints();
    void initialize(int method,const meth_parm &parm);
    double initial_size(const meth_parm &parm,double lo=0.02,double hi=0.1) const;
    bool initialize_simplex(const vecteur &x0,double size);
    void initialize_penalties(const gen &p);
    void randomize_ith_component(vecteur &x,int i,bool use_varbnds) const;
    bool make_random_initial_point(vecteur &x,gen *fval=NULL,int nsamp=10,double sigma_scale=0.01);
    bool restart(vecteur &x,gen *fval=NULL);
    void save_point(const vecteur &pt) { if (initp.empty()) saved_points.push_back(pt); }
    void clear_saved_points() { saved_points.clear(); }
    bool preprocess(); // returns false if infeasibility was detected
    void postprocess_point(gen &x0,const vecteur &oldvars,const fxvars &fx,vecteur &csol);
    void postprocess(gen &x0); // post-process x0
    void reset_timer() { _start_time=clock(); }
    double elapsed_secs() { return double(clock()-_start_time)/CLOCKS_PER_SEC; }
    int var_count() const { return (int)vars.size()-(int)_fixed_vars.size(); }
    int eq_count() const { return (int)eq_cons.size()-(int)_inactive_eq_indices.size(); }
    int ineq_count() const { return (int)ineq_cons.size()-(int)_inactive_ineq_indices.size(); }
    int eval_count() const { return _fevalc; } // return the total number of objective function evaluations so far
    void find_fixed_vars();
    void insert_fixed_vars(vecteur &x) const;
    gen subs_fixed_vars(const gen &g) const;
    vecteur subs_fixed_vars(const vecteur &v) const;
    bool find_inactive_eq(const vecteur &eq);
    bool find_inactive_ineq(const vecteur &ineq);
    vecteur constraints(bool incl_bnds=true,bool all_geq=false) const;
    bool is_smooth() const { return _smooth_obj && _smooth_constr; }
    bool is_convex() const { return _convex_obj && _convex_constr; }
    bool is_unconstrained() const;
    bool has_bounded_variables() const;
    bool is_feasible(const vecteur &x,double feas_tol,bool check_bounds,set<int> *viol_eq=NULL,set<int> *viol_ineq=NULL,set<int> *viol_bnd=NULL) const;
    gen INT(const gen &g) const { return _penalty_scheme==2?_trunc(g,ctx):g; }
    bool is_intvar(int i) const { return _is_intvar[i]; }
    bool has_intvars() const { return std::find(_is_intvar.begin(),_is_intvar.end(),true)!=_is_intvar.end(); }
    int intvar_count() const;
    bool is_intsol(const meth_parm &parm,const vecteur &sol) const;
    bool is_intvars_all_binary() const;
    int minimize_linear(const gen &o,const vecteur &c,vecteur &opt,gen &optval) const;
    int nm_iteration(vecteur &best_point,gen &fval,double &ssize,const meth_parm &parm,bool no_gsl=false);
    bool de_initialize(iterateur &it,const iterateur &itend,iterateur &kt,iterateur &ft,int &best_k);
    bool cobyla(const vecteur &cons,const meth_parm &parm,optima_t &res);
    bool differential_evolution(const meth_parm &parm,optimum_t &res);
    bool nelder_mead(const meth_parm &parm,optima_t &res);
    bool branch_and_bound(const meth_parm &parm,const vecteur &cont_sol,optimum_t &incumbent);
    bool branch_var(const meth_parm &parm,const vecteur &sol,const gen &cost,int strategy,int &var);
    bool bb_solve_subproblem(const meth_parm &parm,vecteur &active_nodes,gen &fval,optima_t &sub_res,optimum_t &incumbent);
    void update_pseudocost(int var,vecteur &pc,vecteur &pc_count,const gen &df,const gen &dx);
    bool outer_approximation(const meth_parm &parm,const vecteur &cont_sol,optimum_t &incumbent);
    int var2intvar(int var) const;
    void get_vars_and_values(const set<pair<int,gen> > &lst,vecteur &v,vecteur &w) const;
    void find_best_solution(const meth_parm &parm,optima_t &res,vecteur &opt,gen &optval,int &ret);
#ifdef HAVE_LIBGSL
    bool gsl_bfgs(const meth_parm &parm,optima_t &res);
    static double gsl_my_f(const gsl_vector *x,void *params);
    static void gsl_my_df(const gsl_vector *x,void *params,gsl_vector *df);
    static void gsl_my_fdf(const gsl_vector *x,void *params,double *f,gsl_vector *df);
    static double gsl_my_uf(double x,void *params);
    gsl_multimin_fminimizer *nm_s;
#endif
public:
    nlp_problem(const gen &f,const vecteur &g,const vecteur &x,GIAC_CONTEXT);
    ~nlp_problem();
    void set_variable_bounds(const vecteur &L,const vecteur &U) { lb=L; ub=U; }
    void set_integer_variable(int i) { _is_intvar[i]=true; }
    void enable_debug_messages(bool yes) { _debug=yes; }
    void collect_iteration_data(bool yes) { _iter_data=true; }
    const vecteur &get_iteration_data() const { return iteration_data; }
    int optimize(int method,const meth_parm &parm,const vecteur &ip,vecteur &opt,gen &optval);
    /* solve lhs*X=rhs using LSQ and write the result in sol; return TRUE if successful and FALSE on error */
    static bool solve_lsq(const matrice &lhs,const vecteur &rhs,vecteur &sol,GIAC_CONTEXT);
    static bool ampl_load(const std::string &filename,vecteur &v,vecteur &o,vecteur &c,string &msg,GIAC_CONTEXT);
};

class cluster_crit {
    const context *ctx;
    const matrice &pts;
    gen distf;
    int N;                                      // number of observations
    int p;                                      // number of variables
    int K;                                      // number of clusters
    std::vector<int> _ci;                       // cluster indices for observations
    std::vector<int> _n;                        // cluster sizes
    std::map<std::pair<int,int>,gen> Delta_map; // cluster-barycenters-distance map
    matrice _mu;                                // observation and cluster barycenters
    vecteur _WGSSk;                             // sums of squared distances of observations from their barycenters
    gen _WGSS;                                  // sum of WGSSk for k=1..K
    matrice _T;     bool have_T;                // total scatter matrix
    matrice _WG;    bool have_WG;               // within-group scatter matrices (for each cluster)
    matrice _BG;    bool have_BG;               // between-group scatter matrix
    matrice _WGsum;                             // sum of matrices WGk for k=1..K
    vecteur _TSS;                               // diagonal of T
    vecteur _BGSS;                              // diagonal of BG
    gen _DB_min;                                // minimal distance between observations from different clusters
    gen _DB_max;                                // maximal distance between observations from different clusters
    gen ET;                                     // depends only on the given data; used for computing PBM index
    /* return the barycenter of k-th cluster or the
     * global barycenter if k is omitted */
    const vecteur &mu(int k=-1) { return *_mu[k+1]._VECTptr; }
    gen WGSS(int k=-1) { return k<0?_WGSS:_WGSSk[k]; }
    const vecteur &T();
    const vecteur &WG(int l=-1);
    const vecteur &BG();
    gen TSS(int k=-1);
    gen BGSS(int k=-1);
    gen DB(int dir);
    int ci(int k) const { return _ci[k]; }
    int n(int k) const { return _n[k]; }
    gen Delta(int k,int l);
    // individual criteria
    gen banfeld_raftery();
    gen calinski_harabasz();
    gen davies_bouldin();
    gen pbm();
    gen ratkowsky_lance();
    gen ray_turi();
    gen scott_symons();
    gen det_ratio();
    gen log_det_ratio();
    gen log_ss_ratio();
    gen ksq_detW();
    gen ball_hall();
    gen trace_W();
    gen trace_WiB();
    // store index value
    bool _updated;
    std::map<int,vecteur> _vmap;
    std::map<int,std::pair<int,double> > _optvals;
    void update(int pos,const gen &val,int dir);
public:
    cluster_crit(const matrice &data,const gen &dist_func,GIAC_CONTEXT);
    void init(int k,const int *indices,const double *bc,const double *wgss);
    bool compute_indices(int crit);
    vecteur get_index_values(int pos=0) { return _vmap[pos]; }
    static int optimal_number_of_clusters(const std::map<int,std::pair<int,double> > &optvals);
    int ncopt() const { return optimal_number_of_clusters(_optvals); }
    static int name2index(const std::string &name);
};

class hclust { // hierarchical agglomerative clustering
public:
    typedef std::pair<double,std::pair<int,int> > drow; // dendrogram row: (delta,a,b)
    typedef std::vector<drow> dendrogram;
private:
    int N;
    const context *ctx;
    std::map<int,vecteur> _vmap;
    vecteur _kvals;
    std::map<int,std::pair<int,double> > _optvals;
    std::vector<int> S,ord;                                     // node labels/order
    std::vector<std::pair<int,double> > slh_a,slh_s;            // silhouette data
    std::vector<std::map<int,std::pair<int,double> > > slh_b;   // silhouette data
    const matrice &points;                                      // input data reference
    gen distf;                                                  // distance function
    double *_dist_cache;                                        // cached distances
    double &dist_cache(int i,int j) { return _dist_cache[(i*(2*N-i-3))/2+j-1]; }
    double dist(int i,int j,bool cache=false);                  // compute d(i,j), use caching if CACHE=true
    bool mst_linkage(dendrogram &dg,bool cache_distances);      // single linkage algorithm
    bool nn_chain_linkage(dendrogram &dg,int meth);             // complete, average, weighted, and Ward linkage algorithm
    double formula(int a,int b,int x,int na,int nb,int nx,double dab,int meth);
    int k_low(const dendrogram &dg);
    void walk_children(const dendrogram &dg,int &pos,int row);
    void reorder_children(dendrogram &dg);
    void swap(dendrogram &d,int row);
public:
    hclust(const matrice &data,const gen &dist_func,GIAC_CONTEXT);
    ~hclust();
    /* linkage with respect to method METH (linkage type), returns
     * false if it fails to compute the distance between a pair of points */
    int linkage(dendrogram &dg,std::vector<int> &ind,int meth,int K,int min_k,int max_k);
    /* return the order of points which is optimal for drawing DG */
    const std::vector<int> &order(dendrogram &dg);
    vecteur get_index_values(int pos=0) { return mtran(makevecteur(_kvals,_vmap[pos])); }
    static int name2index(const std::string &name);
};

class PermuSort { // permutation sorting in ascending order
    const vecteur &v;
    const context *ctx;
    int as;
    class comparator {
        const PermuSort &parent;
    public:
        bool operator()(int i,int j) {
            return is_strictly_greater(parent.v[j-parent.as],parent.v[i-parent.as],parent.ctx);
        }
        comparator(const PermuSort &p) : parent(p) { }
    };
public:
    PermuSort(const vecteur &v_in,GIAC_CONTEXT) : v(v_in), ctx(contextptr) { as=array_start(contextptr); }
    void sort(std::vector<int> &p);
};

class FDWeights { // finite-difference numerical differentiation
    /* implemented according to Sadiq & Viswanath (2014) */
    const vecteur &z;
    vecteur w_lagrange;
    matrice L,R,w;
    int N,M;
    const context *ctx;
    void multbinom(const vecteur &a,vecteur &b,const gen &zeta);
    void convolve(const vecteur &a,const vecteur &b,vecteur &c);
    void lagrange_weights();
    void find_weights(const vecteur &c,int k);
public:
    FDWeights(const vecteur &grid_points,int diff_order,GIAC_CONTEXT);
    void setz0(const gen &zeta=0);
    gen operator()(int k,int m) const { return w[k][m]; } // weight for mth derivative at kth grid-point
};

class siman_config { // a configuration class for simulated annealing
    const context *ctx;
    const gen &energy_func;
    const gen &distance_func;
    const gen &modify_func;
    gen _obj;
    public:
    siman_config(const gen &obj_orig,const gen &efunc,const gen &distfunc,const gen &modfunc,GIAC_CONTEXT)
        : ctx(contextptr), energy_func(efunc), distance_func(distfunc), modify_func(modfunc) { _obj=obj_orig; }
    siman_config(const siman_config &other)
        : ctx(other.ctx), energy_func(other.energy_func), distance_func(other.distance_func), modify_func(other.modify_func) { _obj=other.obj(); }
    void assign(const siman_config &other) { _obj=other.obj(); }
    ~siman_config() {}
    const gen &obj() const { return _obj; }
    bool is_valid() const { return !is_undef(_obj); }
    int step_func_nargs() const { return modify_func._SYMBptr->feuille._VECTptr->front()._VECTptr->size(); }
    gen energy() { return evalf_double(energy_func(_obj,ctx),1,ctx); }
    gen distance(const siman_config &other) { return evalf_double(distance_func(makesequence(_obj,other.obj()),ctx),1,ctx); }
    void modify(double step_size) { _obj=modify_func(step_func_nargs()==1?_obj:makesequence(_obj,step_size),ctx); }
};

gen _implicitdiff(const gen &g,GIAC_CONTEXT);               // implicit differentiation
gen _box_constraints(const gen &g,GIAC_CONTEXT);            // converting a matrix [[L1,U1],[L2,U2],...] to x1=L1..U1,x2=L2..U2,...
gen _minimize(const gen &g,GIAC_CONTEXT);                   // exact global minimization of a differentiable function over a compact domain
gen _maximize(const gen &g,GIAC_CONTEXT);                   // exact gloabl maximization
gen _find_minimum(const gen &g,GIAC_CONTEXT);               // univariate approximation in a segment using Brent's method
gen _extrema(const gen &g,GIAC_CONTEXT);                    // exact local extrema
gen _minimax(const gen &g,GIAC_CONTEXT);                    // minimax polynomial approximation by Remez
gen _tpsolve(const gen &g,GIAC_CONTEXT);                    // transportation problem solver via MODI
gen _nlpsolve(const gen &g,GIAC_CONTEXT);                   // (mixed integer) nonlinear programming problem solver
gen _thiele(const gen &g,GIAC_CONTEXT);                     // rational interpolation using Thiele's algorithm
gen _triginterp(const gen &g,GIAC_CONTEXT);                 // trigonometric interpolation
gen _ratinterp(const gen &g,GIAC_CONTEXT);                  // rational interpolation without poles
gen _kernel_density(const gen &g,GIAC_CONTEXT);             // kernel density estimation
gen _fitdistr(const gen &g,GIAC_CONTEXT);                   // distribution fitting to sample data
gen _fitpoly(const gen &g,GIAC_CONTEXT);                    // fitting a polynomial to a function or data
gen _bvpsolve(const gen &g,GIAC_CONTEXT);                   // iterative boundary-value problem solver
gen _euler_lagrange(const gen &g,GIAC_CONTEXT);             // Euler-Lagrange equations
gen _jacobi_equation(const gen &g,GIAC_CONTEXT);            // Jacobi equation
gen _conjugate_equation(const gen &g,GIAC_CONTEXT);         // conjugate equation
gen _convex(const gen &g,GIAC_CONTEXT);                     // function convexity testing
gen _numdiff(const gen &g,GIAC_CONTEXT);                    // numerical differentiation using Fornberg's algorithm
gen _isolve(const gen &g,GIAC_CONTEXT);                     // Diophantine equation solver
gen _kmeans(const gen &g,GIAC_CONTEXT);                     // k-means clustering using Hartigan-Wong algorithm
gen _cluster(const gen &g,GIAC_CONTEXT);                    // agglomerative hierarchical clustering
gen _levenshtein(const gen &g,GIAC_CONTEXT);                // Levenshtein distance between vectors/strings
gen _isposdef(const gen &g,GIAC_CONTEXT);                   // check whether a symmetric/Hermitian matrix is positive definite
gen _symbol_array(const gen &g,GIAC_CONTEXT);               // Create vectors and matrices of symbols
gen _sortperm(const gen &g,GIAC_CONTEXT);                   // permutation that sorts the given vector in ascending order
gen _siman(const gen &g,GIAC_CONTEXT);                      // simulated annealing optimizer

extern const unary_function_ptr * const at_implicitdiff;
extern const unary_function_ptr * const at_box_constraints;
extern const unary_function_ptr * const at_minimize;
extern const unary_function_ptr * const at_maximize;
extern const unary_function_ptr * const at_find_minimum;
extern const unary_function_ptr * const at_extrema;
extern const unary_function_ptr * const at_minimax;
extern const unary_function_ptr * const at_tpsolve;
extern const unary_function_ptr * const at_nlpsolve;
extern const unary_function_ptr * const at_thiele;
extern const unary_function_ptr * const at_triginterp;
extern const unary_function_ptr * const at_ratinterp;
extern const unary_function_ptr * const at_kernel_density;
extern const unary_function_ptr * const at_fitdistr;
extern const unary_function_ptr * const at_fitpoly;
extern const unary_function_ptr * const at_bvpsolve;
extern const unary_function_ptr * const at_euler_lagrange;
extern const unary_function_ptr * const at_convex;
extern const unary_function_ptr * const at_numdiff;
extern const unary_function_ptr * const at_isolve;
extern const unary_function_ptr * const at_kmeans;
extern const unary_function_ptr * const at_cluster;
extern const unary_function_ptr * const at_levenshtein;
extern const unary_function_ptr * const at_isposdef;
extern const unary_function_ptr * const at_symbol_array;
extern const unary_function_ptr * const at_sortperm;
extern const unary_function_ptr * const at_siman;

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
#endif // __OPTIMIZATION_H
