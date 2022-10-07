/*
 * graphe.h
 *
 * (c) 2018 Luka Marohnić
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

#ifndef __GRAPHE_H
#define __GRAPHE_H
#include <time.h>
#include <stdarg.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "first.h"
#include "gen.h"
#include "unary.h"
#include "moyal.h"
#include <string>
#include <iostream>
#include <fstream>
#include <queue>
#include <stack>
#include <set>
#include <bitset>
#ifdef HAVE_LIBGLPK
#include <glpk.h>
#endif
#ifdef HAVE_LIBGSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

#ifndef DBL_MAX
#define DBL_MAX 1.79769313486e+308
#endif
#define PLASTIC_NUMBER 1.32471795724
#define PLASTIC_NUMBER_2 1.75487766625
#define PLASTIC_NUMBER_3 2.32471795724
#define MARGIN_FACTOR 0.139680581996 // pow(PLASTIC_NUMBER,-7)
#define SIP_NBITS 64

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

typedef unsigned long ulong;

#define _GT_HAAR_LIMIT 60

enum gt_dot_token_type {
    _GT_DOT_TOKEN_TYPE_IDENTIFIER = 1,
    _GT_DOT_TOKEN_TYPE_NUMBER = 2,
    _GT_DOT_TOKEN_TYPE_OPERATOR = 3,
    _GT_DOT_TOKEN_TYPE_STRING = 4,
    _GT_DOT_TOKEN_TYPE_DELIMITER = 5
};

enum gt_attribute {
    _GT_ATTRIB_LABEL,
    _GT_ATTRIB_WEIGHT,
    _GT_ATTRIB_COLOR,
    _GT_ATTRIB_SHAPE,
    _GT_ATTRIB_STYLE,
    _GT_ATTRIB_WIDTH,
    _GT_ATTRIB_DIRECTED,
    _GT_ATTRIB_WEIGHTED,
    _GT_ATTRIB_POSITION,
    _GT_ATTRIB_NAME,
    _GT_ATTRIB_TEMPORARY,
    //add more here
    _GT_ATTRIB_USER  // this one must be the last
};

enum gt_layout_style {
    _GT_STYLE_DEFAULT,
    _GT_STYLE_SPRING,
    _GT_STYLE_PLANAR,
    _GT_STYLE_3D,
    _GT_STYLE_CIRCLE,
    _GT_STYLE_TREE,
    _GT_STYLE_BIPARTITE
};

enum gt_vertex_cover_algorithm {
    _GT_VC_APPROX_CLIQUE,
    _GT_VC_APPROX_ALOM,
    _GT_VC_APPROX_DFS,
    _GT_VC_EXACT
};

enum gt_conn_check {
    _GT_CC_CONNECTED, // current subgraph is connected
    _GT_CC_COMPONENTS_ARE_SUBGRAPHS, // components are subgraphs with indices 1,2,...
    _GT_CC_FIND_COMPONENTS // split graph to connected components and iterate
};

/* special graph constructors */
extern const unary_function_ptr * const at_graph;
extern const unary_function_ptr * const at_cycle_graph;
extern const unary_function_ptr * const at_path_graph;
extern const unary_function_ptr * const at_complete_graph;
extern const unary_function_ptr * const at_complete_binary_tree;
extern const unary_function_ptr * const at_complete_kary_tree;
extern const unary_function_ptr * const at_kneser_graph;
extern const unary_function_ptr * const at_odd_graph;
extern const unary_function_ptr * const at_johnson_graph;
extern const unary_function_ptr * const at_star_graph;
extern const unary_function_ptr * const at_wheel_graph;
extern const unary_function_ptr * const at_web_graph;
extern const unary_function_ptr * const at_prism_graph;
extern const unary_function_ptr * const at_antiprism_graph;
extern const unary_function_ptr * const at_grid_graph;
extern const unary_function_ptr * const at_torus_grid_graph;
extern const unary_function_ptr * const at_sierpinski_graph;
extern const unary_function_ptr * const at_petersen_graph;
extern const unary_function_ptr * const at_flower_snark;
extern const unary_function_ptr * const at_goldberg_snark;
extern const unary_function_ptr * const at_paley_graph;
extern const unary_function_ptr * const at_haar_graph;

class graphe : public gen_user {
public:
    typedef std::vector<int> ivector;
    typedef std::vector<ivector> ivectors;
    typedef std::map<int,gen> attrib;
    typedef std::pair<int,int> ipair;
    typedef std::vector<ipair> ipairs;
    typedef std::pair<double,double> dpair;
    typedef std::vector<dpair> dpairs;
    typedef std::vector<double> point;
    typedef std::vector<point> layout;
    typedef std::map<int,ipair> sparsematrow;
    typedef std::map<int,sparsematrow> sparsemat;
    typedef std::set<ipair> edgeset;
    typedef std::vector<std::map<int,int> > edgemap;
    typedef std::vector<std::map<int,double> > edgemapd;
    typedef std::map<ipair,int> intpoly;
    typedef std::vector<double> dvector;
    typedef std::set<int> iset;
    typedef std::vector<bool> bvector;
    typedef std::vector<bvector> bvectors;
    typedef std::vector<std::bitset<SIP_NBITS> > bitrow;
    typedef std::vector<bitrow> bitmatrix;

    class vertex { // vertex class
        int m_subgraph;
        // used for traversing
        bool m_visited;
        bool m_leaf;
        int m_low;
        int m_disc;
        int m_ancestor;
        int m_color;
        // used for planar embedding
        bool m_embedded;
        int m_number;
        std::map<int,int> m_faces;
        // *
        attrib *m_attributes;
        ivector m_neighbors;
        std::map<int,attrib> *m_neighbor_attributes;
        std::map<int,int> m_multiedges;
        void assign_defaults();
        void assign(const vertex &other);
    public:
        vertex(bool support_attributes=true);
        vertex(const vertex &other);
        vertex(const gen &lab,const attrib &attr);
        ~vertex();
        vertex& operator =(const vertex &other);
        gen label() const;
        bool supports_attributes() const { return m_attributes!=NULL; }
        void unsupport_attributes() { m_attributes=NULL; m_neighbor_attributes=NULL; }
        void set_label(const gen &s) { assert(supports_attributes()); (*m_attributes)[_GT_ATTRIB_LABEL]=s; }
        void set_subgraph(int s) { m_subgraph=s; }
        int subgraph() const { return m_subgraph; }
        void set_embedded(bool yes) { m_embedded=yes; }
        bool is_embedded() const { return m_embedded; }
        void set_number(int n) { m_number=n; }
        int number() const { return m_number; }
        void set_visited(bool yes) { m_visited=yes; }
        bool is_visited() const { return m_visited; }
        void set_low(int l) { m_low=l; }
        int low() const { return m_low; }
        void set_disc(int t) { m_disc=t; }
        int disc() const { return m_disc; }
        void set_ancestor(int i) { m_ancestor=i; }
        void unset_ancestor() { m_ancestor=-1; }
        int ancestor() const { return m_ancestor; }
        void set_color(int c) { m_color=c; }
        int color() const { return m_color; }
        void set_leaf(bool yes) { m_leaf=yes; }
        bool is_leaf() const { return m_leaf; }
        const attrib &attributes() const { assert(supports_attributes()); return *m_attributes; }
        attrib &attributes() { assert(supports_attributes()); return *m_attributes; }
        void set_attribute(int key,const gen &val) { assert(supports_attributes()); (*m_attributes)[key]=val; }
        void set_attributes(const attrib &attr) { assert(supports_attributes()); copy_attributes(attr,*m_attributes); }
        const ivector &neighbors() const { return m_neighbors; }
        int degree() const { return m_neighbors.size(); }
        void add_neighbor(int i,const attrib &attr=attrib());
        bool is_temporary(int i) const;
        attrib &neighbor_attributes(int i);
        const attrib &neighbor_attributes(int i) const;
        bool has_neighbor(int i) const { return binary_search(m_neighbors.begin(),m_neighbors.end(),i); }
        void remove_neighbor(int i);
        void clear_neighbors();
        void map_neighbors(const std::map<int,int> &m);
        void incident_faces(ivector &F) const;
        void add_edge_face(int nb,int f);
        void clear_edge_faces() { m_faces.clear(); }
        int edge_face(int nb) { return m_faces[nb]-1; }
        const std::map<int,int> &edge_faces() const { return m_faces; }
        void set_multiedge(int v,int k);
        int multiedges(int v) const;
        int multiedge_count() const;
        void clear_multiedges() { m_multiedges.clear(); }
        bool has_multiedges() const { return !m_multiedges.empty(); }
    };

    class dotgraph { // temporary structure used in dot parsing
        int m_index;
        attrib vertex_attr,edge_attr,chain_attr;
        ivector m_chain;
        int pos;
    public:
        dotgraph();
        dotgraph(const dotgraph &other);
        dotgraph(int i);
        dotgraph& operator =(const dotgraph &other);
        void assign(const dotgraph &other);
        int index() const { return m_index; }
        void set_index(int i) { m_chain[pos]=i; }
        const attrib &vertex_attributes() const { return vertex_attr; }
        const attrib &edge_attributes() const { return edge_attr; }
        const attrib &chain_attributes() const { return chain_attr; }
        attrib &vertex_attributes() { return vertex_attr; }
        attrib &edge_attributes() { return edge_attr; }
        attrib &chain_attributes() { return chain_attr; }
        const ivector &chain() const { return m_chain; }
        ivector &chain() { return m_chain; }
        int position() const { return pos; }
        void incr() { ++pos; if (int(m_chain.size())<=pos) m_chain.resize(pos+1,0); }
        void clear_chain() { pos=0; m_chain.resize(1); m_chain.front()=0; chain_attr.clear(); }
        bool chain_completed() { return m_chain.back()!=0; }
        bool chain_empty() { return pos==0 && m_chain.front()==0; }
    };

    class walker { // tree node positioner
        graphe *G;
        layout *x;
        double hsep,vsep;
        ivectors levels;
        ivector node_counters,gap_counters,position,gaps,children;
        dvector prelim,modifier;
        std::queue<int> placed;
        int depth;
        void walk(int i,int pass,int level=0,double modsum=0);
        void process_level(int i);
    public:
        walker(graphe *gr,layout *ly,double hs,double vs);
        void positioning(int apex);
    };

    class circ_enum { // circuit enumeration in digraphs
        graphe *G;
        ivector point_stack;
        std::stack<int> marked_stack;
        ivectors A,res;
        bvector mark;
        int s,lb,ub;
        void backtrack(int v,bool &f);
    public:
        circ_enum(graphe *gr,int lo=-1,int hi=-1);
        ivectors find_cycles();
    };

#ifdef HAVE_LIBGLPK
    class painter { // vertex painter
        graphe *G;
        ivectors values;
        ipairs col2ij;
        bvector iscliq,fracts;
        ivector cover_number,initially_colored,branch_candidates,temp_colors,ordering;
        std::set<int> used_colors;
        int lb,ub,maxiter,nxcols;
        bool generate_clique_cuts,select_blb;
        glp_prob *mip;
        double timeout,*heur,*row_coeffs,*best_coeffs;
        int *row_indices,*best_indices;
        void compute_bounds(const ivector &icol,int max_colors);
        void make_values();
        void formulate_mip();
        void fixed_coloring(glp_tree *tree);
        int assign_heur(glp_tree *tree);
        void generate_rows(glp_tree *tree);
        void add_row(glp_prob *prob,int len,int *indices,double *coeffs,double rh);
        static void callback(glp_tree *tree,void *info);
    public:
        painter(graphe *gr) { G=gr; }
        int color_vertices(ivector &colors,const ivector &icol,int max_colors=0,int tm_lim=0,bool verbose=false);
        int select_branching_variable(glp_tree *tree);
        void heur_solution(glp_tree *tree);
    };

    class tsp { // traveling salesman
        struct arc {
            /* arc struct holds only the edge information relevant for TSP */
            int head;
            int tail;
        };
        enum solution_status {
            _GT_TSP_OPTIMAL,
            _GT_TSP_NOT_HAMILTONIAN,
            _GT_TSP_ERROR
        };
        enum heuristic_type {
            _GT_TSP_NO_HEUR                     = 0,
            _GT_TSP_CHRISTOFIDES_SA             = 1,
            _GT_TSP_FARTHEST_INSERTION_HEUR     = 2,
            _GT_TSP_FARTHEST_INSERTION_RANDOM   = 3
        };
        graphe *G;
        glp_prob *mip;
        bool isweighted,*visited,is_symmetric_tsp,verbose,cancellable;
        std::set<ivector> subtours;
        ivector tour,old_sol;
        double *coeff,gap_tol;
        arc *arcs;
        int *indices,nv,ne,heur_type;
        solution_status status;
        std::map<int,std::map<int,double> > weight_map,rlx_sol_map;
        std::map<int,std::map<int,int> > loc_map;
        dvector xev,obj;
        bvector can_branch;
        void formulate_mip();
        bool get_subtours();
        void add_subtours(const ivectors &sv);
        bool find_tours(int k,ivectors &sv,solution_status &status);
        bool subtours_equal(const ivector &st1,const ivector &st2);
        ivector canonical_subtour(const ivector &subtour);
        void append_sce(const ivector &subtour);
        ipair make_edge(int i,int j) const;
        int edge_index(const ipair &e);
        double weight(int i,int j);
        double weight(const ipair &e) { return weight(e.first,e.second); }
        double lower_bound();
        bool is_move_feasible(int k,const ivector &t,const ipairs &x);
        bool lin_kernighan(ivector &hc);
        bool perform_3opt_moves(ivector &hc);
        void straighten(ivector &hc);
        void improve_tour(ivector &hc,bool do_3opt=true);
        void farthest_insertion(int index,ivector &hc);
        bool christofides(ivector &hc,bool show_progress=false);
        static void sample_mean_stddev(const dvector &sample,double &mean,double &stddev);
        bool min_weight_matching_bipartite(const ivector &eind,const dvector &weights,ivector &matched_arcs,bool msg=false);
        void select_branching_variable(glp_tree *tree);
        void rowgen(glp_tree *tree);
        void heur(glp_tree *tree);
        static void callback(glp_tree *tree,void *info);
        void min_wpm_heur(glp_tree *tree,const ivector &eind);
        static void min_wpm_callback(glp_tree *tree,void *info);
        /* min-cut routines, originally written by Andrew Makhorin (<mao@gnu.org>) */
        ivectors mincut_data;
        int max_flow(int nn,int nedg,const ivector &beg,
                     const ivector &end,const ivector &cap,int s,int t,
                     ivector &x);
        int min_st_cut(int nn, int nedg,const ivector &beg,
                       const ivector &end,const ivector &cap,int s,int t,
                       const ivector &x,ivector &cut);
        int minimal_cut(int nn,int nedg,const ivector &beg,
                        const ivector &end,const ivector &cap,ivector &cut);
    public:
        tsp(graphe *gr,double gap_tolerance=0,bool is_verbose=false);
        ~tsp();
        int solve(int k,ivectors &hcv,dvector &costs);
        int solve(ivector &hc,double &cost);
        bool approx(ivector &hc,double &ratio);
        double tour_cost(const ivector &hc);
    };
    
    class atsp { // asymmetric traveling salesman problem
        graphe *G;
        glp_prob *mip;
        ipairs mia; // must include arcs
        ivectors ft; // forbidden tours
        bool isweighted,select_blb,verbose,terminated;
        double gap_tol;
        void formulate_mip();
        static void callback(glp_tree *tree,void *info);
    public:
        atsp(graphe *gr,const ipairs &must_include_arcs,double gap_tolerance=0,bool is_verbose=false);
        ~atsp();
        bool solve(ivector &hc,double &cost); // find the shortest tour
        void ksolve(int k,ivectors &hcv,dvector &costs); // find k shortest tours
    };

    class mvc_solver { // exact and approx minimum vertex cover (MVC) problem solver
        graphe *G;
        glp_prob *ilp;
        static void callback(glp_tree *tree,void *info);
        void preprocess(glp_tree *tree);
        void branch(glp_tree *tree);
        int heuristic(glp_tree *tree);
        int lower_bound(int s);
        bool is_vertex_fixed(glp_prob *p,int j,bool &in_cover);
        void make_vertex_fixed(glp_prob *p,int j,bool in_cover);
        void find_mirrors(int v);
        void packing(glp_tree *tree);
        ivector mirrors,V,V_pos;
        double *heur_sol;
        bool compute_heur,is_k_vc;
        ipairs edges;
        int sg,last_row;
    public:
        mvc_solver(graphe *gr,int s=-1);
        ~mvc_solver();
        int solve(ivector &cover,int k=-1,int tm_lim=0,double gap_tol=0,bool verbose=false);
    };

    typedef struct { double rhs, pi; } mcf_v_data;
    typedef struct { double low, cap, cost, x; } mcf_a_data;
    typedef struct { int set; } wbm_v_data;
    typedef struct { double cost; int x; } wbm_a_data;
#endif

    class rectangle { // simple rectangle class
        double m_x,m_y,m_width,m_height;
        bool m_locked_above,m_locked_right;
        layout *L;
    public:
        struct comparator {
            bool operator()(const rectangle &r1,const rectangle &r2) const {
                return r1.height()>r2.height();
            }
        };
        rectangle();
        rectangle(double X,double Y,double W,double H,layout *ly=NULL);
        rectangle(const rectangle &rect);
        rectangle& operator =(const rectangle &other);
        void assign(const rectangle &other);
        void set_anchor(double x,double y) { m_x=x; m_y=y; }
        double x() const { return m_x; }
        double y() const { return m_y; }
        double width() const { return m_width; }
        double height() const { return m_height; }
        void set_locked_above(bool yes) { m_locked_above=yes; }
        void set_locked_right(bool yes) { m_locked_right=yes; }
        bool is_locked_above() const { return m_locked_above; }
        bool is_locked_right() const { return m_locked_right; }
        bool intersects(const rectangle &other) const;
        bool intersects(const std::vector<rectangle> &rectangles) const;
        bool intersects(std::vector<rectangle>::const_iterator first,std::vector<rectangle>::const_iterator last) const;
        layout *get_layout() const { return L; }
    };

    class cpol {
        void assign(const cpol &other);
    public:
        int nv;
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
        ulong *cg;
        int *col;
#else
        int *adj;
#endif
        size_t sz;
        int frq;
        intpoly poly;
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
        cpol() { nv=0; cg=NULL; col=NULL; sz=0; frq=0; }
#else
        cpol() { nv=0; adj=NULL; sz=0; frq=0; }
#endif
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
        cpol(int n,ulong *g,int *c,size_t s,const intpoly &p);
#else
        cpol(int n,int *a,size_t s,const intpoly &p);
#endif
        cpol(const cpol &other);
        ~cpol() { }
        cpol& operator =(const cpol &other);
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
        bool match(int n,ulong *g,int *c) const;
#else
        bool match(int n,int *a,int adj_sz) const;
#endif
    };

    class ransampl { // random sampling from a given degree distribution
        int n;
        vecteur prob;
        ivector alias;
        const context *ctx;
    public:
        ransampl(const vecteur &p,GIAC_CONTEXT);
        gen data() const;
        int generate() const;
    };

    class bucketsampler { // random sampling from a dynamic distribution
        const context *ctx;
        ivector weights;
        long total_weight;
        std::map<int,int> max_weight,level_weight;
        std::map<int,ivector> levels;
        ipairs positions;
        int nearest_pow2(double a) { return std::floor(0.5+std::log(a)/M_LN2); }
    public:
        bucketsampler(const ivector &W,GIAC_CONTEXT);
        int generate();
        void insert(int w);
        void update(int i,int w);
        void increment(int i) { update(i,weights[i]+1); }
    };

    class unionfind { // disjoint-set data structure
        struct element {
            int id,parent,rank;
        };
        int sz;
        element *elements;
    public:
        unionfind(int n);
        ~unionfind() { delete[] elements; }
        void make_set(int id);
        int find(int id);
        void unite(int id1,int id2);
        void select(int id);
        void clear();
    };

    class ostergard { // clique maximizer
        graphe *G;
        int maxsize;
        bool found,timed_out;
        double timeout; // seconds
        ivector c,incumbent,clique_nodes;
        clock_t start;
        void recurse(ivector &U,int size,ivector &position);
    public:
        ostergard(graphe *gr,double max_time=0) { G=gr; timeout=max_time; }
        int maxclique(ivector &clique);
    };

    class yen { // Yen's k shortest paths algorithm
        typedef struct tree_node {
            int i;
            bool selected;
            tree_node *parent;
            std::vector<tree_node*> children;
        } tree_node;
        struct kspaths_comparator {
            bool operator()(const std::pair<gen,tree_node*> &a,const std::pair<gen,tree_node*> &b) const {
                if (giac::is_zero(a.first-b.first))
                    return a.second<b.second;
                return is_strictly_greater(b.first,a.first,context0);
            }
        };
        graphe *G;
        tree_node *root;
        int src, dest, K;
        std::vector<tree_node*> kspaths;
        tree_node *add_tree_node(tree_node *p);
        tree_node *store_path(const ivector &path,tree_node *r);
        void select_path(tree_node *p);
        void restore_path(tree_node *p,ivector &path);
        void delete_children(tree_node *r);
    public:
        yen(graphe *gr,int s,int d,int k) { G=gr; src=s; dest=d; K=k; root=NULL; }
        ~yen();
        void find_kspaths(ivectors &paths);
    };

    class mm { // An efficient implementation of Edmonds' blossom algorithm
        enum label_t { EVEN=0, ODD=1 };
        graphe *G;
        int *mate,*label,*pred,*bridge,V,s;
        std::queue<int> Q;
        unionfind *ds;
        ivector ap;
        bool alternating_forest();
        bool alternating_tree(int v);
        bool examine(int v,int w);
        void extend_tree(int v,int w);
        void shrink_blossom(int v,int w);
        void shrink_path(int b,int v,int w);
        void augmenting_path(int v,int w);
        ivector find_path(int s,int t);
        int find_base(int v,int w);
        int find_root(int v);
    public:
        mm(graphe *g);
        ~mm();
        void find_maximum_matching(ipairs &matching,int sg=-1);
    };

    class sip { // subgraph isomorphism via Ullmann's algorithm
        int N,n,max_sg,nb;
        bitmatrix A,AT,M;
        bvectors B;
        bvector used_cols;
        ivectors found;
        bool induced,dir,isom;
        ipairs changes;
        std::stack<ipairs> snapshots;
        static void clear_bitrow(bitrow &row);
        static bool get_j(const bitrow &row,int j) { return row[j/SIP_NBITS].test(j%SIP_NBITS); }
        static bool get_ij(const bitmatrix &mat,int i,int j) { return get_j(mat[i],j); }
        static void set_ij(bitmatrix &mat,int i,int j,bool val) { mat[i][j/SIP_NBITS].set(j%SIP_NBITS,val); }
        bool mult_bitrows(const bitrow &r,const bitrow &c);
        bool M_has_empty_row();
        bool recurse(int i=0);
        void prune();
        void edit(int i,int j) { set_ij(M,i,j,false); changes.push_back(std::make_pair(i,j)); }
        void push_changes() { snapshots.push(changes); changes.clear(); }
        void revert_changes();
        bool is_isomorphism();
    public:
        sip(const graphe &G,const graphe &P,int max_sg_in=1);
        int find_subgraphs(bool only_induced=true);
        const ivector &get_subgraph(int i) const { return found[i]; }
    };
    
    struct edges_comparator { // for sorting edges by their weight
        graphe *G;
        bool operator()(const ipair &a,const ipair &b) const {
            return is_strictly_greater(G->weight(b),G->weight(a),G->giac_context());
        }
        edges_comparator(graphe *gr) { G=gr; }
    };

    struct ivectors_comparator { // for sorting ivectors by their length
        bool operator()(const ivector &a,const ivector &b) const {
            return a.size()<b.size();
        }
    };

    struct degree_comparator { // for sorting vertices by their degrees
        graphe *G;
        bool asc;
        bool operator()(int v,int w) const {
            return (asc && G->degree(v)<G->degree(w)) ||
                    (!asc && G->degree(v)>G->degree(w));
        }
        degree_comparator(graphe *gr,bool ascending=true) { G=gr; asc=ascending; }
    };

    struct ivectors_degree_comparator { // for sorting sets of vertices by ascending total degree
        graphe *G;
        bool operator()(const ivector &a,const ivector &b) const {
            int deg_a=0,deg_b=0;
            for (ivector_iter it=a.begin();it!=a.end();++it) {
                deg_a+=G->degree(*it);
            }
            for (ivector_iter it=b.begin();it!=b.end();++it) {
                deg_b+=G->degree(*it);
            }
            return deg_a*b.size()>=deg_b*a.size();
        }
        ivectors_degree_comparator(graphe *gr) { G=gr; }
    };

    /* iterators */
    typedef std::vector<vertex>::const_iterator node_iter;
    typedef std::map<int,attrib>::const_iterator neighbor_iter;
    typedef attrib::const_iterator attrib_iter;
    typedef std::vector<point>::const_iterator layout_iter;
    typedef ivector::const_iterator ivector_iter;
    typedef ivectors::const_iterator ivectors_iter;
    typedef ipairs::const_iterator ipairs_iter;
    typedef edgeset::const_iterator edgeset_iter;
    typedef intpoly::const_iterator intpoly_iter;
    typedef dvector::const_iterator dvector_iter;
    typedef iset::const_iterator iset_iter;

    /* static variables and constants */
    static const gen FAUX;
    static const gen VRAI;
    static bool verbose;
    static int default_vertex_color;
    static int default_edge_color;
    static int default_vertex_label_color;
    static int default_highlighted_edge_color;
    static int default_highlighted_vertex_color;
    static int default_edge_width;
    static int bold_edge_width;
    static std::map<ivector,std::vector<cpol> > cache;
    // special graphs
    static const int clebsch_graph[];
    static const int coxeter_graph[];
    static const char* coxeter_graph_vnames[];
    static const int dyck_graph[];
    static const int grinberg_graph[];
    static const int grotzsch_graph[];
    static const int harries_graph_lcf[];
    static const int harries_wong_graph_lcf[];
    static const int balaban_10cage_lcf[];
    static const int balaban_11cage_lcf[];
    static const int heawood_graph[];
    static const int herschel_graph[];
    static const int mcgee_graph[];
    static const int pappus_graph[];
    static const int robertson_graph[];
    static const int soccer_ball_graph[];
    static const int tetrahedral_graph[];
    static const int octahedral_graph[];
    static const int icosahedral_graph[];
    static const int ljubljana_graph_lcf[];
    static const int foster_graph_lcf[];
    static const int blanusa_graph[];
    static const int bidiakis_cube_graph_lcf[];
    static const int bull_graph[];
    static const int butterfly_graph[];
    static const int diamond_graph[];
    static const int chvatal_graph[];
    static const int franklin_graph_lcf[];
    static const int frucht_graph_lcf[];
    static const int biggs_smith_graph_lcf[];
    static const int moser_spindle_graph[];
    static const int errera_graph[];
    static const int goldner_harary_graph[];
    static const int golomb_graph[];
    static const int hoffman_graph_matrix[8][8];
    static const int poussin_graph[];
    static const int wagner_graph[];
    static const int folkman_graph_lcf[];
    static const int gray_graph_lcf[];
    static const int tutte_12cage_lcf[];
    static const int tutte_8cage_lcf[];
    static const int f26a_graph_lcf[];
    static const int tutte_graph[];
    static const int brinkmann_graph[];
    static const int barnette_bosak_lederberg_graph[];
    static const int double_star_snark[];
    static const int doyle_graph[];
    static const int meringer_graph[];
    static const int robertson_wegner_graph[];
    static const int wong_graph[];
    static const char* const gewirtz_words[];
    static const int harborth_graph[];
    static const int kittell_graph[];
    static const int krackhardt_kite_graph[];
    static const int meredith_graph[];
    static const int perkel_graph[];
    static const int sousselier_graph[];
    static const int walther_graph[];
    static const int watkins_snark[];
    static const int wells_graph[];
    static const int wiener_araya_graph[];
    static const int markstroem_graph[];

    struct spcgraph { const char* name; int nv; int ne; };
    static const spcgraph special_graph[];

    static const std::ios::iostate cout_rdstate;
    static const std::ios::iostate cerr_rdstate;

private:
    const context *ctx;
    std::vector<vertex> nodes;
    attrib attributes;
    std::vector<std::string> user_tags;
    ivector marked_nodes;
    int disc_time;
    ivector disc_nodes;
    std::stack<ipair> edge_stack;
    std::stack<int> node_stack;
    std::queue<int> node_queue;
    std::map<int,iset> visited_edges;
    ivectors maxcliques;
    std::stack<ivector> saved_subgraphs;
    bool m_supports_attributes;
    void clear_node_stack();
    void clear_node_queue();
    void message(const char *str) const;
    void message(const char *format,int a) const;
    void message(const char *format,int a,int b) const;
    void message(const char *format,int a,int b,int c) const;
    void suspend_logging();
    void restore_logging();
    std::string giac_version() const;
    vertex &node(int i) { return nodes[i]; }
    bool dot_parse_attributes(std::ifstream &dotfile,attrib &attr);
    static bool insert_attribute(attrib &attr,int key,const gen &val,bool overwrite=true);
    static bool remove_attribute(attrib &attr,int key);
    static void copy_attributes(const attrib &src,attrib &dest);
    void attrib2genmap(const attrib &attr,gen_map &m,bool keys2tags=false) const;
    void write_attrib(std::ofstream &dotfile,const attrib &attr,bool style=true) const;
    static ivector_iter binsearch(ivector_iter first,ivector_iter last,int a);
    static size_t sets_union(const iset &A,const iset &B,iset &U);
    static size_t sets_intersection(const iset &A,const iset &B,iset &I);
    static size_t sets_intersection(const ivector &A,const iset &B,iset &I);
    static size_t sets_difference(const iset &A,const iset &B,iset &D);
    static size_t sets_difference(const iset &A,const ivector &B,iset &D);
    static size_t sets_difference(const ivector &A,const ivector &B,iset &D);
    static point make_point(double x,double y) { point p(2,x); p.back()=y; return p; }
    static point make_point(double x,double y,double z) { point p(3,x); p[1]=y; p.back()=z; return p; }
    static void add_point(point &a,const point &b);
    static void subtract_point(point &a,const point &b);
    static void scale_point(point &p,double s);
    static double point_vecprod2d(const point &v,const point &w);
    static double point_dotprod(const point &p,const point &q);
    static void clear_point_coords(point &p) { std:fill(p.begin(),p.end(),0); }
    static double point_displacement(const point &p,bool sqroot=true);
    static double point_distance(const point &p,const point &q,point &pq);
    static void point_mirror(double a,double b,double c,const point &src,point &dest);
    static void point_lincomb(const point &p,const point &q,double d1,double d2,point &res);
    static void copy_point(const point &src,point &dest);
    void rand_point(point &p,double radius=1.0);
    static vecteur point2vecteur(const point &p);
    static bool points_coincide(const point &p,const point &q,double tol);
    static void copy_layout(const layout &src,layout &dest);
    static void rotate_layout(layout &x,double phi);
    static double layout_min(const layout &x,int d);
    static double layout_diameter(const layout &x);
    static void point2polar(point &p,double &r,double &phi);
    vecteur draw_edge(int i,int j,const layout &x) const;
    static bool sparse_matrix_element(const sparsemat &A,int i,int j,ipair &val);
    static void multiply_sparse_matrices(const sparsemat &A,const sparsemat &B,sparsemat &P,int ncols,bool symmetric=false);
    static gen sparse_product_element(const sparsemat &A, const sparsemat &B,int i,int j);
    static void transpose_sparsemat(const sparsemat &A,sparsemat &T);
    void multilevel_recursion(layout &x,int d,double R,double K,double tol,int depth=0);
    int mdeg(const ivector &V,int i) const;
    void coarsening(graphe &G,const sparsemat &P,const ivector &V) const;
    void tomita(iset &R,iset &P,iset &X,std::map<int,int> &m,int mode=0);
    int cp_maxclique(ivector &clique);
    void cp_recurse(ivector &C,ivector &P,ivector &incumbent);
    int ost_maxclique(ivector &clique);
    void ost_recursive(ivector &U,int size,int &maxsize,ivector &incumbent,bool &found);
    void find_cut_vertices_dfs(int i,std::set<int> &ap,int sg);
    void find_blocks_dfs(int i,std::vector<ipairs> &blocks,int sg);
    void find_bridges_dfs(int i,ipairs &B,int sg);
    int find_cycle_dfs(int i,int sg);
    bool find_path_dfs(int dest,int i,int sg,bool skip_embedded);
    static void sort_rectangles(std::vector<rectangle> &rectangles);
    static bool embed_rectangles(std::vector<rectangle> &rectangles,double maxheight);
    static bool segments_crossing(const point &p,const point &r,const point &q,const point &s,point &crossing);
    static bool point2segment_projection(const point &p,const point &q,const point &r,point &proj);
    void force_directed_placement(layout &x,double K,double R=DBL_MAX,double tol=0.01,bool ac=true);
    static bool get_node_position(const attrib &attr,point &p,GIAC_CONTEXT);
    void coarsening_mis(const ivector &V,graphe &G,sparsemat &P) const;
    void coarsening_ec(const ipairs &M,graphe &G,sparsemat &P) const;
    static int face_has_edge(const ivector &face,int i,int j);
    int first_neighbor_from_subgraph(const vertex &v,int sg) const;
    void set_nodes_embedded(const ivector &v,bool yes=true);
    void unembed_all_nodes();
    int planar_embedding(ivectors &faces);
    void set_embedding(const ivectors &faces);
    void clear_embedding();
    int choose_embedding_face(const ivectors &faces,int v);
    int choose_outer_face(const ivectors &faces);
    static void make_regular_polygon_layout(layout &x,const ivector &v,double R=1.0,double elongate=0.0);
    bool edges_crossing(const ipair &e1,const ipair &e2,const layout &x,point &crossing) const;
    static void build_block_tree(int i,ivectors &blocks);
    static int common_element(const ivector &v1,const ivector &v2,int offset=0);
    void embed_children_blocks(int i,ivectors &block_tree,std::vector<ivectors> &blocks_faces);
    void periphericity(const ivector &outer_face,ivector &p);
    void tree_height_dfs(int i,int level,int &depth);
    void make_product_nodes(const graphe &G,graphe &P) const;
    static void extract_path_from_cycle(const ivector &cycle,int i,int j,ivector &path);
    static void generate_nk_sets(int n,int k,std::vector<std::bitset<32> > &v);
    void strongconnect_dfs(ivectors &components,bvector &onstack,int i,int sg);
    bool degrees_equal(const ivector &v,int deg=0) const;
    void lca_recursion(int u,const ipairs &p,ivector &lca,unionfind &ds);
    void st_numbering_dfs(int i,ivector &preorder);
    void rdfs(int i,ivector &d,bool rec,int sg,bool skip_embedded);
    bool is_descendant(int v,int anc) const;
    static int pred(int i,int n);
    static int succ(int i,int n);
    static void arc_path(int i,int j,const ivector &cycle,ivector &path);
    void fold_face(const ivector &face,bool subdivide,int &label);
    void find_chords(const ivector &face,ipairs &chords);
    void augment(const ivectors &faces,int outer_face,bool subdivide=false);
    int saturation_degree(const vertex &v,std::set<int> &colors) const;
    int uncolored_degree(const vertex &v) const;
    bool is_partially_colored() const;
    void remove_maximal_clique(iset &V) const;
    bool bipartite_matching_bfs(ivector &dist,int sg=-1);
    bool bipartite_matching_dfs(int u,ivector &dist,int sg=-1);
    static gen make_colon_label(const ivector &v);
    void simplify(graphe &G,bool color_temp_vertices=false) const;
    intpoly tutte_poly_recurse(int vc);
    static void poly_mult(intpoly &a,const intpoly &b);
    static void poly_add(intpoly &a,const intpoly &b);
    static intpoly poly_geom(int var,int k,bool leading_one,bool add_other_var=false);
    static intpoly poly_one();
    static gen intpoly2gen(const intpoly &v,const gen &x,const gen &y);
    void sharc_order();
    static void compute_rutcounts(int n,vecteur &t);
    void ranrut(int n,ivector &tree,const vecteur &pt);
    void ranrut_forest(int m,ivectors &trees,const vecteur &alpha,const vecteur &a);
    void insert_tree(const ivector &tree,int root);
    static ipair rat2ipair(const gen &r);
    static gen ipair2rat(const ipair &p);
    void save_subgraphs();
    void restore_subgraphs();
    int vertex_pair_connectivity(int v,int w);
    static gen harmonic_mean_exact(gen a,gen b,gen c) { return 3*a*b*c/(a*b+b*c+a*c); }
    static double harmonic_mean(double a,double b,double c) { return 3.0*a*b*c/(a*b+b*c+a*c); }
    void strec(int i,int t,int counter,int np,iset &Q,vecteur &timestamp,vecteur &l);
    bool hamcycle_recurse(ivector &path,int pos);
    void grasp_construct(double aplha, ivector &Q,bool cmpl,int sg);
    void grasp_local(ivector &Q,bool cmpl,int sg);
    bool mvc_special(ivector &cover,int sg);
    bool mvc_is_unconfined(int i,int sg=0) const;
    bool mvc_is_dominant(int v,int sg) const;
    bool mvc_reduce_basic(int sg,int c);
    void mvc_half_integral(int sg, ivector &in_cover, ivector &out_cover);
    void mvc_alom(ivector &cover,int sg=-1);
    void mvc_dfs(ivector &cover,int sg=-1);
    void mvc_bipartite(const ivector &U,const ivector &V,ivector &cover,int sg=-1,bool conn=false);
    void alom_candidates(const ivector &V,ivector &cand,int sg);
    int count_edges_in_Nv(int v,int sg=-1) const;
    int count_edges(const ivector &V) const;
    bool is_simplicial(int i,const sparsemat &A,double D=0.0);

public:
    graphe(const context *contextptr=context0,bool support_attributes=true);
    graphe(const graphe &G);
    graphe(const std::string &name,const context *contextptr=context0,bool support_attributes=true);
    static graphe *from_gen(const gen &g);
    virtual ~graphe() { }
    graphe &operator =(const graphe &other);
    bool is_simple() const;
    virtual std::string print (GIAC_CONTEXT) const;
    virtual std::string texprint (GIAC_CONTEXT) const;
    virtual bool operator ==(const gen &) const;
    virtual gen_user *memory_alloc() const { return new graphe(*this); }
    virtual gen giac_constructor(GIAC_CONTEXT) const;
    
    // methods
    int rand_integer(int n) const { assert(n>=0); return n==0?0:giac::giac_rand(ctx)%n; }
    double rand_uniform() const;
    double rand_normal() const { return giac::randNorm(ctx); }
    ivector rand_permu(int n) const;
    static gen to_binary(int number,int chars);
    const context *giac_context() const { return ctx; }
    static gen make_idnt(const char* name,int index=-1,bool intern=true);
    void make_default_labels(vecteur &labels,int n,int n0=0,int offset=-1) const;
    bool labels2iset(const vecteur &labels,iset &s);
    static gen boole(bool b) { return b?VRAI:FAUX; }
    static std::string int2string(int i);
    static gen str2gen(const std::string &str,bool isstring=false);
    static std::string genstring2str(const gen &g);
    static gen plusinf();
    void ivectors2vecteur(const ivectors &v,vecteur &res,bool sort_all=false) const;
    void append_segment(vecteur &drawing, const point &p,const point &q,int color,int width,int style,bool arrow=false) const;
    void append_node(vecteur &drawing,const point &p,int color,int width,int shape) const;
    void append_label(vecteur &drawing,const point &p,const gen &label,int quadrant,int color=_BLACK) const;
    void reserve_nodes(int n) { assert(nodes.empty()); nodes.reserve(n); }
    //bool read_gen(const gen &g);
    void read_special(const int *special_graph);
    void read_special(const char **special_graph);
    void copy(graphe &G) const;
    void copy_nodes(const std::vector<vertex> &V);
    void copy_nodes(graphe &G,std::map<int,int> &vmap,int sg=-1) const;
    bool supports_attributes() const { return m_supports_attributes; }
    void clear();
    void clear_maximal_cliques() { maxcliques.clear(); }
    void find_maximal_cliques();
    const ivectors &maximal_cliques() const { return maxcliques; }
    int tag2index(const std::string &tag);
    std::string index2tag(int index) const;
    int register_user_tag(const std::string &tag);
    void register_user_tags(const std::vector<std::string> &tags);
    const ivector &get_marked_nodes() const { return marked_nodes; }
    void get_marked_nodes(vecteur &V) const;
    void get_marked_nodes_in_subgraph(int s,ivector &m) const;
    void copy_marked_nodes(const ivector &mv) { marked_nodes=ivector(mv.begin(),mv.end()); }
    void mark_node(int v);
    void mark_node(const gen &v) { mark_node(node_index(v)); }
    bool unmark_node(int v);
    bool unmark_node(const gen &v) { return unmark_node(node_index(v)); }
    void unmark_all_nodes() { marked_nodes.clear(); }
    void sort_marked_nodes() { std::sort(marked_nodes.begin(),marked_nodes.end()); }
    void set_edge_visited(int i,int j);
    void set_edge_visited(const ipair &e) { set_edge_visited(e.first,e.second); }
    bool is_edge_visited(int i,int j) const;
    bool is_edge_visited(const ipair &e) const { return is_edge_visited(e.first,e.second); }
    void unvisit_all_edges() { visited_edges.clear(); }
    //gen to_gen();
    int *to_array(int &sz,bool colored,bool reduce=false) const;
    bool write_latex(const std::string &filename,const gen &drawing) const;
    bool write_dot(const std::string &filename,bool style=true) const;
    bool write_lst(const std::string &filename) const;
    bool read_dot(const std::string &filename);
    bool read_lst(const std::string &filename);
    bool is_null() const { return nodes.empty(); }
    bool is_empty(int sg=-1) const;
    void weight_matrix(matrice &W) const;
    gen weight(int i,int j) const;
    gen weight(const ipair &edge) const { return weight(edge.first,edge.second); }
    int edge_count(int sg=-1) const;
    int node_count() const { return nodes.size(); }
    vecteur vertices(int sg=-1) const;
    void unvisit_all_nodes(int sg=-1);
    void unset_all_ancestors(int sg=-1);
    void uncolor_all_nodes(int base_color=0,int sg=-1);
    void set_node_color(int i,int c) { node(i).set_color(c); }
    void dfs(int root,bool rec=true,bool clr=true,ivector *D=NULL,int sg=-1,bool skip_embedded=false);
    void bfs(int root,bool rec=true,bool clr=true,ivector *D=NULL,int sg=-1,bool skip_embedded=false);
    const ivector &get_discovered_nodes() const { return disc_nodes; }
    bool is_connected(int sg=-1);
    bool is_biconnected(int sg=-1);
    bool is_triconnected(int sg=-1);
    bool is_cycle(ipairs &E,int sg=-1);
    void adjacent_nodes(int i,ivector &adj,bool include_temp_edges=true) const;
    void get_edges_as_pairs(ipairs &E,int sg=-1) const;
    vecteur edges(bool include_weights,int sg=-1) const;
    ivector edge_multiplicities() const;
    int sum_of_edge_multiplicities() const;
    int add_node();
    int add_node(const gen &v,const attrib &attr=attrib());
    void add_nodes(const vecteur &v);
    void add_nodes(int n);
    void add_unlabeled_nodes(int n);
    void remove_isolated_nodes(const iset &I);
    void isolate_nodes(const iset &V);
    const vertex &node(int i) const { assert(i>=0 && i<node_count()); return nodes[i]; }
    const gen node_label(int i) const { assert(i>=0 && i<node_count()); return nodes[i].label(); }
    vecteur get_node_labels(const ivector &v=ivector(0)) const;
    int node_index(const gen &v) const;
    int edge_index(const ipair &e) const;
    int largest_integer_label() const;
    void set_subgraph(const ivector &v,int s);
    void set_subgraph(const ipairs &e,int s);
    void get_subgraph(int sg,ivector &v) const;
    int subgraph_size(int sg) const;
    int first_vertex_from_subgraph(int sg) const;
    void merge_subgraphs(int s,int t);
    void unset_subgraphs(int default_sg=-1);
    int max_subgraph_index() const;
    const attrib &graph_attributes() const { return attributes; }
    const attrib &node_attributes(int i) const;
    const attrib &edge_attributes(int i,int j) const;
    attrib &node_attributes(int i);
    attrib &edge_attributes(int i,int j);
    const attrib &edge_attributes(const ipair &e) const { return edge_attributes(e.first,e.second); }
    attrib &edge_attributes(const ipair &e) { return edge_attributes(e.first,e.second); }
    void attrib2vecteurs(const attrib &attr,vecteur &tags,vecteur &values) const;
    bool add_edge(int i,int j,const gen &w=gen(1));
    bool add_edge(int i,int j,const attrib &attr);
    bool add_edge(const ipair &edge) { return add_edge(edge.first,edge.second); }
    bool add_edge(const ipair &edge,const attrib &attr) { return add_edge(edge.first,edge.second,attr); }
    bool add_edge(const gen &v,const gen &w,const gen &weight=gen(1));
    bool add_edge(const gen &v,const gen &w,const attrib &attr);
    void add_temporary_edge(int i,int j);
    bool is_temporary_edge(int i,int j) const;
    void remove_temporary_edges();
    bool remove_edge(int i,int j);
    bool remove_edge(const ipair &p) { return remove_edge(p.first,p.second); }
    bool has_edge(int i,int j) const;
    bool has_edge(const ipair &p) const { return has_edge(p.first,p.second); }
    ipair make_edge(const vecteur &v) const;
    bool edges2ipairs(const vecteur &E,ipairs &ev,bool &notfound) const;
    vecteur ipairs2edges(const ipairs &E) const;
    static void ipairs2edgeset(const ipairs &E,edgeset &Eset);
    bool nodes_are_adjacent(int i,int j) const;
    int in_degree(int index,int sg=-1) const;
    int out_degree(int index,int sg=-1) const;
    int degree(int index,int sg=-1) const;
    int maximum_degree(int sg=-1) const;
    int minimum_degree(int sg=-1) const;
    ivector degree_sequence(int sg=-1) const;
    void sort_by_degrees(graphe &G,ivector &sigma,int sg=-1);
    void adjacency_matrix(matrice &m) const;
    void adjacency_sparse_matrix(sparsemat &sm,bool diag_ones=false,int sg=-1) const;
    void laplacian_matrix(matrice &m,bool normalize=false) const;
    void incidence_matrix(matrice &m) const;
    void set_graph_attribute(int key,const gen &val) { attributes[key]=val; }
    void set_graph_attributes(const attrib &attr) { copy_attributes(attr,attributes); }
    void set_node_attribute(int index,int key,const gen &val);
    void set_edge_attribute(int i,int j,int key,const gen &val);
    bool get_graph_attribute(int key,gen &val) const;
    bool get_node_attribute(int i,int key,gen &val) const;
    bool get_edge_attribute(int i,int j,int key,gen &val) const;
    void discard_graph_attribute(int key);
    void discard_node_attribute(int i,int key);
    void discard_edge_attribute(int i,int j,int key);
    void set_name(const std::string &str) { set_graph_attribute(_GT_ATTRIB_NAME,str2gen(str,true)); }
    std::string name() const { gen s; if (get_graph_attribute(_GT_ATTRIB_NAME,s)) return genstring2str(s); else return ""; }
    bool is_directed() const;
    bool is_weighted() const;
    void set_directed(bool yes) { set_graph_attribute(_GT_ATTRIB_DIRECTED,boole(yes)); }
    void set_weighted(bool yes) { set_graph_attribute(_GT_ATTRIB_WEIGHTED,boole(yes)); }
    void make_weighted(const matrice &m);
    void make_directed() { set_directed(true); }
    void make_unweighted();
    void randomize_edge_weights(double a,double b,bool integral_weights=false);
    void make_lcf_graph(const ivector &jumps,int e);
    int is_regular(int d) const;
    bool is_strongly_regular(ipair &sig);
    bool is_equal(const graphe &G) const;
    void underlying(graphe &G,int sg=-1) const;
    void complement(graphe &G,int sg=-1) const;
    bool isomorphic_copy(graphe &G,const ivector &sigma,bool strip_attributes=false,int sg=-1);
    bool relabel_nodes(const vecteur &labels);
    void induce_subgraph(const ivector &vi,graphe &G) const;
    void extract_subgraph(const ipairs &E,graphe &G) const;
    void subgraph_indices(const graphe &S,const ivector &v,ivector &sv) const;
    bool is_subgraph(const graphe &G) const;
    void renumber_subgraph(int sg,int sg_new);
    void maximal_independent_set(ivector &ind) const;
    void maximum_matching(ipairs &matching,int sg=-1,gt_conn_check cc=_GT_CC_FIND_COMPONENTS);
    void maximal_matching(ipairs &matching,int sg=-1) const;
    bool demoucron(ivectors &faces,int sg=-1);
    void create_random_layout(layout &x,int dim);
    void make_spring_layout(layout &x,int d,double tol=0.001);
    void make_circular_layout(layout &x,const ivector &hull,double A=0,double tol=0.005,double elongate=0.0);
    void make_tutte_layout(layout &x,const ivector &outer_face);
    bool make_planar_layout(layout &x,double *score=NULL);
    bool make_best_planar_layout(layout &x,int iter);
    void make_tree_layout(layout &x,double sep,int apex=0);
    void make_bipartite_layout(layout &x,const ivector &p1,const ivector &p2);
    void layout_best_rotation(layout &x);
    bool is_tree(int sg=-1) { return !is_directed() && edge_count(sg)+1==subgraph_size(sg) && is_connected(sg); }
    bool is_forest(int sg=-1,gt_conn_check cc=_GT_CC_FIND_COMPONENTS);
    bool is_tournament(int sg=-1) const;
    bool is_planar(int sg=-1);
    bool is_clique(int sg=-1) const;
    gen triangle_count(ivectors *dest=NULL,bool ccoeff=false,bool exact=true);
    int tree_height(int root=-1);
    void clique_stats(std::map<int,int> &m,int mode=0);
    int maximum_clique(ivector &clique);
    void greedy_neighborhood_clique_cover_numbers(ivector &cover_numbers);
    bool clique_cover(ivectors &cover,int k=0,int tm_lim=0,bool verbose=false);
    int maximum_independent_set(ivector &v) const;
    int girth(bool odd=false,int sg=-1);
    bool hakimi(const ivector &L);
    void erdos_renyi(double p);
    void preferential_attachment(int d,int o);
    void molloy_reed(const vecteur &p);
    void make_plane_dual(const ivectors &faces);
    void make_lcf_graph(const int *j);
    void make_sierpinski_graph(int n,int k,bool triangle);
    void make_shrikhande_graph();
    void make_complete_graph();
    void make_complete_multipartite_graph(const ivector &partition_sizes,layout *x=NULL);
    void make_petersen_graph(int n,int k,layout *x=NULL);
    bool make_intersection_graph(int n,int k,int sz,bool lab_sets=false);
    void make_path_graph();
    void make_cycle_graph();
    void make_grid_graph(int m,int n,int mode=0);
    void make_web_graph(int n,int m,layout *x=NULL);
    void make_wheel_graph(int n,layout *x=NULL);
    void make_antiprism_graph(int n,layout *x=NULL);
    void make_complete_kary_tree(int k,int d);
    bool make_flower_snark(int n,layout *x=NULL);
    bool make_goldberg_snark(int n);
    bool make_haar_graph(const gen &n);
    bool make_haar_graph(ulong n);
    void make_paley_graph(int p,int k);
    void make_hypercube_graph(int n);
    void make_hoffman_singleton_graph();
    void make_sylvester_graph();
    void make_higman_sims_graph();
    void make_brouwer_haemers_graph();
    void make_gewirtz_graph();
    void make_schlaefli_graph();
    void make_gosset_graph();
    void make_szerekes_snark();
    void make_random_tree(int maxd=0);
    void make_random_rooted_tree();
    void make_random_free_tree();
    void make_random_planar(double p,int connectivity);
    void make_random_sequential(const ivector &d);
    void make_random_bipartite(int a,int b,double p);
    void make_random_regular(int d,bool connected);
    static void translate_layout(layout &x,const point &dx);
    void cartesian_product(const graphe &G,graphe &P) const;
    void tensor_product(const graphe &G,graphe &P) const;
    void connected_components(ivectors &components,int sg=-1,bool skip_embedded=false,int *count=NULL);
    int connected_component_count(int sg=-1);
    int connected_components_to_subgraphs(int sg=-1);
    void biconnected_components(ivectors &components,int sg=-1);
    void strongly_connected_components(ivectors &components,int sg=-1);
    bool has_cut_vertex(int sg=-1,int i=0);
    void find_cut_vertices(ivector &articulation_points,int sg=-1);
    void find_blocks(std::vector<ipairs> &blocks,int sg=-1);
    void find_bridges(ipairs &B,int sg=-1);
    void find_ears(ivectors &ears,int sg=-1);
    bool find_cycle(ivector &cycle,int sg=-1);
    bool find_path(int i,int j,ivector &path,int sg=-1,bool skip_embedded=false);
    bool find_eulerian_trail(ivector &path);
    int eulerian_trail_start(bool &iscycle) const;
    bool fleury(int start,ivector &path);
    void hierholzer(ivector &path);
    bool is_multigraph() const;
    void make_simple();
    int multiedges(const ipair &e) const;
    void set_multiedge(const ipair &e,int k);
    bool weights2multiedges();
    void contract_edge(int i,int j,bool adjust_positions=true);
    void contract_edge(const ipair &e,bool adjust_pos=true) { contract_edge(e.first,e.second,adjust_pos); }
    void subdivide_edge(const ipair &e,int n,int &label);
    void incident_edges(const ivector &V,edgeset &E) const;
    static bool edges_incident(const ipair &e1,const ipair &e2);
    void convex_hull(const layout &x,layout &hull);
    void edge_labels_placement(const layout &x);
    void draw_edges(vecteur &drawing,const layout &x);
    void draw_nodes(vecteur &drawing,const layout &x) const;
    void draw_labels(vecteur &drawing,const layout &x) const;
    void distance(int i,const ivector &J,ivector &dist,ivectors *shortest_paths=NULL);
    void allpairs_distance(matrice &m);
    void dijkstra(int src,const ivector &dest,vecteur &path_weights,ivectors *cheapest_paths=NULL,int sg=-1);
    bool bellman_ford(int src,const ivector &dest,vecteur &path_weights,ivectors *cheapest_paths=NULL);
    bool topologic_sort(ivector &ordering);
    bool is_arborescence() const;
    void reverse(graphe &G) const;
    void spanning_tree(int i,graphe &T,int sg=-1);
    void minimal_spanning_tree(graphe &T,int sg=-1);
    void lowest_common_ancestors(int root,const ipairs &p,ivector &lca_recursion);
    int lowest_common_ancestor(int i,int j,int root);
    void compute_st_numbering(int s,int t);
    vecteur get_st_numbering() const;
    void assign_edge_directions_from_st();
    void parametrized_st_orientation(int s,int t,double p);
    void greedy_vertex_coloring_biggs(ivector &ordering);
    int greedy_vertex_coloring(const ivector &p);
    int exact_vertex_coloring(int max_colors=0,int tm_lim=0,bool verbose=false);
    int exact_edge_coloring(ivector &colors,int *numcol=NULL,int tm_lim=0,bool verbose=false);
    bool edge_coloring_heuristic(ivector &colors);
    int get_node_color(int i) const;
    void get_node_colors(ivector &colors) const;
    bool is_bipartite(ivector &V1,ivector &V2,int sg=-1,gt_conn_check cc=_GT_CC_FIND_COMPONENTS);
    bool is_vertex_colorable(int k);
    void dsatur();
    int color_count() const;
    ipair adjacent_color_count(int i) const;
    bool adjacent_colors(int i,std::set<int> &colors) const;
    ipair chromatic_number_bounds();
    void store_layout(const layout &x);
    bool has_stored_layout(layout &x) const;
    void bipartite_matching(const ivector &p1,const ivector &p2,ipairs &matching,int sg=-1,gt_conn_check cc=_GT_CC_FIND_COMPONENTS);
    bool weighted_bipartite_matching(const ivector &p1,const ivector &p2,ipairs &matching,
        bool minimize=true,double eps=1e-5,int sg=-1,gt_conn_check cc=_GT_CC_FIND_COMPONENTS);
    void line_graph(graphe &G,ipairs &E) const;
    void transitive_closure(graphe &G,bool weighted=false);
    bool is_isomorphic(graphe &other,std::map<int,int> &isom,bool use_colors=true);
    gen aut_generators() const;
    bool canonical_labeling(ivector &lab) const;
    bool bondy_chvatal_closure(graphe &G,ivector &d);
    int hamcond(bool make_closure=true);
    bool is_hamiltonian(ivector &hc);
    bool hamcycle(ivector &path);
    int traveling_salesman(int k,ivectors &hcv,dvector &costs,double gap_tol=0,bool verbose=true);
    bool find_directed_tours(int k,ivectors &hcv,dvector &costs,const ipairs &incl,double gap_tol=0,bool verbose=false);
    bool make_euclidean_distances();
    gen maxflow_edmonds_karp(int s,int t,std::vector<std::map<int,gen> > &flow,const gen &limit=plusinf());
    void minimum_cut(int s,const std::vector<std::map<int,gen> > &flow,ipairs &cut);
    gen tutte_polynomial(const gen &x,const gen &y);
    void fundamental_cycles(ivectors &cycles,int sg=-1,gt_conn_check=_GT_CC_FIND_COMPONENTS);
    void mycielskian(graphe &G) const;
    gen local_clustering_coeff(int i) const;
    gen clustering_coeff(bool approx,bool exact);
    gen transitivity();
    int edge_connectivity();
    int vertex_connectivity();
    void truncate(graphe &dest,const ivectors &faces);
    void condensation(graphe &G);
    void elementary_cycles(ivectors &cyc,int lo,int hi);
    void yen_ksp(int K,int src,int dest,ivectors &paths);
    void compute_in_out_degrees(ivector &ind,ivector &outd) const;
    vecteur distances_from(int k);
    gen betweenness_centrality(int k) const;
    gen communicability_betweenness_centrality(int k) const;
    gen closeness_centrality(int k,bool harmonic=false);
    gen degree_centrality(int k) const;
    vecteur katz_centrality(const gen &att) const;
    gen information_centrality(int k,bool approx=false) const;
    int splittance(int &m,ivector &vseq) const;
    bool is_split_graph(ivector &clq,ivector &indp) const;
    void contract_subgraph(graphe &G,const ivector &sg,const gen &lb) const;
    void grasp_clique(int maxitr,ivector &clq,bool cmpl=false,int sg=-1);
    bool mvc(ivector &cover,int vc_alg,int sg=-1,int tm_lim=0,double gap_tol=0,bool verbose=false);
    int k_vertex_cover(ivector &cover,int k);
    int vertex_cover_number(int sg=-1);
    bool is_reachable(int u,int v);
    void reachable(int u,ivector &r);
    void find_simplicial_vertices(ivector &res);
    void subgraph_isomorphism(graphe &P,int max_sg,bool induced,ivectors &res) const;
    bool pruefer_encode(ivector &code,bool check_tree=false);
    bool pruefer_decode(const ivector &code);
    const char* identify();
    void identify_from_sequences(vecteur &spec,int haar_limit=25);

    // static methods
    static gen customize_display(int options);
    static gen customize_layout(int options);
    static rectangle layout_bounding_rect(layout &ly,double padding=0,int coord_offset=0);
    static void pack_rectangles(std::vector<rectangle> &rectangles);
    static gen point2gen(const point &p,bool vect=false);
    static bool gen2point(const gen &g,point &p,GIAC_CONTEXT);
    static point layout_center(const layout &x);
    static void scale_layout(layout &x,double diam);
    static void scale_layout_2d(layout &x,double sx,double sy);
    static int best_quadrant(const point &p,const layout &x);
    static gen colon_label(int i,int j);
    static gen colon_label(int i,int j,int k);
    static size_t intersect_linear(ivector_iter min1,ivector_iter max1,ivector_iter min2,ivector_iter max2);
    static bool is_graphic_sequence(const ivector &s_orig);
    static ivector_iter insert_sorted(ivector &V,int val);
    static bool erase_sorted(ivector &V,int val);
    static ivector make_ivector(int n,...);
    static double poly_area(const layout &x);
    static int term_hook(void *info,const char *s);
    static bool is_interrupted();
};

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
#endif // __GRAPHE_H
