/*
 * graphe.cc
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

#include "giacPCH.h"
#include "giac.h"
#include "graphe.h"
#include "graphtheory.h"
#include "signalprocessing.h"
#include "optimization.h"
#include <ctype.h>
#include <math.h>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_LIBNAUTY
#include "nautywrapper.h"
#endif
using namespace std;

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

const gen graphe::VRAI=gen(1).change_subtype(_INT_BOOLEAN);
const gen graphe::FAUX=gen(0).change_subtype(_INT_BOOLEAN);
bool graphe::verbose=true;
int graphe::default_edge_color=_BLUE;
int graphe::default_edge_width=_LINE_WIDTH_2;
int graphe::bold_edge_width=_LINE_WIDTH_4;
int graphe::default_highlighted_edge_color=_RED;
int graphe::default_highlighted_vertex_color=_GREEN;
int graphe::default_vertex_color=_YELLOW;
int graphe::default_vertex_label_color=_BLACK;

/* messages and logging */

void graphe::message(const char *str) const {
    if (verbose)
        *logptr(ctx) << gettext(str) << "\n";
}

void graphe::message(const char *format,int a) const {
    char buffer[256];
    sprintf(buffer,gettext(format),a);
    if (verbose)
        *logptr(ctx) << buffer << "\n";
}

void graphe::message(const char *format,int a,int b) const {
    char buffer[256];
    sprintf(buffer,gettext(format),a,b);
    if (verbose)
        *logptr(ctx) << buffer << "\n";
}

void graphe::message(const char *format,int a,int b,int c) const {
    char buffer[256];
    sprintf(buffer,gettext(format),a,b,c);
    if (verbose)
        *logptr(ctx) << buffer << "\n";
}

void graphe::suspend_logging() {
    ctx->globalptr->_logptr_->setstate(ios::ios_base::failbit);
}

void graphe::restore_logging() {
    ios::ios_base::iostate s=ctx->globalptr->_logptr_->rdstate() & ~ios::ios_base::failbit;
    ctx->globalptr->_logptr_->clear(s);
}

string graphe::giac_version() const {
    return genstring2str(_version(change_subtype(vecteur(0),_SEQ__VECT),ctx));
}

/* print functions for debugging */

ostream& operator<<(ostream &os, const graphe::ivector &v)
{
    os << "[";
    for (graphe::ivector_iter it=v.begin();it!=v.end();++it) {
        if (it!=v.begin()) os << ",";
        os << *it;
    }
    os << "]";
    return os;
}

ostream& operator<<(ostream &os, const graphe::ivectors &v)
{
    os << "[";
    for (graphe::ivectors_iter it=v.begin();it!=v.end();++it) {
        if (it!=v.begin()) os << ",";
        os << *it;
    }
    os << "]";
    return os;
}

ostream& operator<<(ostream &os, const graphe::ipairs &v)
{
    os << "[";
    for (graphe::ipairs_iter it=v.begin();it!=v.end();++it) {
        if (it!=v.begin()) os << ",";
        os << "(" << it->first << "," << it->second << ")";
    }
    os << "]";
    return os;
}

/* helper functions */

gen graphe::plusinf() {
    return symbolic(at_plus,_IDNT_infinity());
}

double graphe::rand_uniform() const {
    gen args=change_subtype(vecteur(0),_SEQ__VECT);
    return _rand(args,ctx).to_double(ctx);
}

double graphe::poly_area(const layout &x) {
    double a=0;
    int n=x.size();
    for (int i=0;i<n;++i) {
        const point &p1=x[i],&p2=x[(i+1)%n];
        assert(p1.size()==2 && p2.size()==2);
        a+=p1.front()*p2.back()-p1.back()*p2.front();
    }
    return fabs(a)/2;
}

int graphe::term_hook(void *info,const char *s) {
    *static_cast<std::ostream*>(info) << s;
    return true;
}

bool graphe::is_interrupted() {
    if (interrupted || ctrl_c) {
        interrupted=ctrl_c=false;
        return true;
    }
    return false;
}

/* make a list of integers and return it */
graphe::ivector graphe::make_ivector(int n,...) {
    va_list lst;
    va_start(lst,n);
    ivector res(n);
    for (int i=0;i<n;++i) {
        res[i]=va_arg(lst,int);
    }
    return res;
}

/* convert list of lists of integers into vecteur of vecteurs */
void graphe::ivectors2vecteur(const ivectors &v,vecteur &res,bool sort_all) const {
    res.resize(v.size());
    for (ivectors_iter it=v.begin();it!=v.end();++it) {
        res[it-v.begin()]=sort_all?sort_identifiers(get_node_labels(*it),ctx):get_node_labels(*it);
    }
}

/* vertex class implementation */
void graphe::vertex::assign_defaults() {
    m_subgraph=-1;
    m_visited=false;
    m_ancestor=-1;
    m_low=-1;
    m_disc=-1;
    m_color=0; // white
    m_embedded=false;
    m_number=-1;
}

graphe::vertex::vertex(bool support_attributes) {
    m_attributes=support_attributes?new attrib():NULL;
    m_neighbor_attributes=support_attributes?new map<int,attrib>():NULL;
    assign_defaults();
}

graphe::vertex::vertex(const gen &lab,const attrib &attr) {
    m_attributes=new attrib();
    m_neighbor_attributes=new map<int,attrib>();
    assign_defaults();
    if (!attr.empty())
        set_attributes(attr);
    set_label(lab);
}

void graphe::vertex::assign(const vertex &other) {
    m_subgraph=other.subgraph();
    m_visited=other.is_visited();
    m_ancestor=other.ancestor();
    m_low=other.low();
    m_disc=other.disc();
    m_color=other.color();
    m_embedded=other.is_embedded();
    m_number=other.number();
    m_faces=other.edge_faces();
    if (!other.supports_attributes()) {
        if (supports_attributes()) {
            delete m_attributes;
            delete m_neighbor_attributes;
        }
        unsupport_attributes();
    } else {
        if (!supports_attributes()) {
            m_attributes=new attrib();
            m_neighbor_attributes=new map<int,attrib>();
        }
        m_neighbor_attributes->clear();
        set_attributes(other.attributes());
    }
    m_multiedges.clear();
    m_neighbors.resize(other.degree());
    int k;
    for (ivector_iter it=other.neighbors().begin();it!=other.neighbors().end();++it) {
        m_neighbors[it-other.neighbors().begin()]=*it;
        if (other.supports_attributes())
            copy_attributes(other.neighbor_attributes(*it),(*m_neighbor_attributes)[*it]);
        if ((k=other.multiedges(*it))>0)
            m_multiedges.insert(make_pair(*it,k));
    }
}

graphe::vertex::~vertex() {
    if (supports_attributes()) {
        delete m_attributes;
        delete m_neighbor_attributes;
    }
}

int graphe::vertex::multiedges(int v) const {
    map<int,int>::const_iterator it=m_multiedges.find(v);
    if (it!=m_multiedges.end())
        return it->second;
    return 0;
}

int graphe::vertex::multiedge_count() const {
    int count=0;
    for (map<int,int>::const_iterator it=m_multiedges.begin();it!=m_multiedges.end();++it) {
        count+=it->second;
    }
    return count;
}

void graphe::vertex::set_multiedge(int v,int k) {
    map<int,int>::iterator it=m_multiedges.find(v);
    if (k>0)
        if (it!=m_multiedges.end())
            it->second=k;
        else m_multiedges.insert(make_pair(v,k));
    else {
        if (it!=m_multiedges.end())
            m_multiedges.erase(it);
    }
}

gen graphe::vertex::label() const {
    assert(supports_attributes());
    attrib_iter it=m_attributes->find(_GT_ATTRIB_LABEL);
    if (it==m_attributes->end())
        return undef;
    return it->second;
}

graphe::vertex::vertex(const vertex &other) {
    m_attributes=other.supports_attributes()?new attrib():NULL;
    m_neighbor_attributes=other.supports_attributes()?new map<int,attrib>():NULL;
    assign(other);
}

graphe::vertex& graphe::vertex::operator =(const vertex &other) {
    assign(other);
    return *this;
}

void graphe::vertex::add_neighbor(int i,const attrib &attr) {
    if (m_neighbors.empty())
        m_neighbors.push_back(i);
    else {
        ivector::iterator it=std::lower_bound(m_neighbors.begin(),m_neighbors.end(),i);
        if (it!=m_neighbors.end() && *it==i)
            return;
        m_neighbors.insert(it,i);
    }
    if (supports_attributes())
        copy_attributes(attr,(*m_neighbor_attributes)[i]);
}

bool graphe::vertex::is_temporary(int i) const {
    const attrib &attr=neighbor_attributes(i);
    attrib_iter it=attr.find(_GT_ATTRIB_TEMPORARY);
    return it!=attr.end() && giac::is_one(it->second);
}

graphe::attrib &graphe::vertex::neighbor_attributes(int i) {
    assert(supports_attributes());
    map<int,attrib>::iterator it=m_neighbor_attributes->find(i);
    assert(it!=m_neighbor_attributes->end());
    return it->second;
}

const graphe::attrib &graphe::vertex::neighbor_attributes(int i) const {
    assert(supports_attributes());
    neighbor_iter it=m_neighbor_attributes->find(i);
    assert(it!=m_neighbor_attributes->end());
    return it->second;
}

void graphe::vertex::remove_neighbor(int i) {
    ivector::iterator it;
    if ((it=std::find(m_neighbors.begin(),m_neighbors.end(),i))==m_neighbors.end())
        return;
    m_neighbors.erase(it);
    if (supports_attributes()) {
        map<int,attrib>::iterator jt=m_neighbor_attributes->find(i);
        assert(jt!=m_neighbor_attributes->end());
        m_neighbor_attributes->erase(jt);
    }
    map<int,int>::iterator kt=m_multiedges.find(i);
    if (kt!=m_multiedges.end())
        m_multiedges.erase(kt);
}

void graphe::vertex::clear_neighbors() {
    m_neighbors.clear();
    if (supports_attributes())
        m_neighbor_attributes->clear();
    m_multiedges.clear();
}

void graphe::vertex::map_neighbors(const map<int,int> &m) {
    ivector::iterator it=m_neighbors.begin(),itend=m_neighbors.end();
    map<int,attrib>::iterator jt;
    map<int,int>::iterator mt;
    attrib attr;
    bool sa=supports_attributes();
    for (;it!=itend;++it) {
        if (sa) {
            jt=m_neighbor_attributes->find(*it);
            assert(jt!=m_neighbor_attributes->end());
        }
        mt=m_multiedges.find(*it);
        *it=m.at(*it);
        if (sa) {
            copy_attributes(jt->second,attr);
            m_neighbor_attributes->erase(jt);
            copy_attributes(attr,(*m_neighbor_attributes)[*it]);
        }
        if (mt!=m_multiedges.end()) {
            int me=mt->second;
            m_multiedges.erase(mt);
            m_multiedges[*it]=me;
        }
    }
}

void graphe::vertex::incident_faces(ivector &F) const {
    F.resize(m_faces.size());
    int i=0,f;
    for (map<int,int>::const_iterator it=m_faces.begin();it!=m_faces.end();++it) {
        assert((f=it->second)>0);
        F[i++]=f-1;
    }
}

void graphe::vertex::add_edge_face(int nb,int f) {
    assert(m_faces.find(nb)==m_faces.end());
    m_faces[nb]=f+1;
}

/* set the given planar embedding */
void graphe::set_embedding(const ivectors &faces) {
    int n,f=0,i;
    for (ivectors_iter it=faces.begin();it!=faces.end();++it,++f) {
        const ivector &face=*it;
        n=face.size();
        for (i=0;i<n;++i) {
            vertex &v=node(face[i]);
            v.add_edge_face(face[(i+1)%n],f);
        }
    }
}

/* clear the previously set embedding */
void graphe::clear_embedding() {
    for (std::vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->clear_edge_faces();
    }
}

/* get the first neighbor of v in the subgraph sg */
int graphe::first_neighbor_from_subgraph(const vertex &v,int sg) const {
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        if (node(*it).subgraph()==sg)
            return *it;
    }
    return -1;
}

/* just a safety measure, the stack should be empty before calling this function */
void graphe::clear_node_stack() {
    while (!node_stack.empty())
        node_stack.pop();
}

/* clear the node queue */
void graphe::clear_node_queue() {
    while (!node_queue.empty())
        node_queue.pop();
}

/* return true iff this graph has no edges */
bool graphe::is_empty(int sg) const {
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if ((sg<0 || it->subgraph()==sg) && !it->neighbors().empty())
            return false;
    }
    return true;
}

/* store all subgraph indices */
void graphe::save_subgraphs() {
    ivector sgr(node_count());
    node_iter it=nodes.begin(),itend=nodes.end();
    int i=0;
    for (;it!=itend;++it,++i) {
        sgr[i]=it->subgraph();
    }
    saved_subgraphs.push(sgr);
}

/* restore previously saved subgraph indices */
void graphe::restore_subgraphs() {
    assert(!saved_subgraphs.empty() && int(saved_subgraphs.top().size())==node_count());
    const ivector &sgr=saved_subgraphs.top();
    int i=0;
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it,++i) {
        it->set_subgraph(sgr[i]);
    }
    saved_subgraphs.pop();
}

/* dotgraph class implementation */
graphe::dotgraph::dotgraph() {
    m_index=0;
    pos=0;
    m_chain=ivector(1,0);
}

graphe::dotgraph::dotgraph(int i) {
    m_index=i;
    pos=0;
    m_chain=ivector(1,0);
}

void graphe::dotgraph::assign(const dotgraph &other) {
    m_index=other.index();
    copy_attributes(other.vertex_attributes(),vertex_attr);
    copy_attributes(other.edge_attributes(),edge_attr);
    copy_attributes(other.chain_attributes(),chain_attr);
    m_chain=other.chain();
    pos=other.position();
}

graphe::dotgraph::dotgraph(const dotgraph &other) {
    assign(other);
}

graphe::dotgraph& graphe::dotgraph::operator =(const dotgraph &other) {
    assign(other);
    return *this;
}

/* rectangle class implementation */
graphe::rectangle::rectangle() {
    m_x=m_y=m_width=m_height=0;
    L=NULL;
}

graphe::rectangle::rectangle(double X,double Y,double W,double H,layout *ly) {
    m_x=X;
    m_y=Y;
    m_width=W;
    m_height=H;
    L=ly;
}

void graphe::rectangle::assign(const rectangle &other) {
    m_x=other.x();
    m_y=other.y();
    m_width=other.width();
    m_height=other.height();
    L=other.get_layout();
}

graphe::rectangle::rectangle(const rectangle &rect) {
    assign(rect);
}

graphe::rectangle& graphe::rectangle::operator =(const rectangle &other) {
    assign(other);
    return *this;
}

bool graphe::rectangle::intersects(const rectangle &other) const {
    return x()<other.x()+other.width() &&
            x()+width()>other.x() &&
            y()<other.y()+other.height() &&
            y()+height()>other.y();
}

bool graphe::rectangle::intersects(const vector<rectangle> &rectangles) const {
    vector<rectangle>::const_iterator it=rectangles.begin();
    for (;it!=rectangles.end();++it) {
        if (intersects(*it))
            return true;
    }
    return false;
}

bool graphe::rectangle::intersects(vector<rectangle>::const_iterator first,vector<rectangle>::const_iterator last) const {
    vector<rectangle>::const_iterator it=first;
    for (;it!=last;++it) {
        if (intersects(*it))
            return true;
    }
    return false;
}

/* convert number to binary format and return it as gen of type string */
gen graphe::to_binary(int number,int chars) {
    return str2gen(bitset<32>((unsigned long)number).to_string().substr(32-chars),true);
}

/* make a copy of attr */
void graphe::copy_attributes(const attrib &src,attrib &dest) {
    dest.clear();
    for (attrib_iter it=src.begin();it!=src.end();++it) {
        dest.insert(make_pair(it->first,gen(it->second)));
    }
}

/* fill the vecteur V with first n integers (0- or 1- based, depending on the mode) */
void graphe::make_default_labels(vecteur &labels,int n,int n0,int offset) const {
    int ofs=offset<0?array_start(ctx):offset;
    labels.resize(n);
    for (int i=0;i<n;++i) {
        labels[i]=i+ofs+n0;
    }
}

/* create identifier */
gen graphe::make_idnt(const char* name,int index,bool intern) {
    string id;
    if (intern) id+=" ";
    id+=name;
    if (index>=0) id+=int2string(index);
    return identificateur(id.c_str());
}

/* convert integer to string */
std::string graphe::int2string(int i) {
    return printint(i);
}

/* convert string to gen */
gen graphe::str2gen(const string &str,bool isstring) {
    if (isstring) {
        vector<string> words;
        gen space=_char(32,context0);
        char *cstr=new char[str.size()+1];
        strcpy(cstr,str.c_str());
        char *p=strtok(cstr," ");
        while (p!=NULL) {
            words.push_back(p);
            p=strtok(NULL," ");
        }
        vecteur res;
        for (vector<string>::const_iterator it=words.begin();it!=words.end();++it) {
            res.push_back(string2gen(*it,false));
            if (it+1!=words.end())
                res.push_back(space);
        }
        delete[] cstr;
        return _cat(change_subtype(res,_SEQ__VECT),context0);
    }
    return gen(str,context0);
}

/* convert gen string to std::string (no quotes) */
string graphe::genstring2str(const gen &g) {
    assert(g.type==_STRNG);
    int len=_size(g,context0).val;
    return string(g._STRNGptr->begin(),g._STRNGptr->begin()+len);
}

/* return random permutation of {0,1,..,n-1} */
graphe::ivector graphe::rand_permu(int n) const {
    ivector res=vecteur_2_vector_int(*_randperm(n,ctx)._VECTptr);
    int ofs=array_start(ctx);
    for (ivector::iterator it=res.begin();it!=res.end();++it) {
        *it-=ofs;
    }
    return res;
}

/* insert element val in sorted list V keeping it sorted (in ascending order),
 * return an iterator pointing to the inserted element */
graphe::ivector_iter graphe::insert_sorted(ivector &V,int val) {
    if (V.empty()) {
        V.push_back(val);
        return V.begin()+V.size()-1;
    }
    ivector::iterator it=std::lower_bound(V.begin(),V.end(),val);
    int pos=it-V.begin();
    V.insert(it,val);
    return V.begin()+pos;
}

/* erase element val in sorted list V, return true if val is contained in V,
 * else return false and do nothing */
bool graphe::erase_sorted(ivector &V,int val) {
    if (V.empty())
        return false;
    ivector::iterator it=std::lower_bound(V.begin(),V.end(),val);
    if (it==V.end())
        return false;
    if (*it!=val)
        return false;
    V.erase(it);
    return true;
}

/* compute the union of A and B and store it in U */
size_t graphe::sets_union(const iset &A,const iset &B,iset &U) {
    U.clear();
    std::set_union(A.begin(),A.end(),B.begin(),B.end(),std::inserter(U,U.begin()));
    return U.size();
}

/* compute the intersection of A and B and store it in I */
size_t graphe::sets_intersection(const iset &A,const iset &B,iset &I) {
    I.clear();
    std::set_intersection(A.begin(),A.end(),B.begin(),B.end(),std::inserter(I,I.begin()));
    return I.size();
}

/* compute the intersection of A and B and store it in I */
size_t graphe::sets_intersection(const ivector &A,const iset &B,iset &I) {
    I.clear();
    std::set_intersection(A.begin(),A.end(),B.begin(),B.end(),std::inserter(I,I.begin()));
    return I.size();
}

/* compute the set difference of A and B and store it in D */
size_t graphe::sets_difference(const iset &A,const iset &B,iset &D) {
    D.clear();
    std::set_difference(A.begin(),A.end(),B.begin(),B.end(),std::inserter(D,D.begin()));
    return D.size();
}

/* compute the set difference of A and B and store it in D */
size_t graphe::sets_difference(const iset &A,const ivector &B,iset &D) {
    D.clear();
    std::set_difference(A.begin(),A.end(),B.begin(),B.end(),std::inserter(D,D.begin()));
    return D.size();
}

/* return the size of the difference of A and B */
size_t graphe::sets_difference(const ivector &A,const ivector &B,iset &D) {
    D.clear();
    std::set_difference(A.begin(),A.end(),B.begin(),B.end(),std::inserter(D,D.begin()));
    return D.size();
}

/* binary search for an element a in sorted list first...last */
graphe::ivector_iter graphe::binsearch(ivector_iter first,ivector_iter last,int a) {
    ivector_iter mid;
    while (first!=last) {
        mid=first+int(last-first)/2;
        if (*mid==a) return mid;
        if (*mid>a) last=mid;
        else first=mid+1;
    }
    return first;
}

/* return the intersection size for sorted lists A and B */
size_t graphe::intersect_linear(ivector_iter min1,ivector_iter max1,ivector_iter min2,ivector_iter max2) {
    if (min1==max1 || min2==max2 || *min1>*(max2-1) || *min2>*(max1-1)) return 0;
    if (*min1>*min2) min2=binsearch(min2,max2,*min1);
    else if (*min2>*min1) min1=binsearch(min1,max1,*min2);
    size_t result=0;
    while (min1!=max1 && min2!=max2) {
        if (*min1<*min2) ++min1;
        else if (*min2<*min1) ++min2;
        else {
            result++;
            ++min1; ++min2;
        }
    }
    return result;
}

/* graphe default constructor */
graphe::graphe(GIAC_CONTEXT,bool support_attributes) {
    ctx=contextptr;
    m_supports_attributes=support_attributes;
    set_graph_attribute(_GT_ATTRIB_DIRECTED,FAUX);
    set_graph_attribute(_GT_ATTRIB_WEIGHTED,FAUX);
}

/* graphe constructor, create a copy of G */
graphe::graphe(const graphe &G) {
    m_supports_attributes=G.supports_attributes();
    set_graph_attribute(_GT_ATTRIB_DIRECTED,boole(G.is_directed()));
    set_graph_attribute(_GT_ATTRIB_WEIGHTED,boole(G.is_weighted()));
    ctx=G.giac_context();
    G.copy(*this);
}

/* graphe constructor, create special graph with the specified name */
graphe::graphe(const string &name,GIAC_CONTEXT,bool support_attributes) {
    ctx=contextptr;
    m_supports_attributes=support_attributes;
    set_graph_attribute(_GT_ATTRIB_DIRECTED,FAUX);
    set_graph_attribute(_GT_ATTRIB_WEIGHTED,FAUX);
    ivector hull;
    layout x;
    if (name=="clebsch") {
        read_special(clebsch_graph);
        if (support_attributes) {
            vecteur labels;
            for (int i=0;i<16;++i) {
                labels.push_back(to_binary(i,4));
            }
            relabel_nodes(labels);
        }
    } else if (name=="coxeter") {
        if (support_attributes) {
            read_special(coxeter_graph_vnames);
            for (int i=1;i<=7;++i) {
                string ai=string("a")+int2string(i);
                hull.push_back(node_index(str2gen(ai,true)));
            }
            make_circular_layout(x,hull,3.5);
        } else read_special(coxeter_graph);
    } else if (name=="desargues") {
        make_petersen_graph(10,3,support_attributes?&x:NULL);
    } else if (name=="barnette-bosak-lederberg") {
        read_special(barnette_bosak_lederberg_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="bidiakis") {
        make_lcf_graph(bidiakis_cube_graph_lcf);
        if (support_attributes) {
            for (int i=0;i<node_count();++i) hull.push_back(i);
            make_circular_layout(x,hull);
            layout_best_rotation(x);
        }
    } else if (name=="brinkmann") {
        read_special(brinkmann_graph);
    } else if (name=="bull") {
        read_special(bull_graph);
        if (support_attributes) {
            make_spring_layout(x,2);
            layout_best_rotation(x);
        }
    } else if (name=="butterfly") {
        read_special(butterfly_graph);
        if (support_attributes) {
            make_spring_layout(x,2);
            layout_best_rotation(x);
        }
    } else if (name=="blanusa") {
        read_special(blanusa_graph);
        if (support_attributes) {
            for (int i=0;i<6;++i) hull.push_back(i);
            add_temporary_edge(0,9);
            add_temporary_edge(1,9);
            add_temporary_edge(3,15);
            add_temporary_edge(4,15);
            make_circular_layout(x,hull,2.5,0.005,0.618);
            remove_temporary_edges();
        }
    } else if (name=="brouwer-haemers") {
        make_brouwer_haemers_graph();
    } else if (name=="diamond") {
        read_special(diamond_graph);
        if (support_attributes) {
            make_spring_layout(x,2);
            layout_best_rotation(x);
        }
    } else if (name=="double-star") {
        read_special(double_star_snark);
    } else if (name=="doyle") {
        read_special(doyle_graph);
    } else if (name=="chvatal") {
        read_special(chvatal_graph);
        if (support_attributes) {
            x.resize(12);
            for (layout::iterator it=x.begin();it!=x.end();++it) it->resize(2);
            x[0][0]=0.43; x[0][1]=2.51; x[1][0]=2.16; x[1][1]=2.02;
            x[2][0]=3.89; x[2][1]=2.51; x[3][0]=2.16; x[3][1]=1.12;
            x[4][0]=1.25; x[4][1]=3.67; x[5][0]=0.86; x[5][1]=1.34;
            x[6][0]=3.47; x[6][1]=1.34; x[7][0]=3.07; x[7][1]=3.67;
            x[8][0]=2.83; x[8][1]=2.75; x[9][0]=2.83; x[9][1]=0.43;
            x[10][0]=1.49; x[10][1]=0.43; x[11][0]=1.49; x[11][1]=2.75;
        }
    } else if (name=="dodecahedron") {
        make_petersen_graph(10,2,support_attributes?&x:NULL);
        if (support_attributes) {
            hull=make_ivector(5,0,1,2,12,10);
            make_circular_layout(x,hull);
        }
    } else if (name=="errera") {
        read_special(errera_graph);
        if (support_attributes) {
            for (int i=0;i<3;++i) hull.push_back(i);
            make_circular_layout(x,hull);
            layout_best_rotation(x);
        }
    } else if (name=="poussin") {
        read_special(poussin_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="franklin") {
        make_lcf_graph(franklin_graph_lcf);
        if (support_attributes) {
            for (int i=0;i<node_count();++i) hull.push_back(i);
            make_circular_layout(x,hull);
            layout_best_rotation(x);
        }
    } else if (name=="frucht") {
        make_lcf_graph(frucht_graph_lcf);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="biggs-smith") {
        make_lcf_graph(biggs_smith_graph_lcf);
        if (support_attributes) {
            for (int i=0;i<node_count();++i) hull.push_back(i);
            make_circular_layout(x,hull);
        }
    } else if (name=="duerer") {
        make_petersen_graph(6,2,support_attributes?&x:NULL);
    } else if (name=="dyck") {
        read_special(dyck_graph);
    } else if (name=="folkman") {
        make_lcf_graph(folkman_graph_lcf);
    } else if (name=="gewirtz") {
        make_gewirtz_graph();
    } else if (name=="gray") {
        make_lcf_graph(gray_graph_lcf);
    } else if (name=="grinberg") {
        read_special(grinberg_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="goldner-harary") {
        read_special(goldner_harary_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="golomb") {
        read_special(golomb_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="gosset") {
        make_gosset_graph();
    } else if (name=="groetzsch") {
        read_special(grotzsch_graph);
        if (support_attributes) {
            for (int i=0;i<5;++i) hull.push_back(node_index(i));
            make_circular_layout(x,hull,2.5);
        }
    } else if (name=="f26a") {
        make_lcf_graph(f26a_graph_lcf);
    } else if (name=="harries") {
        make_lcf_graph(harries_graph_lcf);
    } else if (name=="harries-wong") {
        make_lcf_graph(harries_wong_graph_lcf);
    } else if (name=="balaban10") {
        make_lcf_graph(balaban_10cage_lcf);
    } else if (name=="balaban11") {
        make_lcf_graph(balaban_11cage_lcf);
    } else if (name=="harborth") {
        read_special(harborth_graph);
        if (support_attributes) {
            make_spring_layout(x,2);
            layout_best_rotation(x);
        }
    } else if (name=="heawood") {
        read_special(heawood_graph);
    } else if (name=="herschel") {
        read_special(herschel_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="higman-sims") {
        make_higman_sims_graph();
    } else if (name=="hoffman") {
        if (support_attributes) {
            vecteur labels;
            make_default_labels(labels,16);
            add_nodes(labels);
        } else add_nodes(16);
        for (int i=0;i<16;++i) for (int j=0;j<16;++j) {
            if ((i<8 && j<8) || (i>=8 && j>=8)) continue;
            if (hoffman_graph_matrix[i<j?j-8:i-8][i<j?i:j]==1)
                add_edge(i,j);
        }
        if (support_attributes) {
            make_spring_layout(x,2);
            layout_best_rotation(x);
        }
    } else if (name=="hoffman-singleton") {
        make_hoffman_singleton_graph();
    } else if (name=="icosahedron") {
        read_special(icosahedral_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="levi") {
        make_lcf_graph(tutte_8cage_lcf);
    } else if (name=="ljubljana") {
        make_lcf_graph(ljubljana_graph_lcf);
    } else if (name=="foster") {
        make_lcf_graph(foster_graph_lcf);
    } else if (name=="kittell") {
        read_special(kittell_graph);
        if (support_attributes) {
            for (int i=0;i<3;++i) hull.push_back(i);
            make_circular_layout(x,hull);
        }
    } else if (name=="krackhardt") {
        read_special(krackhardt_kite_graph);
        if (support_attributes) {
            make_spring_layout(x,2);
            layout_best_rotation(x);
        }
    } else if (name=="markstroem") {
        read_special(markstroem_graph);
        if (support_attributes) {
            for (int i=0;i<9;++i) hull.push_back(i);
            make_circular_layout(x,hull);
        }
    } else if (name=="mcgee") {
        read_special(mcgee_graph);
    } else if (name=="meredith") {
        read_special(meredith_graph);
    } else if (name=="meringer") {
        read_special(meringer_graph);
    } else if (name=="moebius-kantor") {
        make_petersen_graph(8,3,support_attributes?&x:NULL);
    } else if (name=="moser") {
        read_special(moser_spindle_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="nauru") {
        make_petersen_graph(12,5,support_attributes?&x:NULL);
    } else if (name=="octahedron") {
        read_special(octahedral_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="pappus") {
        read_special(pappus_graph);
    } else if (name=="perkel") {
        read_special(perkel_graph);
    } else if (name=="petersen") {
        make_petersen_graph(5,2,support_attributes?&x:NULL);
    } else if (name=="robertson") {
        read_special(robertson_graph);
        if (support_attributes) {
            for (int i=0;i<node_count();++i) hull.push_back(i);
            make_circular_layout(x,hull);
        }
    } else if (name=="robertson-wegner") {
        read_special(robertson_wegner_graph);
    } else if (name=="schlaefli") {
        make_schlaefli_graph();
    } else if (name=="soccerball") {
        graphe G;
        G.read_special(icosahedral_graph);
        ivectors faces;
        assert(G.demoucron(faces));
        G.truncate(*this,faces);
        if (support_attributes) {
            hull=make_ivector(5,4,0,1,2,3);
            make_circular_layout(x,hull,2.5);
        }
    } else if (name=="sousselier") {
        read_special(sousselier_graph);
        if (support_attributes) {
            hull=make_ivector(15,4,12,5,2,13,9,0,8,7,1,10,14,3,6,11);
            make_circular_layout(x,hull);
        }
    } else if (name=="sylvester") {
        make_sylvester_graph();
    } else if (name=="shrikhande") {
        make_shrikhande_graph();
        if (support_attributes) {
            for (int i=0;i<4;++i) {
                add_temporary_edge(i,i+8);
                hull.push_back(i);
            }
            add_temporary_edge(8,10);
            add_temporary_edge(9,11);
            add_temporary_edge(4,6);
            add_temporary_edge(5,7);
            make_circular_layout(x,hull);
            remove_temporary_edges();
        }
    } else if (name=="szerekes") {
        make_szerekes_snark();
        if (support_attributes) {
            for (int i=0;i<5;++i) hull.push_back(45+(i*3)%5);
            make_circular_layout(x,hull);
        }
    } else if (name=="tetrahedron") {
        read_special(tetrahedral_graph);
        if (support_attributes) {
            make_planar_layout(x);
            layout_best_rotation(x);
        }
    } else if (name=="tietze") {
        make_flower_snark(3,support_attributes?&x:NULL);
    } else if (name=="tutte") {
        read_special(tutte_graph);
        if (support_attributes) {
            hull=make_ivector(10,34,26,23,15,20,28,32,42,3,41);
            make_circular_layout(x,hull);
        }
    } else if (name=="tutte12") {
        make_lcf_graph(tutte_12cage_lcf);
    } else if (name=="wagner") {
        read_special(wagner_graph);
        if (support_attributes) {
            hull=make_ivector(8,0,1,2,3,4,5,6,7);
            make_circular_layout(x,hull);
            layout_best_rotation(x);
        }
    } else if (name=="walther") {
        read_special(walther_graph);
    } else if (name=="watkins") {
        read_special(watkins_snark);
    } else if (name=="wells") {
        read_special(wells_graph);
    } else if (name=="wiener-araya") {
        read_special(wiener_araya_graph);
        if (support_attributes) {
            hull=make_ivector(4,41,38,40,39);
            make_circular_layout(x,hull,0.5);
        }
    } else if (name=="wong") {
        read_special(wong_graph);
    }
    if (support_attributes && !x.empty()) {
        double sep=1.0;
        scale_layout(x,sep*std::sqrt((double)node_count()));
        rectangle rect=layout_bounding_rect(x,sep/PLASTIC_NUMBER_3);
        translate_layout(x,make_point(-rect.x(),-rect.y()));
        store_layout(x);
    }
}

/* construct Hoffman-Singleton graph by joining five pentagons and five pentagrams together */
void graphe::make_hoffman_singleton_graph() {
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,50);
        add_nodes(labels);
    } else add_nodes(50);
    /* make pentagons */
    for (int i=0;i<10;++i) {
        for (int j=0;j<5;++j) {
            add_edge(5*i+j,5*i+((j+1)%5));
        }
    }
    /* join pentagons */
    for (int h=0;h<5;++h) { // pentagons P_h
        for (int i=0;i<5;++i) { // pentagrams Q_i
            for (int j=0;j<5;++j) { // j-th vertex of P_h
                int k=(h*i+j)%5; // corresponding vertex in Q_i
                /* join these two vertices */
                add_edge(5*h+j,25+5*i+((k*2)%5));
            }
        }
    }
}

/* construct Sylvester graph from Hoffman-Singleton graph */
void graphe::make_sylvester_graph() {
    make_hoffman_singleton_graph();
    const vertex &v1=node(0),&v2=node(v1.neighbors().front());
    iset V;
    V.insert(0); V.insert(v1.neighbors().front());
    for (ivector_iter it=v1.neighbors().begin();it!=v1.neighbors().end();++it) {
        V.insert(*it);
    }
    for (ivector_iter it=v2.neighbors().begin();it!=v2.neighbors().end();++it) {
        V.insert(*it);
    }
    assert(V.size()==14);
    isolate_nodes(V);
    remove_isolated_nodes(V);
}

/* construct Higman-Sims graph from the independent sets of size 15 in Hoffman-Singleton graph */
void graphe::make_higman_sims_graph() {
    graphe G(ctx,false),GC(ctx,false);
    G.make_hoffman_singleton_graph();
    G.complement(GC);
    ivectors clq;
    clq.reserve(100);
    GC.find_maximal_cliques();
    for (ivectors_iter it=GC.maximal_cliques().begin();it!=GC.maximal_cliques().end();++it) {
        if (it->size()==15)
            clq.push_back(*it);
    }
    assert(clq.size()==100);
    for (ivectors::iterator it=clq.begin();it!=clq.end();++it) {
        std::sort(it->begin(),it->end());
    }
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,100);
        add_nodes(labels);
    } else add_nodes(100);
    for (int i=0;i<100;++i) {
        for (int j=i+1;j<100;++j) {
            int sz=intersect_linear(clq[i].begin(),clq[i].end(),clq[j].begin(),clq[j].end());
            if (sz==0 || sz==8)
                add_edge(i,j);
        }
    }
}

/* construct Brouwer-Haemers graph from the finite field of order 81 */
void graphe::make_brouwer_haemers_graph() {
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,81);
        add_nodes(labels);
    } else add_nodes(81);
    gen g=identificateur(" g"),k=identificateur(" k"),K=identificateur(" K");
    _purge(makesequence(g,k,K),ctx);
    suspend_logging();
    gen G=_galois_field(makesequence(3,4,makevecteur(k,K,g)),ctx),a,b;
    restore_logging();
    for (int i=-1;i<80;++i) {
        a=i<0?gen(0):_eval(pow(g,i),ctx);
        for (int j=i+1;j<80;++j) {
            b=_eval(pow(g,j),ctx);
            if (giac::is_one(_eval(pow(a-b,20),ctx)))
                add_edge(i+1,j+1);
        }
    }
    _purge(makesequence(g,k,K),ctx);
}

/* construct Szerekes snark */
void graphe::make_szerekes_snark() {
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,50);
        add_nodes(labels);
    } else add_nodes(50);
    for (int i=0;i<=45;++i) {
        add_edge(i,(i+1)%45);
    }
    for (int i=0;i<5;++i) {
        add_edge(9*i+2,9*((i+2)%5)+6);
        for (int j=0;j<3;++j) add_edge(45+i,9*i+1+3*j);
        add_edge(9*i,9*i+5);
        add_edge(9*i+3,9*i+8);
    }
}

/* construct Paley graph from the finite field of order p^k, where p is an odd prime */
void graphe::make_paley_graph(int p,int k) {
    int q=_eval(pow(gen(p),k),ctx).val;
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,q);
        add_nodes(labels);
    } else add_nodes(q);
    bool dir=(q-3)%4==0;
    set_directed(dir);
    if (k==1) {
        for (int i=0;i<q;++i) {
            for (int j=dir?0:i+1;j<q;++j) {
                if (i==j) continue;
                if (giac::is_one(_legendre_symbol(makesequence(j-i,p),ctx)))
                    add_edge(i,j);
            }
        }
    } else {
        gen g=identificateur(" g"),kk=identificateur(" k"),K=identificateur(" K");
        _purge(makesequence(g,kk,K),ctx);
        suspend_logging();
        gen G=_galois_field(makesequence(p,k,makevecteur(kk,K,g)),ctx),a,b;
        restore_logging();
        for (int i=-1;i<q-1;++i) {
            a=i<0?gen(0):_eval(pow(g,i),ctx);
            for (int j=dir?-1:i+1;j<q-1;++j) {
                if (i==j) continue;
                b=j<0?gen(0):_eval(pow(g,j),ctx);
                if (giac::is_one(_eval(pow(a-b,(q-1)/2),ctx)))
                    add_edge(i+1,j+1);
            }
        }
        _purge(makesequence(g,kk,K),ctx);
    }
}

/* make n-dimensional cube graph */
void graphe::make_hypercube_graph(int n) {
    this->clear();
    int N=(1<<n);
    reserve_nodes(N);
    if (supports_attributes()) {
        for (int i=0;i<N;++i) {
            add_node(graphe::to_binary(i,n));
        }
    } else add_nodes(N);
    for (int i=0;i<N;++i) {
        for (int j=i+1;j<N;++j) {
            if (_hamdist(makesequence(i,j),ctx).val==1)
                add_edge(i,j);
        }
    }
}

/* construct Gewirtz graph from the list of 56 words with 6 letters each */
void graphe::make_gewirtz_graph() {
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,56);
        add_nodes(labels);
    } else add_nodes(56);
    for (int i=0;i<56;++i) {
        const char *wi=gewirtz_words[i];
        for (int j=i+1;j<56;++j) {
            const char *wj=gewirtz_words[j];
            bool dis=true;
            for (int k=0;dis && k<6;++k) {
                for (int l=0;l<6;++l) {
                    if (wi[k]==wj[l]) {
                        dis=false;
                        break;
                    }
                }
            }
            if (dis) add_edge(i,j);
        }
    }
}

void graphe::make_schlaefli_graph() {
    vecteur v1=makevecteur(0,0,0,0,0,0,1,0);
    vecteur v2=makevecteur(0,0,0,0,0,0,0,1);
    vecteur v3(8,fraction(1,2));
    vecteur v,w;
    for (int i=0;i<6;++i) {
        w=v1;
        w[i]=1;
        v.push_back(w);
        w=v2;
        w[i]=1;
        v.push_back(w);
    }
    for (int i=0;i<5;++i) {
        for (int j=i+1;j<6;++j) {
            w=v3;
            w[i]=-w[i];
            w[j]=-w[j];
            v.push_back(w);
        }
    }
    assert(v.size()==27);
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,27);
        add_nodes(labels);
    } else add_nodes(27);
    for (int i=0;i<27;++i) {
        for (int j=i+1;j<27;++j) {
            if (giac::is_one(scalarproduct(*v[i]._VECTptr,*v[j]._VECTptr,ctx)))
                add_edge(i,j);
        }
    }
}

void graphe::make_gosset_graph() {
    vecteur v,w,v0(8,-1);
    for (int i=0;i<8;++i) {
        for (int j=i+1;j<8;++j) {
            w=v0;
            w[i]=w[j]=3;
            v.push_back(w);
            v.push_back(-w);
        }
    }
    assert(v.size()==56);
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,56);
        add_nodes(labels);
    } else add_nodes(56);
    for (int i=0;i<56;++i) {
        for (int j=i+1;j<56;++j) {
            if (giac::is_zero(scalarproduct(*v[i]._VECTptr,*v[j]._VECTptr,ctx)-8))
                add_edge(i,j);
        }
    }
}

/* return true iff this graph is simple */
bool graphe::is_simple() const {
    int i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        const ivector &ngh=it->neighbors();
        for (ivector_iter jt=ngh.begin();jt!=ngh.end();++jt) {
            const ivector &ngh2=node(*jt).neighbors();
            if (*jt==i || (std::lower_bound(ngh2.begin(),ngh2.end(),i)==ngh2.end())
                       || (jt+1!=ngh.end() && *jt==*(jt+1)))
                return false;
        }
    }
    return true;
}

/* allocate, initialize and return an integer array of adjacency lists of this graph,
 * in form [c1,a11,a12,..,-1,c2,a21,a22,..,-1,..], where c1,c2,... are vertex colors */
int *graphe::to_array(int &sz,bool colored,bool reduce) const {
    assert(!reduce || !is_directed());
    sz=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        sz+=it->neighbors().size();
    }
    if (reduce) {
        assert((sz%2)==0);
        sz/=2;
    }
    sz+=2*node_count();
    int *res=new int[sz],i=0,j=0;
    for (it=nodes.begin();it!=itend;++it,++j) {
        res[i++]=colored?get_node_color(j):0;
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            if (!reduce || *jt>j)
                res[i++]=*jt;
        }
        res[i++]=-1;
    }
    return res;
}

/* return index associated with the specified attribute tag */
int graphe::tag2index(const string &tag) {
    if (tag=="label")
        return _GT_ATTRIB_LABEL;
    if (tag=="directed")
        return _GT_ATTRIB_DIRECTED;
    if (tag=="weighted")
        return _GT_ATTRIB_WEIGHTED;
    if (tag=="weight")
        return _GT_ATTRIB_WEIGHT;
    if (tag=="color")
        return _GT_ATTRIB_COLOR;
    if (tag=="style")
        return _GT_ATTRIB_STYLE;
    if (tag=="width")
        return _GT_ATTRIB_WIDTH;
    if (tag=="shape")
        return _GT_ATTRIB_SHAPE;
    if (tag=="pos")
        return _GT_ATTRIB_POSITION;
    if (tag=="name")
        return _GT_ATTRIB_NAME;
    if (tag=="temp")
        return _GT_ATTRIB_TEMPORARY;
    return register_user_tag(tag);
}

/* return attribute tag associated with the specified index */
string graphe::index2tag(int index) const {
    int len;
    switch (index) {
    case _GT_ATTRIB_LABEL:
        return "label";
    case _GT_ATTRIB_DIRECTED:
        return "directed";
    case _GT_ATTRIB_WEIGHTED:
        return "weighted";
    case _GT_ATTRIB_COLOR:
        return "color";
    case _GT_ATTRIB_STYLE:
        return "style";
    case _GT_ATTRIB_WIDTH:
        return "width";
    case _GT_ATTRIB_SHAPE:
        return "shape";
    case _GT_ATTRIB_WEIGHT:
        return "weight";
    case _GT_ATTRIB_POSITION:
        return "pos";
    case _GT_ATTRIB_NAME:
        return "name";
    case _GT_ATTRIB_TEMPORARY:
        return "temp";
    case _GT_ATTRIB_USER:
    default:
        len=index-_GT_ATTRIB_USER;
        assert(int(user_tags.size())>len);
        return user_tags[len];
    }
}

/* register custom user attribute tag */
int graphe::register_user_tag(const string &tag) {
    int i=_GT_ATTRIB_USER;
    vector<string>::const_iterator it=user_tags.begin();
    for (;it!=user_tags.end();++it) {
        if (*it==tag)
            return i;
        ++i;
    }
    user_tags.push_back(tag);
    return i;
}

/* register custom user tags */
void graphe::register_user_tags(const vector<string> &tags) {
    for (vector<string>::const_iterator it=tags.begin();it!=tags.end();++it) {
        register_user_tag(*it);
    }
}

/* fill V with marked vertices */
void graphe::get_marked_nodes(vecteur &V) const {
    V.clear();
    ivector_iter it=marked_nodes.begin(),itend=marked_nodes.end();
    for (;it!=itend;++it) {
        V.push_back(node_label(*it));
    }
}

/* fill m with marked nodes belonging to the s-th subgraph */
void graphe::get_marked_nodes_in_subgraph(int s,ivector &m) const {
    m.clear();
    m.reserve(marked_nodes.size());
    ivector_iter it=marked_nodes.begin(),itend=marked_nodes.end();
    for (;it!=itend;++it) {
        if (node(*it).subgraph()==s)
            m.push_back(*it);
    }
}

/* mark vertex with index v */
void graphe::mark_node(int v) {
    if (find(marked_nodes.begin(),marked_nodes.end(),v)==marked_nodes.end())
        marked_nodes.push_back(v);
}

/* remove vertex v from the list of marked vertices */
bool graphe::unmark_node(int v) {
    ivector::iterator it;
    if ((it=find(marked_nodes.begin(),marked_nodes.end(),v))!=marked_nodes.end()) {
        marked_nodes.erase(it);
        return true;
    }
    return false;
}

/* set edge {i,j} visited */
void graphe::set_edge_visited(int i,int j) {
    visited_edges[i<j?i:j].insert(i<j?j:i);
}

/* return true iff edge {i,j} is visited */
bool graphe::is_edge_visited(int i,int j) const {
    std::map<int,iset>::const_iterator ngh=visited_edges.find(i<j?i:j);
    if (ngh==visited_edges.end())
        return false;
    return ngh->second.count(i<j?j:i)!=0;
}

/* return edge as pair of vertex indices */
graphe::ipair graphe::make_edge(const vecteur &v) const {
    assert(v.size()==2);
    int i=node_index(v.front()),j=node_index(v.back());
    return make_pair(i,j);
}

/* convert vecteur E to list of ipairs representing edges, return false iff error occurs,
   if edge is not found set notfound=true */
bool graphe::edges2ipairs(const vecteur &E,ipairs &ev,bool &notfound) const {
    if (E.empty())
        return false;
    int i,j;
    if (ckmatrix(E)) {
        if (E.front()._VECTptr->size()!=2)
            return false;
        for (const_iterateur it=E.begin();it!=E.end();++it) {
            i=node_index(it->_VECTptr->front());
            j=node_index(it->_VECTptr->back());
            if (i<0 || j<0 || !has_edge(i,j)) {
                notfound=true;
                return false;
            }
            ev.push_back(make_pair(i,j));
        }
    } else {
        if (E.size()!=2)
            return false;
        i=node_index(E.front());
        j=node_index(E.back());
        if (i<0 || j<0 || !has_edge(i,j)) {
            notfound=true;
            return false;
        }
        ev.push_back(make_pair(i,j));
    }
    return true;
}

/* convert ipairs to edges */
vecteur graphe::ipairs2edges(const ipairs &E) const {
    vecteur res(E.size());
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        res[it-E.begin()]=makevecteur(node_label(it->first),node_label(it->second));
    }
    return res;
}

/* convert ipairs to ipairs */
void graphe::ipairs2edgeset(const ipairs &E,edgeset &Eset) {
    Eset.clear();
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        Eset.insert(*it);
    }
}

/* return total number of edges/arcs */
int graphe::edge_count(int sg) const {
    int count=0,dup=0;
    bool isdir=is_directed();
    node_iter it=nodes.begin(),itend=nodes.end();
    for(;it!=itend;++it) {
        if (sg>=0 && it->subgraph()!=sg)
            continue;
        if (sg<0)
            count+=it->degree();
        else {
            for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
                if (node(*jt).subgraph()==sg) count++;
            }
        }
        if (!isdir)
            dup+=it->multiedge_count();
    }
    if (!isdir) {
        count/=2;
        count+=dup;
    }
    return count;
}

/* return number of arcs going in the vertex with the specified index (when directed) */
int graphe::in_degree(int index,int sg) const {
    assert(index>=0 && index<node_count());
    int count=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if ((sg<0 || it->subgraph()==sg) && it->has_neighbor(index))
            count++;
    }
    return count;
}

/* return number of arcs going out of vertex with the specified index (when directed) */
int graphe::out_degree(int index,int sg) const {
    assert(index>=0 && index<node_count());
    const vertex &v=node(index);
    if (sg<0)
        return v.degree();
    int cnt=0;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        if (node(*it).subgraph()==sg)
            ++cnt;
    }
    return cnt;
}

/* return degree of vertex with the specified index */
int graphe::degree(int index,int sg) const {
    if (is_directed())
        return in_degree(index,sg)+out_degree(index,sg);
    return out_degree(index,sg);
}

/* compute in and out degrees and store them */
void graphe::compute_in_out_degrees(ivector &ind,ivector &outd) const {
    int n=node_count();
    ind.resize(n,0);
    outd.resize(n);
    for (int i=0;i<n;++i) {
        const vertex &v=node(i);
        outd[i]=v.neighbors().size();
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            ++ind[*it];
        }
    }
}

/* return the maximum vertex degree */
int graphe::maximum_degree(int sg) const {
    int maxdeg=0,d;
    for (int i=0;i<node_count();++i) {
        if ((d=degree(i,sg))>maxdeg)
            maxdeg=d;
    }
    return maxdeg;
}

/* return the minimum vertex degree */
int graphe::minimum_degree(int sg) const {
    int mindeg=rand_max2,d;
    for (int i=0;i<node_count();++i) {
        if ((d=degree(i,sg))<mindeg)
            mindeg=d;
    }
    return mindeg;
}

/* return the degree sequence of this graph */
graphe::ivector graphe::degree_sequence(int sg) const {
    ivector res;
    res.reserve(node_count());
    node_iter it=nodes.begin(),itend=nodes.end();
    int i=0;
    for (;it!=itend;++it,++i) {
        if (sg<0 || it->subgraph()==sg)
            res.push_back(degree(i));
    }
    return res;
}

/* sort the vertices by order of their degrees, starting from the highest */
void graphe::sort_by_degrees(graphe &G,ivector &sigma,int sg) {
    int n=node_count(),i;
    ipairs lst;
    lst.reserve(n);
    node_iter nt=nodes.begin();
    for (i=0;i<n;++i,++nt) {
        if (sg<0 || nt->subgraph()==sg)
            lst.push_back(make_pair(degree(i),i));
    }
    std::sort(lst.rbegin(),lst.rend());
    sigma.clear();
    sigma.reserve(n);
    ipairs_iter it=lst.begin(),itend=lst.end();
    for (;it!=itend;++it)
        sigma.push_back(it->second);
    isomorphic_copy(G,sigma,true,sg);
}

/* create the adjacency matrix of this graph */
void graphe::adjacency_matrix(matrice &m) const {
    int n=node_count();
    m=*_matrix(makesequence(n,n,0),context0)._VECTptr;
    node_iter it=nodes.begin(),itend=nodes.end();
    const_iterateur jt=m.begin();
    for (;it!=itend;++it,++jt) {
        vecteur &row=*jt->_VECTptr;
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            row[*jt]=1;
        }
    }
}

/* create the adjacency matrix as a sparse matrix */
void graphe::adjacency_sparse_matrix(sparsemat &sm,bool diag_ones,int sg) const {
    sm.clear();
    int i,j;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (i=0;it!=itend;++it,++i) {
        if (sg>=0 && it->subgraph()!=sg)
            continue;
        if (diag_ones)
            sm[i][i]=make_pair(1,1);
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            j=*jt;
            if (sg<0 || node(j).subgraph()==sg)
                sm[i][j]=make_pair(1,1);
        }
    }
}

/* create the Laplacian matrix of this graph */
void graphe::laplacian_matrix(matrice &m,bool normalize) const {
    int n=node_count();
    bool isweighted=is_weighted();
    if (isweighted)
        weight_matrix(m);
    else adjacency_matrix(m);
    vecteur D(n);
    for (int i=0;i<n;++i) {
        if (!isweighted)
            D[i]=degree(i);
        else {
            edgeset E;
            ivector V(1,i);
            incident_edges(V,E);
            gen w(0);
            for (edgeset_iter it=E.begin();it!=E.end();++it) {
                w+=weight(*it);
            }
            D[i]=w;
        }
    }
    if (normalize) {
        for (int i=0;i<n;++i) {
            for (int j=0;j<n;++j) {
                m[i]._VECTptr->at(j)=gen(i==j?1:0)-m[i][j]/giac::sqrt(D[i]*D[j],ctx);
            }
        }
    } else m=subvecteur(*_diag(D,ctx)._VECTptr,m);
}

/* return incidence matrix of this graph */
void graphe::incidence_matrix(matrice &m) const {
    ipairs E;
    get_edges_as_pairs(E);
    int nr=node_count(),nc=E.size();
    m=*_matrix(makesequence(nr,nc,0),context0)._VECTptr;
    bool isdir=is_directed();
    for (int i=0;i<nr;++i) {
        for (int j=0;j<nc;++j) {
            gen &mij=m[i]._VECTptr->at(j);
            const ipair &e=E[j];
            if (isdir) {
                if (i==e.first)
                    mij=-1;
                else if (i==e.second)
                    mij=1;
            } else if (i==e.first || i==e.second)
                mij=1;
        }
    }
}

/* return weight associated with edge {i,j} or arc [i,j] */
gen graphe::weight(int i,int j) const {
    const attrib &attr=edge_attributes(i,j);
    attrib_iter it=attr.find(_GT_ATTRIB_WEIGHT);
    if (it==attr.end())
        return undef;
    return it->second;
}

/* return weight matrix */
void graphe::weight_matrix(matrice &W) const {
    assert(is_weighted());
    int n=node_count(),i,j;
    W=*_matrix(makesequence(n,n,0),context0)._VECTptr;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (i=0;it!=itend;++it,++i) {
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            j=*jt;
            W[i]._VECTptr->at(j)=weight(i,j);
        }
    }
}

/* return list of vertices (in the given subgraph) */
vecteur graphe::vertices(int sg) const {
    vecteur V;
    V.reserve(node_count());
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (sg<0 || it->subgraph()==sg)
            V.push_back(it->label());
    }
    return V;
}

/* make all vertices (in the given subgraph) unvisited */
void graphe::unvisit_all_nodes(int sg) {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        if (sg<0 || it->subgraph()==sg)
            it->set_visited(false);
    }
}

/* unset ancestors for all vertices (in the given subgraph) */
void graphe::unset_all_ancestors(int sg) {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        if (sg<0 || it->subgraph()==sg)
            it->unset_ancestor();
    }
}

/* turn the color of each vertex to white */
void graphe::uncolor_all_nodes(int base_color,int sg) {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        if (sg<0 || it->subgraph()==sg)
            it->set_color(base_color);
    }
}

/* fill the list E with edges (in the given subgraph) represented as pairs of vertex indices */
void graphe::get_edges_as_pairs(ipairs &E,int sg) const {
    int i,j;
    bool isdir=is_directed();
    E.clear();
    E.reserve(edge_count());
    node_iter it=nodes.begin(),itend=nodes.end();
    for (i=0;it!=itend;++it,++i) {
        if (sg>=0 && it->subgraph()!=sg)
            continue;
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            j=*jt;
            if ((isdir || j>i) && (sg<0 || node(j).subgraph()==sg))
                E.push_back(make_pair(i,j));
        }
    }
}

/* return list of edges/arcs (in the given subgraph) */
vecteur graphe::edges(bool include_weights,int sg) const {
    ipairs E;
    get_edges_as_pairs(E,sg);
    vecteur edge(2),ret(E.size());
    int i=0;
    bool isdir=is_directed();
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        edge[0]=node(it->first).label();
        edge[1]=node(it->second).label();
        if (!isdir)
            edge=*_sort(edge,ctx)._VECTptr;
        ret[i++]=include_weights?makevecteur(edge,weight(*it)):edge;
    }
    return ret;
}

/* return the sequence of edge multiplicities */
graphe::ivector graphe::edge_multiplicities() const {
    ipairs E;
    get_edges_as_pairs(E);
    int m=E.size();
    ivector res(m);
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        res[it-E.begin()]=multiedges(*it)+1;
    }
    return res;
}

/* return the total number of edges, counting multiple edges too */
int graphe::sum_of_edge_multiplicities() const {
    ivector mult=edge_multiplicities();
    int res=0;
    for (ivector_iter it=mult.begin();it!=mult.end();++it) {
        res+=*it;
    }
    return res;
}

/* write attributes to dot file */
void graphe::write_attrib(ofstream &dotfile,const attrib &attr,bool style) const {
    dotfile << "[";
    bool comma=false;
    for (attrib_iter it=attr.begin();it!=attr.end();++it) {
        if (!style && (it->first==_GT_ATTRIB_COLOR || it->first==_GT_ATTRIB_SHAPE ||
                       it->first==_GT_ATTRIB_STYLE || it->first==_GT_ATTRIB_WIDTH || it->first==_GT_ATTRIB_POSITION))
            continue;
        if (comma) dotfile << ",";
        const gen &val=it->second;
        dotfile << index2tag(it->first) << "=";
        if (it->first==_GT_ATTRIB_POSITION && it->second.is_symb_of_sommet(at_point))
            dotfile << "\"" << it->second._SYMBptr->feuille << "\"";
        else if (val.type==_STRNG || val.type==_IDNT || val.is_integer() || val.type==_DOUBLE_ || val.type==_FLOAT_ || val.type==_REAL)
            dotfile << val;
        else dotfile << "\"" << val << "\"";
        comma=true;
    }
    dotfile << "]";
}

/* export the drawing of this graph in latex format */
bool graphe::write_latex(const string &filename,const gen &drawing) const {
    ofstream texfile;
    texfile.open(filename.c_str());
    if (!texfile.is_open())
        return false;
    show_axes(0,ctx);
    string picture=genstring2str(_latex(drawing,ctx));
    texfile << "% this file was generated by " << giac_version() << "\n";
    texfile << "\\documentclass{standalone}\n" << "\\usepackage[T1]{fontenc}\n"
            << "\\usepackage[utf8]{inputenc}\n" << "\\usepackage{pstricks}\n"
            << "\\begin{document}\n" << picture << "\\end{document}\n";
    texfile.close();
    return true;
}

/* Export this graph to dot file.
 * If style==false, then the color, shape, style, and position attributes are not exported. */
bool graphe::write_dot(const string &filename,bool style) const {
    ofstream dotfile;
    dotfile.open(filename.c_str());
    if (!dotfile.is_open()) {
        message("Error: failed to open file for writing");
        return false;
    }
    dotfile << "# this file was generated by " << giac_version() << "\n";
    ivector u,v;
    string indent("  "),edgeop(is_directed()?" -> ":" -- ");
    dotfile << (is_directed()?"digraph ":"graph ");
    string graph_name=name();
    if (graph_name.empty())
        dotfile << "{\n";
    else
        dotfile << graph_name << " {\n";
    if (!attributes.empty()) {
        dotfile << indent << "graph ";
        write_attrib(dotfile,attributes);
        dotfile << ";\n";
    }
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (!it->attributes().empty()) {
            dotfile << indent << it->label() << " ";
            write_attrib(dotfile,it->attributes(),style);
            dotfile << ";\n";
        }
        u.clear();
        v.clear();
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            if (it->neighbor_attributes(*jt).empty())
                u.push_back(*jt);
            else
                v.push_back(*jt);
        }
        if (!u.empty()) {
            dotfile << indent << it->label() << edgeop << "{ ";
            for (ivector_iter kt=u.begin();kt!=u.end();++kt) {
                dotfile << node(*kt).label() << " ";
            }
            dotfile << "};\n";
        }
        for (ivector_iter kt=v.begin();kt!=v.end();++kt) {
            dotfile << indent << it->label() << edgeop << node(*kt).label() << " ";
            write_attrib(dotfile,it->neighbor_attributes(*kt),style);
            dotfile << ";\n";
        }
    }
    dotfile << "}\n";
    dotfile.close();
    return true;
}

bool graphe::write_lst(const string &filename) const {
    assert(!is_directed());
    ofstream file;
    file.open(filename.c_str());
    if (!file.is_open()) {
        message("Error: failed to open file for writing");
        return false;
    }
    ipairs E;
    get_edges_as_pairs(E);
    std::sort(E.begin(),E.end());
    ipairs_iter et=E.begin();
    stringstream ss;
    iset visited;
    int a,n=node_count();
    for (int i=0;i<n;++i) {
        ss.str("");
        ss << i+1 << ":";
        a=0;
        while (et!=E.end() && et->first==i) {
            visited.insert(et->second);
            ss << " " << et->second+1;
            ++et;
            ++a;
        }
        if (a>0 || visited.find(i)==visited.end())
            file << ss.str() << "\n";
    }
    file.close();
    return true;
}

static int dot_comment_type;
static int dot_subgraph_level;
static bool dot_reading_attributes;
static bool dot_reading_value;
static int dot_newline_count;
static int dot_token_type;

/* returns true iff last read token is an identifier, number or string */
bool dot_token_is_id() {
    return dot_token_type==_GT_DOT_TOKEN_TYPE_IDENTIFIER ||
            dot_token_type==_GT_DOT_TOKEN_TYPE_NUMBER ||
            dot_token_type==_GT_DOT_TOKEN_TYPE_STRING;
}

/* return true iff char c is one of dot format delimiters */
bool dot_is_id_delim(const char &c) {
    return c=='{' || c=='}' || c=='[' || c==']' || c=='-' ||
            c==',' || c==';' || c=='=' || c=='/' || isspace(c);
}

/* read next token in dot file 'dotfile' */
int dot_read_token(ifstream &dotfile,string &token) {
    char c=0,pc,nc;
    bool last=false;
    dot_token_type=0;
    token=string("");
    while (true) {
        if ((dot_token_type==_GT_DOT_TOKEN_TYPE_IDENTIFIER || dot_token_type==_GT_DOT_TOKEN_TYPE_NUMBER)
                && dot_is_id_delim(dotfile.peek()))
            break;
        pc=c;
        if (!dotfile.get(c))
            break;
        if (c=='\n')
            ++dot_newline_count;
        if (dot_token_type==_GT_DOT_TOKEN_TYPE_IDENTIFIER) {
            if (isalnum(c) || c=='_')
                token.push_back(c);
            else return 0;
            continue;
        }
        if (dot_token_type==_GT_DOT_TOKEN_TYPE_NUMBER) {
            if (isdigit(c) || c=='.')
                token.push_back(c);
            else return 0;
            continue;
        }
        nc=dotfile.peek();
        if (nc==EOF)
            last=true;
        if (dot_comment_type!=0) {
            if (((c=='\n' || c=='\r') && dot_comment_type==1) ||
                    ((c=='/' && pc=='*') && dot_comment_type==2)) {
                if (token.empty()) {
                    dot_comment_type=0;
                    continue;
                }
                break;
            }
            continue;
        }
        if (c=='"') {
            if (dot_token_type==_GT_DOT_TOKEN_TYPE_STRING) {
                if (pc=='\\')  {
                    token.erase(token.end()-1);
                    token.push_back(c);
                    continue;
                }
                break;
            }
            dot_token_type=_GT_DOT_TOKEN_TYPE_STRING;
            continue;
        }
        if (dot_token_type==_GT_DOT_TOKEN_TYPE_STRING) {
            if (last)
                return 0;
            token.push_back(c);
            continue;
        }
        if (isblank(c) || (dot_subgraph_level>0 && c==';') ||
                (dot_reading_attributes && c==',') || (!dot_reading_attributes && isspace(c))) {
            if (isspace(c))
                continue;
            token=string(1,c);
            dot_token_type=_GT_DOT_TOKEN_TYPE_DELIMITER;
            return 1;
        }
        if (c=='{' || c=='}') {
            dot_subgraph_level+=c=='{'?1:-1;
            token=string(1,c);
            dot_token_type=_GT_DOT_TOKEN_TYPE_DELIMITER;
            return 1;
        }
        if (c=='[' || c==']') {
            dot_reading_attributes=c=='[';
            token=string(1,c);
            dot_token_type=_GT_DOT_TOKEN_TYPE_DELIMITER;
            return 1;
        }
        if (c=='=') {
            dot_reading_value=true;
            token=string(1,c);
            dot_token_type=_GT_DOT_TOKEN_TYPE_OPERATOR;
            return 1;
        }
        if (c=='#' && (pc==0 || pc=='\n' || pc=='\r')) {
            dot_comment_type=1;
            continue;
        }
        if (c=='/') {
            if (last)
                return 0;
            if (nc=='/')
                dot_comment_type=1;
            if (nc=='*')
                dot_comment_type=2;
            continue;
        }
        if (token.empty() && (isalpha(c) || c=='_')) {
            dot_token_type=_GT_DOT_TOKEN_TYPE_IDENTIFIER;
            token.push_back(c);
            continue;
        }
        if (token.empty() && (isdigit(c) || (c=='-' && isdigit(nc)))) {
            dot_token_type=_GT_DOT_TOKEN_TYPE_NUMBER;
            token.push_back(c);
            continue;
        }
        if ((c=='-' || c=='>') && pc=='-') {
            dot_token_type=_GT_DOT_TOKEN_TYPE_OPERATOR;
            token=string("-")+c;
            return 1;
        }
    }
    if (dot_token_is_id())
        dot_reading_value=false;
    return dot_token_type==0?-1:1;
}

/* insert attribute key=val to attr [and overwrite existing value] */
bool graphe::insert_attribute(attrib &attr,int key,const gen &val,bool overwrite) {
    pair<attrib::iterator,bool> res;
    if ((res=attr.insert(make_pair(key,val))).second==false && overwrite)
        (res.first)->second=val;
    return res.second;
}

/* remove attribute with key 'key' from attr */
bool graphe::remove_attribute(attrib &attr,int key) {
    attrib::iterator it=attr.find(key);
    if (it==attr.end())
        return false;
    attr.erase(it);
    return true;
}

/* parse attributes from dot file */
bool graphe::dot_parse_attributes(ifstream &dotfile,attrib &attr) {
    string token;
    int key;
    while(true) {
        if (dot_read_token(dotfile,token)!=1)
            return false;
        if (token=="]")
            break;
        if (token==";" || token==",")
            continue;
        if (!dot_token_is_id())
            return false;
        key=tag2index(token);
        if (key==-1 || dot_read_token(dotfile,token)!=1 || token!="=" ||
                dot_read_token(dotfile,token)!=1 || dot_reading_value || !dot_token_is_id())
            return false;
        if (key==_GT_ATTRIB_WEIGHT && !is_weighted())
            set_weighted(true);
        insert_attribute(attr,key,str2gen(token,dot_token_type==_GT_DOT_TOKEN_TYPE_STRING));
    }
    return true;
}

/* initialize graph from dot file */
bool graphe::read_dot(const string &filename) {
    ifstream dotfile;
    dotfile.open(filename.c_str());
    if (!dotfile.is_open()) {
        message("Error: failed to open file for reading");
        return false;
    }
    dot_subgraph_level=0;
    dot_reading_attributes=false;
    dot_reading_value=false;
    dot_comment_type=0;
    dot_newline_count=0;
    map<int,attrib> delayed_attributes;
    string token;
    int res,gtype,subgraph_count=0,index;
    bool strict=false,error_raised=false,has_root_graph=false;
    vector<dotgraph> subgraphs;
    while (true) {
        res=dot_read_token(dotfile,token);
        if (res!=1) { error_raised=(bool)(res+1); break; }
        switch (dot_token_type) {
        case _GT_DOT_TOKEN_TYPE_IDENTIFIER:
            if (token=="strict") {
                strict=true;
                continue;
            } else if (token=="graph" || token=="digraph" || token=="subgraph") {
                if (token=="digraph") {
                    if (dot_subgraph_level>0) { error_raised=true; break; }
                    set_directed(true);
                    gtype=1;
                } else if (token=="subgraph") {
                    gtype=2;
                } else gtype=0;
                if (gtype!=2 && strict)
                    strict=false;
                if (dot_read_token(dotfile,token)!=1) { error_raised=true; break; }
                if (dot_token_type!=_GT_DOT_TOKEN_TYPE_DELIMITER) {
                    if (!dot_token_is_id()) { error_raised=true; break; }
                    if (gtype!=2)
                        set_name(token);
                    if (dot_read_token(dotfile,token)!=1 ||
                            dot_token_type!=_GT_DOT_TOKEN_TYPE_DELIMITER) { error_raised=true; break; }
                } else if (token=="[") {
                    if (gtype!=0 || !dot_parse_attributes(dotfile,attributes))
                        error_raised=true;
                    break;
                }
                if ((gtype!=2 && has_root_graph) || token!="{")
                    error_raised=true;
                subgraphs.push_back(dotgraph(subgraph_count++));
                if (gtype!=2)
                    has_root_graph=true;
                break;
            } else if (token=="node" || token=="edge") {
                if (dot_read_token(dotfile,token)!=1 || dot_token_type!=_GT_DOT_TOKEN_TYPE_DELIMITER ||
                        token!="[" || !dot_parse_attributes(dotfile,token=="node"?
                                                            subgraphs.back().vertex_attributes():
                                                            subgraphs.back().edge_attributes()))
                    error_raised=true;
                break;
            }
        case _GT_DOT_TOKEN_TYPE_NUMBER:
        case _GT_DOT_TOKEN_TYPE_STRING:
            index=add_node(str2gen(token,dot_token_type!=_GT_DOT_TOKEN_TYPE_NUMBER));
            if (subgraphs.empty()) { error_raised=true; break; }
            node(index).set_subgraph(subgraphs.back().index());
            subgraphs.back().set_index(index+1);
            break;
        case _GT_DOT_TOKEN_TYPE_DELIMITER:
            if (token==";") {
                if (subgraphs.back().chain_empty())
                    continue;
                if (!subgraphs.back().chain_completed()) { error_raised=true; break; }
                ivector &ch=subgraphs.back().chain();
                attrib &attr=subgraphs.back().chain_attributes();
                if (subgraphs.back().position()==0) {
                    if (ch.front()<=0) { error_raised=true; break; }
                    delayed_attributes[ch.front()-1]=attr;
                }
                int lh,rh;
                for (int i=0;i<=subgraphs.back().position()-1;++i) {
                    lh=ch[i];
                    rh=ch[i+1];
                    if (lh==0 || rh==0) { error_raised=true; break; }
                    node_iter it=nodes.begin(),jt,itend=nodes.end();
                    if (lh>0 && rh>0)
                        add_edge(lh-1,rh-1,attr);
                    else if (lh<0 && rh<0) {
                        int j=0,k=0;
                        for (;it!=itend;++it) {
                            if (it->subgraph()==-lh) {
                                for (jt=nodes.begin();jt!=itend;++jt) {
                                    if (jt->subgraph()==-rh)
                                        add_edge(j,k,attr);
                                    ++k;
                                }
                            }
                            ++j;
                        }
                    } else {
                        int j=0;
                        for (;it!=itend;++it) {
                            if (it->subgraph()==(lh>0?-rh:-lh))
                                add_edge(lh>0?lh-1:j,lh>0?j:rh-1,attr);
                            ++j;
                        }
                    }
                }
                subgraphs.back().clear_chain();
            } else if (token=="[") {
                attrib &attr=subgraphs.back().chain_attributes();
                if (!attr.empty() || !dot_parse_attributes(dotfile,attr))
                    error_raised=true;
            } else if (token=="{") {
                if (dot_subgraph_level<2) { error_raised=true; break; }
                subgraphs.push_back(dotgraph(subgraph_count++));
            } else if (token=="}") {
                if (subgraphs.empty()) { error_raised=true; break; }
                index=subgraphs.back().index();
                subgraphs.pop_back();
                if (!subgraphs.empty())
                    subgraphs.back().set_index(-index);
            }
            break;
        case _GT_DOT_TOKEN_TYPE_OPERATOR:
            if (dot_subgraph_level<1) { error_raised=true; break; }
            if (token=="--" || token=="->") {
                if ((is_directed() && token=="--") || (!is_directed() && token=="->")) { error_raised=true; break; }
                subgraphs.back().incr();
            }
            break;
        default:
            message((string("Error: unknown token '")+token+"'").c_str());
            break;
        }
        if (error_raised || strict)
            break;
    }
    dotfile.close();
    if (dot_subgraph_level!=0 || strict || dot_reading_value || dot_reading_attributes)
        error_raised=true;
    else {
        /* set delayed attributes to nodes */
        for (map<int,attrib>::const_iterator it=delayed_attributes.begin();it!=delayed_attributes.end();++it) {
            vertex &v=node(it->first);
            for (attrib_iter ait=it->second.begin();ait!=it->second.end();++ait) {
                v.set_attribute(ait->first,ait->second);
            }
        }
    }
    return !error_raised;
}

bool graphe::read_lst(const string &filename) {
    ifstream file;
    string word;
    file.open(filename.c_str());
    if (!file.is_open()) {
        message("Error: failed to open file for reading");
        return false;
    }
    char *e;
    map<int,ivector> adj;
    iset vert;
    int v=-1;
    bool failed=false;
    while (file >> word) {
        if (word.empty())
            continue;
        const char *token=word.c_str();
        int k=strtol(token,&e,10);
        if (k<=0) {
            failed=true;
            break;
        }
        vert.insert(k);
        if (*e==':') {
            v=k-1;
        } else {
            if (v<0) {
                failed=true;
                break;
            }
            adj[v].push_back(k-1);
        }
    }
    if (!file.eof()) {
        message("Error: failed to read to the end of the file");
        file.close();
        return false;
    }
    file.close();
    if (failed) {
        message("Error: not a LST file");
        return false;
    }
    int n=*vert.rbegin();
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,n);
        add_nodes(labels);
    } else add_nodes(n);
    for (map<int,ivector>::const_iterator it=adj.begin();it!=adj.end();++it) {
        int i=it->first;
        const ivector &a=it->second;
        for (ivector_iter jt=a.begin();jt!=a.end();++jt) {
            add_edge(i,*jt);
        }
    }
    return true;
}

/* assign weights from matrix m to edges/arcs of this graph */
void graphe::make_weighted(const matrice &m) {
    assert(is_squarematrix(m) && int(m.size())==node_count());
    int i=0;
    bool isdir=is_directed();
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it,++i) {
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            if (isdir || i<*jt)
                insert_attribute(it->neighbor_attributes(*jt),_GT_ATTRIB_WEIGHT,m[i][*jt]);
        }
    }
    set_graph_attribute(_GT_ATTRIB_WEIGHTED,VRAI);
}

/* make this graph unweighted */
void graphe::make_unweighted() {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            remove_attribute(it->neighbor_attributes(*jt),_GT_ATTRIB_WEIGHT);
        }
    }
    set_graph_attribute(_GT_ATTRIB_WEIGHTED,FAUX);
}

/* randomize edge weights, generating them in segment [a,b] */
void graphe::randomize_edge_weights(double a,double b,bool integral_weights) {
    assert(a<=b);
    if (!is_weighted())
        set_graph_attribute(_GT_ATTRIB_WEIGHTED,boole(true));
    int m,n;
    if (integral_weights) {
        m=std::floor(a);
        n=std::floor(b);
    }
    int nc=node_count();
    gen w;
    for (int i=0;i<nc;++i) {
        for (int j=is_directed()?0:i+1;j<nc;++j) {
            if (!has_edge(i,j))
                continue;
            if (integral_weights)
                w=gen(rand_integer(n-m+1)+m);
            else
                w=gen(a+rand_uniform()*(b-a));
            set_edge_attribute(i,j,_GT_ATTRIB_WEIGHT,w);
        }
    }
}

/* store the underlying graph to G (convert arcs to edges and strip all attributes) */
void graphe::underlying(graphe &G,int sg) const {
    assert(supports_attributes() || !G.supports_attributes());
    G.set_directed(false);
    map<int,int> vmap;
    copy_nodes(G,vmap,sg);
    int i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if (sg>=0 && it->subgraph()!=sg)
            continue;
        if (sg<0)
            G.node(sg<0?i:vmap[i]).set_subgraph(it->subgraph());
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            if (sg<0 || node(*jt).subgraph()==sg)
                G.add_edge(sg<0?i:vmap[i],sg<0?*jt:vmap[*jt]);
        }
    }
}

/* store the complement of this graph in G */
void graphe::complement(graphe &G,int sg) const {
    assert(supports_attributes() || !G.supports_attributes());
    map<int,int> vmap;
    copy_nodes(G,vmap,sg);
    bool isdir=is_directed();
    G.set_directed(isdir);
    int n=node_count();
    for (int i=0;i<n;++i) {
        if (sg>=0 && node(i).subgraph()!=sg)
            continue;
        for (int j=G.is_directed()?0:i+1;j<n;++j) {
            if ((sg<0 || node(j).subgraph()==sg) && !has_edge(i,j))
                G.add_edge(sg<0?i:vmap[i],sg<0?j:vmap[j]);
        }
    }
}

/* convert attrib to gen_map */
void graphe::attrib2genmap(const attrib &attr, gen_map &m,bool keys2tags) const {
    m.clear();
    for (attrib_iter it=attr.begin();it!=attr.end();++it) {
        m[keys2tags?string2gen(index2tag(it->first),false):gen(it->first)]=it->second;
    }
}

bool gmap_find(const gen_map &gmap,const gen &key,gen &val) {
    gen_map::const_iterator it=gmap.find(key);
    if (it==gmap.end())
        return false;
    val=it->second;
    return true;
}

/* read special graph from a list of adjacency lists of integers */
void graphe::read_special(const int *special_graph) {
    int state=1;
    int v=-1;
    map<int,ivector> adj;
    iset vert;
    for (const int *p=special_graph;*p!=-1;++p) {
        if (*p==0) {
            state=1;
        } else if (state==1) {
            v=*p-1;
            vert.insert(*p);
            state=2;
        } else if (state==2) {
            assert(v>=0);
            vert.insert(*p);
            adj[v].push_back(*p-1);
        }
    }
    int n=*vert.rbegin();
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,n);
        add_nodes(labels);
    } else add_nodes(n);
    for (map<int,ivector>::const_iterator it=adj.begin();it!=adj.end();++it) {
        int i=it->first;
        const ivector &a=it->second;
        for (ivector_iter jt=a.begin();jt!=a.end();++jt) {
            add_edge(i,*jt);
        }
    }
}

/* read special graph from a list of adjacency lists of strings */
void graphe::read_special(const char **special_graph) {
    const char **p=special_graph;
    int state=1;
    gen v,w;
    string s;
    do {
        s=string(*p);
        if (s.empty()) {
            state=1;
        } else if (state==1) {
            state=2;
            v=str2gen(s,true);
        } else if (state==2) {
            w=str2gen(s,true);
            add_edge(v,w);
        }
        ++p;
    } while (*p!=NULL);
}

graphe &graphe::operator =(const graphe &other) {
    nodes.clear();
    m_supports_attributes=other.supports_attributes();
    other.copy(*this);
    return *this;
}

/* make a copy of this graph and store it in G */
void graphe::copy(graphe &G) const {
    assert(supports_attributes() || !G.supports_attributes());
    G.clear();
    if (G.supports_attributes()) {
        G.set_name(name());
        G.register_user_tags(user_tags);
        G.set_graph_attributes(attributes);
    }
    G.copy_nodes(nodes);
    G.copy_marked_nodes(get_marked_nodes());
}

void graphe::copy_nodes(const vector<vertex> &V) {
    nodes=V;
    if (!supports_attributes()) {
        for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
            it->unsupport_attributes();
        }
    }
}

/* copies the nodes from sg into G, without attributes except node labels */
void graphe::copy_nodes(graphe &G,map<int,int> &vmap,int sg) const {
    int n=subgraph_size(sg);
    vmap.clear();
    G.clear();
    G.reserve_nodes(n);
    ivector V;
    if (sg>=0) {
        get_subgraph(sg,V);
        for (ivector_iter it=V.begin();it!=V.end();++it)
            vmap[*it]=it-V.begin();
    }
    if (G.supports_attributes()) {
        G.add_nodes(sg<0?vertices():get_node_labels(V));
    }
    else G.add_nodes(n);
}

/* returns true iff graph has edge {i,j} */
bool graphe::has_edge(int i,int j) const {
    if (i<0 || i>=node_count() || j<0 || j>=node_count())
        return false;
    return node(i).has_neighbor(j);
}

/* returns true iff i-th and j-th vertices are adjacent */
bool graphe::nodes_are_adjacent(int i,int j) const {
    return node(i).has_neighbor(j) || node(j).has_neighbor(i);
}

/* returns the attributes of the i-th node */
const graphe::attrib &graphe::node_attributes(int i) const {
    assert(i>=0 && i<node_count() && supports_attributes());
    return node(i).attributes();
}

graphe::attrib &graphe::node_attributes(int i) {
    assert(i>=0 && i<node_count() && supports_attributes());
    return node(i).attributes();
}

/* return const reference to the attributes assigned to edge [i,j] */
const graphe::attrib &graphe::edge_attributes(int i,int j) const {
    assert(supports_attributes());
    if (is_directed())
        return node(i).neighbor_attributes(j);
    return node(i<j?i:j).neighbor_attributes(i<j?j:i);
}

/* return the modifiable reference to the attributes assigned to edge [i,j] */
graphe::attrib &graphe::edge_attributes(int i,int j) {
    assert(supports_attributes());
    if (is_directed())
        return node(i).neighbor_attributes(j);
    return node(i<j?i:j).neighbor_attributes(i<j?j:i);
}

/* convert attrib to pair of vecteurs */
void graphe::attrib2vecteurs(const attrib &attr,vecteur &tags,vecteur &values) const {
    for (attrib_iter it=attr.begin();it!=attr.end();++it) {
        tags.push_back(str2gen(index2tag(it->first),true));
        values.push_back(it->second);
    }
}

/* add edge {i,j} or arc [i,j], depending on the type (undirected or directed) */
bool graphe::add_edge(int i,int j,const gen &w) {
    assert(i>=0 && i<node_count() && j>=0 && j<node_count());
    if (has_edge(i,j))
        return false;
    node(i).add_neighbor(j);
    if (!is_directed())
        node(j).add_neighbor(i);
    if (is_weighted()) {
        assert(supports_attributes());
        set_edge_attribute(i,j,_GT_ATTRIB_WEIGHT,w);
    }
    return true;
}

/* add edge {i,j} or arc [i,j] with attributes */
bool graphe::add_edge(int i,int j,const attrib &attr) {
    assert(i>=0 && i<node_count() && j>=0 && j<node_count() && supports_attributes());
    if (has_edge(i,j))
        return false;
    if (is_directed())
        node(i).add_neighbor(j,attr);
    else {
        int v=i<j?i:j,w=i<j?j:i;
        node(v).add_neighbor(w,attr);
        node(w).add_neighbor(v);
    }
    return true;
}

/* add edge {v,w} or arc [v,w], adding vertices v and/or w if necessary */
bool graphe::add_edge(const gen &v,const gen &w,const gen &weight) {
    assert(supports_attributes());
    int i=add_node(v),j=add_node(w);
    return add_edge(i,j,weight);
}

/* add edge {v,w} or arc [v,w], adding vertices v and/or w if necessary */
bool graphe::add_edge(const gen &v,const gen &w,const attrib &attr) {
    assert(supports_attributes());
    int i=add_node(v),j=add_node(w);
    return add_edge(i,j,attr);
}

/* add temporary edge from i-th to j-th vertex */
void graphe::add_temporary_edge(int i,int j) {
    assert(!has_edge(i,j));
    add_edge(i,j);
    if (supports_attributes())
        set_edge_attribute(i,j,_GT_ATTRIB_TEMPORARY,boole(true));
}

/* return true iff the edge {i,j} is temporary */
bool graphe::is_temporary_edge(int i,int j) const {
    assert(supports_attributes());
    assert(has_edge(i,j));
    gen val;
    get_edge_attribute(i,j,_GT_ATTRIB_TEMPORARY,val);
    return giac::is_one(val);
}

/* remove all temporary edges */
void graphe::remove_temporary_edges() {
    assert(supports_attributes());
    stack<ipair> edges;
    int i=0;
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it,++i) {
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            if (is_temporary_edge(i,*jt))
                edges.push(make_pair(i,*jt));
        }
    }
    while (!edges.empty()) {
        remove_edge(edges.top());
        edges.pop();
    }
}

/* remove edge {v,w} or arc [v,w] */
bool graphe::remove_edge(int i,int j) {
    if (!has_edge(i,j))
        return false;
    node(i).remove_neighbor(j);
    if (!is_directed())
        node(j).remove_neighbor(i);
    return true;
}

/* add new vertex to the graph */
int graphe::add_node() {
    assert(!supports_attributes());
    nodes.push_back(vertex(false));
    return node_count()-1;
}

/* add vertex v to the graph */
int graphe::add_node(const gen &v,const attrib &attr) {
    assert(supports_attributes());
    if (v.type==_MAP) {
        attrib a;
        gen lab(undef);
        int i;
        for (gen_map::const_iterator it=v._MAPptr->begin();it!=v._MAPptr->end();++it) {
            if (it->first.type!=_STRNG)
                continue;
            i=tag2index(*it->first._STRNGptr);
            if (i==_GT_ATTRIB_LABEL)
                lab=it->second;
            else a[i]=it->second;
        }
        if (is_undef(lab))
            return -1; // error
        return add_node(lab,a);
    }
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->label()==v)
            return it-nodes.begin();
    }
    nodes.push_back(vertex(v,attr));
    return node_count()-1;
}

/* add vertices from list v to the graph */
void graphe::add_nodes(const vecteur &v) {
    assert(supports_attributes());
    const_iterateur it=v.begin(),itend=v.end();
    for (;it!=itend;++it)
        add_node(*it);
}

/* adds n new nodes to the graph */
void graphe::add_nodes(int n) {
    assert(!supports_attributes());
    for (int i=0;i<n;++i)
        nodes.push_back(vertex(false));
}

/* remove all nodes with indices from V */
void graphe::isolate_nodes(const iset &V) {
    ivector adj;
    bool isdir=is_directed();
    for (iset_iter it=V.begin();it!=V.end();++it) {
        adjacent_nodes(*it,adj);
        for (ivector_iter jt=adj.begin();jt!=adj.end();++jt) {
            remove_edge(*it,*jt);
            if (isdir)
                remove_edge(*jt,*it);
        }
    }
}

/* convert list of labels to set of vertex indices */
bool graphe::labels2iset(const vecteur &labels,iset &s) {
    int i;
    for (const_iterateur it=labels.begin();it!=labels.end();++it) {
        if ((i=node_index(*it))==-1)
            return false;
        s.insert(i);
    }
    return true;
}

/* return vector of node labels */
vecteur graphe::get_node_labels(const ivector &v) const {
    assert(supports_attributes());
    vecteur V;
    V.resize(v.size());
    for (ivector_iter it=v.begin();it!=v.end();++it) {
        V[it-v.begin()]=node_label(*it);
    }
    return V;
}

/* return index of vertex v */
int graphe::node_index(const gen &v) const {
    assert(supports_attributes());
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->label()==v)
            return it-nodes.begin();
    }
    return -1;
}

/* return index of edge (i,j) as returned by get_edges_as_pairs */
int graphe::edge_index(const ipair &e) const {
    int i,j,cnt=0;
    bool isdir=is_directed();
    ipair edge=e;
    if (!isdir && e.first>e.second) {
        int tmp=edge.first;
        edge.first=edge.second;
        edge.second=tmp;
    }
    node_iter it=nodes.begin(),itend=nodes.end();
    for (i=0;it!=itend;++it,++i) {
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            j=*jt;
            if (isdir || j>i) {
                if (i==edge.first && j==edge.second)
                    return cnt;
                ++cnt;
            }
        }
    }
    return -1;
}

/* set 'subgraph' field of all nodes in v to s */
void graphe::set_subgraph(const ivector &v,int s) {
    for (ivector_iter it=v.begin();it!=v.end();++it) {
        node(*it).set_subgraph(s);
    }
}

/* set 'subgraph' field for nodes in subgraph induced by e */
void graphe::set_subgraph(const ipairs &e,int s) {
    for (ipairs_iter it=e.begin();it!=e.end();++it) {
        vertex &v=node(it->first),&w=node(it->second);
        v.set_subgraph(s);
        w.set_subgraph(s);
    }
}

/* collect all vertices within the given subgraph in list v */
void graphe::get_subgraph(int sg,ivector &v) const {
    v.clear();
    node_iter it=nodes.begin(),itend=nodes.end();
    int i=0;
    for (;it!=itend;++it,++i) {
        if (it->subgraph()==sg)
            v.push_back(i);
    }
}

/* change subgraph index */
void graphe::renumber_subgraph(int sg,int sg_new) {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        if (it->subgraph()==sg)
            it->set_subgraph(sg_new);
    }
}

/* return the index of the first vertex with subgraph field equal to sg,
 * or -1 if no such vertices exist */
int graphe::first_vertex_from_subgraph(int sg) const {
    if (sg<0)
        return 0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->subgraph()==sg)
            return it-nodes.begin();
    }
    return -1;
}

/* return the size of subgraph sg */
int graphe::subgraph_size(int sg) const {
    if (sg<0)
        return node_count();
    int cnt=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->subgraph()==sg)
            ++cnt;
    }
    return cnt;
}

/* merge subgraphs s and t, such that the resulting subgraph has index s */
void graphe::merge_subgraphs(int s,int t) {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        if (it->subgraph()==t)
            it->set_subgraph(s);
    }
}

/* set subgraph field to default_sg (by default -1) for every vertex in the graph */
void graphe::unset_subgraphs(int default_sg) {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->set_subgraph(default_sg);
    }
}

/* return maximal subgraph index in this graph */
int graphe::max_subgraph_index() const {
    int mi=-2;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->subgraph()>mi)
            mi=it->subgraph();
    }
    return mi;
}

/* make cycle graph with n vertices */
void graphe::make_cycle_graph() {
    int n=node_count();
    for (int i=0;i<n;++i) {
        add_edge(i,(i+1)%n);
    }
}

/* set vertex attribute key=val */
void graphe::set_node_attribute(int index,int key,const gen &val) {
    assert (index>=0 && index<node_count());
    node(index).set_attribute(key,val);
}

/* set edge {i,j} attribute key=val */
void graphe::set_edge_attribute(int i,int j,int key,const gen &val) {
    attrib &attr=edge_attributes(i,j);
    attr[key]=val;
}

/* return the value corresponding to the graph attribute specified by key */
bool graphe::get_graph_attribute(int key,gen &val) const {
    attrib_iter it=attributes.find(key);
    if (it==attributes.end()) {
        val=undef;
        return false;
    }
    val=it->second;
    return true;
}

/* return the value corresponding to the i-th node attribute specified by key */
bool graphe::get_node_attribute(int i,int key,gen &val) const {
    assert (i>=0 && i<node_count());
    const attrib &attr=node(i).attributes();
    attrib_iter it=attr.find(key);
    if (it==attr.end()) {
        val=undef;
        return false;
    }
    val=it->second;
    return true;
}

/* return the value corresponding to the attribute assigned to edge [i,j] and specified by key */
bool graphe::get_edge_attribute(int i,int j,int key,gen &val) const {
    const attrib &attr=edge_attributes(i,j);
    attrib_iter it=attr.find(key);
    if (it==attr.end()) {
        val=undef;
        return false;
    }
    val=it->second;
    return true;
}

/* discard the graph attribute specified by key */
void graphe::discard_graph_attribute(int key) {
    attrib::iterator it=attributes.find(key);
    if (it!=attributes.end())
        attributes.erase(it);
}

/* discard the attribute assigned to i-th node and specified by key */
void graphe::discard_node_attribute(int i,int key) {
    attrib &attr=node(i).attributes();
    attrib::iterator it=attr.find(key);
    if (it!=attr.end())
        attr.erase(it);
}

/* discard the attribute assigned to edge [i,j] and specified by key */
void graphe::discard_edge_attribute(int i,int j,int key) {
    attrib &attr=edge_attributes(i,j);
    attrib::iterator it=attr.find(key);
    if (it!=attr.end())
        attr.erase(it);
}

/* return true iff the graph is directed */
bool graphe::is_directed() const {
    attrib_iter it=attributes.find(_GT_ATTRIB_DIRECTED);
    assert(it!=attributes.end() && it->second.is_integer());
    return (bool)it->second.val;
}

/* return true iff the graph is weighted */
bool graphe::is_weighted() const {
    attrib_iter it=attributes.find(_GT_ATTRIB_WEIGHTED);
    assert(it!=attributes.end() && it->second.is_integer());
    return (bool)it->second.val;
}

/* create the subgraph defined by vertices from 'vi' and store it in G */
void graphe::induce_subgraph(const ivector &vi,graphe &G) const {
    assert(supports_attributes() || !G.supports_attributes());
    int n=node_count();
    G.clear();
    G.reserve_nodes(vi.size());
    if (!G.supports_attributes())
        G.add_nodes(vi.size());
    ivector sg_pos(n,-1);
    bool isdir=is_directed();
    G.set_directed(isdir);
    G.set_weighted(is_weighted());
    for (ivector_iter it=vi.begin();it!=vi.end();++it) {
        sg_pos[*it]=it-vi.begin();
        if (G.supports_attributes()) {
            gen v_label=node_label(*it);
            const attrib &attri=node(*it).attributes();
            G.add_node(v_label,attri);
        }
    }
    int i,j,k;
    for (ivector_iter it=vi.begin();it!=vi.end();++it) {
        const vertex &v=node(*it);
        i=it-vi.begin();
        for (ivector_iter jt=v.neighbors().begin();jt!=v.neighbors().end();++jt) {
            if ((j=sg_pos[*jt])>=0 && (isdir || i<j)) {
                ipair e=make_pair(i,j);
                ipair Ge=make_pair(*it,*jt);
                if (G.supports_attributes())
                    G.add_edge(e,edge_attributes(Ge));
                else G.add_edge(e);
                if (!isdir && (k=multiedges(Ge))>0) {
                    G.set_multiedge(e,k);
                }
            }
        }
    }
}

/* create the subgraph G defined by a list of edges E */
void graphe::extract_subgraph(const ipairs &E,graphe &G) const {
    assert(supports_attributes() || !G.supports_attributes());
    G.clear();
    bool isdir=is_directed();
    G.set_directed(isdir);
    G.set_weighted(is_weighted());
    iset vset;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        vset.insert(it->first);
        vset.insert(it->second);
    }
    G.reserve_nodes(vset.size());
    map<int,int> index_map;
    for (iset_iter it=vset.begin();it!=vset.end();++it) {
        const vertex &v=node(*it);
        index_map[*it]=G.supports_attributes()?G.add_node(v.label(),v.attributes()):G.add_node();
    }
    ipair e;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        e=make_pair(index_map[it->first],index_map[it->second]);
        if (G.supports_attributes())
            G.add_edge(e,edge_attributes(*it));
        else G.add_edge(e);
        if (!isdir)
            G.set_multiedge(e,multiedges(*it));
    }
}

/* find indices in subgraph S of vertices in v in this graph and store them in sv
 * (in the given order) */
void graphe::subgraph_indices(const graphe &S,const ivector &v,ivector &sv) const {
    assert(supports_attributes() && S.supports_attributes());
    sv.resize(v.size());
    for (ivector_iter it=v.begin();it!=v.end();++it) {
        assert((sv[it-v.begin()]=S.node_index(node_label(*it)))>=0);
    }
}

/* return true iff this graph is subgraph of G */
bool graphe::is_subgraph(const graphe &G) const {
    assert(supports_attributes() && G.supports_attributes());
    if (is_directed()!=G.is_directed() ||
            node_count()>G.node_count() ||
            edge_count()>G.edge_count())
        return false;
    int i,j;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        i=G.node_index(it->label());
        if (i<0)
            return false;
        for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
            j=G.node_index(node(*jt).label());
            if (j<0 || !G.has_edge(i,j))
                return false;
        }
    }
    return true;
}

/* fill the list adj with vertices adjacent to the i-th one */
void graphe::adjacent_nodes(int i,ivector &adj,bool include_temp_edges) const {
    assert(i>=0 && i<node_count());
    const vertex &v=node(i);
    adj.clear();
    adj.reserve(degree(i));
    int j;
    iset s;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        j=*it;
        if (include_temp_edges || !is_temporary_edge(i,j))
            s.insert(j);
    }
    if (is_directed()) {
        node_iter it=nodes.begin(),itend=nodes.end();
        for (j=0;it!=itend;++it,++j) {
            if (i!=j && it->has_neighbor(i) &&
                    (include_temp_edges || !is_temporary_edge(i,j)))
                s.insert(j);
        }
    }
    for (iset::const_iterator it=s.begin();it!=s.end();++it) {
        adj.push_back(*it);
    }
}

/* return a maximal independent set of vertices in undirected graph */
void graphe::maximal_independent_set(ivector &ind) const {
    int n=node_count(),i;
    ivector V(n),gain(n);
    for (i=0;i<n;++i) {
        V[i]=i;
        gain[i]=degree(i);
    }
    ivector::iterator pos;
    ind.clear();
    ind.reserve(n);
    while (!V.empty()) {
        i=V.front();
        pos=V.begin();
        for (ivector::iterator it=V.begin();it!=V.end();++it) {
            if (gain[i]<gain[*it]) {
                i=*it;
                pos=it;
            }
        }
        V.erase(pos);
        ind.push_back(i);
        const vertex &v=node(i);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            if ((pos=find(V.begin(),V.end(),*it))==V.end())
                continue;
            V.erase(pos);
            const vertex &w=node(*it);
            for (ivector_iter jt=w.neighbors().begin();jt!=w.neighbors().end();++jt) {
                if (find(V.begin(),V.end(),*jt)==V.end())
                    continue;
                gain[*jt]++;
            }
        }
    }
    std::sort(ind.begin(),ind.end());
}

/* find a maximum-cardinality matching */
void graphe::maximum_matching(ipairs &matching,int sg,gt_conn_check cc) {
    mm mm_finder(this);
    ipairs cm;
    int s,cmp;
    switch (cc) {
    case _GT_CC_CONNECTED:
        mm_finder.find_maximum_matching(matching,sg);
        break;
    case _GT_CC_FIND_COMPONENTS:
    case _GT_CC_COMPONENTS_ARE_SUBGRAPHS:
        s=(cc==_GT_CC_COMPONENTS_ARE_SUBGRAPHS || sg<0)?1:1+max_subgraph_index();
        cmp=(cc==_GT_CC_FIND_COMPONENTS)?connected_components_to_subgraphs():max_subgraph_index();
        for (;s<=cmp;++s) {
            mm_finder.find_maximum_matching(cm,s);
            matching.insert(matching.end(),cm.begin(),cm.end());
        }
        break;
    default: assert(false);
    }
}

/* find a maximal matching in an undirected graph (fast algorithm) */
void graphe::maximal_matching(ipairs &matching,int sg) const {
    assert(!is_directed());
    int i,j,N=node_count(),n=sg<0?N:subgraph_size(sg);
    ivector match(n,-1),V,V_pos(N,-1);
    bvector node_marked(n,false);
    if (sg>=0)
        get_subgraph(sg,V);
    else {
        V.resize(n);
        for (i=0;i<n;++i) V[i]=i;
    }
    if (sg>=0)
        for (ivector_iter it=V.begin();it!=V.end();++it) V_pos[*it]=it-V.begin();
    while (true) {
        for (i=0;i<n;++i) {
            if (!node_marked[i] && match[i]<0)
                break;
        }
        if (i==n)
            break;
        node_marked[i]=true;
        const vertex &v=node(V[i]);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            if (sg>=0 && node(*it).subgraph()!=sg)
                continue;
            if (match[sg<0?*it:V_pos[*it]]<0) {
                match[sg<0?*it:V_pos[*it]]=i;
                match[i]=sg<0?*it:V_pos[*it];
                break;
            }
        }
    }
    ivector skip;
    matching.clear();
    for (i=0;i<int(match.size());++i) {
        j=match[i];
        if (j<0 || find(skip.begin(),skip.end(),i)!=skip.end())
            continue;
        skip.push_back(j);
        matching.push_back(make_pair(V[i],V[j]));
    }
}

int graphe::pred(int i,int n) {
    if (i>0)
        return i-1;
    return n-1;
}

int graphe::succ(int i,int n) {
    if (i<n-1)
        return i+1;
    return 0;
}

/* return arc path from i-th to j-th vertex (inclusive), going in clockwise direction */
void graphe::arc_path(int i,int j,const ivector &cycle,ivector &path) {
    int n=cycle.size(),k=i,m=j-i,l=0;
    if (m<0)
        m+=n;
    path.resize(m+1);
    path.front()=cycle[i];
    while (k++!=j) {
        if (k==n)
            k=0;
        path[++l]=cycle[k];
    }
}

/* find the main chords of the given face */
void graphe::find_chords(const ivector &face,ipairs &chords) {
    int n=face.size(),h=n%2==0?n/2:(n-1)/2;
    ivector D;
    for (int k=0;k<n;++k) {
        if (degree(face[k])>2)
            D.push_back(k);
    }
    if (D.empty())
        return;
    chords.clear();
    int m=D.size(),k=0,i,j,k_next,l,v,w,r,f;
    map<int,bool> face_switch;
    ivector ifaces;
    while (true) {
        i=D[k];
        v=face[i];
        l=k;
        k_next=succ(k,m);
        ipair c=make_pair(i,-1);
        face_switch.clear();
        node(v).incident_faces(ifaces);
        for (ivector_iter it=ifaces.begin();it!=ifaces.end();++it) {
            face_switch[*it]=true;
        }
        while (true) {
            l=succ(l,m);
            j=D[l];
            r=j-i;
            if (r<0)
                r+=n;
            if (r==0 || r>h || (n%2==0 && r==h && i>=h))
                break;
            w=face[j];
            node(w).incident_faces(ifaces);
            f=-1;
            for (ivector_iter it=ifaces.begin();it!=ifaces.end();++it) {
                if (face_switch[*it])
                    face_switch[f=*it]=false;
            }
            if (r>1 && (f>=0 || has_edge(v,w))) {
                c.second=j;
                k_next=l;
            }
        }
        if (c.second>=0)
            chords.push_back(c);
        if (k_next<k)
            break;
        k=k_next;
    }
    if (chords.size()>1) {
        int end=chords.back().second;
        if (end<chords.back().first) {
            ipairs::iterator it=chords.begin();
            while (it->first<end) ++it;
            chords.erase(chords.begin(),it);
        }
    }
}

/* fold the face along its main chords */
void graphe::fold_face(const ivector &face,bool subdivide,int &label) {
    ipairs chords;
    find_chords(face,chords);
    int n=face.size(),k=chords.size(),i,j,p,q,r,s,t,u;
    if (k==0)
        return;
    if (subdivide) {
        bvector visited(n,false);
        i=add_node(++label);
        for (ipairs_iter it=chords.begin();it!=chords.end();++it) {
            p=it->first;
            q=it->second;
            visited[p]=visited[q]=true;
            for (j=(p+1)%n;j!=q;j=(j+1)%n) {
                add_edge(i,face[j]);
                visited[j]=true;
            }
        }
        for (j=0;j<n;++j) {
            if (!visited[j])
                add_edge(i,face[j]);
        }
    } else if (k==1) {
        i=chords.front().first;
        q=succ(i,n);
        t=pred(i,n);
        add_temporary_edge(face[q],face[t]);
        if (n>4) {
            j=chords.front().second;
            r=pred(j,n);
            u=succ(j,n);
            add_temporary_edge(face[r],face[u]);
        }
    } else if (k>1) {
        q=succ(chords.front().first,n);
        r=succ(chords[1].first,n);
        if (k==2 && q==pred(chords.front().second,n) && r==pred(chords[1].second,n)) {
            add_temporary_edge(face[q],face[r]);
        } else {
            ivector P(0);
            ivector Q;
            for (s=0;s<k;++s) {
                i=s<k-1?chords[s+1].first:chords.front().first;
                j=chords[s].second;
                q=pred(j,n);
                r=succ(i,n);
                add_temporary_edge(face[q],face[r]);
                arc_path(succ(chords[s].first,n),q,face,Q);
                P.insert(P.end(),Q.begin(),Q.end());
            }
            if (P.size()>3)
                fold_face(P,subdivide,label);
        }
    }
}

/* augment the biconnected planar graph for spring drawing */
void graphe::augment(const ivectors &faces,int outer_face,bool subdivide) {
    clear_embedding();
    set_embedding(faces);
    int label=largest_integer_label();
    for (ivectors_iter it=faces.begin();it!=faces.end();++it) {
        if (it->size()<4 || outer_face==int(it-faces.begin()))
            continue;
        fold_face(*it,subdivide,label);
    }
}

/* add point b to point a */
void graphe::add_point(point &a,const point &b) {
    int d=a.size();
    assert(int(b.size())>=d);
    for (int i=0;i<d;++i) {
        a[i]+=b[i];
    }
}

/* subtract point b from point a */
void graphe::subtract_point(point &a,const point &b) {
    int d=a.size();
    assert(int(b.size())>=d);
    for (int i=0;i<d;++i) {
        a[i]-=b[i];
    }
}

/* scale coordinates of point p by factor s */
void graphe::scale_point(point &p,double s) {
    int d=p.size();
    for (int i=0;i<d;++i) {
        p[i]*=s;
    }
}

/* compute v_x*w_y-v_y*w_x */
double graphe::point_vecprod2d(const point &v,const point &w) {
    assert(v.size()==2 && w.size()==2);
    return v[0]*w[1]-v[1]*w[0];
}

/* compute dot product of vectors p and q */
double graphe::point_dotprod(const point &p,const point &q) {
    int n=p.size();
    assert((n==2 || n==3) && n==int(q.size()));
    double res=0;
    for (int i=0;i<n;++i) {
        res+=p[i]*q[i];
    }
    return res;
}

/* convert point to Giac representation (cplx or point(x,y,[z])) */
gen graphe::point2gen(const point &p,bool vect) {
    /*
if (p.size()==2)
    return makecomplex(p[0],p[1]);
*/
    vecteur coords;
    for (point::const_iterator it=p.begin();it!=p.end();++it) {
        coords.push_back(*it);
    }
    if (vect)
        return coords;
    return symbolic(at_point,_feuille(coords,context0));
}

/* compute the square root of sum of squares of coordinates of point p */
double graphe::point_displacement(const point &p,bool sqroot) {
    double norm=0,d;
    for (point::const_iterator it=p.begin();it!=p.end();++it) {
        d=*it;
        norm+=d*d;
    }
    return sqroot?std::sqrt(norm):norm;
}

/* return the distance between points p and q */
double graphe::point_distance(const point &p,const point &q,point &pq) {
    copy_point(q,pq);
    subtract_point(pq,p);
    return point_displacement(pq);
}

/* compute the point dest which is a reflection of src wrt. the line ax+by+c=0 */
void graphe::point_mirror(double a,double b,double c,const point &src,point &dest) {
    double p=src[0],q=src[1],a2=a*a,b2=b*b,r=a2+b2,s=a2-b2;
    dest.resize(2);
    dest[0]=(p*s-2*b*(a*q+c))/r;
    dest[1]=-(q*s+2*a*(b*p+c))/r;
}

/* write d1*p+d2*q to res */
void graphe::point_lincomb(const point &p,const point &q,double d1,double d2,point &res) {
    copy_point(p,res);
    if (d2==0)
        scale_point(res,d1);
    else  {
        scale_point(res,d1/d2);
        add_point(res,q);
        scale_point(res,d2);
    }
}

/* copy layout src to dest (both must be initialized first) */
void graphe::copy_layout(const layout &src,layout &dest) {
    layout_iter st=src.begin();
    layout::iterator dt=dest.begin();
    for (;st!=src.end() && dt!=dest.end();++st,++dt) {
        *dt=*st;
    }
}

/* copy point src to dest (both must be initialized first) */
void graphe::copy_point(const point &src,point &dest) {
    point::const_iterator st=src.begin();
    point::iterator dt=dest.begin();
    for (;st!=src.end() && dt!=dest.end();++st,++dt) {
        *dt=*st;
    }
}

/* generate uniformly random point p on the origin-centered sphere with the specified radius */
void graphe::rand_point(point &p,double radius) {
    double R;
    do {
        R=0;
        for (point::iterator it=p.begin();it!=p.end();++it) {
            *it=rand_normal();
            R+=std::pow(*it,2.0);
        }
    } while (R==0);
    scale_point(p,radius/std::sqrt(R));
}

/* convert point from cartesian to polar representation (2d only) */
void graphe::point2polar(point &p,double &r,double &phi) {
    int d=p.size();
    assert(d==2);
    double x=p[0],y=p[1];
    r=std::sqrt(std::pow(x,2.0)+std::pow(y,2.0));
    phi=atan2(y,x);
}

/* translate layout x by dx */
void graphe::translate_layout(layout &x,const point &dx) {
    for (layout::iterator it=x.begin();it!=x.end();++it) {
        add_point(*it,dx);
    }
}

/* rotate layout x about the point p by the angle phi (2d only) */
void graphe::rotate_layout(layout &x,double phi) {
    if (x.empty())
        return;
    double r,phi0;
    for (layout::iterator it=x.begin();it!=x.end();++it) {
        point &p=*it;
        point2polar(p,r,phi0);
        p[0]=r*std::cos(phi0+phi);
        p[1]=r*std::sin(phi0+phi);
    }
}

/* determine minimum of d-th components among layout points */
double graphe::layout_min(const layout &x,int d) {
    double m=DBL_MAX;
    for (layout_iter it=x.begin();it!=x.end();++it) {
        if (it->at(d)<m)
            m=it->at(d);
    }
    return m;
}

/* determine the diameter of the layout (distance between the pair of farthest vertices) */
double graphe::layout_diameter(const layout &x) {
    point p(2);
    double diam=0,d;
    for (layout_iter it=x.begin();it!=x.end();++it) {
        for (layout_iter jt=it+1;jt!=x.end();++jt) {
            copy_point(*it,p);
            subtract_point(p,*jt);
            d=point_displacement(p,false);
            if (d>diam)
                diam=d;
        }
    }
    return std::sqrt(diam);
}

/* return the bounding rectangle of layout x */
graphe::rectangle graphe::layout_bounding_rect(layout &ly,double padding,int coord_offset) {
    double xmin=DBL_MAX,xmax=-DBL_MAX,ymin=DBL_MAX,ymax=-DBL_MAX,x,y;
    for (layout_iter it=ly.begin();it!=ly.end();++it) {
        x=it->at(0+coord_offset);
        y=it->at(1+coord_offset);
        if (x<xmin)
            xmin=x;
        if (x>xmax)
            xmax=x;
        if (y<ymin)
            ymin=y;
        if (y>ymax)
            ymax=y;
    }
    return rectangle(xmin-padding,ymin-padding,xmax-xmin+2*padding,ymax-ymin+2*padding,&ly);
}

/* compute the center of layout x */
graphe::point graphe::layout_center(const layout &x) {
    assert(!x.empty());
    int d=x.front().size();
    point center(d);
    clear_point_coords(center);
    for (layout_iter it=x.begin();it!=x.end();++it) {
        add_point(center,*it);
    }
    scale_point(center,1.0/x.size());
    return center;
}

/* scale layout x to have the diameter equal to diam */
void graphe::scale_layout(layout &x,double diam) {
    if (x.empty())
        return;
    int d=x.front().size();
    point M(d,-DBL_MAX),m(d,DBL_MAX);
    for (layout_iter it=x.begin();it!=x.end();++it) {
        const point &p=*it;
        for (int i=0;i<d;++i) {
            if (p[i]<m[i])
                m[i]=p[i];
            if (p[i]>M[i])
                M[i]=p[i];
        }
    }
    double D=0;
    for (int i=0;i<d;++i) {
        if (M[i]-m[i]>D)
            D=M[i]-m[i];
    }
    if (D==0)
        return;
    double s=diam/D;
    for (layout::iterator it=x.begin();it!=x.end();++it) {
        scale_point(*it,s);
    }
}

/* scale 2D layout along the x- and y-axis separately */
void graphe::scale_layout_2d(layout &x,double sx,double sy) {
    rectangle br=layout_bounding_rect(x);
    for (graphe::layout::iterator it=x.begin();it!=x.end();++it) {
        for (int i=0;i<2;++i) {
            double c=i==0?br.x():br.y();
            it->at(i)=c+(it->at(i)-c)*(i==0?sx:sy);
        }
    }
}

/* randomize layout x using a Gaussian random generator */
void graphe::create_random_layout(layout &x,int dim) {
    for (layout::iterator it=x.begin();it!=x.end();++it) {
        it->resize(dim);
        it->at(0)=rand_uniform();
        it->at(1)=rand_uniform();
        if (dim==3)
            it->at(2)=rand_uniform();
    }
}

/* lay out the graph using a force-directed algorithm with spring-electrical model */
void graphe::force_directed_placement(layout &x,double K,double R,double tol,bool ac) {
    double step_length=K,shrinking_factor=0.9,eps=K*tol,C=0.01,D=C*K*K;
    double energy=DBL_MAX,energy0,norm,max_displacement;
    int progress=0,n=x.size(),i,j;
    if (n==0)
        return;
    assert (n==node_count() && n>0);
    int d=x.front().size();
    point force(d),p(d),f(d);
    /* keep updating the positions until the system freezes */
    do {
        energy0=energy;
        energy=0;
        max_displacement=0;
        node_iter nt=nodes.begin(),ntend=nodes.end();
        for (i=0;nt!=ntend;++nt,++i) {
            point &xi=x[i];
            clear_point_coords(force);
            /* compute the attractive forces between vertices adjacent to the i-th vertex */
            for (ivector_iter it=nt->neighbors().begin();it!=nt->neighbors().end();++it) {
                j=*it;
                scale_point(p,point_distance(xi,x[j],p)/K);
                add_point(force,p);
            }
            /* compute the repulsive forces for all vertices j!=i which are not too far from the i-th vertex */
            for (layout_iter it=x.begin();it!=x.end();++it) {
                j=it-x.begin();
                if (i!=j) {
                    norm=point_distance(*it,xi,f);
                    if (norm>R)
                        continue;
                    if (norm==0)
                        rand_point(f,norm=shrinking_factor*eps);
                    scale_point(f,D/(norm*norm));
                    add_point(force,f);
                }
            }
            /* move the location of the i-th vertex in the direction of the force f */
            norm=point_displacement(force);
            if (norm==0)
                continue;
            if (step_length<norm) {
                scale_point(force,step_length/norm);
                norm=step_length;
            }
            add_point(xi,force);
            /* update the maximal displacement for this iteration */
            if (norm>max_displacement)
                max_displacement=norm;
            energy+=norm*norm;
        }
        /* update step length */
        if (ac) {
            /* adaptive cooling scheme */
            if (energy<energy0) {
                ++progress;
                if (progress>=5) {
                    progress=0;
                    step_length/=shrinking_factor;
                }
            } else {
                progress=0;
                step_length*=shrinking_factor;
            }
        } else step_length*=shrinking_factor; /* simple cooling scheme */
    } while (max_displacement>eps);
}

/* compute optimal positions of edge labels and store them as "position" attributes of the respective edges */
void graphe::edge_labels_placement(const layout &x) {
    if (x.empty())
        return;
    int dim=x.front().size();
    ipairs E;
    get_edges_as_pairs(E);
    int n=E.size();
    vector<double> D(n);
    bool initialize=n>300;
    point pq(dim);
    bool isdir=is_directed();
    for (int i=0;i<n;++i) {
        ipair &edge=E[i];
        D[i]=point_distance(x[edge.first],x[edge.second],pq);
        if (initialize || D[i]==0)
            set_edge_attribute(edge.first,edge.second,_GT_ATTRIB_POSITION,isdir?0.6:0.5);
    }
    if (!initialize) {
        vector<double> dist;
        int i,k0;
        double maxs,s;
        point c(dim);
        for (ipairs_iter it=E.begin();it!=E.end();++it) {
            i=it-E.begin();
            double &d=D[i];
            if (d==0)
                continue;
            const point &p=x[it->first];
            dist.clear();
            for (ipairs_iter jt=E.begin();jt!=E.end();++jt) {
                if (edges_incident(*it,*jt))
                    continue;
                if (edges_crossing(*it,*jt,x,c))
                    dist.push_back(point_distance(p,c,pq)/d);
            }
            std::sort(dist.begin(),dist.end());
            dist.insert(dist.begin(),MARGIN_FACTOR);
            for (int k=dist.size();k-->1;) {
                if (dist[k]<=dist.front())
                    dist.erase(dist.begin()+k);
            }
            while (dist.size()>1 && dist.back()>=1.0-MARGIN_FACTOR) {
                dist.pop_back();
            }
            dist.push_back(1.0-MARGIN_FACTOR);
            n=dist.size();
            maxs=0;
            k0=-1;
            for (int k=0;k<n-1;++k) {
                s=dist[k+1]-dist[k];
                if (s>maxs) {
                    maxs=s;
                    k0=k;
                }
            }
            assert(k0>=0);
            set_edge_attribute(it->first,it->second,_GT_ATTRIB_POSITION,1.0-(dist[k0]+maxs*(isdir?0.4:0.5)));
        }
    }
}

/* retrieve element in A at position (i,j) and store it in val, return true iff there is such element */
bool graphe::sparse_matrix_element(const sparsemat &A,int i,int j,ipair &val) {
    sparsemat::const_iterator it;
    map<int,ipair>::const_iterator jt;
    if ((it=A.find(i))==A.end() || (jt=it->second.find(j))==it->second.end())
        return false;
    val=jt->second;
    return true;
}

/* convert a rational number r to ipair */
graphe::ipair graphe::rat2ipair(const gen &r) {
    gen num=_numer(r,context0),den=_denom(r,context0);
    assert(num.is_integer() && den.is_integer());
    return make_pair(num.val,den.val);
}

/* convert a pair of integers to rational number */
gen graphe::ipair2rat(const ipair &p) {
    if (p.first==0)
        return 0;
    assert(p.second!=0);
    return fraction(p.first,p.second);
}

/* compute the product A*B and store it in P (ncols is equal to the number of columns of A),
* if A*B=B*A and A and B are symmetric then enabling 'symmetric=true' speeds up the computation */
void graphe::multiply_sparse_matrices(const sparsemat &A,const sparsemat &B,sparsemat &P,int ncols,bool symmetric) {
    int i,isempty;
    ipair p;
    gen elm;
    for (sparsemat::const_iterator it=A.begin();it!=A.end();++it) {
        i=it->first;
        sparsematrow &row=P[i];
        isempty=true;
        for (int j=symmetric?i:0;j<ncols;++j) {
            elm=0;
            for (sparsematrow::const_iterator jt=it->second.begin();jt!=it->second.end();++jt) {
                if (sparse_matrix_element(B,jt->first,j,p))
                    elm+=ipair2rat(jt->second)*ipair2rat(p);
            }
            if (!is_exactly_zero(elm)) {
                isempty=false;
                row[j]=rat2ipair(elm);
                if (symmetric)
                    P[j][i]=rat2ipair(elm);
            }
        }
        if (isempty)
            P.erase(P.find(i));
    }
}

/* return the scalar product of i-th row in A and j-th column in B */
gen graphe::sparse_product_element(const sparsemat &A,const sparsemat &B, int i, int j) {
    sparsemat::const_iterator it=A.find(i),rt;
    sparsematrow::const_iterator jt,bt;
    if (it==A.end())
        return 0;
    const sparsematrow &row=it->second;
    int c;
    gen res(0);
    for (jt=row.begin();jt!=row.end();++jt) {
        c=jt->first;
        if ((rt=B.find(c))!=B.end() && (bt=rt->second.find(j))!=rt->second.end())
            res+=ipair2rat(jt->second)*ipair2rat(bt->second);
    }
    return res;
}

/* store the transposition of A in T */
void graphe::transpose_sparsemat(const sparsemat &A, sparsemat &T) {
    for (sparsemat::const_iterator it=A.begin();it!=A.end();++it) {
        for(map<int,ipair>::const_iterator jt=it->second.begin();jt!=it->second.end();++jt) {
            T[jt->first][it->first]=jt->second;
        }
    }
}

/* return number of vertices in V which are incident to the i-th vertex */
int graphe::mdeg(const ivector &V, int i) const {
    int d=0;
    for (ivector_iter it=V.begin();it!=V.end();++it) {
        if (i==*it)
            return 1;
        if (has_edge(i,*it))
            ++d;
    }
    return d;
}

/* coarsening of the graph with respect to the prolongation matrix P */
void graphe::coarsening(graphe &G,const sparsemat &P,const ivector &V) const {
    sparsemat Q,I,R,IG;
    int n=node_count();
    transpose_sparsemat(P,Q);
    /* create sparse symmetric incidence matrix I of this graph */
    node_iter it=nodes.begin(),itend=nodes.end();
    int i,j;
    for (i=0;it!=itend;++it,++i) {
        const ivector &ngh=it->neighbors();
        if (ngh.empty())
            continue;
        map<int,ipair> &row=I[i];
        for (ivector_iter jt=ngh.begin();jt!=ngh.end();++jt) {
            row[*jt]=make_pair(1,1);
        }
    }
    /* use Galerkin product Q*I*P as the incidence matrix IG for graph G */
    multiply_sparse_matrices(Q,I,R,n);
    multiply_sparse_matrices(R,P,IG,n);
    for (ivector_iter it=V.begin();it!=V.end();++it) {
        G.add_node(node_label(*it));
    }
    sparsemat::const_iterator jt=IG.begin(),jtend=IG.end();
    map<int,ipair>::const_iterator kt,ktend;
    for (;jt!=jtend;++jt) {
        i=jt->first;
        for (kt=jt->second.begin(),ktend=jt->second.end();kt!=ktend;++kt) {
            j=kt->first;
            if (i<j && kt->second.first!=0)
                G.add_edge(i,j);
        }
    }
}

/* make coarser graph G by restricting to maximal independent set */
void graphe::coarsening_mis(const ivector &V,graphe &G,sparsemat &P) const {
    int n=node_count(),m=V.size();
    for (int i=0;i<n;++i) {
        int md=mdeg(V,i);
        assert(md>0);
        for (int j=0;j<m;++j) {
            if (i==V[j])
                P[i][j]=make_pair(1,1);
            else if (has_edge(i,V[j]))
                P[i][j]=make_pair(1,md);
        }
    }
    coarsening(G,P,V);
}

/* make coarser graph by collapsing edges found in maximal matching */
void graphe::coarsening_ec(const ipairs &M,graphe &G,sparsemat &P) const {
    ivector removed_nodes,V;
    for (ipairs_iter it=M.begin();it!=M.end();++it) {
        removed_nodes.push_back(it->second);
    }
    int n=node_count(),m=removed_nodes.size(),J;
    ivector::iterator ivt;
    for (int i=0;i<n;++i) {
        if ((ivt=find(removed_nodes.begin(),removed_nodes.end(),i))==removed_nodes.end())
            V.push_back(i);
        else
            removed_nodes.erase(ivt);
    }
    for (int i=0;i<n;++i) {
        for (int j=0;j<n-m;++j) {
            if (i==(J=V[j]) || find(M.begin(),M.end(),make_pair(J,i))!=M.end())
                P[i][j]=make_pair(1,1);
        }
    }
    coarsening(G,P,V);
}

static int multilevel_depth;
static bool multilevel_mis;

void graphe::multilevel_recursion(layout &x,int d,double R,double K,double tol,int depth) {
    ivector mis;
    ipairs M;
    if (multilevel_mis)
        maximal_independent_set(mis);
    else
        maximal_matching(M);
    int n=node_count(),m=multilevel_mis?mis.size():n-int(M.size());
    x.resize(n);
    if (m>0.75*n) {
        /* coarsening is slow, switch from EC to MIS method since it's faster */
        multilevel_mis=true;
        multilevel_recursion(x,d,R,K,tol,depth);
        return;
    }
    if (m<2) {
        /* the coarsest level, apply force directed algorithm on a random initial layout */
        multilevel_depth=depth;
        create_random_layout(x,d);
        force_directed_placement(x,K,R*(depth+1)*K,tol,false);
    } else {
        /* create coarser graph H and lay it out */
        graphe G(ctx);
        sparsemat P; // prolongation matrix
        if (multilevel_mis)
            coarsening_mis(mis,G,P);
        else
            coarsening_ec(M,G,P);
        layout y;
        G.multilevel_recursion(y,d,R,K,tol,depth+1);
        /* compute x=P*y (layout lifting) */
        ipair pij;
        point yj(d);
        for (int i=0;i<n;++i) {
            x[i]=point(d,0);
            for (int j=0;j<m;++j) {
                if (sparse_matrix_element(P,i,j,pij)) {
                    copy_point(y[j],yj);
                    scale_point(yj,to_real_number(ipair2rat(pij),ctx).to_double(ctx));
                    add_point(x[i],yj);
                }
            }
        }
        /* make the natural spring length K shorter with respect to
     * the current depth level and subsequently refine x */
        double L=K*std::pow(PLASTIC_NUMBER,depth-multilevel_depth);
        force_directed_placement(x,L,R*(depth+1)*L,K*tol/L,false);
    }
}

/* apply multilevel force directed algorithm to layout x */
void graphe::make_spring_layout(layout &x,int d,double tol) {
    int n=node_count();
    if (n==0)
        return;
    if (n==1) {
        x.resize(1);
        x.front().resize(d,0);
    } else if (n<30) {
        x.resize(n);
        create_random_layout(x,d);
        force_directed_placement(x,10.0,DBL_MAX,tol);
    } else {
        multilevel_mis=false;
        multilevel_recursion(x,d,DBL_MAX,10.0,tol);
    }
    if (d==3) {
        /* z-center the layout */
        double minz=DBL_MAX,maxz=-DBL_MAX;
        for (layout_iter it=x.begin();it!=x.end();++it) {
            if (it->back()>maxz)
                maxz=it->back();
            if (it->back()<minz)
                minz=it->back();
        }
        double dz=-(minz+maxz)/2.0;
        for (layout::iterator it=x.begin();it!=x.end();++it) {
            it->back()+=dz;
        }
    }
}

/* layout face as a regular polygon inscribed in circle of radius R */
void graphe::make_regular_polygon_layout(layout &x,const ivector &v,double R,double elongate) {
    int n=v.size(),i;
    double step=2.0*M_PI/(double)n,phi=(n%2)==0?M_PI_2*(1+2.0/n):M_PI_2;
    for (int k=0;k<n;++k) {
        i=v[k];
        point &p=x[i];
        p.resize(2);
        p[0]=R*std::cos(phi);
        p[1]=R*std::sin(phi);
        if (elongate>0) {
            if (p[0]>.001) p[0]+=elongate;
            else if (p[0]<-.001) p[0]-=elongate;
        }
        phi-=step;
    }
}

/* apply force directed algorithm to this graph (must be triconnected), with the specified outer hull,
* using the algorithm by Bor Plestenjak in "An Algorithm for Drawing Planar Graphs" */
void graphe::make_circular_layout(layout &x,const ivector &hull,double A,double tol,double elongate) {
    int n=node_count(),iter_count=0,maxper=0,j;
    /* place facial vertices on the unit circle and all other vertices in the origin */
    x.resize(n);
    make_regular_polygon_layout(x,hull,1.0,elongate);
    if (n==int(hull.size())) // there are no vertices to place inside the circle
        return;
    bvector is_hull_vertex(n,false);
    for (ivector_iter it=hull.begin();it!=hull.end();++it) {
        is_hull_vertex[*it]=true;
    }
    /* compute vertex periphericity */
    ivector per(n);
    periphericity(hull,per);
    for (ivector_iter pt=per.begin();pt!=per.end();++pt) {
        if (*pt>maxper)
            maxper=*pt;
    }
    layout force(n);
    for (int i=0;i<n;++i) {
        force[i].resize(2);
        if (is_hull_vertex[i])
            continue;
        point &p=x[i];
        p.resize(2,0.0);
    }
    /* cool down the system iteratively */
    double d,cool,C=std::sqrt((double)n/M_PI),eps=tol/std::sqrt((double)n);
    do {
        ++iter_count;
        /* compute the resultant force for each non-outer vertex */
        for (int i=0;i<n;++i) {
            if (is_hull_vertex[i])
                continue;
            vertex &v=node(i);
            point &f=force[i],&p=x[i];
            clear_point_coords(f);
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                j=*it;
                point &q=x[j];
                point r(2);
                copy_point(q,r);
                subtract_point(r,p);
                d=C*std::exp(A*(2.0*maxper-per[i]-per[j])/(double)maxper)*point_displacement(r,false);
                scale_point(r,d);
                add_point(f,r);
            }
        }
        /* move each vertex in the direction of the force */
        cool=1.0/(C+std::pow(iter_count,1.5)/C);
        for (int i=0;i<n;++i) {
            if (is_hull_vertex[i])
                continue;
            point &f=force[i],&p=x[i];
            d=point_displacement(f);
            if (d>0) {
                scale_point(f,std::min(d,cool)/d);
                add_point(p,f);
            }
        }
    } while (cool>eps);
}

/* apply Tutte's barycentric method for vertex placement */
void graphe::make_tutte_layout(layout &x,const ivector &outer_face) {
    /* place facial vertices on the unit circle and all other vertices in the origin */
    int n=node_count();
    x.resize(n);
    make_regular_polygon_layout(x,outer_face);
    if (n==int(outer_face.size()))
        return; /* there are no vertices to place inside the circle */
    bvector is_outer(n,false);
    for (ivector_iter it=outer_face.begin();it!=outer_face.end();++it) {
        is_outer[*it]=true;
    }
    for (int i=0;i<n;++i) {
        if (!is_outer[i]) {
            point &p=x[i];
            p.resize(2);
            p[0]=p[1]=0;
        }
    }
    double tol=1e-5;
    point old(2),r(2);
    bool converged;
    int j;
    do {
        converged=true;
        for (int i=0;i<n;++i) {
            if (is_outer[i])
                continue;
            vertex &v=node(i);
            point &p=x[i];
            copy_point(p,old);
            clear_point_coords(p);
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                j=*it;
                add_point(p,x[j]);
            }
            scale_point(p,1.0/(double)v.degree());
            if (point_distance(p,old,r)>tol)
                converged=false;
        }
    } while (!converged);
}

/* Tomita recursive algorithm (a variant of Bron-Kerbosch) for maximal cliques */
void graphe::tomita(iset &R,iset &P,iset &X,map<int,int> &m,int mode) {
    if (P.empty() && X.empty()) {
        if (mode==1) { // store matching
            if (R.size()==2)
                m[*R.begin()]=*R.rbegin();
        } else if (mode==2 || mode==3) { // store the clique
            ivector cliq;
            cliq.reserve(R.size());
            for (iset_iter it=R.begin();it!=R.end();++it) cliq.push_back(*it);
            maxcliques.push_back(cliq);
        }
        if (mode==0 || mode==3)
            ++m[R.size()];
    } else {
        /* choose as the pivot element the node with the highest number of neighbors in P, say i-th */
        iset PUX,PP,XX,S;
        sets_union(P,X,PUX);
        int i=*PUX.begin(),mn=0,n;
        for (iset_iter it=PUX.begin();it!=PUX.end();++it) {
            const vertex &w=node(*it);
            sets_intersection(w.neighbors(),P,S);
            if ((n=int(S.size()))>mn) {
                mn=n;
                i=*it;
            }
        }
        sets_difference(P,node(i).neighbors(),S);
        for (iset_iter it=S.begin();it!=S.end();++it) {
            const ivector &ngh=node(*it).neighbors();
            R.insert(*it);
            sets_intersection(ngh,P,PP);
            sets_intersection(ngh,X,XX);
            tomita(R,PP,XX,m,mode);
            R.erase(*it);
            P.erase(*it);
            X.insert(*it);
        }
    }
}

/* a modified version of the Bron-Kerbosch algorithm for listing all maximal cliques,
* developed by Tomita et al. Number of k-cliques will be stored to m[k] for each k.
* If store_matching is true, store 2-cliques (v,w) as m[v]=w. */
void graphe::clique_stats(map<int,int> &m,int mode) {
    iset R,X,P;
    for (int i=0;i<node_count();++i) P.insert(i);
    tomita(R,P,X,m,mode);
}

/* generate a clique vertex cover of this graph */
void graphe::find_maximal_cliques() {
    map<int,int> m;
    iset R,X,P;
    for (int i=0;i<node_count();++i) P.insert(i);
    clear_maximal_cliques();
    tomita(R,P,X,m,2);
}

/* CP recursive subroutine */
void graphe::cp_recurse(ivector &C,ivector &P,ivector &incumbent) {
    if (C.size()>incumbent.size())
        incumbent=C;
    if (C.size()+P.size()>incumbent.size()) {
        int p,k;
        ivector Q,Cup;
        while (!P.empty()) {
            p=P.back();
            P.pop_back();
            Cup=C;
            if (find(Cup.begin(),Cup.end(),p)==Cup.end())
                Cup.push_back(p);
            const vertex &v=node(p);
            Q.clear();
            Q.resize(std::min(P.size(),v.neighbors().size()));
            k=0;
            for (ivector_iter it=P.begin();it!=P.end();++it) {
                if (v.has_neighbor(*it))
                    Q[k++]=*it;
            }
            Q.resize(k);
            cp_recurse(Cup,Q,incumbent);
        }
    }
}

/* Carraghan-Pardalos algorithm for finding maximum clique in a sparse graph */
int graphe::cp_maxclique(ivector &clique) {
    int n=node_count();
    ivector P(n),C;
    /* initialize P and sort all adjacency lists */
    for (int i=0;i<n;++i)
        P[i]=i;
    /* sort vertices by degree in descending order */
    degree_comparator comp(this);
    std::sort(P.begin(),P.end(),comp);
    std::reverse(P.begin(),P.end());
    cp_recurse(C,P,clique);
    return clique.size();
}

/*
* Östergård class for finding maximum clique
*
* implemented as described in "A fast algorithm for the maximum clique problem",
* Discrete Applied Mathematics 120 (2002) 197–207
*/

void graphe::ostergard::recurse(ivector &U,int size,ivector &position) {
    if (U.empty()) {
        if (size>maxsize) {
            maxsize=size;
            incumbent=clique_nodes;
            found=true;
        }
        return;
    }
    if (timeout>0 && double(clock()-start)/CLOCKS_PER_SEC>timeout && !incumbent.empty()) {
        timed_out=true;
        return;
    }
    int i,j,minpos,p;
    ivector::iterator it;
    ivector V;
    while (!U.empty()) {
        if (size+int(U.size())<=maxsize)
            break;
        i=U.front();
        it=U.begin();
        minpos=-1;
        for (ivector::iterator jt=U.begin();jt!=U.end();++jt) {
            j=*jt;
            p=position[j];
            if (minpos<0 || p<minpos) {
                minpos=p;
                i=j;
                it=jt;
            }
        }
        const vertex &v=G->node(i);
        assert(v.low()>0);
        if (size+v.low()<=maxsize)
            break;
        U.erase(it);
        V.resize(std::min(U.size(),v.neighbors().size()));
        j=0;
        for (it=U.begin();it!=U.end();++it) {
            if (*it<*v.neighbors().begin() || *it>*v.neighbors().rbegin())
                continue;
            if (v.has_neighbor(*it))
                V[j++]=*it;
        }
        V.resize(j);
        clique_nodes.push_back(i);
        recurse(V,size+1,position);
        clique_nodes.pop_back();
        if (found || timed_out)
            break;
    }
}

/* return the size of maximum clique */
int graphe::ostergard::maxclique(ivector &clique) {
    int n=G->node_count();
    ivector S(n),U,position(n,0);
    for (int i=0;i<n;++i) {
        S[i]=i;
        vertex &v=G->node(i);
        v.set_low(0);
    }
    G->greedy_vertex_coloring_biggs(S);
    std::reverse(S.begin(),S.end());
    G->node(S.back()).set_low(1);
    for (ivector_iter it=S.begin();it!=S.end();++it) {
        position[*it]=it-S.begin();
    }
    maxsize=0;
    clique_nodes.clear();
    int k;
    start=clock();
    timed_out=false;
    for (int i=n;i-->0;) {
        found=false;
        k=S[i];
        vertex &v=G->node(k);
        U.clear();
        for (ivector_iter it=S.begin()+i;it!=S.end();++it) {
            if (v.has_neighbor(*it))
                U.push_back(*it);
        }
        clique_nodes.push_back(k);
        recurse(U,1,position);
        if (timed_out)
            break;
        clique_nodes.pop_back();
        v.set_low(maxsize);
    }
    clique=incumbent;
    return clique.size();
}

/*
* End of ostergard class
*/

/* find maximum clique in this graph and return its size */
int graphe::maximum_clique(ivector &clique) {
    assert(!is_directed());
    clique.clear();
    int n=node_count(),m=edge_count();
    if (20*m<=n*(n-1)) // edge density smaller than or equal to 0.1
        return cp_maxclique(clique);
    ostergard ost(this);
    return ost.maxclique(clique);
}

/* remove a maximal clique from vertex set V using a fast greedy algorithm
* (Johnson, J. Comp. Syst. Sci. 1974) */
void graphe::remove_maximal_clique(iset &V) const {
    iset U=V,W,tmp;
    int s,smax,v;
    do {
        smax=0;
        v=*U.begin();
        W.clear();
        for (iset_iter it=U.begin();it!=U.end();++it) {
            if ((s=sets_intersection(node(*it).neighbors(),U,tmp))>smax) {
                smax=s;
                W=tmp;
                v=*it;
            }
        }
        V.erase(v);
    } while (!(U=W).empty());
}

/* approximate neighborhood clique cover number for every vertex */
void graphe::greedy_neighborhood_clique_cover_numbers(ivector &cover_numbers) {
    int n=node_count(),i;
    cover_numbers.resize(n,0);
    node_iter it=nodes.begin(),itend=nodes.end();
    ivector_iter jt,jtend;
    for (i=0;it!=itend;++it,++i) {
        iset V;
        for (jt=it->neighbors().begin(),jtend=it->neighbors().end();jt!=jtend;++jt) {
            V.insert(*jt);
        }
        int &cn=cover_numbers[i];
        while (!V.empty()) {
            remove_maximal_clique(V);
            ++cn;
        }
    }
}

/* return the number of differently colored and uncolored vertices adjacent to the i-th vertex */
graphe::ipair graphe::adjacent_color_count(int i) const {
    std::set<int> colors;
    int c,unc=0;
    const vertex &v=node(i);
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        c=node(*it).color();
        if (c>0)
            colors.insert(c);
        else ++unc;
    }
    return make_pair(colors.size(),unc);
}

/* insert colors adjacent to the i-th vertex to the given set */
bool graphe::adjacent_colors(int i,std::set<int> &colors) const {
    const vertex &v=node(i);
    if (v.color()>0)
        return false;
    int c;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        if ((c=node(*it).color())>0)
            colors.insert(c);
    }
    return true;
}

/*
* painter class implementation
*/

#ifdef HAVE_LIBGLPK
void graphe::painter::compute_bounds(const ivector &icol,int max_colors) {
    G->greedy_neighborhood_clique_cover_numbers(cover_number);
    G->uncolor_all_nodes();
    lb=(initially_colored=icol).size();
    ub=max_colors;
    if (ub==0) {
        for (ivector_iter it=initially_colored.begin();it!=initially_colored.end();++it) {
            G->set_node_color(*it,it-initially_colored.begin()+1);
        }
        G->dsatur();
        ub=G->color_count();
    }
    if ((generate_clique_cuts=G->node_count()<=100))
        G->find_maximal_cliques();
}

void graphe::painter::make_values() {
    int n=G->node_count(),col=0,k;
    ivector_iter pos;
    values.resize(n);
    for (int i=0;i<n;++i) {
        ivector &x=values[i];
        x.resize(ub);
        pos=find(initially_colored.begin(),initially_colored.end(),i);
        k=pos==initially_colored.end()?-1:int(pos-initially_colored.begin());
        for (int j=0;j<ub;++j) {
            if (k<0) {
                if (j<lb && G->node(initially_colored[j]).has_neighbor(i))
                    x[j]=-2;
                else {
                    x[j]=++col;
                    col2ij.push_back(make_pair(i,j));
                }
            } else
                x[j]=j==k?-1:-2;
        }
    }
    nxcols=col;
}

void graphe::painter::formulate_mip() {
    mip=glp_create_prob();
    glp_set_obj_dir(mip,GLP_MIN);
    int n=G->node_count(),i,j,k;
    iscliq.resize(n,false);
    for (ivector_iter it=initially_colored.begin();it!=initially_colored.end();++it) {
        iscliq[*it]=true;
    }
    /* create row variables and count nonzero entries in the constraint matrix */
    int nonzeros=0,cnt=0,rs,cn,val;
    int nrows=n*(ub+2)-3*lb;
    glp_add_rows(mip,nrows);
    for (i=0;i<n;++i) {
        if (!iscliq[i]) {
            glp_set_row_bnds(mip,++cnt,GLP_FX,1.0,1.0);
            ivector &vals=values[i];
            for (ivector_iter it=vals.begin();it!=vals.end();++it) {
                if (*it>0)
                    ++nonzeros;
            }
        }
    }
    for (k=0;k<n;++k) {
        cn=cover_number[k];
        const vertex &v=G->node(k);
        for (j=0;j<ub;++j) {
            if (j<lb && k==initially_colored[j])
                continue;
            nonzeros+=2;
            rs=j<lb?cn:0;
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                i=*it;
                if ((val=values[i][j])>0)
                    ++nonzeros;
                else rs-=val+2;
            }
            glp_set_row_bnds(mip,++cnt,GLP_UP,0.0,rs);
        }
    }
    for (i=0;i<n;++i) {
        if (!iscliq[i]) {
            glp_set_row_bnds(mip,++cnt,GLP_UP,0.0,lb);
            nonzeros+=ub-lb;
            ivector &vals=values[i];
            for (ivector_iter it=vals.begin();it!=vals.end();++it) {
                if (*it>0)
                    ++nonzeros;
            }
        }
    }
    assert(cnt==nrows);
    /* create column variables */
    int ncols=nxcols+ub-lb;
    glp_add_cols(mip,ncols);
    for (i=0;i<ncols;++i) {
        glp_set_col_kind(mip,i+1,GLP_BV);
        if (i>=nxcols)
            glp_set_obj_coef(mip,i+1,1.0);
    }
    /* create the constraint matrix */
    int *ia=new int[nonzeros+1];
    int *ja=new int[nonzeros+1];
    double *ar=new double[nonzeros+1];
    int row=0,col;
    bool ok;
    cnt=0;
    for (i=0;i<n;++i) {
        if (iscliq[i])
            continue;
        ++row;
        ivector &vals=values[i];
        ok=false;
        for (ivector_iter it=vals.begin();it!=vals.end();++it) {
            if ((col=*it)>0) {
                ia[++cnt]=row; ja[cnt]=col;
                ar[cnt]=1.0;
                ok=true;
            } else assert(col==-2);
        }
        assert(ok);
    }
    for (k=0;k<n;++k) {
        cn=cover_number[k];
        const vertex &v=G->node(k);
        for (j=0;j<ub;++j) {
            if (j<lb && k==initially_colored[j])
                continue;
            ++row;
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                col=values[*it][j];
                if (col>0) {
                    ia[++cnt]=row; ja[cnt]=col;
                    ar[cnt]=1.0;
                }
            }
            if ((col=values[k][j])>0) {
                ia[++cnt]=row; ja[cnt]=col;
                ar[cnt]=cn;
            }
            if (j>=lb) {
                col=1+nxcols+j-lb;
                ia[++cnt]=row; ja[cnt]=col;
                ar[cnt]=-cn;
            }
        }
    }
    for (i=0;i<n;++i) {
        if (iscliq[i])
            continue;
        ++row;
        ivector &vals=values[i];
        for (ivector_iter it=vals.begin();it!=vals.end();++it) {
            j=it-vals.begin();
            if ((col=*it)>0) {
                ia[++cnt]=row; ja[cnt]=col;
                ar[cnt]=j+1;
            }
        }
        for (col=nxcols+1;col<=ncols;++col) {
            ia[++cnt]=row; ja[cnt]=col;
            ar[cnt]=-1.0;
        }
    }
    assert(cnt<=nonzeros && row==nrows);
    /* assign the constraint matrix to the MIP */
    glp_load_matrix(mip,cnt,ia,ja,ar);
    delete[] ia; delete[] ja;
    delete[] ar;
    for (i=0;i<nrows;++i) {
        glp_set_row_stat(mip,i+1,GLP_BS);
    }
}

void graphe::painter::callback(glp_tree *tree,void *info) {
    if (is_interrupted()) {
        glp_ios_terminate(tree);
        return;
    }
    painter *pt=static_cast<painter*>(info);
    int j;
    switch (glp_ios_reason(tree)) {
    case GLP_ISELECT:
        if (pt->select_blb)
            glp_ios_select_node(tree,glp_ios_best_node(tree));
        break;
    case GLP_IBRANCH:
        if (//_rand(makesequence(0,1),pt->G->giac_context()).DOUBLE_val()<0.1 &&
                (j=pt->select_branching_variable(tree))>0)
            glp_ios_branch_upon(tree,j,GLP_DN_BRNCH);
        break;
    case GLP_IHEUR:
        pt->heur_solution(tree);
        break;
    case GLP_IROWGEN:
        pt->generate_rows(tree);
        break;
    case GLP_IBINGO:
        pt->select_blb=true;
        break;
    default:
        break;
    }
}

/* generate constraints violated by the current subproblem relaxation solution */
void graphe::painter::generate_rows(glp_tree *tree) {
    glp_prob *prob=glp_ios_get_prob(tree);
    int n=G->node_count(),col,w,len,r=glp_ios_curr_node(tree),best_len=-1;
    int *data=(int*)glp_ios_node_data(tree,r);
    if (*data>maxiter)
        return;
    (*data)+=1;
    double rh,lh,h,maxh=0;
    /* block color cuts */
    for (int c=ub;c-->lb;) {
        w=nxcols+c-lb+1;
        for (int i=0;i<n;++i) {
            if (iscliq[i])
                continue;
            rh=0.0;
            len=0;
            for (int j=c;j<ub;++j) {
                col=values[i][j];
                if (col<0)
                    rh-=col+2;
                else {
                    row_indices[++len]=col;
                    row_coeffs[len]=1.0;
                }
            }
            if (len==0)
                continue;
            row_indices[++len]=w;
            row_coeffs[len]=-1.0;
            lh=0;
            for (int k=0;k<len;++k) {
                lh+=row_coeffs[k+1]*glp_get_col_prim(prob,row_indices[k+1]);
            }
            h=lh-rh;
            //h/=len;
            if (h>maxh) {
                maxh=h;
                memcpy(best_indices+1,row_indices+1,len*sizeof(int));
                memcpy(best_coeffs+1,row_coeffs+1,len*sizeof(double));
                best_len=len;
            }
        }
    }
    if (best_len>0)
        add_row(prob,best_len,best_indices,best_coeffs,rh);
    /* clique cuts */
    if (generate_clique_cuts) {
        const ivectors &cover=G->maximal_cliques();
        double lim=30.0/double(cover.size());
        int nc=0;
        for (int c=ub;c-->lb;) {
            w=nxcols+c-lb+1;
            for (ivectors_iter it=cover.begin();it!=cover.end();++it) {
                if (G->rand_uniform()>=lim)
                    continue;
                const ivector &cliq=*it;
                rh=0.0;
                len=0;
                for (ivector_iter jt=cliq.begin();jt!=cliq.end();++jt) {
                    col=values[*jt][c];
                    if (col<0)
                        rh-=col+2;
                    else {
                        row_indices[++len]=col;
                        row_coeffs[len]=1.0;
                    }
                }
                if (len==0)
                    continue;
                row_indices[++len]=w;
                row_coeffs[len]=-1.0;
                lh=0;
                for (int k=0;k<len;++k) {
                    lh+=row_coeffs[k+1]*glp_get_col_prim(prob,row_indices[k+1]);
                }
                h=lh-rh;
                //h/=len;
                if (h>0) {
                    add_row(prob,len,row_indices,row_coeffs,rh);
                    ++nc;
                }
            }
        }
    }
}

void graphe::painter::add_row(glp_prob *prob,int len,int *indices,double *coeffs,double rh) {
    int r=glp_add_rows(prob,1);
    glp_set_row_bnds(prob,r,GLP_UP,-DBL_MAX,rh);
    glp_set_mat_row(prob,r,len,indices,coeffs);
}

/* color the nodes for which the colors are known */
void graphe::painter::fixed_coloring(glp_tree *tree) {
    glp_prob *prob=glp_ios_get_prob(tree); // current subproblem
    G->uncolor_all_nodes();
    for (ivector_iter it=initially_colored.begin();it!=initially_colored.end();++it) {
        G->set_node_color(*it,int(it-initially_colored.begin())+1);
    }
    for (int j=nxcols;j-->0;) {
        ipair &ij=col2ij[j];
        if (!glp_ios_can_branch(tree,j+1) && glp_get_col_prim(prob,j+1)==1)
            G->set_node_color(ij.first,lb+ij.second+1);
    }
}

/* fast construction of a heuristic solution using greedy coloring */
void graphe::painter::heur_solution(glp_tree *tree) {
    if (maxiter==0) {
        fixed_coloring(tree);
        /* dsatur coloring, quadratic time */
        G->dsatur();
        assign_heur(tree);
    }
    int iter_count=0;
    while (++iter_count<=maxiter) {
        /* greedy coloring, linear time */
        fixed_coloring(tree);
        int n=G->node_count(),i,j;
        ordering=G->rand_permu(n);
        for (ivector_iter ot=ordering.begin();ot!=ordering.end();++ot) {
            i=*ot;
            used_colors.clear();
            if (!G->adjacent_colors(i,used_colors))
                continue;
            j=0;
            std::set<int>::const_iterator it=used_colors.begin();
            for (;it!=used_colors.end();++it) {
                if (*it!=++j) {
                    G->set_node_color(i,j);
                    break;
                }
            }
            if (it==used_colors.end())
                G->set_node_color(i,j+1);
        }
        if (assign_heur(tree)==0)
            break;
    }
}

/* offer heuristic solution to MIP solver */
int graphe::painter::assign_heur(glp_tree *tree) {
    G->get_node_colors(temp_colors);
    int nc=G->color_count(),i,j;
    if (nc>ub)
        return 1;
    nc-=lb;
    for (int col=1;col<=nxcols;++col) {
        const ipair &p=col2ij[col-1];
        i=p.first;
        j=p.second;
        *(heur+col)=temp_colors[i]==j+1?1.0:0.0;
    }
    for (i=0;i<ub-lb;++i) {
        *(heur+nxcols+1+i)=i<nc?1.0:0.0;
    }
    return glp_ios_heur_sol(tree,heur);
}

/* the custom branching rule for MVCP (Brelaz): a fractionable vertex with
* the largest number of different adjacent colors is chosen, ties are broken by
* choosing the vertex with the maximal degree in the uncolored subgraph */
int graphe::painter::select_branching_variable(glp_tree *tree) {
    glp_prob *prob=glp_ios_get_prob(tree); // current subproblem
    std::fill(branch_candidates.begin(),branch_candidates.end(),-1);
    fixed_coloring(tree);
    double fr;
    for (int j=nxcols;j-->0;) {
        ipair &ij=col2ij[j];
        if (glp_ios_can_branch(tree,j+1)) {
            int &bc=branch_candidates[ij.first];
            fr=std::abs(glp_get_col_prim(prob,j+1)-0.5);
            if (bc<0 || fr<fracts[ij.first]) {
                bc=ij.second;
                fracts[ij.first]=fr;
            }
        }
    }
    int i,j=0;
    ipair nc,nc_max=make_pair(-1,-1);
    for (ivector_iter it=branch_candidates.begin();it!=branch_candidates.end();++it) {
        if (*it<0)
            continue;
        i=it-branch_candidates.begin();
        nc=G->adjacent_color_count(i);
        if (nc>nc_max) {
            nc_max=nc;
            j=values[i][*it];
        }
    }
    return j;
}

int graphe::painter::color_vertices(ivector &colors,const ivector &icol,int max_colors,int tm_lim,bool verbose) {
    compute_bounds(icol,max_colors);
    if (ub<lb) {
        if (verbose)
            G->message("Error: there are no feasible solutions");
        return 0;
    }
    int n=G->node_count();
    make_values();
    formulate_mip();
    glp_smcp lparm;
    glp_init_smcp(&lparm);
    lparm.msg_lev=GLP_MSG_OFF;
    glp_simplex(mip,&lparm);
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.msg_lev=GLP_MSG_ON;
    parm.gmi_cuts=GLP_OFF;
    parm.mir_cuts=GLP_OFF;
    parm.clq_cuts=GLP_ON;
    parm.cov_cuts=GLP_ON;
    parm.presolve=GLP_OFF;
    parm.fp_heur=GLP_ON;
    parm.sr_heur=GLP_ON;
    parm.br_tech=GLP_BR_MFV;
    parm.bt_tech=GLP_BT_DFS;
    parm.cb_size=sizeof(int);
    parm.cb_func=&callback;
    parm.cb_info=static_cast<void*>(this);
    if (tm_lim>0) parm.tm_lim=tm_lim;
    branch_candidates.resize(n);
    fracts.resize(n);
    maxiter=std::ceil(15.0*std::exp(-std::pow(n,2)*0.000321887582487));
    heur=new double[nxcols+ub-lb+1];
    row_indices=new int[nxcols+ub-lb+1];
    row_coeffs=new double[nxcols+ub-lb+1];
    best_indices=new int[nxcols+ub-lb+1];
    best_coeffs=new double[nxcols+ub-lb+1];
    select_blb=false;
    int term_old=0;
    if (verbose)
        glp_term_hook(term_hook,static_cast<void*>(logptr(G->giac_context())));
    else term_old=glp_term_out(GLP_OFF);
    int ncolors=0,res=glp_intopt(mip,&parm);
    if (res==0 || res==GLP_ETMLIM || res==GLP_ESTOP) {
        switch (glp_mip_status(mip)) { // solution status:
        case GLP_UNDEF:
            if (verbose)
                G->message("Error: undefined solution");
            break;
        case GLP_FEAS: // feasible but not necessarily optimal
            if (verbose)
                G->message("Warning: the solution is not necessarily optimal");
        case GLP_OPT: // optimal
            ncolors=lb+glp_mip_obj_val(mip);
            /* color the vertices */
            for (int k=0;k<nxcols;++k) {
                ipair &p=col2ij[k];
                if (glp_mip_col_val(mip,k+1)>0)
                    G->set_node_color(p.first,p.second+1);
            }
            G->get_node_colors(colors);
            break;
        case GLP_NOFEAS: // not feasible
            if (verbose)
                G->message("Error: there are no feasible solutions");
            break;
        }
    } else if (verbose)
        G->message("Error: solver failure");
    glp_delete_prob(mip);
    if (verbose)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
    delete[] heur;
    delete[] row_indices;
    delete[] row_coeffs;
    delete[] best_indices;
    delete[] best_coeffs;
    return ncolors;
}
#endif

/*
* end of painter class
*/

/* find optimal vertex coloring using an exact algorithm, requires GLPK */
int graphe::exact_vertex_coloring(int max_colors,int tm_lim,bool verbose) {
    int ncolors=0,n=node_count();
    if (is_clique()) {
        for (int i=0;i<n;++i) {
            set_node_color(i,i+1);
        }
        return n;
    }
    ivector smp;
    find_simplicial_vertices(smp);
    if (!smp.empty()) {
        uncolor_all_nodes();
        ivector nsv,vcols,colors;
        int sz=n-smp.size();
        if (sz>0) {
            nsv.reserve(sz);
            for (int i=0;i<n;++i) {
                if (find(smp.begin(),smp.end(),i)==smp.end())
                    nsv.push_back(i);
            }
            graphe G(ctx,false);
            induce_subgraph(nsv,G);
            ncolors=G.exact_vertex_coloring(max_colors,tm_lim,verbose);
            if (ncolors==0)
                return 0;
            colors.resize(ncolors);
            if (!colors.empty())
                G.get_node_colors(vcols);
            for (int i=0;i<sz;++i) {
                int c=G.node(i).color();
                set_node_color(nsv[i],c);
                colors[c-1]=c;
            }
        }
        for (ivector_iter it=smp.begin();it!=smp.end();++it) {
            bool colored=false;
            const ivector &ngh=node(*it).neighbors();
            for (ivector_iter ct=colors.begin();!colored && ct!=colors.end();++ct) {
                ivector_iter jt=ngh.begin();
                for (;jt!=ngh.end();++jt) {
                    if (node(*jt).color()==*ct)
                        break;
                }
                if (jt==ngh.end()) {
                    set_node_color(*it,*ct);
                    colored=true;
                }
            }
            if (!colored) { // add new color
                colors.push_back(++ncolors);
                set_node_color(*it,ncolors);
            }
        }
        return ncolors;
    }
#ifndef HAVE_LIBGLPK
    message("Error: GLPK library is required for graph coloring");
#else
    painter pt(this);
    ivector colors,clique;
    ostergard ost(this,5.0);
    ost.maxclique(clique);
    ncolors=pt.color_vertices(colors,clique,max_colors,tm_lim,verbose);
    if (ncolors>0 && find(colors.begin(),colors.end(),0)!=colors.end()) {
        uncolor_all_nodes();
        ncolors=0;
    }
#endif
    return ncolors;
}

/* color the edges of this graph using D or D+1 colors,
* return 1 in the former case, 2 in the latter, 0 on error */
int graphe::exact_edge_coloring(ivector &colors,int *numcol,int tm_lim,bool verbose) {
    assert(!is_directed());
    /* apply heuristic */
    if (edge_coloring_heuristic(colors)) {
        if (numcol!=NULL)
            *numcol=maximum_degree();
        if (verbose)
            message("Solution found by heuristic");
        return 1;
    }
    /* exact edge coloring */
    graphe L(ctx,false);
    ipairs E;
    line_graph(L,E);
    if (E.empty()) {
        *numcol=0;
        return 2;
    }
    /* find the vertex with maximum degree in this graph */
    int m=E.size(),maxdeg=0,deg,i=-1,j,k;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if ((deg=it->neighbors().size())>maxdeg) {
            maxdeg=deg;
            i=it-nodes.begin();
        }
    }
    assert(i>=0);
    ivector icol(maxdeg);
    k=0;
    for (j=0;j<m;++j) {
        const ipair &e=E[j];
        if (e.first==i || e.second==i)
            icol[k++]=j;
    }
    assert(k==maxdeg);
#ifdef HAVE_LIBGLPK
    painter pt(&L);
    int ncolors=pt.color_vertices(colors,icol,maxdeg+1,tm_lim,verbose);
#else
    message("Error: GLPK library is required for graph coloring");
    int ncolors=0;
#endif
    if (colors.empty()) // infeasible or time limit reached
        return 0;
    for (k=0;k<maxdeg;++k) {
        colors[icol[k]]=k+1;
    }
    if (ncolors==0 || find(colors.begin(),colors.end(),0)!=colors.end())
        return 0;
    for (j=0;j<m;++j) {
        const ipair &e=E[j];
        set_edge_attribute(e.first,e.second,_GT_ATTRIB_COLOR,colors[j]);
    }
    if (numcol!=NULL)
        *numcol=ncolors;
    return ncolors-maxdeg+1;
}

/* edge-coloring heuristic by finding maximum matching iteratively,
 * returns true iff this graph is in class 1 */
bool graphe::edge_coloring_heuristic(ivector &colors) {
    assert(!is_directed());
    int maxdeg=maximum_degree();
    ivector sigma=vecteur_2_vector_int(*_randperm(node_count(),ctx)._VECTptr);
    graphe G(ctx);
    if (!isomorphic_copy(G,sigma,true))
        return false;
    vector<ipairs> matchings;
    while (int(matchings.size())<maxdeg) {
        matchings.resize(matchings.size()+1);
        ipairs &matching=matchings.back();
        G.maximum_matching(matching);
        for (ipairs_iter it=matching.begin();it!=matching.end();++it) {
            G.remove_edge(*it);
        }
        if (G.is_empty())
            break;
    }
    assert(int(matchings.size())==maxdeg);
    if (!G.is_empty())
        return false;
    colors.resize(edge_count());
    ipairs E;
    get_edges_as_pairs(E);
    int col=1;
    for (vector<ipairs>::iterator it=matchings.begin();it!=matchings.end();++it) {
        for (ipairs::iterator jt=it->begin();jt!=it->end();++jt) {
            jt->first=sigma[jt->first];
            jt->second=sigma[jt->second];
            if (jt->first>jt->second)
                *jt=make_pair(jt->second,jt->first);
            ipairs_iter kt=find(E.begin(),E.end(),*jt);
            assert(kt!=E.end());
            colors[kt-E.begin()]=col;
        }
        col++;
    }
    return true;
}

/* returns true iff there is a clique cover of order not larger than k and finds that cover */
bool graphe::clique_cover(ivectors &cover,int k,int tm_lim,bool verbose) {
    if (triangle_count()==0) {
        /* clique cover consists of matched edges and singleton vertex sets */
        ipairs matching;
        maximum_matching(matching);
        int m=matching.size(),n=node_count(),i=0;
        if (k>0 && n-m>k)
            return false;
        bvector matched(n);
        cover.resize(n-m);
        for (ipairs_iter it=matching.begin();it!=matching.end();++it) {
            ivector &clique=cover[i++];
            clique.resize(2);
            matched[(clique.front()=it->first)]=true;
            matched[(clique.back()=it->second)]=true;
        }
        for (bvector::const_iterator it=matched.begin();it!=matched.end();++it) {
            if (!*it) {
                ivector &clique=cover[i++];
                clique.resize(1);
                clique.front()=it-matched.begin();
            }
        }
        return true;
    }
    if (k>0) { // GRASP clique heuristic
        int sg=1;
        unset_subgraphs(0);
        bool success=true;
        ivector clq;
        while (subgraph_size(0)>0) {
            if (sg>k) {
                success=false;
                break;
            }
            grasp_clique(5,clq,false,0);
            set_subgraph(clq,sg++);
        }
        if (success) {
            cover.clear();
            cover.reserve(sg-1);
            for (int s=1;s<sg;++s) {
                get_subgraph(s,clq);
                cover.push_back(clq);
            }
            unset_subgraphs();
            if (verbose)
                message("Solution found by GRASP heuristic");
            return true;
        } else unset_subgraphs();
    }
    graphe C(ctx,false);
    complement(C);
    int ncliques;
    /* DSATUR heuristic */
    if (k>0)
        C.dsatur();
    if (k<=0 || (ncliques=C.color_count())>k)
        ncliques=C.exact_vertex_coloring(k,tm_lim,verbose);
    else if (verbose)
        message("Solution found by DSATUR heuristic");
    if (ncliques==0)
        return false;
    cover.clear();
    cover.resize(ncliques);
    for (int i=node_count();i-->0;) {
        const vertex &v=C.node(i);
        cover[v.color()-1].push_back(i);
    }
    return true;
}

/* return true iff the graph is complete */
bool graphe::is_clique(int sg) const {
    int n=node_count();
    bool isdir=is_directed();
    int ec_max;
    if (sg<0) {
        ec_max=n*(n-1);
        if (!isdir)
            ec_max/=2;
        return edge_count()==ec_max;
    }
    int m=subgraph_size(sg);
    ec_max=m*(m-1);
    if (!isdir)
        ec_max/=2;
    return edge_count(sg)==ec_max;
}

/* find maximum independent set in this graph and return its size */
int graphe::maximum_independent_set(ivector &v) const {
    graphe C(ctx,false);
    complement(C);
    return C.maximum_clique(v);
}

/* return the number of (directed) triangles in (di)graph */
gen graphe::triangle_count(ivectors *dest,bool ccoeff,bool exact) {
    if (is_directed()) {
        sparsemat M,M2;
        adjacency_sparse_matrix(M);
        multiply_sparse_matrices(M,M,M2,node_count(),false);
        ipair p;
        gen trace=0;
        for (sparsemat::const_iterator it=M.begin();it!=M.end();++it) {
            for (map<int,ipair>::const_iterator jt=it->second.begin();jt!=it->second.end();++jt) {
                if (sparse_matrix_element(M2,jt->first,it->first,p))
                    trace+=fraction(jt->second.first*p.first,jt->second.second*p.second);
            }
        }
        return _ratnormal(trace/gen(3),ctx);
    }
    gen c(0);
    int n=node_count(),i,j,sz1,sz2,sz3,isz,p1,p2,p3;
    double sum=0.0;
    ivector offset(n,0),intersection,trg(3),sigma;
    ivector_iter iter;
    graphe G(ctx,false);
    sort_by_degrees(G,sigma);
    node_iter it=nodes.begin(),itend=nodes.end();
    for (i=0;it!=itend;++it,++i) {
        const ivector &vn=it->neighbors();
        if (ccoeff) { sz1=vn.size(); p1=(sz1*(sz1-1))/2; }
        for (ivector_iter jt=vn.begin();jt!=vn.end();++jt) {
            if ((j=*jt)<i) continue;
            const ivector &wn=node(j).neighbors();
            if (ccoeff) { sz2=wn.size(); p2=(sz2*(sz2-1))/2; }
            if (!ccoeff && dest==NULL)
                c+=intersect_linear(vn.begin(),vn.begin()+offset[i],wn.begin(),wn.begin()+offset[j]);
            else {
                intersection.resize(std::min(offset[i],offset[j]));
                iter=set_intersection(vn.begin(),vn.begin()+offset[i],wn.begin(),wn.begin()+offset[j],intersection.begin());
                intersection.resize(isz=iter-intersection.begin());
                for (ivector_iter kt=intersection.begin();kt!=intersection.end();++kt) {
                    if (ccoeff) { sz3=node(*kt).neighbors().size(); p3=(sz3*(sz3-1))/2; }
                    if (!ccoeff) {
                        trg[0]=i; trg[1]=j; trg[2]=*kt;
                        dest->push_back(trg);
                        c+=1;
                    } else if (exact) c+=gen(1)/harmonic_mean_exact(p1,p2,p3);
                    else sum+=1.0/harmonic_mean(p1,p2,p3);
                }
            }
            ++offset[j];
        }
    }
    if (!ccoeff) return c;
    if (!exact) return 3.0*sum/double(n);
    return _ratnormal(gen(3)*c/gen(n),ctx);
}

/* remove i-th node which is assumed to be isolated */
void graphe::remove_isolated_nodes(const iset &I) {
    int n=node_count(),i;
    ivector sigma(n);
    for (i=0;i<n;++i) {
        sigma[i]=i;
    }
    iset::const_reverse_iterator it=I.rbegin(),itend=I.rend();
    for (;it!=itend;++it) {
        nodes.erase(nodes.begin()+*it);
        sigma.erase(sigma.begin()+*it);
        sigma.push_back(*it);
    }
    /* obtain the inverse of sigma */
    map<int,int> sigma_inv;
    ivector_iter jt=sigma.begin(),jtend=sigma.end();
    for (i=0;jt!=jtend;++jt,++i)
        sigma_inv[*jt]=i;
    /* change vertex indices in edges */
    vector<vertex>::iterator nt=nodes.begin(),ntend=nodes.end();
    for (;nt!=ntend;++nt)
        nt->map_neighbors(sigma_inv);
}

/* return true iff this is a multigraph */
bool graphe::is_multigraph() const {
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->has_multiedges())
            return true;
    }
    return false;
}

/* makes this graph simple by removing multiedges */
void graphe::make_simple() {
    vector<vertex>::iterator it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it)
        it->clear_multiedges();
}

/* return the number of copies of edge e */
int graphe::multiedges(const ipair &e) const {
    int i=e.first,j=e.second;
    return node(i<j?i:j).multiedges(i<j?j:i);
}

/* add k duplicates of e to the graph */
void graphe::set_multiedge(const ipair &e,int k) {
    int i=e.first,j=e.second;
    node(i<j?i:j).set_multiedge(i<j?j:i,k);
}

/* convert integral edge weights to multiedges, returns false if an edge
 * has a weight which is not a positive integer */
bool graphe::weights2multiedges() {
    assert(!is_directed() && is_weighted());
    ipairs E;
    get_edges_as_pairs(E);
    gen w;
    int k;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        const ipair &e=*it;
        w=weight(e);
        if (w.is_integer() && (k=w.val-1)>=0) {
            set_multiedge(e,k);
        } else return false;
    }
    return true;
}

/* contract the edge {i,j}, leaving j-th node isolated (not connected to any other node) */
void graphe::contract_edge(int i,int j,bool adjust_positions) {
    /* assuming that the graph is undirected */
    assert(has_edge(i,j) && multiedges(make_pair(i,j))==0);
    ivector adj;
    ipair e,e_old;
    adjacent_nodes(j,adj);
    int m=sum_of_edge_multiplicities();
    for (ivector_iter it=adj.begin();it!=adj.end();++it) {
        e_old=make_pair(j,*it);
        if (*it!=i) {
            e=make_pair(i,*it);
            if (has_edge(e)) {
                set_multiedge(e,multiedges(e_old)+multiedges(e)+1);
            } else {
                add_edge(e);
                set_multiedge(e,multiedges(e_old));
            }
        }
        remove_edge(e_old);
    }
    assert(m==1+sum_of_edge_multiplicities());
    if (adjust_positions) {
        vertex &v=node(i),&w=node(j);
        point p,q,r;
        if (get_node_position(v.attributes(),p,ctx) && get_node_position(w.attributes(),q,ctx) && p.size()==q.size()) {
            r.resize(p.size());
            copy_point(p,r);
            add_point(r,q);
            scale_point(r,0.5);
            v.set_attribute(_GT_ATTRIB_POSITION,point2gen(r));
        }
    }
}

/* subdivide the edge e with r new vertices */
void graphe::subdivide_edge(const ipair &e,int n,int &label) {
    assert(remove_edge(e));
    int i=e.first,j,nv=node_count();
    for (int k=0;k<n;++k) {
        j=add_node(++label);
        add_edge(i,j);
        i=j;
    }
    add_edge(i,e.second);
    vertex &v=node(e.first),&w=node(e.second);
    point p,q,r;
    if (get_node_position(v.attributes(),p,ctx) && get_node_position(w.attributes(),q,ctx) && p.size()==q.size()) {
        r.resize(p.size());
        copy_point(q,r);
        subtract_point(r,p);
        scale_point(r,1.0/double(n+1));
        for (i=nv;i<node_count();++i) {
            add_point(p,r);
            node(i).set_attribute(_GT_ATTRIB_POSITION,point2gen(p));
        }
    }
}

/* put all edges incident to nodes in V to set E */
void graphe::incident_edges(const ivector &V,edgeset &E) const {
    if (is_directed()) {
        ivector adj;
        for (ivector_iter it=V.begin();it!=V.end();++it) {
            adjacent_nodes(*it,adj);
            for (ivector_iter jt=adj.begin();jt!=adj.end();++jt) {
                if (has_edge(*it,*jt))
                    E.insert(make_pair(*it,*jt));
                if (has_edge(*jt,*it))
                    E.insert(make_pair(*jt,*it));
            }
        }
    } else {
        int i,j;
        for (ivector_iter it=V.begin();it!=V.end();++it) {
            i=*it;
            const vertex &v=node(i);
            for (ivector_iter jt=v.neighbors().begin();jt!=v.neighbors().end();++jt) {
                j=*jt;
                E.insert(make_pair(i<j?i:j,i<j?j:i));
            }
        }
    }
}

/* make a list of all n-tuples of k integers 0,1,...,k-1 */
void ntupk(graphe::ivectors &v,int n,int k,const graphe::ivector &elem,int i) {
    for (int j=0;j<k;++j) {
        graphe::ivector e(elem);
        e[i]=j;
        if (i<n-1)
            ntupk(v,n,k,e,i+1);
        else
            v.push_back(e);
    }
}

/* return true iff the vertices u and v of sierpinski graph are adjacent */
bool is_sierpinski_match(const graphe::ivector &u,const graphe::ivector &v,int n,int h) {
    int t=0;
    for (;t<=h-1;++t) {
        if (u[t]!=v[t])
            return false;
    }
    if (u[h]==v[h])
        return false;
    for (t=h+1;t<n;++t) {
        if (u[t]!=v[h] || v[t]!=u[h])
            return false;
    }
    return true;
}

/* create a graph from decreasing degree sequence L by using Havel-Hakimi algorithm,
* return false iff L is not a graphic sequence */
bool graphe::hakimi(const ivector &L) {
    assert(node_count()==int(L.size()));
    int n=L.size(),d,i,z;
    if (n==0)
        return true;
    ipairs D(n);
    for (int i=0;i<n;++i) {
        D[i]=make_pair(L[i],i);
    }
    do {
        std::sort(D.begin(),D.end());
        d=D.back().first;
        i=D.back().second;
        D.pop_back();
        for (int k=0;k<d;++k) {
            ipair &p=D[D.size()-1-k];
            if (--p.first<0)
                return false;
            add_edge(i,p.second);
        }
        z=0;
        for (ipairs_iter it=D.begin();it!=D.end();++it) {
            if (it->first==0)
                ++z;
        }
    } while (z<int(D.size()));
    return true;
}

/* return the (odd) girth of this graph (the length of the shortest (odd) cycle),
* return -1 if there are no cycles */
int graphe::girth(bool odd,int sg) {
    assert(node_queue.empty());
    int g=rand_max2,h,i,j;
    bool hascycle=false;
    node_iter it=nodes.begin(),itend=nodes.end();
    ivector_iter jt,jtend;
    for (i=0;it!=itend;++it,++i) {
        unvisit_all_nodes(sg);
        vertex &v=node(i);
        if (sg>=0 && v.subgraph()!=sg)
            continue;
        v.unset_ancestor();
        v.set_disc(0);
        node_queue.push(i);
        while (!node_queue.empty()) {
            j=node_queue.front();
            node_queue.pop();
            vertex &w=node(j);
            w.set_visited(true);
            for (jt=w.neighbors().begin(),jtend=w.neighbors().end();jt!=jtend;++jt) {
                vertex &u=node(*jt);
                if (*jt==w.ancestor() || (sg>=0 && u.subgraph()!=sg))
                    continue;
                if (!u.is_visited()) {
                    u.set_ancestor(j);
                    u.set_disc(w.disc()+1);
                    node_queue.push(*jt);
                } else {
                    h=w.disc()+u.disc()+1;
                    if (h<g && (!odd || h%2!=0)) {
                        hascycle=true;
                        g=h;
                    }
                }
            }
        }
    }
    return hascycle?g:-1;
}

/* create a graph from LCF notation [jumps]^e */
void graphe::make_lcf_graph(const ivector &jumps,int e) {
    this->clear();
    int m=jumps.size(),n=m*e;
    reserve_nodes(n);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,n);
        add_nodes(V);
    } else add_nodes(n);
    make_cycle_graph();
    int j=0,k;
    for (int i=0;i<n;++i) {
        k=(i+jumps[j])%n;
        if (k<0)
            k+=n;
        add_edge(i,k);
        j=(j+1)%m;
    }
}

/* create LCF graph from int pointer (list must end with 0) */
void graphe::make_lcf_graph(const int *j) {
    ivector jumps;
    const int *it=j;
    while (*it!=0) {
        jumps.push_back(*it);
        ++it;
    }
    make_lcf_graph(jumps,*(++it));
}

/* create Sierpinski (triangle) graph */
void graphe::make_sierpinski_graph(int n,int k,bool triangle) {
    this->clear();
    ivectors tuples;
    ivector elem(n,0);
    ntupk(tuples,n,k,elem,0);
    int N=std::pow(k,n);
    reserve_nodes(N);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,N,0);
        add_nodes(V);
    } else add_nodes(N);
    for (int i=0;i<N;++i) {
        ivector &u=tuples[i];
        for (int j=i+1;j<N;++j) {
            ivector &v=tuples[j];
            int h=0;
            for (;h<n;++h) {
                if (is_sierpinski_match(u,v,n,h))
                    break;
            }
            if (h<n)
                add_edge(i,j);
        }
    }
    if (triangle && k>1) {
        /* remove all non-clique edges to obtain Sierpinski triangle graph */
        map<int,int> m;
        clique_stats(m,1);
        iset isolated_nodes;
        for (map<int,int>::const_iterator it=m.begin();it!=m.end();++it) {
            int v=it->first,w=it->second;
            contract_edge(v,w,false);
            isolated_nodes.insert(w);
        }
        remove_isolated_nodes(isolated_nodes);
        if (supports_attributes()) {
            vecteur labels;
            make_default_labels(labels,node_count());
            relabel_nodes(labels);
        }
    }
}

/* create the Shrikhande graph */
void graphe::make_shrikhande_graph() {
    this->clear();
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,16);
        add_nodes(labels);
    } else add_nodes(16);
    ipairs v(16);
    int k=0,m,n;
    for (int i=0;i<4;++i) {
        for (int j=0;j<4;++j) {
            v[k++]=make_pair(i,j);
        }
    }
    for (int i=0;i<16;++i) {
        const ipair &vi=v[i];
        for (int j=i+1;j<16;++j) {
            const ipair &vj=v[j];
            m=(vi.first-vj.first+4)%4;
            n=(vi.second-vj.second+4)%4;
            if ((m*n==0 && (m+n)%2!=0) || (m==n && (m*n)%2!=0))
                add_edge(i,j);
        }
    }
}

/* create complete graph with vertices from list V */
void graphe::make_complete_graph() {
    int n=node_count();
    for (int i=0;i<n;++i) {
        for (int j=i+1;j<n;++j) {
            add_edge(i,j);
        }
    }
}

/* create complete k-partite graph with given partition sizes */
void graphe::make_complete_multipartite_graph(const ivector &partition_sizes,layout *x) {
    int k=partition_sizes.size();
    ivectors partitions;
    int ntotal=0,n;
    vecteur v;
    for (ivector_iter it=partition_sizes.begin();it!=partition_sizes.end();++it) {
        n=*it;
        if (supports_attributes()) {
            make_default_labels(v,n,ntotal);
            add_nodes(v);
        } else add_nodes(n);
        ivector iv(n);
        for (int i=0;i<n;++i) iv[i]=ntotal+i;
        partitions.push_back(iv);
        ntotal+=n;
    }
    assert(node_count()==ntotal);
    for (int ip=0;ip<k;++ip) {
        ivector &pi=partitions[ip];
        for (int jp=ip+1;jp<k;++jp) {
            ivector &pj=partitions[jp];
            for (int i=0;i<int(pi.size());++i) {
                for (int j=0;j<int(pj.size());++j) {
                    add_edge(pi[i],pj[j]);
                }
            }
        }
    }
    if (x!=NULL && partition_sizes.size()==2 && partition_sizes.front()==1 && partition_sizes.back()>2) {
        /* create the star layout */
        n=partition_sizes.back();
        ivector hull(n);
        for (int i=0;i<n;++i) hull[i]=i+1;
        make_circular_layout(*x,hull);
    }
}

/* create generalized Petersen graph G(n,k) using Watkins' notation */
void graphe::make_petersen_graph(int n,int k,layout *x) {
    this->clear();
    reserve_nodes(2*n);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,2*n);
        add_nodes(V);
    } else add_nodes(2*n);
    /* add the outer cycle first */
    for (int i=0;i<n;++i) {
        add_edge(i,(i+1)%n);
    }
    /* add the inner edges */
    for (int i=0;i<n;++i) {
        add_edge(i,i+n);
        add_edge(i+n,(i+k)%n+n);
    }
    if (x!=NULL) {
        /* layout the graph */
        ivector hull(n);
        for (int i=0;i<n;++i) hull[i]=i;
        make_circular_layout(*x,hull,2.5);
    }
}

/* generate all k-sets in a n-set */
void graphe::generate_nk_sets(int n,int k,vector<bitset<32> > &v) {
    ulong N=(1<<n);
    int i=0;
    for (ulong m=1;m<N;++m) {
        bitset<32> b(m);
        if (b.count()==(size_t)k)
            v[i++]=b;
    }
}

/* create Kneser or Johnson graph with parameters n (<=20) and k */
bool graphe::make_intersection_graph(int n,int k,int sz,bool lab_sets) {
    this->clear();
    assert(n>1 && n<21 && k>0 && k<n);
    int nchoosek=comb(n,k).val; // number of vertices
    reserve_nodes(nchoosek);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,nchoosek);
        add_nodes(V);
    } else add_nodes(nchoosek);
    vector<bitset<32> > vert(nchoosek);
    vecteur lab;
    if (supports_attributes() && lab_sets)
        lab.reserve(n);
    generate_nk_sets(n,k,vert);
    for (int i=0;i<nchoosek;++i) {
        const bitset<32> &a=vert[i];
        if (supports_attributes() && lab_sets) {
            lab.clear();
            for (int c=0;c<n;++c) {
                if (a.test(c)) lab.push_back(c+1);
            }
            node(i).set_label(lab);
        }
        for (int j=i+1;j<nchoosek;++j) {
            const bitset<32> &b=vert[j];
            if ((a & b).count()==size_t(sz))
                add_edge(i,j);
        }
    }
    return true;
}

/* create path graph with n vertices */
void graphe::make_path_graph() {
    int n=node_count();
    for (int i=0;i<n-1;++i) {
        add_edge(i,i+1);
    }
}

/* Create n times m grid graph.
 * mode=0 -- rectangular grid
 * mode=1 -- triangular grid
 * mode=2 -- toroidal grid */
void graphe::make_grid_graph(int m,int n,int mode) {
    this->clear();
    vecteur V;
    graphe X(ctx,supports_attributes()),Y(ctx,supports_attributes());
    X.reserve_nodes(m);
    Y.reserve_nodes(n);
    if (supports_attributes()) {
        X.make_default_labels(V,m);
        X.add_nodes(V);
        Y.make_default_labels(V,n);
        Y.add_nodes(V);
    } else {
        X.add_nodes(m);
        Y.add_nodes(n);
    }
    if (mode==1) {
        X.make_directed();
        Y.make_directed();
        for (int i=0;i<m-1;++i) X.add_edge(i,i+1);
        for (int i=0;i<n-1;++i) Y.add_edge(i,i+1);
        graphe CP(ctx,supports_attributes()),TP(ctx,supports_attributes());
        X.cartesian_product(Y,CP);
        X.tensor_product(Y,TP);
        ipairs E;
        TP.get_edges_as_pairs(E);
        for (ipairs_iter it=E.begin();it!=E.end();++it) {
            CP.add_edge(*it);
        }
        CP.underlying(*this);
    } else {
        if (mode==2) {
            X.make_cycle_graph();
            Y.make_cycle_graph();
        } else {
            X.make_path_graph();
            Y.make_path_graph();
        }
        X.cartesian_product(Y,*this);
    }
}

/* create n times m web graph as the cartesian product of n-cycle and m-path graphs */
void graphe::make_web_graph(int n,int m,layout *x) {
    this->clear();
    vecteur V;
    graphe C(ctx,supports_attributes());
    C.reserve_nodes(n);
    if (supports_attributes()) {
        C.make_default_labels(V,n);
        C.add_nodes(V);
    } else C.add_nodes(n);
    C.make_cycle_graph();
    graphe P(ctx,supports_attributes());
    P.reserve_nodes(m);
    if (supports_attributes()) {
        P.make_default_labels(V,m);
        P.add_nodes(V);
    } else P.add_nodes(m);
    P.make_path_graph();
    C.cartesian_product(P,*this);
    if (x!=NULL) {
        /* layout the graph */
        ivector hull(n);
        for (int i=0;i<n;++i) hull[i]=m*i;
        make_circular_layout(*x,hull);
    }
}

/* create wheel graph with n+1 vertices */
void graphe::make_wheel_graph(int n,layout *x) {
    this->clear();
    reserve_nodes(n+1);
    int i;
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,n,0,1);
        add_nodes(V);
        i=add_node(0);
    } else {
        add_nodes(n+1);
        i=n;
    }
    for (int j=0;j<n;++j) {
        add_edge(j,(j+1)%n);
    }
    for (int j=0;j<n;++j) {
        add_edge(i,j);
    }
    if (x!=NULL) {
        /* layout the graph */
        ivector hull(n);
        for (i=0;i<n;++i) hull[i]=i;
        make_circular_layout(*x,hull);
    }
}

/* create n-th flower snark, where n>=3 is odd */
bool graphe::make_flower_snark(int n,layout *x) {
    if (n<3 || n%2==0)
        return false;
    this->clear();
    ivector A(n),B(n),C(n),D(n);
    reserve_nodes(4*n);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,4*n);
        add_nodes(V);
    } else add_nodes(4*n);
    for (int i=0;i<4*n;++i) {
        switch (i%4) {
        case 0: A[i/4]=i; break;
        case 1: B[i/4]=i; break;
        case 2: C[i/4]=i; break;
        case 3: D[i/4]=i; break;
        default: assert(false);
        }
    }
    /* create n star graphs on four vertices */
    for (int i=0;i<n;++i) {
        add_edge(A[i],B[i]);
        add_edge(A[i],C[i]);
        add_edge(A[i],D[i]);
    }
    /* create cycles */
    for (int i=0;i<n;++i) add_edge(B[i],B[(i+1)%n]);
    for (int i=0;i<n-1;++i) add_edge(C[i],C[i+1]);
    for (int i=0;i<n-1;++i) add_edge(D[i],D[i+1]);
    add_edge(C[n-1],D[0]);
    add_edge(D[n-1],C[0]);
    /* layout */
    if (x!=NULL) {
        ivector hull(3*n);
        int cc=1,dc=0;
        for (int i=0;i<n;++i) {
            hull[3*i]=A[i];
            hull[3*i+1]=(i%2==0)?D[dc++]:C[cc++];
            hull[3*i+2]=i==n-1?C[0]:((i%2==0)?D[dc++]:C[cc++]);
        }
        make_circular_layout(*x,hull,n==3?1.5:-1.5);
    }
    return true;
}

bool graphe::make_goldberg_snark(int n) {
    if (n<3 || n%2==0)
        return false;
    this->clear();
    reserve_nodes(8*n);
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,8*n);
        add_nodes(labels);
    } else add_nodes(8*n);
    for (int t=0;t<n;++t) {
        int tn=(t+1)%n,tp=(t-1+n)%n;
        add_edge(8*t,8*tp+1);
        add_edge(8*t,8*t+1);
        add_edge(8*t,8*t+6);
        add_edge(8*t+1,8*tn);
        add_edge(8*t+1,8*t+7);
        add_edge(8*t+2,8*tp+3);
        add_edge(8*t+2,8*t+3);
        add_edge(8*t+2,8*t+7);
        add_edge(8*t+3,8*t+6);
        add_edge(8*t+3,8*tn+2);
        add_edge(8*t+4,8*t+5);
        add_edge(8*t+4,8*tp+4);
        add_edge(8*t+4,8*tn+4);
        add_edge(8*t+5,8*t+7);
        add_edge(8*t+5,8*t+6);
    }
    return true;
}

/* construct Haar graph with index n */
bool graphe::make_haar_graph(const gen &n) {
    if (!n.is_integer() || is_strictly_greater(1,n,ctx))
        return false;
    int k=1+_floor(_logb(makesequence(n,2),ctx),ctx).val;
    this->clear();
    vecteur V,be=*_convert(makesequence(n,_BASE,2),ctx)._VECTptr;
    reserve_nodes(2*k);
    if (supports_attributes()) {
        make_default_labels(V,2*k);
        add_nodes(V);
    } else add_nodes(2*k);
    for (int i=0;i<k;++i) {
        for (int j=k;j<2*k;++j) {
            if (giac::is_one(be[(j-i)%k]))
                add_edge(i,j);
        }
    }
    return true;
}

bool graphe::make_haar_graph(ulong n) {
    if (n==0)
        return false;
    int k=1+int(std::floor(std::log(n)/M_LN2));
    this->clear();
    bitset<64> be(n);
    vecteur V;
    reserve_nodes(2*k);
    if (supports_attributes()) {
        make_default_labels(V,2*k);
        add_nodes(V);
    } else add_nodes(2*k);
    for (int i=0;i<k;++i) {
        for (int j=k;j<2*k;++j) {
            if (be.test((j-i)%k))
                add_edge(i,j);
        }
    }
    return true;
}

/* create antiprism graph of order n (with 2n vertices and 4n edges) */
void graphe::make_antiprism_graph(int n,layout *x) {
    this->clear();
    reserve_nodes(2*n);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,2*n,0);
        add_nodes(V);
    } else add_nodes(2*n);
    int j;
    for (int i=0;i<n;++i) {
        j=(i+1)%n;
        add_edge(2*i,2*j);
        add_edge(2*i+1,2*j+1);
        add_edge(2*j,2*i+1);
        add_edge(2*j,2*j+1);
    }
    if (x!=NULL) {
        /* layout the graph */
        ivector hull(n);
        for (int i=0;i<n;++i) hull[i]=2*i;
        make_circular_layout(*x,hull,-1.0);
    }
}

/* create the complete k-ary tree with depth d */
void graphe::make_complete_kary_tree(int k,int d) {
    assert(k>=2);
    this->clear();
    int n=((int)std::pow(k,d+1)-1)/(k-1),v=0,w=1,len=1;
    reserve_nodes(n);
    if (supports_attributes()) {
        vecteur V;
        make_default_labels(V,n);
        add_nodes(V);
    } else add_nodes(n);
    for (int i=0;i<d;++i) {
        for (int j=0;j<len;++j) {
            for (int m=0;m<k;++m) {
                add_edge(v,w++);
            }
            ++v;
        }
        len*=k;
    }
}

/* find all connected components of an undirected graph and store them */
void graphe::connected_components(ivectors &components,int sg,bool skip_embedded,int *count) {
    unvisit_all_nodes(sg);
    unset_all_ancestors(sg);
    disc_time=0;
    if (count==NULL)
        components.resize(node_count());
    int c=0,i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && (!skip_embedded || !it->is_embedded()) && !it->is_visited())
            dfs(i,true,false,&components[c++],sg,skip_embedded);
    }
    if (count==NULL)
        components.resize(c);
    else
        *count=c;
}

/* find all biconnected components as vertex lists */
void graphe::biconnected_components(ivectors &components,int sg) {
    int s=sg<0?1:1+max_subgraph_index(),cnt;
    int cmp=connected_components_to_subgraphs(sg);
    vector<ipairs> blocks;
    set<int> vset;
    components.clear();
    for (;s<=cmp;++s) {
        if (subgraph_size(s)<3) {
            components.resize(components.size()+1);
            get_subgraph(s,components.back());
            continue;
        }
        blocks.clear();
        find_blocks(blocks,s);
        renumber_subgraph(s,sg);
        for (vector<ipairs>::const_iterator it=blocks.begin();it!=blocks.end();++it) {
            const ipairs &b=*it;
            vset.clear();
            for (ipairs_iter jt=b.begin();jt!=b.end();++jt) {
                vset.insert(jt->first);
                vset.insert(jt->second);
            }
            components.resize(components.size()+1);
            ivector &c=components.back();
            c.resize(vset.size());
            cnt=0;
            for (set<int>::const_iterator jt=vset.begin();jt!=vset.end();++jt) {
                c[cnt++]=*jt;
            }
        }
    }
}

/* return the number of connected components in this graph */
int graphe::connected_component_count(int sg) {
    unvisit_all_nodes(sg);
    unset_all_ancestors(sg);
    disc_time=0;
    int count=0,i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited()) {
            dfs(i,false,false,NULL,sg);
            ++count;
        }
    }
    return count;
}

/* get the connected components and set them as subgraphs with indices 1,2,.. */
int graphe::connected_components_to_subgraphs(int sg) {
    if (sg<0)
        unset_subgraphs();
    unvisit_all_nodes();
    unset_all_ancestors();
    disc_time=0;
    int i=0;
    int s=sg<0?0:max_subgraph_index();
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited()) {
            dfs(i,true,false,NULL,sg);
            set_subgraph(disc_nodes,++s);
        }
    }
    return s;
}

void graphe::strongconnect_dfs(ivectors &components,bvector &onstack,int i,int sg) {
    vertex &v=node(i);
    v.set_visited(true);
    v.set_disc(disc_time++);
    v.set_low(v.disc());
    node_stack.push(i);
    onstack[i]=true;
    int j;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        vertex &w=node(j=*it);
        if (sg>=0 && w.subgraph()!=sg)
            continue;
        if (!w.is_visited()) {
            strongconnect_dfs(components,onstack,j,sg);
            v.set_low(std::min(v.low(),w.low()));
        } else if (onstack[j])
            v.set_low(std::min(v.low(),w.disc()));
    }
    if (v.low()==v.disc()) {
        /* output a strongly connected component */
        components.resize(components.size()+1);
        ivector &component=components.back();
        do {
            j=node_stack.top();
            node_stack.pop();
            component.push_back(j);
            onstack[j]=false;
        } while (j!=i);
    }
}

/* find all strongly connected components in directed graph using Tarjan's algorithm */
void graphe::strongly_connected_components(ivectors &components,int sg) {
    assert(node_stack.empty());
    unvisit_all_nodes(sg);
    bvector onstack(node_count(),false);
    disc_time=0;
    int i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited())
            strongconnect_dfs(components,onstack,i,sg);
    }
}

/* create the condensation of tis graph by contracting strongly connected components */
void graphe::condensation(graphe &G) {
    assert(is_directed());
    ivectors comp;
    strongly_connected_components(comp);
    int n=comp.size(),dir;
    G.clear();
    G.set_directed(true);
    if (G.supports_attributes()) {
        vecteur labels;
        G.make_default_labels(labels,n);
        G.add_nodes(labels);
    } else G.add_nodes(n);
    for (ivectors_iter it=comp.begin();it!=comp.end();++it) {
        for (ivectors_iter jt=it+1;jt!=comp.end();++jt) {
            dir=0;
            for (ivector_iter kt=it->begin();kt!=it->end() && dir==0;++kt) {
                for (ivector_iter lt=jt->begin();lt!=jt->end() && dir==0;++lt) {
                    if (has_edge(*kt,*lt)) dir=1;
                    else if (has_edge(*lt,*kt)) dir=-1;
                }
            }
            switch (dir) {
            case 1:
                G.add_edge(it-comp.begin(),jt-comp.begin());
                break;
            case -1:
                G.add_edge(jt-comp.begin(),it-comp.begin());
                break;
            default:
                break;
            }
        }
    }
}

/* return true iff the connected graph is not biconnected (i.e. has an articulation point) */
bool graphe::has_cut_vertex(int sg,int i) {
    vertex &v=node(i);
    if (i==0) {
        unvisit_all_nodes();
        unset_all_ancestors();
        disc_time=0;
    }
    if (sg>=0 && v.subgraph()!=sg)
        return i==node_count()-1?false:has_cut_vertex(sg,i+1);
    v.set_visited(true);
    ++disc_time;
    v.set_disc(disc_time);
    v.set_low(disc_time);
    int children=0;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        int j=*it;
        vertex &w=node(j);
        if (sg>=0 && w.subgraph()!=sg)
            continue;
        if (!w.is_visited()) {
            ++children;
            w.set_ancestor(i);
            if (has_cut_vertex(sg,j))
                return true;
            if (v.ancestor()<0) {
                if (children==2)
                    return true;
            } else {
                v.set_low(std::min(v.low(),w.low()));
                if (w.low()>=v.disc())
                    return true;
            }
        } else if (j!=v.ancestor() && w.disc()<v.disc())
            v.set_low(std::min(v.low(),w.disc()));
    }
    return false;
}

void graphe::find_cut_vertices_dfs(int i,std::set<int> &ap,int sg) {
    vertex &v=node(i);
    v.set_visited(true);
    ++disc_time;
    v.set_disc(disc_time);
    v.set_low(disc_time);
    int children=0;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        int j=*it;
        vertex &w=node(j);
        if (sg>=0 && w.subgraph()!=sg)
            continue;
        if (!w.is_visited()) {
            ++children;
            w.set_ancestor(i);
            find_cut_vertices_dfs(j,ap,sg);
            if (v.ancestor()<0) {
                if (children==2)
                    ap.insert(i);
            } else {
                v.set_low(std::min(v.low(),w.low()));
                if (w.low()>=v.disc())
                    ap.insert(i);
            }
        } else if (j!=v.ancestor() && w.disc()<v.disc())
            v.set_low(std::min(v.low(),w.disc()));
    }
}

/* return list of cut vertices obtained by using depth-first search, time complexity O(n+m) */
void graphe::find_cut_vertices(ivector &articulation_points,int sg) {
    unvisit_all_nodes();
    unset_all_ancestors();
    disc_time=0;
    int i=0;
    std::set<int> ap;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited())
            find_cut_vertices_dfs(i,ap,sg);
    }
    std::set<int>::const_iterator jt=ap.begin(),jtend=ap.end();
    for (;jt!=jtend;++jt) {
        articulation_points.push_back(*jt);
    }
}

void graphe::find_blocks_dfs(int i,vector<ipairs> &blocks,int sg) {
    ++disc_time;
    vertex &v=node(i);
    v.set_disc(disc_time);
    v.set_low(disc_time);
    v.set_visited(true);
    int j;
    ipair edge;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        j=*it;
        vertex &w=node(j);
        if (sg>=0 && w.subgraph()!=sg)
            continue;
        edge=make_pair(i<j?i:j,i<j?j:i);
        if (!w.is_visited()) {
            w.set_ancestor(i);
            edge_stack.push(edge);
            find_blocks_dfs(j,blocks,sg);
            v.set_low(std::min(v.low(),w.low()));
            if (w.low()>=v.disc()) {
                /* output biconnected component to 'blocks' */
                blocks.resize(blocks.size()+1);
                ipairs &block=blocks.back();
                do {
                    block.push_back(edge_stack.top());
                    edge_stack.pop();
                } while (block.back()!=edge);
            }
        } else if (w.disc()<v.disc() && j!=v.ancestor()) {
            edge_stack.push(edge);
            v.set_low(std::min(v.low(),w.disc()));
        }
    }
}

/* create list of all biconnected components (as lists of edges) of the (sub)graph */
void graphe::find_blocks(vector<ipairs> &blocks,int sg) {
    assert(edge_stack.empty());
    unvisit_all_nodes(sg);
    unset_all_ancestors(sg);
    disc_time=0;
    int i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited())
            find_blocks_dfs(i,blocks,sg);
    }
    while (!edge_stack.empty()) edge_stack.pop();
}

/* create list of simple ears in the graph, which must be biconnected and not a cycle */
void graphe::find_ears(ivectors &ears,int sg) {
    assert(!is_directed() && node_stack.empty());
    ears.clear();
    int n=node_count(),root,d,ad,maxd=-1,i,k,p;
    for (i=0;i<n;++i) {
        const vertex &v=node(i);
        if ((sg<0 || v.subgraph()==sg) &&
                (d=v.neighbors().size())>maxd) {
            maxd=d;
            root=i;
        }
    }
    if (maxd<3)
        return;
    dfs(root,true,true,NULL,sg);
    bool rec=false;
    ivector_iter it=disc_nodes.begin(),itend=disc_nodes.end();
    for (;it!=itend;++it) {
        const vertex &v=node(i=*it);
        if (degree(i,sg)==2) {
            assert((p=v.ancestor())>=0); // v is not the root
            ad=degree(p,sg);
            if (multiedges(make_pair(v.ancestor(),i))>0)
                continue;
            if (ad>2) {
                /* start an ear */
                assert(!rec && node_stack.empty());
                node_stack.push(v.ancestor());
                rec=true;
            }
            if (!rec)
                continue;
            node_stack.push(i);
            k=-1;
            for (ivector_iter nt=v.neighbors().begin();nt!=v.neighbors().end();++nt) {
                if ((sg<0 || node(*nt).subgraph()==sg) && *nt!=v.ancestor())
                    k=*nt;
            }
            if (k<0 || multiedges(make_pair(i,k))>0) { // dead end
                clear_node_stack();
                rec=false;
            } else if (degree(k,sg)>2) {
                /* stop the ear */
                node_stack.push(k);
                ivector ear;
                ear.reserve(node_stack.size());
                while (!node_stack.empty()) {
                    ear.push_back(node_stack.top());
                    node_stack.pop();
                }
                if (ear.front()!=ear.back() && !has_edge(ear.front(),ear.back()))
                    ears.push_back(ear);
                rec=false;
            }
        }
    }
}

void graphe::find_bridges_dfs(int i,ipairs &B,int sg) {
    vertex &v=node(i);
    v.set_visited(true);
    ++disc_time;
    v.set_disc(disc_time);
    v.set_low(disc_time);
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        int j=*it;
        vertex &w=node(j);
        if (sg>=0 && w.subgraph()!=sg)
            continue;
        if (!w.is_visited()) {
            w.set_ancestor(i);
            find_bridges_dfs(j,B,sg);
            v.set_low(std::min(v.low(),w.low()));
            if (w.low()>v.disc())
                B.push_back(make_pair(i<j?i:j,i<j?j:i));
        } else if (j!=v.ancestor() && w.disc()<v.disc())
            v.set_low(std::min(v.low(),w.disc()));
    }
}

/* create list B of all bridges in an undirected graph */
void graphe::find_bridges(ipairs &B,int sg) {
    assert(!is_directed());
    unvisit_all_nodes(sg);
    unset_all_ancestors(sg);
    disc_time=0;
    int i=0;
    B.clear();
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited())
            find_bridges_dfs(i,B,sg);
    }
}

/* Fleury's algorithm for eulerian trail, time complexity O(m^2) */
bool graphe::fleury(int start,ivector &path) {
    /* assuming that the graph is undirected and (semi-)eulerian */
    int m=edge_count(),i=start,j;
    path.resize(m+1);
    path.back()=i;
    ipairs bridges;
    ipair edge;
    while (m>0) {
        vertex &v=node(i);
        find_bridges(bridges);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            j=*it;
            edge=make_pair(i<j?i:j,i<j?j:i);
            if (find(bridges.begin(),bridges.end(),edge)==bridges.end())
                break;
        }
        remove_edge(edge);
        path[--m]=j;
        i=j;
    }
    return true;
}

/* Hierholzer algorithm for eulerian cycle, time complexity O(m) */
void graphe::hierholzer(ivector &path) {
    /* assuming that the graph is eulerian */
    int i,j,pos=path.size();
    ivector_iter it;
    while (pos>0) {
        --pos;
        /* find a closed trail starting in path[pos] */
        j=i=path[pos];
        do {
            vertex &v=node(j);
            for (it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                if (!is_edge_visited(j,*it))
                    break;
            }
            if (it==v.neighbors().end())
                break;
            path.insert(path.begin()+(++pos),*it);
            set_edge_visited(j,*it);
            j=*it;
        } while (i!=j);
    }
}

/* return the node in which an eulerian trail may begin, or -1 if the graph has no such trail */
int graphe::eulerian_trail_start(bool &iscycle) const {
    /* assuming that the graph is connected */
    int n=node_count(),count=0,start_node=-1,d;
    for (int i=0;i<n;++i) {
        d=degree(i);
        if (d%2!=0) {
            ++count;
            start_node=i;
        } else if (start_node<0 && d%2==0 && d>0)
            start_node=i;
    }
    if (count!=0 && count!=2)
        return -1;
    iscycle=count==0;
    return start_node;
}

/* attempt to find an eulerian trail in this graph, return true iff one exists */
bool graphe::find_eulerian_trail(ivector &path) {
    int start;
    bool iscycle;
    ivectors components;
    connected_components(components);
    int ic=0;
    for (ivectors_iter it=components.begin();it!=components.end();++it) {
        if (it->size()>1) ++ic;
    }
    if (ic!=1 || (start=eulerian_trail_start(iscycle))<0)
        return false;
    /* use Hierholzer's algorithm */
    assert(visited_edges.empty());
    path.clear();
    path.reserve(edge_count()+1);
    if (!iscycle) {
        /* find a path between two odd-degree nodes and delete it */
        dfs(start,false);
        int j,k=-1;
        for (int i=node_count();i-->0;) {
            if (i!=start && degree(i)%2!=0) {
                k=i;
                break;
            }
        }
        assert(k>=0);
        while (k!=start) {
            path.push_back(k);
            const vertex &v=node(k);
            j=v.ancestor();
            set_edge_visited(j,k);
            k=j;
        }
    }
    path.push_back(start);
    hierholzer(path);
    unvisit_all_edges();
    return true;
}

int graphe::find_cycle_dfs(int i,int sg) {
    vertex &v=node(i);
    v.set_visited(true);
    node_stack.push(i);
    int j,k;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        j=*it;
        vertex &w=node(j);
        if (sg>=0 && w.subgraph()!=sg)
            continue;
        if (w.is_visited()) {
            if (v.ancestor()!=j)
                return j;
            continue;
        }
        w.set_ancestor(i);
        k=find_cycle_dfs(j,sg);
        if (k>=0)
            return k;
    }
    node_stack.pop();
    return -1;
}

bool graphe::find_cycle(ivector &cycle,int sg) {
    assert(node_stack.empty());
    if (is_null())
        return false;
    int n=node_count(),initv,i;
    unvisit_all_nodes(sg);
    unset_all_ancestors(sg);
    cycle.clear();
    cycle.reserve(n);
    if (sg>=0) {
        for (initv=0;initv<n && node(initv).subgraph()!=sg;++initv);
        if (initv==n)
            return false;
    } else initv=rand_integer(n);
    if ((i=find_cycle_dfs(initv,sg))>=0) {
        do {
            cycle.push_back(node_stack.top());
            node_stack.pop();
        } while (cycle.back()!=i);
        clear_node_stack();
        return true;
    }
    return false;
}

bool graphe::find_path_dfs(int dest,int i,int sg,bool skip_embedded) {
    node_stack.push(i);
    if (i==dest)
        return true;
    vertex &v=node(i);
    v.set_visited(true);
    int j;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        j=*it;
        vertex &w=node(j);
        if ((skip_embedded && w.is_embedded()) || w.subgraph()!=sg || w.is_visited())
            continue;
        if (find_path_dfs(dest,j,sg,skip_embedded))
            return true;
    }
    node_stack.pop();
    return false;
}

/* return a path from i-th to j-th vertex using DFS (graph is assumed to be connected) */
bool graphe::find_path(int i,int j,ivector &path,int sg,bool skip_embedded) {
    assert(node_stack.empty());
    unvisit_all_nodes(sg);
    if (find_path_dfs(i,j,sg,skip_embedded)) {
        path.resize(node_stack.size());
        while (!node_stack.empty()) {
            path[node_stack.size()-1]=node_stack.top();
            node_stack.pop();
        }
        return true;
    }
    clear_node_stack();
    return false;
}

/* clear the graph, deleting all nodes and vertices */
void graphe::clear() {
    unmark_all_nodes();
    nodes.clear();
}

/* return true iff the given face contains the edge {i,j} */
int graphe::face_has_edge(const ivector &face,int i,int j) {
    int v,w,k;
    for (ivector_iter it=face.begin();it!=face.end();++it) {
        v=*it;
        k=(it-face.begin()+1)%face.size();
        w=face[k];
        if ((i==v && j==w) || (i==w && j==v))
            return k;
    }
    return -1;
}

/* set 'embedded' field of all given nodes to 'true' */
void graphe::set_nodes_embedded(const ivector &v,bool yes) {
    for (ivector_iter it=v.begin();it!=v.end();++it) {
        node(*it).set_embedded(yes);
    }
}

/* set 'embedded' filed for all nodes to 'false' */
void graphe::unembed_all_nodes() {
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->set_embedded(false);
    }
}

/* finds planar embedding of a biconnected graph as a list of faces,
* returns true iff the graph is planar */
bool graphe::demoucron(ivectors &faces,int sg) {
    ivector cycle,path,updated_face;
    ivectors bridges,components(node_count());
    ipairs admissible_faces;
    int i,j,k,n,f,s=sg,bc=max_subgraph_index(),cnt;
    ivectors_iter ft,ftend;
    ivector_iter ct;
    std::set<int> contact_nodes;
    unembed_all_nodes();
    save_subgraphs();
    /* adding two initial faces, obtained by finding a cycle in graph */
    find_cycle(cycle,s);
    set_nodes_embedded(cycle);
    faces.clear();
    faces.push_back(ivector(cycle.begin(),cycle.end()));
    faces.push_back(ivector(cycle.rbegin(),cycle.rend()));
    while (true) {
        /* find bridges */
        node_iter nt=nodes.begin(),ntend=nodes.end();
        ivector_iter jt,jtend;
        for (i=0;nt!=ntend;++nt,++i) {
            if (!nt->is_embedded() || (s>=0 && nt->subgraph()!=s))
                continue;
            for (jt=nt->neighbors().begin(),jtend=nt->neighbors().end();jt!=jtend;++jt) {
                const vertex &w=node(j=*jt);
                if (!w.is_embedded() || (w.subgraph()==s && j<i))
                    continue;
                for (ft=faces.begin(),ftend=faces.end();ft!=ftend;++ft) {
                    if (face_has_edge(*ft,i,j)>=0)
                        break;
                }
                if (ft==faces.end()) {
                    /* create a bridge consisting only of edge (i,j) */
                    ivector bridge(4);
                    bridge.front()=++bc; // subgraph index
                    bridge[1]=0; // bridge size
                    bridge[2]=i; // contact vertex 1
                    bridge[3]=j; // contact vertex 2
                    bridges.push_back(bridge);
                }
            }
        }
        n=bridges.size();
        connected_components(components,s,true,&cnt);
        bridges.resize(n+cnt);
        for (k=0;k<cnt;++k) {
            ivector comp=components[k];
            ivector &bridge=bridges[n++];
            bridge.clear();
            bridge.push_back(++bc);
            bridge.push_back(comp.size());
            contact_nodes.clear();
            for (ivector_iter it=comp.begin();it!=comp.end();++it) {
                vertex &v=node(*it);
                v.set_subgraph(bc);
                for (ivector_iter jt=v.neighbors().begin();jt!=v.neighbors().end();++jt) {
                    if (node(*jt).is_embedded())
                        contact_nodes.insert(*jt);
                }
            }
            bridge.resize(2+contact_nodes.size());
            i=2;
            for (std::set<int>::const_iterator it=contact_nodes.begin();it!=contact_nodes.end();++it) {
                bridge[i++]=*it;
            }
        }
        if (bridges.empty())
            break;
        /* for each bridge, find all faces in which it can be drawn (admissible faces) */
        admissible_faces.resize(bridges.size());
        for (ivectors_iter it=bridges.begin();it!=bridges.end();++it) {
            ipair &admissible=admissible_faces[it-bridges.begin()];
            admissible.first=0;
            for (ft=faces.begin(),ftend=faces.end();ft!=ftend;++ft) {
                for (ct=it->begin()+2;ct!=it->end();++ct) {
                    if (find(ft->begin(),ft->end(),*ct)==ft->end())
                        break;
                }
                if (ct==it->end()) {
                    admissible.first++;
                    admissible.second=ft-faces.begin();
                }
            }
            if (admissible.first==0) { // the graph is not planar
                restore_subgraphs();
                return false;
            }
        }
        /* select the first bridge with the smallest number of admissible faces */
        n=rand_max2;
        for (ipairs_iter it=admissible_faces.begin();it!=admissible_faces.end();++it) {
            if (it->first<n) {
                k=it-admissible_faces.begin();
                if ((n=it->first)==1)
                    break;
            }
        }
        ivector &bridge=bridges[k];
        /* draw a path between the first two connecting vertices to the admissible face */
        s=bridge.front(); // subgraph index
        if ((n=bridge[1])>0) {
            i=first_neighbor_from_subgraph(node(bridge[2]),s);
            j=first_neighbor_from_subgraph(node(bridge[3]),s);
            assert(i>=0 && j>=0);
            find_path(i,j,path,s,true);
            set_nodes_embedded(path);
        } else path.clear();
        f=admissible_faces[k].second;
        faces.resize(faces.size()+1);
        ivector &face=faces[f],&new_face=faces.back();
        i=find(face.begin(),face.end(),bridge[2])-face.begin();
        j=find(face.begin(),face.end(),bridge[3])-face.begin();
        arc_path(i,j,face,updated_face);
        arc_path(j,i,face,new_face);
        if (path.size()>0) {
            updated_face.insert(updated_face.end(),path.begin(),path.end());
            new_face.insert(new_face.end(),path.rbegin(),path.rend());
        }
        face=updated_face;
        bridges.erase(bridges.begin()+k); // we're done with this bridge
    }
    restore_subgraphs();
    return true;
}

/* return the common element of two sorted lists of integers if there is one, else return -1 */
int graphe::common_element(const ivector &v1,const ivector &v2,int offset) {
    ivector_iter it1=v1.begin()+offset,it2=v2.begin()+offset;
    if (it1==v1.end() || it2==v2.end())
        return -1;
    while (it1!=v1.end() && it2!=v2.end()) {
        if (*it1==*it2)
            return *it1;
        while (*it1<*it2 && it1!=v1.end()) ++it1;
        while (*it2<*it1 && it2!=v2.end()) ++it2;
    }
    return -1;
}

/* recursively build up a tree of blocks */
void graphe::build_block_tree(int i,ivectors &blocks) {
    int n=blocks.size(),c;
    /* block structure: [is_visited (0 or 1), parent_index, connecting_vertex, ap_1, ap_2, ..., ap_n] */
    ivector &parent=blocks[i];
    parent.front()=1; // mark block as visited
    for (int j=0;j<n;++j) {
        /* find children */
        ivector &block=blocks[j];
        if (i==j || block.front()==1 || (c=common_element(parent,block,3))==-1)
            continue;
        /* block is a child connected in the vertex c */
        block.erase(find(block.begin()+3,block.end(),c));
        block[1]=i;
        block[2]=c;
        build_block_tree(j,blocks);
    }
}

/* choose a face to embed among those containing the vertex v */
int graphe::choose_embedding_face(const ivectors &faces,int v) {
    ivector candidates;
    for (ivectors_iter ft=faces.begin();ft!=faces.end();++ft) {
        if (find(ft->begin(),ft->end(),v)!=ft->end())
            candidates.push_back(ft-faces.begin());
    }
    return candidates[rand_integer(candidates.size())];
}

/* embed children blocks to its parent, for each embedding add one temporary edge */
void graphe::embed_children_blocks(int i,ivectors &block_tree,vector<ivectors> &blocks_faces) {
    /* find all children of the i-th block */
    ivector children;
    for (ivectors_iter it=block_tree.begin();it!=block_tree.end();++it) {
        if (it->at(1)==i)
            children.push_back(it-block_tree.begin());
    }
    if (children.empty()) // this is a leaf node in the tree
        return;
    /* embed each child to the largest available face of the parent */
    ivectors &parent_faces=blocks_faces[i];
    int c,pf,cf,k,n;
    ivector::iterator vt;
    ivector::reverse_iterator wt;
    for (ivector_iter it=children.begin();it!=children.end();++it) {
        embed_children_blocks(*it,block_tree,blocks_faces);
        ivectors &child_faces=blocks_faces[*it];
        c=block_tree[*it][2]; // articulation vertex connecting the child with its parent
        /* find the suitable parent and child faces that contain c */
        pf=choose_embedding_face(parent_faces,c);
        cf=choose_embedding_face(child_faces,c);
        assert(pf>=0 && cf>=0);
        n=parent_faces.size();
        parent_faces.resize(n+child_faces.size());
        /* Add all child faces different than child_face to parent_faces,
     * add a new triangular face obtained by connecting neighbors of c in child and parent,
     * and modify the parent_face by adding child_face edges and removing two edges. */
        k=n;
        for (ivectors::iterator ft=child_faces.begin();ft!=child_faces.end();++ft) {
            if (int(ft-child_faces.begin())==cf)
                continue;
            parent_faces[k++].swap(*ft);
        }
        ivector &parent_face=parent_faces[pf],&child_face=child_faces[cf],&new_face=parent_faces.back();
        /* get v, w as neighbors of c in parent resp. child face */
        vt=find(parent_face.begin(),parent_face.end(),c)+1;
        if (vt==parent_face.end())
            vt=parent_face.begin();
        wt=find(child_face.rbegin(),child_face.rend(),c)+1;
        if (wt==child_face.rend())
            wt=child_face.rbegin();
        assert(*vt!=*wt);
        add_temporary_edge(*vt,*wt);
        /* construct the new face */
        new_face.resize(3);
        new_face[0]=*wt;
        new_face[1]=c;
        new_face[2]=*vt;
        /* modify the parent_face */
        ivector path(child_face.size()-1);
        k=0;
        while (*wt!=c) {
            path[k++]=*(wt++);
            if (wt==child_face.rend())
                wt=child_face.rbegin();
        }
        std::reverse(path.begin(),path.end());
        parent_face.insert(vt,path.begin(),path.end());
    }
}

/* return the index of the largest face of the given embedding */
int graphe::choose_outer_face(const ivectors &faces) {
    int f,fsize,maxsize=0;
    for (ivectors_iter it=faces.begin();it!=faces.end();++it) {
        if ((fsize=it->size())>maxsize) {
            maxsize=fsize;
            f=it-faces.begin();
        }
    }
    return f;
}

/* finds planar embedding of a connected graph as a list of faces,
* returns the index of the outer face or -1 if the graph is not planar */
int graphe::planar_embedding(ivectors &faces) {
    /* split graph to blocks */
    vector<ipairs> blocks;
    find_blocks(blocks);
    if (blocks.size()==1)
        return demoucron(faces)?choose_outer_face(faces):-1;
    /* there exist at least one articulation point */
    unset_subgraphs();
    int sg=0;
    vector<ivectors> blocks_faces(blocks.size());
    ivector cv;
    find_cut_vertices(cv);
    int i=0,nf,ssz;
    for (vector<ipairs>::const_iterator it=blocks.begin();it!=blocks.end();++it) {
        /* test each block separately */
        set_subgraph(*it,++sg);
        ivectors &block_faces=blocks_faces[i++];
        if ((ssz=subgraph_size(sg))>2) {
            /* block has three or more vertices */
            if (int(it->size())+6>3*ssz || !demoucron(block_faces,sg))
                return -1;
        } else {
            /* block contains only one edge (a bridge)  */
            assert(it->size()==1);
            ivector bface(2);
            bface.front()=it->front().first;
            bface.back()=it->front().second;
            block_faces.push_back(bface);
        }
        /* push back a vector of articulation points which belong to this component */
        nf=block_faces.size();
        block_faces.resize(block_faces.size()+1);
        ivector &face_cv=block_faces.back();
        for (int k=0;k<nf;++k) {
            ivector &face=block_faces[k];
            for (ivector_iter cvit=cv.begin();cvit!=cv.end();++cvit) {
                if (find(face.begin(),face.end(),*cvit)!=face.end())
                    face_cv.push_back(*cvit);
            }
        }
        assert(!face_cv.empty());
        std::sort(face_cv.begin(),face_cv.end());
    }
    /* Graph is planar and we have a list of faces for each block.
     * Now blocks are embedded into each other, starting from peripheral blocks,
     * by adding temporary edges. */
    /* make a tree of blocks */
    ivectors block_tree(blocks.size());
    for (ivectors::iterator it=block_tree.begin();it!=block_tree.end();++it) {
        it->push_back(0);
        it->push_back(-1);
        it->push_back(-1);
        ivectors &fv=blocks_faces[it-block_tree.begin()];
        it->insert(it->end(),fv.back().begin(),fv.back().end());
    }
    build_block_tree(0,block_tree);
    for (vector<ivectors>::iterator it=blocks_faces.begin();it!=blocks_faces.end();++it) {
        it->pop_back();
    }
    /* embed all blocks to the root block by climbing up the tree recursively */
    for (ivectors_iter bt=block_tree.begin();bt!=block_tree.end();++bt) {
        if (bt->at(1)<0) {
            // root found
            i=bt-block_tree.begin();
            embed_children_blocks(i,block_tree,blocks_faces);
            ivectors &block_faces=blocks_faces[i];
            faces.resize(block_faces.size());
            for (ivectors::iterator it=block_faces.begin();it!=block_faces.end();++it) {
                faces[it-block_faces.begin()].swap(*it);
            }
            break;
        }
    }
    return choose_outer_face(faces);
}

/* compute the periphericity of the vertices with respect to
* the outer face by applying BFS (starting in each of the outer vertices) */
void graphe::periphericity(const ivector &outer_face,ivector &p) {
    assert(node_queue.empty());
    int level,i,j;
    std::fill(p.begin(),p.end(),rand_max2);
    for (ivector_iter jt=outer_face.begin();jt!=outer_face.end();++jt) {
        p[*jt]=0;
    }
    for (ivector_iter it=outer_face.begin();it!=outer_face.end();++it) {
        unvisit_all_nodes();
        node_queue.push(*it);
        level=1;
        while (!node_queue.empty()) {
            i=node_queue.front();
            vertex &v=node(i);
            for (ivector_iter jt=v.neighbors().begin();jt!=v.neighbors().end();++jt) {
                j=*jt;
                vertex &w=node(j);
                if (w.is_visited() || p[j]==0) // skip the outer vertices
                    continue;
                if (level<p[j])
                    p[j]=level;
                node_queue.push(j);
                w.set_visited(true);
            }
            ++level;
            node_queue.pop();
        }
    }
}

/*
* TREE NODE POSITIONING ALGORITHM
*/

/* traverse the tree from top to bottom using DFS */
void graphe::walker::walk(int i,int pass,int level,double modsum) {
    vertex &v=G->node(i);
    v.set_visited(true);
    double m=0;
    if (pass==3) {
        point &p=x->at(i);
        p.resize(2);
        p[0]=prelim[i]+modsum;
        p[1]=-level*vsep;
        m=modifier[i];
    }
    int j;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        vertex &w=G->node(j=*it);
        if (w.is_visited())
            continue;
        if (pass==1) {
            w.set_ancestor(i);
            ++children[i];
        }
        walk(j,pass,level+1,modsum+m);
    }
    if (pass==1) {
        ++node_counters[level];
        depth=std::max(depth,level+1);
    } else if (pass==2) {
        if (node_counters[level]>0 && G->node(levels[level][node_counters[level]-1]).ancestor()!=v.ancestor())
            ++gap_counters[level];
        position[i]=node_counters[level]++;
        levels[level][position[i]]=i;
        if (children[i]>0) {
            gaps[i]=gap_counters[level];
            gap_counters[level]=0;
        }
    }
}

/* set prelim and modifier for each node in V (they are on the same level) */
void graphe::walker::process_level(int i) {
    ivector &L=levels[i];
    int lastp=G->node(L.front()).ancestor(),p,n=0,m=placed.size(),j;
    double ppos=0,xpos=0,dist,shift=0,ssep=0,step=1,min_dist;
    if (m>0) {
        j=placed.front();
        xpos=prelim[j]-position[j]*hsep;
    }
    for (ivector_iter it=L.begin();it!=L.end();++it) {
        const vertex &v=G->node(*it);
        p=v.ancestor();
        if (p!=lastp) {
            prelim[lastp]=ppos/n;
            placed.push(lastp);
            lastp=p;
            ppos=0;
            n=0;
            xpos+=ssep;
        }
        if (children[*it]==0) {
            prelim[*it]=xpos;
        } else {
            prelim[*it]+=(modifier[*it]=shift);
            xpos=prelim[*it];
            step=1;
            ssep=0;
            placed.pop();
            --m;
            if (m>0) {
                vertex &w=G->node(j=placed.front());
                dist=prelim[j]+shift-xpos;
                min_dist=(position[j]-position[*it])*hsep;
                if (dist<min_dist)
                    shift+=min_dist-dist;
                else if (v.ancestor()==w.ancestor())
                    step=dist/min_dist;
                else
                    ssep=(dist-min_dist)/gaps[j];
            }
        }
        ppos+=prelim[*it];
        ++n;
        xpos+=hsep*step;
    }
    assert(m==0);
    prelim[lastp]=ppos/n;
    placed.push(lastp);
}

/* node positioning procedure */
void graphe::walker::positioning(int apex) {
    G->unset_all_ancestors();
    G->unvisit_all_nodes();
    walk(apex,1); // first walk: determine tree depth and set node ancestors
    /* allocate memory for level lists */
    levels.resize(depth);
    gap_counters.resize(depth,0);
    for (int i=0;i<depth;++i) {
        levels[i].resize(node_counters[i]);
        node_counters[i]=0;
    }
    G->unvisit_all_nodes();
    walk(apex,2); // second walk: create levels
    /* do node positioning for each level (except the top one) in a single sweep */
    for (int level=depth;level-->1;) {
        process_level(level);
    }
    G->unvisit_all_nodes();
    walk(apex,3); // third walk: sum up the modifiers (i.e. move subtrees)
}

graphe::walker::walker(graphe *gr,layout *ly,double hs,double vs) {
    G=gr;
    x=ly;
    hsep=hs;
    vsep=vs;
    depth=0;
    int n=G->node_count();
    node_counters.resize(n,0);
    position.resize(n,0);
    prelim.resize(n,0);
    modifier.resize(n,0);
    gaps.resize(n,0);
    children.resize(n,0);
}

/*
* END OF TREE NODE POSITIONING ALGORITHM
*/

/* calculate node positions and store them to the specified layout */
void graphe::make_tree_layout(layout &x,double sep,int apex) {
    int n=node_count();
    x.resize(n);
    if (n==0) return;
    /* layout the tree */
    double hsep=sep,vsep=sep*1.5;
    walker P(this,&x,hsep,vsep);
    P.positioning(apex);
}

/* compute the number of rooted trees on 1,..,n vertices, output as t[1],..,t[n] */
void graphe::compute_rutcounts(int n,vecteur &t) {
    t.resize(n+1);
    t[1]=1;
    gen sum,td;
    int nlast=1,i;
    while (n>nlast) {
        sum=0;
        for (int d=1;d<=nlast;++d) {
            i=nlast+1;
            td=t[d]*gen(d);
            for (int j=1;j<=nlast;++j) {
                i-=d;
                if (i<=0) break;
                sum+=t[i]*td;
            }
        }
        ++nlast;
        t[nlast]=sum/gen(nlast-1);
    }
}

/* Wilf's RANRUT algorithm - the original implementation: creates a random
 * rooted tree on n vertices, the numbers t[k] may be passed as pt if already computed */
void graphe::ranrut(int n,ivector &tree,const vecteur &pt) {
    int l=0,d,j,i,k=n,is1=0,is2=0,m,ll,ls;
    tree.resize(n+1);
    ipairs jd(n+1);
    gen td,z;
    if (k<3) {
        tree[1]=0;
        if (k>1) tree[2]=1;
        return;
    }
    vecteur t;
    if (int(pt.size())>=n+1)
        t=vecteur(pt.begin(),pt.begin()+n+1);
    else
        compute_rutcounts(n,t);
label12:
    if (k<=2) goto label70;
    z=t[k]*exact((k-1)*rand_uniform(),ctx);
    d=0;
label30:
    ++d;
    td=t[d]*gen(d);
    m=k;
    j=0;
label40:
    ++j;
    m-=d;
    if (m<1) goto label30;
    z-=t[m]*td;
    if (is_positive(z,ctx)) goto label40;
    jd[++is1]=make_pair(j,d);
    k=m;
    goto label12;
label70:
    tree[is2+1]=l;
    l=is2+1;
    is2+=k;
    if (k>1) tree[is2]=is2-1;
label80:
    k=jd[is1].second;
    if (k==0) goto label90;
    jd[is1].second=0;
    goto label12;
label90:
    j=jd[is1--].first;
    m=is2-l+1;
    ll=tree[l];
    ls=l+(j-1)*m-1;
    if (j==1) goto label105;
    for (i=l;i<=ls;++i) {
        tree[i+m]=tree[i]+m;
        if ((i-l)%m==0) tree[i+m]=ll;
    }
label105:
    is2=ls+m;
    if (is2==n) return;
    l=ll;
    goto label80;
}

void graphe::ranrut_forest(int m,ivectors &trees,const vecteur &alpha,const vecteur &a) {
    if (m==0)
        return;
    gen z=alpha[m]*exact(m*rand_uniform(),ctx),alphad;
    int d=0,i,j;
label30:
    ++d;
    alphad=a[d]*gen(d);
    i=m;
    j=0;
label40:
    ++j;
    i-=d;
    if (i<0) goto label30;
    z-=alpha[i]*alphad;
    if (is_positive(z,ctx)) goto label40;
    ranrut_forest(m-j*d,trees,alpha,a);
    ivector tree;
    ranrut(d,tree,a);
    tree[0]=j; // number of copies
    trees.push_back(tree);
}

/* insert a tree rooted at vertex 'root' */
void graphe::insert_tree(const ivector &tree,int root) {
    for (ivector_iter it=tree.begin()+2;it!=tree.end();++it) {
        add_edge(it-tree.begin()-1+root,*it-1+root);
    }
}

/* create random tree on n vertices with degree not larger than maxd */
void graphe::make_random_tree(int maxd) {
    /* add one edge at a time to the tree */
    int ofs=array_start(ctx),n=node_count();
    ivector src,labels=vecteur_2_vector_int(*_randperm(n,ctx)._VECTptr);
    for (ivector::iterator it=labels.begin();it!=labels.end();++it) *it-=ofs;
    src.push_back(labels.back());
    labels.pop_back();
    int v,w;
    while (!labels.empty()) {
        v=src[rand_integer(src.size())];
        w=labels.back();
        labels.pop_back();
        add_edge(v,w);
        src.push_back(w);
        if (degree(node_index(v))==maxd) {
            ivector::iterator it=find(src.begin(),src.end(),v);
            assert(it!=src.end());
            src.erase(it);
        }
    }
}

/* create a rooted tree on vertex set V uniformly at random (|V|<=500) */
void graphe::make_random_rooted_tree() {
    ivector tree;
    ranrut(node_count(),tree,vecteur(0));
    insert_tree(tree,0);
}

/* create a free tree on vertex set V uniformly at random */
void graphe::make_random_free_tree() {
    int n=node_count();
    vecteur a;
    compute_rutcounts(n,a);
    /* the following is a correction of Wilf's algorithm:
     * a_n in step (T1) on page 207 should be T_n, the number of unrooted trees
     * on n vertices.
     * This number is correctly computed from a_1,a_2,..,a_n using the formula in:
     * Otter, "Number of trees" (1948), pp. 589 (http://users.math.msu.edu/users/magyar/Math482/Otter-Trees.pdf)
     */
    gen Tn=a[n];
    for (int i=1;i<n;++i) {
        Tn-=a[i]*a[n-i]/2;
    }
    if (n%2==0)
        Tn+=a[n/2]/2;
    if (n%2==0 && is_strictly_greater(a[n/2]*(1+a[n/2])/(2*Tn),exact(rand_uniform(),ctx),ctx)) {
        /* the output tree will have two centroids */
        ivector tree;
        ranrut(n/2,tree,a);
        insert_tree(tree,0);
        if (is_positive(exact(rand_uniform(),ctx)*(a[n/2]+1)-1,ctx))
            ranrut(n/2,tree,a);
        insert_tree(tree,n/2);
        add_edge(0,n/2);
    } else {
        /* the output tree will have one centroid */
        vecteur alpha(n);
        alpha[0]=1;
        gen sum,alphad;
        int m=0,i;
        while (n-1>m) {
            sum=0;
            for (int d=1;2*d<=n-1;++d) {
                i=m+1;
                alphad=a[d]*gen(d);
                for (int j=1;;++j) {
                    i-=d;
                    if (i<0) break;
                    sum+=alpha[i]*alphad;
                }
            }
            ++m;
            alpha[m]=sum/gen(m);
        }
        ivectors trees;
        ivector roots;
        ranrut_forest(n-1,trees,alpha,a);
        m=1;
        for (ivectors_iter it=trees.begin();it!=trees.end();++it) {
            for (int j=0;j<it->front();++j) {
                insert_tree(*it,m);
                roots.push_back(m);
                m+=it->size()-1;
            }
        }
        for (ivector_iter it=roots.begin();it!=roots.end();++it) {
            add_edge(0,*it);
        }
    }
}

/* use DFS traversal to compute tree depth */
void graphe::tree_height_dfs(int i,int level,int &depth) {
    vertex &v=node(i);
    v.set_visited(true);
    depth=std::max(depth,level);
    int j;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        j=*it;
        vertex &w=node(j);
        if (!w.is_visited())
            tree_height_dfs(j,level+1,depth);
    }
}

/* return the height of this tree if root>=0,
 * otherwise (if root<0) return the first root for which the height is minimal */
int graphe::tree_height(int root) {
    /* assuming that this graph is a tree */
    if (node_count()==1)
        return 0;
    if (root>=0) {
        assert(root<node_count());
        unvisit_all_nodes();
        int depth=0;
        tree_height_dfs(root,0,depth);
        return depth;
    }
    /* compute the minimal tree depth (linear time):
     * https://www.geeksforgeeks.org/roots-tree-gives-minimum-height/ */
    queue<int> q;
    ivector d=degree_sequence(),res;
    int n=node_count();
    for (int i=0;i<n;++i) {
        if (d[i]==1)
            q.push(i);
    }
    while (n>2) {
        int n_pops=q.size();
        n-=n_pops;
        for (int i=0;i<n_pops;++i) {
            int t=q.front();
            q.pop();
            const ivector &ngh=node(t).neighbors();
            for (ivector_iter it=ngh.begin();it!=ngh.end();++it) {
                if (--d[*it]==1)
                    q.push(*it);
            }
        }
    }
    /* at most two vertices remain in the queue */
    while (!q.empty()) {
        res.push_back(q.front());
        q.pop();
    }
    assert(!res.empty() && res.size()<3);
    return (res.size()==1 || res.front()<res.back())?res.front():res.back();
}


/* create a random biconnected planar graph with n vertices */
void graphe::make_random_planar(double p,int connectivity) {
    /* connectivity is 0 (any), 1 (connected), 2 (biconnected) or 3 (triconnected) */
    assert(p>=0 && p<1 && connectivity>=0 && connectivity<4);
    set_directed(false);
    /* start with a random maximal planar graph with n vertices */
    int n=node_count();
    if (n==1)
        return;
    if (n==2) {
        add_edge(0,1);
        return;
    }
    ivectors faces;
    ivector face(3),old_face;
    for (int i=0;i<3;++i) {
        face[i]=i;
        add_edge(i,(i+1)%3);
    }
    faces.push_back(face);
    ivector outer_face(face.begin(),face.end());
    std::reverse(outer_face.begin(),outer_face.end());
    faces.push_back(outer_face);
    int k;
    /* create a random triangulated graph on n vertices */
    for (int i=3;i<n;++i) {
        k=rand_integer(faces.size());
        ivector &face_k=faces[k];
        old_face.assign(face_k.begin(),face_k.end());
        face[0]=i;
        for (int j=0;j<3;++j) {
            face[1]=old_face[j];
            face[2]=old_face[(j+1)%3];
            add_edge(i,old_face[j]);
            faces.push_back(face);
        }
        faces.erase(faces.begin()+k);
    }
    if (p==0)
        return;
    int mindeg=std::max(2,connectivity+1);
    ipairs E;
    get_edges_as_pairs(E);
    ivector indices=rand_permu(E.size());
    /* try to remove each edge with probability 0<p<1 */
    for (ivector_iter it=indices.begin();it!=indices.end();++it) {
        ipair &e=E[*it];
        if (degree(e.first)<mindeg || degree(e.second)<mindeg)
            continue;
        if (rand_uniform()<p)
            remove_edge(e);
        switch (connectivity) {
        case 0:
            break;
        case 1:
            if (!is_connected())
                add_edge(e);
            break;
        case 2:
            if (!is_biconnected())
                add_edge(e);
            break;
        case 3:
            if (!is_triconnected())
                add_edge(e);
            break;
        default:
            assert(false); // unreachable
        }
    }
}

/* generate random edges in this graph */
void graphe::erdos_renyi(double p) {
    int n=node_count(),m=std::floor(p),i,j;
    bool isdir=is_directed();
    if (m==0) {
        /* each edge is chosen with probability p */
        try {
            for (int k=0;k<(isdir?2:1);++k) {
                i=1; j=-1;
                while (i<n) {
                    double r=rand_uniform();
                    j+=1+std::floor(std::log(1-r)/std::log(1-p));
                    while (j>=i && i<n) {
                        j-=i;
                        ++i;
                    }
                    if (i<n) add_edge(k==0?i:j,k==0?j:i);
                }
            }
        } catch (const std::runtime_error &e) {
            if (strstr(e.what(),"user interruption")!=NULL || interrupted || ctrl_c)
                throw e;
        } 
    } else {
        int N=(n*(n-1))/(isdir?1:2),i,j,r;
        bool islarge=m>N/2;
        try {
            if (islarge) { // first generate a complete graph
                for (i=0;i<n;++i)
                    for (j=isdir?0:i+1;j<n;++j)
                        if (i!=j) add_edge(i,j);
            }
            for (int k=0;k<(islarge?N-m:m);++k) {
                while (true) {
                    r=rand_integer(N);
                    if (isdir) {
                        i=r/(n-1); j=r%(n-1);
                        if (j>=i) ++j;
                    } else {
                        i=std::floor((1.0+std::sqrt(1.0+8.0*r))/2.0);
                        j=r-i*(i-1)/2;
                    }
                    if (i>=n || j>=n) continue;
                    if (islarge && has_edge(i,j)) {
                        remove_edge(i,j);
                        break;
                    }
                    if (!islarge && !has_edge(i,j)) {
                        add_edge(i,j);
                        break;
                    }
                }
            }
        } catch (const std::runtime_error &e) {
            if (strstr(e.what(),"user interruption")!=NULL || interrupted || ctrl_c)
                throw e;
        }
    }
}

/* generate edges of this graph according to preferential attachment rule */
void graphe::preferential_attachment(int d,int o) {
    assert(!is_directed());
    int n=node_count();
    if (n<2) return;
    add_edge(0,1);
    int j,k;
    bucketsampler sampler(ivector(2,1),ctx);
    for (int i=2;i<n;++i) {
        for (int count=std::min(i,d);count-->0;) {
            do j=sampler.generate(); while (has_edge(i,j));
            add_edge(i,j);
            sampler.increment(j);
        }
        for (int count=0;count<o;++count) {
            const ivector &ngh=node(i).neighbors();
            if (ngh.size()<2) break;
            j=rand_integer(ngh.size());
            do k=rand_integer(ngh.size()); while (k==j);
            j=ngh[j];
            k=ngh[k];
            if (!has_edge(j,k)) {
                add_edge(j,k);
                sampler.increment(j);
                sampler.increment(k);
            }
        }
        sampler.insert(std::min(i,d));
    }
}

/* return true iff s is a graphic sequence */
bool graphe::is_graphic_sequence(const ivector &s_orig) {
    ivector s=s_orig;
    std::sort(s.begin(),s.end());
    std::reverse(s.begin(),s.end());
    long sum=0;
    for (ivector_iter it=s.begin();it!=s.end();++it) sum+=*it;
    if (sum%2!=0) return false;
    int k=1,m,n=s.size();
    sum=0;
    for (ivector_iter it=s.begin();it!=s.end();++it,++k) {
        sum+=*it;
        m=0;
        for (int i=k;i<n;++i) m+=std::min(s[i],k);
        if (sum>k*(k-1)+m)
            return false;
    }
    return true;
}

/* generate edges according to the given degree distribution */
void graphe::molloy_reed(const vecteur &p) {
    int n=node_count(),k,i;
    long sum;
    ivector stubs(n);
    ransampl rs(p,ctx);
    do {
        i=0; sum=0;
        bool done=false;
        while (true) {
            sum+=(stubs[i]=rs.generate());
            if (!done && i<n-1) ++i;
            else if (sum%2!=0) {
                sum-=stubs[(i=rand_integer(n))];
                done=true;
            } else break;
        }
    } while (!is_graphic_sequence(stubs));
    assert(hakimi(stubs));
    ipairs E;
    get_edges_as_pairs(E);
    int m=E.size();
    if (m<2) return;
    int iters=1+std::floor(M_LN2/std::log(double(m)/double(m-1))),tmp,at1,at2;
    ipair *e,*f;
    for (k=0;k<iters;++k) {
        at1=0;
        do {
            ++at1;
            at2=0;
            e=&E[rand_integer(m)];
            do {
                ++at2;
                f=&E[rand_integer(m)];
            } while (at2<10 && (e==f || edges_incident(*e,*f) ||
                                (has_edge(e->first,f->first) && has_edge(e->first,f->second)) ||
                                (has_edge(e->first,f->first) && has_edge(f->first,e->second)) ||
                                (has_edge(f->first,e->second) && has_edge(e->second,f->second)) ||
                                (has_edge(e->first,f->second) && has_edge(e->second,f->second))));
            if (at2<10) break;
        } while (at1<10);
        if (at1>=10) break;
        remove_edge(*e);
        remove_edge(*f);
        if (!has_edge(e->first,f->second) && !has_edge(f->first,e->second)) {
            tmp=e->second;
            e->second=f->second;
            f->second=tmp;
        } else {
            tmp=e->second;
            e->second=f->first;
            f->first=tmp;
        }
        add_edge(*e);
        add_edge(*f);
    }
}

/* create a random graph with the given degree sequence d */
void graphe::make_random_sequential(const ivector &d) {
    /* assuming that d is graphical */
    assert(int(d.size())==node_count() && d.size()>0);
    int s=0;
    for (ivector_iter it=d.begin();it!=d.end();++it) {
        s+=*it;
    }
    assert((s%2)==0);
    int m=s/2,n=d.size(),sz,cnt;
    double p,r;
    map<ipair,bool> used;
    vector<pair<ipair,double> > emap((n*(n-1))/2);
    ivector D;
    ipair e;
    do {
        D=d;
        used.clear();
        cnt=0;
        while(true) {
            p=0;
            sz=0;
            int &i=e.first,&j=e.second;
            for (i=0;i<n;++i) {
                for (j=i+1;j<n;++j) {
                    if (used[e])
                        continue;
                    p+=D[i]*D[j]*(1-d[i]*d[j]/(4*m));
                    emap[sz++]=make_pair(e,p);
                }
            }
            if (p==0)
                break;
            r=rand_uniform()*p;
            for (int k=0;k<sz;++k) {
                pair<ipair,double> &elm=emap[k];
                if (r<=elm.second) {
                    const ipair &edge=elm.first;
                    used[edge]=true;
                    ++cnt;
                    --D[edge.first];
                    --D[edge.second];
                    break;
                }
            }
        }
    } while (cnt<m);
    for (map<ipair,bool>::const_iterator it=used.begin();it!=used.end();++it) {
        if (it->second)
            add_edge(it->first);
    }
}

/* create a random bipartite graph with partition A,B */
void graphe::make_random_bipartite(int a,int b,double p) {
    assert(!is_directed() && a>0 && b>0 && p>0);
    int m=std::floor(p),ec=0,n=node_count(),M=std::min(m,a*b);
    assert(a+b==n);
    for (int i=0;i<a;++i) {
        for (int j=a;j<a+b;++j) {
            if ((m==0 && rand_uniform()<p) ||
                    (m>0 && rand_uniform()<M/double(a*b))) {
                add_edge(i,j);
                ++ec;
            }
        }
    }
    if (m>0) {
        int d=ec-M,i,j;
        if (d>0) { // remove d edges at random
            for (int k=0;k<d;++k) {
                do i=rand_integer(n); while (node(i).degree()==0);
                const vertex &v=node(i);
                j=v.neighbors().at(rand_integer(v.degree()));
                remove_edge(i,j);
            }
        } else if (d<0) { // add d edges at random
            for (int k=0;k<-d;++k) {
                do i=rand_integer(a); while (node(i).degree()>=b);
                do j=a+rand_integer(b); while (has_edge(i,j));
                add_edge(i,j);
            }
        }
    }
}

/* create a random d-regular graph with vertices fromm V */
void graphe::make_random_regular(int d,bool connected) {
    assert(!is_directed());
    ipairs E;
    int n=node_count();
    ivector prob,degrees(n);
    int prob_total,k,dd;
    double r;
    ipair edge;
    do {
        if (connected)
            make_random_tree(d);
        else {
            for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
                it->clear_neighbors();
            }
        }
        for (int i=0;i<n;++i) {
            degrees[i]=degree(i);
        }
        E.clear();
        for (int i=0;i<n;++i) {
            if ((dd=degrees[i])<d) {
                for (int j=i+1;j<n;++j) {
                    if (!nodes_are_adjacent(i,j) && degrees[j]<d)
                        E.push_back(make_pair(i,j));
                }
            }
        }
        while (!E.empty()) {
            prob.resize(E.size());
            prob_total=0;
            for (ipairs_iter it=E.begin();it!=E.end();++it) {
                prob_total+=(d-degrees[it->first])*(d-degrees[it->second]);
                prob[it-E.begin()]=prob_total;
            }
            k=0;
            if (prob_total>0) {
                r=rand_uniform()*prob_total;
                for (ivector_iter it=prob.begin();it!=prob.end();++it) {
                    if (r<=*it)
                        break;
                    ++k;
                }
            }
            edge=E[k];
            E.erase(E.begin()+k);
            ++degrees[edge.first];
            ++degrees[edge.second];
            add_edge(edge);
            for (k=E.size();k-->0;) {
                ipair &e=E[k];
                if (degrees[e.first]==d || degrees[e.second]==d)
                    E.erase(E.begin()+k);
            }
        }
    } while (is_regular(d)<0);
}

/* return -1 iff the graph is not (d-)regular, else return >=0 */
int graphe::is_regular(int d) const {
    int n=node_count();
    int deg=d;
    bool isdir=is_directed();
    for (int i=0;i<n;++i) {
        if (deg<0) {
            deg=degree(i);
            if (isdir && out_degree(i)!=in_degree(i))
                return -1;
        } else {
            if (degree(i)!=deg || (isdir && in_degree(i)!=out_degree(i)))
                return -1;
        }
    }
    return deg;
}

/* return true iff the graph is strongly regular with sig = (lambda,mu) */
bool graphe::is_strongly_regular(ipair &sig) {
    assert(!is_null());
    if (!is_regular(-1))
        return false;
    int n=node_count();
    int lambda=-1,mu=-1,m;
    for (int i=0;i<n;++i) {
        const ivector &Nv=node(i).neighbors();
        for (int j=0;j<n;++j) {
            if (i==j)
                continue;
            const ivector &Nw=node(j).neighbors();
            m=intersect_linear(Nv.begin(),Nv.end(),Nw.begin(),Nw.end());
            if (has_edge(i,j)) {
                if (lambda<0)
                    lambda=m;
                else if (lambda!=m)
                    return false;
            } else {
                if (mu<0)
                    mu=m;
                else if (mu!=m)
                    return false;
            }
        }
    }
    sig.first=lambda;
    sig.second=mu;
    return true;
}

/* return true iff this graph is equal to G */
bool graphe::is_equal(const graphe &G) const {
    if (is_directed()!=G.is_directed() || is_weighted()!=G.is_weighted())
        return false;
    vecteur V1=vertices(),V2=G.vertices();
    if (V1!=V2)
        return false;
    ipairs E1,E2;
    edgeset edg,edG;
    get_edges_as_pairs(E1);
    G.get_edges_as_pairs(E2);
    ipairs2edgeset(E1,edg);
    ipairs2edgeset(E2,edG);
    if (edg!=edG)
        return false;
    if (is_weighted()) {
        for (edgeset_iter it=edg.begin();it!=edg.end();++it) {
            if (weight(*it)!=G.weight(*it))
                return false;
        }
    }
    return true;
}

/* sort rectangles by height */
void graphe::sort_rectangles(vector<rectangle> &rectangles) {
    rectangle::comparator comp;
    std::sort(rectangles.begin(),rectangles.end(),comp);
}

/* packing rectangles (sorted by height) into an enclosing rectangle with specified dimensions,
* returns true if embedding has changed */
bool graphe::embed_rectangles(vector<rectangle> &rectangles,double maxheight) {
    bool embedding_changed=false;
    double old_x,old_y;
    dpair old_anchor;
    rectangle *old_rect;
    int lock_type;
    for (vector<rectangle>::iterator it=rectangles.begin()+1;it!=rectangles.end();++it) {
        rectangle &r=*it;
        old_x=r.x();
        old_y=r.y();
        r.set_anchor(DBL_MAX,-1);
        old_rect=NULL;
        for (vector<rectangle>::iterator jt=rectangles.begin();jt!=it;++jt) {
            rectangle &s=*jt;
            if (s.x()>r.x())
                continue;
            old_anchor.first=r.x();
            old_anchor.second=r.y();
            if (!s.is_locked_above()) {
                r.set_anchor(s.x(),s.y()+s.height());
                if (r.y()+r.height()<=maxheight && !r.intersects(rectangles.begin(),it)) {
                    s.set_locked_above(true);
                    if (old_rect!=NULL) {
                        if (lock_type==1)
                            old_rect->set_locked_above(false);
                        else if (lock_type==2)
                            old_rect->set_locked_right(false);
                    }
                    old_rect=&s;
                    lock_type=1;
                    continue;
                }
            }
            if (!s.is_locked_right()) {
                r.set_anchor(s.x()+s.width(),s.y());
                if (r.x()<old_anchor.first && r.y()+r.height()<=maxheight && !r.intersects(rectangles.begin(),it)) {
                    s.set_locked_right(true);
                    if (old_rect!=NULL) {
                        if (lock_type==1)
                            old_rect->set_locked_above(false);
                        else if (lock_type==2)
                            old_rect->set_locked_right(false);
                    }
                    old_rect=&s;
                    lock_type=2;
                    continue;
                }
            }
            r.set_anchor(old_anchor.first,old_anchor.second);
        }
        assert(r.y()>=0);
        if (old_x!=r.x() || old_y!=r.y())
            embedding_changed=true;
    }
    return embedding_changed;
}

/* pack rectangles (sorted by height) to an enclosing rectangle with minimal perimeter and wasted space */
void graphe::pack_rectangles(vector<rectangle> &rectangles) {
    /* compute total area occupied by the rectangles */
    double total_area=0;
    for (vector<rectangle>::const_iterator it=rectangles.begin();it!=rectangles.end();++it) {
        total_area+=it->width()*it->height();
    }
    /* step = the length of a negligible part of the shorter side of the smallest rectangle */
    double maxwidth=0,minwidth=DBL_MAX,minheight=rectangles.back().height();
    for (vector<rectangle>::const_iterator it=rectangles.begin();it!=rectangles.end();++it) {
        if (it->width()>maxwidth)
            maxwidth=it->width();
        if (it->width()<minwidth)
            minwidth=it->width();
    }
    double step=std::min(minwidth,minheight)*MARGIN_FACTOR;
    double bw=DBL_MAX,bh=rectangles.front().height(); // initial enclosing rectangle has an "unlimited" width
    double perim,best_perim=DBL_MAX,d;
    double xpos=0;
    dpairs best_embedding(rectangles.size());
    for (vector<rectangle>::iterator it=rectangles.begin();it!=rectangles.end();++it) {
        it->set_anchor(xpos,0);
        xpos+=it->width();
    }
    if (rectangles.size()<2)
        return;
    bool firstpass=true;
    clock_t start=clock();
    while (bw>maxwidth+step) { // loop breaks after a stacked embedding is obtained
        for (vector<rectangle>::iterator it=rectangles.begin();it!=rectangles.end();++it) {
            it->set_locked_above(false);
            it->set_locked_right(false);
        }
        if (firstpass || embed_rectangles(rectangles,bh)) {
            bw=bh=0;
            /* find the smallest enclosing rectangle containing the embedding */
            for (vector<rectangle>::const_iterator it=rectangles.begin();it!=rectangles.end();++it) {
                if ((d=it->x()+it->width())>bw)
                    bw=d;
                if ((d=it->y()+it->height())>bh)
                    bh=d;
            }
            /* find the embedding with the smallest perimeter (when scaled by the wasted ratio) */
            if ((perim=(bw+PLASTIC_NUMBER_2*bh)*std::sqrt(bw*bh/total_area))<best_perim) {
                best_perim=perim;
                for (vector<rectangle>::const_iterator it=rectangles.begin();it!=rectangles.end();++it) {
                    dpair &p=best_embedding[it-rectangles.begin()];
                    p.first=it->x();
                    p.second=it->y();
                }
            }
        }
        /* increase enclosing rectangle height (rectangles will end up stacked eventually) */
        bh+=step;
        firstpass=false;
        if (double(clock()-start)/CLOCKS_PER_SEC>1.5)
            break;
    }
    for (vector<rectangle>::iterator it=rectangles.begin();it!=rectangles.end();++it) {
        dpair &p=best_embedding[it-rectangles.begin()];
        it->set_anchor(p.first,p.second);
    }
}

/* return true iff an isomorphic copy with vertices permuted according to sigma is constructed */
bool graphe::isomorphic_copy(graphe &G,const ivector &sigma,bool strip_attributes,int sg) {
    assert(supports_attributes() || !G.supports_attributes());
    int n=sigma.size();
    G.clear();
    G.set_name(name());
    G.register_user_tags(user_tags);
    G.set_graph_attributes(attributes);
    /* add vertices */
    G.reserve_nodes(n);
    if (!G.supports_attributes())
        G.add_nodes(n);
    else {
        for (ivector_iter it=sigma.begin();it!=sigma.end();++it) {
            if (strip_attributes)
                G.add_node(node_label(*it));
            else G.add_node(node_label(*it),node(*it).attributes());
        }
    }
    if (G.node_count()!=n)
        return false;
    /* add edges */
    ipairs E;
    get_edges_as_pairs(E,sg);
    map<int,int> sigma_inv;
    /* obtain the inverse of sigma */
    int i=0;
    for (ivector_iter it=sigma.begin();it!=sigma.end();++it,++i)
        sigma_inv[*it]=i;
    ipair f;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        const ipair &e=*it;
        f=make_pair(sigma_inv[e.first],sigma_inv[e.second]);
        if (G.supports_attributes() && !strip_attributes)
            G.add_edge(f,edge_attributes(e));
        else G.add_edge(f);
        G.set_multiedge(f,multiedges(e));
    }
    return true;
}

/* return true iff the vertices were successfully relabeled with 'labels' */
bool graphe::relabel_nodes(const vecteur &labels) {
    int n=node_count();
    if (int(labels.size())<n)
        return false;
    int i=0;
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->set_label(labels[i++]);
    }
    return true;
}

/* return true iff the segments from p to p+r and from q to q+s intersect (compute the intersection point) */
bool graphe::segments_crossing(const point &p,const point &r,const point &q,const point &s,point &crossing) {
    point qq(2),pp(2),rr(2);
    copy_point(p,pp);
    copy_point(q,qq);
    copy_point(r,rr);
    subtract_point(qq,p);
    double vp1=point_vecprod2d(qq,r),vp2=point_vecprod2d(r,s);
    if (vp1==0 && vp2==0) // collinear edges
        return false;
    if (vp1!=0 && vp2==0) // parallel and non-intersecting edges
        return false;
    if (vp2!=0) {
        double t=point_vecprod2d(qq,s)/vp2,u=vp1/vp2;
        if (t>=0 && u>=0 && t<=1 && u<=1) {
            /* edges intersect each other, obtain the intersection point */
            scale_point(rr,t);
            add_point(pp,rr);
            crossing.resize(2);
            crossing.front()=pp[0];
            crossing.back()=pp[1];
            return true;
        }
    }
    /* edges are non-intersecting */
    return false;
}

/* return true iff the point p can be projected onto a segment [q,r] (also compute the projection) */
bool graphe::point2segment_projection(const point &p,const point &q,const point &r,point &proj) {
    assert(p.size()==2 && q.size()==2 && r.size()==2);
    proj.resize(2);
    point a(2),b(2),c(2);
    copy_point(p,a);
    subtract_point(a,q);
    copy_point(r,b);
    subtract_point(b,q);
    copy_point(r,c);
    subtract_point(c,p);
    double ab;
    if ((ab=point_dotprod(a,b))*point_dotprod(b,c)<=0)
        return false;
    copy_point(b,proj);
    double bn2=point_displacement(b,false);
    if (bn2==0)
        return false;
    scale_point(proj,ab/bn2);
    add_point(proj,q);
    return true;
}

/* return the value of the largest integer node label */
int graphe::largest_integer_label() const {
    assert(supports_attributes());
    int n,m=array_start(ctx)-1;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->label().is_integer() && (n=it->label().val)>m)
            m=n;
    }
    return m;
}

/* return true iff two undirected edges are incident to each other */
bool graphe::edges_incident(const ipair &e1,const ipair &e2) {
    return e1.first==e2.first || e1.first==e2.second || e1.second==e2.first || e1.second==e2.second;
}

/* return true iff two edges cross each other (also compute the crossing point) */
bool graphe::edges_crossing(const ipair &e1,const ipair &e2,const layout &x,point &crossing) const {
    if (edges_incident(e1,e2))
        return false;
    point p(2),q(2),r(2),s(2);
    copy_point(x[e1.first],p);
    copy_point(x[e2.first],q);
    copy_point(x[e1.second],r);
    copy_point(x[e2.second],s);
    subtract_point(r,p);
    subtract_point(s,q);
    return segments_crossing(p,r,q,s,crossing);
}

/* convert point to vecteur */
vecteur graphe::point2vecteur(const point &p) {
    vecteur res;
    for (point::const_iterator it=p.begin();it!=p.end();++it) {
        res.push_back(*it);
    }
    return res;
}

/* return true iff the points p and q coincide with each other within the given tolerance radius */
bool graphe::points_coincide(const point &p,const point &q,double tol) {
    assert(q.size()==p.size());
    point r(p.size());
    return point_distance(p,q,r)<=tol;
}

/*
 * IMPLEMENTATION OF RANDOM SAMPLING CLASSES
 *
 * RANSAMPL:
 * drawing samples from a given probability distribution
 * adapted from: ransampl.c by Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
 * see the original code at: apps.jcns.fz-juelich.de/ransampl
 *
 * Description:
 * Random-number sampling using the Walker-Vose alias
 * method, as described by Keith Schwarz (2011)
 */

graphe::ransampl::ransampl(const vecteur &p,GIAC_CONTEXT) {
    ctx=contextptr;
    n=p.size();
    alias.resize(n);
    prob.resize(n);
    vecteur P(n);
    ivector S(n),L(n);
    gen norm_factor=gen(n)/_sum(p,ctx);
    for (int i=0;i<n;++i) P[i]=_evalf(p[i]*norm_factor,ctx);
    int nS=0,nL=0,a,g;
    for (int i=n-1;i>=0;--i) {
        if (is_strictly_greater(1,P[i],ctx)) S[nS++]=i; else L[nL++]=i;
    }
    while (nS!=0 && nL!=0) {
        a=S[--nS];
        g=L[--nL];
        prob[a]=P[a];
        alias[a]=g;
        P[g]+=P[a]-1;
        if (is_strictly_greater(1,P[g],ctx)) S[nS++]=g; else L[nL++]=g;
    }
    while (nL!=0) prob[L[--nL]]=1;
    while (nS!=0) prob[S[--nS]]=1;
}

gen graphe::ransampl::data() const {
    vecteur ret(1,n);
    ret=mergevecteur(ret,prob);
    return mergevecteur(ret,vector_int_2_vecteur(alias));
}

int graphe::ransampl::generate() const {
    double ran1=giac_rand(ctx)/(rand_max2+1.0);
    double ran2=giac_rand(ctx)/(rand_max2+1.0);
    int i=std::floor(n*ran1);
    return is_strictly_greater(prob[i],ran2,ctx)?i:alias[i];
}

/*
 * BUCKETSAMPLER:
 * draw samples from a given dynamic probability distribution,
 * supports updating the weights
 */

graphe::bucketsampler::bucketsampler(const ivector &W,GIAC_CONTEXT) {
    ctx=contextptr;
    int n=(weights=W).size(),k,w;
    total_weight=0;
    positions.resize(n);
    for (int i=0;i<n;++i) {
        if ((w=weights[i])==0)
            continue;
        k=nearest_pow2(w);
        max_weight[k]=std::max(max_weight[k],w);
        level_weight[k]+=w;
        positions[i]=make_pair(k,levels[k].size());
        levels[k].push_back(i);
        total_weight+=w;
    }
}

int graphe::bucketsampler::generate() {
    int u=giac_rand(ctx)/(rand_max2+1.0)*double(total_weight),i;
    long cweight=0;
    for (map<int,ivector>::const_reverse_iterator it=levels.rbegin();it!=levels.rend();++it) {
        cweight+=level_weight[it->first];
        i=it->first;
        if (cweight>u)
            break;
    }
    ivector &level=levels[i];
    int lmax=max_weight[i];
    double r,frac;
    int idx_in_level,idx,sz=level.size();
    do {
        r=giac_rand(ctx)/(rand_max2+1.0)*double(sz);
        idx_in_level=std::floor(r);
        frac=r-idx_in_level;
        idx=level[idx_in_level];
    } while (weights[idx]<frac*lmax);
    return idx;
}

void graphe::bucketsampler::insert(int w) {
    int k=nearest_pow2(w);
    positions.push_back(make_pair(k,levels[k].size()));
    levels[k].push_back(weights.size());
    weights.push_back(w);
    max_weight[k]=std::max(max_weight[k],w);
    level_weight[k]+=w;
    total_weight+=w;
}

void graphe::bucketsampler::update(int i,int w) {
    int old_w=weights[i];
    ipair &pos=positions[i];
    int old_k=pos.first,index=pos.second,k=nearest_pow2(w),last;
    if (k!=old_k) {
        ivector &old_level=levels[old_k];
        ivector &level=levels[k];
        last=old_level.back();
        ipair &last_pos=positions[last];
        last_pos.second=index;
        old_level[index]=last;
        old_level.pop_back();
        pos=make_pair(k,level.size());
        level.push_back(i);
        if (w>=max_weight[old_k]) {
            int &mw=max_weight[old_k];
            mw=0;
            for (ivector_iter it=old_level.begin();it!=old_level.end();++it) {
                if (weights[*it]>mw)
                    mw=weights[*it];
            }
        }
    }
    weights[i]=w;
    level_weight[old_k]-=old_w;
    level_weight[k]+=w;
    max_weight[k]=std::max(max_weight[k],w);
    total_weight+=w-old_w;
}

/*
 * END OF RANDOM SAMPLING CLASSES
 */

/*
 * IMPLEMENTATION OF THE DISJOINT SET DATA STRUCTURE
 */

/* constructor: reserve n memory blocks for set elements */
graphe::unionfind::unionfind(int n) {
    sz=n;
    elements=new element[sz];
    clear();
}

/* clear all sets */
void graphe::unionfind::clear() {
    for (int i=0;i<sz;++i) {
        element &e=elements[i];
        e.id=e.parent=-1;
        e.rank=0;
    }
}

/* create a singleton set containing id as its
 * only element, if id is already contained in
 * one of sets then do nothing */
void graphe::unionfind::make_set(int id) {
    assert(id>=0 && id<sz);
    element &e=elements[id];
    if (e.id==id)
        return; // the element is already there
    e.id=id;
    e.parent=id;
    e.rank=1;
}

/* return the representative of the set containing id */
int graphe::unionfind::find(int id) {
    assert(id>=0 && id<sz);
    element &e=elements[id];
    if (e.parent!=id)
        e.parent=find(e.parent);
    return e.parent;
}

/* union by rank */
void graphe::unionfind::unite(int id1,int id2) {
    assert(id1>=0 && id1<sz && id2>=0 && id2<sz);
    int p1=find(id1),p2=find(id2);
    if (p1==p2)
        return;
    element &x=elements[p1],&y=elements[p2];
    if (x.rank>y.rank)
        y.parent=x.id;
    else if (x.rank<y.rank)
        x.parent=y.id;
    else {
        y.parent=x.id;
        ++x.rank;
    }
}

/* make id the representative of its set */
void graphe::unionfind::select(int id) {
    assert(id>=0 && id<sz);
    int i=id,p,np=-1;
    while (true) {
        p=elements[i].parent;
        if (np>=0)
            elements[i].parent=np;
        if (p==i)
            break;
        np=i;
        i=p;
    }
    elements[id].parent=id;
}

/*
 * END OF DISJOINT_SET CLASS
 */

/* make planar layout */
bool graphe::make_planar_layout(layout &x,double *score) {
    int n=node_count(),of,m;
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->clear_edge_faces();
    }
    ivectors faces;
    /* create the faces (adding temporary edges if necessary),
 * return the index of the outer face */
    int ec=edge_count();
    of=planar_embedding(faces);
    if (of<0)
        return false; // the graph is not planar
    /* subdivide the concave faces */
    augment(faces,of,false);
    bool is_3conn=edge_count()==ec; // no temporary edges were added
    const ivector &outer_face=faces[of];
    m=outer_face.size();
    /* create a fake outer face */
    ivector fake_outer_face(m);
    int label=largest_integer_label();
    for (ivector_iter it=outer_face.begin();it!=outer_face.end();++it) {
        add_temporary_edge(*it,fake_outer_face[it-outer_face.begin()]=add_node(++label));
    }
    make_tutte_layout(x,fake_outer_face); // create the layout
    /* remove temporary vertices and edges */
    if (score!=NULL && !is_3conn) { // compute the ratio of the largest to smallest face
        double fa_min=DBL_MAX,fa_max=0;
        ivectors_iter it=faces.begin();
        for (;it!=faces.end();++it) {
            if (int(it-faces.begin())==of)
                continue; // skip the outer face since its area is infinitely large
            int fs=it->size(),i=0;
            for (;i<fs && !is_temporary_edge(it->at(i),it->at((i+1)%fs));++i);
            if (i==fs) {
                layout fl(fs);
                for (i=0;i<fs;++i) fl[i]=x[it->at(i)];
                double a=poly_area(fl);
                if (a<fa_min) fa_min=a;
                if (a>fa_max) fa_max=a;
            }
        }
        *score=fa_min==0?DBL_MAX:fa_max/fa_min;
    } else if (score!=NULL) *score=0;
    remove_temporary_edges();
    while (node_count()>n) nodes.pop_back();
    x.resize(n);
    return true;
}

bool graphe::make_best_planar_layout(layout &x,int iter) {
    assert(iter>0);
    double s,min_s=-1;
    for (int i=0;i<iter;++i) {
        layout y;
        if (!make_planar_layout(y,&s))
            return false;
        if (s==0) { // the graph is 3-connected
            x=y;
            break;
        }
        if (min_s<0 || s<min_s) {
            min_s=s;
            x=y;
        }
    }
    return true;
}

/* get the convex hull of 2D graph layout */
void graphe::convex_hull(const layout &x,layout &hull) {
    vecteur v(x.size());
    for (iterateur it=v.begin();it!=v.end();++it) {
        const point &p=x[it-v.begin()];
        *it=makecomplex(p.front(),p.back());
    }
    vecteur h=*_convexhull(v,ctx)._VECTptr;
    hull.resize(h.size());
    for (const_iterateur it=h.begin();it!=h.end();++it) {
        point &p=hull[it-h.begin()];
        p.resize(2);
        p.front()=to_real_number(*it->_CPLXptr,ctx).to_double(ctx);
        p.back()=to_real_number(*(it->_CPLXptr+1),ctx).to_double(ctx);
    }
}

/* rotate layout x such that left to right symmetry is exposed */
void graphe::layout_best_rotation(layout &x) {
    /* center layout in the origin */
    point center=layout_center(x);
    scale_point(center,-1);
    translate_layout(x,center);
    layout hull;
    convex_hull(x,hull);
    int n=hull.size();
    /* compute the hull perimeter */
    double perim=0,score,maxscore=-1;
    point tmp(2),mp(2),rp(2),rq(2),origin(2,0.0),proj(2);
    for (int i=0;i<n;++i) {
        const point &p=hull[i],&q=hull[(i+1)%n];
        perim+=point_distance(p,q,tmp);
    }
    double a,b,c,d,tol,angle=0,half_perim=perim/2.0;
    bvector matched(hull.size());
    for (int i=0;i<n;++i) {
        const point &p=hull[i];
        for (int j=i+1;j<n;++j) {
            const point &q=hull[j];
            if (!point2segment_projection(origin,p,q,proj))
                continue;
            /* obtain the perpendicular bisector ax+by+c=0 of the segment pq */
            copy_point(p,mp);
            add_point(mp,q);
            scale_point(mp,0.5);
            if (p.front()==q.front()) {
                a=1;
                c=-mp.front();
                b=0;
            } else {
                a=(p.back()-q.back())/(p.front()-q.front());
                b=-1;
                c=mp.back()-a*mp.front();
            }
            /* try to match hull edges that are reflections of each other */
            score=0;
            std::fill(matched.begin(),matched.end(),false);
            for (layout_iter it=hull.begin();it!=hull.end();++it) {
                if (matched[it-hull.begin()])
                    continue;
                const point &p1=*it,&q1=it+1==hull.end()?hull.front():*(it+1);
                tol=(d=point_distance(p1,q1,tmp))*MARGIN_FACTOR/2.0;
                point_mirror(a,b,c,p1,rp);
                point_mirror(a,b,c,q1,rq);
                for (layout_iter jt=it;jt!=hull.end();++jt) {
                    if (matched[jt-hull.begin()])
                        continue;
                    const point &p2=*jt,&q2=jt+1==hull.end()?hull.front():*(jt+1);
                    if ((points_coincide(rp,p2,tol) && points_coincide(rq,q2,tol)) ||
                            (points_coincide(rp,q2,tol) && points_coincide(rq,p2,tol))) {
                        score+=d;
                        matched[it-hull.begin()]=true;
                        matched[jt-hull.begin()]=true;
                    }
                }
            }
            score/=half_perim;
            double bigpart=std::pow(PLASTIC_NUMBER,4); // hull symmetry is more important
            score*=bigpart;
            score+=1.0-point_distance(proj,mp,tmp)/point_distance(p,q,tmp);
            score/=bigpart+1;
            if (score>maxscore) {
                maxscore=score;
                angle=a==0?M_PI_2:std::atan(b/a);
            }
        }
    }
    rotate_layout(x,M_PI_2-angle);
    rectangle rect=layout_bounding_rect(x);
    if (rect.y()+rect.height()/2.0<0)
        rotate_layout(x,M_PI);
}

/* return true iff the degrees of vertices in v all coincide */
bool graphe::degrees_equal(const ivector &v,int deg) const {
    int D=deg,d;
    for (ivector_iter it=v.begin();it!=v.end();++it) {
        d=degree(*it);
        if (D==0)
            D=d;
        else if (d!=D)
            return false;
    }

    return true;
}

/* customize the Giac display */
gen graphe::customize_display(int options) {
    return symb_equal(at_display,change_subtype(options,_INT_COLOR));
}

/* append the line segment [p,q] to vecteur v */
void graphe::append_segment(vecteur &drawing,const point &p,const point &q,int color,int width,int style,bool arrow) const {
    vecteur attributs(1,color | width | style);
    vecteur seg=p.size()==2?
                makevecteur(makecomplex(p[0],p[1]),makecomplex(q[0],q[1]))
            : makevecteur(gen(makevecteur(p[0],p[1],p[2]),_POINT__VECT),gen(makevecteur(q[0],q[1],q[2]),_POINT__VECT));
    drawing.push_back(pnt_attrib(gen(seg,arrow?_VECTOR__VECT:_GROUP__VECT),attributs,ctx));
}

/* append the vertex (as a circle) to vecteur v */
void graphe::append_node(vecteur &drawing,const point &p,int color,int width,int shape) const {
    gen P=point2gen(p,true),args;
    args=makesequence(P,customize_display(color | width | shape));
    drawing.push_back(_point(args,ctx));
}

/* append label to vecteur v at the specified quadrant */
void graphe::append_label(vecteur &drawing,const point &p,const gen &label,int quadrant,int color) const {
    gen P=point2gen(p),args=makesequence(P,label,customize_display(quadrant | color));
    drawing.push_back(_legende(args,ctx));
}

/* convert gen to point coordinates */
bool graphe::gen2point(const gen &g,point &p,GIAC_CONTEXT) {
    if (g.type==_VECT || g.is_symb_of_sommet(at_point)) {
        const vecteur &v=g.type==_VECT?*g._VECTptr:*g._SYMBptr->feuille._VECTptr;
        p.resize(v.size());
        for (const_iterateur it=v.begin();it!=v.end();++it) {
            if (!is_real_number(*it,contextptr))
                return false;
            p[it-v.begin()]=to_real_number(*it,contextptr).to_double(contextptr);
        }
    } else { /* assuming that pos is a complex number */
        p.resize(2);
        if (g.type==_CPLX) {
            gen &real=*g._CPLXptr,&imag=*(g._CPLXptr+1);
            if (!is_real_number(real,contextptr) || !is_real_number(imag,contextptr))
                return false;
            p.front()=to_real_number(real,contextptr).to_double(contextptr);
            p.back()=to_real_number(imag,contextptr).to_double(contextptr);
        } else {
            if (!is_real_number(g,contextptr))
                return false;
            p.front()=to_real_number(g,contextptr).to_double(contextptr);
            p.back()=0;
        }
    }
    return true;
}

/* extract position attribute from attr */
bool graphe::get_node_position(const attrib &attr,point &p,GIAC_CONTEXT) {
    attrib_iter it=attr.find(_GT_ATTRIB_POSITION);
    if (it==attr.end())
        return false;
    return gen2point(it->second,p,contextptr);
}

/* add a single edge to the drawing */
vecteur graphe::draw_edge(int i,int j,const layout &x) const {
    int color,width,style;
    const attrib &attr=node(i).neighbor_attributes(j);
    color=default_edge_color;
    attrib_iter ait;
    if ((ait=attr.find(_GT_ATTRIB_COLOR))!=attr.end() && ait->second.is_integer()) {
        color=ait->second.val;
        if (color==7) color=0; /* draw white (invisible) edges as black */
    } else if ((ait=attr.find(_GT_ATTRIB_TEMPORARY))!=attr.end() && giac::is_one(ait->second))
        color=29; // gray
    style=0;
    width=default_edge_width;
    if ((ait=attr.find(_GT_ATTRIB_STYLE))!=attr.end()) {
        const gen &val=ait->second;
        string s=val.type==_STRNG?genstring2str(val):gen2string(val);
        if (s=="dashed")
            style=_DASH_LINE;
        else if (s=="dotted")
            style=_DOT_LINE;
        else if (s=="bold")
            width=bold_edge_width;
    }
    if ((ait=attr.find(_GT_ATTRIB_WIDTH))!=attr.end()) {
        const gen &val=ait->second;
        string s=val.type==_STRNG?genstring2str(val):gen2string(val);
        if (s=="thin" || s=="1")
            width=_LINE_WIDTH_1;
        else if (s=="normal" || s=="2")
            style=_LINE_WIDTH_2;
        else if (s=="thick" || s=="3")
            width=_LINE_WIDTH_3;
        else if (s=="thicker" || s=="4")
            width=_LINE_WIDTH_4;
        else if (s=="very thick" || s=="5")
            width=_LINE_WIDTH_5;
        else if (s=="6")
            width=_LINE_WIDTH_6;
        else if (s=="7")
            width=_LINE_WIDTH_7;
        else if (s=="8")
            width=_LINE_WIDTH_8;
    }
    vecteur res;
    point r(x.front().size());
    const point &p=x[i],&q=x[j];
    if (is_directed()) {
        assert((ait=attr.find(_GT_ATTRIB_POSITION))!=attr.end());
        double d=to_real_number(ait->second,ctx).to_double(ctx);
        point_lincomb(p,q,d,1-d,r);
        append_segment(res,p,r,color,width,style,true);
        append_segment(res,r,q,color,width,style,false);
    } else
        append_segment(res,p,q,color,style,width);
    return res;
}

/* append line segments representing edges (arcs) of the graph to vecteur v */
void graphe::draw_edges(vecteur &drawing,const layout &x) {
    if (x.empty())
        return;
    ipairs E,white,standard,highlighted,other;
    get_edges_as_pairs(E);
    attrib_iter ait;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        const attrib &attr=node(it->first).neighbor_attributes(it->second);
        if ((ait=attr.find(_GT_ATTRIB_COLOR))!=attr.end() && ait->second.is_integer()) {
            if (ait->second.val==default_highlighted_edge_color)
                highlighted.push_back(*it);
            else if (ait->second.val==default_edge_color)
                standard.push_back(*it);
            else if (ait->second.val==_WHITE)
                white.push_back(*it);
            else other.push_back(*it);
        } else if ((ait=attr.find(_GT_ATTRIB_TEMPORARY))!=attr.end() && giac::is_one(ait->second))
            other.push_back(*it);
        else standard.push_back(*it);
    }
    E=white;
    E.insert(E.end(),standard.begin(),standard.end());
    E.insert(E.end(),other.begin(),other.end());
    E.insert(E.end(),highlighted.begin(),highlighted.end());
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        vecteur res=draw_edge(it->first,it->second,x);
        drawing=mergevecteur(drawing,res);
    }
}

/* append points representing vertices of the graph to vecteur v */
void graphe::draw_nodes(vecteur &drawing,const layout &x) const {
    if (x.empty())
        return;
    int color,width,n=node_count();
    if (n<=30)
        width=_POINT_WIDTH_4;
    else if (n<=130)
        width=_POINT_WIDTH_3;
    else if (n<=330)
        width=_POINT_WIDTH_2;
    else
        width=_POINT_WIDTH_1;
    int shape,i=0;
    attrib_iter ait;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        const attrib &attr=it->attributes();
        const point &p=x[i];
        color=default_vertex_color;
        if ((ait=attr.find(_GT_ATTRIB_COLOR))!=attr.end() && ait->second.is_integer()) {
            color=ait->second.val;
            if (color==7)
                color=0; // draw white (invisible) nodes as black
        }
        shape=_POINT_POINT;
        if ((ait=attr.find(_GT_ATTRIB_SHAPE))!=attr.end()) {
            const gen &val=ait->second;
            string s=val.type==_STRNG?genstring2str(val):gen2string(val);
            if (s=="box" || s=="square")
                shape=_POINT_CARRE;
            else if (s=="triangle")
                shape=_POINT_TRIANGLE;
            else if (s=="diamond")
                shape=_POINT_LOSANGE;
            else if (s=="star")
                shape=_POINT_ETOILE;
            else if (s=="plus")
                shape=_POINT_PLUS;
        }
        append_node(drawing,p,color,width,shape);
    }
}

/* return the best quadrant for the placement of the i-th vertex label
* (minimize the collision with the adjacent edges) */
int graphe::best_quadrant(const point &p,const layout &adj) {
    int bestq,n=adj.size();
    if (n==0 || p.size()!=2)
        return _QUADRANT1;
    layout quad(4);
    quad[0]=make_point(0.7071,0.7071);
    quad[1]=make_point(-0.7071,0.7071);
    quad[2]=make_point(-0.7071,-0.7071);
    quad[3]=make_point(0.7071,-0.7071);
    vector<double> min_angular_dist(4,M_PI);
    point u(2);
    double a;
    for (int i=0;i<4;++i) {
        const point &q=quad[i];
        double &mindist=min_angular_dist[i];
        for (layout_iter it=adj.begin();it!=adj.end();++it) {
            scale_point(u,1/point_distance(p,*it,u));
            a=std::acos(point_dotprod(u,q));
            if (a<mindist)
                mindist=a;
        }
    }
    a=-1;
    for (int i=0;i<4;++i) {
        const double &dist=min_angular_dist[i];
        if (a<0 || dist>a) {
            a=dist;
            bestq=i;
        }
    }
    switch (bestq) {
    case 0: return _QUADRANT1;
    case 1: return _QUADRANT2;
    case 2: return _QUADRANT3;
    case 3: return _QUADRANT4;
    default: break;
    }
    /* unreachable */
    return -1;
}

/* append labels of the nodes of this graph to vecteur v */
void graphe::draw_labels(vecteur &drawing,const layout &x) const {
    if (is_null())
        return;
    assert(!x.empty());
    ivector adjv;
    layout adj(2);
    attrib_iter ait;
    double d;
    int color,n=node_count();
    point r(x.front().size());
    bool isdirected=is_directed(),isweighted=is_weighted();
    if (isweighted && x.front().size()==2) {
        /* draw weight labels */
        for (int i=0;i<n;++i) {
            const point &p=x[i];
            const vertex &v=node(i);
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                if (!isdirected && *it<i)
                    continue;
                const point &q=x[*it];
                const attrib &attr=v.neighbor_attributes(*it);
                color=default_edge_color;
                if ((ait=attr.find(_GT_ATTRIB_COLOR))!=attr.end())
                    color=ait->second.val;
                assert((ait=attr.find(_GT_ATTRIB_POSITION))!=attr.end());
                d=to_real_number(ait->second,ctx).to_double(ctx);
                point_lincomb(p,q,d,1-d,r);
                adj.front()=p;
                adj.back()=q;
                ait=attr.find(_GT_ATTRIB_WEIGHT);
                append_label(drawing,r,ait==attr.end()?undef:ait->second,best_quadrant(r,adj),color);
            }
        }
    }
    /* draw vertex labels */
    for (int i=0;i<n;++i) {
        const point &p=x[i];
        adjacent_nodes(i,adjv);
        adj.resize(adjv.size());
        for (ivector_iter it=adjv.begin();it!=adjv.end();++it) {
            adj[it-adjv.begin()]=x[*it];
        }
        append_label(drawing,p,node(i).label(),best_quadrant(p,adj));
    }
}

void graphe::rdfs(int i,ivector &d,bool rec,int sg,bool skip_embedded) {
    vertex &v=node(i);
    v.set_visited(true);
    v.set_disc(++disc_time);
    v.set_low(v.disc());
    if (rec)
        d.push_back(i);
    int j;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        j=*it;
        vertex &w=node(j);
        if ((sg>=0 && w.subgraph()!=sg) || (skip_embedded && w.is_embedded()))
            continue;
        if (!w.is_visited()) {
            w.set_ancestor(i);
            v.set_leaf(false);
            rdfs(j,d,rec,sg,skip_embedded);
            v.set_low(std::min(v.low(),w.low()));
        } else if (j!=v.ancestor())
            v.set_low(std::min(v.low(),w.disc()));
    }
}

/* depth-first graph traversal with O(n+m) time and O(m) space complexity */
void graphe::dfs(int root,bool rec,bool clr,ivector *D,int sg,bool skip_embedded) {
    if (clr) {
        unvisit_all_nodes(sg);
        unset_all_ancestors(sg);
        for (std::vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
            if (sg<0 || it->subgraph()==sg)
                it->set_leaf(true);
        } 
        disc_time=0;
    }
    ivector &d=D!=NULL?*D:disc_nodes;
    if (rec) {
        d.clear();
        d.reserve(node_count());
    }
    assert(sg<0 || node(root).subgraph()==sg);
    rdfs(root,d,rec,sg,skip_embedded);
}

/* breadth-first graph traversal with O(n+m) time and O(m) space complexity */
void graphe::bfs(int root,bool rec,bool clr,ivector *D,int sg,bool skip_embedded) {
    if (clr) {
        unvisit_all_nodes(sg);
        unset_all_ancestors(sg);
        disc_time=0;
    }
    ivector &d=D!=NULL?*D:disc_nodes;
    if (rec) {
        d.clear();
        d.reserve(node_count());
    }
    assert(node_queue.empty());
    node_queue.push(root);
    int i,j;
    while (!node_queue.empty()) {
        i=node_queue.front();
        node_queue.pop();
        vertex &v=node(i);
        if (!v.is_visited()) {
            v.set_disc(disc_time++);
            if (rec)
                d.push_back(i);
            v.set_visited(true);
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                vertex &w=node(j=*it);
                if ((sg>=0 && w.subgraph()!=sg) || (skip_embedded && w.is_embedded()))
                    continue;
                if (!w.is_visited()) {
                    if (w.ancestor()<0)
                        w.set_ancestor(i);
                    node_queue.push(j);
                }
            }
        }
    }
}

/* return true iff the graph is connected */
bool graphe::is_connected(int sg) {
    assert(!is_null());
    if (is_directed()) {
        graphe G(ctx,false);
        underlying(G);
        return G.is_connected(sg);
    }
    int i=sg<0?0:first_vertex_from_subgraph(sg);
    assert(i>=0);
    dfs(i,false,true,NULL,sg);
    node_iter it=nodes.begin()+i,itend=nodes.end();
    for (;it!=itend;++it) {
        if ((sg<0 || it->subgraph()==sg) && !it->is_visited())
            return false;
    }
    return true;
}

/* return true iff the graph is biconnected */
bool graphe::is_biconnected(int sg) {
    assert(!is_null());
    if (is_directed()) {
        graphe G(ctx,false);
        underlying(G);
        return G.is_biconnected(sg);
    }
    return is_connected(sg) && !has_cut_vertex(sg);
}

/* return true iff the graph is triconnected, using a simple O(n*(n+m)) algorithm */
bool graphe::is_triconnected(int sg) {
    assert(!is_null());
    if (is_directed()) {
        graphe G(ctx,false);
        underlying(G);
        return G.is_triconnected(sg);
    }
    int color=max_subgraph_index()+1;
    for (int i=node_count();i-->0;) {
        vertex &v=node(i);
        if (sg>=0 && v.subgraph()!=sg)
            continue;
        if (degree(i)<3)
            return false;
        v.set_subgraph(color);
    }
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        if (it->subgraph()!=color)
            continue;
        it->set_subgraph(color+1);
        if (!is_biconnected(color))
            return false;
        it->set_subgraph(color);
    }
    return true;
}

/* return true iff the (sub)graph is a cycle, linear time */
bool graphe::is_cycle(ipairs &E,int sg) {
    int root=-1,i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        if (sg>=0 && it->subgraph()!=sg)
            continue;
        if (degree(i,sg)>0) {
            root=i;
            break;
        }
    }
    if (root<0)
        return false;
    dfs(root,true,true,NULL,sg);
    int n=disc_nodes.size();
    ipairs ev;
    get_edges_as_pairs(ev,sg);
    if (n!=int(ev.size()))
        return false;
    ipair e;
    E.resize(n);
    for (i=0;i<n;++i) {
        e=make_pair(disc_nodes[i],disc_nodes[(i+1)%n]);
        if (!has_edge(e))
            return false;
        E[i]=e;
    }
    return true;
}

/* return true iff the graph is a forest (check that the connected components are all trees) */
bool graphe::is_forest(int sg,gt_conn_check cc) {
    assert(!is_null() && !is_directed());
    int s,cmp;
    switch (cc) {
    case _GT_CC_CONNECTED:
        return is_tree(sg);
    case _GT_CC_FIND_COMPONENTS:
        s=sg<0?1:1+max_subgraph_index();
        cmp=connected_components_to_subgraphs(sg);
    case _GT_CC_COMPONENTS_ARE_SUBGRAPHS:
        s=1;
        cmp=max_subgraph_index();
        for (;s<=cmp;++s) {
            if (edge_count(s)+1!=subgraph_size(s))
                return false;
        }
        return true;
    default: break;
    }
    assert(false);
}

/* return true iff the graph is a tournament */
bool graphe::is_tournament(int sg) const {
    int n=subgraph_size(sg);
    if (!is_directed() || edge_count(sg)!=n*(n-1)/2)
        return false;
    for (int i=0;i<n;++i) {
        if (sg>=0 && node(i).subgraph()!=sg)
            continue;
        for (int j=0;j<n;++j) {
            if (i==j || (sg>=0 && node(j).subgraph()!=sg))
                continue;
            if (has_edge(i,j)==has_edge(j,i))
                return false;
        }
    }
    return true;
}

/* return true iff the graph is planar */
bool graphe::is_planar(int sg) {
    ivectors comp,faces;
    int m,s=sg<0?1:1+max_subgraph_index();
    int cmp=connected_components_to_subgraphs(sg);
    int tmp_sg=max_subgraph_index()+1;
    for (;s<=cmp;++s) {
        int sz=subgraph_size(s);
        if (sz<5) continue;
        m=edge_count(s);
        if (m>3*sz-6) {
            return false;
        }
        if (m<9) continue;
        vector<ipairs> blocks;
        find_blocks(blocks,s);
        for (vector<ipairs>::const_iterator jt=blocks.begin();jt!=blocks.end();++jt) {
            if (jt->size()<9) continue;
            set_subgraph(*jt,tmp_sg);
            if (subgraph_size(tmp_sg)<5) continue;
            if (!demoucron(faces,tmp_sg))
                return false;
            renumber_subgraph(tmp_sg,s);
        }
        if (sg>=0)
            renumber_subgraph(s,sg);
    }
    return true;
}

/* create set of vertices for product P of this graph and graph G */
void graphe::make_product_nodes(const graphe &G,graphe &P) const {
    int n=node_count(),m=G.node_count();
    for (int i=0;i<n;++i) {
        for (int j=0;j<m;++j) {
            const gen &v=node_label(i),&w=G.node_label(j);
            string str;
            str+=v.type==_STRNG?genstring2str(v):gen2string(v);
            str+=":";
            str+=w.type==_STRNG?genstring2str(w):gen2string(w);
            P.add_node(str2gen(str,true));
        }
    }
}

/* compute the cartesian product of this graph and graph G and store it in P */
void graphe::cartesian_product(const graphe &G,graphe &P) const {
    assert(supports_attributes()==G.supports_attributes());
    assert(supports_attributes()==P.supports_attributes());
    int n=node_count(),m=G.node_count();
    P.clear();
    P.reserve_nodes(n*m);
    if (supports_attributes())
        make_product_nodes(G,P);
    else P.add_nodes(n*m);
    bool G_isdir=G.is_directed(),isdir=is_directed();
    P.set_directed(G_isdir || isdir);
    for (int i=0;i<n;++i) {
        for (int j=0;j<m;++j) {
            const vertex &w=G.node(j);
            for (ivector_iter it=w.neighbors().begin();it!=w.neighbors().end();++it) {
                if (G_isdir || *it>j)
                    P.add_edge(i*m+j,i*m+(*it));
            }
        }
    }
    for (int j=0;j<m;++j) {
        for (int i=0;i<n;++i) {
            const vertex &v=node(i);
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                if (isdir || *it>i)
                    P.add_edge(i*m+j,(*it)*m+j);
            }
        }
    }
}

/* compute tensor product of this graph and graph G and store it in P */
void graphe::tensor_product(const graphe &G,graphe &P) const {
    assert(supports_attributes()==G.supports_attributes());
    assert(supports_attributes()==P.supports_attributes());
    int n=node_count(),m=G.node_count();
    P.clear();
    P.reserve_nodes(m*n);
    if (supports_attributes())
        make_product_nodes(G,P);
    else P.add_nodes(m*n);
    bool G_isdir=G.is_directed(),isdir=is_directed();
    P.set_directed(G_isdir || isdir);
    for (int i=0;i<n;++i) {
        const vertex &v=node(i);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            for (int j=0;j<m;++j) {
                const vertex &w=G.node(j);
                for (ivector_iter jt=w.neighbors().begin();jt!=w.neighbors().end();++jt) {
                    if (G_isdir || *jt>j)
                        P.add_edge(i*m+j,(*it)*m+(*jt));
                }
            }
        }
    }
}

/* compute the distance between i-th and J nodes using BFS and store them to dist
* (also output paths if shortest_paths!=NULL) */
void graphe::distance(int i,const ivector &J,ivector &dist,ivectors *shortest_paths) {
    bfs(i,false);
    int k,len;
    if (shortest_paths!=NULL) {
        shortest_paths->resize(J.size());
    }
    dist.resize(J.size());
    ivector *shortest_path=NULL;
    for (ivector_iter it=J.begin();it!=J.end();++it) {
        k=*it;
        len=0;
        if (shortest_paths!=NULL) {
            shortest_path=&shortest_paths->at(it-J.begin());
            shortest_path->clear();
            shortest_path->push_back(*it);
        }
        while (k!=i) {
            if ((k=node(k).ancestor())<0) {
                len=-1;
                break;
            }
            if (shortest_path!=NULL)
                shortest_path->push_back(k);
            ++len;
        }
        if (shortest_path!=NULL)
            std::reverse(shortest_path->begin(),shortest_path->end());
        dist[it-J.begin()]=len;
    }
}

/* compute the distances between all pairs of vertices */
void graphe::allpairs_distance(matrice &m) {
    int n=node_count(),i,j,k;
    m.reserve(n);
    for (i=0;i<n;++i) m.push_back(vecteur(n,plusinf()));
    gen s;
    if (is_weighted()) { // Floyd & Warshall algorithm, O(n^3)
        node_iter it=nodes.begin(),itend=nodes.end();
        for (i=0;it!=itend;++it,++i) {
            m[i]._VECTptr->at(i)=0;
            for (ivector_iter jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
                j=*jt;
                m[i]._VECTptr->at(j)=weight(i,j);
            }
        }
        for (k=0;k<n;++k) {
            for (i=0;i<n;++i) {
                const gen &mik=m[i][k];
                for (j=0;j<n;++j) {
                    s=mik+m[k][j];
                    gen &mij=m[i]._VECTptr->at(j);
                    if (is_strictly_greater(mij,s,ctx))
                        mij=s;
                }
            }
        }
    } else { // compute shortest paths using BFS
        ivector V(n),dist;
        ivector_iter it,itend;
        for (i=0;i<n;++i) V[i]=i;
        for (i=0;i<n;++i) {
            distance(i,V,dist);
            for (it=dist.begin(),itend=dist.end();it!=itend;++it) {
                j=it-dist.begin();
                if (*it>=0)
                    m[i]._VECTptr->at(j)=*it;
            }
        }
    }
}

/* return the length of the shortest path from src to dest in weighted
* graph (Dijkstra's algorithm), also fill shortest_path with the respective vertices */
void graphe::dijkstra(int src,const ivector &dest,vecteur &path_weights,ivectors *cheapest_paths,int sg) {
    int n=node_count();
    ivector Q;
    vecteur dist(n);
    bool isweighted=is_weighted();
    if (sg>=0) {
        assert(node(src).subgraph()==sg);
        for (ivector_iter it=dest.begin();it!=dest.end();++it) {
            assert(node(*it).subgraph()==sg);
        }
    }
    unset_all_ancestors(sg);
    for (int i=0;i<n;++i) {
        vertex &v=node(i);
        if (sg>=0 && v.subgraph()!=sg) continue;
        Q.push_back(i);
        dist[i]=i==src?gen(0):plusinf();
    }
    unvisit_all_nodes(sg);
    gen min_dist,alt;
    int pos,u;
    while (!Q.empty()) {
        min_dist=plusinf();
        pos=-1;
        for (ivector_iter it=Q.begin();it!=Q.end();++it) {
            if (is_strictly_greater(min_dist,dist[*it],ctx)) {
                min_dist=dist[*it];
                pos=it-Q.begin();
                u=*it;
            }
        }
        if (pos<0) break;
        Q.erase(Q.begin()+pos);
        vertex &v=node(u);
        v.set_visited(true);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            vertex &w=node(*it);
            if ((sg>=0 && w.subgraph()!=sg) || w.is_visited())
                continue;
            alt=dist[u]+(isweighted?weight(u,*it):gen(1));
            if (is_strictly_greater(dist[*it],alt,ctx)) {
                dist[*it]=alt;
                w.set_ancestor(u);
            }
        }
    }
    path_weights.resize(dest.size());
    for (ivector_iter it=dest.begin();it!=dest.end();++it) {
        path_weights[it-dest.begin()]=dist[*it];
    }
    if (cheapest_paths!=NULL) {
        cheapest_paths->resize(dest.size());
        for (ivector_iter it=dest.begin();it!=dest.end();++it) {
            ivector &path=cheapest_paths->at(it-dest.begin());
            path.clear();
            if (is_inf(dist[*it])) continue;
            path.push_back(*it);
            int p=*it;
            while ((p=node(p).ancestor())>=0) path.push_back(p);
            std::reverse(path.begin(),path.end());
        }
    }
}

/* return the length of the shortest path from src to dest in weighted
* graph (Bellman-Ford algorithm), also fill shortest_path with the respective vertices */
bool graphe::bellman_ford(int src,const ivector &dest,vecteur &path_weights,ivectors *cheapest_paths) {
    int n=node_count(),u,v;
    ivector prev(n);
    vecteur dist(n);
    bool isweighted=is_weighted();
    for (int i=0;i<n;++i) {
        dist[i]=i==src?gen(0):plusinf();
        prev[i]=-1;
    }
    ipairs E;
    get_edges_as_pairs(E);
    gen w(1);
    for (int i=1;i<n;++i) {
        for (ipairs_iter it=E.begin();it!=E.end();++it) {
            u=it->first; v=it->second;
            if (isweighted)
                get_edge_attribute(u,v,_GT_ATTRIB_WEIGHT,w);
            if (is_strictly_greater(dist[v],dist[u]+w,ctx)) {
                dist[v]=dist[u]+w;
                prev[v]=u;
            }
        }
    }
    /* check for negative-weight cycles */
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        if (isweighted)
            get_edge_attribute(u,v,_GT_ATTRIB_WEIGHT,w);
        if (is_strictly_greater(dist[v],dist[u]+w,ctx))
            return false; // a negative-weight cycle is found
    }
    if (cheapest_paths!=NULL) {
        cheapest_paths->resize(dest.size());
        for (ivector_iter it=dest.begin();it!=dest.end();++it) {
            ivector &path=cheapest_paths->at(it-dest.begin());
            path.clear();
            path.push_back(*it);
            int p=*it;
            while ((p=prev[p])>=0) path.push_back(p);
            std::reverse(path.begin(),path.end());
        }
    }
    path_weights.resize(dest.size());
    for (ivector_iter it=dest.begin();it!=dest.end();++it) {
        path_weights[it-dest.begin()]=dist[*it];
    }
    return true; // success
}

/* return true iff the graph has a topological ordering and output the indices
* of thus sorted vertices in list v, use Kahn's algorithm with time complexity O(n+m) */
bool graphe::topologic_sort(ivector &ordering) {
    assert(is_directed() && node_stack.empty());
    int n=node_count(),v;
    ordering.clear();
    ordering.reserve(n);
    for (int i=0;i<n;++i) {
        if (in_degree(i)==0)
            node_stack.push(i);
    }
    graphe G(ctx,false);
    copy(G);
    while (!node_stack.empty()) {
        v=node_stack.top();
        node_stack.pop();
        ordering.push_back(v);
        vertex &vert=G.node(v);
        for (ivector_iter it=vert.neighbors().begin();it!=vert.neighbors().end();++it) {
            if (G.in_degree(*it)==1)
                node_stack.push(*it);
        }
        vert.clear_neighbors();
    }
    return G.edge_count()==0;
}

/* return true iff the graph is an arborescence, i.e. a directed tree */
bool graphe::is_arborescence() const {
    assert(is_directed());
    bool has_root=false;
    ivector deg(node_count(),0);
    node_iter it=nodes.begin(),itend=nodes.end();
    ivector_iter jt,jtend;
    for (;it!=itend;++it) {
        for (jt=it->neighbors().begin(),jtend=it->neighbors().end();jt!=jtend;++jt) {
            ++deg[*jt];
        }
    }
    for (jt=deg.begin(),jtend=deg.end();jt!=jtend;++jt) {
        if (*jt==0) {
            if (has_root)
                return false;
            has_root=true;
        } else if (*jt>1)
            return false;
    }
    return has_root;
}

/* write this (directed) graph with arc directions reversed to G */
void graphe::reverse(graphe &G) const {
    assert(is_directed());
    G.set_directed(true);
    G.reserve_nodes(node_count());
    int i;
    node_iter it,itend;
    ivector_iter jt,jtend;
    if (supports_attributes()) {
        G.set_graph_attributes(attributes);
        for (it=nodes.begin(),itend=nodes.end();it!=itend;++it) {
            G.add_node(it->label(),it->attributes());
        }
    } else G.add_nodes(node_count());
    for (it=nodes.begin(),itend=nodes.end(),i=0;it!=itend;++it,++i) {
        for (jt=it->neighbors().begin(),jtend=it->neighbors().end();jt!=jtend;++jt) {
            if (supports_attributes())
                G.add_edge(*jt,i,it->neighbor_attributes(*jt));
            else G.add_edge(*jt,i);
        }
    }
}

/* write a spanning tree of this graph with i-th node as root to T, time complexity O(n+m) */
void graphe::spanning_tree(int i,graphe &T,int sg) {
    T.clear();
    vecteur V=vertices(sg);
    T.reserve_nodes(V.size());
    T.add_nodes(V);
    ivector indices(V.size());
    int v,p;
    if (sg>=0) {
        const_iterateur vt=V.begin(),vtend=V.end();
        for (v=0;vt!=vtend;++vt,++v) {
            indices[v]=node_index(*vt);
        }
    }
    dfs(i,false,true,NULL,sg);
    node_iter it=nodes.begin(),itend=nodes.end();
    for (v=0;it!=itend;++it,++v) {
        if (sg>=0 && it->subgraph()!=sg)
            continue;
        if ((p=it->ancestor())>=0) {
            if (sg>=0)
                T.add_edge(indices[v],indices[p]);
            else
                T.add_edge(v,p);
        }
    }
}

/* write the minimal spanning tree of this graph to T,
* use Kruskal's algorithm with time complexity O(m*log(n)) */
void graphe::minimal_spanning_tree(graphe &T,int sg) {
    assert(!is_directed() && is_weighted());
    ipairs E,res;
    get_edges_as_pairs(E,sg);
    edges_comparator comp(this);
    std::sort(E.begin(),E.end(),comp);
    int v,u;
    unionfind ds(node_count());
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        ds.make_set(it->first);
        ds.make_set(it->second);
    }
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        u=it->first;
        v=it->second;
        if (ds.find(u)!=ds.find(v)) {
            res.push_back(*it);
            ds.unite(u,v);
        }
    }
    extract_subgraph(res,T);
}

/* Tarjan's offline algorithm for the lowest common ancestor, time complexity O(n) */
void graphe::lca_recursion(int u,const ipairs &p,ivector &lca,unionfind &ds) {
    ds.make_set(u);
    vertex &U=node(u);
    U.set_ancestor(u);
    U.set_visited(true);
    int v;
    for (ivector_iter it=U.neighbors().begin();it!=U.neighbors().end();++it) {
        v=*it;
        vertex &V=node(v);
        if (V.is_visited())
            continue;
        lca_recursion(v,p,lca,ds);
        ds.unite(u,v);
        node(ds.find(u)).set_ancestor(u);
    }
    U.set_color(1); // black
    node_iter it=nodes.begin(),itend=nodes.end();
    ipairs_iter jt,jtend;
    for (v=0;it!=itend;++it,++v) {
        if (u==v)
            continue;
        for (jt=p.begin(),jtend=p.end();jt!=jtend;++jt) {
            if ((jt->first==u && jt->second==v) || (jt->first==v && jt->second==u)) {
                if (it->color()==1)
                    lca[jt-p.begin()]=node(ds.find(v)).ancestor();
                break;
            }
        }
    }
}

/* find the lowest common ancestors for all pairs of vertices in p (this must be a tree) */
void graphe::lowest_common_ancestors(int root,const ipairs &p,ivector &lca) {
    unvisit_all_nodes();
    unset_all_ancestors();
    uncolor_all_nodes();
    lca.resize(p.size(),-1);
    unionfind ds(node_count());
    lca_recursion(root,p,lca,ds);
    assert(find(lca.begin(),lca.end(),-1)==lca.end());
}

/* return the lowest common ancestor of i-th node and j-th node of this tree */
int graphe::lowest_common_ancestor(int i,int j,int root) {
    ipairs p;
    p.push_back(make_pair(i,j));
    ivector lca;
    lowest_common_ancestors(root,p,lca);
    return lca.front();
}

/* return true iff v is a descendant of anc wrt the last dfs traversal */
bool graphe::is_descendant(int v,int anc) const {
    int i=v;
    while ((i=node(i).ancestor())>=0) {
        if (i==anc)
            return true;
    }
    return false;
}

void graphe::st_numbering_dfs(int i,ivector &preorder) {
    vertex &v=node(i);
    v.set_disc(++disc_time);
    v.set_low(i);
    int j;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        vertex &w=node(j=*it);
        if (w.disc()==0) {
            st_numbering_dfs(j,preorder);
            w.set_ancestor(i);
            preorder.push_back(j);
            if (node(w.low()).disc()<node(v.low()).disc())
                v.set_low(w.low());
        } else if (w.disc()!=0 && w.disc()<node(v.low()).disc())
            v.set_low(*it);
    }
}

/* find st-numbering using Tarjan's algorithm (1986), time complexity O(n+m) */
void graphe::compute_st_numbering(int s,int t) {
    /* assuming that the graph is biconnected */
    assert(has_edge(s,t) && node_stack.empty());
    int n=node_count(),i;
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->set_disc(0);
    }
    disc_time=1;
    node(s).set_disc(1);
    ivector preorder,L;
    st_numbering_dfs(t,preorder);
    L.push_back(s);
    L.push_back(t);
    bvector sign(n);
    ivector::iterator lit;
    bool sign_plus=true,sign_minus=false;
    sign[s]=sign_minus;
    for (ivector::const_reverse_iterator it=preorder.rbegin();it!=preorder.rend();++it) {
        vertex &v=node(*it);
        lit=find(L.begin(),L.end(),v.ancestor());
        assert(lit!=L.end());
        if (sign[v.low()]==sign_minus) {
            L.insert(lit,*it);
            sign[v.ancestor()]=sign_plus;
        } else if (sign[v.low()]==sign_plus) {
            L.insert(lit+1,*it);
            sign[v.ancestor()]=sign_minus;
        }
    }
    i=0;
    for (ivector_iter it=L.begin();it!=L.end();++it) {
        node(*it).set_number(++i);
    }
}

/* return st-numbering in form of vecteur */
vecteur graphe::get_st_numbering() const {
    /* assuming that st numbering has been computed */
    vecteur st(node_count());
    int ofs=1-array_start(ctx);
    for (iterateur it=st.begin();it!=st.end();++it) {
        *it=node(it-st.begin()).number()-ofs;
    }
    return st;
}

/* assign directions to edges according to the st-numbering */
void graphe::assign_edge_directions_from_st() {
    bool isweighted=is_weighted();
    matrice W;
    if (isweighted)
        weight_matrix(W);
    ipairs E;
    get_edges_as_pairs(E);
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        discard_edge_attribute(it->first,it->second,_GT_ATTRIB_WEIGHT);
    }
    set_weighted(false);
    set_directed(true);
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        vertex &v=*it;
        ivector ngh=v.neighbors();
        for (ivector_iter jt=ngh.begin();jt!=ngh.end();++jt) {
            const vertex &w=node(*jt);
            if (v.number()>w.number())
                v.remove_neighbor(*jt);
        }
    }
    if (isweighted)
        make_weighted(W);
}

void graphe::strec(int i,int t,int counter,int np,iset &Q,vecteur &timestamp,vecteur &l) {
    bool isweighted=is_weighted();
    ++counter;
    vertex &v=node(i);
    v.set_number(counter);
    gen w(1);
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        if (node(*it).subgraph()!=0)
            continue;
        if (*it!=t)
            Q.insert(*it);
        if (isweighted)
            get_edge_attribute(i,*it,_GT_ATTRIB_WEIGHT,w);
        if (is_strictly_greater(l[i]+w,l[*it],ctx))
            l[*it]=l[i]+w;
        timestamp[*it]=isweighted?l[*it]:gen(counter);
    }
    Q.erase(i);
    v.set_subgraph(1);
    if (Q.empty()) {
        node(t).set_number(node_count());
        return;
    }
    /* update blocks */
    vector<iset> B;
    vector<ipairs> blocks;
    find_blocks(blocks,0);
    int nb=blocks.size();
    B.resize(nb);
    for (vector<ipairs>::const_iterator it=blocks.begin();it!=blocks.end();++it) {
        iset &b=B[it-blocks.begin()];
        for (ipairs_iter jt=it->begin();jt!=it->end();++jt) {
            b.insert(jt->first);
            b.insert(jt->second);
        }
    }
    iset df;
    bool has_root=false;
    if (nb>1) { // G-v is not biconnected
        ivector cv;
        find_cut_vertices(cv,0);
        for (int k=nb;k-->0;) {
            iset &b=B[k];
            if (!has_root && b.find(t)!=b.end()) { // root block, delete it
                B.erase(B.begin()+k);
                has_root=true;
            } else {
                sets_intersection(cv,b,df);
                if (df.size()!=1) // non-leaf tree node, delete it
                    B.erase(B.begin()+k);
                else b.erase(*df.begin()); // b is leaf node, detach it
            }
        }
    }
    /* recurse */
    for (vector<iset>::const_iterator it=B.begin();it!=B.end();++it) {
        sets_intersection(*it,Q,df);
        if (df.empty()) continue;
        gen_map m;
        gen ts,max_ts(undef),min_ts(undef);
        for (iset_iter jt=df.begin();jt!=df.end();++jt) {
            ts=timestamp[*jt];
            max_ts=is_undef(max_ts)?ts:_max(makevecteur(max_ts,ts),ctx);
            min_ts=is_undef(min_ts)?ts:_min(makevecteur(min_ts,ts),ctx);
            if (m.find(ts)==m.end())
                m[ts]=vecteur(1,*jt);
            else m[ts]._VECTptr->push_back(*jt);
        }
        const vecteur &cand=*m[counter<=np?max_ts:min_ts]._VECTptr;
        int choice=(cand.size()==1?cand.front():cand[rand_integer(cand.size())]).val;
        strec(choice,t,counter,np,Q,timestamp,l);
    }
}

/* impose parametrized st-orientation on this graph with parameter p in [0,1] */
void graphe::parametrized_st_orientation(int s,int t,double p) {
    assert(p>=0 && p<=1);
    int n=node_count(),counter=0,np=(int)(n*p);
    vecteur timestamp(n,0),l(n,0);
    iset Q;
    Q.insert(s);
    unset_subgraphs(0);
    strec(s,t,counter,np,Q,timestamp,l);
    unset_subgraphs();
}

/* greedy vertex coloring algorithm due to Biggs, also records vertex inclusion time */
void graphe::greedy_vertex_coloring_biggs(ivector &ordering) {
    uncolor_all_nodes();
    int n=node_count(),k=0,i,maxdeg,d,col=0;
    ordering.resize(n);
    ivector_iter jt;
    node_iter it,itend;
    while (k<n) {
        ++col;
        do {
            maxdeg=-1;
            i=-1;
            for (it=nodes.begin(),itend=nodes.end();it!=itend;++it) {
                if (it->color()>0)
                    continue;
                for (jt=it->neighbors().begin();jt!=it->neighbors().end();++jt) {
                    if (node(*jt).color()==col)
                        break;
                }
                if (jt==it->neighbors().end() && (maxdeg<0 || (d=it->neighbors().size())>maxdeg)) {
                    maxdeg=d;
                    i=it-nodes.begin();
                }
            }
            if (i>=0) {
                node(i).set_color(col);
                ordering[k++]=i;
            }
        } while (i>=0);
    }
}

/* classical greedy vertex coloring algorithm, time complexity O(n+m) */
int graphe::greedy_vertex_coloring(const ivector &p) {
    if (is_directed()) {
        graphe G(ctx,false);
        underlying(G);
        return G.greedy_vertex_coloring(p);
    }
    uncolor_all_nodes();
    int c=0,k;
    std::set<int> used;
    for (ivector_iter it=p.begin();it!=p.end();++it) {
        vertex &v=node(*it);
        if (c==0) {
            v.set_color(++c);
            continue;
        }
        used.clear();
        for (ivector_iter jt=v.neighbors().begin();jt!=v.neighbors().end();++jt) {
            vertex &v=node(*jt);
            if (v.color()>0)
                used.insert(v.color());
        }
        k=1;
        for (std::set<int>::const_iterator jt=used.begin();jt!=used.end();++jt) {
            if (*jt>k)
                break;
            ++k;
        }
        v.set_color(k);
        if (k>c)
            c=k;
    }
    return c;
}

/* get color of the i-th node, either as attribute or from vertex */
int graphe::get_node_color(int i) const {
    int c;
    const vertex &v=node(i);
    attrib_iter it;
    if (supports_attributes()) {
        c=default_vertex_color;
        if ((it=v.attributes().find(_GT_ATTRIB_COLOR))!=v.attributes().end()) {
            c=it->second.val;
        }
    } else c=v.color();
    return c;
}

/* extract colors of the vertices and return them in order */
void graphe::get_node_colors(ivector &colors) const {
    colors.resize(node_count());
    node_iter it=nodes.begin(),itend=nodes.end();
    ivector::iterator jt=colors.begin();
    for (;it!=itend;++it,++jt) {
        *jt=it->color();
    }
}

/* return true iff the graph is bipartite, time complexity O(n+m) */
bool graphe::is_bipartite(ivector &V1,ivector &V2,int sg,gt_conn_check cc) {
    assert(node_queue.empty());
    int s,cmp;
    switch (cc) {
    case _GT_CC_CONNECTED:
        break;
    case _GT_CC_COMPONENTS_ARE_SUBGRAPHS:
    case _GT_CC_FIND_COMPONENTS:
        if (is_directed()) {
            graphe G(ctx,false);
            underlying(G);
            return G.is_bipartite(V1,V2,sg,cc);
        }
        V1.clear(); V1.reserve(node_count());
        V2.clear(); V2.reserve(node_count());
        s=(cc==_GT_CC_COMPONENTS_ARE_SUBGRAPHS || sg<0)?1:1+max_subgraph_index();
        cmp=(cc==_GT_CC_FIND_COMPONENTS)?connected_components_to_subgraphs(sg):max_subgraph_index();
        for (;s<=cmp;++s) {
            if (!is_bipartite(V1,V2,s,_GT_CC_CONNECTED))
                return false;
            if (cc==_GT_CC_FIND_COMPONENTS)
                renumber_subgraph(s,sg);
        }
        return true;
    default: assert(false);
    }
    // sg >= 0
    uncolor_all_nodes(-1,sg);
    int i=first_vertex_from_subgraph(sg);
    assert(i>=0);
    node(i).set_color(1);
    node_queue.push(i);
    ivector_iter it,itend;
    while (!node_queue.empty()) {
        i=node_queue.front();
        node_queue.pop();
        vertex &v=node(i);
        for (it=v.neighbors().begin(),itend=v.neighbors().end();it!=itend;++it) {
            vertex &w=node(*it);
            if (sg>=0 && w.subgraph()!=sg)
                continue;
            if (w.color()==-1) {
                w.set_color(1-v.color());
                node_queue.push(*it);
            } else if (w.color()==v.color()) {
                clear_node_queue();
                return false;
            }
        }
    }
    node_iter nt=nodes.begin(),ntend=nodes.end();
    for (i=0;nt!=ntend;++nt,++i) {
        if (sg>=0 && nt->subgraph()!=sg)
            continue;
        if (nt->color()==1)
            V1.push_back(i);
        else V2.push_back(i);
    }
    return true;
}

/* place vertices from partitions v1 and v2 on two parallel lines */
void graphe::make_bipartite_layout(layout &x,const ivector &p1,const ivector &p2) {
    double aspect_ratio;
    int n1=p1.size(),n2=p2.size(),n=std::max(n1,n2);
    assert(n1>1 && n2>1);
    if (n<3)
        aspect_ratio=1.0;
    else if (n<5)
        aspect_ratio=PLASTIC_NUMBER;
    else if (n<8)
        aspect_ratio=PLASTIC_NUMBER_2;
    else
        aspect_ratio=PLASTIC_NUMBER_3;
    double step1=aspect_ratio/(double)(n1-1),step2=aspect_ratio/(double)(n2-1),xpos=0.0,ypos=1.0;
    x.resize(node_count());
    for (ivector_iter it=p1.begin();it!=p1.end();++it) {
        point &p=x[*it];
        p.resize(2);
        p.front()=xpos;
        p.back()=ypos;
        xpos+=step1;
    }
    xpos=ypos=0.0;
    for (ivector_iter it=p2.begin();it!=p2.end();++it) {
        point &p=x[*it];
        p.resize(2);
        p.front()=xpos;
        p.back()=ypos;
        xpos+=step2;
    }
}

/* construct the plane dual of a planar graph with the given faces with time complexity O(n),
* each face must be a list of vertex indices */
void graphe::make_plane_dual(const ivectors &faces) {
    this->clear();
    set_directed(false);
    reserve_nodes(faces.size());
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,faces.size());
        add_nodes(labels);
    } else add_nodes(faces.size());
    int nc=0;
    for (ivectors_iter it=faces.begin();it!=faces.end();++it) {
        for (ivector_iter jt=it->begin();jt!=it->end();++jt) {
            if (*jt>nc)
                nc=*jt;
        }
    }
    edgemap emap(++nc);
    map<int,int>::const_iterator et;
    int i,j,v,w,f,n;
    for (ivectors_iter it=faces.begin();it!=faces.end();++it) {
        f=it-faces.begin();
        const ivector &face=*it;
        n=face.size();
        for (int k=0;k<n;++k) {
            i=face[k];
            j=face[(k+1)%n];
            v=i<j?i:j;
            w=i<j?j:i;
            map<int,int> &m=emap[v];
            if ((et=m.find(w))==m.end())
                m[w]=f;
            else add_edge(f,et->second);
        }
    }
}

/* return the number of different colors adjacent to the i-th vertex in a partially colored graph */
int graphe::saturation_degree(const vertex &v,std::set<int> &colors) const {
    int i,c;
    colors.clear();
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        i=*it;
        const vertex &w=node(i);
        if ((c=w.color())>0)
            colors.insert(c);
    }
    return colors.size();
}

/* return the degree of v in the uncolored subgraph of this graph */
int graphe::uncolored_degree(const vertex &v) const {
    int deg=0;
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        if (node(*it).color()==0)
            ++deg;
    }
    return deg;
}

/* return true iff this graph has at least one uncolored (white) vertex */
bool graphe::is_partially_colored() const {
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        if (it->color()==0)
            return true;
    }
    return false;
}

/* heuristic algorithm by Brelaz for nearly optimal vertex coloring, complexity O(n*m) */
void graphe::dsatur() {
    int col,i,sat,maxsat,deg,maxdeg=0;
    std::set<int> colors,maxcolors;
    ivector indices;
    node_iter nt=nodes.begin(),ntend=nodes.end();
    for (i=0;nt!=ntend;++nt,++i) {
        if (nt->color()==0)
            indices.push_back(i);
    }
    if (indices.empty())
        return;
    ivector::iterator pos,it,itend;
    std::set<int>::const_iterator jt,jtend;
    do {
        maxsat=-1;
        i=-1;
        for (it=indices.begin(),itend=indices.end();it!=itend;++it) {
            const vertex &v=node(*it);
            if ((sat=saturation_degree(v,colors))>maxsat ||
                    (sat==maxsat && (deg=uncolored_degree(v))>maxdeg)) {
                if (sat==maxsat)
                    maxdeg=deg;
                maxsat=sat;
                maxcolors=colors;
                i=*it;
                pos=it;
            }
        }
        if (i>=0) {
            col=1;
            for (jt=maxcolors.begin(),jtend=maxcolors.end();jt!=jtend;++jt) {
                if (*jt==col)
                    ++col;
                else break;
            }
            node(i).set_color(col);
            indices.erase(pos);
        }
    } while (i>=0);
}

/* return the total number of different nonzero vertex colors in this graph */
int graphe::color_count() const {
    std::set<int> colors;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it) {
        colors.insert(it->color());
    }
    return colors.size();
}

/* return true iff the vertices can be colored using at most k different colors,
* the vertices will be colored after the function returns */
bool graphe::is_vertex_colorable(int k) {
    assert(k>=0);
    if (k==0) {
        uncolor_all_nodes();
        return true;
    }
    if (k>node_count()) {
        message("Warning: there are more colors than vertices");
        return false;
    }
    /* try greedy coloring first (linear time), use random order of vertices */
    ivector sigma=rand_permu(node_count());
    if (greedy_vertex_coloring(sigma)<=k)
        return true;
    /* next try dsatur algorithm (quadratic time) */
    uncolor_all_nodes();
    dsatur();
    if (color_count()<=k)
        return true;
    /* finally resort to solving MIP problem */
    return exact_vertex_coloring(k)!=0;
}

/* return the lower and the upper bound for chromatic number (inclusive) */
graphe::ipair graphe::chromatic_number_bounds() {
    /* the lower bound is the size of maximum clique */
    ivector clique;
    ostergard ost(this,3.0); // with timeout
    int lb=ost.maxclique(clique);
    /* the upper bound is given by heuristic dsatur algorithm */
    uncolor_all_nodes();
    for (ivector_iter it=clique.begin();it!=clique.end();++it) {
        node(*it).set_color(it-clique.begin()+1);
    }
    dsatur();
    int ub=color_count();
    return make_pair(lb,ub);
}

/* store custom vertex coordinates */
void graphe::store_layout(const layout &x) {
    assert(int(x.size())>=node_count());
    for (layout_iter it=x.begin();it!=x.end();++it) {
        vertex &v=node(it-x.begin());
        v.set_attribute(_GT_ATTRIB_POSITION,point2gen(*it));
    }
}

/* retrieve the previously stored layout, return true iff successful */
bool graphe::has_stored_layout(layout &x) const {
    x.resize(node_count());
    attrib_iter ait;
    int dim=0,i=0;
    node_iter it=nodes.begin(),itend=nodes.end();
    for (;it!=itend;++it,++i) {
        const attrib &attr=it->attributes();
        point &p=x[i];
        if ((ait=attr.find(_GT_ATTRIB_POSITION))==attr.end() ||
                !gen2point(ait->second,p,ctx) || (dim>0 && int(p.size())!=dim))
            return false;
        if (dim==0)
            dim=p.size();
    }
    return true;
}

bool graphe::bipartite_matching_bfs(ivector &dist,int sg) {
    assert(node_queue.empty() && sg>=0);
    int u,v;
    node_iter it=nodes.begin(),itend=nodes.end();
    ivector_iter jt,jtend;
    for (u=0;it!=itend;++it,++u) {
        if (it->subgraph()!=sg || it->color()!=1)
            continue;
        if (it->number()==0) {
            dist[u+1]=0;
            node_queue.push(u);
        } else dist[u+1]=rand_max2;
    }
    dist.front()=rand_max2;
    while (!node_queue.empty()) {
        u=node_queue.front();
        node_queue.pop();
        if (dist[u+1]<dist.front()) {
            vertex &U=node(u);
            for (jt=U.neighbors().begin(),jtend=U.neighbors().end();jt!=jtend;++jt) {
                v=*jt;
                vertex &V=node(v);
                if (V.subgraph()==sg && dist[V.number()]==rand_max2) {
                    dist[V.number()]=dist[u+1]+1;
                    node_queue.push(V.number()-1);
                }
            }
        }
    }
    return dist.front()!=rand_max2;
}

bool graphe::bipartite_matching_dfs(int u,ivector &dist,int sg) {
    assert(sg>=0);
    if (u>0) {
        int v;
        vertex &U=node(u-1);
        for (ivector_iter it=U.neighbors().begin();it!=U.neighbors().end();++it) {
            v=*it;
            vertex &V=node(v);
            if (V.subgraph()==sg && dist[V.number()]==dist[u]+1) {
                if (bipartite_matching_dfs(V.number(),dist,sg)) {
                    V.set_number(u);
                    U.set_number(v+1);
                    return true;
                }
            }
        }
        dist[u]=rand_max2;
        return false;
    }
    return true;
}

/* find maximum matching in a bipartite graph using Hopcroft-Karp algorithm */
void graphe::bipartite_matching(const ivector &p1,const ivector &p2,ipairs &matching,int sg,gt_conn_check cc) {
    int cmp,s;
    switch (cc) {
    case _GT_CC_CONNECTED:
        break;
    case _GT_CC_COMPONENTS_ARE_SUBGRAPHS:
    case _GT_CC_FIND_COMPONENTS:
        matching.clear();
        s=(cc==_GT_CC_COMPONENTS_ARE_SUBGRAPHS || sg<0)?1:1+max_subgraph_index();
        cmp=(cc==_GT_CC_FIND_COMPONENTS)?connected_components_to_subgraphs(sg):max_subgraph_index();
        for (;s<=cmp;++s) {
            int sz=subgraph_size(s);
            assert(sz!=0);
            if (sz==1)
                ;
            else if (sz<=3) {
                int i=first_vertex_from_subgraph(s);
                assert(!node(i).neighbors().empty());
                int j=node(i).neighbors().front();
                matching.push_back(make_pair(std::min(i,j),std::max(i,j)));
            } else bipartite_matching(p1,p2,matching,s,_GT_CC_CONNECTED);
            if (cc==_GT_CC_FIND_COMPONENTS)
                renumber_subgraph(s,sg);
        }
        return;
    default: assert(false);
    }
    // sg >= 0
    ivector_iter it,itend;
    for (it=p1.begin(),itend=p1.end();it!=itend;++it) {
        if (node(*it).subgraph()==sg)
            set_node_color(*it,1);
    }
    for (it=p2.begin(),itend=p2.end();it!=itend;++it) {
        if (node(*it).subgraph()==sg)
            set_node_color(*it,2);
    }
    for (vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
        it->set_number(0);
    }
    ivector dist(node_count()+1);
    int count=0,u,v;
    node_iter nt,ntend=nodes.end();
    while (bipartite_matching_bfs(dist,sg)) {
        for (u=0,nt=nodes.begin();nt!=ntend;++nt,++u) {
            if (nt->subgraph()!=sg || nt->color()!=1)
                continue;
            if (nt->number()==0 && bipartite_matching_dfs(u+1,dist,sg))
                ++count;
        }
    }
    /* extract matching */
    for (u=0,nt=nodes.begin();nt!=ntend;++nt,++u) {
        if (nt->subgraph()!=sg || nt->color()!=1)
            continue;
        v=nt->number();
        if (v>0) {
            matching.push_back(make_pair(std::min(u,v-1),std::max(u,v-1)));
            --count;
        }
    }
    assert(count==0);
}

/* find a maximum/minimum weighted matching in bipartite graph, return true iff successful */
bool graphe::weighted_bipartite_matching(const ivector &p1,const ivector &p2,ipairs &matching,
        bool minimize,double eps,int sg,gt_conn_check cc) {
    int cmp,s;
    switch (cc) {
    case _GT_CC_CONNECTED:
        break;
    case _GT_CC_COMPONENTS_ARE_SUBGRAPHS:
    case _GT_CC_FIND_COMPONENTS:
        assert(supports_attributes());
        if (is_null() || is_directed() || !is_weighted())
            return false;
        if (is_empty())
            return true;
        matching.clear();
        s=(cc==_GT_CC_COMPONENTS_ARE_SUBGRAPHS || sg<0)?1:1+max_subgraph_index();
        cmp=(cc==_GT_CC_FIND_COMPONENTS)?connected_components_to_subgraphs(sg):max_subgraph_index();
        for (;s<=cmp;++s) {
            int sz=subgraph_size(s);
            assert(sz>0);
            if (sz==1)
                continue;
            if (!weighted_bipartite_matching(p1,p2,matching,minimize,eps,s,_GT_CC_CONNECTED))
                return false;
            if (cc==_GT_CC_FIND_COMPONENTS)
                renumber_subgraph(s,sg);
        }
        return true;
    default: assert(false);
    }
    //sg >= 0
    int n1,n2,i,j;
    matrice W;
    weight_matrix(W);
    gen f;
    ivector P,Q; // filtered partitions
    P.reserve(p1.size());
    Q.reserve(p2.size());
    for (ivector_iter it=p1.begin();it!=p1.end();++it) {
        if (node(*it).subgraph()==sg)
            P.push_back(*it);
    }
    n1=P.size();
    for (ivector_iter it=p2.begin();it!=p2.end();++it) {
        if (node(*it).subgraph()==sg)
            Q.push_back(*it);
    }
    n2=Q.size();
    for (ivector_iter it=P.begin();it!=P.end();++it) {
        for (ivector_iter jt=node(*it).neighbors().begin();jt!=node(*it).neighbors().end();++jt) {
            if (node(*jt).subgraph()!=sg)
                continue;
            if(find(Q.begin(),Q.end(),*jt)==Q.end())
                return false; // not bipartite
            const gen &w=W[*it][*jt];
            if (w.is_integer())
                continue;
            if (!is_real_number(f,ctx)) {
                message("Error: edge weights must be real");
                return false;
            }
            if (!is_rational(w))
                W[*it]._VECTptr->at(*jt)=to_real_number(f,ctx);
        }
    }
#if 0 //def HAVE_LIBGLPK
    /* integralize the weight matrix */
    if (!is_integer_matrice(W,true))
        W=*_apply(makesequence(at_round,divvecteur(W,eps)),ctx)._VECTptr;
    if (minimize) W=multvecteur(-1,W);
    gen minw=_min(_min(W,ctx),ctx);
    /* create the assignment problem as a GLPK graph */
    glp_graph *G;
    G=glp_create_graph(sizeof(wbm_v_data),sizeof(wbm_a_data));
    int start1=glp_add_vertices(G,n1);
    int start2=glp_add_vertices(G,n2);
    glp_arc *a;
    glp_vertex *v;
    for (i=0;i<std::max(n1,n2);++i) {
        if (i<n1) ((wbm_v_data*)(G->v)[start1+i]->data)->set=0;
        if (i<n2) ((wbm_v_data*)(G->v)[start2+i]->data)->set=1;
    }
    for (ivector_iter it=P.begin();it!=P.end();++it) {
        i=it-P.begin();
        for (ivector_iter jt=node(*it).neighbors().begin();jt!=node(*it).neighbors().end();++jt) {
            if (sg>=0 && node(*jt).subgraph()!=sg)
                continue;
            j=find(Q.begin(),Q.end(),*jt)-Q.begin();
            a=glp_add_arc(G,U_start+i,V_start+j);
            ((wbm_a_data*)a->data)->cost=(W[*it][*jt]-minw).val;
        }
    }
    /* solve the problem using out-of-kilter algorithm */
    int res=glp_asnprob_okalg(GLP_ASN_MMP,G,
                offsetof(wbm_v_data,set),offsetof(wbm_a_data,cost),NULL,offsetof(wbm_a_data,x));
    if (res!=0) {
        if (verbose)
            message("Error: GLPK solver error %d",res);
        glp_delete_graph(G);
        return false;
    }
    /* retrieve the matching */
    for (i=0;i<n1;++i) {
        glp_vertex *v=G->v[U_start+i];
        for (a=v->out;a!=NULL;a=a->t_next) {
            if (((wbm_a_data*)a->data)->x==0)
                continue;
            j=a->head->i;
            assert(j-V_start>=0);
            if (has_edge(U[i],V[j-V_start]))
                matching.push_back(make_pair(U[i],Q[j-V_start]));
        }
    }
    glp_delete_graph(G);
    return true;
#else
    /* convert the assignment problem to transportation problem */
    if (!minimize) W=multvecteur(-1,W);
    gen minw=_min(_min(W,ctx),ctx);
    int n=std::max(n1,n2);
    vecteur sd(n,1);
    gen M=identificateur(" M");
    _purge(M,ctx);
    tprob prob(sd,sd,M,ctx);
    matrice C=*_matrix(makesequence(n,n,0),ctx)._VECTptr;
    for (ivector_iter it=P.begin();it!=P.end();++it) {
        i=it-P.begin();
        for (ivector_iter jt=Q.begin();jt!=Q.end();++jt) {
            j=jt-Q.begin();
            C[i]._VECTptr->at(j)=has_edge(*it,*jt)?W[*it][*jt]-minw:M;
        }
    }
    matrice sol;
    prob.solve(C,sol);
    /* retrieve the matching */
    for (ivector_iter it=P.begin();it!=P.end();++it) {
        i=it-P.begin();
        for (ivector_iter jt=Q.begin();jt!=Q.end();++jt) {
            j=jt-Q.begin();
            if (giac::is_one(sol[i][j])) {
                if (!has_edge(P[i],Q[j]))
                    return false;
                matching.push_back(make_pair(P[i],Q[j]));
                break;
            }
        }
    }
    return true;
#endif
}

/* construct the line graph of this graph */
void graphe::line_graph(graphe &G,ipairs &E) const {
    assert(supports_attributes() || !G.supports_attributes());
    if (E.empty())
        get_edges_as_pairs(E);
    G.clear();
    G.reserve_nodes(E.size());
    bool isdir=G.is_directed();
    G.set_directed(isdir);
    if (G.supports_attributes()) {
        vecteur labels;
        gen label;
        for (ipairs_iter it=E.begin();it!=E.end();++it) {
            label=_cat(makesequence(node_label(it->first),str2gen("-",true),node_label(it->second)),ctx);
            labels.push_back(label);
        }
        G.add_nodes(labels);
    } else G.add_nodes(E.size());
    int i,j;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        i=it-E.begin();
        for (ipairs_iter jt=(isdir?E.begin():it+1);jt!=E.end();++jt) {
            j=jt-E.begin();
            if (i==j) continue;
            if ((isdir && it->second==jt->first) || (!isdir && edges_incident(*it,*jt)))
                G.add_edge(i,j);
        }
    }
}

/* construct the transitive closure of this graph */
void graphe::transitive_closure(graphe &G,bool weighted) {
    int n=node_count();
    bool isdir=is_directed(),iswei=is_weighted();
    G.clear();
    G.set_directed(isdir);
    G.set_weighted(weighted);
    G.reserve_nodes(n);
    G.add_nodes(vertices());
    if (weighted) {
        matrice m;
        gen wgh;
        ivector dist,J(n-1);
        int cnt;
        if (iswei)
            allpairs_distance(m);
        for (int i=0;i<n;++i) {
            if (!iswei) {
                cnt=0;
                for (int k=0;k<n;++k) {
                    if (k!=i)
                        J[cnt++]=k;
                }
                distance(i,J,dist);
            }
            for (int j=isdir?0:i+1;j<n;++j) {
                if (i==j)
                    continue;
                if ((iswei && !is_inf(wgh=m[i][j])) ||
                        (!iswei && is_positive(wgh=dist[j<i?j:j-1],ctx)))
                    G.add_edge(i,j,wgh);
            }
        }
    } else {
        for (int i=0;i<n;++i) {
            dfs(i,false);
            for (int j=isdir?0:i+1;j<n;++j) {
                if (i==j)
                    continue;
                if (node(j).is_visited())
                    G.add_edge(i,j);
            }
        }
    }
}

/* return true iff this graph is isomorphic to other, also obtain an isomorphism */
bool graphe::is_isomorphic(graphe &other,map<int,int> &isom,bool use_colors) {
    if (is_directed()!=other.is_directed())
        return false;
    int n=node_count(),sz;
    if (other.node_count()!=n || edge_count()!=other.edge_count())
        return false;
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
    int *adj1=to_array(sz,use_colors);
    int *adj2=other.to_array(sz,use_colors);
    int *sigma=new int[n];
    pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
    bool res=nautywrapper_is_isomorphic(is_directed()?1:0,n,adj1,adj2,sigma)!=0;
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    if (res) {
        /* obtain the isomorphism */
        isom.clear();
        for (int i=0;i<n;++i) {
            isom[i]=sigma[i];
        }
    }
    delete[] adj1;
    delete[] adj2;
    delete[] sigma;
    return res;
#else
    ivector d1=degree_sequence(),d2=other.degree_sequence();
    std::sort(d1.begin(),d1.end());
    std::sort(d2.begin(),d2.end());
    if (d1!=d2 || connected_component_count()!=other.connected_component_count())
        return false;
    /* use Ullmann's algorithm */
    ivectors uv;
    subgraph_isomorphism(other,1,true,uv);
    if (uv.empty())
        return false;
    else {
        isom.clear();
        for (int i=0;i<n;++i) {
            isom[uv.front()[i]]=i;
        }
        return true;
    }
#endif
}

/* convert a sequence of nd digits to nonnegative integer */
int digits2int(char *digits,int nd) {
    int n=1,res=0;
    for (int i=nd;i-->0;) {
        res+=*(digits+i)*n;
        n*=10;
    }
    return res;
}

/* return the set of generators of the automorphism group of this graph */
gen graphe::aut_generators() const {
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
    int n=node_count(),ofs=array_start(ctx),sz;
    vecteur out(0);
    if (n>0) {
        int *adj=to_array(sz,true);
        char buf[L_tmpnam];
        bool has_temp_file=(tmpnam(buf)!=NULL);
        string tmpname(has_temp_file?buf:"nauty-autg");
#ifdef _WIN32
        string tmpdir=getenv("TEMP")?getenv("TEMP"):"c:\\Users\\Public";
        tmpname=tmpdir+tmpname;
#endif
        FILE *f=fopen(tmpname.c_str(),"wb+");
        if (f==NULL) {
            message ("Error: failed to create temporary file");
            return undef;
        }
        pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
        nautywrapper_aut_generators(is_directed()?1:0,n,adj,f);
        pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
        /* parse the generators and output them as a sequence of
     * permutations, each in form of a list of disjoint cycles */
        int sz=ftell(f);
        rewind(f);
        int i=0,nd=0;
        char c;
        vecteur perm,cycle;
        char digits[32];
        for (int cnt=0;cnt<sz;++cnt) {
            c=fgetc(f);
            if (c=='(')
                cycle.clear();
            else if (c==')') {
                cycle.push_back(ofs+digits2int(digits,nd));
                perm.push_back(cycle);
                nd=0;
            } else if (c=='\n') {
                out.push_back(gen(perm,_LIST__VECT));
                perm.clear();
            } else if (c>47 && c<58)
                digits[nd++]=c-48;
            else if (c==32) {
                cycle.push_back(ofs+digits2int(digits,nd));
                nd=0;
            }
            ++i;
        }
        fclose(f);
        if (remove(tmpname.c_str())!=0)
            message("Warning: failed to remove temporary file");
        delete[] adj;
    }
    return gen(out,_LIST__VECT);
#else
    return gensizeerr("nauty library is required for finding graph automorphisms");
#endif
}

/* return the canonical labeling of this graph as a permutation */
bool graphe::canonical_labeling(ivector &lab) const {
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
    int n=node_count(),sz;
    if (n==0)
        return false;
    int *adj=to_array(sz,true);
    int *clab=new int[n];
    pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
    nautywrapper_canonical(is_directed()?1:0,n,adj,clab,NULL,NULL);
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
    lab.resize(n);
    for (int i=0;i<n;++i) {
        lab[i]=clab[i];
    }
    delete[] adj;
    delete[] clab;
    return true;
#else
    return false;
#endif
}

/* construct the closure of this graph and store it to G, complexity O(n^3) */
bool graphe::bondy_chvatal_closure(graphe &G,ivector &d) {
    underlying(G);
    int n=node_count(),di;
    bool yes=false,have_pair;
    do {
        have_pair=false;
        for (int i=0;i<n && !have_pair;++i) {
            di=d[i];
            for (int j=i+1;j<n;++j) {
                if (!G.has_edge(i,j) && di+d[j]>=n) {
                    G.add_edge(i,j);
                    ++d[i]; ++d[j];
                    yes=have_pair=true;
                    break;
                }
            }
        }
    } while (have_pair);
    return yes;
}

/* return 1 if the graph is Hamiltonian, 0 if it is not and -1 if inconclusive */
int graphe::hamcond(bool make_closure) {
    /* test biconnectivity, complexity O(n) */
    if (!is_biconnected())
        return 0;
    int mindeg=rand_max2,deg,n=node_count(),m=edge_count(),i;
    ivector d=degree_sequence();
    /* Dirac criterion, complexity O(n) */
    node_iter it=nodes.begin(),itend=nodes.end();
    for (i=0;it!=itend;++it,++i) {
        if ((deg=d[i]=it->neighbors().size())<mindeg)
            mindeg=deg;
    }
    if (2*mindeg>=n)
        return 1;
    ivector V1,V2,v;
    /* unbalanced bipartite graph is not Hamiltonian, complexity O(n+m) */
    if (is_bipartite(V1,V2) && V1.size()!=V2.size())
        return 0;
    /* Ore criterion, complexity O(n^2) */
    bool yes=true;
    for (int i=0;i<n && yes;++i) {
        deg=d[i];
        for (int j=i+1;j<n;++j) {
            if (!has_edge(i,j) && deg+d[j]<n) {
                yes=false;
                break;
            }
        }
    }
    if (yes) return 1;
    if (make_closure) {
        /* apply the Bondy-Chvátal theorem: the graph is Hamiltonian if
           and only if its closure is Hamiltonian, complexity O(n^2) */
        graphe G(ctx);
        if (bondy_chvatal_closure(G,d))
            return G.hamcond(false);
    }
    if (double(m)/(double(n*n)>1.0-2.0/(1.0+std::exp(.0015*double(n))))) {
        /* Nash-Williams criterion */
        int inum=maximum_independent_set(v);
        if (3*mindeg>=std::max(n+2,3*inum))
            return 1;
    }
    return -1;
}

/* return true iff the graph is Hamiltonian */
bool graphe::is_hamiltonian(ivector &hc) {
    hc.clear();
    if (!is_directed()) {
      switch (hamcond()) {
      case 0: return false;
      case 1: return true;
      default: break;
      }
    } else {
        ivectors components;
        strongly_connected_components(components);  
        if (components.size()>1) return false; // graph is not strongly connected
        ivector dv=degree_sequence();
        bool isham=true;
        int n=node_count();
        for (ivector_iter it=dv.begin();it!=dv.end();++it) {
            if (*it<n) {
               isham=false;
               break;
            }
        }
        if (isham) return true;
        isham=true;
        for (int i=0;isham && i<n;++i) {
            for (int j=0;j<n;++j) {
                if (i==j) continue;
                if (!has_edge(i,j) && !has_edge(j,i) && dv[i]+dv[j]<2*n-1) {
                    isham=false;
                    break;
                }
            }
        }
        if (isham) return true;
    }
    return hamcycle(hc);
}

bool graphe::hamcycle_recurse(ivector &path,int pos) {
    if (pos==node_count()) return has_edge(path[pos-1],path[0]);
    const vertex &v=node(path[pos-1]);
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        vertex &w=node(*it);
        if (w.is_visited()) continue;
        path[pos]=*it;
        w.set_visited(true);
        if (hamcycle_recurse(path,pos+1)) return true;
        w.set_visited(false);
    }
    return false;
}

bool graphe::hamcycle(ivector &path) {
    int n=node_count();
    path.resize(n+1);
    unvisit_all_nodes();
    path[0]=0;
    nodes[0].set_visited(true);
    if (!hamcycle_recurse(path,1)) return false;
    path[n]=path[0];
    return true;
}

/*
* TSP CLASS IMPLEMENTATION
*/

#ifdef HAVE_LIBGLPK

/* MINIMAL CUT implementation, adapted from glpk/examples/tsp (maxflow.c, mincut.c and main.c)
 * Originally written by Andrew Makhorin <mao@gnu.org>
 */

/* compute maximal flow */
int graphe::tsp::max_flow(int nn,int nedg,const ivector &beg,
        const ivector &end,const ivector &cap,int s,int t,
        ivector &x) {
    int k;
    assert(nn>=2);
    assert(nedg>=0);
    assert(1<=s && s<=nn);
    assert(1<=t && t<=nn);
    assert(s!=t);
    for (k=0; k<nedg; k++) {
        assert(1<=beg[k] && beg[k]<end[k] && end[k]<=nn && cap[k]>0);
    }
    glp_prob *lp;
    glp_smcp smcp;
    int i,nz,flow,*rn,*cn;
    double temp,*aa;
    lp=glp_create_prob();
    glp_add_rows(lp,nn);
    for (i=1;i<=nn;++i) {
        glp_set_row_bnds(lp,i,GLP_FX,0.0,0.0);
    }
    glp_add_cols(lp,nedg+1);
    for (k=0; k<nedg; k++) {
        glp_set_col_bnds(lp,k+1,GLP_DB,-cap[k],+cap[k]);
    }
    glp_set_col_bnds(lp,nedg+1,GLP_FR,0.0,0.0);
    rn=new int[1+2*(nedg+1)];
    cn=new int[1+2*(nedg+1)];
    aa=new double[1+2*(nedg+1)];
    nz=0;
    for (k=0;k<nedg;k++) {
        rn[++nz]=beg[k]; cn[nz]=k+1; aa[nz]=-1.0;
        rn[++nz]=end[k]; cn[nz]=k+1; aa[nz]=+1.0;
    }
    rn[++nz]=t; cn[nz]=nedg+1; aa[nz]=-1.0;
    rn[++nz]=s; cn[nz]=nedg+1; aa[nz]=+1.0;
    assert(nz==2*(nedg+1));
    glp_load_matrix(lp,nz,rn,cn,aa);
    delete[] rn;
    delete[] cn;
    delete[] aa;
    glp_set_obj_dir(lp,GLP_MAX);
    glp_set_obj_coef(lp,nedg+1,1.0);
    int old_term=glp_term_out(GLP_OFF);
    glp_adv_basis(lp,0);
    glp_term_out(old_term);
    glp_init_smcp(&smcp);
    smcp.msg_lev=GLP_MSG_OFF;
    assert(glp_simplex(lp,&smcp)==0);
    assert(glp_get_status(lp)==GLP_OPT);
    for (k=0;k<nedg;k++) {
        temp=glp_get_col_prim(lp,k+1);
        x[k]=(int)floor(temp+.5);
        assert(std::abs(x[k]-temp)<=1e-6);
    }
    temp=glp_get_col_prim(lp,nedg+1);
    flow=(int)floor(temp+.5);
    assert(std::abs(flow-temp)<=1e-6);
    glp_delete_prob(lp);
    return flow;
}

/* find min (s,t)-cut for known max flow */
int graphe::tsp::min_st_cut(int nn,int nedg,const ivector &beg,
        const ivector &end,const ivector &cap,int s,int t,
        const ivector &x,ivector &cut) {
    int i,j,k,p,q,temp;
    ivector &head1=mincut_data[18]; if (int(head1.size())<nn) head1.resize(nn);
    ivector &head2=mincut_data[19]; if (int(head2.size())<nn) head2.resize(nn);
    ivector &next1=mincut_data[20]; if (int(next1.size())<nn) next1.resize(nn);
    ivector &next2=mincut_data[21]; if (int(next2.size())<nn) next2.resize(nn);
    for (i=0;i<nn;i++) head1[i]=head2[i]=0;
    for (k=0;k<nedg;k++) {
        i=beg[k]-1; next1[k]=head1[i]; head1[i]=k+1;
        j=end[k]-1; next2[k]=head2[j]; head2[j]=k+1;
    }
    ivector &list=mincut_data[4]; list.resize(nn);
    for (i=0;i<nn;i++) cut[i]=0;
    p=q=0; list[0]=s; cut[s-1]=1;
    while (p<=q) {
        i=list[p++];
        for (k=head1[i-1];k!=0;k=next1[k-1]) {
            j=end[k-1];
            assert(beg[k-1]==i);
            if (cut[j-1]==0 && x[k-1]<+cap[k-1]) {
                list[++q]=j; cut[j-1]=1;
            }
        }
        for (k=head2[i-1];k!=0;k=next2[k-1]) {
            j=beg[k-1];
            assert(end[k-1]==i);
            if (cut[j-1]==0 && x[k-1]>-cap[k-1]) {
                list[++q]=j; cut[j-1]=1;
            }
        }
    }
    assert(!cut[t-1]);
    temp=0;
    for (k=0;k<nedg;k++) {
        i=beg[k],j=end[k];
        if ((cut[i-1] && !cut[j-1]) || (!cut[i-1] && cut[j-1]))
            temp+=cap[k];
    }
    return temp;
}

/* find min cut with Stoer and Wagner algorithm */
int graphe::tsp::minimal_cut(int nn,int nedg,const ivector &beg,
                             const ivector &end,const ivector &cap,ivector &cut) {
    int k;
    assert(nn>=2);
    assert(nedg>=0);
    for (k=0;k<nedg;k++) {
        assert(1<=beg[k] && beg[k]<end[k] && end[k]<=nn);
        assert(cap[k]>0);
    }
    int i,j,min_cut,flow,temp,I,J,K,S,T,DEG,NV,NE;
    ivector &head1=mincut_data[0]; if (int(head1.size())<nn) head1.resize(nn);
    ivector &head2=mincut_data[1]; if (int(head2.size())<nn) head2.resize(nn);
    ivector &next1=mincut_data[2]; if (int(next1.size())<nedg) next1.resize(nedg);
    ivector &next2=mincut_data[3]; if (int(next2.size())<nedg) next2.resize(nedg);
    for (i=0;i<nn;i++) head1[i]=head2[i]=0;
    for (k=0;k<nedg;k++) {
        i=beg[k]; next1[k]=head1[i-1]; head1[i-1]=k+1;
        j=end[k]; next2[k]=head2[j-1]; head2[j-1]=k+1;
    }
    NV=nn;
    ivector &HEAD=mincut_data[4]; if (int(HEAD.size())<nn) HEAD.resize(nn);
    ivector &NEXT=mincut_data[5]; if (int(NEXT.size())<nn) NEXT.resize(nn);
    ivector &NUMB=mincut_data[6]; if (int(NUMB.size())<nn) NUMB.resize(nn);
    for (i=0;i<nn;i++) {
        HEAD[i]=i+1; NEXT[i]=0; NUMB[i]=i+1;
    }
    ivector &BEG=mincut_data[7]; if (int(BEG.size())<nedg) BEG.resize(nedg);
    ivector &END=mincut_data[8]; if (int(END.size())<nedg) END.resize(nedg);
    ivector &CAP=mincut_data[9]; if (int(CAP.size())<nedg) CAP.resize(nedg);
    ivector &X=mincut_data[10]; if (int(X.size())<nedg) X.resize(nedg);
    ivector &ADJ=mincut_data[11]; if (int(ADJ.size())<nn) ADJ.resize(nn);
    ivector &SUM=mincut_data[12]; if (int(SUM.size())<nn) SUM.resize(nn);
    ivector &CUT=mincut_data[13]; if (int(CUT.size())<nn) CUT.resize(nn);
    min_cut=rand_max2;
    while (NV>1) {
        for (I=0;I<NV;I++) SUM[I]=0.0;
        for (I=1;I<=NV;I++) {
            DEG=0;
            for (i=HEAD[I-1];i!=0;i=NEXT[i-1]) {
                for (k=head1[i-1];k!=0;k=next1[k-1]) {
                    j=end[k-1];
                    J=NUMB[j-1];
                    if (I>=J)
                        continue;
                    if (SUM[J-1]==0.0)
                        ADJ[DEG++]=J;
                    assert(cap[k-1]>0.0);
                    SUM[J-1]+=cap[k-1];
                }
                for (k=head2[i-1];k!=0;k=next2[k-1]) {
                    j=beg[k-1];
                    J=NUMB[j-1];
                    if (I>=J)
                        continue;
                    if (SUM[J-1]==0.0)
                        ADJ[DEG++]=J;
                    assert(cap[k-1]>0.0);
                    SUM[J-1]+=cap[k-1];
                }
            }
            NE=0;
            for (K=0;K<DEG;K++) {
                assert(NE<nedg);
                J=ADJ[K];
                BEG[NE]=I; END[NE]=J; CAP[NE]=SUM[J-1];
                SUM[J-1]=0.0;
                ++NE;
            }
        }
        S=1; T=NV;
        flow=max_flow(NV,NE,BEG,END,CAP,S,T,X);
        if (min_cut>flow)
        {  min_cut=flow;
            temp=min_st_cut(NV,NE,BEG,END,CAP,S,T,X,CUT);
            assert(flow==temp);
            for (i=0;i<nn;i++) cut[i]=CUT[NUMB[i]-1];
            if (min_cut==0)
                break;
        }
        assert(T==NV);
        for (i=HEAD[T-1];i!=0;i=NEXT[i-1]) NUMB[i-1]=S;
        i=HEAD[S-1];
        assert(i!=0);
        while (NEXT[i-1]!=0)
            i=NEXT[i-1];
        NEXT[i-1]=HEAD[T-1];
        NV--;
    }
    return min_cut;
}

/* TSP constructor */
graphe::tsp::tsp(graphe *gr,double gap_tolerance,bool is_verbose) {
    assert(!gr->is_directed());
    G=gr;
    verbose=is_verbose;
    gap_tol=gap_tolerance;
    isweighted=G->is_weighted();
    ipairs E;
    G->get_edges_as_pairs(E);
    edges_comparator comp(G);
    std::sort(E.begin(),E.end(),comp); // sort edges by weight (ascending)
    nv=G->node_count();
    ne=E.size();
    tour.reserve(nv);
    obj.resize(ne);
    xev.resize(ne);
    can_branch.resize(ne);
    coeff=new double[ne+1];
    indices=new int[ne+1];
    visited=new bool[nv];
    arcs=new arc[ne];
    cancellable=false;
    int i;
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        i=it-E.begin();
        arc &a=arcs[i];
        a.head=it->second;
        a.tail=it->first;
        loc_map[a.tail][a.head]=i;
        weight_map[a.tail][a.head]=isweighted?to_real_number(G->weight(*it),G->giac_context()).to_double(G->giac_context()):1.0;
    }
}

/* TSP destructor */
graphe::tsp::~tsp() {
    delete[] coeff;
    delete[] indices;
    delete[] visited;
    delete[] arcs;
}

/* create an edge or arc */
graphe::ipair graphe::tsp::make_edge(int i,int j) const {
    return make_pair(i<j?i:j,i<j?j:i);
}

/* return the subgraph index of edge e */
int graphe::tsp::edge_index(const ipair &e) {
    return loc_map[e.first][e.second];
}

/* formulate TSP as MIP, initially without any subtour elimination constraints */
void graphe::tsp::formulate_mip() {
    ivectors rows;
    int nonzeros=0,i,j,k,l,cnt=0,nrows,ncols;
    for (j=0;j<ne;++j) {
        const arc &a=arcs[j];
        obj[j]=weight(a.tail,a.head);
    }
    for (j=0;j<nv;++j) {
        rows.resize(rows.size()+1);
        ivector &row=rows.back();
        row.push_back(j);
        for (i=0;i<nv;++i) {
            if (i!=j && G->has_edge(i,j)) {
                row.push_back(i);
                ++nonzeros;
            }
        }
    }
    ncols=ne;
    nrows=rows.size();
    if (isweighted)
        nonzeros+=ne;
    int *ia=new int[nonzeros+1];
    int *ja=new int[nonzeros+1];
    double *ar=new double[nonzeros+1];
    for (ivectors_iter it=rows.begin();it!=rows.end();++it) {
        const ivector &row=*it;
        k=it-rows.begin();
        i=row.front();
        for (ivector_iter rt=row.begin()+1;rt!=row.end();++rt) {
            j=*rt;
            l=edge_index(make_edge(i,j));
            assert(l<ne);
            ia[++cnt]=k+1; ja[cnt]=l+1; ar[cnt]=1;
        }
    }
    /* build MIP */
    glp_add_rows(mip,nrows);
    double rh=2.0;
    for (i=0;i<nrows;++i) glp_set_row_bnds(mip,i+1,GLP_FX,rh,rh);
    if (isweighted) {
        int r=glp_add_rows(mip,1);
        double lb=lower_bound();
        glp_set_row_bnds(mip,r,GLP_LO,lb,DBL_MAX);
        for (i=0;i<ne;++i) {
            ia[++cnt]=r; ja[cnt]=i+1; ar[cnt]=obj[i];
        }
    }
    glp_add_cols(mip,ncols);
    for (j=0;j<ne;++j) {
        glp_set_col_kind(mip,j+1,GLP_BV);
        glp_set_obj_coef(mip,j+1,obj[j]);
    }
    assert(cnt<=nonzeros);
    glp_load_matrix(mip,cnt,ia,ja,ar);
    delete[] ia; delete[] ja; delete[] ar;
}

/* compute lower bound for the objective function by using a minimal spanning tree */
double graphe::tsp::lower_bound() {
    double cost,maxcost=0;
    for (int i=0;i<nv;++i) {
        vertex &v=G->node(i);
        G->unset_subgraphs(0);
        v.set_subgraph(1);
        graphe T(G->giac_context());
        G->minimal_spanning_tree(T,0);
        ipairs E;
        cost=0;
        T.get_edges_as_pairs(E);
        for (ipairs_iter it=E.begin();it!=E.end();++it) {
            cost+=to_real_number(T.weight(*it),G->giac_context()).to_double(G->giac_context());
        }
        dvector c;
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            c.push_back(weight(i,*it));
        }
        assert(c.size()>1);
        std::sort(c.begin(),c.end());
        cost+=c[0];
        cost+=c[1];
        if (cost>maxcost)
            maxcost=cost;
    }
    return maxcost;
}

/* find the subtours (sv) in the subgraph induced by vertices in v */
bool graphe::tsp::find_tours(int k,ivectors &sv,solution_status &status) {
    sv.clear();
    if (k<=0)
        return true;
    sv.reserve(k);
    glp_erase_prob(mip);
    glp_set_obj_dir(mip,GLP_MIN);
    formulate_mip();
    glp_smcp lparm;
    glp_init_smcp(&lparm);
    lparm.msg_lev=GLP_MSG_OFF;
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.msg_lev=GLP_MSG_ON;
    parm.gmi_cuts=GLP_ON;
    parm.cov_cuts=GLP_ON;
    parm.br_tech=GLP_BR_PCH;
    parm.bt_tech=GLP_BT_BPH;
    parm.sr_heur=GLP_OFF;
    parm.cb_func=&callback;
    parm.cb_info=static_cast<void*>(this);
    parm.mip_gap=gap_tol;
    bool retval=true;
    int res,lres,stat;
    is_symmetric_tsp=isweighted && G->is_clique();
    heur_type=is_symmetric_tsp?_GT_TSP_CHRISTOFIDES_SA:(isweighted?_GT_TSP_FARTHEST_INSERTION_HEUR:_GT_TSP_NO_HEUR);
    int term_old;
    if (verbose)
        glp_term_hook(term_hook,static_cast<void*>(logptr(G->giac_context())));
    else term_old=glp_term_out(GLP_OFF);
    do {
        do {
            /* append subtour elimination constraints */
            for (set<ivector>::const_iterator it=subtours.begin();it!=subtours.end();++it) {
                append_sce(*it);
            }
            subtours.clear();
            /* solve */
            if ((lres=glp_simplex(mip,&lparm))==0 &&
                ((stat=glp_get_status(mip))==GLP_OPT || stat==GLP_FEAS))
                res=glp_intopt(mip,&parm);
            else if (lres!=0 || stat==GLP_UNDEF || stat==GLP_UNBND) {
                status=_GT_TSP_ERROR;
                break;
            } else {
                status=_GT_TSP_NOT_HAMILTONIAN;
                break;
            }
            stat=glp_mip_status(mip);
            if ((res!=0 && res!=GLP_EMIPGAP && res!=GLP_ESTOP) || stat==GLP_UNDEF) {
                status=_GT_TSP_ERROR;
                retval=false;
                break;
            }
            switch (stat) {
            case GLP_FEAS:
                G->message("Warning: the solution is not necessarily optimal");
            case GLP_OPT:
                status=_GT_TSP_OPTIMAL;
                break;
            case GLP_NOFEAS:
                status=_GT_TSP_NOT_HAMILTONIAN;
                break;
            }
            if (status==_GT_TSP_NOT_HAMILTONIAN) break;
            if (!get_subtours()) {
                status=_GT_TSP_NOT_HAMILTONIAN;
                break;
            }
            if (subtours.size()>1 && verbose)
                G->message("Excluded %d subtours",subtours.size());
        } while (subtours.size()>1 && res!=GLP_ESTOP);
        if (status==_GT_TSP_ERROR)
            break;
        if (status==_GT_TSP_NOT_HAMILTONIAN) {
            if (!sv.empty())
                status=_GT_TSP_OPTIMAL;
            break;
        }
        if (subtours.size()>1) // the solver was terminated
            break;
        sv.push_back(*subtours.begin());
        if (verbose)
            G->message("Found tour %d",sv.size());
    } while (int(sv.size())<k);
    if (verbose)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
    return retval;
}

/* append subtour elimination constraint to MIP */
void graphe::tsp::append_sce(const ivector &subtour) {
    int i,j,r=glp_add_rows(mip,1),cnt=0,sz=subtour.size();
    if (3*sz>2*nv+1) {
        for (int k=0;k<ne;++k) {
            if (find(subtour.begin(),subtour.end(),k)!=subtour.end())
                continue;
            i=k;
            const arc &a1=arcs[i];
            for (ivector_iter jt=subtour.begin();jt!=subtour.end();++jt) {
                j=*jt;
                assert(i!=j);
                const arc &a2=arcs[j];
                if (a1.head==a2.head || a1.tail==a2.tail ||
                        a1.head==a2.tail || a1.tail==a2.head) {
                    indices[++cnt]=k+1;
                    coeff[cnt]=1.0;
                    break;
                }
            }
        }
        glp_set_row_bnds(mip,r,GLP_LO,2,DBL_MAX);
    } else {
        for (ivector_iter it=subtour.begin();it!=subtour.end();++it) {
            indices[++cnt]=*it+1;
            coeff[cnt]=1.0;
        }
        glp_set_row_bnds(mip,r,GLP_UP,0.0,sz-1);
    }
    glp_set_mat_row(mip,r,cnt,indices,coeff);
}

/* append subtours from sv */
void graphe::tsp::add_subtours(const ivectors &sv) {
    ivector subtour;
    for (ivectors_iter it=sv.begin();it!=sv.end();++it) {
        subtours.insert(canonical_subtour(*it));
    }
}

/* collect subtours from the current best MIP integer solution */
bool graphe::tsp::get_subtours() {
    ivector sol;
    int m=glp_get_num_cols(mip),i,j;
    for (i=0;i<m;++i) {
        if (glp_mip_col_val(mip,i+1)>0.5)
            sol.push_back(i);
    }
    if (sol==old_sol) return false;
    old_sol=sol;
    bvector vst(sol.size(),false);
    ivector subtour;
    subtour.reserve(ne);
    arc a;
    for (ivector_iter st=sol.begin();st!=sol.end();++st) {
        const arc &fa=arcs[*st];
        i=st-sol.begin();
        if (vst[i])
            continue;
        subtour.clear();
        subtour.push_back(*st);
        vst[i]=true;
        a=fa;
        ivector_iter it;
        do {
            for (it=sol.begin();it!=sol.end();++it) {
                j=it-sol.begin();
                const arc &aa=arcs[*it];
                if (i==j || vst[j])
                    continue;
                if (aa.tail==a.head || aa.tail==a.tail || aa.head==a.head || aa.head==a.tail) {
                    vst[j]=true;
                    subtour.push_back(*it);
                    a.head=aa.head;
                    a.tail=aa.tail;
                    i=j;
                    break;
                }
            }
        } while (it!=sol.end());
        subtours.insert(canonical_subtour(subtour));
    }
    return true;
}

/* return canonical representation of the subtour (starting with the lowest edge index) */
graphe::ivector graphe::tsp::canonical_subtour(const ivector &subtour) {
    int minv=-1,pos;
    for (ivector_iter it=subtour.begin();it!=subtour.end();++it) {
        if (minv<0 || *it<minv) {
            minv=*it;
            pos=it-subtour.begin();
        }
    }
    ivector res(subtour.begin()+pos,subtour.end());
    res.insert(res.end(),subtour.begin(),subtour.begin()+pos);
    if (res.size()>2 && res[1]>res.back())
        std::reverse(res.begin()+1,res.end()); // canonical orientation
    return res;
}

/* find k shortest tours */
int graphe::tsp::solve(int k,ivectors &hcv,dvector &costs) {
    mip=glp_create_prob();
    hcv.clear(); hcv.reserve(k);
    costs.clear(); costs.reserve(k);
    ivectors sv;
    ivector hc; hc.reserve(nv);
    solution_status status;
    if (!find_tours(k,sv,status))
        return -1;
    int retval;
    switch (status) {
    case _GT_TSP_OPTIMAL:
        for (ivectors_iter t=sv.begin();t!=sv.end();++t) {
            assert(int(t->size())==nv);
            hc.clear();
            const arc &first=arcs[t->front()],&last=arcs[t->back()];
            int i=first.tail,j=first.head,k;
            if (i==last.tail || i==last.head) {
                hc.push_back(i);
                hc.push_back(k=j);
            } else {
                hc.push_back(j);
                hc.push_back(k=i);
            }
            for (ivector_iter it=t->begin()+1;it!=t->end();++it) {
                const arc &a=arcs[*it];
                if (a.tail==k) hc.push_back(k=a.head);
                else {
                    assert(a.head==k);
                    hc.push_back(k=a.tail);
                }
            }
            if (gap_tol>0)
                improve_tour(hc);
            hcv.push_back(hc);
            costs.push_back(tour_cost(hc));
        }
        retval=1; // success
        break;
    case _GT_TSP_NOT_HAMILTONIAN:
        retval=0; // fail
        break;
    case _GT_TSP_ERROR:
        retval=-1; // error
        break;
    }
    glp_delete_prob(mip);
    return retval;
}

/* find a shortest tour */
int graphe::tsp::solve(ivector &hc,double &cost) {
    ivectors hcv;
    dvector costs;
    int res=solve(1,hcv,costs);
    if (res==1 && !hcv.empty() && !costs.empty()) {
        hc=hcv.front();
        cost=costs.front();
    }
    return res;
}

/* return the weight of the edge (i,j) */
double graphe::tsp::weight(int i,int j) {
    ipair e=make_edge(i,j);
    return weight_map[e.first][e.second];
}

/* return the cost of the tour hc */
double graphe::tsp::tour_cost(const ivector &hc) {
    int n=hc.size()-1,v,w;
    if (!isweighted)
        return n;
    double cost=0;
    for (int i=0;i<n;++i) {
        v=hc[i];
        w=hc[i+1];
        assert(G->has_edge(v,w));
        cost+=weight(v,w);
    }
    return cost;
}

/* apply heuristics */
void graphe::tsp::heur(glp_tree *tree) {
    if (heur_type==_GT_TSP_NO_HEUR)
        return;
    int i,j,k;
    if (heur_type==_GT_TSP_CHRISTOFIDES_SA) { // symmetric TSP
        heur_type=_GT_TSP_FARTHEST_INSERTION_RANDOM;
        if (!christofides(tour)) {
            heur(tree);
            return;
        }
    } else { // weighted undirected TSP
        /* choose the initial arc a by random such that weight(e) >= median weight,
         * in the first pass try to construct a tour starting from heaviest edge */
        int index=heur_type==_GT_TSP_FARTHEST_INSERTION_RANDOM?(ne+1)/2+G->rand_integer(ne/2):ne-1;
        farthest_insertion(index,tour);
        heur_type=_GT_TSP_FARTHEST_INSERTION_RANDOM;
        if (int(tour.size())<=nv)
            return;
    }
    assert(int(tour.size())==nv+1);
    /* optimize the tour */
    lin_kernighan(tour);
    /* construct the heuristic solution and pass it to the MIP solver */
    for (i=0;i<ne;++i) coeff[i+1]=0.0;
    for (i=0;i<nv;++i) {
        j=i+1;
        k=edge_index(make_edge(tour[i],tour[j]));
        coeff[k+1]=1.0;
    }
    glp_ios_heur_sol(tree,coeff);
}

/* farthest insertion heuristic for obtaining a Hamiltonian tour */
void graphe::tsp::farthest_insertion(int index,ivector &hc) {
    int i,j,k,f;
    const arc &a=arcs[index];
    i=a.tail;
    j=a.head;
    hc.clear();
    hc.push_back(i);
    hc.push_back(j);
    hc.push_back(i);
    ivector_iter jt;
    for (int vis=0;vis<nv;++vis) visited[vis]=false;
    visited[i]=visited[j]=true;
    dvector W;
    double w,maxw,minw;
    /* proceed with farthest insertion heuristic */
    do {
        maxw=0;
        f=-1;
        for (k=0;k<nv;++k) {
            if (visited[k])
                continue;
            W.clear();
            for (ivector_iter it=hc.begin();it!=hc.end();++it) {
                if (G->has_edge(k,*it))
                    W.push_back(weight(k,*it));
            }
            minw=DBL_MAX;
            if (W.size()<2)
                continue;
            else for (dvector_iter wt=W.begin();wt!=W.end();++wt) {
                if (*wt<minw)
                    minw=*wt;
            }
            if (minw>maxw) {
                maxw=minw;
                f=k;
            }
        }
        if ((k=f)<0)
            break;
        i=-1;
        minw=DBL_MAX;
        for (ivector_iter it=hc.begin();it!=hc.end();++it) {
            jt=it+1;
            if (jt!=hc.end() && G->has_edge(*it,k) && G->has_edge(k,*jt) &&
                    minw>(w=weight(*it,k)+weight(k,*jt)-weight(*it,*jt))) {
                minw=w;
                i=it-hc.begin();
            }
        }
        if (i<0)
            break;
        hc.insert(hc.begin()+i+1,k);
        visited[k]=true;
    } while (int(hc.size())<=nv);
}

/* try to lower the cost of the tour hc */
bool graphe::tsp::perform_3opt_moves(ivector &hc) {
    int n=hc.size()-1,b1,e1,b2,e2,b3,e3,i1,j1,i2,j2,i3,j3,i,j,k,var,iter_count=0,moves_count=0;
    double opt_timeout=5.0+25.0*std::exp(-std::pow(std::max(0,1000-n),2)/2e5);
    bvector visited(nv);
    ivectors opt_moves;
    double sw,save,maxsave,wb1e1,wb2e2,wb3e3,wb1b2,wb1e2,we1e2;
    bool hb1e2,hb1b2,he1e2;
    /* improve the tour by making a limited number of 3-opt moves */
    clock_t start_time=clock();
    ivector ijk(3);
    bool triplet_found=false;
    bool timed_out=false;
    do {
        if (cancellable && is_interrupted())
            return false;
        ++iter_count;
        opt_moves.clear();
        std::fill(visited.begin(),visited.end(),false);
        for (i=0;i<n;++i) {
            triplet_found=false;
            b1=hc[i]; e1=hc[i+1];
            if (visited[b1] || visited[e1])
                continue;
            wb1e1=weight(b1,e1);
            for (j=i+2;j<n;++j) {
                if (double(clock()-start_time)/CLOCKS_PER_SEC>opt_timeout) {
                    timed_out=true;
                    break;
                }
                b2=hc[j]; e2=hc[j+1];
                if (visited[b2] || visited[e2])
                    continue;
                wb2e2=weight(b2,e2);
                hb1e2=G->has_edge(b1,e2);
                hb1b2=G->has_edge(b1,b2);
                he1e2=G->has_edge(e1,e2);
                if (hb1b2) wb1b2=weight(b1,b2);
                if (hb1e2) wb1e2=weight(b1,e2);
                if (he1e2) we1e2=weight(e1,e2);
                for (k=j+2;k<n;++k) {
                    b3=hc[k]; e3=hc[k+1];
                    if (visited[b3] || visited[e3])
                        continue;
                    wb3e3=weight(b3,e3);
                    sw=wb1e1+wb2e2+wb3e3;
                    if ((hb1e2 && G->has_edge(b3,e1) && G->has_edge(b2,e3) && sw>wb1e2+weight(b3,e1)+weight(b2,e3)) ||
                            (hb1e2 && G->has_edge(b3,b2) && G->has_edge(e1,e3) && sw>wb1e2+weight(b3,b2)+weight(e1,e3)) ||
                            (hb1b2 && G->has_edge(b3,e1) && G->has_edge(e2,e3) && sw>wb1b2+weight(b3,e1)+weight(e2,e3)) ||
                            (G->has_edge(b1,b3) && he1e2 && G->has_edge(b2,e3) && sw>weight(b1,b3)+we1e2+weight(b2,e3)))
                    {
                        ivector opt_move(6);
                        opt_move[0]=b1; opt_move[1]=b2; opt_move[2]=b3;
                        opt_move[3]=e1; opt_move[4]=e2; opt_move[5]=e3;
                        opt_moves.push_back(opt_move);
                        visited[b1]=visited[b2]=visited[b3]=true;
                        visited[e1]=visited[e2]=visited[e3]=true;
                        triplet_found=true;
                        break;
                    }
                }
                if (triplet_found)
                    break;
            }
            if (timed_out)
                break;
        }
        if (opt_moves.empty()) break;
        for (ivectors_iter it=opt_moves.begin();it!=opt_moves.end();++it) {
            const ivector &opt_move=*it;
            i1=find(hc.begin(),hc.end(),opt_move[0])-hc.begin();
            i2=find(hc.begin(),hc.end(),opt_move[1])-hc.begin();
            i3=find(hc.begin(),hc.end(),opt_move[2])-hc.begin();
            j1=find(hc.begin(),hc.end(),opt_move[3])-hc.begin();
            j2=find(hc.begin(),hc.end(),opt_move[4])-hc.begin();
            j3=find(hc.begin(),hc.end(),opt_move[5])-hc.begin();
            if (std::abs(i1-j1)!=1) {
                if (i1==0) i1=n;
                else if (j1==0) j1=n;
                else assert(false);
            }
            if (std::abs(i2-j2)!=1) {
                if (i2==0) i2=n;
                else if (j2==0) j2=n;
                else assert(false);
            }
            if (std::abs(i3-j3)!=1) {
                if (i3==0) i3=n;
                else if (j3==0) j3=n;
                else assert(false);
            }
            assert(std::abs(i1-j1)==1 && std::abs(i2-j2)==1 && std::abs(i3-j3)==1);
            ijk[0]=std::min(i1,j1);
            ijk[1]=std::min(i2,j2);
            ijk[2]=std::min(i3,j3);
            std::sort(ijk.begin(),ijk.end());
            i=ijk[0]; j=ijk[1]; k=ijk[2];
            b1=hc[i]; b2=hc[j]; b3=hc[k];
            e1=hc[i+1]; e2=hc[j+1]; e3=hc[k+1];
            sw=weight(b1,e1)+weight(b2,e2)+weight(b3,e3);
            maxsave=0.0;
            if (G->has_edge(b1,e2) && G->has_edge(b3,e1) && G->has_edge(b2,e3) &&
                    (save=sw-weight(b1,e2)-weight(b3,e1)-weight(b2,e3))>maxsave) {
                var=0; maxsave=save;
            }
            if (G->has_edge(b1,e2) && G->has_edge(b3,b2) && G->has_edge(e1,e3) &&
                  (save=sw-weight(b1,e2)-weight(b3,b2)-weight(e1,e3))>maxsave) {
                var=1; maxsave=save;
            }
            if (G->has_edge(b1,b2) && G->has_edge(b3,e1) && G->has_edge(e2,e3) &&
                  (save=sw-weight(b1,b2)-weight(b3,e1)-weight(e2,e3))>maxsave) {
                var=2; maxsave=save;
            }
            if (G->has_edge(b1,b3) && G->has_edge(e1,e2) && G->has_edge(b2,e3) &&
                  (save=sw-weight(b1,b3)-weight(e1,e2)-weight(b2,e3))>maxsave) {
                var=3; maxsave=save;
            }
            if (maxsave>.001) {
                ivector part(hc.begin()+j+1,hc.begin()+k+1);
                switch (var) {
                case 1:
                    std::reverse(hc.begin()+i+1,hc.begin()+j+1);
                case 3:
                    if (var==3) std::reverse(part.begin(),part.end());
                case 0:
                    hc.erase(hc.begin()+j+1,hc.begin()+k+1);
                    hc.insert(hc.begin()+i+1,part.begin(),part.end());
                    break;
                case 2:
                    std::reverse(hc.begin()+i+1,hc.begin()+j+1);
                    std::reverse(hc.begin()+j+1,hc.begin()+k+1);
                    break;
                default:
                    assert(false);
                }
                ++moves_count;
            }
        }
    } while (!timed_out);
    return true;
}

void graphe::tsp::straighten(ivector &hc) {
    ivectors opt_moves;
    int i,j,k,l,b1,b2,e1,e2,n=hc.size()-1,iter_count=0,moves_count=0;
    double w0;
    bvector visited(n+1);
    ivector opt_move(4);
    while (true) {
        if (cancellable && is_interrupted())
            break;
        ++iter_count;
        opt_moves.clear();
        std::fill(visited.begin(),visited.end(),false);
        for (i=0;i<n;++i) {
            if (visited[i] && visited[i+1]) continue;
            b1=hc[i]; e1=hc[i+1];
            w0=weight(b1,e1);
            for (j=i+2;j<n;++j) {
                if ((i==0 && j==n-1) || (visited[j] && visited[j+1])) continue;
                b2=hc[j]; e2=hc[j+1];
                if (G->has_edge(b1,b2) && G->has_edge(e1,e2) &&
                        w0+weight(b2,e2)-weight(b1,b2)-weight(e1,e2)>.001) {
                    opt_move[0]=b1; opt_move[1]=b2; opt_move[2]=e1; opt_move[3]=e2;
                    opt_moves.push_back(opt_move);
                    visited[i]=visited[i+1]=visited[j]=visited[j+1]=true;
                    break;
                }
            }
        }
        if (opt_moves.empty()) break;
        for (ivectors_iter it=opt_moves.begin();it!=opt_moves.end();++it) {
            const ivector &opt_move=*it;
            b1=opt_move[0]; b2=opt_move[1]; e1=opt_move[2]; e2=opt_move[3];
            i=find(hc.begin(),hc.end(),b1)-hc.begin(); j=find(hc.begin(),hc.end(),e1)-hc.begin();
            k=find(hc.begin(),hc.end(),b2)-hc.begin(); l=find(hc.begin(),hc.end(),e2)-hc.begin();
            if (std::abs(i-j)!=1) {
                if (i==0) i=n;
                else if (j==0) j=n;
                else assert(false);
            }
            if (std::abs(k-l)!=1) {
                if (k==0) k=n;
                else if (l==0) l=n;
                else assert(false);
            }
            assert(std::abs(i-j)==1 && std::abs(k-l)==1);
            i=std::min(i,j);
            j=std::min(k,l);
            std::reverse(hc.begin()+(i<j?i:j)+1,hc.begin()+(i<j?j:i)+1);
            ++moves_count;
        }
    }
}

/* Helsgaun criterion for k-opt move feasibility */
bool graphe::tsp::is_move_feasible(int k,const ivector &t,const ipairs &x) {
    ivector incl(2*k+1,0);
    int i,j,pi,count;
    for (int c=0;c<k;++c) {
        i=2*c+1;
        j=(i+1)%(2*k);
        incl[i+1]=j+1;
        incl[j+1]=i+1;
    }
    ivector p,q;
    ipairs xs(x.begin()+1,x.begin()+k+1);
    std::sort(xs.begin(),xs.end());
    for (ipairs_iter it=xs.begin();it!=xs.end();++it) {
        p.push_back(find(t.begin()+1,t.begin()+2*k+1,it->first)-t.begin());
    }
    for (i=k;i-->0;) {
        p.insert(p.begin()+i+1,((pi=p[i])%2)==0?pi-1:pi+1);
    }
    q.resize(p.size()+1);
    for (i=p.size();i-->0;) {
        assert((pi=p[i])<=2*k);
        q[pi]=i+1;
    }
    for (i=2*k,count=1;(i=q[incl[p[i-1]]]^1)!=0;count++);
    return count==k;
}

/* Lin-Kernighan backtracking k-opt heuristic for tour improvement */
bool graphe::tsp::lin_kernighan(ivector &hc) {
    int n=hc.size()-1;
    int i,choice,j,j0,j_next,pos,ypos,t1,t2,t3,t4,t5,t2i,t2ip,dir;
    static ipairs x,y,y_tmp;
    static dvector wx,wy;
    static ivector t,t1_alt,y1_alt,y2_alt,cheaper_hc;
    x.resize(n+1); y.resize(n+1);
    wx.resize(n+1); wy.resize(n+1);
    t.resize(2*n+1);
    cheaper_hc.resize(n+1);
    t1_alt.clear();
    double w,wi;
    for (int k=0;k<n;++k) t1_alt.push_back(k); // initialize alternatives for t1
lk2:
    /* Let i:=1. Choose t[1] and let x1:=(t[1],t[2]=succ(t[1])) in hc.
     * Initialize alternatives for y1. */
    i=1;
    if (cancellable && is_interrupted()) return false; // cancelled
    if (t1_alt.empty()) return true; // no more alternatives for t[1]
    choice=t1_alt.size()==1?0:G->rand_integer(t1_alt.size());
    t1=t[1]=t1_alt[choice]; t1_alt.erase(t1_alt.begin()+choice);
    x[1]=make_pair(t1,t2=t[2]=(t1+1)%n);
    w=wx[1]=weight(hc[t1],hc[t2]);
    y1_alt.clear();
    for (int k=0;k<n;++k) {
        if (k!=t1 && k!=t2 && k!=(t1-1+n)%n && k!=(t2+1)%n &&
                G->has_edge(hc[t2],hc[k]) && w>weight(hc[t2],hc[k]))
            y1_alt.push_back(k);
    }
lk4:
    /* Choose y1:=(t[2],t[3]) not in hc such that G1>0. Else, go to lk2. */
    if (y1_alt.empty()) goto lk2;
    choice=y1_alt.size()==1?0:G->rand_integer(y1_alt.size());
    t3=t[3]=y1_alt[choice]; y1_alt.erase(y1_alt.begin()+choice);
    y[1]=make_pair(t2=t[2],t3);
    wy[1]=weight(hc[t2],hc[t3]);
lk5:
    /* Let i:=i+1 and let xi=(t[2i-1],t[2i]) in hc.
     * Update the tour if possible. If i=1 initialize alternatives for y2. */
    ++i;
    t2ip=t[2*i-1];
    x[i]=make_pair(t2ip,t2i=t[2*i]=(t2ip+1)%n);
    wx[i]=weight(hc[t2ip],hc[t2i]);
    w=0.0;
    for (int k=1;k<=i;++k) w+=wx[k]-(k==i?0.0:wy[k]);
    if (t2i!=(t1-1+n)%n && G->has_edge(hc[t2i],hc[t1]) &&
            w-weight(hc[t2i],hc[t1])>.001 && is_move_feasible(i,t,x)) {
        /* update hc with the new cheaper tour starting from t[2] */
        y[i]=make_pair(t2i,t1); // temporary yi
        y_tmp=ipairs(y.begin()+1,y.begin()+i+1);
        j=j0=t[2]; j_next=t[3];
        pos=ypos=dir=0;
        while (true) {
            if (ypos>=0)
                y_tmp.erase(y_tmp.begin()+ypos);
            cheaper_hc[pos++]=hc[j];
            j=j_next;
            if (j==j0) break;
            ypos=-1;
            for (int k=y_tmp.size();k-->0;) {
                const ipair &yk=y_tmp[k];
                if (yk.first==j || yk.second==j) {
                    j_next=yk.first==j?yk.second:yk.first;
                    ypos=k;
                    dir=0;
                    break;
                }
            }
            if (ypos<0) {
                if (dir==0)
                    dir=find(x.begin()+1,x.begin()+i+1,make_pair(j,(j+1)%n))==x.begin()+i+1?1:-1;
                j_next=(j+dir+n)%n;
            }
        }
        assert(pos==n);
        cheaper_hc.back()=hc[j0];
        hc=cheaper_hc;
        t1_alt.clear();
        for (int k=0;k<n;++k) t1_alt.push_back(k); // initialize alternatives for t1
        goto lk2;
    } else if (i==2) { // initialize alternatives for y2
        y2_alt.clear();
        t1=t[1]; t4=t[4];
        for (int k=0;k<n;++k) {
            if (find(t.begin()+1,t.begin()+5,k)==t.begin()+5 &&
                    find(t.begin()+1,t.begin()+5,(k+1)%n)==t.begin()+5 &&
                    G->has_edge(hc[t4],hc[k]) && w-weight(hc[t4],hc[k])>.001)
                y2_alt.push_back(k);
        }
    }
lk7:
    /* Choose yi:=(t[2i],t[2i+1]) not in hc such that Gi>0. If successful, go to lk5. */
    if (i==2) {
        if (y2_alt.empty()) { i=1; goto lk4; }
        choice=y2_alt.size()==1?0:G->rand_integer(y2_alt.size());
        y[2]=make_pair(t4=t[4],t5=t[5]=y2_alt[choice]); y2_alt.erase(y2_alt.begin()+choice);
        wy[2]=weight(hc[t4],hc[t5]);
        goto lk5;
    } else {
        int k;
        t2i=t[2*i]; t1=t[1];
        for (k=0;k<n;++k) {
            if (find(t.begin()+1,t.begin()+2*i+1,k)==t.begin()+2*i+1 &&
                    find(t.begin()+1,t.begin()+2*i+1,(k+1)%n)==t.begin()+2*i+1 &&
                    G->has_edge(hc[t2i],hc[k]) && w-(wi=weight(hc[t2i],hc[k]))>.001) {
                y[i]=make_pair(t2i,t[2*i+1]=k);
                wy[i]=wi;
                break;
            }
        }
        if (k<n) goto lk5;
    }
    i=2;
    goto lk7;
}

/* find minimum weight perfect matching in an undirected bipartite graph with
 * even number of vertices, edge set E and edge weights W */
bool graphe::tsp::min_weight_matching_bipartite(const ivector &eind,const dvector &weights,ivector &matched_arcs,bool msg) {
    set<int> Vset; // set of vertices
    for (ivector_iter it=eind.begin();it!=eind.end();++it) {
        const arc &a=arcs[*it];
        Vset.insert(a.head);
        Vset.insert(a.tail);
    }
    ivector V;
    for (set<int>::const_iterator it=Vset.begin();it!=Vset.end();++it) {
        V.push_back(*it);
    }
    int m=eind.size(),n=V.size(),nz=2*m,cnt=0,v;
    assert((n%2)==0);
    glp_prob *wp=glp_create_prob();
    glp_add_rows(wp,n);
    glp_add_cols(wp,m);
    glp_set_obj_dir(wp,GLP_MIN);
    int *ia=new int[nz+1];
    int *ja=new int[nz+1];
    double *ar=new double[nz+1];
    for (int i=0;i<n;++i) {
        v=V[i];
        for (int j=0;j<m;++j) {
            const arc &a=arcs[eind[j]];
            if (a.head==v || a.tail==v) {
                ia[++cnt]=i+1; ja[cnt]=j+1;
                ar[cnt]=1.0;
            }
        }
    }
    assert(cnt==nz);
    for (int i=0;i<n;++i) glp_set_row_bnds(wp,i+1,GLP_FX,1.0,1.0);
    for (int j=0;j<m;++j) {
        glp_set_col_kind(wp,j+1,GLP_BV);
        glp_set_obj_coef(wp,j+1,weights[j]);
    }
    glp_load_matrix(wp,nz,ia,ja,ar);
    delete[] ia; delete[] ja;
    delete[] ar;
    glp_smcp lparm;
    glp_init_smcp(&lparm);
    lparm.msg_lev=GLP_MSG_OFF;
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.br_tech=GLP_BR_PCH;
    parm.bt_tech=GLP_BT_BPH;
    parm.gmi_cuts=GLP_ON;
    parm.mir_cuts=GLP_ON;
    parm.clq_cuts=GLP_ON;
    parm.cov_cuts=GLP_ON;
    parm.msg_lev=GLP_MSG_ON;
    parm.fp_heur=GLP_ON;
    parm.sr_heur=GLP_ON;
    parm.cb_func=&min_wpm_callback;
    pair<const ivector*,tsp*> info;
    info.first=&eind;
    info.second=this;
    parm.cb_info=static_cast<void*>(&info);
    if (msg) // TSP approximation
        parm.mip_gap=gap_tol;
    int term_old,res;
    if (msg)
        glp_term_hook(term_hook,static_cast<void*>(logptr(G->giac_context())));
    else term_old=glp_term_out(GLP_OFF);
    if (glp_simplex(wp,&lparm)!=0 || ((res=glp_intopt(wp,&parm))!=0 && res!=GLP_EMIPGAP && res!=GLP_ESTOP))
        return false;
    int stat=glp_mip_status(wp);
    bool ret=true;
    if (stat==GLP_OPT || stat==GLP_FEAS) {
        for (int j=0;j<m;++j) {
            if (glp_mip_col_val(wp,j+1)!=0)
                matched_arcs.push_back(j);
        }
    } else ret=false;
    glp_delete_prob(wp);
    if (msg)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
    return ret;
}

void graphe::tsp::min_wpm_heur(glp_tree *tree,const ivector &eind) {
    bvector vertex_matched(nv,false),arc_matched(eind.size(),false);
    iset V;
    int i,j,k,n=eind.size();
    glp_prob *lp=glp_ios_get_prob(tree);
    for (i=0;i<n;++i) {
        if (glp_ios_can_branch(tree,i+1)==0 && glp_get_col_prim(lp,i+1)>0.9) {
            const arc &a=arcs[eind[i]];
            if (vertex_matched[a.head] || vertex_matched[a.tail])
                continue;
            arc_matched[i]=true;
            vertex_matched[a.head]=vertex_matched[a.tail]=true;
        }
    }
    for (ivector_iter it=eind.begin();it!=eind.end();++it) {
        k=it-eind.begin();
        if (arc_matched[k])
            continue;
        const arc &a=arcs[*it];
        i=a.tail; j=a.head;
        V.insert(i); V.insert(j);
        if (vertex_matched[i] || vertex_matched[j])
            continue;
        arc_matched[k]=true;
        vertex_matched[i]=vertex_matched[j]=true;
    }
    int cnt=0;
    for (bvector::const_iterator it=arc_matched.begin();it!=arc_matched.end();++it) {
        if (*it) ++cnt;
        coeff[it-arc_matched.begin()+1]=*it?1.0:0.0;
    }
    if (int(V.size())==2*cnt)
        glp_ios_heur_sol(tree,coeff);
}

/* construct a tour using the method of Christofides for symmetric TSP */
bool graphe::tsp::christofides(ivector &hc,bool show_progress) {
    hc.clear();
    /* create spanning tree T */
    if (show_progress)
        G->message("Finding a minimal spanning tree...");
    graphe T(G->giac_context());
    G->minimal_spanning_tree(T);
    ivector V;
    for (int i=0;i<nv;++i) {
        if ((T.degree(i)%2)!=0)
            V.push_back(G->node_index(T.node_label(i)));
    }
    /* the number of odd degree vertices in T must be even */
    if (V.size()<=0 || (V.size()%2)!=0)
        return false;
    /* get edges in SG\T connecting the vertices from V,
     * these form a bipartite graph B(V,E) */
    ivector eind,matched_arcs;
    dvector weights;
    for (int k=0;k<ne;++k) {
        const arc &a=arcs[k];
        if (find(V.begin(),V.end(),a.head)!=V.end() &&
                find(V.begin(),V.end(),a.tail)!=V.end() &&
                !T.has_edge(T.node_index(G->node_label(a.tail)),
                            T.node_index(G->node_label(a.head)))) {
            eind.push_back(k);
            weights.push_back(weight(a.tail,a.head));
        }
    }
    if (show_progress)
        G->message("Finding a minimum weight matching in a bipartite graph with %d edges...",eind.size());
    bool has_mwm=true;
    if (eind.size()==1)
        matched_arcs.resize(1,0);
    else has_mwm=min_weight_matching_bipartite(eind,weights,matched_arcs,show_progress);
    if (!has_mwm || 2*matched_arcs.size()!=V.size())
        return false;
    /* add matched edges to T and find Eulerian circuit,
     * such exists because T now has all vertex degrees even */
    for (ivector_iter it=matched_arcs.begin();it!=matched_arcs.end();++it) {
        const arc &a=arcs[eind[*it]];
        const attrib &attr=G->edge_attributes(a.tail,a.head);
        T.add_edge(G->node_label(a.tail),G->node_label(a.head),attr);
    }
    ivector etrail; // eulerian trail
    if (!T.find_eulerian_trail(etrail) || etrail.front()!=etrail.back())
        return false;
    bvector visited(nv,false);
    for (ivector_iter it=etrail.begin();it!=etrail.end();++it) {
        if (visited[*it])
            continue;
        hc.push_back(G->node_index(T.node_label(*it)));
        visited[*it]=true;
    }
    hc.push_back(hc.front());
    return true;
}

/* improve the tour hc by making k-opt moves */
void graphe::tsp::improve_tour(ivector &hc,bool do_3opt) {
    cancellable=true;
    if (lin_kernighan(hc))
        if (!do_3opt || perform_3opt_moves(hc))
            straighten(hc);
    cancellable=false;
}

/* construct an approximate tour, G must be a complete graph */
bool graphe::tsp::approx(ivector &hc,double &ratio) {
    assert(isweighted);
    if (!christofides(hc,verbose)) {
        if (verbose)
            G->message("Error: failed to construct Christofides tour");
        return false;
    }
    double old_cost=tour_cost(hc),imp_cost;
    if (verbose)
        *logptr(G->giac_context()) << "Constructed a tour with weight " << old_cost << "\n";
    if (old_cost==0)
        ratio=0;
    else {
        if (verbose)
            G->message("Attempting to improve the tour...");
        improve_tour(hc);
        imp_cost=tour_cost(hc);
        if (verbose)
            *logptr(G->giac_context()) << "Relative improvement: " << 100.0*(old_cost-imp_cost)/old_cost << "%\n";
        ratio=(1+gap_tol)*1.5*imp_cost/old_cost;
    }
    return true;
}

/* compute the mean and the standard deviation of the given sample */
void graphe::tsp::sample_mean_stddev(const dvector &sample,double &mean,double &stddev) {
    assert(!sample.empty());
    /* compute the mean */
    double sum=0;
    for (dvector_iter it=sample.begin();it!=sample.end();++it) {
        sum+=*it;
    }
    mean=sum/double(sample.size());
    if (sample.size()==1)
        stddev=0;
    else {
        double d;
        sum=0;
        for (dvector_iter it=sample.begin();it!=sample.end();++it) {
            d=*it-mean;
            sum+=d*d;
        }
        stddev=std::sqrt(sum/double(sample.size()-1));
    }
}

/* branching variable selection technique proposed by Padberg & Rinaldi */
void graphe::tsp::select_branching_variable(glp_tree *tree) {
    glp_prob *lp=glp_ios_get_prob(tree);
    int n=glp_get_num_cols(lp),b=-1;
    double x0=0.0,x1=1.0,xe;
    double par=0.5,rho=3.5,bctol=0.05; // constants
    for (int i=0;i<n;++i) {
        can_branch[i]=glp_ios_can_branch(tree,i+1)!=0;
        xe=xev[i]=glp_get_col_prim(lp,i+1);
        if (xe-0.5>=0 && xe<x1) x1=xe;
        if (xe-0.5<=0 && xe>x0) x0=xe;
        obj[i]=glp_get_obj_coef(lp,i+1);
    }
    ivector I;
    for (int i=0;i<n;++i) {
        xe=xev[i];
        if (can_branch[i] && (1-par)*x0<=xe && xe<=x1+par*(1-x1))
            I.push_back(i);
    }
    assert(!I.empty());
    if (I.size()==1)
        b=I.front();
    else {
        int N=I.size();
        vector<pair<double,int> > sorted_by_obj_coeff;
        ivector ind(N);
        for (ivector_iter it=I.begin();it!=I.end();++it) {
            sorted_by_obj_coeff.push_back(make_pair(obj[*it],*it));
        }
        std::sort(sorted_by_obj_coeff.begin(),sorted_by_obj_coeff.end());
        dvector data2(N),data1(0);
        for (int i=0;i<N;++i) {
            ind[i]=sorted_by_obj_coeff[i].second;
            data2[i]=sorted_by_obj_coeff[i].first;
        }
        std::reverse(data2.begin(),data2.end());
        double mu1,mu2,sigma1,sigma2;
        int i=1;
        for (;i<N;++i) {
            data1.push_back(data2.back());
            data2.pop_back();
            sample_mean_stddev(data1,mu1,sigma1);
            sample_mean_stddev(data2,mu2,sigma2);
            if (mu1+rho*sigma1<mu2-rho*sigma2)
                break;
        }
        if (i==N) i=0;
        /* the set I2 starts from i in ind */
        ivector I2(ind.begin()+i,ind.end());
        assert(!I2.empty());
        N=I2.size();
        double dist,mindist=DBL_MAX;
        int j0=-1;
        for (int j=0;j<N;++j) {
            if ((dist=std::abs(0.5-xev[I2[j]]))<mindist) {
                mindist=dist;
                j0=j;
            }
        }
        assert(j0>=0);
        x0=xev[I2[j0]];
        double c,maxc=0.0;
        for (ivector_iter it=I2.begin();it!=I2.end();++it) {
            xe=xev[*it];
            if (xe>=x0-bctol && xe<=x0+bctol && (c=obj[*it])>maxc) {
                maxc=c;
                b=*it;
            }
        }
    }
    assert(b>=0);
    /* Padberg & Rinaldi suggest exploring the upper branch first */
    glp_ios_branch_upon(tree,b+1,GLP_UP_BRNCH);
}

/* row generation (originally written by Andrew Makhorin) */
void graphe::tsp::rowgen(glp_tree *tree) {
    if (mincut_data.empty())
        mincut_data.resize(22);
    glp_prob *prob=glp_ios_get_prob(tree);
    int i,j,k,m=ne,n=nv,nedg=0,nz;
    double sum,xk;
    for (k=0;k<m;k++) {
        if (glp_get_col_prim(prob,k+1)>=.001)
            nedg++;
    }
    /* build the capacitated network */
    assert(nedg<=m);
    ivector &beg=mincut_data[14]; if (int(beg.size())<nedg) beg.resize(nedg);
    ivector &end=mincut_data[15]; if (int(end.size())<nedg) end.resize(nedg);
    ivector &cap=mincut_data[16]; if (int(cap.size())<nedg) cap.resize(nedg);
    nz=0;
    for (k=0;k<m;k++) {
        if ((xk=glp_get_col_prim(prob,k+1))>=.001) {
            const arc &a=arcs[k];
            beg[nz]=a.tail+1; end[nz]=a.head+1;
            rlx_sol_map[a.tail][a.head]=xk;
            cap[nz++]=std::ceil(1000*xk); // scale edge capacities to make them integral
        }
    }
    assert(nz==nedg);
    ivector &cut=mincut_data[17]; if (int(cut.size())<n) cut.resize(n);
    minimal_cut(n,nedg,beg,end,cap,cut); // find minimal cut in the capacitated network
    nedg=0; sum=0;
    for (i=0;i<n;i++) {
        for (j=i+1;j<n;j++) {
            if ((cut[i] && !cut[j]) || (!cut[i] && cut[j])) {
                nedg++;
                sum+=rlx_sol_map[i][j];
            }
        }
    }
    /* if the (unscaled) capacity of min cut is less than 2, the
     * corresponding subtour elimination constraint is violated */
    if (sum<=1.999) {
        nz=0;
        for (i=0;i<n;i++) {
            for (j=i+1;j<n;j++) {
                if ((cut[i] && !cut[j]) || (!cut[i] && cut[j])) {
                    nz++;
                    assert(nz<=nedg);
                    indices[nz]=edge_index(make_pair(i,j))+1;
                    coeff[nz]=1.0;
                }
            }
        }
        assert(nz==nedg);
        i=glp_add_rows(prob,1);
        glp_set_row_bnds(prob,i,GLP_LO,2.0,0.0);
        glp_set_mat_row(prob,i,nz,indices,coeff);
    }
}

/* MIP callback routine */
void graphe::tsp::callback(glp_tree *tree,void *info) {
    if (is_interrupted()) {
        glp_ios_terminate(tree);
        return;
    }
    tsp *tsprob=static_cast<tsp*>(info);
    switch (glp_ios_reason(tree)) {
    case GLP_IHEUR:
        tsprob->heur(tree);
        break;
    case GLP_IBRANCH:
        if (tsprob->isweighted && _rand(makesequence(0,1),tsprob->G->giac_context()).DOUBLE_val()<0.14)
            tsprob->select_branching_variable(tree);
        break;
    case GLP_IROWGEN:
        if (tsprob->is_symmetric_tsp)
            tsprob->rowgen(tree);
        break;
    default:
        break;
    }
}

/* min weight perfect matching solver callback routine */
void graphe::tsp::min_wpm_callback(glp_tree *tree,void *info) {
    if (is_interrupted()) {
        glp_ios_terminate(tree);
        return;
    }
    pair<const ivector*,tsp*> *origin=static_cast<pair<const ivector*,tsp*>*>(info);
    const ivector *eind=origin->first;
    tsp *tsprob=origin->second;
    switch (glp_ios_reason(tree)) {
    case GLP_IHEUR:
        tsprob->min_wpm_heur(tree,*eind);
        break;
    default:
        break;
    }
}

/*
* END OF TSP CLASS IMPLEMENTATION
*/

/*
 * ATSP CLASS (asymmetric traveling salesman problem)
 */
 
graphe::atsp::atsp(graphe *gr,const ipairs &must_include_arcs,double gap_tolerance,bool is_verbose) {
    G=gr;
    mia=must_include_arcs;
    mip=NULL;
    isweighted=G->is_weighted();
    gap_tol=gap_tolerance;
    verbose=is_verbose;
}

graphe::atsp::~atsp() {
    if (mip!=NULL) glp_delete_prob(mip);
}

bool graphe::atsp::solve(ivector &hc,double &cost) {
    if (mip!=NULL) glp_erase_prob(mip);
    else mip=glp_create_prob();
    glp_set_obj_dir(mip,GLP_MIN);
    int n=G->node_count(),nz=8+8*n*(n-2)+n*ft.size(),i,j,r,c=0,res;
    int *u=new int[n];
    int **x=new int*[n];
    for (i=0;i<n;++i) x[i]=new int[n];
    int *ia=new int[nz+1];
    int *ja=new int[nz+1];
    double *ar=new double[nz+1];
    double w;
    gen e;
    ipairs E;
    G->get_edges_as_pairs(E);
    /* add variables */
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        i=it->first; j=it->second;
        x[i][j]=glp_add_cols(mip,1);
        if (find(mia.begin(),mia.end(),*it)!=mia.end())
            glp_set_col_bnds(mip,x[i][j],GLP_FX,1,1);
        else glp_set_col_kind(mip,x[i][j],GLP_BV);
        if (isweighted) {
            if (!is_real_number(G->weight(i,j),G->giac_context())) {
                G->message("Warning: arc weight is not numeric, replacing by 1");
                w=1;
            } else w=to_real_number(G->weight(i,j),G->giac_context()).to_double(G->giac_context());
        } else w=1;
        glp_set_obj_coef(mip,x[i][j],w);
    }
    for (i=1;i<n;++i) {
        u[i]=glp_add_cols(mip,1);
        glp_set_col_bnds(mip,u[i],GLP_DB,2,n);
    }
    /* create the basic constraints */
    for (i=0;i<n;++i) {
        r=glp_add_rows(mip,2);
        glp_set_row_bnds(mip,r,GLP_FX,1,1);
        glp_set_row_bnds(mip,r+1,GLP_FX,1,1);
        for (j=0;j<n;++j) {
            if (i==j) continue;
            if (G->has_edge(i,j)) { ia[++c]=x[i][j]; ja[c]=r; ar[c]=1; }
            if (G->has_edge(j,i)) { ia[++c]=x[j][i]; ja[c]=r+1; ar[c]=1; }
        }
    }
    /* SEC: an enhanced formulation of Miller et al. */
    for (i=1;i<n;++i) {
        for (j=1;j<n;++j) {
            if (i==j) continue;
            r=glp_add_rows(mip,1);
            glp_set_row_bnds(mip,r,GLP_UP,0,n-1);
            ia[++c]=u[i]; ja[c]=r; ar[c]=1;
            ia[++c]=u[j]; ja[c]=r; ar[c]=-1;
            if (G->has_edge(i,j)) { ia[++c]=x[i][j]; ja[c]=r; ar[c]=n; }
            if (G->has_edge(j,i)) { ia[++c]=x[j][i]; ja[c]=r; ar[c]=n-2; }
        }
    }
    /* forbid the specified tours */
    for (ivectors_iter it=ft.begin();it!=ft.end();++it) {
        r=glp_add_rows(mip,1);
        glp_set_row_bnds(mip,r,GLP_UP,0,n-1);
        for (int k=0;k<n;++k) {
            i=it->at(k); j=it->at(k+1);
            assert(G->has_edge(i,j));
            ia[++c]=x[i][j]; ja[c]=r; ar[c]=1;
        }
    }
    assert(c<=nz);
    glp_load_matrix(mip,c,ja,ia,ar);
    delete[] ia; delete[] ja; delete[] ar; delete[] u;
    /* solve the problem */
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.msg_lev=GLP_MSG_ON;
    parm.br_tech=GLP_BR_MFV;
    parm.bt_tech=GLP_BT_DFS;
    parm.presolve=GLP_ON;
    parm.gmi_cuts=GLP_ON;
    parm.mir_cuts=GLP_ON;
    parm.cov_cuts=GLP_ON;
    parm.clq_cuts=GLP_ON;
    parm.fp_heur=GLP_ON;
    parm.sr_heur=GLP_ON;
    parm.mip_gap=gap_tol;
    parm.cb_func=&callback;
    parm.cb_info=static_cast<void*>(this);
    select_blb=false;
    bool ret=false;
    if ((res=glp_intopt(mip,&parm))==0 || res==GLP_EMIPGAP || res==GLP_ESTOP) {
        ret=true;
        switch(glp_mip_status(mip)) {
        case GLP_FEAS:
            if (verbose)
                G->message("Warning: the solution is not necessarily optimal");
        case GLP_OPT:
            break;
        case GLP_UNDEF:
            if (verbose)
                G->message("Warning: the solution is undefined");
        case GLP_NOFEAS:
            ret=false;
            break;
        }
        if (ret) {
            cost=glp_mip_obj_val(mip);
            int n=G->node_count(),i,j,pos=0;
            hc.resize(n+1,0);
            while (pos<n) {
                i=hc[pos];
                for (j=0;j<n;++j) {
                    if (i==j || !G->has_edge(i,j)) continue;
                    if (glp_mip_col_val(mip,x[i][j])>0.5) {
                        hc[++pos]=j;
                        break;
                    }
                }
                assert(j<n);
            }
            assert(hc.back()==hc.front());
        }
    }
    terminated=(res==GLP_ESTOP);
    for (i=0;i<n;++i) delete[] x[i];
    delete[] x;
    return ret;
}

/* MIP callback routine */
void graphe::atsp::callback(glp_tree *tree,void *info) {
    if (is_interrupted()) {
        glp_ios_terminate(tree);
        return;
    }
    atsp *atsprob=static_cast<atsp*>(info);
    switch (glp_ios_reason(tree)) {
    case GLP_ISELECT:
        if (atsprob->select_blb)
            glp_ios_select_node(tree,glp_ios_best_node(tree));
        break;
    case GLP_IBINGO:
        atsprob->select_blb=true;
        break;
    default:
        break;
    }
}

void graphe::atsp::ksolve(int k,ivectors &hcv,dvector &costs) {
    hcv.clear();
    costs.clear();
    ft.clear();
    hcv.reserve(k);
    costs.reserve(k);
    ft.reserve(k-1);
    ivector hc;
    double cost;
    int term_old;
    if (verbose)
        glp_term_hook(term_hook,static_cast<void*>(logptr(G->giac_context())));
    else term_old=glp_term_out(GLP_OFF);
    terminated=false;
    for (int cnt=0;!terminated && cnt<k;++cnt) {
        if (!solve(hc,cost)) break;
        if (verbose)
            G->message("Found the tour %d",cnt+1);
        hcv.push_back(hc);
        costs.push_back(cost);
        if (cnt<k-1) ft.push_back(hc);
    }
    if (verbose)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
}

/*
* END OF ATSP CLASS IMPLEMENTATION
*/

#endif

/* Try to find an optimal Hamiltonian cycle.
 * Return 0 if the graph is not Hamiltonian, else store the circuit in h
 * and return 1, if unable to solve return -1 */
int graphe::traveling_salesman(int k,ivectors &hcv,dvector &costs,double gap_tol,bool verbose) {
#ifdef HAVE_LIBGLPK
    tsp t(this,gap_tol,verbose);
    if (k==0) { // tour approximation
        double ratio;
        ivector h;
        if (t.approx(h,ratio)) {
            //message("The tour cost is within %d%% of the optimal value",std::floor((ratio-1.0)*100.0+.5));
            costs.push_back(t.tour_cost(h));
            hcv.push_back(h);
            return 1;
        }
        return -1;
    }
    assert(k>0);
    return t.solve(k,hcv,costs);
#else
    message("Error: GLPK library is required for solving traveling salesman problem");
    return -1;
#endif
}

/* find first k optimal tours in a directed graph */
bool graphe::find_directed_tours(int k,ivectors &hcv,dvector &costs,const ipairs &incl,double gap_tol,bool verbose) {
    if (!is_directed()) {
        message("Error: the graph must be directed");
        return false;
    }
#ifdef HAVE_LIBGLPK
    atsp t(this,incl,gap_tol,verbose);
    t.ksolve(k,hcv,costs);
    return true;
#else
    message("Error: GLPK library is required for solving traveling salesman problem");
    return false;
#endif
}

/* for each edge of the graph, assign the distance between endpoints as a weight */
bool graphe::make_euclidean_distances() {
    layout x;
    if (is_directed() || !has_stored_layout(x))
        return false;
    if (!x.empty()) {
        int n=node_count();
        set_weighted(true);
        point pq(x.front().size());
        for (int i=0;i<n;++i) {
            const vertex &v=node(i);
            const point &p=x[i];
            for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                const point &q=x[*it];
                set_edge_attribute(i,*it,_GT_ATTRIB_WEIGHT,gen(point_distance(p,q,pq)));
            }
        }
    }
    return true;
}

/* find the maximum flow using Edmonds-Karp algorithm (exact computation) */
gen graphe::maxflow_edmonds_karp(int s,int t,vector<map<int,gen> > &flow,const gen &limit) {
    assert(is_directed() && node_queue.empty());
    gen mf(0),df; // the value of maximum flow
    int n=node_count(),i,j;
    vector<map<int,gen> > cap(n);
    /* initialize arc capacities and flows */
    flow.resize(n);
    bool isweighted=is_weighted();
    node_iter nt=nodes.begin(),ntend=nodes.end();
    ivector_iter it,itend;
    for (i=0;nt!=ntend;++nt,++i) {
        map<int,gen> &c=cap[i];
        flow[i].clear();
        for (it=nt->neighbors().begin(),itend=nt->neighbors().end();it!=itend;++it) {
            j=*it;
            c[j]=isweighted?weight(i,j):gen(1); // all capacities are set to 1 in unweighted graphs
        }
    }
    ipairs pred(n);
    ipair nullpair=make_pair(-1,-1);
    while (true) {
        /* run BFS */
        node_queue.push(s);
        std::fill(pred.begin(),pred.end(),nullpair);
        while (!node_queue.empty()) {
            i=node_queue.front();
            node_queue.pop();
            const vertex &v=node(i);
            for (it=v.neighbors().begin(),itend=v.neighbors().end();it!=itend;++it) {
                j=*it;
                if (j!=s && pred[j]==nullpair && is_strictly_greater(cap[i][j],flow[i][j],ctx)) {
                    pred[j]=make_pair(i,j);
                    node_queue.push(j);
                }
            }
        }
        if (pred[t]!=nullpair) {
            /* found an augmenting path */
            df=plusinf();
            /* compute df = how much flow can we send ... */
            for (ipair e=pred[t];e.first>=0;e=pred[e.first]) {
                df=min(df,cap[e.first][e.second]-flow[e.first][e.second],ctx);
            }
            /* ... and increase edge flows by df */
            for (ipair e=pred[t];e.first>=0;e=pred[e.first]) {
                flow[e.first][e.second]+=df;
                flow[e.second][e.first]-=df;
            }
            mf+=df;
            if (is_greater(mf,limit,ctx))
                break;
        } else break;
    }
    return mf; // return the maximum flow
}

/* obtain a minimum cut from maximum flow */
void graphe::minimum_cut(int s,const vector<map<int,gen> > &flow,ipairs &cut) {
    /* create the residual network */
    graphe G(ctx,false);
    G.add_nodes(node_count());
    G.make_directed();
    cut.clear();
    int i,j;
    gen c,f;
    bool isweighted=is_weighted();
    map<int,gen>::const_iterator mit;
    node_iter it,itend;
    ivector_iter jt,jtend;
    for (i=0,it=nodes.begin(),itend=nodes.end();it!=itend;++it,++i) {
        const map<int,gen> &flowi=flow[i];
        for (jt=it->neighbors().begin(),jtend=it->neighbors().end();jt!=jtend;++jt) {
            j=*jt;
            f=max(0,(mit=flowi.find(j))!=flowi.end()?mit->second:gen(0),ctx);
            c=isweighted?weight(i,j):gen(1);
            if (!giac::is_zero(_ratnormal(c-f,ctx),ctx))
                G.add_edge(i,j);
        }
    }
    /* run DFS on the residual network, store all discovered nodes */
    ivector disc;
    G.dfs(s,true,true,&disc);
    unvisit_all_nodes();
    for (jt=disc.begin(),jtend=disc.end();jt!=jtend;++jt) {
        node(*jt).set_visited(true);
    }
    /* if, for an edge (u,v), u is discovered and v is not
     * or vice versa, (u,v) belongs to the cut */
    for (i=0,it=nodes.begin(),itend=nodes.end();it!=itend;++it,++i) {
        for (jt=it->neighbors().begin(),jtend=it->neighbors().end();jt!=jtend;++jt) {
            if (it->is_visited() && !node(*jt).is_visited())
                cut.push_back(make_pair(i,*jt));
        }
    }
}

gen graphe::make_colon_label(const ivector &v) {
    string str;
    for (ivector_iter it=v.begin();it!=v.end();++it) {
        str+=int2string(*it);
        if (it+1!=v.end())
            str+=":";
    }
    return graphe::str2gen(str,true);
}

gen graphe::colon_label(int i, int j) {
    ivector v(2);
    v[0]=i; v[1]=j;
    return make_colon_label(v);
}

gen graphe::colon_label(int i, int j,int k) {
    ivector v(3);
    v[0]=i; v[1]=j; v[2]=k;
    return make_colon_label(v);
}

/* convert this multigraph to simple graph by adding new colored vertices */
void graphe::simplify(graphe &G,bool color_temp_vertices) const {
    assert(!is_directed() && (supports_attributes() || !G.supports_attributes()));
    ipairs E;
    int lab,k,i,hcol=default_highlighted_vertex_color;
    get_edges_as_pairs(E);
    G.clear();
    if (G.supports_attributes()) {
        G.add_nodes(vertices());
        lab=G.largest_integer_label();
    } else G.add_nodes(node_count());
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        const ipair &e=*it;
        k=multiedges(e);
        if (k>0) {
            for (int c=0;c<=k;++c) {
                i=G.supports_attributes()?G.add_node(++lab):G.add_node();
                if (color_temp_vertices) {
                    if (G.supports_attributes())
                        G.set_node_attribute(i,_GT_ATTRIB_COLOR,hcol);
                    else G.set_node_color(i,hcol);
                }
                G.add_edge(e.first,i);
                G.add_edge(e.second,i);
            }
        } else G.add_edge(e);
    }
}

/* permute vertices according to SHARC ordering by Michael Monagan */
void graphe::sharc_order() {
    /* assuming that the graph is connected */
    assert(!is_null() && node_queue.empty());
    int n=node_count(),i,j;
    unset_subgraphs(0);
    node(0).set_subgraph(1); // the set S has subgraph=1
    node(0).set_visited(true);
    ivector order(1,0); // vertices in S
    order.reserve(n);
    set<int> adj;
    bool has_path;
    while (subgraph_size(1)<n) {
        /* collect vertices adjacent to S */
        adj.clear();
        for (ivector_iter it=order.begin();it!=order.end();++it) {
            const vertex &v=node(*it);
            for (ivector_iter jt=v.neighbors().begin();jt!=v.neighbors().end();++jt) {
                if (node(*jt).subgraph()==0)
                    adj.insert(*jt);
            }
        }
        /* add vertices from S to queue */
        for (set<int>::const_iterator it=adj.begin();it!=adj.end();++it) {
            node_queue.push(*it);
        }
        unvisit_all_nodes(0);
        unset_all_ancestors(0);
        has_path=false;
        while (!node_queue.empty()) {
            i=node_queue.front();
            node_queue.pop();
            vertex &v=node(i);
            if (!v.is_visited()) {
                v.set_visited(true);
                for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
                    vertex &w=node(j=*it);
                    if (w.is_visited())
                        continue;
                    if (adj.find(j)!=adj.end()) {
                        /* path found */
                        w.set_ancestor(i);
                        do {
                            order.push_back(j);
                            node(j).set_subgraph(1);
                        } while ((j=node(j).ancestor())>=0);
                        clear_node_queue();
                        has_path=true;
                        break;
                    }
                    if (w.subgraph()==0) {
                        w.set_ancestor(i);
                        node_queue.push(j);
                    }
                }
            }
        }
        if (!has_path) {
            i=*adj.begin();
            node(i).set_subgraph(1);
            order.push_back(i);
        }
    }
    graphe G(*this);
    G.isomorphic_copy(*this,order);
}

map<graphe::ivector,vector<graphe::cpol> > graphe::cache;

/* implementation of cpol class */

void graphe::cpol::assign(const cpol &other) {
    nv=other.nv;
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
    cg=other.cg;
    col=other.col;
#else
    adj=other.adj;
#endif
    sz=other.sz;
    frq=other.frq;
    poly=other.poly;
}

graphe::cpol::cpol(const cpol &other) {
    assign(other);
}

#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
graphe::cpol::cpol(int n,ulong *g,int *c,size_t s,const intpoly &p) {
    nv=n;
    cg=g;
    col=c;
    sz=s;
    frq=0;
    poly=p;
}
#else
graphe::cpol::cpol(int n,int *a,size_t s,const intpoly &p) {
    nv=n;
    adj=a;
    sz=s;
    frq=0;
    poly=p;
}
#endif

graphe::cpol &graphe::cpol::operator =(const cpol &other) {
    assign(other);
    return *this;
}

#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
bool graphe::cpol::match(int n,ulong *g,int *c) const {
    if (nv!=n)
        return false;
    for (int i=0;i<n;++i) {
        if (c[i]!=col[i])
            return false;
    }
    size_t cnt=0;
    for (;cnt<sz;++cnt) {
        if (g[cnt]!=cg[cnt])
            return false;
    }
    return true;
}
#else
bool graphe::cpol::match(int n,int *a,int adj_sz) const {
    if (nv!=n || sz!=adj_sz)
        return false;
    for (int i=0;i<sz;++i) {
        if (adj[i]!=a[i])
            return false;
    }
    return true;
}
#endif

/* end of cpol class implementation */

/* add the polynomial b to a (a is changed in place) */
void graphe::poly_add(intpoly &a,const intpoly &b) {
    for (intpoly_iter it=b.begin();it!=b.end();++it) {
        a[it->first]+=it->second;
    }
}

/* multiply the polynomial a with b (a is changed in place) */
void graphe::poly_mult(intpoly &a,const intpoly &b) {
    intpoly c(a);
    a.clear();
    for (intpoly_iter it=b.begin();it!=b.end();++it) {
        const ipair &p=it->first;
        for (intpoly_iter jt=c.begin();jt!=c.end();++jt) {
            const ipair &q=jt->first;
            a[make_pair(p.first+q.first,p.second+q.second)]+=it->second*jt->second;
        }
    }
}

/* return the geometric sequence in variable var (1 or 2) of length k */
graphe::intpoly graphe::poly_geom(int var,int k,bool leading_one,bool add_other_var) {
    assert((var==1 || var==2) && k>0);
    intpoly res;
    for (int i=leading_one?0:1;i<k;++i) {
        res[make_pair(var==1?i:0,var==2?i:0)]=1;
    }
    if (add_other_var)
        res[make_pair(var==2?1:0,var==1?1:0)]=1;
    return res;
}

/* return 1 as polynom */
graphe::intpoly graphe::poly_one() {
    intpoly res;
    res[make_pair(0,0)]=1;
    return res;
}

/* convert list of integers to a 2-variable polynomial with integral coefficients */
gen graphe::intpoly2gen(const intpoly &v,const gen &x,const gen &y) {
    int pwx,pwy,cf;
    gen p(0);
    for (intpoly_iter it=v.begin();it!=v.end();++it) {
        pwx=it->first.first;
        pwy=it->first.second;
        cf=it->second;
        p+=gen(cf)*(pwx==0?gen(1):pow(x,pwx))*(pwy==0?gen(1):pow(y,pwy));
    }
    return p;
}

static int tutte_iter_count;
static int tutte_hits;
static clock_t tutte_time_start;
static double tutte_matching_time;

/* compute the Tutte polynomial for this graph, using vorder-push heuristic */
graphe::intpoly graphe::tutte_poly_recurse(int vc) {
    ++tutte_iter_count;
    intpoly p=poly_one(),fac;
    int n=node_count(),adj_sz,snv;
    bool isom;
    ipair e;
    vector<ipairs> blocks;
    ipairs E;
    graphe G(ctx,false),Gd(ctx,false),Gc(ctx,false);
    matrice L;
    ivector L_cp;
    map<ivector,vector<cpol> >::iterator ct;
    int *adj,*col;
    ulong *cg;
    size_t cg_sz;
    switch (vc) {
    case 2:
        assert(n>2);
        /* process cycles with multiedges */
        if (is_cycle(E)) {
            p=poly_geom(2,multiedges(E[n-1])+multiedges(E[n-2])+2,false,true);
            for (int i=1;i<=n-2;++i) {
                poly_mult(p,poly_geom(2,1+multiedges(E[i-1]),true));
            }
            for (int i=1;i<=n-2;++i) {
                intpoly q=poly_one();
                q[make_pair(0,0)]=1;
                for (int j=i+1;j<=n;++j) {
                    poly_mult(q,poly_geom(2,1+multiedges(E[j-1]),false,true));
                }
                for (int j=1;j<=i-1;++j) {
                    poly_mult(q,poly_geom(2,1+multiedges(E[j-1]),true));
                }
                poly_add(p,q);
            }
            break;
        }
        /* check for cached isomorphic graph, record the time */
        isom=false;
        tutte_time_start=clock();
        simplify(G,true);
        G.laplacian_matrix(L);
        L_cp=vecteur_2_vector_int(*_eval(symbolic(at_charpoly,L),ctx)._VECTptr); // charpoly of the Laplacian
        adj=G.to_array(adj_sz,false,true);
        snv=G.node_count();
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
        cg_sz=nautywrapper_words_needed(snv)*(size_t)snv;
        cg=new ulong[cg_sz];
        col=new int[snv];
        pthread_setcancelstate(PTHREAD_CANCEL_DISABLE,NULL);
        nautywrapper_canonical(0,snv,adj,NULL,cg,col);
        pthread_setcancelstate(PTHREAD_CANCEL_ENABLE,NULL);
        delete[] adj;
#endif
        if ((ct=cache.find(L_cp))!=cache.end()) {
            for (vector<cpol>::iterator it=ct->second.begin();it!=ct->second.end();++it) {
                cpol &cp=*it;
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
                if (cp.match(snv,cg,col)) {
                    ++cp.frq;
                    ++tutte_hits;
                    delete[] cg;
                    delete[] col;
                    p=cp.poly;
                    isom=true;
                    break;
                }
#else
                if (cp.match(snv,adj,adj_sz)) {
                    ++cp.frq;
                    ++tutte_hits;
                    delete[] adj;
                    p=cp.poly;
                    isom=true;
                    break;
                }
#endif
            }
        }
        tutte_matching_time+=double(clock()-tutte_time_start)/CLOCKS_PER_SEC;
        if (isom) break;
        /* no luck, perform the delete-contract step */
        get_edges_as_pairs(E);
        e=E.front();
        for (ipairs_iter it=E.begin()+1;it!=E.end();++it) {
            if (*it<e) e=*it;
        }
        fac=poly_geom(2,1+multiedges(e),true);
        set_multiedge(e,0);
        copy(Gd); copy(Gc); // copies to perform deletion-contraction on
        Gd.remove_edge(e);
        Gc.contract_edge(e.second,e.first,false);
        p=Gd.tutte_poly_recurse(1);
        poly_mult(fac,Gc.tutte_poly_recurse(1));
        poly_add(p,fac);
        /* cache the graph and its polynomial for future use */
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
        cache[L_cp].push_back(cpol(snv,cg,col,cg_sz,p));
#else
        cache[L_cp].push_back(cpol(snv,adj,adj_sz,p));
#endif
        break;
    case 1:
        find_blocks(blocks);
        for (vector<ipairs>::iterator it=blocks.begin();it!=blocks.end();++it) {
            ipairs &block=*it;
            assert(!block.empty());
            if (block.size()==1) { // bridge
                poly_mult(p,poly_geom(2,1+multiedges(block.front()),false,true));
            } else { // non-trivial biconnected subgraph
                std::sort(block.begin(),block.end());
                extract_subgraph(block,G);
                poly_mult(p,G.tutte_poly_recurse(2));
            }
        }
        break;
    default:
        assert(false);
    }
    return p;
}

/* return the Tutte polynomial of this graph */
gen graphe::tutte_polynomial(const gen &x,const gen &y) {
    assert(cache.empty() && !is_directed());
    tutte_iter_count=tutte_hits=0;
    tutte_matching_time=0;
    intpoly p;
    graphe G(ctx,false);
    ivector sigma;
    p=poly_one();
    connected_components_to_subgraphs();
    int sc=max_subgraph_index(),i;
    for (i=1;i<=sc;++i) {
        if (subgraph_size(i)<2)
            continue;
        sort_by_degrees(G,sigma,i);
        G.sharc_order();
        poly_mult(p,G.tutte_poly_recurse(1));
    }
    /* clear the cache and return the Tutte polynomial p */
    for (map<ivector,vector<cpol> >::iterator it=cache.begin();it!=cache.end();++it) {
        for (vector<cpol>::const_iterator jt=it->second.begin();jt!=it->second.end();++jt) {
#if defined HAVE_LIBNAUTY && defined HAVE_NAUTY_NAUTUTIL_H
            delete[] jt->cg;
            delete[] jt->col;
#else
            delete[] jt->adj;
#endif
        }
        it->second.clear();
    }
    cache.clear();
    return intpoly2gen(p,x,y);
}

/* find all fundamental cycles in this graph */
void graphe::fundamental_cycles(ivectors &cycles,int sg,gt_conn_check cc) {
    assert(!is_directed());
    ivectors components;
    int cmp,s;
    switch (cc) {
    case _GT_CC_CONNECTED:
        break;
    case _GT_CC_COMPONENTS_ARE_SUBGRAPHS:
        cmp=max_subgraph_index();
        for (s=1;s<=cmp;++s) {
            fundamental_cycles(cycles,s,_GT_CC_CONNECTED);
        }
        return;
    case _GT_CC_FIND_COMPONENTS:
        s=sg<0?1:1+max_subgraph_index();
        cmp=connected_components_to_subgraphs(sg);
        for (;s<=cmp;++s) {
            fundamental_cycles(cycles,s,_GT_CC_CONNECTED);
            renumber_subgraph(s,sg);
        }
        return;
    default: assert(false);
    }
    ipairs E,nte;
    get_edges_as_pairs(E,sg);
    int root=sg<0?0:first_vertex_from_subgraph(sg),i,j;
    dfs(root,true,true,NULL,sg);
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        i=it->first; j=it->second;
        if (node(i).ancestor()!=j && node(j).ancestor()!=i)
            nte.push_back(*it);
    }
    for (ipairs_iter it=nte.begin();it!=nte.end();++it) {
        const vertex &v=node(it->first),&w=node(it->second);
        i=v.disc()>w.disc()?it->first:it->second;
        j=v.disc()<w.disc()?it->first:it->second;
        ivector c;
        while (i!=j) {
            c.push_back(i);
            i=node(i).ancestor();
            assert(i>=0);
        }
        c.push_back(j);
        cycles.push_back(c);
    }
}

/* construct Mycielski graph with 2n+1 vertices and 3m+n edges */
void graphe::mycielskian(graphe &G) const {
    assert(!is_directed());
    ipairs E;
    get_edges_as_pairs(E);
    int lab=largest_integer_label(),n=node_count(),m=E.size();
    copy(G);
    ivector u(n);
    for (int i=0;i<n;++i) {
        u[i]=G.add_node(++lab);
    }
    int w=G.add_node(++lab);
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        G.add_edge(it->first,u[it->second]);
        G.add_edge(it->second,u[it->first]);
    }
    for (int i=0;i<n;++i) {
        G.add_edge(w,u[i]);
    }
    assert(G.node_count()==2*n+1 && G.edge_count()==3*m+n);
}

/* return the local clustering coefficient of the i-th vertex */
gen graphe::local_clustering_coeff(int i) const {
    assert(!is_directed());
    const vertex &v=node(i);
    const ivector &vn=v.neighbors();
    int d=vn.size(),cnt=0;
    for (ivector_iter it=vn.begin();it!=vn.end();++it) {
        const ivector &wn=node(*it).neighbors();
        cnt+=intersect_linear(vn.begin(),vn.end(),wn.begin(),wn.end());
    }
    return _ratnormal(fraction(cnt,d*(d-1)),ctx);
}

/* return the clustering coefficient of this graph */
gen graphe::clustering_coeff(bool approx,bool exact) {
    assert(!is_directed());
    int i,j,n=node_count();
    if (n<3) return 0;
    if (approx) {
        /* approximate within +-1/2*1e-2 with probability 0.9999 */
        int k=184207,l=0,u,v;
        for (i=0;i<k;++i) {
            j=rand_integer(n);
            const ivector &ngh=node(j).neighbors();
            u=ngh[rand_integer(ngh.size())];
            while (u==(v=ngh[rand_integer(ngh.size())]));
            if (has_edge(u,v)) ++l;
        }
        return _evalf(fraction(l,k),ctx);
    }
    return triangle_count(NULL,true,exact);
}

/* return the triangle density (transitivity) of this graph */
gen graphe::transitivity() {
    gen num_triangles(0),num_triplets(0);
    if (!is_directed()) {
        num_triangles=gen(3)*triangle_count();
        node_iter it=nodes.begin(),itend=nodes.end();
        for (;it!=itend;++it) {
            num_triplets+=comb(it->neighbors().size(),2);
        }
    } else {
        ipairs E;
        get_edges_as_pairs(E);
        int i,j;
        ipairs_iter it=E.begin(),itend=E.end();
        for (;it!=itend;++it) {
            i=it->first; j=it->second;
            const vertex &v=node(i),&w=node(j);
            const ivector &vn=v.neighbors(),&wn=w.neighbors();
            num_triangles+=intersect_linear(vn.begin(),vn.end(),wn.begin(),wn.end());
            num_triplets+=wn.size()-(w.has_neighbor(i)?1:0);
        }
    }
    if (giac::is_zero(num_triplets))
        return 0;
    return _ratnormal(fraction(num_triangles,num_triplets),ctx);
}

/* return the edge connectivity of this graph using Matula's algorithm */
int graphe::edge_connectivity() {
    int n=node_count();
    assert(n>=2 && !is_directed());
    set<int> D,A;
    vector<map<int,gen> > flow;
    int p,lambda=rand_max2,d,v,w,lambda_vw,maxdeg,i;
    /* set lambda to its upper bound */
    for (i=0;i<n;++i) {
        if ((d=degree(i))<lambda) {
            p=i;
            lambda=d;
        }
    }
    /* construct a dominating set D */
    v=p; maxdeg=lambda;
    for (i=0;i<n;++i) {
        A.insert(i);
        if (i!=v && (d=degree(i))>maxdeg) {
            maxdeg=d;
            v=i;
        }
    }
    set<int>::iterator sit;
    while (true) {
        D.insert(v);
        A.erase(A.find(v));
        const vertex &vv=node(v);
        for (ivector_iter it=vv.neighbors().begin();it!=vv.neighbors().end();++it) {
            sit=A.find(*it);
            if (sit!=A.end())
                A.erase(sit);
        }
        if (A.empty()) break;
        v=*A.begin();
    }
    /* find lambda(G) */
    v=*D.begin();
    D.erase(D.begin());
    set_directed(true);
    for (set<int>::const_iterator it=D.begin();it!=D.end();++it) {
        w=*it;
        lambda_vw=maxflow_edmonds_karp(v,w,flow).val;
        if (lambda_vw<lambda)
            lambda=lambda_vw;
    }
    set_directed(false);
    return lambda;
}

int graphe::vertex_pair_connectivity(int v,int w) {
    int n=node_count();
    graphe D(ctx,false);
    D.set_directed(true);
    D.add_nodes(2*n);
    int u1,u2;
    for (int i=0;i<n;++i) {
        u1=2*i; u2=2*i+1;
        if (i!=v || i!=w)
            D.add_edge(u1,u2);
        const vertex &u=node(i);
        for (ivector_iter it=u.neighbors().begin();it!=u.neighbors().end();++it) {
            if (i!=v) D.add_edge(*it,u1);
            if (i!=w) D.add_edge(u2,*it);
        }
    }
    vector<map<int,gen> > flow;
    return D.maxflow_edmonds_karp(2*v+1,2*w,flow).val;
}

/* return the vertex connectivity of an undirected graph */
int graphe::vertex_connectivity() {
    int n=node_count(),k=rand_max2,mindeg=rand_max2,deg,v;
    for (int i=0;i<n;++i) {
        if ((deg=degree(i))<mindeg) {
            v=i;
            mindeg=deg;
        }
    }
    for (int i=0;i<n;++i) {
        if (i==v || has_edge(i,v))
            continue;
        k=std::min(k,vertex_pair_connectivity(v,i));
    }
    ivector adj;
    adjacent_nodes(v,adj);
    for (ivector_iter it=adj.begin();it!=adj.end();++it) {
        for (ivector_iter jt=it+1;jt!=adj.end();++jt) {
            if (has_edge(*it,*jt)) continue;
            k=std::min(k,vertex_pair_connectivity(*it,*jt));
        }
    }
    return k;
}

/* truncate the graph, assuming that it is planar */
void graphe::truncate(graphe &dest,const ivectors &faces) {
    clear_embedding();
    set_embedding(faces);
    dest.clear();
    int nv=node_count(),ne=edge_count(),n;
    if (dest.supports_attributes()) {
        vecteur V;
        dest.make_default_labels(V,2*ne);
        dest.add_nodes(V);
    } else dest.add_nodes(2*ne);
    int cnt=0,f;
    ivector adj;
    map<ipair,ivector> p;
    for (int i=0;i<nv;++i) {
        vertex &v=node(i);
        n=v.neighbors().size();
        ivector ngh;
        adjacent_nodes(i,adj);
        ngh.push_back(adj.back());
        adj.pop_back();
        while (!adj.empty()) {
            f=v.edge_face(ngh.back());
            ivector::iterator it=adj.begin();
            for (;it!=adj.end() && f!=node(*it).edge_face(i);++it);
            assert(it!=adj.end());
            ngh.push_back(*it);
            adj.erase(it);
        }
        for (int j=0;j<n;++j) {
            p[make_pair(i<ngh[j]?i:ngh[j],i<ngh[j]?ngh[j]:i)].push_back(cnt+j);
            dest.add_edge(cnt+j,cnt+(j+1)%n);
        }
        cnt+=n;
    }
    ipairs E;
    get_edges_as_pairs(E);
    for (ipairs_iter it=E.begin();it!=E.end();++it) {
        const ivector &iv=p[*it];
        assert(iv.size()==2);
        dest.add_edge(iv.front(),iv.back());
    }
}

/*
 * CIRC_ENUM CLASS IMPLEMENTATION
 *
 * Implementation of Tarjan's algorithm published in SIAM J. Comp. 2(3), 1973
 */

graphe::circ_enum::circ_enum(graphe *gr,int lo,int hi) {
    G=gr;
    lb=lo; ub=hi;
    int n=G->node_count();
    A.resize(n);
    for (int i=0;i<n;++i) {
        const vertex &v=G->node(i);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            A[i].push_back(*it);
        }
    }
}

void graphe::circ_enum::backtrack(int v,bool &f) {
    bool g;
    f=false;
    point_stack.push_back(v);
    mark[v]=true;
    marked_stack.push(v);
    int w,u,len;
    for (int i=A[v].size();i-->0;) {
        w=A[v][i];
        if (w<s) A[v].erase(A[v].begin()+i);
        else if (w==s) {
            len=point_stack.size();
            if ((lb<0 || len>=lb) && (ub<0 || len<=ub))
                res.push_back(point_stack);
            f=true;
        } else if (!mark[w]) {
            backtrack(w,g);
            f=f||g;
        }
    }
    /* f is true if an elementary circuit containing
     * the partial path has been found */
    if (f) {
        while (marked_stack.top()!=v) {
            u=marked_stack.top();
            marked_stack.pop();
            mark[u]=false;
        }
        marked_stack.pop();
        mark[v]=false;
    }
    point_stack.pop_back();
}

/* return the list of all elementary cycles in the given digraph */
graphe::ivectors graphe::circ_enum::find_cycles() {
    int n=G->node_count(),u;
    bool flag;
    mark.resize(n,false);
    for (s=0;s<n;++s) {
        backtrack(s,flag);
        while (!marked_stack.empty()) {
            u=marked_stack.top();
            mark[u]=false;
            marked_stack.pop();
        }
    }
    return res;
}

/*
 * END OF CIRC_ENUM CLASS IMPLEMENTATION
 */

void graphe::elementary_cycles(ivectors &cyc,int lo,int hi) {
    assert(is_directed());
    circ_enum ce(this,lo,hi);
    cyc=ce.find_cycles();
}

/*
 * YEN CLASS IMPLEMENTATION
 */

void graphe::yen::delete_children(tree_node *r) {
    for (vector<tree_node*>::const_iterator it=r->children.begin();it!=r->children.end();++it) {
        delete_children(*it);
        delete *it;
    }
}

graphe::yen::tree_node *graphe::yen::add_tree_node(tree_node *p) {
    if (p==NULL) {
        root=new tree_node;
        root->parent=NULL;
        root->i=src;
        return root;
    }
    tree_node *t;
    p->children.push_back(t=new tree_node);
    t->parent=p;
    return t;
}

graphe::yen::~yen() {
    if (root!=NULL) {
        delete_children(root);
        delete root;
    }
}

graphe::yen::tree_node *graphe::yen::store_path(const ivector &path,tree_node *r) {
    tree_node *t=r,*next;
    int i,n=path.size();
    for (i=1;i<n;++i) {
        next=NULL;
        for (vector<tree_node*>::const_iterator it=t->children.begin();it!=t->children.end();++it) {
            if ((*it)->i==path[i]) {
                next=*it;
                break;
            }
        }
        if (next!=NULL) {
            t=next;
            continue;
        }
        t=add_tree_node(t);
        t->i=path[i];
        t->selected=false;
    }
    return t;
}

void graphe::yen::select_path(tree_node *p) {
    kspaths.push_back(p);
    tree_node *t=p;
    while (t!=NULL) {
        if (t->selected) break;
        t->selected=true;
        t=t->parent;
    }
}

void graphe::yen::restore_path(tree_node *p,ivector &path) {
    tree_node *t=p;
    path.clear();
    while (t!=NULL) {
        path.push_back(t->i);
        t=t->parent;
    }
    std::reverse(path.begin(),path.end());
}

void graphe::yen::find_kspaths(ivectors &paths) {
    assert(G->supports_attributes());
    ivectors cp;
    vecteur pw;
    matrice W;
    if (!G->is_weighted()) {
        G->adjacency_matrix(W);
        G->make_weighted(W);
    }
    G->weight_matrix(W);
    paths.clear();
    G->dijkstra(src,ivector(1,dest),pw,&cp);
    if (cp.front().empty()) return;
    kspaths.clear();
    tree_node *p=store_path(cp.front(),add_tree_node(NULL)),*t,*next;
    select_path(p);
    stack<ipair> removed_edges;
    gen wg;
    int i,j,spur_node,tail,head;
    set<pair<gen,tree_node*>,kspaths_comparator> candidates;
    set<pair<gen,tree_node*>,kspaths_comparator>::const_iterator cit;
    ivector path;
    path.reserve(G->node_count());
    G->save_subgraphs();
    G->unset_subgraphs(1);
    for (int k=1;k<K;++k) {
        restore_path(p,path);
        t=root;
        for (i=0;i+1<(int)path.size();++i) {
            spur_node=path[i];
            for (vector<tree_node*>::const_iterator it=t->children.begin();it!=t->children.end();++it) {
                if ((*it)->selected) {
                    head=(*it)->i; tail=(*it)->parent->i;
                    G->set_edge_attribute(tail,head,_GT_ATTRIB_WEIGHT,plusinf());
                    removed_edges.push(make_pair(tail,head));
                    if (head==path[i+1]) next=*it;
                }
            }
            G->dijkstra(spur_node,ivector(1,dest),pw,&cp,1);
            if (!cp.front().empty()) {
                wg=pw.front();
                for (j=0;j<i;++j) wg+=W[path[j]][path[j+1]];
                candidates.insert(make_pair(wg,store_path(cp.front(),t)));
            }
            G->node(spur_node).set_subgraph(0);
            t=next;
        }
        for (ivector_iter it=path.begin();it+1!=path.end();++it) {
            G->node(*it).set_subgraph(1);
        }
        while (!removed_edges.empty()) {
            i=removed_edges.top().first; j=removed_edges.top().second;
            G->set_edge_attribute(i,j,_GT_ATTRIB_WEIGHT,W[i][j]);
            removed_edges.pop();
        }
        if (candidates.empty()) break;
        cit=candidates.begin();
        p=cit->second;
        select_path(p);
        candidates.erase(cit);
    }
    G->restore_subgraphs();
    for (vector<tree_node*>::const_iterator it=kspaths.begin();it!=kspaths.end();++it) {
        restore_path(*it,path);
        paths.push_back(path);
    }
}

/*
 * END OF YEN CLASS
 */

void graphe::yen_ksp(int K,int src,int dest,ivectors &paths) {
    yen Y(this,src,dest,K);
    Y.find_kspaths(paths);
}

/*
 * MM CLASS IMPLEMENTATION
 * Coded as described in lecture notes by Uri Zwick:
 * https://www.cs.tau.ac.il/~zwick/grad-algo-0910/match.pdf
 */

graphe::mm::mm(graphe *g) {
    G=g;
    V=G->node_count();
    mate=new int[V];
    label=new int[V];
    pred=new int[V];
    bridge=new int[2*V];
    ds=new unionfind(V);
    ap.reserve(V);
}

graphe::mm::~mm() {
    delete[] mate;
    delete[] label;
    delete[] pred;
    delete[] bridge;
    delete ds;
}

bool graphe::mm::alternating_forest() {
    ds->clear();
    for (int v=0;v<V;++v) {
        if (s>=0 && G->node(v).subgraph()!=s)
            continue;
        label[v]=pred[v]=-1;
        ds->make_set(v);
    }
    for (int v=0;v<V;++v) {
        if (s>=0 && G->node(v).subgraph()!=s)
            continue;
        if (mate[v]<0 && alternating_tree(v))
            return true;
    }
    return false;
}

bool graphe::mm::alternating_tree(int v) {
    label[v]=EVEN;
    while (!Q.empty()) Q.pop();
    Q.push(v);
    while (!Q.empty()) {
        v=Q.front();
        Q.pop();
        const ivector &ngh=G->node(v).neighbors();
        for (ivector_iter it=ngh.begin();it!=ngh.end();++it) {
            if (s>=0 && G->node(*it).subgraph()!=s)
                continue;
            if (examine(v,*it))
                return true;
        }
    }
    return false;
}

bool graphe::mm::examine(int v,int w) {
    int vb=ds->find(v),wb=ds->find(w);
    if (vb!=wb) {
        if (label[wb]<0) {
            if (mate[w]<0) {
                augmenting_path(v,w);
                return true;
            } else extend_tree(v,w);
        } else if (label[wb]==EVEN)
            shrink_blossom(v,w);
    }
    return false;
}

void graphe::mm::extend_tree(int v,int w) {
    label[w]=ODD;
    pred[w]=v;
    label[mate[w]]=EVEN;
    pred[mate[w]]=w;
    Q.push(mate[w]);
}

void graphe::mm::shrink_blossom(int v,int w) {
    int b=find_base(v,w);
    shrink_path(b,v,w);
    shrink_path(b,w,v);
}

void graphe::mm::shrink_path(int b,int v,int w) {
    int u=ds->find(v);
    while (u!=b) {
        ds->unite(b,u);
        u=mate[u];
        ds->unite(b,u);
        ds->select(b);
        Q.push(u);
        bridge[2*u]=v;
        bridge[2*u+1]=w;
        u=ds->find(pred[u]);
    }
}

void graphe::mm::augmenting_path(int v,int w) {
    int r=find_root(v);
    ap.clear();
    ap.push_back(w);
    ivector p;
    p=find_path(v,r);
    ap.insert(ap.end(),p.begin(),p.end());
}

graphe::ivector graphe::mm::find_path(int s,int t) {
    ivector p(1,s);
    if (s==t)
        ;
    else if (label[s]==EVEN) {
        p.push_back(mate[s]);
        ivector q=find_path(pred[mate[s]],t);
        p.insert(p.end(),q.begin(),q.end());
    } else {
        int v=bridge[2*s],w=bridge[2*s+1];
        ivector q1=find_path(v,mate[s]);
        ivector q2=find_path(w,t);
        p.insert(p.end(),q1.rbegin(),q1.rend());
        p.insert(p.end(),q2.begin(),q2.end());
    }
    return p;
}

int graphe::mm::find_root(int v) {
    int u=v,i;
    while ((i=pred[u])>=0) u=i;
    return ds->find(u);
}

int graphe::mm::find_base(int v,int w) {
    iset vst;
    int i=v,j=w;
    do vst.insert(i); while ((i=pred[i])>=0);
    while (vst.find(j)==vst.end())
        assert((j=pred[j])>=0);
    return ds->find(j);
}

void graphe::mm::find_maximum_matching(ipairs &matching,int sg) {
    s=sg;
    /* start from a maximal matching */
    G->maximal_matching(matching,s);
    for (int i=0;i<V;++i) {
        mate[i]=-1;
    }
    for (ipairs_iter it=matching.begin();it!=matching.end();++it) {
        mate[it->first]=it->second;
        mate[it->second]=it->first;
    }
    /* improve the matching by finding augmenting paths */
    while (alternating_forest()) {
        assert(ap.size()%2==0);
        for (ivector_iter it=ap.begin();it!=ap.end();++it) {
            mate[*it]=(it-ap.begin())%2==0?*(it+1):*(it-1);
        }
    }
    /* output maximum matching */
    matching.clear();
    for (int i=0;i<V;++i) {
        if (mate[i]>i)
            matching.push_back(std::make_pair(i,mate[i]));
    }
}

/*
 * END OF MM CLASS
 */

/* 
 * CENTRALITY MEASURES
 */

/* return the information centrality measure for k-th vertex,
 * if k<0 return the list of IC measures for all vertices */
gen graphe::information_centrality(int k,bool approx) const {
    int n=node_count();
    assert(n>1 && !is_directed());
    matrice A;
    laplacian_matrix(A);
    for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) {
            A[i]._VECTptr->at(j)+=1;
        }
    }
    if (approx)
        A=*_evalf(A,ctx)._VECTptr;
    gen dt=mdet(A,ctx);
    if (approx) {
        if (!is_real_number(dt,ctx)) {
            message("Error: input must be numerical");
            return undef;
        }
        dt=to_real_number(dt,ctx);
    }
    if (giac::is_zero(dt)) {
        message("Error: incremented Laplacian is singular");
        return undef;
    }
    matrice B=minv(A,ctx);
    if (k>=0) {
        assert(k<n);
        gen num(n),den(0);
        for (int i=0;i<n;++i) {
            den+=B[i]._VECTptr->at(i)+B[k]._VECTptr->at(k)-2*B[i]._VECTptr->at(k);
        }
        return num/den;
    }
    vecteur ic(n,0);
    for (int i=0;i<n;++i) {
        for (int j=0;j<n;++j) {
            ic[i]+=B[j]._VECTptr->at(j)+B[i]._VECTptr->at(i)-2*B[j]._VECTptr->at(i);
        }
        ic[i]=gen(n)/ic[i];
    }
    return ic;
}

/* return the normalized degree centrality measure for k-th vertex,
 * if k<0 return the list of NDC measures for all vertices */
gen graphe::degree_centrality(int k) const {
    int n=node_count();
    assert(n>1);
    vecteur dv(n);
    if (k>=0) {
        assert(k<n);
        return gen(degree(k))/gen(n-1);
    }
    for (int i=0;i<n;++i) {
        dv[i]=gen(degree(i))/gen(n-1);
    }
    return dv;
}

/* return the list of shortest distances from
 * k-th vertex to every other vertex */ 
vecteur graphe::distances_from(int k) {
    int n=node_count();
    assert(k>=0 && k<n);
    vecteur res(n);
    if (is_weighted()) {
        ivector dest;
        for (int i=0;i<n;++i) {
            dest.push_back(i);
        }
        dijkstra(k,dest,res);
    } else {
        bfs(k);
        res[k]=0;
        for (int i=0;i<n;++i) {
            gen d(0);
            if (i==k)
                continue;
            if (node(i).ancestor()<0)
                res[i]=plusinf();
            else {
                int j=i;
                while ((j=node(j).ancestor())>=0) {
                    d+=1;
                }
                res[i]=d;
            }
        }
    }
    return res;
}

/* return the closeness centrality measure for k-th vertex,
 * if k<0 return the list of HC measures for all vertices,
 * if harmonic=true, return harmonic centrality */
gen graphe::closeness_centrality(int k,bool harmonic) {
    int n=node_count();
    assert(n>1);
    if (k>=0) assert(k<n);
    vecteur hc(n,0);
    if (is_weighted()) {
        if (k>=0) {
            graphe H(*this);
            if (is_directed())
                reverse(H);
            vecteur pw=H.distances_from(k);
            pw.erase(pw.begin()+k);
            if (harmonic)
                hc[k]=_sum(_apply(makesequence(at_inv,pw),ctx),ctx);
            else
                hc[k]=gen(n-1)/_sum(pw,ctx);
        } else {
            matrice apd;
            allpairs_distance(apd);
            if (harmonic) {
                apd=*_apply(makesequence(at_inv,apd),ctx)._VECTptr;
                hc=*_sum(apd,ctx)._VECTptr;
            } else {
                hc=*_sum(apd,ctx)._VECTptr;
                hc=multvecteur(n-1,*_apply(makesequence(at_inv),ctx)._VECTptr);
            }
        }
    } else {
        graphe H(*this);
        if (is_directed())
            reverse(H);
        for (int i=0;i<n;++i) {
            if (k>=0 && i!=k)
                continue;
            vecteur dv=H.distances_from(i);
            dv.erase(dv.begin()+i);
            if (harmonic)
                hc[i]=_sum(_apply(makesequence(at_inv,dv),ctx),ctx);
            else
                hc[i]=gen(n-1)/_sum(dv,ctx);
        }
    }
    if (k>=0)
        return hc[k];
    return hc;
}

/* return the list of betweenness centrality measures for all vertices
 * using Brandes' algorithm */
gen graphe::betweenness_centrality(int k) const {
    int n=node_count();
    assert(n>1 && (k<0 || k<n));
    vecteur cb(n,0);
    std::stack<int> S;
    std::queue<int> Q;
    std::map<int,ivector> P;
    ivector sigma(n),d(n);
    vecteur delta(n);
    for (int s=0;s<n;++s) {
        P.clear();
        for (int i=0;i<n;++i) {
            sigma[i]=i==s?1:0;
            d[i]=i==s?0:-1;
            delta[i]=0;
        }
        Q.push(s);
        while (!Q.empty()) {
            int v=Q.front();
            Q.pop();
            S.push(v);
            for (ivector_iter it=node(v).neighbors().begin();it!=node(v).neighbors().end();++it) {
                int w=*it;
                if (d[w]<0) {
                    Q.push(w);
                    d[w]=d[v]+1;
                }
                if (d[w]==d[v]+1) {
                    sigma[w]+=sigma[v];
                    P[w].push_back(v);
                }
            }
        }
        while (!S.empty()) {
            int w=S.top();
            S.pop();
            for (ivector_iter it=P[w].begin();it!=P[w].end();++it) {
                int v=*it;
                delta[v]+=gen(sigma[v])/gen(sigma[w])*(gen(1)+delta[w]);
            }
            if (w!=s)
                cb[w]+=delta[w];
        }
    }
    if (!is_directed())
        cb=multvecteur(fraction(1,2),cb);
    if (k>=0)
        return cb[k];
    return cb;
}

/* return the list of communicability betweenness centrality for all vertices */
gen graphe::communicability_betweenness_centrality(int k) const {
    int n=node_count();
    matrice A;
    vecteur row(n),col(n),res(n,0),nullvec(n,0);
    adjacency_matrix(A);
    matrice eA=*exp(_evalf(A,ctx),ctx)._VECTptr;
    for (int r=0;r<n;++r) {
        if (k>=0 && r!=k) continue;
        row=*A[r]._VECTptr;
        A[r]=nullvec;
        matrice eAk=*exp(_evalf(A,ctx),ctx)._VECTptr;
        for (int p=0;p<n;++p) {
            if (p==r) continue;
            for (int q=0;q<n;++q) {
                if (q==p || q==r) continue;
                res[r]+=1-eAk[p][q]/eA[p][q];
            }
        }
        A[r]=row;
    }
    res=divvecteur(res,(n-1)*(n-2));
    if (k>=0)
        return res[k];
    return res;
}
 
/* return the list of Katz centrality measures for all vertices,
 * where att is attenuation factor */
vecteur graphe::katz_centrality(const gen &att) const {
    int n=node_count();
    matrice A;
    adjacency_matrix(A);
    return *_sum(mtran(minv(midn(n)-mscale(att,mtran(A)),ctx)),ctx)._VECTptr;
}

/* returns the splittance of this graph */
int graphe::splittance(int &m,ivector &vseq) const {
    assert(!is_directed() && !is_empty());
    int n=node_count();
    vecteur ds=vector_int_2_vecteur(degree_sequence()),ind(n);
    for (int i=0;i<n;++i) {
        ind[i]=i;
    }
    matrice dsorted=*_tran(_revlist(_sort(_zip(makesequence(ds,ind),ctx),ctx),ctx),ctx)._VECTptr;
    const vecteur &d=*dsorted[0]._VECTptr;
    const vecteur &v=*dsorted[1]._VECTptr;
    assert(is_integer_vecteur(v));
    m=-1;
    vseq.resize(n);
    for (int i=0;i<n;++i) {
        if (is_greater(d[i],i,ctx))
            m=i;
        vseq[i]=v[i].val;
    }
    assert(m>0);
    gen res=gen(m*(m+1));
    for (int i=0;i<n;++i) {
        res+=i>m?d[i]:-d[i];    
    }
    res=res/gen(2);
    assert(res.is_integer());
    return res.val;
}

/* return true iff this is split graph and store the clique and
 * independent set parts in clq and indp in that case */
bool graphe::is_split_graph(ivector &clq,ivector &indp) const {
    int m;
    ivector vseq;
    if (splittance(m,vseq)!=0)
        return false;
    clq=ivector(vseq.begin(),vseq.begin()+m+1);
    indp=ivector(vseq.begin()+m+1,vseq.end());
    return true;
}

/* construct the graph G which is equal to this graph with
 * subgraph sg contracted into a single vertex */
void graphe::contract_subgraph(graphe &G,const ivector &sg,const gen &lb) const {
    assert(supports_attributes());
    bool isdir=is_directed();
    G.set_directed(isdir);
    G.set_graph_attributes(attributes);
    G.reserve_nodes(node_count()-sg.size()+1);
    int cn=-1,n=node_count(),i;
    bvector in_sg(n,false);
    for (ivector_iter it=sg.begin();it!=sg.end();++it) {
        in_sg[*it]=true;
    }
    ivector vmap(n);
    node_iter nt=nodes.begin(),ntend=nodes.end();
    for (i=0;nt!=ntend;++nt,++i) {
        if (!in_sg[i])
            vmap[i]=G.add_node(nt->label(),nt->attributes());
        else {
            if (cn<0)
                cn=G.add_node(is_undef(lb)?get_node_labels(sg):lb);
            vmap[i]=cn;
        }
    }
    ipairs E;
    get_edges_as_pairs(E);
    ipairs_iter it=E.begin(),itend=E.end();
    for (;it!=itend;++it) {
        int u=it->first,v=it->second;
        if (in_sg[u] && in_sg[v])
            continue;
        G.add_edge(vmap[u],vmap[v],node(u).neighbor_attributes(v));
    }
}

/* construct routine for GRASP clique finding algorithm */
void graphe::grasp_construct(double alpha,ivector &Q,bool cmpl,int sg) {
    Q.clear();
    int n=node_count();
    ivector V,deg(n);
    iset D;
    if (sg>=0) {
        get_subgraph(sg,V);
        n=V.size();
    } else {
        V.resize(n);
        for (int i=0;i<n;++i) V[i]=i;
    }
    ivector C(n);
    for (int i=0;i<n;++i) {
        C[i]=sg<0?i:V[i];
    }
    int u,d,dg,d_min,d_max,s=max_subgraph_index();
    ivector rcl;
    rcl.reserve(n);
    while (!C.empty()) {
        set_subgraph(C,++s);
        d_min=0;
        d_max=0;
        for (ivector_iter it=C.begin();it!=C.end();++it) {
            dg=degree(*it,s);
            d=deg[*it]=cmpl?C.size()-dg-1:dg;
            if (d_min==0 || d<d_min)
                d_min=d;
            if (d_max==0 || d>d_max)
                d_max=d;
        }
        rcl.clear();
        for (ivector_iter it=C.begin();it!=C.end();++it) {
            if (deg[*it]>=d_min+alpha*(d_max-d_min))
                rcl.push_back(*it);
        }
        u=rcl[_rand(rcl.size(),ctx).val];
        insert_sorted(Q,u);
        C.clear();
        const vertex &v=node(u);
        const ivector &N=v.neighbors();
        if (cmpl) {
            sets_difference(V,N,D);
            for (iset_iter it=D.begin();it!=D.end();++it) {
                if (*it!=u && node(*it).subgraph()==s)
                    C.push_back(*it);
            }
        } else for (ivector_iter it=N.begin();it!=N.end();++it) {
            if (node(*it).subgraph()==s)
                C.push_back(*it);
        }
    }
    if (sg<0)
        unset_subgraphs();
    else set_subgraph(V,sg);
}

/* local search routine for GRASP clique finding algorithm */
void graphe::grasp_local(ivector &Q,bool cmpl,int sg) {
    iset lv,V,D,res;
    ivectors H;
    int nc;
    node_iter nt,ntend;
    if (cmpl) for (nc=0,nt=nodes.begin(),ntend=nodes.end();nt!=ntend;++nt,++nc)
        if (sg<0 || nt->subgraph()==sg)
            V.insert(nc);
    int tmp_sg=max_subgraph_index()+1;
    save_subgraphs();
    if (sg<0) {
        unset_subgraphs(0);
        set_subgraph(Q,1);
    } else set_subgraph(Q,tmp_sg);
    while (true) {
        lv.clear();
        for (ivector_iter it=Q.begin();it!=Q.end();++it) {
            const vertex &v=node(*it);
            const ivector &N=v.neighbors();
            if (cmpl) {
                sets_difference(V,N,D);
                for (iset_iter jt=D.begin();jt!=D.end();++jt) {
                    if (*it!=*jt && node(*jt).subgraph()==(sg<0?0:sg))
                        lv.insert(*jt);
                }
            } else for (ivector_iter jt=N.begin();jt!=N.end();++jt) {
                if (node(*jt).subgraph()==(sg<0?0:sg))
                    lv.insert(*jt);
            }
        }
        H.clear();
        for (ivector_iter wt=Q.begin();wt!=Q.end();++wt) {
            for (iset_iter it=lv.begin();it!=lv.end();++it) {
                if (cmpl!=nodes_are_adjacent(*it,*wt))
                    continue;
                const ivector &ni=node(*it).neighbors();
                if (cmpl) {
                    sets_difference(V,ni,D);
                    D.erase(*it);
                }
                size_t isz=cmpl?sets_intersection(Q,D,res):intersect_linear(ni.begin(),ni.end(),Q.begin(),Q.end());
                if (isz<Q.size()-1)
                    continue;
                for (iset_iter jt=lv.begin();jt!=lv.end();++jt) {
                    if (*it<*jt && cmpl!=nodes_are_adjacent(*it,*jt) && cmpl==nodes_are_adjacent(*jt,*wt)) {
                        const ivector &nj=node(*jt).neighbors();
                        if (cmpl) {
                            sets_difference(V,nj,D);
                            D.erase(*jt);
                        }
                        isz=cmpl?sets_intersection(Q,D,res):intersect_linear(nj.begin(),nj.end(),Q.begin(),Q.end());
                        if (isz==Q.size()-1) {
                            ivector h(3);
                            h[0]=*it;
                            h[1]=*jt;
                            h[2]=*wt;
                            H.push_back(h);
                        }
                    }
                }
            }
        }
        if (H.empty())
            break;
        const ivector &h=H[_rand(H.size(),ctx).val];
        node(h[0]).set_subgraph(sg<0?1:tmp_sg);
        node(h[1]).set_subgraph(sg<0?1:tmp_sg);
        node(h[2]).set_subgraph(sg<0?0:sg);
        erase_sorted(Q,h[2]);
        insert_sorted(Q,h[0]);
        insert_sorted(Q,h[1]);
    }
    restore_subgraphs();
}

/* GRASP clique finding algorithm (Abello & Pardalos, 1998) with
 * maxitr iterations operating on subgraph sg, the result is
 * written in clq and sorted in ascending order.
 * If cmpl is true, then the search is performed in the complement
 * graph, i.e. a maximal independent set is constructed. */
void graphe::grasp_clique(int maxitr,ivector &clq,bool cmpl,int sg) {
    int clq_size=0;
    assert(maxitr>0);
    ivector Q;
    for (int k=0;k<maxitr;++k) {
        double alpha=_rand(makesequence(0,1),ctx).DOUBLE_val();
        grasp_construct(alpha,Q,cmpl,sg);
        grasp_local(Q,cmpl,sg);
        if (int(Q.size())>clq_size) {
            clq_size=Q.size();
            clq=Q;
        }
    }
}

int graphe::count_edges_in_Nv(int v,int sg) const {
    ivector nv=node(v).neighbors();
    if (sg>=0) for (int i=nv.size();i-->0;) {
        if (node(nv[i]).subgraph()!=sg)
            nv.erase(nv.begin()+i);
    }
    int ec=0;
    for (ivector_iter it=nv.begin();it!=nv.end();++it) {
        const ivector &nw=node(*it).neighbors();
        ec+=intersect_linear(nv.begin(),nv.end(),nw.begin(),nw.end());
    }
    if (!is_directed()) {
        assert(ec%2==0);
        ec/=2;
    }
    return ec;
}

int graphe::count_edges(const ivector &V) const {
    int ec=0;
    for (ivector_iter it=V.begin();it!=V.end();++it) {
        const ivector &nv=node(*it).neighbors();
        ec+=intersect_linear(V.begin(),V.end(),nv.begin(),nv.end());
    }
    if (!is_directed()) {
        assert(ec%2==0);
        ec/=2;
    }
    return ec;
}

/*
 * MVC_SOLVER CLASS IMPLEMENTATION
 */

#ifdef HAVE_LIBGLPK
 
/* constructor */
graphe::mvc_solver::mvc_solver(graphe *gr,int s) {
    G=gr;
    sg=s;
    assert(sg!=0);
    ilp=NULL;
}

/* destructor */
graphe::mvc_solver::~mvc_solver() {
    if (ilp!=NULL)
        glp_delete_prob(ilp);
}

/* exact algorithm for minimal vertex cover */
int graphe::mvc_solver::solve(ivector &cover,int k,int tm_lim,double gap_tol,bool verbose) {
    /* formulate MVCP as ILP */
    G->get_edges_as_pairs(edges,sg);
    if (sg>=0) {
        G->get_subgraph(sg,V);
        V_pos.resize(G->node_count(),-1);
        for (ivector_iter it=V.begin();it!=V.end();++it) {
            V_pos[*it]=it-V.begin();
        }
    }
    int n=sg<0?G->node_count():V.size(),c=0,m=edges.size();
    int *ia=new int[1+2*m+n],*ja=new int[1+2*m+n];
    double *ar=new double[1+2*m+n];
    heur_sol=new double[n+1];
    if (ilp==NULL)
        ilp=glp_create_prob();
    else glp_erase_prob(ilp);
    glp_set_obj_dir(ilp,GLP_MIN);
    glp_add_rows(ilp,m);
    for (int i=1;i<=m;++i) {
        glp_set_row_bnds(ilp,i,GLP_LO,1.0,0.0);
        const ipair &e=edges[i-1];
        ia[++c]=i; ja[c]=(sg<0?e.first:V_pos[e.first])+1; ar[c]=1.0;
        ia[++c]=i; ja[c]=(sg<0?e.second:V_pos[e.second])+1; ar[c]=1.0;
    }
    last_row=glp_add_rows(ilp,1);
    if (k>=0)
        glp_set_row_bnds(ilp,last_row,GLP_FX,k,k);
    else {
        mm mm_finder(G);
        ipairs matching;
        mm_finder.find_maximum_matching(matching,sg);
        glp_set_row_bnds(ilp,last_row,GLP_LO,matching.size(),0);
    }
    for (int j=1;j<=n;++j) {
        ia[++c]=last_row; ja[c]=j; ar[c]=1.0;
    }
    glp_add_cols(ilp,n);
    for (int i=1;i<=n;++i) {
        glp_set_col_kind(ilp,i,GLP_BV);
        if (k<0)
            glp_set_obj_coef(ilp,i,1.0);
    }
    glp_load_matrix(ilp,c,ia,ja,ar);
    /* solver initialization */
    glp_smcp lparm;
    glp_init_smcp(&lparm);
    lparm.msg_lev=GLP_MSG_OFF;
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.msg_lev=GLP_MSG_ON;
    parm.presolve=GLP_OFF;
    parm.sr_heur=GLP_ON;
    parm.fp_heur=GLP_OFF;
    parm.ps_heur=GLP_OFF;
    parm.gmi_cuts=GLP_ON;
    parm.cov_cuts=GLP_OFF;
    parm.mir_cuts=GLP_ON;
    parm.clq_cuts=GLP_ON;
    parm.br_tech=k<0?GLP_BR_PCH:GLP_BR_MFV;
    parm.bt_tech=k<0?GLP_BT_BLB:GLP_BT_DFS;
    parm.cb_func=&callback;
    parm.cb_info=static_cast<void*>(this);
    parm.cb_size=sizeof(int);
    if (tm_lim>0) parm.tm_lim=tm_lim;
    parm.mip_gap=gap_tol;
    compute_heur=k<0;
    is_k_vc=k>=0;
    /* solve */
    int term_old;
    if (verbose)
        glp_term_hook(term_hook,static_cast<void*>(logptr(G->giac_context())));
    else term_old=glp_term_out(GLP_OFF);
    if (parm.presolve==GLP_OFF)
        glp_simplex(ilp,&lparm);
    int res=glp_intopt(ilp,&parm);
    if (verbose)
        glp_term_hook(NULL,NULL);
    else glp_term_out(term_old);
    delete[] ia;
    delete[] ja;
    delete[] ar;
    delete[] heur_sol;
    switch (res) {
    case GLP_EBOUND:
    case GLP_EROOT:
    case GLP_ENOPFS:
    case GLP_ENODFS:
    case GLP_EFAIL:
        if (k<0 || res==GLP_EFAIL) {
            if (verbose)
                G->message("Error: MIP solver failure");
            return -1;
        }
        return 1;
    case GLP_ESTOP:
        return -1;
    default:
        break;
    }
    switch (glp_mip_status(ilp)) {
    case GLP_FEAS:
        if (k<0 && verbose)
            G->message("Warning: the solution is not necessarily optimal");
    case GLP_OPT:
        /* extract solution */
        for (int i=0;i<n;++i) {
            if ((int)round(glp_mip_col_val(ilp,i+1))==1)
                cover.push_back(sg<0?i:V[i]);
        }
        return 0;
    case GLP_UNDEF:
        if (verbose)
            G->message("Error: solution is undefined");
        return -1;
    case GLP_NOFEAS:
        if (k<0) {
            if (verbose)
                G->message("Error: unable to find a feasible solution");
            return -1;
        }
        return 1;
    default:
        break;
    }
    return -1; // should not be reachable
}

void graphe::mvc_solver::find_mirrors(int v) {
    iset n2v,D;
    ivector nv=G->node(v).neighbors();
    if (sg>=0) for (int i=nv.size();i-->0;) {
        if (G->node(nv[i]).subgraph()!=sg)
            nv.erase(nv.begin()+i);
    }
    for (ivector_iter it=nv.begin();it!=nv.end();++it) {
        const ivector &nw=G->node(*it).neighbors();
        for (ivector_iter jt=nw.begin();jt!=nw.end();++jt) {
            if (sg<0 || G->node(*jt).subgraph()==sg)
                n2v.insert(*jt);
        }
    }
    for (iset_iter it=n2v.begin();it!=n2v.end();++it) {
        sets_difference(nv,G->node(*it).neighbors(),D);
        int d=D.size();
        ivector dv(d);
        int c=0;
        for (iset_iter jt=D.begin();jt!=D.end();++jt) dv[c++]=*it;
        if (d<2 || G->count_edges(dv)==(d*(d-1))/2)
            mirrors.push_back(*it);
    }
    mirrors.push_back(-1);
}

void graphe::mvc_solver::branch(glp_tree *tree) {
    int n=sg<0?G->node_count():V.size();
    ipairs can_branch;
    for (int i=0;i<n;++i) {
        if (glp_ios_can_branch(tree,i+1))
            can_branch.push_back(std::make_pair(G->degree(sg<0?i:V[i],sg),i));
    }
    std::sort(can_branch.rbegin(),can_branch.rend());
    ivector cand;
    int maxd=can_branch.begin()->first;
    for (ipairs_iter it=can_branch.begin();it!=can_branch.end();++it) {
        if (it->first==maxd)
            cand.push_back(it->second);
    }
    int ch=-1;
    if (cand.size()>1) {
        int min_me=0;
        for (ivector_iter it=cand.begin();it!=cand.end();++it) {
            int me=G->count_edges_in_Nv(sg<0?*it:V[*it],sg);
            if (ch<0 || me<min_me) {
                ch=*it;
                min_me=me;
            }
        }
    } else ch=cand.front();
    assert(ch>=0);
    *static_cast<int*>(glp_ios_node_data(tree,glp_ios_curr_node(tree)))=mirrors.size();
    mirrors.push_back(ch);
    find_mirrors(sg<0?ch:V[ch]);
    glp_ios_branch_upon(tree,ch+1,GLP_UP_BRNCH);
}

int graphe::mvc_solver::heuristic(glp_tree *tree) {
    glp_prob *p=glp_ios_get_prob(tree);
    if (sg<0)
        G->unset_subgraphs(1);
    int n=sg<0?G->node_count():V.size();
    ivector indp,cov=V;
    for (int i=0;i<n;++i) {
        int j=sg<0?i:V[i];
        if (glp_get_col_type(p,i+1)==GLP_FX) {
            G->node(j).set_subgraph(0);
            if (glp_get_col_ub(p,i+1)==0)
                cov.erase(std::find(cov.begin(),cov.end(),j));
        }
    }
    G->grasp_clique(3,indp,true,sg);
    for (ivector_iter it=indp.begin();it!=indp.end();++it) {
        cov.erase(std::find(cov.begin(),cov.end(),*it));
    }
    for (int i=1;i<=n;++i) {
        heur_sol[i]=0;
    }
    for (ivector_iter it=cov.begin();it!=cov.end();++it) {
        heur_sol[1+V_pos[*it]]=1;
    }
    return cov.size();
}

/* lower bound for MVC according to Luzhi Wang et al. (2019) */
int graphe::mvc_solver::lower_bound(int s) {
    /* clique lower bound */
    ivector clq;
    int sg_size=G->subgraph_size(s),lb_clq=0;
    while (sg_size>0) {
        G->grasp_clique(3,clq,false,s);
        lb_clq+=clq.size()-1;
        G->set_subgraph(clq,0);
        sg_size-=clq.size();
    }
    return lb_clq;
}

bool graphe::mvc_solver::is_vertex_fixed(glp_prob *p,int j,bool &in_cover) {
    int i=(sg<0?j:V_pos[j])+1;
    if (glp_get_col_type(p,i)==GLP_FX) {
        in_cover=glp_get_col_lb(p,i)==1;
        return true;
    }
    return false;
}

void graphe::mvc_solver::make_vertex_fixed(glp_prob *p,int j,bool in_cover) {
    glp_set_col_bnds(p,(sg<0?j:V_pos[j])+1,GLP_FX,in_cover?1:0,0);
}

void graphe::mvc_solver::preprocess(glp_tree *tree) {
    glp_prob *p=glp_ios_get_prob(tree);
    int cur=glp_ios_curr_node(tree),parent;
    if (cur>0 && (parent=glp_ios_up_node(tree,cur))>0) {
        int ms=*static_cast<int*>(glp_ios_node_data(tree,parent));
        int bv=mirrors[ms],pos=ms+1;
        assert(glp_get_col_type(p,bv+1)==GLP_FX);
        /* mirror branching */
        if (glp_get_col_ub(p,bv+1)==0) {
            const ivector &ngh=G->node(sg<0?bv:V[bv]).neighbors();
            for (ivector_iter it=ngh.begin();it!=ngh.end();++it) {
                if (sg<0 || G->node(*it).subgraph()==sg)
                    make_vertex_fixed(p,*it,true);
            }
        } else {
            bool in_cover;
            int cnt=0;
            while (mirrors[pos]>=0) {
                if (!is_vertex_fixed(p,mirrors[pos],in_cover)) {
                    make_vertex_fixed(p,mirrors[pos],true);
                    ++cnt;
                }
                ++pos;
            }
        }
    }
    /* get already fixed vertices */
    int n=sg<0?G->node_count():V.size(),cov_lbound=0;
    if (sg<0)
        G->unset_subgraphs(1);
    bool in_cover;
    for (int i=0;i<n;++i) {
        int j=sg<0?i:V[i];
        if (is_vertex_fixed(p,j,in_cover)) {
            G->node(j).set_subgraph(0);
            if (in_cover)
                ++cov_lbound;
        }
    }
    bool changed;
    int s=sg<0?1:sg;
    /* try to reduce the uncovered part of G */
    do {
        changed=false;
        for (int k=0;k<n;++k) {
            int i=sg<0?k:V[k];
            vertex &v=G->node(i);
            if (v.subgraph()==s) {
                if (G->degree(i,s)==1) {
                    make_vertex_fixed(p,i,false);
                    int j=G->first_neighbor_from_subgraph(G->node(i),s);
                    make_vertex_fixed(p,j,true);
                    ++cov_lbound;
                    G->node(i).set_subgraph(0);
                    G->node(j).set_subgraph(0);
                    changed=true;
                    break;
                } else if (G->mvc_is_dominant(i,s) || G->mvc_is_unconfined(i,s)) {
                    make_vertex_fixed(p,i,true);
                    ++cov_lbound;
                    G->node(i).set_subgraph(0);
                    changed=true;
                    break;
                }
            }
        }
    } while (changed);
    cov_lbound+=lower_bound(s);
    if (cov_lbound>(int)glp_get_row_lb(p,last_row))
        glp_set_row_bnds(p,last_row,GLP_LO,cov_lbound,0);
    if (sg<0)
        G->unset_subgraphs();
    else G->set_subgraph(V,sg);
}

void graphe::mvc_solver::packing(glp_tree *tree) {
    glp_prob *p=glp_ios_get_prob(tree);
    int n=sg<0?G->node_count():V.size();
    bool in_cover;
    int *ind=new int[n+1];
    double *val=new double[n+1];
    int cnt=0;
    for (int i=0;i<n;++i) {
        int v=sg<0?i:V[i];
        if (is_vertex_fixed(p,v,in_cover) && in_cover) {
            ivector nv=G->node(v).neighbors();
            for (int j=nv.size();j-->0;) {
                if (sg>=0 && G->node(nv[j]).subgraph()!=sg)
                    nv.erase(nv.begin()+j);
            }
            double sum=0;
            for (ivector_iter it=nv.begin();it!=nv.end();++it) {
                sum+=glp_get_col_prim(p,(sg<0?*it:V_pos[*it])+1);
            }
            if (sum>nv.size()-0.999) {
                /* add packing constraint */
                int r=glp_add_rows(p,1);
                for (ivector_iter it=nv.begin();it!=nv.end();++it) {
                    ind[it-nv.begin()+1]=(sg<0?*it:V_pos[*it])+1;
                    val[it-nv.begin()+1]=1;
                }
                glp_set_mat_row(p,r,nv.size(),ind,val);
                glp_set_row_bnds(p,r,GLP_UP,0,nv.size()-1);
                ++cnt;
            }
        }
    }
    delete[] ind;
    delete[] val;
}

void graphe::mvc_solver::callback(glp_tree *tree,void *info) {
    if (is_interrupted()) {
        glp_ios_terminate(tree);
        return;
    }
    mvc_solver *mvcs=static_cast<mvc_solver*>(info);
    switch (glp_ios_reason(tree)) {
    case GLP_IBRANCH:
        /* choose a variable to branch upon */
        mvcs->branch(tree);
        break;
    case GLP_IHEUR:
        /* obtain a heuristic solution */
        if (mvcs->compute_heur) {
            mvcs->heuristic(tree);
            glp_ios_heur_sol(tree,mvcs->heur_sol);
            //mvcs->compute_heur=false;
        }
        break;
    case GLP_IROWGEN:
        /* generate packing constraints */
        mvcs->packing(tree);
        break;
    case GLP_IPREPRO:
        /* preprocess subproblem */
        if (!mvcs->is_k_vc)
            mvcs->preprocess(tree);
        break;
    default:
        /* ignore call for other reasons */
        break;
    }
}
#endif

/*
 * END OF MVC SOLVER CLASS IMPLEMENTATION
 */

/* handle connected component indexed by sg if it is a tree, cycle or clique,
 * return false otherwise */
bool graphe::mvc_special(ivector &cover,int sg) {
    assert(sg>0);
    if (is_empty(sg))
        return true;
    int mind=minimum_degree(sg),maxd=maximum_degree(sg);
    int ec=edge_count(sg),sz=subgraph_size(sg);
    ivector s;
    get_subgraph(sg,s);
    if (ec==sz-1) {
        /* tree component detected */
        while (ec>0) {
            for (ivector_iter jt=s.begin();jt!=s.end();++jt) {
                vertex &vl=node(*jt);
                if (vl.subgraph()==0)
                    continue;
                if (degree(*jt,sg)==1) { // leaf node found
                    int p=first_neighbor_from_subgraph(vl,sg);
                    cover.push_back(p);
                    vertex &vp=node(p);
                    ec-=degree(p,sg);
                    vl.set_subgraph(0);
                    vp.set_subgraph(0);
                    break;
                }
            }
        }
        set_subgraph(s,sg);
        return true;
    }
    if (mind==2 && maxd==2) {
        /* cycle component detected */
        int v=s.front();
        ivector V;
        dfs(v,true,true,&V,sg);
        for (int i=0;i<int(V.size());i+=2) {
            cover.push_back(V[i]);
        }
        return true;
    }
    if (mind==sz-1) {
        /* clique component detected */
        cover.insert(cover.end(),s.begin()+1,s.end());
        return true;
    }
    return false;
}

bool graphe::mvc_is_unconfined(int i,int sg) const {
    iset S,res,nu,nS;
    int u,m,dsize,w;
    assert(node(i).subgraph()==sg);
    S.insert(i);
    while (true) {
        nS=S;
        for (iset_iter it=S.begin();it!=S.end();++it) {
            const ivector &ngh=node(*it).neighbors();
            for (ivector_iter jt=ngh.begin();jt!=ngh.end();++jt) {
                if (node(*jt).subgraph()==sg)
                    nS.insert(*jt);
            }
        }
        for (iset_iter it=S.begin();it!=S.end();++it) {
            const ivector &ngh=node(*it).neighbors();
            u=-1;
            for (ivector_iter jt=ngh.begin();jt!=ngh.end();++jt) {
                if (node(*jt).subgraph()!=sg)
                    continue;
                assert(S.find(*jt)==S.end());
                const ivector &ngh2=node(*jt).neighbors();
                nu.clear();
                for (ivector_iter nt=ngh2.begin();nt!=ngh2.end();++nt) {
                    if (node(*nt).subgraph()==sg)
                        nu.insert(*nt);
                }
                if (sets_intersection(nu,S,res)==1) {
                    dsize=sets_difference(nu,nS,res);
                    if (u<0 || dsize<m) {
                        m=dsize;
                        u=*jt;
                        w=*res.begin();
                    }
                }
            }
            if (u<0)
                return false;
            if (m==0)
                return true;
            if (m==1) {
                S.insert(w);
                break;
            }
            return false;
        }
    }
}

bool graphe::mvc_is_dominant(int v,int sg) const {
    int u;
    const ivector &nv=node(v).neighbors();
    for (ivector_iter it=nv.begin();it!=nv.end();++it) {
        u=*it;
        if (sg>=0 && node(u).subgraph()!=sg)
            continue;
        const ivector &nu=node(u).neighbors();
        ivector_iter jt=nu.begin();
        for (;jt!=nu.end();++jt) {
            if (*jt==v || (sg>=0 && node(*jt).subgraph()!=sg))
                continue;
            ivector_iter kt=std::lower_bound(nv.begin(),nv.end(),*jt);
            if (kt==nv.end() || *kt!=*jt)
                break;
        }
        if (jt==nu.end())
            return true;
    }
    return false;
}

/* reduce the graph, set each vertex subgraph field as follows:
 * -2: vertex is in cover
 * -3: vertex is in cover
 * other: undecided */
bool graphe::mvc_reduce_basic(int sg,int c) {
    // sg>=0 is connected
    ivector V1,V2,cover;
    bool bp=false;
    int i;
    vector<vertex>::iterator it;
    if (mvc_special(cover,sg) || (bp=is_bipartite(V1,V2,sg,_GT_CC_CONNECTED))) {
        if (bp)
            mvc_bipartite(V1,V2,cover,sg,true);
        for (i=0,it=nodes.begin();it!=nodes.end();++it,++i) {
            if (it->subgraph()!=sg)
                continue;
            it->set_subgraph(std::find(cover.begin(),cover.end(),i)==cover.end()?c:c+1);
        }
        return false;
    }
    for (i=0,it=nodes.begin();it!=nodes.end();++it,++i) {
        if (it->subgraph()!=sg)
            continue;
        if (degree(i,sg)==1) {
            it->set_subgraph(c);
            node(first_neighbor_from_subgraph(*it,sg)).set_subgraph(c+1);
            break;
        }
        if (mvc_is_unconfined(i,sg)) {
            it->set_subgraph(c+1);
            break;
        }
    }
    return it!=nodes.end();
}

/* compute half-integral solution of MVC LP relaxation
 * and use it to reduce the graph */
void graphe::mvc_half_integral(int sg,ivector &in_cover,ivector &out_cover) {
    ivector S;
    get_subgraph(sg,S);
    graphe B(ctx,false);
    int n=S.size(),i,j;
    B.add_nodes(2*n);
    for (ivector_iter it=S.begin();it!=S.end();++it) {
        i=it-S.begin();
        for (ivector_iter jt=it+1;jt!=S.end();++jt) {
            j=jt-S.begin();
            if (has_edge(*it,*jt)) {
                B.add_edge(std::make_pair(i,n+j));
                B.add_edge(std::make_pair(j,n+i));
            }
        }
    }
    ivector cover,U(n),V(n);
    for (int i=0;i<n;++i) {
        U[i]=i;
        V[i]=n+i;
    }
    B.mvc_bipartite(U,V,cover);
    for (ivector_iter it=S.begin();it!=S.end();++it) {
        int lv=it-S.begin(),rv=n+lv;
        vertex &v=node(*it);
        if (std::find(cover.begin(),cover.end(),lv)==cover.end() &&
                std::find(cover.begin(),cover.end(),rv)==cover.end()) {
            out_cover.push_back(*it);
            const ivector &ngh=v.neighbors();
            for (ivector_iter jt=ngh.begin();jt!=ngh.end();++jt) {
                vertex &w=node(*jt);
                if (w.subgraph()==sg)
                    in_cover.push_back(*jt);
            }
        } else if (std::find(cover.begin(),cover.end(),lv)!=cover.end() &&
                   std::find(cover.begin(),cover.end(),rv)!=cover.end()) {
            in_cover.push_back(*it);
        }
    }
}

/* find a (minimum) vertex cover (of subgraph sg) with respect to
 * the algorithm specification vc_alg */
bool graphe::mvc(ivector &cover,int vc_alg,int sg,int tm_lim,double gap_tol,bool verbose) {
    cover.clear();
    if (is_null())
        return true;
    int c=-1,s,s0,cmp,i;
    bool changed=true;
    while (changed) {
        s0=c<0?(sg<0?1:1+max_subgraph_index()):c+3;
        cmp=connected_components_to_subgraphs(c<0?sg:c+2);
        if (c<0) c=cmp+1;
        for (s=s0;s<=cmp;++s) {
            while (mvc_reduce_basic(s,c));
        }
        changed=false;
        for (int s=s0;s<=cmp;++s) {
            ivector in_cover,out_cover;
            mvc_half_integral(s,in_cover,out_cover);
            for (ivector_iter it=in_cover.begin();it!=in_cover.end();++it) {
                node(*it).set_subgraph(c+1); changed=true;
            }
            for (ivector_iter it=out_cover.begin();it!=out_cover.end();++it) {
                node(*it).set_subgraph(c); changed=true;
            }
        }
        for (std::vector<vertex>::iterator it=nodes.begin();it!=nodes.end();++it) {
            if (it->subgraph()<c || it->subgraph()>c+1)
                it->set_subgraph(c+2);
        }
    }
    for (i=0;i<node_count();++i) {
        if (node(i).subgraph()==c+1)
            cover.push_back(i);
    }
    if (subgraph_size(c+2)>0) {
        cmp=connected_components_to_subgraphs(c+2);
        for (s=c+3;s<=cmp;++s) {
            if (mvc_special(cover,s))
                continue;
            if (vc_alg==_GT_VC_APPROX_CLIQUE) {
                ivector indp;
                grasp_clique(5,indp,true,s);
                node_iter it=nodes.begin(),itend=nodes.end();
                for (i=0;it!=itend;++it,++i) {
                    if (it->subgraph()==s && std::find(indp.begin(),indp.end(),it-nodes.begin())==indp.end())
                        cover.push_back(i);
                }
            } else if (vc_alg==_GT_VC_APPROX_ALOM) {
                mvc_alom(cover,s);
            } else if (vc_alg==_GT_VC_APPROX_DFS) {
                mvc_dfs(cover,s);
            } else if (vc_alg==_GT_VC_EXACT) {
                ivector V1,V2;
                if (is_bipartite(V1,V2,s,_GT_CC_CONNECTED))
                    mvc_bipartite(V1,V2,cover,s,true);
                else {
#ifdef HAVE_LIBGLPK
                    mvc_solver m(this,s);
                    if (m.solve(cover,-1,tm_lim,gap_tol,verbose)!=0)
                        return false;
#else
                    message("Error: GLPK library is required for exact solving");
                    return false;
#endif
                }
            }
        }
    }
    std::sort(cover.begin(),cover.end());
    return true;
}

/* DFS algorithm for vertex cover (fastest) */
void graphe::mvc_dfs(ivector &cover,int sg) {
    ivector V;
    if (sg>=0)
        get_subgraph(sg,V);
    int n=sg<0?node_count():int(V.size());
    int i=_rand(n,ctx).val;
    dfs(sg<0?i:V[i],false,true);
    for (int j=0;j<n;++j) {
        int v=sg<0?j:V[j];
        if (j!=i && !node(v).is_leaf())
            cover.push_back(v);
    }
}

/* return vertex candidates as suggested by M. Alom */
void graphe::alom_candidates(const ivector &V,ivector &cand,int sg) {
    assert(!V.empty());
    ivector maxdv;
    maxdv.reserve(V.size());
    int maxd=0;
    for (ivector_iter it=V.begin();it!=V.end();++it) {
        int d=degree(*it,sg);
        if (d>maxd) {
            maxdv=ivector(1,*it);
            maxd=d;
        } else if (d==maxd)
            maxdv.push_back(*it);
    }
    if (maxdv.size()==1) {
        cand=maxdv;
        return;
    }
    ivector inc(maxdv.size(),0);
    for (ivector_iter it=maxdv.begin();it!=maxdv.end();++it) {
        for (ivector_iter jt=it+1;jt!=maxdv.end();++jt) {
            if (has_edge(*it,*jt)) {
                inc[it-maxdv.begin()]++;
                inc[jt-maxdv.begin()]++;
            }
        }
    }
    int max_inc=0;
    cand.clear();
    for (ivector_iter it=inc.begin();it!=inc.end();++it) {
        if (*it>max_inc) {
            cand=ivector(1,maxdv[it-inc.begin()]);
            max_inc=*it;
        } else if (*it==max_inc)
            cand.push_back(maxdv[it-inc.begin()]);
    }
}

/* Alom's algorithm for vertex cover (destructive) */
void graphe::mvc_alom(ivector &cover,int sg) {
    assert(sg>=0);
    ivector V,cand;
    get_subgraph(sg,V);
    int ne=edge_count(sg);
    cand.reserve(V.size());
    while (ne>0) {
        alom_candidates(V,cand,sg);
        assert(!cand.empty());
        int c=cand[_rand(cand.size(),ctx).val];
        vertex &v=node(c);
        ipairs E;
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            if (node(*it).subgraph()==sg)
                E.push_back(make_pair(c,*it));
        }
        for (ipairs_iter it=E.begin();it!=E.end();++it) {
            assert(remove_edge(*it));
            --ne;
        }
        cover.push_back(c);
    }
}

/* find minimum vertex cover in bipartite (sub)graph,
 * (U and V must be sorted) */
void graphe::mvc_bipartite(const ivector &U,const ivector &V,ivector &cover,int sg,bool conn) {
    ipairs matching;
    ivector s;
    /* find a maximal matching */
    if (sg>=0) get_subgraph(sg,s);
    bipartite_matching(U,V,matching,sg,conn?_GT_CC_CONNECTED:_GT_CC_FIND_COMPONENTS);
    int n=sg<0?node_count():int(s.size());
    bvector matched(node_count(),false);
    for (ipairs_iter it=matching.begin();it!=matching.end();++it) {
        matched[it->first]=matched[it->second]=true;
    }
    /* create the alternating forest */
    for (int i=0;i<n;++i) {
        node(sg<0?i:s[i]).set_number(0);
    }
    ipairs E;
    get_edges_as_pairs(E,sg);
    ivector match(node_count(),-1);
    for (ipairs_iter it=matching.begin();it!=matching.end();++it) {
        match[it->first]=it->second;
        match[it->second]=it->first;
    }
    for (ivector_iter it=U.begin();it!=U.end();++it) {
        if (!matched[*it])
            node(*it).set_number(2);
    }
    while (true) {
        int u,v,w;
        ipairs::iterator it=E.begin();
        for (;it!=E.end();++it) {
            if (node(it->first).number()==2) {
                u=it->first; v=it->second;
                break;
            } else if (node(it->second).number()==2) {
                u=it->second; v=it->first;
                break;
            }
        }
        if (it==E.end())
            break;
        E.erase(it);
        if (matched[v] && node(v).number()==0) {
            assert((w=match[v])>=0);
            it=std::find(E.begin(),E.end(),std::make_pair(v,w));
            if (it==E.end())
                it=std::find(E.begin(),E.end(),std::make_pair(w,v));
            assert(it!=E.end());
            E.erase(it);
            node(v).set_number(1);
            node(w).set_number(2);
        } else if (node(v).number()==1)
            continue;
        else assert(false);
    }
    /* apply Kőnig's theorem */
    ivector Z,D(U.size()),I(V.size());
    for (int i=0;i<n;++i) {
        int k=sg<0?i:s[i];
        if (node(k).number()>0)
            Z.push_back(k);
    }
    std::sort(Z.begin(),Z.end());
    ivector_iter it=std::set_difference(U.begin(),U.end(),Z.begin(),Z.end(),D.begin());
    D.resize(it-D.begin());
    it=std::set_intersection(V.begin(),V.end(),Z.begin(),Z.end(),I.begin());
    I.resize(it-I.begin());
    size_t oldsize=cover.size();
    cover.resize(oldsize+D.size()+I.size());
    it=std::set_union(D.begin(),D.end(),I.begin(),I.end(),cover.begin()+oldsize);
    cover.resize(it-cover.begin());
}

/* find vertex cover of size k, return true on success,
 * if such cover does not exist return false */
int graphe::k_vertex_cover(ivector &cover,int k) {
    int n=node_count();
    if (k<0 || k>n)
        return 1;
    if (k==n)
        return 0;
    if (k==0)
        return is_empty()?0:1;
    /* find initial vertex cover as a complement of a large independent set */
    graphe C(ctx);
    complement(C);
    ivector indp,V(n);
    C.grasp_clique(5,indp);
    for (int i=0;i<n;++i) V[i]=i;
    for (ivector_iter it=indp.begin();it!=indp.end();++it) {
        V.erase(std::find(V.begin(),V.end(),node_index(C.node_label(*it))));
    }
    if (int(V.size())==k) {
        cover=V;
        return 0;
    }
    if (int(V.size())<k) { // add vertices to the obtained cover at random
        cover=V;
        vecteur W(indp.size());
        for (ivector_iter it=indp.begin();it!=indp.end();++it) {
            W[it-indp.begin()]=*it;
        }
        vecteur rnd=*_rand(makesequence(k-V.size(),W),ctx)._VECTptr;
        for (const_iterateur it=rnd.begin();it!=rnd.end();++it) {
            cover.push_back(it->val);
        }
        assert(k==int(cover.size()));
        std::sort(cover.begin(),cover.end());
        return 0;
    }
    /* size of V is larger than k, use exact solver */
#ifdef HAVE_LIBGLPK
    mvc_solver mvcs(this);
    return mvcs.solve(cover,k);
#else
    message("Error: GLPK library is required for exact solving");
    return -1;
#endif
}

/* return vertex cover number of this graph */
int graphe::vertex_cover_number(int sg) {
    ivector V1,V2,V;
    int res=0,s=(sg<0)?1:1+max_subgraph_index();
    int cmp=connected_components_to_subgraphs(sg);
    for (;s<=cmp;++s) {
        get_subgraph(s,V);
        ivector cov;
        if (mvc_special(cov,s)) {
            res+=cov.size();
        } else if (is_bipartite(V1,V2,s,_GT_CC_CONNECTED)) {
            /* apply Kőnig's theorem */
            ipairs m;
            bipartite_matching(V1,V2,m,s,_GT_CC_CONNECTED);
            res+=m.size();
        } else {
            if (!mvc(cov,_GT_VC_EXACT,s))
                return -1;
            res+=cov.size();
        }
        if (sg>=0) renumber_subgraph(s,sg);
    }
    return res;
}

/* return true iff v is reachable from u */
bool graphe::is_reachable(int u,int v) {
    int n=node_count(),w;
    assert(u>=0 && u<n && v>=0 && v<n);
    std::queue<int> q;
    q.push(u);
    unvisit_all_nodes();
    while (!q.empty()) {
        w=q.front();
        q.pop();
        const ivector &ngh=node(w).neighbors();
        for (ivector_iter it=ngh.begin();it!=ngh.end();++it) {
            if (node(*it).is_visited())
                continue;
            if (*it==v)
                return true;
            q.push(*it);
            node(*it).set_visited(true);
        }
    }
    return false;
}

/* fill r with vertices reachable from u */
void graphe::reachable(int u,ivector &r) {
    int n=node_count(),w;
    assert(u>=0 && u<n);
    std::queue<int> q;
    q.push(u);
    unvisit_all_nodes();
    r.clear();
    r.reserve(n-1);
    while (!q.empty()) {
        w=q.front();
        q.pop();
        const ivector &ngh=node(w).neighbors();
        for (ivector_iter it=ngh.begin();it!=ngh.end();++it) {
            if (node(*it).is_visited())
                continue;
            r.push_back(*it);
            q.push(*it);
            node(*it).set_visited(true);
        }
    }
}

/* return true iff the i-th vertex is simplicial and mark
 * high-degree neighbors (deg>D) as visited */
bool graphe::is_simplicial(int i,const sparsemat &A,double D) {
    gen aii=sparse_product_element(A,A,i,i);
    bool yes=true;
    const vertex &v=node(i);
    for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
        if (yes && !giac::is_zero(aii-sparse_product_element(A,A,i,*it)))
            yes=false;
        if (D>0.5 && node(*it).degree()>D)
            node(*it).set_visited(true);
    }
    return yes;
}

/* find all simplicial vertices and store them in res,
 * using the algorithm of Kloks et al. (worst case complexity |E|^3/2) */
void graphe::find_simplicial_vertices(ivector &res) {
    res.clear();    
    res.reserve(node_count());
    int n=node_count(),m=edge_count();
    double alpha=1+std::log(m)/std::log(n);
    double D=std::pow(m,(alpha-1)*(alpha+1));
    sparsemat A;
    adjacency_sparse_matrix(A,true);
    unset_subgraphs(0);
    for (int i=0;i<n;++i) {
        vertex &v=node(i);
        if (v.degree()<=D) { // low degree vertex
            if (is_simplicial(i,A,D))
                res.push_back(i);
        } else v.set_subgraph(1);
    }
    adjacency_sparse_matrix(A,true,1);
    ivector s;
    get_subgraph(1,s);
    for (ivector_iter it=s.begin();it!=s.end();++it) {
        if (!node(*it).is_visited() && is_simplicial(*it,A))
            res.push_back(*it);
    }
}

/* SIP CLASS IMPLEMENTATION */

/* constructor */
graphe::sip::sip(const graphe &G,const graphe &P,int max_sg_in) {
    N=G.node_count();
    n=P.node_count();
    isom=(n==N && G.edge_count()==P.edge_count());
    dir=G.is_directed();
    assert(P.is_directed()==dir);
    max_sg=max_sg_in;
    /* initialize adjacency matrices (A is transposed) */
    matrice PA;
    P.adjacency_matrix(PA);
    B.resize(n);
    for (int i=0;i<n;++i) {
        B[i].resize(n);
        for (int j=0;j<n;++j) {
            B[i][j]=giac::is_one(PA[i][j]);
        }
    }
    used_cols.resize(N,false);
    A.resize(N); AT.resize(N);
    nb=N/SIP_NBITS;
    if (N%SIP_NBITS) ++nb;
    for (int i=0;i<N;++i) {
        A[i].resize(nb);
        AT[i].resize(nb);
    }
    for (int i=0;i<N;++i) {
        const vertex &v=G.node(i);
        for (ivector_iter it=v.neighbors().begin();it!=v.neighbors().end();++it) {
            set_ij(A,*it,i,true);
            set_ij(AT,i,*it,true);
        }
    }
    /* initialize M */
    M.resize(n);
    bool cond;
    for (int i=0;i<n;++i) {
        M[i].resize(nb);
        for (int j=0;j<N;++j) {
            if (isom) {
                if (dir)
                    cond=P.in_degree(i)==G.in_degree(j) && P.out_degree(i)==G.out_degree(j);
                else cond=P.degree(i)==G.degree(j);
            } else {
                if (dir)
                    cond=P.in_degree(i)<=G.in_degree(j) && P.out_degree(i)<=G.out_degree(j);
                else cond=P.degree(i)<=G.degree(j);
            }
            set_ij(M,i,j,cond && P.get_node_color(i)==G.get_node_color(j));
        }
    }
}

/* scalar product of r and c, it is assumed that r has exactly one 1 */
bool graphe::sip::mult_bitrows(const bitrow &r,const bitrow &c) {
    for (int i=0;i<N;++i) {
        if (r[i].any())
            return (r[i] & c[i]).any();
    }
    assert(false);
}

/* fill row with zeros */
void graphe::sip::clear_bitrow(bitrow &row) {
    for (bitrow::iterator it=row.begin();it!=row.end();++it) {
        it->reset();
    }
}

/* check if M defines an isomorphism */
bool graphe::sip::is_isomorphism() {
    bitmatrix prod(n);
    for (int i=0;i<n;++i) {
        prod[i].resize(nb);
        const bitrow &r=M[i];
        for (int j=0;j<N;++j) {
            if (mult_bitrows(r,A[j]))
                set_ij(prod,i,j,true);
        }
    }
    for (int i=0;i<n;++i) {
        const bitrow &r=M[i];
        for (int j=0;j<n;++j) {
            bool p=mult_bitrows(r,prod[j]),q=B[j][i];
            if ((induced && p!=q) || (!induced && q && !p))
                return false;
        }
    }
    return true;
}

/* prune matrix M */
void graphe::sip::prune() {
    bool changed=true;
    while (changed) {
        changed=false;
        for (int i=0;i<n;++i) {
            for (int j=0;j<N;++j) {
                if (get_ij(M,i,j)) {
                    for (int c=0;c<n;++c) {
                        const bitrow &r=AT[j],&s=A[j],&m=M[c];
                        int k=0,sz=r.size();
                        if (B[i][c]) {
                            for (k=0;k<sz;++k) {
                                if ((r[k] & m[k]).any())
                                    break;
                            }
                        }
                        if (dir && k<sz && B[c][i]) {
                            for (k=0;k<sz;++k) {
                                if ((s[k] & m[k]).any())
                                    break;
                            }
                        }
                        if (k==sz) {
                            edit(i,j);
                            changed=true;
                            break;
                        }
                    }
                }
            }
        }
    }
}

/* undo the last set of changes to M */
void graphe::sip::revert_changes() {
    const ipairs &ch=snapshots.top();
    for (ipairs_iter it=ch.begin();it!=ch.end();++it) {
        set_ij(M,it->first,it->second,true);
    }
    snapshots.pop();
}

/* return true iff there is a zero-row in M */
bool graphe::sip::M_has_empty_row() {
    int k=0;
    for (;k<n;++k) {
        const bitrow &mr=M[k];
        bitrow::const_iterator it=mr.begin();
        for (;it!=mr.end();++it) {
            if (it->any())
                break;
        }
        if (it==mr.end())
            break;
    }
    return k<n;
}

/* main recursion procedure */
bool graphe::sip::recurse(int i) {
    if (i==n) {
        if (is_isomorphism()) {
            ivector sg(n);
            for (int k=0;k<n;++k) {
                int l=0;
                for (;l<N;++l) {
                    if (get_ij(M,k,l))
                        break;
                }
                assert(l<N);
                sg[k]=l;
            }
            found.push_back(sg);
        }
    } else {
        prune();
        push_changes();
        if (!M_has_empty_row()) {
            bitrow row=M[i];
            clear_bitrow(M[i]);
            int res;
            for (int j=0;j<N;++j) {
                if (!get_j(row,j) || used_cols[j])
                    continue;
                set_ij(M,i,j,true);
                /* filter */
                for (bvector::const_iterator it=B[i].begin();it!=B[i].end();++it) {
                    if (*it) {
                        for (int a=0;a<N;++a) {
                            if (!get_ij(AT,j,a) && get_ij(M,it-B[i].begin(),a))
                                edit(it-B[i].begin(),a);
                        }
                    }
                }
                if (induced) {
                    for (int a=0;a<N;++a) {
                        if (get_ij(AT,j,a)) {
                            for (bvector::const_iterator it=B[i].begin();it!=B[i].end();++it) {
                                if (!*it && get_ij(M,it-B[i].begin(),a))
                                    edit(it-B[i].begin(),a);
                            }
                        }
                    }
                }
                push_changes();
                /* recurse */
                used_cols[j]=true;
                res=recurse(i+1);
                /* restore */
                revert_changes();
                used_cols[j]=false;
                set_ij(M,i,j,false);
                if (max_sg>0 && res==max_sg)
                    break;
            }
            M[i]=row;
        }
        revert_changes(); // revert pruning
    }
    return found.size();
}

int graphe::sip::find_subgraphs(bool only_induced) {
    found.clear();
    while (!snapshots.empty()) snapshots.pop();
    induced=only_induced;
    if (M_has_empty_row())
        return 0;
    int res=recurse();
    assert(snapshots.empty());
    return res;
}

/* END OF SIP CLASS */

void graphe::subgraph_isomorphism(graphe &P,int max_sg,bool induced,ivectors &res) const {
    int m=P.node_count();
    graphe Q(ctx,P.supports_attributes());
    ivector sigma;
    P.sort_by_degrees(Q,sigma);
    sip finder(*this,Q,max_sg);
    int n=finder.find_subgraphs(induced);
    res.resize(n);
    ivector sg;
    for (int i=0;i<n;++i) {
        sg=finder.get_subgraph(i);
        res[i].resize(m);
        for (int j=0;j<m;++j) {
            res[i][sigma[j]]=sg[j];
        }
    }
}

/* Encode this graph as a Pruefer sequence <code>, return true iff it is a tree.
 * Linear time implementation (https://cp-algorithms.com/graph/pruefer_code.html) */
bool graphe::pruefer_encode(ivector &code,bool check_tree) {
    if (is_directed())
        return false;
    iset v;
    int n=node_count();
    if (n<2) return false; // not a tree
    code.clear();
    if (n==2 && edge_count()==1)
        return true; // sequence is empty
    if (check_tree && !is_tree())
        return false; // not a tree
    code.resize(n-2);
    dfs(n-1,false);
    int ptr=-1;
    ivector d(n);
    for (int i=0;i<n;++i) {
        d[i]=node(i).degree();
        if (d[i]==1 && ptr==-1)
            ptr=i;
    }
    int leaf=ptr,next;
    for (int i=0;i<n-2;++i) {
        next=node(leaf).ancestor();
        code[i]=next;
        if (--d[next]==1 && next<ptr)
            leaf=next;
        else {
            do ptr++; while (d[ptr]!=1);
            leaf=ptr;
        }
    }
    return true;
}

/* Create a tree from Pruefer code.
 * Linear time implementation (https://cp-algorithms.com/graph/pruefer_code.html) */
bool graphe::pruefer_decode(const ivector &code) {
    int n=code.size();
    this->clear();
    if (supports_attributes()) {
        vecteur labels;
        make_default_labels(labels,n+2);
        add_nodes(labels);
    } else add_nodes(n+2);
    ivector d(n+2,1);
    for (ivector_iter it=code.begin();it!=code.end();++it) {
        if (*it<0 || *it>=n+2)
            return false; // not a Pruefer code
        d[*it]++;
    }
    int ptr=0;
    while (d[ptr]!=1) ++ptr;
    int leaf=ptr;
    for (ivector_iter it=code.begin();it!=code.end();++it) {
        add_edge(leaf,*it);
        if (--d[*it]==1 && *it<ptr)
            leaf=*it;
        else {
            do ++ptr; while(d[ptr]!=1);
            leaf=ptr;
        }
    }
    add_edge(leaf,n+1);
    return true;
}

/* return the Giac name of this graph,
 * if no special graph is isomorphic to it then return NULL */ 
const char* graphe::identify() {
    if (is_directed() || !is_connected())
        return NULL;
    int i=0,n=node_count(),m=edge_count();
    map<int,int> isom;
    while (true) {
        const spcgraph &g=special_graph[i++];
        if (g.nv==0)
            break;
        if (g.nv==n && g.ne==m) {
            graphe G(g.name,ctx,false);
            if (is_isomorphic(G,isom,false)) {
                return g.name;
            }
        }
    }
    return NULL;
}

/* Find all special graphs known to Giac which are isomorphic to this graph.
 * Return true iff an isomorphic special graph is found.
 * Parameters of each found isomorphic graph are stored in spec as follows:
 * [SPECIAL_GRAPH_ID, param1, param2, ...] */
void graphe::identify_from_sequences(vecteur &spec,int haar_limit) {
    spec.clear();
    int nv=node_count(),ne=edge_count();
    if (nv<=1)
        return;
    graphe G(ctx,false);
    map<int,int> isom;
    if (is_directed()) {
        /* check Paley graph with p^k vertices where p is congruent to 3 modulo 4 */
        if (ne==nv*(nv-1)/2) {
            gen f=_ifactors(nv,ctx);
            if (f.type==_VECT && f._VECTptr->size()==2 && is_integer_vecteur(*f._VECTptr)) {
                int p=f._VECTptr->front().val,k=f._VECTptr->back().val;
                if ((p-3)%4==0) {
                    G.make_paley_graph(p,k);
                    if (is_isomorphic(G,isom,false))
                        spec.push_back(makevecteur(at_paley_graph,p,k));
                }
            }
        }
    } else { // the graph is undirected
        ivector U,V;
        bool conn=is_connected(),bp=is_bipartite(U,V);
        int maxdeg=maximum_degree(),mindeg=minimum_degree(),n1=U.size(),n2=V.size();
        bool reg=(maxdeg==mindeg);
        graphe GC(ctx,false);
        complement(GC);
        ivectors gc_components;
        GC.connected_components(gc_components);
        /* cycle graphs */
        if (conn && reg && nv>2 && nv==ne && maxdeg==2)
            spec.push_back(makevecteur(at_cycle_graph,nv));
        /* path graphs */
        if (conn && nv>2 && ne==nv-1 && maxdeg==2)
            spec.push_back(makevecteur(at_path_graph,nv));
        /* complete graphs */
        if (conn && nv>2 && 2*ne==nv*(nv-1))
            spec.push_back(makevecteur(at_complete_graph,nv));
        /* complete bipartite graphs */
        if (conn && nv>2 && bp && ne==n1*n2)
            spec.push_back(makevecteur(at_complete_graph,n1,n2));
        /* complete multipartite graphs */
        if (conn && gc_components.size()>2 && int(gc_components.size())<nv) {
            GC.unset_subgraphs();
            int sg=-1;
            ivectors_iter it=gc_components.begin();
            for (;it!=gc_components.end();++it) {
                int nvc=it->size();
                GC.set_subgraph(*it,++sg);
                if (2*GC.edge_count(sg)!=nvc*(nvc-1))
                    break;
            }
            if (it==gc_components.end()) {
                vecteur res;
                for (it=gc_components.begin();it!=gc_components.end();++it) {
                    res.push_back(it->size());
                }
                res=*_SortD(res,ctx)._VECTptr;
                res.insert(res.begin(),at_complete_graph);
                spec.push_back(res);
            }
        }
        /* complete trees */
        if (conn && is_tree()) {
            map<int,int> deg_count;
            node_iter nt=nodes.begin(),ntend=nodes.end();
            for (;nt!=ntend;++nt)
                deg_count[nt->degree()]++;
            if (deg_count.size()>1 && deg_count.size()<=3 && deg_count.begin()->first==1) {
                ivector d(deg_count.size()),cnt(deg_count.size());
                int i=0,k,n;
                for (map<int,int>::const_iterator it=deg_count.begin();it!=deg_count.end();++it) {
                    d[i]=it->first;
                    cnt[i++]=it->second;
                }
                if (cnt[1]==1 && deg_count.size()==2) {
                    spec.push_back(d[1]==2?makevecteur(at_complete_binary_tree,1):
                                           makevecteur(at_complete_kary_tree,d[1],1));
                } else if (cnt[1]==1 && (k=d[1])>1 && d[2]==k+1 && _logb(makesequence(cnt[0],k),ctx).is_integer()) {
                    n=_logb(makesequence(cnt[0],k),ctx).val;
                    if (tree_height(tree_height())==n && cnt[2]==((int)std::pow(k,n)-1)/(k-1)-1)
                        spec.push_back(k==2?makevecteur(at_complete_binary_tree,n):
                                            makevecteur(at_complete_kary_tree,k,n));
                }
            }
        }
        /* Kneser, odd and Johnson graphs */
        if (conn && reg && nv>2) {
            ipairs nk_kneser,nk_johnson;
            for (int n=3;n<=20;++n) {
                for (int k=1;k<=n/2;++k) {
                    int nck=_comb(makesequence(n,k),ctx).val;
                    if (nck==nv) {
                        if (n%2 || k<n/2) {
                            int nmkck=_comb(makesequence(n-k,k),ctx).val;
                            if (2*ne==nck*nmkck && maxdeg==nmkck) // Kneser graphs
                                nk_kneser.push_back(make_pair(n,k));
                        }
                        if (2*ne==k*(n-k)*nck && maxdeg==k*(n-k)) // Johnson graphs
                            nk_johnson.push_back(make_pair(n,k));
                    }
                }
            }
            for (ipairs_iter it=nk_kneser.begin();it!=nk_kneser.end();++it) {
                int n=it->first,k=it->second;
                G.make_intersection_graph(n,k,0);
                assert(nv==G.node_count() && ne==G.edge_count());
                if (is_isomorphic(G,isom,false)) {
                    spec.push_back(makevecteur(at_kneser_graph,n,k));
                    if (n==2*k+1)
                        spec.push_back(makevecteur(at_odd_graph,k+1));
                }
            }
            for (ipairs_iter it=nk_johnson.begin();it!=nk_johnson.end();++it) {
                int n=it->first,k=it->second;
                G.make_intersection_graph(n,k,k-1);
                assert(nv==G.node_count() && ne==G.edge_count());
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_johnson_graph,n,k));
            }
        }
        /* hypercube graphs */
        if (conn && reg && bp && (2*ne)%nv==0) {
            int n=(2*ne)/nv;
            if (nv==(1<<n) && ne==n*(1<<(n-1))) {
                G.make_hypercube_graph(n);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_hypercube_graph,n));
            }
        }
        /* star graphs */
        if (conn && maxdeg==ne && ne==nv-1 && ne>2)
            spec.push_back(makevecteur(at_star_graph,ne));
        /* web graphs */
        if (conn && mindeg==3 && maxdeg<=4 && nv%(2*nv-ne)==0) {
            int m=2*nv-ne,n=nv/m;
            if (m>2 && n>1) {
                G.make_web_graph(m,n);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_web_graph,m,n));
            }
        }
        /* wheel graphs */
        if (conn && nv>3 && maxdeg==nv-1 && mindeg==3 && ne==2*(nv-1))
            spec.push_back(makevecteur(at_wheel_graph,nv-1));
        /* prism graphs */
        if (conn && reg && maxdeg==3 && nv>5 && nv%2==0 && 2*ne==3*nv) {
            G.make_petersen_graph(nv/2,1);
            if (is_isomorphic(G,isom,false))
                spec.push_back(makevecteur(at_prism_graph,nv/2));
        }
        /* antiprism graphs */
        if (conn && reg && maxdeg==4 && nv%2==0 && ne==2*nv) {
            G.make_antiprism_graph(nv/2);
            if (is_isomorphic(G,isom,false))
                spec.push_back(makevecteur(at_antiprism_graph,nv/2));
        }
        /* grid graphs */
        if (conn && mindeg==2 && maxdeg<=4) { // rectangular grid
            int D=ne*ne-4*nv*(ne-nv+1);
            gen n=(2*nv-ne-giac::sqrt(D,ctx))/2;
            gen m=(2*nv-ne+giac::sqrt(D,ctx))/2;
            if (n.is_integer() && m.is_integer() && n.val>1 && m.val>1) {
                G.make_grid_graph(m.val,n.val,0);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_grid_graph,m,n));
            }
        }
        if (conn && mindeg==2 && maxdeg<=6) { // triangular grid
            int D=9*nv*nv+ne*ne-6*nv*ne-10*nv-2*ne+1;
            gen n=(3*nv-ne+1-giac::sqrt(D,ctx))/4;
            gen m=(3*nv-ne+1+giac::sqrt(D,ctx))/4;
            if (n.is_integer() && m.is_integer() && n.val>1 && m.val>1) {
                G.make_grid_graph(m.val,n.val,1);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_grid_graph,m,n,at_triangle));
            }
        }
        if (conn && reg && maxdeg==4 && ne==2*nv) { // toroidal grid
            vecteur d=*_idivis(nv,ctx)._VECTptr;
            for (const_iterateur it=d.begin();it!=d.end();++it) {
                gen m=*it,n=gen(nv)/(*it);
                if (is_greater(m,n,ctx)) {
                    G.make_grid_graph(m.val,n.val,2);
                    if (is_isomorphic(G,isom,false)) {
                        spec.push_back(makevecteur(at_torus_grid_graph,m,n));
                        break;
                    }
                }
            }
        }
        /* Sierpinski graphs */
        if (conn) {
            int n,p;
            gen lgb;
            if ((2*ne)%(nv-1)==0) {
                p=(2*ne)/(nv-1);
                if (p>2 && (lgb=_logb(makesequence(nv,p),ctx)).is_integer() && (n=lgb.val)>=1 &&
                        maxdeg==(n==1?p-1:p) && mindeg==p-1) {
                    G.make_sierpinski_graph(n,p,false);
                    if (is_isomorphic(G,isom,false))
                        spec.push_back(makevecteur(at_sierpinski_graph,n,p));
                }
            }
            int D=4*nv*(nv-1)-8*ne+1;
            gen pg=nv-(giac::sqrt(D,ctx)-1)/2;
            if (pg.is_integer() && (p=pg.val)>2 && (2*ne)%(p-1)==0 &&
                    (lgb=_logb(makesequence((2*ne)/(p-1),p),ctx)).is_integer() && (n=lgb.val)>=1 &&
                    maxdeg==(n==1?p-1:2*(p-1)) && mindeg==p-1) {
                G.make_sierpinski_graph(n,p,true);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_sierpinski_graph,n,p,at_triangle));
            }
        }
        /* generalized Petersen graphs */
        if (conn && reg && maxdeg==3 && nv%2==0 && 2*ne==3*nv && nv/2>2) {
            for (int k=1;k<=(nv/2-1)/2;++k) {
                G.make_petersen_graph(nv/2,k);
                if (is_isomorphic(G,isom,false)) {
                    spec.push_back(makevecteur(at_petersen_graph,nv/2,k));
                    break;
                }
            }
        }
        /* flower and Goldberg snarks */
        if (conn && reg && maxdeg==3 && nv%4==0 && 2*ne==3*nv) {
            int n;
            if ((n=nv/4)>2 && n%2==1) {
                G.make_flower_snark(n);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_flower_snark,n));
            } else if (n%2==0 && (n/2)%2==1) {
                G.make_goldberg_snark(n/2);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_goldberg_snark,n/2));
            }
        }
        /* Paley graphs */
        if (conn && reg && (nv-1)%4==0 && 4*ne==nv*(nv-1) && maxdeg==(nv-1)/2) {
            ipair sig;
            vecteur f=*_ifactors(nv,ctx)._VECTptr;
            if (f.size()==2 && f.front().val>2 && f.back().val>=1 &&
                    is_strongly_regular(sig) && sig.first==(nv-5)/4 && sig.second==(nv-1)/4) {
                G.make_paley_graph(f.front().val,f.back().val);
                if (is_isomorphic(G,isom,false))
                    spec.push_back(makevecteur(at_paley_graph,f.front(),f.back()));
            }
        }
        /* Haar graphs */
        if (reg && nv%2==0 && bp && n1==n2) {
            int k=nv/2,m;
            if (k<=std::min(haar_limit,_GT_HAAR_LIMIT)) {
                ulong lb=(1UL<<(k-1)),ub=(1UL<<k);
                ivector d(2*k),d0=degree_sequence();
                std::sort(d0.begin(),d0.end());
                for (ulong n=lb;n<ub;++n) {
                    bitset<64> be(n);
                    if (maxdeg!=int(be.count()))
                        continue;
                    std::fill(d.begin(),d.end(),0);
                    m=0;
                    for (int i=0;i<k;++i) {
                        for (int j=k;j<2*k;++j) {
                            if (be.test((j-i)%k)) {
                                d[i]++;
                                d[j]++;
                                m++;
                            }
                        }
                    }
                    if (m!=ne)
                        continue;
                    std::sort(d.begin(),d.end());
                    if (d==d0) {
                        G.make_haar_graph(n);
                        if (is_isomorphic(G,isom,false)) {
                            vecteur v(k);
                            for (int i=0;i<k;++i) {
                                v[i]=be.test(i)?1:0;
                            }
                            gen ng=_convert(makesequence(v,_BASE,2),ctx);
                            spec.push_back(makevecteur(at_haar_graph,ng));
                            break;
                        }
                    }
                }
            }
        }
    }
}

/* Return brief description of this graph as a string. */
string graphe::print(GIAC_CONTEXT) const {
    int nv=node_count(),ne=edge_count();
    stringstream ss;
    ss << nv;
    string nvert(ss.str());
    ss.str("");
    ss << ne;
    string nedg(ss.str());
    bool isdir=is_directed();
    string dir_spec=isdir?"directed":"undirected",weight_spec=is_weighted()?"weighted":"unweighted",disp_out;
    nvert=nvert+(nv==1?" vertex":" vertices");
    nedg=nedg+(isdir?(ne==1?" arc":" arcs"):(ne==1?" edge":" edges"));
    string name=this->name();
    if (!name.empty())
        disp_out=name+": ";
    disp_out+=gettext(((isdir?"a ":"an ")+dir_spec+" "+weight_spec+" graph with ").c_str())+nvert+gettext(" and ")+nedg;
    return disp_out;
}
string graphe::texprint(GIAC_CONTEXT) const {
    string ret="\\text{";
    ret+=print(contextptr);
    ret+="}";
    return ret;
}

/* Return the giac command (as string) which creates this graph.
 * This is required for saving sessions in Xcas. */
gen graphe::giac_constructor(GIAC_CONTEXT) const {
    vecteur V,E,A;
    V.reserve(node_count());
    E.reserve(edge_count());
    node_iter it=nodes.begin(),itend=nodes.end();
    ivector_iter jt,jtend;
    for (;it!=itend;++it) {
        gen_map vm;
        attrib2genmap(it->attributes(),vm,true);
        V.push_back(vm);
        jt=it->neighbors().begin(),jtend=it->neighbors().end();
        for (;jt!=jtend;++jt) {
            gen_map em;
            attrib2genmap(it->neighbor_attributes(*jt),em,true);
            E.push_back(makevecteur(makevecteur(it->label(),node(*jt).label()),em));
        }
    }
    vecteur args=makevecteur(V,change_subtype(E,_SET__VECT));
    for (attrib_iter kt=attributes.begin();kt!=attributes.end();++kt)
        args.push_back(symb_equal(string2gen(index2tag(kt->first),false),kt->second));
    return symbolic(at_graph,change_subtype(args,_SEQ__VECT));
}

/* Return true iff the othe graph is equal to this graph. */
bool graphe::operator ==(const gen &g) const {
    graphe *other=from_gen(g);
    if (other==NULL)
        return false;
    return this->is_equal(*other);
}

/* Return the pointer to the graph contained in G.
 * Return true if G contains a graph, otherwise return false. */
graphe *graphe::from_gen(const gen &g) {
    if (g.type!=_USER)
        return NULL;
    return dynamic_cast<graphe*>(g._USERptr);
}

/* SPECIAL GRAPHS
 *
 * Special graphs are initialized either with adjacency lists of integers or
 * strings, or by the corresponding LCF notation (suffix '_lcf').
 * When specified by adjacency lists, each list starts with the corresponding
 * vertex, followed by its neighbors, and ends with 0 resp. an empty string.
 * The graph specification terminates with -1 (integers) or NULL (strings).
 */
const int graphe::clebsch_graph[] = {
    1,2,3,4,5,6,0,2,7,8,9,10,0,3,7,11,12,13,0,4,8,11,14,15,0,5,9,12,14,16,
    0,6,10,13,15,16,0,7,14,15,16,0,8,12,13,16,0,9,11,13,15,0,10,11,12,14,
    0,11,16,0,12,15,0,13,14,0,-1
};
const int graphe::coxeter_graph[] = {
    1,2,3,4,0,2,23,26,0,3,24,28,0,4,25,27,0,5,8,11,19,0,6,9,13,22,0,7,10,
    14,21,0,8,15,18,0,9,16,17,0,10,12,20,0,11,12,23,0,12,25,0,13,14,23,0,
    14,24,0,15,16,24,0,16,25,0,17,20,26,0,18,21,26,0,19,22,28,0,20,28,0,
    21,27,0,22,27,0,-1
};
const char* graphe::coxeter_graph_vnames[] = {
    "a1","a2","a7","z1","","a2","a3","z2","","a3","a4","z3","","a4",
    "a5","z4","","a5","a6","z5","","a6","a7","z6","","a7","z7","",
    "b1","b3","b6","z1","","b2","b4","b7","z2","","b3","b5","z3","",
    "b4","b6","z4","","b5","b7","z5","","b6","z6","","b7","z7","",
    "c1","c4","c5","z1","","c2","c5","c6","z2","","c3","c6","c7","z3",
    "","c4","c7","z4","","c5","z5","","c6","z6","","c7","z7","",NULL
};
const int graphe::dyck_graph[] = {
    1,2,20,32,0,2,3,7,0,3,4,30,0,4,5,17,0,5,6,24,0,6,7,11,0,7,8,0,8,9,21,
    0,9,10,28,0,10,11,15,0,11,12,0,12,13,25,0,13,14,32,0,14,15,19,0,15,
    16,0,16,17,29,0,17,18,0,18,19,23,0,19,20,0,20,21,0,21,22,0,22,23,27,
    0,23,24,0,24,25,0,25,26,0,26,27,31,0,27,28,0,28,29,0,29,30,0,30,31,
    0,31,32,0,-1
};
const int graphe::grinberg_graph[] = {
    1,2,3,4,0,2,5,7,0,3,6,10,0,4,8,9,0,5,6,21,0,6,22,0,7,8,30,0,8,29,0,
    9,10,14,0,10,13,0,11,12,13,14,0,12,15,17,0,13,16,0,14,18,0,15,16,39,
    0,16,46,0,17,18,40,0,18,45,0,19,20,21,22,0,20,23,25,0,21,24,0,22,26,
    0,23,24,36,0,24,44,0,25,26,37,0,26,46,0,27,28,29,30,0,28,31,33,0,29,
    32,0,30,34,0,31,32,42,0,32,45,0,33,34,43,0,34,44,0,35,36,43,44,0,36,
    37,0,37,38,0,38,39,46,0,39,40,0,40,41,0,41,42,45,0,42,43,0,-1
};
const int graphe::grotzsch_graph[] = {
    1,2,5,7,10,0,2,3,6,8,0,3,4,7,9,0,4,5,8,10,0,5,6,9,0,6,11,0,7,11,0,8,
    11,0,9,11,0,10,11,0,-1
};
const int graphe::heawood_graph[] = {
    1,2,6,14,0,2,3,11,0,3,4,8,0,4,5,13,0,5,6,10,0,6,7,0,7,8,12,0,8,9,0,
    9,10,14,0,10,11,0,11,12,0,12,13,0,13,14,0,-1
};
const int graphe::herschel_graph[] = {
    1,2,4,5,0,2,3,6,7,0,3,4,8,0,4,9,10,0,5,6,10,0,6,11,0,7,8,11,0,8,9,0,
    9,11,0,10,11,0,-1
};
const int graphe::mcgee_graph[] = {
    1,2,13,24,0,2,3,9,0,3,4,20,0,4,5,16,0,5,6,12,0,6,7,23,0,7,8,19,0,8,9,
    15,0,9,10,0,10,11,22,0,11,12,18,0,12,13,0,13,14,0,14,15,21,0,15,16,0,
    16,17,0,17,18,24,0,18,19,0,19,20,0,20,21,0,21,22,0,22,23,0,23,24,0,-1
};
const int graphe::pappus_graph[] = {
    1,2,6,18,0,2,3,9,0,3,4,14,0,4,5,11,0,5,6,16,0,6,7,0,7,8,12,0,8,9,15,
    0,9,10,0,10,11,17,0,11,12,0,12,13,0,13,14,18,0,14,15,0,15,16,0,16,
    17,0,17,18,0,-1
};
const int graphe::robertson_graph[] = {
    1,2,5,9,19,0,2,3,7,13,0,3,4,10,15,0,4,5,8,18,0,5,6,12,0,6,7,14,17,0,7,
    8,11,0,8,9,16,0,9,10,14,0,10,11,17,0,11,12,19,0,12,13,16,0,13,14,18,
    0,14,15,0,15,16,19,0,16,17,0,17,18,0,18,19,0,-1
};
const int graphe::tetrahedral_graph[] = {
    1,2,3,4,0,2,3,4,0,3,4,0,-1
};
const int graphe::octahedral_graph[] = {
    1,3,6,5,4,0,2,3,4,5,6,0,3,5,6,0,4,5,6,0,-1
};
const int graphe::icosahedral_graph[] = {
    1,2,3,4,5,9,0,2,3,5,6,7,0,3,7,8,9,0,4,5,9,10,11,0,5,6,10,0,6,7,10,12,
    0,7,8,12,0,8,9,11,12,0,9,11,0,10,11,12,0,11,12,0,-1
};
const int graphe::ljubljana_graph_lcf[] = {
    47,-23,-31,39,25,-21,-31,-41,25,15,29,-41,-19,15,-49,33,39,-35,-21,17,-33,
    49,41,31,-15,-29,41,31,-15,-25,21,31,-51,-25,23,9,-17,51,35,-29,21,-51,-39,
    33,-9,-51,51,-47,-33,19,51,-21,29,21,-31,-39,0,2
};
const int graphe::harries_graph_lcf[] = {
    -35,-27,27,-23,15,-15,-9,-35,23,-27,27,9,15,-15,0,5
};
const int graphe::harries_wong_graph_lcf[] = {
    9,25,31,-17,17,33,9,-29,-15,-9,9,25,-25,29,17,-9,9,-27,35,-9,9,-17,21,27,-29,
    -9,-25,13,19,-9,-33,-17,19,-31,27,11,-25,29,-33,13,-13,21,-29,-21,25,9,-11,
    -19,29,9,-27,-19,-13,-35,-9,9,17,25,-9,9,27,-27,-21,15,-9,29,-29,33,-9,-25,0,1
};
const int graphe::balaban_10cage_lcf[] = {
    -9,-25,-19,29,13,35,-13,-29,19,25,9,-29,29,17,33,21,9,-13,-31,-9,25,17,9,-31,
    27,-9,17,-19,-29,27,-17,-9,-29,33,-25,25,-21,17,-17,29,35,-29,17,-17,21,-25,
    25,-33,29,9,17,-27,29,19,-17,9,-27,31,-9,-17,-25,9,31,13,-9,-21,-33,-17,-29,29,0,1
};
const int graphe::balaban_11cage_lcf[] = {
    44,26,-47,-15,35,-39,11,-27,38,-37,43,14,28,51,-29,-16,41,-11,-26,15,22,-51,
    -35,36,52,-14,-33,-26,-46,52,26,16,43,33,-15,17,-53,23,-42,-35,-28,30,-22,
    45,-44,16,-38,-16,50,-55,20,28,-17,-43,47,34,-26,-41,11,-36,-23,-16,41,17,
    -51,26,-33,47,17,-11,-20,-30,21,29,36,-43,-52,10,39,-28,-17,-52,51,26,37,
    -17,10,-10,-45,-34,17,-26,27,-21,46,53,-10,29,-50,35,15,-47,-29,-41,26,33,
    55,-17,42,-26,-36,16,0,1
};
const int graphe::foster_graph_lcf[] = {
    17,-9,37,-37,9,-17,0,15
};
const int graphe::bidiakis_cube_graph_lcf[] = {
    -6,4,-4,0,4
};
const int graphe::franklin_graph_lcf[] = {
    5,-5,0,6
};
const int graphe::frucht_graph_lcf[] = {
    -5,-2,-4,2,5,-2,2,5,-2,-5,4,2,0,1
};
const int graphe::biggs_smith_graph_lcf[] = {
    16,24,-38,17,34,48,-19,41,-35,47,-20,34,-36,21,14,48,-16,-36,-43,28,-17,21,
    29,-43,46,-24,28,-38,-14,-50,-45,21,8,27,-21,20,-37,39,-34,-44,-8,38,-21,25,
    15,-34,18,-28,-41,36,8,-29,-21,-48,-28,-20,-47,14,-8,-15,-27,38,24,-48,-18,
    25,38,31,-25,24,-46,-14,28,11,21,35,-39,43,36,-38,14,50,43,36,-11,-36,-24,
    45,8,19,-25,38,20,-24,-14,-21,-8,44,-31,-38,-28,37,0,1
};
const int graphe::folkman_graph_lcf[] = {
    5,-7,-7,5,0,5
};
const int graphe::gray_graph_lcf[] = {
    -25,7,-7,13,-13,25,0,9
};
const int graphe::tutte_12cage_lcf[] = {
    17,27,-13,-59,-35,35,-11,13,-53,53,-27,21,57,11,-21,-57, 59,-17,0,7
};
const int graphe::tutte_8cage_lcf[] = {
    -13,-9,7,-7,9,13,0,5
};
const int graphe::f26a_graph_lcf[] = {
    -7,7,0,13
};
const int graphe::blanusa_graph[] = {
    1,2,6,8,0,2,3,12,0,3,4,13,0,4,5,14,0,5,6,18,0,6,7,0,7,9,17,0,8,17,18,
    0,9,10,18,0,10,11,16,0,11,13,14,0,12,14,15,0,13,15,0,15,16,0,16,17,0,-1
};
const int graphe::bull_graph[] = {
    1,2,3,0,2,3,4,0,3,5,0,-1
};
const int graphe::butterfly_graph[] = {
    1,2,3,4,5,0,2,3,0,4,5,0,-1
};
const int graphe::diamond_graph[] = {
    1,2,3,4,0,2,3,0,3,4,0,-1
};
const int graphe::chvatal_graph[] = {
    1,2,4,5,12,0,2,3,6,7,0,3,4,8,9,0,4,10,11,0,5,6,8,9,0,6,11,10,0,7,8,10,
    11,0,8,12,0,9,10,12,0,11,12,0,-1
};
const int graphe::moser_spindle_graph[] = {
    1,2,5,6,7,0,2,3,7,0,3,4,7,0,4,5,6,0,5,6,0,-1
};
const int graphe::errera_graph[] = {
    1,2,3,4,14,15,0,2,3,13,15,17,0,3,13,14,16,0,4,5,6,14,15,0,5,6,7,8,14,16,
    0,6,7,9,15,17,0,7,8,9,10,0,8,10,11,16,0,9,10,12,17,0,10,11,12,0,11,12,
    13,16,0,12,13,17,0,13,16,17,0,14,16,0,15,17,0,-1
};
const int graphe::goldner_harary_graph[] = {
    1,7,9,10,0,2,8,9,11,0,3,8,10,11,0,4,7,9,11,0,5,7,10,11,0,6,8,9,10,0,
    7,9,10,11,0,8,9,10,11,0,9,10,11,0,10,11,0,-1
};
const int graphe::golomb_graph[] = {
    1,2,3,4,0,2,3,6,0,3,8,0,4,5,9,10,0,5,6,10,0,6,7,10,0,7,8,10,0,8,9,10,
    0,9,10,0,-1
};
const int graphe::hoffman_graph_matrix[8][8] = {
    {1,1,1,1,0,0,0,0},
    {1,1,1,0,1,0,0,0},
    {1,0,0,1,0,1,1,0},
    {0,1,0,1,0,1,0,1},
    {0,0,1,1,0,0,1,1},
    {1,0,0,0,1,1,1,0},
    {0,1,0,0,1,1,0,1},
    {0,0,1,0,1,0,1,1}
};
const int graphe::poussin_graph[] = {
    1,2,3,4,5,10,0,2,3,5,7,9,15,0,3,4,6,8,15,0,4,6,10,0,5,7,10,0,6,8,10,11,
    12,0,7,9,10,11,13,0,8,12,14,15,0,9,13,14,15,0,10,11,0,11,12,13,0,12,
    13,14,0,13,14,0,14,15,0,-1
};
const int graphe::wagner_graph[] = {
    1,2,5,8,0,2,3,6,0,3,4,7,0,4,5,8,0,5,6,0,6,7,0,7,8,0,-1
};
const int graphe::tutte_graph[] = {
    1,2,3,4,0,2,45,46,0,3,44,46,0,4,42,43,0,5,6,7,46,0,6,41,45,0,7,40,44,
    0,8,9,12,21,0,9,22,36,0,10,11,17,20,0,11,16,18,0,12,13,19,0,13,22,25,
    0,14,20,23,30,0,15,19,25,31,0,16,21,24,0,17,23,24,0,18,26,37,0,19,29,
    0,20,26,0,21,29,0,22,28,0,23,27,0,24,27,0,25,28,0,26,32,0,27,35,0,28,
    34,0,29,33,0,30,32,35,0,31,33,34,0,32,42,0,33,43,0,34,43,0,35,42,0,36,
    39,44,0,37,38,39,0,38,41,45,0,39,40,0,40,41,-1
};
const int graphe::brinkmann_graph[] = {
    1,2,3,4,5,0,2,12,17,19,0,3,13,16,18,0,4,10,14,21,0,5,11,15,20,0,6,16,17,
    20,21,0,7,9,15,18,21,0,8,9,14,19,20,0,9,12,13,0,10,12,18,20,0,11,13,19,
    21,0,12,16,0,13,17,0,14,15,17,0,15,16,0,18,19,0,-1
};
const int graphe::barnette_bosak_lederberg_graph[] = {
    1,2,3,4,0,2,34,38,0,3,33,37,0,4,35,36,0,5,9,15,19,0,6,10,16,20,0,7,10,
    11,18,0,8,9,12,17,0,9,14,0,10,13,0,11,13,22,0,12,14,21,0,13,25,0,14,
    26,0,15,17,31,0,16,18,32,0,17,24,0,18,23,0,19,26,30,0,20,25,29,0,21,
    24,28,0,22,23,27,0,23,36,0,24,35,0,25,34,0,26,33,0,27,28,38,0,28,37,
    0,29,32,34,0,30,31,33,0,31,35,0,32,36,0,33,30,0,37,38,0,-1
};
const int graphe::double_star_snark[] = {
    1,2,3,4,0,2,25,28,0,3,27,30,0,4,26,29,0,5,8,9,16,0,6,7,12,15,0,7,11,
    16,0,8,10,15,0,9,13,17,0,10,14,18,0,11,13,18,0,12,14,17,0,13,21,0,
    14,22,0,15,27,0,16,26,0,17,28,0,18,28,0,19,22,24,30,0,20,21,23,29,0,
    21,30,0,22,29,0,23,25,27,0,24,25,26,0,-1
};
const int graphe::doyle_graph[] = {
    1,2,3,4,5,0,2,20,22,27,0,3,21,23,26,0,4,16,19,24,0,5,17,18,25,0,6,7,8,
    24,25,0,7,14,19,27,0,8,15,18,26,0,9,12,13,16,17,0,10,12,14,24,26,0,11,
    13,15,25,27,0,12,15,21,0,13,14,20,0,14,22,0,15,23,0,16,23,25,0,17,22,
    24,0,18,21,27,0,19,20,26,0,20,23,0,21,22,0,-1
};
const int graphe::meringer_graph[] = {
    1,2,3,4,5,6,0,2,15,21,24,25,0,3,18,23,26,30,0,4,17,19,22,28,0,5,16,20,
    27,29,0,6,11,12,13,14,0,7,8,15,23,27,28,0,8,19,25,29,30,0,9,10,17,20,
    21,26,0,10,16,18,22,24,0,11,15,18,19,20,0,12,21,22,23,29,0,13,17,24,27,
    30,0,14,16,25,26,28,0,15,16,17,0,16,30,0,17,29,0,18,28,29,0,19,24,26,
    0,20,23,25,0,21,28,30,0,22,25,27,0,23,24,0,26,27,0,-1
};
const int graphe::robertson_wegner_graph[] = {
    1,2,3,4,5,6,0,2,11,22,23,28,0,3,12,17,26,30,0,4,13,16,27,29,0,5,14,19,
    21,24,0,6,15,18,20,25,0,7,16,18,21,28,30,0,8,17,19,20,28,29,0,9,10,23,
    24,26,29,0,10,22,25,27,30,0,11,20,21,26,27,0,12,15,24,27,28,0,13,14,25,
    26,28,0,14,15,23,30,0,15,22,29,0,16,17,22,24,0,17,23,25,0,18,19,23,27,
    0,19,22,26,0,20,24,30,0,21,25,29,0,-1
};
const int graphe::wong_graph[] = {
    1,2,3,4,5,6,0,2,11,22,23,24,0,3,18,28,29,30,0,4,14,17,19,26,0,5,13,16,
    21,25,0,6,12,15,20,27,0,7,11,15,16,17,18,0,8,12,19,22,25,28,0,9,13,20,
    24,26,30,0,10,14,21,23,27,29,0,11,12,13,14,0,12,26,29,0,13,27,28,0,14,
    25,30,0,15,21,24,28,0,16,19,23,30,0,17,20,22,29,0,18,25,26,27,0,19,24,
    27,0,20,23,25,0,21,22,26,0,22,30,0,23,28,0,24,29,0,-1
};
const char* const graphe::gewirtz_words[] = {
    "abcilu", "abdfrs", "abejop", "abgmnq", "acdghp", "acfjnt", "ackmos",
    "ademtu", "adjklq", "aefgik", "aehlns", "afhoqu", "aglort", "ahijmr",
    "aipqst", "aknpru", "bcdekn", "bchjqs", "bcmprt", "bdgijt", "bdhlmo",
    "beflqt", "beghru", "bfhinp", "bfjkmu", "bgklps", "bikoqr", "bnostu",
    "cdfimq", "cdjoru", "cefpsu", "cegjlm", "cehiot", "cfhklr", "cginrs",
    "cgkqtu", "clnopq", "degoqs", "deilpr", "dfglnu", "dfkopt", "dhiksu",
    "dhnqrt", "djmnps", "efmnor", "ehkmpq", "eijnqu", "ejkrst", "fghmst",
    "fgjpqr", "fijlos", "ghjkno", "gimopu", "hjlptu", "iklmnt", "lmqrsu"
};
const int graphe::harborth_graph[] = {
    1,2,3,4,5,0,2,3,47,52,0,3,2,46,52,0,4,5,50,51,0,5,48,49,0,6,10,43,51,
    0,7,10,11,50,0,8,11,12,13,50,0,9,19,27,42,50,0,10,11,43,0,11,12,0,12,
    13,14,0,13,14,42,0,14,41,42,0,15,36,38,39,48,0,16,34,35,37,48,0,17,18,
    19,20,41,0,18,20,21,41,0,19,20,27,0,20,21,0,21,22,23,0,22,23,26,31,0,
    23,25,26,0,24,25,28,29,30,0,25,26,30,0,26,31,0,27,31,34,0,28,29,30,35,
    0,29,32,35,0,30,34,0,31,34,0,32,33,35,37,0,33,36,37,38,0,36,38,39,0,
    37,38,0,39,40,45,0,40,44,45,49,0,41,42,0,43,46,51,0,44,45,47,49,0,45,
    47,0,46,51,52,0,47,52,0,48,49,0,-1
};
const int graphe::kittell_graph[] = {
    1,2,3,5,21,23,0,2,3,4,5,9,0,3,4,21,22,0,4,9,10,22,0,5,9,18,23,0,6,7,11,
    16,21,23,0,7,8,12,16,21,0,8,12,17,21,22,0,9,10,18,0,10,18,20,22,0,11,
    13,16,19,23,0,12,14,16,17,0,13,14,15,16,19,0,14,15,16,17,0,15,17,19,20,
    0,17,20,22,0,18,19,20,23,0,19,20,23,0,20,22,0,21,22,23,0,-1
};
const int graphe::krackhardt_kite_graph[] = {
    1,2,0,2,5,0,3,6,8,10,0,4,7,9,10,0,5,8,9,0,6,7,8,10,0,7,9,10,0,8,9,10,0,9,10,0,-1 
};
const int graphe::meredith_graph[] = {
    1,2,3,4,5,0,2,65,69,70,0,3,66,69,70,0,4,68,69,70,0,5,67,69,70,0,6,7,38,
    39,40,0,7,41,42,43,0,8,13,62,63,64,0,9,12,62,63,64,0,10,15,41,42,43,0,
    11,14,38,39,40,0,12,16,20,21,0,13,17,18,19,0,14,16,20,21,0,15,17,18,19,
    0,16,22,23,0,17,24,25,0,18,24,25,0,19,24,25,0,20,22,23,0,21,22,23,0,
    22,26,0,23,27,0,24,28,0,25,29,0,26,32,33,34,0,27,32,33,34,0,28,35,36,
    37,0,29,35,36,37,0,30,31,32,33,34,0,31,35,36,37,0,32,48,0,33,48,0,34,
    48,0,35,49,0,36,49,0,37,49,0,38,44,45,0,39,44,45,0,40,44,45,0,41,46,
    47,0,42,46,47,0,43,46,47,0,44,52,0,45,53,0,46,50,0,47,51,0,48,55,0,
    49,54,0,50,56,57,58,0,51,56,57,58,0,52,59,60,61,0,53,59,60,61,0,54,59,
    60,61,0,55,56,57,58,0,56,67,0,57,67,0,58,67,0,59,68,0,60,68,0,61,68,
    0,62,65,66,0,63,65,66,0,64,65,66,0,-1   
};
const int graphe::perkel_graph[] = {
    1,2,3,4,5,6,7,0,2,28,32,33,34,35,0,3,31,40,41,44,45,0,4,36,39,47,52,55,
    0,5,37,38,46,53,54,0,6,29,43,49,50,57,0,7,30,42,48,51,56,0,8,23,25,27,
    28,55,57,0,9,22,24,26,28,54,56,0,10,12,14,17,31,54,57,0,11,13,15,16,31,
    55,56,0,12,13,19,35,40,50,0,13,18,34,41,51,0,14,20,24,36,38,44,0,15,21,
    25,37,39,45,0,16,19,22,30,36,49,0,17,18,23,29,37,48,0,18,25,30,33,38,0,
    19,24,29,32,39,0,20,23,26,40,42,47,0,21,22,27,41,43,46,0,22,42,50,53,0,
    23,43,51,52,0,24,34,46,52,0,25,35,47,53,0,26,27,32,45,48,0,27,33,44,49,
    0,28,29,30,31,0,29,45,53,0,30,44,52,0,31,42,43,0,32,37,51,55,0,33,36,
    50,54,0,34,43,45,47,0,35,42,44,46,0,36,37,43,0,37,42,0,38,45,49,56,0,
    39,44,48,57,0,40,49,53,55,0,41,48,52,54,0,46,51,57,0,47,50,56,0,48,50,
    0,49,51,0,52,53,0,54,55,0,56,57,0,-1
};
const int graphe::sousselier_graph[] = {
    1,9,10,16,0,2,8,11,16,0,3,6,14,16,0,4,7,15,16,0,5,12,13,16,0,6,7,13,
    0,7,12,0,8,9,14,0,9,15,0,10,13,14,0,11,12,15,0,12,14,0,13,15,0,-1
};
const int graphe::walther_graph[] = {
    1,25,0,2,10,0,3,9,0,4,23,24,0,5,21,24,0,6,22,23,0,7,19,20,0,8,21,22,
    0,9,18,0,10,17,0,11,16,24,25,0,12,16,20,25,0,13,16,19,21,0,14,15,19,
    22,0,15,17,18,0,17,23,0,18,20,0,-1
};
const int graphe::watkins_snark[] = {
    1,2,3,4,0,2,47,48,0,3,45,50,0,4,46,49,0,5,6,7,44,0,6,37,39,0,7,36,38,
    0,8,9,33,35,0,9,14,15,0,10,23,26,34,0,11,24,25,34,0,12,15,21,40,0,13,
    20,22,40,0,14,21,29,0,15,16,0,16,20,29,0,17,25,31,32,0,18,19,31,33,0,
    19,25,30,0,20,21,0,22,23,27,0,23,28,0,24,26,28,0,26,27,0,27,43,0,28,
    43,0,29,41,0,30,32,42,0,31,42,0,32,41,0,33,39,0,34,38,0,35,36,37,0,
    36,39,0,37,38,0,40,48,0,41,46,0,42,46,0,43,45,0,44,45,47,0,47,50,0,
    48,49,0,49,50,0,-1
};
const int graphe::wells_graph[] = {
    1,2,3,4,5,6,0,2,13,24,25,26,0,3,20,30,31,32,0,4,14,17,21,27,0,5,16,18,
    23,29,0,6,15,19,22,28,0,7,12,13,17,18,19,0,8,12,20,27,28,29,0,9,12,14,
    23,26,30,0,10,12,15,21,24,32,0,11,12,16,22,25,31,0,13,14,15,16,0,14,29,
    32,0,15,27,31,0,16,28,30,0,17,23,24,28,0,18,22,26,27,0,19,21,25,29,0,
    20,21,22,23,0,21,26,0,22,24,0,23,25,0,24,30,0,25,32,0,26,31,0,27,30,
    0,28,32,0,29,31,0,-1
};
const int graphe::wiener_araya_graph[] = {
    1,20,22,30,0,2,19,21,29,0,3,5,15,30,0,4,6,16,29,0,5,24,27,0,6,23,28,0,
    7,12,24,26,0,8,11,23,25,0,9,14,17,26,0,10,13,18,25,0,11,13,28,0,12,14,
    27,0,13,34,0,14,33,0,15,20,35,0,16,19,36,0,17,21,38,0,18,22,37,0,19,
    38,0,20,37,0,21,41,0,22,42,0,23,40,0,24,39,0,25,42,0,26,41,0,27,35,
    0,28,36,0,29,40,0,30,39,0,31,32,35,37,0,32,36,38,0,33,35,38,0,34,36,
    37,0,39,41,42,0,40,41,42,0,-1
};
const int graphe::markstroem_graph[] = {
    1,2,9,10,0,2,3,10,0,3,4,13,0,4,5,24,0,5,6,24,0,6,7,18,0,7,8,23,0,8,9,23,0,
    9,21,0,10,11,0,11,12,13,0,12,13,14,0,14,15,16,0,15,16,17,0,16,20,0,17,18,
    19,0,18,19,0,19,24,0,20,21,22,0,21,22,0,22,23,0,-1
};

/* the list of all special graphs alongside with the dimensions |V| and |E| */
const graphe::spcgraph graphe::special_graph[] = {
    {"balaban10", 70, 105},
    {"balaban11", 112, 168},
    {"barnette-bosak-lederberg", 38, 57},
    {"bidiakis", 12, 18},
    {"biggs-smith", 102, 153},
    {"blanusa", 18, 27},
    {"brinkmann", 21, 42},
    {"brouwer-haemers", 81, 810},
    {"bull", 5, 5},
    {"butterfly", 5, 6},
    {"clebsch", 16, 40},
    {"chvatal", 12, 24},
    {"coxeter", 28, 42},
    {"desargues", 20, 30},
    {"diamond", 4, 5},
    {"dodecahedron", 20, 30},
    {"double-star", 30, 45},
    {"doyle", 27, 54},
    {"duerer", 12, 18},
    {"dyck", 32, 48},
    {"errera", 17, 45},
    {"f26a", 26, 39},
    {"folkman", 20, 40},
    {"foster", 90, 135},
    {"franklin", 12, 18},
    {"frucht", 12, 18},
    {"gewirtz", 56, 280},
    {"goldner-harary", 11, 27},
    {"golomb", 10, 18},
    {"gosset", 56, 756},
    {"gray", 54, 81},
    {"grinberg", 46, 69},
    {"groetzsch", 11, 20},
    {"harborth", 52, 103},
    {"harries", 70, 105},
    {"harries-wong", 70, 105},
    {"heawood", 14, 21},
    {"herschel", 11, 18},
    {"higman-sims", 100, 1100},
    {"hoffman", 16, 32},
    {"hoffman-singleton", 50, 175},
    {"icosahedron", 12, 30},
    {"kittell", 23, 63},
    {"krackhardt", 10, 18},
    {"levi", 30, 45},
    {"ljubljana", 112, 168},
    {"markstroem", 24, 36},
    {"mcgee", 24, 36},
    {"meredith", 70, 140},
    {"meringer", 30, 75},
    {"moser", 7, 11},
    {"moebius-kantor", 16, 24},
    {"nauru", 24, 36},
    {"octahedron", 6, 12},
    {"pappus", 18, 27},
    {"petersen", 10, 15},
    {"perkel", 57, 171},
    {"poussin", 15, 39},
    {"robertson", 19, 38},
    {"robertson-wegner", 30, 75},
    {"shrikhande", 16, 48},
    {"schlaefli", 27, 216},
    {"sousselier", 16, 27},
    {"sylvester", 36, 90},
    {"szerekes", 50, 75},
    {"tehtrahedron", 4, 6},
    {"tietze", 12, 18},
    {"soccerball", 60, 90},
    {"tutte", 46, 69},
    {"tutte12", 126, 189},
    {"wagner", 8, 12},
    {"walther", 25, 31},
    {"watkins", 50, 75},
    {"wells", 32, 80},
    {"wiener-araya", 42, 67},
    {"wong", 30, 75},
    {"", 0, 0} // termination entry
};

#ifndef NO_NAMESPACE_GIAC
}
#endif // ndef NO_NAMESPACE_GIAC
