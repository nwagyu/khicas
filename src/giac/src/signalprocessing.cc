/*
 * signalprocessing.cc
 *
 * (c) 2018-2019 by Luka Marohnić
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
#include "signalprocessing.h"
#include "graphtheory.h"

#ifdef HAVE_LIBSAMPLERATE
#include "samplerate.h"
#endif

using namespace std;

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

#define _SP_BAD_SOUND_DATA "Invalid sound data"
#define _SP_INVALID_RANGE "Invalid range specification"
#define _SP_BAD_WINDOW "Invalid window parameters"

gen generr(const char* msg,bool translate) {
    string m(translate?gettext(msg):msg);
    m.append(".");
    return gensizeerr(m.c_str());
}

gen generrtype(const char* msg,bool translate) {
    string m(translate?gettext(msg):msg);
    m.append(".");
    return gentypeerr(m.c_str());
}

gen generrdim(const char* msg,bool translate) {
    string m(translate?gettext(msg):msg);
    m.append(".");
    return gendimerr(m.c_str());
}

void print_error(const char *msg,GIAC_CONTEXT) {
    *logptr(contextptr) << gettext("Error") << ": " << gettext(msg) << "\n";
}

void print_warning(const char *msg,GIAC_CONTEXT) {
    *logptr(contextptr) << gettext("Warning") << ": " << gettext(msg) << "\n";
}

/*
 * Return TRUE iff g is a real constant.
 */
bool is_real_number(const gen &g,GIAC_CONTEXT) {
    gen eg=_evalf(g,contextptr);
    if (!eg.is_approx())
        return false;
    switch (eg.type) {
    case _INT_: case _DOUBLE_: case _FLOAT_: case _ZINT: case _REAL:
        return true;
    case _CPLX: 
        return (is_zero(*(eg._CPLXptr+1),contextptr));
    default:
        break;
    }
    return false;
}

/*
 * Convert g to inexact real constant.
 */
gen to_real_number(const gen &g,GIAC_CONTEXT) {
    gen eg=_evalf(g,contextptr);
    if (eg.type==_CPLX)
        return *(eg._CPLXptr);
    else return eg;
}

int nextpow2(int n) {
    int m=1;
    while (m<n) m*=2;
    return m;
}

bool is_sound_data(const gen &g,int &nc,int &bd,int &sr,int &len) {
    if (g.type!=_VECT)
        return false;
    const vecteur &gv=*g._VECTptr;
    if (gv.size()<2 || gv.front().type!=_VECT || gv.front()._VECTptr->size()!=4)
        return false;
    vecteur &header=*gv.front()._VECTptr;
    if (!is_integer_vecteur(header))
        return false;
    int n=header[0].val,b=header[1].val,s=header[2].val,l=(8*header[3].val)/(n*b);
    if (n<1 || (b!=8 && b!=16 && b!=24) || s<1 || l<0 ||
            int(gv.size())!=n+1)
        return false;
    for (const_iterateur it=gv.begin()+1;it!=gv.end();++it) {
        if (it->type!=_VECT)
            return false;
        if (int(it->_VECTptr->size())!=l)
            return false;
    }
    nc=n; bd=b; sr=s; len=l;
    return true;
}

vecteur decode_chdata(const vecteur &data,int bd,int start,int len) {
    if (data.empty() || (bd!=8 && bd!=16))
        return vecteur(0);
    vecteur res;
    int l=len<0?data.size():len;
    res.reserve(l);
    int v;
    double denom=std::pow(2,bd-1);
    for (const_iterateur it=data.begin()+start;it!=data.begin()+start+l;++it) {
        if (!it->is_integer())
            return vecteur(0);
        v=it->val;
        if (bd==8) v-=127;
        else if (bd==16 && v>=32768) v-=65536;
        res.push_back(double(v)/denom);
    }
    return res;
}

vecteur encode_chdata(const vecteur &data,int bd,double ratio,GIAC_CONTEXT) {
    if (data.empty() || (bd!=8 && bd!=16))
        return vecteur(0);
    vecteur res;
    res.reserve(data.size());
    double v,fac=std::pow(2.0,bd-1)-1.0;
    int k;
    for (const_iterateur it=data.begin();it!=data.end();++it) {
        if (!is_real_number(*it,contextptr))
            return vecteur(0);
        v=std::max(-1.0,std::min(1.0,ratio*to_real_number(*it,contextptr).to_double(contextptr)));
        k=std::floor(fac*v);
        if (bd==8) k+=127;
        res.push_back(k);
    }
    return res;
}

int read_channel_data(const gen &g,int &nc,vector<vecteur*> &data) {
    int ret;
    if (ckmatrix(g)) {
        const matrice &chns=*g._VECTptr;
        nc=chns.size();
        ret=chns.front()._VECTptr->size();
        const_iterateur jt,jtend;
        for (jt=chns.begin(),jtend=chns.end();jt!=jtend;++jt)
            data.push_back(jt->_VECTptr);
    } else {
        data.push_back(g._VECTptr);
        ret=g._VECTptr->size();
    }
    return ret;
}

gen _createwav(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    int len=0,bd=16,sr=44100,nc=1;
    double ratio=1.0,norm=-1;
    vector<vecteur*> data;
    if (g.is_integer()) {
        if (g.val<1)
            return generr("Expected a length in samples (positive integer)");
        len=g.val;
    } else if (g.type==_VECT) {
        if (g.subtype==_SEQ__VECT) {
            const vecteur &args=*g._VECTptr;
            double secs=-1,plen=-1;
            for (const_iterateur it=args.begin();it!=args.end();++it) {
                if (is_equal(*it)) {
                    const gen &lh=it->_SYMBptr->feuille._VECTptr->front();
                    const gen &rh=it->_SYMBptr->feuille._VECTptr->back();
                    if (lh==at_channels) {
                        if (!rh.is_integer() || (nc=rh.val)<1)
                            return generr("Number of channels must be a positive integer");
                    } else if (lh==at_bit_depth) {
                        if (!rh.is_integer() || (bd=rh.val)<0 || (bd!=8 && bd!=16))
                            return generr("Bit depth must be a nonnegative integer (8 or 16)");
                    } else if (lh==at_samplerate) {
                        if (!rh.is_integer() || (sr=rh.val)<1)
                            return generr("Sample rate must be a positive integer");
                    } else if (lh==at_size) {
                        if (!rh.is_integer() || (plen=rh.val)<1)
                            return generr("Size must be a positive integer");
                    } else if (lh==at_normalize) {
                        if (!is_real_number(rh,contextptr))
                            return generr("Normalization level must be a real number");
                        norm=-to_real_number(rh,contextptr).to_double(contextptr);
                    } else if (lh==at_duration) {
                        if (!is_real_number(rh,contextptr))
                            return generr("Duration must be a real number");
                        secs=to_real_number(rh,contextptr).to_double(contextptr);
                    } else if (lh==at_channel_data) {
                        if (rh.type!=_VECT)
                            return generr("Channel data must be a list");
                        plen=read_channel_data(rh,nc,data);
                    } else return gensizeerr(gettext("Invalid argument ")+print_INT_(it-args.begin()+1)+".");
                } else if (it->type==_VECT)
                    plen=read_channel_data(*it,nc,data);
            }
            if (secs>=0)
                len=1+std::floor(secs*sr);
            else if ((len=plen)<0)
                return gensizeerr(contextptr);
        } else {
            if (is_sound_data(g,nc,bd,sr,len)) {
                for (int i=1;i<=nc;++i) {
                    data.push_back(g._VECTptr->at(i)._VECTptr);
                }
            } else len=read_channel_data(g,nc,data);
        }
    } else return generrtype("Invalid input argument");
    if (!data.empty() && int(data.size())!=nc)
        return generrdim("No data found or invalid number of channels");
    if (len==0) {
        for (vector<vecteur*>::const_iterator it=data.begin();it!=data.end();++it) {
            len=std::max(len,(int)(*it)->size());
        }
    }
    if (norm>=0) {
        double maxamp=0;
        for (vector<vecteur*>::const_iterator it=data.begin();it!=data.end();++it) {
            maxamp=std::max(maxamp,to_real_number(_max(_abs(**it,contextptr),contextptr),contextptr).to_double(contextptr));
        }
        ratio=std::pow(10.0,-norm/20.0)/maxamp;
    }
    /* create wav data */
    vecteur header,res;
    header.resize(4);
    header[0]=nc;
    header[1]=bd;
    header[2]=sr;
    header[3]=(bd*len*nc)/8;
    res.push_back(header);
    for (vector<vecteur*>::const_iterator it=data.begin();it!=data.end();++it) {
        vecteur enc=encode_chdata(**it,bd,ratio,contextptr);
        enc.resize(len,0);
        res.push_back(enc);
    }
    return res;
}
static const char _createwav_s []="createwav";
static define_unary_function_eval (__createwav,&_createwav,_createwav_s);
define_unary_function_ptr5(at_createwav,alias_at_createwav,&__createwav,0,true)

gen _plotwav(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return generrtype("Invalid input argument");
    vecteur opts;
    if (g.subtype==_SEQ__VECT)
        opts=vecteur(g._VECTptr->begin()+1,g._VECTptr->end());
    const gen &wav=g.subtype==_SEQ__VECT?g._VECTptr->front():g;
    int nc,bd,sr,len,slice_start=0;
    if (!is_sound_data(wav,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    for (const_iterateur it=opts.begin();it!=opts.end();++it) {
        if (is_equal(*it)) {
            const gen &lh=it->_SYMBptr->feuille._VECTptr->front();
            const gen &rh=it->_SYMBptr->feuille._VECTptr->back();
            if (lh==at_range) {
                if (rh.is_symb_of_sommet(at_interval)) {
                    gen a=rh._SYMBptr->feuille._VECTptr->front();
                    gen b=rh._SYMBptr->feuille._VECTptr->back();
                    if (!is_real_number(a,contextptr) || !is_real_number(b,contextptr) ||
                            is_positive((a=to_real_number(a,contextptr))-(b=to_real_number(b,contextptr)),contextptr))
                        return generr(_SP_INVALID_RANGE);
                    slice_start=std::max(slice_start,_floor(a*gen(sr),contextptr).val);
                    len=std::min(len-slice_start,_floor(b*gen(sr),contextptr).val-slice_start);
                } else if (rh.type==_VECT) {
                    if (rh._VECTptr->size()!=2 || !rh._VECTptr->front().is_integer() ||
                            !rh._VECTptr->back().is_integer())
                        return generr(_SP_INVALID_RANGE);
                    slice_start=std::max(slice_start,rh._VECTptr->front().val);
                    len=std::min(len-slice_start,rh._VECTptr->back().val-slice_start);
                } else return generr(_SP_INVALID_RANGE);
            }
        }
    }
    if (len==0)
        return vecteur(0);
    double dur=double(len)/double(sr),width=1500.0,sc=double(len)/width,x,y=nc==1?0.0:1.0;
    double fmax,fmin,s,t0=double(slice_start)/double(sr),dt=1.0/sr;
    int step=std::floor(sc),bnd;
    vecteur drawing,tvec;
    if (sc<=1.0) tvec.resize(len);
    for (int c=1;c<=nc;++c) {
        vecteur data=decode_chdata(*wav._VECTptr->at(c)._VECTptr,bd,slice_start,len),unodes,lnodes;
        if (sc<=1.0) {
            for (iterateur it=data.begin();it!=data.end();++it) {
                *it+=gen(y);
                tvec[it-data.begin()]=t0+(it-data.begin())*dt;
            }
            drawing.push_back(_listplot(_zip(makesequence(at_makevector,tvec,data),contextptr),contextptr));
        } else {
            unodes.reserve(std::floor(width));
            lnodes.reserve(std::floor(width));
            for (int i=0;i<len;i+=step) {
                fmax=fmin=0.0;
                bnd=std::min(i+step,len);
                for (int j=i;j<bnd;++j) {
                    s=to_real_number(data[j],contextptr).to_double(contextptr);
                    if (s>fmax) fmax=s;
                    else if (s<fmin) fmin=s;
                }
                x=t0+dur*(i/step)/width;
                unodes.push_back(makecomplex(x,y+fmax));
                lnodes.push_back(makecomplex(x,y+fmin));
            }
            std::reverse(lnodes.begin(),lnodes.end());
            vecteur nodes=mergevecteur(unodes,lnodes);
            nodes.push_back(symbolic(at_equal,makesequence(at_display,_FILL_POLYGON)));
            drawing.push_back(_polygone(change_subtype(nodes,_SEQ__VECT),contextptr));
            drawing.push_back(_segment(makesequence(makecomplex(t0,y),makecomplex(t0+dur,y)),contextptr));
        }
        y-=2.0;
    }
    return drawing;
}
static const char _plotwav_s []="plotwav";
static define_unary_function_eval (__plotwav,&_plotwav,_plotwav_s);
define_unary_function_ptr5(at_plotwav,alias_at_plotwav,&__plotwav,0,true)

gen _stereo2mono(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    int nc,bd,sr,len;
    if (!is_sound_data(g,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    if (nc==1)
        return g;
    vecteur data(len,0);
    for (int i=1;i<=nc;++i) {
        const vecteur &chan=*g._VECTptr->at(i)._VECTptr;
        addvecteur(data,decode_chdata(chan,bd),data);
    }
    multvecteur(fraction(1,nc),data,data);
    vecteur header=*g._VECTptr->front()._VECTptr,edata=encode_chdata(data,bd,1.0,contextptr);
    header[0]=1;
    header[3]=header[3]/gen(nc);
    return makevecteur(header,edata);
}
static const char _stereo2mono_s []="stereo2mono";
static define_unary_function_eval (__stereo2mono,&_stereo2mono,_stereo2mono_s);
define_unary_function_ptr5(at_stereo2mono,alias_at_stereo2mono,&__stereo2mono,0,true)

gen _plotspectrum(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return gentypeerr(contextptr);
    int nc,bd,sr,len,lfreq,ufreq;
    if (is_sound_data(g,nc,bd,sr,len)) {
        gen intrv=symbolic(at_interval,makesequence(0,sr/2));
        return _plotspectrum(makesequence(g,symbolic(at_equal,makesequence(at_range,intrv))),contextptr);
    } else if (g.subtype==_SEQ__VECT) {
        const vecteur &gv=*g._VECTptr;
        if (gv.size()!=2)
            return generr("Expected two input arguments");
        if (!is_sound_data(gv.front(),nc,bd,sr,len))
            return generrtype(_SP_BAD_SOUND_DATA);
        vecteur data;
        if (nc>1)
            data=decode_chdata(*_stereo2mono(gv.front(),contextptr)._VECTptr->at(1)._VECTptr,bd);
        else data=decode_chdata(*gv.front()._VECTptr->at(1)._VECTptr,bd);
        len=data.size();
        if (!is_equal(gv.back()))
            return gensizeerr(contextptr);
        const gen &lh=gv.back()._SYMBptr->feuille._VECTptr->front();
        const gen &rh=gv.back()._SYMBptr->feuille._VECTptr->back();
        if (lh==at_range) {
            gen a,b;
            if (rh.type==_VECT) {
                if (rh._VECTptr->size()!=2 || !is_integer_vecteur(*rh._VECTptr))
                    return generr(_SP_INVALID_RANGE);
                a=rh._VECTptr->front();
                b=rh._VECTptr->back();
            } else if (rh.is_symb_of_sommet(at_interval)) {
                a=rh._SYMBptr->feuille._VECTptr->front();
                b=rh._SYMBptr->feuille._VECTptr->back();
            } else return generr(_SP_INVALID_RANGE);
            if (!is_integer(a) || !is_integer(b))
                return generr(_SP_INVALID_RANGE);
            lfreq=std::max(0,a.val);
            ufreq=std::min(sr/2,b.val);
            if (lfreq>=ufreq)
                return generr(_SP_INVALID_RANGE);
        }
        int n=nextpow2(len);
        data.resize(n,0);
        vecteur spec=*_fft(data,contextptr)._VECTptr;
        vecteur nodes;
        int dfreq=ufreq-lfreq,n1=((long)n*(long)lfreq)/(long)sr,n2=((long)n*(long)ufreq)/(long)sr;
        int width=std::min(dfreq,std::max(1500,dfreq/5));
        nodes.reserve(std::floor(width));
        int step=(n2-n1)/(2*width);
        if (step==0) step=1;
        double f;
        for (int i=n1;i<n2;++i) {
            f=std::max(f,to_real_number(pow(_abs(spec[i],contextptr),2)/gen(n),contextptr).to_double(contextptr));
            if (i%step==0) {
                nodes.push_back(makecomplex(((long)i*(long)sr)/(long)n,f));
                f=0;
            }
        }
        nodes.push_back(ufreq+1);
        nodes.push_back(lfreq);
        nodes.push_back(symbolic(at_equal,makesequence(at_display,_FILL_POLYGON)));
        return _polygone(change_subtype(nodes,_SEQ__VECT),contextptr);
    } else return gentypeerr(contextptr);
}
static const char _plotspectrum_s []="plotspectrum";
static define_unary_function_eval (__plotspectrum,&_plotspectrum,_plotspectrum_s);
define_unary_function_ptr5(at_plotspectrum,alias_at_plotspectrum,&__plotspectrum,0,true)

gen _channels(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    int nc,bd,sr,len;
    if (!is_sound_data(g,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    return nc;
}
static const char _channels_s []="channels";
static define_unary_function_eval (__channels,&_channels,_channels_s);
define_unary_function_ptr5(at_channels,alias_at_channels,&__channels,0,true)

gen _bit_depth(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    int nc,bd,sr,len;
    if (!is_sound_data(g,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    return bd;
}
static const char _bit_depth_s []="bit_depth";
static define_unary_function_eval (__bit_depth,&_bit_depth,_bit_depth_s);
define_unary_function_ptr5(at_bit_depth,alias_at_bit_depth,&__bit_depth,0,true)

gen _samplerate(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    int nc,bd,sr,len;
    if (!is_sound_data(g,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    return sr;
}
static const char _samplerate_s []="samplerate";
static define_unary_function_eval (__samplerate,&_samplerate,_samplerate_s);
define_unary_function_ptr5(at_samplerate,alias_at_samplerate,&__samplerate,0,true)

gen _duration(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    int nc,bd,sr,len;
    if (!is_sound_data(g,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    return _evalf(fraction(len,sr),contextptr);
}
static const char _duration_s []="duration";
static define_unary_function_eval (__duration,&_duration,_duration_s);
define_unary_function_ptr5(at_duration,alias_at_duration,&__duration,0,true)

gen _channel_data(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT)
        return gentypeerr(contextptr);
    vecteur opts;
    if (g.subtype==_SEQ__VECT) {
        if (g._VECTptr->front().type!=_VECT)
            return generrtype(_SP_BAD_SOUND_DATA);
        opts=vecteur(g._VECTptr->begin()+1,g._VECTptr->end());
    }
    const gen &data=g.subtype==_SEQ__VECT?g._VECTptr->front():g;
    int nc,bd,sr,len;
    if (!is_sound_data(data,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    int chan=0,slice_start=0,slice_len=len;
    bool asmatrix=false;
    for (const_iterateur it=opts.begin();it!=opts.end();++it) {
        if (it->is_integer()) {
            if (it->val<1 || it->val>nc)
                return generr("Invalid channel number");
            chan=it->val;
        } else if (*it==at_matrix)
            asmatrix=true;
        else if (is_equal(*it)) {
            const gen &lh=it->_SYMBptr->feuille._VECTptr->front();
            const gen &rh=it->_SYMBptr->feuille._VECTptr->back();
            if (lh==at_range) {
                if (rh.type==_VECT) {
                    if (rh._VECTptr->size()!=2)
                        return generrdim(_SP_INVALID_RANGE);
                    if (!is_integer_vecteur(*rh._VECTptr))
                        return generr(_SP_INVALID_RANGE);
                    int start=rh._VECTptr->front().val,stop=rh._VECTptr->back().val;
                    slice_start=start-1; slice_len=stop-start+1;
                } else if (rh.is_symb_of_sommet(at_interval)) {
                    const gen &a=rh._SYMBptr->feuille._VECTptr->front();
                    const gen &b=rh._SYMBptr->feuille._VECTptr->back();
                    if (!a.is_integer() || !b.is_integer()) {
                        if (!is_real_number(a,contextptr) || !is_real_number(b,contextptr))
                            return generr(_SP_INVALID_RANGE);
                        slice_start=std::floor(sr*to_real_number(a,contextptr).to_double(contextptr));
                        slice_len=std::floor(sr*to_real_number(b-a,contextptr).to_double(contextptr));
                    } else {
                        slice_start=a.val-1;
                        slice_len=(b-a).val+1;
                    }
                } else return generr(_SP_INVALID_RANGE);
            }
        }
    }
    if (slice_start<0 || slice_start>=len || slice_len<0 || slice_len>len || slice_start+slice_len>len)
        return generr("Invalid slice bounds");
    if (chan==0) {
        vecteur ret;
        for (const_iterateur it=data._VECTptr->begin()+1;it!=data._VECTptr->end();++it) {
            ret.push_back(decode_chdata(*it->_VECTptr,bd,slice_start,slice_len));
        }
        if (ret.size()==1)
            return ret.front();
        return asmatrix?ret:change_subtype(ret,_SEQ__VECT);
    }
    return decode_chdata(*data._VECTptr->at(chan)._VECTptr,bd,slice_start,slice_len);
}
static const char _channel_data_s []="channel_data";
static define_unary_function_eval (__channel_data,&_channel_data,_channel_data_s);
define_unary_function_ptr5(at_channel_data,alias_at_channel_data,&__channel_data,0,true)

gen _cross_correlation(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &args=*g._VECTptr;
    if (args.size()!=2 || args.front().type!=_VECT || args.back().type!=_VECT)
        return generr("Expected a sequence of two lists");
    vecteur A=*args.front()._VECTptr,B=*args.back()._VECTptr;
    int m=A.size(),n=B.size(),N=nextpow2(std::max(n,m));
    A.resize(2*N,0);
    B.resize(2*N,0);
    vecteur a=*_fft(A,contextptr)._VECTptr,b=*_fft(B,contextptr)._VECTptr;
    vecteur cc_ffted=*_pointprod(makesequence(a,conj(b,contextptr)),contextptr)._VECTptr;
    vecteur cc=*_apply(makesequence(at_real,_ifft(cc_ffted,contextptr)),contextptr)._VECTptr;
    reverse(cc.begin(),cc.begin()+N);
    reverse(cc.begin()+N,cc.end());
    return vecteur(cc.begin()+N-m,cc.end()-N+n-1);
}
static const char _cross_correlation_s []="cross_correlation";
static define_unary_function_eval (__cross_correlation,&_cross_correlation,_cross_correlation_s);
define_unary_function_ptr5(at_cross_correlation,alias_at_cross_correlation,&__cross_correlation,0,true)

gen _auto_correlation(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG) {
        if (g.subtype==-1) return g;
        return _auto_correlation(_readwav(g,contextptr),contextptr);
    }
    return _cross_correlation(makesequence(g,g),contextptr);
}
static const char _auto_correlation_s []="auto_correlation";
static define_unary_function_eval (__auto_correlation,&_auto_correlation,_auto_correlation_s);
define_unary_function_ptr5(at_auto_correlation,alias_at_auto_correlation,&__auto_correlation,0,true)

void lowpass(vecteur &data,double cutoff,int sr) {
    double rc=0.15915494309/cutoff;
    double dt=1.0/sr;
    gen s=data.front(),alpha(dt/(rc+dt));
    iterateur it=data.begin(),itend=data.end();
    for (;it!=itend;++it) {
        s+=alpha*(*it-s);
        *it=s;
    }
}

void highpass(vecteur &data,double cutoff,int sr) {
    double rc=0.15915494309/cutoff;
    double dt=1.0/sr;
    gen alpha(rc/(rc+dt)),prevdata=data.front(),temp;
    iterateur it=data.begin()+1,itend=data.end();
    for (;it!=itend;++it) {
        temp=*it;
        *it=alpha*(*(it-1)+*it-prevdata);
        prevdata=temp;
    }
}

gen filter(const vecteur &args,filter_type typ,GIAC_CONTEXT) {
    double cutoff=to_real_number(args.at(1),contextptr).to_double(contextptr);
    int nc,bd,sr,len;
    if (is_sound_data(*args.front()._VECTptr,nc,bd,sr,len)) {
        if (cutoff<=0 || cutoff>=sr/2)
            return generr("Invalid cutoff specification");
        gen opt(undef);
        if (args.size()>2) {
            if (!is_equal(args[2]))
                return generr("Third argument must be normalize=<real>");
            if (args[2]._SYMBptr->feuille._VECTptr->front()==at_normalize)
                opt=args[2];
        }
        matrice data;
        for (int i=1;i<=nc;++i) {
            vecteur chdata=decode_chdata(*args.front()._VECTptr->at(i)._VECTptr,bd);
            switch(typ) {
            case _LOWPASS_FILTER: lowpass(chdata,cutoff,sr); break;
            case _HIGHPASS_FILTER: highpass(chdata,cutoff,sr); break;
            }
            data.push_back(chdata);
        }
        gen argv=makesequence(symbolic(at_equal,makesequence(at_samplerate,sr)),
                              symbolic(at_equal,makesequence(at_bit_depth,bd)),data);
        if (!is_undef(opt)) argv._VECTptr->push_back(opt);
        return _createwav(argv,contextptr);
    } else {
        vecteur data=*args.front()._VECTptr;
        sr=44100;
        if (args.size()>2) {
            if (!args.at(2).is_integer() || (sr=args.at(2).val)<=0)
                return generr("Invalid sample rate specification");
        }
        if (cutoff<=0 || cutoff>=sr/2)
            return generr("Invalid cutoff specification");
        switch(typ) {
        case _LOWPASS_FILTER: lowpass(data,cutoff,sr); break;
        case _HIGHPASS_FILTER: highpass(data,cutoff,sr); break;
        }
        return data;
    }
}

gen _lowpass(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()<2 ||
            g._VECTptr->front().type!=_VECT || !is_real_number(g._VECTptr->at(1),contextptr))
        return gentypeerr(contextptr);
    return filter(*g._VECTptr,_LOWPASS_FILTER,contextptr);
}
static const char _lowpass_s []="lowpass";
static define_unary_function_eval (__lowpass,&_lowpass,_lowpass_s);
define_unary_function_ptr5(at_lowpass,alias_at_lowpass,&__lowpass,0,true)

gen _highpass(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()<2 ||
            g._VECTptr->front().type!=_VECT || !is_real_number(g._VECTptr->at(1),contextptr))
        return gentypeerr(contextptr);
    return filter(*g._VECTptr,_HIGHPASS_FILTER,contextptr);
}
static const char _highpass_s []="highpass";
static define_unary_function_eval (__highpass,&_highpass,_highpass_s);
define_unary_function_ptr5(at_highpass,alias_at_highpass,&__highpass,0,true)

gen _moving_average(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    if (gv.size()!=2)
        return generr("Expected two input arguments");
    if (gv.front().type!=_VECT)
        return generr("First argument must be an array");
    if (!gv.back().is_integer() || gv.back().val<=0)
        return generr("Second argument must be a positive integer");
    const vecteur &s=*gv.front()._VECTptr;
    int n=gv.back().val,len=s.size();
    if (n>len)
        return generr("Filter length exceeds array size");
    vecteur res(len-n+1);
    gen acc(0);
    for (int i=0;i<n;++i) acc+=s[i];
    res[0]=acc;
    for (int i=n;i<len;++i) {
        acc-=s[i-n];
        acc+=s[i];
        res[i-n+1]=acc;
    }
    return multvecteur(fraction(1,n),res);
}
static const char _moving_average_s []="moving_average";
static define_unary_function_eval (__moving_average,&_moving_average,_moving_average_s);
define_unary_function_ptr5(at_moving_average,alias_at_moving_average,&__moving_average,0,true)

gen _rms(const gen &g,GIAC_CONTEXT) {
	if (g.type==_STRNG && g.subtype==-1) return g;
	if (g.type!=_VECT)
		return gentypeerr(contextptr);
	const vecteur &gv=*g._VECTptr;
	int n=gv.size();
	if (n==0)
		return generr("Input list is empty");
	gen res(0),rp,ip;
    const_iterateur it=gv.begin(),itend=gv.end();
	for (;it!=itend;++it) {
		reim(*it,rp,ip,contextptr);
		res+=sq(rp)+sq(ip);
	}
	res=res/n;
	return sqrt(res,contextptr);
}
static const char _rms_s []="rms";
static define_unary_function_eval (__rms,&_rms,_rms_s);
define_unary_function_ptr5(at_rms,alias_at_rms,&__rms,0,true)

gen _resample(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
#ifndef HAVE_LIBSAMPLERATE
    print_error("libsamplerate is required for resampling audio",contextptr);
    return vecteur(0);
#else
    if (g.type!=_VECT)
        return gentypeerr(contextptr);
    int nc,bd,sr,len;
    const gen &snd=g.subtype==_SEQ__VECT?g._VECTptr->front():g;
    if (!is_sound_data(snd,nc,bd,sr,len))
        return generrtype(_SP_BAD_SOUND_DATA);
    int nsr=44100;
    int quality=2;
    if (g.subtype==_SEQ__VECT) {
        const vecteur &gv=*g._VECTptr;
        if (gv.size()<2 || !gv[1].is_integer() || (nsr=gv[1].val)<1)
            return generr("Second argument must be a positive integer");
        if (gv.size()>2) {
            if (!gv[2].is_integer() || (quality=gv[2].val)<0 || (quality>4))
                return generr("Invalid quality specification (expected an integer 1-4)");
        }
    }
    matrice chdata;
    for (int i=1;i<=nc;++i) {
        vecteur ch=decode_chdata(*snd._VECTptr->at(i)._VECTptr,bd);
        chdata.push_back(ch);
    }
    SRC_DATA data;
    data.src_ratio=double(nsr)/double(sr);
    int nlen=(int)std::floor(len*data.src_ratio);
    data.input_frames=len;
    data.output_frames=nlen;
    float *indata=new float[len*nc];
    data.data_out=new float[nlen*nc];
    for (int i=0;i<len;++i) {
        for (int j=0;j<nc;++j) {
            indata[i*nc+j]=to_real_number(chdata[j][i],contextptr).to_double(contextptr);
        }
    }
    data.data_in=indata;
    switch(quality) {
    case 0:
        quality=SRC_LINEAR;
        break;
    case 1:
        quality=SRC_ZERO_ORDER_HOLD;
        break;
    case 2:
        quality=SRC_SINC_FASTEST;
        break;
    case 3:
        quality=SRC_SINC_MEDIUM_QUALITY;
        break;
    case 4:
        quality=SRC_SINC_BEST_QUALITY;
        break;
    default:
        assert(false);
    }
    src_simple(&data,quality,nc);
    matrice output=*_matrix(makesequence(nc,nlen,0),contextptr)._VECTptr;
    for (int i=0;i<nlen;++i) {
        for (int j=0;j<nc;++j) {
            output[j]._VECTptr->at(i)=data.data_out[i*nc+j];
        }
    }
    delete[] indata;
    delete[] data.data_out;
    vecteur header=*snd._VECTptr->front()._VECTptr;
    header[2]=nsr;
    header[3]=(bd*nc*nlen)/8;
    vecteur ret;
    ret.push_back(header);
    for (const_iterateur it=output.begin();it!=output.end();++it) {
        ret.push_back(encode_chdata(*it->_VECTptr,bd,1.0,contextptr));
    }
    return ret;
#endif
}
static const char _resample_s []="resample";
static define_unary_function_eval (__resample,&_resample,_resample_s);
define_unary_function_ptr5(at_resample,alias_at_resample,&__resample,0,true)

int varcount=0;

gen _convolution(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &args=*g._VECTptr;
    if (!args.empty() && args.front().type!=_VECT) {
        // convolve real functions
        gen T(0),var=identificateur("x"),tvar=identificateur(" tau"+print_INT_(++varcount));
        int n=args.size(),optstart=2;
        if (n<2)
            return generrdim("Too few input arguments");
        const gen &f1=args[0],&f2=args[1];
        if (n>2) {
            if (args[optstart].type==_IDNT)
                var=args[optstart++];
            // parse options
            for (const_iterateur it=args.begin()+optstart;it!=args.end();++it) {
                if (!is_equal(*it))
                    return generr("Expected option=value");
                vecteur &s=*it->_SYMBptr->feuille._VECTptr;
                if (s.front()==at_shift)
                    T=s.back();
                else return generr("Unrecognized option");
            }
        }
        giac_assume(symb_superieur_egal(tvar,0),contextptr);
        gen c=_integrate(makesequence(f1*_Heaviside(var,contextptr)*
                                      subst(f2*_Heaviside(var,contextptr),var,tvar-var,false,contextptr),
                                      var,minus_inf,plus_inf),contextptr);
        _purge(tvar,contextptr);
        if (contains(_lname(c,contextptr),var))
            return generr("Failed to integrate");
        c=subst(c,tvar,var-T,false,contextptr)*_Heaviside(var-T,contextptr);
        return c;
    }
    // convolve sequences
    if (args.size()!=2 || args.front().type!=_VECT || args.back().type!=_VECT)
        return generr("Expected a pair of lists");
    vecteur A=*args.front()._VECTptr,B=*args.back()._VECTptr;
    int lenA=A.size(),lenB=B.size(),len=2*nextpow2(std::max(lenA,lenB));
    A.resize(len-1,0);
    B.resize(len-1,0);
    vecteur a=*_fft(A,contextptr)._VECTptr,b=*_fft(B,contextptr)._VECTptr;
    vecteur cv=*_apply(makesequence(at_real,_ifft(_pointprod(makesequence(a,b),
                                                             contextptr),
                                                  contextptr)),
                       contextptr)._VECTptr;
    cv.resize(lenA+lenB-1);
    return cv;
}
static const char _convolution_s []="convolution";
static define_unary_function_eval (__convolution,&_convolution,_convolution_s);
define_unary_function_ptr5(at_convolution,alias_at_convolution,&__convolution,0,true)

gen _boxcar(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()!=3)
        return gentypeerr(contextptr);
    const vecteur &args=*g._VECTptr;
    return _Heaviside(args[2]-args[0],contextptr)-_Heaviside(args[2]-args[1],contextptr);
}
static const char _boxcar_s []="boxcar";
static define_unary_function_eval (__boxcar,&_boxcar,_boxcar_s);
define_unary_function_ptr5(at_boxcar,alias_at_boxcar,&__boxcar,0,true)

gen _rect(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type==_VECT && g.subtype==_SEQ__VECT && g._VECTptr->size()==2)
      return g._VECTptr->front()*exp(g._VECTptr->back()*cst_i,contextptr);
    gen hf=fraction(1,2);
    return _boxcar(makesequence(-hf,hf,g),contextptr);
}
static const char _rect_s []="rect";
static define_unary_function_eval (__rect,&_rect,_rect_s);
define_unary_function_ptr5(at_rect,alias_at_rect,&__rect,0,true)

gen _tri(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    return (1-g)*(_Heaviside(g,contextptr)-_Heaviside(g-1,contextptr))+
           (1+g)*(_Heaviside(-g,contextptr)-_Heaviside(-g-1,contextptr));
}
static const char _tri_s []="tri";
static define_unary_function_eval (__tri,&_tri,_tri_s);
define_unary_function_ptr5(at_tri,alias_at_tri,&__tri,0,true)

gen _sinc(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (is_zero(g))
        return gen(1);
    return sin(g,contextptr)/g;
}
static const char _sinc_s []="sinc";
static define_unary_function_eval (__sinc,&_sinc,_sinc_s);
define_unary_function_ptr5(at_sinc,alias_at_sinc,&__sinc,0,true)

gen _logistic(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    gen x(g),x0(0),L(1),k(1); // standard sigmoid by default
    if (g.type==_VECT && g.subtype==_SEQ__VECT) {
        vecteur gv(*g._VECTptr);
        if (gv.empty() || gv.size()>4)
            return gendimerr(contextptr);
        x=gv.front();
        if (gv.size()>=2)
            x0=gv[1];
        if (gv.size()>=3)
            L=gv[2];
        if (gv.size()==4)
            k=gv[3];
    }
    return L/(exp(-k*(x-x0),contextptr)+1);
}
static const char _logistic_s []="logistic";
static define_unary_function_eval (__logistic,&_logistic,_logistic_s);
define_unary_function_ptr5(at_logistic,alias_at_logistic,&__logistic,0,true)

#define MAX_TAILLE 500

/* return true iff g is significantly simpler than h */
bool is_simpler(const gen &g,const gen &h,double scale=1.0) {
    return scale*taille(h,MAX_TAILLE)>taille(g,MAX_TAILLE);
}

bool is_rational_wrt(const gen &e,const identificateur &x) {
    return rlvarx(e,x).size()<=1;
}

gen factorise(const gen &g,GIAC_CONTEXT) {
    return factorcollect(g,false,contextptr);
}

bool ispoly(const gen &e,const identificateur &x,gen &d,GIAC_CONTEXT) {
    if (is_constant_wrt(e,x,contextptr)) {
        d=gen(0);
        return true;
    }
    if (!is_rational_wrt(e,x) || !is_constant_wrt(_denom(e,contextptr),x,contextptr))
        return false;
    d=_degree(makesequence(e,x),contextptr);
    return d.is_integer() && !is_zero(d);
}

void constlin_terms(const gen &g,const identificateur &x,gen &lt,gen &c,gen &rest,GIAC_CONTEXT) {
    gen e=expand(g,contextptr);
    vecteur terms=(e.is_symb_of_sommet(at_plus) && e._SYMBptr->feuille.type==_VECT?*e._SYMBptr->feuille._VECTptr:vecteur(1,e));
    gen a,b;
    rest=c=lt=gen(0);
    for (const_iterateur it=terms.begin();it!=terms.end();++it) {
        if (is_constant_wrt(*it,x,contextptr))
            c+=*it;
        else if (is_linear_wrt(*it,x,a,b,contextptr)) {
            lt+=a;
            c+=b;
        } else rest+=*it;
    }
}

gen logabs_expand(const gen &g,const identificateur &x,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        vecteur res;
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            res.push_back(logabs_expand(*it,x,contextptr));
        }
        return gen(res,g.subtype);
    }
    gen a,b;
    if (g.is_symb_of_sommet(at_ln) && g._SYMBptr->feuille.is_symb_of_sommet(at_abs) &&
            is_linear_wrt(g._SYMBptr->feuille._SYMBptr->feuille,x,a,b,contextptr) &&
            !is_zero(a) && !is_one(a) && is_zero(im(a,contextptr)) && is_zero(im(b,contextptr)))
        return ln(_abs(x+b/a,contextptr),contextptr)+ln(_abs(a,contextptr),contextptr);
    if (g.type==_SYMB)
        return symbolic(g._SYMBptr->sommet,logabs_expand(g._SYMBptr->feuille,x,contextptr));
    return g;
}

vecteur analyze_terms(const gen &g,const identificateur &x,GIAC_CONTEXT) {
    vecteur terms,ret;
    if (g.is_symb_of_sommet(at_plus)) {
        if (g._SYMBptr->feuille.type==_VECT)
            terms=*g._SYMBptr->feuille._VECTptr;
        else terms=vecteur(1,g._SYMBptr->feuille);
    } else terms=vecteur(1,g);
    for (const_iterateur it=terms.begin();it!=terms.end();++it) {
        vecteur factors;
        gen cnst(1),rest(1),deg(0),exprest(0),sh(0),rt=ratnormal(*it,contextptr),st;
        if (!is_constant_wrt(rt,x,contextptr))
            rt=factorise(rt,contextptr);
        if (rt.is_symb_of_sommet(at_neg)) {
            cnst=gen(-1);
            rt=rt._SYMBptr->feuille;
        }
        if (rt.is_symb_of_sommet(at_prod) && rt._SYMBptr->feuille.type==_VECT)
            factors=*rt._SYMBptr->feuille._VECTptr;
        else factors=vecteur(1,rt);
        for (int i=factors.size();i-->0;) {
            gen &fac=factors[i];
            if (fac.is_symb_of_sommet(at_inv)) {
                if (fac._SYMBptr->feuille.is_symb_of_sommet(at_prod) &&
                    fac._SYMBptr->feuille._SYMBptr->feuille.type==_VECT) {
                    const vecteur &d=*fac._SYMBptr->feuille._SYMBptr->feuille._VECTptr;
                    for (const_iterateur jt=d.begin();jt!=d.end();++jt) {
                        factors.push_back(symbolic(at_inv,*jt));
                    }
                    factors.erase(factors.begin()+i);
                }
            }
        }
        for (const_iterateur jt=factors.begin();jt!=factors.end();++jt) {
            if (jt->is_symb_of_sommet(at_exp)) {
                gen lt1,lt2,c1,c2,rest1,rest2;
                constlin_terms(re(jt->_SYMBptr->feuille,contextptr),x,lt1,c1,rest1,contextptr);
                constlin_terms(im(jt->_SYMBptr->feuille,contextptr),x,lt2,c2,rest2,contextptr);
                cnst=cnst*exp(c1+cst_i*c2,contextptr);
                sh+=lt2;
                exprest+=rest1+x*lt1+cst_i*rest2;
            } else if (is_constant_wrt(*jt,x,contextptr))
                cnst=cnst*(*jt);
            else if (ispoly(*jt,x,st,contextptr)) {
                deg+=st;
                cnst=cnst*_lcoeff(makesequence(*jt,x),contextptr);
            } else rest=rest*(*jt);
        }
        rest=rest*exp(exprest,contextptr);
        cnst=recursive_normal(_lin(cnst,contextptr),contextptr);
        ret.push_back(makevecteur(cnst,sh,deg,rest));
    }
    return ret;
}

bool has_rootof(const gen &g) {
    if (g.is_symb_of_sommet(at_rootof)) return true;
    if (g.type==_VECT) {
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (has_rootof(*it)) return true;
        }
    } else if (g.type==_SYMB) return has_rootof(g._SYMBptr->feuille);
    return false;
}

bool has_partial_diff(const gen &g) {
    if (g.is_symb_of_sommet(at_of) &&
            g._SYMBptr->feuille._VECTptr->front().is_symb_of_sommet(at_derive))
        return true;
    if (g.type==_VECT) {
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (has_partial_diff(*it)) return true;
        }
    } else if (g.type==_SYMB) return has_partial_diff(g._SYMBptr->feuille);
    return false;
}

bool is_heavisided(const gen &g,const identificateur &x,gen &rest,gen &a,gen &b,GIAC_CONTEXT) {
    if (g.is_symb_of_sommet(at_neg)) {
        if (!is_heavisided(g._SYMBptr->feuille,x,rest,a,b,contextptr))
            return false;
        rest=-rest;
        return true;
    }
    if (!g.is_symb_of_sommet(at_prod) || g._SYMBptr->feuille.type!=_VECT)
        return false;
    rest=gen(1);
    const vecteur &factors=*g._SYMBptr->feuille._VECTptr;
    bool yes=false;
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (!yes && it->is_symb_of_sommet(at_Heaviside) &&
                is_linear_wrt(it->_SYMBptr->feuille,x,a,b,contextptr) && !is_zero(a))
            yes=true;
        else rest=rest*(*it);
    }
    return yes;
}

gen lin_cplx(const gen &g_orig,GIAC_CONTEXT) {
    gen g=_lin(g_orig,contextptr),ret(1);
    if (g.is_symb_of_sommet(at_neg)) {
        ret=gen(-1);
        g=g._SYMBptr->feuille;
    }
    vecteur factors;
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT)
        factors=*g._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,g);
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (it->is_symb_of_sommet(at_exp)) {
            gen ex=it->_SYMBptr->feuille;
            ret=ret*exp(cst_i*im(ex,contextptr),contextptr)*exp(re(ex,contextptr),contextptr);
        } else ret=ret*(*it);
    }
    return ret;
}

gen hcollect(const gen &g,const gen &h,const identificateur &x,gen &cnst,GIAC_CONTEXT) {
    gen f=factorise(lin_cplx(g,contextptr)+lin_cplx(h,contextptr),contextptr);
    gen rest(1),r1,r2,a1,a2,b1,b2,hvs(undef);
    vecteur factors;
    cnst=gen(1);
    if (f.is_symb_of_sommet(at_neg)) {
        rest=gen(-1);
        f=f._SYMBptr->feuille;
    }
    if (f.is_symb_of_sommet(at_prod) && f._SYMBptr->feuille.type==_VECT)
        factors=*f._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,f);
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (is_undef(hvs) && it->is_symb_of_sommet(at_plus) && it->_SYMBptr->feuille.type==_VECT &&
                it->_SYMBptr->feuille._VECTptr->size()==2 &&
                is_heavisided(it->_SYMBptr->feuille._VECTptr->front(),x,r1,a1,b1,contextptr) &&
                is_heavisided(it->_SYMBptr->feuille._VECTptr->back(),x,r2,a2,b2,contextptr) &&
                is_zero(ratnormal(b1*a2-b2*a1,contextptr)) && is_strictly_positive(-a1*a2,contextptr) &&
                is_zero(recursive_normal(subst(r1,x,-b1/a1+_sign(a1,contextptr)*x,false,contextptr)-
                                         subst(r2,x,-b2/a2+_sign(a2,contextptr)*x,false,contextptr),
                                         contextptr))) {
            hvs=subst(r1,x,-b1/a1+_sign(a1,contextptr)*_abs(x+b1/a1,contextptr),false,contextptr);
        } else if (is_constant_wrt(*it,x,contextptr))
            cnst=cnst*(*it);
        else rest=rest*(*it);
    }
    return rest*hvs;
}

gen Heaviside2abs(const gen &g,const identificateur &x,GIAC_CONTEXT) {
    gen h=expand(g,contextptr),e,c;
    vecteur terms;
    if (h.is_symb_of_sommet(at_plus) && h._SYMBptr->feuille.type==_VECT)
        terms=*h._SYMBptr->feuille._VECTptr;
    else terms=vecteur(1,h);
    for (int i=terms.size();i-->0;) {
        for (int j=terms.size();j-->i+1;) {
            if (!is_undef(e=hcollect(terms[i],terms[j],x,c,contextptr))) {
                terms.erase(terms.begin()+j);
                terms.erase(terms.begin()+i);
                terms.push_back(c*e);
                break;
            }
        }
    }
    return symbolic(at_plus,gen(terms,_SEQ__VECT));
}

bool has_integral(const gen &g) {
    if (g.type==_VECT) {
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            if (has_integral(*it))
                return true;
        }
        return false;
    }
    if (g.type==_SYMB)
        return has_integral(g._SYMBptr->feuille);
    if (g.is_symb_of_sommet(at_integrate) || g.is_symb_of_sommet(at_int) ||
            g.is_symb_of_sommet(at_Int))
        return true;
    return false;
}

bool is_same(const gen &g,const gen &h,GIAC_CONTEXT) {
    return is_zero(ratnormal(g-h,contextptr));
}

vecteur fourier_terms(const gen &g_orig,const identificateur &x,bool do_simp,GIAC_CONTEXT) {
    gen g;
    if (do_simp)
        g=expand(ratnormal(logabs_expand(lnexpand(_lin(_trig2exp(_hyp2exp(_pow2exp(g_orig,contextptr),
                 contextptr),contextptr),contextptr),contextptr),x,contextptr),contextptr),contextptr);
    else g=g_orig;
    vecteur terms=analyze_terms(g,x,contextptr);
    // decompose terms into sums of partial fractions
    for (int i=terms.size();i-->0;) {
        assert(terms[i].type==_VECT);
        const vecteur &v=*terms[i]._VECTptr;
        if (!is_rational_wrt(v[3],x) ||
            is_constant_wrt(_denom(v[3],contextptr),x,contextptr)) continue;
        vecteur rterms=analyze_terms(_cpartfrac(makesequence(v[3],x),contextptr),x,contextptr);
        for (const_iterateur it=rterms.begin();it!=rterms.end();++it) {
            assert(it->type==_VECT);
            const vecteur &w=*(it->_VECTptr);
            terms.push_back(makevecteur(v[0]*w[0],v[1],v[2]+w[2],expand(w[3],contextptr)));
        }
        terms.erase(terms.begin()+i);
    }
    // collect terms to reduce complexity
    for (int i=0;i<int(terms.size());++i) {
        vecteur &v=*terms[i]._VECTptr;
        for (int j=terms.size();j-->i+1;) {
            const vecteur &w=*terms[j]._VECTptr;
            if (is_same(v[3],w[3],contextptr) && is_same(v[1],w[1],contextptr) &&
                    is_same(v[2],w[2],contextptr)) {
                v[0]+=w[0];
                terms.erase(terms.begin()+j);
            }
        }
    }
    // convert: f(-x)*Heaviside(x)+f(x)*Heaviside(-x) -> f(-abs(x))
    for (int i=terms.size();i-->0;) {
        vecteur &v=*terms[i]._VECTptr;
        for (int j=terms.size();j-->i+1;) {
            const vecteur &w=*terms[j]._VECTptr;
            gen e,c;
            if (is_same(v[1],w[1],contextptr) && is_zero(v[2]-w[2]) &&
                    !is_undef(e=hcollect(v[3]*v[0],w[3]*w[0],x,c,contextptr))) {
                e=makevecteur(c,v[1],v[2],e);
                terms.erase(terms.begin()+j);
                terms.erase(terms.begin()+i);
                terms.push_back(e);
                break;
            }
        }
    }
    return terms;
}

static vecteur fourier_table;
static vecteur laplace_table;

gen _addtable(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &args=*g._VECTptr;
    if (args.size()!=5)
        return generrdim("Wrong number of input arguments");
    int ti=0;
    if (args[0]==at_fourier) ti=1;
    if (args[0]==at_laplace) ti=2;
    if (!args[1].is_symb_of_sommet(at_of) || !args[2].is_symb_of_sommet(at_of) ||
            args[1]._SYMBptr->feuille._VECTptr->front().type!=_IDNT ||
            args[2]._SYMBptr->feuille._VECTptr->front().type!=_IDNT ||
            args[3].type!=_IDNT || args[4].type!=_IDNT)
        return generrtype("Invalid input argument");
    const identificateur &f=*args[1]._SYMBptr->feuille._VECTptr->front()._IDNTptr;
    const identificateur &F=*args[2]._SYMBptr->feuille._VECTptr->front()._IDNTptr;
    const identificateur &t=*args[3]._IDNTptr,&s=*args[4]._IDNTptr;
    const vecteur fvars=gen2vecteur(args[1]._SYMBptr->feuille._VECTptr->back());
    const vecteur Fvars=gen2vecteur(args[2]._SYMBptr->feuille._VECTptr->back());
    if (fvars.size()!=Fvars.size())
        return gendimerr(contextptr);
    int k=-1;
    for (int i=fvars.size();i-->0;) {
        if (fvars[i].type!=_IDNT || Fvars[i].type!=_IDNT)
            return generr("Expected an identifier");
        const identificateur &x1=*fvars[i]._IDNTptr,&x2=*Fvars[i]._IDNTptr;
        if (x1==x2) continue;
        if (k>=0 || !(x1==t) || !(x2==s))
            return gensizeerr(contextptr);
        k=i;
    }
    if (k<0)
        return gensizeerr(contextptr);
    if (ti<1 || ti>2) return gen(0);
    vecteur &tbl=(ti==1?fourier_table:laplace_table);
    for (const_iterateur it=tbl.begin();it!=tbl.end();++it) {
        const vecteur &v=*(it->_VECTptr);
        if (v[0]==f || v[0]==F || v[1]==f || v[1]==F) return gen(0);
    }
    tbl.push_back(makevecteur(f,F,fvars.size(),k));
    return gen(1);
}
static const char _addtable_s []="addtable";
static define_unary_function_eval (__addtable,&_addtable,_addtable_s);
define_unary_function_ptr5(at_addtable,alias_at_addtable,&__addtable,0,true)

vecteur fourier_table_find(const identificateur &f,int sz) {
    for (const_iterateur it=fourier_table.begin();it!=fourier_table.end();++it) {
        const vecteur &r=*(it->_VECTptr);
        if (r[2].val!=sz) continue;
        if (r[0]==f) return makevecteur(1,r[1],r[3]);
        if (r[1]==f) return makevecteur(2*cst_pi,r[0],r[3]);
    }
    return vecteur(0);
}

bool is_func(const gen &g,identificateur &f,vecteur &vars) {
    if (!g.is_symb_of_sommet(at_of) || g._SYMBptr->feuille._VECTptr->front().type!=_IDNT)
        return false;
    f=*g._SYMBptr->feuille._VECTptr->front()._IDNTptr;
    vars=gen2vecteur(g._SYMBptr->feuille._VECTptr->back());
    return true;
}

bool is_diff(const gen &g,const identificateur &x,gen &f,int &d) {
    if (!g.is_symb_of_sommet(at_derive) || g._SYMBptr->feuille.type!=_VECT)
        return false;
    const vecteur &args=*g._SYMBptr->feuille._VECTptr;
    int n;
    if (args.size()==2 && args.back()==x) d=1;
    else if (args.size()==3 && args[1]==x && args[2].is_integer() && (n=args[2].val)>0) d=n;
    else return false;
    f=args.front();
    return true;
}

bool is_partialdiff(const gen &g,identificateur &f,vecteur &vars,vecteur &deg,GIAC_CONTEXT) {
    if (!g.is_symb_of_sommet(at_of)) return false;
    gen df=g._SYMBptr->feuille._VECTptr->front();
    const gen &v=g._SYMBptr->feuille._VECTptr->back();
    if (v.type!=_VECT) return false;
    vars=*v._VECTptr;
    int n,sz=vars.size();
    deg=vecteur(sz,0);
    while (df.is_symb_of_sommet(at_derive) && df._SYMBptr->feuille.type==_VECT) {
        const vecteur &args=*df._SYMBptr->feuille._VECTptr;
        if (args.size()!=2 || !args.back().is_integer() || (n=args.back().val)<0 || n>=sz)
            return false;
        deg[n]+=gen(1);
        df=args.front();
    }
    if (df.type!=_IDNT || is_zero(deg,contextptr)) return false;
    f=*df._IDNTptr;
    return true;
}

gen lin_Heaviside(const gen &g_orig,const identificateur &x,GIAC_CONTEXT) {
    gen g=expand(g_orig,contextptr),rest(1),ret(0),a,b;
    if (g.is_symb_of_sommet(at_plus) && g._SYMBptr->feuille.type==_VECT) {
        const vecteur &terms=*g._SYMBptr->feuille._VECTptr;
        for (const_iterateur it=terms.begin();it!=terms.end();++it) {
            ret+=lin_Heaviside(*it,x,contextptr);
        }
        return ret;
    }
    if (g.is_symb_of_sommet(at_neg)) {
        rest=gen(-1);
        g=g._SYMBptr->feuille;
    }
    vecteur factors;
    matrice L;
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT)
        factors=*g._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,g);
    int i=0;
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (i<2 && it->is_symb_of_sommet(at_Heaviside) &&
                is_linear_wrt(it->_SYMBptr->feuille,x,a,b,contextptr)) {
            L.push_back(makevecteur(_sign(a,contextptr),-b/_abs(a,contextptr)));
            ++i;
        } else rest=rest*(*it);
    }
    if (i<2) {
        if (i==1) rest=rest*_Heaviside(L[0][0]*x-L[0][1],contextptr);
        return rest;
    }
    if (is_one(L[0][0]*L[1][0])) {
        if (is_one(L[0][0]))
            rest=rest*_Heaviside(x-(is_greater(L[1][1],L[0][1],contextptr)?L[1][1]:L[0][1]),contextptr);
        else
            rest=rest*_Heaviside((is_greater(L[1][1],L[0][1],contextptr)?L[0][1]:L[1][1])-x,contextptr);
    } else {
        a=is_one(L[0][0])?L[0][1]:L[1][1];
        b=is_one(L[0][0])?L[1][1]:L[0][1];
        if (is_greater(b,a,contextptr))
            rest=rest*(_Heaviside(x-a,contextptr)-_Heaviside(x-b,contextptr));
        else rest=gen(0);
    }
    return lin_Heaviside(rest,x,contextptr);
}

gen cond2Heaviside(const gen &g,const identificateur &x,GIAC_CONTEXT) {
    if ((g.is_symb_of_sommet(at_and) || g.is_symb_of_sommet(at_et)) && g._SYMBptr->feuille.type==_VECT) {
        gen res(1);
        const vecteur &args=*g._SYMBptr->feuille._VECTptr;
        for (const_iterateur it=args.begin();it!=args.end();++it) {
            res=res*cond2Heaviside(*it,x,contextptr);
        }
        return res;
    }
    if (g.is_symb_of_sommet(at_ou) && g._SYMBptr->feuille.type==_VECT) {
        gen res(1);
        const vecteur &args=*g._SYMBptr->feuille._VECTptr;
        for (const_iterateur it=args.begin();it!=args.end();++it) {
            res=res*(1-cond2Heaviside(*it,x,contextptr));
        }
        return 1-res;
    }
    if (g.is_symb_of_sommet(at_inferieur_egal) || g.is_symb_of_sommet(at_inferieur_strict)) {
        if (g._SYMBptr->feuille._VECTptr->front()==x &&
                is_constant_wrt(g._SYMBptr->feuille._VECTptr->back(),x,contextptr))
            return _Heaviside(g._SYMBptr->feuille._VECTptr->back()-x,contextptr);
        if (g._SYMBptr->feuille._VECTptr->back()==x &&
                is_constant_wrt(g._SYMBptr->feuille._VECTptr->front(),x,contextptr))
            return _Heaviside(x-g._SYMBptr->feuille._VECTptr->front(),contextptr);
    }
    if (g.is_symb_of_sommet(at_superieur_egal) || g.is_symb_of_sommet(at_superieur_strict)) {
        if (g._SYMBptr->feuille._VECTptr->front()==x &&
                is_constant_wrt(g._SYMBptr->feuille._VECTptr->back(),x,contextptr))
            return _Heaviside(x-g._SYMBptr->feuille._VECTptr->back(),contextptr);
        if (g._SYMBptr->feuille._VECTptr->back()==x &&
                is_constant_wrt(g._SYMBptr->feuille._VECTptr->front(),x,contextptr))
            return _Heaviside(g._SYMBptr->feuille._VECTptr->front()-x,contextptr);
    }
    return undef;
}

gen sign2Heaviside(const gen &g,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        vecteur res;
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            res.push_back(sign2Heaviside(*it,contextptr));
        }
        return gen(res,g.subtype);
    }
    if (g.is_symb_of_sommet(at_sign))
        return 2*symbolic(at_Heaviside,g._SYMBptr->feuille)-1;
    if (g.type==_SYMB)
        return symbolic(g._SYMBptr->sommet,sign2Heaviside(g._SYMBptr->feuille,contextptr));
    return g;
}

gen Heaviside2sign(const gen &g,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        vecteur res;
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            res.push_back(Heaviside2sign(*it,contextptr));
        }
        return gen(res,g.subtype);
    }
    if (g.is_symb_of_sommet(at_Heaviside))
        return (symbolic(at_sign,g._SYMBptr->feuille)+1)/2;
    if (g.type==_SYMB)
        return symbolic(g._SYMBptr->sommet,Heaviside2sign(g._SYMBptr->feuille,contextptr));
    return g;
}

gen simplify_signs(const gen &g,const identificateur &x,GIAC_CONTEXT) {
    if (g.type==_VECT) {
        vecteur res;
        for (const_iterateur it=g._VECTptr->begin();it!=g._VECTptr->end();++it) {
            res.push_back(simplify_signs(*it,x,contextptr));
        }
        return gen(res,g.subtype);
    }
    gen a,b;
    if ((g.is_symb_of_sommet(at_Heaviside) || g.is_symb_of_sommet(at_sign)) &&
            is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr) && !is_zero(a))
        return symbolic(g._SYMBptr->sommet,g._SYMBptr->feuille/_abs(a,contextptr));
    if (g.is_symb_of_sommet(at_Dirac) && g._SYMBptr->feuille.type!=_VECT &&
            is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr) && !is_zero(a))
        return symbolic(at_Dirac,x+b/a)/_abs(a,contextptr);
    if (g.type==_SYMB)
        return symbolic(g._SYMBptr->sommet,simplify_signs(g._SYMBptr->feuille,x,contextptr));
    return g;
}

vecteur find_Dirac(const gen &g,const identificateur &x,gen &rest,GIAC_CONTEXT) {
    gen a,b,fac;
    rest=gen(1);
    vecteur res,factors;
    if (g.is_symb_of_sommet(at_neg)) {
        res=find_Dirac(g._SYMBptr->feuille,x,rest,contextptr);
        rest=-rest;
        return res;
    }
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT)
        factors=*g._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,g);
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (res.empty() && it->is_symb_of_sommet(at_Dirac)) {
            const gen &arg=(it->_SYMBptr->feuille.type==_VECT?it->_SYMBptr->feuille._VECTptr->front()
                                                             :it->_SYMBptr->feuille);
            if (is_linear_wrt(arg,x,a,b,contextptr) && !is_zero(a)) {
                res=makevecteur(-b/a,it->_SYMBptr->feuille.type==_VECT?
                                     it->_SYMBptr->feuille._VECTptr->back():gen(0));
                continue;
            }
        }
        fac=*it;
        if (fac.is_symb_of_sommet(at_abs))
            fac=fac._SYMBptr->feuille*(2*_Heaviside(fac,contextptr)-1);
        rest=rest*fac;
    }
    rest=ratnormal(rest,contextptr);
    return res;
}

/* expand f(x)*Dirac(x-a,n) to
   sum((-1)^k*binom(n,k)*(d^k f/dx^k)(a)*(d^(n-k) Dirac/dx^(n-k))(x-a),k=0..n) */
gen expand_Dirac_prod(const gen &f_orig,const identificateur &x,const gen &a,int n,GIAC_CONTEXT) {
    gen ret(0),f=f_orig;
    for (int i=0;i<=n;++i) {
        ret+=subst(f,x,a,false,contextptr)*pow(gen(-1),i)*_comb(makesequence(n,i),contextptr)*
             (i==n?_Dirac(x-a,contextptr):_Dirac(makesequence(x-a,n-i),contextptr));
        if (i<n) f=_derive(makesequence(f,x),contextptr);
    }
    return ret;
}

gen simplify_Dirac(const gen &g_orig,const identificateur &x,int depth,GIAC_CONTEXT) {
	if (depth>100)
		return g_orig;
    gen g=expand(g_orig,contextptr),r,ret(0),rest(0);
    vecteur terms;
    if (g.is_symb_of_sommet(at_plus) && g._SYMBptr->feuille.type==_VECT)
        terms=*g._SYMBptr->feuille._VECTptr;
    else terms=vecteur(1,g);
    for (const_iterateur it=terms.begin();it!=terms.end();++it) {
        vecteur d=find_Dirac(*it,x,r,contextptr);
        if (d.empty()) {
            ret+=*it;
        } else {
            assert(d.size()==2 && d.back().is_integer());
            if (!is_zero(ratnormal(subst(_denom(r,contextptr),x,d.front(),false,contextptr),contextptr)))
                ret+=expand_Dirac_prod(r,x,d.front(),d.back().val,contextptr);
            else rest+=r*(is_zero(d.back())?_Dirac(x-d.front(),contextptr):
                                            _Dirac(makesequence(x-d.front(),d.back()),contextptr));
        }
    }
    if (!is_zero(rest))
        rest=simplify_Dirac(ratnormal(rest,contextptr),x,depth+1,contextptr);
    return ret+rest;
}

/* return true iff g = a*(x+b)^2+c */
bool is_poly2(const gen &g,const identificateur &x,gen &a,gen &b,gen &c,GIAC_CONTEXT) {
    gen st;
    if(!ispoly(g,x,st,contextptr) || !is_one(st-1))
        return false;
    gen cf=_canonical_form(makesequence(g,x),contextptr),e,s(undef),d;
    if (cf.is_symb_of_sommet(at_plus) && cf._SYMBptr->feuille.type==_VECT &&
            cf._SYMBptr->feuille._VECTptr->size()==2) {
        e=cf._SYMBptr->feuille._VECTptr->front();
        c=cf._SYMBptr->feuille._VECTptr->back();
    } else {
        e=cf;
        c=0;
    }
    a=gen(1);
    if (e.is_symb_of_sommet(at_neg)) {
        a=gen(-1);
        e=e._SYMBptr->feuille;
    }
    if (e.is_symb_of_sommet(at_prod) && e._SYMBptr->feuille.type==_VECT) {
        const vecteur &factors=*e._SYMBptr->feuille._VECTptr;
        for (const_iterateur it=factors.begin();it!=factors.end();++it) {
            if (is_constant_wrt(*it,x,contextptr))
                a=a*(*it);
            else if (is_undef(s))
                s=*it;
            else return false;
        }
    } else s=e;
    if (s.is_symb_of_sommet(at_pow) && s._SYMBptr->feuille.type==_VECT &&
            s._SYMBptr->feuille._VECTptr->size()==2 && is_one(s._SYMBptr->feuille._VECTptr->back()-1)) {
        gen l=s._SYMBptr->feuille._VECTptr->front();
        if (!is_linear_wrt(l,x,d,b,contextptr) || !is_one(d))
            return false;
    } else return false;
    return true;
}

/* return true iff g = a*(b*x+c)^n */
bool is_monomial(const gen &g,const identificateur &x,gen &n,gen &a,gen &b,gen &c,GIAC_CONTEXT) {
    if (!ispoly(g,x,n,contextptr))
        return false;
    gen gf=factorise(makesequence(g,x),contextptr),s(undef);
    a=gen(1);
    if (gf.is_symb_of_sommet(at_prod) && gf._SYMBptr->feuille.type==_VECT) {
        const vecteur &factors=*gf._SYMBptr->feuille._VECTptr;
        for (const_iterateur it=factors.begin();it!=factors.end();++it) {
            if (is_constant_wrt(*it,x,contextptr))
                a=a*(*it);
            else if (is_undef(s))
                s=*it;
            else return false;
        }
    } else s=gf;
    return s.is_symb_of_sommet(at_pow) && s._SYMBptr->feuille.type==_VECT &&
           is_linear_wrt(s._SYMBptr->feuille._VECTptr->front(),x,b,c,contextptr);
}

/* return true iff g = a*|x+b| */
bool is_linabs(const gen &g_orig,const identificateur &x,gen &a,gen &b,GIAC_CONTEXT) {
    gen g=factorise(g_orig,contextptr),A(0),B;
    a=gen(1);
    if (g.is_symb_of_sommet(at_neg)) {
        a=gen(-1);
        g=g._SYMBptr->feuille;
    }
    vecteur factors;
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT)
        factors=*g._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,g);
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (is_constant_wrt(*it,x,contextptr))
            a=a*(*it);
        else if (is_zero(A) && it->is_symb_of_sommet(at_abs) &&
                    is_linear_wrt(it->_SYMBptr->feuille,x,A,B,contextptr) && !is_zero(A)) {
            a=a*_abs(A,contextptr);
            b=B/A;
        } else return false;
    }
    return true;
}

iterateur find_lin(vecteur &lv,const gen &key,GIAC_CONTEXT) {
    iterateur it=lv.begin();
    for (;it!=lv.end();++it) {
        if (is_same(key,it->_VECTptr->front(),contextptr))
            break;
    }
    return it;
}

void extract_linpow(const gen &g,const identificateur &x,
                    const gen &e,gen &c,vecteur &lv,gen &rest,GIAC_CONTEXT) {
    gen a,b;
    if (is_constant_wrt(g,x,contextptr))
        c=c*_pow(makesequence(g,e),contextptr);
    else if (g.is_symb_of_sommet(at_neg)) {
        c=c*_pow(makesequence(-1,e),contextptr);
        extract_linpow(g._SYMBptr->feuille,x,e,c,lv,rest,contextptr);
    } else if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT) {
        const_iterateur it=g._SYMBptr->feuille._VECTptr->begin();
        for (;it!=g._SYMBptr->feuille._VECTptr->end();++it) {
            extract_linpow(*it,x,e,c,lv,rest,contextptr);
        }
    } else if (g.is_symb_of_sommet(at_pow)) {
        const vecteur &v=*g._SYMBptr->feuille._VECTptr;
        extract_linpow(v.front(),x,e*v.back(),c,lv,rest,contextptr);
    } else if (g.is_symb_of_sommet(at_sqrt)) {
        extract_linpow(g._SYMBptr->feuille,x,e/2,c,lv,rest,contextptr);
    } else if (g.is_symb_of_sommet(at_abs) && is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr)) {
        c=c*_pow(makesequence(_abs(a,contextptr),e),contextptr);
        iterateur it=find_lin(lv,b/a,contextptr);
        if (it==lv.end())
            lv.push_back(makevecteur(b/a,0,e));
        else (*(it->_VECTptr))[2]+=e;
    } else if (is_linear_wrt(g,x,a,b,contextptr)) {
        c=c*_pow(makesequence(a,e),contextptr);
        iterateur it=find_lin(lv,b/a,contextptr);
        if (it==lv.end())
            lv.push_back(makevecteur(b/a,e,0));
        else (*(it->_VECTptr))[1]+=e;
    } else rest=rest*_pow(makesequence(g,e),contextptr);
}

bool is_linpow(const gen &g,const identificateur &x,gen &a,gen &b,gen &c,bool &s,GIAC_CONTEXT) {
    gen c1(1),c2(1),r1(1),r2(1),h=simplify_signs(g,x,contextptr);
    vecteur lv1,lv2;
    gen num=_numer(h,contextptr),den=_denom(h,contextptr);
    if (!is_constant_wrt(num,x,contextptr)) num=factorise(num,contextptr);
    if (!is_constant_wrt(den,x,contextptr)) den=factorise(den,contextptr);
    extract_linpow(num,x,1,c1,lv1,r1,contextptr);
    extract_linpow(den,x,1,c2,lv2,r2,contextptr);
    if (lv1.size()>1 || lv2.size()>1)
        return false;
    if (!lv2.empty()) {
        const vecteur &v=*lv2.front()._VECTptr;
        iterateur it=find_lin(lv1,v.front(),contextptr);
        if (it==lv1.end())
            lv1.push_back(makevecteur(v[0],-v[1],-v[2]));
        else {
            vecteur &w=*(it->_VECTptr);
            w[1]-=v[1];
            w[2]-=v[2];
        }
    }
    if (lv1.size()!=1)
        return false;
    b=ratnormal(lv1.front()._VECTptr->front(),contextptr);
    c=c1/c2;
    gen rest=r1/r2,a1=lv1.front()._VECTptr->at(1),a2=lv1.front()._VECTptr->at(2),A,B;
    if (!is_zero(im(a1,contextptr)) || !is_zero(im(a2,contextptr)) || !a1.is_integer())
        return false;
    a=a1+a2;
    if (is_strictly_greater(a,-1,contextptr)) {
        if (is_zero(_even(a1,contextptr)))
            rest=rest*_sign(x+b,contextptr);
    } else {
        if (!a2.is_integer())
            return false;
        if (is_zero(_even(a2,contextptr)))
            rest=rest*_sign(x+b,contextptr);
    }
    rest=ratnormal(rest,contextptr);
    if (is_constant_wrt(rest,x,contextptr)) {
        c=c*rest;
        s=false;
    } else if (rest.is_symb_of_sommet(at_sign) &&
               is_linear_wrt(rest._SYMBptr->feuille,x,A,B,contextptr) && is_one(A) &&
               is_same(B/A,b,contextptr)) {
        s=true;          
    } else return false;
    return true;
}

bool is_expexp(const gen &g_orig,const identificateur &x,gen &a,gen &b,GIAC_CONTEXT) {
    gen g=expand(g_orig,contextptr);
    if (!g.is_symb_of_sommet(at_plus) || g._SYMBptr->feuille.type!=_VECT)
        return false;
    const vecteur &v=*g._SYMBptr->feuille._VECTptr;
    gen t1(0),t2(undef),a1,b1,a2,b2;
    for (const_iterateur it=v.begin();it!=v.end();++it) {
        if (is_undef(t2) && it->is_symb_of_sommet(at_neg) &&
                 it->_SYMBptr->feuille.is_symb_of_sommet(at_exp) &&
                 is_linear_wrt(it->_SYMBptr->feuille._SYMBptr->feuille,x,a2,b2,contextptr))
            t2=*it;
        else t1+=*it;
    }
    if (!is_linear_wrt(t1,x,a1,b1,contextptr) || !is_zero(b1) || !is_same(a1,a2,contextptr))
        return false;
    a=a1;
    b=b2;
    return true;
}

int is_expsum(const gen &g_orig,const identificateur &x,gen &a,gen &b,gen &c,bool &sq,GIAC_CONTEXT) {
    gen g=factorise(_lin(expand(g_orig,contextptr),contextptr),contextptr);
    c=gen(1);
    if (g.is_symb_of_sommet(at_neg)) {
        c=gen(-1);
        g=g._SYMBptr->feuille;
    }
    vecteur factors;
    vecteur A(2),B(2),S(3);
    gen rpart,ipart,p;
    int i=0;
    sq=false;
    if (g.is_symb_of_sommet(at_prod))
        factors=*g._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,g);
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (is_constant_wrt(*it,x,contextptr))
            c=c*(*it);
        else if (it->is_symb_of_sommet(at_plus) && it->_SYMBptr->feuille.type==_VECT &&
                 it->_SYMBptr->feuille._VECTptr->size()>=2) {
            rpart=re(*it,contextptr);
            ipart=im(*it,contextptr);
            if (is_zero(rpart)) {
                p=ipart;
                c=c*cst_i;
            } else if (is_zero(ipart))
                p=rpart;
            else return 0;
            if (!p.is_symb_of_sommet(at_plus) || p._SYMBptr->feuille.type!=_VECT ||
                    p._SYMBptr->feuille._VECTptr->size()<2)
                return 0;
            const vecteur &terms=*(p._SYMBptr->feuille._VECTptr);
            for (const_iterateur jt=terms.begin();jt!=terms.end();++jt) {
                if (i<2 && jt->is_symb_of_sommet(at_exp) &&
                        is_linear_wrt(jt->_SYMBptr->feuille,x,A[i],B[i],contextptr)) {
                    S[i]=gen(1);
                    ++i;
                } else if (i<2 && jt->is_symb_of_sommet(at_neg) &&
                        jt->_SYMBptr->feuille.is_symb_of_sommet(at_exp) &&
                        is_linear_wrt(jt->_SYMBptr->feuille._SYMBptr->feuille,x,A[i],B[i],contextptr)) {
                    S[i]=gen(-1);
                    ++i;
                } else if (!sq && is_one(_abs(S[2]=*jt/2,contextptr)) && terms.size()==3)
                    sq=true;
                else return 0;
            }
        } else return 0;
    }
    if (i<2 || !is_zero(ratnormal(A[0]+A[1],contextptr)) || !is_zero(ratnormal(B[0]+B[1],contextptr)))
        return 0;
    if (sq) {
        if (!is_one(S[0]) || !is_one(S[1]))
            return false;
        a=A[0]/2;
        b=B[0]/2;
        return 1+(1+S[2]).val/2;
    }
    c=c*S[0];
    a=A[0];
    b=B[0];
    return 1+int(is_one(S[0]*S[1]));
}

bool is_integral(const gen &g,identificateur &x,gen &f,gen &a,gen &b,GIAC_CONTEXT) {
    if ((!g.is_symb_of_sommet(at_Int) && !g.is_symb_of_sommet(at_int) &&
            !g.is_symb_of_sommet(at_integrate)) || g._SYMBptr->feuille.type!=_VECT)
        return false;
    const vecteur &args=*g._SYMBptr->feuille._VECTptr;
    if (args.size()!=2 && args.size()!=4)
        return false;
    f=args.front();
    if (args.size()==2) {
        if (!is_equal(args.back()) ||
                args.back()._SYMBptr->feuille._VECTptr->front().type!=_IDNT ||
                !args.back()._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_interval))
            return false;
        x=*args.back()._SYMBptr->feuille._VECTptr->front()._IDNTptr;
        a=args.back()._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->front();
        b=args.back()._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille._VECTptr->back();
        return true;
    }
    if (args[1].type!=_IDNT)
        return false;
    x=*args[1]._IDNTptr;
    a=args[2];
    b=args[3];
    return true;
}

bool is_convolution(const gen &e,const identificateur &x,gen &f,gen &g,GIAC_CONTEXT) {
    gen a,b,h;
    identificateur t;
    if (!is_integral(e,t,h,a,b,contextptr) || !a.is_symb_of_sommet(at_neg) ||
            !is_inf(a._SYMBptr->feuille) || !b.is_symb_of_sommet(at_plus) ||
            !is_inf(b._SYMBptr->feuille))
        return false;
    vecteur factors;
    if (h.is_symb_of_sommet(at_prod) && h._SYMBptr->feuille.type==_VECT)
        factors=*h._SYMBptr->feuille._VECTptr;
    else factors=vecteur(1,h);
    f=g=gen(1);
    for (const_iterateur it=factors.begin();it!=factors.end();++it) {
        if (is_constant_wrt(*it,x,contextptr))
            f=f*(*it);
        else g=g*(*it);
    }
    g=recursive_normal(subst(g,x,x+t,false,contextptr),contextptr);
    return is_constant_wrt(g,t,contextptr);
}

int fourier_func_type(const gen &g,const identificateur &x,vecteur &params,GIAC_CONTEXT) {
    params.clear();
    if (is_constant_wrt(g,x,contextptr)) {
        params.push_back(g);
        return _FOURIER_FUNCTYPE_ONE;
    }
    gen a,b,c,d,n;
    vecteur tbl,vars,deg;
    identificateur fcn;
    bool sgn,sq;
    int dg,hyp;
    if (is_rational_wrt(g,x)) {
        gen num=_numer(g,contextptr),den=_denom(g,contextptr);
        assert(is_constant_wrt(num,x,contextptr));
        if (is_linear_wrt(den,x,a,b,contextptr)) {
            params=makevecteur(num/a,b/a,1);
        } else if (is_monomial(den,x,n,a,b,c,contextptr)) {
            assert(!is_zero(a) && !is_zero(b));
            params=makevecteur(num/(a*_pow(makesequence(b,n),contextptr)),c/b,n);
        } else return _FOURIER_FUNCTYPE_UNKNOWN;
        return _FOURIER_FUNCTYPE_INV_MONOM;
    }
    if (g.is_symb_of_sommet(at_sign) && is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(a,b);
        return _FOURIER_FUNCTYPE_SGN;
    }
    if (g.is_symb_of_sommet(at_Dirac) && ((g._SYMBptr->feuille.type==_VECT &&
            is_linear_wrt(g._SYMBptr->feuille._VECTptr->front(),x,a,b,contextptr)) ||
            is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr)) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(a,b/a,
            g._SYMBptr->feuille.type==_VECT?g._SYMBptr->feuille._VECTptr->back():gen(0));
        return _FOURIER_FUNCTYPE_DIRAC;
    }
    if (g.is_symb_of_sommet(at_BesselJ) && g._SYMBptr->feuille.type==_VECT &&
            g._SYMBptr->feuille._VECTptr->size()==2 &&
            (c=g._SYMBptr->feuille._VECTptr->front()).is_integer() && is_positive(c,contextptr) &&
            is_linear_wrt(g._SYMBptr->feuille._VECTptr->back(),x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(a,b,c);
        return _FOURIER_FUNCTYPE_BESSELJ;
    }
    if (g.is_symb_of_sommet(at_ln) && g._SYMBptr->feuille.is_symb_of_sommet(at_abs) &&
            is_linear_wrt(g._SYMBptr->feuille._SYMBptr->feuille,x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params.push_back(b/a);
        return _FOURIER_FUNCTYPE_LOGABSX;
    }
    if (g.is_symb_of_sommet(at_inv) && g._SYMBptr->feuille.is_symb_of_sommet(at_abs) &&
            is_linear_wrt(g._SYMBptr->feuille._SYMBptr->feuille,x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(symbolic(at_inv,_abs(a,contextptr)),b/a);
        return _FOURIER_FUNCTYPE_INVABSX;
    }
    if (is_linpow(g,x,a,b,c,sgn,contextptr) && is_zero(im(b,contextptr))) {
        if (!sgn && is_greater(-1,a,contextptr)) {
            params=makevecteur(c,b,-a);
            return _FOURIER_FUNCTYPE_INV_MONOM;
        }
        params=makevecteur(a,b,c,sgn?gen(1):gen(0));
        return _FOURIER_FUNCTYPE_ABSX_ALPHA;
    }
    if (g.is_symb_of_sommet(at_exp) &&
            is_poly2(g._SYMBptr->feuille,x,a,b,c,contextptr) &&
            is_zero(im(b,contextptr)) && is_positive(re(-a,contextptr),contextptr)) {
        params=makevecteur(exp(c,contextptr),b,-a);
        return _FOURIER_FUNCTYPE_GAUSSIAN;
    }
    if (g.is_symb_of_sommet(at_Gamma) && is_one(re(g._SYMBptr->feuille,contextptr)) &&
            is_linear_wrt(im(g._SYMBptr->feuille,contextptr),x,a,b,contextptr)) {
        params=makevecteur(a,b/a);
        return _FOURIER_FUNCTYPE_GAMMA;
    }
    if (g.is_symb_of_sommet(at_prod) && g._SYMBptr->feuille.type==_VECT &&
            g._SYMBptr->feuille._VECTptr->size()==2 &&
            ((g._SYMBptr->feuille._VECTptr->front().is_symb_of_sommet(at_exp) &&
              is_linear_wrt(g._SYMBptr->feuille._VECTptr->front()._SYMBptr->feuille,x,a,d,contextptr) &&
              g._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_Heaviside) &&
              is_linear_wrt(g._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille,x,b,c,contextptr)) ||
             (g._SYMBptr->feuille._VECTptr->back().is_symb_of_sommet(at_exp) &&
              is_linear_wrt(g._SYMBptr->feuille._VECTptr->back()._SYMBptr->feuille,x,a,d,contextptr) &&
              g._SYMBptr->feuille._VECTptr->front().is_symb_of_sommet(at_Heaviside) &&
              is_linear_wrt(g._SYMBptr->feuille._VECTptr->front()._SYMBptr->feuille,x,b,c,contextptr)))
              && is_zero(im(b,contextptr)) && is_zero(im(c,contextptr)) && !is_zero(b) && is_zero(d) &&
              is_strictly_positive(-a/b,contextptr)) {
        params=makevecteur(a,b,c);
        return _FOURIER_FUNCTYPE_EXP_HEAVISIDE;
    }
    if (g.is_symb_of_sommet(at_Heaviside) && is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(a,b/a);
        return _FOURIER_FUNCTYPE_HEAVISIDE;
    }
    if (g.is_symb_of_sommet(at_exp) && is_linabs(g._SYMBptr->feuille,x,a,b,contextptr) &&
            is_strictly_positive(-a,contextptr)) {
        params=makevecteur(-a,b);
        return _FOURIER_FUNCTYPE_EXPABSX;
    }
    if (g.is_symb_of_sommet(at_exp) && is_expexp(g._SYMBptr->feuille,x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(a,b);
        return _FOURIER_FUNCTYPE_EXPEXP;
    }
    if (g.is_symb_of_sommet(at_Airy_Ai) && is_linear_wrt(g._SYMBptr->feuille,x,a,b,contextptr) &&
            is_zero(im(a,contextptr)) && is_zero(im(b,contextptr))) {
        params=makevecteur(a,b/a);
        return _FOURIER_FUNCTYPE_AIRY_AI;
    }
    if (g.is_symb_of_sommet(at_inv) && (hyp=is_expsum(g._SYMBptr->feuille,x,a,b,c,sq,contextptr))!=0) {
        params=makevecteur(c*gen(sq?4:2),a,b/a);
        if (hyp==1 && !sq) return _FOURIER_FUNCTYPE_COSECH;
        if (hyp==2) return sq?_FOURIER_FUNCTYPE_SECH_2:_FOURIER_FUNCTYPE_SECH;
        return _FOURIER_FUNCTYPE_UNKNOWN;
    }
    if (is_convolution(g,x,a,b,contextptr)) {
        params=makevecteur(a,b);
        return _FOURIER_FUNCTYPE_CONVOLUTION;
    }
    if ((g.is_symb_of_sommet(at_piecewise) || g.is_symb_of_sommet(at_when)) &&
            g._SYMBptr->feuille.type==_VECT) {
        const vecteur &args=*g._SYMBptr->feuille._VECTptr;
        bool cnd=true;
        for (const_iterateur it=args.begin();it!=args.end();++it) {
            if (cnd) {
                if (it+1==args.end()) {
                    params.push_back(symbolic(at_plus,unsigned_inf));
                    params.push_back(*it);
                } else if ((it->is_symb_of_sommet(at_inferieur_egal) ||
                            it->is_symb_of_sommet(at_inferieur_strict)) &&
                           it->_SYMBptr->feuille.type==_VECT &&
                           it->_SYMBptr->feuille._VECTptr->front()==x &&
                           is_constant_wrt(it->_SYMBptr->feuille._VECTptr->back(),x,contextptr))
                    params.push_back(it->_SYMBptr->feuille._VECTptr->back());
                else if ((it->is_symb_of_sommet(at_superieur_egal) ||
                          it->is_symb_of_sommet(at_superieur_strict)) &&
                         it->_SYMBptr->feuille.type==_VECT &&
                         it->_SYMBptr->feuille._VECTptr->back()==x &&
                         is_constant_wrt(it->_SYMBptr->feuille._VECTptr->front(),x,contextptr))
                    params.push_back(it->_SYMBptr->feuille._VECTptr->front());
                else return _FOURIER_FUNCTYPE_UNKNOWN;
            } else params.push_back(*it);
            cnd=!cnd;
        }
        assert(params.size()%2==0);
        int i,n=params.size()/2;
        for (i=0;i<n && (i==0 || is_greater(params[2*i],params[2*(i-1)],contextptr));++i);
        return (n>1 && i==n)?_FOURIER_FUNCTYPE_PIECEWISE:_FOURIER_FUNCTYPE_UNKNOWN;
    }
    if (is_diff(g,x,a,dg)) {
        params=makevecteur(a,dg);
        return _FOURIER_FUNCTYPE_DIFF;
    }
    if (is_partialdiff(g,fcn,vars,deg,contextptr) &&
            !(tbl=fourier_table_find(fcn,vars.size())).empty() &&
            is_linear_wrt(vars[tbl[2].val],x,a,b,contextptr)){
        params=mergevecteur(mergevecteur(tbl,makevecteur(a,b)),mergevecteur(vars,deg));
        return _FOURIER_FUNCTYPE_PARTIAL_DIFF;
    }
    if (is_func(g,fcn,vars) && !(tbl=fourier_table_find(fcn,vars.size())).empty() &&
            is_linear_wrt(vars[tbl[2].val],x,a,b,contextptr)) {        
        params=mergevecteur(tbl,mergevecteur(makevecteur(a,b),vars));
        return _FOURIER_FUNCTYPE_FUNC;
    }
    // some functions with denominator equal to ax+b
    gen gnorm=ratnormal(g,contextptr);
    if (is_linear_wrt(_denom(gnorm,contextptr),x,c,d,contextptr) && !is_zero(c)) {
        gen h=_numer(gnorm,contextptr);
        if (h.is_symb_of_sommet(at_atan) && is_linear_wrt(h._SYMBptr->feuille,x,a,b,contextptr) &&
                !is_zero(a) && is_same(b*c,a*d,contextptr)) {
            params=makevecteur(a,c,b/a);
            return _FOURIER_FUNCTYPE_ATAN_OVERX;
        }
        if (h.is_symb_of_sommet(at_exp) && is_linabs(h._SYMBptr->feuille,x,a,b,contextptr) &&
                is_same(b*c,d,contextptr) && is_strictly_positive(-a,contextptr)) {
            params=makevecteur(-a,b,c);
            return _FOURIER_FUNCTYPE_EXPABSX_OVERX;
        }
    }
    return _FOURIER_FUNCTYPE_UNKNOWN;
}

gen derive_FT(const gen &g,const identificateur &x,int n,GIAC_CONTEXT) {
    gen d=g;
    for (int i=0;i<n;++i) {
        d=_derive(makesequence(simplify_signs(sign2Heaviside(d,contextptr),x,contextptr),x),contextptr);
        //d=simplify_Dirac(d,x,contextptr);
    }
    return d;
}

gen fourier(const gen &f_orig,const identificateur &x,const identificateur &s,
            bool neg,bool try_diff,GIAC_CONTEXT) {
    gen fsimp;
    bool do_simp=!has_integral(f_orig) && !has_partial_diff(f_orig);
    if (do_simp) {
        fsimp=lin_Heaviside(sign2Heaviside(simplify_signs(f_orig,x,contextptr),contextptr),x,contextptr);
    } else fsimp=f_orig;
    vecteur terms=fourier_terms(fsimp/(neg?2*cst_pi:gen(1)),x,do_simp,contextptr),p;
    gen ret(0),intgr(0),t,a,b,c,n,minf=symbolic(at_neg,unsigned_inf),pinf=symbolic(at_plus,unsigned_inf);
    gen rest(0);
    vecteur vars,degr;
    int ind,sz;
    for (int ti=0;ti<int(terms.size());++ti) {
        const gen &term=terms[ti];
        assert(term.type==_VECT);
        const vecteur &v=*term._VECTptr;
        gen cnst=v[0],sh=v[1],d=v[2],f=v[3];
        assert(d.is_integer() && is_positive(d,contextptr));
        switch (fourier_func_type(f,x,p,contextptr)) {
        case _FOURIER_FUNCTYPE_PIECEWISE:
            sz=p.size()/2;
            a=p[0];
            b=_Heaviside(a-x,contextptr)*p[1]+_Heaviside(x-p[2*sz-4],contextptr)*p[2*sz-1];
            t=fourier(b,x,s,false,true,contextptr);
            for (ind=1;ind<sz-1;++ind) {
                b=p[2*ind];
                t+=_integrate(makesequence(exp(-cst_i*x*s,contextptr)*p[2*ind+1],x,a,b),contextptr);
                a=b;
            }
            break;
        case _FOURIER_FUNCTYPE_FUNC:
            a=p[3];
            b=p[4];
            c=p[0];
            vars=vecteur(p.begin()+5,p.end());
            vars[p[2].val]=(is_one(c)?s:-s)/a;
            if (vars.size()==1) vars=makevecteur(p[1],vars.front());
            else vars=makevecteur(p[1],gen(vars,_SEQ__VECT));
            t=c*symbolic(at_of,vars)*exp(cst_i*s*b/a,contextptr)/_abs(a,contextptr);
            break;
        case _FOURIER_FUNCTYPE_DIFF:
            t=_pow(makesequence(cst_i*s,p[1]),contextptr)*fourier(p[0],x,s,false,false,contextptr);
            break;
        case _FOURIER_FUNCTYPE_PARTIAL_DIFF:
            a=p[3];
            b=p[4];
            c=p[0];
            vars=vecteur(p.begin()+5,p.end()-(p.size()-5)/2);
            vars[p[2].val]=(is_one(c)?s:-s)/a;
            vars=makevecteur(p[1],gen(vars,_SEQ__VECT));
            degr=vecteur(p.begin()+5+(p.size()-5)/2,p.end());
            t=c*symbolic(at_of,vars)*exp(cst_i*s*b/a,contextptr)/_abs(a,contextptr);
            t=t*_pow(makesequence(cst_i*s/a,degr[p[2].val]),contextptr);
            for (int i=0;i<int(degr.size());++i) {
                if (i==p[2].val) continue;
                t=_derive(makesequence(t,vars[i],degr[i]),contextptr);
            }
            break;
        case _FOURIER_FUNCTYPE_ONE:
            t=2*p[0]*cst_pi*_Dirac(s,contextptr);
            break;
        case _FOURIER_FUNCTYPE_DIRAC:
            t=_pow(makesequence(cst_i*s/p[0],p[2]),contextptr)*exp(cst_i*p[1]*s,contextptr)/
              _abs(p[0],contextptr);
            break;
        case _FOURIER_FUNCTYPE_HEAVISIDE:
            t=exp(cst_i*s*p[1],contextptr)*
              (_sign(p[0],contextptr)/(cst_i*s)+cst_pi*_Dirac(s,contextptr));
            break;
        case _FOURIER_FUNCTYPE_EXP_HEAVISIDE:
            t=exp(p[2]*(cst_i*s-p[0])/p[1],contextptr)*_sign(p[1],contextptr)/(cst_i*s-p[0]);
            break;
        case _FOURIER_FUNCTYPE_BESSELJ:
            t=2*_pow(makesequence(-cst_i,p[2]),contextptr)*exp(cst_i*p[1]*s/p[0],contextptr)*
              _tchebyshev1(makesequence(p[2],s/p[0]),contextptr)*_rect(s/(2*p[0]),contextptr)/
              (_abs(p[0],contextptr)*sqrt(1-pow(s/p[0],2),contextptr));
            break;
        case _FOURIER_FUNCTYPE_LOGABSX:
            t=-exp(cst_i*p[0]*s,contextptr)*cst_pi*
              (_inv(_abs(s,contextptr),contextptr)+2*cst_euler_gamma*_Dirac(s,contextptr));
            break;
        case _FOURIER_FUNCTYPE_INVABSX:
            t=-2*p[0]*exp(cst_i*p[1]*s,contextptr)*(ln(_abs(s,contextptr),contextptr)+cst_euler_gamma);
            break;
        case _FOURIER_FUNCTYPE_GAUSSIAN:
            if (is_strictly_positive(re(p[2],contextptr),contextptr))
                t=p[0]*exp(cst_i*p[1]*s,contextptr)*
                  sqrt(cst_pi/p[2],contextptr)*exp(-pow(s,2)/(4*p[2]),contextptr);
            else {
                b=im(-p[2],contextptr);
                a=_sign(b,contextptr);
                b=_abs(b,contextptr);
                assert(!is_zero(b));
                t=p[0]*sqrt(cst_pi/b,contextptr)*exp(-a*cst_i*(pow(s,2)/(4*b)-cst_pi/4),contextptr);
            }
            break;
        case _FOURIER_FUNCTYPE_GAMMA:
            t=2*cst_pi*exp(cst_i*p[1]*s-s/p[0]-exp(-s/p[0],contextptr),contextptr)/
              _abs(p[0],contextptr);
            break;
        case _FOURIER_FUNCTYPE_EXPABSX:
            t=2*p[0]*exp(cst_i*p[1]*s,contextptr)/(pow(p[0],2)+pow(s,2));
            break;
        case _FOURIER_FUNCTYPE_EXPEXP:
            t=exp(cst_i*p[1]*s/p[0]-p[1],contextptr)*
              _Gamma(cst_i*s/p[0]+1,contextptr)/_abs(p[0],contextptr);
            break;
        case _FOURIER_FUNCTYPE_AIRY_AI:
            t=exp(cst_i*(p[1]*s+pow(s/p[0],3)/3),contextptr)/_abs(p[0],contextptr);
            break;
        case _FOURIER_FUNCTYPE_ATAN_OVERX:
            t=cst_pi*_sign(p[0],contextptr)*exp(cst_i*p[2]*s,contextptr)*
              _upper_incomplete_gamma(makesequence(0,_abs(s/p[0],contextptr)),contextptr)/p[1];
            break;
        case _FOURIER_FUNCTYPE_EXPABSX_OVERX:
            t=-2*cst_i*atan(s/p[0],contextptr)*exp(cst_i*p[1]*s,contextptr)/p[2];
            break;
        case _FOURIER_FUNCTYPE_COSECH:
            t=-cst_i*cst_pi*exp(cst_i*s*p[2],contextptr)*tanh(cst_pi*s/(2*p[1]),contextptr)/
              (_abs(p[1],contextptr)*p[0]);
            break;
        case _FOURIER_FUNCTYPE_SECH:
            t=cst_pi*exp(cst_i*p[2]*s,contextptr)/
              (p[0]*cosh(cst_pi*s/(2*p[1]),contextptr)*_abs(p[1],contextptr));
            break;
        case _FOURIER_FUNCTYPE_SECH_2:
            t=s*cst_pi*exp(cst_i*p[2]*s,contextptr)/
              (p[0]*p[1]*sinh(cst_pi*s/(2*p[1]),contextptr)*_abs(p[1],contextptr));
            break;
        case _FOURIER_FUNCTYPE_CONVOLUTION:
            t=fourier(p[0],x,s,false,true,contextptr)*fourier(p[1],x,s,false,true,contextptr);
            break;
        case _FOURIER_FUNCTYPE_INV_MONOM:
        case _FOURIER_FUNCTYPE_ABSX_ALPHA:
        default:
            rest+=cnst*exp(cst_i*sh*x,contextptr)*pow(x,d.val)*f;
            continue;
        }
        t=cnst*subst(_pow(makesequence(cst_i,d),contextptr)*derive_FT(t,s,d.val,contextptr),
                     s,s-sh,false,contextptr);
        if (neg) t=subst(t,s,-s,false,contextptr);
        ret+=_exp2pow(_lin(recursive_normal(t,contextptr),contextptr),contextptr);
    }
    if (!is_zero(rest)) {
        terms=fourier_terms(Heaviside2sign(rest,contextptr),x,do_simp,contextptr);
        for (const_iterateur it=terms.begin();it!=terms.end();++it) {
            assert(it->type==_VECT);
            const vecteur &v=*(it->_VECTptr);
            gen cnst=v[0],sh=v[1],d=v[2],f=v[3];
            assert(d.is_integer() && is_positive(d,contextptr));
            switch (fourier_func_type(f,x,p,contextptr)) {
            case _FOURIER_FUNCTYPE_SGN:
                t=-2*cst_i*exp(cst_i*p[1]*s/p[0],contextptr)*_sign(p[0],contextptr)/s;
                break;
            case _FOURIER_FUNCTYPE_INV_MONOM:
                n=p[2]-1;
                if (is_zero(im(p[1],contextptr))) {
                    t=_pow(makesequence(-cst_i*s,n),contextptr)*
                      exp(cst_i*p[1]*s,contextptr)*_sign(s,contextptr);
                } else {
                    gen sgn=_sign(im(p[1],contextptr),contextptr);
                    t=2*sgn*_pow(makesequence(-cst_i*s,n),contextptr)*
                      exp(cst_i*p[1]*s,contextptr)*_Heaviside(sgn*s,contextptr);
                }
                t=-t*p[0]*cst_pi*cst_i/_factorial(n,contextptr);
                break;
            case _FOURIER_FUNCTYPE_ABSX_ALPHA:
                t=p[2]*exp(cst_i*p[1]*s,contextptr);
                if (is_greater(-1,p[0],contextptr)) {
                    assert(p[0].is_integer() && is_one(p[3]));
                    n=-p[0]-1;
                    t=-2*t*_pow(makesequence(-cst_i*s,n),contextptr)*
                      (ln(_abs(s,contextptr),contextptr)+cst_euler_gamma)/_factorial(n,contextptr);
                } else {
                    t=-2*t*_Gamma(p[0]+1,contextptr)/
                      _pow(makesequence(_abs(s,contextptr),p[0]+1),contextptr);
                    if (is_one(p[3]))
                        t=t*cst_i*_sign(s,contextptr)*cos(cst_pi*p[0]/2,contextptr);
                    else t=t*sin(cst_pi*p[0]/2,contextptr);
                }
                break;
            default:
                if (!try_diff)
                    return undef;
                f=sign2Heaviside(f,contextptr);
                if (try_diff) {
                    // check whether we can transform the derivative
                    gen df=derive_FT(f,x,1,contextptr);
                    gen F=fourier(df,x,s,false,false,contextptr);
                    if (!is_undef(F) && !is_inf(F)) {
                        t=-cst_i*F/s;
                        break;
                    }
                }
                intgr+=cnst*eval(_Int(makesequence(exp(-cst_i*x*((neg?-s:s)-sh),contextptr)*
                                        pow(x,d.val)*f,x,minf,pinf),contextptr),0,contextptr);
                continue;
            }
            t=cnst*subst(_pow(makesequence(cst_i,d),contextptr)*derive_FT(t,s,d.val,contextptr),
                         s,s-sh,false,contextptr);
            if (neg) t=subst(t,s,-s,false,contextptr);
            ret+=_exp2pow(_lin(recursive_normal(t,contextptr),contextptr),contextptr);
        }
    }
    if (is_zero(recursive_normal(im(ret,contextptr),contextptr)))
        ret=re(ret,contextptr);
    ret=simplify_signs(ret,s,contextptr);
    ret=simplify_Dirac(ret,s,0,contextptr);
    ret=Heaviside2abs(ret,s,contextptr);
    gen rcn=recursive_normal(ret,contextptr);
    if (is_rational_wrt(rcn,s))
        ret=simplify(ret,contextptr);
    else if (is_simpler(rcn,ret))
        ret=rcn;
    if (!is_constant_wrt(ret,s,contextptr)) {
        gen rcn1=ratnormal(sign2Heaviside(ret,contextptr),contextptr);
        gen rcn2=ratnormal(Heaviside2sign(ret,contextptr),contextptr);
        ret=is_simpler(rcn1,rcn2)?rcn1:rcn2;
        rcn=ratnormal(_exp2trig(ret,contextptr),contextptr);
        if (is_simpler(rcn,ret)) ret=rcn;
        rcn=factorise(ret,contextptr);
        if (is_simpler(rcn,ret)) ret=rcn;
    }
    return ret+intgr;
}

gen _fourier(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    gen f_orig;
    identificateur var("x"),tvar(" fourier_s");
    bool has_tvar=false;
    if (g.type==_VECT) {
        if (g.subtype!=_SEQ__VECT)
            return gentypeerr(contextptr);
        const vecteur &args=*g._VECTptr;
        if (args.size()>3 || args.empty())
            return generrdim("Wrong number of input arguments");
        f_orig=args.front();
        if (args.size()>=2) {
            if (args[1].type!=_IDNT)
                return generrtype("Expected an identifier");
            var=*args[1]._IDNTptr;
        }
        if (args.size()==3) {
            if (args[2].type!=_IDNT)
                return generrtype("Expected an identifier");
            if ((tvar=*args[2]._IDNTptr)==var)
                return generr("Original and transform variables must be different");
            has_tvar=true;
        }
    } else f_orig=g;
    gen ret=fourier(f_orig,var,tvar,false,true,contextptr);
    if (!has_tvar)
        ret=subst(ret,tvar,var,false,contextptr);
    return ret;
}
static const char _fourier_s []="fourier";
static define_unary_function_eval (__fourier,&_fourier,_fourier_s);
define_unary_function_ptr5(at_fourier,alias_at_fourier,&__fourier,0,true)

gen _ifourier(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    gen f_orig;
    identificateur var("x"),tvar(" fourier_s");
    bool has_tvar=false;
    if (g.type==_VECT) {
        if (g.subtype!=_SEQ__VECT)
            return gentypeerr(contextptr);
        const vecteur &args=*g._VECTptr;
        if (args.size()>3 || args.empty())
            return generrdim("Wrong number of input arguments");
        f_orig=args.front();
        if (args.size()>=2) {
            if (args[1].type!=_IDNT)
                return generrtype("Expected an identifier");
            var=*args[1]._IDNTptr;
        }
        if (args.size()==3) {
            if (args[2].type!=_IDNT)
                return generrtype("Expected an identifier");
            if ((tvar=*args[2]._IDNTptr)==var)
                return generr("Original and transform variables must be different");
            has_tvar=true;
        }
    } else f_orig=g;
    gen ret=fourier(f_orig,var,tvar,true,true,contextptr);
    if (!has_tvar)
        ret=subst(ret,tvar,var,false,contextptr);
    return ret;
}
static const char _ifourier_s []="ifourier";
static define_unary_function_eval (__ifourier,&_ifourier,_ifourier_s);
define_unary_function_ptr5(at_ifourier,alias_at_ifourier,&__ifourier,0,true)

vecteur apply_window_function(const gen &expr,const identificateur &k,const vecteur &data,int start,int len,GIAC_CONTEXT) {
    vecteur output(len);
    for (int j=0;j<len;++j) {
        output[j]=_evalf(subst(expr,k,gen((double)j),false,contextptr),contextptr)*data[start+j];
    }
    return output;
}

bool nivelate(vecteur &data,int k,const gen &b,const gen &val,const unary_function_ptr *comp,GIAC_CONTEXT) {
    gen r;
    if (has_i(data[k]) && !is_zero((r=_abs(data[k],contextptr)))) {
        if (_eval(symbolic(comp,makesequence(r,b)),contextptr).val!=0) {
            data[k]=val*data[k]/r;
            return true;
        }
    } else if (_eval(symbolic(comp,makesequence(data[k],b)),contextptr).val!=0) {
        data[k]=val;
        return true;
    }
    return false;
}

gen _threshold(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &args=*g._VECTptr;
    if (int(args.size())<2)
        return generrdim("Wrong number of input arguments");
    if (args.front().type!=_VECT)
        return generrtype("Expected a list");
    const vecteur &data=*args.front()._VECTptr;
    gen bnd=args.at(1);
    int n=data.size();
    vecteur output=data;
    if (bnd.type==_VECT) {
        if (int(bnd._VECTptr->size())!=2)
            return generr("Expected a list with two elements");
        gen lb=bnd._VECTptr->front(),ub=bnd._VECTptr->back(),lval,uval;
        if (is_equal(lb)) {
            lval=_rhs(lb,contextptr);
            lb=_lhs(lb,contextptr);
        } else lval=lb;
        if (is_equal(ub)) {
            uval=_rhs(ub,contextptr);
            ub=_lhs(ub,contextptr);
        } else uval=ub;
        for (int k=0;k<n;++k) {
            if (!nivelate(output,k,lb,lval,at_inferieur_strict,contextptr))
                nivelate(output,k,ub,uval,at_superieur_strict,contextptr);
        }
    } else {
        gen val;
        if (is_equal(bnd)) {
            val=_rhs(bnd,contextptr);
            bnd=_lhs(bnd,contextptr);
        } else val=bnd;
        if (!is_real_number(bnd,contextptr))
            return generrtype("Expected a real constant");
        gen comp=at_inferieur_strict,isabs;
        bool absolute=false;
        for (const_iterateur it=args.begin()+2;it!=args.end();++it) {
            if (*it==at_superieur_strict || *it==at_superieur_egal || *it==at_inferieur_egal)
                comp=*it;
            isabs=gen(1);
            if (*it==at_abs || (is_equal(*it) &&
                                it->_SYMBptr->feuille._VECTptr->front()==at_abs &&
                                (isabs=it->_SYMBptr->feuille._VECTptr->back()).type==_INT_ &&
                                isabs.subtype==_INT_BOOLEAN)) {
                if (has_i(data) || !is_strictly_positive(bnd,contextptr))
                    return gensizeerr(contextptr);
                absolute=(bool)isabs.val;
            }
        }
        for (int k=0;k<n;++k) {
            if (absolute) {
                if (_eval(symbolic(comp._FUNCptr,makesequence(_abs(data[k],contextptr),bnd)),contextptr).val!=0)
                    output[k]=is_positive(data[k],contextptr)?val:-val;
            } else nivelate(output,k,bnd,val,comp._FUNCptr,contextptr);
        }
    }
    return _eval(output,contextptr);
}
static const char _threshold_s []="threshold";
static define_unary_function_eval (__threshold,&_threshold,_threshold_s);
define_unary_function_ptr5(at_threshold,alias_at_threshold,&__threshold,0,true)

bool parse_window_parameters(const gen &g,vecteur &data,int &start,int &len,double *alpha,GIAC_CONTEXT) {
    start=0;
    if (g.type==_VECT && g.subtype!=_SEQ__VECT) {
        data=*g._VECTptr;
        len=data.size();
        return true;
    }
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT || g._VECTptr->size()>3 || g._VECTptr->front().type!=_VECT)
        return false;
    vecteur &args=*g._VECTptr;
    data=*args.front()._VECTptr;
    len=data.size();
    bool has_alpha;
    if (is_real_number(args.at(1),contextptr)) {
        has_alpha=true;
        if (!alpha)
            return false;
        *alpha=to_real_number(args.at(1),contextptr).to_double(contextptr);
    } else if (args.size()>2) return false;
    if (args.back().is_symb_of_sommet(at_interval)) {
        gen lh=_lhs(args.back(),contextptr),rh=_rhs(args.back(),contextptr);
        if (!lh.is_integer() || !rh.is_integer() || lh.val<0 || rh.val>=len)
            return false;
        start=lh.val;
        len=rh.val-start+1;
    } else if (args.size()>2 || !has_alpha) return false;
    return true;
}

gen _bartlett_hann_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    double a=0.62,b=0.48,c=0.38;
    gen expr=a-b*_abs(k/(N-1)-fraction(1,2),contextptr)-c*cos(2*k*_IDNT_pi()/(N-1),contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _bartlett_hann_window_s []="bartlett_hann_window";
static define_unary_function_eval (__bartlett_hann_window,&_bartlett_hann_window,_bartlett_hann_window_s);
define_unary_function_ptr5(at_bartlett_hann_window,alias_at_bartlett_hann_window,&__bartlett_hann_window,0,true)

gen _blackman_harris_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen a(0.35875),b(0.48829),c(0.14128),d(0.01168);
    gen K=k*_IDNT_pi()/(N-1),expr=a-b*cos(2*K,contextptr)+c*cos(4*K,contextptr)-d*cos(6*K,contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _blackman_harris_window_s []="blackman_harris_window";
static define_unary_function_eval (__blackman_harris_window,&_blackman_harris_window,_blackman_harris_window_s);
define_unary_function_ptr5(at_blackman_harris_window,alias_at_blackman_harris_window,&__blackman_harris_window,0,true)

gen _blackman_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double alpha=0.16;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<=0)
        return generrtype(_SP_BAD_WINDOW);
    gen K=k*_IDNT_pi()/(N-1),expr=(1-alpha)/2-cos(2*K,contextptr)/2+alpha*cos(4*K,contextptr)/2;
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _blackman_window_s []="blackman_window";
static define_unary_function_eval (__blackman_window,&_blackman_window,_blackman_window_s);
define_unary_function_ptr5(at_blackman_window,alias_at_blackman_window,&__blackman_window,0,true)

gen _bohman_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen K=_abs(2*k/(N-1)-1,contextptr),expr=(1-K)*cos(_IDNT_pi()*K,contextptr)+sin(_IDNT_pi()*K,contextptr)/_IDNT_pi();
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _bohman_window_s []="bohman_window";
static define_unary_function_eval (__bohman_window,&_bohman_window,_bohman_window_s);
define_unary_function_ptr5(at_bohman_window,alias_at_bohman_window,&__bohman_window,0,true)

gen _cosine_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double alpha=1.0;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<=0)
        return generrtype(_SP_BAD_WINDOW);
    gen expr=exp(alpha*ln(sin(k*_IDNT_pi()/(N-1),contextptr),contextptr),contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _cosine_window_s []="cosine_window";
static define_unary_function_eval (__cosine_window,&_cosine_window,_cosine_window_s);
define_unary_function_ptr5(at_cosine_window,alias_at_cosine_window,&__cosine_window,0,true)

gen _gaussian_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double alpha=0.1;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<=0 || alpha>0.5)
        return generrtype(_SP_BAD_WINDOW);
    gen c=(N-1)/2.0,expr=exp(-pow((k-c)/(alpha*c),2)/2,contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _gaussian_window_s []="gaussian_window";
static define_unary_function_eval (__gaussian_window,&_gaussian_window,_gaussian_window_s);
define_unary_function_ptr5(at_gaussian_window,alias_at_gaussian_window,&__gaussian_window,0,true)

gen _hamming_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen a(0.54),b(0.46),expr=a-b*cos(2*_IDNT_pi()*k/(N-1),contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _hamming_window_s []="hamming_window";
static define_unary_function_eval (__hamming_window,&_hamming_window,_hamming_window_s);
define_unary_function_ptr5(at_hamming_window,alias_at_hamming_window,&__hamming_window,0,true)

gen _hann_poisson_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double alpha=1;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&alpha,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen K=2*_IDNT_pi()*k/(N-1);
    gen expr=(1-cos(K,contextptr))*exp(-alpha*_abs(N-1-2*k,contextptr)/(N-1),contextptr)/2;
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _hann_poisson_window_s []="hann_poisson_window";
static define_unary_function_eval (__hann_poisson_window,&_hann_poisson_window,_hann_poisson_window_s);
define_unary_function_ptr5(at_hann_poisson_window,alias_at_hann_poisson_window,&__hann_poisson_window,0,true)

gen _hann_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen expr=pow(sin(_IDNT_pi()*k/(N-1),contextptr),2);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _hann_window_s []="hann_window";
static define_unary_function_eval (__hann_window,&_hann_window,_hann_window_s);
define_unary_function_ptr5(at_hann_window,alias_at_hann_window,&__hann_window,0,true)

gen _parzen_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen K=1-2*k/(N-1),cond=symb_inferieur_egal(symbolic(at_abs,(N-1)/2.0-k),(N-1)/4.0);
    gen f1=1-6*pow(K,2)*(1-_abs(K,contextptr)),f2=2*pow(1-_abs(K,contextptr),3);
    gen expr=symbolic(at_when,makevecteur(cond,f1,f2));
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _parzen_window_s []="parzen_window";
static define_unary_function_eval (__parzen_window,&_parzen_window,_parzen_window_s);
define_unary_function_ptr5(at_parzen_window,alias_at_parzen_window,&__parzen_window,0,true)

gen _poisson_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double alpha=1.0;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&alpha,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen expr=exp(-alpha*_abs(2*k/(N-1)-1,contextptr),contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _poisson_window_s []="poisson_window";
static define_unary_function_eval (__poisson_window,&_poisson_window,_poisson_window_s);
define_unary_function_ptr5(at_poisson_window,alias_at_poisson_window,&__poisson_window,0,true)

gen _riemann_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    gen K=(2*k/(N-1)-1)*_IDNT_pi(),cond=symbolic(at_same,makesequence(k,(N-1)/2.0));
    gen expr=symbolic(at_when,makesequence(cond,1,sin(K,contextptr)/K));
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _riemann_window_s []="riemann_window";
static define_unary_function_eval (__riemann_window,&_riemann_window,_riemann_window_s);
define_unary_function_ptr5(at_riemann_window,alias_at_riemann_window,&__riemann_window,0,true)

gen _triangle_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double L=0;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&L,contextptr) || (L!=1 && L!=-1 && L!=0))
        return generrtype(_SP_BAD_WINDOW);
    gen expr=1-_abs((2*k-N+1)/(N+L),contextptr);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _triangle_window_s []="triangle_window";
static define_unary_function_eval (__triangle_window,&_triangle_window,_triangle_window_s);
define_unary_function_ptr5(at_triangle_window,alias_at_triangle_window,&__triangle_window,0,true)

gen _tukey_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    double alpha=0.5;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,&alpha,contextptr) || alpha<0 || alpha>1)
        return generrtype(_SP_BAD_WINDOW);
    double p=alpha*(N-1)/2.0,q=1-alpha/2;
    gen cond1=symb_inferieur_strict(k,p),cond2=symb_inferieur_egal(k,q*(N-1));
    gen f1=(1+cos(_IDNT_pi()*(k/p-1),contextptr))/2,f2=(1+cos(_IDNT_pi()*(k/p+1-2/alpha),contextptr))/2;
    gen expr=symbolic(at_piecewise,makesequence(cond1,f1,cond2,1,f2));
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _tukey_window_s []="tukey_window";
static define_unary_function_eval (__tukey_window,&_tukey_window,_tukey_window_s);
define_unary_function_ptr5(at_tukey_window,alias_at_tukey_window,&__tukey_window,0,true)

gen _welch_window(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
    vecteur data;
    int start,N;
    identificateur k(" k");
    if (!parse_window_parameters(g,data,start,N,0,contextptr))
        return generrtype(_SP_BAD_WINDOW);
    double p=(N-1)/2.0;
    gen expr=1-pow(1-k/p,2);
    return apply_window_function(expr,k,data,start,N,contextptr);
}
static const char _welch_window_s []="welch_window";
static define_unary_function_eval (__welch_window,&_welch_window,_welch_window_s);
define_unary_function_ptr5(at_welch_window,alias_at_welch_window,&__welch_window,0,true)

/*
 * ANN CLASS IMPLEMENTATION (Artificial Neural Networks)
 * A simple implementation of feedforward neural network (multilayer perceptron) using GSL
 */
#ifdef HAVE_LIBGSL
/* Constructors. */
ann::ann(const vector<int> &topology,int bsize,GIAC_CONTEXT) {
    if (topology.size()<3)
        throw std::runtime_error(gettext("At least three layers are required"));
    ctx=contextptr;
    // block_size: the number of input vectors that can be passed through the network at once
    block_size=bsize;
    // no momentum by default
    momentum=0;
    // default Adam parameters from the paper
    beta1=0.9; beta2=0.999;
    // default initial learning rate and schedule multiplier
    learning_rate=1e-3; schedule=1.0;
    // iteration counter
    iter=accum_count=0;
    // topology
    this->topology=topology;
    // epsilon: used for preventing division by zero (default value from Adam paper)
    eps=1e-8; // epsilon(contextptr);
    // this network is a regression model by default
    _is_classifier=false;
    // construct the network
    create();
    // no regularization by default
    gsl_vector_set_zero(reg_coeff);
    // create temporary symbolic variables
    x=identificateur(" nn_x");
    _purge(x,ctx);
    y.reserve(topology.back());
    t.reserve(topology.back());
    for (int i=0;i<topology.back();++i) {
        y.push_back(identificateur(" y"+print_INT_(i+1)));
        t.push_back(identificateur(" t"+print_INT_(i+1)));
        _purge(makesequence(y.back(),t.back()),ctx);
    }
    fan=makevecteur(identificateur(" fan_in"),identificateur(" fan_out"));
}
ann::ann(const ann &other) {
    assign(other,other.block_size);
}
ann::ann(const ann &other,int bs) {
    assign(other,bs);
}
/* Assign OTHER to this network. */
ann &ann::operator=(const ann &other) {
    deallocate();
    layers.clear();
    cache.clear();
    deltas.clear();
    weights.clear();
    Deltas.clear();
    moment1.clear();
    moment2.clear();
    assign(other,other.block_size);
    return *this;
}
void ann::assign(const ann &other,int bs) {
    ctx=other.ctx;
    topology=other.topology;
    learning_rate=other.learning_rate;
    block_size=bs;
    eps=other.eps;
    momentum=other.momentum;
    beta1=other.beta1;
    beta2=other.beta2;
    iter=other.iter;
    accum_count=other.accum_count;
    schedule=other.schedule;
    y=other.y;
    t=other.t;
    fan=other.fan;
    errfunc=other.errfunc;
    activation=other.activation;
    output_activation=other.output_activation;
    labels=other.labels;
    _is_classifier=other._is_classifier;
    title=other.title;
    create();
    gsl_vector_memcpy(this->reg_coeff,other.reg_coeff);
    int lc=layer_count(),i;
    for (i=0;i<lc;++i) {
        gsl_matrix_memcpy(layers[i],other.layer(i));
        if (i>0) {
            gsl_matrix_memcpy(weights[i-1],other.weight_matrix(i-1));
            gsl_matrix_memcpy(Deltas[i-1],other.weight_updates(i-1));
            gsl_matrix_memcpy(moment1[i-1],other.moment1[i-1]);
            gsl_matrix_memcpy(moment2[i-1],other.moment2[i-1]);
        }
    }
}
/* Create the network parameters and allocate the required space. */
void ann::create() {
    reg_coeff=gsl_vector_alloc(layer_count()-1);
    eout=eout=gsl_matrix_alloc(output_size(),block_size);
    if (reg_coeff==NULL || eout==NULL)
        throw std::bad_alloc();
    int lc=layer_count(),ti,tip,j,sz;
    for (int i=0;i<lc;++i) {
        ti=topology[i];
        sz=ti+(i==lc-1?0:1);
        if (ti<1)
            throw std::runtime_error(gettext("Invalid number of neurons in a layer"));
        layers.push_back(gsl_matrix_alloc(sz,block_size));
        cache.push_back(i==0?NULL:gsl_matrix_alloc(sz,block_size));
        deltas.push_back(i==0?NULL:gsl_matrix_alloc(sz,block_size));
        if (layers.back()==NULL || (i>0 && (deltas.back()==NULL || cache.back()==NULL)))
            throw std::bad_alloc();
        if (i!=lc-1) { // initialize bias neurons and constant cache
            for (int b=0;b<block_size;++b) {
                gsl_matrix_set(layers.back(),ti,b,1.0);
                if (i>0) gsl_matrix_set(cache.back(),ti,b,1.0);
            }
        }
        if (i>0) { // initialize all weights (including bias) to zero
            tip=topology[i-1];
            weights.push_back(gsl_matrix_calloc(tip+1,sz));
            Deltas.push_back(gsl_matrix_calloc(tip+1,sz));
            moment1.push_back(gsl_matrix_calloc(tip+1,sz));
            moment2.push_back(gsl_matrix_calloc(tip+1,sz));
            if (weights.back()==NULL || Deltas.back()==NULL || moment1.back()==NULL || moment2.back()==NULL)
                throw std::bad_alloc();
            gsl_matrix_set(weights.back(),tip,ti,1.0);
        }
    }
}
/* Free all GSL array data in the network. */
void ann::deallocate() {
    if (reg_coeff!=NULL) gsl_vector_free(reg_coeff);
    if (eout!=NULL) gsl_matrix_free(eout);
    vector<gsl_matrix*>::const_iterator jt,jtend;
    for (jt=layers.begin(),jtend=layers.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);
    for (jt=cache.begin(),jtend=cache.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);
    for (jt=deltas.begin(),jtend=deltas.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);
    for (jt=weights.begin(),jtend=weights.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);
    for (jt=Deltas.begin(),jtend=Deltas.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);
    for (jt=moment1.begin(),jtend=moment1.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);
    for (jt=moment2.begin(),jtend=moment2.end();jt!=jtend;++jt)
        if (*jt!=NULL) gsl_matrix_free(*jt);    
}
/* Set momentum mom (for 0<=mom<1 use classical momentum and for mom=1 use Adam). */
void ann::set_momentum(double mom) {
    if (mom<0)
        throw std::runtime_error(gettext("Momentum must be nonnegative"));
    momentum=std::min(1.0,mom);
}
/* Set activation function for LAYER. If LAYER=-1, it is set for all hidden neurons.
 * e.g. set_activation(-1,tanh)
 * or set_activation(2,logistic,[0,2,1.5]) <- assuming layer count of three
 * or set_activation(-1,piecewise(x< -pi/2, -1, x<pi/2, sin(x), 1)). */
void ann::set_activation(bool out,const gen &f,const vecteur &params) {
    gen res,dif;
    vecteur s,act;
    int pf=0;
    if (f==at_id)
        pf=_ANN_LINEAR;
    else if (f==at_logistic)
        pf==_ANN_SIGMOID;
    else if (f==at_tanh)
        pf=_ANN_TANH;
    else if (f.is_integer() && f.subtype==_INT_MAPLECONVERSION)
        pf=f.val;
    if (pf>0) {
        act.push_back(change_subtype(pf,_INT_MAPLECONVERSION));
        if (!is_numericv(params,num_mask_withfrac|num_mask_withint))
            throw std::runtime_error(gettext("Activation function parameters must be numeric"));
        act.push_back(_evalf(params,ctx));
    } else { // custom activation function
        switch (f.type) {
        case _FUNC:
        case _SYMB:
            if (f.type==_FUNC || f.is_symb_of_sommet(at_program)) {
                if (f.type==_SYMB && f._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=params.size()+1)
                    throw std::runtime_error(gettext("Wrong number of arguments to function"));
                res=f(gen(mergevecteur(vecteur(1,x),params),_SEQ__VECT),ctx);
                if (_lname(res,ctx)._VECTptr->size()>1)
                    throw std::runtime_error(gettext("Function contains uninitialized parameters"));
                break;
            } else {
                if (!params.empty())
                    throw std::runtime_error(gettext("No parameters can be passed alongside a symbolic expression"));
                s=*_lname(f,ctx)._VECTptr;
                if (s.size()!=1)
                    throw std::runtime_error(gettext("Expected a univariate symbolic expression"));
                res=subst(f,s.front(),x,false,ctx);
            }
            break;
        default:
            throw std::runtime_error(gettext("Invalid activation function type"));
        }
        dif=_unapply(makesequence(_derive(makesequence(res,x),ctx),x),ctx);
        act=makevecteur(_unapply(makesequence(res,x),ctx),dif);
    }
    if (out) output_activation=act;
    else activation=act;
}
/* Set the error function.
 * e.g. set_error_function(distance2) or, if you have a function f(x,y,a,b,...) where a,b,... are
 * parameters, use: set_error_function(f,[A,B,...]) where A,B,... are the parameter values (fixed). */
void ann::set_error_function(const gen &f,const vecteur &params) {
    if (f.is_integer() && f.subtype==_INT_MAPLECONVERSION) {
        if (!is_numericv(params,num_mask_withfrac|num_mask_withint))
            throw std::runtime_error(gettext("Error function parameters must be numeric"));
        errfunc=makevecteur(f,_evalf(params,ctx));
    } else {
        if (f.type!=_FUNC && !f.is_symb_of_sommet(at_program))
            throw std::runtime_error("Invalid error function type");
        if (f.type==_SYMB && f._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=params.size()+1)
            throw std::runtime_error(gettext("Wrong number of arguments to function"));
        vecteur vars=mergevecteur(t,y);
        gen res=f(gen(mergevecteur(makevecteur(t,y),params),_SEQ__VECT),ctx);
        if ((int)_lname(res,ctx)._VECTptr->size()>2*output_size())
            throw std::runtime_error(gettext("Function contains uninitialized parameters"));
        errfunc=makevecteur(_unapply(gen(mergevecteur(vecteur(1,res),vars),_SEQ__VECT),ctx),
                            _unapply(gen(mergevecteur(vecteur(1,_grad(makesequence(res,y),ctx)),vars),_SEQ__VECT),ctx));
    }
}
/* Set labels for output neurons. */
void ann::set_labels(const vecteur &lab,bool classes) {
    labels=lab;
    if (classes && output_size()==1) {
        if (lab.size()!=2)
            throw std::runtime_error(gettext("Expected 2 labels for binary classification"));
    } else {
        if ((int)lab.size()!=output_size())
            throw std::runtime_error(gettext("Number of labels does not match the output size"));
    }
    _is_classifier=classes;
}
/* Set L2-regularization (weight decay) parameters. */
void ann::set_regularization(const vecteur &reg) const {
    if (reg.size()!=reg_coeff->size)
        throw std::runtime_error(gettext("Invalid number of regularization parameters"));
    if (vecteur2gsl_vector(reg,reg_coeff,ctx)!=GSL_SUCCESS)
        throw std::runtime_error(gettext("Failed to set regularization parameters"));
}
/* Return the network type as a string. */
string ann::network_type() const {
    switch (task()) {
        case 1: return gettext("a classifier");
        case 2: return gettext("a multi-label classifier");
        default: break;
    }
    return gettext("a neural network");
}
/* Apply activation function or its derivative to A. */
void ann::activate(int i,double &a,double *b) const {
    bool ones=_is_classifier && i+1==layer_count();
    const vecteur &f=(i+1==layer_count()?output_activation:activation);
    double e,v;
    if (f.front().is_integer() && f.front().subtype==_INT_MAPLECONVERSION) {
        const vecteur &params=*f.back()._VECTptr;
        bool nop=params.empty();
        switch (f.front().val) {
        case _ANN_LINEAR:
            if (b!=NULL) *b=1.0;
            break;
        case _ANN_SIGMOID:
            e=a>0?std::exp(-a):std::exp(a);
            v=(a>0?1.0:e)/(1.0+e);
            a=v;
            if (b!=NULL) *b=ones?1.0:v*(1-v);
            break;
        case _ANN_TANH:
            v=std::tanh(a);
            a=v;
            if (b!=NULL) *b=ones?1.0:1.0-v*v;
            break;
        case _ANN_RELU: // leaky ReLU if parameter is set
            v=nop?0.0:params.front().to_double(ctx);
            if (b!=NULL) *b=ones||a>0?1.0:v;
            if (a<0) a*=v;
            break;
        default: // not reachable
            a=giac::nan();
            if (b!=NULL) *b=giac::nan();
            break;
        }
    } else {
        if (b!=NULL) *b=ones?1.0:f[1](a,ctx).to_double(ctx);
        a=f[0](a,ctx).to_double(ctx);
    }
}
/* Propagate a block of samples forward. */
void ann::propagate_forward(bool comp_deriv) const {
    int i,j,k,ti,len,lc=layer_count(),os=output_size(),s1;
    gsl_vector_view mc;
    double s,*o,m;
    bool softmax=is_default_classifier() && os>1;
    for (i=1;i<lc;++i) {
        s1=layers[i]->size1;
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,weights[i-1],layers[i-1],0.0,layers[i]);
        if (softmax && i+1==lc) {
            for (k=0;k<block_size;++k) {
                mc=gsl_matrix_column(layers[i],k);
                m=gsl_vector_max(&mc.vector);
                s=0;
                for (j=0;j<os;++j) {
                    o=gsl_vector_ptr(&mc.vector,j);
                    s+=(*o=std::exp(*o-m));
                }
                gsl_vector_scale(&mc.vector,1.0/s);
            }
            if (comp_deriv) gsl_matrix_set_all(cache.back(),1.0);
        } else {
            len=(s1-(i+1==lc?0:1))*block_size;
            for (j=0;j<len;++j)
                activate(i,layers[i]->data[j],comp_deriv?&cache[i]->data[j]:NULL);
        }
    }
} 
/* Return E(t,y). t can be changed! */
gen ann::compute_errfunc(gsl_vector *t,gsl_vector *y) const {
    if (errfunc.size()!=2)
        throw std::runtime_error(gettext("Error function not defined"));
    if (errfunc.front().is_integer() && errfunc.front().subtype==_INT_MAPLECONVERSION) {
        const vecteur &params=*errfunc.back()._VECTptr;
        int i,n=y->size,ef=errfunc.front().val;
        double res,ti,yi;
        switch (ef) {
        case _ANN_HALF_MSE:
            gsl_vector_sub(t,y);
            gsl_blas_ddot(t,t,&res);
            return res/(2.0*n);
        case _ANN_CROSS_ENTROPY:
        case _ANN_LOG_LOSS:
            res=0;
            for (i=0;i<n;++i) {
                ti=gsl_vector_get(t,i);
                yi=gsl_vector_get(y,i);
                if (ti>0.0) res-=ti*std::log(std::max(yi,eps));
                if (ef==_ANN_LOG_LOSS)
                    if (ti<1.0) res-=(1.0-ti)*std::log(std::max(1.0-yi,eps));
            }
            return res;
        default:
            break;
        }
        return undef; // not reachable
    }
    vecteur args=mergevecteur(gsl_vector2vecteur(t),gsl_vector2vecteur(y));
    return errfunc.front()(gen(args,_SEQ__VECT),ctx);
}
/* Compute dE(t,y)/dy and store it to res. */
void ann::compute_errfunc_diff(gsl_matrix *t,gsl_matrix *y,gsl_matrix *res) const {
    if (errfunc.size()!=2)
        throw std::runtime_error(gettext("Error function not defined"));
    int n=y->size1,m=y->size2;
    if (errfunc.front().is_integer() && errfunc.front().subtype==_INT_MAPLECONVERSION) {
        const vecteur &params=*errfunc.back()._VECTptr;
        int i,len=m*n,ef=errfunc.front().val;
        double ti,yi;
        switch (ef) {
        case _ANN_HALF_MSE:
            gsl_matrix_memcpy(res,y);
            gsl_matrix_sub(res,t);
            break;
        case _ANN_CROSS_ENTROPY:
        case _ANN_LOG_LOSS:
            for (i=0;i<len;++i) {
                ti=t->data[i];
                yi=y->data[i];
                if (ti>0.0) res->data[i]=-ti/std::max(yi,eps);
                if (ef==_ANN_LOG_LOSS)
                    if (ti<1.0) res->data[i]-=(1.0-ti)/std::max(1.0-yi,eps);
            }
            break;
        default:
            assert(false); // not reachable
        }
    } else {
        matrice args=mtran(mergevecteur(gsl_matrix2matrice(t),gsl_matrix2matrice(y))),val;
        val.reserve(m);
        const_iterateur it=args.begin(),itend=args.end();
        for (;it!=itend;++it)
            val.push_back(errfunc.back()(gen(*it,_SEQ__VECT),ctx));
        if (matrice2gsl_matrix(val,0,0,0,0,true,res,ctx)!=GSL_SUCCESS)
            throw std::runtime_error(gettext("Failed to compute the error function derivative"));
    }
}
/* Return true iff this network is a classifier with usual error function. */
bool ann::is_default_classifier() const {
    if (!_is_classifier)
        return false;
    const gen &e=errfunc.front();
    if (!e.is_integer() || e.subtype!=_INT_MAPLECONVERSION)
        return false;
    return e.val==(output_size()==1?_ANN_LOG_LOSS:_ANN_CROSS_ENTROPY);
}
/* Calculate the errors made by neurons. */
void ann::calc_deltas() const {
    int i,j,len,lc=layer_count();
    if (is_default_classifier()) {
        gsl_matrix_memcpy(deltas.back(),layers.back());
        gsl_matrix_sub(deltas.back(),eout);
    } else compute_errfunc_diff(eout,layers.back(),deltas.back());
    for (i=lc-2;i>=0;--i) {
        gsl_matrix_mul_elements(deltas[i+1],cache[i+1]);
        if (i>0)
            gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,weights[i],deltas[i+1],0.0,deltas[i]);
    }
}
/* Accumulate weight deltas. */
void ann::update_Deltas() {
    gsl_matrix_view dm,lm;
    int lc=layer_count(),i;
    for (i=0;i<lc-1;++i)
        gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,layers[i],deltas[i+1],1.0,Deltas[i]);
}
/* Update weights. */
void ann::update_weights() {
    if (accum_count==0)
        return;
    ++iter;
    double eta=_evalf(schedule.is_symb_of_sommet(at_program)?schedule(iter,ctx):schedule,ctx).to_double(ctx);
    double lambda,e=std::sqrt(1.0-std::pow(beta2,iter)),c=eps*e,lr=eta*learning_rate,alpha=lr*e/(1.0-std::pow(beta1,iter)),l;
    int i,lc=layer_count();
    size_t j,len,n1,n2,bias_start;
    bool adam=(momentum==1);
    for (i=0;i<lc-1;++i) {
        n1=weights[i]->size1;
        n2=weights[i]->size2;
        len=n1*n2;
        bias_start=n2*(n1-1);
        lambda=gsl_vector_get(reg_coeff,i);
        l=1.0-eta*lambda;
        for (j=0;j<len;++j) {
            if (i<lc-2 && (j+1)%n2==0)
                continue;
            double &d=Deltas[i]->data[j],&w=weights[i]->data[j],&fm=moment1[i]->data[j],&sm=moment2[i]->data[j];
            d/=(double)accum_count;
            if (lambda>0 && j<bias_start) {
                d+=w*lambda;
                w*=l;
            }
            if (adam) {
                fm=beta1*fm+(1.0-beta1)*d;
                sm=beta2*sm+(1.0-beta2)*d*d;
                w-=alpha*fm/(std::sqrt(sm)+c);
            } else w-=(fm=momentum*fm+lr*d);
        }
        gsl_matrix_set_zero(Deltas[i]);
    }
    accum_count=0;
}
/* Backward propagation of the expected output. */
void ann::propagate_backward() {
    calc_deltas();
    update_Deltas();
}
/* Initialize iterative learning. Set ITER to 0 and zero-out the moments. */
void ann::restart_iteration() {
    iter=0;
    int lc=layer_count()-1,i=0;
    for (;i<lc-1;++i) {
        gsl_matrix_set_zero(moment1[i]);
        gsl_matrix_set_zero(moment2[i]);
    }
}
/* Convert output to GSL matrix. */
void ann::output2gsl_matrix(const vecteur &output,int i0,int n,gsl_matrix *m) const {
    int i;
    gsl_matrix_set_zero(m);
    const_iterateur it,jt;
    bool ck;
    switch (this->task()) {
    case 0:
        if (ckmatrix(output))
            assert(matrice2gsl_matrix(output,i0,0,n,0,true,m,ctx)==GSL_SUCCESS);
        else {
            gsl_vector_view m1=gsl_matrix_subrow(m,0,0,n);
            assert(vecteur2gsl_vector(output.begin()+i0,output.begin()+i0+n,&m1.vector,ctx)==GSL_SUCCESS);
        }
        break;
    case 1:
        for (i=0,jt=output.begin()+i0;i<n;++i,++jt) {
            it=std::find(labels.begin(),labels.end(),*jt);
            if (it==labels.end())
                throw std::runtime_error(gettext("Output label not found"));
            if (output_size()==1)
                gsl_matrix_set(m,0,i,it-labels.begin());
            else gsl_matrix_set(m,it-labels.begin(),i,1.0);
        }
        break;
    case 2:
        for (i=0,jt=output.begin()+i0;i<n;++i,++jt) {
            const vecteur &samp=*it->_VECTptr;
            for (it=labels.begin();it!=labels.end();++it) {
                if (contains(samp,*it))
                    gsl_matrix_set(m,it-labels.begin(),i,1.0);
            }
        }
        break;
    }
}
/* Check input-output consistency. */
void ann::ckinput(const matrice &input,const vecteur &output,bool ignore_output) const {
    if (mcols(input)!=input_size())
        throw std::runtime_error(gettext("Invalid input size"));
    if (!is_numericm(input,num_mask_withfrac|num_mask_withint))
        throw std::runtime_error(gettext("Input must be numeric"));
    if (ignore_output) return;
    if (input.size()!=output.size())
        throw std::runtime_error(gettext("Number of samples does not match the number of outputs"));
    if (ckmatrix(output) && mcols(output)!=output_size())
        throw std::runtime_error(gettext("Invalid output dimension"));
    switch (this->task()) {
    case 0:
        if (!is_numericv(output,num_mask_withfrac|num_mask_withint))
            throw std::runtime_error(gettext("Output must be numeric"));
        if (!ckmatrix(output) && output_size()>1)
            throw std::runtime_error(gettext("Invalid output dimension"));
        break;
    case 1:
        break;
    case 2:
        if (!ckmatrix(output))
            throw std::runtime_error(gettext("Output must be a vector of labels"));
        break;
    default:
        throw std::runtime_error(gettext("Invalid network task"));
    }
}
/* Compute the gradient update for the given input and expected output.
 * After batch_size samples are passed forward (or all samples for size 0), update the weights.
 * Setting batch_size to a negative number disables updating the weights. */
void ann::train(const matrice &input,const vecteur &expected_output,int batch_size) {
    int nsamp=input.size(),i,j=0,jprev=0,sz=0,ep;
    if (batch_size>0 && batch_size>nsamp)
        throw std::runtime_error(gettext("Batch size too large"));
    ckinput(input,expected_output,false);
    bool update=false;
    gsl_matrix_view ov;
    while (j<nsamp) {
        j=std::min(j+block_size,nsamp);
        sz+=j-jprev;
        if ((batch_size>0 && sz>=batch_size) || (batch_size>=0 && j==nsamp)) {
            if (batch_size>0 && sz>batch_size)
                j-=sz-batch_size;
            update=true;
            sz=0;
        }
        if (j-jprev<block_size)
            gsl_matrix_set_zero(layers.front());
        assert(matrice2gsl_matrix(input,jprev,0,j-jprev,0,true,layers.front(),ctx)==GSL_SUCCESS);
        propagate_forward(true);
        output2gsl_matrix(expected_output,jprev,j-jprev,eout);
        propagate_backward();
        accum_count+=j-jprev;
        if (update) {
            update_weights();
            update=false;
        }
        jprev=j;
    }
}
/* Convert the output to gen. */
gen ann::gsl_vector2output(gsl_vector *output) const {
    vecteur sel;
    size_t i;
    const_iterateur it;
    switch (this->task()) {
    case 0: // regression
        return gsl_vector2vecteur(output);
    case 1: // classifier
        if (output_size()==1)
            return labels[gsl_vector_get(output,0)>=0.5?1:0];
        return labels[gsl_vector_max_index(output)];
    case 2: // labeller
        sel.reserve(labels.size());
        for (i=0,it=labels.begin();i<output->size;++i,++it) {
            if (gsl_vector_get(output,i)>=0.5)
                sel.push_back(*it);
        }
        return sel;
    default:
        break;
    }
    throw std::runtime_error(gettext("Invalid network task"));
}
/* Feed the input through the network.
 * If expected_output is NULL, return the output in res.
 * If expected_output is not NULL, compute error function values for the obtained
 * and expected output and store these values in res. */
void ann::feed(const matrice &input,vecteur &res,const vecteur &expected_output) const {
    int nsamp=input.size(),i,j=0,jprev=0;
    res.clear();
    if (nsamp==0) return;
    bool comp=!expected_output.empty();
    ckinput(input,expected_output,!comp);
    res.reserve(nsamp);
    gsl_vector_view ov,eov;
    while (j<nsamp) {
        j=std::min(j+block_size,nsamp);
        assert(matrice2gsl_matrix(input,jprev,0,j-jprev,0,true,layers.front(),ctx)==GSL_SUCCESS);
        if (j-jprev<block_size)
            gsl_matrix_set_zero(eout);
        propagate_forward(false);
        if (comp) output2gsl_matrix(expected_output,jprev,j-jprev,eout);
        for (i=0;i<j-jprev;++i) {
            ov=gsl_matrix_column(layers.back(),i);
            if (comp) {
                eov=gsl_matrix_column(eout,i);
                res.push_back(compute_errfunc(&eov.vector,&ov.vector));
            } else res.push_back(gsl_vector2output(&ov.vector));
        }
        jprev=j;
    }
}
/* Feed the argument forward, return the result. */
gen ann::operator()(const gen &g,GIAC_CONTEXT) const {
    if (g.type!=_VECT)
        return generrtype("Input must be a vector or a matrix");
    vecteur res;
    if (g.subtype==_SEQ__VECT) {
        if (g._VECTptr->size()!=2)
            return generrdim("Expected 2 arguments");
        const gen &inp=g._VECTptr->front(),&out=g._VECTptr->back();
        if (inp.type!=_VECT)
            return generrtype("Input must be a vector or a matrix");
        if (ckmatrix(inp)) {
            if (out.type!=_VECT)
                return generrtype("Expected output must be a vector");
            feed(*inp._VECTptr,res,*out._VECTptr);
        } else feed(vecteur(1,inp),res,vecteur(1,out));
    } else if (ckmatrix(g))
        feed(*g._VECTptr,res);
    else feed(vecteur(1,g),res);
    if (ckmatrix(res) && output_size()==1)
        res=*mtran(res).front()._VECTptr;
    return res.size()==1?res.front():res;
}
/* Convert weight matrices to a vector of giac matrices. */
void ann::weights2vecteur(vecteur &res) const {
    int i,lc=layer_count();
    res.clear();
    res.reserve(lc-1);
    gsl_matrix_view mv;
    for (i=0;i<lc-1;++i) {
        mv=gsl_matrix_submatrix(weights[i],0,0,weights[i]->size1,weights[i]->size2-(i+2==lc?0:1));
        res.push_back(gsl_matrix2matrice(&mv.matrix));
    }
}
/* Return the given hyperparameter value or the i-th layer (after feed-forward). */
gen ann::operator[](const gen &g) {
    int i,j,lc=layer_count();
    vecteur ret;
    gsl_matrix_view mv;
    if (!g.is_integer())
        return gentypeerr(ctx);
    i=g.val;
    if (g.subtype==_INT_MAPLECONVERSION) {
        switch (i) {
        case _ANN_LEARNING_RATE:
            if (schedule.is_symb_of_sommet(at_program))
                ret=makevecteur(string2gen("Scheduled",false),learning_rate,schedule);
            else return learning_rate;
            break;
        case _ANN_BLOCK_SIZE:
            return block_size;
        case _ANN_TOPOLOGY:
            ret=vector_int_2_vecteur(topology);
            break;
        case _GT_WEIGHTS:
            weights2vecteur(ret);
            break;
        case _ANN_WEIGHT_DECAY:
            ret=gsl_vector2vecteur(reg_coeff);
            break;
        case _ANN_MOMENTUM:
            if (momentum==1)
                ret=makevecteur(string2gen("Adam",false),beta1,beta2);
            else return momentum;
            break;
        default:
            return generr("Unknown hyperparameter");
        }
    } else if (g.subtype==_INT_PLOT) {
        switch (i) {
        case _LABELS:
            ret=labels;
            break;
        case _TITLE:
            return string2gen(title,false);
        default:
            return generr("Unknown hyperparameter");
        }
    } else {
        if (i<0 || i>=lc)
            return generrdim("Layer index out of range");
        mv=gsl_matrix_submatrix(layers[i],0,0,layers[i]->size1,layers[i]->size2-(i+1==lc?0:1));
        ret=gsl_matrix2matrice(&mv.matrix,true);
    }
    return ret;
}
/* Initialize weights using the initializer f (a random variable). Bias weights are initialized to zero. */
void ann::randomize_initial_weights(const gen &f,const vecteur &fan_in_out) const {
    gen r;
    if (fan_in_out.empty())
        r=f;
    else if (fan_in_out.size()==1)
        r=subst(f,fan_in_out.front(),fan[0],false,ctx);
    else if (fan_in_out.size()==2)
        r=subst(f,fan_in_out,fan,false,ctx);
    else throw std::runtime_error(gettext("Invalid weights initializer specification"));
    int i=0;
    for (;i<layer_count()-1;++i) {
        gen fan_in=topology[i],fan_out=topology[i+1];
        gen args=makesequence(fan_in,fan_out,subst(r,fan,makevecteur(fan_in,fan_out),false,ctx));
        set_initial_weights(*(is_real_number(args._VECTptr->back(),ctx)?_matrix(args,ctx):_ranm(args,ctx))._VECTptr,i);
    }
}
/* Initialize weigths between i-th and (i+1)-th layer directly from a matrix. */
void ann::set_initial_weights(const matrice &w,int i) const {
    if (i<0||i+1>=layer_count())
        throw std::runtime_error(gettext("Invalid layer specification"));
    if (mrows(w)<topology[i] || mrows(w)>topology[i]+1 || mcols(w)!=topology[i+1])
        throw std::runtime_error(gettext("Invalid weight matrix dimensions"));
    if (matrice2gsl_matrix(w,0,0,0,0,false,weights[i],ctx)!=GSL_SUCCESS)
        throw std::runtime_error(gettext("Failed to initialize weights"));
}
/* Return the total number of neurons (including bias neurons). */
int ann::neuron_count() const {
    int s=0;
    for (vector<int>::const_iterator it=topology.begin();it!=topology.end();++it)
        s+=*it;
    return s+layer_count()-1;
}
/* Return an one-line description of this network. */
string ann::print (GIAC_CONTEXT) const {
    stringstream ss;
    if (title.size()>0)
        ss << title << ": ";
    ss << network_type() << " " << gettext("with input of size") << " " << topology.front();
    switch (task()) {
    case 1:
        ss << " " << gettext("and") << " ";
        if (output_size()==1)
            ss << 2;
        else ss << output_size();
        ss << " " << gettext("classes");
        break;
    case 2:
        ss << " " << gettext("and") << " " << output_size() << " " << gettext("labels");
        break;
    default:
        ss << " " << gettext("and") << " " << gettext("output of size") << " " << output_size();
        break;
    }
    return ss.str();
}
/* Return the LaTeX code containing the one-line description of this network. */
string ann::texprint(GIAC_CONTEXT) const {
    string ret="\\text{";
    ret+=print(contextptr);
    ret+="}";
    return ret;
}
/* Return the pointer to the neural network contained in G, or NULL if G is not a neural network. */
ann *ann::from_gen(const gen &g) {
    if (g.type!=_USER)
        return NULL;
    return dynamic_cast<ann*>(g._USERptr);
}
/* Return true iff the other network has the same parameters as this network. */
bool ann::operator ==(const gen &g) const {
    ann *other=from_gen(g);
    if (other==NULL)
        return false;
    if (this->topology!=other->topology || this->labels!=other->labels || this->activation!=other->activation ||
            this->output_activation!=other->output_activation || this->_is_classifier!=other->_is_classifier ||
            this->errfunc!=other->errfunc || this->learning_rate!=other->learning_rate || this->schedule!=other->schedule ||
            this->block_size!=other->block_size || this->momentum!=other->momentum ||
            (momentum==1 && (this->beta1!=other->beta1 || this->beta2!=other->beta2)))
        return false;
    for(int i=0;i<topology.size()-1;++i) {
        if (gsl_vector_get(this->reg_coeff,i)!=gsl_vector_get(other->reg_coeff,i))
            return false;
    }
    return true;
}
/* Implement the giac constructor. */
gen ann::giac_constructor(GIAC_CONTEXT) const {
    vecteur args;
    // topology
    args.push_back(vector_int_2_vecteur(topology));
    // hidden activation
    if (activation.front().is_integer() && activation.front().subtype==_INT_MAPLECONVERSION) {
        vecteur act;
        switch (activation.front().val) {
            case _ANN_LINEAR: act.push_back(at_id); break;
            case _ANN_SIGMOID: act.push_back(at_logistic); break;
            case _ANN_TANH: act.push_back(at_tanh); break;
            default: act.push_back(activation.front()); break;
        }
        act=mergevecteur(act,*activation.back()._VECTptr);
        if (!act.empty())
            args.push_back(symb_equal(change_subtype(_FUNC,_INT_TYPE),act.size()>1?act:act.front()));
    } else args.push_back(symb_equal(change_subtype(_FUNC,_INT_TYPE),activation.front()));
    // output activation
    if (!_is_classifier && output_activation.front().is_integer() && output_activation.front().subtype==_INT_MAPLECONVERSION) {
        vecteur act;
        switch (output_activation.front().val) {
            case _ANN_LINEAR: act.push_back(at_id); break;
            default: break;
        }
        act=mergevecteur(act,*output_activation.back()._VECTptr);
        if (!act.empty())
            args.push_back(symb_equal(at_output,act.size()>1?act:act.front()));
    } else args.push_back(symb_equal(at_output,output_activation.front()));
    // error function
    if (errfunc.front().is_integer() && errfunc.front().subtype==_INT_MAPLECONVERSION) {
        vecteur err(1,errfunc.front());
        err=mergevecteur(err,*errfunc.back()._VECTptr);
        args.push_back(symb_equal(at_erf,err.size()>1?err:err.front()));
    } else args.push_back(symb_equal(at_erf,errfunc.front()));
    // learning rate
    args.push_back(symb_equal(change_subtype(_ANN_LEARNING_RATE,_INT_MAPLECONVERSION),
                              schedule.is_symb_of_sommet(at_program)?makevecteur(learning_rate,schedule):gen(learning_rate)));
    // batch size
    args.push_back(symb_equal(change_subtype(_ANN_BLOCK_SIZE,_INT_MAPLECONVERSION),block_size));
    // weight decay parameters
    args.push_back(symb_equal(change_subtype(_ANN_WEIGHT_DECAY,_INT_MAPLECONVERSION),gsl_vector2vecteur(reg_coeff)));
    // momentum
    args.push_back(symb_equal(change_subtype(_ANN_MOMENTUM,_INT_MAPLECONVERSION),momentum==1?makevecteur(beta1,beta2):momentum));
    // title
    if (title.size()>0)
        args.push_back(symb_equal(change_subtype(_TITLE,_INT_PLOT),string2gen(title,false)));
    // labels
    if (!labels.empty())
        args.push_back(symb_equal(_is_classifier?at_classes:change_subtype(_LABELS,_INT_PLOT),labels));
    // weights
    vecteur w;
    weights2vecteur(w);
    args.push_back(symb_equal(change_subtype(_GT_WEIGHTS,_INT_MAPLECONVERSION),w));
    // finalize
    return symbolic(at_neural_network,change_subtype(args,_SEQ__VECT));
}
#endif // HAVE_LIBGSL
/*
 * END OF ANN CLASS IMPLEMENTATION
 */

gen _neural_network(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
#ifdef HAVE_LIBGSL
    if (g.type!=_VECT)
        return gentypeerr(contextptr);
    vecteur topology,act,out_act,err,lab,winit,reg;
    vector<matrice*> init_weights;
    gen rate(0.001),schedule(1),mom(0.0),fan_in=identificateur(" fan_in"),fan_out=identificateur(" fan_out");
    string name;
    bool cls=false,set_rate=false;
    int i,bs=1;
    ann *src=NULL;
    if (g.subtype==_SEQ__VECT) {
        const vecteur &gv=*g._VECTptr;
        if (gv.empty())
            return generrdim("Sequence of arguments is empty");
        if ((src=ann::from_gen(gv.front()))==NULL) {
            if (gv.front().type!=_VECT)
                return generrtype("Invalid topology specification");
            topology=*gv.front()._VECTptr;
        }
        // parse options
        for (const_iterateur it=gv.begin()+1;it!=gv.end();++it) {
            if (it->is_symb_of_sommet(at_equal)) {
                const gen &prop=it->_SYMBptr->feuille._VECTptr->front();
                const gen &val=it->_SYMBptr->feuille._VECTptr->back();
                if (prop==at_output) {
                    out_act=val.type==_VECT?*val._VECTptr:vecteur(1,val);
                    if (out_act.empty())
                        return generr("Invalid output activation function specification");
                    vecteur params(out_act.begin()+1,out_act.end());
                    out_act=makevecteur(out_act.front(),params);
                } else if (prop==at_classes) {
                    if (val.type!=_VECT)
                        return generrtype("Expected a list of classes");
                    lab=*val._VECTptr;
                    cls=true;
                } else if (prop==at_erf) {
                    err=val.type==_VECT?*val._VECTptr:vecteur(1,val);
                    if (err.empty())
                        return generr("Invalid error function specification");
                    vecteur params(err.begin()+1,err.end());
                    err=makevecteur(err.front(),params);
                } else if (prop.is_integer() && prop.subtype==_INT_MAPLECONVERSION) {
                    switch (prop.val) {
                    case _GT_WEIGHTS:
                        if (val.type==_VECT) {
                            if (!val._VECTptr->empty()) {
                                if (ckmatrix(val._VECTptr->front())) {
                                    if (val._VECTptr->size()!=topology.size()-1)
                                        return generrdim("Invalid number of initial weight matrices");
                                    init_weights.reserve(val._VECTptr->size());
                                    for (const_iterateur wt=val._VECTptr->begin();wt!=val._VECTptr->end();++wt) {
                                        if (wt!=val._VECTptr->begin() && !ckmatrix(*wt))
                                            return generrtype("Expected an initial weight matrix");
                                        init_weights.push_back(wt->_VECTptr);
                                    }
                                } else winit=makevecteur(val._VECTptr->front(),vecteur(val._VECTptr->begin()+1,val._VECTptr->end()));
                            }
                        } else winit=makevecteur(val,vecteur(0));
                        break;
                    case _ANN_LEARNING_RATE:
                        if (val.type==_VECT) {
                            if (val._VECTptr->size()!=2)
                                return generr("Invalid scheduled learning rate specification");
                            if (!val._VECTptr->back().is_symb_of_sommet(at_program) ||
                                    val._VECTptr->back()._SYMBptr->feuille._VECTptr->front()._VECTptr->size()!=1)
                                return generr("Schedule multiplier must be an univariate function");
                            rate=val._VECTptr->front();
                            schedule=val._VECTptr->back();
                        } else rate=val;
                        if (!is_real_number(rate,contextptr) || !is_strictly_positive(rate=to_real_number(rate,contextptr),contextptr))
                            return generr("(Initial) learning rate must be a positive real number");
                        set_rate=true;
                        break;
                    case _ANN_BLOCK_SIZE:
                        if (!val.is_integer() || (bs=val.val)<1)
                            return generr("Block size must be a positive integer");
                        break;
                    case _ANN_MOMENTUM:
                        if (val.is_integer() && val.subtype==_INT_PLOT && val.val==_ADAPTIVE) {
                            mom=1;
                        } else if (val.type==_VECT) {
                            if (val._VECTptr->size()!=2 || !is_numericv(*val._VECTptr,num_mask_withfrac) ||
                                    is_positive(-_min(val,contextptr),contextptr) || is_greater(_max(val,contextptr),1,contextptr))
                                return generr("Invalid specification of adaptive momentum parameters");
                            mom=_evalf(val,contextptr);
                        } else if (is_real_number(val,contextptr)) {
                            mom=to_real_number(val,contextptr);
                            if (!is_positive(mom,contextptr) || is_greater(mom,1,contextptr))
                                return generr("Momentum parameter must be a real number in [0,1)");
                        }
                        break;
                    case _ANN_WEIGHT_DECAY:
                        reg=val.type==_VECT?*val._VECTptr:vecteur(1,val);
                        if (!is_numericv(reg,num_mask_withfrac|num_mask_withint) || is_strictly_positive(-_min(reg,contextptr),contextptr))
                            return generr("Invalid weight decay parameter(s) specification");
                        break;
                    default:
                        return generr("Invalid property specification");
                    }
                } else if (prop.is_integer() && prop.subtype==_INT_TYPE && prop.val==_FUNC) {
                    act=val.type==_VECT?*val._VECTptr:vecteur(1,val);
                    if (act.empty())
                        return generr("Invalid activation function specification");
                    vecteur params(act.begin()+1,act.end());
                    act=makevecteur(act.front(),params);
                } else if (prop.is_integer() && prop.subtype==_INT_PLOT) {
                    if (prop.val==_LABELS) {
                        if (val.type!=_VECT)
                            return generrtype("Expected a list of labels");
                        lab=*val._VECTptr;
                    } else if (prop.val==_TITLE) {
                        if (val.type!=_STRNG)
                            return generrtype("Title must be a string");
                        name=*val._STRNGptr;
                    } else return generr("Invalid property specification");
                } else return generr("Invalid property specification");
            } else return generrtype("Invalid optional argument");
        }
    } else topology=*g._VECTptr;
    vector<int> t;
    if (src==NULL) {
        if (topology.size()<2)
            return generrdim("Too few layers in the network");
        if (!is_integer_vecteur(topology))
            return generrtype("Invalid layer size type");
        if (is_strictly_greater(1,_min(topology,contextptr),contextptr))
            return gensizeerr("Invalid layer size, must be at least 1");
        if (topology.size()==2) { // guess the size of the hidden layer: (input+output)*2/3
            gen hsize=_round(2*_sum(topology,contextptr)/3,contextptr);
            *logptr(contextptr) << gettext("No hidden layers specified, creating 1 hidden layer of size") << " " << hsize << "\n";
            topology.insert(topology.begin()+1,hsize);
        }
        vecteur2vector_int(topology,0,t);
    }
    if (reg.size()==1) {
        gen r=reg.front();
        reg=vecteur((src==NULL?topology.size():src->layer_count())-1,r);
    }
    gen glorot=inv(sqrt(fan_in+fan_out,contextptr),contextptr),lecun=inv(sqrt(fan_in,contextptr),contextptr);
    vecteur fan=makevecteur(fan_in,fan_out);
    vecteur glorot_uniform=makevecteur(symbolic(at_uniformd,makesequence(minus_sqrt6*glorot,plus_sqrt6*glorot)),fan);
    vecteur glorot_normal=makevecteur(symbolic(at_normald,makesequence(0,plus_sqrt2*glorot)),fan);
    vecteur he_uniform=makevecteur(symbolic(at_uniformd,makesequence(minus_sqrt6*lecun,plus_sqrt6*lecun)),fan);
    vecteur he_normal=makevecteur(symbolic(at_normald,makesequence(0,plus_sqrt2*lecun)),fan);
    vecteur lecun_uniform=makevecteur(symbolic(at_uniformd,makesequence(minus_sqrt3*lecun,plus_sqrt3*lecun)),fan);
    vecteur lecun_normal=makevecteur(symbolic(at_normald,makesequence(0,lecun)),fan);
    if (src==NULL && winit.empty() && init_weights.empty()) // set default initialization
        winit=makevecteur(symbolic(at_uniformd,makesequence(-lecun,lecun)),makevecteur(fan_in,fan_out));
    try {
        log_output_redirect lor(contextptr);
        ann network=src==NULL?ann(t,bs,contextptr):ann(*src,bs);
        if (src==NULL || set_rate)
            network.set_learning_rate(rate.to_double(contextptr),schedule);
        if (mom.type==_VECT) {
            network.set_momentum(1);
            network.set_adam_params(mom._VECTptr->front().to_double(contextptr),mom._VECTptr->back().to_double(contextptr));
        } else if (src==NULL) network.set_momentum(mom.to_double(contextptr));
        if (name.length()>0)
            network.set_name(name);
        if (!reg.empty())
            network.set_regularization(reg);
        // set labels
        if (!lab.empty())
            network.set_labels(lab,cls);
        int task=network.task();
        if (task<0)
            return generr("Bad output labels");
        // set activation function for hidden layers
        if (!act.empty())
            network.set_activation(false,act.front(),*act.back()._VECTptr);
        else if (src==NULL) switch (task) {
        case 0:
            network.set_activation(false,change_subtype(_ANN_RELU,_INT_MAPLECONVERSION));
            break;
        case 1:
        case 2:
            network.set_activation(false,at_tanh);
            break;
        default: assert(false);
        }
        // set activation function for the output layer
        if (!out_act.empty()) {
            if (cls)
                *logptr(contextptr) << gettext("Warning") << ": "
                                    << gettext("'output' option is ignored in classifier networks") << "\n";
            network.set_activation(true,out_act.front(),*out_act.back()._VECTptr);
        } else if (src==NULL) switch (task) {
        case 0:
            network.set_activation(true,at_id);
            break;
        case 1:
            if (network.output_size()==1)
                network.set_activation(true,at_logistic);
            else ; // softmax will be used
            break;
        case 2:
            network.set_activation(true,at_logistic);
            break;
        default: assert(false);
        }
        // set error function
        if (!err.empty())
            network.set_error_function(err.front(),*err.back()._VECTptr);
        else if (src==NULL) switch (task) {
        case 0:
            network.set_error_function(change_subtype(_ANN_HALF_MSE,_INT_MAPLECONVERSION));
            break;
        case 1:
            if (network.output_size()>1)
                network.set_error_function(change_subtype(_ANN_CROSS_ENTROPY,_INT_MAPLECONVERSION));
            else network.set_error_function(change_subtype(_ANN_LOG_LOSS,_INT_MAPLECONVERSION));
            break;
        case 2:
            network.set_error_function(change_subtype(_ANN_LOG_LOSS,_INT_MAPLECONVERSION));
            break;
        default: assert(false);
        }
        // initialize weights
        if (!init_weights.empty()) {
            i=0;
            for (vector<matrice*>::const_iterator it=init_weights.begin();it!=init_weights.end();++it,++i)
                network.set_initial_weights(*(*it),i);
        } else if (!winit.empty()) {
            if (winit.front().type==_STRNG) {
                if (!winit.back()._VECTptr->empty())
                    *logptr(contextptr) << gettext("Warning") << ": " << gettext("ignoring parameters for weight initializer") << "\n";
                const string &str=*winit.front()._STRNGptr;
                if (str=="glorot-uniform")
                    winit=glorot_uniform;
                else if (str=="he-uniform")
                    winit=he_uniform;
                else if (str=="lecun-uniform")
                    winit=lecun_uniform;
                else if (str=="glorot-normal")
                    winit=glorot_normal;
                else if (str=="he-normal")
                    winit=he_normal;
                else if (str=="lecun-normal")
                    winit=lecun_normal;
                else return generr("Unknown weight initializer");
            }
            network.randomize_initial_weights(winit.front(),*winit.back()._VECTptr);
        }
        return network;
    } catch (const std::runtime_error &e) {
        return generr(e.what(),false);
    }
#else
    return generr("GNU Scientific Library is required for neural networks");
#endif // HAVE_LIBGSL
}
static const char _neural_network_s []="neural_network";
static define_unary_function_eval (__neural_network,&_neural_network,_neural_network_s);
define_unary_function_ptr5(at_neural_network,alias_at_neural_network,&__neural_network,0,true)

gen _train(const gen &g,GIAC_CONTEXT) {
    if (g.type==_STRNG && g.subtype==-1) return g;
#ifdef HAVE_LIBGSL
    if (g.type!=_VECT || g.subtype!=_SEQ__VECT)
        return gentypeerr(contextptr);
    const vecteur &gv=*g._VECTptr;
    if (gv.size()<3)
        return generrdim("Expected 3 or more arguments");
    ann *net=ann::from_gen(gv.front());
    if (net==NULL)
        return generrtype("First argument must be a neural network");
    const gen &inp=gv[1],&out=gv[2];
    if (!ckmatrix(inp) || out.type!=_VECT)
        return generrtype("Input and output must be a matrix and a vector");
    int batch_size=0;
    if (gv.size()>3 && (!gv[3].is_integer() || (batch_size=gv[3].val)>(int)inp._VECTptr->size()))
        return generr("Invalid batch size specification");
    try {
        net->train(*inp._VECTptr,*out._VECTptr,batch_size);
    } catch (const std::runtime_error &e) {
        return generr(e.what(),false);
    }
    return *net;
#else
    return generr("GNU Scientific Library is required for neural networks");
#endif // HAVE_LIBGSL
}
static const char _train_s []="train";
static define_unary_function_eval (__train,&_train,_train_s);
define_unary_function_ptr5(at_train,alias_at_train,&__train,0,true)

#ifndef NO_NAMESPACE_GIAC
}
#endif // ndef NO_NAMESPACE_GIAC
