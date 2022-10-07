/*
 * signalprocessing.h
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
#ifndef __SIGNALPROCESSING_H
#define __SIGNALPROCESSING_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "first.h"
#include "gen.h"
#include "unary.h"
#ifdef HAVE_LIBGSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

enum filter_type {
    _LOWPASS_FILTER,
    _HIGHPASS_FILTER
};

enum FourierFuncType {
    _FOURIER_FUNCTYPE_UNKNOWN       =0,
    _FOURIER_FUNCTYPE_ONE           =1,
    _FOURIER_FUNCTYPE_DIRAC         =2,
    _FOURIER_FUNCTYPE_INV_MONOM     =3,
    _FOURIER_FUNCTYPE_GAUSSIAN      =4,
    _FOURIER_FUNCTYPE_LOGABSX       =5,
    _FOURIER_FUNCTYPE_BESSELJ       =6,
    _FOURIER_FUNCTYPE_SGN           =7,
    _FOURIER_FUNCTYPE_ABSX_ALPHA    =8,
    _FOURIER_FUNCTYPE_INVABSX       =9,
    _FOURIER_FUNCTYPE_GAMMA         =10,
    _FOURIER_FUNCTYPE_EXPEXP        =11,
    _FOURIER_FUNCTYPE_AIRY_AI       =12,
    _FOURIER_FUNCTYPE_EXP_HEAVISIDE =13,
    _FOURIER_FUNCTYPE_HEAVISIDE     =14,
    _FOURIER_FUNCTYPE_ATAN_OVERX    =15,
    _FOURIER_FUNCTYPE_EXPABSX       =16,
    _FOURIER_FUNCTYPE_EXPABSX_OVERX =17,
    _FOURIER_FUNCTYPE_COSECH        =18,
    _FOURIER_FUNCTYPE_SECH          =19,
    _FOURIER_FUNCTYPE_SECH_2        =20,
    _FOURIER_FUNCTYPE_CONVOLUTION   =21,
    _FOURIER_FUNCTYPE_PIECEWISE     =22,
    _FOURIER_FUNCTYPE_FUNC          =23,
    _FOURIER_FUNCTYPE_DIFF          =24,
    _FOURIER_FUNCTYPE_PARTIAL_DIFF  =25
};

enum ann_activation_function {
    _ANN_LINEAR     =1,
    _ANN_SIGMOID    =2,
    _ANN_TANH       =3
};

bool is_sound_data(const gen &g,int &nc,int &bd,int &sr,int &len);
vecteur decode_chdata(const vecteur &data,int bd,int start=0,int len=-1);
vecteur encode_chdata(const vecteur &data,int bd,double ratio,GIAC_CONTEXT);
gen filter(const vecteur &args,filter_type typ,GIAC_CONTEXT);
bool is_real_number(const gen &g,GIAC_CONTEXT);
gen to_real_number(const gen &g,GIAC_CONTEXT);
gen generr(const char* msg,bool translate=true);
gen generrtype(const char* msg,bool translate=true);
gen generrdim(const char* msg,bool translate=true);
void print_error(const char *msg,GIAC_CONTEXT);
void print_warning(const char *msg,GIAC_CONTEXT);
bool has_rootof(const gen &g);

gen _createwav(const gen &g,GIAC_CONTEXT);
gen _stereo2mono(const gen &g,GIAC_CONTEXT);
gen _plotwav(const gen &g,GIAC_CONTEXT);
gen _plotspectrum(const gen &g,GIAC_CONTEXT);
gen _channels(const gen &g,GIAC_CONTEXT);
gen _bit_depth(const gen &g,GIAC_CONTEXT);
gen _samplerate(const gen &g,GIAC_CONTEXT);
gen _duration(const gen &g,GIAC_CONTEXT);
gen _channel_data(const gen &g,GIAC_CONTEXT);
gen _cross_correlation(const gen &g,GIAC_CONTEXT);
gen _auto_correlation(const gen &g,GIAC_CONTEXT);
gen _lowpass(const gen &g,GIAC_CONTEXT);
gen _highpass(const gen &g,GIAC_CONTEXT);
gen _moving_average(const gen &g,GIAC_CONTEXT);
gen _rms(const gen &g,GIAC_CONTEXT);
gen _resample(const gen &g,GIAC_CONTEXT);
gen _convolution(const gen &g,GIAC_CONTEXT);
gen _threshold(const gen &g,GIAC_CONTEXT);
gen _bartlett_hann_window(const gen &g,GIAC_CONTEXT);
gen _blackman_harris_window(const gen &g,GIAC_CONTEXT);
gen _blackman_window(const gen &g,GIAC_CONTEXT);
gen _bohman_window(const gen &g,GIAC_CONTEXT);
gen _cosine_window(const gen &g,GIAC_CONTEXT);
gen _gaussian_window(const gen &g,GIAC_CONTEXT);
gen _hamming_window(const gen &g,GIAC_CONTEXT);
gen _hann_poisson_window(const gen &g,GIAC_CONTEXT);
gen _hann_window(const gen &g,GIAC_CONTEXT);
gen _parzen_window(const gen &g,GIAC_CONTEXT);
gen _poisson_window(const gen &g,GIAC_CONTEXT);
gen _riemann_window(const gen &g,GIAC_CONTEXT);
gen _triangle_window(const gen &g,GIAC_CONTEXT);
gen _tukey_window(const gen &g,GIAC_CONTEXT);
gen _welch_window(const gen &g,GIAC_CONTEXT);
gen _fourier(const gen &g,GIAC_CONTEXT);
gen _ifourier(const gen &g,GIAC_CONTEXT);
gen _addtable(const gen &g,GIAC_CONTEXT);
gen _rect(const gen &g,GIAC_CONTEXT);
gen _boxcar(const gen &g,GIAC_CONTEXT);
gen _tri(const gen &g,GIAC_CONTEXT);
gen _sinc(const gen &g,GIAC_CONTEXT);
gen _logistic(const gen &g,GIAC_CONTEXT);
gen _neural_network(const gen &g,GIAC_CONTEXT);
gen _train(const gen &g,GIAC_CONTEXT);

extern const unary_function_ptr * const at_createwav;
extern const unary_function_ptr * const at_plotwav;
extern const unary_function_ptr * const at_plotspectrum;
extern const unary_function_ptr * const at_channels;
extern const unary_function_ptr * const at_samplerate;
extern const unary_function_ptr * const at_bit_depth;
extern const unary_function_ptr * const at_duration;
extern const unary_function_ptr * const at_channel_data;
extern const unary_function_ptr * const at_cross_correlation;
extern const unary_function_ptr * const at_auto_correlation;
extern const unary_function_ptr * const at_lowpass;
extern const unary_function_ptr * const at_highpass;
extern const unary_function_ptr * const at_convolution;
extern const unary_function_ptr * const at_threshold;
extern const unary_function_ptr * const at_bartlett_hann_window;
extern const unary_function_ptr * const at_blackman_harris_window;
extern const unary_function_ptr * const at_blackman_window;
extern const unary_function_ptr * const at_bohman_window;
extern const unary_function_ptr * const at_cosine_window;
extern const unary_function_ptr * const at_gaussian_window;
extern const unary_function_ptr * const at_hamming_window;
extern const unary_function_ptr * const at_hann_poisson_window;
extern const unary_function_ptr * const at_hann_window;
extern const unary_function_ptr * const at_parzen_window;
extern const unary_function_ptr * const at_poisson_window;
extern const unary_function_ptr * const at_riemann_window;
extern const unary_function_ptr * const at_triangle_window;
extern const unary_function_ptr * const at_tukey_window;
extern const unary_function_ptr * const at_welch_window;
extern const unary_function_ptr * const at_fourier;
extern const unary_function_ptr * const at_ifourier;
extern const unary_function_ptr * const at_addtable;
extern const unary_function_ptr * const at_rect;
extern const unary_function_ptr * const at_boxcar;
extern const unary_function_ptr * const at_tri;
extern const unary_function_ptr * const at_sinc;
extern const unary_function_ptr * const at_logistic;
extern const unary_function_ptr * const at_neural_network;
extern const unary_function_ptr * const at_train;

/* Neural network class */
#ifdef HAVE_LIBGSL
class ann : public gen_user {
    std::vector<int> topology;
    std::vector<gsl_matrix*> layers;
    std::vector<gsl_matrix*> cache;
    std::vector<gsl_matrix*> deltas;
    std::vector<gsl_matrix*> weights;
    std::vector<gsl_matrix*> Deltas;
    std::vector<gsl_matrix*> moment1;
    std::vector<gsl_matrix*> moment2;
    const context *ctx;
    int block_size,iter,accum_count;
    double learning_rate,eps,momentum,beta1,beta2;
    gen x,schedule;
    vecteur y,t,fan,labels,activation,output_activation,errfunc;
    bool _is_classifier;
    gsl_vector *reg_coeff;
    gsl_matrix *eout;
    std::string title;
    void create();
    void assign(const ann &other,int bs);
    void deallocate();
    void propagate_forward(bool comp_deriv=false) const;
    void propagate_backward();
    void calc_deltas() const;
    void update_Deltas();
    void update_weights();
    void activate(int i,double &a,double *b) const;
    gen compute_errfunc(gsl_vector *t,gsl_vector *y) const;
    void compute_errfunc_diff(gsl_matrix *t,gsl_matrix *y,gsl_matrix *res) const;
    void ckinput(const matrice &input,const vecteur &output,bool ignore_output) const;
    void output2gsl_matrix(const vecteur &output,int i0,int n,gsl_matrix *m) const;
    gen gsl_vector2output(gsl_vector *output) const;
    void feed(const matrice &input,vecteur &res,const matrice &expected_output=vecteur(0)) const;
    std::string network_type() const;
    void weights2vecteur(vecteur &res) const;
public:
    ann(const std::vector<int> &topology,int bsize,GIAC_CONTEXT);
    ann(const ann &other);
    ann(const ann &other,int bs);
    ann &operator=(const ann &other);
    static ann *from_gen(const gen &g);
    // reimplemented from gen_user
    virtual ~ann() { deallocate(); }
    virtual std::string print (GIAC_CONTEXT) const;
    virtual std::string texprint (GIAC_CONTEXT) const;
    virtual gen giac_constructor (GIAC_CONTEXT) const;
    virtual bool operator ==(const gen &) const;
    virtual gen_user *memory_alloc() const { return new ann(*this); }
    virtual gen operator()(const gen &,GIAC_CONTEXT) const;
    virtual gen operator[](const gen &);
    // methods
    gsl_matrix *const layer(int i) const { assert(i>=0&&i<layer_count()); return layers[i]; }
    gsl_matrix *const weight_matrix(int i) const { assert(i>=0&&i<layer_count()-1); return weights[i]; }
    gsl_matrix *const weight_updates(int i) const { assert(i>=0&&i<layer_count()-1); return Deltas[i]; }
    vecteur layout() const { return vector_int_2_vecteur(topology); }
    int layer_count() const { return topology.size(); }
    int layer_size(int i) const { assert(i>=0&&i<layer_count()); return topology[i]; }
    int input_size() const { return layer_size(0); }
    int output_size() const { return layer_size(layer_count()-1); }
    int hidden_layer_count() const { return layer_count()-2; }
    int neuron_count() const;
    const std::string &name() const { return title; }
    bool is_classifier() const { return _is_classifier; }
    bool is_default_classifier() const;
    void set_error_function(const gen &f,const vecteur &params=vecteur(0));
    void set_activation(bool out,const gen &f,const vecteur &params=vecteur(0));
    void set_regularization(const vecteur &reg) const;
    void randomize_initial_weights(const gen &f,const vecteur &fan_in_out) const;
    void set_initial_weights(const matrice &w,int i=0) const;
    void set_labels(const vecteur &lab,bool classes=false);
    void set_name(const std::string &name) { title=name; }
    void set_learning_rate(double lr,const gen &sched) { learning_rate=lr; schedule=sched; }
    void set_momentum(double mom);
    void set_adam_params(double b1,double b2) { beta1=b1; beta2=b2; }
    void restart_iteration();
    int task() const { return labels.empty()?0:(_is_classifier?1:2); }
    const gen &label(int i) const { assert(i>=0&&i<topology.back()); return labels[i]; }
    void train(const matrice &input,const vecteur &expected_output,int batch_size);
};
#endif

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
#endif // __SIGNALPROCESSING_H
