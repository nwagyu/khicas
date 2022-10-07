/*
 * This file is part of the MicroPython project, http://micropython.org/
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2013, 2014 Damien P. George and 2017, 2018 Rami Ali
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "py/compile.h"
#include "py/runtime.h"
#include "py/repl.h"
#include "py/gc.h"
#include "py/mperrno.h"
#include "lib/utils/pyexec.h"

#include "library.h"
#if defined EMCC && !defined NO_QSTR
#include <emscripten.h>
#endif

#if MICROPY_ENABLE_COMPILER
int do_str(const char *src, mp_parse_input_kind_t input_kind,bool is_repl) {
    int ret = 0;
    nlr_buf_t nlr;
    if (nlr_push(&nlr) == 0) {
        mp_lexer_t *lex = mp_lexer_new_from_str_len(MP_QSTR__lt_stdin_gt_, src, strlen(src), 0);
        qstr source_name = lex->source_name;
        mp_parse_tree_t parse_tree = mp_parse(lex, input_kind);
        mp_obj_t module_fun = mp_compile(&parse_tree, source_name, is_repl);
        mp_call_function_0(module_fun);
        nlr_pop();
    } else {
        // uncaught exception
        if (mp_obj_is_subclass_fast(mp_obj_get_type((mp_obj_t)nlr.ret_val), &mp_type_SystemExit)) {
            mp_obj_t exit_val = mp_obj_exception_get_value(MP_OBJ_FROM_PTR(nlr.ret_val));
            if (exit_val != mp_const_none) {
                mp_int_t int_val;
                if (mp_obj_get_int_maybe(exit_val, &int_val)) {
                    ret = int_val & 255;
                } else {
                    ret = 1;
                }
            }
        } else {
            mp_obj_print_exception(&mp_plat_print, (mp_obj_t)nlr.ret_val);
            ret = 1;
        }
    }
    return ret;
}
#endif

static char *stack_top;

int mp_js_do_str(const char *code) {
  EM_ASM_ARGS({
      var msg = UTF8ToString($0);
      console.log(msg);
    },code);
  if (strncmp(code,"show",4)==0)
    return do_str(code, MP_PARSE_FILE_INPUT,false);
  const char * s=code;
  for (;*s;++s){
    if (*s=='\n')
      return do_str(code, MP_PARSE_FILE_INPUT,false);
  }
  return do_str(code, MP_PARSE_SINGLE_INPUT,true);
}

int mp_js_process_char(int c) {
    return pyexec_event_repl_process_char(c);
}

char * mp_js_init(int heap_size) {
    int stack_dummy;
    stack_top = (char*)&stack_dummy;
    char * ptr=0;
    #if MICROPY_ENABLE_GC
    char *heap = (char*)malloc(heap_size * sizeof(char));
    ptr=heap;
    gc_init(heap, heap + heap_size);
    #endif

    #if MICROPY_ENABLE_PYSTACK
    static mp_obj_t pystack[1024];
    mp_pystack_init(pystack, &pystack[MP_ARRAY_SIZE(pystack)]);
    #endif

    mp_init();
    
    mp_obj_list_init(mp_sys_path, 0);
    mp_obj_list_append(mp_sys_path, MP_OBJ_NEW_QSTR(MP_QSTR_));
    mp_obj_list_init(mp_sys_argv, 0);
    return ptr;
}

void mp_js_init_repl() {
    pyexec_event_repl_init();
}

void gc_collect(void) {
    // WARNING: This gc_collect implementation doesn't try to get root
    // pointers from CPU registers, and thus may function incorrectly.
    jmp_buf dummy;
    if (setjmp(dummy) == 0) {
        longjmp(dummy, 1);
    }
    gc_collect_start();
    gc_collect_root((void*)stack_top, ((mp_uint_t)(void*)(&dummy + 1) - (mp_uint_t)stack_top) / sizeof(mp_uint_t));
    gc_collect_end();
}

#if 0 // in py/lexer.c
const char numpy_script[]="(import linalg\nimport math\nclass array: \n    def __init__(self, a): \n        self.a = a \n  \n    def __add__(self, other): \n        return array(linalg.add(self.a , other.a))\n  \n    def __sub__(self, other): \n        return array(linalg.sub(self.a , other.a))\n  \n    def __mul__(self, other):\n        if type(self)==array:\n            if type(other)==array:\n                return array(linalg.mul(self.a , other.a))\n            return array(linalg.mul(self.a,other))\n        return array(linalg.mul(self,other.a))\n    \n    def __rmul__(self, other): \n        if type(self)==array:\n            if type(other)==array:\n                return array(linalg.mul(self.a , other.a))\n            return array(linalg.mul(self.a,other))\n        return array(linalg.mul(self,other.a))\n\n    def __matmul__(self, other):\n        return __mul(self,other)\n\n    def __getitem__(self,key):\n        r=(self.a)[key]\n        if type(r)==list or type(r)==tuple:\n            return array(r)\n        return r\n\n    def __setitem__(self, key, value):\n        if (type(value)==array):\n            (self.a)[key]=value.a\n        else:\n            (self.a)[key]=value\n        return None\n\n    def __len__(self):\n        return len(self.a)\n    \n    def __str__(self): \n        return 'array('+str(self.a)+')'\n  \n    def __repr__(self): \n        return 'array('+str(self.a)+')'\n  \n    def __neg__(self):\n        return array(-self.a)\n\n    def __pos__(self):\n        return self\n    \n    def __abs__(self):\n        return array(linalg.abs(self.a))\n\n    def __round__(self):\n        return array(linalg.apply(round,self.a,linalg.matrix))\n\n    def __trunc__(self):\n        return array(linalg.apply(trunc,self.a,linalg.matrix))\n\n    def __floor__(self):\n        return array(linalg.apply(floor,self.a,linalg.matrix))\n\n    def __ceil__(self):\n        return array(linalg.apply(ceil,self.a,linalg.matrix))\n\n    def T(self):\n        return array(linalg.transpose(self.a))\n            \ndef real(x):\n    if type(x)==array:\n        return array(linalg.re(x.a))\n    return x.real\n\ndef imag(x):\n    if type(x)==array:\n        return array(linalg.im(x.a))\n    return x.imag\n\ndef conj(x):\n    if type(x)==array:\n        return array(linalg.conj(x.a))\n    return linalg.conj(x)\n\ndef sin(x):\n    if type(x)==array:\n        return array(linalg.apply(math.sin,x.a,linalg.matrix))\n    return math.sin(x)\n\ndef cos(x):\n    if type(x)==array:\n        return array(linalg.apply(math.cos,x.a,linalg.matrix))\n    return math.cos(x)\n\ndef tan(x):\n    if type(x)==array:\n        return array(linalg.apply(math.tan,x.a,linalg.matrix))\n    return math.tan(x)\n\ndef asin(x):\n    if type(x)==array:\n        return array(linalg.apply(math.asin,x.a,linalg.matrix))\n    return math.asin(x)\n\ndef acos(x):\n    if type(x)==array:\n        return array(linalg.apply(math.acos,x.a,linalg.matrix))\n    return math.acos(x)\n\ndef atan(x):\n    if type(x)==array:\n        return array(linalg.apply(math.atan,x.a,linalg.matrix))\n    return math.atan(x)\n\ndef sinh(x):\n    if type(x)==array:\n        return array(linalg.apply(math.sinh,x.a,linalg.matrix))\n    return math.sinh(x)\n\ndef cosh(x):\n    if type(x)==array:\n        return array(linalg.apply(math.cosh,x.a,linalg.matrix))\n    return math.cosh(x)\n\ndef tanh(x):\n    if type(x)==array:\n        return array(linalg.apply(math.tanh,x.a,linalg.matrix))\n    return math.tanh(x)\n\ndef exp(x):\n    if type(x)==array:\n        return array(linalg.apply(math.exp,x.a,linalg.matrix))\n    return math.exp(x)\n\ndef log(x):\n    if type(x)==array:\n        return array(linalg.apply(math.log,x.a,linalg.matrix))\n    return math.log(x)\n\ndef size(x):\n    if type(x)==array:\n        return linalg.size(x.a)\n    return linalg.size(x)\n\ndef shape(x):\n    if type(x)==array:\n        return linalg.shape(x.a)\n\ndef dot(a,b):\n    return a*b\n\ndef transpose(a):\n    if type(x)==array:\n        return array(linalg.transpose(x.a))\n\ndef trn(a):\n    if type(x)==array:\n        return array(linalg.conj(linalg.transpose(x.a)))\n    return linalg.conj(linalg.transpose(x.a))\n\ndef zeros(n,m=0):\n    return array(linalg.zeros(n,m))\n\ndef ones(n,m=0):\n    return array(linalg.ones(n,m))\n\ndef eye(n):\n    return array(linalg.eye(n))\n\ndef det(x):\n    if type(x)==array:\n        return linalg.det(x.a)\n    return linalg.det(x)\n\ndef inv(x):\n    if type(x)==array:\n        return array(linalg.inv(x.a))\n    return linalg.inv(x)\n\ndef solve(a,b):\n    if type(a)==array:\n        if type(b)==array:\n            return array(linalg.solve(a.a,b.a))\n        return array(linalg.solve(a.a,b))\n    if type(b)==array:\n        return array(linalg.solve(a,b.a))\n    return linalg.solve(a,b)\n\ndef eig(a):\n    if type(a)==array:\n        r=linalg.eig(a.a)\n        return array(r[0]),array(r[1])\n    return linalg.eig(a)\n\ndef linspace(a,b,c):\n    return array(linalg.linspace(a,b,c))\n\ndef arange(a,b,c=1):\n    return array(linalg.arange(a,b,c))\n\ndef reshape(a,n,m):\n    if type(n)==tuple:\n        m=n[1]\n        n=n[0]\n    if type(a)==array:\n        return array(linalg.matrix(n,m,a.a))\n    return linalg.matrix(n,m,a)\n)";

mp_lexer_t * mp_lexer_new_from_file(const char * filename) {
  if (strcmp(filename,"numpy.py")==0)
    return mp_lexer_new_from_str_len(qstr_from_str(filename), numpy_script, strlen(numpy_script), 0 /* size_t free_len*/);
  else
    mp_raise_OSError(MP_ENOENT);
}
#endif

mp_import_stat_t mp_import_stat(const char *path) {
  if (strcmp(path,"numpy.py")==0) {
    return MP_IMPORT_STAT_FILE;
  }
  return MP_IMPORT_STAT_NO_EXIST;
}

mp_obj_t mp_builtin_open(size_t n_args, const mp_obj_t *args, mp_map_t *kwargs) {
    return mp_const_none;
}
MP_DEFINE_CONST_FUN_OBJ_KW(mp_builtin_open_obj, 1, mp_builtin_open);

void nlr_jump_fail(void *val) {
    while (1);
}

void NORETURN __fatal_error(const char *msg) {
    while (1);
}

#ifndef NDEBUG
void MP_WEAK __assert_func(const char *file, int line, const char *func, const char *expr) {
    printf("Assertion '%s' failed, at file %s:%d\n", expr, file, line);
    __fatal_error("Assertion failed");
}
#endif
