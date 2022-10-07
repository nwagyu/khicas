/*
 * This file is part of the MicroPython project, http://micropython.org/
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2013, 2014 Damien P. George
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
#include <stdbool.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "py/mpstate.h"
#include "py/nlr.h"
#include "py/compile.h"
#include "py/runtime.h"
#include "py/builtin.h"
#include "py/repl.h"
#include "py/gc.h"
#include "py/stackctrl.h"
#include "py/mphal.h"
#include "py/mpthread.h"
#include "py/mperrno.h"
#include "extmod/misc.h"
#include "genhdr/mpversion.h"
#include "input.h"
#ifndef PATH_MAX
#define PATH_MAX 128
#endif

void console_output(const char *,int);
const char * read_file(const char * filename);
bool file_exists(const char * filename);
int ctrl_c_interrupted();

int micropython_port_vm_hook_loop() {
  static int c = 0;

  c++;
  if (c & 0x1f) {
    return 0;
  }
  if (ctrl_c_interrupted()){
    nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Keyboard interrupt"));
    return 1;
  }
  return 0;
}

#if 0
const char numpy_script[]=R"(import linalg
import math
class array: 
    def __init__(self, a): 
        self.a = a 
  
    def __add__(self, other): 
        return array(linalg.add(self.a , other.a))
  
    def __sub__(self, other): 
        return array(linalg.sub(self.a , other.a))
  
    def __mul__(self, other):
        if type(self)==array:
            if type(other)==array:
                return array(linalg.mul(self.a , other.a))
            return array(linalg.mul(self.a,other))
        return array(linalg.mul(self,other.a))
    
    def __rmul__(self, other): 
        if type(self)==array:
            if type(other)==array:
                return array(linalg.mul(self.a , other.a))
            return array(linalg.mul(self.a,other))
        return array(linalg.mul(self,other.a))

    def __matmul__(self, other):
        return __mul(self,other)

    def __getitem__(self,key):
        r=(self.a)[key]
        if type(r)==list or type(r)==tuple:
            return array(r)
        return r

    def __setitem__(self, key, value):
        if (type(value)==array):
            (self.a)[key]=value.a
        else:
            (self.a)[key]=value
        return None

    def __len__(self):
        return len(self.a)
    
    def __str__(self): 
        return 'array('+str(self.a)+')'
  
    def __repr__(self): 
        return 'array('+str(self.a)+')'
  
    def __neg__(self):
        return array(-self.a)

    def __pos__(self):
        return self
    
    def __abs__(self):
        return array(linalg.abs(self.a))

    def __round__(self):
        return array(linalg.apply(round,self.a,linalg.matrix))

    def __trunc__(self):
        return array(linalg.apply(trunc,self.a,linalg.matrix))

    def __floor__(self):
        return array(linalg.apply(floor,self.a,linalg.matrix))

    def __ceil__(self):
        return array(linalg.apply(ceil,self.a,linalg.matrix))

    def T(self):
        return array(linalg.transpose(self.a))
            
def real(x):
    if type(x)==array:
        return array(linalg.re(x.a))
    return x.real

def imag(x):
    if type(x)==array:
        return array(linalg.im(x.a))
    return x.imag

def conj(x):
    if type(x)==array:
        return array(linalg.conj(x.a))
    return linalg.conj(x)

def sin(x):
    if type(x)==array:
        return array(linalg.apply(math.sin,x.a,linalg.matrix))
    return math.sin(x)

def cos(x):
    if type(x)==array:
        return array(linalg.apply(math.cos,x.a,linalg.matrix))
    return math.cos(x)

def tan(x):
    if type(x)==array:
        return array(linalg.apply(math.tan,x.a,linalg.matrix))
    return math.tan(x)

def asin(x):
    if type(x)==array:
        return array(linalg.apply(math.asin,x.a,linalg.matrix))
    return math.asin(x)

def acos(x):
    if type(x)==array:
        return array(linalg.apply(math.acos,x.a,linalg.matrix))
    return math.acos(x)

def atan(x):
    if type(x)==array:
        return array(linalg.apply(math.atan,x.a,linalg.matrix))
    return math.atan(x)

def sinh(x):
    if type(x)==array:
        return array(linalg.apply(math.sinh,x.a,linalg.matrix))
    return math.sinh(x)

def cosh(x):
    if type(x)==array:
        return array(linalg.apply(math.cosh,x.a,linalg.matrix))
    return math.cosh(x)

def tanh(x):
    if type(x)==array:
        return array(linalg.apply(math.tanh,x.a,linalg.matrix))
    return math.tanh(x)

def exp(x):
    if type(x)==array:
        return array(linalg.apply(math.exp,x.a,linalg.matrix))
    return math.exp(x)

def log(x):
    if type(x)==array:
        return array(linalg.apply(math.log,x.a,linalg.matrix))
    return math.log(x)

def size(x):
    if type(x)==array:
        return linalg.size(x.a)
    return linalg.size(x)

def shape(x):
    if type(x)==array:
        return linalg.shape(x.a)

def dot(a,b):
    return a*b

def transpose(a):
    if type(x)==array:
        return array(linalg.transpose(x.a))

def trn(a):
    if type(x)==array:
        return array(linalg.conj(linalg.transpose(x.a)))
    return linalg.conj(linalg.transpose(x.a))

def zeros(n,m=0):
    return array(linalg.zeros(n,m))

def ones(n,m=0):
    return array(linalg.ones(n,m))

def eye(n):
    return array(linalg.eye(n))

def det(x):
    if type(x)==array:
        return linalg.det(x.a)
    return linalg.det(x)

def inv(x):
    if type(x)==array:
        return array(linalg.inv(x.a))
    return linalg.inv(x)

def solve(a,b):
    if type(a)==array:
        if type(b)==array:
            return array(linalg.solve(a.a,b.a))
        return array(linalg.solve(a.a,b))
    if type(b)==array:
        return array(linalg.solve(a,b.a))
    return linalg.solve(a,b)

def eig(a):
    if type(a)==array:
        r=linalg.eig(a.a)
        return array(r[0]),array(r[1])
    return linalg.eig(a)

def linspace(a,b,c):
    return array(linalg.linspace(a,b,c))

def arange(a,b,c=1):
    return array(linalg.arange(a,b,c))

def reshape(a,n,m=0):
    if type(n)==tuple:
        m=n[1]
        n=n[0]
    if type(a)==array:
        return array(linalg.matrix(n,m,a.a))
    return linalg.matrix(n,m,a)
)";

mp_lexer_t * mp_lexer_new_from_file(const char * filename) {
  const char * script=read_file(filename);
  if (!script && strcmp(filename,"numpy.py")==0)
    script=numpy_script;
  if (script) 
    return mp_lexer_new_from_str_len(qstr_from_str(filename), script, strlen(script), 0 /* size_t free_len*/);
  else
    mp_raise_OSError(MP_ENOENT);
}
#endif

mp_import_stat_t mp_import_stat(const char *path) {
  if (strcmp(path,"numpy.py")==0 || file_exists(path)) {
    return MP_IMPORT_STAT_FILE;
  }
  return MP_IMPORT_STAT_NO_EXIST;
}


#if MICROPY_KBD_EXCEPTION
int mp_interrupt_char=5;

void mp_hal_set_interrupt_char(int c) {
  if (c != -1) {
    mp_obj_exception_clear_traceback(MP_OBJ_FROM_PTR(&MP_STATE_VM(mp_kbd_exception)));
  }
  mp_interrupt_char = c;
}

void mp_keyboard_interrupt(void) {
  MP_STATE_VM(mp_pending_exception) = MP_OBJ_FROM_PTR(&MP_STATE_VM(mp_kbd_exception));
}

bool back_key_pressed();
int getkey(int allow_suspend);




#endif

// Command line options, with their defaults
// STATIC uint emit_opt = MP_EMIT_OPT_NONE;

#if MICROPY_ENABLE_GC
// Heap size of GC heap (if enabled)
// Make it larger on a 64 bit machine, because pointers are larger.
long heap_size = 1024*1024 * (sizeof(mp_uint_t) / 4);
#endif

STATIC void stderr_print_strn(void *env, const char *str, size_t len) {
  console_output(str,len); return;
    (void)env;
    ssize_t dummy = write(STDERR_FILENO, str, len);
    mp_uos_dupterm_tx_strn(str, len);
    (void)dummy;
}

const mp_print_t mp_stderr_print = {NULL, stderr_print_strn};

#define FORCED_EXIT (0x100)
// If exc is SystemExit, return value where FORCED_EXIT bit set,
// and lower 8 bits are SystemExit value. For all other exceptions,
// return 1.
STATIC int handle_uncaught_exception(mp_obj_base_t *exc) {
    // check for SystemExit
    if (mp_obj_is_subclass_fast(MP_OBJ_FROM_PTR(exc->type), MP_OBJ_FROM_PTR(&mp_type_SystemExit))) {
        // None is an exit value of 0; an int is its value; anything else is 1
        mp_obj_t exit_val = mp_obj_exception_get_value(MP_OBJ_FROM_PTR(exc));
        mp_int_t val = 0;
        if (exit_val != mp_const_none && !mp_obj_get_int_maybe(exit_val, &val)) {
            val = 1;
        }
        return FORCED_EXIT | (val & 255);
    }

    // Report all other exceptions
    mp_obj_print_exception(&mp_stderr_print, MP_OBJ_FROM_PTR(exc));
    return 1;
}

#define LEX_SRC_STR (1)
#define LEX_SRC_VSTR (2)
#define LEX_SRC_FILENAME (3)
#define LEX_SRC_STDIN (4)

// Returns standard error codes: 0 for success, 1 for all other errors,
// except if FORCED_EXIT bit is set then script raised SystemExit and the
// value of the exit is in the lower 8 bits of the return value
int execute_from_lexer(int source_kind, const void *source, mp_parse_input_kind_t input_kind, bool is_repl) {

    nlr_buf_t nlr;
    if (nlr_push(&nlr) == 0) {
        // create lexer based on source kind
        mp_lexer_t *lex;
        if (source_kind == LEX_SRC_STR) {
            const char *line = source;
            lex = mp_lexer_new_from_str_len(MP_QSTR__lt_stdin_gt_, line, strlen(line), false);
        } else if (source_kind == LEX_SRC_VSTR) {
            const vstr_t *vstr = source;
            lex = mp_lexer_new_from_str_len(MP_QSTR__lt_stdin_gt_, vstr->buf, vstr->len, false);
        } else if (source_kind == LEX_SRC_FILENAME) {
            lex = mp_lexer_new_from_file((const char*)source);
        } else { // LEX_SRC_STDIN
            lex = mp_lexer_new_from_fd(MP_QSTR__lt_stdin_gt_, 0, false);
        }

        qstr source_name = lex->source_name;

        #if MICROPY_PY___FILE__
        if (input_kind == MP_PARSE_FILE_INPUT) {
            mp_store_global(MP_QSTR___file__, MP_OBJ_NEW_QSTR(source_name));
        }
        #endif
	parser_errorline=0; parser_errorcol=0;
        mp_parse_tree_t parse_tree = mp_parse(lex, input_kind);

        mp_obj_t module_fun = mp_compile(&parse_tree, source_name, is_repl);

        // execute it
        mp_call_function_0(module_fun);
        // check for pending exception
        if (MP_STATE_VM(mp_pending_exception) != MP_OBJ_NULL) {
            mp_obj_t obj = MP_STATE_VM(mp_pending_exception);
            MP_STATE_VM(mp_pending_exception) = MP_OBJ_NULL;
            nlr_raise(obj);
        }

        nlr_pop();
        return 0;
    } else {
        // uncaught exception
        return handle_uncaught_exception(nlr.ret_val);
    }
}

STATIC char *strjoin(const char *s1, int sep_char, const char *s2) {
    int l1 = strlen(s1);
    int l2 = strlen(s2);
    char *s = malloc(l1 + l2 + 2);
    memcpy(s, s1, l1);
    if (sep_char != 0) {
        s[l1] = sep_char;
        l1 += 1;
    }
    memcpy(s + l1, s2, l2);
    s[l1 + l2] = 0;
    return s;
}

int do_repl(void) {
    mp_hal_stdout_tx_str("MicroPython " MICROPY_GIT_TAG " on " MICROPY_BUILD_DATE "\n");

    // use simple readline

    for (;;) {
        char *line = prompt(">>> ");
        if (line == NULL) {
            // EOF
            return 0;
        }
        while (mp_repl_continue_with_input(line)) {
            char *line2 = prompt("... ");
            if (line2 == NULL) {
                break;
            }
            char *line3 = strjoin(line, '\n', line2);
            free(line);
            free(line2);
            line = line3;
        }

        int ret = execute_from_lexer(LEX_SRC_STR, line, MP_PARSE_SINGLE_INPUT, true);
        if (ret & FORCED_EXIT) {
            return ret;
        }

        free(line);
    }
}

#include "main.h"
int do_file(const char *file) {
    return execute_from_lexer(LEX_SRC_FILENAME, file, MP_PARSE_FILE_INPUT, false);
}

/*
STATIC int do_str(const char *str) {
    return execute_from_lexer(LEX_SRC_STR, str, MP_PARSE_FILE_INPUT, false);
}
*/

STATIC void set_sys_argv(char *argv[], int argc, int start_arg) {
    for (int i = start_arg; i < argc; i++) {
        mp_obj_list_append(mp_sys_argv, MP_OBJ_NEW_QSTR(qstr_from_str(argv[i])));
    }
}

MP_NOINLINE int main_(int argc, char **argv);


int micropy_eval(const char * str){
  nlr_buf_t nlr;
  if (nlr_push(&nlr) == 0) {
    parser_errorline=0;
    mp_lexer_t *lex = mp_lexer_new_from_str_len(0, str, strlen(str), false);
    mp_parse_tree_t pt = mp_parse(lex, MP_PARSE_FILE_INPUT);
    mp_obj_t module_fun = mp_compile(&pt, lex->source_name,true);
    // mp_hal_set_interrupt_char((int)Ion::Keyboard::Key::Back);
    mp_call_function_0(module_fun);
    // mp_hal_set_interrupt_char(-1); // Disable interrupt
    nlr_pop();
  } else { // Uncaught exception
#if 1
    handle_uncaught_exception(nlr.ret_val);
#else
    /* mp_obj_print_exception is supposed to handle error printing. However,
     * because we want to print custom information, we copied and modified the
     * content of mp_obj_print_exception instead of calling it. */
    if (mp_obj_is_exception_instance((mp_obj_t)nlr.ret_val)) {
        size_t n, *values;
        mp_obj_exception_get_traceback((mp_obj_t)nlr.ret_val, &n, &values);
        if (n > 0) {
            assert(n % 3 == 0);
            for (int i = n - 3; i >= 0; i -= 3) {
              if (values[i] != 0 || i == 0) {
                if (values[i] == 0) {
                  mp_printf(&mp_plat_print, "  Last command\n");
                } else {
#if MICROPY_ENABLE_SOURCE_LINE
                  mp_printf(&mp_plat_print, "  File \"%q\", line %d", values[i], (int)values[i + 1]);
#else
                  mp_printf(&mp_plat_print, "  File \"%q\"", values[i]);
#endif
                  // the block name can be NULL if it's unknown
                  qstr block = values[i + 2];
                  if (block == MP_QSTR_NULL) {
                    mp_print_str(&mp_plat_print, "\n");
                  } else {
                    mp_printf(&mp_plat_print, ", in %q\n", block);
                  }
                }
              }
            }
        }
    }
    mp_obj_print_helper(&mp_plat_print, (mp_obj_t)nlr.ret_val, PRINT_EXC);
    mp_print_str(&mp_plat_print, "\n");
    /* End of mp_obj_print_exception. */
#endif
  }
  return 0;
}

#ifndef MICROPY_LIB
int main(int argc, char **argv) {
    #if MICROPY_PY_THREAD
    mp_thread_init();
    #endif
    // We should capture stack top ASAP after start, and it should be
    // captured guaranteedly before any other stack variables are allocated.
    // For this, actual main (renamed main_) should not be inlined into
    // this function. main_() itself may have other functions inlined (with
    // their own stack variables), that's why we need this main/main_ split.
    mp_stack_ctrl_init();
    return main_(argc, argv);
}
#endif

char * micropy_init(int stack_size,int heap_size){
  mp_stack_ctrl_init();
  mp_stack_set_limit(stack_size);// (32768);

  if (heap_size){
#if MICROPY_ENABLE_GC
    char *heap = malloc(heap_size);
    if (!heap)
      return 0;
    gc_init(heap, heap + heap_size - 1);
#endif
    mp_init();
    return (char *) 1;
  }
  return 0;
}

MP_NOINLINE int main_(int argc, char **argv) {
  
  char * heap=micropy_init(16*1024,32*1024);
#if MICROPY_ENABLE_GC
  if (!heap)
    return 1;
#endif
    const int NOTHING_EXECUTED = -2;
    int ret = do_repl();

    mp_deinit();

    if (heap) free(heap);

    return ret & 0xff;
}

void nlr_jump_fail(void *val) {
    printf("FATAL: uncaught NLR %p\n", val);
    exit(1);
}
