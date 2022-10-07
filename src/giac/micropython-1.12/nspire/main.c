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
#include <libndls.h>

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
#include "extmod/misc.h"
#include "genhdr/mpversion.h"
#include "input.h"

// Command line options, with their defaults
// STATIC uint emit_opt = MP_EMIT_OPT_NONE;

#if MICROPY_ENABLE_GC
// Heap size of GC heap (if enabled)
// Make it larger on a 64 bit machine, because pointers are larger.
long heap_size = 1024*1024 * (sizeof(mp_uint_t) / 4);
#endif

int micropython_port_vm_hook_loop() {
  static int c = 0;

  c++;
  if (c & 0x1f) {
    return 0;
  }
  if (on_key_pressed()){
    nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Keyboard interrupt"));
    return 1;
  }
  return 0;
}

void console_output(const char *,int);
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
STATIC int execute_from_lexer(int source_kind, const void *source, mp_parse_input_kind_t input_kind, bool is_repl) {

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

STATIC int do_repl(void) {
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


int micropy_eval(const char * line){
  return execute_from_lexer(LEX_SRC_STR, line, MP_PARSE_FILE_INPUT, true);
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
  // mp_stack_set_limit(40000 * (BYTES_PER_WORD / 4));
  mp_stack_set_limit(stack_size);
  
  //Disable output buffering, otherwise interactive mode becomes useless
  setbuf(stdout, NULL);
  
  mp_stack_set_limit(stack_size);
  
#if MICROPY_ENABLE_GC
    char *heap = malloc(heap_size);
    if(!heap)
    {
	_show_msgbox("Micropython", "Heap allocation failed. Please reboot.", 0);
	return 0;
    }
    gc_init(heap, heap + heap_size - 1);
#endif

    mp_init();

    return heap;
}

MP_NOINLINE int main_(int argc, char **argv) {
    cfg_register_fileext("py", "micropython");
  
    char * heap=micropy_init(32768,heap_size);
#if MICROPY_ENABLE_GC
  if (!heap)
    return 1;
#endif
    uint path_num = 2;
    mp_obj_list_init(MP_OBJ_TO_PTR(mp_sys_path), path_num);
    mp_obj_t *path_items;
    mp_obj_list_get(mp_sys_path, &path_num, &path_items);
    mp_obj_list_init(MP_OBJ_TO_PTR(mp_sys_argv), 0);

    path_items[0] = MP_OBJ_NEW_QSTR(MP_QSTR_);
    path_items[1] = MP_OBJ_NEW_QSTR(qstr_from_str("/documents/ndless"));

    mp_obj_list_init(mp_sys_argv, 0);
  
    const int NOTHING_EXECUTED = -2;
    int ret = NOTHING_EXECUTED;

    for (int a = 1; a < argc; a++) {
        char *pathbuf = malloc(PATH_MAX);
        char *basedir = realpath(argv[a], pathbuf);
        if (basedir == NULL) {
            mp_printf(&mp_stderr_print, "%s: can't open file '%s': [Errno %d] %s\n", argv[0], argv[a], errno, strerror(errno));
            // CPython exits with 2 in such case
            ret = 2;
            break;
        }

        // Set base dir of the script as first entry in sys.path
        char *p = strrchr(basedir, '/');
        path_items[0] = MP_OBJ_NEW_QSTR(qstr_from_strn(basedir, p - basedir));
        free(pathbuf);

        set_sys_argv(argv, argc, a);
        ret = do_file(argv[a]);
        break;
    }

    if (ret == NOTHING_EXECUTED) {
        ret = do_repl();
    }
    else
    {
        puts("Press any key to exit.");
        wait_key_pressed();
    }

    mp_deinit();

    if (heap) free(heap);

    return ret & 0xff;
}

mp_import_stat_t mp_import_stat(const char *path) {
  if (strcmp(path,"numpy.py.tns")==0)
    return MP_IMPORT_STAT_FILE;
  struct stat st;
  if (stat(path, &st) == 0) {
    if (S_ISDIR(st.st_mode)) {
      return MP_IMPORT_STAT_DIR;
    } else if (S_ISREG(st.st_mode)) {
      return MP_IMPORT_STAT_FILE;
    }
  }
  return MP_IMPORT_STAT_NO_EXIST;
}

void nlr_jump_fail(void *val) {
    printf("FATAL: uncaught NLR %p\n", val);
    exit(1);
}
