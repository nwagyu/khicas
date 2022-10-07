/*
 * This file is part of the Micro Python project, http://micropython.org/
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

int micropython_port_vm_hook_loop() ;
// options to control how Micro Python is built
#define MICROPY_VM_HOOK_LOOP micropython_port_vm_hook_loop();

#define MICROPY_ALLOC_PATH_MAX      (PATH_MAX)
#define MICROPY_EMIT_ARM            (1)
#define MICROPY_ENABLE_GC           (1)
#define MICROPY_ENABLE_FINALISER    (1)
#define MICROPY_MEM_STATS           (0)
#define MICROPY_DEBUG_PRINTERS      (1)
#define MICROPY_HELPER_REPL         (1)
#define MICROPY_HELPER_LEXER_UNIX   (1)
#define MICROPY_ENABLE_SOURCE_LINE  (1)
#define MICROPY_FLOAT_IMPL          (MICROPY_FLOAT_IMPL_DOUBLE)
#define MICROPY_LONGINT_IMPL        (MICROPY_LONGINT_IMPL_MPZ)
#define MICROPY_STREAMS_NON_BLOCK   (1)
#define MICROPY_OPT_COMPUTED_GOTO   (1)
#define MICROPY_PY_BUILTINS_STR_UNICODE (1)
#define MICROPY_PY_BUILTINS_FROZENSET (1)
#define MICROPY_PY_BUILTINS_COMPILE (1)
#define MICROPY_PY_SYS_EXIT         (1)
#define MICROPY_PY_SYS_PLATFORM     "nspire"
#define MICROPY_PY_SYS_MAXSIZE      (1)
#define MICROPY_PY_SYS_STDFILES     (1)
#define MICROPY_PY_CMATH            (1)
#define MICROPY_PY_IO_FILEIO        (0)
#define MICROPY_PY_GC_COLLECT_RETVAL (1)
#define MICROPY_COMP_MODULE_CONST   (1)
#define MICROPY_COMP_TRIPLE_TUPLE_ASSIGN (1)
#define MICROPY_STACK_CHECK         (1)
#define MICROPY_MALLOC_USES_ALLOCATED_SIZE (1)
#define MICROPY_HELPER_REPL         (1)
#define MICROPY_READER_POSIX        (1)
#define MICROPY_REPL_AUTO_INDENT    (0)
#define MICROPY_OPT_CACHE_MAP_LOOKUP_IN_BYTECODE (1)
#define MICROPY_CAN_OVERRIDE_BUILTINS (1)
#define MICROPY_USE_INTERNAL_PRINTF (0)
#define MICROPY_PY_FUNCTION_ATTRS   (1)
#define MICROPY_PY_DESCRIPTORS      (1)
#define MICROPY_PY_BUILTINS_STR_SPLITLINES (1)
#define MICROPY_PY_BUILTINS_MEMORYVIEW (1)
#define MICROPY_PY_BUILTINS_NOTIMPLEMENTED (1)
#define MICROPY_PY_MICROPYTHON_MEM_INFO (1)
#define MICROPY_PY_ALL_SPECIAL_METHODS (1)
#define MICROPY_PY_ARRAY_SLICE_ASSIGN (1)
#define MICROPY_PY_BUILTINS_SLICE_ATTRS (1)
#define MICROPY_PY_SYS_EXC_INFO     (1)
#define MICROPY_PY_COLLECTIONS_ORDEREDDICT (1)
#define MICROPY_PY_MATH_SPECIAL_FUNCTIONS (1)
#define MICROPY_MODULE_FROZEN       (0)
#define MICROPY_STACKLESS           (0)
#define MICROPY_STACKLESS_STRICT    (0)
#define MICROPY_PY_UZLIB            (1)
#define MICROPY_PY_UJSON            (1)
#define MICROPY_PY_URE              (1)
#define MICROPY_PY_UHEAPQ           (1)
#define MICROPY_PY_UHASHLIB         (1)
#define MICROPY_PY_UBINASCII        (1)
#define MICROPY_PY_MACHINE          (1)
#define MICROPY_WARNINGS            (1)

#define MICROPY_PY_UCTYPES          (1)

// Define to MICROPY_ERROR_REPORTING_DETAILED to get function, etc.
// names in exception messages (may require more RAM).
#define MICROPY_ERROR_REPORTING     (MICROPY_ERROR_REPORTING_DETAILED)
// Define to 1 to use untested inefficient GC helper implementation
// (if more efficient arch-specific one is not available).
#ifndef MICROPY_GCREGS_SETJMP
#define MICROPY_GCREGS_SETJMP       (0)
#endif

#define MICROPY_ENABLE_EMERGENCY_EXCEPTION_BUF   (1)
#define MICROPY_EMERGENCY_EXCEPTION_BUF_SIZE  (128)

extern const struct _mp_obj_module_t mp_module_os;
extern const struct _mp_obj_module_t mp_module_nsp;
extern const struct _mp_obj_module_t mp_module_graphic;
extern const struct _mp_obj_module_t mp_module_cas;
extern const struct _mp_obj_module_t mp_module_turtle;
extern const struct _mp_obj_module_t mp_module_matplotl;
extern const struct _mp_obj_module_t mp_module_linalg;
extern const struct _mp_obj_module_t mp_module_arit;
extern const struct _mp_obj_module_t ulab_user_cmodule;

#define MICROPY_PORT_BUILTIN_MODULES \
	{ MP_OBJ_NEW_QSTR(MP_QSTR__os), (mp_obj_t) &mp_module_os }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_nsp), (mp_obj_t) &mp_module_nsp }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_time), (mp_obj_t) &mp_module_nsp }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_cas), (mp_obj_t) &mp_module_cas }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_xcas), (mp_obj_t) &mp_module_cas }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_arit), (mp_obj_t) &mp_module_arit }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_linalg), (mp_obj_t) &mp_module_linalg }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_matplotl), (mp_obj_t) &mp_module_matplotl }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_matplotlib), (mp_obj_t) &mp_module_matplotl }, \
        { MP_OBJ_NEW_QSTR(MP_QSTR_matplotlib_dot_pyplot), (mp_obj_t) &mp_module_matplotl }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_pylab), (mp_obj_t) &mp_module_matplotl }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_turtle), (mp_obj_t) &mp_module_turtle }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_casioplot), (mp_obj_t) &mp_module_graphic }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_kandinsky), (mp_obj_t) &mp_module_graphic }, \
	{ MP_OBJ_NEW_QSTR(MP_QSTR_graphic), (mp_obj_t) &mp_module_graphic },\
        { MP_ROM_QSTR(MP_QSTR_ulab), (mp_obj_t) &ulab_user_cmodule },

typedef int mp_int_t;
typedef unsigned int mp_uint_t;
typedef long mp_off_t;

#define BYTES_PER_WORD sizeof(mp_int_t)
#define MP_SSIZE_MAX INT_MAX

typedef void *machine_ptr_t; // must be of pointer size
typedef const void *machine_const_ptr_t; // must be of pointer size

#define MICROPY_PORT_BUILTINS \
    { MP_OBJ_NEW_QSTR(MP_QSTR_input), (mp_obj_t)&mp_builtin_input_obj }, \
    { MP_OBJ_NEW_QSTR(MP_QSTR_open), (mp_obj_t)&mp_builtin_open_obj },

#include <alloca.h>
