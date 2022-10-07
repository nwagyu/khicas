#ifdef QUICKJS
/*
 * QuickJS stand alone interpreter
 * 
 * Copyright (c) 2017-2020 Fabrice Bellard
 * Copyright (c) 2017-2020 Charlie Gordon
 * Changes made by B. Parisse for Giac/Xcas/KhiCAS interaction
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
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */ 
#include "cutils.h"
#include "quickjs-libc.h"
#include <qjsgiac.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <time.h>
#if defined(__APPLE__)
#include <malloc/malloc.h>
#elif defined(__linux__)
#include <malloc.h>
#endif

#ifdef NSPIRE_NEWLIB
#include "os.h"
#else
#include "js.h"
#endif

#ifdef EMCC
#include <emscripten.h>
#endif


extern size_t pythonjs_heap_size,pythonjs_stack_size; // in global.h
const char * caseval(const char *);
extern void * bf_ctx_ptr;

JSRuntimeContext global_js_context={0,0};
char * js_ck_eval(const char * line,JSRuntimeContext * js_contextptr){
  if (!js_contextptr->ctx){
    *js_contextptr=js_init(pythonjs_heap_size,pythonjs_stack_size,0,0);
    bf_ctx_ptr=get_bf_context(js_contextptr->rt);
    //caseval("python_compat(-1)");
  }
  return js_eval(line,js_contextptr);
}

char * quickjs_ck_eval(const char * line){
  return js_ck_eval(line,&global_js_context);
}

// extern const uint8_t qjsc_repl[];
// extern const uint32_t qjsc_repl_size;
#ifdef CONFIG_BIGNUM
extern const uint8_t qjsc_qjscalc[];
extern const uint32_t qjsc_qjscalc_size;
static int bignum_ext;
#endif

static int eval_buf(JSContext *ctx, const void *buf, int buf_len,
                    const char *filename, int eval_flags)
{
  JSValue val;
  int ret;

  if ((eval_flags & JS_EVAL_TYPE_MASK) == JS_EVAL_TYPE_MODULE) {
    /* for the modules, we compile then run to be able to set
       import.meta */
    val = JS_Eval(ctx, (const char *) buf, buf_len, filename,
		  eval_flags | JS_EVAL_FLAG_COMPILE_ONLY);
    if (!JS_IsException(val)) {
      js_module_set_import_meta(ctx, val, TRUE, TRUE);
      val = JS_EvalFunction(ctx, val);
    }
  } else {
    val = JS_Eval(ctx, (const char *) buf, buf_len, filename, eval_flags);
  }
  if (JS_IsException(val)) {
    js_std_dump_error(ctx);
    ret = -1;
  } else {
    ret = 0;
  }
  JS_FreeValue(ctx, val);
  return ret;
}

static int eval_file(JSContext *ctx, const char *filename, int module)
{
  uint8_t *buf;
  int ret, eval_flags;
  size_t buf_len;
    
  buf = js_load_file(ctx, &buf_len, filename);
  if (!buf) {
    perror(filename);
    exit(1);
  }

  if (module < 0) {
    module = (has_suffix(filename, ".mjs") ||
	      JS_DetectModule((const char *)buf, buf_len));
  }
  if (module)
    eval_flags = JS_EVAL_TYPE_MODULE;
  else
    eval_flags = JS_EVAL_TYPE_GLOBAL;
  ret = eval_buf(ctx, buf, buf_len, filename, eval_flags);
  js_free(ctx, buf);
  return ret;
}

/* also used to initialize the worker context */
static JSContext *JS_NewCustomContext(JSRuntime *rt)
{
  JSContext *ctx;
  ctx = JS_NewContext(rt);
  if (!ctx)
    return NULL;
#ifdef CONFIG_BIGNUM
  if (bignum_ext) {
    JS_AddIntrinsicBigFloat(ctx);
    JS_AddIntrinsicBigDecimal(ctx);
    JS_AddIntrinsicOperators(ctx);
    JS_EnableBignumExt(ctx, TRUE);
  }
#endif
  /* system modules */
  js_init_module_std(ctx, "std");
  js_init_module_os(ctx, "os");
  return ctx;
}

void js_add_cas( JSContext *ctx);
void js_add_graphic( JSContext *ctx); 
void (*qjs_interrupt_callback)(void)=0;
#ifdef NSPIRE_NEWLIB
void process_freeze(); void os_fill_rect(int,int,int,int,int);
int interrupt_handler(JSRuntime *rt, void *opaque){
  if (isKeyPressed(KEY_NSPIRE_ESC))
    return 1;
  return 0;
}
#else
int ctrl_c_interrupted(int exception);
int interrupt_handler(JSRuntime *rt, void *opaque){
  if (qjs_interrupt_callback)
    qjs_interrupt_callback();
  if (ctrl_c_interrupted(0))
    return 1;
  return 0;
}
void process_freeze(){
}
#endif


JSRuntimeContext js_init(size_t memory_limit,size_t stack_size,int argc,char ** argv){
  JSRuntime *rt;
  JSContext *ctx;
  int load_std = 0;
  if (!memory_limit) memory_limit = 2*1024*1024;
  if (!stack_size) stack_size = 1*1024*1024;
  if (stack_size<64*1024) stack_size=32*1024;
  //if (memory_limit<512*1024) memory_limit = 512*1024;
  //stack_size=0;
  memory_limit=0;
    
#ifdef CONFIG_BIGNUM
  int load_jscalc=1;
  if (load_jscalc)
    bignum_ext = 1;
#endif

  rt = JS_NewRuntime();

  if (!rt) {
    fprintf(stderr, "qjs: cannot allocate JS runtime\n");
    JSRuntimeContext res={0,0};
    return res;
  }
  if (memory_limit != 0)
    JS_SetMemoryLimit(rt, memory_limit); // DO IT for nspire
  if (stack_size != 0)
    JS_SetMaxStackSize(rt, stack_size); // DO IT for nspire
  js_std_set_worker_new_context_func(JS_NewCustomContext);
  js_std_init_handlers(rt);
  ctx = JS_NewCustomContext(rt);
  if (!ctx) {
    fprintf(stderr, "qjs: cannot allocate JS context\n");
    JSRuntimeContext res={0,0};
    return res;
  }

  /* loader for ES6 modules */
  JS_SetModuleLoaderFunc(rt, NULL, js_module_loader, NULL);

  int dump_unhandled_promise_rejection = 0;
  if (dump_unhandled_promise_rejection) {
    JS_SetHostPromiseRejectionTracker(rt, js_std_promise_rejection_tracker,
				      NULL);
  }
      
#ifdef CONFIG_BIGNUM
  if (load_jscalc) {
#if 1 // def NSPIRE_NEWLIB
#ifndef NSPIRE_NEWLIB
    const char * ptr=qjscalcjs;
#if 0
    EM_ASM_ARGS({
	var msg = Pointer_stringify($0); // Convert message to JS string
	console.log('qjscalcjs',msg);                      
      }, qjscalcjs);
#endif
#else
    extern char _binary_qjscalcjs_js_start,_binary_qjscalcjs_js_end;
    // JSValue val=JS_Eval(ctx,&_binary_qjscalcjs_js_start,&_binary_qjscalcjs_js_end-&_binary_qjscalcjs_js_start-2,0,0);
    const char * ptr=&_binary_qjscalcjs_js_start;
#endif // NSPIRE_NEWLIB
    size_t len=strlen(ptr);
    JSValue val=JS_Eval(ctx,ptr,len,"qjscalcjs",0);
    if (JS_IsException(val)){
      const char * s=JS_ToCString(ctx,JS_GetException(ctx));
      if (s)
	printf("Loading qjscalcjs: %s\n",s);
    }
    JS_FreeValue(ctx,val);
#else
    js_std_eval_binary(ctx, qjsc_qjscalc, qjsc_qjscalc_size, 0);
#endif // 1
  }
#endif // CONFIG_BIGNUM

  js_std_add_helpers(ctx, argc , argv );
  js_add_cas(ctx);
  js_add_graphic(ctx);
  JS_SetInterruptHandler(JS_GetRuntime(ctx), interrupt_handler, NULL);
#ifndef NSPIRE_NEWLIB
  if (1){
#if 1 // def EMCC
    const char * ptr=xcasjs;
#if 0
    EM_ASM_ARGS({
	var msg = Pointer_stringify($0); // Convert message to JS string
	console.log('xcasjs',msg);                      
      }, xcasjs);
#endif
#else
    extern char _binary_xcasjs_js_start,_binary_xcasjs_js_end;
    const char * ptr=&_binary_xcasjs_js_start;
#endif
    size_t len=strlen(ptr);
    JSValue val=JS_Eval(ctx,ptr,len,"xcasjs",0);
    if (JS_IsException(val)){
      const char * s=JS_ToCString(ctx,JS_GetException(ctx));
      if (s)
	printf("Loading xcasjs: %s\n",s);
    }
    JS_FreeValue(ctx,val);
  }
#endif

#ifdef NSPIRE_NEWLIB
  extern char _binary_tijs_js_start,_binary_tijs_js_end;
  //JSValue val=JS_Eval(ctx,&_binary_tijs_js_start,&_binary_tijs_js_end-&_binary_tijs_js_start-2,0,0);
  JSValue val=JS_Eval(ctx,&_binary_tijs_js_start,strlen(&_binary_tijs_js_start),"<cmdline>",0);
  if (JS_IsException(val)){
    const char * s=JS_ToCString(ctx,JS_GetException(ctx));
    if (s)
      printf("Loading ti: %s\n",s);
  }
  JS_FreeValue(ctx,val);
#else // FIXME
#endif
    
  /* make 'std' and 'os' visible to non module code */
  if (load_std) {
    const char *str = "import * as std from 'std';\n"
      "import * as os from 'os';\n"
      "globalThis.std = std;\n"
      "globalThis.os = os;\n";
    eval_buf(ctx, str, strlen(str), "<input>", JS_EVAL_TYPE_MODULE);
  }
  JSRuntimeContext rtctx={rt,ctx};
  return rtctx;
}

int js_loop(JSRuntimeContext rtctx){
  // LOOP
  process_freeze();
  printf("*** QuickJS interpreter, type 0 to quit ***\n");
  while (1){
    char line[1024];
    scanf("%s",line);
    if (strlen(line)==1 && line[0]=='0')
      break;
    JSValue val=JS_Eval(rtctx.ctx,line,strlen(line),"<cmdline>",0);
    process_freeze();
    if (JS_IsException(val)){
      const char * s=JS_ToCString(rtctx.ctx,JS_GetException(rtctx.ctx));
      if (s)
	printf("  %s\n",s);
    }
    else {
      const char * s=JS_ToCString(rtctx.ctx,val);
      if (s)
	printf("  %s\n",s);
      else
	printf("  Error\n");
    }
    JS_FreeValue(rtctx.ctx,val);
    js_std_loop(rtctx.ctx);
  }

  return 0;

 fail:
  return 1;
}

int quickjs_parse_error_col=-1;
char * js_eval(const char * line,JSRuntimeContext * js_contextptr){
  quickjs_parse_error_col=quickjs_parse_error_line=-1;
  JSValue val=JS_Eval(js_contextptr->ctx,line,strlen(line),"<cmdline>",0);
  const char * s=0;
  if (JS_IsException(val)){
    /* this assumes that js_parse_error in quickjs.c has been modified : 
...
int quickjs_parse_error_line=-1;
const char * quickjs_parse_ptr=0;
int __attribute__((format(printf, 2, 3))) js_parse_error(JSParseState *s, const char *fmt, ...){
  quickjs_parse_error_line=s->line_num;
  quickjs_parse_ptr=s->buf_ptr;
  ...
     */ 
    if (quickjs_parse_error_line>0){
      const char * ptr=line;int l=0;
      for (;*ptr && l<quickjs_parse_error_line;++ptr){
	if (*ptr=='\n')
	  ++l;
      }
      quickjs_parse_error_col=1+(ptr-quickjs_parse_ptr);
    }
    //fprintf(stderr,"Line %i Column %i\n",quickjs_parse_error_line,quickjs_parse_error_col);
    s=JS_ToCString(js_contextptr->ctx,JS_GetException(js_contextptr->ctx));
  }
  else {
    s=JS_ToCString(js_contextptr->ctx,val);
    JS_FreeValue(js_contextptr->ctx,val);
  }
  char * ptr=(char *) malloc(s?strlen(s)+128:16);
  strcpy(ptr,s?s:" Error");
  if (s) JS_FreeCString(js_contextptr->ctx,s);
  if (quickjs_parse_error_col>=0 && quickjs_parse_error_line>=0){
    sprintf(ptr+strlen(ptr)," line %i col %i",quickjs_parse_error_line,quickjs_parse_error_col);
  }
  return ptr;
}

void js_end(JSRuntimeContext rtctx){
  js_std_free_handlers(rtctx.rt);
  JS_FreeContext(rtctx.ctx);
  JS_FreeRuntime(rtctx.rt);
  rtctx.ctx=0; rtctx.rt=0;
}

#endif // QUICKJS
