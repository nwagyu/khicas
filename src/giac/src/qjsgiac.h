#ifdef QUICKJS
#ifndef QJSGIAC_H
#define QJSGIAC_H

#ifdef __cplusplus
extern "C" {
#endif
  extern int quickjs_parse_error_col,quickjs_parse_error_line;  
  extern void (*qjs_interrupt_callback)(void);
  struct JSRuntime;
  struct JSContext;
  typedef struct  {
    JSRuntime * rt;
    JSContext * ctx;
  } JSRuntimeContext;
  char * js_ck_eval(const char * line,JSRuntimeContext * js_contextptr);
  extern JSRuntimeContext global_js_context;
  char * quickjs_ck_eval(const char * line);
  
  char * js_eval(const char * line,JSRuntimeContext * js_contextptr);
  JSRuntimeContext js_init(size_t memory_limit,size_t stack_size,int argc,char ** argv);
  int js_loop(JSRuntimeContext rtctx);
  void js_end(JSRuntimeContext rtctx);
  
#ifdef __cplusplus
}
#endif

#endif // QJSGIAC_H
#endif // QUICKJS
