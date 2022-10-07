#ifdef QUICKJS
#include "quickjs.h"

#define countof(x) (sizeof(x) / sizeof((x)[0]))

#if 1
const char * caseval(const char *);
#else
const char * caseval(const char * s){
  return "CAS not available";
}
#endif
const char * console_prompt(const char * s);

static JSValue js_caseval(JSContext *ctx, JSValueConst this_val,int argc, JSValueConst *argv){
  if (argc==0)
    return JS_UNDEFINED;
  const char * input=JS_ToCString(ctx,argv[0]);
  const char * output=caseval(input);
  JS_FreeCString(ctx,input);
  return JS_NewString(ctx,output);
}

static JSValue js_prompt(JSContext *ctx, JSValueConst this_val,int argc, JSValueConst *argv){
  const char * input=argc?JS_ToCString(ctx,argv[0]):0;
  const char * output=console_prompt(input?input:"?");
  if (input) JS_FreeCString(ctx,input);
  return JS_NewString(ctx,output);
}

void js_add_cas( JSContext *ctx){
  JSValue global_obj, cas;
  global_obj = JS_GetGlobalObject(ctx);
  cas = JS_NewObject(ctx);
  JS_SetPropertyStr(ctx, cas, "prompt",
		    JS_NewCFunction(ctx, js_prompt, "prompt", 1));
  JS_SetPropertyStr(ctx, cas, "caseval",
		    JS_NewCFunction(ctx, js_caseval, "caseval", 1));
  JS_SetPropertyStr(ctx, global_obj, "cas", cas);
  JS_FreeValue(ctx, global_obj);
}
#endif
