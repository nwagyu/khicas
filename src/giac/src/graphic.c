#ifdef QUICKJS
#if 1 //def EMCC
void sync_screen(); //in mp.h
#endif

#ifdef KHICAS
#include <libndls.h>
#include "os.h" // Ndless/ndless-sdk/include/os.h
#include <unistd.h>
#include <stdio.h>
#include <dirent.h>
#include <ngc.h>
#include "k_defs.h"
#include "k_csdk.h"
#else
void set_pixel(int x,int y,int c);
void os_set_pixel(int x,int y,int c){
  return set_pixel(x,y,c);
}
void os_fill_rect(int x,int y,int w,int h,int c){
  int i=0,j=0;
  for (i=0;i<w;++i){
    for (j=0;j<h;++j){
      os_set_pixel(x+i,y+j,c);
    }
  }
}
#endif
void console_freeze();
#include "quickjs.h"

#define countof(x) (sizeof(x) / sizeof((x)[0]))


static JSValue js_set_pixel(JSContext *ctx, JSValueConst this_val,
                      int argc, JSValueConst *argv)
{
  if (argc==0){
    sync_screen();
    return JS_UNDEFINED;
  }
  if (argc<2 || argc>3){
    sync_screen();
    return JS_EXCEPTION;
  }
  int32_t i,j;
  if (JS_ToInt32(ctx, &i, argv[0]) || JS_ToInt32(ctx, &j, argv[1]))
    return JS_EXCEPTION;
  int32_t c=0;
  if (argc==3 && JS_ToInt32(ctx, &c, argv[2]))
    return JS_EXCEPTION;
  console_freeze();
  os_set_pixel(i,j,c);
  return JS_UNDEFINED;
}

static JSValue js_fill_rect(JSContext *ctx, JSValueConst this_val,
                      int argc, JSValueConst *argv)
{
  if (argc<4 || argc>5){
    return JS_EXCEPTION;
  }
  int32_t x,y,w,h;
  if (JS_ToInt32(ctx, &x, argv[0]) || JS_ToInt32(ctx, &y, argv[1]) ||
      JS_ToInt32(ctx, &w, argv[2]) || JS_ToInt32(ctx, &h, argv[3])
      )
    return JS_EXCEPTION;
  int32_t c=0;
  if (argc==5 && JS_ToInt32(ctx, &c, argv[4]))
    return JS_EXCEPTION;
  console_freeze();
  os_fill_rect(x,y,w,h,c);
  return JS_UNDEFINED;
}

void js_add_graphic( JSContext *ctx){
  JSValue global_obj, graphic;
  global_obj = JS_GetGlobalObject(ctx);
  graphic = JS_NewObject(ctx);
  JS_SetPropertyStr(ctx, graphic, "setPixel",
		    JS_NewCFunction(ctx, js_set_pixel, "setPixel", 1));
  JS_SetPropertyStr(ctx, graphic, "fillRect",
		    JS_NewCFunction(ctx, js_fill_rect, "fillRect", 1));
  JS_SetPropertyStr(ctx, global_obj, "graphic", graphic);
  JS_FreeValue(ctx, global_obj);
}
#endif
