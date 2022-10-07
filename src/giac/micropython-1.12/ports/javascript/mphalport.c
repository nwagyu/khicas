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

#include "library.h"
#include "mphalport.h"
#if defined EMCC && !defined NO_QSTR
#include <emscripten.h>
#endif

void console_output(const char *s,int i){
  EM_ASM_ARGS({
      var text = UTF8ToString($0); // Convert message to JS string
      text=text.substr(0,$1);
      UI.add_python_output(text);
      text=UI.clean_for_html(text);
      var tmp=document.getElementById('consolediv');
      if (tmp!=null && tmp.style.display != 'block') {
	tmp.style.display = 'block';
	UI.set_config_width();
      }
      var element=document.getElementById('output');
      if (element!=null){
	element.style.display = 'inherit';
	element.innerHTML += text; // element.value += text + "\n";
	element.scrollTop = 99999; // focus on bottom      
	// Module.print(text);
      }
    },s,i);
}

void mp_hal_stdout_tx_strn(const char *str, size_t len) {
  console_output(str,len);
  //mp_js_write(str, len);
}

void mp_hal_delay_ms(mp_uint_t ms) {
    uint32_t start = mp_hal_ticks_ms();
    while (mp_hal_ticks_ms() - start < ms) {
    }
}

void mp_hal_delay_us(mp_uint_t us) {
    uint32_t start = mp_hal_ticks_us();
    while (mp_hal_ticks_us() - start < us) {
    }
}

mp_uint_t mp_hal_ticks_us(void) {
    return mp_js_ticks_ms() * 1000;
}

mp_uint_t mp_hal_ticks_ms(void) {
    return mp_js_ticks_ms();
}

mp_uint_t mp_hal_ticks_cpu(void) {
    return 0;
}

extern int mp_interrupt_char;

int mp_hal_get_interrupt_char(void) {
    return mp_interrupt_char;
}
