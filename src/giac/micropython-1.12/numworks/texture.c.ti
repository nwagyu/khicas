#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "mpconfig.h"
#include "nlr.h"
#include "misc.h"
#include "gc.h"
#include "qstr.h"
#include "obj.h"
#include "objstr.h"
#include "runtime.h"
#include "texture.h"

#include <libndls.h>
#include <nucleus.h>

/* 
 * Small example:
 * 
 * from nsp import Texture
 * t = Texture(320, 240, None)
 * t.fill(0x0000)
 * t.display()
 * t.delete()
 * 
 * Textures are quite big, so call delete() as early as possible to free heap space.
 * 
 * Available functions:
 * fill(color): Fills the entire texture with color.
 * display(): The texture is shown on screen. It must have been created with "nsp.Texture(320, 240, None)"!
 * getPx(x, y): Returns the color of the pixel at (x/y). Throws exception if out of bounds.
 * setPx(x, y, color): Sets color of the pixel at (x/y) to color. Throws exception if out of bounds.
 * setData(str): Writes the data of the base64 string str to the texture, has to be correct size.
 * drawOnto(dest, src_x = 0, src_y = 0, src_w = self.width, src_h = self.height, dest_x = 0, dest_y = 0, dest_w = src_w, dest_h = src_h): Draws part of the texture onto dest.
 * delete(): Frees the allocated memory. Should be done manually.
 */

static mp_obj_t nsp_texture_make_new(const mp_obj_type_t *type_in, size_t n_args, size_t n_kw, const mp_obj_t *args)
{
	mp_arg_check_num(n_args, n_kw, 3, MP_OBJ_FUN_ARGS_MAX, true);

	nsp_texture_obj_t *self = m_new_obj(nsp_texture_obj_t);
	self->base.type = &nsp_texture_type;

	self->width = mp_obj_get_int(args[0]);
	self->height = mp_obj_get_int(args[1]);
	if(args[2] == mp_const_none)
		self->has_transparency = false;
	else
		self->transparent_color = mp_obj_get_int(args[2]);

	self->bitmap = gc_alloc(self->width * self->height * 2, false);

	if (!self->bitmap)
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Allocation of texture buffer failed!"));

	return self;
}

static void nsp_texture_print(const mp_print_t *print, mp_obj_t self_in, mp_print_kind_t kind)
{
	if(mp_obj_get_type(self_in) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return;
	}
		
	nsp_texture_obj_t *self = self_in;

	mp_printf(print, "Texture (w=%u, h=%u, transparent=", self->width, self->height);
	if(!self->has_transparency)
		mp_printf(print, "false, ptr=%p)", self->bitmap);
	else
		mp_printf(print, "%u, ptr=%p)", self->transparent_color, self->bitmap);
}

static mp_obj_t nsp_texture_display(mp_obj_t self_in)
{
	if(mp_obj_get_type(self_in) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return mp_const_none;
	}

	nsp_texture_obj_t *self = self_in;
	
	if(self->width != 320 || self->height != 240 || self->has_transparency)
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "The texture must have the dimensions 320x240 without transparency!"));
	else if(has_colors)
		lcd_blit(self->bitmap, SCR_320x240_565);
	else
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Not implemented in this version"));
	
	return mp_const_none;
}
static MP_DEFINE_CONST_FUN_OBJ_1(nsp_texture_display_obj, nsp_texture_display);

static mp_obj_t nsp_texture_fill(mp_obj_t self_in, mp_obj_t arg)
{
	if(mp_obj_get_type(self_in) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return mp_const_none;
	}

	nsp_texture_obj_t *self = self_in;
	uint16_t color = mp_obj_get_int(arg);
	uint16_t *start = self->bitmap, *end = self->bitmap + (self->width * self->height);
	while(start < end)
		*start++ = color;
	
	return mp_const_none;
}
static MP_DEFINE_CONST_FUN_OBJ_2(nsp_texture_fill_obj, nsp_texture_fill);

static mp_obj_t nsp_texture_setPx(uint n_args, const mp_obj_t *args) 
{
	if(mp_obj_get_type(args[0]) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return mp_const_none;
	}
	
	nsp_texture_obj_t *self = args[0];
	uint16_t x = mp_obj_get_int(args[1]), y = mp_obj_get_int(args[2]), color = mp_obj_get_int(args[3]);
	if(x < self->width && y < self->height)
		self->bitmap[x + y * self->width] = color;
	else
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Texture coordinates out of range!"));

	return mp_const_none;
}
MP_DEFINE_CONST_FUN_OBJ_VAR_BETWEEN(nsp_texture_setPx_obj, 4, 4, nsp_texture_setPx);

static mp_obj_t nsp_texture_getPx(mp_obj_t self_in, mp_obj_t x_in, mp_obj_t y_in)
{
	if(mp_obj_get_type(self_in) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return mp_const_none;
	}
	
	nsp_texture_obj_t *self = self_in;
	uint16_t x = mp_obj_get_int(x_in), y = mp_obj_get_int(y_in);
	if(x < self->width && y < self->height)
		return MP_OBJ_NEW_SMALL_INT(self->bitmap[x + y * self->width]);
	else
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Texture coordinates out of range!"));
	
	return mp_const_none;
}
static MP_DEFINE_CONST_FUN_OBJ_3(nsp_texture_getPx_obj, nsp_texture_getPx);

static const mp_arg_t nsp_drawOnto_args[] = {
	{ MP_QSTR_src, MP_ARG_REQUIRED | MP_ARG_OBJ, {} },
	{ MP_QSTR_dest,  MP_ARG_REQUIRED | MP_ARG_OBJ, {} },
	{ MP_QSTR_src_x, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_src_y, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_src_w, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_src_h, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_dest_x, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_dest_y, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_dest_w, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
	{ MP_QSTR_dest_h, MP_ARG_KW_ONLY | MP_ARG_INT, { .u_int = -1 } },
};

static mp_obj_t nsp_texture_drawOnto(uint n_args, const mp_obj_t *args, mp_map_t *kw_args)
{
	nsp_texture_obj_t *src, *dest;
	uint16_t src_x = 0, src_y = 0, src_w, src_h, dest_x = 0, dest_y = 0, dest_w, dest_h;
	
	if(mp_obj_get_type(args[0]) != &nsp_texture_type || mp_obj_get_type(args[1]) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return mp_const_none;
	}
	
	mp_arg_val_t args_val[sizeof(nsp_drawOnto_args)/sizeof(*nsp_drawOnto_args)];
	mp_arg_parse_all(n_args, args, kw_args, sizeof(nsp_drawOnto_args)/sizeof(*nsp_drawOnto_args), nsp_drawOnto_args, args_val);
	
	src = args_val[0].u_obj;
	dest = args_val[1].u_obj;
	
	if(args_val[2].u_int != -1)
		src_x = args_val[2].u_int;
	
	if(args_val[3].u_int != -1)
		src_y = args_val[3].u_int;
	
	if(args_val[4].u_int == -1)
		src_w = src->width;
	else
		src_w = args_val[4].u_int;
	
	if(args_val[5].u_int == -1)
		src_h = src->height;
	else
		src_h = args_val[5].u_int;
	
	if(args_val[6].u_int != -1)
		dest_x = args_val[6].u_int;
	
	if(args_val[7].u_int != -1)
		dest_y = args_val[7].u_int;
	
	if(args_val[8].u_int == -1)
		dest_w = src->width;
	else
		dest_w = args_val[8].u_int;
	
	if(args_val[9].u_int == -1)
		dest_h = src->height;
	else
		dest_h = args_val[9].u_int;
	
	if(src_x + src_w > src->width || src_y + src_h > src->height || dest_x + dest_w > dest->width || dest_y + dest_h > dest->height)
		return mp_const_none;
	
	uint16_t *dest_ptr = dest->bitmap + dest_x + dest_y * dest->width;
	const unsigned int dest_nextline = dest->width - dest_w;
	
	//Special cases, for better performance
	if(src_w == dest_w && src_h == dest_h)
	{
		dest_w = MIN(src_w, dest->width - dest_x);
		dest_h = MIN(src_h, dest->height - dest_y);

		const uint16_t *src_ptr = src->bitmap + src_x + src_y * src->width;
		const unsigned int src_nextline = src->width - dest_w;
		
		if(!src->has_transparency)
		{
			for(unsigned int i = dest_h; i--;)
			{
				for(unsigned int j = dest_w; j--;)
					*dest_ptr++ = *src_ptr++;
				
				dest_ptr += dest_nextline;
				src_ptr += src_nextline;
			}
		}
		else
		{
			for(unsigned int i = dest_h; i--;)
			{
				for(unsigned int j = dest_w; j--;)
				{
					uint16_t c = *src_ptr++;
					if(c != src->transparent_color)
						*dest_ptr = c;
					
					++dest_ptr;
				}
				
				dest_ptr += dest_nextline;
				src_ptr += src_nextline;
			}
		}
		
		return mp_const_none;
	}
	
	const float dx_src = ((float)src_w) / dest_w, dy_src = ((float)src_h) / dest_h;
	uint16_t *dest_line_end = MIN(dest_ptr + dest_w, dest->bitmap + dest->width * dest->height);
	const uint16_t *dest_ptr_end = MIN(dest_ptr + dest->width * dest_h, dest->bitmap + dest->height * dest->width);
	
	float src_fx = src_x, src_fy = src_y;
	
	if(!src->has_transparency)
	{
		while(dest_ptr < dest_ptr_end)
		{
			src_fx = src_x;

			while(dest_ptr < dest_line_end)
			{	
				*dest_ptr++ = src->bitmap[((unsigned int)src_fy * src->width) + (unsigned int)src_fx];
				
				src_fx += dx_src;
			}

			dest_ptr += dest_nextline;
			dest_line_end += dest->width;
			src_fy += dy_src;
		}
	}
	else
	{
		while(dest_ptr < dest_ptr_end)
		{
			src_fx = src_x;
			
			while(dest_ptr < dest_line_end)
			{
				uint16_t c = src->bitmap[((unsigned int)src_fy * src->width) + (unsigned int)src_fx];
				if(c != src->transparent_color)
					*dest_ptr = c;
				
				src_fx += dx_src;
				++dest_ptr;
			}
			
			dest_ptr += dest_nextline;
			dest_line_end += dest->width;
			src_fy += dy_src;
		}
	}
	
	return mp_const_none;
}
static MP_DEFINE_CONST_FUN_OBJ_KW(nsp_texture_drawOnto_obj, 1, nsp_texture_drawOnto);

/* Base64 decoder from wikipedia */

#define WHITESPACE 64
#define EQUALS     65
#define INVALID    66 
static const unsigned char d[] = {
    66,66,66,66,66,66,66,66,66,64,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,
    66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,62,66,66,66,63,52,53,
    54,55,56,57,58,59,60,61,66,66,66,65,66,66,66, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,66,66,66,66,66,66,26,27,28,
    29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,66,66,
    66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,
    66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,
    66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,
    66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,
    66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,66,
    66,66,66,66,66,66
};
 
int base64decode(const char *in, unsigned char *out, size_t outLen) { 
    const char *end = in + strlen(in);
    size_t buf = 1, len = 0;
 
    while (in < end) {
        unsigned char c = d[(int) *in++];
 
        switch (c) {
        case WHITESPACE: continue;   /* skip whitespace */
        case INVALID:    return 1;   /* invalid input, return error */
        case EQUALS:                 /* pad character, end of data */
            in = end;
            continue;
        default:
            buf = buf << 6 | c;
 
            /* If the buffer is full, split it into bytes */
            if (buf & 0x1000000) {
                if ((len += 3) > outLen) return 1; /* buffer overflow */
                *out++ = buf >> 16;
                *out++ = buf >> 8;
                *out++ = buf;
                buf = 1;
            }   
        }
    }
 
    if (buf & 0x40000) {
        if ((len += 2) > outLen) return 1; /* buffer overflow */
        *out++ = buf >> 10;
        *out++ = buf >> 2;
    }
    else if (buf & 0x1000) {
        if (++len > outLen) return 1; /* buffer overflow */
        *out++ = buf >> 4;
    }
    return 0;
}

static mp_obj_t nsp_texture_setData(mp_obj_t self_in, mp_obj_t str)
{
        if(mp_obj_get_type(self_in) != &nsp_texture_type)
        {
                nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
                return mp_const_none;
        }

	nsp_texture_obj_t *self = self_in;

	GET_STR_DATA_LEN(str, str_data, str_len)

	base64decode((const char*) str_data, (unsigned char*)self->bitmap, self->width * self->height * 2);

	return mp_const_none;
}

static MP_DEFINE_CONST_FUN_OBJ_2(nsp_texture_setData_obj, nsp_texture_setData);

static mp_obj_t nsp_texture_delete(mp_obj_t self_in)
{
	if(mp_obj_get_type(self_in) != &nsp_texture_type)
	{
		nlr_raise(mp_obj_new_exception_msg(&mp_type_ValueError, "Wrong type of argument."));
		return mp_const_none;
	}
	
	nsp_texture_obj_t *self = self_in;
	
	if(!self->bitmap)
		return mp_const_none;

	self->width = 0;
	self->height = 0;
	gc_free(self->bitmap);
	self->bitmap = (uint16_t*) 0;
	
	return mp_const_none;
}
static MP_DEFINE_CONST_FUN_OBJ_1(nsp_texture_delete_obj, nsp_texture_delete);

static const mp_map_elem_t nsp_texture_locals_dict_table[] = {
	{ MP_OBJ_NEW_QSTR(MP_QSTR_display), (mp_obj_t) &nsp_texture_display_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_fill), (mp_obj_t) &nsp_texture_fill_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_setPx), (mp_obj_t) &nsp_texture_setPx_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_getPx), (mp_obj_t) &nsp_texture_getPx_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_drawOnto), (mp_obj_t) &nsp_texture_drawOnto_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_setData), (mp_obj_t) &nsp_texture_setData_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_delete), (mp_obj_t) &nsp_texture_delete_obj },
};

static MP_DEFINE_CONST_DICT(nsp_texture_locals_dict, nsp_texture_locals_dict_table);

const mp_obj_type_t nsp_texture_type = {
    { &mp_type_type },
    .name = MP_QSTR_Texture,
    .print = nsp_texture_print,
    .make_new = nsp_texture_make_new,
    .locals_dict = (mp_obj_t)&nsp_texture_locals_dict
};
