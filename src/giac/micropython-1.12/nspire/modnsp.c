#include <libndls.h>

#include "mpconfig.h"
#include "misc.h"
#include "nlr.h"
#include "qstr.h"
#include "obj.h"
#include "runtime.h"
#include "objtuple.h"
#include "texture.h"

static mp_obj_t nsp_readRTC()
{
	return mp_obj_new_int(*(unsigned int*)0x90090000);
}
static MP_DEFINE_CONST_FUN_OBJ_0(nsp_readRTC_obj, nsp_readRTC);

static mp_obj_t nsp_waitKeypress()
{
	wait_key_pressed();
	return mp_const_none;
}
static MP_DEFINE_CONST_FUN_OBJ_0(nsp_waitKeypress_obj, nsp_waitKeypress);

STATIC const mp_map_elem_t mp_module_nsp_globals_table[] = {
        { MP_ROM_QSTR(MP_QSTR___name__), MP_ROM_QSTR(MP_QSTR_nsp) },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_Texture), (mp_obj_t) &nsp_texture_type },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_waitKeypress), (mp_obj_t) &nsp_waitKeypress_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_readRTC), (mp_obj_t) &nsp_readRTC_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_monotonic), (mp_obj_t) &nsp_readRTC_obj },
	{ MP_OBJ_NEW_QSTR(MP_QSTR_time), (mp_obj_t) &nsp_readRTC_obj },
};

STATIC const mp_obj_dict_t mp_module_nsp_globals = {
    .base = {&mp_type_dict},
    .map = {
        .all_keys_are_qstrs = 1,
        .is_fixed = 1,
        .used = MP_ARRAY_SIZE(mp_module_nsp_globals_table),
        .alloc = MP_ARRAY_SIZE(mp_module_nsp_globals_table),
        .table = (mp_map_elem_t*) mp_module_nsp_globals_table,
    },
};

const mp_obj_module_t mp_module_nsp = {
    .base = { &mp_type_module },
    .globals = (mp_obj_dict_t*) &mp_module_nsp_globals,
};
