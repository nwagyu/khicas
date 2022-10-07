#include <stdint.h>

typedef struct nsp_texture_obj_t {
    mp_obj_base_t base;
    uint16_t width;
    uint16_t height;
    uint16_t *bitmap;
    bool has_transparency;
    uint16_t transparent_color;
} nsp_texture_obj_t;

extern const mp_obj_type_t nsp_texture_type;
