include ../../py/mkenv.mk

CROSS = 0

QSTR_DEFS = qstrdefsport.h

include $(TOP)/py/py.mk

CC = emcc -g4
LD = emcc -g4

INC += -I.
INC += -I../../py
INC += -I$(TOP)
INC += -I$(BUILD)

CPP = clang -E

ifdef EMSCRIPTEN
    CPP += -isystem $(EMSCRIPTEN)/system/include/libc -cxx-isystem $(EMSCRIPTEN)/system/include/libcxx
else
	EMCC=/home/parisse/emsdk/fastcomp/emscripten
	CPP += -I$(EMCC)/system/include/libc -I$(EMCC)/system/include/libcxx
endif

CFLAGS = -m32 -Wall -Werror $(INC)  $(COPT)
LDFLAGS = -m32 -Wl,-Map=$@.map,--cref -Wl,--gc-sections -s ERROR_ON_UNDEFINED_SYMBOLS=0

CFLAGS += -O0 -DNDEBUG -DMICROPY_LIB -DEMCC
CFLAGS += -fdata-sections -ffunction-sections

ifneq ($(FROZEN_MPY_DIR),)
# To use frozen bytecode, put your .py files in a subdirectory (eg frozen/) and
# then invoke make with FROZEN_MPY_DIR=frozen (be sure to build from scratch).
CFLAGS += -DMICROPY_QSTR_EXTRA_POOL=mp_qstr_frozen_const_pool
CFLAGS += -DMICROPY_MODULE_FROZEN_MPY
endif

SRC_LIB = $(addprefix lib/,\
	utils/interrupt_char.c \
	utils/stdout_helpers.c \
	utils/pyexec.c \
	mp-readline/readline.c \
	)

SRC_C = \
	main.c \
	mphalport.c \
	modutime.c \
	graphic.c file.c \

USERMODULES_DIR := ulab

CFLAGS_USERMOD += -I$(USERMODULES_DIR)

override CFLAGS_EXTRA += -DMODULE_ULAB_ENABLED=1

SRC_C += $(USERMODULES_DIR)/scipy/linalg/linalg.c
SRC_C += $(USERMODULES_DIR)/scipy/optimize/optimize.c
SRC_C += $(USERMODULES_DIR)/scipy/signal/signal.c
SRC_C += $(USERMODULES_DIR)/scipy/special/special.c
SRC_C += $(USERMODULES_DIR)/ndarray_operators.c
SRC_C += $(USERMODULES_DIR)/ulab_tools.c
SRC_C += $(USERMODULES_DIR)/ndarray.c
SRC_C += $(USERMODULES_DIR)/numpy/ndarray/ndarray_iter.c
SRC_C += $(USERMODULES_DIR)/ndarray_properties.c
SRC_C += $(USERMODULES_DIR)/numpy/approx.c
SRC_C += $(USERMODULES_DIR)/numpy/compare.c
SRC_C += $(USERMODULES_DIR)/ulab_create.c
SRC_C += $(USERMODULES_DIR)/numpy/fft/fft.c
SRC_C += $(USERMODULES_DIR)/numpy/fft/fft_tools.c
SRC_C += $(USERMODULES_DIR)/numpy/filter.c
SRC_C += $(USERMODULES_DIR)/numpy/linalg/linalg.c
SRC_C += $(USERMODULES_DIR)/numpy/linalg/linalg_tools.c
SRC_C += $(USERMODULES_DIR)/numpy/numerical.c
SRC_C += $(USERMODULES_DIR)/numpy/poly.c
SRC_C += $(USERMODULES_DIR)/numpy/stats.c
SRC_C += $(USERMODULES_DIR)/numpy/transform.c
SRC_C += $(USERMODULES_DIR)/numpy/vector.c
SRC_C += $(USERMODULES_DIR)/numpy/numpy.c
SRC_C += $(USERMODULES_DIR)/scipy/scipy.c
SRC_C += $(USERMODULES_DIR)/user/user.c
SRC_C += $(USERMODULES_DIR)/utils/utils.c
SRC_C += $(USERMODULES_DIR)/ulab.c

SRC_QSTR += $(SRC_C)

OBJ = 
OBJ = $(PY_O) 
OBJ += $(addprefix $(BUILD)/, $(SRC_LIB:.c=.o))
OBJ += $(addprefix $(BUILD)/, $(SRC_C:.c=.o))

JSFLAGS = -O0 -s EXPORTED_FUNCTIONS="['_mp_js_init', '_mp_js_init_repl', '_mp_js_do_str', '_mp_js_process_char', '_mp_hal_get_interrupt_char', '_mp_keyboard_interrupt']" -s EXTRA_EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" -s "BINARYEN_TRAP_MODE='clamp'" --memory-init-file 0 --js-library library.js

all: $(BUILD)/micropython.js

$(BUILD)/micropython.js: $(OBJ) library.js wrapper.js
	$(ECHO) "LINK $(BUILD)/firmware.js"
	$(Q)emcc $(LDFLAGS) -o $(BUILD)/firmware.js $(OBJ) $(JSFLAGS)
	cat wrapper.js $(BUILD)/firmware.js > $@

distdir:
	/bin/cp *.c *.h *.mk Makefile mklib $(distdir)

min: $(BUILD)/micropython.js
	uglifyjs $< -c -o $(BUILD)/micropython.min.js

test: $(BUILD)/micropython.js $(TOP)/tests/run-tests
	$(eval DIRNAME=ports/$(notdir $(CURDIR)))
	cd $(TOP)/tests && MICROPY_MICROPYTHON=../ports/javascript/node_run.sh ./run-tests

include $(TOP)/py/mkrules.mk
