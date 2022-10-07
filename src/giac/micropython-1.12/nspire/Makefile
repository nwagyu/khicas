include ../py/mkenv.mk
-include mpconfigport.mk

CC := nspire-gcc
LD := arm-none-eabi-ld.gold

ZEHNFLAGS = --name "micropython" --author "micropython team, ported by Fabian Vogt" --notice "Interactive python console" --version 09 --compress

# define main target
PROG = micropython

# qstr definitions (must come before including py.mk)
QSTR_DEFS = qstrdefsport.h

# OS name, for simple autoconfig
UNAME_S := $(shell uname -s)

# Disable stripping of the main executable, relocation info is needed by genzehn
NOSTRIP := 1
# DEBUG := 1

# include py core make definitions
include ../py/py.mk

INC =  -I.
INC += -I..
INC += -I../py
INC += -I$(BUILD)

# Debugging/Optimization
ifdef DEBUG
COPT = -O0 -g
else
COPT = -Os -g #-fomit-frame-pointer -flto
endif

# compiler settings
CWARN = -Wall -Werror -Wno-error=cpp -Wno-error=format
CFLAGS = -DMICROPY_LIB $(INC) $(CWARN) -ansi -std=gnu99 $(CFLAGS_MOD) $(COPT) $(CFLAGS_EXTRA)

LDFLAGS = $(LDFLAGS_MOD) -lm $(LDFLAGS_EXTRA) -Wl,--nspireio

# source files
SRC_C = $(shell find . -name \*.c)

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

OBJ = $(PY_O) $(addprefix $(BUILD)/, $(SRC_C:.c=.o))

include ../py/mkrules.mk

all: $(PROG).tns
	ln -sf ../ulab .

distdir:
	/bin/cp *.c *.h *.mk Makefile mklib $(distdir)

$(PROG).tns: $(PROG)
	+genzehn --input $^ --output $@.zehn $(ZEHNFLAGS)
	+make-prg $@.zehn $@
	+rm $@.zehn
