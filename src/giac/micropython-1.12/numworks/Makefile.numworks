include ../py/mkenv.mk
-include mpconfigport.mk

PREFIX=arm-none-eabi-
CC = $(PREFIX)gcc
CXX = $(PREFIX)g++
LD = $(PREFIX)ld 

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
COPT = -Ofast -g #-fomit-frame-pointer -flto
endif

# compiler settings
CWARN = -Wall -Werror -Wno-error=cpp -Wno-error=format
CFLAGS = -DNUMWORKS -DMICROPY_LIB $(INC) $(CWARN) -ansi -std=gnu99 $(CFLAGS_MOD) $(COPT) $(CFLAGS_EXTRA)
CFLAGS += -mthumb -march=armv7e-m -mfloat-abi=hard  -mcpu=cortex-m7 -mfpu=fpv5-sp-d16  -fdata-sections -ffunction-sections -fno-strict-aliasing -fno-exceptions # -fpermissive


LDFLAGS = $(LDFLAGS_MOD) -lm $(LDFLAGS_EXTRA) 

# source files
SRC_C = $(shell find . -name \*.c)

SRC_QSTR += $(SRC_C)

OBJ = $(PY_O) $(addprefix $(BUILD)/, $(SRC_C:.c=.o))

include ../py/mkrules.mk

all: $(PROG)

