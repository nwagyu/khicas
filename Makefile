Q ?= @
CC = arm-none-eabi-gcc
CXX = arm-none-eabi-g++
AR = arm-none-eabi-ar
RANLIB = arm-none-eabi-ranlib
STRIP = arm-none-eabi-strip
NWLINK = nwlink

DEBUG = 0
LINK_GC = 1
LTO = 1

LIBS_PATH = $(shell pwd)/output/libs

CFLAGS = $(shell $(NWLINK) eadk-cflags)
CFLAGS += -fno-strict-aliasing -fno-threadsafe-statics
CFLAGS += -fno-exceptions -fno-use-cxa-atexit

LDFLAGS = -Wl,--relocatable
LDFLAGS += -nostartfiles
LDFLAGS += --specs=nosys.specs
LDFLAGS += -L$(LIBS_PATH)/lib -lc_nano -lstdc++_nano -lm -lgmp -lmpfr -lmpfi

ifeq ($(DEBUG),1)
CFLAGS += -ggdb
else
CFLAGS += -Os -DNDEBUG
endif

ifeq ($(LINK_GC),1)
CFLAGS += -fdata-sections -ffunction-sections
LDFLAGS += -Wl,-e,main -Wl,-u,eadk_app_name -Wl,-u,eadk_app_icon -Wl,-u,eadk_api_level
LDFLAGS += -Wl,--gc-sections
endif

ifeq ($(LTO),1)
AR = arm-none-eabi-gcc-ar
RANLIB = arm-none-eabi-gcc-ranlib
CFLAGS += -flto -fno-fat-lto-objects
CFLAGS += -flto-partition=none
CFLAGS += -fwhole-program
CFLAGS += -fvisibility=internal
LDFLAGS += -flinker-output=nolto-rel
endif

.PHONY: build
build: output/khicas.nwa

###
# External libs
###

LIBS_CFLAGS = $(CFLAGS) -I$(LIBS_PATH)/include --specs=nosys.specs

AUTOCONF_FLAGS = --host=arm-none-eabi
AUTOCONF_FLAGS += --prefix=$(LIBS_PATH)
AUTOCONF_FLAGS += CFLAGS="$(LIBS_CFLAGS)"
ifeq ($(LTO),1)
AUTOCONF_FLAGS += AR="arm-none-eabi-gcc-ar" RANLIB="arm-none-eabi-gcc-ranlib"
endif

libs = $(addprefix output/libs/lib/,libgmp.a libmpfr.a libmpfi.a)

define autoconf
output/libs/lib/lib$(1).a:
	@mkdir -p output/autoconf/$(1)
	@echo "AUTOCNF $(1)"
	$(Q) cd output/autoconf/$(1) && ../../../src/$(1)/configure $(2) $(AUTOCONF_FLAGS) > configure.log 2>&1
	@echo "MAKE    $(1)"
	$(Q) cd output/autoconf/$(1) && make install > make.log 2>&1
endef

$(eval $(call autoconf,gmp,--disable-assembly))
$(eval $(call autoconf,mpfr,--with-gmp=$(LIBS_PATH) --disable-shared))
$(eval $(call autoconf,mpfi,--with-gmp=$(LIBS_PATH) --with-mpfr=$(LIBS_PATH)))

###
# Giac
###

giac_objs = $(addprefix output/giac/,sym2poly.o gausspol.o threaded.o maple.o ti89.o mathml.o moyal.o misc.o permu.o desolve.o input_parser.o symbolic.o index.o modpoly.o modfactor.o ezgcd.o derive.o solve.o intg.o intgab.o risch.o lin.o series.o subst.o vecteur.o sparse.o csturm.o tex.o global.o ifactor.o alg_ext.o gauss.o isom.o help.o plot.o plot3d.o rpn.o prog.o pari.o cocoa.o unary.o usual.o identificateur.o gen.o input_lexer.o tinymt32.o first.o quater.o kdisplay.o kadd.o k_csdk.o)

GIAC_CFLAGS = -w -Wno-psabi -Wno-deprecated-declarations -Wno-c++20-extensions] # Silence warnings
GIAC_CFLAGS += -DDEVICE -DNUMWORKS -DHAVE_CONFIG_H -DIN_GIAC -DTIMEOUT -DNO_STDEXCEPT -DGIAC_BINARY_ARCHIVE -DNO_UNARY_FUNCTION_COMPOSE -DNO_TEMPLATE_MULTGCD -DGIAC_NO_OPTIMIZATIONS -DSTATIC_BUILTIN_LEXER_FUNCTIONS # -DMICROPY_LIB

.PHONY: patch_kdisplay
patch_kdisplay:
	@echo "PATCH   src/giac/src/kdisplay.cc"
	$(Q) patch --forward --strip=1 --reject-file=- < src/giac_src_kdisplay.patch || true

output/giac/kdisplay.o: patch_kdisplay

output/giac/global.o output/giac/kdisplay.o: GIAC_CFLAGS := $(filter-out -DDEVICE, $(GIAC_CFLAGS))

# We need to do a weird header danse because a "config.h" file is versionned in giac's source code...
output/giac/%.o: src/giac/src/%.cc $(libs)
	@mkdir -p output/giac
	$(Q) mv src/giac/src/config.h src/giac/src/config.h.original
	$(Q) cp src/giac/src/config.h.numworks.gmp src/giac/src/config.h
	@echo "CXX     $<"
	$(Q) $(CXX) $(LIBS_CFLAGS) $(GIAC_CFLAGS) -c $< -o $@
	$(Q) mv src/giac/src/config.h.original src/giac/src/config.h
output/giac/%.o: src/giac/src/%.c $(libs)
	@mkdir -p output/giac
	$(Q) mv src/giac/src/config.h src/giac/src/config.h.original
	$(Q) cp src/giac/src/config.h.numworks.gmp src/giac/src/config.h
	@echo "CC      $<"
	$(Q) $(CC) $(LIBS_CFLAGS) $(GIAC_CFLAGS) -c $< -o $@
	$(Q) mv src/giac/src/config.h.original src/giac/src/config.h

output/main.o: src/main.cpp
	@mkdir -p output
	@echo "CXX     $<"
	$(Q) $(CXX) $(CFLAGS)  -c $< -o $@

.PHONY: build
build: output/khicas.nwa

.PHONY: check
check: output/khicas.bin

.PHONY: run
run: output/khicas.nwa
	@echo "INSTALL $<"
	$(Q) $(NWLINK) install-nwa $<

output/%.bin: output/%.nwa
	@echo "BIN     $@"
	$(Q) $(NWLINK) nwa-bin $< $@

output/%.elf: output/%.nwa
	@echo "ELF     $@"
	$(Q) $(NWLINK) nwa-elf $< $@

output/khicas.nwa: output/main.o output/icon.o $(giac_objs)
	@echo "LD      $@"
	$(Q) $(CC) $(CFLAGS) $(LDFLAGS) $^ -lc_nano -lstdc++_nano -lm -lmpfi -lmpfr -lgmp -o $@
	@echo "STRIP   $@"
	$(Q) $(STRIP) --strip-debug $@

output/icon.o: src/icon.png
	@mkdir -p output
	@echo "ICON    $<"
	$(Q) $(NWLINK) png-icon-o $< $@

.PHONY: clean
clean:
	@echo "CLEAN"
	$(Q) rm -rf output
