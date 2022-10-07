include ../py/mkenv.mk
-include mpconfigport.mk

#PREFIX=/usr/bin/
CC = $(PREFIX)gcc # gcc-4.7 works if gcc is gcc-4.4 
CXX = $(PREFIX)g++
LD = $(PREFIX)ld 

# define main target
PROG = libmicropython.a

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
COPT =  -g
else
COPT =  -g #-fomit-frame-pointer -flto
endif

# compiler settings
CWARN = -Wall #-Werror -Wno-error=cpp -Wno-error=format
CFLAGS = -fPIC -DMICROPY_LIB $(INC) $(CWARN) -ansi -std=gnu99 $(CFLAGS_MOD) $(COPT) $(CFLAGS_EXTRA)



LDFLAGS = $(LDFLAGS_MOD) -lm $(LDFLAGS_EXTRA) 

# source files
SRC_C = $(shell find . -name \*.c)

SRC_QSTR += $(SRC_C)

OBJ = $(PY_O) $(addprefix $(BUILD)/, $(SRC_C:.c=.o))

include ../py/mkrules.mk

all: $(PROG)
	rm -f libmicropython.a
	ar cru libmicropython.a build/py/mpstate.o build/py/nlr.o build/py/nlrx86.o build/py/nlrx64.o build/py/nlrthumb.o build/py/nlrpowerpc.o build/py/nlrxtensa.o build/py/nlrsetjmp.o build/py/malloc.o build/py/gc.o build/py/pystack.o build/py/qstr.o build/py/vstr.o build/py/mpprint.o build/py/unicode.o build/py/mpz.o build/py/reader.o build/py/lexer.o build/py/parse.o build/py/scope.o build/py/compile.o build/py/emitcommon.o build/py/emitbc.o build/py/asmbase.o build/py/asmx64.o build/py/emitnx64.o build/py/asmx86.o build/py/emitnx86.o build/py/asmthumb.o build/py/emitnthumb.o build/py/emitinlinethumb.o build/py/asmarm.o build/py/emitnarm.o build/py/asmxtensa.o build/py/emitnxtensa.o build/py/emitinlinextensa.o build/py/emitnxtensawin.o build/py/formatfloat.o build/py/parsenumbase.o build/py/parsenum.o build/py/emitglue.o build/py/persistentcode.o build/py/runtime.o build/py/runtime_utils.o build/py/scheduler.o build/py/nativeglue.o build/py/ringbuf.o build/py/stackctrl.o build/py/argcheck.o build/py/warning.o build/py/profile.o build/py/map.o build/py/obj.o build/py/objarray.o build/py/objattrtuple.o build/py/objbool.o build/py/objboundmeth.o build/py/objcell.o build/py/objclosure.o build/py/objcomplex.o build/py/objdeque.o build/py/objdict.o build/py/objenumerate.o build/py/objexcept.o build/py/objfilter.o build/py/objfloat.o build/py/objfun.o build/py/objgenerator.o build/py/objgetitemiter.o build/py/objint.o build/py/objint_longlong.o build/py/objint_mpz.o build/py/objlist.o build/py/objmap.o build/py/objmodule.o build/py/objobject.o build/py/objpolyiter.o build/py/objproperty.o build/py/objnone.o build/py/objnamedtuple.o build/py/objrange.o build/py/objreversed.o build/py/objset.o build/py/objsingleton.o build/py/objslice.o build/py/objstr.o build/py/objstrunicode.o build/py/objstringio.o build/py/objtuple.o build/py/objtype.o build/py/objzip.o build/py/opmethods.o build/py/sequence.o build/py/stream.o build/py/binary.o build/py/builtinimport.o build/py/builtinevex.o build/py/builtinhelp.o build/py/modarray.o build/py/modbuiltins.o build/py/modcollections.o build/py/modgc.o build/py/modio.o build/py/modmath.o build/py/modcmath.o build/py/modmicropython.o build/py/modstruct.o build/py/modsys.o build/py/moduerrno.o build/py/modthread.o build/py/vm.o build/py/bc.o build/py/showbc.o build/py/repl.o build/py/smallint.o build/py/frozenmod.o build/extmod/moductypes.o build/extmod/modujson.o build/extmod/modure.o build/extmod/moduzlib.o build/extmod/moduheapq.o build/extmod/modutimeq.o build/extmod/moduhashlib.o build/extmod/moducryptolib.o build/extmod/modubinascii.o build/extmod/virtpin.o build/extmod/machine_mem.o build/extmod/machine_pinbase.o build/extmod/machine_signal.o build/extmod/machine_pulse.o build/extmod/machine_i2c.o build/extmod/machine_spi.o build/extmod/modbluetooth.o build/extmod/modussl_axtls.o build/extmod/modussl_mbedtls.o build/extmod/modurandom.o build/extmod/moduselect.o build/extmod/moduwebsocket.o build/extmod/modwebrepl.o build/extmod/modframebuf.o build/extmod/vfs.o build/extmod/vfs_blockdev.o build/extmod/vfs_reader.o build/extmod/vfs_posix.o build/extmod/vfs_posix_file.o build/extmod/vfs_fat.o build/extmod/vfs_fat_diskio.o build/extmod/vfs_fat_file.o build/extmod/vfs_lfs.o build/extmod/utime_mphal.o build/extmod/uos_dupterm.o build/lib/embed/abort_.o build/lib/utils/printf.o build/./unix_mphal.o build/./modos.o build/./gccollect.o  build/./input.o build/./file.o build/./graphic.o build/./main.o build/./nop.o build/modtime.o
	/bin/cp libmicropython.a ../..
distdir:
	/bin/cp *.c *.h *.mk Makefile $(distdir)
install: 
	/bin/cp libmicropython.a ../..

EMPTY_AUTOMAKE_TARGETS = dvi pdf ps info html tags ctags install-data install-exec uninstall install-dvi install-html install-info install-ps install-pdf installdirs check installcheck mostlyclean distclean maintainer-clean
.PHONY: $(EMPTY_AUTOMAKE_TARGETS)
$(EMPTY_AUTOMAKE_TARGETS):
