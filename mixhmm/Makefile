#SHELL = /bin/sh
#.SUFFIXES: .c .o
COMMON_DEPS = util.h


CC = gcc
LD = gcc
CCFLAGS = -Wall -g
#CFLAGS = `perl -MExtUtils::Embed -e ccopts` -fPIC
LDFLAGS = `perl -MExtUtils::Embed -e ldopts` -lm
LDFLAGS = -lm
#ARCH = `perl -MConfig -e 'print $$Config{archname}'`

all: utest

utest: runner
	./runner
#include/arrays.o
runner: test_biot.o kc.o include/nrutil.o
	$(LD) $(LDFLAGS) -o $@ $^

%.o: %.c $(COMMON_DEPS)
	$(CC) $(CCFLAGS) -c $<

#khmm.so: khmm_wrap.o khmm.o kc.o include/nrutil.o include/arrays.o
#    $(LD) $(LDFLAGS) -shared -o $@ $^
#    mkdir -p $(ARCH)/
#    mkdir -p $(ARCH)/auto/
#    mv $@ $(ARCH)/

