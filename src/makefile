###############################################################################
#
# Installation DIRECTORIES
#
# Change these if you want to install somewhere else.

INSTALL_PLUGINS_DIR	=	/usr/lib/ladspa/

###############################################################################
#
# GENERAL
#

INCLUDES	=	-I.
LIBRARIES	=	-ldl -lm
CFLAGS		=	$(INCLUDES) -Wall -O3 -fPIC
CXXFLAGS	=	$(CFLAGS)
PLUGINS		=	../plugins/biquad-eq-t5.so
CC		=	cc
CPP		=	c++

###############################################################################
#
# RULES TO BUILD PLUGINS FROM C OR C++ CODE
biquad-eq:	plugins/biquad-eq-t5.c
	$(CC) $(CFLAGS) -o plugins/biquad-eq-t5.o -c plugins/biquad-eq-t5.c
	$(LD) -o ../plugins/biquad-eq-t5.so plugins/biquad-eq-t5.o -shared

###############################################################################
#
# TARGETS
#

install:	targets
	cp ../plugins/* $(INSTALL_PLUGINS_DIR)

targets:	biquad-eq

###############################################################################
#
# UTILITIES
#

always:	

clean:
	-rm -f `find . -name "*.o"` ../bin/* ../plugins/*
	-rm -f `find .. -name "*~"`
	-rm -f *.bak core score.srt
	-rm -f *.bb *.bbg *.da *-ann gmon.out bb.out
	-rm -f `find .. -name "*.class"`

###############################################################################
