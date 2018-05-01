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
PLUGINS		=	../plugins/t5_3band_parameq_with_shelves.so
CC		=	cc

###############################################################################
#
# RULES TO BUILD PLUGINS FROM C OR C++ CODE
t5_3band_parameq_with_shelves:	plugins/t5_3band_parameq_with_shelves.c
	$(CC) $(CFLAGS) -o plugins/t5_3band_parameq_with_shelves.o -c plugins/t5_3band_parameq_with_shelves.c
	$(LD) -o ../plugins/t5_3band_parameq_with_shelves.so plugins/t5_3band_parameq_with_shelves.o -shared

###############################################################################
#
# TARGETS
#

install:	targets
	cp ../plugins/* $(INSTALL_PLUGINS_DIR)

targets:	t5_3band_parameq_with_shelves

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
