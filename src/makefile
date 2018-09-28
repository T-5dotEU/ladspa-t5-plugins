INSTALL_PLUGINS_DIR	=	/usr/lib/ladspa/

INCLUDES	=	-I.
LIBRARIES	=	-ldl -lm
CFLAGS		=	$(INCLUDES) -Wall -O3 -fPIC
CXXFLAGS	=	$(CFLAGS)
CC			=	cc

targets: t5_lr4_lowpass t5_lr4_highpass t5_3band_parameq_with_shelves

install:	targets
	cp ../plugins/* $(INSTALL_PLUGINS_DIR)

t5_3band_parameq_with_shelves:	plugins/t5_3band_parameq_with_shelves.c
	$(CC) $(CFLAGS) -o plugins/t5_3band_parameq_with_shelves.o -c plugins/t5_3band_parameq_with_shelves.c
	$(LD) -o ../plugins/t5_3band_parameq_with_shelves.so plugins/t5_3band_parameq_with_shelves.o -shared

t5_lr4_lowpass:
	$(CC) $(CFLAGS) -o plugins/t5_lr4_lowpass.o -c plugins/t5_lr4_lowpass.c
	$(LD) -o ../plugins/t5_lr4_lowpass.so plugins/t5_lr4_lowpass.o -shared

t5_lr4_highpass:
	$(CC) $(CFLAGS) -o plugins/t5_lr4_highpass.o -c plugins/t5_lr4_highpass.c
	$(LD) -o ../plugins/t5_lr4_highpass.so plugins/t5_lr4_highpass.o -shared

always:	

clean:
	-rm -f `find . -name "*.o"` ../bin/* ../plugins/*
