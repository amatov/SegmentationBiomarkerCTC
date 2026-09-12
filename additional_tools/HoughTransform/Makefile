###############################################################################
# Houghtool - The Software Package for efficiency measuring and visualization
# of the HT and it's variants for line detection
#
#
# Contents of the package:
#   Standard Hough Transform (SHT):				sht_line
#   Randomized Hough Transform (RHT):				rht_line
#   Dynamic Randomized Hough Transform (DRHT):			drht_line
#   Window Randomized Hough Transform (WRHT):			wrht_line
#   Random Window Randomized Hough Transform (RWRHT):		rwrht_line
#   Curve Fitting Hough Transform (CFHT):			cfht_line
#   Probabilistic Hough Transform of Kiryati et al. (ProbHT):	probht_line
#   Adaptive Hough Transform (AHT):				aht_line
#   Combinatorial Hough Transform (CHT):			cht_line
#   Dynamic Combinatorial Hough Transform (DCHT):		dcht_line
#
#   Graphical X Window System User Interface:			xhoughtool
#
# Makefile usage: make [programname|all|clean|help]
#                   programname is one of the names listed above, e.g. make rht_line. 
#                   Without programname or with the option "all" all programs are compiled. 
#                   The option "clean" deletes core and all object files. Help is 
#                   available with the option "help".
#
#
# Address: Lappeenranta University of Technology
#          Department of Information Technology
#          Laboratory of Information Processing
#          P.O. Box 20
#          FIN-53851 Lappeenranta
#          Finland
#
# Authors: Jouni Ikonen (Jouni.Ikonen@lut.fi), Petri Hirvonen 
#          Pekka Kultanen, Heikki Kalviainen (Heikki.Kalviainen@lut.fi)
#           
#
# File:    	Makefile
# Date:    	Jun 1, 1993
# Last change: 	Oct 10,1995
#
###############################################################################
CC		= cc
CFLAGS		= -g
#
# Define BIG_IMAGES and you can use images up to size of 512x512 (def 256x256)
#
#CFLAGS		= -O4 -DBIG_IMAGES
#CFLAGS		= -O4 -temp=temp
#
# Use following CFLAGS when the linker doesn't find nint() and/or rint()
#
#CFLAGS		= -O -DNON_ANSI_RECOVERY

XINCCFLAGS	= -DVISUAL_PACK
XCFLAGS		= -I/usr/local/packages/OpenLook_3.0/share/include
LIBS		= -lm
XLIBS		= -lxview -lolgx -lX11
PROGRAMS	= sht_line rht_line drht_line wrht_line rwrht_line\
		  cfht_line probht_line aht_line cht_line dcht_line\
		  xhoughtool
X_OBJECTS	= xhoughtool.o graphics.o drawaccu.o
XHT_OBJECTS	= x_ht_fndmaxline.o\
		  x_sht_linedetect.o x_sht_fndmaxline.o\
		  x_rht_linedetect.o rht_infmat.o\
		  rht_PSfuncs.o rht_dist.o rht_imgutils.o\
		  x_drht_linedetect.o drht_blockpixels.o x_wrht_linedetect.o\
		  x_rwrht_linedetect.o x_rwrht_imgutils.o\
		  x_cfht_linedetect.o cfht_infmat.o cfht_fndmaxline.o\
		  probht_imgutils.o x_aht_linedetect.o aht_fndmaxline.o\
		  x_cht_linedetect.o x_cht_fndmaxline.o x_dcht_linedetect.o
FILT_OBJECTS	= ht_fndmaxline.o ht_filtutils.o
HT_OBJECTS	= ht_ldaimg.o ht_stoimg.o ht_imgutils.o ht_random.o\
		  ht_statics.o ht_nonansi.o
SHT_OBJECTS	= sht_lineserver.o sht_linedetect.o sht_fndmaxline.o
RHT_OBJECTS	= rht_lineserver.o rht_linedetect.o rht_infmat.o\
		  rht_PSfuncs.o rht_dist.o rht_imgutils.o
DRHT_OBJECTS	= drht_lineserver.o drht_linedetect.o drht_blockpixels.o\
		  rht_infmat.o rht_PSfuncs.o rht_dist.o rht_imgutils.o
WRHT_OBJECTS	= wrht_lineserver.o wrht_linedetect.o rht_infmat.o\
		  rht_PSfuncs.o rht_imgutils.o
RWRHT_OBJECTS	= rwrht_lineserver.o rwrht_linedetect.o rwrht_imgutils.o\
		  rht_infmat.o rht_PSfuncs.o rht_dist.o rht_imgutils.o
CFHT_OBJECTS	= cfht_lineserver.o cfht_linedetect.o cfht_infmat.o\
		  cfht_fndmaxline.o rht_infmat.o rht_imgutils.o\
		  rht_PSfuncs.o
ProbHT_OBJECTS	= probht_lineserver.o probht_imgutils.o sht_linedetect.o\
		  sht_fndmaxline.o
AHT_OBJECTS	= aht_lineserver.o aht_linedetect.o aht_fndmaxline.o
CHT_OBJECTS	= cht_lineserver.o cht_linedetect.o cht_fndmaxline.o
DCHT_OBJECTS	= dcht_lineserver.o dcht_linedetect.o rht_imgutils.o
X_HFILES	= xhoughtool.h
HT_HFILES	= ht_hough.h ht_formats.h ht_visilog.h ht_cvl.h\
		  ht_graphmacros.h
RHT_HFILES	= rht_infmat.h

all:			$(PROGRAMS)

sht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(SHT_OBJECTS)
			$(CC) -o sht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
				       $(SHT_OBJECTS) $(LIBS)

rht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(RHT_OBJECTS)
			$(CC) -o rht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
				       $(RHT_OBJECTS) $(LIBS)

drht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(DRHT_OBJECTS)
			$(CC) -o drht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
					$(DRHT_OBJECTS) $(LIBS)

wrht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(WRHT_OBJECTS)
			$(CC) -o wrht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
					$(WRHT_OBJECTS) $(LIBS)

rwrht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(RWRHT_OBJECTS)
			$(CC) -o rwrht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
					 $(RWRHT_OBJECTS) $(LIBS)

cfht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(CFHT_OBJECTS)
			$(CC) -o cfht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
					$(CFHT_OBJECTS) $(LIBS)

probht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(ProbHT_OBJECTS)
			$(CC) -o probht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
					  $(ProbHT_OBJECTS) $(LIBS)

aht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(AHT_OBJECTS)
			$(CC) -o aht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
				       $(AHT_OBJECTS) $(LIBS)

cht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(CHT_OBJECTS)
			$(CC) -o cht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
				       $(CHT_OBJECTS) $(LIBS)

dcht_line:		$(FILT_OBJECTS) $(HT_OBJECTS) $(DCHT_OBJECTS)
			$(CC) -o dcht_line $(FILT_OBJECTS) $(HT_OBJECTS)\
					$(DCHT_OBJECTS) $(LIBS)

xhoughtool:		$(XHT_OBJECTS) $(X_OBJECTS) $(HT_OBJECTS)
			$(CC) -o xhoughtool $(XHT_OBJECTS) $(X_OBJECTS)\
					 $(HT_OBJECTS) $(LIBS) $(XLIBS)

xhoughtool.o:		xhoughtool.c xhoughtool.icon
			$(CC) $(CFLAGS) $(XCFLAGS) -c  xhoughtool.c

graphics.o:		graphics.c
			$(CC) $(CFLAGS) $(XCFLAGS) -c  graphics.c

drawaccu.o:		drawaccu.c
			$(CC) $(CFLAGS) $(XCFLAGS) -c  drawaccu.c

ht_filtutils.o:		ht_filtutils.c
			$(CC) $(CFLAGS) -c  ht_filtutils.c

ht_ldaimg.o:		ht_ldaimg.c
			$(CC) $(CFLAGS) -c  ht_ldaimg.c

ht_stoimg.o:		ht_stoimg.c
			$(CC) $(CFLAGS) -c  ht_stoimg.c

ht_imgutils.o:		ht_imgutils.c
			$(CC) $(CFLAGS) -c  ht_imgutils.c

ht_fndmaxline.o:	ht_fndmaxline.c
			$(CC) $(CFLAGS) -c  ht_fndmaxline.c

x_ht_fndmaxline.o:	ht_fndmaxline.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  ht_fndmaxline.c -o x_ht_fndmaxline.o

ht_random.o:		ht_random.c
			$(CC) $(CFLAGS) -c  ht_random.c

ht_statics.o:		ht_statics.c
			$(CC) $(CFLAGS) -c  ht_statics.c

ht_nonansi.o:		ht_nonansi.c
			$(CC) $(CFLAGS) -c  ht_nonansi.c

sht_lineserver.o:	sht_lineserver.c
			$(CC) $(CFLAGS) -c  sht_lineserver.c

sht_linedetect.o:	sht_linedetect.c
			$(CC) $(CFLAGS) -c  sht_linedetect.c

x_sht_linedetect.o:	sht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  sht_linedetect.c -o x_sht_linedetect.o

sht_fndmaxline.o:	sht_fndmaxline.c
			$(CC) $(CFLAGS) -c  sht_fndmaxline.c

x_sht_fndmaxline.o:	sht_fndmaxline.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  sht_fndmaxline.c -o x_sht_fndmaxline.o

rht_lineserver.o:	rht_lineserver.c
			$(CC) $(CFLAGS) -c  rht_lineserver.c

rht_linedetect.o:	rht_linedetect.c
			$(CC) $(CFLAGS) -c  rht_linedetect.c

x_rht_linedetect.o:	rht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  rht_linedetect.c -o x_rht_linedetect.o

rht_infmat.o:		rht_infmat.c
			$(CC) $(CFLAGS) -c  rht_infmat.c

rht_PSfuncs.o:		rht_PSfuncs.c
			$(CC) $(CFLAGS) -c  rht_PSfuncs.c

rht_dist.o:		rht_dist.c
			$(CC) $(CFLAGS) -c  rht_dist.c

rht_imgutils.o:		rht_imgutils.c
			$(CC) $(CFLAGS) -c  rht_imgutils.c

drht_lineserver.o:	drht_lineserver.c
			$(CC) $(CFLAGS) -c  drht_lineserver.c

drht_linedetect.o:	drht_linedetect.c
			$(CC) $(CFLAGS) -c  drht_linedetect.c

x_drht_linedetect.o:	drht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  drht_linedetect.c -o x_drht_linedetect.o

drht_blockpixels.o:	drht_blockpixels.c
			$(CC) $(CFLAGS) -c  drht_blockpixels.c

wrht_lineserver.o:	wrht_lineserver.c
			$(CC) $(CFLAGS) -c  wrht_lineserver.c

wrht_linedetect.o:	wrht_linedetect.c
			$(CC) $(CFLAGS) -c  wrht_linedetect.c

x_wrht_linedetect.o:	wrht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  wrht_linedetect.c -o x_wrht_linedetect.o

rwrht_lineserver.o:	rwrht_lineserver.c
			$(CC) $(CFLAGS) -c  rwrht_lineserver.c

rwrht_linedetect.o:	rwrht_linedetect.c
			$(CC) $(CFLAGS) -c  rwrht_linedetect.c

x_rwrht_linedetect.o:	rwrht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  rwrht_linedetect.c -o x_rwrht_linedetect.o

rwrht_imgutils.o:	rwrht_imgutils.c
			$(CC) $(CFLAGS) -c  rwrht_imgutils.c

x_rwrht_imgutils.o:	rwrht_imgutils.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  rwrht_imgutils.c -o x_rwrht_imgutils.o

cfht_lineserver.o:	cfht_lineserver.c
			$(CC) $(CFLAGS) -c  cfht_lineserver.c

cfht_linedetect.o:	cfht_linedetect.c
			$(CC) $(CFLAGS) -c  cfht_linedetect.c

x_cfht_linedetect.o:	cfht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  cfht_linedetect.c -o x_cfht_linedetect.o

cfht_infmat.o:		cfht_infmat.c
			$(CC) $(CFLAGS) -c  cfht_infmat.c

cfht_fndmaxline.o:	cfht_fndmaxline.c
			$(CC) $(CFLAGS) -c  cfht_fndmaxline.c

probht_lineserver.o:	probht_lineserver.c
			$(CC) $(CFLAGS) -c  probht_lineserver.c

probht_imgutils.o:	probht_imgutils.c
			$(CC) $(CFLAGS) -c  probht_imgutils.c

aht_lineserver.o:	aht_lineserver.c
			$(CC) $(CFLAGS) -c  aht_lineserver.c

aht_linedetect.o:	aht_linedetect.c
			$(CC) $(CFLAGS) -c  aht_linedetect.c

x_aht_linedetect.o:	aht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  aht_linedetect.c -o x_aht_linedetect.o

aht_fndmaxline.o:	aht_fndmaxline.c
			$(CC) $(CFLAGS) -c  aht_fndmaxline.c

cht_lineserver.o:	cht_lineserver.c
			$(CC) $(CFLAGS) -c  cht_lineserver.c

cht_linedetect.o:	cht_linedetect.c
			$(CC) $(CFLAGS) -c  cht_linedetect.c

x_cht_linedetect.o:	cht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  cht_linedetect.c -o x_cht_linedetect.o

cht_fndmaxline.o:	cht_fndmaxline.c
			$(CC) $(CFLAGS) -c  cht_fndmaxline.c

x_cht_fndmaxline.o:	cht_fndmaxline.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  cht_fndmaxline.c -o x_cht_fndmaxline.o

dcht_lineserver.o:	dcht_lineserver.c
			$(CC) $(CFLAGS) -c  dcht_lineserver.c

dcht_linedetect.o:	dcht_linedetect.c
			$(CC) $(CFLAGS) -c  dcht_linedetect.c

x_dcht_linedetect.o:	dcht_linedetect.c
			$(CC) $(CFLAGS) $(XINCCFLAGS) $(XCFLAGS)\
			 -c  dcht_linedetect.c -o x_dcht_linedetect.o


$(FILT_OBJECTS):	$(HT_HFILES)

$(XHT_OBJECTS):		$(HT_HFILES) $(RHT_HFILES)

$(X_OBJECTS):		$(HT_HFILES) $(X_HFILES)

$(HT_OBJECTS):		$(HT_HFILES)

$(SHT_OBJECTS):		$(HT_HFILES)

$(RHT_OBJECTS):		$(HT_HFILES) $(RHT_HFILES)

$(DRHT_OBJECTS):	$(HT_HFILES) $(RHT_HFILES)

$(WRHT_OBJECTS):	$(HT_HFILES) $(RHT_HFILES)

$(RWRHT_OBJECTS):	$(HT_HFILES) $(RHT_HFILES)

$(CFHT_OBJECTS):	$(HT_HFILES) $(RHT_HFILES)

$(ProbHT_OBJECTS):	$(HT_HFILES)

$(AHT_OBJECTS):		$(HT_HFILES)

$(CHT_OBJECTS):		$(HT_HFILES)

$(DCHT_OBJECTS):	$(HT_HFILES)

clean:
			$(RM) $(PROGRAMS) *.o a.out core

help:
			@echo "Makefile usage: make [programname|all|clean|help]"
			@echo
			@echo "   programname is a name of the method, e.g. make rht_line."
			@echo "   Without programname or with the option "all" all programs are compiled."
			@echo "   The option "clean" deletes core and all object files."
			@echo 
			@echo "   Options to programname: "
			@echo "   Standard Hough Transform (SHT):			sht_line"
			@echo "   Randomized Hough Transform (RHT):			rht_line"
			@echo "   Dynamic Randomized Hough Transform (DRHT):		drht_line"
			@echo "   Window Randomized Hough Transform (WRHT):		wrht_line"
			@echo "   Random Window Randomized Hough Transform (RWRHT):	rwrht_line"
			@echo "   Curve Fitting Hough Transform (CFHT):		cfht_line
			@echo "   Probabilistic Hough Transform by Kiryati et al. (ProbHT): probht_line"
			@echo "   Adaptive Hough Transform (AHT):			aht_line"
			@echo "   Combinatorial Hough Transform (CHT):			cht_line"
			@echo "   Dynamic Combinatorial Hough Transform (DCHT):	dcht_line"
			@echo "   Graphical X Window System User Interface:		xhoughtool"
			@echo
			@echo "Objects needed to link the programs:"
			@echo
			@echo "ALL: $(HT_OBJECTS)"
			@echo
			@echo "Filters: $(FILT_OBJECTS)"
			@echo
			@echo "sht: $(SHT_OBJECTS)"
			@echo
			@echo "rht: $(RHT_OBJECTS)"
			@echo
			@echo "drht: $(DRHT_OBJECTS)"
			@echo
			@echo "wrht: $(WRHT_OBJECTS)"
			@echo
			@echo "rwrht: $(RWRHT_OBJECTS)"
			@echo
			@echo "cfht: $(CFHT_OBJECTS)"
			@echo
			@echo "probht: $(ProbHT_OBJECTS)"
			@echo
			@echo "aht: $(AHT_OBJECTS)"
			@echo
			@echo "cht: $(CHT_OBJECTS)"
			@echo
			@echo "dcht: $(DCHT_OBJECTS)"
			@echo
			@echo "xhoughtool: $(XHT_OBJECTS) $(X_OBJECTS)"
