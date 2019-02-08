##############################################################################

MK_TOP = ../../../..

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmmaskstat
LIB_FILES = 
PAR_FILES         = dmmaskstat.par
INC_FILES         = 
XML_FILES         = dmmaskstat.xml

SRCS	=           dmmaskstat.c t_dmmaskstat.c


OBJS = $(SRCS:.c=.o)
LOCAL_INC = -I../dmimgfilt
MAKECLEAN += filters.o


MAKETEST_SCRIPT   = dmmaskstat.t



include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------


$(EXEC): $(OBJS) filters.o
	$(LINK) filters.o
	@echo

filters.o: ../dmimgfilt/filters.c
	rm -f filters.o
	(cd ../dmimgfilt; $(MAKE) $(MKMACROS) filters.o )
	cp -f ../dmimgfilt/filters.o .


announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |                Building dmmaskstat                        | "
	@echo "   \----------------------------------------------------------/ "

