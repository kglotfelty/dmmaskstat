tool = dmmaskstat

binexedir = $(prefix)/binexe
paramdir = $(prefix)/param
ahelpdir = $(prefix)/share/doc/xml

binexe_PROGRAMS = $(tool)

dmmaskstat_SOURCES = dmmaskstat.c t_dmmaskstat.c
dmmaskstat_CPPFLAGS = $(CIAO_CFLAGS)
dmmaskstat_LDADD = $(CIAO_LIBS) -ldmimgio -ldmfilters
dmmaskstat_LINK = $(CXX) -o $@ -Wl,-rpath,$(prefix)/lib -Wl,-rpath,$(prefix)/ots/lib 

# problems with libstdc++ and cxcparam on Linux :(
if LINUX
dmmaskstat_LDADD += -L$(prefix)/ots/lib -lstdc++
endif


dist_param_DATA = $(tool).par
dist_ahelp_DATA = $(tool).xml

install-data-hook:
	chmod a-w $(paramdir)/$(dist_param_DATA)
	chmod a-w $(ahelpdir)/$(dist_ahelp_DATA)
