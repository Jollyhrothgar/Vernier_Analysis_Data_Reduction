AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libVernierDSTReduction.la 

noinst_HEADERS = \
  VernierDSTReductionLinkDef.h 

libVernierDSTReduction_la_SOURCES = \
  VernierDSTReduction.cc \
  VernierDSTReduction_Dict.C

libVernierDSTReduction_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lSubsysReco

# Rule for generating CINT dictionaries from class headers.
VernierDSTReduction_Dict.C: VernierDSTReduction.h VernierDSTReductionLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*


