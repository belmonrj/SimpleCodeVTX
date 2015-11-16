 AUTOMAKE_OPTIONS = foreign

 INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include 

 lib_LTLIBRARIES = libtrinn_SimpleTreeSVXCNT.la

noinst_HEADERS = \
   SimpleTreeSVXCNTLinkDef.h
 
 include_HEADERS = \
   SimpleTreeSVXCNT.h

 AM_LDFLAGS = \
   -L$(libdir) \
   -L$(OFFLINE_MAIN)/lib

 libtrinn_SimpleTreeSVXCNT_la_LIBADD = \
   -L$(libdir) \
   -lCNT

 libtrinn_SimpleTreeSVXCNT_la_SOURCES = \
   SimpleTreeSVXCNT.C \
   SimpleTreeSVXCNT_Dict.C
 
 BUILT_SOURCES = \
   testexternals.C

 noinst_PROGRAMS = \
   testexternals

 testexternals_LDADD = \
   libtrinn_SimpleTreeSVXCNT.la

 testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

 SimpleTreeSVXCNT_Dict.C: SimpleTreeSVXCNT.h SimpleTreeSVXCNTLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
