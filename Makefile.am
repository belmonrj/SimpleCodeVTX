 AUTOMAKE_OPTIONS = foreign

 INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

 lib_LTLIBRARIES = libSimpleTreeSVXCNT.la

 noinst_HEADERS = \
   SimpleTreeSVXCNTLinkDef.h

 include_HEADERS = \
   SimpleTreeSVXCNT.h

 AM_LDFLAGS = \
   -L$(libdir) \
   -L$(OFFLINE_MAIN)/lib

 libSimpleTreeSVXCNT_la_LIBADD = \
   -L$(libdir) \
   -lCNT

 libSimpleTreeSVXCNT_la_SOURCES = \
   SimpleTreeSVXCNT.C \
   SimpleTreeSVXCNT_Dict.C

 BUILT_SOURCES = \
   testexternals.C

 noinst_PROGRAMS = \
   testexternals

 testexternals_LDADD = \
   libSimpleTreeSVXCNT.la

 testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

 SimpleTreeSVXCNT_Dict.C: SimpleTreeSVXCNT.h SimpleTreeSVXCNTLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^
