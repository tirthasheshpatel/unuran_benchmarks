## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

SUBDIRS = \
 utils \
 methods \
 specfunct \
 distr \
 distributions \
 parser \
 tests \
 uniform \
 urng \
 .

DIST_SUBDIRS = $(SUBDIRS)

SUBLIBS = \
 methods/libmethods.la \
 distr/libdistr.la \
 distributions/libdistributions.la \
 specfunct/libspecfunct.la \
 parser/libparser.la \
 tests/libtests.la \
 uniform/libuniform.la \
 urng/liburng.la \
 utils/libutils.la

lib_LTLIBRARIES = libunuran.la
libunuran_la_SOURCES = 
libunuran_la_LIBADD = $(SUBLIBS)
libunuran_la_LDFLAGS = -version-info $(UNURAN_LT_VERSION)

include_HEADERS = unuran.h

noinst_HEADERS = \
 unuran_config.h \
 unuran.h.in \
 unur_cookies.h \
 unur_source.h \
 unur_struct.h \
 unur_typedefs.h


# make unuran.h (process file unuran.h.in)
sinclude .dep-unuran_h
unuran.h: unuran.h.in 
	$(top_srcdir)/scripts/merge_h.pl $(top_srcdir)/src/unuran.h.in > unuran.h
	sed -e "s#^unuran\.h\:#stringparser\.c\:#" .dep-unuran_h | sed -e "s#\./#\.\./#g" > ./parser/.dep-stringparser_c

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 .dep-unuran_h \
 unuran.h
