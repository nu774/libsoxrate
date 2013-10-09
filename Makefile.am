ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

lib_LTLIBRARIES = libsoxrate.la

libsoxrate_la_SOURCES = \
    alloc.c \
    dft_filter.c \
    effects_i_dsp.c \
    fft4g.c \
    lsxfft.c \
    rate.c \
    threaded_module.c

if HAVE_WINDRES
    libsoxrate_la_SOURCES += libsoxrate.rc
endif

libsoxrate_la_LDFLAGS = -no-undefined -export-symbols-regex '^lsx_.*$$'
libsoxrate_la_LIBADD = -lm

include_HEADERS = libsoxrate.h

.rc.lo: 
	$(LIBTOOL) $(AM_LIBATOOLFLAGS) $(LIBTOOLFLAGS) --tag=RC --mode=compile $(RC) $(RCFLAGS) $< -o $@

