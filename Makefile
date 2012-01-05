OBJS=alloc.o dft_filter.o effects_i_dsp.o fft4g.o lsxfft.o rate.o threaded_module.o

all: libsoxrate.so

CFLAGS = -fPIC -O

libsoxrate.so: $(OBJS)
	$(CC) -shared -o $@ $(OBJS) -lpthread

clean:
	$(RM) -f *.o
