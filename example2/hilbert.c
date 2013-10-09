#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <sndfile.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/time.h>
#endif
#include "libsoxrate.h"

typedef struct file_t {
    SNDFILE *ifile, *ofile;
    SF_INFO iinfo, oinfo;
} file_t;

#ifdef _WIN32
SNDFILE* sf_openx(const char *path, int mode, SF_INFO *sfinfo)
{
    int size;
    wchar_t *wp;
    size = MultiByteToWideChar(65001, 0, path, -1, 0, 0);
    wp = _alloca(size * sizeof(wchar_t));
    size = MultiByteToWideChar(65001, 0, path, -1, wp, size);
    return sf_wchar_open(wp, mode, sfinfo);
}
#else
#define sf_openx sf_open
#endif

void usage(void)
{
    fprintf(stderr, "usage: hilbert infile outfile\n");
    exit(1);
}

SNDFILE *open_input(const char *file, SF_INFO *info)
{
    SNDFILE *fp;
    if (!strcmp(file, "-")) {
	if (!(fp = sf_open_fd(0, SFM_READ, info, 0)))
	    fputs("ERROR: sf_open_fd()\n", stderr);
    } else {
	if (!(fp = sf_openx(file, SFM_READ, info)))
	    fprintf(stderr, "ERROR: sf_open()\n");
    }
    return fp;
}

SNDFILE *open_output(const char *file, SF_INFO *info)
{
    SNDFILE *fp;
    if (!file || !strcmp(file, "-")) {
	if (!(fp = sf_open_fd(1, SFM_WRITE, info, 0)))
	    fputs("ERROR: sf_open_fd()\n", stderr);
    } else {
	if (!(fp = sf_openx(file, SFM_WRITE, info)))
	    fprintf(stderr, "ERROR: sf_open()\n");
    }
    return fp;
}

static void hilbert(double *coefs, size_t numcoefs)
{
    size_t i, origin = (numcoefs - 1) >> 1;
    coefs[origin] = 0.0;
    for (i = 1; i <= origin; ++i) {
	double x = (i & 1) ? 1.0 / i : 0.0;
	coefs[origin + i] = -x;
	coefs[origin - i] = x;
    }
}

static void apply_hamming(double *coefs, size_t numcoefs)
{
    size_t i, origin = (numcoefs - 1) >> 1;
    for (i = 1; i <= origin; ++i) {
	double w = 0.54 - 0.46 *
	    cos(2.0 * M_PI * (origin + i) / (numcoefs - 1));
	coefs[origin + i] *= w;
	coefs[origin - i] *= w;
    }
}

static double calc_gain(double *coefs, size_t numcoefs)
{
    double gain = 0.0;
    size_t i, origin = (numcoefs - 1) >> 1;
    int odd = 0;
    for (i = origin + 1; i < numcoefs; i += 2) {
	if (odd) gain += coefs[i];
	else gain -= coefs[i];
	odd ^= 1;
    }
    if (odd) gain *= -1;
    gain *= 2.0;
    return gain;
}

double get_time()
{
#ifdef _WIN32
    return GetTickCount() / 1000.0;
#else
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_usec / 1000000.0 + tv.tv_sec;
#endif
}

void process_file(const file_t *file)
{
#define BS 0x4000
    float *ibuf = 0, *obuf = 0, *ip;
    sf_count_t count = 0, total = 0;
    size_t nch = file->iinfo.channels;
    int i, eof = 0;
    lsx_fir_t *conv = 0;
    double start, last, now;
    double *ht = 0;
    size_t numtaps = file->iinfo.samplerate / 12.0;
    double gain;
    if (!(numtaps & 1)) ++numtaps;

    if (!(ht = malloc(sizeof(double) * numtaps))) {
	fputs("ERROR: malloc()\n", stderr);
	goto done;
    }
    if (!(ibuf = malloc(sizeof(float) * BS * file->iinfo.channels))) {
	fputs("ERROR: malloc()\n", stderr);
	goto done;
    }
    if (!(obuf = malloc(sizeof(float) * BS * file->iinfo.channels))) {
	fputs("ERROR: malloc()\n", stderr);
	goto done;
    }
    hilbert(ht, numtaps);
    apply_hamming(ht, numtaps);
    gain = calc_gain(ht, numtaps);

    ip = ibuf;
    
    if (!(conv = lsx_fir_create(nch, ht, numtaps, (numtaps-1)/ 2, 1))) {
	fputs("ERROR: lsx_fir_create()\n", stderr);
	goto done;
    }
    if (lsx_fir_start(conv) < 0) {
	fputs("ERROR: lsx_fir_start()\n", stderr);
	goto done;
    }
    start = last = get_time();
    for (;;) {
	size_t ilen, olen;

	if (!eof && !count) {
	    count = sf_readf_float(file->ifile, ibuf, BS);
	    for (i = 0; i < count * nch; ++i)
		ibuf[i] /= gain;
	    total += count;
	    ip = ibuf;
	    if (!count) eof = 1;
	    now = get_time();
	    if (now - last > 0.1 || eof) {
		double ellapsed = now - start;
		double percent = 100.0 * total / file->iinfo.frames;
		double seconds = (double)total / file->iinfo.samplerate;
		last = now;
		fprintf(stderr, "\r%.1f%% (%.1fx) ",
		    percent, seconds / ellapsed);
		fflush(stderr);
	    }
	}
	ilen = count;
	olen = BS;
	if (lsx_fir_process(conv, ip, obuf, &ilen, &olen) < 0) {
	    fputs("\nERROR: lsx_rate_process()\n", stderr);
	    break;
	}
	if (!count && !olen) {
	    fputs("...done\n", stderr);
	    break;
	}
	ip += ilen * nch;
	count -= ilen;
	sf_writef_float(file->ofile, obuf, olen);
    }
done:
    if (conv) lsx_fir_close(conv);
    if (obuf) free(obuf);
    if (ibuf) free(ibuf);
    if (ht) free(ht);
#undef BS
}

#if defined(_MSC_VER) || defined(__MINGW32__)
int utf8_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
    file_t file = { 0 };
    int rc = 2;

    if (argc < 3) usage();

    if (!(file.ifile = open_input(argv[1], &file.iinfo)))
	goto done;

    file.oinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    file.oinfo.samplerate = file.iinfo.samplerate;
    file.oinfo.channels = file.iinfo.channels;

    if (!(file.ofile = open_output(argv[2], &file.oinfo)))
	goto done;

    process_file(&file);
    rc = 0;
done:
    if (file.ifile) sf_close(file.ifile);
    if (file.ofile) sf_close(file.ofile);
    return rc;
}
