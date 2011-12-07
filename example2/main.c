#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <windows.h>
#include "sndfile.h"
#include "libsoxrate.h"

typedef struct file_t {
    SNDFILE *ifile, *ofile;
    SF_INFO iinfo, oinfo;
} file_t;

struct {
    HMODULE handle;
    SNDFILE *(*open_fd)(int, int, SF_INFO *, int);
    const char *(*strerror)(SNDFILE *);
    int (*command)(SNDFILE *, int, void *, int);
    sf_count_t (*readf_float)(SNDFILE *, float *, sf_count_t);
    sf_count_t (*writef_float)(SNDFILE *, const float *, sf_count_t);
    int (*close)(SNDFILE *);
    SNDFILE *(*wchar_open)(const wchar_t *, int, SF_INFO *);
} sndfile;

void load_libsndfile(void)
{
#define CHECK(expr) do { if (!(expr)) goto error_end_1; } while (0)
#define FETCH(x) \
    CHECK(sndfile.x = (void*)(GetProcAddress(sndfile.handle,"sf_" #x)))

    if (!(sndfile.handle = LoadLibraryW(L"libsndfile-1.dll")))
	goto error_end;
    FETCH(open_fd);
    FETCH(strerror);
    FETCH(command);
    FETCH(readf_float);
    FETCH(writef_float);
    FETCH(close);
    FETCH(wchar_open);
    return;
#undef CHECK
#undef FETCH
error_end_1:
    FreeLibrary(sndfile.handle);
error_end:
    fputs("can't load libsndfile-1.dll\n", stderr);
    exit(2);
}

void usage(void)
{
    fprintf(stderr, "usage: hilbert infile outfile\n");
    exit(1);
}

SNDFILE *open_input(const wchar_t *file, SF_INFO *info)
{
    SNDFILE *fp;
    if (!wcscmp(file, L"-")) {
	if (!(fp = sndfile.open_fd(0, SFM_READ, info, 0)))
	    fputs("ERROR: sf_open_fd()\n", stderr);
    } else {
	if (!(fp = sndfile.wchar_open(file, SFM_READ, info)))
	    fprintf(stderr, "ERROR: sf_wchar_open(): %ls\n", file);
    }
    return fp;
}

SNDFILE *open_output(const wchar_t *file, SF_INFO *info)
{
    SNDFILE *fp;
    if (!file || !wcscmp(file, L"-")) {
	if (!(fp = sndfile.open_fd(1, SFM_WRITE, info, 0)))
	    fputs("ERROR: sf_open_fd()\n", stderr);
    } else {
	if (!(fp = sndfile.wchar_open(file, SFM_WRITE, info)))
	    fprintf(stderr, "ERROR: sf_wchar_open(): %ls\n", file);
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

void process_file(const file_t *file)
{
#define BS 0x4000
    float *ibuf = 0, *obuf = 0, *ip;
    sf_count_t count = 0, total = 0;
    size_t nch = file->iinfo.channels;
    int i, eof = 0;
    lsx_fir_t *conv = 0;
    DWORD start, last, now;
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
    start = last = GetTickCount();
    for (;;) {
	size_t ilen, olen;

	if (!eof && !count) {
	    count = sndfile.readf_float(file->ifile, ibuf, BS);
	    for (i = 0; i < count * nch; ++i)
		ibuf[i] /= gain;
	    total += count;
	    ip = ibuf;
	    if (!count) eof = 1;
	    now = GetTickCount();
	    if (now - last > 100 || eof) {
		double ellapsed = (double)(now - start) / 1000.0;
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
	sndfile.writef_float(file->ofile, obuf, olen);
    }
done:
    if (conv) lsx_fir_close(conv);
    if (obuf) free(obuf);
    if (ibuf) free(ibuf);
    if (ht) free(ht);
#undef BS
}

int wmain(int argc, wchar_t **argv)
{
    file_t file = { 0 };

    if (argc < 3) usage();
    load_libsndfile();

    if (!(file.ifile = open_input(argv[1], &file.iinfo)))
	goto done;

    file.oinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    file.oinfo.samplerate = file.iinfo.samplerate;
    file.oinfo.channels = file.iinfo.channels;

    if (!(file.ofile = open_output(argv[2], &file.oinfo)))
	goto done;

    process_file(&file);
done:
    if (file.ifile) sndfile.close(file.ifile);
    if (file.ofile) sndfile.close(file.ofile);
    if (sndfile.handle) FreeLibrary(sndfile.handle);
}
