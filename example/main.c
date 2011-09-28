#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include "sndfile.h"
#include "getopt.h"
#include "libsoxrate.h"

typedef struct option_t {
    wchar_t *ofilename;
    unsigned rate;
    int quality;
    int phase;
    double band;
    int alias;
    int threads;
} option_t;

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
    fprintf(stderr,
"usage: soxrate [options] infile\n"
"options:\n"
"-r <n>         samplerate [required]\n"
"-q <n>         quality, 0:mid 1:high 2:very high\n"
"-p <n>         phase response, 0:minimum, 1:intermediate, 2:linear\n"
"-b <n>         bandwidth, 80 - 99.7\n"
"-a             allow aliasing\n"
"-t             use threading\n"
"-o <filename>  output filename\n"
    );
    exit(1);
}

void parse_option(option_t *opts, int *argc, wchar_t ***argv)
{
    int c;
    option_t tmp = { 0, 0, -1, -1, 0, 0, 0 };
    while ((c = getopt(*argc, *argv, L"o:r:q:p:b:at")) != EOF) {
	if (c == 'o')
	    tmp.ofilename = optarg;
	else if (c == 'r') {
	    if (swscanf(optarg, L"%u", &tmp.rate) != 1)
		usage();
	} else if (c == 'q') {
	    if (swscanf(optarg, L"%u", &tmp.quality) != 1)
		usage();
	} else if (c == 'p') {
	    if (swscanf(optarg, L"%u", &tmp.phase) != 1)
		usage();
	} else if (c == 'b') {
	    if (swscanf(optarg, L"%lf", &tmp.band) != 1)
		usage();
	} else if (c == 'a')
	    tmp.alias = 1;
	else if (c == 't')
	    tmp.threads = 1;
    }
    *argc += optind;
    *argv += optind;
    if (!tmp.rate || *argc == 0)
	usage();
    *opts = tmp;
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

void config_rate_converter(lsx_rate_t *conv, const option_t *opts)
{
    if (opts->quality >= 0) {
	if (lsx_rate_config(conv, SOX_RATE_QUALITY, opts->quality) < 0)
	    fputs("WARNING: lsx_rate_config(): quality\n", stderr);
    }
    if (opts->phase >= 0) {
	if (lsx_rate_config(conv, SOX_RATE_PHASE_RESPONSE, opts->phase) < 0)
	    fputs("WARNING: lsx_rate_config(): phase\n", stderr);
    }
    if (opts->band > 0) {
	if (lsx_rate_config(conv, SOX_RATE_PHASE_RESPONSE, opts->band) < 0)
	    fputs("WARNING: lsx_rate_config(): band\n", stderr);
    }
    if (lsx_rate_config(conv, SOX_RATE_ALLOW_ALIASING, opts->alias) < 0)
	fputs("WARNING: lsx_rate_config(): alias\n", stderr);
    if (lsx_rate_config(conv, SOX_RATE_USE_THREADS, opts->threads) < 0)
	fputs("WARNING: lsx_rate_config(): thread\n", stderr);
}

void process_file(const file_t *file, const option_t *opts)
{
#define BS 0x4000
    float *ibuf = 0, *obuf = 0, *ip;
    sf_count_t count = 0, total = 0;
    size_t nch = file->iinfo.channels;
    int eof = 0;
    lsx_rate_t *conv = 0;
    DWORD start, last, now;

    if (!(ibuf = malloc(sizeof(float) * BS * file->iinfo.channels))) {
	fputs("ERROR: malloc()\n", stderr);
	goto done;
    }
    if (!(obuf = malloc(sizeof(float) * BS * file->iinfo.channels))) {
	fputs("ERROR: malloc()\n", stderr);
	goto done;
    }
    ip = ibuf;
    
    if (!(conv = lsx_rate_create(nch, file->iinfo.samplerate, opts->rate))) {
	fputs("ERROR: lsx_rate_create()\n", stderr);
	goto done;
    }
    config_rate_converter(conv, opts);
    if (lsx_rate_start(conv) < 0) {
	fputs("ERROR: lsx_rate_start()\n", stderr);
	goto done;
    }
    start = last = GetTickCount();
    for (;;) {
	size_t ilen, olen;

	if (!eof && !count) {
	    count = sndfile.readf_float(file->ifile, ibuf, BS);
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
	if (lsx_rate_process(conv, ip, obuf, &ilen, &olen) < 0) {
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
    if (conv) lsx_rate_close(conv);
    if (ibuf) free(ibuf);
    if (obuf) free(obuf);
#undef BS
}

int wmain(int argc, wchar_t **argv)
{
    option_t opts;
    file_t file = { 0 };

    parse_option(&opts, &argc, &argv);
    load_libsndfile();

    if (!(file.ifile = open_input(argv[0], &file.iinfo)))
	goto done;

    file.oinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    file.oinfo.samplerate = opts.rate;
    file.oinfo.channels = file.iinfo.channels;

    if (!(file.ofile = open_output(opts.ofilename, &file.oinfo)))
	goto done;

    process_file(&file, &opts);
done:
    if (file.ifile) sndfile.close(file.ifile);
    if (file.ofile) sndfile.close(file.ofile);
    if (sndfile.handle) FreeLibrary(sndfile.handle);
}
