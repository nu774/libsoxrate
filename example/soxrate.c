#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sndfile.h>
#ifdef _WIN32
#include <windows.h>
#include "getopt.h"
#else
#include <unistd.h>
#include <sys/time.h>
#endif
#include "libsoxrate.h"

typedef struct option_t {
    char *ofilename;
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

void parse_option(option_t *opts, int *argc, char ***argv)
{
    int c;
    option_t tmp = { 0, 0, -1, -1, 0, 0, 0 };
    while ((c = getopt(*argc, *argv, "o:r:q:p:b:at")) != EOF) {
	if (c == 'o')
	    tmp.ofilename = optarg;
	else if (c == 'r') {
	    if (sscanf(optarg, "%u", &tmp.rate) != 1)
		usage();
	} else if (c == 'q') {
	    if (sscanf(optarg, "%u", &tmp.quality) != 1)
		usage();
	} else if (c == 'p') {
	    if (sscanf(optarg, "%u", &tmp.phase) != 1)
		usage();
	} else if (c == 'b') {
	    if (sscanf(optarg, "%lf", &tmp.band) != 1)
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

void process_file(const file_t *file, const option_t *opts)
{
#define BS 0x4000
    double *ibuf = 0, *obuf = 0, *ip;
    sf_count_t count = 0, total = 0;
    size_t nch = file->iinfo.channels;
    int eof = 0;
    lsx_rate_t *conv = 0;
    double start, last, now;

    if (!(ibuf = malloc(sizeof(double) * BS * file->iinfo.channels))) {
	fputs("ERROR: malloc()\n", stderr);
	goto done;
    }
    if (!(obuf = malloc(sizeof(double) * BS * file->iinfo.channels))) {
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
    start = last = get_time();
    for (;;) {
	size_t ilen, olen;

	if (!eof && !count) {
	    count = sf_readf_double(file->ifile, ibuf, BS);
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
	if (lsx_rate_process_double(conv, ip, obuf, &ilen, &olen) < 0) {
	    fputs("\nERROR: lsx_rate_process()\n", stderr);
	    break;
	}
#if 0
	{
	    FILE *fp = fopen("soxrate.log", "a");
	    fprintf(fp, "in %d\tout %d\n", ilen, olen);
	    fclose(fp);
	}
#endif
	if (!count && !olen) {
	    fputs("...done\n", stderr);
	    break;
	}
	ip += ilen * nch;
	count -= ilen;
	sf_writef_double(file->ofile, obuf, olen);
    }
done:
    if (conv) lsx_rate_close(conv);
    if (ibuf) free(ibuf);
    if (obuf) free(obuf);
#undef BS
}

#ifdef _WIN32
int utf8_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
    option_t opts = { 0 };
    file_t file = { 0 };
    int rc = 2;

    parse_option(&opts, &argc, &argv);

    if (!(file.ifile = open_input(argv[0], &file.iinfo)))
	goto done;

    file.oinfo.format = SF_FORMAT_WAV | SF_FORMAT_FLOAT;
    file.oinfo.samplerate = opts.rate;
    file.oinfo.channels = file.iinfo.channels;

    if (!(file.ofile = open_output(opts.ofilename, &file.oinfo)))
	goto done;

    process_file(&file, &opts);
    rc = 0;
done:
    if (file.ifile) sf_close(file.ifile);
    if (file.ofile) sf_close(file.ofile);
    return rc;
}
