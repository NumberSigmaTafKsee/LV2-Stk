#include "SoundPipe.h"
#include "sndfile.h"
#include "dr_wav.h"
#include "dr_wav.c"

#define MAX_SECTION 50
#define MAX_NAME 50

#ifndef dB
/* if below -100dB, set to -100dB to prevent taking log of zero */
#define dB(x) 20.0 * ((x) > 0.00001 ? log10(x) : log10(0.00001))
#endif

#ifndef dB2lin
#define dB2lin(x)           pow( 10.0, (x) / 20.0 )
#endif

#define max(a,b) ((a < b) ? b : a)
#define min(a,b) ((a < b) ? a : b)

#define ROOT2 (1.4142135623730950488)

#define WUTR_SOUND_DECAY 0.95
#define WUTR_SYSTEM_DECAY 0.996
#define WUTR_GAIN 1.0
#define WUTR_NUM_SOURCES 10.0
#define WUTR_CENTER_FREQ0 450.0
#define WUTR_CENTER_FREQ1 600.0
#define WUTR_CENTER_FREQ2 750.0
#define WUTR_RESON 0.9985
#define WUTR_FREQ_SWEEP 1.0001
#define MAX_SHAKE 2000

#define PFRAC1(x)   ((SPFLOAT)((x) & ftp1->lomask) * ftp1->lodiv)

#define MPIDSR -M_PI/sp->sr

#define log001 (-(SPFLOAT)6.9078)    /* log(.001) */
#define tpd360  0.0174532925199433

#define OPENLPC_FRAMESIZE_1_8	    250
#define OPENLPC_FRAMESIZE_1_4	    320
#define OPENLPC_ENCODED_FRAME_SIZE  7

#define PREEMPH

#define bcopy(a, b, n)	  memmove(b, a, n)

#define LPC_FILTORDER		10
#define FS		my_fs /* Sampling rate */
#define MAXWINDOW	1000	/* Max analysis window length */

#define FC		200.0	/* Pitch analyzer filter cutoff */
#define DOWN		5	/* Decimation for pitch analyzer */
#define MINPIT		40.0	/* Minimum pitch (observed: 74) */
#define MAXPIT		320.0	/* Maximum pitch (observed: 250) */

#define MINPER		(int)(FS/(DOWN*MAXPIT)+.5)	/* Minimum period  */
#define MAXPER		(int)(FS/(DOWN*MINPIT)+.5)	/* Maximum period  */

#define REAL_MINPER	 (DOWN*MINPER) /* converted to samples units */

#define WSCALE		1.5863	/* Energy loss due to windowing */

#define BITS_FOR_LPC 38

#define ARCSIN_Q /* provides better quantization of first two k[] at low bitrates */

#define FP_BITS(fp) (*(int *)&(fp))
#define FIST_FLOAT_MAGIC_S (float)(7.0f * 2097152.0f)


#define PLUKMIN 64

const double SQRT2 = sqrt(2.0);

/* Strip whitespace chars off end of given string, in place. Return s. */
static char* rstrip(char* s)
{
    char* p = s + strlen(s);
    while (p > s && isspace((unsigned char)(*--p)))
        *p = '\0';
    return s;
}

/* Return pointer to first non-whitespace char in given string. */
static char* lskip(const char* s)
{
    while (*s && isspace((unsigned char)(*s)))
        s++;
    return (char*)s;
}

/* Return pointer to first char c or ';' comment in given string, or pointer to
   null at end of string if neither found. ';' must be prefixed by a whitespace
   character to register as a comment. */
static char* find_char_or_comment(const char* s, char c)
{
    int was_whitespace = 0;
    while (*s && *s != c && !(was_whitespace && *s == ';')) {
        was_whitespace = isspace((unsigned char)(*s));
        s++;
    }
    return (char*)s;
}

/* Version of strncpy that ensures dest (size bytes) is null-terminated. */
static char* strncpy0(char* dest, const char* src, size_t size)
{
    strncpy(dest, src, size);
    dest[size - 1] = '\0';
    return dest;
}

/* See documentation in header file. */
int ini_parse_file(FILE* file,
                   int (*handler)(void*, const char*, const char*,
                                  const char*),
                   void* user)
{
    /* Uses a fair bit of stack (use heap instead if you need to) */
#if INI_USE_STACK
    char line[INI_MAX_LINE];
#else
    char* line;
#endif
    char section[MAX_SECTION] = "";
    char prev_name[MAX_NAME] = "";

    char* start;
    char* end;
    char* name;
    char* value;
    int lineno = 0;
    int error = 0;

#if !INI_USE_STACK
    line = (char*)malloc(INI_MAX_LINE);
    if (!line) {
        return -2;
    }
#endif

    /* Scan through file line by line */
    while (fgets(line, INI_MAX_LINE, file) != NULL) {
        lineno++;

        start = line;
#if INI_ALLOW_BOM
        if (lineno == 1 && (unsigned char)start[0] == 0xEF &&
                           (unsigned char)start[1] == 0xBB &&
                           (unsigned char)start[2] == 0xBF) {
            start += 3;
        }
#endif
        start = lskip(rstrip(start));

        if (*start == ';' || *start == '#') {
            /* Per Python ConfigParser, allow '#' comments at start of line */
        }
#if INI_ALLOW_MULTILINE
        else if (*prev_name && *start && start > line) {
            /* Non-black line with leading whitespace, treat as continuation
               of previous name's value (as per Python ConfigParser). */
            if (!handler(user, section, prev_name, start) && !error)
                error = lineno;
        }
#endif
        else if (*start == '[') {
            /* A "[section]" line */
            end = find_char_or_comment(start + 1, ']');
            if (*end == ']') {
                *end = '\0';
                strncpy0(section, start + 1, sizeof(section));
                *prev_name = '\0';
            }
            else if (!error) {
                /* No ']' found on section line */
                error = lineno;
            }
        }
        else if (*start && *start != ';') {
            /* Not a comment, must be a name[=:]value pair */
            end = find_char_or_comment(start, '=');
            if (*end != '=') {
                end = find_char_or_comment(start, ':');
            }
            if (*end == '=' || *end == ':') {
                *end = '\0';
                name = rstrip(start);
                value = lskip(end + 1);
                end = find_char_or_comment(value, '\0');
                if (*end == ';')
                    *end = '\0';
                rstrip(value);

                /* Valid name[=:]value pair found, call handler */
                strncpy0(prev_name, name, sizeof(prev_name));
                if (!handler(user, section, name, value) && !error)
                    error = lineno;
            }
            else if (!error) {
                /* No '=' or ':' found on name[=:]value line */
                error = lineno;
            }
        }

#if INI_STOP_ON_FIRST_ERROR
        if (error)
            break;
#endif
    }

#if !INI_USE_STACK
    free(line);
#endif

    return error;
}

/* See documentation in header file. */
int ini_parse(const char* filename,
              int (*handler)(void*, const char*, const char*, const char*),
              void* user)
{
    FILE* file;
    int error;

    file = fopen(filename, "r");
    if (!file)
        return -1;
    error = ini_parse_file(file, handler, user);
    fclose(file);
    return error;
}




enum { CLEAR, ATTACK, DECAY, SUSTAIN, RELEASE, KEY_ON, KEY_OFF };

int sp_adsr_create(sp_adsr **p)
{
    *p = malloc(sizeof(sp_adsr));
    return SP_OK;
}

int sp_adsr_destroy(sp_adsr **p)
{
    free(*p);
    return SP_OK;
}

int sp_adsr_init(sp_data *sp, sp_adsr *p)
{
    p->atk = 0.1;
    p->dec = 0.1;
    p->sus = 0.5;
    p->rel = 0.3;
    p->timer = 0;
    p->a = 0;
    p->b = 0;
    p->y = 0;
    p->x = 0;
    p->prev = 0;
    p->atk_time = p->atk * sp->sr;
    p->mode = CLEAR;
    return SP_OK;
}

static SPFLOAT tau2pole(sp_data *sp, sp_adsr *p, SPFLOAT tau)
{
    return exp(-1.0 / (tau * sp->sr));
}

static SPFLOAT adsr_filter(sp_data *sp, sp_adsr *p)
{
    p->y = p->b * p->x  + p->a * p->y;
    return p->y;
}


static void addCheckButton (void* ui_interface, const char* label, FAUSTFLOAT* zone)
{
    sp_phaser *p = ui_interface;
    p->args[p->argpos] = zone;
    p->argpos++;
}
static void addVerticalSlider(void* ui_interface, const char* label, FAUSTFLOAT* zone, FAUSTFLOAT init, FAUSTFLOAT min, FAUSTFLOAT max, FAUSTFLOAT step)
{
    sp_autowah *p = ui_interface;
    p->args[p->argpos] = zone;
    p->argpos++;
}
static void addHorizontalSlider(void* ui_interface, const char* label, FAUSTFLOAT* zone, FAUSTFLOAT init, FAUSTFLOAT min, FAUSTFLOAT max, FAUSTFLOAT step)
{
    sp_compressor *p = ui_interface;
    p->args[p->argpos] = zone;
    p->argpos++;
}

int sp_adsr_compute(sp_data *sp, sp_adsr *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT pole;
    if(p->prev < *in && p->mode != DECAY) {
        p->mode = ATTACK;
        p->timer = 0;
        /* quick fix: uncomment if broken */
        /* pole = tau2pole(sp, p, p->atk * 0.75); */
        /* p->atk_time = p->atk * sp->sr * 1.5; */
        pole = tau2pole(sp, p, p->atk * 0.6);
        p->atk_time = p->atk * sp->sr;
        p->a = pole;
        p->b = 1 - pole;
    } else if(p->prev > *in) {
        p->mode = RELEASE;
        pole = tau2pole(sp, p, p->rel);
        p->a = pole;
        p->b = 1 - pole;
    }

    p->x = *in;
    p->prev = *in;

    switch(p->mode) {
        case CLEAR:
            *out = 0;
            break;
        case ATTACK:
            p->timer++;
            *out = adsr_filter(sp, p);
            /* quick fix: uncomment if broken */
            /* if(p->timer > p->atk_time) { */
            if(*out > 0.99) {
                p->mode = DECAY;
                pole = tau2pole(sp, p, p->dec);
                p->a = pole;
                p->b = 1 - pole;
            }
            break;
        case DECAY:
        case RELEASE:
            p->x *= p->sus;
            *out = adsr_filter(sp, p);
        default:
            break;        
    }

    return SP_OK;
}


int sp_allpass_create(sp_allpass **p)
{
    *p = malloc(sizeof(sp_allpass));
    return SP_OK;
}

int sp_allpass_destroy(sp_allpass **p)
{
    sp_allpass *pp = *p;
    sp_auxdata_free(&pp->aux);
    free(*p);
    return SP_OK;
}

int sp_allpass_init(sp_data *sp, sp_allpass *p, SPFLOAT looptime)
{
    p->revtime = 3.5;
    p->looptime = looptime;
    p->bufsize = 0.5 + looptime * sp->sr;
    sp_auxdata_alloc(&p->aux, p->bufsize * sizeof(SPFLOAT));
    p->prvt = 0.0;
    p->coef = 0.0;
    p->bufpos = 0;
    return SP_OK;
}

int sp_allpass_compute(sp_data *sp, sp_allpass *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT y, z;
    SPFLOAT coef = p->coef;
    SPFLOAT *buf = (SPFLOAT *)p->aux.ptr;
    if(p->prvt != p->revtime) {
        p->prvt = p->revtime;
        coef = p->coef = exp(-6.9078 * p->looptime / p->prvt);
    }
    y = buf[p->bufpos];
    z = coef * y + *in; 
    buf[p->bufpos] = z;
    *out = y - coef * z;

    p->bufpos++;
    p->bufpos %= p->bufsize; 
    return SP_OK;
}


int sp_atone_create(sp_atone **p)
{
    *p = malloc(sizeof(sp_atone));
    return SP_OK;
}

int sp_atone_destroy(sp_atone **p)
{
    free(*p);
    return SP_OK;
}

int sp_atone_init(sp_data *sp, sp_atone *p)
{
    p->hp = 1000;
    SPFLOAT b;
    p->tpidsr = (2.0 * M_PI) / sp->sr * 1.0;
    p->prvhp = (SPFLOAT)p->hp;
    b = 2.0 - cos((SPFLOAT)(p->prvhp * p->tpidsr));
    p->c2 = b - sqrt(b * b - 1.0);
    p->c1 = 1.0 - p->c2;
    p->yt1 = 0.0;
    return SP_OK;
}

int sp_atone_compute(sp_data *sp, sp_atone *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT c2 = p->c2, yt1 = p->yt1;
    SPFLOAT x;

    if (p->hp != p->prvhp) {
      SPFLOAT b;
      p->prvhp = p->hp;
      b = 2.0 - cos((SPFLOAT)(p->hp * p->tpidsr));
      p->c2 = c2 = b - sqrt(b * b - 1.0);
    }

    x = yt1 = c2 * (yt1 + *in);
    *out = x;
    yt1 -= *in;
    p->yt1 = yt1;
    return SP_OK;
}



static float faustpower2_f(float value) {
	return (value * value);
}

typedef struct {
	
	float fRec0[3];
	float fRec3[2];
	float fRec2[2];
	float fRec1[2];
	float fRec4[2];
	float fRec5[2];
	FAUSTFLOAT fVslider0;
	int fSamplingFreq;
	int iConst0;
	float fConst1;
	float fConst2;
	float fConst3;
	float fConst4;
	float fConst5;
	float fConst6;
	FAUSTFLOAT fVslider1;
	FAUSTFLOAT fVslider2;
	
} autowah;

autowah* newautowah() { 
	autowah* dsp = (autowah*)malloc(sizeof(autowah));
	return dsp;
}

void deleteautowah(autowah* dsp) { 
	free(dsp);
}

void instanceInitautowah(autowah* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->fVslider0 = (FAUSTFLOAT)0.;
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fConst1 = (1413.72f / (float)dsp->iConst0);
	dsp->fConst2 = exp((0.f - (100.f / (float)dsp->iConst0)));
	dsp->fConst3 = (1.f - dsp->fConst2);
	dsp->fConst4 = exp((0.f - (10.f / (float)dsp->iConst0)));
	dsp->fConst5 = (1.f - dsp->fConst4);
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->fRec3[i0] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec2[i1] = 0.f;
			
		}
		
	}
	dsp->fConst6 = (2827.43f / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fRec1[i2] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i3;
		for (i3 = 0; (i3 < 2); i3 = (i3 + 1)) {
			dsp->fRec4[i3] = 0.f;
			
		}
		
	}
	dsp->fVslider1 = (FAUSTFLOAT)100.;
	dsp->fVslider2 = (FAUSTFLOAT)0.1;
	/* C99 loop */
	{
		int i4;
		for (i4 = 0; (i4 < 2); i4 = (i4 + 1)) {
			dsp->fRec5[i4] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i5;
		for (i5 = 0; (i5 < 3); i5 = (i5 + 1)) {
			dsp->fRec0[i5] = 0.f;
			
		}
		
	}
	
}


void initautowah(autowah* dsp, int samplingFreq) {
	instanceInitautowah(dsp, samplingFreq);
}

void buildUserInterfaceautowah(autowah* dsp, UIGlue* interface) {
	interface->addVerticalSlider(interface->uiInterface, "level", &dsp->fVslider2, 0.1f, 0.f, 1.f, 0.01f);
	interface->addVerticalSlider(interface->uiInterface, "wah", &dsp->fVslider0, 0.f, 0.f, 1.f, 0.01f);
	interface->addVerticalSlider(interface->uiInterface, "wet_dry", &dsp->fVslider1, 100.f, 0.f, 100.f, 1.f);
}

void computeautowah(autowah* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* input0 = inputs[0];
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fVslider0;
	float fSlow1 = (float)dsp->fVslider1;
	float fSlow2 = (0.01f * (fSlow1 * (float)dsp->fVslider2));
	float fSlow3 = ((1.f - (0.01f * fSlow1)) + (1.f - fSlow0));
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			float fTemp0 = (float)input0[i];
			float fTemp1 = fabs(fTemp0);
			dsp->fRec3[0] = max(fTemp1, ((dsp->fConst4 * dsp->fRec3[1]) + (dsp->fConst5 * fTemp1)));
			dsp->fRec2[0] = ((dsp->fConst2 * dsp->fRec2[1]) + (dsp->fConst3 * dsp->fRec3[0]));
			float fTemp2 = min(1.f, dsp->fRec2[0]);
			float fTemp3 = pow(2.f, (2.3f * fTemp2));
			float fTemp4 = (1.f - (dsp->fConst1 * (fTemp3 / pow(2.f, (1.f + (2.f * (1.f - fTemp2)))))));
			dsp->fRec1[0] = ((0.999f * dsp->fRec1[1]) + (0.001f * (0.f - (2.f * (fTemp4 * cos((dsp->fConst6 * fTemp3)))))));
			dsp->fRec4[0] = ((0.999f * dsp->fRec4[1]) + (0.001f * faustpower2_f(fTemp4)));
			dsp->fRec5[0] = ((0.999f * dsp->fRec5[1]) + (0.0001f * pow(4.f, fTemp2)));
			dsp->fRec0[0] = (0.f - (((dsp->fRec1[0] * dsp->fRec0[1]) + (dsp->fRec4[0] * dsp->fRec0[2])) - (fSlow2 * (dsp->fRec5[0] * fTemp0))));
			output0[i] = (FAUSTFLOAT)((fSlow0 * (dsp->fRec0[0] - dsp->fRec0[1])) + (fSlow3 * fTemp0));
			dsp->fRec3[1] = dsp->fRec3[0];
			dsp->fRec2[1] = dsp->fRec2[0];
			dsp->fRec1[1] = dsp->fRec1[0];
			dsp->fRec4[1] = dsp->fRec4[0];
			dsp->fRec5[1] = dsp->fRec5[0];
			dsp->fRec0[2] = dsp->fRec0[1];
			dsp->fRec0[1] = dsp->fRec0[0];
		}
		
	}
	
}


int sp_autowah_create(sp_autowah **p)
{
    *p = malloc(sizeof(sp_autowah));
    return SP_OK;
}

int sp_autowah_destroy(sp_autowah **p)
{
    sp_autowah *pp = *p;
    autowah *dsp = pp->faust;
    deleteautowah (dsp);
    free(*p);
    return SP_OK;
}




int sp_autowah_init(sp_data *sp, sp_autowah *p)
{
    autowah *dsp = newautowah(); 
    UIGlue UI;
    p->argpos = 0;
    UI.addVerticalSlider= addVerticalSlider;
    UI.uiInterface = p;
    buildUserInterfaceautowah(dsp, &UI);
    initautowah(dsp, sp->sr);
    
    p->level = p->args[0]; 
    p->wah = p->args[1]; 
    p->mix = p->args[2];

    p->faust = dsp;
    return SP_OK;
}

int sp_autowah_compute(sp_data *sp, sp_autowah *p, SPFLOAT *in, SPFLOAT *out) 
{
    autowah *dsp = p->faust;
    SPFLOAT *faust_out[] = {out};
    SPFLOAT *faust_in[] = {in};
    computeautowah(dsp, 1, faust_in, faust_out);
    return SP_OK;
}

int sp_bal_create(sp_bal **p)
{
    *p = malloc(sizeof(sp_bal));
    return SP_OK;
}

int sp_bal_destroy(sp_bal **p)
{
    free(*p);
    return SP_OK;
}

int sp_bal_init(sp_data *sp, sp_bal *p)
{

    SPFLOAT b;
    p->ihp = 10;
    b = 2.0 - cos((SPFLOAT)(p->ihp * (2.0 * M_PI / sp->sr)));
    p->c2 = b - sqrt(b*b - 1.0);
    p->c1 = 1.0 - p->c2;
    p->prvq = p->prvr = p->prva = 0.0;

    return SP_OK;
}

int sp_bal_compute(sp_data *sp, sp_bal *p, SPFLOAT *sig, SPFLOAT *comp, SPFLOAT *out)
{
    SPFLOAT q, r, a, diff;
    SPFLOAT c1 = p->c1, c2 = p->c2;

    q = p->prvq;
    r = p->prvr;
    SPFLOAT as = *sig;
    SPFLOAT cs = *comp;

    q = c1 * as * as + c2 * q;
    r = c1 * cs * cs + c2 * r;

    p->prvq = q;
    p->prvr = r;

    if (q != 0.0) {
        a = sqrt(r/q);
    } else {
        a = sqrt(r);
    }

    if((diff = a - p->prva) != 0.0) {
        *out = *sig * p->prva;
        p->prva = a;
    } else {
        *out = *sig * a;
    }

    return SP_OK;
}


int sp_bar_create(sp_bar **p)
{
    *p = malloc(sizeof(sp_bar));
    return SP_OK;
}

int sp_bar_destroy(sp_bar **p)
{
    sp_bar *pp = *p;
    sp_auxdata_free(&pp->w_aux);
    free(*p);
    return SP_OK;
}

int sp_bar_init(sp_data *sp, sp_bar *p, SPFLOAT iK, SPFLOAT ib)
{
    p->bcL = 1;
    p->bcR = 1;
    p->iK = iK;
    p->ib = ib;
    p->scan = 0.23;
    p->T30 = 3;
    p->pos = 0.2;
    p->vel = 500;
    p->wid = 0.05;

    SPFLOAT K = p->iK;       /* ~=3.0  stiffness parameter, dimensionless */
    SPFLOAT T30 = p->T30;   /* ~=5.0; 30 db decay time (s) */
    SPFLOAT b = p->ib;       /* ~=0.001 high-frequency loss parameter
                               (keep small) */

    /* derived parameters */
    SPFLOAT dt = 1.0 / sp->sr;
    SPFLOAT sig = (2.0 * sp->sr) * (pow(10.0, 3.0 * dt / T30) - 1.0);
    SPFLOAT dxmin = sqrt(dt * (b+hypot(b, K+K) ));
    int N = (int) (1.0/dxmin);
    SPFLOAT dx = 1.0/N;

    /* scheme coefficients */
    p->s0 = (2.0-6.0*K*K*dt*dt/(dx*dx*dx*dx)-2.0*b*dt/(dx*dx))/(1.0+sig*dt*0.5);
    p->s1 = (4.0*K*K*dt*dt/(dx*dx*dx*dx)+b*dt/(dx*dx))/(1.0+sig*dt*0.5);
    p->s2 = -K*K*dt*dt/((dx*dx*dx*dx)*(1.0+sig*dt*0.5));
    p->t0 = (-1.0+2.0*b*dt/(dx*dx)+sig*dt*0.5)/(1.0+sig*dt*0.5);
    p->t1 = (-b*dt)/(dx*dx*(1.0+sig*dt*0.5));

    sp_auxdata_alloc(&p->w_aux, (size_t) 3 * ((N + 5) * sizeof(SPFLOAT)));
    p->w = (SPFLOAT *) p->w_aux.ptr;
    p->w1 = &(p->w[N + 5]);
    p->w2 = &(p->w1[N + 5]);
    p->step = p->first = 0;
    p->N = N;
    p->first = 0;
    return SP_OK;
}

int sp_bar_compute(sp_data *sp, sp_bar *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT xofreq = 2 * M_PI * (p->scan)/sp->sr;
    SPFLOAT xo, xofrac;
    int xoint;
    int step = p->step;
    int first = p->first;
    int N = p->N, rr;
    SPFLOAT *w = p->w, *w1 = p->w1, *w2 = p->w2;
    SPFLOAT s0 = p->s0, s1 = p->s1, s2 = p->s2, t0 = p->t0, t1 = p->t1;
    int bcL = (int)lrintf((SPFLOAT)p->bcL);
    int bcR = (int)lrintf((SPFLOAT)p->bcR);
    SPFLOAT SINNW = sin(xofreq*step);
    SPFLOAT COSNW = cos(xofreq*step);
    SPFLOAT SIN1W = sin(xofreq);
    SPFLOAT COS1W = cos(xofreq);

    if(*in) {
        p->first = 0;
        SPFLOAT K = p->iK;
        SPFLOAT T30 = p->T30;
        SPFLOAT b = p->ib;

        SPFLOAT dt = 1.0 / sp->sr;
        SPFLOAT sig = (2.0 * sp->sr) * (pow(10.0, 3.0 * dt / T30) - 1.0);
        SPFLOAT dxmin = sqrt(dt * (b+hypot(b, K+K) ));
        int N = (int) (1.0/dxmin);
        SPFLOAT dx = 1.0/N;

        p->s0 = (2.0-6.0*K*K*dt*dt/(dx*dx*dx*dx)-2.0*b*dt/(dx*dx))/(1.0+sig*dt*0.5);
        p->s1 = (4.0*K*K*dt*dt/(dx*dx*dx*dx)+b*dt/(dx*dx))/(1.0+sig*dt*0.5);
        p->s2 = -K*K*dt*dt/((dx*dx*dx*dx)*(1.0+sig*dt*0.5));
        p->t0 = (-1.0+2.0*b*dt/(dx*dx)+sig*dt*0.5)/(1.0+sig*dt*0.5);
        p->t1 = (-b*dt)/(dx*dx*(1.0+sig*dt*0.5));

        s0 = p->s0; 
        s1 = p->s1; 
        s2 = p->s2; 
        t0 = p->t0; 
        t1 = p->t1;
    }

    if ((bcL|bcR)&(~3) && (bcL|bcR)!=0) {
        fprintf(stderr,
                "sp_bar: Ends must be clamped(1), pivoting(2), or free(3)\n");
        return SP_NOT_OK;
    }

    if (bcL == 3) {
        w1[1] = 2.0*w1[2]-w1[3];
        w1[0] = 3.0*w1[1]-3.0*w1[2]+w1[3];
    }
    else if (bcL == 1) {
        w1[2] = 0.0;
        w1[3] = 0.0;
    }
    else if (bcL == 2) {
        w1[2] = 0.0;
        w1[1] = -w1[3];
    }

    if (bcR == 3) {
        w1[N+3] = 2.0*w1[N+2]-w1[N+1];
        w1[N+4] = 3.0*w1[N+3]-3.0*w1[N+2]+w1[N+1];
    } else if (bcR == 1) {
        w1[N+1] = 0.0;
        w1[N+2] = 0.0;
    } else if (bcR == 2) {
        w1[N+2] = 0.0;
        w1[N+3] = -w1[N+1];
    }

    /* Iterate model */
    for (rr = 0; rr < N+1; rr++) {
        w[rr+2] = s0*w1[rr+2] + s1*(w1[rr+3]+w1[rr+1]) + s2*(w1[rr+4]+w1[rr]) +
                  t0*w2[rr+2] + t1*(w2[rr+3]+w2[rr+1]);
    }

    /*  strike inputs */

    if (first == 0) {
        p->first = first = 1;
        for (rr = 0; rr < N; rr++) {
            if (fabs(rr/(SPFLOAT)N - p->pos) <= p->wid) {
                w[rr+2] += (1.0/sp->sr)*(p->vel)*0.5*
                    (1.0+cos(M_PI*fabs(rr/(SPFLOAT)N-(p->pos))/(p->wid)));
            }
        }
    }
    {
        SPFLOAT xx = SINNW*COS1W + COSNW*SIN1W;
        SPFLOAT yy = COSNW*COS1W - SINNW*SIN1W;

        SINNW = xx;
        COSNW = yy;
    }
    xo = 0.5 + 0.5*SINNW;
    xoint = (int) (xo*N) + 2;
    xofrac = xo*N - (int)(xo*N);

    *out = ((1.0-xofrac)*w[xoint] + xofrac*w[xoint+1]);
    step++;
    {
        SPFLOAT *ww = w2;

        w2 = w1;
        w1 = w;
        w = ww;
    }
    p->step = step;
    p->w = w;
    p->w1 = w1;
    p->w2 = w2;
    return SP_OK;
}


int sp_create(sp_data **spp)
{
    *spp = (sp_data *) malloc(sizeof(sp_data));
    sp_data *sp = *spp;
    sprintf(sp->filename, "test.wav");
    sp->nchan = 1;
    SPFLOAT *out = malloc(sizeof(SPFLOAT) * sp->nchan);
    *out = 0;
    sp->out = out;
    sp->sr = 44100;
    sp->len = 5 * sp->sr;
    sp->pos = 0;
    sp->rand = 0;
    return 0;
}

int sp_createn(sp_data **spp, int nchan)
{
    *spp = (sp_data *) malloc(sizeof(sp_data));
    sp_data *sp = *spp;
    sprintf(sp->filename, "test.wav");
    sp->nchan = nchan;
    SPFLOAT *out = malloc(sizeof(SPFLOAT) * sp->nchan);
    *out = 0;
    sp->out = out;
    sp->sr = 44100;
    sp->len = 5 * sp->sr;
    sp->pos = 0;
    sp->rand = 0;
    return 0;
}

int sp_destroy(sp_data **spp)
{
    sp_data *sp = *spp;
    free(sp->out);
    free(*spp);
    return 0;
}

#ifndef NO_LIBSNDFILE

int sp_process(sp_data *sp, void *ud, void (*callback)(sp_data *, void *))
{
    SNDFILE *sf[sp->nchan];
    char tmp[140];
    SF_INFO info;
    memset(&info, 0, sizeof(SF_INFO));
    SPFLOAT buf[sp->nchan][SP_BUFSIZE];
    info.samplerate = sp->sr;
    info.channels = 1;
    info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;
    int numsamps, i, chan;
    if(sp->nchan == 1) {
        sf[0] = sf_open(sp->filename, SFM_WRITE, &info);
    } else {
        for(chan = 0; chan < sp->nchan; chan++) {
            sprintf(tmp, "%02d_%s", chan, sp->filename);
            sf[chan] = sf_open(tmp, SFM_WRITE, &info);
        }
    }

    while(sp->len > 0){
        if(sp->len < SP_BUFSIZE) {
            numsamps = sp->len;
        }else{
            numsamps = SP_BUFSIZE;
        }
        for(i = 0; i < numsamps; i++){
            callback(sp, ud);
            for(chan = 0; chan < sp->nchan; chan++) {
                buf[chan][i] = sp->out[chan];
            }
            sp->pos++;
        }
        for(chan = 0; chan < sp->nchan; chan++) {
#ifdef USE_DOUBLE
            sf_write_double(sf[chan], buf[chan], numsamps);
#else
            sf_write_float(sf[chan], buf[chan], numsamps);
#endif
        }
        sp->len -= numsamps;
    }
    for(i = 0; i < sp->nchan; i++) {
        sf_close(sf[i]);
    }
    return 0;
}

#endif

int sp_process_raw(sp_data *sp, void *ud, void (*callback)(sp_data *, void *))
{
    int chan;
    if(sp->len == 0) {
        while(1) {
            callback(sp, ud);
            for (chan = 0; chan < sp->nchan; chan++) {
                fwrite(&sp->out[chan], sizeof(SPFLOAT), 1, stdout);
            }
            sp->len--;
        }
    } else {
        while(sp->len > 0) {
            callback(sp, ud);
            for (chan = 0; chan < sp->nchan; chan++) {
                fwrite(&sp->out[chan], sizeof(SPFLOAT), 1, stdout);
            }
            sp->len--;
            sp->pos++;
        }
    }
    return SP_OK;
}

int sp_process_plot(sp_data *sp, void *ud, void (*callback)(sp_data *, void *))
{
    int chan;
    fprintf(stdout, "sp_out =  [ ... \n");
    while(sp->len > 0) {
        callback(sp, ud);
        for (chan = 0; chan < sp->nchan; chan++) {
            /* fwrite(&sp->out[chan], sizeof(SPFLOAT), 1, stdout); */
            fprintf(stdout, "%g ", sp->out[chan]);
        }
        fprintf(stdout, "; ...\n");
        sp->len--;
        sp->pos++;
    }
    fprintf(stdout, "];\n");

    fprintf(stdout, "plot(sp_out);\n");
    fprintf(stdout, "title('Plot generated by Soundpipe');\n");
    fprintf(stdout, "xlabel('Time (samples)');\n");
    fprintf(stdout, "ylabel('Amplitude');\n");
    return SP_OK;
}

int sp_auxdata_alloc(sp_auxdata *aux, size_t size)
{
    aux->ptr = malloc(size);
    aux->size = size;
    memset(aux->ptr, 0, size);
    return SP_OK;
}

int sp_auxdata_free(sp_auxdata *aux)
{
    free(aux->ptr);
    return SP_OK;
}


SPFLOAT sp_midi2cps(SPFLOAT nn)
{
    return pow(2, (nn - 69.0) / 12.0) * 440.0;
}

int sp_set(sp_param *p, SPFLOAT val) {
    p->state = 1;
    p->val = val;
    return SP_OK;
}

int sp_out(sp_data *sp, uint32_t chan, SPFLOAT val)
{
    if(chan > sp->nchan - 1) {
        fprintf(stderr, "sp_out: Invalid channel\n");
        return SP_NOT_OK;
    }
    sp->out[chan] = val;
    return SP_OK;
}

uint32_t sp_rand(sp_data *sp)
{
    uint32_t val = (1103515245 * sp->rand + 12345) % SP_RANDMAX;
    sp->rand = val;
    return val;
}

void sp_srand(sp_data *sp, uint32_t val)
{
    sp->rand = val;
}


int sp_biquad_create(sp_biquad **p)
{
    *p = malloc(sizeof(sp_biquad));
    return SP_OK;
}

int sp_biquad_destroy(sp_biquad **p)
{
    free(*p);
    return SP_OK;
}

int sp_biquad_init(sp_data *sp, sp_biquad *p)
{
    p->tpidsr = 2.0*M_PI / sp->sr;
    p->sr = sp->sr;

    p->cutoff = 500;
    p->res = 0.7;
    p->reinit = 0.0;

    SPFLOAT fcon = p->cutoff * p->tpidsr;
    SPFLOAT alpha = 1-2*p->res*cos(fcon)*cos(fcon)+p->res*p->res*cos(2*fcon);
    SPFLOAT beta = p->res*p->res*sin(2*fcon)-2*p->res*cos(fcon)*sin(fcon);
    SPFLOAT gamma = 1+cos(fcon);
    SPFLOAT m1 = alpha*gamma+beta*sin(fcon);
    SPFLOAT m2 = alpha*gamma-beta*sin(fcon);
    SPFLOAT den = sqrt(m1*m1+m2*m2);

    p->b0 = 1.5*(alpha*alpha+beta*beta)/den;
    p->b1 = p->b0;
    p->b2 = 0.0;
    p->a0 = 1.0;
    p->a1 = -2.0*p->res*cos(fcon);
    p->a2 = p->res*p->res;


   if(p->reinit == 0.0){
      p->xnm1 = p->xnm2 = p->ynm1 = p->ynm2 = 0.0;
   }
   return SP_OK; 
}

int sp_biquad_compute(sp_data *sp, sp_biquad *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT xn, yn;
    SPFLOAT a0 = p->a0, a1 = p->a1, a2 = p->a2;
    SPFLOAT b0 = p->b0, b1 = p->b1, b2 = p->b2;
    SPFLOAT xnm1 = p->xnm1, xnm2 = p->xnm2, ynm1 = p->ynm1, ynm2 = p->ynm2;

    xn = *in;
    yn = ( b0 * xn + b1 * xnm1 + b2 * xnm2 -
             a1 * ynm1 - a2 * ynm2) / a0;
    xnm2 = xnm1;
    xnm1 = xn;
    ynm2 = ynm1;
    ynm1 = yn;
    *out = yn;
    
    p->xnm1 = xnm1; p->xnm2 = xnm2; p->ynm1 = ynm1; p->ynm2 = ynm2;
    return SP_OK;
}


int sp_biscale_create(sp_biscale **p)
{
    *p = malloc(sizeof(sp_biscale));
    return SP_OK;
}

int sp_biscale_destroy(sp_biscale **p)
{
    free(*p);
    return SP_OK;
}

int sp_biscale_init(sp_data *sp, sp_biscale *p)
{
    p->min = 0;
    p->max = 1;
    return SP_OK;
}

int sp_biscale_compute(sp_data *sp, sp_biscale *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = p->min + (*in + 1.0) / 2.0 * (p->max - p->min);
    return SP_OK;
}


int sp_bitcrush_create(sp_bitcrush **p)
{
    *p = malloc(sizeof(sp_bitcrush));
    return SP_OK;
}

int sp_bitcrush_destroy(sp_bitcrush **p)
{
    sp_bitcrush *pp = *p;
    sp_fold_destroy(&pp->fold);
    free(*p);
    return SP_OK;
}

int sp_bitcrush_init(sp_data *sp, sp_bitcrush *p)
{
    p->bitdepth = 8;
    p->srate = 10000;
    sp_fold_create(&p->fold);
    sp_fold_init(sp, p->fold);
    return SP_OK;
}

int sp_bitcrush_compute(sp_data *sp, sp_bitcrush *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT bits = pow(2, floor(p->bitdepth));
    SPFLOAT foldamt = sp->sr / p->srate;
    SPFLOAT sig;
    *out = *in * 65536.0;
    *out += 32768;
    *out *= (bits / 65536.0);
    *out = floor(*out);
    *out = *out * (65536.0 / bits) - 32768;
    sig = *out;
    p->fold->incr = foldamt;
    sp_fold_compute(sp, p->fold, &sig, out);
    *out /= 65536.0;
    return SP_OK;
}

typedef struct {

	float fRec0[2];
	float fVec0[2];
	float fVec1[2];
	int fSamplingFreq;
	int iConst0;
	FAUSTFLOAT fHslider0;
	FAUSTFLOAT fHslider1;
	float fConst1;
	float fConst2;
} blsaw;

blsaw* newblsaw() {
	blsaw* dsp = (blsaw*)malloc(sizeof(blsaw));
	return dsp;
}

void deleteblsaw(blsaw* dsp) {
	free(dsp);
}

void instanceInitblsaw(blsaw* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fHslider0 = (FAUSTFLOAT)1.;
	dsp->fHslider1 = (FAUSTFLOAT)440.;
	dsp->fConst1 = (float)dsp->iConst0;
	dsp->fConst2 = (2.f / dsp->fConst1);
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->fRec0[i0] = 0.f;

		}

	}
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fVec0[i1] = 0.f;

		}

	}
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fVec1[i2] = 0.f;

		}

	}

}

void initblsaw(blsaw* dsp, int samplingFreq) {
	instanceInitblsaw(dsp, samplingFreq);
}

void buildUserInterfaceblsaw(blsaw* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "freq", &dsp->fHslider1, 440.f, 0.f, 20000.f, 0.0001f);
	interface->addHorizontalSlider(interface->uiInterface, "amp", &dsp->fHslider0, 1.f, 0.f, 1.f, 0.0001f);
}

void computeblsaw(blsaw* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fHslider1;
	float fSlow1 = ((float)dsp->iConst0 * ((float)dsp->fHslider0 / fSlow0));
	float fSlow2 = (dsp->fConst2 * fSlow0);
	float fSlow3 = (dsp->fConst1 / fSlow0);
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			dsp->fRec0[0] = fmod((1.f + dsp->fRec0[1]), fSlow3);
			float fTemp0 = faustpower2_f(((fSlow2 * dsp->fRec0[0]) - 1.f));
			dsp->fVec0[0] = fTemp0;
			dsp->fVec1[0] = 0.25f;
			output0[i] = (FAUSTFLOAT)(fSlow1 * ((fTemp0 - dsp->fVec0[1]) * dsp->fVec1[1]));
			dsp->fRec0[1] = dsp->fRec0[0];
			dsp->fVec0[1] = dsp->fVec0[0];
			dsp->fVec1[1] = dsp->fVec1[0];

		}

	}

}


int sp_blsaw_create(sp_blsaw **p)
{
    *p = malloc(sizeof(sp_blsaw));
    return SP_OK;
}

int sp_blsaw_destroy(sp_blsaw **p)
{
    sp_blsaw *pp = *p;
    blsaw *dsp = pp->ud;
    deleteblsaw (dsp);
    free(*p);
    return SP_OK;
}

int sp_blsaw_init(sp_data *sp, sp_blsaw *p)
{
    blsaw *dsp = newblsaw();
    UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfaceblsaw(dsp, &UI);
    initblsaw(dsp, sp->sr);


    p->freq = p->args[0];
    p->amp = p->args[1];

    p->ud = dsp;
    return SP_OK;
}

int sp_blsaw_compute(sp_data *sp, sp_blsaw *p, SPFLOAT *in, SPFLOAT *out)
{

    blsaw *dsp = p->ud;
    SPFLOAT out1 = 0;
    SPFLOAT *faust_out[] = {&out1};
    SPFLOAT *faust_in[] = {in};
    computeblsaw(dsp, 1, faust_in, faust_out);

    *out = out1;
    return SP_OK;
}

typedef struct {
	float fVec2[4096];
	int iVec0[2];
	float fRec0[2];
	float fVec1[2];
	FAUSTFLOAT fHslider0;
	int fSamplingFreq;
	float fConst0;
	FAUSTFLOAT fHslider1;
	FAUSTFLOAT fHslider2;
	float fConst1;
	float fConst2;
	int IOTA;
} blsquare;

blsquare* newblsquare() {
	blsquare* dsp = (blsquare*)malloc(sizeof(blsquare));
	return dsp;
}

void deleteblsquare(blsquare* dsp) {
	free(dsp);
}


void instanceInitblsquare(blsquare* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->fHslider0 = (FAUSTFLOAT)1.;
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->iVec0[i0] = 0;

		}

	}
	dsp->fConst0 = (float)min(192000, max(1, dsp->fSamplingFreq));
	dsp->fHslider1 = (FAUSTFLOAT)0.5;
	dsp->fHslider2 = (FAUSTFLOAT)440.;
	dsp->fConst1 = (0.25f * dsp->fConst0);
	dsp->fConst2 = (1.f / dsp->fConst0);
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec0[i1] = 0.f;

		}

	}
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fVec1[i2] = 0.f;

		}

	}
	dsp->IOTA = 0;
	/* C99 loop */
	{
		int i3;
		for (i3 = 0; (i3 < 4096); i3 = (i3 + 1)) {
			dsp->fVec2[i3] = 0.f;

		}

	}

}

void initblsquare(blsquare* dsp, int samplingFreq) {
	instanceInitblsquare(dsp, samplingFreq);
}

void buildUserInterfaceblsquare(blsquare* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "frequency", &dsp->fHslider2, 440.f, 0.f, 20000.f, 0.0001f);
	interface->addHorizontalSlider(interface->uiInterface, "amp", &dsp->fHslider0, 1.f, 0.f, 1.f, 1e-05f);
	interface->addHorizontalSlider(interface->uiInterface, "width", &dsp->fHslider1, 0.5f, 0.f, 1.f, 0.f);
}

void computeblsquare(blsquare* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fHslider0;
	float fSlow1 = max((float)dsp->fHslider2, 23.4489f);
	float fSlow2 = max(0.f, min(2047.f, (dsp->fConst0 * ((float)dsp->fHslider1 / fSlow1))));
	int iSlow3 = (int)fSlow2;
	int iSlow4 = (1 + iSlow3);
	float fSlow5 = ((float)iSlow4 - fSlow2);
	float fSlow6 = (dsp->fConst1 / fSlow1);
	float fSlow7 = (dsp->fConst2 * fSlow1);
	float fSlow8 = (fSlow2 - (float)iSlow3);
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			dsp->iVec0[0] = 1;
			dsp->fRec0[0] = fmodf((dsp->fRec0[1] + fSlow7), 1.f);
			float fTemp0 = faustpower2_f(((2.f * dsp->fRec0[0]) - 1.f));
			dsp->fVec1[0] = fTemp0;
			float fTemp1 = (fSlow6 * ((fTemp0 - dsp->fVec1[1]) * (float)dsp->iVec0[1]));
			dsp->fVec2[(dsp->IOTA & 4095)] = fTemp1;
			output0[i] = (FAUSTFLOAT)(fSlow0 * (0.f - (((fSlow5 * dsp->fVec2[((dsp->IOTA - iSlow3) & 4095)]) + (fSlow8 * dsp->fVec2[((dsp->IOTA - iSlow4) & 4095)])) - fTemp1)));
			dsp->iVec0[1] = dsp->iVec0[0];
			dsp->fRec0[1] = dsp->fRec0[0];
			dsp->fVec1[1] = dsp->fVec1[0];
			dsp->IOTA = (dsp->IOTA + 1);

		}

	}

}


int sp_blsquare_create(sp_blsquare **p)
{
    *p = malloc(sizeof(sp_blsquare));
    return SP_OK;
}

int sp_blsquare_destroy(sp_blsquare **p)
{
    sp_blsquare *pp = *p;
    blsquare *dsp = pp->ud;
    deleteblsquare (dsp);
    free(*p);
    return SP_OK;
}

int sp_blsquare_init(sp_data *sp, sp_blsquare *p)
{
    blsquare *dsp = newblsquare(); UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfaceblsquare(dsp, &UI);
    initblsquare(dsp, sp->sr);


    p->freq = p->args[0];
    p->amp = p->args[1];
    p->width = p->args[2];

    p->ud = dsp;
    return SP_OK;
}

int sp_blsquare_compute(sp_data *sp, sp_blsquare *p, SPFLOAT *in, SPFLOAT *out)
{

    blsquare *dsp = p->ud;
    SPFLOAT out1 = 0;
    SPFLOAT *faust_out[] = {&out1};
    SPFLOAT *faust_in[] = {in};
    computeblsquare(dsp, 1, faust_in, faust_out);

    *out = out1;
    return SP_OK;
}


typedef struct {

	float fVec2[4096];
	int iVec0[2];
	float fRec1[2];
	float fVec1[2];
	float fRec0[2];
	int fSamplingFreq;
	int iConst0;
	float fConst1;
	FAUSTFLOAT fHslider0;
	FAUSTFLOAT fHslider1;
	float fConst2;
	float fConst3;
	float fConst4;
	float fConst5;
	int IOTA;

} bltriangle;

bltriangle* newbltriangle() {
	bltriangle* dsp = (bltriangle*)malloc(sizeof(bltriangle));
	return dsp;
}

void deletebltriangle(bltriangle* dsp) {
	free(dsp);
}

void instanceInitbltriangle(bltriangle* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->iVec0[i0] = 0;

		}

	}
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fConst1 = (4.f / (float)dsp->iConst0);
	dsp->fHslider0 = (FAUSTFLOAT)440.;
	dsp->fHslider1 = (FAUSTFLOAT)1.;
	dsp->fConst2 = (float)dsp->iConst0;
	dsp->fConst3 = (0.5f * dsp->fConst2);
	dsp->fConst4 = (0.25f * dsp->fConst2);
	dsp->fConst5 = (1.f / dsp->fConst2);
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec1[i1] = 0.f;

		}

	}
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fVec1[i2] = 0.f;

		}

	}
	dsp->IOTA = 0;
	/* C99 loop */
	{
		int i3;
		for (i3 = 0; (i3 < 4096); i3 = (i3 + 1)) {
			dsp->fVec2[i3] = 0.f;

		}

	}
	/* C99 loop */
	{
		int i4;
		for (i4 = 0; (i4 < 2); i4 = (i4 + 1)) {
			dsp->fRec0[i4] = 0.f;

		}

	}

}

void initbltriangle(bltriangle* dsp, int samplingFreq) {
	instanceInitbltriangle(dsp, samplingFreq);
}

void buildUserInterfacebltriangle(bltriangle* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "freq", &dsp->fHslider0, 440.f, 0.f, 20000.f, 0.0001f);
	interface->addHorizontalSlider(interface->uiInterface, "amp", &dsp->fHslider1, 1.f, 0.f, 1.f, 1e-05f);
}

void computebltriangle(bltriangle* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fHslider0;
	float fSlow1 = (dsp->fConst1 * (fSlow0 * (float)dsp->fHslider1));
	float fSlow2 = max(fSlow0, 23.4489f);
	float fSlow3 = max(0.f, min(2047.f, (dsp->fConst3 / fSlow2)));
	int iSlow4 = (int)fSlow3;
	int iSlow5 = (1 + iSlow4);
	float fSlow6 = ((float)iSlow5 - fSlow3);
	float fSlow7 = (dsp->fConst4 / fSlow2);
	float fSlow8 = (dsp->fConst5 * fSlow2);
	float fSlow9 = (fSlow3 - (float)iSlow4);
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			dsp->iVec0[0] = 1;
			dsp->fRec1[0] = fmodf((dsp->fRec1[1] + fSlow8), 1.f);
			float fTemp0 = faustpower2_f(((2.f * dsp->fRec1[0]) - 1.f));
			dsp->fVec1[0] = fTemp0;
			float fTemp1 = (fSlow7 * ((fTemp0 - dsp->fVec1[1]) * (float)dsp->iVec0[1]));
			dsp->fVec2[(dsp->IOTA & 4095)] = fTemp1;
			dsp->fRec0[0] = (0.f - (((fSlow6 * dsp->fVec2[((dsp->IOTA - iSlow4) & 4095)]) + (fSlow9 * dsp->fVec2[((dsp->IOTA - iSlow5) & 4095)])) - ((0.999f * dsp->fRec0[1]) + fTemp1)));
			output0[i] = (FAUSTFLOAT)(fSlow1 * dsp->fRec0[0]);
			dsp->iVec0[1] = dsp->iVec0[0];
			dsp->fRec1[1] = dsp->fRec1[0];
			dsp->fVec1[1] = dsp->fVec1[0];
			dsp->IOTA = (dsp->IOTA + 1);
			dsp->fRec0[1] = dsp->fRec0[0];

		}

	}

}


int sp_bltriangle_create(sp_bltriangle **p)
{
    *p = malloc(sizeof(sp_bltriangle));
    return SP_OK;
}

int sp_bltriangle_destroy(sp_bltriangle **p)
{
    sp_bltriangle *pp = *p;
    bltriangle *dsp = pp->ud;
    deletebltriangle (dsp);
    free(*p);
    return SP_OK;
}

int sp_bltriangle_init(sp_data *sp, sp_bltriangle *p)
{
    bltriangle *dsp = newbltriangle(); UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfacebltriangle(dsp, &UI);
    initbltriangle(dsp, sp->sr);


    p->freq = p->args[0];
    p->amp = p->args[1];

    p->ud = dsp;
    return SP_OK;
}

int sp_bltriangle_compute(sp_data *sp, sp_bltriangle *p, SPFLOAT *in, SPFLOAT *out)
{

    bltriangle *dsp = p->ud;
    SPFLOAT out1 = 0;
    SPFLOAT *faust_out[] = {&out1};
    SPFLOAT *faust_in[] = {in};
    computebltriangle(dsp, 1, faust_in, faust_out);

    *out = out1;
    return SP_OK;
}


int sp_brown_create(sp_brown **p)
{
    *p = malloc(sizeof(sp_brown));
    return SP_OK;
}

int sp_brown_destroy(sp_brown **p)
{
    free(*p);
    return SP_OK;
}

int sp_brown_init(sp_data *sp, sp_brown *p)
{
    p->brown = 0.0;
    return SP_OK;
}

int sp_brown_compute(sp_data *sp, sp_brown *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT r;
    while(1) {
        r = (sp_rand(sp) % SP_RANDMAX) / (SPFLOAT)(SP_RANDMAX);
        r = ((r * 2) - 1) * 0.5;
        p->brown += r;
        if(p->brown < -8.0f || p->brown > 8.0f) {
            p->brown -= r;
        } else {
            break;
        }
    }

    *out = p->brown * 0.0625;
    return SP_OK;
}



int sp_butbp_create(sp_butbp **p)
{
    *p = malloc(sizeof(sp_butbp));
    return SP_OK;
}

int sp_butbp_destroy(sp_butbp **p)
{
    free(*p);
    return SP_OK;
}

int sp_butbp_init(sp_data *sp, sp_butbp *p)
{
    p->istor = 0.0;
    p->sr = sp->sr;
    p->freq = 1000;
    p->bw = 10;
    p->pidsr = M_PI / sp->sr * 1.0;
    p->tpidsr = 2 * M_PI / sp->sr * 1.0;
    p->a[6] = p->a[7] = 0.0;
    p->lkf = 0.0;
    p->lkb = 0.0;
    return SP_OK;
}

int sp_butbp_compute(sp_data *sp, sp_butbp *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT *a = p->a;
    SPFLOAT t, y;

    if (p->bw <= 0.0) {
       *out = 0;
       return SP_OK;
    }

    SPFLOAT bw, fr;
    bw = p->bw;
    fr = p->freq;

    if (bw != p->lkb || fr != p->lkf) {
        SPFLOAT c, d;
        p->lkf = fr;
        p->lkb = bw;
        c = 1.0 / tan((SPFLOAT)(p->pidsr * bw));
        d = 2.0 * cos((SPFLOAT)(p->tpidsr * fr));
        a[1] = 1.0 / (1.0 + c);
        a[2] = 0.0;
        a[3] = -a[1];
        a[4] = - c * d * a[1];
        a[5] = (c - 1.0) * a[1];
    }
    t = *in - a[4] * a[6] - a[5] * a[7];
    y = t * a[1] + a[2] * a[6] + a[3] * a[7];
    a[7] = a[6];
    a[6] = t;
    *out = y;
    return SP_OK;
}

int sp_butbr_create(sp_butbr **p)
{
    *p = malloc(sizeof(sp_butbr));
    return SP_OK;
}

int sp_butbr_destroy(sp_butbr **p)
{
    free(*p);
    return SP_OK;
}

int sp_butbr_init(sp_data *sp, sp_butbr *p)
{
    p->istor = 0.0;
    p->sr = sp->sr;
    p->freq = 1000;
    p->bw = 1000;
    p->pidsr = M_PI / sp->sr * 1.0;
    p->tpidsr = 2 * M_PI / sp->sr * 1.0;
    p->a[6] = p->a[7] = 0.0;
    p->lkf = 0.0;
    p->lkb = 0.0;
    return SP_OK;
}

int sp_butbr_compute(sp_data *sp, sp_butbr *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT *a = p->a;
    SPFLOAT t, y;

    if (p->bw <= 0.0) {
      *out = 0;
      return SP_OK;
    }

    SPFLOAT bw, fr;
    bw = p->bw;
    fr = p->freq;
    if (bw != p->lkb || fr != p->lkf) {
        SPFLOAT c, d;
        p->lkf = fr;
        p->lkb = bw;
        c = tan((SPFLOAT)(p->pidsr * bw));
        d = 2.0 * cos((SPFLOAT)(p->tpidsr * fr));
        a[1] = 1.0 / (1.0 + c);
        a[2] = - d * a[1];
        a[3] = a[1];
        a[4] = a[2];
        a[5] = (1.0 - c) * a[1];
    }
    t = (SPFLOAT)*in - a[4] * a[6] - a[5] * a[7];
    y = t * a[1] + a[2] * a[6] + a[3] * a[7];
    a[7] = a[6];
    a[6] = t;
    *out = y;
    return SP_OK;
}


/* Filter loop */

static int sp_butter_filter(SPFLOAT *in, SPFLOAT *out, SPFLOAT *a)
{
    SPFLOAT t, y;
    t = *in - a[4] * a[6] - a[5] * a[7];
    y = t * a[1] + a[2] * a[6] + a[3] * a[7];
    a[7] = a[6];
    a[6] = t;
    *out = y;
    return SP_OK;
}


int sp_buthp_create(sp_buthp **p)
{
    *p = malloc(sizeof(sp_buthp));
    return SP_OK;
}

int sp_buthp_destroy(sp_buthp **p)
{
    free(*p);
    return SP_OK;
}

int sp_buthp_init(sp_data *sp, sp_buthp *p)
{
    p->istor = 0.0;
    p->sr = sp->sr;
    p->freq = 1000;
    p->pidsr = M_PI / sp->sr * 1.0;
    if (p->istor==0.0) {
        p->a[6] = p->a[7] = 0.0;
        p->lkf = 0.0;
    }
    return SP_OK;
}

int sp_buthp_compute(sp_data *sp, sp_buthp *p, SPFLOAT *in, SPFLOAT *out)
{
    if (p->freq <= 0.0)     {
      *out = 0;
      return SP_OK;
    }

    if (p->freq != p->lkf)      {
      SPFLOAT *a, c;
      a = p->a;
      p->lkf = p->freq;
      c = tan((SPFLOAT)(p->pidsr * p->lkf));

      a[1] = 1.0 / ( 1.0 + ROOT2 * c + c * c);
      a[2] = -(a[1] + a[1]);
      a[3] = a[1];
      a[4] = 2.0 * ( c*c - 1.0) * a[1];
      a[5] = ( 1.0 - ROOT2 * c + c * c) * a[1];
    }
    sp_butter_filter(in, out, p->a);
    return SP_OK;
}



int sp_butlp_create(sp_butlp **p)
{
    *p = malloc(sizeof(sp_butlp));
    return SP_OK;
}

int sp_butlp_destroy(sp_butlp **p)
{
    free(*p);
    return SP_OK;
}

int sp_butlp_init(sp_data *sp, sp_butlp *p)
{
    p->istor = 0.0;
    p->sr = sp->sr;
    p->freq = 1000;
    p->pidsr = M_PI / sp->sr * 1.0;
    if (p->istor==0.0) {
        p->a[6] = p->a[7] = 0.0;
        p->lkf = 0.0;
    }
    return SP_OK;
}

int sp_butlp_compute(sp_data *sp, sp_butlp *p, SPFLOAT *in, SPFLOAT *out)
{
    if (p->freq <= 0.0){
      *out = 0;
      return SP_OK;
    }

    if (p->freq != p->lkf){
        SPFLOAT *a, c;
        a = p->a;
        p->lkf = p->freq;
        c = 1.0 / tan((SPFLOAT)(p->pidsr * p->lkf));
        a[1] = 1.0 / ( 1.0 + ROOT2 * c + c * c);
        a[2] = a[1] + a[1];
        a[3] = a[1];
        a[4] = 2.0 * ( 1.0 - c*c) * a[1];
        a[5] = ( 1.0 - ROOT2 * c + c * c) * a[1];
    }

    sp_butter_filter(in, out, p->a);
    return SP_OK;
}


int sp_clip_create(sp_clip **p)
{
    *p = malloc(sizeof(sp_clip));
    return SP_OK;
}

int sp_clip_destroy(sp_clip **p)
{
    free(*p);
    return SP_OK;
}

int sp_clip_init(sp_data *sp, sp_clip *p)
{
    p->lim = 1;
    p->k1 = M_PI / (2.0 * p->lim);
    return SP_OK;
}

int sp_clip_compute(sp_data *sp, sp_clip *p, SPFLOAT *in, SPFLOAT *out)
{
    p->k1 = M_PI / (2.0 * p->lim);
    SPFLOAT k1 = p->k1;
    SPFLOAT limit = p->lim;
    SPFLOAT x;

    x = *in;
    if (x >= limit) {
        x = limit;
    } else if (x <= -limit) {
        x = -limit;
    } else {
        x = limit * sin(k1 * x);
    }
    *out = x;

    return SP_OK;
}


int sp_clock_create(sp_clock **p)
{
    *p = malloc(sizeof(sp_clock));
    return SP_OK;
}

int sp_clock_destroy(sp_clock **p)
{
    free(*p);
    return SP_OK;
}

int sp_clock_init(sp_data *sp, sp_clock *p)
{
    p->subdiv = 1.0;
    p->bpm = 120;
    p->counter = 0;
    return SP_OK;
}

int sp_clock_compute(sp_data *sp, sp_clock *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0.0;
    if(p->counter == 0 || *in != 0) {
        *out = 1.0;
        p->counter = (int)(sp->sr * (60.0 / (p->bpm * p->subdiv))) + 1;
    }
    p->counter--; 
    return SP_OK;
}




int sp_comb_create(sp_comb **p)
{
    *p = malloc(sizeof(sp_comb));
    return SP_OK;
}

int sp_comb_destroy(sp_comb **p)
{
    sp_comb *pp = *p;
    sp_auxdata_free(&pp->aux);
    free(*p);
    return SP_OK;
}

int sp_comb_init(sp_data *sp, sp_comb *p, SPFLOAT looptime)
{
    p->revtime = 3.5;
    p->looptime = looptime;
    p->bufsize = (uint32_t) (0.5 + looptime * sp->sr);
    sp_auxdata_alloc(&p->aux, p->bufsize * sizeof(SPFLOAT));
    p->prvt = 0.0;
    p->coef = 0.0;
    p->bufpos = 0;
    return SP_OK;
}

int sp_comb_compute(sp_data *sp, sp_comb *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT tmp = 0;
    SPFLOAT coef = p->coef;
    SPFLOAT outsamp = 0;
    SPFLOAT *buf = (SPFLOAT *)p->aux.ptr;

    if(p->prvt != p->revtime) {
        p->prvt = p->revtime;
        SPFLOAT exp_arg = (SPFLOAT) (log001 * p->looptime / p->prvt);
        if(exp_arg < -36.8413615) {
            coef = p->coef = 0;
        } else {
            coef = p->coef = exp(exp_arg);
        }
    }
    outsamp = buf[p->bufpos];
    tmp = outsamp;
    tmp *= coef;
    tmp += *in;
    buf[p->bufpos] = tmp;
    *out = outsamp;

    p->bufpos++;
    p->bufpos %= p->bufsize; 
    return SP_OK;
}


typedef struct {
	float fRec2[2];
	float fRec1[2];
	float fRec0[2];
	int fSamplingFreq;
	int iConst0;
	float fConst1;
	FAUSTFLOAT fHslider0;
	FAUSTFLOAT fHslider1;
	float fConst2;
	FAUSTFLOAT fHslider2;
	FAUSTFLOAT fHslider3;
} compressor;

static compressor* newcompressor() { 
	compressor* dsp = (compressor*)malloc(sizeof(compressor));
	return dsp;
}

static void deletecompressor(compressor* dsp) { 
	free(dsp);
}

static void instanceInitcompressor(compressor* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fConst1 = (2.f / (float)dsp->iConst0);
	dsp->fHslider0 = (FAUSTFLOAT)0.1;
	dsp->fHslider1 = (FAUSTFLOAT)1.;
	dsp->fConst2 = (1.f / (float)dsp->iConst0);
	dsp->fHslider2 = (FAUSTFLOAT)0.1;
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->fRec2[i0] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec1[i1] = 0.f;
			
		}
		
	}
	dsp->fHslider3 = (FAUSTFLOAT)0.;
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fRec0[i2] = 0.f;
			
		}
		
	}
	
}

static void initcompressor(compressor* dsp, int samplingFreq) {
	instanceInitcompressor(dsp, samplingFreq);
}

static void buildUserInterfacecompressor(compressor* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "ratio", &dsp->fHslider1, 1.f, 1.f, 40.f, 0.001f);
	interface->addHorizontalSlider(interface->uiInterface, "thresh", &dsp->fHslider3, 0.f, -80.f, 0.f, 0.001f);
	interface->addHorizontalSlider(interface->uiInterface, "atk", &dsp->fHslider0, 0.1f, 0.f, 10.f, 0.001f);
	interface->addHorizontalSlider(interface->uiInterface, "rel", &dsp->fHslider2, 0.1f, 0.f, 10.f, 0.001f);
}

static void computecompressor(compressor* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* input0 = inputs[0];
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fHslider0;
	float fSlow1 = exp((0.f - (dsp->fConst1 / fSlow0)));
	float fSlow2 = ((1.f - fSlow1) * ((1.f / (float)dsp->fHslider1) - 1.f));
	float fSlow3 = exp((0.f - (dsp->fConst2 / fSlow0)));
	float fSlow4 = exp((0.f - (dsp->fConst2 / (float)dsp->fHslider2)));
	float fSlow5 = (float)dsp->fHslider3;
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			float fTemp0 = (float)input0[i];
			float fTemp1 = fabs(fTemp0);
			float fTemp2 = ((dsp->fRec1[1] > fTemp1)?fSlow4:fSlow3);
			dsp->fRec2[0] = ((dsp->fRec2[1] * fTemp2) + ((1.f - fTemp2) * fTemp1));
			dsp->fRec1[0] = dsp->fRec2[0];
			dsp->fRec0[0] = ((fSlow1 * dsp->fRec0[1]) + (fSlow2 * max(((20.f * log10(dsp->fRec1[0])) - fSlow5), 0.f)));
			output0[i] = (FAUSTFLOAT)(pow(10.f, (0.05f * dsp->fRec0[0])) * fTemp0);
			dsp->fRec2[1] = dsp->fRec2[0];
			dsp->fRec1[1] = dsp->fRec1[0];
			dsp->fRec0[1] = dsp->fRec0[0];
			
		}
		
	}
	
}


int sp_compressor_create(sp_compressor **p)
{
    *p = malloc(sizeof(sp_compressor));
    return SP_OK;
}

int sp_compressor_destroy(sp_compressor **p)
{
    sp_compressor *pp = *p;
    compressor *dsp = pp->faust;
    deletecompressor (dsp);
    free(*p);
    return SP_OK;
}

int sp_compressor_init(sp_data *sp, sp_compressor *p)
{
    compressor *dsp = newcompressor(); 
    UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfacecompressor(dsp, &UI);
    initcompressor(dsp, sp->sr);

     
    p->ratio = p->args[0]; 
    p->thresh = p->args[1]; 
    p->atk = p->args[2]; 
    p->rel = p->args[3];

    p->faust = dsp;
    return SP_OK;
}

int sp_compressor_compute(sp_data *sp, sp_compressor *p, SPFLOAT *in, SPFLOAT *out) 
{

    compressor *dsp = p->faust;
    SPFLOAT *faust_out[] = {out};
    SPFLOAT *faust_in[] = {in};
    computecompressor(dsp, 1, faust_in, faust_out);
    return SP_OK;
}

static void multiply_fft_buffers(SPFLOAT *outBuf, SPFLOAT *ringBuf,
                                 SPFLOAT *IR_Data, int partSize, int nPartitions,
                                 int ringBuf_startPos)
{
    SPFLOAT   re, im, re1, re2, im1, im2;
    SPFLOAT   *rbPtr, *irPtr, *outBufPtr, *outBufEndPm2, *rbEndP;

    /* note: partSize must be at least 2 samples */
    partSize <<= 1;
    outBufEndPm2 = (SPFLOAT*) outBuf + (int) (partSize - 2);
    rbEndP = (SPFLOAT*) ringBuf + (int) (partSize * nPartitions);
    rbPtr = &(ringBuf[ringBuf_startPos]);
    irPtr = IR_Data;
    outBufPtr = outBuf;
    memset(outBuf, 0, sizeof(SPFLOAT)*(partSize));
    do {
      /* wrap ring buffer position */
      if (rbPtr >= rbEndP)
        rbPtr = ringBuf;
      outBufPtr = outBuf;
      *(outBufPtr++) += *(rbPtr++) * *(irPtr++);    /* convolve DC */
      *(outBufPtr++) += *(rbPtr++) * *(irPtr++);    /* convolve Nyquist */
      re1 = *(rbPtr++);
      im1 = *(rbPtr++);
      re2 = *(irPtr++);
      im2 = *(irPtr++);
      re = re1 * re2 - im1 * im2;
      im = re1 * im2 + re2 * im1;
      while (outBufPtr < outBufEndPm2) {
        /* complex multiply */
        re1 = rbPtr[0];
        im1 = rbPtr[1];
        re2 = irPtr[0];
        im2 = irPtr[1];
        outBufPtr[0] += re;
        outBufPtr[1] += im;
        re = re1 * re2 - im1 * im2;
        im = re1 * im2 + re2 * im1;
        re1 = rbPtr[2];
        im1 = rbPtr[3];
        re2 = irPtr[2];
        im2 = irPtr[3];
        outBufPtr[2] += re;
        outBufPtr[3] += im;
        re = re1 * re2 - im1 * im2;
        im = re1 * im2 + re2 * im1;
        outBufPtr += 4;
        rbPtr += 4;
        irPtr += 4;
      }
      outBufPtr[0] += re;
      outBufPtr[1] += im;
    } while (--nPartitions);
}

static int buf_bytes_alloc(int nChannels, int partSize, int nPartitions)
{
    int nSmps;

    nSmps = (partSize << 1);                                /* tmpBuf     */
    nSmps += ((partSize << 1) * nPartitions);               /* ringBuf    */
    nSmps += ((partSize << 1) * nChannels * nPartitions);   /* IR_Data    */
    nSmps += ((partSize << 1) * nChannels);                 /* outBuffers */

    return ((int) sizeof(SPFLOAT) * nSmps);
}

static void set_buf_pointers(sp_conv *p,
                             int nChannels, int partSize, int nPartitions)
{
    SPFLOAT *ptr;
    int   i;

    ptr = (SPFLOAT *) (p->auxData.ptr);
    p->tmpBuf = ptr;
    ptr += (partSize << 1);
    p->ringBuf = ptr;
    ptr += ((partSize << 1) * nPartitions);
    for (i = 0; i < nChannels; i++) {
      p->IR_Data[i] = ptr;
      ptr += ((partSize << 1) * nPartitions);
    }
    for (i = 0; i < nChannels; i++) {
      p->outBuffers[i] = ptr;
      ptr += (partSize << 1);
    }
}

int sp_conv_create(sp_conv **p)
{
    *p = malloc(sizeof(sp_conv));
    return SP_OK;
}

int sp_conv_destroy(sp_conv **p)
{
    sp_conv *pp = *p;
    sp_auxdata_free(&pp->auxData);
    sp_fft_destroy(&pp->fft);
    free(*p);
    return SP_OK;
}

int sp_conv_init(sp_data *sp, sp_conv *p, sp_ftbl *ft, SPFLOAT iPartLen)
{
    int     i, j, k, n, nBytes, skipSamples;
    SPFLOAT FFTscale;

    p->iTotLen = ft->size;
    p->iSkipSamples = 0;
    p->iPartLen = iPartLen;

    p->nChannels = 1;
    /* partition length */
    p->partSize = (int)lrintf(p->iPartLen);
    if (p->partSize < 4 || (p->partSize & (p->partSize - 1)) != 0) {
        fprintf(stderr, "conv: invalid partition size.\n");
        return SP_NOT_OK;  
    }

    sp_fft_init(&p->fft, (int)log2(p->partSize << 1));
    n = (int) ft->size / p->nChannels;
    skipSamples = (int)lrintf(p->iSkipSamples);
    n -= skipSamples;

    if (lrintf(p->iTotLen) > 0 && n > lrintf(p->iTotLen)) {
        n = (int)lrintf(p->iTotLen);
    }

    if (n <= 0) {
        fprintf(stderr, "uh oh.\n");
        return SP_NOT_OK;
    }

    p->nPartitions = (n + (p->partSize - 1)) / p->partSize;
    /* calculate the amount of aux space to allocate (in bytes) */
    nBytes = buf_bytes_alloc(p->nChannels, p->partSize, p->nPartitions);
    sp_auxdata_alloc(&p->auxData, nBytes);
    /* if skipping samples: check for possible truncation of IR */
    /* initialise buffer pointers */
    set_buf_pointers(p, p->nChannels, p->partSize, p->nPartitions);
    /* clear ring buffer to zero */
    n = (p->partSize << 1) * p->nPartitions;
    memset(p->ringBuf, 0, n*sizeof(SPFLOAT));
    p->cnt = 0;
    p->rbCnt = 0;
    FFTscale = 1.0;
    for (j = 0; j < p->nChannels; j++) {
        /* table read position */
        i = (skipSamples * p->nChannels) + j; 
        /* IR write position */
        n = (p->partSize << 1) * (p->nPartitions - 1); 
        do { 
            for (k = 0; k < p->partSize; k++) {
                if (i >= 0 && i < (int) ft->size) {
                    p->IR_Data[j][n + k] = ft->tbl[i] * FFTscale;
                } else {
                    p->IR_Data[j][n + k] = 0.0;
                } 
                i += p->nChannels;
            }
        /* pad second half of IR to zero */
            for (k = p->partSize; k < (p->partSize << 1); k++) {
                p->IR_Data[j][n + k] = 0.0;
            }
            /* calculate FFT */
            sp_fftr(&p->fft, &(p->IR_Data[j][n]), (p->partSize << 1));
            n -= (p->partSize << 1);
        } while (n >= 0);
    }
    /* clear output buffers to zero */
    for (j = 0; j < p->nChannels; j++) {
        for (i = 0; i < (p->partSize << 1); i++)
        p->outBuffers[j][i] = 0.0;
    }
    p->initDone = 1;

    return SP_OK;
}

int sp_conv_compute(sp_data *sp, sp_conv *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT *x, *rBuf;
    int i, n, nSamples, rBufPos;

    nSamples = p->partSize;
    rBuf = &(p->ringBuf[p->rbCnt * (nSamples << 1)]);
    /* store input signal in buffer */
    rBuf[p->cnt] = *in;
    /* copy output signals from buffer */
    *out = p->outBuffers[0][p->cnt]; 

    /* is input buffer full ? */
    if (++p->cnt < nSamples) {
        return SP_OK;                   
    }
    /* reset buffer position */
    p->cnt = 0;
    /* calculate FFT of input */
    for (i = nSamples; i < (nSamples << 1); i++) {
        /* Zero padding */
        rBuf[i] = 0.0;
    }
    sp_fftr(&p->fft, rBuf, (nSamples << 1));
    /* update ring buffer position */
    p->rbCnt++;

    if (p->rbCnt >= p->nPartitions){ 
        p->rbCnt = 0; 
    }

    rBufPos = p->rbCnt * (nSamples << 1);
    rBuf = &(p->ringBuf[rBufPos]);
    /* PB: will only loop once since nChannels == 1*/
    for (n = 0; n < p->nChannels; n++) {
        /* multiply complex arrays */
        multiply_fft_buffers(p->tmpBuf, p->ringBuf, p->IR_Data[n],
                     nSamples, p->nPartitions, rBufPos);
        /* inverse FFT */
        sp_ifftr(&p->fft, p->tmpBuf, (nSamples << 1));
        /* copy to output buffer, overlap with "tail" of previous block */
        x = &(p->outBuffers[n][0]);
        for (i = 0; i < nSamples; i++) {
            x[i] = p->tmpBuf[i] + x[i + nSamples];
            x[i + nSamples] = p->tmpBuf[i + nSamples];
        }
    }
    return SP_OK;
}


int sp_count_create(sp_count **p)
{
    *p = malloc(sizeof(sp_count));
    return SP_OK;
}

int sp_count_destroy(sp_count **p)
{
    free(*p);
    return SP_OK;
}

int sp_count_init(sp_data *sp, sp_count *p)
{
    p->count = 4;
    p->curcount = -1;
    p->mode = 0;
    return SP_OK;
}

int sp_count_compute(sp_data *sp, sp_count *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in){
        if(p->mode == 0) {
            p->curcount = (p->curcount + 1) % p->count;
        } else {
            if(p->curcount == -2) {
                *out = -2;
                return SP_OK;
            }
            if(p->curcount >= p->count - 1) {
                p->curcount = -2;
            } else {
                if(p->curcount == -1) p->curcount = 0;
                else p->curcount++;
            }
        }
    }
    *out = p->curcount;
    return SP_OK;
}

int sp_crossfade_create(sp_crossfade **p)
{
    *p = malloc(sizeof(sp_crossfade));
    return SP_OK;
}

int sp_crossfade_destroy(sp_crossfade **p)
{
    free(*p);
    return SP_OK;
}

int sp_crossfade_init(sp_data *sp, sp_crossfade *p)
{
    p->pos = 0.5;
    return SP_OK;
}

int sp_crossfade_compute(sp_data *sp, sp_crossfade *p, SPFLOAT *in1, SPFLOAT *in2, SPFLOAT *out)
{
    *out = *in2 * p->pos + *in1 * (1 - p->pos);
    return SP_OK;
}

int sp_dcblock_create(sp_dcblock **p)
{
    *p = malloc(sizeof(sp_dcblock));
    return SP_OK;
}

int sp_dcblock_destroy(sp_dcblock **p)
{
    free(*p);
    return SP_OK;
}

int sp_dcblock_init(sp_data *sp, sp_dcblock *p)
{
    p->outputs = 0.0;
    p->inputs = 0.0;
    p->gain = 0.99;
    if (p->gain == 0.0 || p->gain>=1.0 || p->gain<=-1.0)
      p->gain = 0.99;
    return SP_OK;
}

int sp_dcblock_compute(sp_data *sp, sp_dcblock *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT gain = p->gain;
    SPFLOAT outputs = p->outputs;
    SPFLOAT inputs = p->inputs;

    SPFLOAT sample = *in;
    outputs = sample - inputs + (gain * outputs);
    inputs = sample;
    *out = outputs;
    p->outputs = outputs;
    p->inputs = inputs;
    return SP_OK;
}

int sp_delay_create(sp_delay **p)
{
    *p = malloc(sizeof(sp_delay));
    return SP_OK;
}

int sp_delay_destroy(sp_delay **p)
{
    sp_delay *pp = *p;
    sp_auxdata_free(&pp->buf);
    free(*p);
    return SP_OK;
}

int sp_delay_init(sp_data *sp, sp_delay *p, SPFLOAT time)
{
    p->time = time;
    p->bufsize = round(time * sp->sr);
    sp_auxdata_alloc(&p->buf, p->bufsize * sizeof(SPFLOAT));
    p->bufpos = 0;
    p->feedback = 0;
    p->last = 0;
    return SP_OK;
}

int sp_delay_compute(sp_data *sp, sp_delay *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT delay = 0, sig = 0;
    SPFLOAT *buf = (SPFLOAT *)p->buf.ptr; 
    delay = buf[p->bufpos];
    sig = (delay * p->feedback) + *in;
    buf[p->bufpos] = sig;
    p->bufpos = (p->bufpos + 1) % p->bufsize;
    *out = delay;
    return SP_OK;
}


int sp_diode_create(sp_diode **p)
{
    *p = malloc(sizeof(sp_diode));
    return SP_OK;
}

int sp_diode_destroy(sp_diode **p)
{
    free(*p);
    return SP_OK;
}

static SPFLOAT sp_diode_opva_fdbk_out(sp_data *sp, sp_diode *p, int filt)
{
    return p->opva_beta[filt] * 
        (p->opva_z1[filt] + p->opva_fdbk[filt] * p->opva_delta[filt]);
}

static SPFLOAT sp_diode_opva_compute(sp_data *sp, sp_diode *p, SPFLOAT in, int filt)
{
    /*
	double x_in = (xn*m_dGamma + m_dFeedback + m_dEpsilon*getFeedbackOutput());
	
	double vn = (m_da0*x_in - m_dZ1)*m_dAlpha;
	double out = vn + m_dZ1;
	m_dZ1 = vn + out;
    */

	/* m_dBeta*(m_dZ1 + m_dFeedback*m_dDelta); */

	SPFLOAT x_in = (in*p->opva_gamma[filt]
        + p->opva_fdbk[filt]
        + p->opva_eps[filt] * sp_diode_opva_fdbk_out(sp, p, filt));
	SPFLOAT vn = (p->opva_a0[filt]*x_in - 
        p->opva_z1[filt])*p->opva_alpha[filt];
	SPFLOAT out = vn + p->opva_z1[filt];
	p->opva_z1[filt] = vn + out;
    return out;
}

static void sp_diode_update(sp_data *sp, sp_diode *p)
{
	/* calculate alphas */
	SPFLOAT G1, G2, G3, G4;
	SPFLOAT wd = 2*M_PI*p->freq;          
	SPFLOAT T  = 1/(SPFLOAT)sp->sr;             
	SPFLOAT wa = (2/T)*tan(wd*T/2); 
	SPFLOAT g = wa*T/2;  
    int i;

    /* Big G's */
	G4 = 0.5*g/(1.0 + g);
	G3 = 0.5*g/(1.0 + g - 0.5*g*G4);
	G2 = 0.5*g/(1.0 + g - 0.5*g*G3);
	G1 = g/(1.0 + g - g*G2);

    /* big G value gamma */
    p->gamma  = G4*G3*G2*G1;
	
    p->SG[0] =  G4*G3*G2; 
	p->SG[1] =  G4*G3; 
	p->SG[2] =  G4; 
	p->SG[3] =  1.0; 

	/* set alphas */
    for(i = 0; i < 4; i++) p->opva_alpha[i] = g/(1.0 + g);
	
    /* set betas */
    p->opva_beta[0] = 1.0/(1.0 + g - g*G2);
	p->opva_beta[1] = 1.0/(1.0 + g - 0.5*g*G3);
	p->opva_beta[2] = 1.0/(1.0 + g - 0.5*g*G4);
	p->opva_beta[3] = 1.0/(1.0 + g);

	/* set gammas */
	p->opva_gamma[0] = 1.0 + G1*G2;
	p->opva_gamma[1] = 1.0 + G2*G3;
	p->opva_gamma[2] = 1.0 + G3*G4;

    /* set deltas */
	p->opva_delta[0] = g;
	p->opva_delta[1] = 0.5*g;
	p->opva_delta[2] = 0.5*g;

    /* set epsilons */
	p->opva_eps[0] = G2;
	p->opva_eps[1] = G3;
	p->opva_eps[2] = G4;
}

int sp_diode_init(sp_data *sp, sp_diode *p)
{
    int i;
    /* initialize the 4 one-pole VA filters */

    for(i = 0; i < 4; i++) {
        p->opva_alpha[i] = 1.0;		
        p->opva_beta[i] = -1.0;		
        p->opva_gamma[i] = 1.0;
        p->opva_delta[i] = 0.0;
        p->opva_eps[i] = 1.0;
        p->opva_fdbk[i] = 0.0;
        p->opva_a0[i] = 1.0;
        p->opva_z1[i] = 0.0;

        p->SG[i] = 0.0;
    }

	p->gamma = 0.0;
	p->K = 0.0;

	/* Filter coeffs that are constant */
	/* set a0s */
    p->opva_a0[0] = 1.0;
    p->opva_a0[1] = 0.5;
    p->opva_a0[2] = 0.5;
    p->opva_a0[3] = 0.5;

	/* last LPF has no feedback path */
    p->opva_gamma[3] = 1.0;
    p->opva_delta[3] = 0.0;
    p->opva_eps[3] = 0.0;
    p->opva_fdbk[3] = 0.0;

    /* default cutoff to 1000hz */
    p->freq = 1000;
    p->res = 0;
    /* update filter coefs */

    sp_diode_update(sp, p);
    return SP_OK;
}

int sp_diode_compute(sp_data *sp, sp_diode *p, SPFLOAT *in, SPFLOAT *out)
{
    int i;
    SPFLOAT sigma;
    SPFLOAT un;
    SPFLOAT tmp = 0.0;

    /* update filter coefficients */
    p->K = p->res * 17;
    sp_diode_update(sp, p);

    p->opva_fdbk[2] = sp_diode_opva_fdbk_out(sp, p, 3);
    p->opva_fdbk[1] = sp_diode_opva_fdbk_out(sp, p, 2);
    p->opva_fdbk[0] = sp_diode_opva_fdbk_out(sp, p, 1);

    sigma = 
        p->SG[0] * sp_diode_opva_fdbk_out(sp, p, 0) +
        p->SG[1] * sp_diode_opva_fdbk_out(sp, p, 1) +
        p->SG[2] * sp_diode_opva_fdbk_out(sp, p, 2) +
        p->SG[3] * sp_diode_opva_fdbk_out(sp, p, 3);

    un = (*in - p->K * sigma) / (1 + p->K * p->gamma);
    tmp = un;
    for(i = 0; i < 4; i++) {
        tmp = sp_diode_opva_compute(sp, p, tmp, i);
    }
    *out = tmp;
    return SP_OK;
}

struct sp_diskin {
    SNDFILE *file;
    SF_INFO info;
    SPFLOAT buffer[1024];
    int bufpos;
    int loaded;
    int count;
};

int sp_diskin_create(sp_diskin **p)
{
    *p = malloc(sizeof(sp_diskin));
    return SP_OK;
}

int sp_diskin_destroy(sp_diskin **p)
{
    sp_diskin *pp = *p;
    if(pp->loaded) sf_close(pp->file);
    free(*p);
    return SP_OK;
}

int sp_diskin_init(sp_data *sp, sp_diskin *p, const char *filename)
{
    p->info.format = 0;
    memset(&p->info, 0, sizeof(SF_INFO));
    p->file = sf_open(filename, SFM_READ, &p->info);
    p->loaded = 0;
    p->bufpos = 0;

    if(p->file == NULL) {
        fprintf(stderr, "Error: could not open file \"%s\"\n", filename);
        exit(1);
    }

    if(p->info.channels != 1) {
        fprintf(stderr, "Warning: file \"%s\" has %d channels," 
                "when it is expecting only 1\n", filename, p->info.channels);
    }

    p->loaded = 1;

    if(p->info.frames < 1024) {
        p->count = p->info.frames;
    } else {
        p->count = 1024;
    }
    memset(p->buffer, 0, sizeof(SPFLOAT) * 1024);
    return SP_OK;
}

int sp_diskin_compute(sp_data *sp, sp_diskin *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->bufpos == 0 && p->loaded && p->count > 0) {
#ifdef USE_DOUBLE
        p->count = sf_read_double(p->file, p->buffer, p->count);
#else
        p->count = sf_read_float(p->file, p->buffer, p->count);
#endif
    } 

    if(p->count <= 0) {
        *out = 0;
        return SP_OK;
    }

    *out = p->buffer[p->bufpos++];
    p->bufpos %= 1024; 
    return SP_OK;
}

int sp_dist_create(sp_dist **p)
{
    *p = malloc(sizeof(sp_dist));
    return SP_OK;
}

int sp_dist_destroy(sp_dist **p)
{
    free(*p);
    return SP_OK;
}

int sp_dist_init(sp_data *sp, sp_dist *p)
{
    p->pregain = 2.0;
    p->postgain = 0.5;
    p->shape1 = 0;
    p->shape2 = 0;
    return SP_OK;
}

int sp_dist_compute(sp_data *sp, sp_dist *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT pregain = p->pregain, postgain  = p->postgain;
    SPFLOAT shape1 = p->shape1, shape2 = p->shape2;
    SPFLOAT sig;
    
    pregain   *=  6.5536;
    postgain  *=  0.61035156;
    shape1    *=  4.096;
    shape2    *=  4.096;

    /* IV - Dec 28 2002 */
    shape1 += pregain;
    shape2 -= pregain;
    postgain *= 0.5;
    sig = *in;
    /* Generate tanh distortion and output the result */
    *out =                          
    ((exp(sig * shape1) - exp(sig * shape2))
             / cosh(sig * pregain))
    * postgain;
    return SP_OK;
}


int sp_dmetro_create(sp_dmetro **p)
{
    *p = malloc(sizeof(sp_dmetro));
    return SP_OK;
}

int sp_dmetro_destroy(sp_dmetro **p)
{
    free(*p);
    return SP_OK;
}

int sp_dmetro_init(sp_data *sp, sp_dmetro *p)
{
    p->counter = 0;
    p->time = 1.0;
    return SP_OK;
}

int sp_dmetro_compute(sp_data *sp, sp_dmetro *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0; 

    if(p->counter == 0) {
        *out = 1.0;
        p->counter = (int)(sp->sr * p->time) + 1;
    }

    p->counter--; 

    return SP_OK;
}




static int my_random(sp_data *sp, int max)
{                      
    return (sp_rand(sp) % (max + 1));
}

static SPFLOAT noise_tick(sp_data *sp)                                        
{                       
    SPFLOAT temp;                                                                
    temp = 1.0 * sp_rand(sp) - 1073741823.5;
    return temp * (1.0 / 1073741823.0);
}                                                                              

int sp_drip_create(sp_drip **p)
{
    *p = malloc(sizeof(sp_drip));
    return SP_OK;
}

int sp_drip_destroy(sp_drip **p)
{
    free(*p);
    return SP_OK;
}

int sp_drip_init(sp_data *sp, sp_drip *p, SPFLOAT dettack)
{

    SPFLOAT temp;
    p->dettack = dettack;
    p->num_tubes = 10;
    p->damp = 0.2;
    p->shake_max = 0;
    p->freq = 450.0;
    p->freq1 = 600.0;
    p->freq2 = 720.0;
    p->amp = 0.3;

    p->sndLevel = 0.0;
    SPFLOAT tpidsr = 2.0 * M_PI / sp->sr;

    p->kloop = (sp->sr * p->dettack);
    p->outputs00 = 0.0;
    p->outputs01 = 0.0;
    p->outputs10 = 0.0;
    p->outputs11 = 0.0;
    p->outputs20 = 0.0;
    p->outputs21 = 0.0;

    p->totalEnergy = 0.0;

    p->center_freqs0 = p->res_freq0 = WUTR_CENTER_FREQ0;
    p->center_freqs1 = p->res_freq1 = WUTR_CENTER_FREQ1;
    p->center_freqs2 = p->res_freq2 = WUTR_CENTER_FREQ2;
    p->num_objectsSave = p->num_objects = WUTR_NUM_SOURCES;
    p->soundDecay = WUTR_SOUND_DECAY;
    p->systemDecay = WUTR_SYSTEM_DECAY;
    temp = log(WUTR_NUM_SOURCES) * WUTR_GAIN / WUTR_NUM_SOURCES;
    p->gains0 = p->gains1 = p->gains2 = temp;
    p->coeffs01 = WUTR_RESON * WUTR_RESON;
    p->coeffs00 = -WUTR_RESON * 2.0 *
      cos(WUTR_CENTER_FREQ0 * tpidsr);
    p->coeffs11 = WUTR_RESON * WUTR_RESON;
    p->coeffs10 = -WUTR_RESON * 2.0 *
      cos(WUTR_CENTER_FREQ1 * tpidsr);
    p->coeffs21 = WUTR_RESON * WUTR_RESON;
    p->coeffs20 = -WUTR_RESON * 2.0 *
      cos(WUTR_CENTER_FREQ2 * tpidsr);
                                
    p->shakeEnergy = p->amp * 1.0 * MAX_SHAKE * 0.1;
    p->shake_damp = 0.0;
    if (p->shakeEnergy > MAX_SHAKE) p->shakeEnergy = MAX_SHAKE;
    p->shake_maxSave = 0.0;
    p->num_objects = 10;        
    p->finalZ0 = p->finalZ1 = p->finalZ2 = 0.0;
    return SP_OK;
}

int sp_drip_compute(sp_data *sp, sp_drip *p, SPFLOAT *trig, SPFLOAT *out)
{
    SPFLOAT data;
    SPFLOAT lastOutput;

    SPFLOAT tpidsr = 2.0 * M_PI / sp->sr;

    if(*trig) {
        sp_drip_init(sp, p, p->dettack);
    } 
    if (p->num_tubes != 0.0 && p->num_tubes != p->num_objects) {
        p->num_objects = p->num_tubes;
        if (p->num_objects < 1.0) p->num_objects = 1.0;
    }
    if (p->freq != 0.0 && p->freq != p->res_freq0) {
        p->res_freq0 = p->freq;
        p->coeffs00 = -WUTR_RESON * 2.0 *
        cos(p->res_freq0 * tpidsr);
    }
    if (p->damp != 0.0 && p->damp != p->shake_damp) {
        p->shake_damp = p->damp;
        p->systemDecay = WUTR_SYSTEM_DECAY + (p->shake_damp * 0.002);
    }
    if (p->shake_max != 0.0 && p->shake_max != p->shake_maxSave) {
        p->shake_maxSave = p->shake_max;
        p->shakeEnergy += p->shake_maxSave * MAX_SHAKE * 0.1;
        if (p->shakeEnergy > MAX_SHAKE) p->shakeEnergy = MAX_SHAKE;
    }
    if (p->freq1 != 0.0 && p->freq1 != p->res_freq1) {
        p->res_freq1 = p->freq1;
        p->coeffs10 = -WUTR_RESON * 2.0 *
        cos(p->res_freq1 * tpidsr);
    }
    if (p->freq2 != 0.0 && p->freq2 != p->res_freq2) {
        p->res_freq2 = p->freq2;
        p->coeffs20 = -WUTR_RESON * 2.0 *
        cos(p->res_freq2 * tpidsr);
    }
    if ((--p->kloop) == 0) {
        p->shakeEnergy = 0.0;
    }

    SPFLOAT shakeEnergy = p->shakeEnergy;
    SPFLOAT systemDecay = p->systemDecay;
    SPFLOAT sndLevel = p->sndLevel;
    SPFLOAT num_objects = p->num_objects;
    SPFLOAT soundDecay = p->soundDecay;
    SPFLOAT inputs0, inputs1, inputs2;

    shakeEnergy *= systemDecay; /* Exponential system decay */

    sndLevel = shakeEnergy;
    if (my_random(sp, 32767) < num_objects) {
        int j;
        j = my_random(sp, 3);
        if (j == 0) {
            p->center_freqs0 = p->res_freq1 *
            (0.75 + (0.25 * noise_tick(sp)));
            p->gains0 = fabs(noise_tick(sp));
        } else if (j == 1) {
            p->center_freqs1 = p->res_freq1 *
            (1.0 + (0.25 * noise_tick(sp)));
            p->gains1 = fabs(noise_tick(sp));
        } else  {
            p->center_freqs2 = p->res_freq1 *
            (1.25 + (0.25 * noise_tick(sp)));
            p->gains2 = fabs(noise_tick(sp));
        }
    }

    p->gains0 *= WUTR_RESON;
    if (p->gains0 > 0.001) {
        p->center_freqs0  *= WUTR_FREQ_SWEEP;
        p->coeffs00 = -WUTR_RESON * 2.0 *
        cos(p->center_freqs0 * tpidsr);
    }
    p->gains1 *= WUTR_RESON;
    if (p->gains1 > 0.00) {
        p->center_freqs1 *= WUTR_FREQ_SWEEP;
        p->coeffs10 = -WUTR_RESON * 2.0 *
        cos(p->center_freqs1 * tpidsr);
    }
    p->gains2 *= WUTR_RESON;
    if (p->gains2 > 0.001) {
        p->center_freqs2 *= WUTR_FREQ_SWEEP;
        p->coeffs20 = -WUTR_RESON * 2.0 *
        cos(p->center_freqs2 * tpidsr);
    }

    sndLevel *= soundDecay;   
    inputs0 = sndLevel;
    inputs0 *= noise_tick(sp); 
    inputs1 = inputs0 * p->gains1;
    inputs2 = inputs0 * p->gains2;
    inputs0 *= p->gains0;
    inputs0 -= p->outputs00*p->coeffs00;
    inputs0 -= p->outputs01*p->coeffs01;
    p->outputs01 = p->outputs00;
    p->outputs00 = inputs0;
    data = p->gains0*p->outputs00;
    inputs1 -= p->outputs10*p->coeffs10;
    inputs1 -= p->outputs11*p->coeffs11;
    p->outputs11 = p->outputs10;
    p->outputs10 = inputs1;
    data += p->gains1*p->outputs10;
    inputs2-= p->outputs20*p->coeffs20;
    inputs2 -= p->outputs21*p->coeffs21;
    p->outputs21 = p->outputs20;
    p->outputs20 = inputs2;
    data += p->gains2*p->outputs20;

    p->finalZ2 = p->finalZ1;
    p->finalZ1 = p->finalZ0;
    p->finalZ0 = data * 4.0;

    lastOutput = p->finalZ2 - p->finalZ0;
    lastOutput *= 0.005;
    *out = lastOutput;
    p->shakeEnergy = shakeEnergy;
    p->sndLevel = sndLevel;
    return SP_OK;
}


int sp_dtrig_create(sp_dtrig **p)
{
    *p = malloc(sizeof(sp_dtrig));
    return SP_OK;
}

int sp_dtrig_destroy(sp_dtrig **p)
{
    free(*p);
    return SP_OK;
}

int sp_dtrig_init(sp_data *sp, sp_dtrig *p, sp_ftbl *ft)
{
    p->ft = ft;
    p->counter = 0;
    p->pos = 0; 
    p->running = 0;
    p->loop = 0;
    p->delay = 0;
    p->scale = 1;
    return SP_OK;
}

int sp_dtrig_compute(sp_data *sp, sp_dtrig *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in == 1.0){
        p->running = 1.0;
        p->pos = 0;
        p->counter = p->delay * sp->sr;
    } 
    if((p->pos < p->ft->size) && p->running){
        if(p->counter == 0){
            p->counter = (uint32_t)(p->scale * p->ft->tbl[p->pos] * sp->sr - 1);
            *out = 1.0;
            p->pos++; 
            if(p->loop){
                p->pos %= p->ft->size;
            }
            return SP_OK;
        }else{
            *out = 0;
            p->counter--;
            return SP_OK;
        }
    }    
    *out = 0;
    return SP_NOT_OK;
}


int sp_dust_create(sp_dust **p)
{
    *p = malloc(sizeof(sp_dust));
    return SP_OK;
}

int sp_dust_destroy(sp_dust **p) 
{
    free(*p);
    return SP_OK;
}

int sp_dust_init(sp_data *sp, sp_dust *p) 
{
    p->density = 10;
    p->amp = 0.4;
    p->density0 = 0.0;
    p->thresh = 0.0;
    p->scale = 0.0;
    p->rand = sp_rand(sp);
    p->onedsr = 1.0 / sp->sr;
    p->bipolar = 0;
    return SP_OK;
}

int sp_dust_compute(sp_data *sp, sp_dust *p, SPFLOAT *in, SPFLOAT *out) 
{
    SPFLOAT density, thresh, scale;
    const SPFLOAT dv2_31 = 4.656612873077392578125e-10;
    density = p->density;
    if (density != p->density0) {
        thresh = p->thresh = density * p->onedsr;
        if(p->bipolar) {
            scale  = p->scale  = (thresh > 0.0 ? 2.0 / thresh : 0.0);
        } else {
            scale  = p->scale  = (thresh > 0.0 ? 1.0 / thresh : 0.0);
        }
        p->density0 = density;
    } else {
        thresh = p->thresh;
        scale  = p->scale;
    }
    *out = 0;
    SPFLOAT r;
    p->rand = sp_rand(sp);
    r = (SPFLOAT)p->rand * dv2_31;

    if(p->bipolar) {
        *out = p->amp * (r < thresh ? r*scale - 1.0 : 0.0);
    } else {
        *out = p->amp * (r < thresh ? r*scale : 0.0);
    }

    return SP_OK;
}


int sp_eqfil_create(sp_eqfil **p)
{
    *p = malloc(sizeof(sp_eqfil));
    return SP_OK;
}

int sp_eqfil_destroy(sp_eqfil **p)
{
    free(*p);
    return SP_OK;
}

int sp_eqfil_init(sp_data *sp, sp_eqfil *p)
{
    p->sr = sp->sr;
    p->z1 = p->z2 = 0.0;
    p->freq = 1000;
    p->bw = 125;
    p->gain = 2;

    p->frv = p->freq; p->bwv = p->bw;
    p->d = cos(2 * M_PI * p->frv /p->sr);
    p->c = tan(M_PI * p->bwv / p->sr);
    return SP_OK;
}

int sp_eqfil_compute(sp_data *sp, sp_eqfil *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT z1 = p->z1, z2 = p->z2, c, d, w, a, y;
    SPFLOAT g;

    if(p->bw != p->bwv || p->freq != p->frv) {
        SPFLOAT sr = sp->sr;
        p->frv = p->freq; p->bwv = p->bw;
        p->d = cos(2 * M_PI * p->frv / sr);
        p->c = tan(M_PI * p->bwv / sr);
    }

    c = p->c;
    d = p->d;
    a = (1.0 - c) / (1.0 + c);
    g = p->gain;

    w = *in + d * (1.0 + a) * z1 - a * z2;
    y = w * a - d * (1.0 + a) * z1 + z2;
    z2 = z1;
    z1 = w;
    *out = (0.5 * (y + *in + g * (*in - y)));

    p->z1 = z1;
    p->z2 = z2;

return SP_OK;
}


static SPFLOAT intpow1(SPFLOAT x, int32_t n)
{
    SPFLOAT ans = 1.0;
    while (n!=0) {
      if (n&1) ans = ans * x;
      n >>= 1;
      x = x*x;
    }
    return ans;
}

static SPFLOAT intpow(SPFLOAT x, int32_t n)
{
    if (n<0) {
      n = -n;
      x = 1.0/x;
    }
    return intpow1(x, n);
}


static int newpulse(sp_data *sp,
                    sp_fof *p, sp_fof_overlap *ovp, SPFLOAT amp, SPFLOAT fund, SPFLOAT form)
{
    SPFLOAT   octamp = amp, oct;
    int32_t   rismps, newexp = 0;

    ovp->timrem = p->dur * sp->sr;

    if ((oct = p->oct) > 0.0) {
        int32_t ioct = (int32_t)oct, bitpat = ~(-1L << ioct);
        if (bitpat & ++p->fofcount) return(0);
        if ((bitpat += 1) & p->fofcount) octamp *= (1.0 + ioct - oct);
    }
    if (fund == 0.0) ovp->formphs = 0;
    else ovp->formphs = (int32_t)(p->fundphs * form / fund) & SP_FT_PHMASK;

    ovp->forminc = (int32_t)(form * p->ftp1->sicvt);

    if (p->band != p->prvband) {
        p->prvband = p->band;
        p->expamp = exp(p->band * MPIDSR);
        newexp = 1;
    }
    /* Init grain rise ftable phase. Negative kform values make
    the kris (ifnb) initial index go negative and crash csound.
    So insert another if-test with compensating code. */
    if (p->ris >= (1.0 / sp->sr) && form != 0.0) {
        if (form < 0.0 && ovp->formphs != 0)
            ovp->risphs = (int32_t)((SP_FT_MAXLEN - ovp->formphs) / -form / p->ris);
        else ovp->risphs = (int32_t)(ovp->formphs / form / p->ris);

        ovp->risinc = (int32_t)(p->ftp1->sicvt / p->ris);
        rismps = SP_FT_MAXLEN / ovp->risinc;
    } else {
        ovp->risphs = SP_FT_MAXLEN;
        rismps = 0;
    }

    if (newexp || rismps != p->prvsmps) {
        if ((p->prvsmps = rismps))
        p->preamp = intpow(p->expamp, -rismps);
        else p->preamp = 1.0;
    }

    ovp->curamp = octamp * p->preamp;
    ovp->expamp = p->expamp;
    if ((ovp->dectim = (int32_t)(p->dec * sp->sr)) > 0)
        ovp->decinc = (int32_t)(p->ftp1->sicvt / p->dec);
    ovp->decphs = SP_FT_PHMASK;
    return 1;
}

int sp_fof_create(sp_fof **p)
{
    *p = malloc(sizeof(sp_fof));
    return SP_OK;
}

int sp_fof_destroy(sp_fof **p)
{
    sp_fof *pp = *p;
    sp_auxdata_free(&pp->auxch);
    free(*p);
    return SP_OK;
}

int sp_fof_init(sp_data *sp, sp_fof *p, sp_ftbl *sine, sp_ftbl *win, int iolaps, SPFLOAT iphs)
{
    p->amp = 0.5;
    p->fund = 100;
    p->form = 500;
    p->oct = 0;
    p->band = 50;
    p->ris = 0.003;
    p->dec = 0.0007;
    p->dur = 0.02;
    p->iolaps = iolaps;
    p->iphs = iphs;
    p->ftp1 = sine;
    p->ftp2 = win;

    sp_fof_overlap *ovp, *nxtovp;
    int32_t olaps;

    if (p->iphs == 0.0) p->fundphs = SP_FT_MAXLEN;                  
    else p->fundphs = (int32_t)(p->iphs * SP_FT_MAXLEN) & SP_FT_PHMASK;

    olaps = (int32_t)p->iolaps;

    if (p->iphs >= 0.0) {
        sp_auxdata_alloc(&p->auxch, (size_t)olaps * sizeof(sp_fof_overlap));
    }

    ovp = &p->basovrlap;
    nxtovp = (sp_fof_overlap *) p->auxch.ptr;

    do {
        ovp->nxtact = NULL;
        ovp->nxtfree = nxtovp;
        ovp = nxtovp++;
    } while (--olaps);
    ovp->nxtact = NULL;
    ovp->nxtfree = NULL;
    p->fofcount = -1;
    p->prvband = 0.0;
    p->expamp = 1.0;
    p->prvsmps = (int32_t)0;
    p->preamp = 1.0;

    p->ampcod   = 1;
    p->fundcod  = 1;
    p->formcod  = 1;
    p->xincod   = p->ampcod || p->fundcod || p->formcod;
    p->fmtmod = 0;
    p->foftype = 1;
    return SP_OK;
}

int sp_fof_compute(sp_data *sp, sp_fof *p, SPFLOAT *in, SPFLOAT *out)
{
    sp_fof_overlap *ovp;
    sp_ftbl *ftp1, *ftp2;
    SPFLOAT amp, fund, form;
    int32_t fund_inc, form_inc;
    SPFLOAT v1, fract ,*ftab;

    amp = p->amp;
    fund = p->fund;
    form = p->form;
    ftp1 = p->ftp1;
    ftp2 = p->ftp2;
    fund_inc = (int32_t)(fund * ftp1->sicvt);
    form_inc = (int32_t)(form * ftp1->sicvt);

    if (p->fundphs & SP_FT_MAXLEN) {
        p->fundphs &= SP_FT_PHMASK;
        ovp = p->basovrlap.nxtfree;
        if (newpulse(sp, p, ovp, amp, fund, form)) {
            ovp->nxtact = p->basovrlap.nxtact;
            p->basovrlap.nxtact = ovp;
            p->basovrlap.nxtfree = ovp->nxtfree;
        }
    }
    *out = 0.0;
    ovp = &p->basovrlap;
    while (ovp->nxtact != NULL) {
        SPFLOAT  result;
        sp_fof_overlap *prvact = ovp;
        ovp = ovp->nxtact;
        fract = PFRAC1(ovp->formphs);
        ftab = ftp1->tbl + (ovp->formphs >> ftp1->lobits);
        v1 = *ftab++;
        result = v1 + (*ftab - v1) * fract;
        if (p->foftype) {
            if (p->fmtmod)
            ovp->formphs += form_inc;           
            else ovp->formphs += ovp->forminc;
        }
        else {
            /* SPFLOAT ovp->glissbas = kgliss / grain length. ovp->sampct is
             incremented each sample. We add glissbas * sampct to the
             pitch of grain at each a-rate pass (ovp->formphs is the
             index into ifna; ovp->forminc is the stepping factor that
             decides pitch) */
            ovp->formphs += (int32_t)(ovp->forminc + ovp->glissbas * ovp->sampct++);
        }
        ovp->formphs &= SP_FT_PHMASK;
        if (ovp->risphs < SP_FT_MAXLEN) {
            result *= *(ftp2->tbl + (ovp->risphs >> ftp2->lobits) );
            ovp->risphs += ovp->risinc;
        }
        if (ovp->timrem <= ovp->dectim) {
            result *= *(ftp2->tbl + (ovp->decphs >> ftp2->lobits) );
            if ((ovp->decphs -= ovp->decinc) < 0) ovp->decphs = 0;
        }
        *out += (result * ovp->curamp);
        if (--ovp->timrem) ovp->curamp *= ovp->expamp;
        else {
            prvact->nxtact = ovp->nxtact;
            ovp->nxtfree = p->basovrlap.nxtfree;
            p->basovrlap.nxtfree = ovp;
            ovp = prvact;
        }
    }
    p->fundphs += fund_inc;
    return SP_OK;
}



int sp_fofilt_create(sp_fofilt **p)
{
    *p = malloc(sizeof(sp_fofilt));
    return SP_OK;
}

int sp_fofilt_destroy(sp_fofilt **p)
{
    free(*p);
    return SP_OK;
}

int sp_fofilt_init(sp_data *sp, sp_fofilt *p)
{
   p->tpidsr = 2.0*M_PI / sp->sr;
   p->sr = sp->sr;

   p->freq = 1000;
   p->atk = 0.007;
   p->dec = 0.04;
   p->istor = 0.0;

   int i;
   if (p->istor==0.0){
        for (i=0; i<4; i++)
         p->delay[i] = 0.0;
   }
   return SP_OK;
}

int sp_fofilt_compute(sp_data *sp, sp_fofilt *p, SPFLOAT *in, SPFLOAT *out)
{

    SPFLOAT freq = p->freq;
    SPFLOAT ris = p->atk;
    SPFLOAT dec = p->dec;
    SPFLOAT *delay = p->delay,ang=0,fsc,rrad1=0,rrad2=0;
    SPFLOAT w1,y1,w2,y2;
    SPFLOAT lfrq = -1.0, lrs = -1.0, ldc = -1.0;

    SPFLOAT frq = freq;
    SPFLOAT rs = ris;
    SPFLOAT dc = dec;
    if (frq != lfrq || rs != lrs || dc != ldc) {
        lfrq = frq; lrs = rs; ldc = dc;
        ang = (SPFLOAT)p->tpidsr*frq;
        fsc = sin(ang) - 3.0;

        rrad1 =  pow(10.0, fsc/(dc*sp->sr));
        rrad2 =  pow(10.0, fsc/(rs*sp->sr));
    }

    w1  = *in + 2.0*rrad1*cos(ang)*delay[0] - rrad1*rrad1*delay[1];
    y1 =  w1 - delay[1];
    delay[1] = delay[0];
    delay[0] = w1;

    w2  = *in + 2.0*rrad2*cos(ang)*delay[2] - rrad2*rrad2*delay[3];
    y2 =  w2 - delay[3];
    delay[3] = delay[2];
    delay[2] = w2;

    *out = (SPFLOAT) (y1 - y2);

    return SP_OK;
}



static int newpulse2(sp_data *sp, sp_fog *p, sp_fog_overlap *ovp, SPFLOAT amp,
                    SPFLOAT fund, SPFLOAT ptch)
{
    SPFLOAT octamp = amp, oct;
    SPFLOAT form = ptch / p->ftp1->sicvt, fogcvt = p->fogcvt;
    int32_t rismps, newexp = 0;
    ovp->timrem = (int32_t)(p->dur * sp->sr);

    if ((oct = p->oct) > 0.0) {
        int32_t ioct = (int32_t)oct, bitpat = ~(-1L << ioct);
        if (bitpat & ++p->fofcount) return(0);
        if ((bitpat += 1) & p->fofcount) octamp *= (1.0) + ioct - oct;
    }

    if (fund == 0.0) ovp->formphs = 0;
    else ovp->formphs = (int32_t)(p->fundphs * form / fund) & SP_FT_PHMASK;

    ovp->forminc = (int32_t)(ptch * fogcvt);

    if (p->band != p->prvband) {
        p->prvband = p->band;
        p->expamp = exp(p->band * MPIDSR);
        newexp = 1;
    }

    if (p->ris >= (1.0 / sp->sr)  && form != 0.0) {
        ovp->risphs = (uint32_t)(ovp->formphs / (fabs(form))
                                    / p->ris);
        ovp->risinc = (int32_t)(p->ftp1->sicvt / p->ris);
        rismps = SP_FT_MAXLEN / ovp->risinc;
    } else {
        ovp->risphs = SP_FT_MAXLEN;
        rismps = 0;
    }
    ovp->formphs = (ovp->formphs + p->spdphs) & SP_FT_PHMASK;

    if (newexp || rismps != p->prvsmps) {
        if ((p->prvsmps = rismps)) p->preamp = intpow(p->expamp, -rismps);
        else p->preamp = 1.0;
    }

    ovp->curamp = octamp * p->preamp;
    ovp->expamp = p->expamp;

    if ((ovp->dectim = (int32_t)(p->dec * sp->sr )) > 0) {
        ovp->decinc = (int32_t)(p->ftp1->sicvt / p->dec);
    }

    ovp->decphs = SP_FT_PHMASK;

    ovp->pos = p->spd * p->ftp1->size;
    ovp->inc = p->trans;

    return 1;
}

int sp_fog_create(sp_fog **p)
{
    *p = malloc(sizeof(sp_fog));
    return SP_OK;
}

int sp_fog_destroy(sp_fog **p)
{
    sp_fog *pp = *p;
    sp_auxdata_free(&pp->auxch);
    free(*p);
    return SP_OK;
}

int sp_fog_init(sp_data *sp, sp_fog *p, sp_ftbl *wav, sp_ftbl *win, int iolaps, SPFLOAT iphs)
{
    p->amp = 0.5;
    p->dens = 80;
    p->trans = 1;
    p->spd = 0;
    p->oct = 0;
    p->band = 50;
    p->ris = 0.01;
    p->dec = 0.07;
    p->dur = 0.1;
    p->iolaps = iolaps;
    p->iphs = iphs;
    p->ftp1 = wav;
    p->ftp2 = win;

    sp_fog_overlap *ovp, *nxtovp;
    int32_t olaps;
    p->fogcvt = SP_FT_MAXLEN/(p->ftp1)->size;
    p->spdphs = 0L;
    if (p->iphs == 0.0) p->fundphs = SP_FT_MAXLEN;
    else p->fundphs = (int32_t)(p->iphs * SP_FT_MAXLEN) & SP_FT_PHMASK;

    olaps = (int32_t)p->iolaps;

    sp_auxdata_alloc(&p->auxch, (size_t)olaps * sizeof(sp_fog_overlap));
    ovp = &p->basovrlap;
    nxtovp = (sp_fog_overlap *) p->auxch.ptr;

    do {
        ovp->nxtact = NULL;
        ovp->nxtfree = nxtovp;
        ovp = nxtovp++;
    } while (--olaps);

    ovp->nxtact  = NULL;
    ovp->nxtfree = NULL;
    p->fofcount = -1;
    p->prvband = 0.0;
    p->expamp = 1.0;
    p->prvsmps = 0;
    p->preamp = 1.0;
    p->fmtmod  = 0;
    return SP_OK;
}

int sp_fog_compute(sp_data *sp, sp_fog *p, SPFLOAT *in, SPFLOAT *out)
{
    sp_fog_overlap *ovp;
    sp_ftbl *ftp1,  *ftp2;
    SPFLOAT  amp, fund, ptch, speed;
    SPFLOAT fract;
    int32_t fund_inc;
    int32_t incr;

    int32_t ndx;
    SPFLOAT x1, x2;


    amp = p->amp;
    fund = p->dens;
    ptch = p->trans;
    speed = p->spd;
    ftp1 = p->ftp1;
    ftp2 = p->ftp2;
    fund_inc = (int32_t)(fund * ftp1->sicvt);

    if (p->fundphs & SP_FT_MAXLEN) {
        p->fundphs &= SP_FT_PHMASK;
        ovp = p->basovrlap.nxtfree;
        if (newpulse2(sp, p, ovp, amp, fund, ptch)) {
            ovp->nxtact = p->basovrlap.nxtact;
            p->basovrlap.nxtact = ovp;
            p->basovrlap.nxtfree = ovp->nxtfree;
        }
    }
    *out = 0.0;
    ovp = &p->basovrlap;
    while (ovp->nxtact != NULL) {
        SPFLOAT result;
        sp_fog_overlap *prvact = ovp;
        ovp = ovp->nxtact;
        ndx = floor(ovp->pos);
        fract = ovp->pos - ndx;

        while(ndx >= ftp1->size) {
            ndx -= ftp1->size;
        }

        while(ndx < 0) ndx += ftp1->size;

        x1 = ftp1->tbl[ndx];
        x2 = ftp1->tbl[ndx + 1];

        result = x1 + (x2 - x1) * fract;

        ovp->pos += ovp->inc;

        if (ovp->risphs < SP_FT_MAXLEN) {
            /* bounds checking so it doesn't segfault */
            incr = (ovp->risphs >> ftp2->lobits);
            if(incr <= ftp2->size) {
                result *= *(ftp2->tbl + incr );
            } else {
                result = 0;
            }
            ovp->risphs += ovp->risinc;
        }
        if (ovp->timrem <= ovp->dectim) {
            incr = (ovp->decphs >> ftp2->lobits);
            if(incr <= ftp2->size) {
                result *= *(ftp2->tbl + incr );
            } else {
                result = 0;
            }
            if ((ovp->decphs -= ovp->decinc) < 0)
            ovp->decphs = 0;
        }
        *out += (result * ovp->curamp);
        if (--ovp->timrem) ovp->curamp *= ovp->expamp;
        else {
            prvact->nxtact = ovp->nxtact;
            ovp->nxtfree = p->basovrlap.nxtfree;
            p->basovrlap.nxtfree = ovp;
            ovp = prvact;
        }
    }

    p->fundphs += fund_inc;
    p->spdphs = (int32_t)(speed * SP_FT_MAXLEN);
    p->spdphs &= SP_FT_PHMASK;
    return SP_OK;
}

int sp_fold_create(sp_fold **p)
{
    *p = malloc(sizeof(sp_fold));
    return SP_OK;
}

int sp_fold_destroy(sp_fold **p)
{
    free(*p);
    return SP_OK;
}

int sp_fold_init(sp_data *sp, sp_fold *p)
{
    p->incr = 1000;
    p->sample_index = 0;
    p->index = 0.0;
    p->value = 0.0; 
    return SP_OK;
}

int sp_fold_compute(sp_data *sp, sp_fold *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT index = p->index;
    int32_t sample_index = p->sample_index;
    SPFLOAT value = p->value;
    if (index < (SPFLOAT)sample_index) {
        index += p->incr;
        *out = value = *in;
    } else {
        *out = value;
    }
    sample_index++;
    p->index = index;
    p->sample_index = sample_index;
    p->value = value;
    return SP_OK;
}


int sp_foo_create(sp_foo **p)
{
    *p = malloc(sizeof(sp_foo));
    return SP_OK;
}

int sp_foo_destroy(sp_foo **p)
{
    free(*p);
    return SP_OK;
}

int sp_foo_init(sp_data *sp, sp_foo *p)
{
    /* Initalize variables here. */
    p->bar = 123;
    return SP_OK;
}

int sp_foo_compute(sp_data *sp, sp_foo *p, SPFLOAT *in, SPFLOAT *out)
{
    /* Send the signal's input to the output */
    *out = *in;
    return SP_OK;
}


int sp_fosc_create(sp_fosc **p)
{
    *p = malloc(sizeof(sp_fosc));
    return SP_OK;
}

int sp_fosc_destroy(sp_fosc **p)
{
    free(*p);
    return SP_OK;
}

int sp_fosc_init(sp_data *sp, sp_fosc *p, sp_ftbl *ft)
{
    p->freq = 440;
    p->amp = 0.4;
    p->iphs = 0.0;
    p->ft = ft;

    p->mod = 1.0;
    p->car = 1.0;
    p->indx = 1.0;

    p->cphs = p->mphs = (int32_t)(p->iphs * SP_FT_MAXLEN);

    return SP_OK;
}

int sp_fosc_compute(sp_data *sp, sp_fosc *p, SPFLOAT *in, SPFLOAT *out)
{

    sp_ftbl *ftp;

    SPFLOAT  amp, cps, fract, v1, v2, car, fmod, cfreq, mod;
    SPFLOAT  xcar, xmod, ndx, *ftab;
    int32_t  mphs, cphs, minc, cinc, lobits;
    SPFLOAT  sicvt = p->ft->sicvt;
    SPFLOAT  *ft;

    ftp = p->ft;
    ft = ftp->tbl;
    lobits = ftp->lobits;
    mphs = p->mphs;
    cphs = p->cphs;
    cps  = p->freq;
    amp  = p->amp;
    xcar = p->car;
    xmod = p->mod;

    car = cps * xcar;
    mod = cps * xmod;
    ndx = p->indx * mod;
    minc = (int32_t)(mod * sicvt);
    mphs &= SP_FT_PHMASK;
    fract = ((mphs) & ftp->lomask) * ftp->lodiv;
    ftab = ft + (mphs >> lobits);
    v1 = ftab[0];

    if(ftab[0] == p->ft->tbl[p->ft->size - 1]) {
        v2 = p->ft->tbl[0];
    } else {
        v2 = ftab[1];
    }

    fmod = (v1 + (v2 - v1) * fract) * ndx;
    mphs += minc;
    cfreq = car + fmod;
    cinc = (int32_t)(cfreq * sicvt);
    cphs &= SP_FT_PHMASK;
    fract = ((cphs) & ftp->lomask) * ftp->lodiv;
    ftab = ft + (cphs >>lobits);
    v1 = ftab[0];

    if(ftab[0] == p->ft->tbl[p->ft->size - 1]) {
        v2 = p->ft->tbl[0];
    } else {
        v2 = ftab[1];
    }

    *out = (v1 + (v2 - v1) * fract) * amp;
    cphs += cinc;
    p->mphs = mphs;
    p->cphs = cphs;

    return SP_OK;
}



/* initialize constants in ftable */
int sp_ftbl_init(sp_data *sp, sp_ftbl *ft, size_t size)
{
    ft->size = size;
    ft->sicvt = 1.0 * SP_FT_MAXLEN / sp->sr;
    ft->lobits = log2(SP_FT_MAXLEN / size);
    ft->lomask = (1<<ft->lobits) - 1;
    ft->lodiv = 1.0 / (1<<ft->lobits);
    ft->del = 1;
    return SP_OK;
}

int sp_ftbl_create(sp_data *sp, sp_ftbl **ft, size_t size)
{
    *ft = malloc(sizeof(sp_ftbl));
    sp_ftbl *ftp = *ft;
    ftp->tbl = malloc(sizeof(SPFLOAT) * (size + 1));
    memset(ftp->tbl, 0, sizeof(SPFLOAT) * (size + 1));
   
    sp_ftbl_init(sp, ftp, size);
    return SP_OK;
}

int sp_ftbl_bind(sp_data *sp, sp_ftbl **ft, SPFLOAT *tbl, size_t size)
{
    *ft = malloc(sizeof(sp_ftbl));
    sp_ftbl *ftp = *ft;
    ftp->tbl = tbl;
    sp_ftbl_init(sp, ftp, size);
    ftp->del = 0;
    return SP_OK;
}

int sp_ftbl_destroy(sp_ftbl **ft)
{
    sp_ftbl *ftp = *ft;
    if(ftp->del) free(ftp->tbl);
    free(*ft);
    return SP_OK;
}

/* TODO: handle spaces at beginning of string */
static char * tokenize(char **next, int *size)
{
    if(*size <= 0) return NULL;
    char *token = *next;
    char *str = *next;

    char *peak = str + 1;

    while((*size)--) {
        if(*str == ' ') {
            *str = 0;
            if(*peak != ' ') break;
        }
        str = str + 1;
        peak = str + 1;
    }
    *next = peak;
    return token;
}

int sp_gen_vals(sp_data *sp, sp_ftbl *ft, const char *string)
{
    int size = strlen(string);
    char *str = malloc(sizeof(char) * size + 1);
    strcpy(str, string);
    char *out; 
    char *ptr = str;
    int j = 0;
    while(size > 0) {
        out = tokenize(&str, &size);
        if(ft->size < j + 1){
            ft->tbl = realloc(ft->tbl, sizeof(SPFLOAT) * (ft->size + 2));
            /* zero out new tables */
            ft->tbl[ft->size] = 0;
            ft->tbl[ft->size + 1] = 0;
            ft->size++;
        }
        ft->tbl[j] = atof(out);
        j++;
    }
  
    sp_ftbl_init(sp, ft, ft->size);
    free(ptr); 
    return SP_OK;
}

int sp_gen_sine(sp_data *sp, sp_ftbl *ft)
{
    unsigned long i;
    SPFLOAT step = 2 * M_PI / ft->size;
    for(i = 0; i < ft->size; i++){
        ft->tbl[i] = sin(i * step);
    }
    return SP_OK;
}

#ifndef NO_LIBSNDFILE
/*TODO: add error checking, make tests */
int sp_gen_file(sp_data *sp, sp_ftbl *ft, const char *filename)
{
    SF_INFO info;
    memset(&info, 0, sizeof(SF_INFO));
    info.format = 0;
    SNDFILE *snd = sf_open(filename, SFM_READ, &info);
#ifdef USE_DOUBLE
    sf_readf_double(snd, ft->tbl, ft->size);
#else
    sf_readf_float(snd, ft->tbl, ft->size);
#endif
    sf_close(snd);
    return SP_OK;
}

int sp_ftbl_loadfile(sp_data *sp, sp_ftbl **ft, const char *filename)
{
    *ft = malloc(sizeof(sp_ftbl));
    sp_ftbl *ftp = *ft;
    SF_INFO info;
    memset(&info, 0, sizeof(SF_INFO));
    info.format = 0;
    SNDFILE *snd = sf_open(filename, SFM_READ, &info);
    if(snd == NULL) {
        return SP_NOT_OK;
    }
    size_t size = info.frames * info.channels;

    ftp->tbl = malloc(sizeof(SPFLOAT) * (size + 1));

    sp_ftbl_init(sp, ftp, size);

#ifdef USE_DOUBLE
    sf_readf_double(snd, ftp->tbl, ftp->size);
#else
    sf_readf_float(snd, ftp->tbl, ftp->size);
#endif
    sf_close(snd);
    return SP_OK;
}
#endif

/* port of GEN10 from Csound */
int sp_gen_sinesum(sp_data *sp, sp_ftbl *ft, const char *argstring)
{
    sp_ftbl *args;
    sp_ftbl_create(sp, &args, 1);
    sp_gen_vals(sp, args, argstring);

    int32_t phs;
    SPFLOAT amp;
    int32_t flen = (int32_t)ft->size;
    SPFLOAT tpdlen = 2.0 * M_PI / (SPFLOAT) flen;

    int32_t i, n;

    for(i = (int32_t)args->size; i > 0; i--){
        amp = args->tbl[args->size - i];
        if(amp != 0) {
            for(phs = 0, n = 0; n < ft->size; n++){
                ft->tbl[n] += sin(phs * tpdlen) * amp;
                phs += i;
                phs %= flen;
            }
        }
    }
    sp_ftbl_destroy(&args);
    return SP_OK;
}

int sp_gen_line(sp_data *sp, sp_ftbl *ft, const char *argstring)
{
    uint16_t i, n = 0, seglen;
    SPFLOAT incr, amp = 0;
    SPFLOAT x1, x2, y1, y2;
    sp_ftbl *args;
    sp_ftbl_create(sp, &args, 1);
    sp_gen_vals(sp, args, argstring);

    if((args->size % 2) == 1 || args->size == 1) {
        fprintf(stderr, "Error: not enough arguments for gen_line.\n");
        sp_ftbl_destroy(&args);
        return SP_NOT_OK;
    } else if(args->size == 2) {
        for(i = 0; i < ft->size; i++) {
            ft->tbl[i] = args->tbl[1];
        }
        return SP_OK;
    }

    x1 = args->tbl[0];
    y1 = args->tbl[1];
    for(i = 2; i < args->size; i += 2) {
        x2 = args->tbl[i];
        y2 = args->tbl[i + 1];

        if(x2 < x1) {
            fprintf(stderr, "Error: x coordiates must be sequential!\n");
            break;
        }

        seglen = (x2 - x1);
        incr = (SPFLOAT)(y2 - y1) / (seglen - 1);
        amp = y1;

        while(seglen != 0){
            if(n < ft->size) {
                ft->tbl[n] = amp;
                amp += incr;
                seglen--;
                n++;
            } else {
                break;
            }
        }
        y1 = y2;
        x1 = x2;
    }

    sp_ftbl_destroy(&args);
    return SP_OK;
}

int sp_gen_xline(sp_data *sp, sp_ftbl *ft, const char *argstring)
{
    uint16_t i, n = 0, seglen;
    SPFLOAT mult, amp = 0;
    SPFLOAT x1, x2, y1, y2;
    sp_ftbl *args;
    sp_ftbl_create(sp, &args, 1);
    sp_gen_vals(sp, args, argstring);

    if((args->size % 2) == 1 || args->size == 1) {
        fprintf(stderr, "Error: not enough arguments for gen_line.\n");
        sp_ftbl_destroy(&args);
        return SP_NOT_OK;
    } else if(args->size == 2) {
        for(i = 0; i < ft->size; i++) {
            ft->tbl[i] = args->tbl[1];
        }
        return SP_OK;
    }

    x1 = args->tbl[0];
    y1 = args->tbl[1];
    for(i = 2; i < args->size; i += 2) {
        x2 = args->tbl[i];
        y2 = args->tbl[i + 1];

        if(x2 < x1) {
            fprintf(stderr, "Error: x coordiates must be sequential!\n");
            break;
        }

        if(y1 == 0) {
            y1 = 0.000001;
        }

        if(y2 == 0) {
            y2 = 0.000001;
        }

        seglen = (uint32_t)(x2 - x1);
        mult = (y2 / y1);
        mult = pow(mult, (SPFLOAT)1.0 / seglen);
        amp = y1;

        while(seglen != 0){
            if(n < ft->size) {
                ft->tbl[n] = amp;
                amp *= mult;
                seglen--;
                n++;
            } else {
                break;
            }
        }
        y1 = y2;
        x1 = x2;
    }

    sp_ftbl_destroy(&args);
    return SP_OK;

}


static SPFLOAT gaussrand(sp_randmt *p, SPFLOAT scale)
{
    int64_t r1 = -((int64_t)0xFFFFFFFFU * 6);
    int n = 12;
    SPFLOAT x;

    do {
      r1 += (int64_t)sp_randmt_compute(p);
    } while (--n);

    x = (SPFLOAT)r1;
    return (SPFLOAT)(x * ((SPFLOAT)scale * (1.0 / (3.83 * 4294967295.03125))));
}

int sp_gen_gauss(sp_data *sp, sp_ftbl *ft, SPFLOAT scale, uint32_t seed)
{
    int n;

    sp_randmt rand;

    sp_randmt_seed(&rand, NULL, seed);

    for(n = 0; n < ft->size; n++) {
        ft->tbl[n] = gaussrand(&rand, scale);
    }

    return SP_OK;
}

/* based off of GEN 19 */
int sp_gen_composite(sp_data *sp, sp_ftbl *ft, const char *argstring)
{
    SPFLOAT phs, inc, amp, dc, tpdlen = 2 * M_PI/ (SPFLOAT) ft->size;
    int i, n;
    
    sp_ftbl *args;
    sp_ftbl_create(sp, &args, 1);
    sp_gen_vals(sp, args, argstring);

    for(n = 0; n < args->size; n += 4) {
        inc = args->tbl[n] * tpdlen;
        amp = args->tbl[n + 1];
        phs = args->tbl[n + 2] * tpd360;
        dc = args->tbl[n + 3];

        for (i = 0; i <ft->size ; i++) {
            ft->tbl[i] += (SPFLOAT) (sin(phs) * amp + dc);
            if ((phs += inc) >= 2 * M_PI) phs -= 2 * M_PI;
        }
    }

    sp_ftbl_destroy(&args);
    return SP_OK;
}

int sp_gen_rand(sp_data *sp, sp_ftbl *ft, const char *argstring)
{
    sp_ftbl *args;
    sp_ftbl_create(sp, &args, 1);
    sp_gen_vals(sp, args, argstring);
    int n, pos = 0, i, size = 0;

    for(n = 0; n < args->size; n += 2) {
        size = round(ft->size * args->tbl[n + 1]);
        for(i = 0; i < size; i++) {
            if(pos < ft->size) {
                ft->tbl[pos] = args->tbl[n];
                pos++;
            }
        }
    }
    if(pos <= ft->size) {
        ft->size = pos;
    }
    sp_ftbl_destroy(&args);
    return SP_OK;
}

int sp_gen_triangle(sp_data *sp, sp_ftbl *ft)
{
    unsigned int i;
    unsigned int counter;
    SPFLOAT incr;
    int step;

    incr = 1.0f / (SPFLOAT)ft->size;
    incr *= 2;

    step = 1;

    counter = 0;

    for(i = 0; i < ft->size; i++) {
        if(i == ft->size / 2) {
            step = -1;
        }
        ft->tbl[i] = (2.f*(counter * incr) - 1.f);

        counter += step;
    }

    return SP_OK;
}




int sp_gbuzz_create(sp_gbuzz **p)
{
    *p = malloc(sizeof(sp_gbuzz));
    return SP_OK;
}

int sp_gbuzz_destroy(sp_gbuzz **p)
{
    free(*p);
    return SP_OK;
}

int sp_gbuzz_init(sp_data *sp, sp_gbuzz *p, sp_ftbl *ft, SPFLOAT iphs)
{
    p->freq = 440;
    p->amp = 0.4;
    p->nharm = 4;
    p->lharm = 1;
    p->mul = 0.1;
    p->ft = ft;
    p->iphs = iphs; 
    
    if (p->iphs >= 0) {
        p->lphs = (int32_t)(p->iphs * SP_FT_MAXLEN);
        p->prvr = 0.0;
    }
    p->last = 1.0;
    return SP_OK;
}

int sp_gbuzz_compute(sp_data *sp, sp_gbuzz *p, SPFLOAT *in, SPFLOAT *out)
{
    sp_ftbl *ftp;
    SPFLOAT *ftbl;
    int32_t phs, inc, lobits, lenmask, k, km1, kpn, kpnm1;
    SPFLOAT r, absr, num, denom, scal, last = p->last;
    int32_t nn, lphs = p->lphs;
    
    ftp = p->ft;
    ftbl = ftp->tbl;
    lobits = ftp->lobits;
    lenmask = (int32_t) ftp->size - 1;
    k = (int32_t)p->lharm;
    
    if ((nn = (int32_t)p->nharm)<0) nn = -nn;
    
    if (nn == 0) {
        nn = 1;
    }
    km1 = k - 1;
    kpn = k + nn;
    kpnm1 = kpn - 1;
    
    if ((r = p->mul) != p->prvr || nn != p->prvn) {
        p->twor = r + r;
        p->rsqp1 = r * r + 1.0;
        p->rtn = intpow1(r, nn);
        p->rtnp1 = p->rtn * r;
        
        if ((absr = fabs(r)) > 0.999 && absr < 1.001) {
            p->rsumr = 1.0 / nn;
        } else {
            p->rsumr = (1.0 - absr) / (1.0 - fabs(p->rtn));
        }
        
        p->prvr = r;
        p->prvn = (int16_t)nn;
    }
    
    scal =  p->amp * p->rsumr;
    inc = (int32_t)(p->freq * ftp->sicvt);
    phs = lphs >>lobits;
    denom = p->rsqp1 - p->twor * ftbl[phs];
    num = ftbl[phs * k & lenmask]
        - r * ftbl[phs * km1 & lenmask]
        - p->rtn * ftbl[phs * kpn & lenmask]
        + p->rtnp1 * ftbl[phs * kpnm1 & lenmask];
    
    if (denom > 0.0002 || denom < -0.0002) {
        *out = last = num / denom * scal;
    } else if (last<0) {
        *out = last = - *out;
    } else {
        *out = last = *out;
    }
    
    lphs += inc;
    lphs &= SP_FT_PHMASK;
    p->last = last;
    p->lphs = lphs;
    return SP_OK;
}


int sp_hilbert_create(sp_hilbert **p)
{
    *p = malloc(sizeof(sp_hilbert));
    return SP_OK;
}

int sp_hilbert_destroy(sp_hilbert **p)
{
    free(*p);
    return SP_OK;
}

int sp_hilbert_init(sp_data *sp, sp_hilbert *p)
{
    int j; 
    SPFLOAT onedsr = 1.0 / sp->sr;
    /* pole values taken from Bernie Hutchins, "Musical Engineer's Handbook" */
    SPFLOAT poles[12] = {0.3609, 2.7412, 11.1573, 44.7581, 179.6242, 798.4578,
                        1.2524, 5.5671, 22.3423, 89.6271, 364.7914, 2770.1114};
    SPFLOAT polefreq, rc, alpha, beta;
    /* calculate coefficients for allpass filters, based on sampling rate */
    for (j=0; j<12; j++) {
        polefreq = poles[j] * 15.0;
        rc = 1.0 / (2.0 * M_PI * polefreq);
        alpha = 1.0 / rc;
        alpha = alpha * 0.5 * onedsr;
        beta = (1.0 - alpha) / (1.0 + alpha);
        p->xnm1[j] = p->ynm1[j] = 0.0;
        p->coef[j] = -(SPFLOAT)beta;
    }
    return SP_OK;
}

int sp_hilbert_compute(sp_data *sp, sp_hilbert *p, SPFLOAT *in, SPFLOAT *out1, SPFLOAT *out2)
{
    SPFLOAT xn1 = 0;
    SPFLOAT yn1 = 0; 
    SPFLOAT xn2 = 0;
    SPFLOAT yn2 = 0;
    SPFLOAT *coef;
    int j;

    coef = p->coef;

    xn1 = *in;
    /* 6th order allpass filter for sine output. Structure is
    * 6 first-order allpass sections in series. Coefficients
    * taken from arrays calculated at i-time.
    */
    for (j=0; j < 6; j++) {
        yn1 = coef[j] * (xn1 - p->ynm1[j]) + p->xnm1[j];
        p->xnm1[j] = xn1;
        p->ynm1[j] = yn1;
        xn1 = yn1;
    }
    xn2 = *in;
    /* 6th order allpass filter for cosine output. Structure is
    * 6 first-order allpass sections in series. Coefficients
    * taken from arrays calculated at i-time.
    */
    for (j=6; j < 12; j++) {
        yn2 = coef[j] * (xn2 - p->ynm1[j]) + p->xnm1[j];
        p->xnm1[j] = xn2;
        p->ynm1[j] = yn2;
        xn2 = yn2;
    }
    *out1 = yn2;
    *out2 = yn1;
    return SP_OK;
}


int sp_in_create(sp_in **p)
{
    *p = malloc(sizeof(sp_in));
    return SP_OK;
}

int sp_in_destroy(sp_in **p)
{
    sp_in *pp = *p;
    fclose(pp->fp);
    free(*p);
    return SP_OK;
}

int sp_in_init(sp_data *sp, sp_in *p)
{
    p->fp = stdin; 
    return SP_OK;
}

int sp_in_compute(sp_data *sp, sp_in *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0;
    fread(out, sizeof(SPFLOAT), 1, p->fp);
    return SP_OK;
}

int sp_incr_create(sp_incr **p)
{
    *p = malloc(sizeof(sp_incr));
    return SP_OK;
}

int sp_incr_destroy(sp_incr **p)
{
    free(*p);
    return SP_OK;
}

int sp_incr_init(sp_data *sp, sp_incr *p, SPFLOAT val)
{
    p->min = 0;
    p->max = 1;
    p->step = 0.1;
    p->val = val;
    return SP_OK;
}

int sp_incr_compute(sp_data *sp, sp_incr *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in > 0 ) {
        p->val += p->step;
        p->val = max(min(p->val, p->max), p->min);
    } else if (*in < 0) {
        p->val -= p->step;
        p->val = max(min(p->val, p->max), p->min);
    }
    *out = p->val;
    return SP_OK;
}


/* the randgabs are essentially magic incantations from Csound */

static SPFLOAT sp_jitter_randgab(sp_data *sp) 
{
    SPFLOAT out = (SPFLOAT) ((sp_rand(sp) >> 1) & 0x7fffffff) *
    (4.656612875245796924105750827168e-10);
    return out;
}

static SPFLOAT sp_jitter_birandgab(sp_data *sp) 
{
    SPFLOAT out = (SPFLOAT) (sp_rand(sp) & 0x7fffffff) *
    (4.656612875245796924105750827168e-10);
    return out;
}

int sp_jitter_create(sp_jitter **p)
{
    *p = malloc(sizeof(sp_jitter));
    return SP_OK;
}

int sp_jitter_destroy(sp_jitter **p)
{
    free(*p);
    return SP_OK;
}

int sp_jitter_init(sp_data *sp, sp_jitter *p)
{
    p->amp = 0.5;
    p->cpsMin = 0.5;
    p->cpsMax = 4; 
    p->num2 = sp_jitter_birandgab(sp);
    p->initflag = 1;
    p->phs=0;
    return SP_OK;
}

int sp_jitter_compute(sp_data *sp, sp_jitter *p, SPFLOAT *in, SPFLOAT *out)
{
    if (p->initflag) {
      p->initflag = 0;
      *out = p->num2 * p->amp;
      p->cps = sp_jitter_randgab(sp) * (p->cpsMax - p->cpsMin) + p->cpsMin;
      p->phs &= SP_FT_PHMASK;
      p->num1 = p->num2;
      p->num2 = sp_jitter_birandgab(sp);
      p->dfdmax = 1.0 * (p->num2 - p->num1) / SP_FT_MAXLEN;
      return SP_OK;
    }
    
    *out = (p->num1 + (SPFLOAT)p->phs * p->dfdmax) * p->amp;
    p->phs += (int32_t)(p->cps * (SPFLOAT)(SP_FT_MAXLEN / sp->sr));

    if (p->phs >= SP_FT_MAXLEN) {
      p->cps   = sp_jitter_randgab(sp) * (p->cpsMax - p->cpsMin) + p->cpsMin;
      p->phs   &= SP_FT_PHMASK;
      p->num1   = p->num2;
      p->num2 =  sp_jitter_birandgab(sp);
      p->dfdmax = 1.0 * (p->num2 - p->num1) / SP_FT_MAXLEN;
    }
    return SP_OK;
}


int sp_line_create(sp_line **p)
{
    *p = malloc(sizeof(sp_line));
    return SP_OK;
}

int sp_line_destroy(sp_line **p)
{
    free(*p);
    return SP_OK;
}

static void line_reinit(sp_data *sp, sp_line *p)
{
    SPFLOAT onedsr = 1.0 / sp->sr;
    p->incr = (SPFLOAT)((p->b - p->a) / (p->dur)) * onedsr;
    p->val = p->a;
    p->stime = 0;
    p->sdur = sp->sr * p->dur;
}

int sp_line_init(sp_data *sp, sp_line *p)
{
    p->a = 0;
    p->dur = 0.5;
    p->b = 1;
    line_reinit(sp, p);
    p->init = 1;
    return SP_OK;
}

int sp_line_compute(sp_data *sp, sp_line *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in != 0 ) {
        line_reinit(sp, p);
        p->init = 0;
    }

    if(p->init) {
        *out = 0;
        return SP_OK;
    }

    if(p->stime < p->sdur) {
        SPFLOAT val = p->val;
        p->val += p->incr;
        p->stime++;
        *out = val;
    } else {
        *out = p->b;
    }

    return SP_OK;
}



typedef struct openlpc_e_state openlpc_encoder_state;
typedef struct openlpc_d_state openlpc_decoder_state;

openlpc_encoder_state *create_openlpc_encoder_state(void);
void init_openlpc_encoder_state(openlpc_encoder_state *st, int framelen);
int  openlpc_encode(const short *in, unsigned char *out, openlpc_encoder_state *st);
void destroy_openlpc_encoder_state(openlpc_encoder_state *st);

openlpc_decoder_state *create_openlpc_decoder_state(void);
void init_openlpc_decoder_state(openlpc_decoder_state *st, int framelen);
int  openlpc_decode(sp_data *sp, unsigned char *in, short *out, openlpc_decoder_state *st);
void destroy_openlpc_decoder_state(openlpc_decoder_state *st);

void openlpc_sr(float sr);

size_t openlpc_get_encoder_state_size(void);
size_t openlpc_get_decoder_state_size(void);

int sp_lpc_create(sp_lpc **lpc)
{
    *lpc = malloc(sizeof(sp_lpc));
    return SP_OK;
}

int sp_lpc_destroy(sp_lpc **lpc)
{
    sp_lpc *plpc;
    plpc = *lpc;
    sp_auxdata_free(&plpc->m_e);
    sp_auxdata_free(&plpc->m_d);
    sp_auxdata_free(&plpc->m_out);
    sp_auxdata_free(&plpc->m_in);
    free(*lpc);
    return SP_OK;
}

int sp_lpc_init(sp_data *sp, sp_lpc *lpc, int framesize)
{
    int i;
    lpc->counter = 0;
    lpc->clock = 0;
    lpc->block = 4;
    lpc->samp = 0;
    lpc->mode = 0;
    lpc->framesize = framesize;
    openlpc_sr(sp->sr / lpc->block);

    sp_auxdata_alloc(&lpc->m_d, openlpc_get_decoder_state_size());
    sp_auxdata_alloc(&lpc->m_e, openlpc_get_encoder_state_size());
    lpc->d = lpc->m_d.ptr;
    lpc->e = lpc->m_e.ptr;

    sp_auxdata_alloc(&lpc->m_in, sizeof(short) * framesize);
    sp_auxdata_alloc(&lpc->m_out, sizeof(short) * framesize);

    lpc->out = lpc->m_out.ptr;
    lpc->in = lpc->m_in.ptr;

    init_openlpc_decoder_state(lpc->d, framesize);
    init_openlpc_encoder_state(lpc->e, framesize);

    for(i = 0; i < framesize; i++) {
        lpc->in[i] = 0;
        lpc->out[i] = 0;
        if(i < 7) lpc->data[i] = 0;
    }
    return SP_OK;
}

int sp_lpc_compute(sp_data *sp, sp_lpc *lpc, SPFLOAT *in, SPFLOAT *out)
{
    int i;


    if(lpc->clock == 0) {
        if(lpc->counter == 0) {
            if(lpc->mode == 0) { 
                openlpc_encode(lpc->in, lpc->data, lpc->e);
            } else {
                for(i = 0; i < 7; i++) {
                    lpc->y[i] = 
                        lpc->smooth*lpc->y[i] + 
                        (1-lpc->smooth)*lpc->ft->tbl[i];
                    lpc->data[i] = 255 * lpc->y[i];
                }
            }
            openlpc_decode(sp, lpc->data, lpc->out, lpc->d);
        }

        if(lpc->mode == 0) lpc->in[lpc->counter] = *in * 32767; 
        lpc->samp = lpc->out[lpc->counter] / 32767.0;

        lpc->counter = (lpc->counter + 1) % lpc->framesize;
    }


    lpc->clock = (lpc->clock + 1) % lpc->block;
    *out = lpc->samp;

    return SP_OK;
}

int sp_lpc_synth(sp_data *sp, sp_lpc *lpc, sp_ftbl *ft)
{
    int i;
    int sr;
    sr = sp->sr;

    sr = sr / 4;
    sr = sr / lpc->framesize;
    lpc->ft = ft;
    lpc->mode = 1;
    for(i = 0; i < 7; i++) lpc->y[i] = 0;
    lpc->smooth = exp(-1.0 / (0.01 * sr));
    return SP_OK;
}


static float my_fs = 11025.0;

typedef struct openlpc_e_state{
	float   s[MAXWINDOW], y[MAXWINDOW], h[MAXWINDOW];
    int     framelen, buflen;
    float   xv1[3], yv1[3], 
            xv2[2], yv2[2], 
			xv3[1], yv3[3], 
			xv4[2], yv4[2];
    float   w[MAXWINDOW], r[LPC_FILTORDER+1];
} openlpc_e_state_t;

typedef struct openlpc_d_state{
		float Oldper, OldG, Oldk[LPC_FILTORDER + 1];
        float bp[LPC_FILTORDER+1];
        float exc;
		int pitchctr, framelen, buflen;
} openlpc_d_state_t;


#if BITS_FOR_LPC == 38
/* (38 bit LPC-10, 2.7 Kbit/s @ 20ms, 2.4 Kbit/s @ 22.5 ms */
static int parambits[LPC_FILTORDER] = {6,5,5,4,4,3,3,3,3,2};
#elif BITS_FOR_LPC == 32
/* (32 bit LPC-10, 2.4 Kbit/s, not so good */
static int parambits[LPC_FILTORDER] = {5,5,5,4,3,3,2,2,2,1};
#else /* BITS_FOR_LPC == 80	*/
/* 80-bit LPC10, 4.8 Kbit/s */
static int parambits[LPC_FILTORDER] = {8,8,8,8,8,8,8,8,8,8};
#endif

static float logmaxminper;
static int sizeofparm;	/* computed by lpc_init */

/*
static int lrintf(float inval)
{
    float tmp = FIST_FLOAT_MAGIC_S + inval;
    int res = ((FP_BITS(tmp)<<10)-0x80000000);
    return res>>10;
}
*/
static void auto_correl1(float *w, int n, float *r)
{
    int i, k;
    
    for (k=0; k <= MAXPER; k++, n--) {
        r[k] = 0.0;
        for (i=0; i < n; i++) {
            r[k] += (w[i] *  w[i+k]);
        }
    }
}

static void auto_correl2(float *w, int n, float *r)
{
    int i, k;
    
    for (k=0; k <= LPC_FILTORDER; k++, n--) {
        r[k] = 0.0;
        for (i=0; i < n; i++) {
            r[k] += (w[i] *  w[i+k]);
        }
    }
}

static void durbin(float r[], int p, float k[], float *g)
{
    int i, j;
    float a[LPC_FILTORDER+1], at[LPC_FILTORDER+1], e;
    
    for (i=0; i <= p; i++) a[i] = at[i] = 0.0;
    
    e = r[0];
    for (i=1; i <= p; i++) {
        k[i] = -r[i];
        for (j=1; j < i; j++) {
            at[j] = a[j];
            k[i] -= a[j] * r[i-j];
        }
        if (e == 0) {  /* fix by John Walker */
            *g = 0;
            return;
        }
        k[i] /= e;
        a[i] = k[i];
        for (j=1; j < i; j++) a[j] = at[j] + k[i] * at[i-j];
        e *= 1.0f - k[i]*k[i];
    }
    if (e < 0) {
        e = 0; /* fix by John Walker */
    }
    *g = (float)sqrt(e);
}

/* Enzo's streamlined pitch extractor - on the signal, not the residue */

static void calc_pitch(float *w, int len, float *per)
{
    int i, j, rpos;
    float d[MAXWINDOW/DOWN], r[MAXPER+1], rmax;
    float rval, rm, rp;
    float x, y;
    float vthresh;

    /* decimation */
    for (i=0, j=0; i < len; i+=DOWN) 
        d[j++] = w[i];
    
    auto_correl1(d, len/DOWN, r); 
    
    /* find peak between MINPER and MAXPER */
    x = 1;
    rpos = 0;
    rmax = 0.0;
    y = 0;
    rm = 0;
    rp = 0;

    vthresh = 0.;
     
    for (i = 1; i < MAXPER; i++) {
        rm = r[i-1];
        rp = r[i+1];
        y = rm+r[i]+rp; /* find max of integral from i-1 to i+1 */
        
        if ((y > rmax) && (r[i] > rm) && (r[i] > rp) && (i > MINPER)) {
            rmax = y;
            rpos = i;
        }
    }
    
    /* consider adjacent values */
    rm = r[rpos-1];
    rp = r[rpos+1];
    
#if 0
    {
        float a, b, c, x, y;
        /* parabolic interpolation */
        a = 0.5f * rm - rmax + 0.5f * rp;
        b = -0.5f*rm*(2.0f*rpos+1.0f) + 2.0f*rpos*rmax + 0.5f*rp*(1.0f-2.0f*rpos);
        c = 0.5f*rm*(rpos*rpos+rpos) + rmax*(1.0f-rpos*rpos) + 0.5f*rp*(rpos*rpos-rpos);
        
        /* find max of interpolating parabole */
        x = -b / (2.0f * a);
        y = a*x*x + b*x + c;
        
        rmax = y;
        /* normalize, so that 0. < rval < 1. */ 
        rval = (r[0] == 0 ? 1.0f : rmax / r[0]);
    }
#else
    if(rpos > 0) {
        x = ((rpos-1)*rm + rpos*r[rpos] + (rpos+1)*rp)/(rm+r[rpos]+rp); 
    }
    /* normalize, so that 0. < rval < 1. */ 
    rval = (r[0] == 0 ? 0 : r[rpos] / r[0]);
#endif
    
    /* periods near the low boundary and at low volumes
    are usually spurious and 
    manifest themselves as annoying mosquito buzzes */
    
    *per = 0;	/* default: unvoiced */
    if ( x > MINPER &&  /* x could be < MINPER or even < 0 if pos == MINPER */
        x < (MAXPER+1) /* same story */
        ) {
        
        vthresh = 0.6f; 
        if(r[0] > 0.002)	   /* at low volumes (< 0.002), prefer unvoiced */ 
            vthresh = 0.25;       /* drop threshold at high volumes */
        
        if(rval > vthresh)
            *per = x * DOWN;
    }
}

/* Initialization of various parameters */

openlpc_encoder_state *create_openlpc_encoder_state(void)
{
    openlpc_encoder_state *state;
    
    state = (openlpc_encoder_state *)calloc(1, sizeof(openlpc_encoder_state));
    
    return state;
}


void init_openlpc_encoder_state(openlpc_encoder_state *st, int framelen)
{
    int i, j;
    
    st->framelen = framelen;
    
    st->buflen = framelen*3/2;
    /*  (st->buflen > MAXWINDOW) return -1;*/
    
    for(i=0, j=0; i<sizeof(parambits)/sizeof(parambits[0]); i++) {
        j += parambits[i];
    }
    sizeofparm = (j+7)/8 + 2;
    for (i = 0; i < st->buflen; i++) {
        st->s[i] = 0.0;
        st->h[i] = (float)(WSCALE*(0.54 - 0.46 * cos(2 * M_PI * i / (st->buflen-1.0))));
    }
    /* init the filters */
    st->xv1[0] = st->xv1[1] = st->xv1[2] = st->yv1[0] = st->yv1[1] = st->yv1[2] = 0.0f;
    st->xv2[0] = st->xv2[1] = st->yv2[0] = st->yv2[1] = 0.0f;
    st->xv3[0] = st->yv3[0] = st->yv3[1] = st->yv3[2] = 0.0f;
    st->xv4[0] = st->xv4[1] = st->yv4[0] = st->yv4[1] = 0.0f;
    
    logmaxminper = (float)log((float)MAXPER/MINPER);
    
}

void destroy_openlpc_encoder_state(openlpc_encoder_state *st)
{
    if(st != NULL)
    {
        free(st);
        st = NULL;
    }
}

/* LPC Analysis (compression) */

int openlpc_encode(const short *buf, unsigned char *parm, openlpc_encoder_state *st)
{
    int i, j;
    float per, gain, k[LPC_FILTORDER+1];
    float per1, per2;
    float xv10, xv11, xv12, yv10, yv11, yv12,
        xv20, xv21, yv20, yv21,
        xv30, yv30, yv31, yv32,
        xv40, xv41, yv40, yv41;
    
    xv10 = st->xv1[0];
    xv11 = st->xv1[1];
    xv12 = st->xv1[2];
    yv10 = st->yv1[0];
    yv11 = st->yv1[1];
    yv12 = st->yv1[2];
    xv30 = st->xv3[0];
    yv30 = st->yv3[0];
    yv31 = st->yv3[1];
    yv32 = st->yv3[2];
    for(i = 0; i < LPC_FILTORDER + 1; i++) k[i] = 0;
    /* convert short data in buf[] to signed lin. data in s[] and prefilter */
    for (i=0, j=st->buflen - st->framelen; i < st->framelen; i++, j++) {
        
        float u = (float)(buf[i]/32768.0f);
        
        /* Anti-hum 2nd order Butterworth high-pass, 100 Hz corner frequency */
        /* Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
        mkfilter -Bu -Hp -o 2 -a 0.0125 -l -z */
        
        xv10 = xv11;
        xv11 = xv12; 
        xv12 = (float)(u * 0.94597831f); /* /GAIN */
        
        yv10 = yv11;
        yv11 = yv12; 
        yv12 = (float)((xv10 + xv12) - 2 * xv11
            + ( -0.8948742499f * yv10) + ( 1.8890389823f * yv11));
        
        u = st->s[j] = yv12;	/* also affects input of next stage, to the LPC filter synth */
        
        /* low-pass filter s[] -> y[] before computing pitch */
        /* second-order Butterworth low-pass filter, corner at 300 Hz */
        /* Digital filter designed by mkfilter/mkshape/gencode   A.J. Fisher
        MKFILTER.EXE -Bu -Lp -o 2 -a 0.0375 -l -z */
        
        /*st->xv3[0] = (float)(u / 2.127814584e+001);*/ /* GAIN */
        xv30 = (float)(u * 0.04699658f); /* GAIN */
        yv30 = yv31;
        yv31 = yv32; 
        yv32 = xv30 + (float)(( -0.7166152306f * yv30) + (1.6696186545f * yv31));
        st->y[j] = yv32;
    }
    st->xv1[0] = xv10;
    st->xv1[1] = xv11;
    st->xv1[2] = xv12;
    st->yv1[0] = yv10;
    st->yv1[1] = yv11;
    st->yv1[2] = yv12;
    st->xv3[0] = xv30;
    st->yv3[0] = yv30;
    st->yv3[1] = yv31;
    st->yv3[2] = yv32;
#ifdef PREEMPH
    /* operate optional preemphasis s[] -> s[] on the newly arrived frame */
    xv20 = st->xv2[0];
    xv21 = st->xv2[1];
    yv20 = st->yv2[0];
    yv21 = st->yv2[1];
    xv40 = st->xv4[0];
    xv41 = st->xv4[1];
    yv40 = st->yv4[0];
    yv41 = st->yv4[1];
    for (j=st->buflen - st->framelen; j < st->buflen; j++) {
        float u = st->s[j];
        
        /* handcoded filter: 1 zero at 640 Hz, 1 pole at 3200 */
#define TAU (FS/3200.f)
#define RHO (0.1f)
        xv20 = xv21; 	/* e(n-1) */
        xv21 = (float)(u * 1.584f);		/* e(n)	, add 4 dB to compensate attenuation */
        yv20 = yv21;
        yv21 = (float)(TAU/(1.f+RHO+TAU) * yv20 	 /* u(n) */
            + (RHO+TAU)/(1.f+RHO+TAU) * xv21
            - TAU/(1.f+RHO+TAU) * xv20);
        u = yv21;
        
        /* cascaded copy of handcoded filter: 1 zero at 640 Hz, 1 pole at 3200 */
        xv40 = xv41;
        xv41 = (float)(u * 1.584f);
        yv40 = yv41;
        yv41 = (float)(TAU/(1.f+RHO+TAU) * yv40
            + (RHO+TAU)/(1.f+RHO+TAU) * xv41
            - TAU/(1.f+RHO+TAU) * xv40);
        u = yv41;
        
        st->s[j] = u;
    }
    st->xv2[0] = xv20;
    st->xv2[1] = xv21;
    st->yv2[0] = yv20;
    st->yv2[1] = yv21;
    st->xv4[0] = xv40;
    st->xv4[1] = xv41;
    st->yv4[0] = yv40;
    st->yv4[1] = yv41;
#endif
    
    /* operate windowing s[] -> w[] */
    
    for (i=0; i < st->buflen; i++)
        st->w[i] = st->s[i] * st->h[i];
    
    /* compute LPC coeff. from autocorrelation (first 11 values) of windowed data */
    auto_correl2(st->w, st->buflen, st->r);
    durbin(st->r, LPC_FILTORDER, k, &gain);
    
    /* calculate pitch */
    calc_pitch(st->y, st->framelen, &per1);                 /* first 2/3 of buffer */
    calc_pitch(st->y + st->buflen - st->framelen, st->framelen, &per2); /* last 2/3 of buffer */
    if(per1 > 0 && per2 >0)
        per = (per1+per2)/2;
    else if(per1 > 0)
        per = per1;
    else if(per2 > 0)
        per = per2;
    else
        per = 0;
    
    /* logarithmic q.: 0 = MINPER, 256 = MAXPER */
    parm[0] = (unsigned char)(per == 0? 0 : (unsigned char)(log(per/(REAL_MINPER)) / logmaxminper * (1<<8)));
    
#ifdef LINEAR_G_Q
    i = gain * (1<<7);
    if(i > 255) 	/* bug fix by EM */
        i = 255;
#else
    i = (int)(float)(256.0f * log(1 + (2.718-1.f)/10.f*gain)); /* deriv = 5.82 allowing to reserve 2 bits */
    if(i > 255) i = 255;	 /* reached when gain = 10 */
    i = (i+2) & 0xfc;
#endif
    
    parm[1] = (unsigned char)i;
    
    if(per1 > 0)
        parm[1] |= 1;
    if(per2 > 0)
        parm[1] |= 2;
    
    for(j=2; j < sizeofparm; j++)
        parm[j] = 0;
    
    for (i=0; i < LPC_FILTORDER; i++) {
        int bitamount = parambits[i];
        int bitc8 = 8-bitamount;
        int q = (1 << bitc8);  /* quantum: 1, 2, 4... */
        float u = k[i+1];
        int iu;
#ifdef ARCSIN_Q
        if(i < 2) u = (float)(asin(u)*2.f/M_PI);
#endif
        u *= 127;
        if(u < 0)
            u += (0.6f * q);
        else
            u += (0.4f * q); /* highly empirical! */
        
        iu = lrintf(u);
        iu = iu & 0xff; /* keep only 8 bits */
        
        /* make room at the left of parm array shifting left */
        for(j=sizeofparm-1; j >= 3; j--) {
            parm[j] = (unsigned char)((parm[j] << bitamount) | (parm[j-1] >> bitc8));
        }
        parm[2] = (unsigned char)((parm[2] << bitamount) | (iu >> bitc8)); /* parm[2] */
    }
    
    bcopy(st->s + st->framelen, st->s, (st->buflen - st->framelen)*sizeof(st->s[0]));
    bcopy(st->y + st->framelen, st->y, (st->buflen - st->framelen)*sizeof(st->y[0]));
    
    return sizeofparm;
}

openlpc_decoder_state *create_openlpc_decoder_state(void)
{
    openlpc_decoder_state *state;
    
    state = (openlpc_decoder_state *)calloc(1, sizeof(openlpc_decoder_state));
    
    return state;
}

void init_openlpc_decoder_state(openlpc_decoder_state *st, int framelen)
{
    int i, j;
    
    st->Oldper = 0.0f;
    st->OldG = 0.0f;
    for (i = 0; i <= LPC_FILTORDER; i++) {
        st->Oldk[i] = 0.0f;
        st->bp[i] = 0.0f;
    }
    st->pitchctr = 0;
    st->exc = 0.0f;
    logmaxminper = (float)log((float)MAXPER/MINPER);
    
    for(i=0, j=0; i<sizeof(parambits)/sizeof(parambits[0]); i++) {
        j += parambits[i];
    }
    sizeofparm = (j+7)/8 + 2;

    /* test for a valid frame len? */
    st->framelen = framelen;
    st->buflen = framelen*3/2;
}

/* LPC Synthesis (decoding) */

int openlpc_decode(sp_data *sp, unsigned char *parm, short *buf, openlpc_decoder_state *st)
{
    int i, j, flen=st->framelen;
    float per, gain, k[LPC_FILTORDER+1];
    float f, u, newgain, Ginc, Newper, perinc;
    float Newk[LPC_FILTORDER+1], kinc[LPC_FILTORDER+1];
    float gainadj;
    int hframe;
    float hper[2];
    int ii;
    float bp0, bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8, bp9, bp10;
            float kj;
    
    bp0 = st->bp[0];
    bp1 = st->bp[1];
    bp2 = st->bp[2];
    bp3 = st->bp[3];
    bp4 = st->bp[4];
    bp5 = st->bp[5];
    bp6 = st->bp[6];
    bp7 = st->bp[7];
    bp8 = st->bp[8];
    bp9 = st->bp[9];
    bp10 = st->bp[10];
    
    per = (float)(parm[0]);
    
    per = (float)(per == 0? 0: REAL_MINPER * exp(per/(1<<8) * logmaxminper));

    hper[0] = hper[1] = per;

    if((parm[1] & 0x1) == 0) hper[0] = 0;
    if((parm[1] & 0x2) == 0) hper[1] = 0;
    
#ifdef LINEAR_G_Q
    gain = (float)parm[1] / (1<<7);
#else
    gain = (float)parm[1] / 256.f;
    gain = (float)((exp(gain) - 1)/((2.718-1.f)/10));
#endif
    
    k[0] = 0.0;
    
    for (i=LPC_FILTORDER-1; i >= 0; i--) {
        int bitamount = parambits[i];
        int bitc8 = 8-bitamount;
        /* casting to char should set the sign properly */
        char c = (char)(parm[2] << bitc8);
        
        for(j=2; j<sizeofparm; j++)
            parm[j] = (unsigned char)((parm[j] >> bitamount) | (parm[j+1] << bitc8)); 
        
        k[i+1] = ((float)c / (1<<7));
#ifdef ARCSIN_Q
        if(i<2) k[i+1] = (float)sin(M_PI/2*k[i+1]);
#endif
    }
    
    /* k[] are the same in the two subframes */
    for (i=1; i <= LPC_FILTORDER; i++) {
        Newk[i] = st->Oldk[i];
        kinc[i] = (k[i] - st->Oldk[i]) / flen;
    }
    
    /* Loop on two half frames */
    
    for(hframe=0, ii=0; hframe<2; hframe++) {
        
        Newper = st->Oldper;
        newgain = st->OldG;
        
        Ginc = (gain - st->OldG) / (flen/2);
        per = hper[hframe];
        
        if (per == 0.0) {			 /* if unvoiced */
            gainadj = /* 1.5874 * */ (float)sqrt(3.0f/st->buflen);
        } else {
            gainadj = (float)sqrt(per/st->buflen); 
        }
        
        /* Interpolate period ONLY if both old and new subframes are voiced, gain and K always */ 
        
        if (st->Oldper != 0 && per != 0) {
            perinc = (per - st->Oldper) / (flen/2);
        } else {
            perinc = 0.0f; 
            Newper = per;
        }
        
        if (Newper == 0.f) st->pitchctr = 0;
        
        for (i=0; i < flen/2; i++, ii++) {
            if (Newper == 0.f) {
                u = (float)(((sp_rand(sp)*(1/(1.0f+RAND_MAX))) - 0.5f ) * newgain * gainadj); 
            } else {			/* voiced: send a delta every per samples */
                /* triangular excitation */
                if (st->pitchctr == 0) {
                    st->exc = newgain * 0.25f * gainadj;
                    st->pitchctr = (int) Newper;
                } else {
                    st->exc -= newgain/Newper * 0.5f * gainadj;
                    st->pitchctr--;
                }
                u = st->exc;
            }
            f = u;
	        /* excitation */
            kj = Newk[10];
            f -= kj * bp9;
            bp10 = bp9 + kj * f;

            kj = Newk[9];
            f -= kj * bp8;
            bp9 = bp8 + kj * f;

            kj = Newk[8];
            f -= kj * bp7;
            bp8 = bp7 + kj * f;

            kj = Newk[7];
            f -= kj * bp6;
            bp7 = bp6 + kj * f;

            kj = Newk[6];
            f -= kj * bp5;
            bp6 = bp5 + kj * f;

            kj = Newk[5];
            f -= kj * bp4;
            bp5 = bp4 + kj * f;

            kj = Newk[4];
            f -= kj * bp3;
            bp4 = bp3 + kj * f;

            kj = Newk[3];
            f -= kj * bp2;
            bp3 = bp2 + kj * f;

            kj = Newk[2];
            f -= kj * bp1;
            bp2 = bp1 + kj * f;

            kj = Newk[1];
            f -= kj * bp0;
            bp1 = bp0 + kj * f;

            bp0 = f;
            u = f;

            if (u  < -0.9999f) {
                u = -0.9999f;
            } else if (u > 0.9999f) {
                u = 0.9999f;
            }

            buf[ii] = (short)lrintf(u * 32767.0f);

            Newper += perinc;
            newgain += Ginc;
            for (j=1; j <= LPC_FILTORDER; j++) Newk[j] += kinc[j];

        }

        st->Oldper = per;
        st->OldG = gain;
    }
    st->bp[0] = bp0;
    st->bp[1] = bp1;
    st->bp[2] = bp2;
    st->bp[3] = bp3;
    st->bp[4] = bp4;
    st->bp[5] = bp5;
    st->bp[6] = bp6;
    st->bp[7] = bp7;
    st->bp[8] = bp8;
    st->bp[9] = bp9;
    st->bp[10] = bp10;

    for (j=1; j <= LPC_FILTORDER; j++) st->Oldk[j] = k[j];

    return flen;
}

void destroy_openlpc_decoder_state(openlpc_decoder_state *st)
{
    if(st != NULL)
    {
        free(st);
        st = NULL;
    }
}

void openlpc_sr(float sr)
{
    my_fs = sr;
}

size_t openlpc_get_encoder_state_size(void)
{
    return sizeof(openlpc_encoder_state);
}

size_t openlpc_get_decoder_state_size(void)
{
    return sizeof(openlpc_decoder_state);
}


int sp_lpf18_create(sp_lpf18 **p)
{
    *p = malloc(sizeof(sp_lpf18));
    return SP_OK;
}

int sp_lpf18_destroy(sp_lpf18 **p)
{
    free(*p);
    return SP_OK;
}

int sp_lpf18_init(sp_data *sp, sp_lpf18 *p)
{
    p->cutoff = 1000;
    p->res = 0.8;
    p->dist = 2;

    p->ay1 = 0.0;
    p->ay2 = 0.0;
    p->aout = 0.0;
    p->lastin = 0.0;
    p->onedsr = 1.0 / sp->sr;
    return SP_OK;
}

int sp_lpf18_compute(sp_data *sp, sp_lpf18 *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT ay1 = p->ay1;
    SPFLOAT ay2 = p->ay2;
    SPFLOAT aout = p->aout;
    SPFLOAT lastin = p->lastin;
    double value = 0.0;
    int   flag = 1;
    SPFLOAT lfc=0, lrs=0, kres=0, kfcn=0, kp=0, kp1=0,  kp1h=0;
    double lds = 0.0;

    SPFLOAT fco, res, dist;
    SPFLOAT ax1  = lastin;
    SPFLOAT ay11 = ay1;
    SPFLOAT ay31 = ay2;
    fco = p->cutoff;
    res = p->res;
    dist = p->dist;

    if (fco != lfc || flag) {
        lfc = fco;
        kfcn = 2.0 * fco * p->onedsr;
        kp = ((-2.7528 * kfcn + 3.0429) * kfcn +
                1.718) * kfcn - 0.9984;
        kp1 = kp + 1.0;
        kp1h = 0.5 * kp1;
        flag = 1;
    }

    if (res != lrs || flag) {
        lrs = res;
        kres = res * (((-2.7079 * kp1 + 10.963) * kp1
                           - 14.934) * kp1 + 8.4974);
        flag = 1;
    }

    if (dist != lds || flag) {
        lds = dist;
        value = 1.0 + (dist * (1.5 + 2.0 * res * (1.0 - kfcn)));
    }

    flag = 0;
    lastin = *in - tanh(kres*aout);
    ay1 = kp1h * (lastin + ax1) - kp * ay1;
    ay2 = kp1h * (ay1 + ay11) - kp * ay2;
    aout = kp1h * (ay2 + ay31) - kp * aout;

    *out = tanh(aout * value);

    p->ay1 = ay1;
    p->ay2 = ay2;
    p->aout = aout;
    p->lastin = lastin;
    return SP_OK;
}


int sp_maygate_create(sp_maygate **p)
{
    *p = malloc(sizeof(sp_maygate));
    return SP_OK;
}

int sp_maygate_destroy(sp_maygate **p)
{
    free(*p);
    return SP_OK;
}

int sp_maygate_init(sp_data *sp, sp_maygate *p)
{
    p->prob = 0.0;
    p->gate = 0;
    p->mode = 0;
    return SP_OK;
}

int sp_maygate_compute(sp_data *sp, sp_maygate *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in == 0) {
        if(p->mode) {
            *out = 0;
        } else {
            *out = p->gate;
        }
        return SP_OK;
    }

    if((1.0 * sp_rand(sp) / SP_RANDMAX) <= p->prob) {
        *out = 1;
        p->gate = 1;
    } else {
        *out = 0;
        p->gate = 0;
    }
    return SP_OK;
}

int sp_metro_create(sp_metro **p)
{
    *p = malloc(sizeof(sp_metro));
    return SP_OK;
}

int sp_metro_destroy(sp_metro **p)
{
    free(*p);
    return SP_OK;
}

int sp_metro_init(sp_data *sp, sp_metro *p)
{
    p->iphs = 0;
    p->freq= 2.0;
    SPFLOAT phs = p->iphs;
    int32_t  longphs = phs;
    if (phs >= 0.0){
      p->curphs = (SPFLOAT)phs - (SPFLOAT)longphs;
    }
    p->flag=1;
    p->onedsr = 1.0 / sp->sr;
    return SP_OK;
}

int sp_metro_compute(sp_data *sp, sp_metro *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT phs= p->curphs;
    if (phs == 0.0 && p->flag){
        *out = 1.0;
        p->flag = 0;
    }else if ((phs += p->freq * p->onedsr) >= 1.0){
        *out = 1.0;
        phs -= 1.0;
        p->flag = 0;
    }else{
        *out = 0.0;
    }
    p->curphs = phs;
    return SP_OK;
}

int sp_mincer_create(sp_mincer **p)
{
    *p = malloc(sizeof(sp_mincer));
    return SP_OK;
}

int sp_mincer_destroy(sp_mincer **p)
{
    sp_mincer *pp = *p;
    sp_fft_destroy(&pp->fft);
    sp_auxdata_free(&pp->fwin);
    sp_auxdata_free(&pp->bwin);
    sp_auxdata_free(&pp->prev);
    sp_auxdata_free(&pp->framecount);
    sp_auxdata_free(&pp->outframe);
    sp_auxdata_free(&pp->win);
    free(*p);
    return SP_OK;
}

static int find_power(int n) {
    int pow = -1;
    while(n > 0) {
        n >>= 1;
        pow++;
    }
    return pow;
}

int sp_mincer_init(sp_data *sp, sp_mincer *p, sp_ftbl *ft, int winsize)
{
    p->ft = ft;
    p->idecim = 4;
    p->iN = winsize;
    p->lock = 1;
    p->pitch = 1;
    p->amp = 1;
    p->time = 0;
    int N =  p->iN, ui;
    unsigned int size;
    int decim = p->idecim;
    int pow;

    /* find power to use for fft */
    
    pow = find_power(winsize);
    /* 2^11 = 2048, the default fftsize, will probably not change */
    sp_fft_init(&p->fft, pow);


    if (decim == 0) decim = 4;

    p->hsize = N/decim;
    p->cnt = p->hsize;
    p->curframe = 0;
    p->pos = 0;

    size = (N+2)*sizeof(SPFLOAT);
    sp_auxdata_alloc(&p->fwin, size);
    sp_auxdata_alloc(&p->bwin, size);
    sp_auxdata_alloc(&p->prev, size);
    size = decim*sizeof(int);
    sp_auxdata_alloc(&p->framecount, size);
    {
      int k=0;
        for (k=0; k < decim; k++) {
            ((int *)(p->framecount.ptr))[k] = k*N;
        }
    }
    size = decim*sizeof(SPFLOAT)*N;
    sp_auxdata_alloc(&p->outframe, size);
    
    size = N*sizeof(SPFLOAT);
    sp_auxdata_alloc(&p->win, size);
    {
        SPFLOAT x = 2.0 * M_PI/N;
        for (ui=0; ui < N; ui++)
        ((SPFLOAT *)p->win.ptr)[ui] = 0.5 - 0.5 * cos((SPFLOAT)ui*x);
    }

    p->N = N;
    p->decim = decim;

    return SP_OK;
}

int sp_mincer_compute(sp_data *sp, sp_mincer *p, SPFLOAT *in2, SPFLOAT *out)
{
    SPFLOAT pitch = p->pitch, time = p->time, lock = p->lock, amp =p->amp;
    SPFLOAT *tab, frac;
    sp_ftbl *ft = p->ft;
    int N = p->N, hsize = p->hsize, cnt = p->cnt;
    int sizefrs, size, post, i;
    long spos = p->pos;
    SPFLOAT pos;
    SPFLOAT *fwin, *bwin, insig = 0,
    *prev, *win = (SPFLOAT *) p->win.ptr;
    SPFLOAT *outframe;
    SPFLOAT ph_real, ph_im, tmp_real, tmp_im, divi;
    int *framecnt;
    int curframe = p->curframe, decim = p->decim;
    SPFLOAT scaling = (8./decim)/3.;

    if (cnt == hsize) {
        tab = ft->tbl;
        size = (int)ft->size;

        /* spos is the reading position in samples, hsize is hopsize,
        time[n] is current read position in secs
        esr is sampling rate
        */
        spos  = hsize*(long)((time)*sp->sr/hsize);
        sizefrs = size;
        while(spos > sizefrs) spos -= sizefrs;
        while(spos <= 0)  spos += sizefrs;


        pos = spos;
        bwin = (SPFLOAT *) p->bwin.ptr;
        fwin = (SPFLOAT *) p->fwin.ptr;
        prev = (SPFLOAT *)p->prev.ptr;
        framecnt  = (int *)p->framecount.ptr;
        outframe= (SPFLOAT *) p->outframe.ptr;
        /* this loop fills two frames/windows with samples from table,
        reading is linearly-interpolated,
        frames are separated by 1 hopsize
        */
        for (i=0; i < N; i++) {
            /* front window, fwin */
            post = (int) pos;
            frac = pos  - post;
            while (post < 0) post += size;
            while (post >= size) post -= size;
            if(post + 1 <  size)
            insig = tab[post] + frac*(tab[post+ 1] - tab[post]);
            else insig = tab[post];

            /* window it */
            fwin[i] = insig * win[i]; 
            /* back windo, bwin */
            post = (int) (pos - hsize*pitch);
            post *= 1;
            post += 0;
            while(post < 0) post += size;
            while(post >= size) post -= size;
            if(post + 1<  size)
            insig = tab[post] + frac*(tab[post + 1] - tab[post]);
            else insig = tab[post];
            bwin[i] = insig * win[i];  /* window it */
            /* increment read pos according to pitch transposition */
            pos += pitch;
        }

        /* take the FFT of both frames
        re-order Nyquist bin from pos 1 to N
        */
        sp_fftr(&p->fft, bwin, N);
        bwin[N] = bwin[1];
        bwin[N+1] = 0.0;
        sp_fftr(&p->fft, fwin, N);
        fwin[N] = fwin[1];
        fwin[N+1] = 0.0;

        /* phase vocoder processing */

        for (i=0; i < N + 2; i+=2) {
            /* phases of previous output frame in exponential format,
            obtained by dividing by magnitude */
            divi =  1.0/(hypot(prev[i], prev[i+1]) + 1e-20);
            ph_real  =    prev[i]*divi;
            ph_im =       prev[i+1]*divi;

            /* back window magnitudes, phase differences between
            prev and back windows */
            tmp_real =   bwin[i] * ph_real + bwin[i+1] * ph_im;
            tmp_im =   bwin[i] * ph_im - bwin[i+1] * ph_real;
            bwin[i] = tmp_real;
            bwin[i+1] = tmp_im;
        }

        for (i=0; i < N + 2; i+=2) {
            if (lock) {  /* phase-locking */
                if (i > 0) {
                    if (i < N){
                        tmp_real = bwin[i] + bwin[i-2] + bwin[i+2];
                        tmp_im = bwin[i+1] + bwin[i-1] + bwin[i+3];
                    } else { /* Nyquist */
                        tmp_real = bwin[i] + bwin[i-2];
                        tmp_im = 0.0;
                    } 
                } else { /* 0 Hz */
                    tmp_real = bwin[i] + bwin[i+2];
                    tmp_im = 0.0;
                }
            } else { /* no locking */
                tmp_real = bwin[i];
                tmp_im = bwin[i+1];
            }

            tmp_real += 1e-15;
            divi =  1.0/(hypot(tmp_real, tmp_im));

            /* phases of tmp frame */
            ph_real = tmp_real*divi;
            ph_im = tmp_im*divi;

            /* front window mags, phase sum of
            tmp and front windows */
            tmp_real =   fwin[i] * ph_real - fwin[i+1] * ph_im;
            tmp_im =   fwin[i] * ph_im + fwin[i+1] * ph_real;

            /* phase vocoder output */
            prev[i] = fwin[i] = tmp_real;
            prev[i+1] = fwin[i+1] = tmp_im;
        }
        /* re-order bins and take inverse FFT */
        fwin[1] = fwin[N];
        sp_ifftr(&p->fft, fwin, N);
        /* frame counter */
        framecnt[curframe] = curframe*N;
        /* write to overlapped output frames */
        for (i=0;i<N;i++) outframe[framecnt[curframe]+i] = win[i]*fwin[i];

        cnt=0;
        curframe++;
        if (curframe == decim) curframe = 0;
    }

    framecnt  = (int *) p->framecount.ptr;
    outframe  = (SPFLOAT *) p->outframe.ptr;
    *out = (SPFLOAT)0;
    /* write output */
    for (i = 0; i < decim; i++) {
        *out += outframe[framecnt[i]];
        framecnt[i]++;
    }
    /* scale output */
    *out *= amp*scaling;
    cnt++;

    p->cnt = cnt;
    p->curframe = curframe;

    return SP_OK;
}

int sp_mode_create(sp_mode **p)
{
    *p = malloc(sizeof(sp_mode));
    return SP_OK;
}

int sp_mode_destroy(sp_mode **p)
{
    free(*p);
    return SP_OK;
}

int sp_mode_init(sp_data *sp, sp_mode *p)
{
    p->freq = 500.0;
    p->q = 50;

    p->xnm1 = p->ynm1 = p->ynm2 = 0.0;
    p->a0 = p->a1 = p->a2 = p->d = 0.0;
    p->lfq = -1.0;
    p->lq = -1.0;

    p->sr = sp->sr;

    return SP_OK;
}

int sp_mode_compute(sp_data *sp, sp_mode *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT lfq = p->lfq, lq = p->lq;

    SPFLOAT xn, yn, a0=p->a0, a1=p->a1, a2=p->a2,d=p->d;
    SPFLOAT xnm1 = p->xnm1, ynm1 = p->ynm1, ynm2 = p->ynm2;

    SPFLOAT kfq = p->freq;
    SPFLOAT kq  = p->q;
    if (lfq != kfq || lq != kq) {
        SPFLOAT kfreq  = kfq*(2.0 * M_PI);
        SPFLOAT kalpha = (p->sr/kfreq);
        SPFLOAT kbeta  = kalpha*kalpha;
               d      = 0.5*kalpha;

        lq = kq; lfq = kfq;
        a0     = 1.0/ (kbeta+d/kq);
        a1     = a0 * (1.0-2.0*kbeta);
        a2     = a0 * (kbeta-d/kq);
     }
     xn = *in;

     yn = a0*xnm1 - a1*ynm1 - a2*ynm2;

     xnm1 = xn;
     ynm2 = ynm1;
     ynm1 = yn;

     yn = yn*d;

     *out  = yn;

    p->xnm1 = xnm1;  p->ynm1 = ynm1;  p->ynm2 = ynm2;
    p->lfq = lfq;    p->lq = lq;      p->d = d;
    p->a0 = a0;      p->a1 = a1;      p->a2 = a2;
    return SP_OK;
}


#define SPFLOAT2LONG(x) lrintf(x)


/* John ffitch tanh function to speed up inner loop */ 

static double my_tanh(double x)
{
    /* use the fact that tanh(-x) = - tanh(x)
       and if x>~4 tanh is approx constant 1
       and for small x tanh(x) =~ x
       So giving a cheap approximation */
    int sign = 1;
    if (x<0) { 
        sign=-1; 
        x= -x;
    }
    if (x>=4.0) {
      return sign;
    }
    if (x<0.5) return x*sign;
    return sign*tanh(x);
}

int sp_moogladder_create(sp_moogladder **t){
    *t = malloc(sizeof(sp_moogladder));
    return SP_OK;
}

int sp_moogladder_destroy(sp_moogladder **t){
    free(*t);
    return SP_OK;
}

int sp_moogladder_init(sp_data *sp, sp_moogladder *p){
    p->istor = 0.0;
    p->res = 0.4;
    p->freq = 1000;
    int i;

    if (p->istor == 0.0) {
      for (i = 0; i < 6; i++)
        p->delay[i] = 0.0;
      for (i = 0; i < 3; i++)
        p->tanhstg[i] = 0.0;
      p->oldfreq = 0.0;
      p->oldres = -1.0;     /* ensure calculation on first cycle */
    }
    return SP_OK;
}

int sp_moogladder_compute(sp_data *sp, sp_moogladder *p, SPFLOAT *in, SPFLOAT *out){
    SPFLOAT freq = p->freq;
    SPFLOAT res = p->res;
    SPFLOAT res4;
    SPFLOAT *delay = p->delay;
    SPFLOAT *tanhstg = p->tanhstg;
    SPFLOAT stg[4], input;
    SPFLOAT acr, tune;
#define THERMAL (0.000025) /* (1.0 / 40000.0) transistor thermal voltage  */
    int     j, k;

    if (res < 0) res = 0;

    if (p->oldfreq != freq || p->oldres != res) {
        SPFLOAT f, fc, fc2, fc3, fcr;
        p->oldfreq = freq;
        /* sr is half the actual filter sampling rate  */
        fc =  (SPFLOAT)(freq/sp->sr);
        f  =  0.5*fc;
        fc2 = fc*fc;
        fc3 = fc2*fc;
        /* frequency & amplitude correction  */
        fcr = 1.8730*fc3 + 0.4955*fc2 - 0.6490*fc + 0.9988;
        acr = -3.9364*fc2 + 1.8409*fc + 0.9968;
        tune = (1.0 - exp(-((2 * M_PI)*f*fcr))) / THERMAL;   /* filter tuning  */
        p->oldres = res;
        p->oldacr = acr;
        p->oldtune = tune;
    } else {
        res = p->oldres;
        acr = p->oldacr;
        tune = p->oldtune;
    }
    res4 = 4.0*(SPFLOAT)res*acr;

    /* oversampling  */
    for (j = 0; j < 2; j++) {
        /* filter stages  */
        input = *in - res4 /*4.0*res*acr*/ *delay[5];
        delay[0] = stg[0] = delay[0] + tune*(my_tanh(input*THERMAL) - tanhstg[0]);
        for (k = 1; k < 4; k++) {
          input = stg[k-1];
          stg[k] = delay[k]
            + tune*((tanhstg[k-1] = my_tanh(input*THERMAL))
                    - (k != 3 ? tanhstg[k] : my_tanh(delay[k]*THERMAL)));
          delay[k] = stg[k];
        }
        /* 1/2-sample delay for phase compensation  */
        delay[5] = (stg[3] + delay[4])*0.5;
        delay[4] = stg[3];
    }
    *out = (SPFLOAT) delay[5];
    return SP_OK;
}


int sp_noise_create(sp_noise **ns)
{
    *ns = malloc(sizeof(sp_noise));
    return SP_OK;
}

int sp_noise_init(sp_data *sp, sp_noise *ns)
{
    ns->amp = 1.0;
    return SP_OK;
}

int sp_noise_compute(sp_data *sp, sp_noise *ns, SPFLOAT *in, SPFLOAT *out)
{
    *out = ((sp_rand(sp) % SP_RANDMAX) / (SP_RANDMAX * 1.0));
    *out = (*out * 2) - 1;
    *out *= ns->amp;
    return SP_OK;
}

int sp_noise_destroy(sp_noise **ns)
{
    free(*ns);
    return SP_OK;
}


int nano_dict_add(nano_dict *dict, const char *name)
{
    nano_entry *entry = malloc(sizeof(nano_entry));
    entry->size = 0;
    entry->speed = 1;
    entry->pos = 0;
    strcpy(entry->name, name);
    dict->last->next = entry;
    dict->last = entry;
    dict->nval++;
    return SP_OK;
}

int nano_ini_handler(void *user, const char *section, const char *name,
        const char *value)
{
    nanosamp *ss = user;
    nano_dict *dict = &ss->dict;
    const char *entry_name = dict->last->name;

    if(dict->init) {
        nano_dict_add(dict, section);
        dict->init = 0;
    } else if(strncmp(entry_name, section, 50) != 0) {
        nano_dict_add(dict, section);
    }

    dict->last->speed = 1.0;

    if(strcmp(name, "pos") == 0) {
        dict->last->pos = (uint32_t)(atof(value) * ss->sr);
    } else if(strcmp(name, "size") == 0) {
        dict->last->size = (uint32_t)(atof(value) * ss->sr);
    } else if(strcmp(name, "speed") == 0) {
        dict->last->speed = atof(value);
    }

    return SP_OK;
}

int nano_create(nanosamp **smp, const char *ini, int sr)
{
    *smp = malloc(sizeof(nanosamp));
    nanosamp *psmp = *smp;
    strcpy(psmp->ini, ini);
    psmp->dict.last = &psmp->dict.root;
    psmp->dict.nval = 0;
    psmp->dict.init = 1;
    psmp->selected = 0;
    psmp->curpos = 0;
    psmp->sr = sr;
    if(ini_parse(psmp->ini, nano_ini_handler, psmp) < 0) {
        printf("Can't load file %s\n", psmp->ini);
        return SP_NOT_OK;
    }

    return SP_OK;
}

int nano_select_from_index(nanosamp *smp, uint32_t pos)
{
    pos %= smp->dict.nval;
    smp->selected = 1;
    smp->sample = smp->index[pos];
    smp->curpos = 0;
    return SP_OK;
}

uint32_t nano_keyword_to_index(nanosamp *smp, const char *keyword)
{
    uint32_t i;
    for (i = 0; i < smp->dict.nval; i++) {
        if(strcmp(keyword, smp->index[i]->name)) {
            return i;
        }
    }
    return 0;
}

int nano_select(nanosamp *smp, const char *keyword)
{
    uint32_t i;
    nano_dict *dict = &smp->dict;
    nano_entry *entry = dict->root.next;
    smp->curpos = 0;
    smp->selected = 0;
    for(i = 0; i < dict->nval; i++) {
        if(strncmp(keyword, entry->name, 50) == 0) {
            smp->selected = 1;
            smp->sample = entry;
            smp->curpos = 0;
            break;
        } else {
            entry = entry->next;
        }
    }

    if(smp->selected == 1) return SP_OK;
    else return SP_NOT_OK;
}


int nano_compute(sp_data *sp, nanosamp *smp, SPFLOAT *out)
{
    if(!smp->selected) {
        *out = 0;
        return SP_NOT_OK;
    }

    if(smp->curpos < (SPFLOAT)smp->sample->size) {
        SPFLOAT x1 = 0 , x2 = 0, frac = 0, tmp = 0;
        uint32_t index = 0;
        SPFLOAT *tbl = smp->ft->tbl;
        tmp = (smp->curpos + smp->sample->pos);
        index = floorf(tmp);
        frac = fabs(tmp - index);

        if(index >= smp->ft->size) {
            index = smp->ft->size - 1;
        }

        x1 = tbl[index];
        x2 = tbl[index + 1];
        *out = x1 + (x2 - x1) * frac;
        smp->curpos += smp->sample->speed;
    } else {
        smp->selected = 0;
        *out = 0;
    }

    return SP_OK;
}

int nano_dict_destroy(nano_dict *dict)
{
    int i;
    nano_entry *entry, *next;
    entry = dict->root.next;

    for(i = 0; i < dict->nval; i++) {
        next = entry->next;
        free(entry);
        entry = next;
    }
    return SP_OK;
}

int nano_destroy(nanosamp **smp)
{
    nanosamp *psmp = *smp;
    nano_dict_destroy(&psmp->dict);
    free(*smp);
    return SP_OK;
}



int nano_create_index(nanosamp *smp)
{
    nano_dict *dict = &smp->dict;
    smp->index = malloc(dict->nval * sizeof(nano_entry *));
    int i;
    nano_entry *entry, *next;
    entry = dict->root.next;

    for(i = 0; i < dict->nval; i++) {
        next = entry->next;
        smp->index[i] = entry;
        entry = next;
    }
    return SP_OK;
}

int nano_destroy_index(nanosamp *smp)
{
    free(smp->index);
    return SP_OK;
}

int sp_nsmp_create(sp_nsmp **p)
{
    *p = malloc(sizeof(sp_nsmp));
    return SP_OK;
}

int sp_nsmp_destroy(sp_nsmp **p)
{
    sp_nsmp *pp = *p;
    nano_destroy_index(pp->smp);
    nano_destroy(&pp->smp);
    free(*p);
    return SP_OK;
}

int sp_nsmp_init(sp_data *sp, sp_nsmp *p, sp_ftbl *ft, int sr, const char *ini)
{
    if (nano_create(&p->smp, ini, sr) == SP_NOT_OK) {
        nano_destroy(&p->smp);
        return SP_NOT_OK;
    }
    nano_create_index(p->smp);
    p->smp->sr = sr;
    p->index= 0;
    p->triggered = 0;
    p->smp->ft = ft;
    return SP_OK;
}

int sp_nsmp_compute(sp_data *sp, sp_nsmp *p, SPFLOAT *trig, SPFLOAT *out)
{
    if (*trig != 0) {
       p->triggered = 1;
       nano_select_from_index(p->smp, p->index);
    }

    if(p->triggered == 1) {
        nano_compute(sp, p->smp, out);
    } else {
        *out = 0;
    }

    return SP_OK;
}

int sp_nsmp_print_index(sp_data *sp, sp_nsmp *p)
{
    uint32_t i;
    for(i = 0; i < p->smp->dict.nval; i++) {
        printf("%d: key = %s\n", i, p->smp->index[i]->name);
    }
    return SP_OK;
}


int sp_osc_create(sp_osc **osc)
{
    *osc = malloc(sizeof(sp_osc));
    return SP_OK;
}

int sp_osc_destroy(sp_osc **osc)
{
    free(*osc);
    return SP_NOT_OK;
}

int sp_osc_init(sp_data *sp, sp_osc *osc, sp_ftbl *ft, SPFLOAT iphs)
{
    osc->freq = 440.0;
    osc->amp = 0.2;
    osc->tbl = ft;
    osc->iphs = fabs(iphs);
    osc->inc = 0;
    if (osc->iphs >= 0){
        osc->lphs = ((int32_t)(osc->iphs * SP_FT_MAXLEN)) & SP_FT_PHMASK;
    }

    return SP_OK;
}

int sp_osc_compute(sp_data *sp, sp_osc *osc, SPFLOAT *in, SPFLOAT *out)
{
    sp_ftbl *ftp;
    SPFLOAT amp, cps, fract, v1, v2, *ft;
    int32_t phs, lobits;
    int32_t pos;
    SPFLOAT sicvt = osc->tbl->sicvt;

    ftp = osc->tbl;
    lobits = osc->tbl->lobits;
    amp = osc->amp;
    cps = osc->freq;
    phs = osc->lphs;
    ft = osc->tbl->tbl;
    
    osc->inc = (int32_t)lrintf(cps * sicvt);

    fract = ((phs) & ftp->lomask) * ftp->lodiv;
    pos = phs>>lobits;
    v1 = *(ft + pos);
    v2 = *(ft + ((pos + 1) % ftp->size));
    *out = (v1 + (v2 - v1) * fract) * amp;
    phs += osc->inc;
    phs &= SP_FT_PHMASK;

    osc->lphs = phs;
    return SP_OK;
}


int sp_oscmorph_create(sp_oscmorph **p)
{
    *p = malloc(sizeof(sp_oscmorph));
    return SP_OK;
}

int sp_oscmorph_destroy(sp_oscmorph **p)
{
    free(*p);
    return SP_OK;
}

int sp_oscmorph_init(sp_data *sp, sp_oscmorph *osc, sp_ftbl **ft, int nft, SPFLOAT iphs)
{
    int i;
    osc->freq = 440.0;
    osc->amp = 0.2;
    osc->tbl = ft;
    osc->iphs = fabs(iphs);
    osc->inc = 0;
    osc->lphs = ((int32_t)(osc->iphs * SP_FT_MAXLEN)) & SP_FT_PHMASK;
    osc->wtpos = 0.0;
    osc->nft = nft;
    uint32_t prev = (uint32_t)ft[0]->size;
    for(i = 0; i < nft; i++) {
        if(prev != ft[i]->size) {
            fprintf(stderr, "sp_oscmorph: size mismatch\n");
            return SP_NOT_OK;
        }
        prev = (uint32_t)ft[i]->size;
    }
    return SP_OK;
}

int sp_oscmorph_compute(sp_data *sp, sp_oscmorph *osc, SPFLOAT *in, SPFLOAT *out)
{
    sp_ftbl *ftp1;
    SPFLOAT amp, cps, fract, v1, v2;
    SPFLOAT *ft1, *ft2;
    int32_t phs, lobits, pos;
    SPFLOAT sicvt = osc->tbl[0]->sicvt;

    /* Use only the fractional part of the position or 1 */
    if (osc->wtpos > 1.0) {
        osc->wtpos -= (int)osc->wtpos;
    }
    SPFLOAT findex = osc->wtpos * (osc->nft - 1);
    int index = floor(findex);
    SPFLOAT wtfrac = findex - index;

    lobits = osc->tbl[0]->lobits;
    amp = osc->amp;
    cps = osc->freq;
    phs = osc->lphs;
    ftp1 = osc->tbl[index];
    ft1 = osc->tbl[index]->tbl;

    if(index >= osc->nft - 1) {
        ft2 = ft1;
    } else {
        ft2 = osc->tbl[index + 1]->tbl;
    }
    
    osc->inc = (int32_t)lrintf(cps * sicvt);

    fract = ((phs) & ftp1->lomask) * ftp1->lodiv;

    pos = phs >> lobits;

    v1 = (1 - wtfrac) * 
        *(ft1 + pos) + 
        wtfrac * 
        *(ft2 + pos);
    v2 = (1 - wtfrac) * 
        *(ft1 + ((pos + 1) % ftp1->size))+ 
        wtfrac * 
        *(ft2 + ((pos + 1) % ftp1->size));

    *out = (v1 + (v2 - v1) * fract) * amp;

    phs += osc->inc;
    phs &= SP_FT_PHMASK;

    osc->lphs = phs;
    return SP_OK;
}


int sp_gen_padsynth(sp_data *sp, sp_ftbl *ps, sp_ftbl *amps, 
        SPFLOAT f, SPFLOAT bw) 
{

    int i, nh;
    int N = (int) ps->size;
    int number_harmonics = (int) amps->size;
    SPFLOAT *A = amps->tbl;
    SPFLOAT *smp = ps->tbl;

    SPFLOAT *freq_amp = malloc((N / 2) * sizeof(SPFLOAT));
    SPFLOAT *freq_phase = malloc((N / 2) * sizeof(SPFLOAT));

    for (i=0;i<N/2;i++) freq_amp[i]=0.0;

    for (nh=1;nh<number_harmonics;nh++) {
        SPFLOAT bw_Hz;
        SPFLOAT bwi;
        SPFLOAT fi;
        bw_Hz = (pow(2.0, bw/1200.0) - 1.0) * f * nh;
        bwi = bw_Hz/(2.0*ps->size);
        fi = f*nh/ps->size;
        for (i = 0; i < N/2 ; i++) {
            SPFLOAT hprofile;
            hprofile = sp_padsynth_profile((i / (SPFLOAT) N) - fi, bwi);
            freq_amp[i] += hprofile*A[nh];
        }
    }

    for (i=0;i<N/2;i++) {
        freq_phase[i]= (sp_rand(sp) / (SP_RANDMAX + 1.0)) * 2.0 * M_PI;
    };

    sp_padsynth_ifft(N,freq_amp,freq_phase,smp);
    sp_padsynth_normalize(N,smp);

    free(freq_amp);
    free(freq_phase);
    return SP_OK;
}

/* This is the profile of one harmonic
   In this case is a Gaussian distribution (e^(-x^2))
   The amplitude is divided by the bandwidth to ensure that the harmonic
   keeps the same amplitude regardless of the bandwidth */

SPFLOAT sp_padsynth_profile(SPFLOAT fi, SPFLOAT bwi) 
{
    SPFLOAT x =fi/bwi;
    x *= x;

/* 
 * this avoids computing the e^(-x^2) where 
 * it's results are very close to zero
 */
    if (x>14.71280603) return 0.0;

    return exp(-x)/bwi;
}

int sp_padsynth_ifft(int N, SPFLOAT *freq_amp, 
        SPFLOAT *freq_phase, SPFLOAT *smp) 
{
    int i;
    FFTwrapper *fft;
    FFTwrapper_create(&fft, N);
    FFTFREQS fftfreqs;
    newFFTFREQS(&fftfreqs,N/2);

    for (i=0; i<N/2; i++){
        fftfreqs.c[i]=freq_amp[i]*cos(freq_phase[i]);
        fftfreqs.s[i]=freq_amp[i]*sin(freq_phase[i]);
    };
    freqs2smps(fft, &fftfreqs,smp);
    deleteFFTFREQS(&fftfreqs);
    FFTwrapper_destroy(&fft);
    return SP_OK;
}

/*
    Simple normalization function. It normalizes the sound to 1/sqrt(2)
*/

int sp_padsynth_normalize(int N, SPFLOAT *smp) 
{
    int i;
    SPFLOAT max=0.0;
    for (i=0;i<N;i++) if (fabs(smp[i])>max) max=fabs(smp[i]);
    if (max<1e-5) max=1e-5;
    for (i=0;i<N;i++) smp[i]/=max*1.4142;
    return SP_OK;
}


int sp_pan2_create(sp_pan2 **p)
{
    *p = malloc(sizeof(sp_pan2));
    return SP_OK;
}

int sp_pan2_destroy(sp_pan2 **p)
{
    free(*p);
    return SP_OK;
}

int sp_pan2_init(sp_data *sp, sp_pan2 *p)
{
    p->type = 0;
    p->pan = 0;
    return SP_OK;
}

int sp_pan2_compute(sp_data *sp, sp_pan2 *p, SPFLOAT *in, SPFLOAT *out1, SPFLOAT *out2)
{
    /* Send the signal's input to the output */
    uint32_t type = p->type;
    SPFLOAT pan = (1 + p->pan) * 0.5;
    SPFLOAT cc, ss, l, r;

    type %= 4;

    switch (type) {
        /* Equal power */
        case 0:
        pan = M_PI * 0.5 * pan;
        *out1 = *in * cos(pan);
        *out2 = *in * sin(pan);
        break;

        /* Square root */
        case 1:
        *out1 = *in * sqrt(pan);
        *out2 = *in * sqrt(1.0 - pan);
        break;

        /* simple linear */
        case 2:
        *out1 = *in * (1.0 - pan);
        *out2 = *in * pan;
        break;

        /* Equal power (alternative) */
        case 3:

        cc = cos(M_PI * pan * 0.5);
        ss = sin(M_PI * pan * 0.5);
        l = SQRT2 * (cc + ss) * 0.5;
        r = SQRT2 * (cc - ss) * 0.5;
        *out1 = *in * l;
        *out2 = *in * r;
        break;
    }

    return SP_OK;
}

int sp_panst_create(sp_panst **p)
{
    *p = malloc(sizeof(sp_panst));
    return SP_OK;
}

int sp_panst_destroy(sp_panst **p)
{
    free(*p);
    return SP_OK;
}

int sp_panst_init(sp_data *sp, sp_panst *p)
{
    p->type = 0;
    p->pan = 0;
    return SP_OK;
}

int sp_panst_compute(sp_data *sp, sp_panst *p, SPFLOAT *in1, SPFLOAT *in2, SPFLOAT *out1, SPFLOAT *out2)
{
    /* Send the signal's input to the output */
    uint32_t type = p->type;
    SPFLOAT pan = (p->pan + 1.0) * 0.5;
    SPFLOAT cc, ss, l, r;

    type %= 4;

    switch (type) {
        /* Equal power */
        case 0:
        pan = M_PI * 0.5 * pan;
        *out1 = *in1 * cos(pan);
        *out2 = *in2 * sin(pan);
        break;

        /* Square root */
        case 1:
        *out1 = *in1 * sqrt(pan);
        *out2 = *in2 * sqrt(1.0 - pan);
        break;

        /* simple linear */
        case 2:
        *out1 = *in1 * (1.0 - pan);
        *out2 = *in2 * pan;
        break;

        /* Equal power (alternative) */
        case 3:

        cc = cos(M_PI * pan * 0.5);
        ss = sin(M_PI * pan * 0.5);
        l = SQRT2 * (cc + ss) * 0.5;
        r = SQRT2 * (cc - ss) * 0.5;
        *out1 = *in1 * l;
        *out2 = *in2 * r;
        break;
    }

    return SP_OK;
}


int sp_pareq_create(sp_pareq **p)
{
    *p = malloc(sizeof(sp_pareq));
    return SP_OK;
}

int sp_pareq_destroy(sp_pareq **p)
{
    free(*p);
    return SP_OK;
}

int sp_pareq_init(sp_data *sp, sp_pareq *p)
{
    p->q = 0.707;
    p->v = 1;
    p->mode = 0;
    p->fc = 1000;

    p->xnm1 = p->xnm2 = p->ynm1 = p->ynm2 = 0.0;
    p->prv_fc = p->prv_v = p->prv_q = -1.0;
    p->tpidsr = (2 * M_PI) / sp->sr;
    return SP_OK;
}

int sp_pareq_compute(sp_data *sp, sp_pareq *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT xn, yn;
    SPFLOAT sq;

    if (p->fc != p->prv_fc || p->v != p->prv_v || p->q != p->prv_q) {
        SPFLOAT omega = (SPFLOAT)(p->tpidsr * p->fc), k, kk, vkk, vk, vkdq, a0;
        p->prv_fc = p->fc; p->prv_v = p->v; p->prv_q = p->q;
        switch ((int)p->mode) {
            /* Low Shelf */
            case 1: 
                sq = sqrt(2.0 * (SPFLOAT) p->prv_v);
                k = tan(omega * 0.5);
                kk = k * k;
                vkk = (SPFLOAT)p->prv_v * kk;
                p->b0 = 1.0 + sq * k + vkk;
                p->b1 = 2.0 * (vkk - 1.0);
                p->b2 = 1.0 - sq * k + vkk;
                a0 = 1.0 + k / (SPFLOAT)p->prv_q + kk;
                p->a1 = 2.0 * (kk - 1.0);
                p->a2 = 1.0 - k / (SPFLOAT)p->prv_q + kk;
                break;

            /* High Shelf */
            case 2: 
                sq = sqrt(2.0 * (SPFLOAT) p->prv_v);
                k = tan((M_PI - omega) * 0.5);
                kk = k * k;
                vkk = (SPFLOAT)p->prv_v * kk;
                p->b0 = 1.0 + sq * k + vkk;
                p->b1 = -2.0 * (vkk - 1.0);
                p->b2 = 1.0 - sq * k + vkk;
                a0 = 1.0 + k / (SPFLOAT)p->prv_q + kk;
                p->a1 = -2.0 * (kk - 1.0);
                p->a2 = 1.0 - k / (SPFLOAT)p->prv_q + kk;
                break;

            /* Peaking EQ */
            default: 
                k = tan(omega * 0.5);
                kk = k * k;
                vk = (SPFLOAT)p->prv_v * k;
                vkdq = vk / (SPFLOAT)p->prv_q;
                p->b0 = 1.0 + vkdq + kk;
                p->b1 = 2.0 * (kk - 1.0);
                p->b2 = 1.0 - vkdq + kk;
                a0 = 1.0 + k / (SPFLOAT)p->prv_q + kk;
                p->a1 = 2.0 * (kk - 1.0);
                p->a2 = 1.0 - k / (SPFLOAT)p->prv_q + kk;
        }
        a0 = 1.0 / a0;
        p->a1 *= a0; p->a2 *= a0; p->b0 *= a0; p->b1 *= a0; p->b2 *= a0;
    }
    {
        SPFLOAT a1 = p->a1, a2 = p->a2;
        SPFLOAT b0 = p->b0, b1 = p->b1, b2 = p->b2;
        SPFLOAT xnm1 = p->xnm1, xnm2 = p->xnm2, ynm1 = p->ynm1, ynm2 = p->ynm2;
        xn = *in;
        yn = b0 * xn + b1 * xnm1 + b2 * xnm2 - a1 * ynm1 - a2 * ynm2;
        xnm2 = xnm1;
        xnm1 = xn;
        ynm2 = ynm1;
        ynm1 = yn;
        *out = yn;
        p->xnm1 = xnm1; p->xnm2 = xnm2; p->ynm1 = ynm1; p->ynm2 = ynm2;
    }
    return SP_OK;
}



static void compute_block(sp_data *sp, sp_paulstretch *p) {
    uint32_t istart_pos = floor(p->start_pos);
    uint32_t pos; 
    uint32_t i;
    uint32_t windowsize = p->windowsize;
    uint32_t half_windowsize = p->half_windowsize;
    SPFLOAT *buf = p->buf;
    SPFLOAT *hinv_buf = p->hinv_buf;
    SPFLOAT *old_windowed_buf= p->old_windowed_buf;
    SPFLOAT *tbl = p->ft->tbl;
    SPFLOAT *window = p->window;
    SPFLOAT *output= p->output;
    for(i = 0; i < windowsize; i++) {
        /* Loop through buffer */
        pos = (istart_pos + i) % p->ft->size;

        if(p->wrap) {
            pos %= p->ft->size;
        } 

        if(pos < p->ft->size) {
            buf[i] = tbl[pos] * window[i];
        } else {
            buf[i] = 0;
        }
    }
    kiss_fftr(p->fft, buf, p->tmp1);
    for(i = 0; i < windowsize / 2; i++) {
        SPFLOAT mag = sqrt(p->tmp1[i].r*p->tmp1[i].r + p->tmp1[i].i*p->tmp1[i].i);
        SPFLOAT ph = ((SPFLOAT)sp_rand(sp) / SP_RANDMAX) * 2 * M_PI;
        p->tmp1[i].r = mag * cos(ph); 
        p->tmp1[i].i = mag * sin(ph); 
    }
    kiss_fftri(p->ifft, p->tmp1, buf);
    for(i = 0; i < windowsize; i++) {
        buf[i] *= window[i];
        if(i < half_windowsize) {
            output[i] = (SPFLOAT)(buf[i] + old_windowed_buf[half_windowsize + i]) / windowsize;
            output[i] *= hinv_buf[i];
        }
        old_windowed_buf[i] = buf[i];
    }
    p->start_pos += p->displace_pos;
}

int sp_paulstretch_create(sp_paulstretch **p)
{
    *p = malloc(sizeof(sp_paulstretch));
    return SP_OK;
}

int sp_paulstretch_destroy(sp_paulstretch **p)
{
    sp_paulstretch *pp = *p;
    sp_auxdata_free(&pp->m_window);
    sp_auxdata_free(&pp->m_old_windowed_buf);
    sp_auxdata_free(&pp->m_hinv_buf);
    sp_auxdata_free(&pp->m_buf);
    sp_auxdata_free(&pp->m_output);
    kiss_fftr_free(pp->fft);
    kiss_fftr_free(pp->ifft);
    KISS_FFT_FREE(pp->tmp1);
    free(*p);
    return SP_OK;
}

int sp_paulstretch_init(sp_data *sp, sp_paulstretch *p, sp_ftbl *ft, SPFLOAT windowsize, SPFLOAT stretch)
{
    uint32_t i;
    p->ft = ft;
    p->windowsize = (uint32_t)(sp->sr * windowsize);
    p->stretch = stretch;
    if(p->windowsize < 16) {
        p->windowsize = 16;
    }
    p->half_windowsize = p->windowsize / 2;
    p->displace_pos = (p->windowsize * 0.5) / p->stretch;

    sp_auxdata_alloc(&p->m_window, sizeof(SPFLOAT) * p->windowsize);
    p->window = p->m_window.ptr;

    sp_auxdata_alloc(&p->m_old_windowed_buf, sizeof(SPFLOAT) * p->windowsize);
    p->old_windowed_buf = p->m_old_windowed_buf.ptr;

    sp_auxdata_alloc(&p->m_hinv_buf, sizeof(SPFLOAT) * p->half_windowsize);
    p->hinv_buf = p->m_hinv_buf.ptr;

    sp_auxdata_alloc(&p->m_buf, sizeof(SPFLOAT) * p->windowsize);
    p->buf = p->m_buf.ptr;

    sp_auxdata_alloc(&p->m_output, sizeof(SPFLOAT) * p->half_windowsize);
    p->output = p->m_output.ptr;

    /* Create Hann window */
    for(i = 0; i < p->windowsize; i++) {
        p->window[i] = 0.5 - cos(i * 2.0 * M_PI / (p->windowsize - 1)) * 0.5;
    }
    /* creatve inverse hann window */
    SPFLOAT hinv_sqrt2 = (1 + sqrt(0.5)) * 0.5;
    for(i = 0; i < p->half_windowsize; i++) {
        p->hinv_buf[i] = hinv_sqrt2 - (1.0 - hinv_sqrt2) * cos(i * 2.0 * M_PI / p->half_windowsize);
    }

    p->start_pos = 0.0;
    p->counter = 0;

    /* set up kissfft */
    p->fft = kiss_fftr_alloc(p->windowsize, 0, NULL, NULL);
    p->ifft = kiss_fftr_alloc(p->windowsize, 1, NULL, NULL);
    kiss_fft_cpx *tmp1 = malloc(sizeof(kiss_fft_cpx) * p->windowsize);
    memset(tmp1, 0, sizeof(SPFLOAT) * p->windowsize);
    p->tmp1 = tmp1;

    /* turn on wrap mode by default */
    p->wrap = 1;
    return SP_OK;
}

int sp_paulstretch_compute(sp_data *sp, sp_paulstretch *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->counter == 0) {
        compute_block(sp, p);
    }
    *out = p->output[p->counter];
    p->counter = (p->counter + 1) % p->half_windowsize;
    
    return SP_OK;
}


int sp_pdhalf_create(sp_pdhalf **p)
{
    *p = malloc(sizeof(sp_pdhalf));
    return SP_OK;
}

int sp_pdhalf_destroy(sp_pdhalf **p)
{
    free(*p);
    return SP_OK;
}

int sp_pdhalf_init(sp_data *sp, sp_pdhalf *p)
{
    p->ibipolar = 0;
    p->ifullscale = 1.0;
    p->amount = 0;
    return SP_OK;
}

int sp_pdhalf_compute(sp_data *sp, sp_pdhalf *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT cur, maxampl, midpoint, leftslope, rightslope;

    maxampl = p->ifullscale;
    if (maxampl == 0.0)  maxampl = 1.0;

    if (p->ibipolar != 0.0) {
        midpoint =  (p->amount >= 1.0 ? maxampl :
                    (p->amount <= -1.0 ? -maxampl :
                    (p->amount * maxampl)));

    if (midpoint != -maxampl) 
        leftslope  = maxampl / (midpoint + maxampl);
    else leftslope  = 0.0;
    if (midpoint != maxampl)  
        rightslope = maxampl / (maxampl - midpoint);
    else rightslope = 0.0;

    cur = *in;
    if (cur < midpoint) *out = leftslope * (cur - midpoint);
    else *out = rightslope * (cur - midpoint);
    } else {
        SPFLOAT halfmaxampl = 0.5 * maxampl;
        midpoint =  (p->amount >= 1.0 ? maxampl :
                    (p->amount <= -1.0 ? 0.0 :
                    ((p->amount + 1.0) * halfmaxampl)));

        if (midpoint != 0.0) 
            leftslope = halfmaxampl / midpoint;
        else leftslope  = 0.0;
        if (midpoint != maxampl) 
            rightslope = halfmaxampl / (maxampl - midpoint);
        else rightslope = 0.0;

        cur = *in;
        if (cur < midpoint) { 
            *out = leftslope * cur;
        } else { 
            *out = rightslope * (cur - midpoint) + halfmaxampl;
        }
    }

    return SP_OK;
}



int sp_peaklim_create(sp_peaklim **p)
{
    *p = malloc(sizeof(sp_peaklim));
    return SP_OK;
}

int sp_peaklim_destroy(sp_peaklim **p)
{
    free(*p);
    return SP_OK;
}

int sp_peaklim_init(sp_data *sp, sp_peaklim *p)
{
    p->a1_r = 0;
    p->b0_r = 1;
    p->a1_a = 0;
    p->b0_a = 1;
    p->atk = 0.1;
    p->rel = 0.1;
    p->patk = -100;
    p->prel = -100;
    p->level = 0;
    return SP_OK;
}

int sp_peaklim_compute(sp_data *sp, sp_peaklim *p, SPFLOAT *in, SPFLOAT *out)
{

    SPFLOAT db_gain = 0;
    SPFLOAT gain = 0;

    /* change coefficients, if needed */

    if(p->patk != p->atk) {
        p->patk = p->atk;
		p->a1_a = exp( -1.0 / ( p->rel * sp->sr ) );
		p->b0_a = 1 - p->a1_a;
    }

    if(p->prel != p->rel) {
        p->prel = p->rel;
		p->a1_r = exp( -1.0 / ( p->rel * sp->sr ) );
		p->b0_r = 1 - p->a1_r;
    }

    
    if ( fabs(*in) > p->level)
        p->level += p->b0_a * ( fabs(*in) - p->level);
    else
        p->level += p->b0_r * ( fabs(*in) - p->level);

    db_gain = min(0.0, dB(dB2lin(p->thresh)/p->level));
    gain = dB2lin(db_gain);		

    *out = *in * gain;

    return SP_OK;
}

static float faustpower3_f(float value) {
	return ((value * value) * value);
	
}
static float faustpower4_f(float value) {
	return (((value * value) * value) * value);
	
}

typedef struct {
	
	float fRec4[3];
	float fRec3[3];
	float fRec2[3];
	float fRec1[3];
	float fRec11[3];
	float fRec10[3];
	float fRec9[3];
	float fRec8[3];
	int iVec0[2];
	float fRec5[2];
	float fRec6[2];
	float fRec0[2];
	float fRec7[2];
	FAUSTFLOAT fHslider0;
	FAUSTFLOAT fCheckbox0;
	FAUSTFLOAT fHslider1;
	int fSamplingFreq;
	int iConst0;
	float fConst1;
	FAUSTFLOAT fHslider2;
	FAUSTFLOAT fHslider3;
	FAUSTFLOAT fHslider4;
	FAUSTFLOAT fHslider5;
	float fConst2;
	FAUSTFLOAT fHslider6;
	FAUSTFLOAT fHslider7;
	FAUSTFLOAT fCheckbox1;
	
} phaser;

phaser* newphaser() { 
	phaser* dsp = (phaser*)malloc(sizeof(phaser));
	return dsp;
}

void deletephaser(phaser* dsp) { 
	free(dsp);
}

void instanceInitphaser(phaser* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->fHslider0 = (FAUSTFLOAT)0.;
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->iVec0[i0] = 0;
			
		}
		
	}
	dsp->fCheckbox0 = (FAUSTFLOAT)0.;
	dsp->fHslider1 = (FAUSTFLOAT)1.;
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fConst1 = (1.f / (float)dsp->iConst0);
	dsp->fHslider2 = (FAUSTFLOAT)1000.;
	dsp->fHslider3 = (FAUSTFLOAT)1.5;
	dsp->fHslider4 = (FAUSTFLOAT)100.;
	dsp->fHslider5 = (FAUSTFLOAT)800.;
	dsp->fConst2 = (0.10472f / (float)dsp->iConst0);
	dsp->fHslider6 = (FAUSTFLOAT)30.;
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec5[i1] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fRec6[i2] = 0.f;
			
		}
		
	}
	dsp->fHslider7 = (FAUSTFLOAT)0.;
	/* C99 loop */
	{
		int i3;
		for (i3 = 0; (i3 < 3); i3 = (i3 + 1)) {
			dsp->fRec4[i3] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i4;
		for (i4 = 0; (i4 < 3); i4 = (i4 + 1)) {
			dsp->fRec3[i4] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i5;
		for (i5 = 0; (i5 < 3); i5 = (i5 + 1)) {
			dsp->fRec2[i5] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i6;
		for (i6 = 0; (i6 < 3); i6 = (i6 + 1)) {
			dsp->fRec1[i6] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i7;
		for (i7 = 0; (i7 < 2); i7 = (i7 + 1)) {
			dsp->fRec0[i7] = 0.f;
			
		}
		
	}
	dsp->fCheckbox1 = (FAUSTFLOAT)0.;
	/* C99 loop */
	{
		int i8;
		for (i8 = 0; (i8 < 3); i8 = (i8 + 1)) {
			dsp->fRec11[i8] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i9;
		for (i9 = 0; (i9 < 3); i9 = (i9 + 1)) {
			dsp->fRec10[i9] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i10;
		for (i10 = 0; (i10 < 3); i10 = (i10 + 1)) {
			dsp->fRec9[i10] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i11;
		for (i11 = 0; (i11 < 3); i11 = (i11 + 1)) {
			dsp->fRec8[i11] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i12;
		for (i12 = 0; (i12 < 2); i12 = (i12 + 1)) {
			dsp->fRec7[i12] = 0.f;
			
		}
		
	}
	
}

void initphaser(phaser* dsp, int samplingFreq) {
	instanceInitphaser(dsp, samplingFreq);
}

void buildUserInterfacephaser(phaser* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "MaxNotch1Freq", &dsp->fHslider5, 800.f, 20.f, 10000.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "MinNotch1Freq", &dsp->fHslider4, 100.f, 20.f, 5000.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "Notch width", &dsp->fHslider2, 1000.f, 10.f, 5000.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "NotchFreq", &dsp->fHslider3, 1.5f, 1.1f, 4.f, 0.01f);
	interface->addCheckButton(interface->uiInterface, "VibratoMode", &dsp->fCheckbox0);
	interface->addHorizontalSlider(interface->uiInterface, "depth", &dsp->fHslider1, 1.f, 0.f, 1.f, 0.01f);
	interface->addHorizontalSlider(interface->uiInterface, "feedback gain", &dsp->fHslider7, 0.f, 0.f, 1.f, 0.01f);
	interface->addCheckButton(interface->uiInterface, "invert", &dsp->fCheckbox1);
	interface->addHorizontalSlider(interface->uiInterface, "level", &dsp->fHslider0, 0.f, -60.f, 10.f, 0.1f);
	interface->addHorizontalSlider(interface->uiInterface, "lfobpm", &dsp->fHslider6, 30.f, 24.f, 360.f, 1.f);
}

void computephaser(phaser* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* input0 = inputs[0];
	FAUSTFLOAT* input1 = inputs[1];
	FAUSTFLOAT* output0 = outputs[0];
	FAUSTFLOAT* output1 = outputs[1];
	float fSlow0 = pow(10.f, (0.05f * (float)dsp->fHslider0));
	float fSlow1 = (0.5f * ((int)(float)dsp->fCheckbox0?2.f:(float)dsp->fHslider1));
	float fSlow2 = (1.f - fSlow1);
	float fSlow3 = exp((dsp->fConst1 * (0.f - (3.14159f * (float)dsp->fHslider2))));
	float fSlow4 = faustpower2_f(fSlow3);
	float fSlow5 = (0.f - (2.f * fSlow3));
	float fSlow6 = (float)dsp->fHslider3;
	float fSlow7 = (dsp->fConst1 * fSlow6);
	float fSlow8 = (float)dsp->fHslider4;
	float fSlow9 = (6.28319f * fSlow8);
	float fSlow10 = (0.5f * ((6.28319f * max(fSlow8, (float)dsp->fHslider5)) - fSlow9));
	float fSlow11 = (dsp->fConst2 * (float)dsp->fHslider6);
	float fSlow12 = sin(fSlow11);
	float fSlow13 = cos(fSlow11);
	float fSlow14 = (0.f - fSlow12);
	float fSlow15 = (float)dsp->fHslider7;
	float fSlow16 = (dsp->fConst1 * faustpower2_f(fSlow6));
	float fSlow17 = (dsp->fConst1 * faustpower3_f(fSlow6));
	float fSlow18 = (dsp->fConst1 * faustpower4_f(fSlow6));
	float fSlow19 = ((int)(float)dsp->fCheckbox1?(0.f - fSlow1):fSlow1);
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			dsp->iVec0[0] = 1;
			float fTemp0 = (float)input0[i];
			dsp->fRec5[0] = ((fSlow12 * dsp->fRec6[1]) + (fSlow13 * dsp->fRec5[1]));
			dsp->fRec6[0] = ((1.f + ((fSlow13 * dsp->fRec6[1]) + (fSlow14 * dsp->fRec5[1]))) - (float)dsp->iVec0[1]);
			float fTemp1 = ((fSlow10 * (1.f - dsp->fRec5[0])) + fSlow9);
			float fTemp2 = (dsp->fRec4[1] * cos((fSlow7 * fTemp1)));
			dsp->fRec4[0] = (0.f - (((fSlow5 * fTemp2) + (fSlow4 * dsp->fRec4[2])) - ((fSlow0 * fTemp0) + (fSlow15 * dsp->fRec0[1]))));
			float fTemp3 = (dsp->fRec3[1] * cos((fSlow16 * fTemp1)));
			dsp->fRec3[0] = ((fSlow5 * (fTemp2 - fTemp3)) + (dsp->fRec4[2] + (fSlow4 * (dsp->fRec4[0] - dsp->fRec3[2]))));
			float fTemp4 = (dsp->fRec2[1] * cos((fSlow17 * fTemp1)));
			dsp->fRec2[0] = ((fSlow5 * (fTemp3 - fTemp4)) + (dsp->fRec3[2] + (fSlow4 * (dsp->fRec3[0] - dsp->fRec2[2]))));
			float fTemp5 = (dsp->fRec1[1] * cos((fSlow18 * fTemp1)));
			dsp->fRec1[0] = ((fSlow5 * (fTemp4 - fTemp5)) + (dsp->fRec2[2] + (fSlow4 * (dsp->fRec2[0] - dsp->fRec1[2]))));
			dsp->fRec0[0] = ((fSlow4 * dsp->fRec1[0]) + ((fSlow5 * fTemp5) + dsp->fRec1[2]));
			output0[i] = (FAUSTFLOAT)((fSlow0 * (fSlow2 * fTemp0)) + (dsp->fRec0[0] * fSlow19));
			float fTemp6 = (float)input1[i];
			float fTemp7 = ((fSlow10 * (1.f - dsp->fRec6[0])) + fSlow9);
			float fTemp8 = (dsp->fRec11[1] * cos((fSlow7 * fTemp7)));
			dsp->fRec11[0] = (0.f - (((fSlow5 * fTemp8) + (fSlow4 * dsp->fRec11[2])) - ((fSlow0 * fTemp6) + (fSlow15 * dsp->fRec7[1]))));
			float fTemp9 = (dsp->fRec10[1] * cos((fSlow16 * fTemp7)));
			dsp->fRec10[0] = ((fSlow5 * (fTemp8 - fTemp9)) + (dsp->fRec11[2] + (fSlow4 * (dsp->fRec11[0] - dsp->fRec10[2]))));
			float fTemp10 = (dsp->fRec9[1] * cos((fSlow17 * fTemp7)));
			dsp->fRec9[0] = ((fSlow5 * (fTemp9 - fTemp10)) + (dsp->fRec10[2] + (fSlow4 * (dsp->fRec10[0] - dsp->fRec9[2]))));
			float fTemp11 = (dsp->fRec8[1] * cos((fSlow18 * fTemp7)));
			dsp->fRec8[0] = ((fSlow5 * (fTemp10 - fTemp11)) + (dsp->fRec9[2] + (fSlow4 * (dsp->fRec9[0] - dsp->fRec8[2]))));
			dsp->fRec7[0] = ((fSlow4 * dsp->fRec8[0]) + ((fSlow5 * fTemp11) + dsp->fRec8[2]));
			output1[i] = (FAUSTFLOAT)((fSlow0 * (fSlow2 * fTemp6)) + (dsp->fRec7[0] * fSlow19));
			dsp->iVec0[1] = dsp->iVec0[0];
			dsp->fRec5[1] = dsp->fRec5[0];
			dsp->fRec6[1] = dsp->fRec6[0];
			dsp->fRec4[2] = dsp->fRec4[1];
			dsp->fRec4[1] = dsp->fRec4[0];
			dsp->fRec3[2] = dsp->fRec3[1];
			dsp->fRec3[1] = dsp->fRec3[0];
			dsp->fRec2[2] = dsp->fRec2[1];
			dsp->fRec2[1] = dsp->fRec2[0];
			dsp->fRec1[2] = dsp->fRec1[1];
			dsp->fRec1[1] = dsp->fRec1[0];
			dsp->fRec0[1] = dsp->fRec0[0];
			dsp->fRec11[2] = dsp->fRec11[1];
			dsp->fRec11[1] = dsp->fRec11[0];
			dsp->fRec10[2] = dsp->fRec10[1];
			dsp->fRec10[1] = dsp->fRec10[0];
			dsp->fRec9[2] = dsp->fRec9[1];
			dsp->fRec9[1] = dsp->fRec9[0];
			dsp->fRec8[2] = dsp->fRec8[1];
			dsp->fRec8[1] = dsp->fRec8[0];
			dsp->fRec7[1] = dsp->fRec7[0];
			
		}
		
	}
	
}



int sp_phaser_create(sp_phaser **p)
{
    *p = malloc(sizeof(sp_phaser));
    return SP_OK;
}

int sp_phaser_destroy(sp_phaser **p)
{
    sp_phaser *pp = *p;
    phaser *dsp = pp->faust;
    deletephaser (dsp);
    free(*p);
    return SP_OK;
}

int sp_phaser_init(sp_data *sp, sp_phaser *p)
{
    phaser *dsp = newphaser(); 
    UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.addCheckButton = addCheckButton;
    UI.uiInterface = p;
    buildUserInterfacephaser(dsp, &UI);
    initphaser(dsp, sp->sr);

     
    p->MaxNotch1Freq = p->args[0]; 
    p->MinNotch1Freq = p->args[1]; 
    p->Notch_width = p->args[2]; 
    p->NotchFreq = p->args[3]; 
    p->VibratoMode = p->args[4]; 
    p->depth = p->args[5]; 
    p->feedback_gain = p->args[6]; 
    p->invert = p->args[7]; 
    p->level = p->args[8]; 
    p->lfobpm = p->args[9];

    p->faust = dsp;
    return SP_OK;
}

int sp_phaser_compute(sp_data *sp, sp_phaser *p, 
	SPFLOAT *in1, SPFLOAT *in2, SPFLOAT *out1, SPFLOAT *out2) 
{
    phaser *dsp = p->faust;
    SPFLOAT *faust_out[] = {out1, out2};
    SPFLOAT *faust_in[] = {in1, in2};
    computephaser(dsp, 1, faust_in, faust_out);
    return SP_OK;
}


int sp_phasor_create(sp_phasor **p)
{
    *p = malloc(sizeof(sp_phasor));
    return SP_OK;
}

int sp_phasor_destroy(sp_phasor **p)
{
    free(*p);
    return SP_OK;
}

int sp_phasor_init(sp_data *sp, sp_phasor *p, SPFLOAT iphs)
{
    p->freq = 440;
    p->phs = iphs;
    p->curphs = iphs;
    p->onedsr = 1.0 / sp->sr;
    return SP_OK;
}

int sp_phasor_compute(sp_data *sp, sp_phasor *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT phase;
    SPFLOAT incr;

    phase = p->curphs;
    incr = p->freq * p->onedsr;
    *out = phase;
    phase += incr;
    if (phase >= 1.0) {
        phase -= 1.0;
    } else if (phase < 0.0) {
        phase += 1.0;
    }
    p->curphs = phase;
    return SP_OK;
}


static uint32_t ctz[64] =
{
    6, 0, 1, 0, 2, 0, 1, 0,
    3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0,
    3, 0, 1, 0, 2, 0, 1, 0,
    5, 0, 1, 0, 2, 0, 1, 0,
    3, 0, 1, 0, 2, 0, 1, 0,
    4, 0, 1, 0, 2, 0, 1, 0,
    3, 0, 1, 0, 2, 0, 1, 0,
};

int sp_pinknoise_create(sp_pinknoise **p)
{
    *p = malloc(sizeof(sp_pinknoise));
    return SP_OK;
}

int sp_pinknoise_destroy(sp_pinknoise **p)
{
    free(*p);
    return SP_OK;
}

int sp_pinknoise_init(sp_data *sp, sp_pinknoise *p)
{
    int i;
    p->amp = 1.0;
    p->seed = sp_rand(sp);
    p->total = 0;
    p->counter = 0;
    for(i = 0; i < 7; i++) {
        p->dice[i] = 0;
    }
    return SP_OK;
}

int sp_pinknoise_compute(sp_data *sp, sp_pinknoise *p, SPFLOAT *in, SPFLOAT *out) 
{
    uint32_t k = ctz[p->counter & 63];
    p->prevrand = p->dice[k];
    p->seed = 1664525 * p->seed + 1013904223;
    p->newrand = p->seed >> 3;
    p->dice[k] = p->newrand;
    p->total += (p->newrand - p->prevrand);
    p->seed = 1103515245 * p->seed + 12345;
    p->newrand = p->seed >> 3;
    short tmp = (short) ((((p->total + p->newrand) * (1.0f / (3 << 29)) - 1) - .25f) * 16384.0f);
    
    *out = ((SPFLOAT) tmp / 32767) * p->amp;
    p->counter = (p->counter + 1) % 0xFFFFFFFF;
    return SP_OK;
}


/* #define lrintf(x) lrintf(x) */

int sp_pitchamdf_create(sp_pitchamdf **p)
{
    *p = malloc(sizeof(sp_pitchamdf));
    return SP_OK;
}

int sp_pitchamdf_destroy(sp_pitchamdf **p)
{
    sp_pitchamdf *pp = *p;
    sp_auxdata_free(&pp->median);
/* This mirrors the original code */
    if(pp->rmsmedisize) {
        sp_auxdata_free(&pp->rmsmedian);
    }
    sp_auxdata_free(&pp->buffer);
    free(*p);
    return SP_OK;
}

int sp_pitchamdf_init(sp_data *sp, sp_pitchamdf *p, SPFLOAT imincps, SPFLOAT imaxcps)
{
    SPFLOAT srate, downs;
    int32_t size, minperi, maxperi, downsamp, upsamp, msize, bufsize;
    uint32_t interval;

    p->imincps = imincps;
    p->imaxcps = imaxcps;

    /* TODO: should we expose these variables? */
    p->icps = 0;
    p->imedi = 1;
    p->idowns = 1;
    p->iexcps = 0;
    p->irmsmedi = 0;

    p->inerr = 0;
    downs = p->idowns;

    if (downs < (-1.9)) {
        upsamp = (int)lrintf((-downs));
        downsamp = 0;
        srate = sp->sr * (SPFLOAT)upsamp;
    } else {
        downsamp = (int)lrintf(downs);
        if (downsamp < 1) downsamp = 1;
        srate = sp->sr / (SPFLOAT)downsamp;
        upsamp = 0;
    }

    minperi = (int32_t)(srate / p->imaxcps);
    maxperi = (int32_t)(0.5 + srate / p->imincps);
    if (maxperi <= minperi) {
        p->inerr = 1;
        return SP_NOT_OK;
    }

    if (p->iexcps < 1)
        interval = maxperi;
    else
        interval = (uint32_t)(srate / p->iexcps);

    size = maxperi + interval;
    bufsize = sizeof(SPFLOAT)*(size + maxperi + 2);

    p->srate = srate;
    p->downsamp = downsamp;
    p->upsamp = upsamp;
    p->minperi = minperi;
    p->maxperi = maxperi;
    p->size = size;
    p->readp = 0;
    p->index = 0;
    p->lastval = 0.0;

    if (p->icps < 1) {
        p->peri = (minperi + maxperi) / 2;
    } else {
        p->peri = (int)(srate / p->icps);
    }

    if (p->irmsmedi < 1) {
        p->rmsmedisize = 0;
    } else {
        p->rmsmedisize = ((int)lrintf(p->irmsmedi))*2+1;
    }

    p->rmsmediptr = 0;

    if (p->rmsmedisize) {
        msize = p->rmsmedisize * 3 * sizeof(SPFLOAT);
        sp_auxdata_alloc(&p->rmsmedian, msize);
    }

    if (p->imedi < 1) {
        p->medisize = 0;
    } else {
        p->medisize = (int)lrintf(p->imedi) * 2 + 1;
    }

    p->mediptr = 0;

    if (p->medisize) {
        msize = p->medisize * 3 * sizeof(SPFLOAT);
        sp_auxdata_alloc(&p->median, msize);
    }

    sp_auxdata_alloc(&p->buffer, bufsize);
    return SP_OK;
}


#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp

static SPFLOAT medianvalue(uint32_t n, SPFLOAT *vals)
{   
    /* vals must point to 1 below relevant data! */
    uint32_t i, ir, j, l, mid;
    uint32_t k = (n + 1) / 2;
    SPFLOAT a, temp;

    l = 1;
    ir = n;
    while (1) {
        if (ir <= l+1) {
            if (ir == l+1 && vals[ir] < vals[l]) {
                SWAP(vals[l], vals[ir]);
            }
            return vals[k];
        } else {
            mid = (l+ir) >> 1;
            SWAP(vals[mid], vals[l+1]);
            if (vals[l+1] > vals[ir]) {
                SWAP(vals[l+1], vals[ir]);
            }
            if (vals[l] > vals[ir]) {
                SWAP(vals[l], vals[ir]);
            }
            if (vals[l+1] > vals[l]) {
                SWAP(vals[l+1], vals[l]);
            }
            i = l + 1;
            j = ir;
            a = vals[l];
            while (1) {
                do i++; while (vals[i] < a);
                do j--; while (vals[j] > a);
                if (j < i) break;
                SWAP(vals[i], vals[j]);
            }
            vals[l] = vals[j];
            vals[j] = a;
            if (j >= k) ir = j-1;
            if (j <= k) l = i;
        }
    }
}
#undef SWAP

int sp_pitchamdf_compute(sp_data *sp, sp_pitchamdf *p, SPFLOAT *in, 
    SPFLOAT *cps, SPFLOAT *rms_out)
{
    SPFLOAT *buffer = (SPFLOAT*)p->buffer.ptr;
    SPFLOAT *rmsmedian = (SPFLOAT*)p->rmsmedian.ptr;
    int32_t rmsmedisize = p->rmsmedisize;
    int32_t rmsmediptr = p->rmsmediptr;
    SPFLOAT *median = (SPFLOAT*)p->median.ptr;
    int32_t medisize = p->medisize;
    int32_t mediptr = p->mediptr;
    int32_t size = p->size;
    int32_t index = p->index;
    int32_t minperi = p->minperi;
    int32_t maxperi = p->maxperi;
    SPFLOAT srate = p->srate;
    int32_t peri = p->peri;
    int32_t upsamp = p->upsamp;
    SPFLOAT upsmp = (SPFLOAT)upsamp;
    SPFLOAT lastval = p->lastval;
    SPFLOAT newval, delta;
    int32_t readp = p->readp;
    int32_t interval = size - maxperi;
    int i;
    int32_t i1, i2;
    SPFLOAT val, rms;
    SPFLOAT sum;
    SPFLOAT acc, accmin, diff;

    if (upsamp) {
        newval = *in;
        delta = (newval-lastval) / upsmp;
        lastval = newval;

        for (i=0; i<upsamp; i++) {
            newval += delta;
            buffer[index++] = newval;

            if (index == size) {
                peri = minperi;
                accmin = 0.0;
                for (i2 = 0; i2 < size; ++i2) {
                    diff = buffer[i2+minperi] - buffer[i2];
                    if (diff > 0) accmin += diff;
                    else accmin -= diff;
                }
                for (i1 = minperi + 1; i1 <= maxperi; ++i1) {
                    acc = 0.0;
                    for (i2 = 0; i2 < size; ++i2) {
                        diff = buffer[i1+i2] - buffer[i2];
                        if (diff > 0) acc += diff;
                        else acc -= diff;
                    }
                    if (acc < accmin) {
                        accmin = acc;
                        peri = i1;
                    }
                }

                for (i1 = 0; i1 < interval; i1++) { 
                    buffer[i1] = buffer[i1+interval]; 
                }

                index = maxperi;

                if (medisize) {
                    median[mediptr] = (SPFLOAT)peri;
                    for (i1 = 0; i1 < medisize; i1++) {
                        median[medisize+i1] = median[i1];
                    }

                    median[medisize*2+mediptr] =
                    medianvalue(medisize, &median[medisize-1]);
                    peri = (int32_t)median[medisize*2 +
                        ((mediptr+medisize/2+1) % medisize)];

                    mediptr = (mediptr + 1) % medisize;
                    p->mediptr = mediptr;
                }
            }
        }
        p->lastval = lastval;
    } else {
        int32_t  downsamp = p->downsamp;
        buffer[index++] = *in;
        readp += downsamp;

        if (index == size) {
            peri = minperi;
            accmin = 0.0;

            for (i2 = 0; i2 < size; ++i2) {
                diff = buffer[i2+minperi] - buffer[i2];
                if (diff > 0.0) accmin += diff;
                else accmin -= diff;
            }

            for (i1 = minperi + 1; i1 <= maxperi; ++i1) {
                acc = 0.0;
                for (i2 = 0; i2 < size; ++i2) {
                    diff = buffer[i1+i2] - buffer[i2];
                    if (diff > 0.0) acc += diff;
                    else acc -= diff;
                }
                if (acc < accmin) {
                    accmin = acc;
                    peri = i1;
                }
            }

            for (i1 = 0; i1 < interval; i1++) {
                buffer[i1] = buffer[i1+interval];
            }

            index = maxperi;

            if (medisize) {
                median[mediptr] = (SPFLOAT)peri;

                for (i1 = 0; i1 < medisize; i1++) {
                    median[medisize+i1] = median[i1];
                }

                median[medisize*2+mediptr] =
                medianvalue(medisize, &median[medisize-1]);
                peri = (int32_t)median[medisize*2 +
                    ((mediptr+medisize/2+1) % medisize)];

                mediptr = (mediptr + 1) % medisize;
                p->mediptr = mediptr;
            }
        }
    }
    buffer = &buffer[(index + size - peri) % size];
    sum = 0.0;
    for (i1=0; i1<peri; i1++) {
        val = buffer[i1];
        sum += (SPFLOAT)(val * val);
    }
    if (peri==0)      
        rms = 0.0;
    else
        rms = (SPFLOAT)sqrt(sum / (SPFLOAT)peri);
    if (rmsmedisize) {
        rmsmedian[rmsmediptr] = rms;
        for (i1 = 0; i1 < rmsmedisize; i1++) {
            rmsmedian[rmsmedisize+i1] = rmsmedian[i1];
        }

        rmsmedian[rmsmedisize*2+rmsmediptr] =
            medianvalue(rmsmedisize, &rmsmedian[rmsmedisize-1]);
        rms = rmsmedian[rmsmedisize*2 +
            ((rmsmediptr+rmsmedisize/2+1) % rmsmedisize)];

        rmsmediptr = (rmsmediptr + 1) % rmsmedisize;
        p->rmsmediptr = rmsmediptr;
    }

    if (peri==0) {
        *cps = 0.0;
    } else {
        *cps = srate / (SPFLOAT)peri;
    }

    *rms_out = rms;
    p->index = index;
    p->peri = peri;
    p->readp = readp;
    return SP_OK;
}




 
int sp_pluck_create(sp_pluck **p)
{
    *p = malloc(sizeof(sp_pluck));
    return SP_OK;
}

int sp_pluck_destroy(sp_pluck **p)
{
    sp_pluck *pp = *p;
    sp_auxdata_free(&pp->auxch);
    free(*p);
    return SP_OK;
}

static void sp_pluck_reinit(sp_data *sp, sp_pluck *p)
{
    int n;
    SPFLOAT val = 0;
    SPFLOAT *ap = (SPFLOAT *)p->auxch.ptr;
    for (n=p->npts; n--; ) {   
        val = (SPFLOAT) ((SPFLOAT) sp_rand(sp) / SP_RANDMAX);
        *ap++ = (val * 2) - 1;
    }
    p->phs256 = 0;
}

int sp_pluck_init(sp_data *sp, sp_pluck *p, SPFLOAT ifreq)
{
    int32_t npts;

    p->amp = 0.5;
    p->ifreq = ifreq;
    p->freq = ifreq;

    if ((npts = (int32_t)(sp->sr / p->ifreq)) < PLUKMIN) {
        npts = PLUKMIN;                  
    }
    
    sp_auxdata_alloc(&p->auxch, (npts + 1) * sizeof(SPFLOAT));
    p->maxpts = npts;
    p->npts = npts;

    sp_pluck_reinit(sp, p);
    /* tuned pitch convt */
    p->sicps = (npts * 256.0 + 128.0) * (1.0 / sp->sr);
    p->init = 1;
    return SP_OK;
}

int sp_pluck_compute(sp_data *sp, sp_pluck *p, SPFLOAT *trig, SPFLOAT *out)
{
    SPFLOAT *fp;
    int32_t phs256, phsinc, ltwopi, offset;
    SPFLOAT frac, diff;


    if(*trig != 0) {
        p->init = 0;
        sp_pluck_reinit(sp, p);
    }

    if(p->init) {
        *out = 0;
        return SP_OK;
    }

    phsinc = (int32_t)(p->freq * p->sicps);
    phs256 = p->phs256;
    ltwopi = p->npts << 8;
    offset = phs256 >> 8;
    fp = (SPFLOAT *)p->auxch.ptr + offset;     /* lookup position   */
    diff = fp[1] - fp[0];
    frac = (SPFLOAT)(phs256 & 255) / 256.0; /*  w. interpolation */
    *out = (fp[0] + diff*frac) * p->amp; /*  gives output val */
    if ((phs256 += phsinc) >= ltwopi) {
        int nn;
        SPFLOAT preval;
        phs256 -= ltwopi;               
        fp=(SPFLOAT *)p->auxch.ptr;
        preval = fp[0];                
        fp[0] = fp[p->npts];
        fp++;
        nn = p->npts;
        do {          
            /* 1st order recursive filter*/
            preval = (*fp + preval) * 0.5;
            *fp++ = preval;
        } while (--nn);
    }
    p->phs256 = phs256;
    return SP_OK;
}


int sp_port_create(sp_port **p)
{
    *p = malloc(sizeof(sp_port));
    return SP_OK;
}

int sp_port_destroy(sp_port **p)
{
    free(*p);
    return SP_OK;
}

int sp_port_init(sp_data *sp, sp_port *p, SPFLOAT htime)
{
    p->yt1 = 0;
    p->prvhtim = -100.0;
    p->htime = htime;

    p->sr = sp->sr;
    p->onedsr = 1.0/p->sr;
    return SP_OK;
}

int sp_port_compute(sp_data *sp, sp_port *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->prvhtim != p->htime){
        p->c2 = pow(0.5, p->onedsr / p->htime);
        p->c1 = 1.0 - p->c2;
        p->prvhtim = p->htime;
    }

    *out = p->yt1 = p->c1 * *in + p->c2 * p->yt1;
    return SP_OK;
}

int sp_port_reset(sp_data *sp, sp_port *p, SPFLOAT *in)
{
    p->yt1 = *in;
    return SP_OK;
}


int sp_posc3_create(sp_posc3 **posc3)
{
    *posc3 = malloc(sizeof(sp_posc3));
    return SP_OK;
}

int sp_posc3_destroy(sp_posc3 **posc3)
{
    free(*posc3);
    return SP_NOT_OK;
}

int sp_posc3_init(sp_data *sp, sp_posc3 *posc3, sp_ftbl *ft)
{

    posc3->amp = 0.2;
    posc3->freq = 440.0;
    posc3->iphs = 0.0;
    posc3->onedsr = 1.0 / sp->sr;

    posc3->tbl = ft;
    posc3->tablen = (int32_t) ft->size;
    posc3->tablenUPsr = posc3->tablen * posc3->onedsr;
    posc3->phs = posc3->iphs * posc3->tablen;
    return SP_OK;
}

int sp_posc3_compute(sp_data *sp, sp_posc3 *posc3, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT *ftab;
    SPFLOAT fract;
    SPFLOAT phs  = posc3->phs;
    SPFLOAT si   = posc3->freq * posc3->tablen * posc3->onedsr;
    SPFLOAT amp = posc3->amp;
    int x0;
    SPFLOAT y0, y1, ym1, y2;

    ftab = posc3->tbl->tbl;

    x0    = (int32_t )phs;
    fract = (SPFLOAT)(phs - (SPFLOAT)x0);
    x0--;

    if (x0<0) {
        ym1 = ftab[posc3->tablen-1]; x0 = 0;
    }
    else ym1 = ftab[x0++];
    y0    = ftab[x0++];
    y1    = ftab[x0++];
    if (x0>posc3->tablen) y2 = ftab[1];
    else y2 = ftab[x0];
    {
        SPFLOAT frsq = fract*fract;
        SPFLOAT frcu = frsq*ym1;
        SPFLOAT t1   = y2 + y0+y0+y0;
        *out     = amp * (y0 + 0.5 *frcu +
        fract*(y1 - frcu/6.0 - t1/6.0
        - ym1/3.0) +
        frsq*fract*(t1/6.0 - 0.5*y1) +
        frsq*(0.5* y1 - y0));
    }
    phs += si;
    while (phs >= posc3->tablen) {
        phs -= posc3->tablen;
    }
    while (phs < 0.0) {
        phs += posc3->tablen;
        posc3->phs = phs;
    }
    posc3->phs = phs;
    return SP_OK;
}


int sp_progress_create(sp_progress **p)
{
    *p = malloc(sizeof(sp_progress));
    return SP_OK;
}

int sp_progress_destroy(sp_progress **p)
{
    free(*p);
    return SP_OK;
}

int sp_progress_init(sp_data *sp, sp_progress *p)
{
    p->nbars = 40;
    p->skip = 1000;
    p->counter = 0;
    p->len = (uint32_t) sp->len;
    return SP_OK;
}

int sp_progress_compute(sp_data *sp, sp_progress *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->counter == 0 || sp->pos == p->len - 1) {
        int n;
        SPFLOAT slope = 1.0 / p->nbars;
        if(sp->pos == 0) fprintf(stderr, "\e[?25l");
        SPFLOAT percent = ((SPFLOAT)sp->pos / p->len);
        fprintf(stderr, "[");
        for(n = 0; n < p->nbars; n++) {
            if(n * slope <= percent) {
                fprintf(stderr, "#");
            }else {
                fprintf(stderr, " ");
            }
        }
        fprintf(stderr, "] %.2f%%\t\r", 100 * percent);

    }
    if(sp->pos == p->len - 1) fprintf(stderr, "\n\e[?25h");
    fflush(stderr);
    p->counter++;
    p->counter %= p->skip;
    return SP_OK;
}


typedef struct {
    uint32_t size;
    prop_list **ar;
} prop_slice;

static int prop_create(prop_data **pd);
static int prop_parse(prop_data *pd, const char *str);
static prop_event prop_next(sp_data *sp, prop_data *pd);
static float prop_time(prop_data *pd, prop_event evt);
static int prop_destroy(prop_data **pd);

static int prop_val_free(prop_val val);
static int prop_list_init(prop_list *lst);
static int prop_list_destroy(prop_list *lst);
static int prop_list_append(prop_list *lst, prop_val val);
static void prop_list_reset(prop_list *lst);
static int prop_list_copy(prop_list *src, prop_list **dst);

static void mode_insert_event(prop_data *pd, char type);
static void mode_insert_slice(prop_data *pd);
static void mode_list_start(prop_data *pd);
static void mode_list_end(prop_data *pd);
static void prop_slice_encap(prop_data *pd);
static void prop_slice_append(prop_data *pd);
static void reset(prop_data *pd);
static void back_to_top(prop_data *pd);

enum {
PTYPE_SLICE,
PTYPE_LIST,
PTYPE_EVENT,
PTYPE_OFF,
PTYPE_ON,
PTYPE_MAYBE,
PMODE_INSERT,
PMODE_SETDIV,
PMODE_SETMUL,
PMODE_UNSETMUL,
PMODE_INIT,
PSTATUS_NOTOK,
PSTATUS_OK,
PTYPE_NULL
};

int sp_prop_create(sp_prop **p)
{
    *p = malloc(sizeof(sp_prop));
    return SP_OK;
}

int sp_prop_destroy(sp_prop **p)
{
    sp_prop *pp = *p;
    prop_destroy(&pp->prp);
    free(*p);
    return SP_OK;
}

int sp_prop_init(sp_data *sp, sp_prop *p, const char *str)
{
    p->count = 0;

    prop_create(&p->prp);
    if(prop_parse(p->prp, str) == PSTATUS_NOTOK) {
        fprintf(stderr,"There was an error parsing the string.\n");
        return SP_NOT_OK;
    }
    p->bpm = 60;
    p->lbpm = 60;
    return SP_OK;
}

int sp_prop_compute(sp_data *sp, sp_prop *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->count == 0) {
        if(p->bpm != p->lbpm) {
            p->prp->scale = (SPFLOAT) 60.0 / p->bpm;
            p->lbpm = p->bpm;
        }
        p->evt = prop_next(sp, p->prp);
        p->count = prop_time(p->prp, p->evt) * sp->sr;
        switch(p->evt.type) {
            case PTYPE_ON: 
                *out = 1.0;
                break;
            case PTYPE_MAYBE: 
                if( ((SPFLOAT) sp_rand(sp) / SP_RANDMAX) > 0.5) *out = 1.0;
                else *out = 0.0;
                break;
            default:
                *out = 0.0;
                break;
        }
        return SP_OK;
    }
    *out = 0;
    p->count--;

    return SP_OK;
}

static int stack_push(prop_stack *ps, uint32_t val)
{
    if(ps->pos++ < 16) {
        ps->stack[ps->pos] = val;
    }
    return SP_OK;
}

static void stack_init(prop_stack *ps)
{
    ps->pos = -1;
    int n;
    for(n = 0; n < 16; n++) ps->stack[n] = 1;
}

static uint32_t stack_pop(prop_stack *ps)
{
    if(ps->pos >= 0) {
        return ps->stack[ps->pos--];
    }
    return 1;
}

static void mode_insert_event(prop_data *pd, char type)
{
#ifdef DEBUG_PROP
    if(type == PTYPE_ON) {
        printf("mode_insert: PTYPE_ON\n");
    } else {
        printf("mode_insert: PTYPE_OFF\n");
    }
    printf("\tval/mul = %d, pos = %d, cons = %d, div = %d\n", 
            pd->mul, pd->num, pd->cons_mul, pd->div);
#endif

    prop_val val;
    val.type = PTYPE_EVENT;
    prop_event *evt = malloc(sizeof(prop_event));
    evt->type = type;
    evt->val = pd->mul;
    evt->cons = pd->cons_mul;
    val.ud = evt;
    prop_list_append(pd->main, val);
}

static void mode_setdiv(prop_data *pd, char n)
{
    if(pd->tmp == 0 && n == 0) n = 1;
    pd->tmp *= 10;
    pd->tmp += n;
}

static void mode_setmul(prop_data *pd)
{
    pd->mul *= pd->tmp;
    pd->div = pd->tmp;
    stack_push(&pd->mstack, pd->tmp);
    pd->tmp = 0;
}

static void mode_unsetmul(prop_data *pd)
{
    uint32_t div = stack_pop(&pd->mstack);
#ifdef DEBUG_PROP
    printf("mul / div = %d / %d\n", pd->mul, div);
#endif
    pd->mul /= div;
}

static void mode_setcons(prop_data *pd)
{
    pd->cons_mul *= pd->tmp;
    pd->cons_div = pd->tmp;
    stack_push(&pd->cstack, pd->tmp);
    pd->tmp = 0;
}

static void mode_unsetcons(prop_data *pd)
{
    uint32_t div = stack_pop(&pd->cstack);
#ifdef DEBUG_PROP
    printf("mul / div = %d / %d\n", pd->cons_mul, div);
#endif
    pd->cons_mul /= div;
}

static int prop_create(prop_data **pd)
{
    *pd = malloc(sizeof(prop_data));
    prop_data *pdp = *pd;

    pdp->mul = 1;
    pdp->div = 0;
    pdp->scale = 1;
    pdp->cons_mul = 1;
    pdp->cons_div = 0;
    pdp->mode = PMODE_INIT;
    pdp->pos = 1;
    pdp->main = &pdp->top;
    pdp->main->lvl = 0;
    pdp->tmp = 0;

    stack_init(&pdp->mstack);
    stack_init(&pdp->cstack);
    prop_list_init(pdp->main);

    return PSTATUS_OK;
}

static int prop_parse(prop_data *pd, const char *str)
{
    char c;
    while(*str != 0) {
        c = str[0];

        switch(c) {
            case '+':
                mode_insert_event(pd, PTYPE_ON);
                break;
            case '?':
                mode_insert_event(pd, PTYPE_MAYBE);
                break;
            case '-':
                mode_insert_event(pd, PTYPE_OFF);
                break;

            case '0':
                mode_setdiv(pd, 0);
                break;
            case '1':
                mode_setdiv(pd, 1);
                break;
            case '2':
                mode_setdiv(pd, 2);
                break;
            case '3':
                mode_setdiv(pd, 3);
                break;
            case '4':
                mode_setdiv(pd, 4);
                break;
            case '5':
                mode_setdiv(pd, 5);
                break;
            case '6':
                mode_setdiv(pd, 6);
                break;
            case '7':
                mode_setdiv(pd, 7);
                break;
            case '8':
                mode_setdiv(pd, 8);
                break;
            case '9':
                mode_setdiv(pd, 9);
                break;
            case '(':
                mode_setmul(pd);
                break;
            case ')':
                mode_unsetmul(pd);
                break;
            case '[':
                mode_setcons(pd);
                break;
            case ']':
                mode_unsetcons(pd);
                break;
            case '|':
                mode_insert_slice(pd);
                break;
            case '{':
                mode_list_start(pd);
                break;
            case '}':
                mode_list_end(pd);
                break;
            case ' ': break;
            case '\n': break;
            case '\t': break;

            default:
                return PSTATUS_NOTOK;
        }
        pd->pos++;
        str++;
    }
    prop_list_reset(&pd->top);
    pd->main = &pd->top;
    return PSTATUS_OK;
}

prop_val prop_list_iterate(prop_list *lst)
{
    if(lst->pos >= lst->size) {
        prop_list_reset(lst);
    }
    prop_val val = lst->last->val;
    lst->last = lst->last->next;
    lst->pos++;
    return val; 
}

static void back_to_top(prop_data *pd)
{
    prop_list *lst = pd->main;
    prop_list_reset(lst);
    pd->main = lst->top;
    reset(pd);
}

static void reset(prop_data *pd)
{
    prop_list *lst = pd->main;
    if(lst->pos >= lst->size) {
        back_to_top(pd);
    }
}

prop_event prop_next(sp_data *sp, prop_data *pd)
{
/*
    prop_list *lst = pd->main;

    if(lst->pos >= lst->size) {
        //prop_list_reset(lst);
        pd->main = lst->top;
    }
*/
    reset(pd); 
    prop_list *lst = pd->main;

    prop_val val = lst->last->val;
    lst->last = lst->last->next;
    lst->pos++;

    switch(val.type) {
        case PTYPE_SLICE: {
            prop_slice *slice = (prop_slice *)val.ud;

            uint32_t pos = floor(
                ((SPFLOAT)sp_rand(sp) / SP_RANDMAX) 
                * slice->size);

            pd->main = slice->ar[pos];
            prop_list_reset(pd->main);
            return prop_next(sp, pd);
            break;
        }
        case PTYPE_LIST: {
            prop_list *lst = (prop_list *)val.ud;
            pd->main = lst;
            prop_list_reset(pd->main);
            return prop_next(sp, pd);
            break;
        }
        default:
            break;
    }
    prop_event *p = (prop_event *)val.ud;
    return *p;
}

static float prop_time(prop_data *pd, prop_event evt)
{
    float val = evt.cons * (pd->scale / evt.val);
    return val;
}

static int prop_destroy(prop_data **pd)
{
    prop_data *pdp = *pd;

    prop_list_destroy(&pdp->top);

    free(*pd);
    return PSTATUS_OK;
}

static int prop_list_init(prop_list *lst)
{
    lst->last = &lst->root;
    lst->size = 0;
    lst->pos = 0;
    lst->root.val.type = PTYPE_NULL;
    lst->top = lst;
    return PSTATUS_OK;
}

static int prop_list_append(prop_list *lst, prop_val val)
{
    prop_entry *new_ = malloc(sizeof(prop_entry));
    new_->val = val;
    lst->last->next = new_;
    lst->last = new_;
    lst->size++;
    return PSTATUS_OK;
}

static int prop_slice_free(prop_slice *slice)
{
    uint32_t i;
    for(i = 0; i < slice->size; i++) {
        prop_list_destroy(slice->ar[i]);   
        free(slice->ar[i]); 
    }
    free(slice->ar);
    return PSTATUS_OK;
}

static int prop_val_free(prop_val val)
{
    switch(val.type) {
        case PTYPE_SLICE:
            prop_slice_free((prop_slice *)val.ud);
            free(val.ud);
            break;
        case PTYPE_LIST:
            prop_list_destroy((prop_list *)val.ud);
            free(val.ud);
            break;
        default:
            free(val.ud);
            break;
    }
    return PSTATUS_OK;
}

static int prop_list_destroy(prop_list *lst) 
{
    prop_entry *entry = lst->root.next;
    prop_entry *next;
    uint32_t i;

    for(i = 0; i < lst->size; i++) {
        next = entry->next;
        prop_val_free(entry->val);
        free(entry);
        entry = next;
    }
    return PSTATUS_OK;
}

static void prop_list_reset(prop_list *lst)
{
    lst->last = lst->root.next;
    lst->pos = 0;
}   

static void mode_insert_slice(prop_data *pd)
{
    prop_entry *entry = pd->main->top->last;
    if(entry->val.type != PTYPE_SLICE) {
        prop_slice_encap(pd);
    } else {
        prop_slice_append(pd);
    }
}

static void prop_slice_encap(prop_data *pd)
{
    prop_val val;
    prop_list *top = pd->main->top;
    val.type = PTYPE_SLICE;
    prop_slice *slice = malloc(sizeof(prop_slice));
    val.ud = slice;
    prop_list *lst, *new_;
    prop_list_copy(pd->main, &lst);
    new_ = malloc(sizeof(prop_list));
    new_->lvl = pd->main->lvl;
    slice->size = 2;
    slice->ar = 
        (prop_list **)malloc(sizeof(prop_list *) * slice->size);
    slice->ar[0] = lst;
    /* reinit main list */
    prop_list_init(pd->main);
    prop_list_append(pd->main, val);
    slice->ar[1] = new_;
    prop_list_init(slice->ar[1]);
    pd->main = slice->ar[1];

    slice->ar[0]->top = top;
    slice->ar[1]->top = top;
}

static void prop_slice_append(prop_data *pd)
{
    prop_entry *entry = pd->main->top->last;
    prop_slice *slice = (prop_slice *)entry->val.ud;
    
    prop_list *new_ = malloc(sizeof(prop_list));
    prop_list_init(new_);
    slice->size++;
    slice->ar = (prop_list **)
        realloc(slice->ar, sizeof(prop_list *) * slice->size);
    slice->ar[slice->size - 1] = new_;
    new_->top = pd->main->top;
    pd->main = new_;
}

static int prop_list_copy(prop_list *src, prop_list **dst)
{
    *dst = malloc(sizeof(prop_list));
    prop_list *pdst = *dst;
    pdst->root = src->root;
    pdst->last = src->last;
    pdst->size = src->size;
    pdst->pos = src->pos;
    pdst->lvl = src->lvl;
    return PSTATUS_OK;
}

static void mode_list_start(prop_data *pd)
{
    prop_val val;
    val.type = PTYPE_LIST;
    prop_list *new_ = malloc(sizeof(prop_list));
    prop_list_init(new_);
    new_->lvl = pd->main->lvl + 1;
    val.ud = new_;
    prop_list_append(pd->main, val);
    new_->top = pd->main;
    pd->main = new_;
}

static void mode_list_end(prop_data *pd)
{
    pd->main = pd->main->top;
}

int sp_prop_reset(sp_data *sp, sp_prop *p)
{
    back_to_top(p->prp);
    p->count = 0;
    return SP_OK;
}


typedef struct {
	float fVec0[65536];
	float fRec0[2];
	int IOTA;
	FAUSTFLOAT fHslider0;
	FAUSTFLOAT fHslider1;
	FAUSTFLOAT fHslider2;
	int fSamplingFreq;
} pshift;

static pshift* newpshift() { 
	pshift* dsp = (pshift*)malloc(sizeof(pshift));
	return dsp;
}

static void deletepshift(pshift* dsp) { 
	free(dsp);
}

static void instanceInitpshift(pshift* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->IOTA = 0;
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 65536); i0 = (i0 + 1)) {
			dsp->fVec0[i0] = 0.f;
			
		}
		
	}
	dsp->fHslider0 = (FAUSTFLOAT)1000.;
	dsp->fHslider1 = (FAUSTFLOAT)0.;
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec0[i1] = 0.f;
			
		}
		
	}
	dsp->fHslider2 = (FAUSTFLOAT)10.;
}

static void initpshift(pshift* dsp, int samplingFreq) {
	instanceInitpshift(dsp, samplingFreq);
}

static void buildUserInterfacepshift(pshift* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "shift", &dsp->fHslider1, 0.f, -24.f, 24.f, 0.1f);
	interface->addHorizontalSlider(interface->uiInterface, "window", &dsp->fHslider0, 1000.f, 50.f, 10000.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "xfade", &dsp->fHslider2, 10.f, 1.f, 10000.f, 1.f);
}

static void computepshift(pshift* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* input0 = inputs[0];
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fHslider0;
	float fSlow1 = ((1.f + fSlow0) - powf(2.f, (0.0833333f * (float)dsp->fHslider1)));
	float fSlow2 = (1.f / (float)dsp->fHslider2);
	float fSlow3 = (fSlow0 - 1.f);
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			float fTemp0 = (float)input0[i];
			dsp->fVec0[(dsp->IOTA & 65535)] = fTemp0;
			dsp->fRec0[0] = fmodf((dsp->fRec0[1] + fSlow1), fSlow0);
			int iTemp1 = (int)dsp->fRec0[0];
			int iTemp2 = (1 + iTemp1);
			float fTemp3 = min((fSlow2 * dsp->fRec0[0]), 1.f);
			float fTemp4 = (dsp->fRec0[0] + fSlow0);
			int iTemp5 = (int)fTemp4;
			output0[i] = (FAUSTFLOAT)((((dsp->fVec0[((dsp->IOTA - (iTemp1 & 65535)) & 65535)] * ((float)iTemp2 - dsp->fRec0[0])) + ((dsp->fRec0[0] - (float)iTemp1) * dsp->fVec0[((dsp->IOTA - (iTemp2 & 65535)) & 65535)])) * fTemp3) + (((dsp->fVec0[((dsp->IOTA - (iTemp5 & 65535)) & 65535)] * (0.f - ((dsp->fRec0[0] + fSlow3) - (float)iTemp5))) + ((fTemp4 - (float)iTemp5) * dsp->fVec0[((dsp->IOTA - ((1 + iTemp5) & 65535)) & 65535)])) * (1.f - fTemp3)));
			dsp->IOTA = (dsp->IOTA + 1);
			dsp->fRec0[1] = dsp->fRec0[0];
			
		}
		
	}
	
}


int sp_pshift_create(sp_pshift **p)
{
    *p = malloc(sizeof(sp_pshift));
    return SP_OK;
}

int sp_pshift_destroy(sp_pshift **p)
{
    sp_pshift *pp = *p;
    pshift *dsp = pp->faust;
    deletepshift (dsp);
    free(*p);
    return SP_OK;
}

int sp_pshift_init(sp_data *sp, sp_pshift *p)
{
    pshift *dsp = newpshift(); 
    UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfacepshift(dsp, &UI);
    initpshift(dsp, sp->sr);

     
    p->shift = p->args[0]; 
    p->window = p->args[1]; 
    p->xfade = p->args[2];

    p->faust = dsp;
    return SP_OK;
}

int sp_pshift_compute(sp_data *sp, sp_pshift *p, SPFLOAT *in, SPFLOAT *out) 
{

    pshift *dsp = p->faust;
    SPFLOAT out1 = 0;
    SPFLOAT *faust_out[] = {&out1};
    SPFLOAT *faust_in[] = {in};
    computepshift(dsp, 1, faust_in, faust_out);

    *out = out1;
    return SP_OK;
}


#define MINFREQINBINS 5
#define MAXHIST 3
#define MAXWINSIZ 8192
#define MINWINSIZ 128
#define DEFAULTWINSIZ 1024
#define NPREV 20
#define MAXPEAKNOS 100
#define DEFAULTPEAKNOS 20
#define MINBW 0.03
#define BINPEROCT 48
#define BPEROOVERLOG2 69.24936196
#define FACTORTOBINS 4/0.0145453
#define BINGUARD 10
#define PARTIALDEVIANCE 0.023
#define DBSCAL 3.333
#define DBOFFSET -92.3
#define MINBIN 3
#define MINAMPS 40
#define MAXAMPS 50


#define THRSH 10.

#define COEF1 ((SPFLOAT)(.5 * 1.227054))
#define COEF2 ((SPFLOAT)(.5 * -0.302385))
#define COEF3 ((SPFLOAT)(.5 * 0.095326))
#define COEF4 ((SPFLOAT)(.5 * -0.022748))
#define COEF5 ((SPFLOAT)(.5 * 0.002533))
#define FLTLEN 5

#define NPARTIALONSET ((int)(sizeof(partialonset)/sizeof(SPFLOAT)))

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

static const SPFLOAT partialonset[] =
{
    0.0,
    48.0,
    76.0782000346154967102,
    96.0,
    111.45254855459339269887,
    124.07820003461549671089,
    134.75303625876499715823,
    144.0,
    152.15640006923099342109,
    159.45254855459339269887,
    166.05271769459026829915,
    172.07820003461549671088,
    177.62110647077242370064,
    182.75303625876499715892,
    187.53074858920888940907,
    192.0,
};

/*TODO: rename these structs */

typedef struct histopeak
{
  SPFLOAT hpitch;
  SPFLOAT hvalue;
  SPFLOAT hloud;
  int hindex;
  int hused;
} HISTOPEAK;

typedef struct peak
{
  SPFLOAT pfreq;
  SPFLOAT pwidth;
  SPFLOAT ppow;
  SPFLOAT ploudness;
} PEAK;

int sp_ptrack_create(sp_ptrack **p)
{
    *p = malloc(sizeof(sp_ptrack));
    return SP_OK;
}

int sp_ptrack_destroy(sp_ptrack **p)
{
    sp_ptrack *pp = *p;
    sp_auxdata_free(&pp->signal);
    sp_auxdata_free(&pp->prev);
    sp_auxdata_free(&pp->sin);
    sp_auxdata_free(&pp->spec2);
    sp_auxdata_free(&pp->spec1);
    sp_auxdata_free(&pp->peakarray);
    sp_fft_destroy(&pp->fft);
    free(*p);
    return SP_OK;
}

int sp_ptrack_init(sp_data *sp, sp_ptrack *p, int ihopsize, int ipeaks)
{
    p->size = ihopsize;

    int i, winsize = p->size*2, powtwo, tmp;
    SPFLOAT *tmpb;


    /* TODO: fix this warning */
    if (winsize < MINWINSIZ || winsize > MAXWINSIZ) {
      fprintf(stderr, "Woops\n");
      return SP_NOT_OK;
    }

    tmp = winsize;

    powtwo = -1;

    while (tmp) {
      tmp >>= 1;
      powtwo++;
    }

    /* 3 days of debugging later... I found this off by one error */
    /* powtwo needs to be powtwo - 1 for fft_init */
    sp_fft_init(&p->fft, powtwo - 1) ;

    /* TODO: make this error better */
    if (winsize != (1 << powtwo)) {
        fprintf(stderr, "Woops\n");
        return SP_NOT_OK;
    }

    p->hopsize = p->size;

    sp_auxdata_alloc(&p->signal, p->hopsize * sizeof(SPFLOAT));
    sp_auxdata_alloc(&p->prev, (p->hopsize*2 + 4*FLTLEN)*sizeof(SPFLOAT));
    sp_auxdata_alloc(&p->sin, (p->hopsize*2)*sizeof(SPFLOAT));
    sp_auxdata_alloc(&p->spec2, (winsize*4 + 4*FLTLEN)*sizeof(SPFLOAT));
    sp_auxdata_alloc(&p->spec1, (winsize*4)*sizeof(SPFLOAT));

    for (i = 0, tmpb = (SPFLOAT *)p->signal.ptr; i < p->hopsize; i++)
        tmpb[i] = 0.0;
    for (i = 0, tmpb = (SPFLOAT *)p->prev.ptr; i < winsize + 4 * FLTLEN; i++)
        tmpb[i] = 0.0;
    for (i = 0, tmpb = (SPFLOAT *)p->sin.ptr; i < p->hopsize; i++) {
        tmpb[2*i] =   (SPFLOAT) cos((M_PI*i)/(winsize));
        tmpb[2*i+1] = -(SPFLOAT)sin((M_PI*i)/(winsize));
    }

    p->cnt = 0;
    p->numpks = ipeaks;

    sp_auxdata_alloc(&p->peakarray, (p->numpks+1)*sizeof(PEAK));

    p->cnt = 0;
    p->histcnt = 0;
    p->sr = sp->sr;
    for (i = 0; i < NPREV; i++) p->dbs[i] = -144.0;
    p->amplo = MINAMPS;
    p->amphi = MAXAMPS;
    p->npartial = 7;
    p->dbfs = 32768.0;
    p->prevf = p->cps = 100.0;

    return SP_OK;
}

static void ptrack(sp_data *sp, sp_ptrack *p)
{
    SPFLOAT *spec = (SPFLOAT *)p->spec1.ptr;
    SPFLOAT *spectmp = (SPFLOAT *)p->spec2.ptr;
    SPFLOAT *sig = (SPFLOAT *)p->signal.ptr;
    SPFLOAT *sinus  = (SPFLOAT *)p->sin.ptr;
    SPFLOAT *prev  = (SPFLOAT *)p->prev.ptr;
    PEAK  *peaklist = (PEAK *)p->peakarray.ptr;
    HISTOPEAK histpeak;
    int i, j, k, hop = p->hopsize, n = 2*hop, npeak = 0, logn = -1, count, tmp;
    SPFLOAT totalpower = 0, totalloudness = 0, totaldb = 0;
    SPFLOAT maxbin,  *histogram = spectmp + BINGUARD;
    SPFLOAT hzperbin = (SPFLOAT) p->sr / (n + n);
    int numpks = p->numpks;
    int indx, halfhop = hop>>1;
    SPFLOAT best;
    SPFLOAT cumpow = 0, cumstrength = 0, freqnum = 0, freqden = 0;
    int npartials = 0,  nbelow8 = 0;
    SPFLOAT putfreq;

    count = p->histcnt + 1;
    if (count == NPREV) count = 0;
    p->histcnt = count;

    tmp = n;
    while (tmp) {
        tmp >>= 1;
        logn++;
    }
    maxbin = BINPEROCT * (logn-2);
    for (i = 0, k = 0; i < hop; i++, k += 2) {
        spec[k]   = sig[i] * sinus[k];
        spec[k+1] = sig[i] * sinus[k+1];
    }

    sp_fft_cpx(&p->fft, spec, hop);

    for (i = 0, k = 2*FLTLEN; i < hop; i+=2, k += 4) {
        spectmp[k]   = spec[i];
        spectmp[k+1] = spec[i+1];
    }

    for (i = n - 2, k = 2*FLTLEN+2; i >= 0; i-=2, k += 4) {
        spectmp[k]   = spec[i];
        spectmp[k+1] = -spec[i+1];
    }

    for (i = (2*FLTLEN), k = (2*FLTLEN-2);i<FLTLEN*4; i+=2, k-=2) {
        spectmp[k]   = spectmp[i];
        spectmp[k+1] = -spectmp[i+1];
    }

    for (i = (2*FLTLEN+n-2), k =(2*FLTLEN+n); i>=0; i-=2, k+=2) {
        spectmp[k]   = spectmp[i];
        spectmp[k+1] = -spectmp[k+1];
    }

    for (i = j = 0, k = 2*FLTLEN; i < halfhop; i++, j+=8, k+=2) {
        SPFLOAT re,  im;

        re= COEF1 * ( prev[k-2] - prev[k+1]  + spectmp[k-2] - prev[k+1]) +
            COEF2 * ( prev[k-3] - prev[k+2]  + spectmp[k-3]  - spectmp[ 2]) +
            COEF3 * (-prev[k-6] +prev[k+5]  -spectmp[k-6] +spectmp[k+5]) +
            COEF4 * (-prev[k-7] +prev[k+6]  -spectmp[k-7] +spectmp[k+6]) +
            COEF5 * ( prev[k-10] -prev[k+9]  +spectmp[k-10] -spectmp[k+9]);

        im= COEF1 * ( prev[k-1] +prev[k]  +spectmp[k-1] +spectmp[k]) +
            COEF2 * (-prev[k-4] -prev[k+3]  -spectmp[k-4] -spectmp[k+3]) +
            COEF3 * (-prev[k-5] -prev[k+4]  -spectmp[k-5] -spectmp[k+4]) +
            COEF4 * ( prev[k-8] +prev[k+7]  +spectmp[k-8] +spectmp[k+7]) +
            COEF5 * ( prev[k-9] +prev[k+8]  +spectmp[k-9] +spectmp[k+8]);

        spec[j]   = 0.707106781186547524400844362104849 * (re + im);
        spec[j+1] = 0.707106781186547524400844362104849 * (im - re);
        spec[j+4] = prev[k] + spectmp[k+1];
        spec[j+5] = prev[k+1] - spectmp[k];

        j += 8;
        k += 2;

        re= COEF1 * ( prev[k-2] -prev[k+1]  -spectmp[k-2] +spectmp[k+1]) +
            COEF2 * ( prev[k-3] -prev[k+2]  -spectmp[k-3] +spectmp[k+2]) +
            COEF3 * (-prev[k-6] +prev[k+5]  +spectmp[k-6] -spectmp[k+5]) +
            COEF4 * (-prev[k-7] +prev[k+6]  +spectmp[k-7] -spectmp[k+6]) +
            COEF5 * ( prev[k-10] -prev[k+9]  -spectmp[k-10] +spectmp[k+9]);

        im= COEF1 * ( prev[k-1] +prev[k]  -spectmp[k-1] -spectmp[k]) +
            COEF2 * (-prev[k-4] -prev[k+3]  +spectmp[k-4] +spectmp[k+3]) +
            COEF3 * (-prev[k-5] -prev[k+4]  +spectmp[k-5] +spectmp[k+4]) +
            COEF4 * ( prev[k-8] +prev[k+7]  -spectmp[k-8] -spectmp[k+7]) +
            COEF5 * ( prev[k-9] +prev[k+8]  -spectmp[k-9] -spectmp[k+8]);

        spec[j]   = 0.707106781186547524400844362104849 * (re + im);
        spec[j+1] = 0.707106781186547524400844362104849 * (im - re);
        spec[j+4] = prev[k] - spectmp[k+1];
        spec[j+5] = prev[k+1] + spectmp[k];

    }


    for (i = 0; i < n + 4*FLTLEN; i++) prev[i] = spectmp[i];

    for (i = 0; i < MINBIN; i++) spec[4*i + 2] = spec[4*i + 3] =0.0;

    for (i = 4*MINBIN, totalpower = 0; i < (n-2)*4; i += 4) {
        SPFLOAT re = spec[i] - 0.5 * (spec[i-8] + spec[i+8]);
        SPFLOAT im = spec[i+1] - 0.5 * (spec[i-7] + spec[i+9]);
        spec[i+3] = (totalpower += (spec[i+2] = re * re + im * im));
    }

    if (totalpower > 1.0e-9) {
        totaldb = (SPFLOAT)DBSCAL * logf(totalpower/n);
        totalloudness = (SPFLOAT)sqrtf((SPFLOAT)sqrtf(totalpower));
        if (totaldb < 0) totaldb = 0;
    }
    else totaldb = totalloudness = 0.0;

    p->dbs[count] = totaldb + DBOFFSET;

    if (totaldb >= p->amplo) {
        npeak = 0;

        for (i = 4*MINBIN;i < (4*(n-2)) && npeak < numpks; i+=4) {
            SPFLOAT height = spec[i+2], h1 = spec[i-2], h2 = spec[i+6];
            SPFLOAT totalfreq, peakfr, tmpfr1, tmpfr2, m, var, stdev;

            if (height < h1 || height < h2 ||
            h1 < 0.00001*totalpower ||
            h2 < 0.00001*totalpower) continue;

            peakfr= ((spec[i-8] - spec[i+8]) * (2.0 * spec[i] -
                                        spec[i+8] - spec[i-8]) +
             (spec[i-7] - spec[i+9]) * (2.0 * spec[i+1] -
                                        spec[i+9] - spec[i-7])) /
            (height + height);
            tmpfr1=  ((spec[i-12] - spec[i+4]) *
              (2.0 * spec[i-4] - spec[i+4] - spec[i-12]) +
              (spec[i-11] - spec[i+5]) * (2.0 * spec[i-3] -
                                          spec[i+5] - spec[i-11])) /
            (2.0 * h1) - 1;
            tmpfr2= ((spec[i-4] - spec[i+12]) * (2.0 * spec[i+4] -
                                         spec[i+12] - spec[i-4]) +
             (spec[i-3] - spec[i+13]) * (2.0 * spec[i+5] -
                                         spec[i+13] - spec[i-3])) /
            (2.0 * h2) + 1;


            m = 0.333333333333 * (peakfr + tmpfr1 + tmpfr2);
            var = 0.5 * ((peakfr-m)*(peakfr-m) +
                     (tmpfr1-m)*(tmpfr1-m) + (tmpfr2-m)*(tmpfr2-m));

            totalfreq = (i>>2) + m;
            if (var * totalpower > THRSH * height
            || var < 1.0e-30) continue;

            stdev = (SPFLOAT)sqrt((SPFLOAT)var);
            if (totalfreq < 4) totalfreq = 4;


            peaklist[npeak].pwidth = stdev;
            peaklist[npeak].ppow = height;
            peaklist[npeak].ploudness = sqrt(sqrt(height));
            peaklist[npeak].pfreq = totalfreq;
            npeak++;
        }

          if (npeak > numpks) npeak = numpks;
          for (i = 0; i < maxbin; i++) histogram[i] = 0;
          for (i = 0; i < npeak; i++) {
            SPFLOAT pit = (SPFLOAT)(BPEROOVERLOG2 * logf(peaklist[i].pfreq) - 96.0);
            SPFLOAT binbandwidth = FACTORTOBINS * peaklist[i].pwidth/peaklist[i].pfreq;
            SPFLOAT putbandwidth = (binbandwidth < 2.0 ? 2.0 : binbandwidth);
            SPFLOAT weightbandwidth = (binbandwidth < 1.0 ? 1.0 : binbandwidth);
            SPFLOAT weightamp = 4.0 * peaklist[i].ploudness / totalloudness;
            for (j = 0; j < NPARTIALONSET; j++) {
              SPFLOAT bin = pit - partialonset[j];
              if (bin < maxbin) {
                SPFLOAT para, pphase, score = 30.0 * weightamp /
                  ((j+p->npartial) * weightbandwidth);
                int firstbin = bin + 0.5 - 0.5 * putbandwidth;
                int lastbin = bin + 0.5 + 0.5 * putbandwidth;
                int ibw = lastbin - firstbin;
                if (firstbin < -BINGUARD) break;
                para = 1.0 / (putbandwidth * putbandwidth);
                for (k = 0, pphase = firstbin-bin; k <= ibw;
                     k++,pphase += 1.0)
                  histogram[k+firstbin] += score * (1.0 - para * pphase * pphase);

              }
            }
          }


        for (best = 0, indx = -1, j=0; j < maxbin; j++) {
            if (histogram[j] > best) {
                indx = j;  
                best = histogram[j];
            }
        }

        histpeak.hvalue = best;
        histpeak.hindex = indx;

        putfreq = expf((1.0 / BPEROOVERLOG2) * (histpeak.hindex + 96.0));

        for (j = 0; j < npeak; j++) {
            SPFLOAT fpnum = peaklist[j].pfreq/putfreq;
            int pnum = (int)(fpnum + 0.5);
            SPFLOAT fipnum = pnum;
            SPFLOAT deviation;
            if (pnum > 16 || pnum < 1) continue;
            deviation = 1.0 - fpnum/fipnum;
            if (deviation > -PARTIALDEVIANCE && deviation < PARTIALDEVIANCE) {
                SPFLOAT stdev, weight;
                npartials++;
                if (pnum < 8) nbelow8++;
                cumpow += peaklist[j].ppow;
                cumstrength += sqrt(sqrt(peaklist[j].ppow));
                stdev = (peaklist[j].pwidth > MINBW ?
                   peaklist[j].pwidth : MINBW);
                weight = 1.0 / ((stdev*fipnum) * (stdev*fipnum));
                freqden += weight;
                freqnum += weight * peaklist[j].pfreq/fipnum;
            }
        }
        if ((nbelow8 < 4 || npartials < 7) && cumpow < 0.01 * totalpower) {
            histpeak.hvalue = 0;
        } else {
            SPFLOAT pitchpow = (cumstrength * cumstrength);
            SPFLOAT freqinbins = freqnum/freqden;
            pitchpow = pitchpow * pitchpow;

            if (freqinbins < MINFREQINBINS) {
                histpeak.hvalue = 0;
            } else {
                p->cps = histpeak.hpitch = hzperbin * freqnum/freqden;
                histpeak.hloud = DBSCAL * logf(pitchpow/n);
            }
        }
    }
}

int sp_ptrack_compute(sp_data *sp, sp_ptrack *p, SPFLOAT *in, SPFLOAT *freq, SPFLOAT *amp)
{
    SPFLOAT *buf = (SPFLOAT *)p->signal.ptr;
    int pos = p->cnt, h = p->hopsize;
    SPFLOAT scale = p->dbfs;

    if (pos == h) {
        ptrack(sp,p);
        pos = 0;
    }
    buf[pos] = *in * scale;
    pos++;

    *freq = p->cps;
    *amp =  exp(p->dbs[p->histcnt] / 20.0 * log(10.0));
    
    p->cnt = pos;

    return SP_OK;
}


int sp_randh_create(sp_randh **p)
{
    *p = malloc(sizeof(sp_randh));
    return SP_OK;
}

int sp_randh_destroy(sp_randh **p)
{
    free(*p);
    return SP_OK;
}

int sp_randh_init(sp_data *sp, sp_randh *p)
{
    p->counter = 0;
    p->freq = 10;
    p->dur = (sp->sr / p->freq);
    p->min = 0;
    p->max = 1;
    p->val = 0;
    return SP_OK;
}

int sp_randh_compute(sp_data *sp, sp_randh *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->counter == 0) {
        p->val = p->min + ((SPFLOAT) sp_rand(sp) / SP_RANDMAX) * (p->max - p->min);
        
        if(p->freq == 0) {
            p->dur = 1;
        } else {
            p->dur = (sp->sr / p->freq) + 1; 
        }

        *out = p->val;
    } else {
        *out = p->val;
    }
    p->counter = (p->counter + 1) % p->dur;
    return SP_OK;
}


#define sp_oneUp31Bit      (4.656612875245796924105750827168e-10)

#define sp_randGab   ((SPFLOAT)     \
        (((p->holdrand = p->holdrand * 214013 + 2531011) >> 1)  \
         & 0x7fffffff) * sp_oneUp31Bit)


int sp_randi_create(sp_randi **p)
{
    *p = malloc(sizeof(sp_randi));
    return SP_OK;
}

int sp_randi_destroy(sp_randi **p)
{
    free(*p);
    return SP_OK;
}

int sp_randi_init(sp_data *sp, sp_randi *p)
{
    p->sicvt = 1.0 * SP_FT_MAXLEN / sp->sr;
    p->phs = 0;
    p->min = 0;
    p->max = 1;
    p->cps = 3;
    p->mode = 3;
    p->holdrand = sp_rand(sp);
    p->fstval = 0;

    int mode = (int)(p->mode);
    switch (mode) {
    case 1: /* immediate interpolation between kmin and 1st random number */
        p->num1 = 0.0;
        p->num2 = sp_randGab;
        p->dfdmax = (p->num2 - p->num1) / SP_FT_MAXLEN * 1.0;
        break;
    case 2: /* immediate interpolation between ifirstval and 1st random number */
        p->num1 = (p->max - p->min) ?
          (p->fstval - p->min) / (p->max - p->min) : 0.0;
        p->num2 = sp_randGab;
        p->dfdmax = (p->num2 - p->num1) / SP_FT_MAXLEN * 1.0;
        break;
    case 3: /* immediate interpolation between 1st and 2nd random number */
        p->num1 = sp_randGab;
        p->num2 = sp_randGab;
        p->dfdmax = (p->num2 - p->num1) / SP_FT_MAXLEN * 1.0;
        break;
    default: /* old behaviour as developped by Gabriel */
        p->num1 = p->num2 = 0.0;
        p->dfdmax = 0.0;
    }
    return SP_OK;
}

int sp_randi_compute(sp_data *sp, sp_randi *p, SPFLOAT *in, SPFLOAT *out)
{
    int32_t phs = p->phs, inc;
    SPFLOAT cpsp;
    SPFLOAT amp, min;

    cpsp = p->cps;
    min = p->min;
    amp =  (p->max - min);
    inc = (int32_t)(cpsp * p->sicvt);

    *out = (p->num1 + (SPFLOAT)phs * p->dfdmax) * amp + min;
    phs += inc;
    if (phs >= SP_FT_MAXLEN) {
        phs &= SP_FT_PHMASK;
        p->num1 = p->num2;
        p->num2 = sp_randGab;
        p->dfdmax = 1.0 * (p->num2 - p->num1) / SP_FT_MAXLEN;
    }
    p->phs = phs;

    return SP_OK;
}


#define N           (624)
#define M           (397)
#define MATRIX_A    0x9908B0DFU     /* constant vector a */
#define UPPER_MASK  0x80000000U     /* most significant w-r bits */
#define LOWER_MASK  0x7FFFFFFFU     /* least significant r bits */

static void MT_update_state(uint32_t *mt)
{
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    const uint32_t  mag01[2] = { (uint32_t) 0, (uint32_t) MATRIX_A };
    int       i;
    uint32_t  y;

    for (i = 0; i < (N - M); i++) {
      y = (mt[i] & UPPER_MASK) | (mt[i + 1] & LOWER_MASK);
      mt[i] = mt[i + M] ^ (y >> 1) ^ mag01[y & (uint32_t) 1];
    }
    for ( ; i < (N - 1); i++) {
      y = (mt[i] & UPPER_MASK) | (mt[i + 1] & LOWER_MASK);
      mt[i] = mt[i + (M - N)] ^ (y >> 1) ^ mag01[y & (uint32_t) 1];
    }
    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & (uint32_t) 1];
}

/* generates a random number on [0,0xffffffff]-interval */

uint32_t sp_randmt_compute(sp_randmt *p)
{
    int       i = p->mti;
    uint32_t  y;

    if (i >= N) {                   /* generate N words at one time */
      MT_update_state(&(p->mt[0]));
      i = 0;
    }
    y = p->mt[i];
    p->mti = i + 1;
    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & (uint32_t) 0x9D2C5680U;
    y ^= (y << 15) & (uint32_t) 0xEFC60000U;
    y ^= (y >> 18);

    return y;
}

void sp_randmt_seed(sp_randmt *p,
    const uint32_t *initKey, uint32_t keyLength)
{
    int       i, j, k;
    uint32_t  x;

    /* if array is NULL, use length parameter as simple 32 bit seed */
    x = (initKey == NULL ? keyLength : (uint32_t) 19650218);
    p->mt[0] = x;
    for (i = 1; i < N; i++) {
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      x = ((uint32_t) 1812433253 * (x ^ (x >> 30)) + (uint32_t) i);
      p->mt[i] = x;
    }
    p->mti = N;
    if (initKey == NULL)
      return;
    i = 0; j = 0;
    k = (N > (int) keyLength ? N : (int) keyLength);
    for ( ; k; k--) {
      x = p->mt[i++];
      p->mt[i] = (p->mt[i] ^ ((x ^ (x >> 30)) * (uint32_t) 1664525))
                 + initKey[j] + (uint32_t) j;   /* non linear */
      if (i == (N - 1)) {
        p->mt[0] = p->mt[N - 1];
        i = 0;
      }
      if (++j >= (int) keyLength)
        j = 0;
    }
    for (k = (N - 1); k; k--) {
      x = p->mt[i++];
      p->mt[i] = (p->mt[i] ^ ((x ^ (x >> 30)) * (uint32_t) 1566083941))
                 - (uint32_t) i;                /* non linear */
      if (i == (N - 1)) {
        p->mt[0] = p->mt[N - 1];
        i = 0;
      }
    }
    /* MSB is 1; assuring non-zero initial array */
    p->mt[0] = (uint32_t) 0x80000000U;
}

#undef N
#undef M

int sp_reson_create(sp_reson **p)
{
    *p = malloc(sizeof(sp_reson));
    return SP_OK;
}

int sp_reson_destroy(sp_reson **p)
{
    free(*p);
    return SP_OK;
}

int sp_reson_init(sp_data *sp, sp_reson *p)
{
    p->scale = 0;
    p->freq = 4000;
    p->bw = 1000;
    p->prvfreq = p->prvbw = -100.0;
    p->tpidsr = (2.0 * M_PI) / sp->sr * 1.0;
    p->yt1 = p->yt2 = 0.0;
    return SP_OK;
}


int sp_reson_compute(sp_data *sp, sp_reson *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT c3p1, c3t4;
    SPFLOAT yt1, yt2, c1 = p->c1, c2 = p->c2, c3 = p->c3;
    int flag = 0;

    yt1 = p->yt1; 
    yt2 = p->yt2;
    
    SPFLOAT yt0;
    SPFLOAT cf = p->freq;
    
    /* bw needs to stay positive so it doesn't blow the filter up */
    SPFLOAT bw = fabs(p->bw);
    
    if (cf != p->prvfreq ) {
        p->prvfreq = cf;
        p->cosf = cos(cf * (p->tpidsr));
        flag = 1;
    }
    
    if (bw != p->prvbw) {
        p->prvbw = bw;
        c3 = p->c3 = exp(bw * (-1.0 * p->tpidsr));
        flag = 1;
    }
    
    if (flag) {
        c3p1 = c3 + 1.0;
        c3t4 = c3 * 4.0;
        c2 = p->c2 = c3t4 * p->cosf / c3p1;
        c1 = p->c1 = 1.0;
        flag = 0;
    }
    
    yt0 = c1 * *in  + c2 * yt1 - c3 * yt2;
    *out = yt0;
    yt2 = yt1;
    yt1 = yt0;
    p->yt1 = yt1; p->yt2 = yt2;
    return SP_OK;
}


int sp_reverse_create(sp_reverse **p)
{
    *p = malloc(sizeof(sp_reverse));
    return SP_OK;
}

int sp_reverse_destroy(sp_reverse **p)
{
    sp_reverse *pp = *p;
    sp_auxdata_free(&pp->buf);
    free(*p);
    return SP_OK;
}

int sp_reverse_init(sp_data *sp, sp_reverse *p, SPFLOAT delay)
{
    size_t size = delay * sp->sr * sizeof(SPFLOAT) * 2;
    p->bufpos = 0;
    sp_auxdata_alloc(&p->buf, size);
    p->bufsize = (uint32_t)p->buf.size / sizeof(SPFLOAT);
    return SP_OK;
}

int sp_reverse_compute(sp_data *sp, sp_reverse *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT *buf = (SPFLOAT *)p->buf.ptr;
    *out = buf[p->bufpos];
    buf[(p->bufsize - 1) - p->bufpos] = *in;
    p->bufpos = (p->bufpos + 1) % p->bufsize;
    return SP_OK;
}

#define DEFAULT_SRATE   44100.0
#define MIN_SRATE       5000.0
#define MAX_SRATE       1000000.0
#define MAX_PITCHMOD    20.0
#define DELAYPOS_SHIFT  28
#define DELAYPOS_SCALE  0x10000000
#define DELAYPOS_MASK   0x0FFFFFFF

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

/* reverbParams[n][0] = delay time (in seconds)                     */
/* reverbParams[n][1] = random variation in delay time (in seconds) */
/* reverbParams[n][2] = random variation frequency (in 1/sec)       */
/* reverbParams[n][3] = random seed (0 - 32767)                     */

static const SPFLOAT reverbParams[8][4] = {
    { (2473.0 / DEFAULT_SRATE), 0.0010, 3.100,  1966.0 },
    { (2767.0 / DEFAULT_SRATE), 0.0011, 3.500, 29491.0 },
    { (3217.0 / DEFAULT_SRATE), 0.0017, 1.110, 22937.0 },
    { (3557.0 / DEFAULT_SRATE), 0.0006, 3.973,  9830.0 },
    { (3907.0 / DEFAULT_SRATE), 0.0010, 2.341, 20643.0 },
    { (4127.0 / DEFAULT_SRATE), 0.0011, 1.897, 22937.0 },
    { (2143.0 / DEFAULT_SRATE), 0.0017, 0.891, 29491.0 },
    { (1933.0 / DEFAULT_SRATE), 0.0006, 3.221, 14417.0 }
};

static int delay_line_max_samples(SPFLOAT sr, SPFLOAT iPitchMod, int n);
static int init_delay_line(sp_revsc *p, sp_revsc_dl *lp, int n);
static int delay_line_bytes_alloc(SPFLOAT sr, SPFLOAT iPitchMod, int n);
static const SPFLOAT outputGain  = 0.35;
static const SPFLOAT jpScale     = 0.25;
int sp_revsc_create(sp_revsc **p){
    *p = malloc(sizeof(sp_revsc));
    return SP_OK;
}

int sp_revsc_init(sp_data *sp, sp_revsc *p)
{
    p->iSampleRate = sp->sr;
    p->sampleRate = sp->sr;
    p->feedback = 0.97;
    p->lpfreq = 10000;
    p->iPitchMod = 1;
    p->iSkipInit = 0;
    p->dampFact = 1.0;
    p->prv_LPFreq = 0.0;
    p->initDone = 1;
    int i, nBytes = 0;
    for(i = 0; i < 8; i++){
        nBytes += delay_line_bytes_alloc(sp->sr, 1, i);
    }
    sp_auxdata_alloc(&p->aux, nBytes);
    nBytes = 0;
    for (i = 0; i < 8; i++) {
        p->delayLines[i].buf = (p->aux.ptr) + nBytes;
        init_delay_line(p, &p->delayLines[i], i);
        nBytes += delay_line_bytes_alloc(sp->sr, 1, i);
    }

    return SP_OK;
}


int sp_revsc_destroy(sp_revsc **p)
{
    sp_revsc *pp = *p;
    sp_auxdata_free(&pp->aux);
    free(*p);
    return SP_OK;
}

static int delay_line_max_samples(SPFLOAT sr, SPFLOAT iPitchMod, int n)
{
    SPFLOAT maxDel;

    maxDel = reverbParams[n][0];
    maxDel += (reverbParams[n][1] * (SPFLOAT) iPitchMod * 1.125);
    return (int) (maxDel * sr + 16.5);
}

static int delay_line_bytes_alloc(SPFLOAT sr, SPFLOAT iPitchMod, int n)
{
    int nBytes = 0;

    nBytes += (delay_line_max_samples(sr, iPitchMod, n) * (int) sizeof(SPFLOAT));
    return nBytes;
}

static void next_random_lineseg(sp_revsc *p, sp_revsc_dl *lp, int n)
{
    SPFLOAT prvDel, nxtDel, phs_incVal;

    /* update random seed */
    if (lp->seedVal < 0)
      lp->seedVal += 0x10000;
    lp->seedVal = (lp->seedVal * 15625 + 1) & 0xFFFF;
    if (lp->seedVal >= 0x8000)
      lp->seedVal -= 0x10000;
    /* length of next segment in samples */
    lp->randLine_cnt = (int) ((p->sampleRate / reverbParams[n][2]) + 0.5);
    prvDel = (SPFLOAT) lp->writePos;
    prvDel -= ((SPFLOAT) lp->readPos
               + ((SPFLOAT) lp->readPosFrac / (SPFLOAT) DELAYPOS_SCALE));
    while (prvDel < 0.0)
      prvDel += lp->bufferSize;
    prvDel = prvDel / p->sampleRate;    /* previous delay time in seconds */
    nxtDel = (SPFLOAT) lp->seedVal * reverbParams[n][1] / 32768.0;
    /* next delay time in seconds */
    nxtDel = reverbParams[n][0] + (nxtDel * (SPFLOAT) p->iPitchMod);
    /* calculate phase increment per sample */
    phs_incVal = (prvDel - nxtDel) / (SPFLOAT) lp->randLine_cnt;
    phs_incVal = phs_incVal * p->sampleRate + 1.0;
    lp->readPosFrac_inc = (int) (phs_incVal * DELAYPOS_SCALE + 0.5);
}

static int init_delay_line(sp_revsc *p, sp_revsc_dl *lp, int n)
{
    SPFLOAT readPos;
    /* int     i; */

    /* calculate length of delay line */
    lp->bufferSize = delay_line_max_samples(p->sampleRate, 1, n);
    lp->dummy = 0;
    lp->writePos = 0;
    /* set random seed */
    lp->seedVal = (int) (reverbParams[n][3] + 0.5);
    /* set initial delay time */
    readPos = (SPFLOAT) lp->seedVal * reverbParams[n][1] / 32768;
    readPos = reverbParams[n][0] + (readPos * (SPFLOAT) p->iPitchMod);
    readPos = (SPFLOAT) lp->bufferSize - (readPos * p->sampleRate);
    lp->readPos = (int) readPos;
    readPos = (readPos - (SPFLOAT) lp->readPos) * (SPFLOAT) DELAYPOS_SCALE;
    lp->readPosFrac = (int) (readPos + 0.5);
    /* initialise first random line segment */
    next_random_lineseg(p, lp, n);
    /* clear delay line to zero */
    lp->filterState = 0.0;
    memset(lp->buf, 0, sizeof(SPFLOAT) * lp->bufferSize);
    return SP_OK;
}


int sp_revsc_compute(sp_data *sp, sp_revsc *p, SPFLOAT *in1, SPFLOAT *in2, SPFLOAT *out1, SPFLOAT *out2)
{
    SPFLOAT ainL, ainR, aoutL, aoutR;
    SPFLOAT vm1, v0, v1, v2, am1, a0, a1, a2, frac;
    sp_revsc_dl *lp;
    int readPos;
    uint32_t n;
    int bufferSize; /* Local copy */
    SPFLOAT dampFact = p->dampFact;

    if (p->initDone <= 0) return SP_NOT_OK;

    /* calculate tone filter coefficient if frequency changed */

    if (p->lpfreq != p->prv_LPFreq) {
        p->prv_LPFreq = p->lpfreq;
        dampFact = 2.0 - cos(p->prv_LPFreq * (2 * M_PI) / p->sampleRate);
        dampFact = p->dampFact = dampFact - sqrt(dampFact * dampFact - 1.0);
    }

    /* calculate "resultant junction pressure" and mix to input signals */

    ainL = aoutL = aoutR = 0.0;
    for (n = 0; n < 8; n++) {
        ainL += p->delayLines[n].filterState;
    }
    ainL *= jpScale;
    ainR = ainL + *in2;
    ainL = ainL + *in1;

    /* loop through all delay lines */

    for (n = 0; n < 8; n++) {
        lp = &p->delayLines[n];
        bufferSize = lp->bufferSize;

        /* send input signal and feedback to delay line */

        lp->buf[lp->writePos] = (SPFLOAT) ((n & 1 ? ainR : ainL)
                                 - lp->filterState);
        if (++lp->writePos >= bufferSize) {
            lp->writePos -= bufferSize;
        }

        /* read from delay line with cubic interpolation */

        if (lp->readPosFrac >= DELAYPOS_SCALE) {
            lp->readPos += (lp->readPosFrac >> DELAYPOS_SHIFT);
            lp->readPosFrac &= DELAYPOS_MASK;
        }
        if (lp->readPos >= bufferSize)
        lp->readPos -= bufferSize;
        readPos = lp->readPos;
        frac = (SPFLOAT) lp->readPosFrac * (1.0 / (SPFLOAT) DELAYPOS_SCALE);

        /* calculate interpolation coefficients */

        a2 = frac * frac; a2 -= 1.0; a2 *= (1.0 / 6.0);
        a1 = frac; a1 += 1.0; a1 *= 0.5; am1 = a1 - 1.0;
        a0 = 3.0 * a2; a1 -= a0; am1 -= a2; a0 -= frac;

        /* read four samples for interpolation */

        if (readPos > 0 && readPos < (bufferSize - 2)) {
            vm1 = (SPFLOAT) (lp->buf[readPos - 1]);
            v0  = (SPFLOAT) (lp->buf[readPos]);
            v1  = (SPFLOAT) (lp->buf[readPos + 1]);
            v2  = (SPFLOAT) (lp->buf[readPos + 2]);
        }
        else {

        /* at buffer wrap-around, need to check index */

        if (--readPos < 0) readPos += bufferSize;
            vm1 = (SPFLOAT) lp->buf[readPos];
        if (++readPos >= bufferSize) readPos -= bufferSize;
            v0 = (SPFLOAT) lp->buf[readPos];
        if (++readPos >= bufferSize) readPos -= bufferSize;
            v1 = (SPFLOAT) lp->buf[readPos];
        if (++readPos >= bufferSize) readPos -= bufferSize;
            v2 = (SPFLOAT) lp->buf[readPos];
        }
        v0 = (am1 * vm1 + a0 * v0 + a1 * v1 + a2 * v2) * frac + v0;

        /* update buffer read position */

        lp->readPosFrac += lp->readPosFrac_inc;

        /* apply feedback gain and lowpass filter */

        v0 *= (SPFLOAT) p->feedback;
        v0 = (lp->filterState - v0) * dampFact + v0;
        lp->filterState = v0;

        /* mix to output */

        if (n & 1) {
            aoutR += v0;
        }else{
            aoutL += v0;
        }

        /* start next random line segment if current one has reached endpoint */

        if (--(lp->randLine_cnt) <= 0) {
            next_random_lineseg(p, lp, n);
        }
    }
    /* someday, use aoutR for multimono out */

    *out1  = aoutL * outputGain;
    *out2 = aoutR * outputGain;
    return SP_OK;
}


int sp_rms_create(sp_rms **p)
{
    *p = malloc(sizeof(sp_rms));
    return SP_OK;
}

int sp_rms_destroy(sp_rms **p)
{
    free(*p);
    return SP_OK;
}

int sp_rms_init(sp_data *sp, sp_rms *p)
{
    p->ihp = 10;
    p->istor = 0;

    SPFLOAT b;

    b = 2.0 - cos((SPFLOAT)(p->ihp * (2 * M_PI / sp->sr)));
    p->c2 = b - sqrt(b*b - 1.0);
    p->c1 = 1.0 - p->c2;
    if (!p->istor) p->prvq = 0.0;
    return SP_OK;
}

int sp_rms_compute(sp_data *sp, sp_rms *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT q;
    SPFLOAT c1 = p->c1, c2 = p->c2;

    q = p->prvq;
    
    SPFLOAT as = *in;
    q = c1 * as * as + c2 * q;
    
    p->prvq = q;
    *out = sqrt(q);
    return SP_OK;
}


static int sp_rpt_set(sp_rpt *p, SPFLOAT bpm, int div, int rep);

int sp_rpt_create(sp_rpt **p)
{
    *p = malloc(sizeof(sp_rpt));
    return SP_OK;
}

int sp_rpt_destroy(sp_rpt **p)
{
    sp_rpt *pp = *p;
    sp_auxdata_free(&pp->aux);
    free(*p);
    return SP_OK;
}

int sp_rpt_init(sp_data *sp, sp_rpt *p, SPFLOAT maxdur)
{
    sp_auxdata_alloc(&p->aux, sizeof(SPFLOAT) * (uint32_t)maxdur * sp->sr);
    p->playpos = 0;
    p->bufpos = 0;
    p->running = 0;
    p->reps = 4;
    p->count = p->reps;
    p->size = (int)p->aux.size;
    p->sr = sp->sr;
    p->bpm = 130;
    p->div = 4;
    p->rep = 4;
    p->rc = SP_OK;
    return SP_OK;
}

int sp_rpt_compute(sp_data *sp, sp_rpt *p, SPFLOAT *trig,
        SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT *buf = (SPFLOAT *)p->aux.ptr;
    if(p->rc == SP_NOT_OK) {
        *out = 0;
        return SP_NOT_OK;
    }
    if(*trig > 0){
        p->rc = sp_rpt_set(p, p->bpm, p->div, p->rep);
        p->running = 1;
        p->playpos = 0;
        p->bufpos = 0;
        p->count = p->reps + 1;
    }
    if(p->bufpos * sizeof(SPFLOAT) < p->aux.size){
        p->rc = sp_rpt_set(p, p->bpm, p->div, p->rep);
        buf[p->bufpos] = *in;
        p->bufpos++;
    }else{
        p->running = 0;
    }
    if(p->running && p->count > 0){
        if(p->playpos == 0){
            p->count--;
        }
        *out = buf[p->playpos];
        p->playpos = (p->playpos + 1) % p->size;
    }else{
        *out = *in;
    }
    return SP_OK;
}

static int sp_rpt_set(sp_rpt *p, SPFLOAT bpm, int div, int rep)
{
    uint32_t size = (p->sr * (60.0 / bpm)) / (SPFLOAT) div;
    p->reps = rep;
    if(size * sizeof(SPFLOAT) > p->aux.size){
        fprintf(stderr, "Error: not enough memory allocated for buffer.\n");
        return SP_NOT_OK;
    }else if(size <= 0){
        fprintf(stderr, "Error: Size cannot be zero.\n");
        return SP_NOT_OK;
    }else{
        p->size = size;
    }
    return SP_OK;
}


#define sp_oneUp31Bit      (4.656612875245796924105750827168e-10)
#define sp_randGab   ((SPFLOAT)     \
        (((p->holdrand = p->holdrand * 214013 + 2531011) >> 1)  \
         & 0x7fffffff) * sp_oneUp31Bit)
#define sp_BiRandGab ((SPFLOAT)     \
        (p->holdrand = p->holdrand * -214013 + 2531011) * sp_oneUp31Bit)

int sp_rspline_create(sp_rspline **p)
{
    *p = malloc(sizeof(sp_rspline));
    return SP_OK;
}

int sp_rspline_destroy(sp_rspline **p)
{
    free(*p);
    return SP_OK;
}

int sp_rspline_init(sp_data *sp, sp_rspline *p)
{
    p->holdrand = sp_rand(sp);
    p->num1 = sp_randGab;
    p->num2 = sp_randGab;
    p->df1 = 0.0;
    p->phs = 0.0;
    p->init = 1;
    p->cps_max = 3;
    p->cps_min = 0.1;
    p->min = 0;
    p->max = 1;
    p->onedsr = 1.0 / (SPFLOAT)sp->sr;
    return SP_OK;
}

int sp_rspline_compute(sp_data *sp, sp_rspline *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT x;
    SPFLOAT f0 = p->num0, df0 = p->df0;
    SPFLOAT phs = p->phs;
    SPFLOAT slope, resd1, resd0, f2, f1;

    /* the original source code used an init flag alongside a goto */
    /* I didn't want to use gotos, so I opted to use this approach */

    if(p->init) {
        p->si = (sp_randGab * (p->cps_max - p->cps_min) + p->cps_min)*p->onedsr;
        p->init = 1;
    }

    phs += p->si;
    if (phs >= 1.0 || p->init) {
        p->si = (sp_randGab * (p->cps_max - p->cps_min)+p->cps_min)*p->onedsr;
        while (phs > 1.0) phs -= 1.0;
        f0 = p->num0 = p->num1;
        f1 = p->num1 = p->num2;
        f2 = p->num2 = sp_BiRandGab;
        df0 = p->df0 = p->df1;
        p->df1 = (f2 - f0) * 0.5;
        slope = f1 - f0;
        resd0 = df0 - slope;
        resd1 = p->df1 - slope;
        p->c3 = resd0 + resd1;
        p->c2 = - (resd1 + 2.0 * resd0);
    }

    x = (SPFLOAT) phs;
    *out = (((p->c3 * x + p->c2) * x + df0) * x + f0) * 
        (p->max - p->min) + p->min;
    p->phs = phs;
    p->init = 0;
    return SP_OK;
}


int sp_samphold_create(sp_samphold **p)
{
    *p = malloc(sizeof(sp_samphold));
    return SP_OK;
}

int sp_samphold_destroy(sp_samphold **p)
{
    free(*p);
    return SP_OK;
}

int sp_samphold_init(sp_data *sp, sp_samphold *p)
{
    p->val = 0;
    return SP_OK;
}

int sp_samphold_compute(sp_data *sp, sp_samphold *p, SPFLOAT *trig, SPFLOAT *in, SPFLOAT *out)
{
    if(*trig != 0) {
        p->val = *in;
    }
    *out = p->val;
    return SP_OK;
}


static void bilinear_transform(SPFLOAT acoefs[], SPFLOAT dcoefs[], SPFLOAT fs)
{
    SPFLOAT b0, b1, b2, a0, a1, a2;
    SPFLOAT bz0, bz1, bz2, az0, az1, az2;

    b0 = acoefs[0]; b1 = acoefs[1]; b2 = acoefs[2];
    a0 = acoefs[3]; a1 = acoefs[4]; a2 = acoefs[5];

    bz0 = 1.0; bz1 = 0.0; bz2 = 0.0;
    az0 = 1.0; az1 = 0.0; az2 = 0.0;

    az0 = a2*4*fs*fs + a1*2*fs + a0;

    bz2 = (b2*4*fs*fs - b1*2*fs + b0) / az0;
    bz1 = (-b2*8*fs*fs + 2*b0) / az0;
    bz0 = (b2*4*fs*fs+ b1*2*fs + b0) / az0;
    az2 = (a2*4*fs*fs - a1*2*fs + a0) / az0;
    az1 = (-a2*8*fs*fs + 2*a0) / az0;

    dcoefs[0] = bz0; dcoefs[1] = bz1; dcoefs[2] = bz2;
    dcoefs[3] = az1; dcoefs[4] = az2;
}

int sp_saturator_create(sp_saturator **p)
{
    *p = malloc(sizeof(sp_saturator));
    return SP_OK;
}

int sp_saturator_destroy(sp_saturator **p)
{
    free(*p);
    return SP_OK;
}

int sp_saturator_init(sp_data *sp, sp_saturator *p)
{
    int i, j;
    const SPFLOAT aacoefs[6][7] =
    {
        {2.60687e-05, 2.98697e-05, 2.60687e-05, -1.31885, 0.437162, 0.0, 0.0},
        {1, -0.800256, 1, -1.38301, 0.496576, 0.0, 0.0},
        {1, -1.42083, 1, -1.48787, 0.594413, 0.0, 0.0},
        {1, -1.6374, 1, -1.60688, 0.707142, 0.0, 0.0},
        {1, -1.7261, 1, -1.7253, 0.822156, 0.0, 0.0},
        {1, -1.75999, 1, -1.84111, 0.938811, 0.0, 0.0}
    };

    SPFLOAT wc_dc = 5*2*M_PI;
    SPFLOAT scoeffs[6] = {  0, 1, 0, wc_dc, 1, 0 };
    SPFLOAT zcoeffs[5];
    p->drive = 1;
    p->dcoffset = 0;

    for(i = 0; i < 6; i++){
        for(j = 0; j < 7; j++){
            p->aa[i][j] =  aacoefs[i][j];
            p->ai[i][j] =  aacoefs[i][j];
        }
    }
    bilinear_transform(scoeffs, zcoeffs, sp->sr*8);
    for(i = 0; i < 2; i++){
        for(j = 0; j < 5; j++)
            p->dcblocker[i][j] = zcoeffs[j];
        p->dcblocker[i][5] = 0.0;
        p->dcblocker[i][6] = 0.0;
    }
        return SP_OK;
}

static int quad_compute(SPFLOAT p[7],  SPFLOAT *input, SPFLOAT* output)
{
    SPFLOAT in = *input;
    *output = p[5] + in * p[0];
    p[5] = p[6] + in * p[1] - *output*p[3];
    p[6] = in * p[2] - *output*p[4];
    return SP_OK;
}


int sp_saturator_compute(sp_data *sp, sp_saturator *p, SPFLOAT *in, SPFLOAT *out)
{
    int i, j;
    SPFLOAT fsignal, usignal, dsignal;

    fsignal = p->drive * *in;
    for(i = 0; i < 8; i++){
        usignal = (i == 0) ? 8 *fsignal : 0.0;
        for(j = 0; j < 6; j++)
            quad_compute(p->ai[j], &usignal, &usignal);

        dsignal = (usignal + p->dcoffset) / (1.0 + fabs(usignal + p->dcoffset));

        quad_compute(p->dcblocker[0], &dsignal, &dsignal);
        quad_compute(p->dcblocker[1], &dsignal, &dsignal);

        for(j = 0; j < 6; j++)
            quad_compute(p->aa[j], &dsignal, out);
    }
    return SP_OK;
}


int sp_scale_create(sp_scale **p)
{
    *p = malloc(sizeof(sp_scale));
    return SP_OK;
}

int sp_scale_destroy(sp_scale **p)
{
    free(*p);
    return SP_OK;
}

int sp_scale_init(sp_data *sp, sp_scale *p)
{
    p->min = -1;
    p->max = 1;
    return SP_OK;
}

int sp_scale_compute(sp_data *sp, sp_scale *p, SPFLOAT *in, SPFLOAT *out)
{
    *out =  *in * (p->max - p->min) + p->min;
    return SP_OK;
}


int sp_gen_scrambler(sp_data *sp, sp_ftbl *src, sp_ftbl **dest)
{

    uint32_t size = (src->size % 2 == 0) ? src->size : src->size - 1;
    sp_ftbl *dst;
    sp_ftbl_create(sp, &dst, size);
    kiss_fftr_cfg fft, ifft;
    kiss_fft_cpx *tmp;

    /* set up kissfft */
    fft = kiss_fftr_alloc(size, 0, NULL, NULL);
    ifft = kiss_fftr_alloc(size, 1, NULL, NULL);
    tmp = malloc(sizeof(kiss_fft_cpx) * size);
    memset(tmp, 0, sizeof(SPFLOAT) * size);
    kiss_fftr(fft, src->tbl, tmp);

    uint32_t i;
    SPFLOAT mag, phs;
    for(i = 0; i < size / 2; i++) {
        mag = sqrt(tmp[i].r * tmp[i].r + tmp[i].i * tmp[i].i) / size;
        phs = ((SPFLOAT)sp_rand(sp) / SP_RANDMAX) * 2 * M_PI;
        tmp[i].r = mag * cos(phs);
        tmp[i].i = mag * sin(phs);
    }

    tmp[0].r = 0;
    tmp[0].i = 0;
    tmp[size / 2 - 1].r = 0;
    tmp[size / 2 - 1].i = 0;

    kiss_fftri(ifft, tmp, dst->tbl);
    SPFLOAT max = -1;
    SPFLOAT val = 0;
    for(i = 0; i < size; i++) {
        val = fabs(dst->tbl[i]); 
        if(val > max) {
            max = val;
        }
    }

    for(i = 0; i < size; i++) {
       dst->tbl[i] /= max;
    }

    kiss_fftr_free(fft);
    kiss_fftr_free(ifft);
    KISS_FFT_FREE(tmp);
    
    *dest = dst;
    return SP_OK;
}

int sp_sdelay_create(sp_sdelay **p)
{
    *p = malloc(sizeof(sp_sdelay));
    sp_sdelay *pp = *p;
    pp->size = 0;
    return SP_OK;
}

int sp_sdelay_destroy(sp_sdelay **p)
{
    sp_sdelay *pp = *p;

    if(pp->size > 0) {
        free(pp->buf);
    }

    free(*p);
    return SP_OK;
}

int sp_sdelay_init(sp_data *sp, sp_sdelay *p, int size)
{
    int n;
    p->size = size;
    p->buf = malloc(size * sizeof(SPFLOAT));
    for(n = 0; n < p->size; n++) p->buf[n] = 0;
    p->pos = 0;
    return SP_OK;
}

int sp_sdelay_compute(sp_data *sp, sp_sdelay *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = p->buf[p->pos];
    p->buf[p->pos] = *in;
    p->pos = (p->pos + 1) % p->size;
    return SP_OK;
}


int sp_slice_create(sp_slice **p)
{
    *p = malloc(sizeof(sp_slice));
    return SP_OK;
}

int sp_slice_destroy(sp_slice **p)
{
    free(*p);
    return SP_OK;
}

int sp_slice_init(sp_data *sp, sp_slice *p, sp_ftbl *vals, sp_ftbl *buf)
{
    p->vals = vals;
    p->buf = buf;
    p->pos = 0;
    p->nextpos = 0;
    p->id = 0;
    return SP_OK;
}

int sp_slice_compute(sp_data *sp, sp_slice *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0;
    if(*in != 0) {
        if(p->id < p->vals->size) {
            p->pos = p->vals->tbl[p->id];
            if(p->id == p->vals->size - 1) {
                p->nextpos = p->buf->size;
            } else {
                p->nextpos = p->vals->tbl[p->id + 1];
            }
        }
    }

    if(p->pos < p->nextpos) {
        *out = p->buf->tbl[p->pos];
        p->pos++;
    }

    return SP_OK;
}


int sp_smoothdelay_create(sp_smoothdelay **p)
{
    *p = malloc(sizeof(sp_smoothdelay));
    return SP_OK;
}

int sp_smoothdelay_destroy(sp_smoothdelay **p)
{
    sp_smoothdelay *pp = *p;
    sp_auxdata_free(&pp->buf1);
    sp_auxdata_free(&pp->buf2);
    free(*p);
    return SP_OK;
}

int sp_smoothdelay_init(sp_data *sp, sp_smoothdelay *p, 
        SPFLOAT maxdel, uint32_t interp)
{
    uint32_t n = (int32_t)(maxdel * sp->sr)+1;
    p->sr = sp->sr;
    p->del = maxdel * 0.5;
    p->pdel = -1;
    p->maxdel = maxdel;
    p->feedback = 0;
    p->maxbuf = n - 1;
    p->maxcount = interp;

    sp_auxdata_alloc(&p->buf1, n * sizeof(SPFLOAT));
    p->bufpos1 = 0;
    p->deltime1 = (uint32_t) (p->del * sp->sr);

    sp_auxdata_alloc(&p->buf2, n * sizeof(SPFLOAT));
    p->bufpos2 = 0;
    p->deltime2 = p->deltime1;

    p->counter = 0;
    p->curbuf = 0;
    return SP_OK;
}

static SPFLOAT delay_sig(SPFLOAT *buf, 
        uint32_t *bufpos, 
        uint32_t deltime, 
        SPFLOAT fdbk, 
        SPFLOAT in)
{
    SPFLOAT delay = buf[*bufpos];
    SPFLOAT sig = (delay * fdbk) + in;
    buf[*bufpos] = sig;
    *bufpos = (*bufpos + 1) % deltime;
    return delay;
}

int sp_smoothdelay_compute(sp_data *sp, sp_smoothdelay *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0;
    if(p->del != p->pdel && p->counter == 0) {
        uint32_t dels = min((uint32_t)(p->del * sp->sr), p->maxbuf);

        /* initial delay time sets time for both buffers */

        if(p->pdel < 0) {
            p->deltime1 = dels;
            p->deltime2 = dels;
        }

        p->pdel = p->del;

        if(dels == 0) dels = 1;

        if(p->curbuf == 0) {
            p->curbuf = 1;
            p->deltime2 = dels;
        } else {
            p->curbuf = 0;
            p->deltime1 = dels;
        }
        p->counter = p->maxcount;
    }



    SPFLOAT *buf1 = (SPFLOAT *)p->buf1.ptr; 
    SPFLOAT *buf2 = (SPFLOAT *)p->buf2.ptr; 
    SPFLOAT it = (SPFLOAT)p->counter / p->maxcount;
    if(p->counter != 0) p->counter--;
  
    SPFLOAT del1 = delay_sig(buf1, &p->bufpos1, 
            p->deltime1, p->feedback, *in);

    SPFLOAT del2 = delay_sig(buf2, &p->bufpos2, 
            p->deltime2, p->feedback, *in);

    if(p->curbuf == 0) {
        /* 1 to 2 */
        *out = (del1 * it) + (del2 * (1 - it));
    } else {
        /* 2 to 1 */
        *out = (del2 * it) + (del1 * (1 - it));
    }
    return SP_OK;
}


#define SPA_BUFSIZE 4096

int sp_spa_create(sp_spa **p)
{
    *p = malloc(sizeof(sp_spa));
    return SP_OK;
}

int sp_spa_destroy(sp_spa **p)
{
    sp_spa *pp = *p;
    sp_auxdata_free(&pp->aux);
    spa_close(&pp->spa);
    free(*p);
    return SP_OK;
}

int sp_spa_init(sp_data *sp, sp_spa *p, const char *filename)
{
    if(spa_open(sp, &p->spa, filename, SPA_READ) != SP_OK) {
        return SP_NOT_OK;
    }
    
    p->pos = 0;

    p->bufsize = SPA_BUFSIZE;
    sp_auxdata_alloc(&p->aux, sizeof(SPFLOAT) * p->bufsize);

    p->buf = p->aux.ptr;

    return SP_OK;
}

int sp_spa_compute(sp_data *sp, sp_spa *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->bufsize == 0) {
        *out = 0.0;
        return SP_OK;
    }

    if(p->pos == 0) {
        p->bufsize = spa_read_buf(sp, &p->spa, p->buf, SPA_BUFSIZE);
        if(p->bufsize == 0) {
            *out = 0.0;
            return SP_OK;
        }
    }

    *out = p->buf[p->pos];
    p->pos = (p->pos + 1) % p->bufsize;
    return SP_OK;
}

int sp_sparec_create(sp_sparec **p)
{
    *p = malloc(sizeof(sp_sparec));
    return SP_OK;
}

int sp_sparec_destroy(sp_sparec **p)
{
    sp_sparec *pp = *p;
    sp_auxdata_free(&pp->aux);
    spa_close(&pp->spa);
    free(*p);
    return SP_OK;
}

int sp_sparec_init(sp_data *sp, sp_sparec *p, const char *filename)
{
    if(spa_open(sp, &p->spa, filename, SPA_WRITE) != SP_OK) {
        return SP_NOT_OK;
    }
    
    p->pos = SPA_BUFSIZE;

    p->bufsize = SPA_BUFSIZE;
    sp_auxdata_alloc(&p->aux, sizeof(SPFLOAT) * p->bufsize);

    p->buf = p->aux.ptr;
    return SP_OK;
}

int sp_sparec_compute(sp_data *sp, sp_sparec *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->pos == 0) {
        p->pos = p->bufsize;
        spa_write_buf(sp, &p->spa, p->buf, p->bufsize);
    }

    p->buf[p->bufsize - p->pos] = *in;

    p->pos--;
    *out = *in;
    return SP_OK;
}

/* call this to close sparec. will write the rest of the buffer */
int sp_sparec_close(sp_data *sp, sp_sparec *p)
{
    if(p->pos < p->bufsize - 1) {
        spa_write_buf(sp, &p->spa, p->buf, p->bufsize - p->pos);
    }
    return SP_OK;
}


int sp_streson_create(sp_streson **p) 
{
    *p = malloc(sizeof(sp_streson));
    return SP_OK;
}

int sp_streson_destroy(sp_streson **p) 
{
    sp_streson *pp = *p;
    sp_auxdata_free(&pp->buf);
    free(*p);
    return SP_OK;
}

int sp_streson_init(sp_data *sp, sp_streson *p) 
{
    int n;
    p->freq = 440.0;
    p->fdbgain = 0.8;
    p->size = (int) (sp->sr/20);   /* size of delay line */
    sp_auxdata_alloc(&p->buf, p->size * sizeof(SPFLOAT));
    p->Cdelay = (SPFLOAT*) p->buf.ptr; /* delay line */
    p->LPdelay = p->APdelay = 0.0; /* reset the All-pass and Low-pass delays */
    p->wpointer = p->rpointer = 0; /* reset the read/write pointers */
    for (n = 0; n < p->size; n++){
      p->Cdelay[n] = 0.0;
    }
    return SP_OK;
}

int sp_streson_compute(sp_data *sp, sp_streson *p, SPFLOAT *in, SPFLOAT *out) 
{
    SPFLOAT g = p->fdbgain;
    SPFLOAT freq;
    SPFLOAT a, s, w, sample, tdelay, fracdelay;
    int delay;
    int rp = p->rpointer, wp = p->wpointer;
    int size = p->size;
    SPFLOAT APdelay = p->APdelay;
    SPFLOAT LPdelay = p->LPdelay;
    int vdt;

    freq = p->freq;
    if (freq < 20.0) freq = 20.0;   /* lowest freq is 20 Hz */
    tdelay = sp->sr/freq;
    delay = (int) (tdelay - 0.5); /* comb delay */
    fracdelay = tdelay - (delay + 0.5); /* fractional delay */
    vdt = size - delay;       /* set the var delay */
    a = (1.0-fracdelay)/(1.0+fracdelay);   /* set the all-pass gain */
    
    SPFLOAT tmpo;
    rp = (vdt + wp);
    if (rp >= size) rp -= size;
    tmpo = p->Cdelay[rp];
    w = *in + tmpo;
    s = (LPdelay + w)*0.5;
    LPdelay = w;
    *out = sample = APdelay + s*a;
    APdelay = s - (sample*a);
    p->Cdelay[wp] = sample*g;
    wp++;
    if (wp == size) wp=0;
    p->rpointer = rp; p->wpointer = wp;
    p->LPdelay = LPdelay; p->APdelay = APdelay;
    return SP_OK;
}


int sp_switch_create(sp_switch **p)
{
    *p = malloc(sizeof(sp_switch));
    return SP_OK;
}

int sp_switch_destroy(sp_switch **p)
{
    free(*p);
    return SP_OK;
}

int sp_switch_init(sp_data *sp, sp_switch *p)
{
    p->mode = 0;
    return SP_OK;
}

int sp_switch_compute(sp_data *sp, sp_switch *p, SPFLOAT *trig,
    SPFLOAT *in1, SPFLOAT *in2, SPFLOAT *out)
{
    if (*trig) {
        p->mode = p->mode == 0 ? 1 : 0;
    }

    if(p->mode == 0) {
        *out = *in1;
    } else {
        *out = *in2;
    }

    return SP_OK;
}



int sp_tabread_create(sp_tabread **p)
{
    *p = malloc(sizeof(sp_tabread));
    return SP_OK;
}

int sp_tabread_destroy(sp_tabread **p)
{
    free(*p);
    return SP_OK;
}

int sp_tabread_init(sp_data *sp, sp_tabread *p, sp_ftbl *ft, int mode)
{
    p->ft = ft;
    p->mode = (SPFLOAT) mode;
    p->offset = 0;
    p->wrap = 0;
    return SP_OK;
}

int sp_tabread_compute(sp_data *sp, sp_tabread *p, SPFLOAT *in, SPFLOAT *out)
{
    int ndx, len = (int)p->ft->size;
    int32_t mask = (int)p->ft->size - 1;
    SPFLOAT index = p->index;
    SPFLOAT *tbl = p->ft->tbl;
    SPFLOAT offset = p->offset;
    SPFLOAT mul = 1, tmp, frac;

    if (p->mode) {
        mul = p->ft->size;
    }else {
        p->mul = 1;
    }

    int32_t iwrap = (int32_t)p->wrap;

    SPFLOAT x1, x2;
    tmp = (index + offset) * mul;
    ndx = floor(tmp);
    frac = tmp - ndx;
    if (iwrap) {
        if ((mask ? 0 : 1)) {
            while(ndx >= len) ndx -= len;
            while(ndx < 0)  ndx += len;
        }
        else ndx &= mask;
    } else {
        if (ndx >= len) ndx = len - 1;
        else if (ndx < 0) ndx = 0;
    }

    x1 = tbl[ndx];
    x2 = tbl[ndx+1];
    *out = x1 + (x2 - x1) * frac;
    return SP_OK;
}


static void make_Envelope(sp_tadsr *e)
{
    e->target = 0.0;
    e->value = 0.0;
    e->rate = 0.001;
    e->state = 1;
}

static void make_ADSR(sp_tadsr *a)
{
    make_Envelope(a);
    a->target = 0.0;
    a->value = 0.0;
    a->attackRate = 0.001;
    a->decayRate = 0.001;
    a->sustainLevel = 0.0;
    a->releaseRate = 0.01;
    a->state = CLEAR;
}

static void ADSR_keyOn(sp_tadsr *a)
{
    a->target = 1.0;
    a->rate = a->attackRate;
    a->state = ATTACK;
}

static void ADSR_keyOff(sp_tadsr *a)
{
    a->target = 0.0;
    a->rate = a->releaseRate;
    a->state = RELEASE;
}

static void ADSR_setSustainLevel(sp_data *sp, sp_tadsr *a, SPFLOAT aLevel)
{
   a->sustainLevel = aLevel;
}

static void ADSR_setAttackTime(sp_data *sp, sp_tadsr *a, SPFLOAT aTime)
{
    a->attackRate = 1.0 / (aTime * sp->sr);
}

static void ADSR_setDecayTime(sp_data *sp, sp_tadsr *a, SPFLOAT aTime)
{
    a->decayRate = 1.0 / (aTime * sp->sr);
}

static void ADSR_setReleaseTime(sp_data *sp, sp_tadsr *a, SPFLOAT aTime)
{
    a->releaseRate = 1.0 / (aTime * sp->sr);
}

static void ADSR_setAllTimes(sp_data *sp, sp_tadsr *a, SPFLOAT attTime, SPFLOAT decTime,
                      SPFLOAT susLevel, SPFLOAT relTime)
{
    ADSR_setAttackTime(sp, a, attTime);
    ADSR_setDecayTime(sp, a, decTime);
    ADSR_setSustainLevel(sp, a, susLevel);
    ADSR_setReleaseTime(sp, a, relTime);
}

static SPFLOAT ADSR_tick(sp_tadsr *a)
{
    SPFLOAT out = 0;
    if (a->state == ATTACK) {
        a->value += a->rate;
        if (a->value >= a->target) {
            a->value = a->target;
            a->rate = a->decayRate;
            a->target = a->sustainLevel;
            a->state = DECAY;
        }
        out = a->value;
    } else if (a->state == DECAY) {
        a->value -= a->decayRate;
        out = a->value;
        if (a->value <= a->sustainLevel) {
            a->value = a->sustainLevel;
            out = a->sustainLevel;
            a->rate = 0.0;
            a->state = SUSTAIN;
        }
    } else if (a->state == RELEASE)  {
        a->value -= a->releaseRate;
        if (a->value <= 0.0) {
            a->value = 0.0;
            a->state = CLEAR;
        }
        out = a->value;
    } else if (a->state == SUSTAIN)  {
        out = a->sustainLevel;
    }
    return out;
}

int sp_tadsr_create(sp_tadsr **p)
{
    *p = malloc(sizeof(sp_tadsr));
    return SP_OK;
}

int sp_tadsr_destroy(sp_tadsr **p)
{
    free(*p);
    return SP_OK;
}

int sp_tadsr_init(sp_data *sp, sp_tadsr *p)
{
    make_ADSR(p);
    p->atk = 0.5;
    p->dec = 0.5;
    p->sus = 0.0;
    p->rel = 0.5;
    p->mode = KEY_OFF;
    return SP_OK;
}

int sp_tadsr_compute(sp_data *sp, sp_tadsr *p, SPFLOAT *trig, SPFLOAT *out)
{
    if(*trig != 0) {

        if(*trig == 2) {
            ADSR_keyOff(p);
            p->mode = KEY_OFF;
        }else if(p->mode == KEY_OFF) {
            ADSR_setAllTimes(sp, p, p->atk, p->dec, p->sus, p->rel);
            ADSR_keyOn(p);
            p->mode = KEY_ON;
        } else {
            ADSR_keyOff(p);
            p->mode = KEY_OFF;
        }
    }
    *out = ADSR_tick(p);
    return SP_OK;
}



#ifndef TWO_PI
#define TWO_PI 6.28318530717958647692528676655901
#endif

#define ORD_MAX 50

static void lpc_durbin(SPFLOAT *r, int p, float *k, float *g)
{
  int i, j;
  SPFLOAT a[ORD_MAX], at[ORD_MAX], e=r[0];
    
  for(i=0; i<=p; i++) a[i] = at[i] = 0.0f;

  for(i=1; i<=p; i++) 
  {
    k[i] = -r[i];

    for(j=1; j<i; j++) 
    { 
      at[j] = a[j];
      k[i] -= a[j] * r[i-j]; 
    }
    if(fabs(e) < 1.0e-20f) { e = 0.0f;  break; }
    k[i] /= e;
    
    a[i] = k[i];
    for(j=1; j<i; j++) a[j] = at[j] + k[i] * at[i-j];
    
    e *= 1.0f - k[i] * k[i];
  }

  if(e < 1.0e-20f) e = 0.0f;
  *g = (float)sqrt(e);
}

static void lpc(float *buf, float *car, uint32_t n, uint32_t o)
{
    SPFLOAT z[ORD_MAX], r[ORD_MAX], k[ORD_MAX], G, x;
    uint32_t i, j, nn=n;
    SPFLOAT min;

    /* buf[] is already emphasized and windowed */
    for(j=0; j<=o; j++, nn--) {
        z[j] = r[j] = 0.0f;
        for(i=0; i<nn; i++) r[j] += buf[i] * buf[i+j]; /* autocorrelation */
    }

    r[0] *= 1.001f;  /* stability fix */
    
    min = 0.00001f;
    if(r[0] < min) { for(i=0; i<n; i++) buf[i] = 0.0f; return; } 

    lpc_durbin(r, o, k, &G);  /* calc reflection coeffs */

    for(i=1; i<=o; i++) {
        if(k[i] > 0.995f) k[i] = 0.995f; else if(k[i] < -0.995f) k[i] = -.995f;
    }

    for(i=0; i<n; i++) {
        x = G * car[i];
        /* lattice filter */
        for(j=o; j>0; j--) { 
            x -= k[j] * z[j-1];
            z[j] = z[j-1] + k[j] * x;
        }
        buf[i] = z[0] = x;  /* output buf[] will be windowed elsewhere */
    }
}

int sp_talkbox_create(sp_talkbox **p)
{
    *p = malloc(sizeof(sp_talkbox));
    return SP_OK;
}

int sp_talkbox_destroy(sp_talkbox **p)
{
    free(*p);
    return SP_OK;
}

int sp_talkbox_init(sp_data *sp, sp_talkbox *p)
{
    uint32_t n;
    p->quality = 1.f;
    p->N = 1;
    p->K = 0;

    n = (uint32_t)(0.01633f * sp->sr);
    if(n > SP_TALKBOX_BUFMAX) n = SP_TALKBOX_BUFMAX;

    /* calculate hanning window */
    if(n != p->N) {
        p->N = n;
        SPFLOAT dp = TWO_PI / (SPFLOAT)p->N;
        SPFLOAT pos = 0.0f;
        for(n=0; n < p->N; n++) {
            p->window[n] = 0.5f - 0.5f * (SPFLOAT)cos(pos);
            pos += dp;
        }
    }
 
    /* zero out variables and buffers */
    p->pos = p->K = 0;
    p->emphasis = 0.0f;
    p->FX = 0;

    p->u0 = p->u1 = p->u2 = p->u3 = p->u4 = 0.0f;
    p->d0 = p->d1 = p->d2 = p->d3 = p->d4 = 0.0f;

    memset(p->buf0, 0, SP_TALKBOX_BUFMAX * sizeof(SPFLOAT));
    memset(p->buf1, 0, SP_TALKBOX_BUFMAX * sizeof(SPFLOAT));
    memset(p->car0, 0, SP_TALKBOX_BUFMAX * sizeof(SPFLOAT));
    memset(p->car1, 0, SP_TALKBOX_BUFMAX * sizeof(SPFLOAT));
    return SP_OK;
}

int sp_talkbox_compute(sp_data *sp, sp_talkbox *t, SPFLOAT *src, SPFLOAT *exc, SPFLOAT *out)
{
    int32_t p0=t->pos, p1 = (t->pos + t->N/2) % t->N;
    SPFLOAT e=t->emphasis, w, o, x, fx=t->FX;
    SPFLOAT p, q, h0=0.3f, h1=0.77f;
    SPFLOAT den;

    t->O = (uint32_t)((0.0001f + 0.0004f * t->quality) * sp->sr);

    o = *src;
    x = *exc;

    p = t->d0 + h0 * x; 
    t->d0 = t->d1;  
    t->d1 = x - h0 * p;

    q = t->d2 + h1 * t->d4;
    t->d2 = t->d3;
    t->d3 = t->d4 - h1 * q;  

    t->d4 = x;

    x = p + q;

    if(t->K++) {
        t->K = 0;

        /* carrier input */
        t->car0[p0] = t->car1[p1] = x;

        /* 6dB/oct pre-emphasis */
        x = o - e;  e = o;  

        /* 50% overlapping hanning windows */
        w = t->window[p0]; fx = t->buf0[p0] * w;  t->buf0[p0] = x * w;  
        if(++p0 >= t->N) { lpc(t->buf0, t->car0, t->N, t->O);  p0 = 0; }

        w = 1.0f - w;  fx += t->buf1[p1] * w;  t->buf1[p1] = x * w;
        if(++p1 >= t->N) { lpc(t->buf1, t->car1, t->N, t->O);  p1 = 0; }
    }

    p = t->u0 + h0 * fx; 
    t->u0 = t->u1;  
    t->u1 = fx - h0 * p;

    q = t->u2 + h1 * t->u4;
    t->u2 = t->u3;
    t->u3 = t->u4 - h1 * q;  

    t->u4 = fx;
    x = p + q;

    o = x * 0.5;
    *out = o;
    t->emphasis = e;
    t->pos = p0;
    t->FX = fx;

    den = 1.0e-10f; 
    /* anti-denormal */
    if(fabs(t->d0) < den) t->d0 = 0.0f; 
    if(fabs(t->d1) < den) t->d1 = 0.0f;
    if(fabs(t->d2) < den) t->d2 = 0.0f;
    if(fabs(t->d3) < den) t->d3 = 0.0f;
    if(fabs(t->u0) < den) t->u0 = 0.0f;
    if(fabs(t->u1) < den) t->u1 = 0.0f;
    if(fabs(t->u2) < den) t->u2 = 0.0f;
    if(fabs(t->u3) < den) t->u3 = 0.0f;
    return SP_OK;
}


int sp_tblrec_create(sp_tblrec **p)
{
    *p = malloc(sizeof(sp_tblrec));
    return SP_OK;
}

int sp_tblrec_destroy(sp_tblrec **p)
{
    free(*p);
    return SP_OK;
}

int sp_tblrec_init(sp_data *sp, sp_tblrec *p, sp_ftbl *ft)
{
    p->index = 0;
    p->record = 0;
    p->ft = ft;
    return SP_OK;
}

int sp_tblrec_compute(sp_data *sp, sp_tblrec *p, SPFLOAT *in, SPFLOAT *trig, SPFLOAT *out)
{
    if(*trig != 0) {
        if(p->record == 1) {
            p->record = 0;
        } else {
            p->record = 1;
            p->index = 0;
            memset(p->ft->tbl, 0, sizeof(SPFLOAT) * p->ft->size);
        }
    }

    if(p->record) {
        p->ft->tbl[p->index] = *in;
        p->index = (p->index + 1) % p->ft->size;
    }

    *out = *in;
    return SP_OK;
}


int sp_tbvcf_create(sp_tbvcf **p)
{
    *p = malloc(sizeof(sp_tbvcf));
    return SP_OK;
}

int sp_tbvcf_destroy(sp_tbvcf **p)
{
    free(*p);
    return SP_OK;
}

int sp_tbvcf_init(sp_data *sp, sp_tbvcf *p)
{
    p->fco = 500.0;
    p->res = 0.8;
    p->dist = 2.0;
    p->asym = 0.5;

    p->sr = sp->sr;
    p->onedsr = 1.0 / p->sr;

    p->iskip = 0.0;
    if(p->iskip == 0.0){
        p->y = p->y1 = p->y2 = 0.0;
    }
    p->fcocod = p->fco;
    p->rezcod = p->res;


    return SP_OK;
}

/* TODO: clean up code here. */
int sp_tbvcf_compute(sp_data *sp, sp_tbvcf *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT x;
    SPFLOAT fco, res, dist, asym;
    SPFLOAT y = p->y, y1 = p->y1, y2 = p->y2;
    /* The initialisations are fake to fool compiler warnings */
    SPFLOAT ih, fdbk, d, ad;
    SPFLOAT fc=0.0, fco1=0.0, q=0.0, q1=0.0;

    ih  = 0.001; /* ih is the incremental factor */

 /* Set up the pointers
    fcoptr  = p->fco;
    resptr  = p->res;
    distptr = p->dist;
    asymptr = p->asym; */

 /* Get the values for the k-rate variables
    fco  = (SPFLOAT)*fcoptr;
    res  = (SPFLOAT)*resptr;
    dist = (SPFLOAT)*distptr;
    asym = (SPFLOAT)*asymptr; */

    /* This should work in sp world */
    fco = p->fco;
    res = p->res;
    dist = p->dist;
    asym = p->asym;

 /* Try to decouple the variables */
    if ((p->rezcod==0) && (p->fcocod==0)) { /* Calc once only */
        q1   = res/(1.0 + sqrt(dist));
        fco1 = pow(fco*260.0/(1.0+q1*0.5),0.58);
        q    = q1*fco1*fco1*0.0005;
        fc   = fco1*p->onedsr*(p->sr/8.0);
    }
    if ((p->rezcod!=0) || (p->fcocod!=0)) {
        q1  = res/(1.0 + sqrt(dist));
        fco1 = pow(fco*260.0/(1.0+q1*0.5),0.58);
        q  = q1*fco1*fco1*0.0005;
        fc  = fco1*p->onedsr*(p->sr/8.0);
    }
    x  = *in;
    fdbk = q*y/(1.0 + exp(-3.0*y)*asym);
    y1  = y1 + ih*((x - y1)*fc - fdbk);
    d  = -0.1*y*20.0;
    ad  = (d*d*d + y2)*100.0*dist;
    y2  = y2 + ih*((y1 - y2)*fc + ad);
    y  = y + ih*((y2 - y)*fc);
    *out = (y*fc/1000.0*(1.0 + q1)*3.2);

    p->y = y; p->y1 = y1; p->y2 = y2;
    return SP_OK;
}


int sp_tdiv_create(sp_tdiv **p)
{
    *p = malloc(sizeof(sp_tdiv));
    return SP_OK;
}

int sp_tdiv_destroy(sp_tdiv **p)
{
    free(*p);
    return SP_OK;
}

int sp_tdiv_init(sp_data *sp, sp_tdiv *p)
{
    p->num = 2;
    p->counter = 0;
    p->offset = 0;
    return SP_OK;
}

int sp_tdiv_compute(sp_data *sp, sp_tdiv *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0.0;
    if(*in != 0) {
        if(p->counter == p->offset) *out = 1.0;
        else *out = 0.0;
        p->counter = (p->counter + 1) % p->num;
    }
    return SP_OK;
}



int sp_tenv_create(sp_tenv **p)
{
    *p = malloc(sizeof(sp_tenv));
    sp_tenv *pp = *p;
    sp_tevent_create(&pp->te);
    return SP_OK;
}

int sp_tenv_destroy(sp_tenv **p)
{
    sp_tenv *pp = *p;
    sp_tevent_destroy(&pp->te);
    free(*p);
    return SP_OK;
}

static void sp_tenv_reinit(void *ud)
{
    sp_tenv *env = ud;
    env->pos = 0;
    env->atk_end = env->sr * env->atk;
    env->rel_start = env->sr * (env->atk + env->hold);
    env->atk_slp = 1.0 / env->atk_end;
    env->rel_slp = -1.0 / (env->sr * env->rel);
    env->totaldur = env->sr * (env->atk + env->hold + env->rel);
}

static void sp_tenv_comp(void *ud, SPFLOAT *out)
{
    sp_tenv *env = ud;
    SPFLOAT sig = 0;
    uint32_t pos = env->pos;
    *out = 0.0;
    if(pos < env->atk_end){
        sig = env->last + env->atk_slp;
    }else if (pos < env->rel_start){
        sig = 1.0;
    }else if (pos < env->totaldur){
        sig = env->last + env->rel_slp;
    }else{
        sig = 0.0;
    }
    sig = (sig > 1.0) ? 1.0 : sig;
    sig = (sig < 0.0) ? 0.0 : sig;

    /* Internal input signal mode */
    if(env->sigmode) {
        *out = env->input * sig;
    } else {
        *out = sig;
    }


    env->pos++;
    env->last = sig;
}

int sp_tenv_init(sp_data *sp, sp_tenv *p)
{
    p->pos = 0;
    p->last = 0;
    p->atk = 0.1;
    p->hold = 0.3;
    p->rel = 0.2;
    p->sigmode = 0;
    p->input = 0;

    p->sr = sp->sr;
    p->atk_end = p->sr * p->atk;
    p->rel_start = p->sr * (p->atk + p->hold);
    p->atk_slp = 1.0 / p->atk_end;
    p->rel_slp = -1.0 / (p->sr * p->rel);
    p->totaldur = p->sr * (p->atk + p->hold + p->rel);
    sp_tevent_init(sp, p->te, sp_tenv_reinit, sp_tenv_comp, p);
    return SP_OK;
}

int sp_tenv_compute(sp_data *sp, sp_tenv *p, SPFLOAT *in, SPFLOAT *out)
{
    sp_tevent_compute(sp, p->te, in, out);
    return SP_OK;
}


enum {
    T_ON,
    T_OFF,
    T_INIT
};

int sp_tenv2_create(sp_tenv2 **p)
{
    *p = malloc(sizeof(sp_tenv2));
    return SP_OK;
}

int sp_tenv2_destroy(sp_tenv2 **p)
{
    free(*p);
    return SP_OK;
}

int sp_tenv2_init(sp_data *sp, sp_tenv2 *p)
{
    p->state = T_INIT;
    p->atk = 0.1;
    p->rel = 0.1;
    p->slope = 0;
    p->last = 0;
    return SP_OK;
}

int sp_tenv2_compute(sp_data *sp, sp_tenv2 *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in != 0) {
        if(p->state == T_INIT || p->state == T_OFF) {
            p->state = T_ON;
            p->timer = (uint32_t)(sp->sr * p->atk);
            p->totaltime = p->timer;
            p->slope = 1.0 / p->totaltime;
        } else if (p->state == T_ON) { 
            p->state = T_OFF;
            p->timer = (uint32_t)(sp->sr * p->rel);
            p->totaltime = p->timer;
            p->slope = 1.0 / p->totaltime;
        }
    }
    if(p->timer == 0) {
        if(p->state == T_ON) *out = 1;
        else *out = 0;
        return SP_OK;
    } else {
        p->timer--;
        if(p->state == T_ON)  {
            *out = p->last + p->slope;
        } else if (p->state == T_OFF) {
            *out = p->last - p->slope;
        }

        if(*out > 1) *out = 1;
        if(*out < 0) *out = 0;

        p->last = *out;

        return SP_OK;
    }
    return SP_OK;
}

int sp_tenvx_create(sp_tenvx **p)
{
    *p = malloc(sizeof(sp_tenvx));
    return SP_OK;
}

int sp_tenvx_destroy(sp_tenvx **p)
{
    free(*p);
    return SP_OK;
}


int sp_tenvx_init(sp_data *sp, sp_tenvx *p)
{
    p->hold = 0.5;
    p->atk = 0.5;
    p->rel = 0.5;
    p->a_a = p->b_a = 0;
    p->a_r = p->b_r = 0;
    p->y = 0;
    p->count = (uint32_t) (p->hold * sp->sr);
    return SP_OK;
}


int sp_tenvx_compute(sp_data *sp, sp_tenvx *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0;

    if(*in != 0) {
        p->a_a = exp(-1.0/(p->atk * sp->sr));
        p->b_a = 1.0 - p->a_a;
        p->a_r = exp(-1.0/(p->rel * sp->sr));
        p->b_r = 1.0 - p->a_r;
        p->count = (uint32_t) (p->hold * sp->sr);
    }

    if(p->count > 0) {
        *out = p->b_a + p->a_a * p->y;
        p->y = *out;
        p->count--;
    } else {
        *out = p->a_r * p->y;
        p->y = *out;
    }

    return SP_OK;
}


int sp_tevent_create(sp_tevent **te)
{
    *te = malloc(sizeof(sp_tevent));
    return SP_NOT_OK;
}

int sp_tevent_destroy(sp_tevent **te)
{
    free(*te);
    return SP_NOT_OK;
}

int sp_tevent_init(sp_data *sp, sp_tevent *te, 
        void (*reinit)(void*), void (*compute)(void *, SPFLOAT *out), void *ud)
{
    te->reinit = reinit;
    te->compute = compute;
    te->ud = ud;
    te->started = 0;
    return SP_OK;
}

int sp_tevent_compute(sp_data *sp, sp_tevent *te, SPFLOAT *in, SPFLOAT *out)
{
    if(*in){
        te->reinit(te->ud);
        te->started = 1;
    }
    if(te->started) {
        te->compute(te->ud, out);
    }
    else {
        *out = 0;
    }

    return SP_OK;
}


int sp_tgate_create(sp_tgate **p)
{
    *p = malloc(sizeof(sp_tgate));
    return SP_OK;
}

int sp_tgate_destroy(sp_tgate **p)
{
    free(*p);
    return SP_OK;
}

int sp_tgate_init(sp_data *sp, sp_tgate *p)
{
    p->time = 0;
    p->timer = 0;
    return SP_OK;
}

int sp_tgate_compute(sp_data *sp, sp_tgate *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = 0;
    if(*in != 0) {
        p->timer = p->time * sp->sr;
    }
    if(p->timer != 0) {
        *out = 1;
        p->timer--;
    }
    return SP_OK;
}


int sp_thresh_create(sp_thresh **p)
{
    *p = malloc(sizeof(sp_thresh));
    return SP_OK;
}

int sp_thresh_destroy(sp_thresh **p)
{
    free(*p);
    return SP_OK;
}

int sp_thresh_init(sp_data *sp, sp_thresh *p)
{
    /* Initalize variables here. */
    p->init = 1;
    p->mode = 0;
    p->prev = 0;
    p->thresh = 0;
    return SP_OK;
}

int sp_thresh_compute(sp_data *sp, sp_thresh *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->init) {
        *out = 0;
        p->prev = *in;
        p->init = 0;
        return SP_OK;
    }

    switch((int)p->mode) {
        /* input signal goes above threshold */
        case 0:
            *out = (*in > p->thresh && p->prev <= p->thresh);
            break;

        /* input signal goes below threshold */
        case 1:
            *out = (*in < p->thresh && p->prev >= p->thresh);
            break;

        /* input signal goes below or above threshold */
        case 2:
            *out = (*in < p->thresh && p->prev >= p->thresh) ||
                (*in > p->thresh && p->prev <= p->thresh);
            break;

        default:
            return SP_NOT_OK;
    }

    p->prev = *in;
    
    return SP_OK;
}


int sp_timer_create(sp_timer **p)
{
    *p = malloc(sizeof(sp_timer));
    return SP_OK;
}

int sp_timer_destroy(sp_timer **p)
{
    free(*p);
    return SP_OK;
}

int sp_timer_init(sp_data *sp, sp_timer *p)
{
    p->mode = 0;
    p->pos = 0;
    p->time = 0;
    return SP_OK;
}

int sp_timer_compute(sp_data *sp, sp_timer *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in != 0) {
        if(p->mode == 0) {
            p->pos = 0;
            p->mode = 1;
        } else if(p->mode == 1) {
            p->time = (SPFLOAT) p->pos / sp->sr;
            p->mode = 0;
        }
    }

    if(p->mode == 1) {
        p->pos++;
    }

    *out = p->time;
    return SP_OK;
}


int sp_tin_create(sp_tin **p)
{
    *p = malloc(sizeof(sp_tin));
    return SP_OK;
}

int sp_tin_destroy(sp_tin **p)
{
    free(*p);
    return SP_OK;
}

int sp_tin_init(sp_data *sp, sp_tin *p)
{
    p->fp = stdin; 
    p->val = 0;
    return SP_OK;
}

int sp_tin_compute(sp_data *sp, sp_tin *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in) {
        fread(&p->val, sizeof(SPFLOAT), 1, p->fp);
    }

    *out = p->val;
    return SP_OK;
}


int sp_tone_create(sp_tone **t)
{
    *t = malloc(sizeof(sp_tone));
    return SP_OK;
}

int sp_tone_destroy(sp_tone **t)
{
    free(*t);
    return SP_OK;
}

int sp_tone_init(sp_data *sp, sp_tone *p)
{
    p->hp = 1000;
    SPFLOAT b;
    p->tpidsr = (2.0 * M_PI) / sp->sr * 1.0;
    p->prvhp = (SPFLOAT)p->hp;
    b = 2.0 - cos((SPFLOAT)(p->prvhp * p->tpidsr));
    p->c2 = b - sqrt(b * b - 1.0);
    p->c1 = 1.0 - p->c2;
    p->yt1 = 0.0;

    return SP_OK;
}

int sp_tone_compute(sp_data *sp, sp_tone *p, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT c1 = p->c1, c2 = p->c2;
    SPFLOAT yt1 = p->yt1;

    if (p->hp != p->prvhp) {
      SPFLOAT b;
      p->prvhp = p->hp;
      b = 2.0 - cos((p->prvhp * p->tpidsr));
      p->c2 = c2 = b - sqrt(b * b - 1.0);
      p->c1 = c1 = 1.0 - c2;
    }

    yt1 = c1 * (*in) + c2 * yt1;
    *out = yt1;

    p->yt1 = yt1;
    return SP_OK;
}


int sp_trand_create(sp_trand **p)
{
    *p = malloc(sizeof(sp_trand));
    return SP_OK;
}

int sp_trand_destroy(sp_trand **p)
{
    free(*p);
    return SP_OK;
}

int sp_trand_init(sp_data *sp, sp_trand *p)
{
    p->min = 0;
    p->max = 1;
    p->val = 0;
    return SP_OK;
}

int sp_trand_compute(sp_data *sp, sp_trand *p, SPFLOAT *in, SPFLOAT *out)
{
    if(*in != 0) {
        p->val = p->min + ((SPFLOAT) sp_rand(sp) / SP_RANDMAX) * (p->max - p->min);
        *out = p->val;
    } else {
        *out = p->val;
    }
    return SP_OK;
}


int sp_tseg_create(sp_tseg **p)
{
    *p = malloc(sizeof(sp_tseg));
    return SP_OK;
}

int sp_tseg_destroy(sp_tseg **p)
{
    free(*p);
    return SP_OK;
}

int sp_tseg_init(sp_data *sp, sp_tseg *p, SPFLOAT ibeg)
{
    p->beg = ibeg;
    p->end = 1.0;
    p->val = ibeg;
    p->type = -3;
    p->slope = 1.0;
    p->dur = 1.0;
    p->tdivnsteps = 0.0;
    p->count = 0;
    p->steps = p->dur * sp->sr;
    return SP_OK;
}

int sp_tseg_compute(sp_data *sp, sp_tseg *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = p->val;
    if(*in != 0) {
        p->slope = 1.0 / (1.0 - exp(p->type));
        p->beg = p->val;
        p->count = 0;
        p->steps = p->dur * sp->sr;
        p->tdivnsteps = (SPFLOAT)p->type / (p->steps - 1);
    }
    if(p->count < p->steps) {
        *out = p->beg + (p->end - p->beg) * 
            ((1 - exp(p->count * p->tdivnsteps)) * p->slope);
        p->val = *out;
        p->count++;
    }
    return SP_OK;
}


int sp_tseq_create(sp_tseq **p)
{
    *p = malloc(sizeof(sp_tseq));
    return SP_OK;
}

int sp_tseq_destroy(sp_tseq **p)
{
    free(*p);
    return SP_OK;
}

int sp_tseq_init(sp_data *sp, sp_tseq *p, sp_ftbl *ft)
{
    p->ft = ft;
    p->pos = -1;
    p->val = 0;
    p->shuf = 0;
    return SP_OK;
}

int sp_tseq_compute(sp_data *sp, sp_tseq *p, SPFLOAT *trig, SPFLOAT *val)
{    
    if(*trig != 0){
        if(p->shuf) {
            p->pos = sp_rand(sp) % p->ft->size;
        } else {
            p->pos = (p->pos + 1) % p->ft->size;
        }
        p->val = p->ft->tbl[p->pos];
    }
    *val = p->val;
    return SP_OK;
}


int sp_vdelay_create(sp_vdelay **p)
{
    *p = malloc(sizeof(sp_vdelay));
    return SP_OK;
}

int sp_vdelay_destroy(sp_vdelay **p)
{
    sp_vdelay *pp = *p;
    sp_auxdata_free(&pp->buf);
    free(*p);
    return SP_OK;
}

int sp_vdelay_init(sp_data *sp, sp_vdelay *p, SPFLOAT maxdel)
{
    uint32_t n = (int32_t)(maxdel * sp->sr)+1;
    p->sr = sp->sr;
    p->del = maxdel * 0.5;
    p->maxdel = maxdel;
    sp_auxdata_alloc(&p->buf, n * sizeof(SPFLOAT));
    p->left = 0;
    p->feedback = 0;
    p->prev = 0;
    return SP_OK;
}

int sp_vdelay_compute(sp_data *sp, sp_vdelay *p, SPFLOAT *in, SPFLOAT *out)
{
    int32_t maxd, indx;
    *out = p->sr;
    SPFLOAT del = p->del;
    SPFLOAT b0, b1, b2, b3;
    int32_t v0, v1, v2, v3;
    SPFLOAT fv1;
    indx = p->left;
    SPFLOAT *buf = (SPFLOAT *)p->buf.ptr;

    buf[indx] = *in + p->prev*p->feedback;

    fv1 = del * (-1.0 * sp->sr);
    v1 = (int32_t)fv1;
    fv1 -= (SPFLOAT) v1;
    v1 += (int32_t)indx;
    maxd = (uint32_t) (p->maxdel * p->sr);
    /* Make sure we're inside the buffer */
    if ((v1 < 0) || (fv1 < 0.0)) {
        fv1++; v1--;
        while (v1 < 0) {
            v1 += (int32_t)maxd;
        }
    } else {
        while (v1 >= (int32_t)maxd) {
        v1 -= (int32_t)maxd;
        }
    }
    /* Find next sample for interpolation      */
    v2 = (v1 == (int32_t)(maxd - 1UL) ? 0L : v1 + 1L);

    if (maxd<4) {
        b1 = buf[v1];
        b2 = buf[v2];
        *out = b1 + fv1 * (b2 - b1);
    } else {
        v0 = (v1==0 ? maxd-1 : v1-1);
        v3 = (v2==(int32_t)maxd-1 ? 0 : v2+1);
        {
            SPFLOAT w, x, y, z;
            z = fv1 * fv1; z--;
            z *= 0.1666666667;
            y = fv1;
            y++; w = (y *= 0.5); w--;
            x = 3.0 * z; y -= x; w -= z; x -= fv1;
            b0 = buf[v0];
            b1 = buf[v1];
            b2 = buf[v2];
            b3 = buf[v3];
            *out = (w*b0 + x*b1 + y*b2 + z*b3)
            * fv1 + b1;
        }
    }
    if (++indx == maxd) indx = 0;
    p->left = indx;
    p->prev = *out;
    return SP_OK;
}

int sp_vdelay_reset(sp_data *sp, sp_vdelay *p)
{
    SPFLOAT *buf;
    uint32_t n;
    int32_t maxd;
    
    buf = (SPFLOAT *)p->buf.ptr;
    p->left = 0;
    maxd = (uint32_t) (p->maxdel * p->sr);

    for(n = 0; n < maxd; n++) buf[n] = 0;

    return SP_OK;
}




#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

#ifndef MAX
#define MAX(A,B) ((A) > (B) ? (A) : (B))
#endif

#define EPSILON 1.0e-38
#define MAX_TRANSIENTS 4


typedef struct {
    SPFLOAT  freq;
    SPFLOAT  tenseness;
    SPFLOAT  Rd;
    SPFLOAT  waveform_length;
    SPFLOAT  time_in_waveform;

    SPFLOAT  alpha;
    SPFLOAT  E0;
    SPFLOAT  epsilon;
    SPFLOAT  shift;
    SPFLOAT  delta;
    SPFLOAT  Te;
    SPFLOAT  omega;

    SPFLOAT  T;
} glottis;



typedef struct transient {
    int  position;
    SPFLOAT  time_alive;
    SPFLOAT  lifetime;
    SPFLOAT  strength;
    SPFLOAT  exponent;
    char is_free;
    unsigned int id;
    struct transient *next;
} transient;


typedef struct {
    transient pool[MAX_TRANSIENTS];
    transient *root;
    int size;
    int next_free;
} transient_pool;



typedef struct {
    int n;
     
    SPFLOAT  diameter[44];
    SPFLOAT  rest_diameter[44];
    SPFLOAT  target_diameter[44];
    SPFLOAT  new_diameter[44];
    SPFLOAT  R[44]; 
    SPFLOAT  L[44]; 
    SPFLOAT  reflection[45];
    SPFLOAT  new_reflection[45];
    SPFLOAT  junction_outL[45];
    SPFLOAT  junction_outR[45];
    SPFLOAT  A[44];

    int nose_length;


    int nose_start; 

    int tip_start;
    SPFLOAT  noseL[28];
    SPFLOAT  noseR[28];
    SPFLOAT  nose_junc_outL[29];
    SPFLOAT  nose_junc_outR[29];
    SPFLOAT  nose_reflection[29];
    SPFLOAT  nose_diameter[28];
    SPFLOAT  noseA[28];

    SPFLOAT  reflection_left;
    SPFLOAT  reflection_right;
    SPFLOAT  reflection_nose;

    SPFLOAT  new_reflection_left;
    SPFLOAT  new_reflection_right;
    SPFLOAT  new_reflection_nose;

    SPFLOAT  velum_target;

    SPFLOAT  glottal_reflection;
    SPFLOAT  lip_reflection;
    int  last_obstruction;
    SPFLOAT  fade;
    SPFLOAT  movement_speed; 
    SPFLOAT  lip_output;
    SPFLOAT  nose_output;
    SPFLOAT  block_time;

    transient_pool tpool;
    SPFLOAT  T;
} tract;


struct sp_voc {
    glottis  glot; /*The Glottis*/
    tract  tr; /*The Vocal Tract */
    SPFLOAT  buf[512];
    int counter;
};





static void glottis_setup_waveform(glottis *glot, SPFLOAT lambda)
{
    
    SPFLOAT Rd;
    SPFLOAT Ra;
    SPFLOAT Rk;
    SPFLOAT Rg;
    
    SPFLOAT Ta;
    SPFLOAT Tp;
    SPFLOAT Te;
    
    SPFLOAT epsilon;
    SPFLOAT shift;
    SPFLOAT delta;
    SPFLOAT rhs_integral;
    
    SPFLOAT lower_integral;
    SPFLOAT upper_integral;
    
    SPFLOAT omega;
    SPFLOAT s;
    SPFLOAT y;
    SPFLOAT z;
    
    SPFLOAT alpha;
    SPFLOAT E0;

    
    glot->Rd = 3 * (1 - glot->tenseness);
    glot->waveform_length = 1.0 / glot->freq;
    
    Rd = glot->Rd;
    if(Rd < 0.5) Rd = 0.5;
    if(Rd > 2.7) Rd = 2.7;

    
    Ra = -0.01 + 0.048*Rd;
    Rk = 0.224 + 0.118*Rd;
    Rg = (Rk/4)*(0.5 + 1.2*Rk)/(0.11*Rd-Ra*(0.5+1.2*Rk));

    
    Ta = Ra;
    Tp = (SPFLOAT)1.0 / (2*Rg);
    Te = Tp + Tp*Rk;


    
    epsilon = (SPFLOAT)1.0 / Ta;
    shift = exp(-epsilon * (1 - Te));
    delta = 1 - shift;


    
    rhs_integral = (SPFLOAT)(1.0/epsilon) * (shift-1) + (1-Te)*shift;
    rhs_integral = rhs_integral / delta;
    lower_integral = - (Te - Tp) / 2 + rhs_integral;
    upper_integral = -lower_integral;

    
    omega = M_PI / Tp;
    s = sin(omega * Te);
    
    y = -M_PI * s * upper_integral / (Tp*2);
    z = log(y);
    alpha = z / (Tp/2 - Te);
    E0 = -1 / (s * exp(alpha*Te));


    
    glot->alpha = alpha;
    glot->E0 = E0;
    glot->epsilon = epsilon;
    glot->shift = shift;
    glot->delta = delta;
    glot->Te = Te;
    glot->omega = omega;
}


static void glottis_init(glottis *glot, SPFLOAT sr)
{
    glot->freq = 140; /* 140Hz frequency by default */
    glot->tenseness = 0.6; /* value between 0 and 1 */
    glot->T = 1.0/sr; /* big T */
    glot->time_in_waveform = 0;
    glottis_setup_waveform(glot, 0);
}


static SPFLOAT glottis_compute(sp_data *sp, glottis *glot, SPFLOAT lambda)
{
    SPFLOAT out;
    SPFLOAT aspiration;
    SPFLOAT noise;
    SPFLOAT t;
    SPFLOAT intensity;

    out = 0;
    intensity = 1.0;
    glot->time_in_waveform += glot->T;

    if(glot->time_in_waveform > glot->waveform_length) {
        glot->time_in_waveform -= glot->waveform_length;
        glottis_setup_waveform(glot, lambda);

    }

    t = (glot->time_in_waveform / glot->waveform_length);

    if(t > glot->Te) {
        out = (-exp(-glot->epsilon * (t-glot->Te)) + glot->shift) / glot->delta;
    } else {
        out = glot->E0 * exp(glot->alpha * t) * sin(glot->omega * t);
    }




    noise = 2.0 * ((SPFLOAT) sp_rand(sp) / SP_RANDMAX) - 1;


    aspiration = intensity * (1 - sqrt(glot->tenseness)) * 0.3 * noise;

    aspiration *= 0.2;

    out += aspiration;

    return out;
}




static void tract_calculate_reflections(tract *tr)
{
    int i;
    SPFLOAT  sum; 

    for(i = 0; i < tr->n; i++) {
        tr->A[i] = tr->diameter[i] * tr->diameter[i];
        /* Calculate area from diameter squared*/
    }

    for(i = 1; i < tr->n; i++) {
        tr->reflection[i] = tr->new_reflection[i];
        if(tr->A[i] == 0) {
            tr->new_reflection[i] = 0.999; /* to prevent bad behavior if 0 */
        } else {
            tr->new_reflection[i] =
                (tr->A[i - 1] - tr->A[i]) / (tr->A[i - 1] + tr->A[i]);
        }
    }

    tr->reflection_left = tr->new_reflection_left;
    tr->reflection_right = tr->new_reflection_right;
    tr->reflection_nose = tr->new_reflection_nose;

    sum = tr->A[tr->nose_start] + tr->A[tr->nose_start + 1] + tr->noseA[0];
    tr->new_reflection_left = (SPFLOAT)(2 * tr->A[tr->nose_start] - sum) / sum;
    tr->new_reflection_right = (SPFLOAT)(2 * tr->A[tr->nose_start + 1] - sum) / sum;
    tr->new_reflection_nose = (SPFLOAT)(2 * tr->noseA[0] - sum) / sum;
}


static void tract_calculate_nose_reflections(tract *tr)
{
    int i;

    for(i = 0; i < tr->nose_length; i++) {
        tr->noseA[i] = tr->nose_diameter[i] * tr->nose_diameter[i];
    }

    for(i = 1; i < tr->nose_length; i++) {
        tr->nose_reflection[i] = (tr->noseA[i - 1] - tr->noseA[i]) /
            (tr->noseA[i-1] + tr->noseA[i]);
    }
}



static int append_transient(transient_pool *pool, int position)
{
    int i;
    int free_id;
    transient *t;

    free_id = pool->next_free;
    if(pool->size == MAX_TRANSIENTS) return 0;

    if(free_id == -1) {
        for(i = 0; i < MAX_TRANSIENTS; i++) {
            if(pool->pool[i].is_free) {
                free_id = i;
                break;
            }
        }
    }

    if(free_id == -1) return 0;

    t = &pool->pool[free_id];
    t->next = pool->root;
    pool->root = t;
    pool->size++;
    t->is_free = 0;
    t->time_alive = 0;
    t->lifetime = 0.2;
    t->strength = 0.3;
    t->exponent = 200;
    t->position = position;
    pool->next_free = -1;
    return 0;
}


static void remove_transient(transient_pool *pool, unsigned int id)
{
    int i;
    transient *n;

    pool->next_free = id;
    n = pool->root;
    if(id == n->id) {
        pool->root = n->next;
        pool->size--;
        return;
    }

    for(i = 0; i < pool->size; i++) {
        if(n->next->id == id) {
            pool->size--;
            n->next->is_free = 1;
            n->next = n->next->next;
            break;
        }
        n = n->next;
    }
}



static SPFLOAT move_towards(SPFLOAT current, SPFLOAT target,
        SPFLOAT amt_up, SPFLOAT amt_down)
{
    SPFLOAT tmp;
    if(current < target) {
        tmp = current + amt_up;
        return MIN(tmp, target);
    } else {
        tmp = current - amt_down;
        return MAX(tmp, target);
    }
    return 0.0;
}

static void tract_reshape(tract *tr)
{
    SPFLOAT amount;
    SPFLOAT slow_return;
    SPFLOAT diameter;
    SPFLOAT target_diameter;
    int i;
    int current_obstruction;

    current_obstruction = -1;
    amount = tr->block_time * tr->movement_speed;

    for(i = 0; i < tr->n; i++) {
        slow_return = 0;
        diameter = tr->diameter[i];
        target_diameter = tr->target_diameter[i];

        if(diameter < 0.001) current_obstruction = i;

        if(i < tr->nose_start) slow_return = 0.6;
        else if(i >= tr->tip_start) slow_return = 1.0;
        else {
            slow_return =
                0.6+0.4*(i - tr->nose_start)/(tr->tip_start - tr->nose_start);
        }

        tr->diameter[i] = move_towards(diameter, target_diameter,
                slow_return * amount, 2 * amount);

    }

    if(tr->last_obstruction > -1 && current_obstruction == -1 &&
            tr->noseA[0] < 0.05) {
        append_transient(&tr->tpool, tr->last_obstruction);
    }
    tr->last_obstruction = current_obstruction;

    tr->nose_diameter[0] = move_towards(tr->nose_diameter[0], tr->velum_target,
            amount * 0.25, amount * 0.1);
    tr->noseA[0] = tr->nose_diameter[0] * tr->nose_diameter[0];
}


static void tract_init(sp_data *sp, tract *tr)
{
    int i;
    SPFLOAT diameter, d; /* needed to set up diameter arrays */
    
    tr->n = 44;
    tr->nose_length = 28;
    tr->nose_start = 17;
    
    tr->reflection_left = 0.0;
    tr->reflection_right = 0.0;
    tr->reflection_nose = 0.0;
    tr->new_reflection_left = 0.0;
    tr->new_reflection_right= 0.0;
    tr->new_reflection_nose = 0.0;
    tr->velum_target = 0.01;
    tr->glottal_reflection = 0.75;
    tr->lip_reflection = -0.85;
    tr->last_obstruction = -1;
    tr->movement_speed = 15;
    tr->lip_output = 0;
    tr->nose_output = 0;
    tr->tip_start = 32;

    
    memset(tr->diameter, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->rest_diameter, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->target_diameter, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->new_diameter, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->L, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->R, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->reflection, 0, (tr->n + 1) * sizeof(SPFLOAT));
    memset(tr->new_reflection, 0, (tr->n + 1) * sizeof(SPFLOAT));
    memset(tr->junction_outL, 0, (tr->n + 1) * sizeof(SPFLOAT));
    memset(tr->junction_outR, 0, (tr->n + 1) * sizeof(SPFLOAT));
    memset(tr->A, 0, tr->n * sizeof(SPFLOAT));
    memset(tr->noseL, 0, tr->nose_length * sizeof(SPFLOAT));
    memset(tr->noseR, 0, tr->nose_length * sizeof(SPFLOAT));
    memset(tr->nose_junc_outL, 0, (tr->nose_length + 1) * sizeof(SPFLOAT));
    memset(tr->nose_junc_outR, 0, (tr->nose_length + 1) * sizeof(SPFLOAT));
    memset(tr->nose_diameter, 0, tr->nose_length * sizeof(SPFLOAT));
    memset(tr->noseA, 0, tr->nose_length * sizeof(SPFLOAT));

    
    for(i = 0; i < tr->n; i++) {
        diameter = 0;
        if(i < 7 * (SPFLOAT)tr->n / 44 - 0.5) {
            diameter = 0.6;
        } else if( i < 12 * (SPFLOAT)tr->n / 44) {
            diameter = 1.1;
        } else {
            diameter = 1.5;
        }
    
        tr->diameter[i] =
            tr->rest_diameter[i] =
            tr->target_diameter[i] =
            tr->new_diameter[i] = diameter;
    
    }

    
        for(i = 0; i < tr->nose_length; i++) {
            d = 2 * ((SPFLOAT)i / tr->nose_length);
            if(d < 1) {
                diameter = 0.4 + 1.6 * d;
            } else {
                diameter = 0.5 + 1.5*(2-d);
            }
            diameter = MIN(diameter, 1.9);
            tr->nose_diameter[i] = diameter;
        }


    tract_calculate_reflections(tr);
    tract_calculate_nose_reflections(tr);
    tr->nose_diameter[0] = tr->velum_target;

    tr->block_time = 512.0 / (SPFLOAT)sp->sr;
    tr->T = 1.0 / (SPFLOAT)sp->sr;
    
    tr->tpool.size = 0;
    tr->tpool.next_free = 0;
    for(i = 0; i < MAX_TRANSIENTS; i++) {
        tr->tpool.pool[i].is_free = 1;
        tr->tpool.pool[i].id = i;
        tr->tpool.pool[i].position = 0;
        tr->tpool.pool[i].time_alive = 0;
        tr->tpool.pool[i].strength = 0;
        tr->tpool.pool[i].exponent = 0;
    }
}


static void tract_compute(sp_data *sp, tract *tr,
    SPFLOAT  in,
    SPFLOAT  lambda)
{
     SPFLOAT  r, w;
    int i;
    SPFLOAT  amp;
    int current_size;
    transient_pool *pool;
    transient *n;

    
    
        pool = &tr->tpool;
        current_size = pool->size;
        n = pool->root;
        for(i = 0; i < current_size; i++) {
            amp = n->strength * pow(2, -1.0 * n->exponent * n->time_alive);
            tr->L[n->position] += amp * 0.5;
            tr->R[n->position] += amp * 0.5;
            n->time_alive += tr->T * 0.5;
            if(n->time_alive > n->lifetime) {
                 remove_transient(pool, n->id);
            }
            n = n->next;
        }

    
    tr->junction_outR[0] = tr->L[0] * tr->glottal_reflection + in;
    tr->junction_outL[tr->n] = tr->R[tr->n - 1] * tr->lip_reflection;
    
    for(i = 1; i < tr->n; i++) {
        r = tr->reflection[i] * (1 - lambda) + tr->new_reflection[i] * lambda;
        w = r * (tr->R[i - 1] + tr->L[i]);
        tr->junction_outR[i] = tr->R[i - 1] - w;
        tr->junction_outL[i] = tr->L[i] + w;
    }

    
    i = tr->nose_start;
    r = tr->new_reflection_left * (1-lambda) + tr->reflection_left*lambda;
    tr->junction_outL[i] = r*tr->R[i-1] + (1+r)*(tr->noseL[0]+tr->L[i]);
    r = tr->new_reflection_right * (1 - lambda) + tr->reflection_right * lambda;
    tr->junction_outR[i] = r*tr->L[i] + (1+r)*(tr->R[i-1]+tr->noseL[0]);
    r = tr->new_reflection_nose * (1 - lambda) + tr->reflection_nose * lambda;
    tr->nose_junc_outR[0] = r * tr->noseL[0]+(1+r)*(tr->L[i]+tr->R[i-1]);

    
    for(i = 0; i < tr->n; i++) {
        tr->R[i] = tr->junction_outR[i]*0.999;
        tr->L[i] = tr->junction_outL[i + 1]*0.999;
    }
    tr->lip_output = tr->R[tr->n - 1];

    
    tr->nose_junc_outL[tr->nose_length] =
        tr->noseR[tr->nose_length-1] * tr->lip_reflection;
    
    for(i = 1; i < tr->nose_length; i++) {
        w = tr->nose_reflection[i] * (tr->noseR[i-1] + tr->noseL[i]);
        tr->nose_junc_outR[i] = tr->noseR[i - 1] - w;
        tr->nose_junc_outL[i] = tr->noseL[i] + w;
    }

    
    for(i = 0; i < tr->nose_length; i++) {
        tr->noseR[i] = tr->nose_junc_outR[i];
        tr->noseL[i] = tr->nose_junc_outL[i + 1];
    }
    tr->nose_output = tr->noseR[tr->nose_length - 1];

}




int sp_voc_create(sp_voc **voc)
{
    *voc = malloc(sizeof(sp_voc));
    return SP_OK;
}


int sp_voc_destroy(sp_voc **voc)
{
    free(*voc);
    return SP_OK;
}


int sp_voc_init(sp_data *sp, sp_voc *voc)
{
    glottis_init(&voc->glot, sp->sr); /* initialize glottis */
    tract_init(sp, &voc->tr); /* initialize vocal tract */
    voc->counter = 0;
    return SP_OK;
}


int sp_voc_compute(sp_data *sp, sp_voc *voc, SPFLOAT *out)
{
    SPFLOAT vocal_output, glot;
    SPFLOAT lambda1, lambda2;
    int i;

    

    if(voc->counter == 0) {
        tract_reshape(&voc->tr);
        tract_calculate_reflections(&voc->tr);
        for(i = 0; i < 512; i++) {
            vocal_output = 0;
            lambda1 = (SPFLOAT) i / 512;
            lambda2 = (SPFLOAT) (i + 0.5) / 512;
            glot = glottis_compute(sp, &voc->glot, lambda1);

            tract_compute(sp, &voc->tr, glot, lambda1);
            vocal_output += voc->tr.lip_output + voc->tr.nose_output;

            tract_compute(sp, &voc->tr, glot, lambda2);
            vocal_output += voc->tr.lip_output + voc->tr.nose_output;
            voc->buf[i] = vocal_output * 0.125;
        }
    }


    *out = voc->buf[voc->counter];
    voc->counter = (voc->counter + 1) % 512;
    return SP_OK;
}


int sp_voc_tract_compute(sp_data *sp, sp_voc *voc, SPFLOAT *in, SPFLOAT *out)
{
    SPFLOAT vocal_output;
    SPFLOAT lambda1, lambda2;

    if(voc->counter == 0) {
        tract_reshape(&voc->tr);
        tract_calculate_reflections(&voc->tr);
    }

    vocal_output = 0;
    lambda1 = (SPFLOAT) voc->counter / 512;
    lambda2 = (SPFLOAT) (voc->counter + 0.5) / 512;

    tract_compute(sp, &voc->tr, *in, lambda1);
    vocal_output += voc->tr.lip_output + voc->tr.nose_output;
    tract_compute(sp, &voc->tr, *in, lambda2);
    vocal_output += voc->tr.lip_output + voc->tr.nose_output;


    *out = vocal_output * 0.125;
    voc->counter = (voc->counter + 1) % 512;
    return SP_OK;
}


void sp_voc_set_frequency(sp_voc *voc, SPFLOAT freq)
{
    voc->glot.freq = freq;
}


SPFLOAT * sp_voc_get_frequency_ptr(sp_voc *voc)
{
    return &voc->glot.freq;
}


SPFLOAT* sp_voc_get_tract_diameters(sp_voc *voc)
{
    return voc->tr.target_diameter;
}


SPFLOAT* sp_voc_get_current_tract_diameters(sp_voc *voc)
{
    return voc->tr.diameter;
}


int sp_voc_get_tract_size(sp_voc *voc)
{
    return voc->tr.n;
}


SPFLOAT* sp_voc_get_nose_diameters(sp_voc *voc)
{
    return voc->tr.nose_diameter;
}


int sp_voc_get_nose_size(sp_voc *voc)
{
    return voc->tr.nose_length;
}


void sp_voc_set_diameters(sp_voc *voc, 
    int blade_start, 
    int lip_start, 
    int tip_start, 
    SPFLOAT tongue_index,
    SPFLOAT tongue_diameter, 
    SPFLOAT *diameters) {

    int i;
    SPFLOAT t;
    SPFLOAT fixed_tongue_diameter;
    SPFLOAT curve;
    int grid_offset = 0;

    for(i = blade_start; i < lip_start; i++) {
        t = 1.1 * M_PI *
            (SPFLOAT)(tongue_index - i)/(tip_start - blade_start);
        fixed_tongue_diameter = 2+(tongue_diameter-2)/1.5;
        curve = (1.5 - fixed_tongue_diameter + grid_offset) * cos(t);
        if(i == blade_start - 2 || i == lip_start - 1) curve *= 0.8;
        if(i == blade_start || i == lip_start - 2) curve *= 0.94;
        diameters[i] = 1.5 - curve;
    }
}


void sp_voc_set_tongue_shape(sp_voc *voc,
    SPFLOAT tongue_index, 
    SPFLOAT tongue_diameter) {
    SPFLOAT *diameters;
    diameters = sp_voc_get_tract_diameters(voc);
    sp_voc_set_diameters(voc, 10, 39, 32,
            tongue_index, tongue_diameter, diameters);
}


int sp_voc_get_counter(sp_voc *voc)
{
    return voc->counter;
}


void sp_voc_set_tenseness(sp_voc *voc, SPFLOAT tenseness)
{
    voc->glot.tenseness = tenseness;
}


SPFLOAT * sp_voc_get_tenseness_ptr(sp_voc *voc)
{
    return &voc->glot.tenseness;
}


void sp_voc_set_velum(sp_voc *voc, SPFLOAT velum)
{
    voc->tr.velum_target = velum;
}


SPFLOAT *sp_voc_get_velum_ptr(sp_voc *voc)
{
    return &voc->tr.velum_target;
}



#define CHECK_NULL_FILE(fp) if(fp == NULL) return SP_NOT_OK

int spa_open(sp_data *sp, sp_audio *spa, const char *name, int mode)
{
    spa->mode = SPA_NULL;
    spa_header *header = &spa->header;
    spa->offset = sizeof(spa_header);
    if(mode == SPA_READ) {
        spa->fp = fopen(name, "rb");
        CHECK_NULL_FILE(spa->fp);
        fread(header, spa->offset, 1, spa->fp);
    } else if(mode == SPA_WRITE) {
        spa->fp = fopen(name, "wb");
        CHECK_NULL_FILE(spa->fp);
        header->magic = 100;
        header->nchan = sp->nchan;
        header->len = sp->len;
        header->sr = sp->sr;
        fwrite(header, spa->offset, 1, spa->fp);
    } else {
        return SP_NOT_OK;
    }

    spa->mode = mode;

    return SP_OK;
}

size_t spa_write_buf(sp_data *sp, sp_audio *spa, SPFLOAT *buf, uint32_t size)
{
    if(spa->mode != SPA_WRITE) {
        return 0;
    }
    return fwrite(buf, sizeof(SPFLOAT), size, spa->fp);
}

size_t spa_read_buf(sp_data *sp, sp_audio *spa, SPFLOAT *buf, uint32_t size)
{
    if(spa->mode != SPA_READ) {
        return 0;
    }
    return fread(buf, sizeof(SPFLOAT), size, spa->fp);
}

int spa_close(sp_audio *spa)
{
    if(spa->fp != NULL) fclose(spa->fp);
    return SP_OK;
}

int sp_process_spa(sp_data *sp, void *ud, void (*callback)(sp_data *, void *))
{
    sp_audio spa;
    if(spa_open(sp, &spa, sp->filename, SPA_WRITE) == SP_NOT_OK) {
        fprintf(stderr, "Error: could not open file %s.\n", sp->filename);    
    }
    while(sp->len > 0) {
        callback(sp, ud);
        spa_write_buf(sp, &spa, sp->out, sp->nchan);
        sp->len--;
        sp->pos++;
    }
    spa_close(&spa);
    return SP_OK;
}

int sp_ftbl_loadspa(sp_data *sp, sp_ftbl **ft, const char *filename)
{
    *ft = malloc(sizeof(sp_ftbl));
    sp_ftbl *ftp = *ft;

    sp_audio spa;

    spa_open(sp, &spa, filename, SPA_READ);

    size_t size = spa.header.len;

    ftp->tbl = malloc(sizeof(SPFLOAT) * (size + 1));
    sp_ftbl_init(sp, ftp, size);

    spa_read_buf(sp, &spa, ftp->tbl, ftp->size);
    spa_close(&spa);
    return SP_OK;
}











typedef struct {
	
	float fRec0[3];
	float fRec3[3];
	float fRec4[3];
	float fRec7[3];
	float fRec8[3];
	float fRec11[3];
	float fRec12[3];
	float fRec15[3];
	float fRec16[3];
	float fRec19[3];
	float fRec20[3];
	float fRec23[3];
	float fRec24[3];
	float fRec27[3];
	float fRec28[3];
	float fRec31[3];
	float fRec32[3];
	float fRec35[3];
	float fRec36[3];
	float fRec39[3];
	float fRec40[3];
	float fRec43[3];
	float fRec44[3];
	float fRec47[3];
	float fRec48[3];
	float fRec51[3];
	float fRec52[3];
	float fRec55[3];
	float fRec56[3];
	float fRec59[3];
	float fRec60[3];
	float fRec63[3];
	float fRec2[2];
	float fRec1[2];
	float fRec6[2];
	float fRec5[2];
	float fRec10[2];
	float fRec9[2];
	float fRec14[2];
	float fRec13[2];
	float fRec18[2];
	float fRec17[2];
	float fRec22[2];
	float fRec21[2];
	float fRec26[2];
	float fRec25[2];
	float fRec30[2];
	float fRec29[2];
	float fRec34[2];
	float fRec33[2];
	float fRec38[2];
	float fRec37[2];
	float fRec42[2];
	float fRec41[2];
	float fRec46[2];
	float fRec45[2];
	float fRec50[2];
	float fRec49[2];
	float fRec54[2];
	float fRec53[2];
	float fRec58[2];
	float fRec57[2];
	float fRec62[2];
	float fRec61[2];
	int fSamplingFreq;
	int iConst0;
	float fConst1;
	float fConst2;
	FAUSTFLOAT fHslider0;
	float fConst3;
	float fConst4;
	float fConst5;
	FAUSTFLOAT fHslider1;
	FAUSTFLOAT fHslider2;
	float fConst6;
	float fConst7;
	float fConst8;
	float fConst9;
	float fConst10;
	float fConst11;
	float fConst12;
	float fConst13;
	float fConst14;
	float fConst15;
	float fConst16;
	float fConst17;
	float fConst18;
	float fConst19;
	float fConst20;
	float fConst21;
	float fConst22;
	float fConst23;
	float fConst24;
	float fConst25;
	float fConst26;
	float fConst27;
	float fConst28;
	float fConst29;
	float fConst30;
	float fConst31;
	float fConst32;
	float fConst33;
	float fConst34;
	float fConst35;
	float fConst36;
	float fConst37;
	float fConst38;
	float fConst39;
	float fConst40;
	float fConst41;
	float fConst42;
	float fConst43;
	float fConst44;
	float fConst45;
	float fConst46;
	float fConst47;
	float fConst48;
	float fConst49;
	float fConst50;
	float fConst51;
	float fConst52;
	float fConst53;
	float fConst54;
	float fConst55;
	float fConst56;
	float fConst57;
	float fConst58;
	float fConst59;
	float fConst60;
	float fConst61;
	float fConst62;
	float fConst63;
	float fConst64;
	float fConst65;
	
} vocoder;

static vocoder* newvocoder() { 
	vocoder* dsp = (vocoder*)malloc(sizeof(vocoder));
	return dsp;
}

static void deletevocoder(vocoder* dsp) { 
	free(dsp);
}


static void instanceInitvocoder(vocoder* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fConst1 = tan((115.99f / (float)dsp->iConst0));
	dsp->fConst2 = (1.f / dsp->fConst1);
	dsp->fHslider0 = (FAUSTFLOAT)0.5;
	dsp->fConst3 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst1))));
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 3); i0 = (i0 + 1)) {
			dsp->fRec0[i0] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 3); i1 = (i1 + 1)) {
			dsp->fRec3[i1] = 0.f;
			
		}
		
	}
	dsp->fConst4 = (0.f - dsp->fConst2);
	dsp->fConst5 = (1.f / (float)dsp->iConst0);
	dsp->fHslider1 = (FAUSTFLOAT)0.01;
	dsp->fHslider2 = (FAUSTFLOAT)0.01;
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 2); i2 = (i2 + 1)) {
			dsp->fRec2[i2] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i3;
		for (i3 = 0; (i3 < 2); i3 = (i3 + 1)) {
			dsp->fRec1[i3] = 0.f;
			
		}
		
	}
	dsp->fConst6 = tan((171.297f / (float)dsp->iConst0));
	dsp->fConst7 = (1.f / dsp->fConst6);
	dsp->fConst8 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst6))));
	/* C99 loop */
	{
		int i4;
		for (i4 = 0; (i4 < 3); i4 = (i4 + 1)) {
			dsp->fRec4[i4] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i5;
		for (i5 = 0; (i5 < 3); i5 = (i5 + 1)) {
			dsp->fRec7[i5] = 0.f;
			
		}
		
	}
	dsp->fConst9 = (0.f - dsp->fConst7);
	/* C99 loop */
	{
		int i6;
		for (i6 = 0; (i6 < 2); i6 = (i6 + 1)) {
			dsp->fRec6[i6] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i7;
		for (i7 = 0; (i7 < 2); i7 = (i7 + 1)) {
			dsp->fRec5[i7] = 0.f;
			
		}
		
	}
	dsp->fConst10 = tan((252.975f / (float)dsp->iConst0));
	dsp->fConst11 = (1.f / dsp->fConst10);
	dsp->fConst12 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst10))));
	/* C99 loop */
	{
		int i8;
		for (i8 = 0; (i8 < 3); i8 = (i8 + 1)) {
			dsp->fRec8[i8] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i9;
		for (i9 = 0; (i9 < 3); i9 = (i9 + 1)) {
			dsp->fRec11[i9] = 0.f;
			
		}
		
	}
	dsp->fConst13 = (0.f - dsp->fConst11);
	/* C99 loop */
	{
		int i10;
		for (i10 = 0; (i10 < 2); i10 = (i10 + 1)) {
			dsp->fRec10[i10] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i11;
		for (i11 = 0; (i11 < 2); i11 = (i11 + 1)) {
			dsp->fRec9[i11] = 0.f;
			
		}
		
	}
	dsp->fConst14 = tan((373.6f / (float)dsp->iConst0));
	dsp->fConst15 = (1.f / dsp->fConst14);
	dsp->fConst16 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst14))));
	/* C99 loop */
	{
		int i12;
		for (i12 = 0; (i12 < 3); i12 = (i12 + 1)) {
			dsp->fRec12[i12] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i13;
		for (i13 = 0; (i13 < 3); i13 = (i13 + 1)) {
			dsp->fRec15[i13] = 0.f;
			
		}
		
	}
	dsp->fConst17 = (0.f - dsp->fConst15);
	/* C99 loop */
	{
		int i14;
		for (i14 = 0; (i14 < 2); i14 = (i14 + 1)) {
			dsp->fRec14[i14] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i15;
		for (i15 = 0; (i15 < 2); i15 = (i15 + 1)) {
			dsp->fRec13[i15] = 0.f;
			
		}
		
	}
	dsp->fConst18 = tan((551.743f / (float)dsp->iConst0));
	dsp->fConst19 = (1.f / dsp->fConst18);
	dsp->fConst20 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst18))));
	/* C99 loop */
	{
		int i16;
		for (i16 = 0; (i16 < 3); i16 = (i16 + 1)) {
			dsp->fRec16[i16] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i17;
		for (i17 = 0; (i17 < 3); i17 = (i17 + 1)) {
			dsp->fRec19[i17] = 0.f;
			
		}
		
	}
	dsp->fConst21 = (0.f - dsp->fConst19);
	/* C99 loop */
	{
		int i18;
		for (i18 = 0; (i18 < 2); i18 = (i18 + 1)) {
			dsp->fRec18[i18] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i19;
		for (i19 = 0; (i19 < 2); i19 = (i19 + 1)) {
			dsp->fRec17[i19] = 0.f;
			
		}
		
	}
	dsp->fConst22 = tan((814.828f / (float)dsp->iConst0));
	dsp->fConst23 = (1.f / dsp->fConst22);
	dsp->fConst24 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst22))));
	/* C99 loop */
	{
		int i20;
		for (i20 = 0; (i20 < 3); i20 = (i20 + 1)) {
			dsp->fRec20[i20] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i21;
		for (i21 = 0; (i21 < 3); i21 = (i21 + 1)) {
			dsp->fRec23[i21] = 0.f;
			
		}
		
	}
	dsp->fConst25 = (0.f - dsp->fConst23);
	/* C99 loop */
	{
		int i22;
		for (i22 = 0; (i22 < 2); i22 = (i22 + 1)) {
			dsp->fRec22[i22] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i23;
		for (i23 = 0; (i23 < 2); i23 = (i23 + 1)) {
			dsp->fRec21[i23] = 0.f;
			
		}
		
	}
	dsp->fConst26 = tan((1203.36f / (float)dsp->iConst0));
	dsp->fConst27 = (1.f / dsp->fConst26);
	dsp->fConst28 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst26))));
	/* C99 loop */
	{
		int i24;
		for (i24 = 0; (i24 < 3); i24 = (i24 + 1)) {
			dsp->fRec24[i24] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i25;
		for (i25 = 0; (i25 < 3); i25 = (i25 + 1)) {
			dsp->fRec27[i25] = 0.f;
			
		}
		
	}
	dsp->fConst29 = (0.f - dsp->fConst27);
	/* C99 loop */
	{
		int i26;
		for (i26 = 0; (i26 < 2); i26 = (i26 + 1)) {
			dsp->fRec26[i26] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i27;
		for (i27 = 0; (i27 < 2); i27 = (i27 + 1)) {
			dsp->fRec25[i27] = 0.f;
			
		}
		
	}
	dsp->fConst30 = tan((1777.15f / (float)dsp->iConst0));
	dsp->fConst31 = (1.f / dsp->fConst30);
	dsp->fConst32 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst30))));
	/* C99 loop */
	{
		int i28;
		for (i28 = 0; (i28 < 3); i28 = (i28 + 1)) {
			dsp->fRec28[i28] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i29;
		for (i29 = 0; (i29 < 3); i29 = (i29 + 1)) {
			dsp->fRec31[i29] = 0.f;
			
		}
		
	}
	dsp->fConst33 = (0.f - dsp->fConst31);
	/* C99 loop */
	{
		int i30;
		for (i30 = 0; (i30 < 2); i30 = (i30 + 1)) {
			dsp->fRec30[i30] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i31;
		for (i31 = 0; (i31 < 2); i31 = (i31 + 1)) {
			dsp->fRec29[i31] = 0.f;
			
		}
		
	}
	dsp->fConst34 = tan((2624.55f / (float)dsp->iConst0));
	dsp->fConst35 = (1.f / dsp->fConst34);
	dsp->fConst36 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst34))));
	/* C99 loop */
	{
		int i32;
		for (i32 = 0; (i32 < 3); i32 = (i32 + 1)) {
			dsp->fRec32[i32] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i33;
		for (i33 = 0; (i33 < 3); i33 = (i33 + 1)) {
			dsp->fRec35[i33] = 0.f;
			
		}
		
	}
	dsp->fConst37 = (0.f - dsp->fConst35);
	/* C99 loop */
	{
		int i34;
		for (i34 = 0; (i34 < 2); i34 = (i34 + 1)) {
			dsp->fRec34[i34] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i35;
		for (i35 = 0; (i35 < 2); i35 = (i35 + 1)) {
			dsp->fRec33[i35] = 0.f;
			
		}
		
	}
	dsp->fConst38 = tan((3876.f / (float)dsp->iConst0));
	dsp->fConst39 = (1.f / dsp->fConst38);
	dsp->fConst40 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst38))));
	/* C99 loop */
	{
		int i36;
		for (i36 = 0; (i36 < 3); i36 = (i36 + 1)) {
			dsp->fRec36[i36] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i37;
		for (i37 = 0; (i37 < 3); i37 = (i37 + 1)) {
			dsp->fRec39[i37] = 0.f;
			
		}
		
	}
	dsp->fConst41 = (0.f - dsp->fConst39);
	/* C99 loop */
	{
		int i38;
		for (i38 = 0; (i38 < 2); i38 = (i38 + 1)) {
			dsp->fRec38[i38] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i39;
		for (i39 = 0; (i39 < 2); i39 = (i39 + 1)) {
			dsp->fRec37[i39] = 0.f;
			
		}
		
	}
	dsp->fConst42 = tan((5724.18f / (float)dsp->iConst0));
	dsp->fConst43 = (1.f / dsp->fConst42);
	dsp->fConst44 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst42))));
	/* C99 loop */
	{
		int i40;
		for (i40 = 0; (i40 < 3); i40 = (i40 + 1)) {
			dsp->fRec40[i40] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i41;
		for (i41 = 0; (i41 < 3); i41 = (i41 + 1)) {
			dsp->fRec43[i41] = 0.f;
			
		}
		
	}
	dsp->fConst45 = (0.f - dsp->fConst43);
	/* C99 loop */
	{
		int i42;
		for (i42 = 0; (i42 < 2); i42 = (i42 + 1)) {
			dsp->fRec42[i42] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i43;
		for (i43 = 0; (i43 < 2); i43 = (i43 + 1)) {
			dsp->fRec41[i43] = 0.f;
			
		}
		
	}
	dsp->fConst46 = tan((8453.61f / (float)dsp->iConst0));
	dsp->fConst47 = (1.f / dsp->fConst46);
	dsp->fConst48 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst46))));
	/* C99 loop */
	{
		int i44;
		for (i44 = 0; (i44 < 3); i44 = (i44 + 1)) {
			dsp->fRec44[i44] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i45;
		for (i45 = 0; (i45 < 3); i45 = (i45 + 1)) {
			dsp->fRec47[i45] = 0.f;
			
		}
		
	}
	dsp->fConst49 = (0.f - dsp->fConst47);
	/* C99 loop */
	{
		int i46;
		for (i46 = 0; (i46 < 2); i46 = (i46 + 1)) {
			dsp->fRec46[i46] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i47;
		for (i47 = 0; (i47 < 2); i47 = (i47 + 1)) {
			dsp->fRec45[i47] = 0.f;
			
		}
		
	}
	dsp->fConst50 = tan((12484.5f / (float)dsp->iConst0));
	dsp->fConst51 = (1.f / dsp->fConst50);
	dsp->fConst52 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst50))));
	/* C99 loop */
	{
		int i48;
		for (i48 = 0; (i48 < 3); i48 = (i48 + 1)) {
			dsp->fRec48[i48] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i49;
		for (i49 = 0; (i49 < 3); i49 = (i49 + 1)) {
			dsp->fRec51[i49] = 0.f;
			
		}
		
	}
	dsp->fConst53 = (0.f - dsp->fConst51);
	/* C99 loop */
	{
		int i50;
		for (i50 = 0; (i50 < 2); i50 = (i50 + 1)) {
			dsp->fRec50[i50] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i51;
		for (i51 = 0; (i51 < 2); i51 = (i51 + 1)) {
			dsp->fRec49[i51] = 0.f;
			
		}
		
	}
	dsp->fConst54 = tan((18437.5f / (float)dsp->iConst0));
	dsp->fConst55 = (1.f / dsp->fConst54);
	dsp->fConst56 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst54))));
	/* C99 loop */
	{
		int i52;
		for (i52 = 0; (i52 < 3); i52 = (i52 + 1)) {
			dsp->fRec52[i52] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i53;
		for (i53 = 0; (i53 < 3); i53 = (i53 + 1)) {
			dsp->fRec55[i53] = 0.f;
			
		}
		
	}
	dsp->fConst57 = (0.f - dsp->fConst55);
	/* C99 loop */
	{
		int i54;
		for (i54 = 0; (i54 < 2); i54 = (i54 + 1)) {
			dsp->fRec54[i54] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i55;
		for (i55 = 0; (i55 < 2); i55 = (i55 + 1)) {
			dsp->fRec53[i55] = 0.f;
			
		}
		
	}
	dsp->fConst58 = tan((27228.9f / (float)dsp->iConst0));
	dsp->fConst59 = (1.f / dsp->fConst58);
	dsp->fConst60 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst58))));
	/* C99 loop */
	{
		int i56;
		for (i56 = 0; (i56 < 3); i56 = (i56 + 1)) {
			dsp->fRec56[i56] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i57;
		for (i57 = 0; (i57 < 3); i57 = (i57 + 1)) {
			dsp->fRec59[i57] = 0.f;
			
		}
		
	}
	dsp->fConst61 = (0.f - dsp->fConst59);
	/* C99 loop */
	{
		int i58;
		for (i58 = 0; (i58 < 2); i58 = (i58 + 1)) {
			dsp->fRec58[i58] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i59;
		for (i59 = 0; (i59 < 2); i59 = (i59 + 1)) {
			dsp->fRec57[i59] = 0.f;
			
		}
		
	}
	dsp->fConst62 = tan((40212.4f / (float)dsp->iConst0));
	dsp->fConst63 = (1.f / dsp->fConst62);
	dsp->fConst64 = (2.f * (1.f - (1.f / faustpower2_f(dsp->fConst62))));
	/* C99 loop */
	{
		int i60;
		for (i60 = 0; (i60 < 3); i60 = (i60 + 1)) {
			dsp->fRec60[i60] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i61;
		for (i61 = 0; (i61 < 3); i61 = (i61 + 1)) {
			dsp->fRec63[i61] = 0.f;
			
		}
		
	}
	dsp->fConst65 = (0.f - dsp->fConst63);
	/* C99 loop */
	{
		int i62;
		for (i62 = 0; (i62 < 2); i62 = (i62 + 1)) {
			dsp->fRec62[i62] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i63;
		for (i63 = 0; (i63 < 2); i63 = (i63 + 1)) {
			dsp->fRec61[i63] = 0.f;
			
		}
		
	}
	
}

static void initvocoder(vocoder* dsp, int samplingFreq) {
	instanceInitvocoder(dsp, samplingFreq);
}

static void buildUserInterfacevocoder(vocoder* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "atk", &dsp->fHslider1, 0.01f, 0.0001f, 0.5f, 1e-05f);
	interface->addHorizontalSlider(interface->uiInterface, "rel", &dsp->fHslider2, 0.01f, 0.0001f, 0.5f, 1e-05f);
	interface->addHorizontalSlider(interface->uiInterface, "bwratio", &dsp->fHslider0, 0.5f, 0.1f, 2.f, 0.001f);
}

static void computevocoder(vocoder* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* input0 = inputs[0];
	FAUSTFLOAT* input1 = inputs[1];
	FAUSTFLOAT* output0 = outputs[0];
	float fSlow0 = (float)dsp->fHslider0;
	float fSlow1 = (0.645744f * fSlow0);
	float fSlow2 = (1.f / (1.f + (dsp->fConst2 * (fSlow1 + dsp->fConst2))));
	float fSlow3 = (1.f + (dsp->fConst2 * (dsp->fConst2 - fSlow1)));
	float fSlow4 = exp((0.f - (dsp->fConst5 / (float)dsp->fHslider1)));
	float fSlow5 = exp((0.f - (dsp->fConst5 / (float)dsp->fHslider2)));
	float fSlow6 = (0.645744f * fSlow0);
	float fSlow7 = (1.f / (1.f + (dsp->fConst7 * (fSlow6 + dsp->fConst7))));
	float fSlow8 = (1.f + (dsp->fConst7 * (dsp->fConst7 - fSlow6)));
	float fSlow9 = (1.f / (1.f + (dsp->fConst11 * (fSlow1 + dsp->fConst11))));
	float fSlow10 = (1.f + (dsp->fConst11 * (dsp->fConst11 - fSlow1)));
	float fSlow11 = (1.f / (1.f + (dsp->fConst15 * (fSlow6 + dsp->fConst15))));
	float fSlow12 = (1.f + (dsp->fConst15 * (dsp->fConst15 - fSlow6)));
	float fSlow13 = (0.645744f * fSlow0);
	float fSlow14 = (1.f / (1.f + (dsp->fConst19 * (fSlow13 + dsp->fConst19))));
	float fSlow15 = (1.f + (dsp->fConst19 * (dsp->fConst19 - fSlow13)));
	float fSlow16 = (1.f / (1.f + (dsp->fConst23 * (fSlow1 + dsp->fConst23))));
	float fSlow17 = (1.f + (dsp->fConst23 * (dsp->fConst23 - fSlow1)));
	float fSlow18 = (0.645744f * fSlow0);
	float fSlow19 = (1.f / (1.f + (dsp->fConst27 * (fSlow18 + dsp->fConst27))));
	float fSlow20 = (1.f + (dsp->fConst27 * (dsp->fConst27 - fSlow18)));
	float fSlow21 = (1.f / (1.f + (dsp->fConst31 * (fSlow13 + dsp->fConst31))));
	float fSlow22 = (1.f + (dsp->fConst31 * (dsp->fConst31 - fSlow13)));
	float fSlow23 = (0.645744f * fSlow0);
	float fSlow24 = (1.f / (1.f + (dsp->fConst35 * (fSlow23 + dsp->fConst35))));
	float fSlow25 = (1.f + (dsp->fConst35 * (dsp->fConst35 - fSlow23)));
	float fSlow26 = (1.f / (1.f + (dsp->fConst39 * (fSlow6 + dsp->fConst39))));
	float fSlow27 = (1.f + (dsp->fConst39 * (dsp->fConst39 - fSlow6)));
	float fSlow28 = (1.f / (1.f + (dsp->fConst43 * (fSlow6 + dsp->fConst43))));
	float fSlow29 = (1.f + (dsp->fConst43 * (dsp->fConst43 - fSlow6)));
	float fSlow30 = (1.f / (1.f + (dsp->fConst47 * (fSlow1 + dsp->fConst47))));
	float fSlow31 = (1.f + (dsp->fConst47 * (dsp->fConst47 - fSlow1)));
	float fSlow32 = (1.f / (1.f + (dsp->fConst51 * (fSlow13 + dsp->fConst51))));
	float fSlow33 = (1.f + (dsp->fConst51 * (dsp->fConst51 - fSlow13)));
	float fSlow34 = (1.f / (1.f + (dsp->fConst55 * (fSlow1 + dsp->fConst55))));
	float fSlow35 = (1.f + (dsp->fConst55 * (dsp->fConst55 - fSlow1)));
	float fSlow36 = (1.f / (1.f + (dsp->fConst59 * (fSlow13 + dsp->fConst59))));
	float fSlow37 = (1.f + (dsp->fConst59 * (dsp->fConst59 - fSlow13)));
	float fSlow38 = (1.f / (1.f + (dsp->fConst63 * (fSlow6 + dsp->fConst63))));
	float fSlow39 = (1.f + (dsp->fConst63 * (dsp->fConst63 - fSlow6)));
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			float fTemp0 = (float)input1[i];
			dsp->fRec0[0] = (fTemp0 - (fSlow2 * ((fSlow3 * dsp->fRec0[2]) + (dsp->fConst3 * dsp->fRec0[1]))));
			float fTemp1 = (float)input0[i];
			dsp->fRec3[0] = (fTemp1 - (fSlow2 * ((fSlow3 * dsp->fRec3[2]) + (dsp->fConst3 * dsp->fRec3[1]))));
			float fTemp2 = fabs((fSlow2 * ((dsp->fConst2 * dsp->fRec3[0]) + (dsp->fConst4 * dsp->fRec3[2]))));
			float fTemp3 = ((dsp->fRec1[1] > fTemp2)?fSlow5:fSlow4);
			dsp->fRec2[0] = ((dsp->fRec2[1] * fTemp3) + (fTemp2 * (1.f - fTemp3)));
			dsp->fRec1[0] = dsp->fRec2[0];
			dsp->fRec4[0] = (fTemp0 - (fSlow7 * ((fSlow8 * dsp->fRec4[2]) + (dsp->fConst8 * dsp->fRec4[1]))));
			dsp->fRec7[0] = (fTemp1 - (fSlow7 * ((fSlow8 * dsp->fRec7[2]) + (dsp->fConst8 * dsp->fRec7[1]))));
			float fTemp4 = fabs((fSlow7 * ((dsp->fConst7 * dsp->fRec7[0]) + (dsp->fConst9 * dsp->fRec7[2]))));
			float fTemp5 = ((dsp->fRec5[1] > fTemp4)?fSlow5:fSlow4);
			dsp->fRec6[0] = ((dsp->fRec6[1] * fTemp5) + (fTemp4 * (1.f - fTemp5)));
			dsp->fRec5[0] = dsp->fRec6[0];
			dsp->fRec8[0] = (fTemp0 - (fSlow9 * ((fSlow10 * dsp->fRec8[2]) + (dsp->fConst12 * dsp->fRec8[1]))));
			dsp->fRec11[0] = (fTemp1 - (fSlow9 * ((fSlow10 * dsp->fRec11[2]) + (dsp->fConst12 * dsp->fRec11[1]))));
			float fTemp6 = fabs((fSlow9 * ((dsp->fConst11 * dsp->fRec11[0]) + (dsp->fConst13 * dsp->fRec11[2]))));
			float fTemp7 = ((dsp->fRec9[1] > fTemp6)?fSlow5:fSlow4);
			dsp->fRec10[0] = ((dsp->fRec10[1] * fTemp7) + (fTemp6 * (1.f - fTemp7)));
			dsp->fRec9[0] = dsp->fRec10[0];
			dsp->fRec12[0] = (fTemp0 - (fSlow11 * ((fSlow12 * dsp->fRec12[2]) + (dsp->fConst16 * dsp->fRec12[1]))));
			dsp->fRec15[0] = (fTemp1 - (fSlow11 * ((fSlow12 * dsp->fRec15[2]) + (dsp->fConst16 * dsp->fRec15[1]))));
			float fTemp8 = fabs((fSlow11 * ((dsp->fConst15 * dsp->fRec15[0]) + (dsp->fConst17 * dsp->fRec15[2]))));
			float fTemp9 = ((dsp->fRec13[1] > fTemp8)?fSlow5:fSlow4);
			dsp->fRec14[0] = ((dsp->fRec14[1] * fTemp9) + (fTemp8 * (1.f - fTemp9)));
			dsp->fRec13[0] = dsp->fRec14[0];
			dsp->fRec16[0] = (fTemp0 - (fSlow14 * ((fSlow15 * dsp->fRec16[2]) + (dsp->fConst20 * dsp->fRec16[1]))));
			dsp->fRec19[0] = (fTemp1 - (fSlow14 * ((fSlow15 * dsp->fRec19[2]) + (dsp->fConst20 * dsp->fRec19[1]))));
			float fTemp10 = fabs((fSlow14 * ((dsp->fConst19 * dsp->fRec19[0]) + (dsp->fConst21 * dsp->fRec19[2]))));
			float fTemp11 = ((dsp->fRec17[1] > fTemp10)?fSlow5:fSlow4);
			dsp->fRec18[0] = ((dsp->fRec18[1] * fTemp11) + (fTemp10 * (1.f - fTemp11)));
			dsp->fRec17[0] = dsp->fRec18[0];
			dsp->fRec20[0] = (fTemp0 - (fSlow16 * ((fSlow17 * dsp->fRec20[2]) + (dsp->fConst24 * dsp->fRec20[1]))));
			dsp->fRec23[0] = (fTemp1 - (fSlow16 * ((fSlow17 * dsp->fRec23[2]) + (dsp->fConst24 * dsp->fRec23[1]))));
			float fTemp12 = fabs((fSlow16 * ((dsp->fConst23 * dsp->fRec23[0]) + (dsp->fConst25 * dsp->fRec23[2]))));
			float fTemp13 = ((dsp->fRec21[1] > fTemp12)?fSlow5:fSlow4);
			dsp->fRec22[0] = ((dsp->fRec22[1] * fTemp13) + (fTemp12 * (1.f - fTemp13)));
			dsp->fRec21[0] = dsp->fRec22[0];
			dsp->fRec24[0] = (fTemp0 - (fSlow19 * ((fSlow20 * dsp->fRec24[2]) + (dsp->fConst28 * dsp->fRec24[1]))));
			dsp->fRec27[0] = (fTemp1 - (fSlow19 * ((fSlow20 * dsp->fRec27[2]) + (dsp->fConst28 * dsp->fRec27[1]))));
			float fTemp14 = fabs((fSlow19 * ((dsp->fConst27 * dsp->fRec27[0]) + (dsp->fConst29 * dsp->fRec27[2]))));
			float fTemp15 = ((dsp->fRec25[1] > fTemp14)?fSlow5:fSlow4);
			dsp->fRec26[0] = ((dsp->fRec26[1] * fTemp15) + (fTemp14 * (1.f - fTemp15)));
			dsp->fRec25[0] = dsp->fRec26[0];
			dsp->fRec28[0] = (fTemp0 - (fSlow21 * ((fSlow22 * dsp->fRec28[2]) + (dsp->fConst32 * dsp->fRec28[1]))));
			dsp->fRec31[0] = (fTemp1 - (fSlow21 * ((fSlow22 * dsp->fRec31[2]) + (dsp->fConst32 * dsp->fRec31[1]))));
			float fTemp16 = fabs((fSlow21 * ((dsp->fConst31 * dsp->fRec31[0]) + (dsp->fConst33 * dsp->fRec31[2]))));
			float fTemp17 = ((dsp->fRec29[1] > fTemp16)?fSlow5:fSlow4);
			dsp->fRec30[0] = ((dsp->fRec30[1] * fTemp17) + (fTemp16 * (1.f - fTemp17)));
			dsp->fRec29[0] = dsp->fRec30[0];
			dsp->fRec32[0] = (fTemp0 - (fSlow24 * ((fSlow25 * dsp->fRec32[2]) + (dsp->fConst36 * dsp->fRec32[1]))));
			dsp->fRec35[0] = (fTemp1 - (fSlow24 * ((fSlow25 * dsp->fRec35[2]) + (dsp->fConst36 * dsp->fRec35[1]))));
			float fTemp18 = fabs((fSlow24 * ((dsp->fConst35 * dsp->fRec35[0]) + (dsp->fConst37 * dsp->fRec35[2]))));
			float fTemp19 = ((dsp->fRec33[1] > fTemp18)?fSlow5:fSlow4);
			dsp->fRec34[0] = ((dsp->fRec34[1] * fTemp19) + (fTemp18 * (1.f - fTemp19)));
			dsp->fRec33[0] = dsp->fRec34[0];
			dsp->fRec36[0] = (fTemp0 - (fSlow26 * ((fSlow27 * dsp->fRec36[2]) + (dsp->fConst40 * dsp->fRec36[1]))));
			dsp->fRec39[0] = (fTemp1 - (fSlow26 * ((fSlow27 * dsp->fRec39[2]) + (dsp->fConst40 * dsp->fRec39[1]))));
			float fTemp20 = fabs((fSlow26 * ((dsp->fConst39 * dsp->fRec39[0]) + (dsp->fConst41 * dsp->fRec39[2]))));
			float fTemp21 = ((dsp->fRec37[1] > fTemp20)?fSlow5:fSlow4);
			dsp->fRec38[0] = ((dsp->fRec38[1] * fTemp21) + (fTemp20 * (1.f - fTemp21)));
			dsp->fRec37[0] = dsp->fRec38[0];
			dsp->fRec40[0] = (fTemp0 - (fSlow28 * ((fSlow29 * dsp->fRec40[2]) + (dsp->fConst44 * dsp->fRec40[1]))));
			dsp->fRec43[0] = (fTemp1 - (fSlow28 * ((fSlow29 * dsp->fRec43[2]) + (dsp->fConst44 * dsp->fRec43[1]))));
			float fTemp22 = fabs((fSlow28 * ((dsp->fConst43 * dsp->fRec43[0]) + (dsp->fConst45 * dsp->fRec43[2]))));
			float fTemp23 = ((dsp->fRec41[1] > fTemp22)?fSlow5:fSlow4);
			dsp->fRec42[0] = ((dsp->fRec42[1] * fTemp23) + (fTemp22 * (1.f - fTemp23)));
			dsp->fRec41[0] = dsp->fRec42[0];
			dsp->fRec44[0] = (fTemp0 - (fSlow30 * ((fSlow31 * dsp->fRec44[2]) + (dsp->fConst48 * dsp->fRec44[1]))));
			dsp->fRec47[0] = (fTemp1 - (fSlow30 * ((fSlow31 * dsp->fRec47[2]) + (dsp->fConst48 * dsp->fRec47[1]))));
			float fTemp24 = fabs((fSlow30 * ((dsp->fConst47 * dsp->fRec47[0]) + (dsp->fConst49 * dsp->fRec47[2]))));
			float fTemp25 = ((dsp->fRec45[1] > fTemp24)?fSlow5:fSlow4);
			dsp->fRec46[0] = ((dsp->fRec46[1] * fTemp25) + (fTemp24 * (1.f - fTemp25)));
			dsp->fRec45[0] = dsp->fRec46[0];
			dsp->fRec48[0] = (fTemp0 - (fSlow32 * ((fSlow33 * dsp->fRec48[2]) + (dsp->fConst52 * dsp->fRec48[1]))));
			dsp->fRec51[0] = (fTemp1 - (fSlow32 * ((fSlow33 * dsp->fRec51[2]) + (dsp->fConst52 * dsp->fRec51[1]))));
			float fTemp26 = fabs((fSlow32 * ((dsp->fConst51 * dsp->fRec51[0]) + (dsp->fConst53 * dsp->fRec51[2]))));
			float fTemp27 = ((dsp->fRec49[1] > fTemp26)?fSlow5:fSlow4);
			dsp->fRec50[0] = ((dsp->fRec50[1] * fTemp27) + (fTemp26 * (1.f - fTemp27)));
			dsp->fRec49[0] = dsp->fRec50[0];
			dsp->fRec52[0] = (fTemp0 - (fSlow34 * ((fSlow35 * dsp->fRec52[2]) + (dsp->fConst56 * dsp->fRec52[1]))));
			dsp->fRec55[0] = (fTemp1 - (fSlow34 * ((fSlow35 * dsp->fRec55[2]) + (dsp->fConst56 * dsp->fRec55[1]))));
			float fTemp28 = fabs((fSlow34 * ((dsp->fConst55 * dsp->fRec55[0]) + (dsp->fConst57 * dsp->fRec55[2]))));
			float fTemp29 = ((dsp->fRec53[1] > fTemp28)?fSlow5:fSlow4);
			dsp->fRec54[0] = ((dsp->fRec54[1] * fTemp29) + (fTemp28 * (1.f - fTemp29)));
			dsp->fRec53[0] = dsp->fRec54[0];
			dsp->fRec56[0] = (fTemp0 - (fSlow36 * ((fSlow37 * dsp->fRec56[2]) + (dsp->fConst60 * dsp->fRec56[1]))));
			dsp->fRec59[0] = (fTemp1 - (fSlow36 * ((fSlow37 * dsp->fRec59[2]) + (dsp->fConst60 * dsp->fRec59[1]))));
			float fTemp30 = fabs((fSlow36 * ((dsp->fConst59 * dsp->fRec59[0]) + (dsp->fConst61 * dsp->fRec59[2]))));
			float fTemp31 = ((dsp->fRec57[1] > fTemp30)?fSlow5:fSlow4);
			dsp->fRec58[0] = ((dsp->fRec58[1] * fTemp31) + (fTemp30 * (1.f - fTemp31)));
			dsp->fRec57[0] = dsp->fRec58[0];
			dsp->fRec60[0] = (fTemp0 - (fSlow38 * ((fSlow39 * dsp->fRec60[2]) + (dsp->fConst64 * dsp->fRec60[1]))));
			dsp->fRec63[0] = (fTemp1 - (fSlow38 * ((fSlow39 * dsp->fRec63[2]) + (dsp->fConst64 * dsp->fRec63[1]))));
			float fTemp32 = fabs((fSlow38 * ((dsp->fConst63 * dsp->fRec63[0]) + (dsp->fConst65 * dsp->fRec63[2]))));
			float fTemp33 = ((dsp->fRec61[1] > fTemp32)?fSlow5:fSlow4);
			dsp->fRec62[0] = ((dsp->fRec62[1] * fTemp33) + (fTemp32 * (1.f - fTemp33)));
			dsp->fRec61[0] = dsp->fRec62[0];
			output0[i] = (FAUSTFLOAT)((((((((((((((((fSlow2 * ((dsp->fRec0[2] * (0.f - (dsp->fConst2 * dsp->fRec1[0]))) + (dsp->fConst2 * (dsp->fRec0[0] * dsp->fRec1[0])))) + (fSlow7 * ((dsp->fRec4[2] * (0.f - (dsp->fConst7 * dsp->fRec5[0]))) + (dsp->fConst7 * (dsp->fRec4[0] * dsp->fRec5[0]))))) + (fSlow9 * ((dsp->fRec8[2] * (0.f - (dsp->fConst11 * dsp->fRec9[0]))) + (dsp->fConst11 * (dsp->fRec8[0] * dsp->fRec9[0]))))) + (fSlow11 * ((dsp->fRec12[2] * (0.f - (dsp->fConst15 * dsp->fRec13[0]))) + (dsp->fConst15 * (dsp->fRec12[0] * dsp->fRec13[0]))))) + (fSlow14 * ((dsp->fRec16[2] * (0.f - (dsp->fConst19 * dsp->fRec17[0]))) + (dsp->fConst19 * (dsp->fRec16[0] * dsp->fRec17[0]))))) + (fSlow16 * ((dsp->fRec20[2] * (0.f - (dsp->fConst23 * dsp->fRec21[0]))) + (dsp->fConst23 * (dsp->fRec20[0] * dsp->fRec21[0]))))) + (fSlow19 * ((dsp->fRec24[2] * (0.f - (dsp->fConst27 * dsp->fRec25[0]))) + (dsp->fConst27 * (dsp->fRec24[0] * dsp->fRec25[0]))))) + (fSlow21 * ((dsp->fRec28[2] * (0.f - (dsp->fConst31 * dsp->fRec29[0]))) + (dsp->fConst31 * (dsp->fRec28[0] * dsp->fRec29[0]))))) + (fSlow24 * ((dsp->fRec32[2] * (0.f - (dsp->fConst35 * dsp->fRec33[0]))) + (dsp->fConst35 * (dsp->fRec32[0] * dsp->fRec33[0]))))) + (fSlow26 * ((dsp->fRec36[2] * (0.f - (dsp->fConst39 * dsp->fRec37[0]))) + (dsp->fConst39 * (dsp->fRec36[0] * dsp->fRec37[0]))))) + (fSlow28 * ((dsp->fRec40[2] * (0.f - (dsp->fConst43 * dsp->fRec41[0]))) + (dsp->fConst43 * (dsp->fRec40[0] * dsp->fRec41[0]))))) + (fSlow30 * ((dsp->fRec44[2] * (0.f - (dsp->fConst47 * dsp->fRec45[0]))) + (dsp->fConst47 * (dsp->fRec44[0] * dsp->fRec45[0]))))) + (fSlow32 * ((dsp->fRec48[2] * (0.f - (dsp->fConst51 * dsp->fRec49[0]))) + (dsp->fConst51 * (dsp->fRec48[0] * dsp->fRec49[0]))))) + (fSlow34 * ((dsp->fRec52[2] * (0.f - (dsp->fConst55 * dsp->fRec53[0]))) + (dsp->fConst55 * (dsp->fRec52[0] * dsp->fRec53[0]))))) + (fSlow36 * ((dsp->fRec56[2] * (0.f - (dsp->fConst59 * dsp->fRec57[0]))) + (dsp->fConst59 * (dsp->fRec56[0] * dsp->fRec57[0]))))) + (fSlow38 * ((dsp->fRec60[2] * (0.f - (dsp->fConst63 * dsp->fRec61[0]))) + (dsp->fConst63 * (dsp->fRec60[0] * dsp->fRec61[0])))));
			dsp->fRec0[2] = dsp->fRec0[1];
			dsp->fRec0[1] = dsp->fRec0[0];
			dsp->fRec3[2] = dsp->fRec3[1];
			dsp->fRec3[1] = dsp->fRec3[0];
			dsp->fRec2[1] = dsp->fRec2[0];
			dsp->fRec1[1] = dsp->fRec1[0];
			dsp->fRec4[2] = dsp->fRec4[1];
			dsp->fRec4[1] = dsp->fRec4[0];
			dsp->fRec7[2] = dsp->fRec7[1];
			dsp->fRec7[1] = dsp->fRec7[0];
			dsp->fRec6[1] = dsp->fRec6[0];
			dsp->fRec5[1] = dsp->fRec5[0];
			dsp->fRec8[2] = dsp->fRec8[1];
			dsp->fRec8[1] = dsp->fRec8[0];
			dsp->fRec11[2] = dsp->fRec11[1];
			dsp->fRec11[1] = dsp->fRec11[0];
			dsp->fRec10[1] = dsp->fRec10[0];
			dsp->fRec9[1] = dsp->fRec9[0];
			dsp->fRec12[2] = dsp->fRec12[1];
			dsp->fRec12[1] = dsp->fRec12[0];
			dsp->fRec15[2] = dsp->fRec15[1];
			dsp->fRec15[1] = dsp->fRec15[0];
			dsp->fRec14[1] = dsp->fRec14[0];
			dsp->fRec13[1] = dsp->fRec13[0];
			dsp->fRec16[2] = dsp->fRec16[1];
			dsp->fRec16[1] = dsp->fRec16[0];
			dsp->fRec19[2] = dsp->fRec19[1];
			dsp->fRec19[1] = dsp->fRec19[0];
			dsp->fRec18[1] = dsp->fRec18[0];
			dsp->fRec17[1] = dsp->fRec17[0];
			dsp->fRec20[2] = dsp->fRec20[1];
			dsp->fRec20[1] = dsp->fRec20[0];
			dsp->fRec23[2] = dsp->fRec23[1];
			dsp->fRec23[1] = dsp->fRec23[0];
			dsp->fRec22[1] = dsp->fRec22[0];
			dsp->fRec21[1] = dsp->fRec21[0];
			dsp->fRec24[2] = dsp->fRec24[1];
			dsp->fRec24[1] = dsp->fRec24[0];
			dsp->fRec27[2] = dsp->fRec27[1];
			dsp->fRec27[1] = dsp->fRec27[0];
			dsp->fRec26[1] = dsp->fRec26[0];
			dsp->fRec25[1] = dsp->fRec25[0];
			dsp->fRec28[2] = dsp->fRec28[1];
			dsp->fRec28[1] = dsp->fRec28[0];
			dsp->fRec31[2] = dsp->fRec31[1];
			dsp->fRec31[1] = dsp->fRec31[0];
			dsp->fRec30[1] = dsp->fRec30[0];
			dsp->fRec29[1] = dsp->fRec29[0];
			dsp->fRec32[2] = dsp->fRec32[1];
			dsp->fRec32[1] = dsp->fRec32[0];
			dsp->fRec35[2] = dsp->fRec35[1];
			dsp->fRec35[1] = dsp->fRec35[0];
			dsp->fRec34[1] = dsp->fRec34[0];
			dsp->fRec33[1] = dsp->fRec33[0];
			dsp->fRec36[2] = dsp->fRec36[1];
			dsp->fRec36[1] = dsp->fRec36[0];
			dsp->fRec39[2] = dsp->fRec39[1];
			dsp->fRec39[1] = dsp->fRec39[0];
			dsp->fRec38[1] = dsp->fRec38[0];
			dsp->fRec37[1] = dsp->fRec37[0];
			dsp->fRec40[2] = dsp->fRec40[1];
			dsp->fRec40[1] = dsp->fRec40[0];
			dsp->fRec43[2] = dsp->fRec43[1];
			dsp->fRec43[1] = dsp->fRec43[0];
			dsp->fRec42[1] = dsp->fRec42[0];
			dsp->fRec41[1] = dsp->fRec41[0];
			dsp->fRec44[2] = dsp->fRec44[1];
			dsp->fRec44[1] = dsp->fRec44[0];
			dsp->fRec47[2] = dsp->fRec47[1];
			dsp->fRec47[1] = dsp->fRec47[0];
			dsp->fRec46[1] = dsp->fRec46[0];
			dsp->fRec45[1] = dsp->fRec45[0];
			dsp->fRec48[2] = dsp->fRec48[1];
			dsp->fRec48[1] = dsp->fRec48[0];
			dsp->fRec51[2] = dsp->fRec51[1];
			dsp->fRec51[1] = dsp->fRec51[0];
			dsp->fRec50[1] = dsp->fRec50[0];
			dsp->fRec49[1] = dsp->fRec49[0];
			dsp->fRec52[2] = dsp->fRec52[1];
			dsp->fRec52[1] = dsp->fRec52[0];
			dsp->fRec55[2] = dsp->fRec55[1];
			dsp->fRec55[1] = dsp->fRec55[0];
			dsp->fRec54[1] = dsp->fRec54[0];
			dsp->fRec53[1] = dsp->fRec53[0];
			dsp->fRec56[2] = dsp->fRec56[1];
			dsp->fRec56[1] = dsp->fRec56[0];
			dsp->fRec59[2] = dsp->fRec59[1];
			dsp->fRec59[1] = dsp->fRec59[0];
			dsp->fRec58[1] = dsp->fRec58[0];
			dsp->fRec57[1] = dsp->fRec57[0];
			dsp->fRec60[2] = dsp->fRec60[1];
			dsp->fRec60[1] = dsp->fRec60[0];
			dsp->fRec63[2] = dsp->fRec63[1];
			dsp->fRec63[1] = dsp->fRec63[0];
			dsp->fRec62[1] = dsp->fRec62[0];
			dsp->fRec61[1] = dsp->fRec61[0];
			
		}
		
	}
	
}


int sp_vocoder_create(sp_vocoder **p)
{
    *p = malloc(sizeof(sp_vocoder));
    return SP_OK;
}

int sp_vocoder_destroy(sp_vocoder **p)
{
    sp_vocoder *pp = *p;
    vocoder *dsp = pp->faust;
    deletevocoder (dsp);
    free(*p);
    return SP_OK;
}

int sp_vocoder_init(sp_data *sp, sp_vocoder *p)
{
    vocoder *dsp = newvocoder(); 
    UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfacevocoder(dsp, &UI);
    initvocoder(dsp, sp->sr);

     
    p->atk = p->args[0]; 
    p->rel = p->args[1]; 
    p->bwratio = p->args[2];

    p->faust = dsp;
    return SP_OK;
}

int sp_vocoder_compute(sp_data *sp, sp_vocoder *p, SPFLOAT *source, SPFLOAT *excite, SPFLOAT *out)
{

    vocoder *dsp = p->faust;
    SPFLOAT *faust_out[] = {out};
    SPFLOAT *faust_in[] = {source, excite};
    computevocoder(dsp, 1, faust_in, faust_out);
    return SP_OK;
}


int sp_waveset_create(sp_waveset **p)
{
    *p = malloc(sizeof(sp_waveset));
    return SP_OK;
}

int sp_waveset_destroy(sp_waveset **p)
{
    sp_waveset *pp = *p;
    sp_auxdata_free(&pp->auxch);
    free(*p);
    return SP_OK;
}

int sp_waveset_init(sp_data *sp, sp_waveset *p, SPFLOAT ilen)
{
    p->length = 1 + (sp->sr * ilen);

    sp_auxdata_alloc(&p->auxch, p->length * sizeof(SPFLOAT));
    p->cnt = 1;
    p->start = 0;
    p->current = 0;
    p->end = 0;
    p->direction = 1;
    p->lastsamp = 1.0;
    p->noinsert = 0;
    return SP_OK;
}

int sp_waveset_compute(sp_data *sp, sp_waveset *p, SPFLOAT *in, SPFLOAT *out)
{
    int index = p->end;
    SPFLOAT *insert = (SPFLOAT*)(p->auxch.ptr) + index;

    if (p->noinsert) goto output;
    *insert++ = *in;
    if (++index ==  p->start) {
        p->noinsert = 1;
    }
    if (index==p->length) {  
        index = 0;
        insert = (SPFLOAT*)(p->auxch.ptr);
    }

    output:

    p->end = index;
    index = p->current;
    insert = (SPFLOAT*)(p->auxch.ptr) + index;
    SPFLOAT samp = *insert++;
    index++;

    if (index==p->length) {
        index = 0;
        insert = (SPFLOAT*)(p->auxch.ptr);
        p->noinsert = 0;
    }

    if (samp != 0.0 && p->lastsamp*samp < 0.0) {
        if (p->direction == 1) {
            p->direction = -1;
        } else {
            p->direction = 1;
            if (++p->cnt > p->rep) {
                p->cnt = 1;
                p->start = index;
                p->noinsert = 0;
            } else { index = p->start;
                insert = (SPFLOAT*)(p->auxch.ptr) + index;
            }
        }
    }

    if (samp != 0.0) p->lastsamp = samp;
    *out = samp;
    p->current = index;

    return SP_OK;
}



#define WAVIN_BUFSIZE 1024

struct sp_wavin {
    SPFLOAT buf[WAVIN_BUFSIZE];
    int count;
    drwav wav;
    drwav_uint64 pos;
};

int sp_wavin_create(sp_wavin **p)
{
    *p = malloc(sizeof(sp_wavin));
    return SP_OK;
}

int sp_wavin_destroy(sp_wavin **p)
{
    drwav_uninit(&(*p)->wav);
    free(*p);
    return SP_OK;
}

int sp_wavin_init(sp_data *sp, sp_wavin *p, const char *filename)
{
    p->count = 0;
    p->pos = 0;
    drwav_init_file(&p->wav, filename);
    return SP_OK;
}

int sp_wavin_compute(sp_data *sp, sp_wavin *p, SPFLOAT *in, SPFLOAT *out)
{
    if(p->pos > p->wav.totalSampleCount) {
        *out = 0;
        return SP_OK;
    }
    if(p->count == 0) {
        drwav_read_f32(&p->wav, WAVIN_BUFSIZE, p->buf);
    }

    *out = p->buf[p->count];
    p->count = (p->count + 1) % WAVIN_BUFSIZE;
    p->pos++;
    return SP_OK;
}


#define WAVOUT_BUFSIZE 1024

struct sp_wavout {
    drwav *wav;
    drwav_data_format format;
    SPFLOAT buf[WAVOUT_BUFSIZE];
    int count;
};

int sp_wavout_create(sp_wavout **p)
{
    *p = malloc(sizeof(sp_wavout));
    return SP_OK;
}

int sp_wavout_destroy(sp_wavout **p)
{
    /* write any remaining samples */
    if((*p)->count != 0) {
        drwav_write((*p)->wav, (*p)->count, (*p)->buf);
    }
    drwav_close((*p)->wav);
    free(*p);
    return SP_OK;
}

int sp_wavout_init(sp_data *sp, sp_wavout *p, const char *filename)
{
    p->count = 0;
    p->format.container = drwav_container_riff;
    p->format.format = DR_WAVE_FORMAT_IEEE_FLOAT;
    p->format.channels = 1;
    p->format.sampleRate = sp->sr;
    p->format.bitsPerSample = 32;
    p->wav = drwav_open_file_write(filename, &p->format);
    return SP_OK;
}

int sp_wavout_compute(sp_data *sp, sp_wavout *p, SPFLOAT *in, SPFLOAT *out)
{
    *out = *in;
    if(p->count == WAVOUT_BUFSIZE) {
        drwav_write(p->wav, WAVOUT_BUFSIZE, p->buf);
        p->count = 0;
    }
    p->buf[p->count] = *in;
    p->count++;
    return SP_OK;
}


static void update(sp_data *sp, sp_wpkorg35 *wpk)
{
	/* prewarp for BZT */
	SPFLOAT wd = 2*M_PI*wpk->cutoff;          
	SPFLOAT T  = 1.0/(SPFLOAT)sp->sr;             
	SPFLOAT wa = (2/T)*tan(wd*T/2); 
	SPFLOAT g  = wa*T/2.0;    

	/* the feedforward coeff in the VA One Pole */
	SPFLOAT G = g/(1.0 + g);

    /* set alphas */
    wpk->lpf1_a = G;
    wpk->lpf2_a = G;
    wpk->hpf_a = G;

    /* set betas */
	wpk->lpf2_b = (wpk->res - wpk->res*G)/(1.0 + g);
	wpk->hpf_b = -1.0/(1.0 + g);

	wpk->alpha = 1.0/(1.0 - wpk->res*G + wpk->res*G*G); ;
}

SPFLOAT wpk_doFilter(sp_wpkorg35 *wpk)
{
    return 0.0;
}

int sp_wpkorg35_create(sp_wpkorg35 **p)
{
    *p = malloc(sizeof(sp_wpkorg35));
    return SP_OK;
}

int sp_wpkorg35_destroy(sp_wpkorg35 **p)
{
    free(*p);
    return SP_OK;
}

int sp_wpkorg35_init(sp_data *sp, sp_wpkorg35 *p)
{
    p->alpha = 0.0;
    p->pcutoff = p->cutoff = 1000; 
    p->pres = p->res = 1.0; 

    /* reset memory for filters */
    p->lpf1_z = 0;
    p->lpf2_z = 0;
    p->hpf_z = 0;


    /* initialize LPF1 */

    p->lpf1_a = 1.0;
    p->lpf1_z = 0.0;
    
    /* initialize LPF2 */

    p->lpf2_a = 1.0;
    p->lpf2_b = 1.0;
    p->lpf2_z = 0.0;

    p->nonlinear = 0;
    p->saturation = 0;

    /* update filters */
    update(sp, p);
    return SP_OK;
}

int sp_wpkorg35_compute(sp_data *sp, sp_wpkorg35 *p, SPFLOAT *in, SPFLOAT *out)
{
    /* TODO: add previous values */

    if(p->pcutoff != p->cutoff || p->pres != p->res) update(sp, p);

    /* initialize variables */
    SPFLOAT y1 = 0.0;
    SPFLOAT S35 = 0.0;
    SPFLOAT u = 0.0;
    SPFLOAT y = 0.0;
    SPFLOAT vn = 0.0;

    /* process input through LPF1 */
    vn = (*in - p->lpf1_z) * p->lpf1_a;
    y1 = vn + p->lpf1_z;
    p->lpf1_z = y1 + vn;

    /* form feedback value */
    
    S35 = (p->hpf_z * p->hpf_b) + (p->lpf2_z * p->lpf2_b); 

    /* Calculate u */
    u = p->alpha * (y1 + S35);

    /* Naive NLP */

    if(p->saturation > 0) {
        u = tanh(p->saturation * u);
    }

    /* Feed it to LPF2 */
    vn = (u - p->lpf2_z) * p->lpf2_a;
    y = (vn + p->lpf2_z);
    p->lpf2_z = y + vn;
    y *= p->res;

    /* Feed y to HPF2 */

    vn = (y - p->hpf_z) * p->hpf_a;
    p->hpf_z = vn + (vn + p->hpf_z); 

    /* Auto-normalize */

    if(p->res > 0) {
        y *= 1.0 / p->res;
    }

    *out = y;

    p->pcutoff = p->cutoff;
    p->pres = p->res;
    return SP_OK;
}

typedef struct {
	
	float fVec1[32768];
	float fVec4[32768];
	float fVec8[32768];
	float fVec6[16384];
	float fVec10[16384];
	float fVec12[16384];
	float fVec14[16384];
	float fVec16[16384];
	float fVec0[8192];
	float fVec2[8192];
	float fVec5[4096];
	float fVec7[4096];
	float fVec9[4096];
	float fVec13[4096];
	float fVec15[4096];
	float fVec3[2048];
	float fVec11[2048];
	float fVec17[2048];
	float fRec4[3];
	float fRec5[3];
	float fRec6[3];
	float fRec7[3];
	float fRec8[3];
	float fRec9[3];
	float fRec10[3];
	float fRec11[3];
	float fRec3[3];
	float fRec2[3];
	float fRec45[3];
	float fRec44[3];
	float fRec0[2];
	float fRec1[2];
	float fRec15[2];
	float fRec14[2];
	float fRec12[2];
	float fRec19[2];
	float fRec18[2];
	float fRec16[2];
	float fRec23[2];
	float fRec22[2];
	float fRec20[2];
	float fRec27[2];
	float fRec26[2];
	float fRec24[2];
	float fRec31[2];
	float fRec30[2];
	float fRec28[2];
	float fRec35[2];
	float fRec34[2];
	float fRec32[2];
	float fRec39[2];
	float fRec38[2];
	float fRec36[2];
	float fRec43[2];
	float fRec42[2];
	float fRec40[2];
	FAUSTFLOAT fHslider0;
	FAUSTFLOAT fHslider1;
	int IOTA;
	int fSamplingFreq;
	int iConst0;
	float fConst1;
	FAUSTFLOAT fHslider2;
	FAUSTFLOAT fHslider3;
	FAUSTFLOAT fHslider4;
	FAUSTFLOAT fHslider5;
	float fConst2;
	float fConst3;
	FAUSTFLOAT fHslider6;
	float fConst4;
	FAUSTFLOAT fHslider7;
	FAUSTFLOAT fHslider8;
	float fConst5;
	FAUSTFLOAT fHslider9;
	float fConst6;
	int iConst7;
	float fConst8;
	FAUSTFLOAT fHslider10;
	int iConst9;
	float fConst10;
	float fConst11;
	float fConst12;
	int iConst13;
	int iConst14;
	float fConst15;
	float fConst16;
	float fConst17;
	int iConst18;
	int iConst19;
	float fConst20;
	float fConst21;
	float fConst22;
	int iConst23;
	int iConst24;
	float fConst25;
	float fConst26;
	float fConst27;
	int iConst28;
	int iConst29;
	float fConst30;
	float fConst31;
	float fConst32;
	int iConst33;
	int iConst34;
	float fConst35;
	float fConst36;
	float fConst37;
	int iConst38;
	int iConst39;
	float fConst40;
	float fConst41;
	float fConst42;
	int iConst43;
	int iConst44;
	
} zitarev;

static zitarev* newzitarev() { 
	zitarev* dsp = (zitarev*)malloc(sizeof(zitarev));
	return dsp;
}

static void deletezitarev(zitarev* dsp) { 
	free(dsp);
}

static void instanceInitzitarev(zitarev* dsp, int samplingFreq) {
	dsp->fSamplingFreq = samplingFreq;
	dsp->fHslider0 = (FAUSTFLOAT)-20.;
	/* C99 loop */
	{
		int i0;
		for (i0 = 0; (i0 < 2); i0 = (i0 + 1)) {
			dsp->fRec0[i0] = 0.f;
			
		}
		
	}
	dsp->fHslider1 = (FAUSTFLOAT)1.;
	/* C99 loop */
	{
		int i1;
		for (i1 = 0; (i1 < 2); i1 = (i1 + 1)) {
			dsp->fRec1[i1] = 0.f;
			
		}
		
	}
	dsp->IOTA = 0;
	/* C99 loop */
	{
		int i2;
		for (i2 = 0; (i2 < 8192); i2 = (i2 + 1)) {
			dsp->fVec0[i2] = 0.f;
			
		}
		
	}
	dsp->iConst0 = min(192000, max(1, dsp->fSamplingFreq));
	dsp->fConst1 = (6.28319f / (float)dsp->iConst0);
	dsp->fHslider2 = (FAUSTFLOAT)1500.;
	dsp->fHslider3 = (FAUSTFLOAT)0.;
	dsp->fHslider4 = (FAUSTFLOAT)315.;
	dsp->fHslider5 = (FAUSTFLOAT)0.;
	dsp->fConst2 = floorf((0.5f + (0.219991f * (float)dsp->iConst0)));
	dsp->fConst3 = ((0.f - (6.90776f * dsp->fConst2)) / (float)dsp->iConst0);
	dsp->fHslider6 = (FAUSTFLOAT)2.;
	dsp->fConst4 = (6.28319f / (float)dsp->iConst0);
	dsp->fHslider7 = (FAUSTFLOAT)6000.;
	dsp->fHslider8 = (FAUSTFLOAT)3.;
	dsp->fConst5 = (3.14159f / (float)dsp->iConst0);
	dsp->fHslider9 = (FAUSTFLOAT)200.;
	/* C99 loop */
	{
		int i3;
		for (i3 = 0; (i3 < 2); i3 = (i3 + 1)) {
			dsp->fRec15[i3] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i4;
		for (i4 = 0; (i4 < 2); i4 = (i4 + 1)) {
			dsp->fRec14[i4] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i5;
		for (i5 = 0; (i5 < 32768); i5 = (i5 + 1)) {
			dsp->fVec1[i5] = 0.f;
			
		}
		
	}
	dsp->fConst6 = floorf((0.5f + (0.019123f * (float)dsp->iConst0)));
	dsp->iConst7 = (int)((int)(dsp->fConst2 - dsp->fConst6) & 32767);
	/* C99 loop */
	{
		int i6;
		for (i6 = 0; (i6 < 8192); i6 = (i6 + 1)) {
			dsp->fVec2[i6] = 0.f;
			
		}
		
	}
	dsp->fConst8 = (0.001f * (float)dsp->iConst0);
	dsp->fHslider10 = (FAUSTFLOAT)60.;
	/* C99 loop */
	{
		int i7;
		for (i7 = 0; (i7 < 2048); i7 = (i7 + 1)) {
			dsp->fVec3[i7] = 0.f;
			
		}
		
	}
	dsp->iConst9 = (int)((int)(dsp->fConst6 - 1.f) & 2047);
	/* C99 loop */
	{
		int i8;
		for (i8 = 0; (i8 < 2); i8 = (i8 + 1)) {
			dsp->fRec12[i8] = 0.f;
			
		}
		
	}
	dsp->fConst10 = floorf((0.5f + (0.256891f * (float)dsp->iConst0)));
	dsp->fConst11 = ((0.f - (6.90776f * dsp->fConst10)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i9;
		for (i9 = 0; (i9 < 2); i9 = (i9 + 1)) {
			dsp->fRec19[i9] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i10;
		for (i10 = 0; (i10 < 2); i10 = (i10 + 1)) {
			dsp->fRec18[i10] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i11;
		for (i11 = 0; (i11 < 32768); i11 = (i11 + 1)) {
			dsp->fVec4[i11] = 0.f;
			
		}
		
	}
	dsp->fConst12 = floorf((0.5f + (0.027333f * (float)dsp->iConst0)));
	dsp->iConst13 = (int)((int)(dsp->fConst10 - dsp->fConst12) & 32767);
	/* C99 loop */
	{
		int i12;
		for (i12 = 0; (i12 < 4096); i12 = (i12 + 1)) {
			dsp->fVec5[i12] = 0.f;
			
		}
		
	}
	dsp->iConst14 = (int)((int)(dsp->fConst12 - 1.f) & 4095);
	/* C99 loop */
	{
		int i13;
		for (i13 = 0; (i13 < 2); i13 = (i13 + 1)) {
			dsp->fRec16[i13] = 0.f;
			
		}
		
	}
	dsp->fConst15 = floorf((0.5f + (0.192303f * (float)dsp->iConst0)));
	dsp->fConst16 = ((0.f - (6.90776f * dsp->fConst15)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i14;
		for (i14 = 0; (i14 < 2); i14 = (i14 + 1)) {
			dsp->fRec23[i14] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i15;
		for (i15 = 0; (i15 < 2); i15 = (i15 + 1)) {
			dsp->fRec22[i15] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i16;
		for (i16 = 0; (i16 < 16384); i16 = (i16 + 1)) {
			dsp->fVec6[i16] = 0.f;
			
		}
		
	}
	dsp->fConst17 = floorf((0.5f + (0.029291f * (float)dsp->iConst0)));
	dsp->iConst18 = (int)((int)(dsp->fConst15 - dsp->fConst17) & 16383);
	/* C99 loop */
	{
		int i17;
		for (i17 = 0; (i17 < 4096); i17 = (i17 + 1)) {
			dsp->fVec7[i17] = 0.f;
			
		}
		
	}
	dsp->iConst19 = (int)((int)(dsp->fConst17 - 1.f) & 4095);
	/* C99 loop */
	{
		int i18;
		for (i18 = 0; (i18 < 2); i18 = (i18 + 1)) {
			dsp->fRec20[i18] = 0.f;
			
		}
		
	}
	dsp->fConst20 = floorf((0.5f + (0.210389f * (float)dsp->iConst0)));
	dsp->fConst21 = ((0.f - (6.90776f * dsp->fConst20)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i19;
		for (i19 = 0; (i19 < 2); i19 = (i19 + 1)) {
			dsp->fRec27[i19] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i20;
		for (i20 = 0; (i20 < 2); i20 = (i20 + 1)) {
			dsp->fRec26[i20] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i21;
		for (i21 = 0; (i21 < 32768); i21 = (i21 + 1)) {
			dsp->fVec8[i21] = 0.f;
			
		}
		
	}
	dsp->fConst22 = floorf((0.5f + (0.024421f * (float)dsp->iConst0)));
	dsp->iConst23 = (int)((int)(dsp->fConst20 - dsp->fConst22) & 32767);
	/* C99 loop */
	{
		int i22;
		for (i22 = 0; (i22 < 4096); i22 = (i22 + 1)) {
			dsp->fVec9[i22] = 0.f;
			
		}
		
	}
	dsp->iConst24 = (int)((int)(dsp->fConst22 - 1.f) & 4095);
	/* C99 loop */
	{
		int i23;
		for (i23 = 0; (i23 < 2); i23 = (i23 + 1)) {
			dsp->fRec24[i23] = 0.f;
			
		}
		
	}
	dsp->fConst25 = floorf((0.5f + (0.125f * (float)dsp->iConst0)));
	dsp->fConst26 = ((0.f - (6.90776f * dsp->fConst25)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i24;
		for (i24 = 0; (i24 < 2); i24 = (i24 + 1)) {
			dsp->fRec31[i24] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i25;
		for (i25 = 0; (i25 < 2); i25 = (i25 + 1)) {
			dsp->fRec30[i25] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i26;
		for (i26 = 0; (i26 < 16384); i26 = (i26 + 1)) {
			dsp->fVec10[i26] = 0.f;
			
		}
		
	}
	dsp->fConst27 = floorf((0.5f + (0.013458f * (float)dsp->iConst0)));
	dsp->iConst28 = (int)((int)(dsp->fConst25 - dsp->fConst27) & 16383);
	/* C99 loop */
	{
		int i27;
		for (i27 = 0; (i27 < 2048); i27 = (i27 + 1)) {
			dsp->fVec11[i27] = 0.f;
			
		}
		
	}
	dsp->iConst29 = (int)((int)(dsp->fConst27 - 1.f) & 2047);
	/* C99 loop */
	{
		int i28;
		for (i28 = 0; (i28 < 2); i28 = (i28 + 1)) {
			dsp->fRec28[i28] = 0.f;
			
		}
		
	}
	dsp->fConst30 = floorf((0.5f + (0.127837f * (float)dsp->iConst0)));
	dsp->fConst31 = ((0.f - (6.90776f * dsp->fConst30)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i29;
		for (i29 = 0; (i29 < 2); i29 = (i29 + 1)) {
			dsp->fRec35[i29] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i30;
		for (i30 = 0; (i30 < 2); i30 = (i30 + 1)) {
			dsp->fRec34[i30] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i31;
		for (i31 = 0; (i31 < 16384); i31 = (i31 + 1)) {
			dsp->fVec12[i31] = 0.f;
			
		}
		
	}
	dsp->fConst32 = floorf((0.5f + (0.031604f * (float)dsp->iConst0)));
	dsp->iConst33 = (int)((int)(dsp->fConst30 - dsp->fConst32) & 16383);
	/* C99 loop */
	{
		int i32;
		for (i32 = 0; (i32 < 4096); i32 = (i32 + 1)) {
			dsp->fVec13[i32] = 0.f;
			
		}
		
	}
	dsp->iConst34 = (int)((int)(dsp->fConst32 - 1.f) & 4095);
	/* C99 loop */
	{
		int i33;
		for (i33 = 0; (i33 < 2); i33 = (i33 + 1)) {
			dsp->fRec32[i33] = 0.f;
			
		}
		
	}
	dsp->fConst35 = floorf((0.5f + (0.174713f * (float)dsp->iConst0)));
	dsp->fConst36 = ((0.f - (6.90776f * dsp->fConst35)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i34;
		for (i34 = 0; (i34 < 2); i34 = (i34 + 1)) {
			dsp->fRec39[i34] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i35;
		for (i35 = 0; (i35 < 2); i35 = (i35 + 1)) {
			dsp->fRec38[i35] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i36;
		for (i36 = 0; (i36 < 16384); i36 = (i36 + 1)) {
			dsp->fVec14[i36] = 0.f;
			
		}
		
	}
	dsp->fConst37 = floorf((0.5f + (0.022904f * (float)dsp->iConst0)));
	dsp->iConst38 = (int)((int)(dsp->fConst35 - dsp->fConst37) & 16383);
	/* C99 loop */
	{
		int i37;
		for (i37 = 0; (i37 < 4096); i37 = (i37 + 1)) {
			dsp->fVec15[i37] = 0.f;
			
		}
		
	}
	dsp->iConst39 = (int)((int)(dsp->fConst37 - 1.f) & 4095);
	/* C99 loop */
	{
		int i38;
		for (i38 = 0; (i38 < 2); i38 = (i38 + 1)) {
			dsp->fRec36[i38] = 0.f;
			
		}
		
	}
	dsp->fConst40 = floorf((0.5f + (0.153129f * (float)dsp->iConst0)));
	dsp->fConst41 = ((0.f - (6.90776f * dsp->fConst40)) / (float)dsp->iConst0);
	/* C99 loop */
	{
		int i39;
		for (i39 = 0; (i39 < 2); i39 = (i39 + 1)) {
			dsp->fRec43[i39] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i40;
		for (i40 = 0; (i40 < 2); i40 = (i40 + 1)) {
			dsp->fRec42[i40] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i41;
		for (i41 = 0; (i41 < 16384); i41 = (i41 + 1)) {
			dsp->fVec16[i41] = 0.f;
			
		}
		
	}
	dsp->fConst42 = floorf((0.5f + (0.020346f * (float)dsp->iConst0)));
	dsp->iConst43 = (int)((int)(dsp->fConst40 - dsp->fConst42) & 16383);
	/* C99 loop */
	{
		int i42;
		for (i42 = 0; (i42 < 2048); i42 = (i42 + 1)) {
			dsp->fVec17[i42] = 0.f;
			
		}
		
	}
	dsp->iConst44 = (int)((int)(dsp->fConst42 - 1.f) & 2047);
	/* C99 loop */
	{
		int i43;
		for (i43 = 0; (i43 < 2); i43 = (i43 + 1)) {
			dsp->fRec40[i43] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i44;
		for (i44 = 0; (i44 < 3); i44 = (i44 + 1)) {
			dsp->fRec4[i44] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i45;
		for (i45 = 0; (i45 < 3); i45 = (i45 + 1)) {
			dsp->fRec5[i45] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i46;
		for (i46 = 0; (i46 < 3); i46 = (i46 + 1)) {
			dsp->fRec6[i46] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i47;
		for (i47 = 0; (i47 < 3); i47 = (i47 + 1)) {
			dsp->fRec7[i47] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i48;
		for (i48 = 0; (i48 < 3); i48 = (i48 + 1)) {
			dsp->fRec8[i48] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i49;
		for (i49 = 0; (i49 < 3); i49 = (i49 + 1)) {
			dsp->fRec9[i49] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i50;
		for (i50 = 0; (i50 < 3); i50 = (i50 + 1)) {
			dsp->fRec10[i50] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i51;
		for (i51 = 0; (i51 < 3); i51 = (i51 + 1)) {
			dsp->fRec11[i51] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i52;
		for (i52 = 0; (i52 < 3); i52 = (i52 + 1)) {
			dsp->fRec3[i52] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i53;
		for (i53 = 0; (i53 < 3); i53 = (i53 + 1)) {
			dsp->fRec2[i53] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i54;
		for (i54 = 0; (i54 < 3); i54 = (i54 + 1)) {
			dsp->fRec45[i54] = 0.f;
			
		}
		
	}
	/* C99 loop */
	{
		int i55;
		for (i55 = 0; (i55 < 3); i55 = (i55 + 1)) {
			dsp->fRec44[i55] = 0.f;
			
		}
		
	}
	
}

static void initzitarev(zitarev* dsp, int samplingFreq) {
	instanceInitzitarev(dsp, samplingFreq);
}

static void buildUserInterfacezitarev(zitarev* dsp, UIGlue* interface) {
	interface->addHorizontalSlider(interface->uiInterface, "in_delay", &dsp->fHslider10, 60.f, 10.f, 100.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "lf_x", &dsp->fHslider9, 200.f, 50.f, 1000.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "rt60_low", &dsp->fHslider8, 3.f, 1.f, 8.f, 0.1f);
	interface->addHorizontalSlider(interface->uiInterface, "rt60_mid", &dsp->fHslider6, 2.f, 1.f, 8.f, 0.1f);
	interface->addHorizontalSlider(interface->uiInterface, "hf_damping", &dsp->fHslider7, 6000.f, 1500.f, 47040.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "eq1_freq", &dsp->fHslider4, 315.f, 40.f, 2500.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "eq1_level", &dsp->fHslider5, 0.f, -15.f, 15.f, 0.1f);
	interface->addHorizontalSlider(interface->uiInterface, "eq2_freq", &dsp->fHslider2, 1500.f, 160.f, 10000.f, 1.f);
	interface->addHorizontalSlider(interface->uiInterface, "eq2_level", &dsp->fHslider3, 0.f, -15.f, 15.f, 0.1f);
	interface->addHorizontalSlider(interface->uiInterface, "mix", &dsp->fHslider1, 1.f, 0.f, 1.f, 0.001f);
	interface->addHorizontalSlider(interface->uiInterface, "level", &dsp->fHslider0, -20.f, -70.f, 40.f, 0.1f);
}

static void computezitarev(zitarev* dsp, int count, FAUSTFLOAT** inputs, FAUSTFLOAT** outputs) {
	FAUSTFLOAT* input0 = inputs[0];
	FAUSTFLOAT* input1 = inputs[1];
	FAUSTFLOAT* output0 = outputs[0];
	FAUSTFLOAT* output1 = outputs[1];
	float fSlow0 = (0.001f * powf(10.f, (0.05f * (float)dsp->fHslider0)));
	float fSlow1 = (0.001f * (float)dsp->fHslider1);
	float fSlow2 = (float)dsp->fHslider2;
	float fSlow3 = powf(10.f, (0.05f * (float)dsp->fHslider3));
	float fSlow4 = (dsp->fConst1 * (fSlow2 / sqrtf(max(0.f, fSlow3))));
	float fSlow5 = ((1.f - fSlow4) / (1.f + fSlow4));
	float fSlow6 = ((0.f - cosf((dsp->fConst1 * fSlow2))) * (1.f + fSlow5));
	float fSlow7 = (float)dsp->fHslider4;
	float fSlow8 = powf(10.f, (0.05f * (float)dsp->fHslider5));
	float fSlow9 = (dsp->fConst1 * (fSlow7 / sqrtf(max(0.f, fSlow8))));
	float fSlow10 = ((1.f - fSlow9) / (1.f + fSlow9));
	float fSlow11 = ((0.f - cosf((dsp->fConst1 * fSlow7))) * (1.f + fSlow10));
	float fSlow12 = (float)dsp->fHslider6;
	float fSlow13 = expf((dsp->fConst3 / fSlow12));
	float fSlow14 = faustpower2_f(fSlow13);
	float fSlow15 = cosf((dsp->fConst4 * (float)dsp->fHslider7));
	float fSlow16 = (1.f - (fSlow14 * fSlow15));
	float fSlow17 = (1.f - fSlow14);
	float fSlow18 = (fSlow16 / fSlow17);
	float fSlow19 = sqrtf(max(0.f, ((faustpower2_f(fSlow16) / faustpower2_f(fSlow17)) - 1.f)));
	float fSlow20 = (fSlow18 - fSlow19);
	float fSlow21 = (((1.f + fSlow19) - fSlow18) * fSlow13);
	float fSlow22 = (float)dsp->fHslider8;
	float fSlow23 = ((expf((dsp->fConst3 / fSlow22)) / fSlow13) - 1.f);
	float fSlow24 = (1.f / tanf((dsp->fConst5 * (float)dsp->fHslider9)));
	float fSlow25 = (1.f + fSlow24);
	float fSlow26 = (0.f - ((1.f - fSlow24) / fSlow25));
	float fSlow27 = (1.f / fSlow25);
	int iSlow28 = (int)((int)(dsp->fConst8 * (float)dsp->fHslider10) & 8191);
	float fSlow29 = expf((dsp->fConst11 / fSlow12));
	float fSlow30 = faustpower2_f(fSlow29);
	float fSlow31 = (1.f - (fSlow15 * fSlow30));
	float fSlow32 = (1.f - fSlow30);
	float fSlow33 = (fSlow31 / fSlow32);
	float fSlow34 = sqrtf(max(0.f, ((faustpower2_f(fSlow31) / faustpower2_f(fSlow32)) - 1.f)));
	float fSlow35 = (fSlow33 - fSlow34);
	float fSlow36 = (((1.f + fSlow34) - fSlow33) * fSlow29);
	float fSlow37 = ((expf((dsp->fConst11 / fSlow22)) / fSlow29) - 1.f);
	float fSlow38 = expf((dsp->fConst16 / fSlow12));
	float fSlow39 = faustpower2_f(fSlow38);
	float fSlow40 = (1.f - (fSlow15 * fSlow39));
	float fSlow41 = (1.f - fSlow39);
	float fSlow42 = (fSlow40 / fSlow41);
	float fSlow43 = sqrtf(max(0.f, ((faustpower2_f(fSlow40) / faustpower2_f(fSlow41)) - 1.f)));
	float fSlow44 = (fSlow42 - fSlow43);
	float fSlow45 = (((1.f + fSlow43) - fSlow42) * fSlow38);
	float fSlow46 = ((expf((dsp->fConst16 / fSlow22)) / fSlow38) - 1.f);
	float fSlow47 = expf((dsp->fConst21 / fSlow12));
	float fSlow48 = faustpower2_f(fSlow47);
	float fSlow49 = (1.f - (fSlow15 * fSlow48));
	float fSlow50 = (1.f - fSlow48);
	float fSlow51 = (fSlow49 / fSlow50);
	float fSlow52 = sqrtf(max(0.f, ((faustpower2_f(fSlow49) / faustpower2_f(fSlow50)) - 1.f)));
	float fSlow53 = (fSlow51 - fSlow52);
	float fSlow54 = (((1.f + fSlow52) - fSlow51) * fSlow47);
	float fSlow55 = ((expf((dsp->fConst21 / fSlow22)) / fSlow47) - 1.f);
	float fSlow56 = expf((dsp->fConst26 / fSlow12));
	float fSlow57 = faustpower2_f(fSlow56);
	float fSlow58 = (1.f - (fSlow15 * fSlow57));
	float fSlow59 = (1.f - fSlow57);
	float fSlow60 = (fSlow58 / fSlow59);
	float fSlow61 = sqrtf(max(0.f, ((faustpower2_f(fSlow58) / faustpower2_f(fSlow59)) - 1.f)));
	float fSlow62 = (fSlow60 - fSlow61);
	float fSlow63 = (((1.f + fSlow61) - fSlow60) * fSlow56);
	float fSlow64 = ((expf((dsp->fConst26 / fSlow22)) / fSlow56) - 1.f);
	float fSlow65 = expf((dsp->fConst31 / fSlow12));
	float fSlow66 = faustpower2_f(fSlow65);
	float fSlow67 = (1.f - (fSlow15 * fSlow66));
	float fSlow68 = (1.f - fSlow66);
	float fSlow69 = (fSlow67 / fSlow68);
	float fSlow70 = sqrtf(max(0.f, ((faustpower2_f(fSlow67) / faustpower2_f(fSlow68)) - 1.f)));
	float fSlow71 = (fSlow69 - fSlow70);
	float fSlow72 = (((1.f + fSlow70) - fSlow69) * fSlow65);
	float fSlow73 = ((expf((dsp->fConst31 / fSlow22)) / fSlow65) - 1.f);
	float fSlow74 = expf((dsp->fConst36 / fSlow12));
	float fSlow75 = faustpower2_f(fSlow74);
	float fSlow76 = (1.f - (fSlow15 * fSlow75));
	float fSlow77 = (1.f - fSlow75);
	float fSlow78 = (fSlow76 / fSlow77);
	float fSlow79 = sqrtf(max(0.f, ((faustpower2_f(fSlow76) / faustpower2_f(fSlow77)) - 1.f)));
	float fSlow80 = (fSlow78 - fSlow79);
	float fSlow81 = (((1.f + fSlow79) - fSlow78) * fSlow74);
	float fSlow82 = ((expf((dsp->fConst36 / fSlow22)) / fSlow74) - 1.f);
	float fSlow83 = expf((dsp->fConst41 / fSlow12));
	float fSlow84 = faustpower2_f(fSlow83);
	float fSlow85 = (1.f - (fSlow15 * fSlow84));
	float fSlow86 = (1.f - fSlow84);
	float fSlow87 = (fSlow85 / fSlow86);
	float fSlow88 = sqrtf(max(0.f, ((faustpower2_f(fSlow85) / faustpower2_f(fSlow86)) - 1.f)));
	float fSlow89 = (fSlow87 - fSlow88);
	float fSlow90 = (((1.f + fSlow88) - fSlow87) * fSlow83);
	float fSlow91 = ((expf((dsp->fConst41 / fSlow22)) / fSlow83) - 1.f);
	/* C99 loop */
	{
		int i;
		for (i = 0; (i < count); i = (i + 1)) {
			dsp->fRec0[0] = ((0.999f * dsp->fRec0[1]) + fSlow0);
			dsp->fRec1[0] = ((0.999f * dsp->fRec1[1]) + fSlow1);
			float fTemp0 = (1.f - dsp->fRec1[0]);
			float fTemp1 = (float)input0[i];
			dsp->fVec0[(dsp->IOTA & 8191)] = fTemp1;
			float fTemp2 = (fSlow6 * dsp->fRec2[1]);
			float fTemp3 = (fSlow11 * dsp->fRec3[1]);
			dsp->fRec15[0] = ((fSlow26 * dsp->fRec15[1]) + (fSlow27 * (dsp->fRec11[1] + dsp->fRec11[2])));
			dsp->fRec14[0] = ((fSlow20 * dsp->fRec14[1]) + (fSlow21 * (dsp->fRec11[1] + (fSlow23 * dsp->fRec15[0]))));
			dsp->fVec1[(dsp->IOTA & 32767)] = ((0.353553f * dsp->fRec14[0]) + 1e-20f);
			float fTemp4 = (float)input1[i];
			dsp->fVec2[(dsp->IOTA & 8191)] = fTemp4;
			float fTemp5 = (0.3f * dsp->fVec2[((dsp->IOTA - iSlow28) & 8191)]);
			float fTemp6 = (((0.6f * dsp->fRec12[1]) + dsp->fVec1[((dsp->IOTA - dsp->iConst7) & 32767)]) - fTemp5);
			dsp->fVec3[(dsp->IOTA & 2047)] = fTemp6;
			dsp->fRec12[0] = dsp->fVec3[((dsp->IOTA - dsp->iConst9) & 2047)];
			float fRec13 = (0.f - (0.6f * fTemp6));
			dsp->fRec19[0] = ((fSlow26 * dsp->fRec19[1]) + (fSlow27 * (dsp->fRec7[1] + dsp->fRec7[2])));
			dsp->fRec18[0] = ((fSlow35 * dsp->fRec18[1]) + (fSlow36 * (dsp->fRec7[1] + (fSlow37 * dsp->fRec19[0]))));
			dsp->fVec4[(dsp->IOTA & 32767)] = ((0.353553f * dsp->fRec18[0]) + 1e-20f);
			float fTemp7 = (((0.6f * dsp->fRec16[1]) + dsp->fVec4[((dsp->IOTA - dsp->iConst13) & 32767)]) - fTemp5);
			dsp->fVec5[(dsp->IOTA & 4095)] = fTemp7;
			dsp->fRec16[0] = dsp->fVec5[((dsp->IOTA - dsp->iConst14) & 4095)];
			float fRec17 = (0.f - (0.6f * fTemp7));
			dsp->fRec23[0] = ((fSlow26 * dsp->fRec23[1]) + (fSlow27 * (dsp->fRec9[1] + dsp->fRec9[2])));
			dsp->fRec22[0] = ((fSlow44 * dsp->fRec22[1]) + (fSlow45 * (dsp->fRec9[1] + (fSlow46 * dsp->fRec23[0]))));
			dsp->fVec6[(dsp->IOTA & 16383)] = ((0.353553f * dsp->fRec22[0]) + 1e-20f);
			float fTemp8 = (dsp->fVec6[((dsp->IOTA - dsp->iConst18) & 16383)] + (fTemp5 + (0.6f * dsp->fRec20[1])));
			dsp->fVec7[(dsp->IOTA & 4095)] = fTemp8;
			dsp->fRec20[0] = dsp->fVec7[((dsp->IOTA - dsp->iConst19) & 4095)];
			float fRec21 = (0.f - (0.6f * fTemp8));
			dsp->fRec27[0] = ((fSlow26 * dsp->fRec27[1]) + (fSlow27 * (dsp->fRec5[1] + dsp->fRec5[2])));
			dsp->fRec26[0] = ((fSlow53 * dsp->fRec26[1]) + (fSlow54 * (dsp->fRec5[1] + (fSlow55 * dsp->fRec27[0]))));
			dsp->fVec8[(dsp->IOTA & 32767)] = ((0.353553f * dsp->fRec26[0]) + 1e-20f);
			float fTemp9 = (fTemp5 + ((0.6f * dsp->fRec24[1]) + dsp->fVec8[((dsp->IOTA - dsp->iConst23) & 32767)]));
			dsp->fVec9[(dsp->IOTA & 4095)] = fTemp9;
			dsp->fRec24[0] = dsp->fVec9[((dsp->IOTA - dsp->iConst24) & 4095)];
			float fRec25 = (0.f - (0.6f * fTemp9));
			dsp->fRec31[0] = ((fSlow26 * dsp->fRec31[1]) + (fSlow27 * (dsp->fRec10[1] + dsp->fRec10[2])));
			dsp->fRec30[0] = ((fSlow62 * dsp->fRec30[1]) + (fSlow63 * (dsp->fRec10[1] + (fSlow64 * dsp->fRec31[0]))));
			dsp->fVec10[(dsp->IOTA & 16383)] = ((0.353553f * dsp->fRec30[0]) + 1e-20f);
			float fTemp10 = (0.3f * dsp->fVec0[((dsp->IOTA - iSlow28) & 8191)]);
			float fTemp11 = (dsp->fVec10[((dsp->IOTA - dsp->iConst28) & 16383)] - (fTemp10 + (0.6f * dsp->fRec28[1])));
			dsp->fVec11[(dsp->IOTA & 2047)] = fTemp11;
			dsp->fRec28[0] = dsp->fVec11[((dsp->IOTA - dsp->iConst29) & 2047)];
			float fRec29 = (0.6f * fTemp11);
			dsp->fRec35[0] = ((fSlow26 * dsp->fRec35[1]) + (fSlow27 * (dsp->fRec6[1] + dsp->fRec6[2])));
			dsp->fRec34[0] = ((fSlow71 * dsp->fRec34[1]) + (fSlow72 * (dsp->fRec6[1] + (fSlow73 * dsp->fRec35[0]))));
			dsp->fVec12[(dsp->IOTA & 16383)] = ((0.353553f * dsp->fRec34[0]) + 1e-20f);
			float fTemp12 = (dsp->fVec12[((dsp->IOTA - dsp->iConst33) & 16383)] - (fTemp10 + (0.6f * dsp->fRec32[1])));
			dsp->fVec13[(dsp->IOTA & 4095)] = fTemp12;
			dsp->fRec32[0] = dsp->fVec13[((dsp->IOTA - dsp->iConst34) & 4095)];
			float fRec33 = (0.6f * fTemp12);
			dsp->fRec39[0] = ((fSlow26 * dsp->fRec39[1]) + (fSlow27 * (dsp->fRec8[1] + dsp->fRec8[2])));
			dsp->fRec38[0] = ((fSlow80 * dsp->fRec38[1]) + (fSlow81 * (dsp->fRec8[1] + (fSlow82 * dsp->fRec39[0]))));
			dsp->fVec14[(dsp->IOTA & 16383)] = ((0.353553f * dsp->fRec38[0]) + 1e-20f);
			float fTemp13 = ((fTemp10 + dsp->fVec14[((dsp->IOTA - dsp->iConst38) & 16383)]) - (0.6f * dsp->fRec36[1]));
			dsp->fVec15[(dsp->IOTA & 4095)] = fTemp13;
			dsp->fRec36[0] = dsp->fVec15[((dsp->IOTA - dsp->iConst39) & 4095)];
			float fRec37 = (0.6f * fTemp13);
			dsp->fRec43[0] = ((fSlow26 * dsp->fRec43[1]) + (fSlow27 * (dsp->fRec4[1] + dsp->fRec4[2])));
			dsp->fRec42[0] = ((fSlow89 * dsp->fRec42[1]) + (fSlow90 * (dsp->fRec4[1] + (fSlow91 * dsp->fRec43[0]))));
			dsp->fVec16[(dsp->IOTA & 16383)] = ((0.353553f * dsp->fRec42[0]) + 1e-20f);
			float fTemp14 = ((dsp->fVec16[((dsp->IOTA - dsp->iConst43) & 16383)] + fTemp10) - (0.6f * dsp->fRec40[1]));
			dsp->fVec17[(dsp->IOTA & 2047)] = fTemp14;
			dsp->fRec40[0] = dsp->fVec17[((dsp->IOTA - dsp->iConst44) & 2047)];
			float fRec41 = (0.6f * fTemp14);
			float fTemp15 = (fRec41 + fRec37);
			float fTemp16 = (fRec29 + (fRec33 + fTemp15));
			dsp->fRec4[0] = (dsp->fRec12[1] + (dsp->fRec16[1] + (dsp->fRec20[1] + (dsp->fRec24[1] + (dsp->fRec28[1] + (dsp->fRec32[1] + (dsp->fRec36[1] + (dsp->fRec40[1] + (fRec13 + (fRec17 + (fRec21 + (fRec25 + fTemp16))))))))))));
			dsp->fRec5[0] = (0.f - ((dsp->fRec12[1] + (dsp->fRec16[1] + (dsp->fRec20[1] + (dsp->fRec24[1] + (fRec13 + (fRec17 + (fRec25 + fRec21))))))) - (dsp->fRec28[1] + (dsp->fRec32[1] + (dsp->fRec36[1] + (dsp->fRec40[1] + fTemp16))))));
			float fTemp17 = (fRec33 + fRec29);
			dsp->fRec6[0] = (0.f - ((dsp->fRec12[1] + (dsp->fRec16[1] + (dsp->fRec28[1] + (dsp->fRec32[1] + (fRec13 + (fRec17 + fTemp17)))))) - (dsp->fRec20[1] + (dsp->fRec24[1] + (dsp->fRec36[1] + (dsp->fRec40[1] + (fRec21 + (fRec25 + fTemp15))))))));
			dsp->fRec7[0] = (0.f - ((dsp->fRec20[1] + (dsp->fRec24[1] + (dsp->fRec28[1] + (dsp->fRec32[1] + (fRec21 + (fRec25 + fTemp17)))))) - (dsp->fRec12[1] + (dsp->fRec16[1] + (dsp->fRec36[1] + (dsp->fRec40[1] + (fRec13 + (fRec17 + fTemp15))))))));
			float fTemp18 = (fRec37 + fRec29);
			float fTemp19 = (fRec41 + fRec33);
			dsp->fRec8[0] = (0.f - ((dsp->fRec12[1] + (dsp->fRec20[1] + (dsp->fRec28[1] + (dsp->fRec36[1] + (fRec13 + (fRec21 + fTemp18)))))) - (dsp->fRec16[1] + (dsp->fRec24[1] + (dsp->fRec32[1] + (dsp->fRec40[1] + (fRec17 + (fRec25 + fTemp19))))))));
			dsp->fRec9[0] = (0.f - ((dsp->fRec16[1] + (dsp->fRec24[1] + (dsp->fRec28[1] + (dsp->fRec36[1] + (fRec17 + (fRec25 + fTemp18)))))) - (dsp->fRec12[1] + (dsp->fRec20[1] + (dsp->fRec32[1] + (dsp->fRec40[1] + (fRec13 + (fRec21 + fTemp19))))))));
			float fTemp20 = (fRec37 + fRec33);
			float fTemp21 = (fRec41 + fRec29);
			dsp->fRec10[0] = (0.f - ((dsp->fRec16[1] + (dsp->fRec20[1] + (dsp->fRec32[1] + (dsp->fRec36[1] + (fRec17 + (fRec21 + fTemp20)))))) - (dsp->fRec12[1] + (dsp->fRec24[1] + (dsp->fRec28[1] + (dsp->fRec40[1] + (fRec13 + (fRec25 + fTemp21))))))));
			dsp->fRec11[0] = (0.f - ((dsp->fRec12[1] + (dsp->fRec24[1] + (dsp->fRec32[1] + (dsp->fRec36[1] + (fRec13 + (fRec25 + fTemp20)))))) - (dsp->fRec16[1] + (dsp->fRec20[1] + (dsp->fRec28[1] + (dsp->fRec40[1] + (fRec17 + (fRec21 + fTemp21))))))));
			float fTemp22 = (0.37f * (dsp->fRec5[0] + dsp->fRec6[0]));
			dsp->fRec3[0] = (0.f - ((fTemp3 + (fSlow10 * dsp->fRec3[2])) - fTemp22));
			float fTemp23 = (fSlow10 * dsp->fRec3[0]);
			float fTemp24 = (0.5f * ((fTemp23 + (dsp->fRec3[2] + (fTemp22 + fTemp3))) + (fSlow8 * ((fTemp23 + (fTemp3 + dsp->fRec3[2])) - fTemp22))));
			dsp->fRec2[0] = (0.f - ((fTemp2 + (fSlow5 * dsp->fRec2[2])) - fTemp24));
			float fTemp25 = (fSlow5 * dsp->fRec2[0]);
			output0[i] = (FAUSTFLOAT)(dsp->fRec0[0] * ((fTemp0 * fTemp1) + (0.5f * (dsp->fRec1[0] * ((fTemp25 + (dsp->fRec2[2] + (fTemp24 + fTemp2))) + (fSlow3 * ((fTemp25 + (fTemp2 + dsp->fRec2[2])) - fTemp24)))))));
			float fTemp26 = (fSlow6 * dsp->fRec44[1]);
			float fTemp27 = (fSlow11 * dsp->fRec45[1]);
			float fTemp28 = (0.37f * (dsp->fRec5[0] - dsp->fRec6[0]));
			dsp->fRec45[0] = (0.f - ((fTemp27 + (fSlow10 * dsp->fRec45[2])) - fTemp28));
			float fTemp29 = (fSlow10 * dsp->fRec45[0]);
			float fTemp30 = (0.5f * ((fTemp29 + (dsp->fRec45[2] + (fTemp28 + fTemp27))) + (fSlow8 * ((fTemp29 + (fTemp27 + dsp->fRec45[2])) - fTemp28))));
			dsp->fRec44[0] = (0.f - ((fTemp26 + (fSlow5 * dsp->fRec44[2])) - fTemp30));
			float fTemp31 = (fSlow5 * dsp->fRec44[0]);
			output1[i] = (FAUSTFLOAT)(dsp->fRec0[0] * ((fTemp0 * fTemp4) + (0.5f * (dsp->fRec1[0] * ((fTemp31 + (dsp->fRec44[2] + (fTemp30 + fTemp26))) + (fSlow3 * ((fTemp31 + (fTemp26 + dsp->fRec44[2])) - fTemp30)))))));
			dsp->fRec0[1] = dsp->fRec0[0];
			dsp->fRec1[1] = dsp->fRec1[0];
			dsp->IOTA = (dsp->IOTA + 1);
			dsp->fRec15[1] = dsp->fRec15[0];
			dsp->fRec14[1] = dsp->fRec14[0];
			dsp->fRec12[1] = dsp->fRec12[0];
			dsp->fRec19[1] = dsp->fRec19[0];
			dsp->fRec18[1] = dsp->fRec18[0];
			dsp->fRec16[1] = dsp->fRec16[0];
			dsp->fRec23[1] = dsp->fRec23[0];
			dsp->fRec22[1] = dsp->fRec22[0];
			dsp->fRec20[1] = dsp->fRec20[0];
			dsp->fRec27[1] = dsp->fRec27[0];
			dsp->fRec26[1] = dsp->fRec26[0];
			dsp->fRec24[1] = dsp->fRec24[0];
			dsp->fRec31[1] = dsp->fRec31[0];
			dsp->fRec30[1] = dsp->fRec30[0];
			dsp->fRec28[1] = dsp->fRec28[0];
			dsp->fRec35[1] = dsp->fRec35[0];
			dsp->fRec34[1] = dsp->fRec34[0];
			dsp->fRec32[1] = dsp->fRec32[0];
			dsp->fRec39[1] = dsp->fRec39[0];
			dsp->fRec38[1] = dsp->fRec38[0];
			dsp->fRec36[1] = dsp->fRec36[0];
			dsp->fRec43[1] = dsp->fRec43[0];
			dsp->fRec42[1] = dsp->fRec42[0];
			dsp->fRec40[1] = dsp->fRec40[0];
			dsp->fRec4[2] = dsp->fRec4[1];
			dsp->fRec4[1] = dsp->fRec4[0];
			dsp->fRec5[2] = dsp->fRec5[1];
			dsp->fRec5[1] = dsp->fRec5[0];
			dsp->fRec6[2] = dsp->fRec6[1];
			dsp->fRec6[1] = dsp->fRec6[0];
			dsp->fRec7[2] = dsp->fRec7[1];
			dsp->fRec7[1] = dsp->fRec7[0];
			dsp->fRec8[2] = dsp->fRec8[1];
			dsp->fRec8[1] = dsp->fRec8[0];
			dsp->fRec9[2] = dsp->fRec9[1];
			dsp->fRec9[1] = dsp->fRec9[0];
			dsp->fRec10[2] = dsp->fRec10[1];
			dsp->fRec10[1] = dsp->fRec10[0];
			dsp->fRec11[2] = dsp->fRec11[1];
			dsp->fRec11[1] = dsp->fRec11[0];
			dsp->fRec3[2] = dsp->fRec3[1];
			dsp->fRec3[1] = dsp->fRec3[0];
			dsp->fRec2[2] = dsp->fRec2[1];
			dsp->fRec2[1] = dsp->fRec2[0];
			dsp->fRec45[2] = dsp->fRec45[1];
			dsp->fRec45[1] = dsp->fRec45[0];
			dsp->fRec44[2] = dsp->fRec44[1];
			dsp->fRec44[1] = dsp->fRec44[0];
			
		}
		
	}
	
}

int sp_zitarev_create(sp_zitarev **p)
{
    *p = malloc(sizeof(sp_zitarev));
    return SP_OK;
}

int sp_zitarev_destroy(sp_zitarev **p)
{
    sp_zitarev *pp = *p;
    zitarev *dsp = pp->faust;
    deletezitarev (dsp);
    free(*p);
    return SP_OK;
}

int sp_zitarev_init(sp_data *sp, sp_zitarev *p)
{
    zitarev *dsp = newzitarev(); 
    UIGlue UI;
    p->argpos = 0;
    UI.addHorizontalSlider= addHorizontalSlider;
    UI.uiInterface = p;
    buildUserInterfacezitarev(dsp, &UI);
    initzitarev(dsp, sp->sr);

    p->in_delay = p->args[0]; 
    p->lf_x = p->args[1]; 
    p->rt60_low = p->args[2]; 
    p->rt60_mid = p->args[3]; 
    p->hf_damping = p->args[4]; 
    p->eq1_freq = p->args[5]; 
    p->eq1_level = p->args[6]; 
    p->eq2_freq = p->args[7]; 
    p->eq2_level = p->args[8]; 
    p->mix = p->args[9]; 
    p->level = p->args[10];

    p->faust = dsp;
    return SP_OK;
}

int sp_zitarev_compute(sp_data *sp, sp_zitarev *p, SPFLOAT *in1, SPFLOAT *in2, SPFLOAT *out1, SPFLOAT *out2) 
{

    zitarev *dsp = p->faust;
    SPFLOAT *faust_out[] = {out1, out2};
    SPFLOAT *faust_in[] = {in1, in2};
    computezitarev(dsp, 1, faust_in, faust_out);
    return SP_OK;
}