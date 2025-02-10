/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_modular.h>
#include <flint/acb_elliptic.h>
#include <flint/acb_mat.h>
#include <flint/acb_theta.h>
#include <flint/profiler.h>
#include <flint/thread_support.h>

/* some useful color operations */
#define CLAMP(y) FLINT_MAX(0.0, FLINT_MIN((y), 1.0))
#define BLEND(x,y) (0.5*(x) + 0.5*(y))
#define DODGE(a,b) ((a) / ((1.0 - (b)) + 1/256.0))

/* HLS algorithm from python's colorsys module */
static double
vv(double m1, double m2, double hue)
{
    hue = hue - floor(hue);

    if (hue < 1/6.)
        return m1 + (m2-m1)*hue*6.0;
    if (hue < 0.5)
        return m2;
    if (hue < 2/3.)
        return m1 + (m2-m1)*(2/3.-hue)*6.0;
    return m1;
}

static void
hls_to_rgb(double *R, double *G, double *B, double h, double l, double s)
{
    double m1, m2;

    if (s == 0.0)
    {
        *R = *G = *B = l;
        return;
    }

    if (l <= 0.5)
        m2 = l * (1.0+s);
    else
        m2 = l+s-(l*s);

    m1 = 2.0*l - m2;

    *R = vv(m1, m2, h + 1/3.);
    *G = vv(m1, m2, h);
    *B = vv(m1, m2, h - 1/3.);
}

static void
rgb_to_hls(double *H, double *L, double *S, double R, double G, double B)
{
    double h, l, s, hi, lo, d;

    hi = FLINT_MAX(FLINT_MAX(R, G), B);
    lo = FLINT_MIN(FLINT_MIN(R, G), B);

    l = 0.5 * (lo + hi);
    d = hi - lo;

    if (hi == lo)
    {
        s = 0.0;
        h = 0.0;
    }
    else
    {
        if (l <= 0.5)
            s = d / (hi + lo);
        else
            s = d / (2.0 - hi - lo);

        if (d == 0.0)
            d = 1.0;

        if (R == hi)
            h = (G - B) / d;
        else if (G == hi)
            h = (B - R) / d + 2.0;
        else
            h = (R - G) / d + 4.0;

        h = h / 6.0;
        if (h < 0.0)
            h += 1.0;
    }

    *H = h;
    *L = l;
    *S = s;
}

/* color balance algorithm from gimp */
static double balance_channel(double value, double l,
    double shadows, double midtones, double highlights)
{
    double a = 0.25, b = 0.333, scale = 0.7;

    shadows    *= CLAMP((l - b) / (-a) + 0.5) * scale;
    midtones   *= CLAMP((l - b) / ( a) + 0.5) *
                  CLAMP((l + b - 1.0) / (-a) + 0.5) * scale;
    highlights *= CLAMP((l + b - 1.0) / ( a) + 0.5) * scale;

    value += shadows;
    value += midtones;
    value += highlights;
    return CLAMP(value);
}

static void balance(double * R, double * G, double * B,
    double ra, double rb, double rc,  /* red shadows, midtones, highlights */
    double ga, double gb, double gc,  /* green */
    double ba, double bb, double bc)  /* blue */
{
    double h, l, s;
    double h2, l2, s2;

    rgb_to_hls(&h, &l, &s, *R, *G, *B);

    *R = balance_channel(*R, *R, ra, rb, rc);
    *G = balance_channel(*G, *G, ga, gb, gc);
    *B = balance_channel(*B, *B, ba, bb, bc);

    /* preserve lightness */
    rgb_to_hls(&h2, &l2, &s2, *R, *G, *B);
    hls_to_rgb(R, G, B, h2, l, s2);
}

#define PI 3.1415926535898

const double blue_orange_colors[][4] = {
  {-1.0,  0.0, 0.0, 0.0},
  {-0.95, 0.1, 0.2, 0.5},
  {-0.5,  0.0, 0.5, 1.0},
  {-0.05, 0.4, 0.8, 0.8},
  { 0.0,  1.0, 1.0, 1.0},
  { 0.05, 1.0, 0.9, 0.3},
  { 0.5,  0.9, 0.5, 0.0},
  { 0.95, 0.7, 0.1, 0.0},
  { 1.0,  0.0, 0.0, 0.0},
  { 2.0,  0.0, 0.0, 0.0},
};

void
color_function(double * R, double * G, double * B, const acb_t z, int mode)
{
    double H, L, S;
    slong prec, i;
    arb_t t, u;

    if (!acb_is_finite(z) || acb_rel_accuracy_bits(z) < 4)
    {
        *R = *G = *B = 0.5;
        return;
    }

    if (mode >= 2)
    {
        double R1, G1, B1;
        double R2, G2, B2;

        /* combine both color functions */
        color_function(&R1, &G1, &B1, z, 0);
        color_function(&R2, &G2, &B2, z, 1);

        *R = BLEND(R1, CLAMP(DODGE(R1, R2)));
        *G = BLEND(G1, CLAMP(DODGE(G1, G2)));
        *B = BLEND(B1, CLAMP(DODGE(B1, B2)));

        /* then play with the levels */
        if (mode == 3)
            balance(R, G, B, 0.0, -0.5, 0.2, 0.0, 0.0, -0.1, 0.0, -1.0, -0.2);
        else if (mode == 4)
            balance(R, G, B, 0.0, -0.5, 0.2, 0.0, 0.5, -0.1, 0.0, -0.3, -1.0);
        else if (mode == 5)
            balance(R, G, B, 0.0, -0.5, -1.0, 0.0, -0.1, -0.67, 0.0, -0.55, -0.12);
        else if (mode == 6)
            balance(R, G, B, 0.86, 0.0, 0.13, 0.57, 0.19, -0.52, 0.31, -0.30, -0.94);

        return;
    }

    arb_init(t);
    arb_init(u);

    prec = 32;

    arf_set_round(arb_midref(t), arb_midref(acb_realref(z)), prec, ARF_RND_DOWN);
    arf_set_round(arb_midref(u), arb_midref(acb_imagref(z)), prec, ARF_RND_DOWN);

    arb_atan2(t, u, t, prec);

    H = arf_get_d(arb_midref(t), ARF_RND_DOWN);

    if (mode == 0)
    {
        H = (H + PI) / (2 * PI) + 0.5;
        H = H - floor(H);

        acb_abs(t, z, prec);

        if (arf_cmpabs_2exp_si(arb_midref(t), 200) > 0)
        {
            L = 1.0;
        }
        else if (arf_cmpabs_2exp_si(arb_midref(t), -200) < 0)
        {
            L = 0.0;
        }
        else
        {
            L = arf_get_d(arb_midref(t), ARF_RND_DOWN);
            L = 1.0 - 1.0/(1.0 + pow(L, 0.2));
        }

        S = 0.8;

        hls_to_rgb(R, G, B, H, L, S);
    }
    else
    {
        H = H / PI;
        H = FLINT_MAX(FLINT_MIN(H, 1.0), -1.0);

        for (i = 1; ; i++)
        {
            if (blue_orange_colors[i][0] > H)
            {
                double a, ra, ga, ba, b, rb, gb, bb, s;

                a  = blue_orange_colors[i-1][0];
                ra = blue_orange_colors[i-1][1];
                ga = blue_orange_colors[i-1][2];
                ba = blue_orange_colors[i-1][3];
                b  = blue_orange_colors[i][0];
                rb = blue_orange_colors[i][1];
                gb = blue_orange_colors[i][2];
                bb = blue_orange_colors[i][3];

                s = (H - a) / (b - a);
                *R = ra + (rb - ra) * s;
                *G = ga + (gb - ga) * s;
                *B = ba + (bb - ba) * s;
                break;
            }
        }
    }

    arb_clear(t);
    arb_clear(u);
}

typedef void (*func_ptr)(acb_t, const acb_t, slong);

void
ai(acb_t res, const acb_t z, slong prec)
{
    acb_hypgeom_airy(res, NULL, NULL, NULL, z, prec);
}

void
bi(acb_t res, const acb_t z, slong prec)
{
    acb_hypgeom_airy(NULL, NULL, res, NULL, z, prec);
}

void
besselj(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_j(res, nu, z, prec);
    acb_clear(nu);
}

void
bessely(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_y(res, nu, z, prec);
    acb_clear(nu);
}

void
besseli(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_i(res, nu, z, prec);
    acb_clear(nu);
}

void
besselk(acb_t res, const acb_t z, slong prec)
{
    acb_t nu;
    acb_init(nu);
    acb_hypgeom_bessel_k(res, nu, z, prec);
    acb_clear(nu);
}

/* this function looks better when scaled */
void
modj(acb_t res, const acb_t z, slong prec)
{
    acb_modular_j(res, z, prec);
    acb_div_ui(res, res, 1728, prec);
}

void
modjq(acb_t res, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_log(res, z, prec);
    acb_const_pi(t, prec);
    acb_div(res, res, t, prec);
    acb_mul_2exp_si(res, res, -1);
    acb_div_onei(res, res);
    acb_modular_j(res, res, prec);
    acb_div_ui(res, res, 1728, prec);
    acb_clear(t);
}

void
modetaq(acb_t res, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_log(res, z, prec);
    acb_const_pi(t, prec);
    acb_div(res, res, t, prec);
    acb_mul_2exp_si(res, res, -1);
    acb_div_onei(res, res);
    acb_modular_eta(res, res, prec);
    acb_clear(t);
}

void
modlambdaq(acb_t res, const acb_t z, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_log(res, z, prec);
    acb_const_pi(t, prec);
    acb_div(res, res, t, prec);
    acb_mul_2exp_si(res, res, -1);
    acb_div_onei(res, res);
    acb_modular_lambda(res, res, prec);
    acb_clear(t);
}

void
ellipp(acb_t res, const acb_t z, slong prec)
{
    acb_onei(res);
    acb_elliptic_p(res, z, res, prec);
}

void
ellipzeta(acb_t res, const acb_t z, slong prec)
{
    acb_onei(res);
    acb_elliptic_zeta(res, z, res, prec);
}

void
ellipsigma(acb_t res, const acb_t z, slong prec)
{
    acb_onei(res);
    acb_elliptic_sigma(res, z, res, prec);
}

void
fresnels(acb_t res, const acb_t z, slong prec)
{
    acb_hypgeom_fresnel(res, NULL, z, 0, prec);
}

void
fresnelc(acb_t res, const acb_t z, slong prec)
{
    acb_hypgeom_fresnel(NULL, res, z, 0, prec);
}

void
riemanntheta(acb_t res, const acb_t z, slong prec)
{
    acb_mat_t tau;
    acb_ptr th, t;

    acb_mat_init(tau, 2, 2);

    acb_set_d_d(acb_mat_entry(tau, 0, 0), 1.0, 2.0);
    acb_set_d_d(acb_mat_entry(tau, 0, 1), -1.0, -1.0);
    acb_set_d_d(acb_mat_entry(tau, 1, 0), -1.0, -1.0);
    acb_set_d_d(acb_mat_entry(tau, 1, 1), 1.0, 2.0);

    t = _acb_vec_init(2);
    th = _acb_vec_init(16);

    acb_set(t, z);
    acb_div_onei(t + 1, z);
    acb_mul_2exp_si(t + 1, t + 1, -1);
    acb_theta_all(th, t, tau, 0, prec);
    acb_set(res, th + 9);
    acb_mul_2exp_si(res, res, -2);

    _acb_vec_clear(th, 16);
    _acb_vec_clear(t, 2);

    acb_mat_clear(tau);
}

typedef struct
{
    arf_ptr xa;
    arf_ptr xb;
    arf_ptr ya;
    arf_ptr yb;
    slong xnum;
    slong ynum;
    slong y;
    func_ptr func;
    unsigned char * buf;
    int color_mode;
}
work_t;

void worker(slong x, work_t * work)
{
    slong prec;
    acb_t z, w;
    arf_ptr xa, xb, ya, yb;
    slong xnum, ynum, y;
    double R, G, B;

    acb_init(z);
    acb_init(w);

    xa = work->xa;
    xb = work->xb;
    ya = work->ya;
    yb = work->yb;

    xnum = work->xnum;
    ynum = work->ynum;
    y = work->y;

    for (prec = 30; prec < 500; prec *= 2)
    {
        arf_sub(arb_midref(acb_imagref(z)), yb, ya, prec, ARF_RND_DOWN);
        arf_mul_ui(arb_midref(acb_imagref(z)),
            arb_midref(acb_imagref(z)), y, prec, ARF_RND_DOWN);
        arf_div_ui(arb_midref(acb_imagref(z)),
            arb_midref(acb_imagref(z)), ynum - 1, prec, ARF_RND_DOWN);
        arf_add(arb_midref(acb_imagref(z)),
            arb_midref(acb_imagref(z)), ya, prec, ARF_RND_DOWN);

        arf_sub(arb_midref(acb_realref(z)), xb, xa, prec, ARF_RND_DOWN);
        arf_mul_ui(arb_midref(acb_realref(z)),
            arb_midref(acb_realref(z)), x, prec, ARF_RND_DOWN);
        arf_div_ui(arb_midref(acb_realref(z)),
            arb_midref(acb_realref(z)), xnum - 1, prec, ARF_RND_DOWN);
        arf_add(arb_midref(acb_realref(z)),
            arb_midref(acb_realref(z)), xa, prec, ARF_RND_DOWN);

        work->func(w, z, prec);

        if (acb_rel_accuracy_bits(w) > 4)
            break;
    }

    color_function(&R, &G, &B, w, work->color_mode);

    work->buf[3 * (y * xnum + x) + 0] = FLINT_MIN(255, floor(R * 255));
    work->buf[3 * (y * xnum + x) + 1] = FLINT_MIN(255, floor(G * 255));
    work->buf[3 * (y * xnum + x) + 2] = FLINT_MIN(255, floor(B * 255));

    acb_clear(z);
    acb_clear(w);
}

int main(int argc, char *argv[])
{
    slong x, y, xnum, ynum, i;
    double dxa, dxb, dya, dyb;
    FILE * fp;
    arf_t xa, xb, ya, yb;
    acb_t z, w;
    func_ptr func;
    int color_mode;
    unsigned char * buf;
    slong num_threads;

    if (argc < 2)
    {
        printf("complex_plot [-range xa xb ya yb] [-size xn yn] [-color n] [-threads n] <func>\n\n");

        printf("Plots one of the predefined functions on [xa,xb] + [ya,yb]i\n");
        printf("using domain coloring, at a resolution of xn by yn pixels.\n\n");

        printf("Defaults parameters are [-10,10] + [-10,10]i and xn = yn = 512.\n\n");

        printf("A color function can be selected with -color. The choices are:\n");
        printf("0   phase=hue, magnitude=brightness\n");
        printf("1   phase only, white-gold-black-blue-white counterclockwise\n");
        printf("2   0+1 (dodge filter)\n");
        printf("3   0+1, shiny\n");
        printf("4   0+1, warm\n");
        printf("5   0+1, cold\n");
        printf("6   0+1, tomato\n\n");

        printf("The output is written to arbplot.ppm. If you have ImageMagick,\n");
        printf("run [convert arbplot.ppm arbplot.png] to get a PNG.\n\n");

        printf("Function codes <func> are:\n");
        printf("  sin        - Sine\n");
        printf("  gamma      - Gamma function\n");
        printf("  digamma    - Digamma function\n");
        printf("  lgamma     - Logarithmic gamma function\n");
        printf("  zeta       - Riemann zeta function\n");
        printf("  erf        - Error function\n");
        printf("  ai         - Airy function Ai\n");
        printf("  bi         - Airy function Bi\n");
        printf("  besselj    - Bessel function J_0\n");
        printf("  bessely    - Bessel function Y_0\n");
        printf("  besseli    - Bessel function I_0\n");
        printf("  besselk    - Bessel function K_0\n");
        printf("  modj       - Modular j-function\n");
        printf("  modjq      - Modular j-function (as function of q)\n");
        printf("  modeta     - Dedekind eta function\n");
        printf("  modetaq    - Dedekind eta function (as function of q)\n");
        printf("  modlambda  - Modular lambda function\n");
        printf("  modlambdaq - Modular lambda function (as function of q)\n");
        printf("  ellipp     - Weierstrass elliptic function (on square lattice)\n");
        printf("  ellipzeta  - Weierstrass elliptic function (on square lattice)\n");
        printf("  ellipsigma - Weierstrass elliptic function (on square lattice)\n");
        printf("  barnesg    - Barnes G-function\n");
        printf("  agm        - Arithmetic geometric mean\n");
        printf("  fresnels   - Fresnel integral S\n");
        printf("  fresnelc   - Fresnel integral C\n\n");

        return 1;
    }

    xnum = 512;
    ynum = 512;
    dxa = dya = -10;
    dxb = dyb = 10;
    func = acb_gamma;
    color_mode = 0;
    num_threads = 1;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-size"))
        {
            xnum = atol(argv[i+1]);
            ynum = atol(argv[i+2]);
            i += 2;
        }
        else if (!strcmp(argv[i], "-range"))
        {
            dxa = atof(argv[i+1]);
            dxb = atof(argv[i+2]);
            dya = atof(argv[i+3]);
            dyb = atof(argv[i+4]);
            i += 4;
        }
        else if (!strcmp(argv[i], "-color"))
        {
            color_mode = atoi(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-threads"))
        {
            num_threads = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "sin"))
            func = acb_sin;
        else if (!strcmp(argv[i], "gamma"))
            func = acb_gamma;
        else if (!strcmp(argv[i], "digamma"))
            func = acb_digamma;
        else if (!strcmp(argv[i], "lgamma"))
            func = acb_lgamma;
        else if (!strcmp(argv[i], "zeta"))
            func = acb_zeta;
        else if (!strcmp(argv[i], "erf"))
            func = acb_hypgeom_erf;
        else if (!strcmp(argv[i], "ai"))
            func = ai;
        else if (!strcmp(argv[i], "bi"))
            func = bi;
        else if (!strcmp(argv[i], "besselj"))
            func = besselj;
        else if (!strcmp(argv[i], "bessely"))
            func = bessely;
        else if (!strcmp(argv[i], "besseli"))
            func = besseli;
        else if (!strcmp(argv[i], "besselk"))
            func = besselk;
        else if (!strcmp(argv[i], "modj"))
            func = modj;
        else if (!strcmp(argv[i], "modjq"))
            func = modjq;
        else if (!strcmp(argv[i], "modeta"))
            func = acb_modular_eta;
        else if (!strcmp(argv[i], "modetaq"))
            func = modetaq;
        else if (!strcmp(argv[i], "modlambda"))
            func = acb_modular_lambda;
        else if (!strcmp(argv[i], "modlambdaq"))
            func = modlambdaq;
        else if (!strcmp(argv[i], "ellipp"))
            func = ellipp;
        else if (!strcmp(argv[i], "ellipzeta"))
            func = ellipzeta;
        else if (!strcmp(argv[i], "ellipsigma"))
            func = ellipsigma;
        else if (!strcmp(argv[i], "barnesg"))
            func = acb_barnes_g;
        else if (!strcmp(argv[i], "agm"))
            func = acb_agm1;
        else if (!strcmp(argv[i], "fresnels"))
            func = fresnels;
        else if (!strcmp(argv[i], "fresnelc"))
            func = fresnelc;
        else if (!strcmp(argv[i], "riemanntheta"))
            func = riemanntheta;
        else
        {
            printf("unknown option: %s\n", argv[i]);
            return 1;
        }
    }

    acb_init(z);
    acb_init(w);

    arf_init(xa);
    arf_init(xb);
    arf_init(ya);
    arf_init(yb);

    arf_set_d(xa, dxa);
    arf_set_d(xb, dxb);
    arf_set_d(ya, dya);
    arf_set_d(yb, dyb);

    buf = flint_malloc(3 * xnum * ynum);

    flint_set_num_threads(num_threads);

    TIMEIT_ONCE_START

    for (y = ynum - 1; y >= 0; y--)
    {
        work_t work;

        if (y % (ynum / 16) == 0)
            printf("row %ld\n", y);

        work.xa = xa;
        work.xb = xb;
        work.ya = ya;
        work.yb = yb;
        work.xnum = xnum;
        work.ynum = ynum;
        work.y = y;
        work.func = func;
        work.buf = buf;
        work.color_mode = color_mode;

        flint_parallel_do((do_func_t) worker, &work, xnum, num_threads, FLINT_PARALLEL_STRIDED);
    }

    TIMEIT_ONCE_STOP

    fp = fopen("arbplot.ppm", "w");
    fprintf(fp, "P6\n%ld %ld 255\n", xnum, ynum);

    for (y = ynum - 1; y >= 0; y--)
    {
        for (x = 0; x < xnum; x++)
        {
            fputc(buf[3 * (y * xnum + x) + 0], fp);
            fputc(buf[3 * (y * xnum + x) + 1], fp);
            fputc(buf[3 * (y * xnum + x) + 2], fp);
        }
    }

    flint_free(buf);

    fclose(fp);

    arf_clear(xa);
    arf_clear(xb);
    arf_clear(ya);
    arf_clear(yb);

    acb_clear(z);
    acb_clear(w);

    flint_cleanup_master();
    return 0;
}

