/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    The implementations of the middle product over Z/nZ, side by side.

    For `1 <= fn <= gn` and `outlen = gn - fn + 1`, the middle product of
    `(f, fn)` and `(g, gn)` is the range `[fn - 1, gn)` of `f*g`, the
    `outlen` coefficients that are sums of the full `fn` terms. Every
    timed function computes exactly that, so they are directly
    comparable; they are all called with the longer operand first, as
    the fft_small entry points document (their heuristics read the
    second length only, so getting this wrong sends one of them down the
    single-prime direct-fft branch while another does not, which is
    worth tens of per cent and swamps what is being measured).

      #0 mulmid       _nmod_poly_mulmid, the dispatcher: what a caller
                      actually gets today;
      #1 classical    _nmod_poly_mulmid_classical, the O(fn*outlen)
                      dot-product loop;
      #2 KS           _nmod_poly_mulmid_KS, Kronecker substitution;

    and, where fft_small is available,

      #3 fft_small    _nmod_poly_mulmid_fft_small, i.e. the 16-bit
                      repacking when it applies and #4 otherwise;
      #4 window       _nmod_poly_mul_mid_default_mpn_ctx, the window of
                      the full fn by gn product;
      #5 mul          the same entry point as #4 but asked for the full
                      product of lengths fn and outlen: not a middle
                      product, but the cost floor.

    The ratio column `win/mul` = #4/#5 is what the middle product
    currently costs above that floor. The floor is a real one: a middle
    product of this shape is a bilinear map of the same size as a
    product of lengths fn and outlen, and by the transposition principle
    it can be computed by the transpose of that multiplication, whereas
    taking a window of the full fn by gn product convolves to length
    `2*fn + outlen - 2` instead of `fn + outlen - 1`.

    Every implementation is checked against #4 on each shape (against #2
    without fft_small) before anything is timed, and #1 and #2 drop out
    of the sweep once they are hopeless, so that the large shapes are
    not spent re-establishing that dot products are quadratic.

    The sweep walks fn over a geometric grid of ratio 1.5, so that both
    the `2^d` and the `1.5*2^d` transform lengths are sampled, and at
    each fn takes outlen = r*fn for several r on both sides of 1. The
    balanced shapes r = 1 are the ones Newton iteration and power series
    division spend their time in; the unbalanced ones matter because the
    window implementation convolves to length `2*fn + outlen - 2`, which
    is not symmetric in (fn, outlen) while the middle product itself is.

    Usage:
        p-mulmid
            this help, with the function numbering;
        p-mulmid -s [nbits [minfn [maxfn]]]
            the sweep, as a table (defaults: 60 bits, 16, 65536);
        p-mulmid nbits fun fn outlen
            a single measurement, for scripting one point at a time.

    nbits is the bit length of the modulus, which is taken to be the
    first prime at least 2^(nbits-1); nbits = 0 selects one of the
    fft_small context's own fft primes, for which the plan needs a
    single prime and no CRT.
*/

#include <stdlib.h>
#include <string.h>
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fft_small.h"
#include "profiler.h"

/* ------------------------------------------------------------------ */
/* the timed functions, under one signature                            */
/* ------------------------------------------------------------------ */

/* writes to z the outlen = gn - fn + 1 coefficients of the range
   [fn - 1, gn) of the product of (f, fn) and (g, gn), 1 <= fn <= gn */
typedef void (*mulmid_fun) (nn_ptr z, nn_srcptr f, slong fn,
                            nn_srcptr g, slong gn, nmod_t mod);

static void
mid_dispatch(nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn, nmod_t mod)
{
    _nmod_poly_mulmid(z, g, gn, f, fn, fn - 1, gn, mod);
}

static void
mid_classical(nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn, nmod_t mod)
{
    _nmod_poly_mulmid_classical(z, g, gn, f, fn, fn - 1, gn, mod);
}

static void
mid_KS(nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn, nmod_t mod)
{
    _nmod_poly_mulmid_KS(z, g, gn, f, fn, fn - 1, gn, mod);
}

#if FLINT_HAVE_FFT_SMALL

static void
mid_fft_small(nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn, nmod_t mod)
{
    _nmod_poly_mulmid_fft_small(z, g, gn, f, fn, fn - 1, gn, mod);
}

static void
mid_window(nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn, nmod_t mod)
{
    _nmod_poly_mul_mid_default_mpn_ctx(z, fn - 1, gn, g, gn, f, fn, mod);
}

/* not a middle product: the fft_small product of lengths fn and outlen,
   as a cost floor. It is taken through the same entry point as #4 so
   that the two differ in the shape alone. Its output has length gn,
   like the others' input g, so it fits in the same buffer */
static void
mul_floor(nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn, nmod_t mod)
{
    slong outlen = gn - fn + 1;

    if (fn >= outlen)
        _nmod_poly_mul_mid_default_mpn_ctx(z, 0, gn, f, fn, g, outlen, mod);
    else
        _nmod_poly_mul_mid_default_mpn_ctx(z, 0, gn, g, outlen, f, fn, mod);
}

# define NFUNS 6
# define NMIDS 5        /* the last one is the cost floor, not a mulmid */
# define IFUN_REF 4     /* results are checked against the window */

#else

# define NFUNS 3
# define NMIDS 3
# define IFUN_REF 2     /* no fft_small: KS is the only fast reference */

#endif

static const mulmid_fun funs[NFUNS] = {
    mid_dispatch,       /* 0 */
    mid_classical,      /* 1 */
    mid_KS,             /* 2 */
#if FLINT_HAVE_FFT_SMALL
    mid_fft_small,      /* 3 */
    mid_window,         /* 4 */
    mul_floor,          /* 5 */
#endif
};

static const char * const shortname[NFUNS] = {
    "mulmid", "classical", "KS",
#if FLINT_HAVE_FFT_SMALL
    "fft_small", "window", "mul",
#endif
};

static const char * const collabel[NFUNS] = {
    "#0", "#1", "#2",
#if FLINT_HAVE_FFT_SMALL
    "#3", "#4", "#5",
#endif
};

static const char * const description[NFUNS] = {
    "#0  --> _nmod_poly_mulmid                    (the dispatcher)",
    "#1  --> _nmod_poly_mulmid_classical          (dot products)",
    "#2  --> _nmod_poly_mulmid_KS                 (Kronecker substitution)",
#if FLINT_HAVE_FFT_SMALL
    "#3  --> _nmod_poly_mulmid_fft_small          (repacking, else #4)",
    "#4  --> _nmod_poly_mul_mid_default_mpn_ctx   (window of the full product)",
    "#5  --> _nmod_poly_mul_mid_default_mpn_ctx   (full product of lengths fn, outlen)",
#endif
};

/* the classical one is quadratic; above this many coefficient products
   its timing is a foregone conclusion and only costs wall time */
#define CLASSICAL_MAX_WORK 2.0e7

/* #1 and #2 are dropped from the rest of the sweep once they are this
   much slower than the reference and slow in absolute terms */
#define SKIP_FACTOR 25.0
#define SKIP_MIN_TIME 5.0e-2

/* ------------------------------------------------------------------ */
/* timing                                                              */
/* ------------------------------------------------------------------ */

/* wall clock seconds for one call, on already-allocated operands */
static double
time_mulmid(int ifun, nn_ptr z, nn_srcptr f, slong fn, nn_srcptr g, slong gn,
            nmod_t mod)
{
    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START;
    funs[ifun](z, f, fn, g, gn, mod);
    TIMEIT_STOP_VALUES(tcpu, twall);

    return twall;
}

/* the first prime at least 2^(nbits-1), or an fft prime of the
   fft_small default context for nbits = 0 */
static void
select_modulus(nmod_t * mod, ulong nbits)
{
#if FLINT_HAVE_FFT_SMALL
    if (nbits == 0)
    {
        mpn_ctx_struct * R = get_default_mpn_ctx();

        nmod_init(mod, R->ffts[1].mod.n);
        return;
    }
#endif

    if (nbits < 2)
        nbits = 60;

    nmod_init(mod, n_nextprime(UWORD(1) << (nbits - 1), 1));
}

/* ------------------------------------------------------------------ */
/* the sweep                                                           */
/* ------------------------------------------------------------------ */

/* outlen/fn, on both sides of the balanced case */
static const double ratios[] = {0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
#define NRATIOS (sizeof(ratios)/sizeof(ratios[0]))

static void
run_shape(slong fn, slong outlen, nmod_t mod, flint_rand_t state,
          nn_ptr f, nn_ptr g, nn_ptr z, nn_ptr zref, int * enabled)
{
    slong gn = fn + outlen - 1;
    double t[NFUNS];
    int run[NFUNS];
    int j;

    _nmod_vec_randtest(f, state, fn, mod);
    _nmod_vec_randtest(g, state, gn, mod);

    for (j = 0; j < NFUNS; j++)
        run[j] = enabled[j];

    if ((double) fn * (double) outlen > CLASSICAL_MAX_WORK)
        run[1] = 0;

    /* a wrong answer computed quickly is not a data point */
    funs[IFUN_REF](zref, f, fn, g, gn, mod);

    for (j = 0; j < NMIDS; j++)
    {
        if (!run[j] || j == IFUN_REF)
            continue;

        _nmod_vec_zero(z, outlen);
        funs[j](z, f, fn, g, gn, mod);

        if (!_nmod_vec_equal(z, zref, outlen))
        {
            flint_printf("\nFAIL: %s disagrees at fn = %wd, outlen = %wd, "
                         "mod.n = %wu\n", shortname[j], fn, outlen, mod.n);
            flint_abort();
        }
    }

    for (j = 0; j < NFUNS; j++)
        t[j] = run[j] ? time_mulmid(j, z, f, fn, g, gn, mod) : 0.0;

    flint_printf("%9wd%9wd |", fn, outlen);

    for (j = 0; j < NFUNS; j++)
    {
        if (j == NMIDS)
            flint_printf(" |");

        if (run[j])
            flint_printf("%10.2e", t[j]);
        else
            flint_printf("%10s", "-");
    }

#if FLINT_HAVE_FFT_SMALL
    flint_printf(" |%9.2f", t[4]/t[5]);
#endif
    flint_printf("\n");
    fflush(stdout);

    /* the sweep grows, and so does the gap: once one of the two
       sub-quadratic-at-best implementations is hopeless it is dropped
       for good rather than spending the rest of the run confirming it */
    for (j = 1; j <= 2; j++)
        if (run[j] && t[j] > SKIP_FACTOR*t[IFUN_REF] && t[j] > SKIP_MIN_TIME)
            enabled[j] = 0;
}

static void
sweep(ulong nbits, slong minfn, slong maxfn, flint_rand_t state)
{
    nmod_t mod;
    nn_ptr f, g, z, zref;
    slong fn, maxgn;
    int enabled[NFUNS];
    int j;

    for (j = 0; j < NFUNS; j++)
        enabled[j] = 1;

    select_modulus(&mod, nbits);

    flint_printf("# middle product: the range [fn - 1, gn) of f*g, with\n"
                 "# fn = len(f), gn = len(g) = fn + outlen - 1, "
                 "mod.n = %wu (%wu bits), %wd thread(s)\n",
                 mod.n, FLINT_BIT_COUNT(mod.n), flint_get_num_threads());

    for (j = 0; j < NFUNS; j++)
        flint_printf("# %s\n", description[j]);

    flint_printf("# times are wall clock seconds for one call; "
                 "'-' is a shape not timed\n");
#if FLINT_HAVE_FFT_SMALL
    flint_printf("# win/mul = #4/#5\n");
#endif

    flint_printf("%9s%9s |", "fn", "outlen");

    for (j = 0; j < NFUNS; j++)
    {
        if (j == NMIDS)
            flint_printf(" |");

        flint_printf("%10s", collabel[j]);
    }

#if FLINT_HAVE_FFT_SMALL
    flint_printf(" |%9s", "win/mul");
#endif
    flint_printf("\n");

    maxgn = maxfn + (slong) (ratios[NRATIOS - 1]*maxfn);

    f    = _nmod_vec_init(maxfn);
    g    = _nmod_vec_init(maxgn);
    z    = _nmod_vec_init(maxgn);
    zref = _nmod_vec_init(maxgn);

    for (fn = minfn; fn <= maxfn; fn += 1 + fn/2)
    {
        for (ulong k = 0; k < NRATIOS; k++)
        {
            slong outlen = (slong) (ratios[k]*fn);

            outlen = FLINT_MAX(outlen, WORD(1));

            run_shape(fn, outlen, mod, state, f, g, z, zref, enabled);
        }
    }

    _nmod_vec_clear(f);
    _nmod_vec_clear(g);
    _nmod_vec_clear(z);
    _nmod_vec_clear(zref);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(int argc, char ** argv)
{
    flint_rand_t state;
    int j;

    flint_rand_init(state);

    if (argc == 1)
    {
        flint_printf("Usage:\n");
        flint_printf("   `%s`\n", argv[0]);
        flint_printf("      shows this help;\n");
        flint_printf("   `%s -s [nbits [minfn [maxfn]]]`\n", argv[0]);
        flint_printf("      runs the sweep, fn from minfn to maxfn by "
                     "steps of 1.5x and outlen/fn\n"
                     "      from 1/8 to 8 "
                     "(defaults: 60 bits, minfn = 16, maxfn = 65536);\n");
        flint_printf("   `%s [nbits] [fun] [fn] [outlen]`\n", argv[0]);
        flint_printf("      times one function on one shape.\n");
        flint_printf("\n");
        flint_printf("   - nbits: bit length of the modulus, which is the "
                     "first prime >= 2^(nbits-1);\n"
                     "     nbits = 0 selects an fft prime of the fft_small "
                     "default context\n"
                     "     (single-prime plan), and falls back to 60 bits "
                     "without fft_small\n");
        flint_printf("   - fun: id number of the timed function (see below)\n");
        flint_printf("   - fn, outlen: len(f) and the number of output "
                     "coefficients, so that\n"
                     "     len(g) = fn + outlen - 1\n");
        flint_printf("\nAvailable functions:\n");

        for (j = 0; j < NFUNS; j++)
            flint_printf("   %s\n", description[j]);

#if !FLINT_HAVE_FFT_SMALL
        flint_printf("\n   (built without fft_small, so the fft_small "
                     "implementations are absent)\n");
#endif

        flint_rand_clear(state);
        return 0;
    }

    /* warm up: the first fft_small call of the process builds the
       default context, which is not what is being measured */
    {
        nmod_t mod;
        nn_ptr f, g, z;

        nmod_init(&mod, n_nextprime(UWORD(1) << 59, 1));

        f = _nmod_vec_init(200);
        g = _nmod_vec_init(399);
        z = _nmod_vec_init(399);
        _nmod_vec_randtest(f, state, 200, mod);
        _nmod_vec_randtest(g, state, 399, mod);

        for (j = 0; j < NFUNS; j++)
            funs[j](z, f, 200, g, 399, mod);

        _nmod_vec_clear(f);
        _nmod_vec_clear(g);
        _nmod_vec_clear(z);
    }

    if (!strcmp(argv[1], "-s"))
    {
        ulong nbits = (argc > 2) ? (ulong) atoi(argv[2]) : 60;
        slong minfn = (argc > 3) ? (slong) atol(argv[3]) : 16;
        slong maxfn = (argc > 4) ? (slong) atol(argv[4]) : 65536;

        sweep(nbits, FLINT_MAX(minfn, WORD(1)), FLINT_MAX(maxfn, minfn), state);
    }
    else if (argc == 5)
    {
        const ulong nbits = (ulong) atoi(argv[1]);
        const int ifun = atoi(argv[2]);
        const slong fn = (slong) atol(argv[3]);
        const slong outlen = (slong) atol(argv[4]);
        const slong gn = fn + outlen - 1;
        nmod_t mod;
        nn_ptr f, g, z;

        if (ifun < 0 || ifun >= NFUNS || fn < 1 || outlen < 1)
        {
            flint_printf("bad arguments; run with no argument for help\n");
            flint_rand_clear(state);
            return 1;
        }

        select_modulus(&mod, nbits);

        f = _nmod_vec_init(fn);
        g = _nmod_vec_init(gn);
        z = _nmod_vec_init(gn);
        _nmod_vec_randtest(f, state, fn, mod);
        _nmod_vec_randtest(g, state, gn, mod);

        flint_printf("bits fun fn        outlen    \n");
        flint_printf("%-4wu %-3d %-10wd%-10wd", FLINT_BIT_COUNT(mod.n), ifun,
                     fn, outlen);
        flint_printf("%.2e", time_mulmid(ifun, z, f, fn, g, gn, mod));
        flint_printf("\n");

        _nmod_vec_clear(f);
        _nmod_vec_clear(g);
        _nmod_vec_clear(z);
    }
    else
    {
        flint_printf("bad arguments; run with no argument for help\n");
        flint_rand_clear(state);
        return 1;
    }

    flint_rand_clear(state);
    return 0;
}
