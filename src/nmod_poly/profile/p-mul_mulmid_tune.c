/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    For deciding tuning thresholds, the three following runs can be useful
    (a few minutes each):

    ./build/nmod_poly/profile/p-mul_mulmid_tune -q -n 3 -b 0,20,30,50,60,64 -f 16:1024:x1.25 -r 0.1,1,10
    ./build/nmod_poly/profile/p-mul_mulmid_tune -q -n 3 -b 0,12,16,20,24,28,31,32,36,40,44,48,50,52,56,60,62,64 -f 8:400:x1.1 -r 1
    ./build/nmod_poly/profile/p-mul_mulmid_tune -q -n 3 -b 0,16,20,24,28,31,32,40,50,60,64 -f 4:256:x1.2 -o 1024,4096
    P=./build/nmod_poly/profile/p-mul_mulmid_tune; RANDOM=2660; B=(20 30 50 60 64 0); for i in $(seq 120); do b=${B[RANDOM%6]}; fn=$((2+RANDOM%1500)); gn=$((2+RANDOM%1500)); z=$((fn+gn-1)); lo=$((RANDOM%z)); hi=$((lo+1+RANDOM%(z-lo))); $P $b $fn $gn $lo $hi | awk -v s="$b $fn $gn $lo $hi" '$1=="mulmid"||$1=="mul"{printf "%s %s", (++k==1?s" ":" "), $2} END{print ""}'; done
*/

/*
    The multiplication and middle product algorithms of nmod_poly, side by
    side, on shapes related by transposition, for tuning their dispatchers.

    For `fn >= 1` and `outlen >= 1`, let `gn = fn + outlen - 1`. The middle
    product of `(f, fn)` and `(g, gn)` is the range `[fn - 1, gn)` of `f*g`,
    the `outlen` coefficients that are sums of the full `fn` terms. It is
    the transpose of the product of lengths `fn` and `outlen`, and costs
    about the same with each algorithm that transposes (classical,
    fft_small; not KS, which computes the whole `fn` by `gn` product). So
    every shape `(fn, outlen)` is timed as a middle product and as the
    product of `(f, fn)` and `(g, outlen)`:

      middle product                  product of lengths fn and outlen

      #0 mulmid   _nmod_poly_mulmid   #4 mul      _nmod_poly_mul
      #1 cl       .._mulmid_classical #5 cl       _nmod_poly_mul_classical
      #2 KS       .._mulmid_KS        #6 KS       _nmod_poly_mul_KS
      #3 fft      .._mulmid_fft_small #7 KS2      _nmod_poly_mul_KS2
                                      #8 KS4      _nmod_poly_mul_KS4
                                      #9 fft      _nmod_poly_mullow_fft_small

    #0 and #4 are the dispatchers, the others the algorithms they choose
    from (#3 and #9 repack the coefficients of moduli up to 23 when they
    can, and otherwise call fft_small's window product; #3 and #9 are
    absent without fft_small). All are called with the longer operand
    first.

    Ratio columns:

      mid best    which of #1, #2, #3 is fastest;
      #0/best     what the choice of the mulmid dispatcher costs (1.00: the
                  right choice);
      mul best    which of #5..#9 is fastest;
      #4/best     the same for the mul dispatcher;
      #0/#4       middle product against product: ideally about 1;
      #3/#9       the same for fft_small alone;

    and, where fft_small is available, the plan fft_small uses for #3:
    `np` the number of primes (`1d` when it transforms directly modulo
    mod.n) and `ztrunc` the transform length. The cost of fft_small is
    mostly a function of these two, and ztrunc is a staircase in fn and
    outlen (a power of two below 512, a multiple of 256 above).

    Before timing, every middle product is checked against #3 (#2 without
    fft_small) and every product against #9 (#6). Classical algorithms are
    not timed above 2e7 coefficient products, and unless -a is given, the
    classical and KS algorithms of each family drop out of the rest of the
    run once they are hopeless (25x slower than the reference and above
    50 ms).

    Usage:
        p-mul_mulmid_tune [options]

        -b B1,B2,...    bit lengths of the modulus, one table each
                        (default 60). The modulus is the first prime
                        >= 2^(B-1); B = 0 selects one of the fft_small
                        context's own primes (single-prime plan, no CRT).
                        "=N" uses the modulus N itself (N >= 2).
        -f LIST         the lengths fn (default 16:4096, see below);
        -o LIST         absolute values of outlen;
        -r R1,R2,...    outlen = round(R*fn), for each real R.
                        -o and -r can be combined: each fn is timed with
                        the union of both. Neither: -r 1 (balanced).
        -F LIST         time only these functions (default: all). The
                        references are still computed for the checks.
        -t T            number of threads (default 1).
        -n K            each time is the minimum over K timings (default
                        1); useful on noisy machines.
        -q              quick: aim at 20 ms per timing rather than 100 ms.
        -a              never drop classical and KS from the run.

        A LIST is comma separated, each item one of
            a           the single value a;
            a:b         a to b by steps x -> x + 1 + x/2;
            a:b:xR      a to b by factors R > 1 (at least +1 each step);
            a:b:+S      a to b by steps of S.

    Examples:

        p-mul_mulmid_tune -b 60 -f 100:220:+8 -r 1
            the balanced crossovers at 60 bits;
        p-mul_mulmid_tune -b 20,30,50,60,64 -f 2:256:x1.25 -o 4096
            one short operand against many coefficients;
        p-mul_mulmid_tune -q -n 3 -b 0,20,30,50,60,64 -f 16:1024:x1.25 -r 0.1,1,10
            the grid used to compare machines.

    Single measurements, for scripting (nbits as in -b; each prints one
    line per function, name and seconds for one call):

        p-mul_mulmid_tune nbits fn outlen
            the shape (fn, outlen) above, all functions;
        p-mul_mulmid_tune nbits fn gn nlo nhi
            the middle products #0..#3 on a general window: coefficients
            [nlo, nhi) of the product of f of length fn and g of length
            gn, 0 <= nlo < nhi <= fn + gn - 1; the products #4..#9 on the
            full product of f and g (nlo and nhi do not apply to them).
            Results are checked first.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#if FLINT_HAVE_FFT_SMALL
# include "fft_small.h"
#endif
#include "profiler.h"

/* ------------------------------------------------------------------ */
/* the timed functions                                                 */
/* ------------------------------------------------------------------ */

/* All take the operands (a, an), (b, bn), an >= bn. The middle products
   write the coefficients [nlo, nhi) of a*b, which is exactly the
   signature of _nmod_poly_mulmid and its variants; the products write
   a*b and ignore nlo, nhi. */
typedef void (*timed_fun) (nn_ptr z, nn_srcptr a, slong an,
                           nn_srcptr b, slong bn, slong nlo, slong nhi,
                           nmod_t mod);

#define PRODUCT(name, call) \
static void \
name(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn, \
     slong nlo, slong nhi, nmod_t mod) \
{ \
    call; \
}

PRODUCT(mul, _nmod_poly_mul(z, a, an, b, bn, mod))
PRODUCT(mul_classical, _nmod_poly_mul_classical(z, a, an, b, bn, mod))
PRODUCT(mul_KS, _nmod_poly_mul_KS(z, a, an, b, bn, mod))
PRODUCT(mul_KS2, _nmod_poly_mul_KS2(z, a, an, b, bn, mod))
PRODUCT(mul_KS4, _nmod_poly_mul_KS4(z, a, an, b, bn, mod))
#if FLINT_HAVE_FFT_SMALL
PRODUCT(mul_fft_small, _nmod_poly_mullow_fft_small(z, a, an, b, bn, an + bn - 1, mod))
#endif

typedef struct
{
    timed_fun f;        /* NULL if not available */
    const char * name;
    const char * col;   /* for the "best" columns */
    int skippable;      /* classical (1) or KS (2): may drop out of a run */
}
fun_info;

#if FLINT_HAVE_FFT_SMALL
# define IF_FFT(f) (f)
#else
# define IF_FFT(f) NULL
#endif

#define NFUNS 10
#define NMID 4                  /* #0..#3 are middle products */
#define IFUN_MID 0
#define IFUN_MUL 4
#if FLINT_HAVE_FFT_SMALL
# define IFUN_MID_REF 3         /* the references for the checks */
# define IFUN_MUL_REF 9
#else
# define IFUN_MID_REF 2
# define IFUN_MUL_REF 6
#endif

static const fun_info funs[NFUNS] = {
    { _nmod_poly_mulmid,                   "mulmid",           "-",   0 },
    { _nmod_poly_mulmid_classical,         "mulmid_classical", "cl",  1 },
    { _nmod_poly_mulmid_KS,                "mulmid_KS",        "KS",  2 },
    { IF_FFT(_nmod_poly_mulmid_fft_small), "mulmid_fft_small", "fft", 0 },
    { mul,                                 "mul",              "-",   0 },
    { mul_classical,                       "mul_classical",    "cl",  1 },
    { mul_KS,                              "mul_KS",           "KS",  2 },
    { mul_KS2,                             "mul_KS2",          "KS2", 2 },
    { mul_KS4,                             "mul_KS4",          "KS4", 2 },
    { IF_FFT(mul_fft_small),               "mullow_fft_small", "fft", 0 },
};

/* the classical algorithms are quadratic; above this many coefficient
   products their timing is a foregone conclusion */
#define CLASSICAL_MAX_WORK 2.0e7

/* classical and KS are dropped from the rest of the run once they are
   this much slower than the reference and slow in absolute terms */
#define SKIP_FACTOR 25.0
#define SKIP_MIN_TIME 5.0e-2

/* ------------------------------------------------------------------ */
/* options                                                             */
/* ------------------------------------------------------------------ */

typedef struct
{
    slong * v;
    slong len;
    slong alloc;
}
slong_list_struct;

typedef struct
{
    /* moduli: either a bit length (is_n = 0; 0 for the fft prime) or an
       explicit modulus (is_n = 1) */
    ulong mods[64];
    int mod_is_n[64];
    slong nmods;

    slong_list_struct fns;
    slong_list_struct outlens;
    slong_list_struct funs;
    double ratios[64];
    slong nratios;

    slong threads;
    slong nrep;
    slong target_ms;
    int no_skip;
}
options_struct;

static void
list_push(slong_list_struct * L, slong x)
{
    if (L->len == L->alloc)
    {
        L->alloc = FLINT_MAX(16, 2 * L->alloc);
        L->v = flint_realloc(L->v, L->alloc * sizeof(slong));
    }

    L->v[L->len++] = x;
}

static int
cmp_slong(const void * a, const void * b)
{
    slong x = *(const slong *) a, y = *(const slong *) b;
    return (x > y) - (x < y);
}

/* sorts and removes duplicates */
static void
list_normalise(slong_list_struct * L)
{
    slong i, j;

    if (L->len == 0)
        return;

    qsort(L->v, L->len, sizeof(slong), cmp_slong);

    for (i = j = 1; i < L->len; i++)
        if (L->v[i] != L->v[j - 1])
            L->v[j++] = L->v[i];

    L->len = j;
}

/* Each parser reads a comma separated list; the item parsers stop at the
   comma, which the loop then checks. They return 0 on a syntax error. */

/* a LIST of integers >= vmin into L (see the usage) */
static int
parse_list(slong_list_struct * L, const char * s, slong vmin)
{
    const char * p = s;
    char * end;

    for (;;)
    {
        slong a, b, x;

        a = strtol(p, &end, 10);
        if (end == p || a < vmin)
            return 0;

        if (*end != ':')
        {
            list_push(L, a);
        }
        else
        {
            p = end + 1;
            b = strtol(p, &end, 10);
            if (end == p || b < a)
                return 0;

            if (end[0] == ':' && end[1] == '+')
            {
                slong step = strtol(end + 2, &end, 10);

                if (step < 1)
                    return 0;

                for (x = a; x <= b; x += step)
                    list_push(L, x);
            }
            else if (end[0] == ':' && end[1] == 'x')
            {
                double r = strtod(end + 2, &end), y;

                if (!(r > 1.0))
                    return 0;

                /* a*r^k rounded, with no repetition */
                for (x = a, y = a; x <= b; )
                {
                    list_push(L, x);
                    y *= r;
                    x = FLINT_MAX(x + 1, (slong) floor(y + 0.5));
                }
            }
            else
            {
                for (x = a; x <= b; x += 1 + x / 2)
                    list_push(L, x);
            }
        }

        if (*end == '\0')
            return 1;
        if (*end != ',')
            return 0;
        p = end + 1;
    }
}

static int
parse_mods(options_struct * O, const char * s)
{
    const char * p = s;
    char * end;

    for (O->nmods = 0; ; p = end + 1)
    {
        int is_n = (p[0] == '=');
        ulong x = strtoul(p + is_n, &end, 10);

        if (end == p + is_n || O->nmods >= 64
            || (is_n && x < 2) || (!is_n && x > FLINT_BITS))
            return 0;

        O->mods[O->nmods] = x;
        O->mod_is_n[O->nmods] = is_n;
        O->nmods++;

        if (*end == '\0')
            return 1;
        if (*end != ',')
            return 0;
    }
}

static int
parse_ratios(options_struct * O, const char * s)
{
    const char * p = s;
    char * end;

    for (O->nratios = 0; ; p = end + 1)
    {
        double r = strtod(p, &end);

        if (end == p || !(r > 0.0) || O->nratios >= 64)
            return 0;

        O->ratios[O->nratios++] = r;

        if (*end == '\0')
            return 1;
        if (*end != ',')
            return 0;
    }
}

/* ------------------------------------------------------------------ */
/* timing, moduli, plans                                               */
/* ------------------------------------------------------------------ */

/* Seconds for one call of function j, with the microsecond wall clock of
   profiler.h: the repetition count grows until one batch takes target_ms,
   then the minimum over nrep batches is returned. (TIMEIT_START/STOP
   would do the first part, with a fixed target of 100 ms of cpu time and
   a millisecond clock, which does not allow -q and -n.) */
static double
time_fun(int j, nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn,
         slong nlo, slong nhi, nmod_t mod, const options_struct * O)
{
    timeit_t T;
    slong target = 1000 * O->target_ms, reps, k, r;
    double best;

    for (reps = 1; ; )
    {
        timeit_start_us(T);
        for (k = 0; k < reps; k++)
            funs[j].f(z, a, an, b, bn, nlo, nhi, mod);
        timeit_stop_us(T);

        if (T->wall >= target)
            break;

        /* aim straight at the target once the batch is measurable */
        if (T->wall >= target / 20)
            reps = (slong) (reps * (1.05 * target / T->wall)) + 1;
        else
            reps *= 10;
    }

    best = (double) T->wall / reps;

    for (r = 1; r < O->nrep; r++)
    {
        timeit_start_us(T);
        for (k = 0; k < reps; k++)
            funs[j].f(z, a, an, b, bn, nlo, nhi, mod);
        timeit_stop_us(T);

        best = FLINT_MIN(best, (double) T->wall / reps);
    }

    return 1e-6 * best;
}

/* the first prime at least 2^(nbits-1), or a prime of the fft_small
   default context for nbits = 0, or n itself */
static void
select_modulus(nmod_t * mod, ulong x, int is_n)
{
    if (is_n)
        nmod_init(mod, x);
#if FLINT_HAVE_FFT_SMALL
    else if (x == 0)
        nmod_init(mod, get_default_mpn_ctx()->ffts[1].mod.n);
#endif
    else if (x == FLINT_BITS)
        nmod_init(mod, n_nextprime(UWORD_MAX - 1000, 1));
    else
        nmod_init(mod, n_nextprime(UWORD(1) << (FLINT_MAX(x, 2) - 1), 1));
}

#if FLINT_HAVE_FFT_SMALL

/* as in fft_small/nmod_poly_mul.c */
#define _len_trunc(x) \
    ((n_clog2(x) < LG_BLK_SZ) ? n_pow2(n_max((ulong) 4, n_clog2(x))) \
                              : n_round_up((x), BLK_SZ))

/* prints the plan of the window product behind #3, with the arguments of
   _nmod_poly_mul_mid_mpn_ctx for the middle product of (f, fn), (g, gn) */
static void
print_plan(slong fn, slong gn, nmod_t mod)
{
    fft_small_plan_t P;
    ulong xt = n_max(_len_trunc((ulong) gn), _len_trunc((ulong) fn));

    if (!fft_small_plan_init_nmod(P, get_default_mpn_ctx(), fn - 1, gn,
                fn + gn - 1, xt, fn, 2 * NMOD_BITS(mod), mod, fn))
    {
        flint_printf(" |  -      -");
        return;
    }

    flint_printf(" | %2wu%s %6wu", P->np, P->use_direct_fft ? "d" : " ",
                 P->ztrunc);
    fft_small_plan_clear(P);
}

#endif

/* ------------------------------------------------------------------ */
/* measuring                                                           */
/* ------------------------------------------------------------------ */

/* Checks and times the functions with run[j] set: the middle products on
   the window [nlo, nhi) of (a, an) times (b, bn), the products on
   (pa, pan) times (pb, pbn); both pairs longer first. z, zref have room
   for pan + pbn - 1 and nhi - nlo coefficients. */
static void
measure(double * t, const int * run,
        nn_srcptr a, slong an, nn_srcptr b, slong bn, slong nlo, slong nhi,
        nn_srcptr pa, slong pan, nn_srcptr pb, slong pbn,
        nn_ptr z, nn_ptr zref, nmod_t mod, const options_struct * O)
{
    int j, ref = -1;

    for (j = 0; j < NFUNS; j++)
    {
        int mid = (j < NMID);
        slong len = mid ? nhi - nlo : pan + pbn - 1;
        nn_srcptr x = mid ? a : pa, y = mid ? b : pb;
        slong xn = mid ? an : pan, yn = mid ? bn : pbn;

        t[j] = 0.0;
        if (!run[j])
            continue;

        /* a wrong answer computed quickly is not a data point */
        if (ref != (mid ? IFUN_MID_REF : IFUN_MUL_REF))
        {
            ref = mid ? IFUN_MID_REF : IFUN_MUL_REF;
            funs[ref].f(zref, x, xn, y, yn, nlo, nhi, mod);
        }

        _nmod_vec_zero(z, len);
        funs[j].f(z, x, xn, y, yn, nlo, nhi, mod);

        if (!_nmod_vec_equal(z, zref, len))
        {
            flint_printf("\nFAIL: %s disagrees with %s, mod.n = %wu, "
                         "lengths %wd, %wd, window [%wd, %wd)\n", funs[j].name,
                         funs[ref].name, mod.n, xn, yn, nlo, nhi);
            flint_abort();
        }

        t[j] = time_fun(j, z, x, xn, y, yn, nlo, nhi, mod, O);
    }
}

/* ------------------------------------------------------------------ */
/* the table                                                           */
/* ------------------------------------------------------------------ */

static void
print_header(nmod_t mod, const options_struct * O)
{
    int j;

    flint_printf("# middle product: coefficients [fn - 1, gn) of f*g, fn = len(f), "
                 "gn = len(g) = fn + outlen - 1;\n"
                 "# product: f*g, len(f) = fn, len(g) = outlen\n"
                 "# mod.n = %wu (%wu bits), %wd thread(s), times: wall clock "
                 "seconds for one call (min of %wd), '-' not timed\n",
                 mod.n, FLINT_BIT_COUNT(mod.n), flint_get_num_threads(), O->nrep);

    for (j = 0; j < NFUNS; j++)
        if (funs[j].f != NULL)
            flint_printf("# #%d --> _nmod_poly_%s%s\n", j, funs[j].name,
                         (j == IFUN_MID || j == IFUN_MUL) ? "  (dispatcher)" : "");

    flint_printf("# mid, mul: fastest of #1..#3, of #5..#9"
#if FLINT_HAVE_FFT_SMALL
                 "; np, ztrunc: fft_small plan of #3 (d = direct)"
#endif
                 "\n%6s %6s |", "fn", "outlen");

    for (j = 0; j < NFUNS; j++)
    {
        if (j == IFUN_MUL)
            flint_printf(" |");
        if (funs[j].f != NULL)
            flint_printf("       #%d", j);
    }

    flint_printf(" | %4s %7s | %4s %7s | %6s", "mid", "#0/best", "mul",
                 "#4/best", "#0/#4");
#if FLINT_HAVE_FFT_SMALL
    flint_printf(" %6s | %3s %6s", "#3/#9", "np", "ztrunc");
#endif
    flint_printf("\n");
}

static void
print_ratio(double num, double den, int width)
{
    if (num > 0.0 && den > 0.0)
        flint_printf(" %*.2f", width, num / den);
    else
        flint_printf(" %*s", width, "-");
}

/* the fastest timed function among lo..hi, printed with the dispatcher
   j0 over it */
static void
print_best(const double * t, int j0, int lo, int hi)
{
    int j, best = -1;

    for (j = lo; j <= hi; j++)
        if (t[j] > 0.0 && (best == -1 || t[j] < t[best]))
            best = j;

    flint_printf(" | %4s", best == -1 ? "-" : funs[best].col);
    print_ratio(t[j0], best == -1 ? 0.0 : t[best], 7);
}

static void
run_shape(slong fn, slong outlen, nmod_t mod, flint_rand_t state,
          nn_ptr f, nn_ptr g, nn_ptr z, nn_ptr zref,
          int * enabled, const options_struct * O)
{
    slong gn = fn + outlen - 1;
    double t[NFUNS];
    int run[NFUNS];
    int j;

    _nmod_vec_randtest(f, state, fn, mod);
    _nmod_vec_randtest(g, state, gn, mod);

    for (j = 0; j < NFUNS; j++)
        run[j] = enabled[j] && funs[j].f != NULL
                 && !(funs[j].skippable == 1 && (double) fn * outlen > CLASSICAL_MAX_WORK);

    /* the middle product of (f, fn), (g, gn); the product of (f, fn),
       (g, outlen) */
    if (fn >= outlen)
        measure(t, run, g, gn, f, fn, fn - 1, gn, f, fn, g, outlen, z, zref, mod, O);
    else
        measure(t, run, g, gn, f, fn, fn - 1, gn, g, outlen, f, fn, z, zref, mod, O);

    flint_printf("%6wd %6wd |", fn, outlen);

    for (j = 0; j < NFUNS; j++)
    {
        if (j == IFUN_MUL)
            flint_printf(" |");

        if (funs[j].f == NULL)
            continue;

        if (t[j] > 0.0)
            flint_printf(" %8.2e", t[j]);
        else
            flint_printf(" %8s", "-");
    }

    print_best(t, IFUN_MID, 1, NMID - 1);
    print_best(t, IFUN_MUL, IFUN_MUL + 1, NFUNS - 1);
    flint_printf(" |");
    print_ratio(t[IFUN_MID], t[IFUN_MUL], 6);

#if FLINT_HAVE_FFT_SMALL
    print_ratio(t[3], t[9], 6);
    print_plan(fn, gn, mod);
#endif

    flint_printf("\n");
    fflush(stdout);

    /* the run grows, and so does the gap: once classical or KS is
       hopeless for a family it is dropped for good */
    if (!O->no_skip)
    {
        for (j = 0; j < NFUNS; j++)
        {
            int ref = (j < NMID) ? IFUN_MID_REF : IFUN_MUL_REF;

            if (funs[j].skippable && t[ref] > 0.0
                    && t[j] > SKIP_FACTOR * t[ref] && t[j] > SKIP_MIN_TIME)
                enabled[j] = 0;
        }
    }
}

/* the outlen values for one fn: the absolute ones, and the ratios */
static void
outlens_for(slong_list_struct * L, slong fn, const options_struct * O)
{
    slong i;

    L->len = 0;

    for (i = 0; i < O->outlens.len; i++)
        list_push(L, O->outlens.v[i]);

    for (i = 0; i < O->nratios; i++)
        list_push(L, FLINT_MAX(WORD(1), (slong) floor(O->ratios[i] * fn + 0.5)));

    list_normalise(L);
}

static void
run_table(ulong modx, int is_n, const options_struct * O, flint_rand_t state)
{
    nmod_t mod;
    nn_ptr f, g, z, zref;
    slong maxfn = 1, maxgn = 1, i, k;
    slong_list_struct L = {NULL, 0, 0};
    int enabled[NFUNS];
    int j;

    for (j = 0; j < NFUNS; j++)
        enabled[j] = (O->funs.len == 0);
    for (i = 0; i < O->funs.len; i++)
        enabled[O->funs.v[i]] = 1;

    select_modulus(&mod, modx, is_n);
    print_header(mod, O);

    for (i = 0; i < O->fns.len; i++)
    {
        outlens_for(&L, O->fns.v[i], O);
        maxfn = FLINT_MAX(maxfn, O->fns.v[i]);
        for (k = 0; k < L.len; k++)
            maxgn = FLINT_MAX(maxgn, O->fns.v[i] + L.v[k] - 1);
    }

    f    = _nmod_vec_init(maxfn);
    g    = _nmod_vec_init(maxgn);
    z    = _nmod_vec_init(maxgn);
    zref = _nmod_vec_init(maxgn);

    for (i = 0; i < O->fns.len; i++)
    {
        outlens_for(&L, O->fns.v[i], O);
        for (k = 0; k < L.len; k++)
            run_shape(O->fns.v[i], L.v[k], mod, state, f, g, z, zref, enabled, O);
    }

    _nmod_vec_clear(f);
    _nmod_vec_clear(g);
    _nmod_vec_clear(z);
    _nmod_vec_clear(zref);
    flint_free(L.v);
}

/* ------------------------------------------------------------------ */
/* single measurements                                                 */
/* ------------------------------------------------------------------ */

/* p-mul_mulmid_tune nbits fn outlen
   p-mul_mulmid_tune nbits fn gn nlo nhi */
static int
single(int argc, char ** argv, const options_struct * O, flint_rand_t state)
{
    options_struct M;
    nmod_t mod;
    slong fn, gn, outlen, nlo, nhi;
    nn_ptr f, g, z, zref;
    nn_srcptr a, b;
    slong an, bn;
    double t[NFUNS];
    int run[NFUNS];
    int j, window = (argc == 6);

    memset(&M, 0, sizeof(M));
    fn = atol(argv[2]);
    gn = window ? atol(argv[3]) : fn + atol(argv[3]) - 1;
    outlen = gn - fn + 1;
    nlo = window ? atol(argv[4]) : fn - 1;
    nhi = window ? atol(argv[5]) : gn;

    if (!parse_mods(&M, argv[1]) || M.nmods != 1 || fn < 1 || gn < 1
            || (!window && outlen < 1)
            || nlo < 0 || nhi <= nlo || nhi > fn + gn - 1)
    {
        flint_printf("bad arguments; run with -h for help\n");
        return 1;
    }

    select_modulus(&mod, M.mods[0], M.mod_is_n[0]);

    f = _nmod_vec_init(fn);
    g = _nmod_vec_init(gn);
    z = _nmod_vec_init(fn + gn);
    zref = _nmod_vec_init(fn + gn);
    _nmod_vec_randtest(f, state, fn, mod);
    _nmod_vec_randtest(g, state, gn, mod);

    if (fn >= gn)
        a = f, an = fn, b = g, bn = gn;
    else
        a = g, an = gn, b = f, bn = fn;

    for (j = 0; j < NFUNS; j++)
        run[j] = (funs[j].f != NULL);

    if (window)
    {
        flint_printf("# %wu bits, mod.n = %wu: len(f) = %wd, len(g) = %wd; "
                     "#0..#3: coefficients [%wd, %wd) of f*g; #4..#9: f*g\n",
                     FLINT_BIT_COUNT(mod.n), mod.n, fn, gn, nlo, nhi);
        measure(t, run, a, an, b, bn, nlo, nhi, a, an, b, bn, z, zref, mod, O);
    }
    else
    {
        flint_printf("# %wu bits, mod.n = %wu: fn = %wd, outlen = %wd; "
                     "#0..#3: middle product of lengths fn, gn = %wd; "
                     "#4..#9: product of lengths fn, outlen\n",
                     FLINT_BIT_COUNT(mod.n), mod.n, fn, outlen, gn);
        if (fn >= outlen)
            measure(t, run, g, gn, f, fn, nlo, nhi, f, fn, g, outlen, z, zref, mod, O);
        else
            measure(t, run, g, gn, f, fn, nlo, nhi, g, outlen, f, fn, z, zref, mod, O);
    }

    for (j = 0; j < NFUNS; j++)
        if (run[j])
            flint_printf("%-18s %.3e\n", funs[j].name, t[j]);

    _nmod_vec_clear(f);
    _nmod_vec_clear(g);
    _nmod_vec_clear(z);
    _nmod_vec_clear(zref);
    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

static void
usage(const char * name)
{
    int j;

    flint_printf("Usage: %s [options]\n\n", name);
    flint_printf(
"  -b B1,B2,...  modulus bit lengths, one table each (default 60); the modulus\n"
"                is the first prime >= 2^(B-1); B = 0: an fft_small prime\n"
"                (single-prime plan); =N: the modulus N itself\n"
"  -f LIST       the lengths fn = len(f) (default 16:4096)\n"
"  -o LIST       absolute values of outlen\n"
"  -r R1,R2,...  outlen = round(R*fn); -o and -r combine (union);\n"
"                neither: -r 1\n"
"  -F LIST       time only these functions (default all)\n"
"  -t T          threads (default 1)\n"
"  -n K          report the min over K timings (default 1)\n"
"  -q            quick: 20 ms per timing instead of 100 ms\n"
"  -a            never drop classical and KS from the run once hopeless\n"
"\n"
"  LIST: comma separated items a | a:b (steps x -> x+1+x/2) | a:b:xR\n"
"        (factor R) | a:b:+S (step S)\n"
"\n"
"  Each shape (fn, outlen) is timed as the middle product of lengths fn and\n"
"  gn = fn + outlen - 1 (#0..#3) and as the product of lengths fn and outlen\n"
"  (#4..#9), its transpose.\n"
"\n"
"  e.g. %s -q -n 3 -b 0,20,30,50,60,64 -f 16:1024:x1.25 -r 0.1,1,10\n"
"\n"
"  single timings, all functions, one line each:\n"
"       %s nbits fn outlen\n"
"       %s nbits fn gn nlo nhi   (#0..#3: coefficients [nlo, nhi) of f*g,\n"
"                                  len(f) = fn, len(g) = gn; #4..#9: f*g)\n"
"\nFunctions:\n", name, name, name, name);

    for (j = 0; j < NFUNS; j++)
        if (funs[j].f != NULL)
            flint_printf("   #%d  _nmod_poly_%s\n", j, funs[j].name);

    flint_printf("\n"
"  Ratio columns: mid, mul = fastest of #1..#3, of #5..#9 (what each\n"
"  dispatcher chooses from; #3 and #9 call the fft_small window product when\n"
"  they do not repack); #0/best, #4/best; #0/#4 (middle product against\n"
"  the product it is the transpose of); #3/#9.\n");

#if !FLINT_HAVE_FFT_SMALL
    flint_printf("\n  (built without fft_small: #3 and #9 are absent)\n");
#endif
}

int main(int argc, char ** argv)
{
    flint_rand_t state;
    options_struct O;
    int i, ok = 1;

    memset(&O, 0, sizeof(O));
    O.mods[0] = 60;
    O.nmods = 1;
    O.threads = 1;
    O.nrep = 1;
    O.target_ms = 100;

    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
    {
        usage(argv[0]);
        return 0;
    }

    flint_rand_init(state);

    /* warm up: the first fft_small call of the process builds the default
       context, which is not what is being measured */
    {
        options_struct W;
        double t[NFUNS];
        int run[NFUNS];
        nmod_t mod;
        nn_ptr f = _nmod_vec_init(400), z = _nmod_vec_init(800);

        memset(&W, 0, sizeof(W));
        W.nrep = 1;
        W.target_ms = 1;
        for (i = 0; i < NFUNS; i++)
            run[i] = (funs[i].f != NULL);

        nmod_init(&mod, n_nextprime(UWORD(1) << 59, 1));
        _nmod_vec_randtest(f, state, 400, mod);
        measure(t, run, f, 200, f + 200, 200, 199, 200, f, 200, f + 200, 200,
                z, z + 400, mod, &W);

        _nmod_vec_clear(f);
        _nmod_vec_clear(z);
    }

    if ((argc == 4 || argc == 6) && argv[1][0] != '-')
    {
        ok = !single(argc, argv, &O, state);
        flint_rand_clear(state);
        flint_cleanup_master();
        return !ok;
    }

    for (i = 1; i < argc && ok; i++)
    {
        const char * a = argv[i];
        const char * v = (i + 1 < argc) ? argv[i + 1] : NULL;

        if (!strcmp(a, "-q"))
            O.target_ms = 20;
        else if (!strcmp(a, "-a"))
            O.no_skip = 1;
        else if (v == NULL)
            ok = 0;
        else if (!strcmp(a, "-b"))
            ok = parse_mods(&O, v), i++;
        else if (!strcmp(a, "-f"))
            ok = parse_list(&O.fns, v, 1), i++;
        else if (!strcmp(a, "-o"))
            ok = parse_list(&O.outlens, v, 1), i++;
        else if (!strcmp(a, "-r"))
            ok = parse_ratios(&O, v), i++;
        else if (!strcmp(a, "-F"))
            ok = parse_list(&O.funs, v, 0), i++;
        else if (!strcmp(a, "-t"))
            O.threads = atol(v), ok = (O.threads >= 1), i++;
        else if (!strcmp(a, "-n"))
            O.nrep = atol(v), ok = (O.nrep >= 1), i++;
        else
            ok = 0;
    }

    for (i = 0; i < O.funs.len; i++)
        ok = ok && O.funs.v[i] < NFUNS;

    if (ok && O.fns.len == 0)
        ok = parse_list(&O.fns, "16:4096", 1);

    if (ok && O.outlens.len == 0 && O.nratios == 0)
        ok = parse_ratios(&O, "1");

    if (ok)
    {
        list_normalise(&O.fns);
        flint_set_num_threads(O.threads);

        for (i = 0; i < O.nmods; i++)
        {
            if (i > 0)
                flint_printf("\n");
            run_table(O.mods[i], O.mod_is_n[i], &O, state);
        }
    }
    else
        flint_printf("bad arguments; run with -h for help\n");

    flint_free(O.fns.v);
    flint_free(O.outlens.v);
    flint_free(O.funs.v);
    flint_rand_clear(state);
    flint_cleanup_master();
    return !ok;
}
