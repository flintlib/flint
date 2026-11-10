/*
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <string.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "thread_support.h"
#include "ulong_extras.h"
#include "fmpq.h"
#include "arith.h"
#include "bernoulli.h"
#include "mp_real.h"
#include "impl.h"

/* Khinchin's constant K from

       S = log K log 2 = sum_{k >= 2} f(1/k),
       f(x) = -log(1 - x) log(1 + x) = sum_{n >= 1} c_n x^(2n),
       c_n = A_n / n,  A_n = 1 - 1/2 + 1/3 - ... + 1/(2n - 1),

   the cells (k, n) of the double sum being covered in three ways.

   ROWS.  For 2 <= k < K0, f(1/k) = (l_k - l_(k-1)) (l_(k+1) - l_k) with
   l_m = log m, an integer combination of the logarithms lambda_j of the
   primes p_j <= K0 (_mp_real_log_primes_vec: one set of Machin-type
   series).  The sum of the rows is then a quadratic form
   sum_ij C_ij lambda_i lambda_j with small integer C, evaluated as
   sum_i lambda_i (sum_(j>=i) C'_ij lambda_j) with one full product per
   prime instead of one per row.

   COLUMNS.  The rest is sum_n c_n Z_n with Z_n = zeta(2n, K0)
   = sum_(k >= K0) k^(-2n), which has K0^(-2n) as its leading term, so
   the columns n <= M with K0^(-2M) ~ 2^-p are needed.  The weights are
   applied by summation by parts: with a_1 = 1,
   a_i = -1 / ((2i - 2)(2i - 1)) (so that A_n = a_1 + ... + a_n) and the
   suffix sums R_i = sum_(lo <= n <= hi, n >= i) X_n / n,

       sum_(n=lo)^hi c_n X_n = A_(lo-1) R_lo + sum_(i=lo)^hi a_i R_i,

   two single-limb divisions and two additions per column, with no
   multiplication by the real A_n.

   - Large n (n1 <= n <= M): Z_n by direct summation of k^(-2n) over
     K0 <= k <= K(n) ~ B^(W / 2n), in fixed point with W fraction limbs
     (an absolute error counted in ulps B^-W), with recycled powers.
     Only k in [K0, 2 K0) and the odd k >= 2 K0 are summed: the even
     k >= 2 K0 contribute 4^-n Z_n, so Z_n = X_n / (1 - 4^-n) for the
     partial sum X_n, by a few shifted additions.  Each power is kept as
     an integer X = k^(-2n) B^(W + e) with e extended limbs below the
     frame, created at the first column where it matters (binary
     powering and a division) and stepped down by X <- k^2 X
     (mpn_mul_1, much faster than the division by k^2 of an upward
     recurrence).  A downward step multiplies the absolute error by k^2
     but keeps the relative accuracy, so the power is created with the
     relative accuracy needed at the bottom of its range and its bottom
     limb is dropped each time the error count grows past B; the range
     is cut into up to three segments, each created afresh, so that the
     length of the multiplications tapers off with the accuracy needed.
     The column sums, and with them the suffix sums, taper off by
     themselves: k^(-2n) only adds its limbs above B^-W, and R_n is of
     the size of Z_n.  The columns are taken in batches of 16, each
     block of 256 limbs of a power being stepped through the batch while
     in cache (_khinchin_tiled_batch), since streaming every power once
     per column is bound by the memory bandwidth at high precision.

   - Small n (1 <= n < n1): sum_n c_n Z_n = sum_n c_n zeta(2n)
     - sum_n c_n (1 + sum_(k=2)^(K0-1) k^(-2n)), the second sum in fixed
     point like the direct columns (one unit limb), the first from
     zeta(2n) = |B_2n| t_n / 2 with the exact Bernoulli numbers of the
     reverse iterator (bernoulli_rev: sums of odd powers, recycled
     downward, at the precision of the Bernoulli numbers, which is
     below the working precision for these n) and
     t_n = (2 pi)^(2n) / (2n)!.  The columns are taken in blocks of b:
     within a block [lo, hi] everything is computed in units of t_hi,
     where zeta(2m) / t_hi = (|B_2m| c) (2 pi)^(-2j) / 2 with
     j = hi - m and the integer c = (2m + 1) ... (2 hi), so with the
     b + 1 powers (2 pi)^(-2j) computed once a column costs one product
     of the Bernoulli number by a full-precision power instead of two
     full products, and each block four full products.

   The direct method costs about (K(n) - K0) (W - 2n log_B K0) limb
   operations in column n and the Bernoulli column a product of the
   size of the Bernoulli number by a full-precision number plus the
   iterator, so the switch n1 is where K(n1) is a few hundred to a few
   thousand (_khinchin_params).  Finally K = exp(S / log 2).

   THREADS.  The sums being linear, the direct columns split into groups
   of k with separate summations by parts, and the Bernoulli columns
   into ranges of n with separate iterators (chunk c of [lo, hi]
   contributing res_c + A_(lo-1) R_lo); these, the subtracted powers
   and the rows run as independent tasks. */

/* ---- parameters ---- */

/* 2 pi e: |B_2m| ~ 4 sqrt(pi m) (m / (pi e))^(2m), so B_2m has about
   2m log2(2m / (2 pi e)) bits */
#define KHINCHIN_2PIE 17.079468445347132

static double
_log2_d(double x)
{
    return log(x) * 1.4426950408889634074;
}

/* the largest n with k^(-2n) >= 2^-F, up to the rounding of doubles
   (which may err either way; the truncation bounds account for what is
   left out) */
static slong
_khinchin_n_hi(ulong k, double F)
{
    return (slong) (F / (2.0 * _log2_d((double) k)));
}

static void
_khinchin_params(slong * K0, slong * n1, slong W)
{
    double x, p = (double) FLINT_BITS * W;
    slong k0;

    /* x = log2 K(n1), K(n1) the largest k summed directly (measured
       optimum on 1 to 8 threads, on Sapphire Rapids and Zen 3; flat
       within some 10% around it, as is K0 between 6 and 12, larger K0
       being slower with many threads; see p-const_khinchin_tune) */
    if (p < 300)
        x = 3.0;
    else if (p < 700)
        x = 5.0;
    else if (p < 1500)
        x = 6.0;
    else if (p < 6000)
        x = 7.5;
    else if (p < 30000)
        x = 9.0;
    else if (p < 60000)
        x = 10.0;
    else if (p < 150000)
        x = 11.0;
    else if (p < 600000)
        x = 11.5;
    else if (p < 3000000)
        x = 13.0;
    else
    {
        /* beyond, where B_2n1 has about p bits (B_2n has about
           2n log2(2n / (2 pi e)) bits): the direct columns at n1 then
           take about as many powers k^(-2n), K(n1) = 2^x, as the
           iterator takes j^(-2n), j <= 2n / (2 pi e), with the same
           precision; this matches the tuned values below (12.2 at 10^6
           bits, 13.8 at 3.3 10^6) and grows to 16.8 at 3.3 10^7 bits
           (10^7 digits), where the Bernoulli columns, which grow like
           p^3, would otherwise take nearly all the time */
        double lo = 1000.0, hi = 1e12, m;
        int i;
        for (i = 0; i < 100; i++)
        {
            m = 0.5 * (lo + hi);
            if (m * _log2_d(m / KHINCHIN_2PIE) < p)
                lo = m;
            else
                hi = m;
        }
        x = FLINT_MAX(14.0, p / lo);
    }

    if (p < 20000)
        k0 = 8;
    else if (p < 200000)
        k0 = 6;
    else if (p < 600000)
        k0 = 10;
    else
        k0 = 12;

    *K0 = k0;
    *n1 = FLINT_MAX((slong) (p / (2.0 * x)), 1);
}

/* the most powers (2 pi)^(-2j) held at once */
#define KHINCHIN_BMAX 128

/* the most bits of the power tables of the Bernoulli iterators running
   at once (1 GiB), unless a single one needs more */
#define KHINCHIN_TABLE_MAX 8589934592.0

/* direct groups per thread, and the fewest columns per range of the
   subtracted powers */
#define KHINCHIN_DGROUPS 2
#define KHINCHIN_SUB_MIN 2000

/* the direct columns: columns per batch, limbs per block, powers per
   pass of the tiled sums */
#define KHINCHIN_DCOLS 16
#define KHINCHIN_DTILE 256
#define KHINCHIN_DGRAB 8

/* the divisions of a group of columns in parallel from this precision */
#define KHINCHIN_PAR_DIV 1000

/* ---- the Bernoulli chunks ---- */

/* the bits of B_2m */
static double
_khinchin_bernoulli_bits(slong m)
{
    double x = 2.0 * m;
    return FLINT_MAX(x * _log2_d(x / KHINCHIN_2PIE) + 64.0, 64.0);
}

/* cum[j] = the estimated cost of the Bernoulli columns j <= m < n1 (the
   iterator and the product by a full-precision power, at p bits),
   0 < j <= n1 */
static double *
_khinchin_bernoulli_costs(slong n1, double p)
{
    double * cum = flint_malloc((n1 + 1) * sizeof(double)), bits;
    slong j;

    cum[n1] = 0.0;
    for (j = n1 - 1; j >= 1; j--)
    {
        bits = _khinchin_bernoulli_bits(j);
        cum[j] = cum[j + 1] + pow(bits, 1.5) + p * pow(bits, 0.5);
    }
    return cum;
}

/* the chunks [zlo[c], zhi[c]] of the columns 1 <= m < n1, chunk 0 on top,
   with costs in proportion to budgets of threads bud[c] (nz of them,
   summing to nt); returns the number of chunks (fewer if the columns
   run out) */
static slong
_khinchin_bernoulli_split(slong * zlo, slong * zhi, const slong * bud,
    slong nz, slong n1, slong nt, const double * cum)
{
    slong g = 0, j, bsum = bud[0];

    zhi[0] = n1 - 1;
    for (j = n1 - 1; j > 1 && g < nz - 1; j--)
    {
        if (cum[j] >= cum[1] * bsum / nt)
        {
            zlo[g] = j;
            g++;
            bsum += bud[g];
            zhi[g] = j - 1;
        }
    }
    zlo[g] = 1;
    return g + 1;
}

static void
_khinchin_bernoulli_budgets(slong * bud, slong nz, slong nt)
{
    slong c;
    for (c = 0; c < nz; c++)
        bud[c] = nt / nz + (c < nt % nz);
}

/* the power table of the iterator from B_N holds about N^2 / 24 bits,
   about 0.55 of that from N = 10^5 (where it is sieved by 3 and 5) */
static double
_khinchin_table_bits(slong N)
{
    return (double) N * (double) N / 24.0 * ((N >= 100000) ? 0.55 : 1.0);
}

/* The number of Bernoulli chunks run at once on nt threads: one per
   thread, the fastest, unless the power tables of their iterators
   would exceed KHINCHIN_TABLE_MAX bits in all, in which case the most
   chunks that stay below it (or below the table of a single iterator)
   among ceil(nt / T), T = 2, 3, ..., each with about T threads for the
   parallelism inside it. */
static slong
_khinchin_lanes(slong n1, slong nt, const double * cum)
{
    slong nb, c, nz, T, * zlo, * zhi, * bud;
    double lim, tot;

    if (nt <= 1 || n1 <= 2)
        return 1;

    lim = FLINT_MAX(KHINCHIN_TABLE_MAX, _khinchin_table_bits(2 * (n1 - 1)));
    zlo = flint_malloc(3 * (nt + 1) * sizeof(slong));
    zhi = zlo + nt + 1;
    bud = zhi + nt + 1;

    /* nb = ceil(nt / T) for T = 1, 2, ... threads per chunk */
    for (T = 1; ; T++)
    {
        nb = FLINT_MIN((nt + T - 1) / T, n1 - 1);
        if (nb <= 1)
            break;
        _khinchin_bernoulli_budgets(bud, nb, nt);
        nz = _khinchin_bernoulli_split(zlo, zhi, bud, nb, n1, nt, cum);
        tot = 0.0;
        for (c = 0; c < nz; c++)
            tot += _khinchin_table_bits(2 * zhi[c]);
        if (tot <= lim)
            break;
    }

    flint_free(zlo);
    return FLINT_MAX(nb, 1);
}


void
_mp_real_const_khinchin_params(slong * K0, slong * n1, slong * nb, slong n)
{
    slong nt = (n < 64) ? 1 : flint_get_num_threads();
    double * cum;

    _khinchin_params(K0, n1, n);
    cum = _khinchin_bernoulli_costs(*n1, (double) FLINT_BITS * n);
    *nb = _khinchin_lanes(*n1, nt, cum);
    flint_free(cum);
}

/* ---- powers ---- */

/* a power k^(-2n) as (X, len) B^-(W + e), the integer X within eps of
   k^(-2n) B^(W + e); the limbs live at buf[off ...] */
typedef struct
{
    nn_ptr buf;
    slong cap;
    slong off;
    slong len;
    slong e;
    double eps;
    ulong k2;
}
_khinchin_pow_struct;

/* create k^(-2m) (k < 2^(FLINT_BITS / 2)) with a relative accuracy of
   about r limbs: (X, len) = floor(k^(-2m) B^(W + e)) to within a few
   ulps, e chosen so that X has r + 2 limbs (e >= 1) */
static void
_khinchin_pow_init(_khinchin_pow_struct * P, ulong k, ulong m, slong W, slong r)
{
    mp_real_t v, one;
    ulong err;
    slong vtop, e, h;

    /* about the number of limbs of k^(-2m) B^W (at most one more) */
    vtop = (slong) ceil(((double) FLINT_BITS * W - 2.0 * m * _log2_d((double) k))
        / FLINT_BITS + 1e-9);
    vtop = FLINT_MAX(vtop, 0);
    e = FLINT_MAX(r + 2 - vtop, 1);
    h = vtop + e + 1;

    mp_real_init(v);
    mp_real_init(one);

    P->cap = h + 32;
    P->buf = flint_malloc(P->cap * sizeof(ulong));
    P->off = 0;
    P->e = e;
    P->k2 = k * k;

    mp_real_set_ui(one, k * k);
    mp_real_pow_ui(v, one, m, h + 2);
    mp_real_set_ui(one, 1);
    mp_real_div(v, one, v, h + 2);
    /* k^(-2m) B^(W + e - h) < 1, the floor of which times B^h is
       floor(k^(-2m) B^(W + e)) */
    mp_real_mul_2exp_si(v, v, FLINT_BITS * (W + e - h));
    _mp_real_get_fixed(P->buf, &err, v, h);
    if (err == UWORD_MAX)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    P->len = h;
    while (P->len > 0 && P->buf[P->len - 1] == 0)
        P->len--;
    P->eps = (double) err;

    mp_real_clear(v);
    mp_real_clear(one);
}

static void
_khinchin_pow_clear(_khinchin_pow_struct * P)
{
    flint_free(P->buf);
}

/* ---- summation by parts in fixed point ---- */

/* (z, xn) = floor(x / (a b)), by two divisions if a b does not fit a
   limb (the floor of the floor being the floor) */
static void
_khinchin_divrem_ab(nn_ptr z, nn_srcptr x, slong xn, ulong a, ulong b)
{
    ulong hi, lo;
    umul_ppmm(hi, lo, a, b);
    if (hi == 0)
        mpn_divrem_1(z, 0, x, xn, lo);
    else
    {
        mpn_divrem_1(z, 0, x, xn, a);
        mpn_divrem_1(z, 0, z, xn, b);
    }
}

/* the columns of a range [lo, hi] taken downward: R = R_n and
   T = sum_(i >= n, i >= 2) R_i / ((2i - 2)(2i - 1)), with W fraction
   limbs and one unit limb (and one spare) */
typedef struct
{
    nn_ptr R;
    nn_ptr T;
    nn_ptr tmp;
    slong rlen;
    slong W;
    double err;     /* bound for the error of every R_i, in ulps */
    double ndiv;    /* the number of divisions into T */
}
_khinchin_abel_struct;

static void
_khinchin_abel_init(_khinchin_abel_struct * A, slong W)
{
    A->R = flint_calloc(W + 2, sizeof(ulong));
    A->T = flint_calloc(W + 2, sizeof(ulong));
    A->tmp = flint_malloc((W + 2) * sizeof(ulong));
    A->rlen = 0;
    A->W = W;
    A->err = 0.0;
    A->ndiv = 0.0;
}

static void
_khinchin_abel_clear(_khinchin_abel_struct * A)
{
    flint_free(A->R);
    flint_free(A->T);
    flint_free(A->tmp);
}

/* column n with (Z, zlen) = X_n (zlen <= W + 1) to within zerr ulps;
   Z is destroyed */
static void
_khinchin_abel_step(_khinchin_abel_struct * A, ulong n, nn_ptr Z, slong zlen,
    double zerr)
{
    slong W = A->W, rl;
    ulong cy;

    /* the error of X_n / n and the truncation of the division */
    A->err += zerr / n + 1.0;

    while (zlen > 0 && Z[zlen - 1] == 0)
        zlen--;

    if (zlen > 0)
    {
        /* R += floor(Z / n) */
        mpn_divrem_1(Z, 0, Z, zlen, n);
        while (zlen > 0 && Z[zlen - 1] == 0)
            zlen--;
        if (zlen > A->rlen)
        {
            flint_mpn_zero(A->R + A->rlen, zlen - A->rlen);
            A->rlen = zlen;
        }
        if (zlen > 0)
        {
            cy = mpn_add(A->R, A->R, A->rlen, Z, zlen);
            if (cy)
                A->R[A->rlen++] = cy;
        }
    }

    rl = A->rlen;
    FLINT_ASSERT(rl <= W + 1);

    if (n >= 2 && rl > 0)
    {
        /* T += floor(R / ((2n - 2)(2n - 1))) */
        _khinchin_divrem_ab(A->tmp, A->R, rl, 2 * n - 2, 2 * n - 1);
        cy = mpn_add_n(A->T, A->T, A->tmp, rl);
        if (cy)
            mpn_add_1(A->T + rl, A->T + rl, W + 2 - rl, cy);
        A->ndiv += 1.0;
    }
}

/* res = A_(lo-1) R_lo - T (Alo = NULL for lo = 1: R_1 - T), with the
   error bound (1 + sum_(i>=2) |a_i|) err + ndiv, 1 + sum |a_i|
   = 2 - log 2 < 1.31 */
static void
_khinchin_abel_get(mp_real_t res, const _khinchin_abel_struct * A,
    const mp_real_t Alo, slong prec)
{
    slong W = A->W;
    mp_real_t t;

    mp_real_init(t);
    _mp_real_set_mpn_2exp(res, A->R, W + 2, -FLINT_BITS * W);
    if (Alo != NULL)
        mp_real_mul(res, res, Alo, prec);
    _mp_real_set_mpn_2exp(t, A->T, W + 2, -FLINT_BITS * W);
    mp_real_sub(res, res, t, prec);
    _mp_real_add_error_ulps_at(res, (1.31 * A->err + A->ndiv) * MP_REAL_EPS, -W);
    mp_real_clear(t);
}

/* ---- the rows ---- */

/* sum_(k=2)^(K0-1) f(1/k), and log 2, at prec limbs */
static void
_khinchin_rows(mp_real_t res, mp_real_t log2, slong K0, slong prec)
{
    slong num, i, j, k, a;
    ulong * primes;
    slong * C, * va, * vb;
    mp_real_struct * lam;
    mp_real_t s, t;
    n_primes_t iter;

    num = n_prime_pi(K0);
    if (num > MP_REAL_LOG_PRIMES_MAX)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    primes = flint_malloc(num * sizeof(ulong));
    n_primes_init(iter);
    for (i = 0; i < num; i++)
        primes[i] = n_primes_next(iter);
    n_primes_clear(iter);

    C = flint_calloc(num * num, sizeof(slong));
    va = flint_malloc(2 * num * sizeof(slong));
    vb = va + num;

    for (k = 2; k < K0; k++)
    {
        /* the exponents of log(k / (k - 1)) and log((k + 1) / k) */
        for (j = 0; j < num; j++)
        {
            slong e0 = 0, e1 = 0, e2 = 0;
            ulong p = primes[j], m;
            for (m = k - 1; m % p == 0; m /= p) e0++;
            for (m = k; m % p == 0; m /= p) e1++;
            for (m = k + 1; m % p == 0; m /= p) e2++;
            va[j] = e1 - e0;
            vb[j] = e2 - e1;
        }
        for (i = 0; i < num; i++)
            if (va[i] != 0)
                for (j = 0; j < num; j++)
                    C[i * num + j] += va[i] * vb[j];
    }

    lam = flint_malloc(num * sizeof(mp_real_struct));
    for (i = 0; i < num; i++)
        mp_real_init(lam + i);
    mp_real_init(s);
    mp_real_init(t);

    _mp_real_log_primes_vec(lam, num, prec);

    mp_real_zero(res);
    for (i = 0; i < num; i++)
    {
        mp_real_zero(s);
        for (j = i; j < num; j++)
        {
            a = (j == i) ? C[i * num + i] : C[i * num + j] + C[j * num + i];
            if (a > 0)
                mp_real_addmul_ui(s, s, lam + j, a, prec);
            else if (a < 0)
                mp_real_submul_ui(s, s, lam + j, -a, prec);
        }
        if (!mp_real_is_zero(s))
        {
            mp_real_mul(t, lam + i, s, prec);
            mp_real_add(res, res, t, prec);
        }
    }

    mp_real_set(log2, lam + 0);

    for (i = 0; i < num; i++)
        mp_real_clear(lam + i);
    flint_free(lam);
    mp_real_clear(s);
    mp_real_clear(t);
    flint_free(C);
    flint_free(va);
    flint_free(primes);
}

/* ---- the direct columns ---- */

/* the k of the direct columns: those in [K0, 2 K0) and the odd k >= 2 K0
   with n_hi(k) >= n1 (num of them, ks[num] the next one), and n_hi */
static slong
_khinchin_direct_ks(ulong ** ks, slong ** nhi, slong K0, slong n1, double F)
{
    slong num, i;
    ulong k;

    num = 0;
    for (k = K0; _khinchin_n_hi(k, F) >= n1; k += (k < 2 * K0 - 1) ? 1 : 2)
        num++;

    *ks = flint_malloc((num + 1) * sizeof(ulong));
    *nhi = flint_malloc((num + 1) * sizeof(slong));
    for (i = 0, k = K0; i <= num; i++, k += (k < 2 * K0 - 1) ? 1 : 2)
    {
        (*ks)[i] = k;
        (*nhi)[i] = _khinchin_n_hi(k, F);
    }

    return num;
}

/* limbs of k^(-2 n1) above B^-W (rounded up, plus one for the
   doubles) */
static slong
_khinchin_direct_E(ulong k, slong n1, double F)
{
    slong E = (slong) ((F - 2.0 * n1 * _log2_d((double) k)) / FLINT_BITS) + 2;
    return FLINT_MAX(E, 1);
}

/* the column n from the partial sum (X, xlen) with error zerr ulps
   (omitted terms included): Z_n = X (1 + 4^-n + 4^-2n + ...), each
   shifted term truncated, into A; Z (W + 2 limbs) is scratch */
static void
_khinchin_direct_column(_khinchin_abel_struct * A, nn_srcptr X, slong xlen,
    slong n, double zerr, nn_ptr Z)
{
    slong zlen, j;

    flint_mpn_copyi(Z, X, xlen);
    zlen = xlen;
    for (j = 1; ; j++)
    {
        ulong sh = 2 * (ulong) n * j, sl = sh / FLINT_BITS, sb = sh % FLINT_BITS;
        ulong cy;
        slong tl;

        if ((slong) sl >= xlen)
            break;
        tl = xlen - sl;
        if (sb != 0)
            mpn_rshift(A->tmp, X + sl, tl, sb);
        else
            flint_mpn_copyi(A->tmp, X + sl, tl);
        cy = mpn_add(Z, Z, zlen, A->tmp, tl);
        if (cy)
            Z[zlen++] = cy;
        zerr += 1.0;
    }
    /* the shifted terms left out are below one ulp each, below 4/3 ulps
       together */
    zerr += 2.0;
    zerr *= 4.0 / 3.0;

    _khinchin_abel_step(A, n, Z, zlen, zerr);
}

/* the omitted terms of column n, each below B^-W, in ulps:
   sum_(k >= k1) k^(-2n) <= k1^(-2n) (1 + k1 / (2n - 1)) */
static double
_khinchin_direct_omitted(ulong k1, slong n, double F)
{
    double t = F - 2.0 * n * _log2_d((double) k1);
    return ldexp(1.0, (int) ceil(t + 1e-6)) * (1.0 + (double) k1 / (2.0 * n - 1.0));
}

/* the segments of a power: the range [n1, n_hi] of its columns is cut
   into nseg segments of equal length, each created afresh with the
   relative accuracy needed at its bottom, so that the multiplications
   taper off with the accuracy actually needed */
typedef struct
{
    _khinchin_pow_struct P;
    slong seg, nseg, segb;
    /* in a batch: the first active step, whether created in the batch,
       the limb at which X index 0 sits, the top limb (exclusive) once
       all steps are done, and the drops before each step */
    slong s, e0, top;
    int fresh;
    unsigned char d[KHINCHIN_DCOLS + 1];
}
_khinchin_dpow_struct;

/* (re)create the power of k at column m, in the segment that contains
   the column nL (the bottom of the batch) */
static void
_khinchin_dpow_create(_khinchin_dpow_struct * D, ulong k, slong nh, slong m,
    slong nL, slong n1, slong W, double F, int first)
{
    slong R = nh - n1 + 1;

    if (first)
    {
        D->nseg = (R >= 64) ? 3 : (R >= 24) ? 2 : 1;
        D->seg = 0;
        D->segb = nh + 1;
    }
    else
        _khinchin_pow_clear(&D->P);

    while (D->segb > nL && D->seg < D->nseg)
    {
        D->seg++;
        D->segb = nh - (D->seg * R) / D->nseg;
    }

    _khinchin_pow_init(&D->P, k, m, W, _khinchin_direct_E(k, D->segb, F) + 1);
    D->fresh = 1;
}

/* One batch of the columns n0 - t, t < cnt, of the powers D[0 ... np)
   (each from its step s, with no product at that step if fresh): the
   partial sums X[t] (L limbs each) of the parts of the powers above
   B^-W, and their errors added to zerr[t].  The columns are taken in
   batches with one partial sum per column, and the powers
   KHINCHIN_DGRAB at a time by blocks of KHINCHIN_DTILE limbs of the X
   positions: each block of a power is stepped through the columns of
   the batch while in cache, and the blocks of the X[t] serve all the
   powers of the pass, instead of streaming each power and its column sum
   once per column, which is bound by the memory bandwidth once the
   powers outgrow the L2 cache (about twice slower on one thread for
   powers of 2 10^5 limbs, and more so on several threads).

   A step X <- k^2 X (k^2 < B, so one more limb at most) multiplies the
   error of the power, in units of its bottom limb, by k^2; once that
   error reaches 2^20 B the bottom limb is noise and is dropped after the
   step, which adds one unit of the next limb (kept negligible against
   the error, so that the losses do not accumulate over the drops).  The
   drops depend only on the error counts, so they are scheduled before
   the products; the carries of the products between blocks are kept per
   power and step.  The part above B^-W of a power with e extended limbs
   is within eps (e = 0), 1 + eps / B (e = 1) or 2 ulps (eps < 2^20 B k^2
   < B^2 when e >= 2). */
static void
_khinchin_tiled_batch(_khinchin_dpow_struct * D, slong np, slong cnt,
    nn_ptr * X, slong L, double * zerr)
{
    ulong car[KHINCHIN_DGRAB][KHINCHIN_DCOLS];
    slong i, t, g;

    /* the schedule of each power: room, drops, errors */
    for (i = 0; i < np; i++)
    {
        _khinchin_dpow_struct * Di = D + i;
        _khinchin_pow_struct * P = &Di->P;
        double eps;
        slong e, d, nmul;

        /* room for cnt more limbs */
        if (P->off + P->len + cnt + 1 > P->cap)
        {
            if (P->len + cnt + 1 > P->cap)
            {
                nn_ptr nb;
                P->cap = P->len + cnt + 33;
                nb = flint_malloc(P->cap * sizeof(ulong));
                flint_mpn_copyi(nb, P->buf + P->off, P->len);
                flint_free(P->buf);
                P->buf = nb;
            }
            else if (P->off != 0)
                memmove(P->buf, P->buf + P->off, P->len * sizeof(ulong));
            P->off = 0;
        }

        eps = P->eps;
        e = P->e;
        d = 0;
        nmul = 0;
        for (t = 0; t < cnt; t++)
        {
            Di->d[t] = d;
            if (t < Di->s)
                continue;
            if (t > Di->s || !Di->fresh)
            {
                nmul++;
                eps = eps * (double) P->k2 * MP_REAL_EPS;
                if (eps >= 0x1p20 * MP_REAL_D_B && e > 0)
                {
                    d++;
                    e--;
                    eps = eps * MP_REAL_D_BINV * MP_REAL_EPS + 1.0;
                }
            }
            if (e == 0)
                zerr[t] += eps;
            else if (e == 1)
                zerr[t] += 1.0 + eps * MP_REAL_D_BINV * MP_REAL_EPS;
            else
                zerr[t] += 2.0;
        }
        Di->d[cnt] = d;
        Di->e0 = P->e;
        /* the power grows by at most a limb per step, but stays below
           B^(W + e0) in the frame of the batch (k^(-2m) < 1, and the
           truncations only lower it): its part above B^-W has at most
           W limbs, within the L = W + 2 limbs of the sums */
        Di->top = FLINT_MIN(P->len + nmul, P->e + (L - 2));
        FLINT_ASSERT(P->len <= Di->top);
        flint_mpn_zero(P->buf + P->off + P->len, nmul + 1);
        P->eps = eps;
    }

    for (t = 0; t < cnt; t++)
        flint_mpn_zero(X[t], L);

    /* the products and sums, KHINCHIN_DGRAB powers at a time */
    for (g = 0; g < np; g += KHINCHIN_DGRAB)
    {
        slong gn = FLINT_MIN(KHINCHIN_DGRAB, np - g), xlo, xhi, o, u;

        xlo = WORD_MAX;
        xhi = WORD_MIN;
        for (u = 0; u < gn; u++)
        {
            _khinchin_dpow_struct * Di = D + g + u;
            xlo = FLINT_MIN(xlo, -Di->e0);
            xhi = FLINT_MAX(xhi, Di->top - Di->e0);
            for (t = 0; t < cnt; t++)
                car[u][t] = 0;
        }
        /* blocks aligned to multiples of the tile in X positions */
        xlo = (xlo >= 0) ? (xlo / KHINCHIN_DTILE) * KHINCHIN_DTILE
            : -(((-xlo) + KHINCHIN_DTILE - 1) / KHINCHIN_DTILE) * KHINCHIN_DTILE;

        for (o = xlo; o < xhi; o += KHINCHIN_DTILE)
        {
            for (u = 0; u < gn; u++)
            {
                _khinchin_dpow_struct * Di = D + g + u;
                nn_ptr p = Di->P.buf + Di->P.off;
                ulong k2 = Di->P.k2;
                slong jb = FLINT_MAX(o + Di->e0, 0);
                slong je = FLINT_MIN(o + KHINCHIN_DTILE + Di->e0, Di->top);
                slong ab = FLINT_MAX(jb, Di->e0);

                if (jb >= je)
                    continue;

                for (t = Di->s; t < cnt; t++)
                {
                    /* the step: the limbs from the bottom one left */
                    if (t > Di->s || !Di->fresh)
                    {
                        slong lo = FLINT_MAX(jb, Di->d[t]);
                        if (lo < je)
                        {
                            ulong cin = (lo == jb) ? car[u][t] : 0, hi;
                            hi = mpn_mul_1(p + lo, p + lo, je - lo, k2);
                            hi += mpn_add_1(p + lo, p + lo, je - lo, cin);
                            car[u][t] = hi;
                        }
                        else
                            car[u][t] = 0;
                    }

                    /* the part above B^-W into X[t] */
                    if (ab < je)
                    {
                        slong x = ab - Di->e0;
                        ulong cy = mpn_add_n(X[t] + x, X[t] + x, p + ab, je - ab);
                        FLINT_ASSERT(x + (je - ab) <= L - 2);
                        if (cy)
                        {
                            mpn_add_1(X[t] + x + (je - ab), X[t] + x + (je - ab),
                                L - x - (je - ab), cy);
                        }
                    }
                }
            }
        }
    }

    /* the new state of the powers */
    for (i = 0; i < np; i++)
    {
        _khinchin_dpow_struct * Di = D + i;
        _khinchin_pow_struct * P = &Di->P;
        slong d = Di->d[cnt];

        P->off += d;
        P->e -= d;
        P->len = Di->top - d;
        while (P->len > 0 && P->buf[P->off + P->len - 1] == 0)
            P->len--;
        Di->fresh = 0;
        Di->s = 0;
    }
}

/* the partial sums X[t] (L limbs) and the scratch of the batches */
static nn_ptr
_khinchin_tiled_sums(nn_ptr * X, slong L)
{
    nn_ptr Xbuf = flint_malloc(KHINCHIN_DCOLS * L * sizeof(ulong));
    slong t;
    for (t = 0; t < KHINCHIN_DCOLS; t++)
        X[t] = Xbuf + t * L;
    return Xbuf;
}

/* The part of the direct columns n1 <= n <= M from the k = ks[i],
   a <= i < b, into A, in batches of KHINCHIN_DCOLS columns.  Only k in
   [K0, 2 K0) and the odd k >= 2 K0 are summed: the even k >= 2 K0
   contribute 4^-n Z_n, so the column is X_n / (1 - 4^-n) for the
   partial sum X_n (the sums being linear, this holds for every part
   separately).  A power enters at its column n_hi(k); the powers that
   enter in a batch, and those whose segment ends in it, are created for
   the column where they enter, resp. the top of the batch, with the
   accuracy of the segment containing the bottom of the batch.  The part
   containing the first k not yet entered in column n bounds the omitted
   terms. */
static void
_khinchin_direct(_khinchin_abel_struct * A, const ulong * ks, const slong * nhi,
    slong a, slong b, slong n1, slong W)
{
    double F = (double) FLINT_BITS * W, zerr[KHINCHIN_DCOLS];
    slong n0, nL, cnt, act, q, i, t, L = W + 2;
    _khinchin_dpow_struct * D;
    nn_ptr Xbuf, X[KHINCHIN_DCOLS], Z;

    if (a >= b)
        return;

    D = flint_malloc((b - a) * sizeof(_khinchin_dpow_struct));
    Xbuf = _khinchin_tiled_sums(X, L);
    Z = flint_malloc(L * sizeof(ulong));
    act = a;
    q = a;

    for (n0 = nhi[a]; n0 >= n1; n0 -= cnt)
    {
        nL = FLINT_MAX(n0 - KHINCHIN_DCOLS + 1, n1);
        cnt = n0 - nL + 1;

        /* the powers entering in this batch (fresh), and those whose
           segment ends in it */
        for (i = 0; i < act - a; i++)
            if (D[i].segb > nL)
                _khinchin_dpow_create(D + i, ks[a + i], nhi[a + i], n0, nL, n1, W, F, 0);
        while (act < b && nhi[act] >= nL)
        {
            _khinchin_dpow_struct * Di = D + (act - a);
            _khinchin_dpow_create(Di, ks[act], nhi[act], nhi[act], nL, n1, W, F, 1);
            Di->s = n0 - nhi[act];
            act++;
        }

        for (t = 0; t < cnt; t++)
            zerr[t] = 0.0;
        _khinchin_tiled_batch(D, act - a, cnt, X, L, zerr);

        for (t = 0; t < cnt; t++)
        {
            slong n = n0 - t, xlen = L;

            while (q < b && nhi[q] >= n)
                q++;
            if (nhi[b] < n)
                zerr[t] += _khinchin_direct_omitted(ks[q], n, F);

            while (xlen > 0 && X[t][xlen - 1] == 0)
                xlen--;
            _khinchin_direct_column(A, X[t], xlen, n, zerr[t], Z);
        }
    }

    for (i = 0; i < act - a; i++)
        _khinchin_pow_clear(&D[i].P);
    flint_free(D);
    flint_free(Xbuf);
    flint_free(Z);
}

/* ---- the Bernoulli columns ---- */

/* sum_(n=lo)^hi c_n (1 + sum_(k=2)^(K0-1) k^(-2n)), into A by summation
   by parts (the part A_(lo-1) R_lo being applied by the caller), the
   powers in batches of columns like the direct columns */
static void
_khinchin_bernoulli_sub(_khinchin_abel_struct * A, slong K0, slong lo, slong hi, slong W)
{
    double zerr[KHINCHIN_DCOLS];
    slong np = K0 - 2, n0, cnt, k, t, L = W + 2;
    _khinchin_dpow_struct * D;
    nn_ptr Xbuf, X[KHINCHIN_DCOLS];

    D = flint_malloc(FLINT_MAX(np, 1) * sizeof(_khinchin_dpow_struct));
    Xbuf = _khinchin_tiled_sums(X, L);

    /* needed down to n = lo (n = 1: k^-2 has W limbs) */
    for (k = 2; k < K0; k++)
    {
        _khinchin_dpow_struct * Dk = D + (k - 2);
        _khinchin_pow_init(&Dk->P, k, hi, W,
            FLINT_MIN(_khinchin_direct_E(k, lo, (double) FLINT_BITS * W) + 1, W + 2));
        Dk->s = 0;
        Dk->fresh = 1;
    }

    for (n0 = hi; n0 >= lo; n0 -= cnt)
    {
        cnt = FLINT_MIN(KHINCHIN_DCOLS, n0 - lo + 1);

        for (t = 0; t < cnt; t++)
            zerr[t] = 0.0;
        _khinchin_tiled_batch(D, np, cnt, X, L, zerr);

        for (t = 0; t < cnt; t++)
        {
            slong xlen = L;

            /* the 1 */
            mpn_add_1(X[t] + W, X[t] + W, L - W, 1);
            while (xlen > 0 && X[t][xlen - 1] == 0)
                xlen--;
            _khinchin_abel_step(A, n0 - t, X[t], xlen, zerr[t]);
        }
    }

    for (k = 2; k < K0; k++)
        _khinchin_pow_clear(&D[k - 2].P);
    flint_free(D);
    flint_free(Xbuf);
}

/* A_n = 1 - 1/2 + ... + 1/(2n - 1) = H_(2n-1) - H_(n-1), n >= 1 */
static void
_khinchin_A(mp_real_t res, slong n, slong prec)
{
    fmpq_t h, g;
    mp_real_t t;

    fmpq_init(h);
    fmpq_init(g);
    mp_real_init(t);
    arith_harmonic_number(h, 2 * n - 1);
    arith_harmonic_number(g, n - 1);
    fmpq_sub(h, h, g);
    mp_real_set_fmpz(res, fmpq_numref(h));
    mp_real_set_fmpz(t, fmpq_denref(h));
    mp_real_div(res, res, t, prec);
    fmpq_clear(h);
    fmpq_clear(g);
    mp_real_clear(t);
}

/* the columns m = top - i of a group, i < len (j = hi - m <= b):
   x_i = zeta(2m) / (m t_hi) = |B_2m| c_m Y_j / (2 m denom) */
typedef struct
{
    const fmpz * const * num;
    const fmpz * const * den;
    const fmpz * c;
    mp_real_struct * x;
    const mp_real_struct * Y;
    slong top, hi, prec;
}
_khinchin_col_struct;

static void
_khinchin_col(slong i, void * args)
{
    const _khinchin_col_struct * C = args;
    slong m = C->top - i, prec = C->prec;
    const fmpz * den = C->den[i];
    mp_real_struct * x = C->x + i;
    fmpz_t Nc;
    mp_real_t u;

    fmpz_init(Nc);
    mp_real_init(u);

    fmpz_mul(Nc, C->num[i], C->c + i);
    fmpz_abs(Nc, Nc);
    mp_real_set_fmpz(u, Nc);
    mp_real_mul(x, u, C->Y + (C->hi - m), prec);

    if (fmpz_abs_fits_ui(den) && fmpz_bits(den) + FLINT_BIT_COUNT(2 * m) <= FLINT_BITS)
    {
        mp_real_div_ui(x, x, 2 * m * fmpz_get_ui(den), prec);
    }
    else
    {
        mp_real_set_fmpz(u, den);
        mp_real_mul_ui(u, u, 2 * m, prec);
        mp_real_div(x, x, u, prec);
    }

    fmpz_clear(Nc);
    mp_real_clear(u);
}

/* x_i = Rx_i / ((2m - 2)(2m - 1)), m = top - i >= 2 */
static void
_khinchin_col_div(slong i, void * args)
{
    const _khinchin_col_struct * C = args;
    slong m = C->top - i;

    if (m >= 2)
    {
        ulong a = 2 * m - 2, b = 2 * m - 1, hi, lo;
        umul_ppmm(hi, lo, a, b);
        if (hi == 0)
            mp_real_div_ui(C->x + i, C->x + i, lo, C->prec);
        else
        {
            mp_real_div_ui(C->x + i, C->x + i, a, C->prec);
            mp_real_div_ui(C->x + i, C->x + i, b, C->prec);
        }
    }
}

/* For the columns lo0 <= m <= hi0: res = sum_(i=lo0)^hi0 a_i R_i,
   Rabs = R_lo0 and sig = sum_(i=lo0)^hi0 a_i, where
   R_i = sum_(n=i)^hi0 zeta(2n) / n, so that
   sum_(n=lo0)^hi0 c_n zeta(2n) = res + A_(lo0-1) Rabs.  y = (2 pi)^2 and
   Y_j = y^-j, j <= b, the block size.  With threads (the budget of the
   calling thread), the columns are taken in groups whose products by
   the powers Y_j run in parallel, as do the outputs of the iterator. */
static void
_khinchin_bernoulli_zeta(mp_real_t res, mp_real_t Rabs, mp_real_t sig0,
    slong lo0, slong hi0, const mp_real_t y, const mp_real_struct * Y, slong b,
    int use_cache, slong prec)
{
    slong m, hi, lo, g, gn, i, nt;
    mp_real_t t, x, u, Rx, Ux, sig;
    fmpz_t c, sn, sd;
    fmpz * num, * den, * cs;
    const fmpz ** nump, ** denp;
    mp_real_struct * xs;
    _khinchin_col_struct C;
    bernoulli_rev_t iter;

    nt = flint_get_num_threads();
    g = (nt == 1) ? 1 : nt * ((16 + nt - 1) / nt);

    mp_real_init(t);
    mp_real_init(x);
    mp_real_init(u);
    mp_real_init(Rx);
    mp_real_init(Ux);
    mp_real_init(sig);
    fmpz_init(c);
    fmpz_init(sn);
    fmpz_init(sd);
    num = _fmpz_vec_init(3 * g);
    den = num + g;
    cs = den + g;
    nump = flint_malloc(2 * g * sizeof(fmpz *));
    denp = nump + g;
    xs = flint_malloc(g * sizeof(mp_real_struct));
    for (i = 0; i < g; i++)
        mp_real_init(xs + i);

    C.num = nump;
    C.den = denp;
    C.c = cs;
    C.x = xs;
    C.Y = Y;
    C.prec = prec;

    /* t = t_hi = y^hi / (2 hi)! */
    hi = hi0;
    mp_real_pow_ui(t, y, hi, prec);
    fmpz_fac_ui(c, 2 * hi);
    mp_real_set_fmpz(u, c);
    mp_real_div(t, t, u, prec);

    mp_real_zero(res);
    mp_real_zero(Rabs);
    mp_real_zero(sig0);

    if (!use_cache)
        bernoulli_rev_init(iter, 2 * hi);

    while (hi >= lo0)
    {
        lo = FLINT_MAX(hi - b + 1, lo0);

        mp_real_zero(Rx);
        mp_real_zero(Ux);
        fmpz_zero(sn);
        fmpz_one(sd);
        fmpz_one(c);

        for (m = hi; m >= lo; m -= gn)
        {
            gn = FLINT_MIN(g, m - lo + 1);

            /* c = (2m + 1) ... (2 hi) and B_2m for the group */
            for (i = 0; i < gn; i++)
            {
                slong mm = m - i;

                if (mm != hi)
                    fmpz_mul2_uiui(c, c, 2 * mm + 1, 2 * mm + 2);
                fmpz_set(cs + i, c);

                if (use_cache)
                {
                    nump[i] = fmpq_numref(bernoulli_cache + 2 * mm);
                    denp[i] = fmpq_denref(bernoulli_cache + 2 * mm);
                }
                else
                {
                    bernoulli_rev_next(num + i, den + i, iter);
                    nump[i] = num + i;
                    denp[i] = den + i;
                }
            }

            C.top = m;
            C.hi = hi;
            if (gn >= 2 && nt >= 2)
                flint_parallel_do(_khinchin_col, &C, gn, nt, FLINT_PARALLEL_STRIDED);
            else
                for (i = 0; i < gn; i++)
                    _khinchin_col(i, &C);

            /* Rx = sum_(n >= m) zeta(2n) / (n t_hi) (the prefix sums of
               the group in place), Ux = sum a_i Rx_i (the divisions in
               parallel), sig = sum a_i (exactly, as sn / sd) */
            mp_real_add(xs + 0, xs + 0, Rx, prec);
            for (i = 1; i < gn; i++)
                mp_real_add(xs + i, xs + i, xs + i - 1, prec);
            mp_real_set(Rx, xs + gn - 1);

            if (gn >= 2 && nt >= 2 && prec >= KHINCHIN_PAR_DIV)
                flint_parallel_do(_khinchin_col_div, &C, gn, nt, FLINT_PARALLEL_STRIDED);
            else
                for (i = 0; i < gn; i++)
                    _khinchin_col_div(i, &C);

            for (i = 0; i < gn; i++)
            {
                slong mm = m - i;

                if (mm >= 2)
                {
                    mp_real_sub(Ux, Ux, xs + i, prec);
                    fmpz_mul2_uiui(sn, sn, 2 * mm - 2, 2 * mm - 1);
                    fmpz_sub(sn, sn, sd);
                    fmpz_mul2_uiui(sd, sd, 2 * mm - 2, 2 * mm - 1);
                }
                else
                {
                    mp_real_add(Ux, Ux, xs + i, prec);
                    fmpz_add(sn, sn, sd);
                }
            }
        }

        mp_real_set_fmpz(sig, sn);
        mp_real_set_fmpz(u, sd);
        mp_real_div(sig, sig, u, prec);

        /* sum_(i in block) a_i R_i = sig Rabs + t_hi Ux, where
           R_i = Rabs + t_hi Rx_i */
        mp_real_mul(x, sig, Rabs, prec);
        mp_real_add(res, res, x, prec);
        mp_real_mul(x, t, Ux, prec);
        mp_real_add(res, res, x, prec);
        mp_real_mul(x, t, Rx, prec);
        mp_real_add(Rabs, Rabs, x, prec);
        mp_real_add(sig0, sig0, sig, prec);

        /* t_(lo-1) = t_hi y^-(hi - lo + 1) (2 lo - 1)(2 lo) ... (2 hi) */
        if (lo > lo0)
        {
            fmpz_mul2_uiui(c, c, 2 * lo - 1, 2 * lo);
            mp_real_mul(t, t, Y + (hi - lo + 1), prec);
            mp_real_set_fmpz(u, c);
            mp_real_mul(t, t, u, prec);
        }

        hi = lo - 1;
    }

    if (!use_cache)
        bernoulli_rev_clear(iter);

    mp_real_clear(t);
    mp_real_clear(x);
    mp_real_clear(u);
    mp_real_clear(Rx);
    mp_real_clear(Ux);
    mp_real_clear(sig);
    fmpz_clear(c);
    fmpz_clear(sn);
    fmpz_clear(sd);
    _fmpz_vec_clear(num, 3 * g);
    flint_free(nump);
    for (i = 0; i < g; i++)
        mp_real_clear(xs + i);
    flint_free(xs);
}

/* ---- the tasks ---- */

/* the work is split into tasks: the Bernoulli columns in nz chunks
   (each with its own iterator and its own t_hi), the direct columns in
   nd groups of k (each with its own summation by parts), the
   subtracted powers of the Bernoulli columns in ns ranges of columns
   (likewise), and the rows */
typedef struct
{
    slong K0, n1, W, prec;
    int use_cache;
    /* Bernoulli chunks [zlo[c], zhi[c]], with budgets of zbud[c] threads */
    slong nz;
    slong * zlo, * zhi, * zbud;
    mp_real_struct * zres, * zR, * zsig;
    mp_real_struct y;
    mp_real_struct * Y;
    slong b;
    /* direct groups [da[g], da[g + 1]) of ks */
    slong nd;
    slong * da;
    const ulong * ks;
    const slong * nhi;
    _khinchin_abel_struct * dA;
    /* subtracted powers in the columns [slo[c], shi[c]], with
       sAlo[c] = A_(slo[c] - 1) */
    slong ns;
    slong * slo, * shi;
    _khinchin_abel_struct * sA;
    mp_real_struct * sAlo;
    /* rows */
    mp_real_struct rows, log2;
}
_khinchin_work_struct;

static void
_khinchin_task(slong i, void * args)
{
    _khinchin_work_struct * w = args;

    if (i < w->nz)
    {
        _khinchin_bernoulli_zeta(w->zres + i, w->zR + i, w->zsig + i,
            w->zlo[i], w->zhi[i], &w->y, w->Y, w->b, w->use_cache, w->prec);
        return;
    }
    i -= w->nz;

    if (i < w->nd)
    {
        _khinchin_direct(w->dA + i, w->ks, w->nhi, w->da[i], w->da[i + 1],
            w->n1, w->W);
        return;
    }
    i -= w->nd;

    if (i < w->ns)
    {
        _khinchin_bernoulli_sub(w->sA + i, w->K0, w->slo[i], w->shi[i], w->W);
        if (w->slo[i] > 1)
            _khinchin_A(w->sAlo + i, w->slo[i] - 1, w->prec);
        return;
    }

    _khinchin_rows(&w->rows, &w->log2, w->K0, w->prec);
}

/* the tasks after the Bernoulli chunks */
static void
_khinchin_task_rest(slong i, void * args)
{
    _khinchin_work_struct * w = args;
    _khinchin_task(i + w->nz, args);
}

void
_mp_real_const_khinchin_tune(mp_real_t res, slong n, slong K0, slong n1, slong nb)
{
    slong W = n, prec = n + 2, M, num, i, j, nt, g;
    double F = (double) FLINT_BITS * W, tail, tot, acc, * cum;
    _khinchin_work_struct w;
    ulong * ks;
    slong * nhi;
    double * cost;
    mp_real_t s, c, Alo;

    if (K0 == 0 || n1 == 0)
    {
        slong a, b;
        _khinchin_params(&a, &b, W);
        if (K0 == 0) K0 = a;
        if (n1 == 0) n1 = b;
    }

    /* at most 64 primes up to K0 (for the rows), and the k summed
       directly k < 2^(F / 2 n1) <= 2^20 (n1 > W also ensures
       k < 2^(FLINT_BITS / 2), so that k^2 fits a limb) */
    K0 = FLINT_MAX(K0, 3);
    K0 = FLINT_MIN(K0, 300);
    M = _khinchin_n_hi(K0, F);
    n1 = FLINT_MAX(n1, W + 1);
    n1 = FLINT_MAX(n1, (slong) (F / 40.0) + 1);
    n1 = FLINT_MIN(n1, M + 1);

    nt = flint_get_num_threads();
    if (W < 64)
        nt = 1;

    w.K0 = K0;
    w.n1 = n1;
    w.W = W;
    w.prec = prec;

    mp_real_init(s);
    mp_real_init(c);
    mp_real_init(Alo);
    mp_real_init(&w.rows);
    mp_real_init(&w.log2);

    /* the Bernoulli chunks: nb of them (at most nt), run at once, with
       budgets of zbud[c] threads summing to nt and costs in proportion */
    cum = (n1 > 1) ? _khinchin_bernoulli_costs(n1, F) : NULL;
    if (nb <= 0)
        nb = (n1 > 1) ? _khinchin_lanes(n1, nt, cum) : 1;
    nb = FLINT_MAX(FLINT_MIN(nb, nt), 1);
    w.nz = (n1 > 1) ? FLINT_MIN(nb, n1 - 1) : 0;
    w.zlo = flint_malloc(3 * (w.nz + 1) * sizeof(slong));
    w.zhi = w.zlo + w.nz + 1;
    w.zbud = w.zhi + w.nz + 1;
    w.zres = flint_malloc(3 * (w.nz + 1) * sizeof(mp_real_struct));
    w.zR = w.zres + w.nz + 1;
    w.zsig = w.zR + w.nz + 1;
    for (i = 0; i < 3 * (w.nz + 1); i++)
        mp_real_init(w.zres + i);
    mp_real_init(&w.y);
    w.Y = NULL;
    w.b = 0;
    w.use_cache = 0;

    if (w.nz > 0)
    {
        _khinchin_bernoulli_budgets(w.zbud, w.nz, nt);
        w.nz = _khinchin_bernoulli_split(w.zlo, w.zhi, w.zbud, w.nz, n1, nt, cum);

        /* the cache (thread-local) only on one thread */
        w.use_cache = (nt == 1) && ((2 * (n1 - 1) < bernoulli_cache_num) || (n1 <= 65));
        if (w.use_cache)
        {
            BERNOULLI_ENSURE_CACHED(2 * (n1 - 1) + 1)
        }

        /* y = (2 pi)^2 = 64 (pi/4)^2, Y_j = y^-j; b + 1 powers and four
           products per block of b columns (at most KHINCHIN_BMAX, which
           bounds the memory of the powers at a cost of a few percent) */
        w.b = (slong) (2.0 * sqrt((double) (n1 - 1) / w.nz) + 1.0);
        w.b = FLINT_MIN(w.b, KHINCHIN_BMAX);
        w.b = FLINT_MAX(FLINT_MIN(w.b, n1 - 1), 1);
        mp_real_const_pi4(&w.y, prec, 1);
        mp_real_mul(&w.y, &w.y, &w.y, prec);
        mp_real_mul_2exp_si(&w.y, &w.y, 6);
        w.Y = flint_malloc((w.b + 1) * sizeof(mp_real_struct));
        for (j = 0; j <= w.b; j++)
            mp_real_init(w.Y + j);
        mp_real_set_ui(w.Y + 0, 1);
        mp_real_div(w.Y + 1, w.Y + 0, &w.y, prec);
        for (j = 2; j <= w.b; j++)
            mp_real_mul(w.Y + j, w.Y + j - 1, w.Y + 1, prec);
    }

    /* the direct groups, of about equal cost */
    num = _khinchin_direct_ks(&ks, &nhi, K0, n1, F);
    w.ks = ks;
    w.nhi = nhi;
    w.nd = (num > 0) ? FLINT_MIN(KHINCHIN_DGROUPS * nt, num) : 0;
    w.da = flint_malloc((w.nd + 1) * sizeof(slong));
    w.dA = flint_malloc(FLINT_MAX(w.nd, 1) * sizeof(_khinchin_abel_struct));
    for (g = 0; g < w.nd; g++)
        _khinchin_abel_init(w.dA + g, W);
    if (w.nd > 0)
    {
        cost = flint_malloc(num * sizeof(double));
        tot = 0.0;
        for (i = 0; i < num; i++)
        {
            double E = _khinchin_direct_E(ks[i], n1, F);
            cost[i] = (nhi[i] - n1 + 1.0) * (E + 4.0) + 20.0 * E;
            tot += cost[i];
        }
        w.da[0] = 0;
        acc = 0.0;
        g = 1;
        for (i = 0; i < num && g < w.nd; i++)
        {
            acc += cost[i];
            if (acc >= tot * g / w.nd)
                w.da[g++] = i + 1;
        }
        w.nd = g;
        w.da[g] = num;
        flint_free(cost);
    }

    /* the subtracted powers in ranges of about equal numbers of columns
       (at least KHINCHIN_SUB_MIN each) */
    w.ns = (n1 > 1) ? FLINT_MIN(nt, (n1 - 1 + KHINCHIN_SUB_MIN - 1) / KHINCHIN_SUB_MIN) : 0;
    w.slo = flint_malloc(2 * (w.ns + 1) * sizeof(slong));
    w.shi = w.slo + w.ns + 1;
    w.sA = flint_malloc((w.ns + 1) * sizeof(_khinchin_abel_struct));
    w.sAlo = flint_malloc((w.ns + 1) * sizeof(mp_real_struct));
    for (g = 0; g < w.ns; g++)
    {
        w.shi[g] = (n1 - 1) - ((n1 - 1) * g) / w.ns;
        w.slo[g] = (n1 - 1) - ((n1 - 1) * (g + 1)) / w.ns + 1;
        _khinchin_abel_init(w.sA + g, W);
        mp_real_init(w.sAlo + g);
    }

    /* serially on this thread when nt = 1 (which the use of the
       thread-local Bernoulli cache relies on); with one thread per
       Bernoulli chunk, all the tasks in one queue; otherwise the
       Bernoulli chunks with their budgets, then the rest */
    if (nt == 1)
    {
        for (i = 0; i < w.nz + w.nd + w.ns + 1; i++)
            _khinchin_task(i, &w);
    }
    else if (w.nz == 0 || w.zbud[0] == 1)
    {
        _mp_real_parallel_tasks(_khinchin_task, &w, w.nz + w.nd + w.ns + 1);
    }
    else
    {
        /* the chunks at once, chunk c with a budget of zbud[c] threads
           (if the pool grants fewer threads, the remaining chunks go to
           the threads that free up first) */
        double * zcost = flint_malloc(w.nz * sizeof(double));
        for (i = 0; i < w.nz; i++)
            zcost[i] = (double) w.zbud[i];
        _mp_real_parallel_tasks_cost(_khinchin_task, &w, w.nz, w.nz, zcost);
        flint_free(zcost);
        _mp_real_parallel_tasks(_khinchin_task_rest, &w, w.nd + w.ns + 1);
    }

    /* A_(n1-1) = sum of all a_i, i < n1, and the Bernoulli columns:
       chunk c contributes res_c + A_(lo_c - 1) R_c */
    mp_real_zero(s);
    mp_real_zero(Alo);
    for (g = w.nz - 1; g >= 0; g--)
    {
        mp_real_add(s, s, w.zres + g, prec);
        if (!mp_real_is_zero(Alo))
        {
            mp_real_mul(c, Alo, w.zR + g, prec);
            mp_real_add(s, s, c, prec);
        }
        mp_real_add(Alo, Alo, w.zsig + g, prec);
    }

    for (g = 0; g < w.ns; g++)
    {
        _khinchin_abel_get(c, w.sA + g, (w.slo[g] > 1) ? w.sAlo + g : NULL, prec);
        mp_real_sub(s, s, c, prec);
    }

    /* the direct groups; the columns n > M:
       sum c_n Z_n <= sum_(n > M) K0^(-2n) (1 + K0 / (2n - 1)) / n,
       with K0^(-2(M+1)) < 2^-F */
    tail = ldexp(1.0, (int) ceil(F - 2.0 * (M + 1) * _log2_d((double) K0) + 1e-6));
    tail *= (1.0 + (double) K0 / (2.0 * M + 1.0)) / ((M + 1.0) * (1.0 - 1.0 / (K0 * K0)));
    if (w.nd == 0)
    {
        mp_real_zero(c);
        _mp_real_add_error_ulps_at(c, 1.31 * tail * MP_REAL_EPS, -W);
        mp_real_add(s, s, c, prec);
    }
    for (g = 0; g < w.nd; g++)
    {
        if (g == 0)
            w.dA[g].err += tail;
        _khinchin_abel_get(c, w.dA + g, (n1 > 1) ? Alo : NULL, prec);
        mp_real_add(s, s, c, prec);
    }

    mp_real_add(s, s, &w.rows, prec);
    mp_real_div(s, s, &w.log2, prec);
    mp_real_exp_bits(res, s, FLINT_BITS * (n + 1));

    for (g = 0; g < w.nd; g++)
        _khinchin_abel_clear(w.dA + g);
    flint_free(w.dA);
    flint_free(w.da);
    flint_free(ks);
    flint_free(nhi);
    if (w.Y != NULL)
    {
        for (j = 0; j <= w.b; j++)
            mp_real_clear(w.Y + j);
        flint_free(w.Y);
    }
    mp_real_clear(&w.y);
    for (i = 0; i < 3 * (w.nz + 1); i++)
        mp_real_clear(w.zres + i);
    flint_free(w.zres);
    flint_free(w.zlo);
    for (g = 0; g < w.ns; g++)
    {
        _khinchin_abel_clear(w.sA + g);
        mp_real_clear(w.sAlo + g);
    }
    flint_free(w.sA);
    flint_free(w.sAlo);
    flint_free(w.slo);
    mp_real_clear(&w.rows);
    mp_real_clear(&w.log2);
    flint_free(cum);
    mp_real_clear(s);
    mp_real_clear(c);
    mp_real_clear(Alo);
}

/* the theta splitting (const_khinchin_afe.c) costs p^(3/2 + o(1)) against
   the p^(2 + o(1)) or more of the zeta series above; it is faster from
   somewhere between 1 and 2 million digits (depending on the machine and
   the thread count) */
void
_mp_real_const_khinchin_compute(mp_real_t res, slong n)
{
    if (n >= MP_REAL_KHINCHIN_AFE_CUTOFF)
        _mp_real_const_khinchin_afe(res, n, 0, 0);
    else
        _mp_real_const_khinchin_tune(res, n, 0, 0, 0);
}
