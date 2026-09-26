/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "mp_real.h"
#include "impl.h"

/* k-th roots and reciprocal k-th roots (k >= 4; the square and cube
   roots have their own fixed-point code) by a high-order iteration
   written in mp_real arithmetic, which handles the exponents and the
   propagation of the arithmetic errors: with z any approximation of
   x^(-1/k) and u = 1 - x z^k,

       x^(-1/k) = z (1 - u)^(-1/k),
       x^(1/k)  = S (1 - u)^(-(k-1)/k),    S = x z^(k-1),

   exactly, and the binomial series (1 - u)^(-b) = sum_j c_j u^j with
   c_0 = 1, c_j = c_(j-1) (b + j - 1) / j (decreasing, since b < 1)
   truncated after u^(r-1) has a tail below 2 c_r |u|^r for |u| <= 1/2.
   So z accurate to a fraction 1/r of the precision, obtained by a
   recursive call whose radius is then discarded (only the size of the
   final u matters, and that is computed rigorously), gives the root
   to the full precision with one evaluation of S, u and r - 1 series
   terms, each at the precision it contributes at. */

/* res = x^m by binary powering at precision n, m >= 1 */
static void
_pow_ui(mp_real_t res, const mp_real_t x, ulong m, slong n)
{
    if (m == 1)
        mp_real_set(res, x);
    else
    {
        mp_real_t t;
        mp_real_init(t);
        _pow_ui(t, x, m / 2, n);
        mp_real_mul(t, t, t, n);
        if (m % 2)
            mp_real_mul(t, t, x, n);
        mp_real_swap(res, t);
        mp_real_clear(t);
    }
}

/* the highest order of a step (the tuned orders are 3 and 4) */
#define ROOT_MAX_ORDER 16

/* the coefficients of the series of a step of order r for the
   exponent -bk/k: c_j = prod_{i<j} (bk + i k) / (k^j j!) over a common
   denominator D -- the lcm of the denominators of the c_j in lowest
   terms when the numerators, the denominators and the lcm are words
   (reduced = 1; for the square roots D is then a power of two, 8 at
   order 3 and 16 at order 4, and the division a shift), else
   D = k^(r-1) (r-1)! (Dfits = 1 when it is a word; else (r-1)! and
   r - 1 factors k are divided out one at a time), with D c_j =
   k^(r-1-j) ((r-1)!/j!) prod_{i<j} (bk + i k) as words when they fit
   (cfits[j]) and computed as exact balls otherwise.  Built once per
   root and shared by the steps of the recursion (the reductions cost
   a few hundred cycles, visible at one limb). */
typedef struct
{
    int r, reduced, Dfits;
    ulong k, bk, D, fac;
    ulong c[ROOT_MAX_ORDER];
    int cfits[ROOT_MAX_ORDER];
} root_coeffs_struct;

/* the last tables built, per thread: a program takes many roots of
   the same degree, and the reductions cost a few hundred cycles */
#define ROOT_COEFFS_CACHE 4
static FLINT_TLS_PREFIX root_coeffs_struct root_coeffs_cache[ROOT_COEFFS_CACHE];
static FLINT_TLS_PREFIX int root_coeffs_cache_next = 0;

static void
_root_coeffs(root_coeffs_struct * rc, ulong k, ulong bk, int r)
{
    slong i, j;
    ulong hi, D, g;
    int allfit;

    for (i = 0; i < ROOT_COEFFS_CACHE; i++)
    {
        if (root_coeffs_cache[i].r == r && root_coeffs_cache[i].k == k
            && root_coeffs_cache[i].bk == bk)
        {
            *rc = root_coeffs_cache[i];
            return;
        }
    }

    rc->r = r;
    rc->k = k;
    rc->bk = bk;

    /* the generic common denominator and coefficients, as words when
       they fit */
    rc->fac = 1;
    for (j = 2; j < r; j++)
        rc->fac *= (ulong) j;        /* (r-1)!, r small */
    D = rc->fac;
    rc->Dfits = 1;
    for (j = 1; j < r; j++)
    {
        umul_ppmm(hi, D, D, k);
        if (hi != 0)
            rc->Dfits = 0;
    }
    rc->D = D;
    allfit = rc->Dfits;
    for (j = 1; j < r; j++)
    {
        ulong c = 1;
        int cfits = 1;
        for (i = j + 1; i < r && cfits; i++)
        {
            umul_ppmm(hi, c, c, (ulong) i);
            cfits = (hi == 0);
        }
        for (i = 0; i < r - 1 - j && cfits; i++)
        {
            umul_ppmm(hi, c, c, k);
            cfits = (hi == 0);
        }
        for (i = 0; i < j && cfits; i++)
        {
            umul_ppmm(hi, c, c, bk + (ulong) i * k);
            cfits = (hi == 0);
        }
        rc->c[j] = c;
        rc->cfits[j] = cfits;
        allfit = allfit && cfits;
    }

    /* in lowest terms: D and the coefficients over their common gcd
       (D becomes the lcm of the denominators of the c_j) */
    rc->reduced = allfit;
    if (allfit)
    {
        g = D;
        for (j = 1; j < r && g > 1; j++)
            g = n_gcd(g, rc->c[j]);
        if (g > 1)
        {
            rc->D = D / g;
            for (j = 1; j < r; j++)
                rc->c[j] /= g;
        }
    }

    root_coeffs_cache[root_coeffs_cache_next] = *rc;
    root_coeffs_cache_next = (root_coeffs_cache_next + 1) % ROOT_COEFFS_CACHE;
}

/* one step of order r from z0 (taken as exact): res = the root
   (recip = 0) or the reciprocal root (recip = 1) of v to n limbs.
   rz0 is a bound on the relative error of z0, as a power of two
   (from the radius of the ball it came from), or unknown (a large
   value): known, |1 - v z0^k| < 2^(kb + 1 + rz0) and only the window
   of the product S z0 that reaches the residual is computed. */
static void
_root_step(mp_real_t res, const mp_real_t v, const mp_real_t z0, ulong k,
    slong n, int r, int recip, slong rz0, const root_coeffs_struct * rc)
{
    mp_real_t S, u, T;
    slong p = n + 2, uexp, j, cp, kb = FLINT_BIT_COUNT(k);
    ulong bk = recip ? 1 : k - 1;   /* (1 - u)^(-bk/k) */
    int zk_first;

    mp_real_init(S);
    mp_real_init(u);
    mp_real_init(T);

    /* the residual u = 1 - v z0^k, and the base S = v z0^(k-1) of the
       root.  For the root: u = 1 - S z0, one more product by the
       full-size S (a window of it when the size of the residual is
       known from z0's radius).  For the reciprocal root S is not
       needed: z0^k by binary powering first (z0 is short, a fraction
       1/r of the precision, so its powers are short squarings) and
       one product by v -- against (v z0^(k-1)) z0, one product by the
       full-size v instead of two, and for a short v (rsqrt(2)) no
       full-size product at all. */
    mp_real_set_ui(T, 1);
    if (!recip)
    {
        _pow_ui(S, z0, k - 1, p);
        mp_real_mul(S, S, v, p);
    }
    /* a short v (sqrt(2)): the root too takes its residual from z0^k
       (a squaring of the short z0 for k = 2, then a short product),
       where S z0 would be a product of two full-size operands */
    zk_first = recip || v->size <= p / 16;
    if (zk_first)
        _pow_ui(T, z0, k, p);
    if (rz0 < -8)
    {
        /* |1 - v z0^k| below 2k 2^rz0, plus the rounding of the
           factor taken at its midpoint */
        slong ub = FLINT_MAX(kb + 1 + rz0, mp_real_rel_radius_lt_2exp_si(zk_first ? T : S) + 1) + 1;
        slong E = (ub >= 0) ? ub / FLINT_BITS + 1 : -((-ub) / FLINT_BITS) + 1;
        /* u is needed to p limbs of the base, that is to p - |ub|/B
           limbs of its own: only that much of the window is computed */
        slong pu = FLINT_MAX(2, p + ub / FLINT_BITS);
        mp_real_set_ui(u, 1);
        if (zk_first)
            _mp_real_submul_bounded(u, u, v, T, E, pu);
        else
            _mp_real_submul_bounded(u, u, S, z0, E, pu);
    }
    else
    {
        if (zk_first)
            mp_real_mul(u, T, v, p);
        else
            mp_real_mul(u, S, z0, p);
        mp_real_set_ui(T, 1);
        mp_real_sub(u, T, u, p);
    }
    uexp = mp_real_abs_bound_lt_2exp_si(u);

    if (uexp > -2)
        flint_throw(FLINT_ERROR, "mp_real_root_ui: iteration failed\n");

    /* the series: res = base (1 + sum_{j<r} c_j u^j) as
       base + base W / D with W = sum_{j<r} (D c_j) u^j over the common
       denominator D of the table rc (see _root_coeffs): the powers of
       u come one product each by squaring (u^(2m) = (u^m)^2,
       u^(2m+1) = u^(2m) u; rectangular splitting would pay only for
       orders far above the tuned 3 and 4) at the precision they
       contribute at (the term u^j is 2^(j uexp) the size of the
       base), the division by D acts on W (a fraction 1 - 1/r of the
       precision) as a shift, one mp_real_div_ui or the chain of the
       generic denominator, and one full product applies the base. */
    {
        const mp_real_struct * base = recip ? z0 : S;
        mp_real_struct pw[ROOT_MAX_ORDER];
        mp_real_t W;
        slong cw = FLINT_MAX(2, p + uexp / FLINT_BITS + 2), i;
        root_coeffs_struct rc_local;

        FLINT_ASSERT(rc == NULL || (rc->r == r && rc->k == k && rc->bk == bk));
        if (rc == NULL)
        {
            _root_coeffs(&rc_local, k, bk, r);
            rc = &rc_local;
        }

        mp_real_init(W);
        for (j = 1; j < r; j++)
        {
            const mp_real_struct * uj;

            cp = FLINT_MAX(2, p + (j * uexp) / FLINT_BITS + 2);
            if (j == 1)
                uj = u;
            else
            {
                mp_real_init(pw + j);
                if (j % 2 == 0)
                    mp_real_mul(pw + j, (j / 2 == 1) ? u : pw + j / 2,
                        (j / 2 == 1) ? u : pw + j / 2, cp);
                else
                    mp_real_mul(pw + j, pw + j - 1, u, cp);
                uj = pw + j;
            }

            if (rc->cfits[j])
            {
                if (j == 1)
                    mp_real_mul_ui(W, uj, rc->c[j], cp);
                else
                    mp_real_addmul_ui(W, W, uj, rc->c[j], cw);
                continue;
            }
            else
            {
                /* D c_j = k^(r-1-j) ((r-1)!/j!) prod_{i<j} (bk + i k)
                   exactly (up to 20 limbs: k < MP_REAL_ROOT_K_MAX, r < 16) */
                mp_real_t coef;
                mp_real_init(coef);
                mp_real_set_ui(coef, 1);
                for (i = j + 1; i < r; i++)
                    mp_real_mul_ui(coef, coef, (ulong) i, 32);
                for (i = 0; i < r - 1 - j; i++)
                    mp_real_mul_ui(coef, coef, k, 32);
                for (i = 0; i < j; i++)
                    mp_real_mul_ui(coef, coef, bk + (ulong) i * k, 32);
                mp_real_mul(T, uj, coef, cp);
                mp_real_clear(coef);
            }
            if (j == 1)
                mp_real_swap(W, T);
            else
                mp_real_add(W, W, T, cw);
        }
        for (j = 2; j < r; j++)
            mp_real_clear(pw + j);

        if (rc->Dfits && rc->D == 1)
            ;
        else if (rc->Dfits && (rc->D & (rc->D - 1)) == 0)
            mp_real_mul_2exp_si(W, W, -(slong) FLINT_BIT_COUNT(rc->D) + 1);
        else if (rc->Dfits)
            mp_real_div_ui(W, W, rc->D, cw);
        else
        {
            mp_real_div_ui(W, W, rc->fac, cw);
            for (j = 1; j < r; j++)
                mp_real_div_ui(W, W, k, cw);
        }
        mp_real_mul(T, base, W, p);
        mp_real_add(res, base, T, p);
        mp_real_clear(W);

        /* the tail: below 2 c_r |u|^r |base| with c_r <= 1 */
        mp_real_add_error_2exp_si(res, r * uexp + mp_real_abs_bound_lt_2exp_si(base) + 1);
    }

    mp_real_clear(S);
    mp_real_clear(u);
    mp_real_clear(T);
}

/* the seed: v^(-1/k) for v = md 2^ep, md in [1/2, 1), 0 <= ep < k,
   from a double (through the logarithm, as v itself reaches 2^k),
   accurate to 2^-49, as an exact ball of 64 bits */
static void
_seed(mp_real_t z, double md, slong ep, ulong k)
{
    double zd = exp2(-((double) ep + log2(md)) / (double) k);
    int e0;

    zd = frexp(zd, &e0);
#if FLINT_BITS == 64
    {
        ulong mm = (ulong) d_mul_2exp_inrange(zd, FLINT_BITS);
        _mp_real_set_mpn_2exp(z, &mm, 1, (slong) e0 - FLINT_BITS);
    }
#else
    {
        /* all 53 bits of the double: two limbs (zd 2^64 is an integer
           below 2^64, split exactly) */
        ulong mm[2];
        double t = d_mul_2exp_inrange(zd, 2 * FLINT_BITS);
        mm[1] = (ulong) d_mul_2exp_inrange(t, -FLINT_BITS);
        mm[0] = (ulong) (t - d_mul_2exp_inrange((double) mm[1], FLINT_BITS));
        _mp_real_set_mpn_2exp(z, mm, 2, (slong) e0 - 2 * FLINT_BITS);
    }
#endif
}

/* the accuracy in bits needed from the approximation before a step of
   order r to n limbs: the step leaves about r (a - kb - 2) bits */
static slong
_need(slong n, ulong kb, int r)
{
    return (FLINT_BITS * n + 16 + r * (kb + 2) + r - 1) / r;
}

/* z ~= v^(-1/k) taken as exact (radius zeroed), accurate to need bits,
   v in [1/2, 2^k): the seed when that suffices, otherwise a step of
   order r from the recursion (or from the seed with the order that
   reaches need bits when the recursion would not shorten it) */
static void
_rroot_rec(mp_real_t z, slong * rz, const mp_real_t v, double md, slong ep,
    ulong k, slong kb, slong need, int r, const root_coeffs_struct * rc)
{
    slong n = (need + 8 + FLINT_BITS - 1) / FLINT_BITS;   /* limbs */
    slong need0 = _need(n, kb, r), rz0;
    mp_real_t z0;

    if (need <= 49)
    {
        /* the seed's accuracy is not a rigorous bound: unknown */
        _seed(z, md, ep, k);
        *rz = 0;
        return;
    }

    mp_real_init(z0);

    if (need0 <= 49 || need0 >= need)
    {
        /* from the seed, with the order that reaches need bits */
        int r0 = (int) ((need + (49 - kb - 3) - 1) / (49 - kb - 3));
        r0 = FLINT_MAX(r0, 2);
        _seed(z0, md, ep, k);
        _root_step(z, v, z0, k, n, r0, 1, 0, NULL);
    }
    else
    {
        _rroot_rec(z0, &rz0, v, md, ep, k, kb, need0, r, rc);
        _root_step(z, v, z0, k, n, r, 1, rz0, rc);
    }

    /* the radius, as a relative bound, outlives the ball's */
    *rz = mp_real_rel_radius_lt_2exp_si(z);
    z->err = 0;
    mp_real_clear(z0);
}

/* x^(1/k) (recip = 0) or x^(-1/k) (recip = 1) to n limbs by steps of
   order r (r = 0: the tuned default, 3 for k < 8 and 4 above --
   measured: within 5% of each other for k = 3, the fourth order 5-15%
   faster for k = 7 at all sizes) */
void
_mp_real_root_ui_order(mp_real_t res, const mp_real_t x, ulong k, slong n,
    int r, int recip)
{
    mp_real_t v, z;
    slong E, q, kb, need, rz;
    double md;

    if (x->size == 0 || x->negative)
        flint_throw(FLINT_ERROR, "mp_real_root_ui: need x > 0\n");
    if (k < 2 || k >= MP_REAL_ROOT_K_MAX)
        flint_throw(FLINT_ERROR, "mp_real_root_ui: k out of range\n");

    kb = FLINT_BIT_COUNT(k);
    if (r == 0)
        r = (kb <= 3) ? 3 : 4;
    if (r < 2 || r >= ROOT_MAX_ORDER)
        flint_throw(FLINT_ERROR, "mp_real_root_ui: order out of range\n");

    mp_real_init(v);
    mp_real_init(z);

    /* x = md 2^E with md in [1/2, 1); v = x 2^(-k q) with q = floor(E/k)
       lies in [1/2, 2^k) and the root scales back by 2^q */
    E = FLINT_BITS * x->exp;
    md = (double) x->d[x->size - 1];
    if (x->size >= 2)
        md += (double) x->d[x->size - 2] * MP_REAL_D_BINV;
    md = md * MP_REAL_D_BINV;
    q = (E >= 0) ? E / (slong) k : -(((slong) k - 1 - E) / (slong) k);

    mp_real_set(v, x);
    mp_real_mul_2exp_si(v, v, -(slong) k * q);

    /* the reciprocal root to a fraction of the precision, its radius
       discarded, then the final step */
    need = _need(n, kb, r);
    {
        root_coeffs_struct rc1, rc2;
        _root_coeffs(&rc1, k, 1, r);
        _rroot_rec(z, &rz, v, md, E - (slong) k * q, k, kb, need, r, &rc1);
        if (recip)
            _root_step(res, v, z, k, n, r, recip, rz, &rc1);
        else
        {
            _root_coeffs(&rc2, k, k - 1, r);
            _root_step(res, v, z, k, n, r, recip, rz, &rc2);
        }
    }

    mp_real_mul_2exp_si(res, res, recip ? -q : q);

    mp_real_clear(v);
    mp_real_clear(z);
}

void
mp_real_root_ui(mp_real_t res, const mp_real_t x, ulong k, slong n)
{
    if (k == 0)
        flint_throw(FLINT_ERROR, "mp_real_root_ui: k = 0\n");

    if (k == 1)
        mp_real_set(res, x);
    else if (k == 2)
    {
        if (res == x)
        {
            mp_real_t t;
            mp_real_init(t);
            mp_real_set(t, x);
            mp_real_sqrt(res, t, n);
            mp_real_clear(t);
        }
        else
            mp_real_sqrt(res, x, n);
    }
    else
        _mp_real_root_ui_order(res, x, k, n, 0, 0);
}

void
mp_real_rroot_ui(mp_real_t res, const mp_real_t x, ulong k, slong n)
{
    if (k == 0)
        flint_throw(FLINT_ERROR, "mp_real_rroot_ui: k = 0\n");

    if (k == 2)
    {
        if (res == x)
        {
            mp_real_t t;
            mp_real_init(t);
            mp_real_set(t, x);
            mp_real_rsqrt(res, t, n);
            mp_real_clear(t);
        }
        else
            mp_real_rsqrt(res, x, n);
    }
    else if (k == 1)
    {
        mp_real_t one;
        mp_real_init(one);
        mp_real_set_ui(one, 1);
        mp_real_div(res, one, x, n);
        mp_real_clear(one);
    }
    else
        _mp_real_root_ui_order(res, x, k, n, 0, 1);
}
