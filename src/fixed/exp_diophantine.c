/*
    Copyright (C) 2022, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fixed.h"

/* exp by diophantine (multi-prime) argument reduction (port of
   arb_exp_arf_log_reduction, src/arb/log_reduce.c):

       exp(x) = 2^c_0 (p / q) exp(t),    t = x - sum_j c_j log(p_j),

   with p_j the primes 2, 3, ..., 41 (p_0 = 2 is free: its power is
   a shift), integer coefficients c_j chosen by the reduction below,
   p = prod_{c_j > 0} p_j^(c_j), q = prod_{c_j < 0} p_j^(-c_j), and
   the reduced argument t made tiny (some 2^-40 at a couple of
   thousand bits, down to about 2^-186 at a million bits and beyond,
   the depth of the relation table).  The work beyond exp(t) --
   whose series now needs a small fraction of the terms of a
   bitwise-reduced argument -- is one long-by-short multiplication
   by p, one long-by-short division by q, and the products p and q
   themselves.

   THE REDUCTION (_fixed_log_reduce) is a greedy descent through a
   precomputed table of integer relations: row i holds coefficients
   d_i with sum_j d_ij log(p_j) = epsilon_i, |epsilon_i| decreasing
   down the table (output of _arb_log_precompute_reductions, coefficient
   bound 2^8 per row; rel_tab.c).
   After the free log 2 step, each row subtracts the nearest integer
   multiple of epsilon_i from the residual, keeping a running tally
   of the coefficients and stopping when the total "weight"
   sum_{j > 0} |c_j| log2(p_j) / log(2) -- a proxy for the size of
   p q -- would exceed max_weight (the working precision in bits by
   default, as in arb).  The residual is tracked in double precision
   and resynchronized every eight rows from a short (256 to 768
   bit) fixed-point evaluation of x - sum c_j log(p_j), whose
   two's complement sign limb carries the sign.

   SIGN OF THE RESIDUAL.  Rounding to nearest leaves t of either
   sign, but fixed_exp_reduced takes t >= 0.  Rather than paying for
   a reciprocal, a negative t is lifted by adding one relation with
   a positive epsilon -- the last such row the descent reached --
   which leaves 0 < t <= epsilon_j, at most a factor two (a few bits
   when the last positive row sits a step or two back) worse than
   the nearest-rounding residual.

   PRIME POWERS (_fixed_prime_powers).  p and q are built by vector
   exponentiation on the bits of the exponent vector,

       P(c) = P(floor(c / 2))^2 * prod_{j : c_j odd} p_j,

   one squaring per exponent bit with the odd-subset product a
   single limb (the twelve odd primes multiply to 1.5 10^14 < 2^48;
   on 32-bit limbs the accumulator is flushed as needed), so the
   whole product costs about one squaring at the final size (the
   sizes double up the chain) instead of the several full
   multiplications of a product tree over separately computed
   powers.  Since p / q 2^c_0 = exp(x) up to the residual, p and q
   have the same size, sum_{j>0} |c_j| log2 p_j / 2 bits each, at
   most 0.35 max_weight: 0.15-0.3 times the precision at arb's
   budget, up to 1.4 at four times it and 9 at 32 times.  Beyond
   the working precision the products are carried TRUNCATED: the
   chain caps its mantissa at wn + 3 limbs (the top levels then
   cost one sqrhigh at the cap each), E p is a windowed middle
   product, and the division sees at most wn + 3 divisor limbs.

   WORKING PRECISION AND ERROR.  Everything after the reduction runs
   at wn = n + G limbs.  The cached logarithms are one-sided floors,
   so the dot product is short of the exact t by less than
   sum |c_j| + 1 working ulps, a relative error of exp(t) of the same
   size; fixed_exp_reduced (or fixed_exp_notab when the residual is
   too large for it) adds its documented budget; the multiplication
   is exact, the division is a floor; the factor 2^c_0 p / q < e
   amplifies all of it by less than 3.  G is chosen (one limb in
   practice) so that this total stays below one output ulp, and the
   final truncation adds one more:
   FIXED_EXP_DIOPHANTINE_MAX_ERR = 3.

   SCOPE.  This version takes x in [0, 1) and returns exp(x) in
   [1, e) with n fraction limbs and a unit limb; the reduction of a
   general argument by log 2 (arb_exp_arf_huge) is left for the
   caller for now.

   STATUS (measured on the development VM, single core; see
   tune/tune-exp-diophantine.c).  The reduction is bit-identical to
   arb's (0 differing coefficient vectors in 3600 random inputs
   across 1k-4M bits) and reaches the same residual depth to within
   a bit; the fixed-point steps around it are 2-3x faster than their
   arb counterparts, but the reduced exponential is 93-95% of the
   time in both, so the whole is only 1.05-1.4x faster than
   arb_exp_arf_log_reduction.  Against fixed_exp_bitwise_rs it is
   1.1-1.6x SLOWER per call at every size measured (n = 40 .. 10000
   limbs): the relation table buys r = 60-150 leading zero bits
   (the 13-prime table bottoms out at 2^-186), the bitwise table
   r = 768 at a cost linear in r n.  Its place is the other side of
   the trade: the precomputation is one logarithm per prime, 3-10x
   cheaper than the bitwise table at the same precision (2.8 ms vs
   30 ms at 32k bits, 30 ms vs 310 ms at 171k bits with 13 primes),
   so the bitwise scheme only amortizes its table after some 50-400
   evaluations.  Tuning: 20 primes are consistently the fastest per
   call (r ~ 120-150; 10-15% ahead of 13) but cost 2.7x the
   precomputation through arb's size-classed logarithm vector; the
   weight budget is flat between 0.5 and 4 times the precision (arb's
   choice, the precision, is fine) and 32x is always a loss to the
   division by q. */

/* ---- signed fixed-point helpers ------------------------------------ */

/* acc (len limbs, two's complement, radix point below limb len - 1)
   = base - sum_j rel_j alpha_j, with alpha_j the len-limb angle at
   alpha + j stride (len - 1 fraction limbs, unit limb on top).  Arithmetic is modulo B^len: the intermediate sums may wrap,
   but the final value, of magnitude far below B/2, is recovered
   exactly by the two's complement reading. */
void
_fixed_log_dot(nn_ptr acc, nn_srcptr base, slong len, const slong * rel,
    slong num, nn_srcptr alpha, slong stride)
{
    slong j;

    flint_mpn_copyi(acc, base, len);

    for (j = 0; j < num; j++)
    {
        if (rel[j] > 0)
            mpn_submul_1(acc, alpha + j * stride, len, (ulong) rel[j]);
        else if (rel[j] < 0)
            mpn_addmul_1(acc, alpha + j * stride, len, (ulong) (-rel[j]));
    }
}

/* the value of (a, len) read as a two's complement fixed-point
   number with len - 1 fraction limbs, to double precision:
   the magnitude TRUNCATED to 53 bits, exactly as fmpz_get_d would
   read it (arb's reduction compares against this rounding, so the
   descent takes bit-identical steps); tmp holds len limbs */
double
_fixed_signed_get_d(nn_srcptr a, slong len, nn_ptr tmp)
{
    int neg = (a[len - 1] >> (FLINT_BITS - 1)) != 0;
    nn_srcptr m = a;
    slong k, bits;
    ulong hi;
    unsigned int bc;
    double d;

    if (neg)
    {
        mpn_neg(tmp, a, len);
        m = tmp;
    }

    for (k = len - 1; k >= 0 && m[k] == 0; k--)
        ;

    if (k < 0)
        return 0.0;

    /* the top FLINT_BITS bits of the magnitude, leading bit on top */
    bc = FLINT_BIT_COUNT(m[k]);
    if (bc == FLINT_BITS)
        hi = m[k];
    else
        hi = (m[k] << (FLINT_BITS - bc))
            | ((k > 0) ? (m[k - 1] >> bc) : UWORD(0));
    bits = FLINT_BITS * k + (slong) bc;   /* leading bit position + 1 */

#if FLINT_BITS == 64
    d = (double) (hi >> 11);              /* 53 bits, exact */
    d = ldexp(d, (int) (bits - 53 - FLINT_BITS * (len - 1)));
#else
    d = (double) hi;
    if (k > 0)
    {
        /* the top 53 bits as an integer: hi (32 bits) followed by the
           top 21 of the next limb, so hi scales by 2^21, not 2^32 */
        ulong lo = (bc == FLINT_BITS) ? m[k - 1]
            : ((m[k - 1] << (FLINT_BITS - bc))
                | ((k > 1) ? (m[k - 2] >> bc) : UWORD(0)));
        d = ldexp(d, 21) + (double) (lo >> 11);
        d = ldexp(d, (int) (bits - 53 - FLINT_BITS * (len - 1)));
    }
    else
        d = ldexp(d, (int) (bits - 32 - FLINT_BITS * (len - 1)));
#endif

    return neg ? -d : d;
}

/* ---- the reduction --------------------------------------------------- */

/* Port of _arb_log_reduce_fixed.  Chooses integers rel[0..num-1]
   making x - sum_j rel_j alpha_j small, alpha_0 assumed free (its
   coefficient carries no weight).  x is a fixed-point fraction with
   wr limbs; the alpha_j are read from the logarithm cache at wr
   fraction limbs: entry j is the wr + 1 limbs at alpha + j stride
   (unit limb on top).  tab is the relation table (rel_tab.c):
   row i has sum_j d_ij alpha_j = epsilon_i.

   The descent also stops at the first row with |epsilon_i| below
   eps_min: rows finer than the precision at which the residual will
   be evaluated cannot reduce the argument, and would leave the
   residual's sign to the evaluation noise.

   Returns the index of the last row with a positive epsilon among
   those the descent passed (0 if it broke out before any), i.e. the
   row to add should the exact residual turn out negative. */
slong
_fixed_log_reduce(slong * rel, const fixed_rel_struct * tab,
    nn_srcptr x, slong wr, double max_weight, double eps_min,
    nn_srcptr alpha, slong stride)
{
    const short * d = tab->d;
    const double * epsilon = tab->epsilon;
    const double * epsilon_inv = tab->epsilon_inv;
    const float * weights = tab->weights;
    slong num = tab->num;
    slong i, j, n, last = -1, jfix;
    slong * new_rel;
    const short * d_row;
    double dalpha, weight, dx;
    nn_ptr acc, base, tmp;
    TMP_INIT;

    TMP_START;
    acc = TMP_ALLOC(3 * (wr + 1) * sizeof(ulong));
    base = acc + (wr + 1);
    tmp = base + (wr + 1);
    new_rel = TMP_ALLOC(num * sizeof(slong));

    for (j = 0; j < num; j++)
        rel[j] = 0;

    flint_mpn_copyi(base, x, wr);
    base[wr] = 0;

    dx = _fixed_signed_get_d(base, wr + 1, tmp);

    /* reduce by the first alpha, which is assumed to be free */
    dalpha = _fixed_signed_get_d(alpha, wr + 1, tmp);
    n = (slong) floor(dx / dalpha + 0.5);
    dx -= dalpha * n;
    rel[0] = n;

    /* recompute accurately if there is significant cancellation */
    if (FLINT_ABS(n) > 10)
    {
        _fixed_log_dot(acc, base, wr + 1, rel, num, alpha, stride);
        dx = _fixed_signed_get_d(acc, wr + 1, tmp);
    }

    for (i = 0; ; i++)
    {
        d_row = d + i * num;

        if (i >= tab->rows || fabs(epsilon[i]) < eps_min)
            break;

        for (j = 0; j < num; j++)
            new_rel[j] = d_row[j];

        /* a multiple beyond any sensible weight budget (row entries
           are below 2^16, so with the weight clamp this keeps rel and
           the products n d within the slong range on 32- and 64-bit
           limbs): with a deep table and a tiny argument dx / epsilon_i
           can be astronomically large -- arb never reduces such
           arguments -- and the rows below would need even larger
           multiples, so stop here */
        if (fabs(dx * epsilon_inv[i]) > ldexp(1.0, FLINT_BITS / 2 - 2))
            break;

        n = (slong) floor(dx * epsilon_inv[i] + 0.5);

        if (n != 0)
        {
            weight = 0.0;
            for (j = 0; j < num; j++)
            {
                new_rel[j] = rel[j] + n * new_rel[j];
                if (j != 0)
                    weight += FLINT_ABS(new_rel[j]) * weights[j] * 1.442695;
            }

            if (weight > max_weight)
                break;

            for (j = 0; j < num; j++)
                rel[j] = new_rel[j];

            dx -= n * epsilon[i];
        }

        last = i;

        if (i % 8 == 7)
        {
            _fixed_log_dot(acc, base, wr + 1, rel, num, alpha, stride);
            dx = _fixed_signed_get_d(acc, wr + 1, tmp);
        }
    }

    for (jfix = last; jfix > 0 && epsilon[jfix] <= 0.0; jfix--)
        ;
    if (jfix < 0)
        jfix = 0;

    /* Sign fix.  The descent above is arb's, step for step; its
       nearest rounding leaves a residual of either sign with
       |t| <= |epsilon_last| / 2.  A negative one is lifted by the
       last positive row reached (t in (0, epsilon_jfix]) and then
       driven down again through the rows below it with FLOOR
       rounding -- n = sign(epsilon_i) floor(t / |epsilon_i|) keeps
       t in [0, |epsilon_i|) -- ending below |epsilon_last|, within
       a bit of the nearest-rounding residual, while the lift alone
       would cost several bits whenever the last positive row sits a
       few rows back. */
    _fixed_log_dot(acc, base, wr + 1, rel, num, alpha, stride);
    dx = _fixed_signed_get_d(acc, wr + 1, tmp);

    if (dx < 0.0 && last >= 0)
    {
        d_row = d + jfix * num;
        for (j = 0; j < num; j++)
            rel[j] -= d_row[j];
        dx += epsilon[jfix];

        for (i = jfix + 1; i <= last; i++)
        {
            double ae = fabs(epsilon[i]);

            d_row = d + i * num;
            if (fabs(dx / ae) > ldexp(1.0, FLINT_BITS / 2 - 2))
                break;
            n = (slong) floor(dx / ae);
            if (epsilon[i] < 0.0)
                n = -n;

            if (n != 0)
            {
                weight = 0.0;
                for (j = 0; j < num; j++)
                {
                    new_rel[j] = rel[j] + n * d_row[j];
                    if (j != 0)
                        weight += FLINT_ABS(new_rel[j]) * weights[j] * 1.442695;
                }

                if (weight > max_weight)
                    break;

                for (j = 0; j < num; j++)
                    rel[j] = new_rel[j];

                dx -= n * epsilon[i];
            }

            if (i % 8 == 7)
            {
                _fixed_log_dot(acc, base, wr + 1, rel, num, alpha, stride);
                dx = _fixed_signed_get_d(acc, wr + 1, tmp);
            }
        }
    }

    TMP_END;
    return jfix;
}

/* ---- prime powers ----------------------------------------------------- */

/* res B^(*e_out) ~= prod_j primes[j]^e[j] over the j with e[j] > 0
   (an empty product is 1), by vector exponentiation on the exponent
   bits, with the mantissa capped at cap limbs: once a squaring would
   exceed the cap, the operand is truncated to cap limbs (the dropped
   limbs going into the exponent) and only the high limbs of the
   square are computed (flint_mpn_sqrhigh: a lower bound within
   2 cap ulps of its returned limb), so the top levels of the chain
   cost a cap-limb squaring however large the true product is.  The
   result is then a lower bound within about (levels 2 cap) B^-cap
   relative, which the caller's cap = wn + 3 puts far below one
   working ulp.  res and tmp each need 2 cap + 4 limbs.  Returns the
   mantissa size (top limb nonzero). */
static slong
_fixed_prime_powers(nn_ptr res, nn_ptr tmp, slong * e_out,
    const ulong * primes, const slong * e, slong len, slong cap)
{
    slong emax = 0, an, j, k, L, ex = 0;
    ulong m;
    nn_ptr a = res, b = tmp;

    for (j = 0; j < len; j++)
        emax = FLINT_MAX(emax, e[j]);

    if (emax == 0)
    {
        res[0] = 1;
        *e_out = 0;
        return 1;
    }

    L = FLINT_BIT_COUNT((ulong) emax);

    a[0] = 1;
    an = 1;

    for (k = L - 1; k >= 0; k--)
    {
        if (k != L - 1)
        {
            if (an > cap)
            {
                /* truncate to cap limbs and square the top limbs
                   only: cap limbs of the 2 cap-limb square */
                ulong lo;

                ex += an - cap;
                ex *= 2;   /* (a B^ex)^2 = a^2 B^(2 ex) */
                lo = flint_mpn_sqrhigh(b, a + (an - cap), cap);
                if (b[cap - 1] == 0)
                {
                    /* a zero top limb: shift in the returned limb
                       below the window to keep cap limbs */
                    flint_mpn_copyd(b + 1, b, cap - 1);
                    b[0] = lo;
                    ex += cap - 1;
                }
                else
                    ex += cap;
                an = cap;
            }
            else
            {
                /* b = a^2 in full, into the other buffer (the generic
                   squaring forbids aliasing), then truncated to cap
                   limbs */
                ex *= 2;
                flint_mpn_sqr(b, a, an);
                an = 2 * an - (b[2 * an - 1] == 0);
                if (an > cap)
                {
                    flint_mpn_copyi(b, b + (an - cap), cap);
                    ex += an - cap;
                    an = cap;
                }
            }
            a = (b == tmp) ? tmp : res;
            b = (a == tmp) ? res : tmp;
        }

        /* multiply in the primes whose exponent has bit k set,
           accumulated in a limb and flushed when it would overflow
           (it never does on 64-bit limbs) */
        m = 1;
        for (j = 0; j < len; j++)
        {
            if (e[j] > 0 && ((e[j] >> k) & 1))
            {
                ulong pj = primes[j];

                if (m > UWORD_MAX / pj)
                {
                    ulong cy = mpn_mul_1(a, a, an, m);
                    if (cy)
                        a[an++] = cy;
                    m = 1;
                }
                m *= pj;
            }
        }

        if (m > 1)
        {
            ulong cy = mpn_mul_1(a, a, an, m);
            if (cy)
                a[an++] = cy;
        }
        FLINT_ASSERT(an <= 2 * cap + 4);
    }

    if (a != res)
        flint_mpn_copyi(res, a, an);

    *e_out = ex;
    return an;
}

/* ---- the exponential ---------------------------------------------------- */

void
_fixed_exp_diophantine_tune(nn_ptr y, nn_srcptr x, slong n, slong num,
    double max_weight)
{
    const fixed_rel_struct * tab;
    slong * rel, * neg;
    slong xn = n, wr, wn, G, G2, j, jfix, r, pn, qn, cap, Nn, Qn, i;
    slong ep, eq, eN;
    ulong relsum;
    double eps_min;
    nn_ptr base, t, E, p, q, tmp, N, Q;
    TMP_INIT;

    FLINT_ASSERT(num >= 2 && num <= FIXED_REL_MAX);

    /* the weight bounds sum |rel_j| log2 p_j: keep the coefficients
       (and their products in the descent) within the slong range
       whatever budget is requested */
    max_weight = FLINT_MIN(max_weight, ldexp(1.0, FLINT_BITS - 11));

    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    if (xn == 0)
    {
        flint_mpn_zero(y, n);
        y[n] = 1;
        return;
    }

    tab = fixed_rel_table(0, num);

    /* precision of the reduction search: arb's ladder (256, 512,
       768 bits by the target precision), and in any case enough to
       resolve the deepest row of the table with margin */
    if (FLINT_BITS * n <= 10000)
        wr = 256 / FLINT_BITS;
    else if (FLINT_BITS * n <= 100000)
        wr = 512 / FLINT_BITS;
    else
        wr = 768 / FLINT_BITS;
    wr = FLINT_MAX(wr, (slong) ((-log2(tab->epsilon_min) + 80) / FLINT_BITS) + 1);

    _fixed_log_primes_ensure(num, FLINT_MAX(wr, n + 1));

    TMP_START;
    rel = TMP_ALLOC(2 * num * sizeof(slong));
    neg = rel + num;

    /* the top wr limbs of x (zero-padded below when x is shorter) */
    base = TMP_ALLOC(wr * sizeof(ulong));
    {
        slong pad = FLINT_MAX(wr - n, 0);
        for (i = 0; i < pad; i++)
            base[i] = 0;
        flint_mpn_copyi(base + pad, x + FLINT_MAX(n - wr, 0), wr - pad);
    }

    /* rows with epsilon below 2^-(64 n + 32) are useless: the residual
       is evaluated with (at least) one guard limb, whose dot-product
       error of sum |rel_j| ulps must stay well below the lifting
       relation (matters only for n <= 4; underflows to 0 beyond);
       likewise rows the search precision cannot resolve */
    eps_min = ldexp(1.0, -(int) FLINT_MIN(FLINT_BITS * n + 32, 2000));
    eps_min = FLINT_MAX(eps_min, ldexp(1.0, -(int) (FLINT_BITS * wr - 48)));

    jfix = _fixed_log_reduce(rel, tab, base, wr, max_weight, eps_min,
        _fixed_log_primes_entry(0, wr), _fixed_log_primes_n);

    /* guard limbs: the residual is short of the truth by at most
       sum |rel_j| + 1 working ulps of the cached floors (plus the
       lifting relation, at most twice that), the exponential adds
       up to 128, and 2^rel_0 p / q < e amplifies by under 3 */
    relsum = 0;
    for (j = 0; j < num; j++)
        relsum += (ulong) FLINT_ABS(rel[j]);
    G = (FLINT_BIT_COUNT(4 * relsum + 256) + 2 + 8 + FLINT_BITS - 1)
        / FLINT_BITS;
    wn = n + G;

    _fixed_log_primes_ensure(num, wn);   /* a no-op unless G > 1 */

    /* t = x B^G - sum_j rel_j log(p_j), two's complement in wn + 1
       limbs.  The reduction returns a nonnegative residual at its
       search precision; as a safety net against the evaluation
       error at the working precision, a negative t is lifted by the
       positive relation jfix (moving on to the previous positive
       row should that prove too small) */
    t = TMP_ALLOC(2 * (wn + 1) * sizeof(ulong));
    tmp = t + (wn + 1);
    flint_mpn_zero(tmp, G);
    flint_mpn_copyi(tmp + G, x, n);
    tmp[wn] = 0;

    _fixed_log_dot(t, tmp, wn + 1, rel, num,
        _fixed_log_primes_entry(0, wn), _fixed_log_primes_n);

    for (i = 0; (t[wn] >> (FLINT_BITS - 1)) != 0; i++)
    {
        const short * d_row;

        if (i >= 4 && jfix > 0)
        {
            for (jfix--; jfix > 0 && tab->epsilon[jfix] <= 0.0; jfix--)
                ;
            i = 0;
        }

        d_row = tab->d + jfix * num;
        for (j = 0; j < num; j++)
            rel[j] -= d_row[j];   /* t += epsilon_jfix */
        _fixed_log_dot(t, tmp, wn + 1, rel, num,
            _fixed_log_primes_entry(0, wn), _fixed_log_primes_n);
    }
    FLINT_ASSERT(t[wn] == 0);   /* 0 <= t < 1 */

    /* E = exp(t), wn fraction limbs and a unit limb */
    E = TMP_ALLOC((wn + 1) * sizeof(ulong));
    for (i = wn - 1; i >= 0 && t[i] == 0; i--)
        ;
    if (i < 0)
    {
        flint_mpn_zero(E, wn);
        E[wn] = 1;
    }
    else
    {
        r = FLINT_BITS * (wn - 1 - i)
            + (FLINT_BITS - FLINT_BIT_COUNT(t[i]));
        if (r >= 16)
            fixed_exp_reduced(E, t, wn, (flint_bitcnt_t) r, 0);
        else
            fixed_exp_notab(E, t, wn);
    }

    /* p = pm B^ep and q = qm B^eq, mantissas capped at cap = wn + 3
       limbs: beyond that the products only carry truncated top
       limbs, the exponents the rest */
    cap = wn + 3;
    p = TMP_ALLOC(4 * (2 * cap + 4) * sizeof(ulong));
    q = p + 2 * (2 * cap + 4);
    for (j = 0; j < num; j++)
        neg[j] = -rel[j];
    pn = _fixed_prime_powers(p, p + (2 * cap + 4), &ep, tab->primes + 1,
        rel + 1, num - 1, cap);
    qn = _fixed_prime_powers(q, q + (2 * cap + 4), &eq, tab->primes + 1,
        neg + 1, num - 1, cap);

    /* N = E p 2^rel_0 as Nm B^eN with wn fraction limbs: the full
       product when it fits within cap + 2 limbs, otherwise its top
       cap limbs by the windowed middle product (a lower bound short
       by at most min(wn + 1, pn) ulps of its second limb).  The
       power of two is a bit shift of the mantissa (the right shift
       floors, which composes exactly with the floor of the division)
       and a limb offset in the exponent. */
    {
        slong sh = FLINT_ABS(rel[0]), shq = sh / FLINT_BITS, ps;
        int shb = (int) (sh % FLINT_BITS);

        N = TMP_ALLOC((cap + 4 + wn + 2 + FLINT_MAX(pn, 0)) * sizeof(ulong));
        eN = ep;

        if (pn == 1 && p[0] == 1)
        {
            flint_mpn_copyi(N, E, wn + 1);
            Nn = wn + 1;
        }
        else if ((ps = wn + 1 + pn) <= cap + 2)
        {
            if (wn + 1 >= pn)
                flint_mpn_mul(N, E, wn + 1, p, pn);
            else
                flint_mpn_mul(N, p, pn, E, wn + 1);
            Nn = ps;
        }
        else
        {
            flint_mpn_mulmid(N, E, wn + 1, p, pn, ps - cap, ps);
            Nn = cap;
            eN += ps - cap;
        }
        while (Nn > 0 && N[Nn - 1] == 0)
            Nn--;

        if (rel[0] > 0)
        {
            if (shb)
            {
                N[Nn] = mpn_lshift(N, N, Nn, shb);
                Nn++;
            }
            eN += shq;
        }
        else if (rel[0] < 0)
        {
            if (shb && Nn > 0)
                mpn_rshift(N, N, Nn, shb);
            eN -= shq;
        }
        while (Nn > 0 && N[Nn - 1] == 0)
            Nn--;
    }

    /* exp(x) B^wn = (Nm / qm) B^(eN - eq): shift the numerator by
       the limb offset (a right shift floors, composing exactly with
       the division's floor) and divide.  Balanced shapes go through
       fixed_div_newton (2.4-2.9 M(n) with fft_small at every size
       measured, 1000 .. 150000 limbs); short divisors through the
       quotient-only flint_mpn_tdiv_q, which is ahead there. */
    {
        slong sft = eN - eq, numn, k;
        nn_ptr num;

        if (sft >= 0)
        {
            num = TMP_ALLOC((Nn + sft + 2) * sizeof(ulong));
            flint_mpn_zero(num, sft);
            flint_mpn_copyi(num + sft, N, Nn);
            numn = Nn + sft;
        }
        else
        {
            num = N + FLINT_MIN(-sft, Nn);
            numn = FLINT_MAX(Nn + sft, 0);
        }

        if (qn == 1 && q[0] == 1)
        {
            Q = num;
            Qn = numn;
            G2 = G;
        }
        else if (numn < qn)
        {
            Q = num;
            Qn = 0;
            G2 = G;
        }
        else if (3 * qn >= wn)
        {
            /* fixed_div_newton on the fractions num / B^numn and
               qm / B^qn: their quotient is exp(x) B^k with
               k = wn + qn - numn, which is 0 or -1 (num < 4 qm B^wn),
               so with wn + 3 fraction limbs exp(x)'s own B^-wn limb
               sits at index 3 + k; the error 4 B^-(wn+3) / a with
               a >= 1/B is under 4 B^-(wn+1) in exp(x) */
            k = wn + qn - numn;
            FLINT_ASSERT(k == 0 || k == -1);
            Q = TMP_ALLOC((wn + 6) * sizeof(ulong));
            fixed_div_newton(Q, num, numn, q, qn, wn + 3);
            Qn = wn + 5;
            G2 = G + 3 + k;
        }
        else
        {
            Q = TMP_ALLOC((numn - qn + 2) * sizeof(ulong));
            flint_mpn_tdiv_q(Q, num, numn, q, qn);
            Qn = numn - qn + 1;
            G2 = G;
        }
    }

    if (Qn >= wn + 1)
    {
        flint_mpn_copyi(y, Q + G2, n + 1);
    }
    else
    {
        /* the value rounded below B^wn with high limbs missing:
           only possible for x = 0 + epsilon with the tiny error
           pointing down */
        flint_mpn_zero(y, n + 1);
        if (Qn > G2)
            flint_mpn_copyi(y, Q + G2, Qn - G2);
    }

    TMP_END;
}

void
fixed_exp_diophantine(nn_ptr y, nn_srcptr x, slong n)
{
    _fixed_exp_diophantine_tune(y, x, n, 13, (double) (FLINT_BITS * n));
}
