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

/* sin and cos by diophantine (multi-prime) argument reduction (the trigonometric
   counterpart of exp_diophantine.c; port of
   arb_sin_cos_arf_atan_reduction with a different value side).

   THE REDUCTION.  With theta_0 = pi/2 and theta_j = 2 arg(pi_j) for
   the nonreal Gaussian primes pi_j = a_j + b_j i, the descent of
   _fixed_log_reduce over the relation table (rel_tab.c) picks
   integers c_j with

       x = c_0 pi/2 + sum_{j>0} c_j theta_j + t,

   t tiny.  The angles are DOUBLE arguments so that each rotation is
   a Gaussian rational with a rational integer denominator,
   e^(i theta_j) = pi_j / conj(pi_j) = pi_j^2 / p_j, p_j = N(pi_j).

   THE VALUE SIDE.  Let A = prod_{c_j > 0} pi_j^(c_j) prod_{c_j < 0}
   conj(pi_j)^(-c_j), a Gaussian integer; then

       e^(ix) = i^c_0 e^(it) A / conj(A) = i^c_0 e^(it) A^2 / N,
       N = |A|^2 = prod p_j^|c_j|,

   and the normalization is a division by the RATIONAL INTEGER N,
   which is moreover A_r^2 + A_i^2 exactly (two short squarings, no
   separate prime-power product).  arb divides by the complex number
   conj(P) Q instead (an acb_div); this is the same quantity, but a
   real divisor costs a fraction of a complex division.

   WHY NOT THE TANGENT.  The bitwise routine (tan_bitwise_rs.c) wins
   by iterating on the tangent: its rotation W = prod (1 + i 2^-k)
   has an irrational modulus, so normalizing e^(it) W would need
   |W|^2 and a square root, while the ratio Im/Re needs neither --
   at the price of three products TT DE, DE^2, TT^2 and a reciprocal
   of their sum in the reconstruction.  Here the modulus is known
   exactly and rational: normalizing costs one reciprocal of an
   integer and two multiplications, less than the ratio route's
   reconstruction, and the reduced argument's sine and cosine are
   needed anyway.  So the direct form is the right one for the
   diophantine (multi-prime) reduction; the tangent, if wanted, is the single
   division Z_i / Z_r of the same Z = e^(it) A^2, with no N at all.

   THE COMPUTATION.  t = x - sum c_j theta_j is formed at wn = n + G
   limbs in two's complement from the cached floors of the angles
   (atan_gauss.c); its sign only flips the sine of the residual, so
   no lifting is needed (unlike the exponential).  (s, g) =
   (sin |t|, 1 - cos |t|) come from fixed_sin_cos_reduced.  A is
   built by vector exponentiation on the exponent bits with
   flint_mpn_sqr_complex, the odd-subset product of Gaussian primes
   accumulated in a limb pair; A^2 is one more complex squaring;
   Z = A^2 ((1 - g) + i sigma s) one complex product
   (flint_mpn_mul_complex, three real products by Karatsuba); the
   two quotients Z_r / N and Z_i / N -- one reciprocal of N and a
   product each when both are wanted (15-28% faster than two
   divisions, measured 1024 .. 40000 limbs), one division for a
   single output -- are then rotated by i^c_0 (a permutation and
   sign change of the pair, the results being nonnegative for x in
   [0, 1)).

   SIZES.  With a weight budget of w times the precision, N has up
   to 0.72 w prec bits and A half of that.  Beyond cap = wn + 3
   limbs the Gaussian quantities are carried TRUNCATED with a limb
   exponent: the power chain reads the top cap limbs of both parts
   and squares them with flint_mpn_sqrhigh_n_complex, N and A^2
   likewise come from the high halves of the squares of A's top
   limbs, Z from flint_mpn_mulhigh_n_complex on the two pairs
   windowed to a common length, and the divisions shift by the
   accumulated exponents.  Each high product is within 3 ulps of its
   lowest returned limb, a relative 2^(-64 (cap - 1)) that the guard
   limb absorbs over the chain's ~20 levels.

   ERROR.  The cached angles are floors (sum |c_j| + 1 working ulps),
   the reduced pair is within 96 ulps each, Z is exact given (s, g),
   and the divisions floor; one guard limb puts the total below one
   output ulp.  FIXED_SIN_COS_DIOPHANTINE_MAX_ERR = 4.

   MEASURED (development VM, fft_small).  The reduced sine and cosine
   are 90% of the time: at 2677 limbs fixed_sin_cos_reduced costs
   11.8 ms at r = 105 against 7.5 ms for fixed_exp_reduced, and the
   two only meet from r ~ 300 -- so, unlike the exponential, the
   trigonometric reduction wants to be deep, and the value side
   (Gaussian powers, N and A^2, Z, the two divisions: 1-2 ms) is a
   small fraction.  Sweeping primes x weight budget with the capped
   products (us per call, min of 3; the bitwise routine for
   reference):

       n = 1024, bitwise 1656:  13 primes 2600-2900; 20: 2350-2600;
                                32: 1970-2380 (r 128..263);
                                48: 1810-2030
       n = 2677, bitwise 5917:  13: 10500-15000; 20: 8700-10300;
                                32: 6830 (w = 8, r = 295) .. 7830;
                                48: 6830-9000
       n = 10000, bitwise 38.0 ms:  32 primes 52-59 ms, flat in w

   hence 32 primes at four times the precision as weight (arb's
   choice, 13 at half the precision, is 1.5-2x slower here).  Per
   call this is 1.1-1.4x the bitwise routine and 1.3-2.2x faster
   than arb_sin_cos_arf_atan_reduction.  Precomputation: the 32
   angles come from the 32-term Machin-type set of atan_gauss.c
   (see machin_tab.c):
   20 ms at 1024 limbs, 66 at 2677, 380 at 10000, against 50, 227
   and 1710 ms for the bitwise atan table (13 angles: 10, 31, 145),
   so the bitwise routine amortizes its table after some 100-165
   evaluations.  The alg choice of fixed_sin_cos_reduced at r ~ 100
   for a few thousand limbs picks the sqrt variant where the
   one-burst-step variant is 20% faster; that tuning is a separate
   matter.  */

/* Gaussian quantities carry a limb exponent: (r, rn, i, in) with
   signed lengths as flint_mpn_mul_complex reports, value
   (r + i im) B^e.  Parts are magnitudes; a zero part is one zero
   limb. */

/* the complex products of mpn_extras normalize their outputs, so a
   vanishing component comes back with length 0 (A is real whenever
   every Gaussian exponent is zero, as for a tiny argument); the mpn
   primitives require at least one limb, so restore that */
static void
gz_fixlen(nn_ptr z, slong * zn)
{
    if (*zn == 0)
    {
        z[0] = 0;
        *zn = 1;
    }
}

/* window the two parts of (ar, arn, ai, ain) to the top m limbs at a
   common scale: wr, wi get m limbs each (zero padded), the dropped
   limb count is returned */
static slong
gz_window(nn_ptr wr, nn_ptr wi, nn_srcptr ar, slong arn, nn_srcptr ai,
    slong ain, slong m)
{
    slong L = FLINT_MAX(FLINT_ABS(arn), FLINT_ABS(ain)), drop = L - m, k;
    FLINT_ASSERT(drop >= 0);
    for (k = 0; k < m; k++)
    {
        slong idx = drop + k;
        wr[k] = (idx < FLINT_ABS(arn)) ? ar[idx] : 0;
        wi[k] = (idx < FLINT_ABS(ain)) ? ai[idx] : 0;
    }
    return drop;
}

/* Z = prod_j pi_j^(c_j) over the primes (a, b) at primes + 2j
   (conjugates for negative c_j), a Gaussian integer with signed
   lengths and a limb exponent *ze: parts are kept at most cap limbs
   long, a squaring whose operand exceeds the cap reading the top cap
   limbs of both parts and computing only the high half of the
   square (flint_mpn_sqrhigh_n_complex, within 3 ulps of its lowest
   returned limb).  zr, zi must have room for 2 cap + 4 limbs. */
static void
_fixed_gauss_powers(nn_ptr zr, slong * zrn, nn_ptr zi, slong * zin,
    slong * ze, const signed char * primes, const slong * c, slong num,
    slong cap)
{
    slong emax = 0, j, k, L, rn, in, e = 0;
    nn_ptr ar, ai, br, bi, wr, wi;
    slong arn, ain, brn, bin;
    ulong m[2];
    slong mr, mi;   /* accumulator a + b i, small */

    for (j = 0; j < num; j++)
        emax = FLINT_MAX(emax, FLINT_ABS(c[j]));

    if (emax == 0)
    {
        zr[0] = 1; *zrn = 1;
        zi[0] = 0; *zin = 1;
        *ze = 0;
        return;
    }

    L = FLINT_BIT_COUNT((ulong) emax);

    /* each part gets 2 cap + 4 limbs: an exact squaring of parts at
       most cap long writes at most 2 cap + 1 */
    ar = flint_malloc(4 * (cap + 2) * sizeof(ulong));
    ai = ar + 2 * (cap + 2);
    br = flint_malloc(4 * (cap + 2) * sizeof(ulong));
    bi = br + 2 * (cap + 2);
    wr = flint_malloc(2 * (cap + 1) * sizeof(ulong));
    wi = wr + (cap + 1);

    ar[0] = 1; arn = 1;
    ai[0] = 0; ain = 1;

    for (k = L - 1; k >= 0; k--)
    {
        if (k != L - 1)
        {
            slong len = FLINT_MAX(FLINT_ABS(arn), FLINT_ABS(ain));

            if (len <= cap)
            {
                /* (b) = (a)^2 exactly, then truncated to cap limbs
                   when it grew beyond that */
                slong len2;
                flint_mpn_sqr_complex(br, &brn, bi, &bin,
                    ar, FLINT_ABS(arn), arn < 0, ai, FLINT_ABS(ain), ain < 0);
                gz_fixlen(br, &brn);
                gz_fixlen(bi, &bin);
                e *= 2;
                len2 = FLINT_MAX(FLINT_ABS(brn), FLINT_ABS(bin));
                if (len2 > cap)
                {
                    slong drop = len2 - cap, t;
                    nn_ptr src;
                    for (t = 0; t < 2; t++)
                    {
                        nn_ptr dst = t ? bi : br;
                        slong * ln = t ? &bin : &brn;
                        slong al = FLINT_ABS(*ln), nl;
                        int sg = *ln < 0;
                        src = dst;
                        nl = FLINT_MAX(al - drop, 0);
                        if (nl > 0)
                            flint_mpn_copyi(dst, src + drop, nl);
                        else
                        {
                            dst[0] = 0;
                            nl = 1;
                        }
                        *ln = sg ? -nl : nl;
                    }
                    e += drop;
                }
            }
            else
            {
                /* the high half of the square of the top cap limbs */
                int sr, si;
                slong drop = gz_window(wr, wi, ar, arn, ai, ain, cap);
                flint_mpn_sqrhigh_n_complex(br, &sr, bi, &si,
                    wr, arn < 0, wi, ain < 0, cap);
                /* cap + 1 limbs each: limbs [cap, 2 cap] of the square */
                brn = cap + 1; bin = cap + 1;
                while (brn > 1 && br[brn - 1] == 0) brn--;
                while (bin > 1 && bi[bin - 1] == 0) bin--;
                if (sr) brn = -brn;
                if (si) bin = -bin;
                e = 2 * (e + drop) + cap;
            }
            { nn_ptr t; slong tn;
              t = ar; ar = br; br = t; t = ai; ai = bi; bi = t;
              tn = arn; arn = brn; brn = tn; tn = ain; ain = bin; bin = tn; }
        }

        /* multiply in the primes whose exponent has bit k set, as
           conjugates for negative exponents, accumulated in a small
           Gaussian integer and flushed before it could overflow
           (never for the 13-prime table: the norms' product is below
           2^70, so the parts stay below 2^36) */
        mr = 1; mi = 0;
        for (j = 0; j < num; j++)
        {
            slong ex = FLINT_ABS(c[j]);
            if (ex > 0 && ((ex >> k) & 1))
            {
                slong a = primes[2 * j];
                slong b = primes[2 * j + 1];
                slong nr, ni;
                if (c[j] < 0)
                    b = -b;
                if (FLINT_ABS(mr) > (WORD(1) << (FLINT_BITS - 8))
                    || FLINT_ABS(mi) > (WORD(1) << (FLINT_BITS - 8)))
                {
                    /* flush */
                    m[0] = FLINT_ABS(mr); m[1] = FLINT_ABS(mi);
                    flint_mpn_mul_complex(br, &brn, bi, &bin,
                        ar, FLINT_ABS(arn), arn < 0, ai, FLINT_ABS(ain), ain < 0,
                        m, 1, mr < 0, m + 1, 1, mi < 0);
                    gz_fixlen(br, &brn);
                    gz_fixlen(bi, &bin);
                    { nn_ptr t; slong tn;
                      t = ar; ar = br; br = t; t = ai; ai = bi; bi = t;
                      tn = arn; arn = brn; brn = tn; tn = ain; ain = bin; bin = tn; }
                    mr = 1; mi = 0;
                }
                nr = mr * a - mi * b;
                ni = mr * b + mi * a;
                mr = nr; mi = ni;
            }
        }
        if (mr != 1 || mi != 0)
        {
            m[0] = FLINT_ABS(mr); m[1] = FLINT_ABS(mi);
            flint_mpn_mul_complex(br, &brn, bi, &bin,
                ar, FLINT_ABS(arn), arn < 0, ai, FLINT_ABS(ain), ain < 0,
                m, 1, mr < 0, m + 1, 1, mi < 0);
            gz_fixlen(br, &brn);
            gz_fixlen(bi, &bin);
            { nn_ptr t; slong tn;
              t = ar; ar = br; br = t; t = ai; ai = bi; bi = t;
              tn = arn; arn = brn; brn = tn; tn = ain; ain = bin; bin = tn; }
        }
    }

    rn = FLINT_ABS(arn); in = FLINT_ABS(ain);
    flint_mpn_copyi(zr, ar, rn);
    flint_mpn_copyi(zi, ai, in);
    *zrn = arn; *zin = ain;
    *ze = e;

    flint_free(ar < br ? ar : br);
    flint_free(ar < br ? br : ar);
    flint_free(wr);
}

/* res = floor(|z| B^sh / N) -> (res, wn + 1): the quotient with wn
   fraction limbs once the scales are accounted for (sh limbs up, or
   down: a right shift floors, composing exactly with the division) */
static void
_div_by_N(nn_ptr res, nn_srcptr z, slong zn, slong sh, nn_srcptr N,
    slong Nn, slong wn)
{
    slong qn, nn;
    nn_ptr q, num;
    TMP_INIT;

    while (zn > 0 && z[zn - 1] == 0)
        zn--;
    flint_mpn_zero(res, wn + 1);
    TMP_START;
    if (sh >= 0)
    {
        num = TMP_ALLOC((zn + sh + 1) * sizeof(ulong));
        flint_mpn_zero(num, sh);
        flint_mpn_copyi(num + sh, z, zn);
        nn = zn + sh;
    }
    else
    {
        num = (nn_ptr) z + FLINT_MIN(-sh, zn);
        nn = FLINT_MAX(zn + sh, 0);
    }
    if (nn < Nn)
    {
        TMP_END;
        return;
    }
    qn = nn - Nn + 1;
    q = TMP_ALLOC((qn + 1) * sizeof(ulong));
    flint_mpn_tdiv_q(q, num, nn, N, Nn);
    while (qn > 0 && q[qn - 1] == 0)
        qn--;
    FLINT_ASSERT(qn <= wn + 1);
    flint_mpn_copyi(res, q, qn);
    TMP_END;
}

/* tan = |zi| / |zr| from the two fixed-point parts (same scale, wn
   fraction limbs), into (res, wn + 1) with wn fraction limbs and a
   unit limb.  Both are read in the window of m = wn + 3 limbs below
   the top of the larger one (absolute precision is what the quotient
   needs), and fixed_div_newton divides the two fractions; when the
   numerator is the longer one (tan > 1) the denominator's top window
   limb is zero and it is passed one limb shorter, which makes the
   quotient tan / B and moves the read-off up a limb -- the same
   offset either way since the fraction count grows by one too.
   |error| <= 4 B^-(wn+2) / a < 4 B^-(wn+1), far below one ulp. */
static void
_tan_from_parts(nn_ptr res, nn_srcptr zi, slong zin, nn_srcptr zr,
    slong zrn, slong wn)
{
    slong m = wn + 3, top, k;
    nn_ptr a, b, q;
    TMP_INIT;

    while (zrn > 0 && zr[zrn - 1] == 0) zrn--;
    while (zin > 0 && zi[zin - 1] == 0) zin--;
    FLINT_ASSERT(zrn > 0 && zin <= zrn + 1);

    TMP_START;
    a = TMP_ALLOC((3 * m + 6) * sizeof(ulong));
    b = a + m;
    q = b + m;

    top = FLINT_MAX(zrn, zin);
    for (k = 0; k < m; k++)
    {
        slong idx = top - m + k;
        a[k] = (idx >= 0 && idx < zrn) ? zr[idx] : 0;
        b[k] = (idx >= 0 && idx < zin) ? zi[idx] : 0;
    }

    if (a[m - 1] != 0)
        fixed_div_newton(q, b, m, a, m, wn + 2);       /* q = tan */
    else
        fixed_div_newton(q, b, m, a, m - 1, wn + 3);   /* q = tan / B */

    /* tan B^wn = sum q[k] B^(k-2) in both cases */
    flint_mpn_copyi(res, q + 2, wn + 1);
    TMP_END;
}

/* C = floor(|zr| B^sh / N), S = floor(|zi| B^sh / N) for both parts
   at once (sh as in _div_by_N): one reciprocal of N (fixed_inv_newton
   on N's top m = wn + 3 limbs) and one product per part, against two
   divisions.  Each part is read in its top m limbs; with q the
   reciprocal to wn + 3 fraction limbs, part B^sh / N = (part_window
   q_int) B^(e) for the limb offset e worked out below, and the
   quotient is the product's window at that offset.  Errors: the
   reciprocal within 4 B^-(wn+3) / a <= 4 B^-(wn+2) relatively, the
   windows within B^-(wn+2) relatively, the final floor one working
   ulp -- all inside the guard limb. */
static void
_div2_by_N(nn_ptr C, nn_ptr S, nn_srcptr zr, slong zrn, nn_srcptr zi,
    slong zin, slong sh, nn_srcptr N, slong Nn, slong wn)
{
    slong m = wn + 3, dN, k, t;
    nn_ptr Nw, q, zw, P;
    TMP_INIT;

    while (Nn > 0 && N[Nn - 1] == 0) Nn--;
    FLINT_ASSERT(Nn > 0);

    TMP_START;
    Nw = TMP_ALLOC(m * sizeof(ulong));
    q = TMP_ALLOC((m + 3) * sizeof(ulong));
    zw = TMP_ALLOC(m * sizeof(ulong));
    P = TMP_ALLOC((2 * m + 4) * sizeof(ulong));

    /* N = Nw_int B^dN, Nw_int / B^m in [1/B, 1) */
    dN = Nn - m;
    for (k = 0; k < m; k++)
    {
        slong idx = dN + k;
        Nw[k] = (idx >= 0 && idx < Nn) ? N[idx] : 0;
    }
    fixed_inv_newton(q, Nw, m, m);      /* q_int = B^(2m) / Nw_int, m + 2 limbs */

    for (t = 0; t < 2; t++)
    {
        nn_srcptr z = t ? zi : zr;
        slong zn = t ? zin : zrn, dZ, e, pn;
        nn_ptr res = t ? S : C;

        while (zn > 0 && z[zn - 1] == 0) zn--;
        flint_mpn_zero(res, wn + 1);
        if (zn == 0)
            continue;

        /* z = zw_int B^dZ */
        dZ = zn - m;
        for (k = 0; k < m; k++)
        {
            slong idx = dZ + k;
            zw[k] = (idx >= 0 && idx < zn) ? z[idx] : 0;
        }
        /* z B^sh / N = zw_int B^(dZ + sh) / (Nw_int B^dN)
                      = zw_int q_int B^(dZ + sh - dN - 2m) */
        flint_mpn_mul(P, q, m + 2, zw, m);   /* the longer operand first */
        pn = 2 * m + 2;
        e = dZ + sh - dN - 2 * m;
        /* the quotient is floor(P B^e): drop -e limbs (e <= 0 here:
           the quotient is at most one unit) */
        FLINT_ASSERT(e <= 0);
        if (-e < pn)
        {
            slong cnt = FLINT_MIN(pn + e, wn + 1);
            flint_mpn_copyi(res, P - e, cnt);
        }
    }
    TMP_END;
}

/* the worker behind the three entry points: any of ysin, ycos, ytan
   may be NULL; ytan is tan x with n fraction limbs and a unit limb */
static void
_fixed_trig_diophantine(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x,
    slong n, slong num, double max_weight)
{
    const fixed_rel_struct * tab;
    slong * rel;
    slong xn = n, wr, wn, G, j, i, r, cap, Nn;
    slong arn, ain, Xn, Yn, Zrn, Zin, eA, eN, eX, eZ;
    ulong relsum;
    double eps_min;
    nn_ptr base, t, tmp, s, g, cr, Ar, Ai, X, Y, N, Zr, Zi, C, S;
    int neg, rot;
    TMP_INIT;

    FLINT_ASSERT(num >= 2 && num <= FIXED_ATAN_GAUSS_MAX);

    /* keep the coefficients within the slong range whatever budget
       is requested (see exp_diophantine.c) */
    max_weight = FLINT_MIN(max_weight, ldexp(1.0, FLINT_BITS - 11));

    while (xn > 0 && x[xn - 1] == 0)
        xn--;
    if (xn == 0)
    {
        if (ysin != NULL) flint_mpn_zero(ysin, n + 1);
        if (ycos != NULL) { flint_mpn_zero(ycos, n); ycos[n] = 1; }
        if (ytan != NULL) flint_mpn_zero(ytan, n + 1);
        return;
    }

    tab = fixed_rel_table(1, num);

    if (FLINT_BITS * n <= 10000)
        wr = 256 / FLINT_BITS;
    else if (FLINT_BITS * n <= 100000)
        wr = 512 / FLINT_BITS;
    else
        wr = 768 / FLINT_BITS;
    wr = FLINT_MAX(wr, (slong) ((-log2(tab->epsilon_min) + 80) / FLINT_BITS) + 1);

    _fixed_atan_gauss_ensure(num, FLINT_MAX(wr, n + 1));

    TMP_START;
    rel = TMP_ALLOC(num * sizeof(slong));

    base = TMP_ALLOC(wr * sizeof(ulong));
    {
        slong pad = FLINT_MAX(wr - n, 0);
        for (i = 0; i < pad; i++)
            base[i] = 0;
        flint_mpn_copyi(base + pad, x + FLINT_MAX(n - wr, 0), wr - pad);
    }

    eps_min = ldexp(1.0, -(int) FLINT_MIN(FLINT_BITS * n + 32, 2000));
    eps_min = FLINT_MAX(eps_min, ldexp(1.0, -(int) (FLINT_BITS * wr - 48)));

    _fixed_log_reduce(rel, tab, base, wr, max_weight, eps_min,
        _fixed_atan_gauss_entry(0, wr), _fixed_atan_gauss_n);

    relsum = 0;
    for (j = 0; j < num; j++)
        relsum += (ulong) FLINT_ABS(rel[j]);
    G = (FLINT_BIT_COUNT(4 * relsum + 256) + 2 + 8 + FLINT_BITS - 1)
        / FLINT_BITS;
    wn = n + G;

    _fixed_atan_gauss_ensure(num, wn);

    /* t = x B^G - sum c_j theta_j, two's complement; |t| < 1 and its
       sign only flips the sine of the residual */
    t = TMP_ALLOC(2 * (wn + 1) * sizeof(ulong));
    tmp = t + (wn + 1);
    flint_mpn_zero(tmp, G);
    flint_mpn_copyi(tmp + G, x, n);
    tmp[wn] = 0;
    _fixed_log_dot(t, tmp, wn + 1, rel, num,
        _fixed_atan_gauss_entry(0, wn), _fixed_atan_gauss_n);
    neg = (t[wn] >> (FLINT_BITS - 1)) != 0;
    if (neg)
        mpn_neg(t, t, wn + 1);
    FLINT_ASSERT(t[wn] == 0);

    /* (s, g) = (sin |t|, 1 - cos |t|) */
    s = TMP_ALLOC(3 * (wn + 1) * sizeof(ulong));
    g = s + (wn + 1);
    cr = g + (wn + 1);
    for (i = wn - 1; i >= 0 && t[i] == 0; i--)
        ;
    if (i < 0)
    {
        flint_mpn_zero(s, wn);
        flint_mpn_zero(g, wn);
    }
    else
    {
        r = FLINT_BITS * (wn - 1 - i) + (FLINT_BITS - FLINT_BIT_COUNT(t[i]));
        if (r >= 16)
            fixed_sin_cos_reduced(s, g, t, wn, (flint_bitcnt_t) r, 0);
        else
        {
            /* a large residual (tiny weight budgets): the general
               routine, 1 - cos from the cosine */
            nn_ptr c2 = TMP_ALLOC((wn + 1) * sizeof(ulong));
            fixed_sin_cos_notab(s, c2, t, wn);
            /* g = 1 - cos in wn fraction limbs (cos in [0, 1]) */
            flint_mpn_zero(g, wn);
            if (c2[wn] == 0)
            {
                mpn_neg(g, c2, wn);   /* B^wn - cos */
            }
        }
    }
    /* cr = 1 - g = cos |t|, wn fraction limbs and a unit limb */
    {
        slong k;
        for (k = 0; k < wn && g[k] == 0; k++)
            ;
        if (k == wn)
        {
            flint_mpn_zero(cr, wn);
            cr[wn] = 1;
        }
        else
        {
            mpn_neg(cr, g, wn);
            cr[wn] = 0;
        }
    }

    /* A = prod pi_j^(c_j) (over j >= 1: the (1+i) factor of index 0
       is the quarter-turn i^c_0 applied at the end), parts capped at
       cap = wn + 3 limbs with a limb exponent eA */
    cap = wn + 3;
    Ar = TMP_ALLOC(2 * (2 * cap + 4) * sizeof(ulong));
    Ai = Ar + (2 * cap + 4);
    _fixed_gauss_powers(Ar, &arn, Ai, &ain, &eA, _fixed_gaussian_primes + 2,
        rel + 1, num - 1, cap);

    /* N = A_r^2 + A_i^2 (the norm, a rational integer -- up to the
       truncations) and A^2 = X + iY, exact when A fits the cap and
       otherwise the high halves of the squares of A's top cap limbs */
    N = TMP_ALLOC((2 * cap + 4) * sizeof(ulong));
    X = TMP_ALLOC((2 * cap + 4) * sizeof(ulong));
    Y = TMP_ALLOC((2 * cap + 4) * sizeof(ulong));
    {
        slong rn = FLINT_ABS(arn), in = FLINT_ABS(ain);
        slong len = FLINT_MAX(rn, in);
        nn_ptr sq = TMP_ALLOC((2 * cap + 4) * sizeof(ulong));

        if (len <= cap)
        {
            flint_mpn_sqr(N, Ar, rn);
            Nn = 2 * rn;
            flint_mpn_sqr(sq, Ai, in);
            if (2 * in > Nn)
            {
                flint_mpn_zero(N + Nn, 2 * in - Nn);
                Nn = 2 * in;
            }
            N[Nn] = mpn_add(N, N, Nn, sq, 2 * in);
            Nn += (N[Nn] != 0);
            while (Nn > 0 && N[Nn - 1] == 0)
                Nn--;
            eN = 2 * eA;
            flint_mpn_sqr_complex(X, &Xn, Y, &Yn, Ar, rn, arn < 0, Ai, in, ain < 0);
            gz_fixlen(X, &Xn);
            gz_fixlen(Y, &Yn);
            eX = 2 * eA;
        }
        else
        {
            nn_ptr wr = TMP_ALLOC(2 * (cap + 1) * sizeof(ulong));
            nn_ptr wi = wr + (cap + 1);
            slong drop = gz_window(wr, wi, Ar, arn, Ai, ain, cap);
            int sx, sy;
            /* high halves: limbs [cap, 2 cap) of each square, the
               returned limb below them dropped (within 3 ulps) */
            flint_mpn_sqrhigh(N, wr, cap);
            flint_mpn_sqrhigh(sq, wi, cap);
            N[cap] = mpn_add_n(N, N, sq, cap);
            Nn = cap + 1;
            while (Nn > 0 && N[Nn - 1] == 0)
                Nn--;
            eN = 2 * (eA + drop) + cap;
            flint_mpn_sqrhigh_n_complex(X, &sx, Y, &sy, wr, arn < 0, wi, ain < 0, cap);
            Xn = cap + 1; Yn = cap + 1;
            while (Xn > 1 && X[Xn - 1] == 0) Xn--;
            while (Yn > 1 && Y[Yn - 1] == 0) Yn--;
            if (sx) Xn = -Xn;
            if (sy) Yn = -Yn;
            eX = 2 * (eA + drop) + cap;
        }
    }

    /* Z = A^2 (cos t + i sin t); sin t has the residual's sign.  Exact
       (wn fraction limbs: eZ = eX - wn) when A^2 fits the cap, else
       the high half of the product of the two pairs windowed to a
       common length m (the shorter pair padded below with zeros,
       which its exponent absorbs) */
    Zr = TMP_ALLOC((2 * cap + wn + 6) * sizeof(ulong));
    Zi = TMP_ALLOC((2 * cap + wn + 6) * sizeof(ulong));
    {
        slong lx = FLINT_MAX(FLINT_ABS(Xn), FLINT_ABS(Yn));

        if (lx <= cap)
        {
            flint_mpn_mul_complex(Zr, &Zrn, Zi, &Zin,
                X, FLINT_ABS(Xn), Xn < 0, Y, FLINT_ABS(Yn), Yn < 0,
                cr, wn + 1, 0, s, wn, neg);
            gz_fixlen(Zr, &Zrn);
            gz_fixlen(Zi, &Zin);
            eZ = eX - wn;
        }
        else
        {
            slong m2 = FLINT_MAX(lx, wn + 1), dropx, k;
            nn_ptr wr = TMP_ALLOC(4 * (m2 + 1) * sizeof(ulong));
            nn_ptr wi = wr + (m2 + 1), vr = wi + (m2 + 1), vi = vr + (m2 + 1);
            int sr, si;
            dropx = gz_window(wr, wi, X, Xn, Y, Yn, m2);
            /* cr (wn + 1 limbs), s (wn) padded below to m2 limbs */
            for (k = 0; k < m2; k++)
            {
                slong idx = k - (m2 - (wn + 1));
                vr[k] = (idx >= 0 && idx < wn + 1) ? cr[idx] : 0;
                vi[k] = (idx >= 0 && idx < wn) ? s[idx] : 0;
            }
            flint_mpn_mulhigh_n_complex(Zr, &sr, Zi, &si,
                wr, Xn < 0, wi, Yn < 0, vr, 0, vi, neg, m2);
            Zrn = m2 + 1; Zin = m2 + 1;
            while (Zrn > 1 && Zr[Zrn - 1] == 0) Zrn--;
            while (Zin > 1 && Zi[Zin - 1] == 0) Zin--;
            if (sr) Zrn = -Zrn;
            if (si) Zin = -Zin;
            /* X B^(eX + dropx) times e^(it) B^(-wn - (m2 - wn - 1)),
               the product's high half at B^m2 */
            eZ = (eX + dropx) + (-wn - (m2 - (wn + 1))) + m2;
        }
    }

    /* rotate by i^c_0: e^(ix) = i^c_0 (Z_r + i Z_i) / N; which part
       and sign feed which output */
    rot = (int) (rel[0] & 3);
    {
        int sr = Zrn < 0, si = Zin < 0;
        nn_srcptr cosp, sinp;
        slong cosn, sinn;
        int cos_neg, sin_neg;
        switch (rot)
        {
            case 0: cosp = Zr; cosn = FLINT_ABS(Zrn); cos_neg = sr;
                    sinp = Zi; sinn = FLINT_ABS(Zin); sin_neg = si; break;
            case 1: cosp = Zi; cosn = FLINT_ABS(Zin); cos_neg = !si;
                    sinp = Zr; sinn = FLINT_ABS(Zrn); sin_neg = sr; break;
            case 2: cosp = Zr; cosn = FLINT_ABS(Zrn); cos_neg = !sr;
                    sinp = Zi; sinn = FLINT_ABS(Zin); sin_neg = !si; break;
            default: cosp = Zi; cosn = FLINT_ABS(Zin); cos_neg = si;
                     sinp = Zr; sinn = FLINT_ABS(Zrn); sin_neg = !sr; break;
        }
        /* x in [0, 1): both nonnegative; a negative sign can only
           accompany a zero magnitude */
        (void) cos_neg; (void) sin_neg;

        C = TMP_ALLOC(2 * (wn + 2) * sizeof(ulong));
        S = C + (wn + 2);

        /* the quotients with wn fraction limbs: Z B^eZ / (N B^eN) B^wn;
           both through one reciprocal of N when both are wanted */
        if (ycos != NULL && ysin != NULL)
        {
            _div2_by_N(C, S, cosp, cosn, sinp, sinn, eZ - eN + wn, N, Nn, wn);
            flint_mpn_copyi(ycos, C + G, n + 1);
            flint_mpn_copyi(ysin, S + G, n + 1);
        }
        else if (ycos != NULL)
        {
            _div_by_N(C, cosp, cosn, eZ - eN + wn, N, Nn, wn);
            flint_mpn_copyi(ycos, C + G, n + 1);
        }
        else if (ysin != NULL)
        {
            _div_by_N(S, sinp, sinn, eZ - eN + wn, N, Nn, wn);
            flint_mpn_copyi(ysin, S + G, n + 1);
        }
        if (ytan != NULL)
        {
            /* tan x = sin / cos = the ratio of the parts: no N */
            nn_ptr T = TMP_ALLOC((wn + 2) * sizeof(ulong));
            _tan_from_parts(T, sinp, sinn, cosp, cosn, wn);
            flint_mpn_copyi(ytan, T + G, n + 1);
        }
    }

    TMP_END;
}

void
_fixed_sin_cos_diophantine_tune(nn_ptr ysin, nn_ptr ycos, nn_srcptr x,
    slong n, slong num, double max_weight)
{
    _fixed_trig_diophantine(ysin, ycos, NULL, x, n, num, max_weight);
}

/* 32 primes at four times the precision as weight: the reduced sine
   and cosine only reach the exponential's efficiency from r ~ 300,
   so a deeper reduction pays more here than in exp_diophantine, and
   with the products capped a large budget costs nothing (see the
   tuning notes above) */
#define FIXED_TRIG_DIOPHANTINE_NUM 32
#define FIXED_TRIG_DIOPHANTINE_WEIGHT 4.0

void
fixed_sin_cos_diophantine(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)
{
    _fixed_trig_diophantine(ysin, ycos, NULL, x, n, FIXED_TRIG_DIOPHANTINE_NUM,
        FIXED_TRIG_DIOPHANTINE_WEIGHT * (double) (FLINT_BITS * n));
}

void
fixed_tan_diophantine(nn_ptr res, nn_srcptr x, slong n)
{
    _fixed_trig_diophantine(NULL, NULL, res, x, n, FIXED_TRIG_DIOPHANTINE_NUM,
        FIXED_TRIG_DIOPHANTINE_WEIGHT * (double) (FLINT_BITS * n));
}

void
_fixed_tan_diophantine_tune(nn_ptr res, nn_srcptr x, slong n, slong num,
    double max_weight)
{
    _fixed_trig_diophantine(NULL, NULL, res, x, n, num, max_weight);
}
