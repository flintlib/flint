/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "padic_radix.h"
#include "gr.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "nmod_poly.h"

/*
    Port of padic_exp_rectangular (src/padic/exp_rectangular.c) to the
    padic_radix module.

    The series

        exp(x) = sum_{i=0}^{n-1} x^i / i!,   x = p^v u,

    is evaluated by rectangular splitting.  All of the heavy arithmetic is
    done in fixed point with radix_integer_t, working modulo (p^e)^l = B^l,
    where l is the minimal limb precision required to reach the working digit
    precision.  Reciprocals of the factorials are deferred: the splitting
    builds an integer numerator S together with a denominator F (a running
    product of the integers being divided out), both modulo B^l, and a single
    Hensel inversion of F recovers the unit at the end.

    The working digit precision is N + k, where k bounds v_p((n-1)!) (the
    p-adic valuation of the denominator that is divided out).  Carrying these k
    guard digits is purely internal: the answer is returned modulo p^N, and the
    caller's input only needs to be known to N digits.
*/

/* exp_factab[m][i] = m!/i! for 0 <= i <= m for factorials that fit a
   limb: 20! < 2^64 and 12! < 2^32. */

#if FLINT_BITS == 64
#define EXP_FACTAB_MAXM 20
static const ulong exp_factab[EXP_FACTAB_MAXM + 1][EXP_FACTAB_MAXM + 1] = {
    { UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(1), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(2), UWORD(2), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(6), UWORD(6), UWORD(3), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(24), UWORD(24), UWORD(12), UWORD(4), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(120), UWORD(120), UWORD(60), UWORD(20), UWORD(5), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(720), UWORD(720), UWORD(360), UWORD(120), UWORD(30), UWORD(6), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(5040), UWORD(5040), UWORD(2520), UWORD(840), UWORD(210), UWORD(42), UWORD(7), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(40320), UWORD(40320), UWORD(20160), UWORD(6720), UWORD(1680), UWORD(336), UWORD(56), UWORD(8), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(362880), UWORD(362880), UWORD(181440), UWORD(60480), UWORD(15120), UWORD(3024), UWORD(504), UWORD(72), UWORD(9), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(3628800), UWORD(3628800), UWORD(1814400), UWORD(604800), UWORD(151200), UWORD(30240), UWORD(5040), UWORD(720), UWORD(90), UWORD(10), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(39916800), UWORD(39916800), UWORD(19958400), UWORD(6652800), UWORD(1663200), UWORD(332640), UWORD(55440), UWORD(7920), UWORD(990), UWORD(110), UWORD(11), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(479001600), UWORD(479001600), UWORD(239500800), UWORD(79833600), UWORD(19958400), UWORD(3991680), UWORD(665280), UWORD(95040), UWORD(11880), UWORD(1320), UWORD(132), UWORD(12), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(6227020800), UWORD(6227020800), UWORD(3113510400), UWORD(1037836800), UWORD(259459200), UWORD(51891840), UWORD(8648640), UWORD(1235520), UWORD(154440), UWORD(17160), UWORD(1716), UWORD(156), UWORD(13), UWORD(1), 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(87178291200), UWORD(87178291200), UWORD(43589145600), UWORD(14529715200), UWORD(3632428800), UWORD(726485760), UWORD(121080960), UWORD(17297280), UWORD(2162160), UWORD(240240), UWORD(24024), UWORD(2184), UWORD(182), UWORD(14), UWORD(1), 0, 0, 0, 0, 0, 0 },
    { UWORD(1307674368000), UWORD(1307674368000), UWORD(653837184000), UWORD(217945728000), UWORD(54486432000), UWORD(10897286400), UWORD(1816214400), UWORD(259459200), UWORD(32432400), UWORD(3603600), UWORD(360360), UWORD(32760), UWORD(2730), UWORD(210), UWORD(15), UWORD(1), 0, 0, 0, 0, 0 },
    { UWORD(20922789888000), UWORD(20922789888000), UWORD(10461394944000), UWORD(3487131648000), UWORD(871782912000), UWORD(174356582400), UWORD(29059430400), UWORD(4151347200), UWORD(518918400), UWORD(57657600), UWORD(5765760), UWORD(524160), UWORD(43680), UWORD(3360), UWORD(240), UWORD(16), UWORD(1), 0, 0, 0, 0 },
    { UWORD(355687428096000), UWORD(355687428096000), UWORD(177843714048000), UWORD(59281238016000), UWORD(14820309504000), UWORD(2964061900800), UWORD(494010316800), UWORD(70572902400), UWORD(8821612800), UWORD(980179200), UWORD(98017920), UWORD(8910720), UWORD(742560), UWORD(57120), UWORD(4080), UWORD(272), UWORD(17), UWORD(1), 0, 0, 0 },
    { UWORD(6402373705728000), UWORD(6402373705728000), UWORD(3201186852864000), UWORD(1067062284288000), UWORD(266765571072000), UWORD(53353114214400), UWORD(8892185702400), UWORD(1270312243200), UWORD(158789030400), UWORD(17643225600), UWORD(1764322560), UWORD(160392960), UWORD(13366080), UWORD(1028160), UWORD(73440), UWORD(4896), UWORD(306), UWORD(18), UWORD(1), 0, 0 },
    { UWORD(121645100408832000), UWORD(121645100408832000), UWORD(60822550204416000), UWORD(20274183401472000), UWORD(5068545850368000), UWORD(1013709170073600), UWORD(168951528345600), UWORD(24135932620800), UWORD(3016991577600), UWORD(335221286400), UWORD(33522128640), UWORD(3047466240), UWORD(253955520), UWORD(19535040), UWORD(1395360), UWORD(93024), UWORD(5814), UWORD(342), UWORD(19), UWORD(1), 0 },
    { UWORD(2432902008176640000), UWORD(2432902008176640000), UWORD(1216451004088320000), UWORD(405483668029440000), UWORD(101370917007360000), UWORD(20274183401472000), UWORD(3379030566912000), UWORD(482718652416000), UWORD(60339831552000), UWORD(6704425728000), UWORD(670442572800), UWORD(60949324800), UWORD(5079110400), UWORD(390700800), UWORD(27907200), UWORD(1860480), UWORD(116280), UWORD(6840), UWORD(380), UWORD(20), UWORD(1) },
};
#else
#define EXP_FACTAB_MAXM 12
static const ulong exp_factab[EXP_FACTAB_MAXM + 1][EXP_FACTAB_MAXM + 1] = {
    { UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(1), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(2), UWORD(2), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(6), UWORD(6), UWORD(3), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(24), UWORD(24), UWORD(12), UWORD(4), UWORD(1), 0, 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(120), UWORD(120), UWORD(60), UWORD(20), UWORD(5), UWORD(1), 0, 0, 0, 0, 0, 0, 0 },
    { UWORD(720), UWORD(720), UWORD(360), UWORD(120), UWORD(30), UWORD(6), UWORD(1), 0, 0, 0, 0, 0, 0 },
    { UWORD(5040), UWORD(5040), UWORD(2520), UWORD(840), UWORD(210), UWORD(42), UWORD(7), UWORD(1), 0, 0, 0, 0, 0 },
    { UWORD(40320), UWORD(40320), UWORD(20160), UWORD(6720), UWORD(1680), UWORD(336), UWORD(56), UWORD(8), UWORD(1), 0, 0, 0, 0 },
    { UWORD(362880), UWORD(362880), UWORD(181440), UWORD(60480), UWORD(15120), UWORD(3024), UWORD(504), UWORD(72), UWORD(9), UWORD(1), 0, 0, 0 },
    { UWORD(3628800), UWORD(3628800), UWORD(1814400), UWORD(604800), UWORD(151200), UWORD(30240), UWORD(5040), UWORD(720), UWORD(90), UWORD(10), UWORD(1), 0, 0 },
    { UWORD(39916800), UWORD(39916800), UWORD(19958400), UWORD(6652800), UWORD(1663200), UWORD(332640), UWORD(55440), UWORD(7920), UWORD(990), UWORD(110), UWORD(11), UWORD(1), 0 },
    { UWORD(479001600), UWORD(479001600), UWORD(239500800), UWORD(79833600), UWORD(19958400), UWORD(3991680), UWORD(665280), UWORD(95040), UWORD(11880), UWORD(1320), UWORD(132), UWORD(12), UWORD(1) },
};
#endif

/* p-adic valuation in digits of a nonnegative radix array; 0 when a == 0. */
static slong
_radix_array_valuation_digits(nn_srcptr a, slong n, const radix_t radix)
{
    slong j = 0;

    while (j < n && a[j] == 0)
        j++;

    if (j == n)
        return 0;

    return j * radix->exp + (slong) _radix_valuation_digits_1(a[j], radix);
}

/*
    rop = unit of exp(p^v u) modulo p^N, with u treated as an exact integer.

    The result is a 1-unit (valuation 0); rop is its nonnegative residue modulo
    p^N.  Requires v >= 1 (odd p) or v >= 2 (p = 2) and N >= 1.  Returns
    GR_SUCCESS, or GR_UNABLE only on the (degenerate) failure of the final
    inversion.
*/
int
_padic_radix_exp_rectangular(radix_integer_t rop, const radix_integer_t u,
    slong v, slong N, const radix_t radix)
{
    ulong p = DIGIT_RADIX(radix);
    ulong B = LIMB_RADIX(radix);
    nmod_t Bmod = radix->B;
    slong e = radix->exp;
    slong n, k, M, l, lN, npows, npows_fit, nsums, i, j;
    slong fn, w;
    ulong nm1, c;
    nn_ptr scratch, x, s, summ, tmp, f, acc;
    nn_ptr * pows;
    radix_integer_t xi, summ_ri, f_ri;

    if (N <= 0)
    {
        radix_integer_zero(rop, radix);
        return GR_SUCCESS;
    }

    if (v >= N || radix_integer_is_zero(u, radix))
    {
        radix_integer_one(rop, radix);
        return GR_SUCCESS;
    }

    n = _padic_radix_exp_bound(v, N, p);

    /* k >= v_p((n-1)!): guard digits divided out at the end.  Using the bound
       (n-2)/(p-1) >= ((n-1) - s_p(n-1))/(p-1) = v_p((n-1)!). */
    k = (n >= 2) ? (n - 2) / (slong) (p - 1) : 0;

    M = N + k;                       /* working digit precision */
    l = (M + e - 1) / e;             /* minimal limb precision: l = ceil(M/e) */
    lN = (N + e - 1) / e;            /* limbs covering the target precision N */

    /* Fast path for the smallest precisions: a single limb of working precision
       and a denominator (n-1)! that fits the limb radix. Then the whole series
       numerator is one Horner evaluation of the table row against x mod B, the
       row's first entry is (n-1)!, and one inversion + multiplication finishes. */
    if (l == 1 && n - 1 <= EXP_FACTAB_MAXM)
    {
        const ulong * row = exp_factab[n - 1];     /* coefficients (n-1)!/i! */
        ulong den = row[0];                        /* (n-1)! */

        if (den < B)
        {
            ulong xl, num, den_u, num_u, dinv, res;
            ulong mag0, ulo;

            /* x = p^v u  (mod B), nonnegative single limb */
            mag0 = (u->size == 0) ? 0 : u->d[0];           /* |u| mod B */
            ulo = (u->size >= 0) ? mag0 : (mag0 == 0 ? 0 : B - mag0);
            xl = (v >= e) ? 0 : nmod_mul(ulo, radix->bpow[v], Bmod);

            /* numerator = sum_i (n-1)!/i! * x^i   (mod B) */
            num = _nmod_poly_evaluate_nmod(row, n, xl, Bmod);

            /* Strip the common p-power (= v_p((n-1)!) = v_p(num)) from both. */
            if (p == 2)
            {
                w = flint_ctz(num);
                num_u = num >> w;
                den_u = den >> w;
            }
            else
            {
                w = (slong) _radix_valuation_digits_1(num, radix);
                num_u = num * radix->bpow_oddinv[w].a;
                den_u = den * radix->bpow_oddinv[w].a;
            }

            if (!radix_invmod_bn(&dinv, &den_u, 1, 1, radix))
                flint_abort();

            res = nmod_mul(num_u, dinv, Bmod);   /* exp unit (mod B) */
            res = res % radix->bpow[N];          /* reduce to N digits (N <= e) */

            radix_integer_set_ui(rop, res, radix);
            return GR_SUCCESS;
        }
    }

    /* Choose npows so that a whole block factorial (a product of up to npows of
       the series indices, each < n) fits in a single radix limb: then the
       coefficient recurrence is plain ulong arithmetic and every multiplication
       by c is a radix_mul_1.  npows_fit is the largest m with (n-1)^m < B.  If
       even one index can reach B (n - 1 >= B) no such m >= 1 exists; we decline
       and let the balanced algorithm handle it. */
    nm1 = (ulong) (n - 1);
    if (nm1 <= 1)
    {
        npows_fit = n;               /* factors are 0/1: never overflow */
    }
    else
    {
        ulong prod = 1, hh, ll;
        npows_fit = 0;
        for (;;)
        {
            umul_ppmm(hh, ll, prod, nm1);
            if (hh != 0 || ll >= B)
                break;
            prod = ll;
            npows_fit++;
        }
    }

    if (npows_fit < 1)
        return GR_UNABLE;

    npows = (slong) n_sqrt((ulong) n);
    if (npows < 1)
        npows = 1;
    if (npows > npows_fit)
        npows = npows_fit;
    nsums = (n + npows - 1) / npows;

    /* One allocation for every multi-limb temporary: x, s, summ, tmp, f, the
       power table pows[0..npows] (each l limbs), and the 3-words-per-limb
       delayed-carry accumulator acc (3*l limbs). */
    scratch = flint_malloc(sizeof(ulong) * (npows + 9) * (size_t) l);
    x    = scratch + (slong) 0 * l;
    s    = scratch + (slong) 1 * l;
    summ = scratch + (slong) 2 * l;
    tmp  = scratch + (slong) 3 * l;
    f    = scratch + (slong) 4 * l;
    pows = flint_malloc(sizeof(nn_ptr) * (npows + 1));
    for (i = 0; i <= npows; i++)
        pows[i] = scratch + (slong) (5 + i) * l;
    acc  = scratch + (slong) (npows + 6) * l;

    /* x = p^v u  (mod B^l), as a nonnegative residue, copied into the array. */
    radix_integer_init(xi, radix);
    radix_integer_mod_digits(xi, u, e * l, radix);
    radix_integer_lshift_digits(xi, xi, v, radix);
    radix_integer_mod_limbs(xi, xi, l, radix);
    for (i = 0; i < l; i++)
        x[i] = (i < xi->size) ? xi->d[i] : 0;
    radix_integer_clear(xi, radix);

    /* pows[i] = x^i  (mod B^l),  i = 0, ..., npows */
    for (i = 0; i < l; i++)
        pows[0][i] = 0;
    pows[0][0] = 1;
    flint_mpn_copyi(pows[1], x, l);
    for (i = 2; i <= npows; i++)
        radix_mulmid(pows[i], pows[i - 1], l, pows[1], l, 0, l, radix);

    /* Horner over blocks of npows terms; summ is the numerator, f the running
       denominator, both modulo B^l.  f is kept to its used length fn. */
    for (i = 0; i < l; i++)
        summ[i] = 0;
    for (i = 0; i < l; i++)
        f[i] = 0;
    f[0] = 1;
    fn = 1;

    for (i = nsums - 1; i >= 0; i--)
    {
        slong lo = i * npows;
        slong hi = FLINT_MIN(n - 1, lo + npows - 1);
        ulong carry;
        ulong cy[3];

        /* s = sum over the block of pows[hi-lo] * c, accumulated with delayed
           carries: limbs 3*t, 3*t+1, 3*t+2 of acc hold the unreduced
           contribution to limb t of s.  Each term adds pows[hi-lo][t]*c (a two
           word product) into the three accumulator words for slot t. */
        for (j = 0; j < 3 * l; j++)
        {
            acc[j] = 0;
        }
        c = 1;

        for ( ; hi >= lo; hi--)
        {
            nn_srcptr pj = pows[hi - lo];

            for (j = 0; j < l; j++)
            {
                ulong phi, plo;
                umul_ppmm(phi, plo, pj[j], c);
                add_sssaaaaaa(acc[3 * j + 2], acc[3 * j + 1], acc[3 * j],
                              acc[3 * j + 2], acc[3 * j + 1], acc[3 * j],
                              (ulong) 0, phi, plo);
            }

            if (hi != 0)
                c *= (ulong) hi;       /* stays < B by the npows constraint */
        }

        /* Reconstruct s by carry-propagating acc and reducing mod B, keeping the
           low l limbs (mod B^l); the carry past limb l is discarded. */
        cy[0] = cy[1] = cy[2] = 0;
        if (Bmod.norm == 0)
        {
            for (j = 0; j < l; j++)
            {
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                              acc[3 * j + 2], acc[3 * j + 1], acc[3 * j]);
                s[j] = flint_mpn_divrem_3_1_preinv_norm(cy, cy, Bmod.n, Bmod.ninv);
            }
        }
        else
        {
            for (j = 0; j < l; j++)
            {
                add_sssaaaaaa(cy[2], cy[1], cy[0], cy[2], cy[1], cy[0],
                              acc[3 * j + 2], acc[3 * j + 1], acc[3 * j]);
                s[j] = flint_mpn_divrem_3_1_preinv_unnorm(cy, cy, Bmod.n, Bmod.ninv, Bmod.norm);
            }
        }

        /* summ = s * f + pows[npows] * summ   (mod B^l) */
        radix_mulmid(tmp, pows[npows], l, summ, l, 0, l, radix);
        radix_mulmid(summ, s, l, f, fn, 0, l, radix);
        radix_add(summ, summ, l, tmp, l, radix);

        /* f *= c   (in place; grow the used length by the carry limb) */
        carry = radix_mul_1(f, f, fn, c, radix);
        if (fn < l && carry != 0)
        {
            f[fn] = carry;
            fn++;
        }
    }

    flint_free(pows);

    /* Divide out the factorial. Since exp is a unit, v_p(summ) = v_p(f); strip
       the common power of p from each, then form summ / f modulo p^N. */
    {
        nn_ptr sp = summ, fp = f;
        slong slen = l, flen = l;

        w = _radix_array_valuation_digits(summ, l, radix);
        if (w > 0)
        {
            slong wf = _radix_array_valuation_digits(f, l, radix);
            slong wl, wr;

            wl = w / e;  wr = w % e;
            sp = summ + wl;  slen = l - wl;
            if (wr > 0)
                radix_rshift_digits(sp, sp, slen, (unsigned int) wr, radix);

            wl = wf / e; wr = wf % e;
            fp = f + wl;  flen = l - wl;
            if (wr > 0)
                radix_rshift_digits(fp, fp, flen, (unsigned int) wr, radix);
        }

        while (slen > 0 && sp[slen - 1] == 0)
            slen--;
        while (flen > 0 && fp[flen - 1] == 0)
            flen--;

        summ_ri->d = sp;  summ_ri->size = slen;  summ_ri->alloc = l - (slong) (sp - summ);
        f_ri->d    = fp;  f_ri->size    = flen;  f_ri->alloc    = l - (slong) (fp - f);

        /* rop = summ / f  (mod p^N) */
        if (!radix_integer_divmod_limbs(rop, summ_ri, f_ri, lN, radix))
            flint_throw(FLINT_ERROR, "_padic_radix_exp_rectangular: "
                "factorial denominator is not invertible after stripping\n");
    }

    radix_integer_mod_digits(rop, rop, N, radix);

    flint_free(scratch);

    return GR_SUCCESS;
}

