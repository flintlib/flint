/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz_vec.h"
#include "arb_hypgeom.h"

#ifdef __GNUC__
# define ldexp __builtin_ldexp
# define log __builtin_log
#else
# include <math.h>
#endif

static double
d_log2_fac(double n)
{
    return (n * (log(n) - 1.0)) * 1.4426950408889634;
}

static slong
tail_precision(slong k, double min_k, slong alen, slong blen, double log2z, double log2max, slong prec)
{
    double term_magnitude;
    slong new_prec;

    if (prec <= 128 || k <= 5 || k <= min_k)
        return prec;

    term_magnitude = k * log2z;
    if (alen != blen)
        term_magnitude += (alen - blen) * d_log2_fac(k);

    new_prec = prec - (log2max - term_magnitude) + 10;
    new_prec = FLINT_MIN(new_prec, prec);
    new_prec = FLINT_MAX(new_prec, 32);

    /* printf("term %ld, max %f, log2x %f, magn %f   new_prec %ld\n", k, log2z, log2max, term_magnitude, new_prec); */

    return new_prec;
}

/* Return approximation of log2(|x|), approximately clamped between COEFF_MIN and COEFF_MAX. */
double
arf_get_d_log2_abs_approx_clamped(const arf_t x)
{
    if (arf_is_zero(x))
    {
        return (double) COEFF_MIN;
    }
    else if (!arf_is_finite(x))
    {
        return (double) COEFF_MAX;
    }
    else if (COEFF_IS_MPZ(ARF_EXP(x)))
    {
        if (fmpz_sgn(ARF_EXPREF(x)) < 0)
            return (double) COEFF_MIN;
        else
            return (double) COEFF_MAX;
    }
    else
    {
        mp_srcptr tp;
        mp_size_t tn;
        double v;
        slong e = ARF_EXP(x);

        ARF_GET_MPN_READONLY(tp, tn, x);

        if (tn == 1)
            v = (double)(tp[0]);
        else
            v = (double)(tp[tn - 1]) + (double)(tp[tn - 2]) * ldexp(1.0, -FLINT_BITS);

        v *= ldexp(1.0, -FLINT_BITS);

        return 1.4426950408889634074 * mag_d_log_upper_bound(v) + (double) e;
    }
}

void
arb_hypgeom_sum_fmpq_arb_rs(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
{
    slong m, i, j, k, l, jlen, jbot, jtop, wp;
    double log2z, log2max, adaptive_min_k;
    int want_adaptive_precision;
    arb_t s, t;
    arb_ptr zpow;
    fmpz_t c, den;
    fmpz * cs;
    slong Nbits, acbits, bcbits, numbits, denbits;

    if (N <= 1)
    {
        if (N == 1)
            arb_one(res);
        else
            arb_zero(res);
        return;
    }

    m = n_sqrt(N);
    m = FLINT_MAX(m, 2);
    k = N - 1;
    j = k % m;
    jlen = 0;
    jbot = j;

    fmpz_init(c);
    fmpz_init(den);
    arb_init(s);
    arb_init(t);
    zpow = _arb_vec_init(m + 1);
    cs = _fmpz_vec_init(m + 1);

    fmpz_one(c);
    fmpz_one(den);
    for (i = 0; i < alen; i++)
        fmpz_mul(den, den, fmpq_denref(a + i));
    for (i = 0; i < blen; i++)
        fmpz_mul(c, c, fmpq_denref(b + i));

    if (reciprocal)
    {
        arb_mul_fmpz(zpow + m, z, den, prec);
        arb_set_fmpz(t, c);
        arb_div(zpow + m, t, zpow + m, prec);
    }
    else
    {
        arb_mul_fmpz(zpow + m, z, c, prec);
        arb_div_fmpz(zpow + m, zpow + m, den, prec);
    }

    want_adaptive_precision = N > 5;

    Nbits = FLINT_BIT_COUNT(N);
    acbits = 0;
    for (i = 0; i < alen; i++)
    {
        numbits = fmpz_bits(fmpq_numref(a + i));
        denbits = fmpz_bits(fmpq_denref(a + i));
        want_adaptive_precision = want_adaptive_precision && (FLINT_ABS(numbits - denbits) < 2);
        acbits += FLINT_MAX(denbits + Nbits, numbits) + 1;
    }

    bcbits = 0;
    for (i = 0; i < blen; i++)
    {
        numbits = fmpz_bits(fmpq_numref(b + i));
        denbits = fmpz_bits(fmpq_denref(b + i));
        want_adaptive_precision = want_adaptive_precision && (FLINT_ABS(numbits - denbits) < 2);
        bcbits += FLINT_MAX(denbits + Nbits, numbits) + 1;
    }

    log2max = 0.0;
    log2z = 0.0;
    adaptive_min_k = 0.0;
    if (want_adaptive_precision)
    {
        log2z = arf_get_d_log2_abs_approx_clamped(arb_midref(z));
        if (reciprocal)
            log2z = -log2z;

        /* Terms are increasing, so don't change the precision. */
        if (alen >= blen && log2z >= 0.0)
        {
            want_adaptive_precision = 0;
        }
        else
        {
            if (alen >= blen)
            {
                log2max = 0.0;
            }
            else
            {
                slong r = blen - alen;

                /* r = 1 -> exp(z) */
                /* r = 2 -> exp(2*z^(1/2)) */
                /* r = 3 -> exp(3*z^(1/3)) */
                /* ... */
#ifdef __GNUC__
                log2max = r * __builtin_exp(log2z * 0.693147180559945 / r) * 1.44269504088896;

                /* fixme */
                if (r == 1)
                    adaptive_min_k = __builtin_exp(log2z * log(2));
                else
                    adaptive_min_k = __builtin_exp(0.5 * log2z * log(2));
#else
                log2max = r * exp(log2z * 0.693147180559945 / r) * 1.44269504088896;

                /* fixme */
                if (r == 1)
                    adaptive_min_k = exp(log2z * log(2));
                else
                    adaptive_min_k = exp(0.5 * log2z * log(2));
#endif
            }
        }
    }

    _arb_vec_set_powers(zpow, zpow + m, m + 1, prec);

    while (k >= 0)
    {
        jlen = 1;
        jtop = jbot = k;

        if (jtop > 0)
        {
            while (jlen <= j && jlen <= 8 && jbot >= 2)
            {
                jbot--;
                jlen++;
            }
        }

        if (jbot != jtop - jlen + 1)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        /* Factors between jbot and jtop inclusive */
        if (jbot == 0 || alen == 0)
        {
            fmpz_one(cs + 0);
        }
        else
        {
            if (acbits <= FLINT_BITS - 1)
            {
                slong ac;

                ac = 1;
                for (l = 0; l < alen; l++)
                    ac *= *fmpq_denref(a + l) * (jbot - 1) + *fmpq_numref(a + l);

                fmpz_set_si(cs + 0, ac);
            }
            else
            {
                fmpz_mul_ui(cs + 0, fmpq_denref(a + 0), jbot - 1);
                fmpz_add(cs + 0, cs + 0, fmpq_numref(a + 0));

                for (l = 1; l < alen; l++)
                {
                    fmpz_mul_ui(c, fmpq_denref(a + l), jbot - 1);
                    fmpz_add(c, c, fmpq_numref(a + l));
                    fmpz_mul(cs + 0, cs + 0, c);
                }
            }
        }

        for (i = 1; i < jlen; i++)
        {
            if (alen == 0)
            {
                fmpz_set(cs + i, cs + i - 1);
            }
            else
            {
                if (acbits <= FLINT_BITS - 1)
                {
                    slong ac;

                    ac = 1;
                    for (l = 0; l < alen; l++)
                        ac *= *fmpq_denref(a + l) * (jbot + i - 1) + *fmpq_numref(a + l);

                    fmpz_mul_si(cs + i, cs + i - 1, ac);
                }
                else
                {
                    fmpz_mul_ui(cs + i, fmpq_denref(a + 0), jbot + i - 1);
                    fmpz_add(cs + i, cs + i, fmpq_numref(a + 0));

                    for (l = 1; l < alen; l++)
                    {
                        fmpz_mul_ui(c, fmpq_denref(a + l), jbot + i - 1);
                        fmpz_add(c, c, fmpq_numref(a + l));
                        fmpz_mul(cs + i, cs + i, c);
                    }

                    fmpz_mul(cs + i, cs + i, cs + i - 1);
                }
            }
        }

        if (blen != 0)
        {
            fmpz_one(den);
            for (i = jlen - 1; i >= 0; i--)
            {
                if (i != jlen - 1)
                    fmpz_mul(cs + i, cs + i, den);

                if (i != 0 || jbot != 0)
                {
                    if (bcbits <= FLINT_BITS - 1)
                    {
                        slong bc;

                        bc = 1;
                        for (l = 0; l < blen; l++)
                            bc *= *fmpq_denref(b + l) * (jbot + i - 1) + *fmpq_numref(b + l);

                        fmpz_mul_si(den, den, bc);
                    }
                    else
                    {
                        for (l = 0; l < blen; l++)
                        {
                            fmpz_mul_ui(c, fmpq_denref(b + l), jbot + i - 1);
                            fmpz_add(c, c, fmpq_numref(b + l));
                            fmpz_mul(den, den, c);
                        }
                    }
                }
            }
        }

        if (want_adaptive_precision)
            wp = tail_precision(k - jlen, adaptive_min_k, alen, blen, log2z, log2max, prec);
        else
            wp = prec;

        arb_add(t, s, zpow + j, wp);
        arb_swap(zpow + j, t);
        arb_dot_fmpz(s, NULL, 0, zpow + j - jlen + 1, 1, cs, 1, jlen, wp);
        arb_swap(zpow + j, t);

        if (blen != 0)
            arb_div_fmpz(s, s, den, wp);

        k -= jlen;
        j -= (jlen - 1);

        if (j == 0 && k >= 1)
        {
            arb_mul(s, s, zpow + m, wp);
            j = m - 1;
        }
        else
        {
            j--;
        }
    }

    arb_swap(res, s);

    _arb_vec_clear(zpow, m + 1);
    _fmpz_vec_clear(cs, m + 1);
    arb_clear(s);
    arb_clear(t);
    fmpz_clear(c);
    fmpz_clear(den);
}

