/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <mpfr.h>
#include "fmpz.h"
#include "mpn_helpers.h"

#define CMAXN 12

/* the mantissa {d, n} 2^(exp - n FLINT_BITS) as an mpfr number (exactly) */
static void
_cmpn_get_mpfr(mpfr_t f, nn_srcptr d, slong n, slong exp)
{
    fmpz_t t;
    mpz_t m;
    fmpz_init(t);
    mpz_init(m);
    fmpz_set_ui_array(t, d, n);
    fmpz_get_mpz(m, t);
    mpfr_set_prec(f, n * FLINT_BITS);
    mpfr_set_z_2exp(f, m, exp - n * FLINT_BITS, MPFR_RNDN);
    fmpz_clear(t);
    mpz_clear(m);
}

/* |r - v| in ulps of the rn-limb result at exponent exp, into d */
static void
_ulp_error(mpfr_t d, const mpfr_t r, const mpfr_t v, slong rn, slong exp)
{
    mpfr_sub(d, r, v, MPFR_RNDN);
    mpfr_abs(d, d, MPFR_RNDN);
    mpfr_mul_2si(d, d, -(exp - rn * FLINT_BITS), MPFR_RNDN);
}

/* the complex kernels against mpfr formulas at high precision: each part
   within about an ulp of the exact value (the intermediates are held
   with a guard limb; the bound here leaves a margin over the observed
   maximum, and matches the 2^(2 - prec) of acf's tests) */
TEST_FUNCTION_START(_arf_mpn_complex, state)
{
    slong iter;

    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        ulong a[CMAXN], b[CMAXN], c[CMAXN], d[CMAXN], re[CMAXN], im[CMAXN];
        slong an, bn, cn, dn, ea, eb, ec, ed, prec, rn, reexp, imexp, hp;
        int op, sa, sb, sc, sd, reneg, imneg;
        mpfr_t fa, fb, fc, fd, fs, fu, fv, fr, fi, fe, tol;

        an = 1 + n_randint(state, CMAXN);
        bn = 1 + n_randint(state, CMAXN);
        cn = 1 + n_randint(state, CMAXN);
        dn = 1 + n_randint(state, CMAXN);
        _arf_mpn_test_randtest_mantissa(a, an, state);
        _arf_mpn_test_randtest_mantissa(b, bn, state);
        _arf_mpn_test_randtest_mantissa(c, cn, state);
        _arf_mpn_test_randtest_mantissa(d, dn, state);
        if (n_randint(state, 8) == 0)
        {
            /* powers of two and equal parts: cancellations, roots just
               below a power of two */
            a[an - 1] = UWORD(1) << (FLINT_BITS - 1);
            flint_mpn_zero(a, an - 1);
            if (n_randint(state, 2))
            {
                flint_mpn_copyi(b, a, an);
                bn = an;
            }
            if (n_randint(state, 2))
            {
                flint_mpn_copyi(c, a, an);
                cn = an;
                flint_mpn_copyi(d, b, bn);
                dn = bn;
            }
        }
        ea = (slong) n_randint(state, 100) - 50;
        eb = (slong) n_randint(state, 100) - 50;
        ec = (slong) n_randint(state, 100) - 50;
        ed = (slong) n_randint(state, 100) - 50;
        if (n_randint(state, 2))
        {
            eb = ea + (slong) n_randint(state, 7) - 3;
            ed = ec + (slong) n_randint(state, 7) - 3;
        }
        sa = n_randint(state, 2);
        sb = n_randint(state, 2);
        sc = n_randint(state, 2);
        sd = n_randint(state, 2);
        prec = 1 + n_randint(state, CMAXN * FLINT_BITS);
        rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        op = n_randint(state, 4);   /* 0 inv, 1 div, 2 sqrt, 3 rsqrt */

        hp = 4 * CMAXN * FLINT_BITS + 256;
        mpfr_init2(fa, hp); mpfr_init2(fb, hp); mpfr_init2(fc, hp); mpfr_init2(fd, hp);
        mpfr_init2(fs, hp); mpfr_init2(fu, hp); mpfr_init2(fv, hp);
        mpfr_init2(fr, hp); mpfr_init2(fi, hp); mpfr_init2(fe, hp); mpfr_init2(tol, 64);
        _cmpn_get_mpfr(fa, a, an, ea);
        _cmpn_get_mpfr(fb, b, bn, eb);
        _cmpn_get_mpfr(fc, c, cn, ec);
        _cmpn_get_mpfr(fd, d, dn, ed);
        if (sa) mpfr_neg(fa, fa, MPFR_RNDN);
        if (sb) mpfr_neg(fb, fb, MPFR_RNDN);
        if (sc) mpfr_neg(fc, fc, MPFR_RNDN);
        if (sd) mpfr_neg(fd, fd, MPFR_RNDN);

        reneg = imneg = 0;
        reexp = imexp = 0;

        if (op == 0)
        {
            /* (a - b i) / (a^2 + b^2): magnitudes */
            _arf_mpn_cinv_approx(re, &reexp, im, &imexp, a, an, ea, b, bn, eb, prec);
            mpfr_sqr(fs, fa, MPFR_RNDN);
            mpfr_sqr(fu, fb, MPFR_RNDN);
            mpfr_add(fs, fs, fu, MPFR_RNDN);
            mpfr_div(fu, fa, fs, MPFR_RNDN);
            mpfr_div(fv, fb, fs, MPFR_RNDN);
            mpfr_abs(fu, fu, MPFR_RNDN);
            mpfr_abs(fv, fv, MPFR_RNDN);
        }
        else if (op == 1)
        {
            _arf_mpn_cdiv_approx(re, &reexp, &reneg, im, &imexp, &imneg, a, an, ea, sa, b, bn, eb, sb, c, cn, ec, sc, d, dn, ed, sd, prec);
            /* (a c + b d, b c - a d) / (c^2 + d^2), signed */
            mpfr_sqr(fs, fc, MPFR_RNDN);
            mpfr_sqr(fe, fd, MPFR_RNDN);
            mpfr_add(fs, fs, fe, MPFR_RNDN);
            mpfr_mul(fu, fa, fc, MPFR_RNDN);
            mpfr_mul(fe, fb, fd, MPFR_RNDN);
            mpfr_add(fu, fu, fe, MPFR_RNDN);
            mpfr_mul(fv, fb, fc, MPFR_RNDN);
            mpfr_mul(fe, fa, fd, MPFR_RNDN);
            mpfr_sub(fv, fv, fe, MPFR_RNDN);
            mpfr_div(fu, fu, fs, MPFR_RNDN);
            mpfr_div(fv, fv, fs, MPFR_RNDN);
        }
        else
        {
            /* r = |z|, q = sqrt((r + |a|) / 2): sqrt(z) = (q, b / (2 q))
               for a >= 0, (|b| / (2 q), sgn(b) q) for a < 0; 1 / sqrt(z)
               the conjugate divided by r */
            if (op == 2)
                _arf_mpn_csqrt_approx(re, &reexp, im, &imexp, a, an, ea, sa, b, bn, eb, prec);
            else
                _arf_mpn_crsqrt_approx(re, &reexp, im, &imexp, a, an, ea, sa, b, bn, eb, prec);
            mpfr_hypot(fs, fa, fb, MPFR_RNDN);
            mpfr_abs(fe, fa, MPFR_RNDN);
            mpfr_add(fe, fs, fe, MPFR_RNDN);
            mpfr_div_2ui(fe, fe, 1, MPFR_RNDN);
            mpfr_sqrt(fe, fe, MPFR_RNDN);     /* q */
            mpfr_abs(fu, fb, MPFR_RNDN);
            mpfr_div(fu, fu, fe, MPFR_RNDN);
            mpfr_div_2ui(fu, fu, 1, MPFR_RNDN);   /* |b| / (2 q) */
            if (!sa)
                mpfr_swap(fu, fe);    /* fu = re = q, fe = im = |b| / 2q */
            else
                ;                     /* fu = re = |b| / 2q, fe = im = q */
            mpfr_set(fv, fe, MPFR_RNDN);
            if (op == 3)
            {
                mpfr_div(fu, fu, fs, MPFR_RNDN);
                mpfr_div(fv, fv, fs, MPFR_RNDN);
            }
        }

        /* the parts as computed */
        if (reexp == WORD_MIN)
            mpfr_set_ui(fr, 0, MPFR_RNDN);
        else
            _cmpn_get_mpfr(fr, re, rn, reexp);
        if (imexp == WORD_MIN)
            mpfr_set_ui(fi, 0, MPFR_RNDN);
        else
            _cmpn_get_mpfr(fi, im, rn, imexp);
        mpfr_set_prec(fe, hp);
        if (reneg) mpfr_neg(fr, fr, MPFR_RNDN);
        if (imneg) mpfr_neg(fi, fi, MPFR_RNDN);

        if (reexp != WORD_MIN && !(re[rn - 1] >> (FLINT_BITS - 1)))
            TEST_FUNCTION_FAIL("real part not normalised (op %d)\n", op);
        if (imexp != WORD_MIN && !(im[rn - 1] >> (FLINT_BITS - 1)))
            TEST_FUNCTION_FAIL("imaginary part not normalised (op %d)\n", op);

        mpfr_set_d(tol, 4.0, MPFR_RNDN);

        /* a cancelled part must be exactly zero */
        if ((reexp == WORD_MIN && !mpfr_zero_p(fu)) || (imexp == WORD_MIN && !mpfr_zero_p(fv)))
        {
            flint_printf("FAIL: spurious zero (op %d, an %wd, bn %wd, cn %wd, dn %wd, prec %wd, reexp %wd imexp %wd)\n", op, an, bn, cn, dn, prec, reexp, imexp);
            mpfr_printf("fa = %.20Re fb = %.20Re\nfc = %.20Re fd = %.20Re\nfu = %.20Re\nfv = %.20Re\n", fa, fb, fc, fd, fu, fv);
            flint_abort();
        }

        /* the quotient's numerator parts may cancel: its error is
           relative to the magnitude of the result, as nfloat's */
        if (op == 1)
        {
            slong m = FLINT_MAX(reexp, imexp);
            if (reexp != WORD_MIN) reexp = m;
            if (imexp != WORD_MIN) imexp = m;
        }

        if (reexp != WORD_MIN)
        {
            _ulp_error(fe, fr, fu, rn, reexp);
            if (mpfr_cmp(fe, tol) > 0)
            {
                flint_printf("FAIL: real part (op %d, an %wd, bn %wd, cn %wd, dn %wd, prec %wd)\n", op, an, bn, cn, dn, prec);
                mpfr_printf("fr = %.30Re\nfu = %.30Re\nerr = %.10Re ulp\n", fr, fu, fe);
                flint_abort();
            }
        }
        if (imexp != WORD_MIN)
        {
            _ulp_error(fe, fi, fv, rn, imexp);
            if (mpfr_cmp(fe, tol) > 0)
            {
                flint_printf("FAIL: imaginary part (op %d, an %wd, bn %wd, cn %wd, dn %wd, prec %wd)\n", op, an, bn, cn, dn, prec);
                mpfr_printf("fi = %.30Re\nfv = %.30Re\nerr = %.10Re ulp\n", fi, fv, fe);
                flint_abort();
            }
        }

        mpfr_clear(fa); mpfr_clear(fb); mpfr_clear(fc); mpfr_clear(fd);
        mpfr_clear(fs); mpfr_clear(fu); mpfr_clear(fv); mpfr_clear(fr); mpfr_clear(fi); mpfr_clear(fe); mpfr_clear(tol);
    }

    mpfr_free_cache();
    TEST_FUNCTION_END(state);
}
