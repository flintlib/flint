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

#define MAXN 24

/* the mantissa {d, n} 2^(exp - n FLINT_BITS) as an mpfr number (exactly) */
static void
_mpn_get_mpfr(mpfr_t f, nn_srcptr d, slong n, slong exp)
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

/* the mpfr rounding mode of a sign-free mode */
static mpfr_rnd_t
_mpfr_rnd(arf_rnd_t rnd)
{
    if (rnd == ARF_RND_DOWN) return MPFR_RNDZ;
    if (rnd == ARF_RND_UP) return MPFR_RNDA;
    return MPFR_RNDN;
}

/* quotients, reciprocals and roots (of positive mantissas, the signs the
   caller's) against mpfr: the rounded results equal to mpfr's correct
   rounding with the same inexact flag; the _approx results within an
   ulp of the value */
TEST_FUNCTION_START(_arf_mpn_div, state)
{
    slong iter;

    for (iter = 0; iter < 5000 * flint_test_multiplier(); iter++)
    {
        ulong x[MAXN], y[MAXN], r[MAXN + 2];
        mpfr_t fx, fy, fv, fr, fd, fe;
        slong xn, yn, prec, rn, i, xexp, yexp, zexp, exp, base;
        int op, inexact, inexact_ref, t, use_approx;
        arf_rnd_t rnd;

        int small = n_randint(state, 2), tie = 0;
        slong tie_prec = 0, alias, yn_stripped;
        nn_ptr z;

        /* the register and chain paths: short operands half of the time */
        xn = 1 + n_randint(state, small ? 3 : MAXN);
        yn = 1 + n_randint(state, small ? 3 : MAXN);
        _arf_mpn_test_randtest_mantissa(x, xn, state);
        _arf_mpn_test_randtest_mantissa(y, yn, state);
        if (n_randint(state, 4) == 0)
        {
            for (i = 0; i < xn - 1; i++) x[i] = n_randint(state, 2) ? ~UWORD(0) : 0;
            for (i = 0; i < yn - 1; i++) y[i] = n_randint(state, 2) ? ~UWORD(0) : 0;
        }
        /* padded mantissas (zero low limbs, as nfloat's) */
        if (n_randint(state, 4) == 0)
        {
            slong k = n_randint(state, xn);
            for (i = 0; i < k; i++) x[i] = 0;
            k = n_randint(state, yn);
            for (i = 0; i < k; i++) y[i] = 0;
        }
        /* powers of two, squares and exact quotients now and then */
        if (n_randint(state, 8) == 0)
        {
            if (n_randint(state, 2))
                yn = 1;
            y[yn - 1] = (UWORD(1) << (FLINT_BITS - 1));
            for (i = 0; i < yn - 1; i++) y[i] = 0;
        }
        if (n_randint(state, 8) == 0)
        {
            if (n_randint(state, 2))
                xn = 1;
            x[xn - 1] = (UWORD(1) << (FLINT_BITS - 1));
            for (i = 0; i < xn - 1; i++) x[i] = 0;
        }
        if (n_randint(state, 8) == 0 && xn <= MAXN / 2)
        {
            /* x = a square: y^2 for a random y of half the length */
            slong hn = (xn + 1) / 2;
            ulong tt[MAXN];
            _arf_mpn_test_randtest_mantissa(tt, hn, state);
            if (hn >= 2 && n_randint(state, 2))
            {
                /* a tie of the root at (hn - 1) limbs: the half bit alone
                   below the result */
                tt[0] = UWORD(1) << (FLINT_BITS - 1);
                tie = 1;
                tie_prec = (hn - 1) * FLINT_BITS;
            }
            flint_mpn_sqr(x, tt, hn);
            xn = 2 * hn;
            if (x[xn - 1] == 0) xn--;
            if (!(x[xn - 1] >> (FLINT_BITS - 1))) { mpn_lshift(x, x, xn, 1); }
        }

        xexp = (slong) n_randint(state, 200) - 100;
        yexp = (slong) n_randint(state, 200) - 100;
        op = n_randint(state, 4);   /* 0 div, 1 inv, 2 sqrt, 3 rsqrt */
        if (op == 0 && n_randint(state, 8) == 0)
        {
            /* an exact quotient: x = y * q */
            slong qn = 1 + n_randint(state, FLINT_MAX(1, MAXN - yn));
            ulong q[MAXN];
            _arf_mpn_test_randtest_mantissa(q, qn, state);
            if (qn + yn > MAXN) qn = MAXN - yn;
            if (qn >= 2 && n_randint(state, 2))
            {
                /* a tie of the quotient at (qn - 1) limbs */
                q[0] = UWORD(1) << (FLINT_BITS - 1);
                tie = 1;
                tie_prec = (qn - 1) * FLINT_BITS;
            }
            if (qn >= 1)
            {
                if (qn >= yn) flint_mpn_mul(x, q, qn, y, yn); else flint_mpn_mul(x, y, yn, q, qn);
                xn = qn + yn;
                if (x[xn - 1] == 0) xn--;
                if (!(x[xn - 1] >> (FLINT_BITS - 1))) mpn_lshift(x, x, xn, 1);
            }
        }
        t = xexp & 1;

        mpfr_init2(fx, 64); mpfr_init2(fy, 64); mpfr_init2(fv, 64); mpfr_init2(fr, 64); mpfr_init2(fd, 64); mpfr_init2(fe, 64);
        _mpn_get_mpfr(fx, x, xn, xexp);
        _mpn_get_mpfr(fy, y, yn, yexp);

        prec = 1 + n_randint(state, n_randint(state, 4) ? 3 * FLINT_BITS : MAXN * FLINT_BITS);
        if (tie && (op == 0 || op == 2))
            prec = tie_prec;
        else if (n_randint(state, 4) == 0)
            prec = FLINT_BITS * (1 + n_randint(state, 3));   /* limb boundaries */
        rn = (prec + FLINT_BITS - 1) / FLINT_BITS;
        use_approx = n_randint(state, 3) == 0;
        rnd = _arf_mpn_test_rnd_signfree(state);
        r[rn] = 0xdeadbeef;

        /* the divisor's nonzero limbs (the kernels strip padding) */
        yn_stripped = yn;
        for (i = 0; y[i] == 0; i++)
            yn_stripped--;

        /* the output aliasing an operand (the kernels allow it) */
        alias = n_randint(state, 3);
        z = (alias == 1) ? x : (alias == 2) ? y : r;

        /* the natural exponent of the result */
        switch (op)
        {
            case 0: base = xexp - yexp; break;
            case 1: base = -yexp; break;
            case 2: base = (xexp + t) >> 1; break;
            default: base = (1 - t) - ((xexp - t) >> 1); break;
        }

        inexact = 0;
        if (use_approx)
        {
            switch (op)
            {
                case 0: _arf_mpn_div_approx(z, &zexp, x, xn, y, yn, prec); break;
                case 1: _arf_mpn_inv_approx(z, &zexp, y, yn, prec); break;
                case 2: _arf_mpn_sqrt_approx(z, &zexp, x, xn, t, prec); break;
                default: _arf_mpn_rsqrt_approx(z, &zexp, x, xn, t, prec); break;
            }
        }
        else
        {
            switch (op)
            {
                case 0: inexact = _arf_mpn_div(z, &zexp, x, xn, y, yn, prec, rnd); break;
                case 1: inexact = _arf_mpn_inv(z, &zexp, y, yn, prec, rnd); break;
                case 2: inexact = _arf_mpn_sqrt(z, &zexp, x, xn, t, prec, rnd); break;
                default: inexact = _arf_mpn_rsqrt(z, &zexp, x, xn, t, prec, rnd); break;
            }
        }
        exp = base + zexp;
        if (z != r)
            flint_mpn_copyi(r, z, rn);

        if (r[rn] != 0xdeadbeef)
        {
            flint_printf("FAIL: wrote beyond the buffer (op %d)\n", op);
            flint_abort();
        }

        _mpn_get_mpfr(fr, r, rn, exp);

        if (!(r[rn - 1] >> (FLINT_BITS - 1)) || (!use_approx &&
            (r[0] & ((rn * FLINT_BITS - prec == 0) ? 0 : ((UWORD(1) << (rn * FLINT_BITS - prec)) - 1))) != 0))
        {
            flint_printf("FAIL: form (op %d, xn %wd, yn %wd, prec %wd, approx %d, rnd %d, alias %wd)\n", op, xn, yn, prec, use_approx, (int) rnd, alias);
            flint_printf("r ="); for (i = rn - 1; i >= 0; i--) flint_printf(" %wx", r[i]); flint_printf("\n");
            mpfr_printf("fx = %.40Re\nfy = %.40Re\n", fx, fy);
            flint_abort();
        }

        if (use_approx)
        {
            /* the value to many bits: |r - v| below an ulp of the rn
               limbs (plus the truncation allowance), the truncation
               itself for the square root and for a divisor of a single
               limb, within (1 + 2^-37) ulp for the reciprocal root */
            int exact_trunc = (op == 2 || (op <= 1 && yn_stripped == 1));
            mpfr_set_prec(fv, (rn + xn + yn) * FLINT_BITS + 128);
            switch (op)
            {
                case 0: mpfr_div(fv, fx, fy, MPFR_RNDN); break;
                case 1: mpfr_ui_div(fv, 1, fy, MPFR_RNDN); break;
                case 2: mpfr_sqrt(fv, fx, MPFR_RNDN); break;
                default: mpfr_rec_sqrt(fv, fx, MPFR_RNDN); break;
            }
            mpfr_set_prec(fd, (rn + xn + yn) * FLINT_BITS + 128);
            mpfr_sub(fd, fr, fv, MPFR_RNDN);
            mpfr_mul_2si(fd, fd, -(exp - rn * FLINT_BITS), MPFR_RNDN);  /* in ulps */
            if (op <= 1)
                mpfr_set_ui_2exp(fe, 1, -40, MPFR_RNDN);
            else
                mpfr_set_ui_2exp(fe, 1, -37, MPFR_RNDN);
            mpfr_add_ui(fe, fe, 1, MPFR_RNDN);
            /* the reference's rounding: below 2^-100 of an ulp */
            mpfr_set_ui_2exp(fv, 1, -100, MPFR_RNDN);
            mpfr_add(fe, fe, fv, MPFR_RNDN);
            /* -fe < fd <= fe, and fd <= 0 (+ the reference's rounding) for
               a truncation */
            mpfr_neg(fe, fe, MPFR_RNDN);
            if (mpfr_cmp(fd, fe) <= 0)
                exact_trunc = -1;
            mpfr_neg(fe, fe, MPFR_RNDN);
            if (mpfr_cmp(fd, fe) > 0)
                exact_trunc = -1;
            if (exact_trunc == 1 && mpfr_cmp(fd, fv) > 0)
                exact_trunc = -1;
            if (exact_trunc == -1)
            {
                flint_printf("FAIL: approx (op %d, xn %wd, yn %wd, prec %wd)\n", op, xn, yn, prec);
                mpfr_printf("fd = %.20Re\nfr = %.40Re\n", fd, fr);
                flint_abort();
            }
        }
        else
        {
            mpfr_set_prec(fv, prec);
            switch (op)
            {
                case 0: inexact_ref = mpfr_div(fv, fx, fy, _mpfr_rnd(rnd)); break;
                case 1: inexact_ref = mpfr_ui_div(fv, 1, fy, _mpfr_rnd(rnd)); break;
                case 2: inexact_ref = mpfr_sqrt(fv, fx, _mpfr_rnd(rnd)); break;
                default: inexact_ref = mpfr_rec_sqrt(fv, fx, _mpfr_rnd(rnd)); break;
            }

            if ((inexact != 0) != (inexact_ref != 0) || !mpfr_equal_p(fr, fv))
            {
                flint_printf("FAIL: rounded (op %d, xn %wd, yn %wd, prec %wd, rnd %d)\n", op, xn, yn, prec, (int) rnd);
                mpfr_printf("fv = %.40Re (inexact %d)\nfr = %.40Re (inexact %d)\n", fv, inexact_ref, fr, inexact);
                flint_printf("r ="); for (i = rn - 1; i >= 0; i--) flint_printf(" %wx", r[i]); flint_printf("\n");
                flint_abort();
            }
        }
        mpfr_clear(fx); mpfr_clear(fy); mpfr_clear(fv); mpfr_clear(fr); mpfr_clear(fd); mpfr_clear(fe);
    }

    mpfr_free_cache();
    TEST_FUNCTION_END(state);
}
