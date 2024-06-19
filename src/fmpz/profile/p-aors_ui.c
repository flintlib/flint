/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "gmpcompat.h"

#define ntests 30

#define OLD_ALBIN 1

#if OLD_ALBIN

static void
_fmpz_add_mpn_1(fmpz_t f, const ulong * glimbs, slong gsz, ulong x);

static void
_fmpz_sub_mpn_1(fmpz_t f, const ulong * glimbs, slong gsz, ulong x);

void
fmpz_add_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    mpz_ptr mf;
    slong g1 = *g;
    slong f1 = *f;

    if (!COEFF_IS_MPZ(g1))  /* g is small */
    {
        slong sz = 2;
        if (g1 >= 0)
        {
            {   /* add with jump if carry */
                ulong tmp = g1;
                g1 += x;
                if (((ulong) g1) < tmp)
                    goto carry;
            }
            if (((ulong) g1) <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = g1;
                return;
            }
nocarry:    sz = 1; /* No carry, but result does not fit in small fmpz */
carry:      if (COEFF_IS_MPZ(f1))
                mf = COEFF_TO_PTR(f1);
            else
            {
                mf = _fmpz_new_mpz();
                *f = PTR_TO_COEFF(mf);
            }
            mf->_mp_size = sz;
            mf->_mp_d[0] = g1;
            mf->_mp_d[1] = 1; /* Set carry (not used if sz = 1) */
        }
        else /* g < 0 */
        {
            g1 += x;
            if (((slong) x) >= 0 && g1 <= COEFF_MAX)
            {
                /* If x > 0 does not have its top bit set
                 * and COEFF_MIN <= g < 0, we can interpret x + g as a slong.
                 * So if the result in g1 is smaller than COEFF_MAX, it is a
                 * small fmpz. */
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = g1;
                return;
            }
            else
            {
                /* 1) If top bit is set in x, the result is going to be positive
                 *    but will be larger than COEFF_MAX since
                 *
                 *          x + g  >=  (2^63) - (2^62 - 1)  =  2^62 + 1.
                 *
                 *    However, it will be contained in one limb since g < 0.
                 *
                 * 2) If top bit is not set, then result is larger than
                 *    COEFF_MAX, and so it cannot be a small fmpz. However, it
                 *    must fit in one limb since
                 *
                 *          x + g  <=  (2^63 - 1) + (-1)  =  2^63 - 2,
                 *
                 *    which is contained in one limb. */

                goto nocarry;
            }
        }
    }
    else
    {
        mpz_ptr mg = COEFF_TO_PTR(g1);
        slong gsz = mg->_mp_size;
        ulong * glimbs = mg->_mp_d;

        if (gsz > 0)
            _fmpz_add_mpn_1(f, glimbs, gsz, x);
        else
            _fmpz_sub_mpn_1(f, glimbs, gsz, x);
    }
}

/* "Add" two number with same sign. Decide sign from g. */
static void
_fmpz_add_mpn_1(fmpz_t f, const ulong * glimbs, slong gsz, ulong x)
{
    mpz_ptr mf;
    ulong * flimbs;
    slong gabssz = FLINT_ABS(gsz);

    /* Promote f as it is guaranteed to be large */
    if (COEFF_IS_MPZ(*f))
        mf = COEFF_TO_PTR(*f);
    else
    {
        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
    }
    flimbs = mf->_mp_d;

    if (mf->_mp_alloc < (gabssz + 1)) /* Ensure result fits */
    {
        ulong * tmp = flimbs;
        flimbs = _mpz_realloc(mf, gabssz + 1);

        /* If f and g are aliased, then we need to change glimbs as well. */
        if (tmp == glimbs)
            glimbs = flimbs;
    }

    /* Use GMP to calculate result */
    flimbs[gabssz] = mpn_add_1(flimbs, glimbs, gabssz, x);

    /* flimbs[gabssz] is the carry from mpn_add_1,
     * and so gabssz + flimbs[gabssz] is valid to determine the size of f */
    mf->_mp_size = gabssz + flimbs[gabssz];
    if (gsz < 0)
    {
        /* g and x has same sign. If g is negative, we negate the result */
        mf->_mp_size = -mf->_mp_size;
    }
}

/* Subtract two limbs (they have different sign) and decide the sign via g. */
static void
_fmpz_sub_mpn_1(fmpz_t f, const ulong * glimbs, slong gsz, ulong x)
{
    mpz_ptr mf;
    ulong * flimbs;
    slong gabssz = FLINT_ABS(gsz);

    /* If size of g is 1, we have a higher probability of the result being
     * small. */
    if (gabssz == 1)
    {
        if (x <= glimbs[0]) /* Result is zero or has the same sign as g */
        {
            x = glimbs[0] - x;
L1:         if (x <= COEFF_MAX) /* Fits in small fmpz */
            {
                if (COEFF_IS_MPZ(*f))
                    _fmpz_clear_mpz(*f);
                *f = (gsz > 0) ? x : -x; /* With consideration of sign */
            }
            else /* Does not fit in small fmpz */
            {
                if (COEFF_IS_MPZ(*f))
                    mf = COEFF_TO_PTR(*f);
                else
                {
                    mf = _fmpz_new_mpz();
                    *f = PTR_TO_COEFF(mf);
                }
                mf->_mp_d[0] = x;
                mf->_mp_size = gsz; /* Sign of f is the same as for g */
            }
        }
        else /* |x| > |g|, which implies f has opposite sign of g */
        {
            /* Set x to the absolute value of |g - x|. By switching sign of
             * gsz, we can reuse the code above. */
            x -= glimbs[0];
            gsz = -gsz;
            goto L1;
        }
        return;
    }

    /* As g has more than one limb, it is a very high probability that result
     * does not fit inside small fmpz. */
    if (COEFF_IS_MPZ(*f))
        mf = COEFF_TO_PTR(*f);
    else
    {
        mf = _fmpz_new_mpz();
        *f = PTR_TO_COEFF(mf);
    }
    flimbs = mf->_mp_d;

    if (gabssz == 2)
    {
        /* Special case. Can result in a small fmpz, but as |g| > |x| the sign
         * cannot change. */
        sub_ddmmss(flimbs[1], flimbs[0], glimbs[1], glimbs[0], 0, x);
        if (flimbs[1] != 0)
        {
            /* Most likely: upper limb not zero, so we just have set the sign
             * of f to g's. */
            mf->_mp_size = gsz;
        }
        else if (flimbs[0] > COEFF_MAX)
        {
            /* Still very likely: Upper limb is zero but lower limb does not
             * fit inside a small fmpz. Sign is the same as for g, but the
             * absolute value of the size is one. */
            mf->_mp_size = (gsz > 0) ? 1 : -1;
        }
        else
        {
            /* Upper limb is zero and lower limb fits inside a small fmpz.
             * Therefore we set f to +/- flimbs[0] and clear the mpz associated
             * to f. */
            slong tmp = flimbs[0]; /* We will clear this mpz, so save first. */
            _fmpz_clear_mpz(*f);
            *f = (gsz > 0) ? tmp : -tmp;
        }
    }
    else
    {
        /* As the absolute value of g's size is larger than 2, the result won't
         * fit inside a small fmpz. */
        if (mf->_mp_alloc < gabssz) /* Ensure result fits */
        {
            /* The allocation size of g is always larger than the absolute value
             * of g. Therefore, if f's allocation size is smaller than g's
             * size, they cannot be aliased. */
            flimbs = _mpz_realloc(mf, gabssz);
        }

        mpn_sub_1(flimbs, glimbs, gabssz, x); /* Subtract via GMP */

        /* If last limb is zero, we have to set f's absolute size to one less
         * than g's. */
        mf->_mp_size = gabssz - (flimbs[gabssz - 1] == 0);
        if (gsz < 0) /* If g is negative, then f is negative as well. */
            mf->_mp_size = -mf->_mp_size;
    }
}

void
fmpz_sub_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    mpz_ptr mf;
    slong g1 = *g;
    slong f1 = *f;

    if (!COEFF_IS_MPZ(g1))  /* g is small */
    {
        slong sz = -2;
        if (g1 <= 0)
        {
            /* "add" with jump if carry */
            g1 = x - g1; /* g1 = x + |g| */
            if (((ulong) g1) < x)
                goto carry;

            if (((ulong) g1) <= COEFF_MAX)
            {
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = -g1;
                return;
            }
nocarry:    sz = -1; /* No carry, but result is not a small fmpz */
carry:      if (COEFF_IS_MPZ(f1))
                mf = COEFF_TO_PTR(f1);
            else
            {
                mf = _fmpz_new_mpz();
                *f = PTR_TO_COEFF(mf);
            }
            mf->_mp_size = sz;
            mf->_mp_d[0] = g1;
            mf->_mp_d[1] = 1; /* Set carry (not used if sz = -1) */
        }
        else
        {
            g1 = x - g1; /* -(g - x) */
            if (((slong) x) >= 0 && g1 <= COEFF_MAX)
            {
                /* If x > 0 does not have its top bit set
                 * and 0 < g <= COEFF_MAX, we can interpret x - g as a slong.
                 * So if the result in g1 is smaller than COEFF_MAX, it is a
                 * small fmpz. */
                if (COEFF_IS_MPZ(f1))
                    _fmpz_clear_mpz(f1);
                *f = -g1; /* g - x = -(x - g) */
                return;
            }
            else
            {
                /* 1) If top bit is set in x, the result is going to be negative
                 *    but will be larger than COEFF_MAX since
                 *
                 *          x - g  >=  2^63 - (2^62 - 1)  =  2^62 + 1.
                 *
                 *    However, it will be contained in one limb since g > 0.
                 *
                 * 2) If top bit is not set, then result is smaller than
                 *    COEFF_MIN, and so it cannot be a small fmpz. However, it
                 *    must fit in one limb since
                 *
                 *          x - g  <=  (2^63 - 1) - 1  =  2^63 - 2,
                 *
                 *    which is contained in one limb. */

                goto nocarry;
            }
        }
    }
    else
    {
        mpz_ptr mg = COEFF_TO_PTR(g1);
        slong gsz = mg->_mp_size;
        ulong * glimbs = mg->_mp_d;

        if (gsz > 0)
            _fmpz_sub_mpn_1(f, glimbs, gsz, x);
        else
            _fmpz_add_mpn_1(f, glimbs, gsz, x);
    }
}

#else

void fmpz_add_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c = *g;

    if (!COEFF_IS_MPZ(c))  /* g is small */
    {
        ulong sum[2];
        if (c >= WORD(0))  /* both operands non-negative */
        {
            add_ssaaaa(sum[1], sum[0], 0, c, 0, x);
            fmpz_set_uiui(f, sum[1], sum[0]);
        }
        else  /* coeff is negative, x positive */
        {
            if (-c > x)
                fmpz_set_si(f, x + c); /* can't overflow as g is small and x smaller */
            else
                fmpz_set_ui(f, x + c);  /* won't be negative and has to be less than x */
        }
    }
    else
    {
        mpz_ptr mf = _fmpz_promote(f);  /* g is already large */
        mpz_ptr mc = COEFF_TO_PTR(c);
        flint_mpz_add_ui(mf, mc, x);
        _fmpz_demote_val(f);  /* cancellation may have occurred */
    }
}

void
fmpz_sub_ui_old(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz c = *g;

    if (!COEFF_IS_MPZ(c))       /* coeff is small */
    {
        ulong sum[2];
        if (c < WORD(0))             /* g negative, x positive, so difference is negative */
        {
            add_ssaaaa(sum[1], sum[0], 0, -c, 0, x);
            fmpz_neg_uiui(f, sum[1], sum[0]);
        }
        else                    /* coeff is non-negative, x non-negative */
        {
            if (x < c)
                fmpz_set_ui(f, c - x);  /* won't be negative and is smaller than c */
            else
                fmpz_neg_ui(f, x - c);  /* positive or zero */
        }
    }
    else
    {
        mpz_ptr mc, mf;
        mf = _fmpz_promote(f);    /* g is already large */
        mc = COEFF_TO_PTR(c);
        flint_mpz_sub_ui(mf, mc, x);
        _fmpz_demote_val(f);    /* cancellation may have occurred */
    }
}

#endif

void
sample_add_new(void * arg, ulong count)
{
    fmpz *res, *a;
    ulong *b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    res = _fmpz_vec_init(ntests);
    a = _fmpz_vec_init(ntests);
    b = flint_malloc(sizeof(ulong) * ntests);

    for (ix = 0; ix < 10 * count; ix++)
    {
        for (jx = 0; jx < ntests; jx++)
        {
            fmpz_randtest(a + jx, state, bits);
            fmpz_randtest(res + jx, state, bits);
            b[jx] = n_randtest(state);
        }

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_add_ui(res + jx, a + jx, b[jx]);
        prof_stop();
    }

    _fmpz_vec_clear(res, ntests);
    _fmpz_vec_clear(a, ntests);
    flint_free(b);
    flint_rand_clear(state);
}

void
sample_add_old(void * arg, ulong count)
{
    fmpz *res, *a;
    ulong *b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    res = _fmpz_vec_init(ntests);
    a = _fmpz_vec_init(ntests);
    b = flint_malloc(sizeof(ulong) * ntests);

    for (ix = 0; ix < 10 * count; ix++)
    {
        for (jx = 0; jx < ntests; jx++)
        {
            fmpz_randtest(a + jx, state, bits);
            fmpz_randtest(res + jx, state, bits);
            b[jx] = n_randtest(state);
        }

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_add_ui_old(res + jx, a + jx, b[jx]);
        prof_stop();
    }

    _fmpz_vec_clear(res, ntests);
    _fmpz_vec_clear(a, ntests);
    flint_free(b);
    flint_rand_clear(state);
}

void
sample_sub_new(void * arg, ulong count)
{
    fmpz *res, *a;
    ulong *b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    res = _fmpz_vec_init(ntests);
    a = _fmpz_vec_init(ntests);
    b = flint_malloc(sizeof(ulong) * ntests);

    for (ix = 0; ix < 10 * count; ix++)
    {
        for (jx = 0; jx < ntests; jx++)
        {
            fmpz_randtest(a + jx, state, bits);
            fmpz_randtest(res + jx, state, bits);
            b[jx] = n_randtest(state);
        }

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_sub_ui(res + jx, a + jx, b[jx]);
        prof_stop();
    }

    _fmpz_vec_clear(res, ntests);
    _fmpz_vec_clear(a, ntests);
    flint_free(b);
    flint_rand_clear(state);
}

void
sample_sub_old(void * arg, ulong count)
{
    fmpz *res, *a;
    ulong *b;
    ulong ix, jx;
    int bits = *((int *) arg);

    FLINT_TEST_INIT(state);

    res = _fmpz_vec_init(ntests);
    a = _fmpz_vec_init(ntests);
    b = flint_malloc(sizeof(ulong) * ntests);

    for (ix = 0; ix < 10 * count; ix++)
    {
        for (jx = 0; jx < ntests; jx++)
        {
            fmpz_randtest(a + jx, state, bits);
            fmpz_randtest(res + jx, state, bits);
            b[jx] = n_randtest(state);
        }

        prof_start();
        for (jx = 0; jx < ntests; jx++)
            fmpz_sub_ui_old(res + jx, a + jx, b[jx]);
        prof_stop();
    }

    _fmpz_vec_clear(res, ntests);
    _fmpz_vec_clear(a, ntests);
    flint_free(b);
    flint_rand_clear(state);
}


slong sizes[] = { 10, 30, 60, 62, 64, 66, 80, 128, 160, 256, 512, 1024, 4096, 0 };

int
main(void)
{
    double minnew, maxnew, minold, maxold;
    int i, bits;

    flint_printf("ADD\n");
    for (i = 0; (bits = sizes[i]) != 0; i++)
    {
        prof_repeat(&minnew, &maxnew, sample_add_new, &bits);
        prof_repeat(&minold, &maxold, sample_add_old, &bits);

        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    flint_printf("\nSUB\n");
    for (i = 0; (bits = sizes[i]) != 0; i++)
    {
        prof_repeat(&minnew, &maxnew, sample_sub_new, &bits);
        prof_repeat(&minold, &maxold, sample_sub_old, &bits);

        flint_printf("%d bits:      min %.2fx,    max %.2fx\n",
                bits, minold / minnew, maxold / maxnew);
    }

    return 0;
}
