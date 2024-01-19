/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2009 Andy Novocin
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2015 Kushagra Singh
    Copyright (C) 2018, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#if defined(_WIN64) || defined(__mips64)
# include <stdint.h> /* to enable mpfr_set_sj in mpfr.h */
#endif
#include <mpfr.h>
#include "mpn_extras.h"
#include "gmpcompat.h"
#include "nmod.h"
#include "fmpz.h"

double
fmpz_get_d_2exp(slong * exp, const fmpz_t f)
{
    fmpz d = *f;

    if (!COEFF_IS_MPZ(d))
    {
        ulong d_abs;

        if (d == WORD(0))
        {
            (*exp) = WORD(0);
            return 0.0;
        }

        d_abs = FLINT_ABS(d);
        *exp = FLINT_BIT_COUNT(d_abs);

        if (d < WORD(0))
            return flint_mpn_get_d((mp_limb_t *) &d_abs, WORD(1), WORD(-1), -*exp);
        else
            return flint_mpn_get_d((mp_limb_t *) &d, WORD(1), WORD(1), -*exp);
    }
    else
    {
       long exp2;
       double m = mpz_get_d_2exp(&exp2, COEFF_TO_PTR(d));
       *exp = exp2;
       return m;
    }
}

#if FLINT64   /* 2^53 */
#define DOUBLE_MAX WORD(9007199254740992)
#define DOUBLE_MIN WORD(-9007199254740992)
#else
#define DOUBLE_MAX COEFF_MAX
#define DOUBLE_MIN COEFF_MIN
#endif

double
fmpz_get_d(const fmpz_t f)
{
    fmpz c = *f;

    if (c >= DOUBLE_MIN && c <= DOUBLE_MAX)
    {
        return (double) c;
    }
    else if (!COEFF_IS_MPZ(c))
    {
        mp_limb_t d;

        if (c > 0)
        {
            d = c;
            return flint_mpn_get_d(&d, 1, 1, 0);
        }
        else
        {
            d = -c;
            return flint_mpn_get_d(&d, 1, -1, 0);
        }
    }
    else
        return mpz_get_d(COEFF_TO_PTR(c));
}

void
fmpz_get_mpf(mpf_t x, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
        flint_mpf_set_si(x, *f);      /* set x to small value */
    else
        mpf_set_z(x, COEFF_TO_PTR(*f)); /* set x to large value */
}

void
fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd)
{
    if (!COEFF_IS_MPZ(*f))
#if defined(_WIN64) || defined(__mips64)
        mpfr_set_sj(x, *f, rnd);
#else
        mpfr_set_si(x, *f, rnd);    /* set x to small value */
#endif
    else
        mpfr_set_z(x, COEFF_TO_PTR(*f), rnd);   /* set x to large value */
}

int
fmpz_get_mpn(mp_ptr *n, fmpz_t n_in)
{
    mp_limb_t n_size;
    mp_ptr temp;

    n_size = fmpz_size(n_in);
    *n = flint_malloc(n_size * sizeof(mp_limb_t));

    if (n_size <= 1)
    {
        (*n)[0] = fmpz_get_ui(n_in);
        return 1;
    }
    else
    {
        temp = COEFF_TO_PTR(*n_in)->_mp_d;
        flint_mpn_copyi(*n, temp, n_size);
        return n_size;
    }
}

void
fmpz_get_mpz(mpz_t x, const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))
        flint_mpz_set_si(x, *f);      /* set x to small value */
    else
        mpz_set(x, COEFF_TO_PTR(*f));   /* set x to large value */
}

mp_limb_t fmpz_get_nmod(const fmpz_t aa, nmod_t mod)
{
    fmpz A = *aa;
    mp_limb_t r, SA, UA;

    if (!COEFF_IS_MPZ(A))
    {
        SA = FLINT_SIGN_EXT(A);
        UA = FLINT_ABS(A);
        NMOD_RED(r, UA, mod);
    }
    else
    {
        mpz_srcptr a = COEFF_TO_PTR(A);
        mp_srcptr ad = a->_mp_d;
        slong an = a->_mp_size;

        if (an < 0)
        {
            SA = -UWORD(1);
            an = -an;
        }
        else
        {
            SA = 0;
        }

        if (an < 5)
        {
            r = 0;
            while (an > 0)
            {
                NMOD_RED2(r, r, ad[an - 1], mod);
                an--;
            }
        }
        else
        {
            r = mpn_mod_1(ad, an, mod.n);
        }
    }

    return (SA == 0 || r == 0) ? r : (mod.n - r);
}

slong
fmpz_get_si(const fmpz_t f)
{
    return (!COEFF_IS_MPZ(*f) ? *f : flint_mpz_get_si(COEFF_TO_PTR(*f)));
}

void fmpz_get_signed_ui_array(mp_limb_t * r, slong n, const fmpz_t x)
{
    int neg;
    slong i, sz;

    FLINT_ASSERT(n > 0);

    if (!COEFF_IS_MPZ(*x))
    {
        neg = *x < 0;
        r[0] = FLINT_ABS(*x);
        i = 1;
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        neg = p->_mp_size < 0;
        sz = FLINT_ABS(p->_mp_size);

        for (i = 0; i < n && i < sz; i++)
            r[i] = p->_mp_d[i];
    }

    for ( ; i < n; i++)
        r[i] = 0;

    if (neg)
        mpn_neg(r, r, n);
}

void fmpz_get_signed_uiui(mp_limb_t * hi, mp_limb_t * lo, const fmpz_t x)
{
    ulong r0, r1, s;

    if (!COEFF_IS_MPZ(*x))
    {
        r0 = *x;
        r1 = FLINT_SIGN_EXT(r0);
    }
    else
    {
        __mpz_struct * p = COEFF_TO_PTR(*x);
        s = -(ulong)(p->_mp_size < 0);
        r0 = p->_mp_d[0];
        if (p->_mp_size > 1 || p->_mp_size < -1)
            r1 = p->_mp_d[1];
        else
            r1 = 0;

        sub_ddmmss(r1, r0, r1^s, r0^s, s, s);
    }

    *lo = r0;
    *hi = r1;
}

/*
    Assuming that "in" is non negative and has a limb count <= out_len,
    write the limbs to "out" and zero extend to "out_len" limbs.
*/
void fmpz_get_ui_array(ulong * out, slong out_len, const fmpz_t in)
{
    slong size = 0;
    FLINT_ASSERT(out_len > 0);

    /* copy limbs */
    if (fmpz_abs_fits_ui(in))
    {
        *out++ = fmpz_get_ui(in);
        size++;
    } else
    {
        __mpz_struct * mpz = COEFF_TO_PTR(*in);
        FLINT_ASSERT(mpz->_mp_size <= out_len);
        while (size < mpz->_mp_size)
            *out++ = mpz->_mp_d[size++];
    }

    /* zero extend to out_len */
    while (size++ < out_len)
        *out++ = UWORD(0);
}

ulong
fmpz_get_ui(const fmpz_t f)
{
    if (!COEFF_IS_MPZ(*f))      /*value is small */
        return (*f < WORD(0) ? -*f : *f);
    else                        /* value is large */
        return flint_mpz_get_ui(COEFF_TO_PTR(*f));
}
