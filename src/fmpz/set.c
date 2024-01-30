/*
    Copyright (C) 2009, 2016 William Hart
    Copyright (C) 2011, 2017 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2018, 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "double_extras.h"
#include "gmpcompat.h"
#include "fmpz.h"

void
fmpz_set(fmpz_t f, const fmpz_t g)
{
    if (f == g)
        return;                 /* aliased inputs */

    if (!COEFF_IS_MPZ(*g))      /* g is small */
    {
        _fmpz_demote(f);
        *f = *g;
    }
    else                        /* g is large */
    {
        __mpz_struct * mf = _fmpz_promote(f);
        mpz_set(mf, COEFF_TO_PTR(*g));
    }
}

#if FLINT64   /* 2^53 */
#define DOUBLE_MAX 9007199254740992.0
#define DOUBLE_MIN -9007199254740992.0
#else
#define DOUBLE_MAX COEFF_MAX
#define DOUBLE_MIN COEFF_MIN
#endif

void
fmpz_set_d(fmpz_t f, double c)
{
    if (c >= DOUBLE_MIN && c <= DOUBLE_MAX)
    {
        _fmpz_demote(f);
        /* guaranteed to fit, since c gets truncated */
        *f = (slong) c;
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(f);
        mpz_set_d(z, c);
        _fmpz_demote_val(f);
    }
}

void fmpz_set_d_2exp(fmpz_t f, double m, slong exp)
{
   int exp2;

   m = frexp(m, &exp2);
   exp += exp2;

   if (exp >= 53)
   {
      fmpz_set_d(f, m * ldexp(1.0, 53));
      fmpz_mul_2exp(f, f, exp - 53);
   } else if (exp < 0)
      fmpz_set_ui(f, 0);
   else
      fmpz_set_d(f, d_mul_2exp_inrange(m, exp));
}

void
fmpz_set_mpf(fmpz_t f, const mpf_t x)
{
    int check;

    check = flint_mpf_fits_slong_p(x);

    if (check)
    {
        slong cx = flint_mpf_get_si(x);
        fmpz_set_si(f, cx);
    }
    else
    {
        __mpz_struct *z = _fmpz_promote(f);
        mpz_set_f(z, x);
    }
}

void
fmpz_set_mpz(fmpz_t f, const mpz_t x)
{
    int size = (slong) x->_mp_size;

    if (size == 0)             /* x is zero */
    {
        fmpz_zero(f);
    }
    else if (size == 1)        /* x is positive and 1 limb */
    {
        fmpz_set_ui(f, flint_mpz_get_ui(x));
    }
    else if (size == -1)       /* x is negative and 1 limb */
    {
        ulong uval = flint_mpz_get_ui(x);
        if (uval <= COEFF_MAX)  /* x is small */
        {
            _fmpz_demote(f);
            *f = -uval;
        }
        else                    /* x is large but one limb */
        {
            __mpz_struct * mf = _fmpz_promote(f);
            flint_mpz_set_ui(mf, uval);
            mpz_neg(mf, mf);
        }
    }
    else                        /* x is more than one limb */
    {
        __mpz_struct * mf = _fmpz_promote(f);
        mpz_set(mf, x);
    }
}

/*
    Given an array of limbs "c" representing a integer mod 2^(FLINT_BITS*n),
    set "f" to the symmetric remainder with the halfway point
    2^(FLINT_BITS*n/2) mapping to -2^(FLINT_BITS*n/2)
*/
void fmpz_set_signed_ui_array(fmpz_t f, const ulong * c, slong n)
{
    ulong csign;

    FLINT_ASSERT(n > 0);

    csign = FLINT_SIGN_EXT(c[n - 1]);

    while (n > 0 && c[n - 1] == csign)
        n--;

    if (n < 2)
    {
        if (csign == 0)
            fmpz_set_ui(f, c[0]);
        else if (c[0] != 0)
            fmpz_neg_ui(f, -c[0]);
        else
            fmpz_neg_uiui(f, 1, 0);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(f);
        mp_limb_t * zd = FLINT_MPZ_REALLOC(z, n);

        if (csign == 0)
        {
            flint_mpn_copyi(zd, c, n);
            z->_mp_size = n;
        }
        else
        {
            if (mpn_neg(zd, c, n))
            {
                FLINT_ASSERT(zd[n - 1] != 0);
                z->_mp_size = -n;
            }
            else
            {
                zd = FLINT_MPZ_REALLOC(z, n + 1);
                zd[n] = 1;
                z->_mp_size = -(n + 1);
            }
        }
    }
}

void
fmpz_set_signed_uiuiui(fmpz_t r, ulong hi, ulong mid, ulong lo)
{
    int negate = 0;

    if ((slong) hi < 0)
    {
        hi = -hi - ((lo != 0) || (mid != 0));
        mid = -mid - (lo != 0);
        lo = -lo;
        negate = 1;
    }

    if (hi == 0)
    {
        if (negate)
            fmpz_neg_uiui(r, mid, lo);
        else
            fmpz_set_uiui(r, mid, lo);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(r);
        if (z->_mp_alloc < 3)
            mpz_realloc2(z, 3 * FLINT_BITS);
        z->_mp_d[0] = lo;
        z->_mp_d[1] = mid;
        z->_mp_d[2] = hi;
        z->_mp_size = negate ? -3 : 3;
    }
}

/*
    Given an array of limbs "in" representing a non negative integer,
    set "out" to this integer.
*/
void fmpz_set_ui_array(fmpz_t out, const ulong * in, slong in_len)
{
    slong size = in_len;
    FLINT_ASSERT(in_len > 0);

    /* find end of zero extension */
    while (size > WORD(1) && in[size - 1] == UWORD(0))
        size--;

    /* copy limbs */
    if (size == WORD(1))
    {
        fmpz_set_ui(out, in[0]);
    }
    else
    {
        __mpz_struct * mpz = _fmpz_promote(out);
        if (mpz->_mp_alloc < size)
            mpz_realloc2(mpz, FLINT_BITS * size);
        mpz->_mp_size = size;
        flint_mpn_copyi(mpz->_mp_d, in, size);
    }
}
