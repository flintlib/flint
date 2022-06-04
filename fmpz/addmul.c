/*
    Copyright (C) 2009 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz-impl.h"


#define MPZ_FIT_SIZE(z, nlimbs) \
    do { \
        if (z->_mp_alloc < nlimbs) \
            _mpz_realloc(z, nlimbs); \
    } while (0)

/* Will not get called with x or y small. */
void
_flint_mpz_addmul_large(mpz_ptr z, mpz_srcptr x, mpz_srcptr y, int negate)
{
    mp_size_t xn, yn, tn, zn, zn_signed, zn_new, x_sgn, y_sgn, sgn, alloc;
    mp_srcptr xd, yd;
    mp_ptr zd;
    mp_ptr td;
    mp_limb_t top;
    TMP_INIT;

    xn = x->_mp_size;
    yn = y->_mp_size;
    x_sgn = xn;
    y_sgn = yn;
    xn = FLINT_ABS(xn);
    yn = FLINT_ABS(yn);

    if (xn < yn)
    {
        mpz_srcptr t;
        mp_size_t tn;

        t = x; x = y; y = t;
        tn = xn; xn = yn; yn = tn;
        tn = x_sgn; x_sgn = y_sgn; y_sgn = tn;
    }

    if (negate)
        y_sgn = -y_sgn;

    xd = x->_mp_d;
    yd = y->_mp_d;

    /* todo: could inline code for (zn <= 5) + (xn <= 2) x (yn <= 2) */

    if (yn == 1)
    {
        if (y_sgn >= 0)
            flint_mpz_addmul_ui(z, x, yd[0]);
        else
            flint_mpz_submul_ui(z, x, yd[0]);
        return;
    }

    sgn = x_sgn ^ y_sgn;

    zn_signed = z->_mp_size;
    zn = FLINT_ABS(zn_signed);
    sgn ^= zn_signed;
    tn = xn + yn;

    if (zn == 0)
    {
        /* Cannot have aliasing here, because x and y are not
           small and therefore not 0. We can therefore resize z
           or write to zd without invalidating the inputs. */
        FLINT_ASSERT(xn >= 2 || xd[0] > COEFF_MAX);
        FLINT_ASSERT(yn >= 2 || yd[0] > COEFF_MAX);

        MPZ_FIT_SIZE(z, tn + 1);
        zd = z->_mp_d;
        zn = tn;

        if (x == y)
        {
            flint_mpn_sqr(zd, xd, xn);
            top = zd[zn - 1];
        }
        else
        {
            top = flint_mpn_mul(zd, xd, xn, yd, yn);
        }

        zn -= (top == 0);
        z->_mp_size = (sgn >= 0) ? zn : -zn;
        return;
    }

    TMP_START;
    td = TMP_ALLOC(tn * sizeof(mp_limb_t));

    if (x == y)
    {
        flint_mpn_sqr(td, xd, xn);
        top = td[tn - 1];
    }
    else
    {
        top = flint_mpn_mul(td, xd, xn, yd, yn);
    }

    tn -= (top == 0);
    alloc = FLINT_MAX(tn, zn) + 1;
    MPZ_FIT_SIZE(z, alloc);
    zd = z->_mp_d;

    if (sgn >= 0)
    {
        if (zn >= tn)
        {
            top = mpn_add(zd, zd, zn, td, tn);
            zn_new = zn;
        }
        else
        {
            top = mpn_add(zd, td, tn, zd, zn);
            zn_new = tn;
        }

        zd[zn_new] = top;
        zn_new += (top != 0);
    }
    else
    {
        if (zn > tn || (zn == tn && mpn_cmp(zd, td, zn) >= 0))
        {
            mpn_sub(zd, zd, zn, td, tn);
            zn_new = zn;
        }
        else
        {
            mpn_sub(zd, td, tn, zd, zn);
            zn_new = tn;
            zn_signed = -zn_signed;
        }

        while (zn_new >= 1 && zd[zn_new - 1] == 0)
            zn_new--;
    }

    z->_mp_size = (zn_signed >= 0) ? zn_new : -zn_new;

    TMP_END;
}

void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
{
    fmpz c1, c2;
    __mpz_struct * mf;
	
    c1 = *g;
	
	if (!COEFF_IS_MPZ(c1))  /* g is small */
	{
		if (c1 < WORD(0)) fmpz_submul_ui(f, h, -c1);
		else fmpz_addmul_ui(f, h, c1);
		return;
	} 

	c2 = *h;
   
	if (!COEFF_IS_MPZ(c2))  /* h is small */
	{
		if (c2 < WORD(0)) fmpz_submul_ui(f, g, -c2);
		else fmpz_addmul_ui(f, g, c2);
		return;
	}

	/* both g and h are large */
    mf = _fmpz_promote_val(f);
    _flint_mpz_addmul_large(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2), 0);
    _fmpz_demote_val(f);  /* cancellation may have occurred	*/
}
