/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2021 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz-impl.h"

/* Will not get called with x or y small. */
void
_flint_mpz_addmul_large(mpz_ptr z, mpz_srcptr x, mpz_srcptr y, int negate)
{
    slong xn, yn, tn, zn, zn_signed, zn_new, x_sgn, y_sgn, sgn, alloc;
    nn_srcptr xd, yd;
    nn_ptr zd;
    nn_ptr td;
    ulong top;
    TMP_INIT;

    xn = x->_mp_size;
    yn = y->_mp_size;
    x_sgn = xn;
    y_sgn = yn;
    xn = FLINT_ABS(xn);
    yn = FLINT_ABS(yn);

    if (xn < yn)
    {
        FLINT_SWAP(mpz_srcptr, x, y);
        FLINT_SWAP(slong, xn, yn);
        FLINT_SWAP(slong, x_sgn, y_sgn);
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

    /* Currently unreachable because fmpz_addmul and fmpz_submul both
       handle z == 0 specially. */
#if 1
    FLINT_ASSERT(zn != 0);
#else
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
#endif

    TMP_START;
    td = TMP_ALLOC(tn * sizeof(ulong));

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
    zd = FLINT_MPZ_REALLOC(z, alloc);

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
    fmpz c1, c2, c3;
    mpz_ptr mf;

    c1 = *g;
	c2 = *h;
    c3 = *f;

    /* todo: are the zero checks worth it for small input? */

    if (c1 == 0 || c2 == 0)
        return;

    if (c3 == 0)
    {
        fmpz_mul(f, g, h);
        return;
    }

    if (!COEFF_IS_MPZ(c1))  /* g is small */
    {
        if (!COEFF_IS_MPZ(c2))  /* both are small */
        {
            ulong p1, p0;
            smul_ppmm(p1, p0, c1, c2);

            if (!COEFF_IS_MPZ(c3))
            {
                ulong F1 = FLINT_SIGN_EXT(c3);
                add_ssaaaa(p1, p0, p1, p0, F1, c3);

                fmpz_set_signed_uiui(f, p1, p0);
            }
            else
            {
                mpz_ptr pF = COEFF_TO_PTR(c3);
                flint_mpz_add_signed_uiui(pF, pF, p1, p0);
                _fmpz_demote_val(f);  /* cancellation may have occurred	*/
            }
        }
        else
        {
            fmpz_addmul_si(f, h, c1);
        }
    }
    else if (!COEFF_IS_MPZ(c2))  /* h is small */
	{
		fmpz_addmul_si(f, g, c2);
	}
    else
    {
        mf = _fmpz_promote_val(f);
        _flint_mpz_addmul_large(mf, COEFF_TO_PTR(c1), COEFF_TO_PTR(c2), 0);
        _fmpz_demote_val(f);  /* cancellation may have occurred	*/
    }
}

void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slong x)
{
    fmpz F, G;

    G = *g;
    if (x == 0 || G == 0)
        return;

    F = *f;
    if (F == 0)
    {
        fmpz_mul_si(f, g, x);
        return;
    }

    if (!COEFF_IS_MPZ(G))
    {
        ulong p1, p0;
        smul_ppmm(p1, p0, G, x);

        if (!COEFF_IS_MPZ(F))
        {
            ulong F1 = FLINT_SIGN_EXT(F);
            add_ssaaaa(p1, p0, p1, p0, F1, F);
            fmpz_set_signed_uiui(f, p1, p0);
        }
        else
        {
            mpz_ptr pF = COEFF_TO_PTR(F);
            flint_mpz_add_signed_uiui(pF, pF, p1, p0);
            _fmpz_demote_val(f);
        }
    }
    else
    {
        mpz_ptr pG = COEFF_TO_PTR(G);
        mpz_ptr pF = _fmpz_promote_val(f);

        if (x < 0)
            flint_mpz_submul_ui(pF, pG, -x);
        else
            flint_mpz_addmul_ui(pF, pG, x);

        _fmpz_demote_val(f);
    }
}

void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong x)
{
    fmpz F, G;

    G = *g;
    if (x == 0 || G == 0)
        return;

    F = *f;
    if (F == 0)
    {
        fmpz_mul_ui(f, g, x);
        return;
    }

    if (!COEFF_IS_MPZ(G))
    {
        ulong p1, p0;

        if (x <= WORD_MAX)
        {
            smul_ppmm(p1, p0, G, x);
        }
        else
        {
            umul_ppmm(p1, p0, FLINT_ABS(G), x);
            if (G < 0)
            {
                p1 = -p1 - (p0 != 0);
                p0 = -p0;
            }
        }

        if (!COEFF_IS_MPZ(F))
        {
            ulong F1 = FLINT_SIGN_EXT(F);
            add_ssaaaa(p1, p0, p1, p0, F1, F);
            fmpz_set_signed_uiui(f, p1, p0);
        }
        else
        {
            mpz_ptr pF = COEFF_TO_PTR(F);
            flint_mpz_add_signed_uiui(pF, pF, p1, p0);
            _fmpz_demote_val(f);
        }
    }
    else
    {
        mpz_ptr pG = COEFF_TO_PTR(G);
        mpz_ptr pF = _fmpz_promote_val(f);

        flint_mpz_addmul_ui(pF, pG, x);
        _fmpz_demote_val(f);
    }
}
