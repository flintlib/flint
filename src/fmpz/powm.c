/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "gmpcompat.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "gr.h"

/* assumes e is large */
static void
_fmpz_powm(fmpz_t res, const fmpz_t x, const fmpz_t e, const fmpz_t m)
{
    if (!COEFF_IS_MPZ(*m))  /* m is small */
    {
        /* todo: make n_powmod2_fmpz_preinv faster than mpz_powm and use here */
        ulong c1, c2;
        mpz_t zx, zm;
        __mpz_struct * zres;

        c2 = *m;
        c1 = fmpz_fdiv_ui(x, c2);

        zx->_mp_d = &c1;
        zx->_mp_size = (c1 != 0);
        zx->_mp_alloc = 1;

        zm->_mp_d = &c2;
        zm->_mp_size = 1;
        zm->_mp_alloc = 1;

        zres = _fmpz_promote(res);
        mpz_powm(zres, zx, COEFF_TO_PTR(*e), zm);
        _fmpz_demote_val(res);
    }
    else if (fmpz_is_zero(x) || fmpz_is_one(x))
    {
        fmpz_set(res, x);
    }
#ifdef FLINT_HAVE_FFT_SMALL
    else if (fmpz_bits(m) >= 70000)
    {
        gr_ctx_t gctx;
        fmpz_t t;

        gr_ctx_init_fmpz_mod(gctx, m);
        fmpz_init(t);

        /* fmpz_mod input must be reduced */
        GR_MUST_SUCCEED(gr_set_fmpz(t, x, gctx));

        if (!COEFF_IS_MPZ(*x))
            GR_MUST_SUCCEED(gr_generic_pow_fmpz_binexp(res, t, e, gctx));
        else
            GR_MUST_SUCCEED(gr_generic_pow_fmpz_sliding(res, t, e, gctx));

        fmpz_clear(t);
        gr_ctx_clear(gctx);
    }
#endif
    else
    {
        if (!COEFF_IS_MPZ(*x))  /* x is small */
        {
            mpz_t zx;
            __mpz_struct * zres;
            ulong c1;

            c1 = FLINT_ABS(*x);

            zx->_mp_d = &c1;
            zx->_mp_size = (*x == 0) ? 0 : ((*x > 0) ? 1 : -1);
            zx->_mp_alloc = 1;

            zres = _fmpz_promote(res);
            mpz_powm(zres, zx, COEFF_TO_PTR(*e), COEFF_TO_PTR(*m));
            _fmpz_demote_val(res);
        }
        else  /* x is large */
        {
            __mpz_struct * zres = _fmpz_promote(res);
            mpz_powm(zres, COEFF_TO_PTR(*x), COEFF_TO_PTR(*e), COEFF_TO_PTR(*m));
            _fmpz_demote_val(res);
        }
    }
}

void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)
{
    if (fmpz_sgn(m) <= 0)
    {
        flint_throw(FLINT_ERROR, "Exception in fmpz_powm: "
                                                  "Modulus is less than 1.\n");
    }
    else if (!COEFF_IS_MPZ(*e))  /* e is small */
    {
        if (*e >= 0)
        {
            fmpz_powm_ui(f, g, *e, m);
        }
        else
        {
            fmpz_t g_inv;
            fmpz_init(g_inv);
            if (!fmpz_invmod(g_inv, g, m))
            {
                fmpz_clear(g_inv);
                flint_throw(FLINT_ERROR, "Exception in fmpz_powm: "
                                                 "Base is not invertible.\n");
            }
            else
            {
                fmpz_powm_ui(f, g_inv, -*e, m);
                fmpz_clear(g_inv);
            }
        }
    }
    else  /* e is large */
    {
        _fmpz_powm(f, g, e, m);
    }
}
