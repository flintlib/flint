/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr.h"
#include "gr_generic.h"

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
#if FLINT_HAVE_FFT_SMALL
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

/* assumes m is large */
static void
_fmpz_powm_ui(fmpz_t res, const fmpz_t x, ulong e, const fmpz_t m)
{
    if (e == 1)
    {
        fmpz_mod(res, x, m);
    }
    /* don't bother handling aliasing for the extremely unusual case where res == m */
    else if (e == 2 && res != m)
    {
        fmpz_mul(res, x, x);
        fmpz_mod(res, res, m);
    }
    else if (e == 3 && res != m)
    {
        if (res == x)
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_mul(t, x, x);
            fmpz_mod(t, t, m);
            fmpz_mul(t, t, x);
            fmpz_mod(res, t, m);
            fmpz_clear(t);
        }
        else
        {
            fmpz_mul(res, x, x);
            fmpz_mod(res, res, m);
            fmpz_mul(res, res, x);
            fmpz_mod(res, res, m);
        }
    }
    else if (e == 4 && res != m)
    {
        fmpz_mul(res, x, x);
        fmpz_mod(res, res, m);
        fmpz_mul(res, res, res);
        fmpz_mod(res, res, m);
    }
    else if (fmpz_is_zero(x) || fmpz_is_one(x))
    {
        fmpz_set(res, x);
    }
#if FLINT_HAVE_FFT_SMALL
    else if (fmpz_bits(m) >= 70000)
    {
        gr_ctx_t gctx;
        fmpz_t t;

        gr_ctx_init_fmpz_mod(gctx, m);
        fmpz_init(t);

        /* fmpz_mod input must be reduced */
        GR_MUST_SUCCEED(gr_set_fmpz(t, x, gctx));

        if (!COEFF_IS_MPZ(*x) || FLINT_BIT_COUNT(e) < 20)
            GR_MUST_SUCCEED(gr_generic_pow_ui_binexp(res, t, e, gctx));
        else
            GR_MUST_SUCCEED(gr_generic_pow_ui_sliding(res, t, e, gctx));

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
            flint_mpz_powm_ui(zres, zx, e, COEFF_TO_PTR(*m));
            _fmpz_demote_val(res);
        }
        else  /* x is large */
        {
            __mpz_struct * zres = _fmpz_promote(res);
            flint_mpz_powm_ui(zres, COEFF_TO_PTR(*x), e, COEFF_TO_PTR(*m));
            _fmpz_demote_val(res);
        }
    }
}

void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong e, const fmpz_t m)
{
    if (fmpz_sgn(m) <= 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_powm_ui). Modulus is less than 1.\n");
    }

    if (fmpz_is_one(m))
    {
        fmpz_zero(f);
    }
    else if (e == 0)
    {
        fmpz_one(f);
    }
    else  /* e != 0, m > 0 */
    {
        fmpz g2 = *g;
        fmpz m2 = *m;

        if (!COEFF_IS_MPZ(m2))  /* m is small */
        {
            if (!COEFF_IS_MPZ(g2))  /* g is small */
            {
                mp_limb_t minv = n_preinvert_limb(m2);

                _fmpz_demote(f);

                if (g2 >= 0)
                {
                    g2 = n_mod2_preinv(g2, m2, minv);
                    *f = n_powmod2_ui_preinv(g2, e, m2, minv);
                }
                else
                {
                    g2 = n_mod2_preinv(-g2, m2, minv);
                    *f = n_powmod2_ui_preinv(g2, e, m2, minv);
                    if ((e & UWORD(1)))
                        *f = n_negmod(*f, m2);
                }
            }
            else  /* g is large */
            {
                ulong c3;
                mpz_t m3;
                __mpz_struct *ptr = _fmpz_promote(f);

                c3 = m2;
                m3->_mp_d = &c3;
                m3->_mp_size = 1;
                m3->_mp_alloc = 1;

                flint_mpz_powm_ui(ptr, COEFF_TO_PTR(g2), e, m3);
                _fmpz_demote_val(f);
            }
        }
        else  /* m is large */
        {
            _fmpz_powm_ui(f, g, e, m);
        }
    }
}
