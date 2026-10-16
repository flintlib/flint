/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "aprcl.h"

/* -------------------------------------------------------------------------- */
/* from unity_zp_init.c */

void
unity_zp_init(unity_zp f, ulong p, ulong exp, const fmpz_t n)
{
    f->p = p;
    f->exp = exp;
    fmpz_mod_ctx_init(f->ctx, n);
    fmpz_mod_poly_init(f->poly, f->ctx);
}

void
unity_zp_clear(unity_zp f)
{
    fmpz_mod_poly_clear(f->poly, f->ctx);
    fmpz_mod_ctx_clear(f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_copy.c */

#if FLINT_WANT_ASSERT
# include "fmpz.h"
# include "fmpz_mod.h"
#endif

void
unity_zp_copy(unity_zp f, const unity_zp g)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_set(f->poly, g->poly, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_swap.c */

#if FLINT_WANT_ASSERT
# include "fmpz.h"
# include "fmpz_mod.h"
#endif

void
unity_zp_swap(unity_zp f, unity_zp g)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_swap(f->poly, g->poly, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_set_zero.c */

void
unity_zp_set_zero(unity_zp f)
{
    fmpz_mod_poly_zero(f->poly, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_coeff.c */

void
unity_zp_coeff_add_fmpz(unity_zp f, ulong ind, const fmpz_t x)
{
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mod_poly_get_coeff_fmpz(coeff, f->poly, ind, f->ctx);
    if (fmpz_is_zero(coeff))
    {
        fmpz_clear(coeff);
        fmpz_mod_poly_set_coeff_fmpz(f->poly, ind, x, f->ctx);
        return;
    }
    fmpz_clear(coeff);
    fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind, x);
    if (fmpz_cmp(f->poly->coeffs + ind, fmpz_mod_ctx_modulus(f->ctx)) >= 0)
        fmpz_sub(f->poly->coeffs + ind, f->poly->coeffs + ind,
                                                 fmpz_mod_ctx_modulus(f->ctx));
}

void
unity_zp_coeff_add_ui(unity_zp f, ulong ind, ulong x)
{
    fmpz_t coeff;
    fmpz_init(coeff);
    fmpz_mod_poly_get_coeff_fmpz(coeff, f->poly, ind, f->ctx);
    if (fmpz_is_zero(coeff))
    {
        fmpz_clear(coeff);
        fmpz_mod_poly_set_coeff_ui(f->poly, ind, x, f->ctx);
        return;
    }
    fmpz_clear(coeff);
    fmpz_add_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, x);
    if (fmpz_cmp(f->poly->coeffs + ind, fmpz_mod_ctx_modulus(f->ctx)) >= 0)
        fmpz_sub(f->poly->coeffs + ind, f->poly->coeffs + ind,
                                                 fmpz_mod_ctx_modulus(f->ctx));
}

void
unity_zp_coeff_inc(unity_zp f, ulong ind)
{
    if (ind >= f->poly->length)
    {
        fmpz_mod_poly_set_coeff_ui(f->poly, ind, 1, f->ctx);
        return;
    }

    fmpz_add_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, 1);
    if (fmpz_equal(f->poly->coeffs + ind, fmpz_mod_ctx_modulus(f->ctx)))
        fmpz_set_ui(f->poly->coeffs + ind, 0);
}

void
unity_zp_coeff_dec(unity_zp f, ulong ind)
{
    if (ind >= f->poly->length)
    {
        fmpz_mod_poly_set_coeff_si(f->poly, ind, -1, f->ctx);
        return;
    }

    fmpz_sub_ui(f->poly->coeffs + ind, f->poly->coeffs + ind, 1);
    if (fmpz_equal_si(f->poly->coeffs + ind, -1))
        fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind,
                                                 fmpz_mod_ctx_modulus(f->ctx));
}

void
unity_zp_coeff_set_fmpz(unity_zp f, ulong ind, const fmpz_t x)
{
    fmpz_mod_poly_set_coeff_fmpz(f->poly, ind, x, f->ctx);
}

void
unity_zp_coeff_set_ui(unity_zp f, ulong ind, ulong x)
{
    fmpz_mod_poly_set_coeff_ui(f->poly, ind, x, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_add.c */

#if FLINT_WANT_ASSERT
# include "fmpz.h"
# include "fmpz_mod.h"
#endif

void
unity_zp_add(unity_zp f, const unity_zp g, const unity_zp h)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(h->ctx)));
    fmpz_mod_poly_add(f->poly, g->poly, h->poly, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_equal.c */

int
unity_zp_equal(unity_zp f, unity_zp g)
{
    /*
        f and g can be reduced only by modylo x^{p^k} - 1,
        so reduce by cyclotomic polynomial
    */
    _unity_zp_reduce_cyclotomic(f);
    _unity_zp_reduce_cyclotomic(g);

    return fmpz_mod_poly_equal(f->poly, g->poly, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_is_unity.c */

slong
unity_zp_is_unity(unity_zp f)
{
    ulong result;
    ulong i, p_pow;
    unity_zp unity;

    p_pow = n_pow(f->p, f->exp);
    unity_zp_init(unity, f->p, f->exp, fmpz_mod_ctx_modulus(f->ctx));

    /* if the power was not found returns -1 */
    result = -1;
    for (i = 0; i < p_pow; i++)
    {
        /* set unity = \zeta_{p^k}^i */
        unity_zp_set_zero(unity);
        unity_zp_coeff_set_ui(unity, i, 1);

        /* check if f = zeta_{p^k}^i */
        if (unity_zp_equal(unity, f) == 1)
        {
            /* if so, returns \zeta_{p^k} power */
            result = i;
            break;
        }
    }

    unity_zp_clear(unity);
    return result;
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_reduce_cyclotomic.c */

void
_unity_zp_reduce_cyclotomic_divmod(unity_zp f)
{
    ulong i, j, ppow1, ppow2, cycl_pow;

    ppow2 = n_pow(f->p, f->exp - 1);
    ppow1 = ppow2 * f->p;
    cycl_pow = (f->p - 1) * ppow2;

    for (i = f->poly->length - 1; i >= ppow1; i--)
    {
        fmpz_add(f->poly->coeffs + i - ppow1,
                f->poly->coeffs + i - ppow1, f->poly->coeffs + i);

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    for (i = f->poly->length - 1; i >= cycl_pow; i--)
    {
        if (fmpz_is_zero(f->poly->coeffs + i))
            continue;

        for (j = 0; j < f->p - 1; j++)
        {
            ulong ind = i - cycl_pow + j * ppow2;
            fmpz_sub(f->poly->coeffs + ind,
                    f->poly->coeffs + ind, f->poly->coeffs + i);
        }

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    _fmpz_mod_poly_normalise(f->poly);
    _fmpz_vec_scalar_mod_fmpz(f->poly->coeffs,
               f->poly->coeffs, f->poly->length, fmpz_mod_ctx_modulus(f->ctx));
    _fmpz_mod_poly_normalise(f->poly);
}

void
_unity_zp_reduce_cyclotomic(unity_zp f)
{
    ulong i, j, ppow, cycl_pow;

    if (f->poly->length == 0)
        return;

    ppow = n_pow(f->p, f->exp - 1);
    cycl_pow = (f->p - 1) * ppow;

    for (i = f->poly->length - 1; i >= cycl_pow; i--)
    {
        if (fmpz_is_zero(f->poly->coeffs + i))
            continue;

        for (j = 0; j < f->p - 1; j++)
        {
            ulong ind = i - cycl_pow + j * ppow;
            fmpz_sub(f->poly->coeffs + ind,
                    f->poly->coeffs + ind, f->poly->coeffs + i);

            if (fmpz_cmp_ui(f->poly->coeffs + ind, 0) < 0)
                fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind,
                                                 fmpz_mod_ctx_modulus(f->ctx));
        }

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    _fmpz_mod_poly_normalise(f->poly);
}

void
unity_zp_reduce_cyclotomic(unity_zp f, const unity_zp g)
{
    unity_zp_copy(f, g);
    _unity_zp_reduce_cyclotomic(f);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_aut.c */

/*
    Computes f such that \sigma_x(g) = f.
*/
void
unity_zp_aut(unity_zp f, const unity_zp g, ulong x)
{
    ulong i, p_pow, p_pow_preinv;
    fmpz_t coeff;
    fmpz_init(coeff);

    p_pow = n_pow(f->p, f->exp);
    p_pow_preinv = n_preinvert_limb(p_pow);

    unity_zp_set_zero(f);

    /* for i = 0, 1,..., p^k set f[i * x mod p^k] = g[i] */
    for (i = 0; i < p_pow; i++)
    {
        /* compute x * i mod p^k */
        ulong ind = n_mulmod2_preinv(x, i, p_pow, p_pow_preinv);

        /* set f[ind] = g[i] */
        fmpz_mod_poly_get_coeff_fmpz(coeff, g->poly, i, g->ctx);
        unity_zp_coeff_add_fmpz(f, ind, coeff);
    }

    _unity_zp_reduce_cyclotomic(f);
    fmpz_clear(coeff);
}

/* -------------------------------------------------------------------------- */
/* from unity_zp_aut_inv.c */

/*
    Computes f such that \sigma_x(f) = g.
*/
void
unity_zp_aut_inv(unity_zp f, const unity_zp g, ulong x)
{
    ulong i, j, p_pow1, p_pow2, m, p_pow_preinv;
    fmpz_t f_coeff, g_coeff;

    fmpz_init(f_coeff);
    fmpz_init(g_coeff);
    p_pow1 = n_pow(f->p, f->exp - 1);   /* p_pow1 = p^{k - 1}       */
    p_pow2 = p_pow1 * f->p;             /* p_pow2 = p^k             */
    m = (f->p - 1) * p_pow1;            /* m = (p - 1) * p^{k - 1}  */
    p_pow_preinv = n_preinvert_limb(p_pow2);
    unity_zp_set_zero(f);

    /* for i = 0, 1,..., m - 1 set f[i] = g[xi mod p^k] */
    for (i = 0; i < m; i++)
    {
        /* set g_ind = x * i mod p^k */
        ulong g_ind = n_mulmod2_preinv(x, i, p_pow2, p_pow_preinv);

        /* set g_coeff to g[g_ind] */
        fmpz_mod_poly_get_coeff_fmpz(g_coeff, g->poly, g_ind, g->ctx);

        /* set f[i] = g[x * i mod p^k] */
        unity_zp_coeff_set_fmpz(f, i, g_coeff);
    }

    /*
        for i = m, m + 1,..., p^k - 1
        for j = 1, 2,..., p - 1
        set f[i - j * p^{k - 1}] =
        (f[i - j * p^{k - 1}] - g[x * i mod p^k]) mod n
    */
    for (i = m; i < p_pow2; i++)
    {
        /* set g_ind = x * i mod p^k */
        ulong g_ind = n_mulmod2_preinv(x, i, p_pow2, p_pow_preinv);

        for (j = 1; j < f->p; j++)
        {
            /* set f_ind = i - j * p^{k - 1} */
            ulong f_ind = i - j * p_pow1;

            /* set g_coeff = g[x * i mod p^k] */
            fmpz_mod_poly_get_coeff_fmpz(g_coeff, g->poly, g_ind, g->ctx);

            /* set f_coeff = f[i - j * p^{k - 1}] */
            fmpz_mod_poly_get_coeff_fmpz(f_coeff, f->poly, f_ind, f->ctx);

            /* set f_coeff = f[i - j * p^{k - 1}] - g[x * i mod p^k] */
            fmpz_sub(f_coeff, f_coeff, g_coeff);

            /* set f[i - j * p^{k - 1}] = f_coeff */
            unity_zp_coeff_set_fmpz(f, f_ind, f_coeff);
        }
    }

    fmpz_clear(f_coeff);
    fmpz_clear(g_coeff);
}
