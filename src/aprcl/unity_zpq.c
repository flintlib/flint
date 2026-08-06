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
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "aprcl.h"

/* -------------------------------------------------------------------------- */
/* from unity_zpq_init.c */

void
unity_zpq_init(unity_zpq f, ulong q, ulong p, const fmpz_t n)
{
    slong i;

    f->p = p;
    f->q = q;

    fmpz_mod_ctx_init(f->ctx, n);
    f->polys = (fmpz_mod_poly_t *) flint_malloc(p * sizeof(fmpz_mod_poly_t));

    for (i = 0; i < p; i++)
    {
        fmpz_mod_poly_init(f->polys[i], f->ctx);
    }
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_clear.c */

void
unity_zpq_clear(unity_zpq f)
{
    slong i;

    for (i = 0; i < f->p; i++)
    {
        fmpz_mod_poly_clear(f->polys[i], f->ctx);
    }

    f->p = 0;
    f->q = 0;

    fmpz_mod_ctx_clear(f->ctx);
    flint_free(f->polys);
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_copy.c */

void
unity_zpq_copy(unity_zpq f, const unity_zpq g)
{
    slong i;

    for (i = 0; i < f->p; i++)
    {
        fmpz_mod_poly_set(f->polys[i], g->polys[i], g->ctx);
    }
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_swap.c */

void
unity_zpq_swap(unity_zpq f, unity_zpq g)
{
    fmpz_mod_poly_t *temp = f->polys;
    f->polys = g->polys;
    g->polys = temp;
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_equal.c */

int
unity_zpq_equal(const unity_zpq f, const unity_zpq g)
{
    slong i;

    if (f->p != g->p)
        return 0;

    if (f->q != g->q)
        return 0;

    if (!fmpz_equal(fmpz_mod_ctx_modulus(f->ctx), fmpz_mod_ctx_modulus(g->ctx)))
        return 0;

    for (i = 0; i < f->p; i++)
        if (!fmpz_mod_poly_equal(f->polys[i], g->polys[i], g->ctx))
            return 0;

    return 1;
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_coeff_set.c */

void
unity_zpq_coeff_set_fmpz(unity_zpq f, slong i, slong j, const fmpz_t x)
{
    fmpz_mod_poly_set_coeff_fmpz(f->polys[j], i, x, f->ctx);
}

void
unity_zpq_coeff_set_ui(unity_zpq f, slong i, slong j, ulong x)
{
    fmpz_mod_poly_set_coeff_ui(f->polys[j], i, x, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_coeff_add.c */

void
unity_zpq_coeff_add(unity_zpq f, slong i, slong j, const fmpz_t x)
{
    if (i >= fmpz_mod_poly_length(f->polys[j], f->ctx))
    {
        fmpz_mod_poly_set_coeff_fmpz(f->polys[j], i, x, f->ctx);
        return;
    }

    fmpz_add(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i, x);
    if (fmpz_cmp(f->polys[j]->coeffs + i, fmpz_mod_ctx_modulus(f->ctx)) >= 0)
        fmpz_sub(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i,
                                                 fmpz_mod_ctx_modulus(f->ctx));
}

void
unity_zpq_coeff_add_ui(unity_zpq f, slong i, slong j, ulong x)
{
    if (i >= fmpz_mod_poly_length(f->polys[j], f->ctx))
    {
        fmpz_mod_poly_set_coeff_ui(f->polys[j], i, x, f->ctx);
        return;
    }

    fmpz_add_ui(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i, x);
    if (fmpz_cmp(f->polys[j]->coeffs + i, fmpz_mod_ctx_modulus(f->ctx)) >= 0)
        fmpz_sub(f->polys[j]->coeffs + i, f->polys[j]->coeffs + i,
                                                 fmpz_mod_ctx_modulus(f->ctx));
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_add.c */

#if FLINT_WANT_ASSERT
# include "fmpz.h"
# include "fmpz_mod.h"
#endif

void
unity_zpq_add(unity_zpq f, const unity_zpq g, const unity_zpq h)
{
    ulong i;

    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(h->ctx)));

    for (i = 0; i < f->p; i++)
    {
        fmpz_mod_poly_add(f->polys[i], g->polys[i], h->polys[i], f->ctx);
    }
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_mul.c */

void
unity_zpq_mul(unity_zpq f, const unity_zpq g, const unity_zpq h)
{
    slong i, j, k;
    ulong p, q;
    fmpz_mod_poly_t temp;

    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(h->ctx)));

    q = f->q;
    p = f->p;
    fmpz_mod_poly_init(temp, f->ctx);

    for (i = 0; i < p; i++)
    {
        fmpz_mod_poly_zero(f->polys[i], f->ctx);
    }

    for (i = 0; i < p; i++)
    {
        for (j = 0; j < p; j++)
        {
            ulong qpow;

            qpow = n_addmod(i, j, p);
            fmpz_mod_poly_mul(temp, g->polys[i], h->polys[j], f->ctx);

            if (temp->length == 0)
                continue;

            for (k = temp->length - 1; k >= q; k--)
            {
                fmpz_add(temp->coeffs + k - q,
                        temp->coeffs + k - q, temp->coeffs + k);
                fmpz_set_ui(temp->coeffs + k, 0);
                fmpz_mod(temp->coeffs + k - q, temp->coeffs + k - q,
                                                 fmpz_mod_ctx_modulus(f->ctx));
            }
            _fmpz_mod_poly_normalise(temp);

            fmpz_mod_poly_add(f->polys[qpow], f->polys[qpow], temp, f->ctx);
        }
    }

    fmpz_mod_poly_clear(temp, f->ctx);
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_mul_unity_p.c */

void
_unity_zpq_mul_unity_p(unity_zpq f)
{
    slong i;

    for (i = f->p - 1; i > 0; i--)
        fmpz_mod_poly_swap(f->polys[i], f->polys[i - 1], f->ctx);
}

/*
    Computes unity_zpq * \zeta_p by swapping poly coeffs.
*/
void
unity_zpq_mul_unity_p_pow(unity_zpq f, const unity_zpq g, slong k)
{
    slong i;

    unity_zpq_copy(f, g);

    for (i = 0; i < k; i++)
        _unity_zpq_mul_unity_p(f);
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_pow.c */

void
unity_zpq_pow(unity_zpq f, const unity_zpq g, const fmpz_t pow)
{
    unity_zpq value;
    fmpz_t power, rem;

    unity_zpq_init(value, f->q, f->p, fmpz_mod_ctx_modulus(f->ctx));
    fmpz_init_set(power, pow);
    fmpz_init(rem);

    unity_zpq_coeff_set_ui(f, 0, 0, 1);

    unity_zpq_copy(value, g);

    while (fmpz_is_zero(power) == 0)
    {
        unity_zpq temp_pow;

        fmpz_fdiv_r_2exp(rem, power, 1);
        if (fmpz_is_zero(rem) == 0)
        {
            unity_zpq temp;
            unity_zpq_init(temp, f->q, f->p, fmpz_mod_ctx_modulus(f->ctx));

            unity_zpq_mul(temp, f, value);
            unity_zpq_swap(f, temp);

            unity_zpq_clear(temp);
        }

        unity_zpq_init(temp_pow, f->q, f->p, fmpz_mod_ctx_modulus(f->ctx));
        unity_zpq_mul(temp_pow, value, value);
        unity_zpq_swap(value, temp_pow);
        fmpz_fdiv_q_2exp(power, power, 1);

        unity_zpq_clear(temp_pow);
    }


    fmpz_clear(power);
    fmpz_clear(rem);
    unity_zpq_clear(value);
}

void
unity_zpq_pow_ui(unity_zpq f, const unity_zpq g, ulong pow)
{
    fmpz_t p;

    fmpz_init_set_ui(p, pow);

    unity_zpq_pow(f, g, p);

    fmpz_clear(p);
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_gauss_sum.c */

/*
    Computes gauss sum for character \chi corresponding (q, p).
*/
void
unity_zpq_gauss_sum(unity_zpq f, ulong q, ulong p)
{
    slong i, qinv, qpow, ppow, g;

    g = n_primitive_root_prime(q);
    qinv = n_preinvert_limb(q);
    qpow = 1;
    ppow = 0;

    for (i = 1; i < q; i++)
    {
        qpow = n_mulmod2_preinv(qpow, g, q, qinv);
        ppow = n_addmod(ppow, 1, p);
        unity_zpq_coeff_add_ui(f, qpow, ppow, 1);
    }
}

/* -------------------------------------------------------------------------- */
/* from unity_zpq_gauss_sum_character_pow.c */

/*
    Computes gauss sum for character \chi^n corresponding (q, p).
*/
void
unity_zpq_gauss_sum_character_pow(unity_zpq f, ulong q, ulong p, ulong pow)
{
    slong i, qinv, pinv, qpow, ppow, g;

    g = n_primitive_root_prime(q);
    qinv = n_preinvert_limb(q);
    pinv = n_preinvert_limb(p);
    qpow = 1;

    for (i = 1; i < q; i++)
    {
        qpow = n_mulmod2_preinv(qpow, g, q, qinv);
        ppow = n_mulmod2_preinv(i, pow, p, pinv);
        unity_zpq_coeff_add_ui(f, qpow, ppow, 1);
    }
}

/*
    Computes gauss sum for character \chi^n corresponding (q, p).
*/
void
unity_zpq_gauss_sum_sigma_pow(unity_zpq f, ulong q, ulong p)
{
    ulong n;

    n = fmpz_fdiv_ui(fmpz_mod_ctx_modulus(f->ctx), p);
    unity_zpq_gauss_sum_character_pow(f, q, p, n);
}
