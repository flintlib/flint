/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include "aprcl.h"

void
unity_zp_mont_ninv(fmpz_t ninv, const unity_zp_mont f)
{
    fmpz_t d, n, r, rinv;

    fmpz_init(d);
    fmpz_init(n);
    fmpz_init(r);
    fmpz_init(rinv);

    fmpz_set(n, f->n);
    fmpz_set_ui(d, 1);
    fmpz_setbit(r, fmpz_bits(n));
    fmpz_neg(n, n);


    fmpz_xgcd(d, rinv, ninv, r, n);
    fmpz_mod(ninv, ninv, r);

    fmpz_clear(d);
    fmpz_clear(n);
    fmpz_clear(r);
    fmpz_clear(rinv);
}

/* Init Montgomery form structure */
void
unity_zp_mont_init(unity_zp_mont f, ulong p, ulong exp, const fmpz_t n, const fmpz_t ninv)
{
    f->p = p;
    f->exp = exp;
    fmpz_init_set(f->n, n);
    fmpz_init(f->nr);

    fmpz_init_set(f->ninv, ninv);
    f->r = fmpz_bits(n);

    fmpz_poly_init(f->poly);
}

void
unity_zp_mont_clear(unity_zp_mont f)
{
    fmpz_clear(f->nr);
    fmpz_clear(f->n);
    fmpz_clear(f->ninv);
    fmpz_poly_clear(f->poly);
}

/* 
    Reduce Montgomery form poly by cyclotomic polynomial; 
    same as unity_zp_reduce_cyclotomic.
 */
void _unity_zp_mont_reduce_cyclotomic(unity_zp_mont f)
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

            while (fmpz_cmp_ui(f->poly->coeffs + ind, 0) < 0)
                fmpz_add(f->poly->coeffs + ind, f->poly->coeffs + ind, f->nr);
        }

        fmpz_set_ui(f->poly->coeffs + i, 0);
    }

    _fmpz_poly_normalise(f->poly);
}

void
unity_zp_mont_set_zero(unity_zp_mont f)
{
    fmpz_poly_zero(f->poly);
}

void
unity_zp_mont_sqr(unity_zp_mont f, const unity_zp_mont g)
{
    ulong i, p;

    fmpz_poly_sqr(f->poly, g->poly);

    if (f->poly->length == 0)
        return;

    p = n_pow(f->p, f->exp);
    for (i = f->poly->length - 1; i >= p; i--)
    {
        fmpz_add(f->poly->coeffs + i - p,
                f->poly->coeffs + i - p, f->poly->coeffs + i);

        fmpz_set_ui(f->poly->coeffs + i, 0);
        while (fmpz_cmp(f->poly->coeffs + i - p, f->nr) >= 0)
            fmpz_sub(f->poly->coeffs + i - p, f->poly->coeffs + i - p, f->nr);
    }

    _unity_zp_mont_reduce_cyclotomic(f);
    unity_zp_mont_reduction(f);
}

void
unity_zp_mont_mul(unity_zp_mont f, const unity_zp_mont g, const unity_zp_mont h)
{
    ulong i, p;

    fmpz_poly_mul(f->poly, g->poly, h->poly);

    if (f->poly->length == 0)
        return;

    p = n_pow(f->p, f->exp);
    for (i = f->poly->length - 1; i >= p; i--)
    {
        fmpz_add(f->poly->coeffs + i - p,
                f->poly->coeffs + i - p, f->poly->coeffs + i);

        fmpz_set_ui(f->poly->coeffs + i, 0);
        if (fmpz_cmp(f->poly->coeffs + i - p, f->nr) >= 0)
            fmpz_sub(f->poly->coeffs + i - p, f->poly->coeffs + i - p, f->nr);
    }

    _unity_zp_mont_reduce_cyclotomic(f);
    unity_zp_mont_reduction(f);
}

void
unity_zp_mont_swap(unity_zp_mont f, unity_zp_mont g)
{
    fmpz_poly_swap(f->poly, g->poly);
}

void
unity_zp_mont_copy(unity_zp_mont f, const unity_zp_mont g)
{
    fmpz_poly_set(f->poly, g->poly);
}

/* Convert unity_zp into Montgomery form */
void
unity_zp_to_mont(unity_zp_mont f, const unity_zp g, const fmpz_t r)
{
    fmpz_mod_poly_get_fmpz_poly(f->poly, g->poly);
    _fmpz_vec_scalar_mul_fmpz(f->poly->coeffs,
            f->poly->coeffs, f->poly->length, r);
    _fmpz_vec_scalar_mod_fmpz(f->poly->coeffs,
            f->poly->coeffs, f->poly->length, f->n);
}

/* Convert from Montgomery form into unity_zp */
void
unity_zp_from_mont(unity_zp f, unity_zp_mont g)
{
    unity_zp_mont_reduction(g);
    fmpz_mod_poly_set_fmpz_poly(f->poly, g->poly);
}

/* Montgomery REDC for all elements of f->poly */
void
unity_zp_mont_reduction(unity_zp_mont f)
{
    slong i;
    fmpz_t m, t;
    fmpz_init(m);
    fmpz_init(t);

    for (i = 0; i < f->poly->length; i++)
    {
        /* 
            set f[i] = f[i] * r^(-1);
            f[i] = (f[i] + (f[i] * ninv mod r) * n) / r
        */
        fmpz_fdiv_r_2exp(m, f->poly->coeffs + i, f->r);
        fmpz_mul(t, m, f->ninv);
        fmpz_fdiv_r_2exp(m, t, f->r);
        fmpz_mul(t, f->n, m);
        fmpz_add(f->poly->coeffs + i, f->poly->coeffs + i, t);
        fmpz_fdiv_q_2exp(f->poly->coeffs + i, f->poly->coeffs + i, f->r);
        if (fmpz_cmp(f->n, f->poly->coeffs + i) <= 0)
            fmpz_sub(f->poly->coeffs + i, f->poly->coeffs + i, f->n);
    }

    fmpz_clear(m);
    fmpz_clear(t);
}

