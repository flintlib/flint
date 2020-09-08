/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <math.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "fmpz_mod_poly.h"

void
fmpz_mod_poly_factor_distinct_deg(fmpz_mod_poly_factor_t res,
                             const fmpz_mod_poly_t poly, slong * const *degs,
                                                      const fmpz_mod_ctx_t ctx)
{
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t f, g, v, vinv, tmp;
    fmpz_mod_poly_t *h, *H, *I;
    slong i, j, l, m, n, index, d;
    fmpz_mat_t HH, HHH;
    double beta;

    fmpz_mod_poly_init(v, ctx);

    fmpz_mod_poly_make_monic(v, poly, ctx);
    
    n = fmpz_mod_poly_degree(poly, ctx);
    if (n == 1)
    {
        fmpz_mod_poly_factor_insert(res, v, 1, ctx);
        (*degs)[0] = 1;
        fmpz_mod_poly_clear(v, ctx);
        return;
    }
    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    /* initialization */
    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_init(g, ctx);
    fmpz_mod_poly_init(vinv, ctx);
    fmpz_mod_poly_init(tmp, ctx);

    if (!(h = flint_malloc((2 * m + l + 1) * sizeof(fmpz_mod_poly_struct))))
    {
        flint_printf("Exception (fmpz_mod_poly_factor_distinct_deg):\n");
        flint_printf("Not enough memory.\n");
        flint_abort();
    }
    H = h + (l + 1);
    I = H + m;

    for (i = 0; i < 2*m + l + 1; i++)
        fmpz_mod_poly_init(h[i], ctx);

    fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
    fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length, ctx);

    /* compute baby steps: h[i]=x^{p^i}mod v */
    fmpz_mod_poly_set_coeff_ui(h[0], 1, 1, ctx);
    fmpz_mod_poly_powmod_x_fmpz_preinv(h[1], p, v, vinv, ctx);
    if (fmpz_sizeinbase(p, 2) > ((n_sqrt(v->length - 1) + 1) * 3) / 4)
    {
        for (i= 1; i < FLINT_BIT_COUNT (l); i++)
            fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(*(h + 1 +
                                                             (1 << (i - 1))),
                                                             *(h + 1),
                                                             (1 << (i - 1)),
                                                             (1 << (i - 1)),
                                                             *(h + (1 << (i - 1))),
                                                             v, vinv, ctx);
        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv(*(h + 1 +
                                                         (1 << (i - 1))),
                                                         *(h + 1),
                                                         (1 << (i - 1)),
                                                         l - (1 << (i - 1)),
                                                         *(h + (1 << (i - 1))),
                                                         v, vinv, ctx);
    }
    else
    {
        for (i = 2; i < l + 1; i++)
        {
            fmpz_mod_poly_init(h[i], ctx);
            fmpz_mod_poly_powmod_fmpz_binexp_preinv(h[i], h[i - 1], p,
                                              v, vinv, ctx);
        }
    }

    /* compute coarse distinct-degree factorisation */
    index= 0;
    fmpz_mod_poly_set(H[0], h[l], ctx);
    fmpz_mat_init(HH, n_sqrt(v->length - 1) + 1, v->length - 1);
    fmpz_mod_poly_precompute_matrix(HH, H[0], v, vinv, ctx);
    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[i]=x^{p^(li)}mod v */
        if (j > 0)
        {
            if (I[j - 1]->length > 1)
            {
                _fmpz_mod_poly_reduce_matrix_mod_poly(HHH, HH, v, ctx);
                fmpz_mat_clear(HH);
                fmpz_mat_init_set(HH, HHH);
                fmpz_mat_clear(HHH);
                fmpz_mod_poly_rem(tmp, H[j - 1], v, ctx);
                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H[j], tmp,
                                                             HH, v, vinv, ctx);
            }
            else
                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H[j],
                                                   H[j - 1], HH, v, vinv, ctx);
        }
        /* compute interval polynomials */
        fmpz_mod_poly_set_coeff_ui(I[j], 0, 1, ctx);
        for (i = l - 1; (i >= 0) && (2 * d <= v->length - 1); i--, d++)
        {
            fmpz_mod_poly_rem(tmp, h[i], v, ctx);
            fmpz_mod_poly_sub(tmp, H[j], tmp, ctx);
            fmpz_mod_poly_mulmod_preinv(I[j], tmp, I[j], v, vinv, ctx);
        }

        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        fmpz_mod_poly_gcd(I[j], v, I[j], ctx);
        if (I[j]->length > 1)
        {
            fmpz_mod_poly_remove(v, I[j], ctx);
            fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
            fmpz_mod_poly_inv_series_newton(vinv, vinv, v->length, ctx);
        }
        if (v->length - 1 < 2 * d)
        {
            break;
        }
    }
    if (v->length > 1)
    {
        fmpz_mod_poly_factor_insert(res, v, 1, ctx);
        (*degs)[index++] = v->length - 1;
    }

    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j]->length - 1 > (j + 1)*l || j == 0)
        {
            fmpz_mod_poly_set(g, I[j], ctx);
            for (i = l - 1; i >= 0 && (g->length > 1); i--)
            {
                /* compute f^{[l*(j+1)-i]} */
                fmpz_mod_poly_sub(tmp, H[j], h[i], ctx);
                fmpz_mod_poly_gcd(f, g, tmp, ctx);
                if (f->length > 1)
                {
                    /* insert f^{[l*(j+1)-i]} into res */
                    fmpz_mod_poly_make_monic(f, f, ctx);
                    fmpz_mod_poly_factor_insert(res, f, 1, ctx);
                    (*degs)[index++] = l * (j + 1) - i;

                    fmpz_mod_poly_remove(g, f, ctx);
                }
            }
        }
        else if (I[j]->length > 1)
        {
            fmpz_mod_poly_make_monic(I[j], I[j], ctx);
            fmpz_mod_poly_factor_insert(res, I[j], 1, ctx);
            (*degs)[index++] = I[j]->length - 1;
        }
    }

    /* cleanup */
    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_poly_clear(g, ctx);
    fmpz_mod_poly_clear(v, ctx);
    fmpz_mod_poly_clear(vinv, ctx);
    fmpz_mod_poly_clear(tmp, ctx);

    fmpz_mat_clear(HH);

    for (i = 0; i < l + 1; i++)
        fmpz_mod_poly_clear(h[i], ctx);
    for (i = 0; i < m; i++)
    {
        fmpz_mod_poly_clear(H[i], ctx);
        fmpz_mod_poly_clear(I[i], ctx);
    }
    flint_free(h);
}
