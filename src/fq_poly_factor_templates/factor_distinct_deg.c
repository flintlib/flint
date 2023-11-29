/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

#include <math.h>

void
TEMPLATE(T, poly_factor_distinct_deg) (TEMPLATE(T, poly_factor_t) res,
                                       const TEMPLATE(T, poly_t) poly,
                                       slong * const *degs,
                                       const TEMPLATE(T, ctx_t) ctx)
{
    TEMPLATE(T, poly_t) f, g, s, reducedH0, v, vinv, tmp;
    TEMPLATE(T, poly_t) * h, *H, *I;
    fmpz_t q;
    slong i, j, l, m, n, index, d;
    double beta;
    TEMPLATE(T, mat_t) HH, HHH;

    TEMPLATE(T, poly_init) (v, ctx);
    TEMPLATE(T, poly_make_monic) (v, poly, ctx);

    n = TEMPLATE(T, poly_degree) (poly, ctx);
    if (n == 1)
    {
        TEMPLATE(T, poly_factor_insert) (res, poly, 1, ctx);
        (*degs)[0] = 1;
        TEMPLATE(T, poly_clear) (v, ctx);
        return;
    }

    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    fmpz_init(q);
    TEMPLATE(T, ctx_order) (q, ctx);

    TEMPLATE(T, poly_init) (f, ctx);
    TEMPLATE(T, poly_init) (g, ctx);
    TEMPLATE(T, poly_init) (s, ctx);
    TEMPLATE(T, poly_init) (reducedH0, ctx);
    TEMPLATE(T, poly_init) (vinv, ctx);
    TEMPLATE(T, poly_init) (tmp, ctx);

    h = flint_malloc((2 * m + l + 1) * sizeof(TEMPLATE(T, poly_struct)));
    H = h + (l + 1);
    I = H + m;
    for (i = 0; i < l + 1; i++)
        TEMPLATE(T, poly_init) (h[i], ctx);
    for (i = 0; i < m; i++)
    {
        TEMPLATE(T, poly_init) (H[i], ctx);
        TEMPLATE(T, poly_init) (I[i], ctx);
    }

    TEMPLATE(T, poly_make_monic) (v, poly, ctx);

    TEMPLATE(T, poly_reverse) (vinv, v, v->length, ctx);
    TEMPLATE(T, poly_inv_series_newton) (vinv, vinv, v->length, ctx);

    /* compute baby steps: h[i]=x^{q^i}mod v */
    /*     h[0] = x */

    TEMPLATE(T, poly_iterated_frobenius_preinv) (h, l + 1, v, vinv, ctx);

    /* compute coarse distinct-degree factorisation */
    index = 0;
    TEMPLATE(T, poly_set) (s, v, ctx);
    TEMPLATE(T, poly_set) (H[0], h[l], ctx);
    TEMPLATE(T, poly_set) (reducedH0, H[0], ctx);
    TEMPLATE(T, mat_init) (HH, n_sqrt(v->length - 1) + 1, v->length - 1, ctx);
    TEMPLATE(T, poly_precompute_matrix) (HH, reducedH0, s, vinv, ctx);

    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[j]=x^{q^(lj)}mod s */
        if (j > 0)
        {
            if (I[j - 1]->length > 1)
            {
                _TEMPLATE(T, poly_reduce_matrix_mod_poly) (HHH, HH, s, ctx);
                TEMPLATE(T, mat_clear) (HH, ctx);
                TEMPLATE(T, mat_init_set) (HH, HHH, ctx);
                TEMPLATE(T, mat_clear) (HHH, ctx);
                TEMPLATE(T, poly_rem) (reducedH0, reducedH0, s, ctx);
                TEMPLATE(T, poly_rem) (tmp, H[j - 1], s, ctx);
                TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)
                    (H[j], tmp, HH, s, vinv, ctx);
            }
            else
            {
                TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)
                    (H[j], H[j - 1], HH, s, vinv, ctx);

            }
        }

        /* compute interval polynomials */
        TEMPLATE(T, poly_one) (I[j], ctx);
        for (i = l - 1; (i >= 0) && (2 * d <= s->length - 1); i--, d++)
        {
            TEMPLATE(T, poly_rem) (tmp, h[i], s, ctx);
            TEMPLATE(T, poly_sub) (tmp, H[j], tmp, ctx);
            TEMPLATE(T, poly_mulmod_preinv) (I[j], tmp, I[j], s, vinv, ctx);
        }

        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        TEMPLATE(T, poly_gcd) (I[j], s, I[j], ctx);
        if (I[j]->length > 1)
        {
            TEMPLATE(T, poly_remove) (s, I[j], ctx);
            TEMPLATE(T, poly_reverse) (vinv, s, s->length, ctx);
            TEMPLATE(T, poly_inv_series_newton) (vinv, vinv, s->length, ctx);
        }
        if (s->length - 1 < 2 * d)
        {
            break;
        }
    }
    if (s->length > 1)
    {
        TEMPLATE(T, poly_factor_insert) (res, s, 1, ctx);
        (*degs)[index++] = s->length - 1;
    }


    /* compute fine distinct-degree factorisation */
    for (j = 0; j < m; j++)
    {
        if (I[j]->length - 1 > (j + 1) * l || j == 0)
        {
            TEMPLATE(T, poly_set) (g, I[j], ctx);
            for (i = l - 1; i >= 0 && (g->length > 1); i--)
            {
                /* compute f^{[l*(j+1)-i]} */
                TEMPLATE(T, poly_sub) (tmp, H[j], h[i], ctx);
                TEMPLATE(T, poly_gcd) (f, g, tmp, ctx);
                if (f->length > 1)
                {
                    /* insert f^{[l*(j+1)-i]} into res */
                    TEMPLATE(T, poly_make_monic) (f, f, ctx);
                    TEMPLATE(T, poly_factor_insert) (res, f, 1, ctx);
                    (*degs)[index++] = l * (j + 1) - i;

                    TEMPLATE(T, poly_remove) (g, f, ctx);
                }
            }
        }
        else if (I[j]->length > 1)
        {
            TEMPLATE(T, poly_make_monic) (I[j], I[j], ctx);
            TEMPLATE(T, poly_factor_insert) (res, I[j], 1, ctx);
            (*degs)[index++] = I[j]->length - 1;
        }
    }

    /* cleanup */
    fmpz_clear(q);
    TEMPLATE(T, poly_clear) (f, ctx);
    TEMPLATE(T, poly_clear) (g, ctx);
    TEMPLATE(T, poly_clear) (s, ctx);
    TEMPLATE(T, poly_clear) (reducedH0, ctx);
    TEMPLATE(T, poly_clear) (v, ctx);
    TEMPLATE(T, poly_clear) (vinv, ctx);
    TEMPLATE(T, poly_clear) (tmp, ctx);
    TEMPLATE(T, mat_clear) (HH, ctx);

    for (i = 0; i < l + 1; i++)
        TEMPLATE(T, poly_clear) (h[i], ctx);
    for (i = 0; i < m; i++)
    {
        TEMPLATE(T, poly_clear) (H[i], ctx);
        TEMPLATE(T, poly_clear) (I[i], ctx);
    }
    flint_free(h);
}


#endif
