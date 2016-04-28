/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013, 2014 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#undef ulong
#define ulong ulongxx/* interferes with system includes */

#include <math.h>

#undef ulong

#include <gmp.h>

#define ulong mp_limb_t

#include "fmpz_mod_poly.h"

int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t poly)
{
    fmpz_mod_poly_t f, v, vinv, tmp;
    fmpz_mod_poly_t *h, *H, *I;
    fmpz_mat_t HH;
    slong i, j, l, m, n, d;
    fmpz_t p;
    double beta;
    int result = 1;

    n = fmpz_mod_poly_degree(poly);

    if (n < 2)
        return 1;

    if (!fmpz_mod_poly_is_squarefree(poly))
        return 0;

    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    /* initialization */
    fmpz_init(p);
    fmpz_set(p, &poly->p);

    fmpz_mod_poly_init(f, p);
    fmpz_mod_poly_init(v, p);
    fmpz_mod_poly_init(vinv, p);
    fmpz_mod_poly_init(tmp, p);

    if (!(h = flint_malloc((2 * m + l + 1) * sizeof(fmpz_mod_poly_struct))))
    {
        flint_printf("Exception (fmpz_mod_poly_is_irreducible_ddf): \n");
        flint_printf("Not enough memory.\n");
        flint_abort();
    }
    H = h + (l + 1);
    I = H + m;
    fmpz_mod_poly_init(h[0], p);
    fmpz_mod_poly_init(h[1], p);
    for (i = 0; i < m; i++)
    {
        fmpz_mod_poly_init(H[i], p);
        fmpz_mod_poly_init(I[i], p);
    }

    fmpz_mod_poly_make_monic(v, poly);

    fmpz_mod_poly_reverse (vinv, v, v->length);
    fmpz_mod_poly_inv_series_newton (vinv, vinv, v->length);
    /* compute baby steps: h[i]=x^{p^i}mod v */
    fmpz_mod_poly_set_coeff_ui(h[0], 1, 1);
    fmpz_mod_poly_powmod_x_fmpz_preinv(h[1], p, v, vinv);
    if (fmpz_sizeinbase(p, 2) > ((n_sqrt(v->length - 1) + 1) * 3) / 4)
    {
        for (i= 1; i < FLINT_BIT_COUNT (l); i++)
            fmpz_mod_poly_compose_mod_brent_kung_vec_preinv (*(h + 1 +
                                                             (1 << (i - 1))),
                                                             *(h + 1),
                                                             (1 << (i - 1)),
                                                             (1 << (i - 1)), v,
                                                             vinv);
        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv (*(h + 1 +
                                                         (1 << (i - 1))),
                                                         *(h + 1),
                                                         (1 << (i - 1)),
                                                         l - (1 << (i - 1)), v,
                                                         vinv);
    }
    else
    {
        for (i = 2; i < l + 1; i++)
        {
            fmpz_mod_poly_init(h[i], p);
            fmpz_mod_poly_powmod_fmpz_binexp_preinv(h[i], h[i - 1], p,
                                              v, vinv);
        }
    }

    /* compute coarse distinct-degree factorisation */
    fmpz_mod_poly_set(H[0], h[l]);
    fmpz_mat_init(HH, n_sqrt(v->length - 1) + 1, v->length - 1);
    fmpz_mod_poly_precompute_matrix(HH, H[0], v, vinv);
    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[i]=x^{p^(li)}mod v */
        if (j > 0)
            fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H[j],
                                                        H[j - 1], HH, v, vinv);
        /* compute interval polynomials */
        fmpz_mod_poly_set_coeff_ui(I[j], 0, 1);
        for (i = l - 1; (i >= 0) && (2*d <= v->length - 1); i--, d++)
        {
            fmpz_mod_poly_rem(tmp, h[i], v);
            fmpz_mod_poly_sub(tmp, H[j], tmp);
            fmpz_mod_poly_mulmod_preinv (I[j], tmp, I[j], v, vinv);
        }

        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        fmpz_mod_poly_gcd(I[j], v, I[j]);
        if (I[j]->length > 1)
        {
            result = 0;
            break;
        }
    }

    fmpz_clear(p);
    fmpz_mod_poly_clear(f);
    fmpz_mod_poly_clear(v);
    fmpz_mod_poly_clear(vinv);
    fmpz_mod_poly_clear(tmp);

    fmpz_mat_clear(HH);

    for (i = 0; i < l + 1; i++)
        fmpz_mod_poly_clear(h[i]);
    for (i = 0; i < m; i++)
    {
        fmpz_mod_poly_clear(H[i]);
        fmpz_mod_poly_clear(I[i]);
    }
    flint_free(h);

    return result;
}
