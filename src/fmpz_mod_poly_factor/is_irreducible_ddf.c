/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013, 2014 Martin Lee
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_mod_poly_factor.h"

#ifdef __GNUC__
# define ceil __builtin_ceil
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

int fmpz_mod_poly_is_irreducible_ddf(const fmpz_mod_poly_t poly,
                                                      const fmpz_mod_ctx_t ctx)
{
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_mod_poly_t f, v, vinv, tmp;
    fmpz_mod_poly_t *h, *H, *I;
    fmpz_mat_t HH;
    slong i, j, l, m, n, d;
    double beta;
    int result = 1;

    n = fmpz_mod_poly_degree(poly, ctx);

    if (n < 2)
        return 1;

    if (!fmpz_mod_poly_is_squarefree(poly, ctx))
        return 0;

    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    /* initialization */

    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_init(v, ctx);
    fmpz_mod_poly_init(vinv, ctx);
    fmpz_mod_poly_init(tmp, ctx);

    h = flint_malloc((2 * m + l + 1) * sizeof(fmpz_mod_poly_struct));
    H = h + (l + 1);
    I = H + m;

    for (i = 0; i < 2*m + l + 1; i++)
        fmpz_mod_poly_init(h[i], ctx);

    fmpz_mod_poly_make_monic(v, poly, ctx);

    fmpz_mod_poly_reverse(vinv, v, v->length, ctx);
    fmpz_mod_poly_inv_series(vinv, vinv, v->length, ctx);
    /* compute baby steps: h[i]=x^{p^i}mod v */
    fmpz_mod_poly_set_coeff_ui(h[0], 1, 1, ctx);
    fmpz_mod_poly_powmod_x_fmpz_preinv(h[1], p, v, vinv, ctx);
    if (fmpz_sizeinbase(p, 2) > ((n_sqrt(v->length - 1) + 1) * 3) / 4)
    {
        for (i= 1; i < FLINT_BIT_COUNT (l); i++)
            fmpz_mod_poly_compose_mod_brent_kung_vec_preinv (*(h + 1 +
                                                             (1 << (i - 1))),
                                                             *(h + 1),
                                                             (1 << (i - 1)),
                                                             (1 << (i - 1)),
                                                             *(h + (1 << (i - 1))),
                                                             v, vinv, ctx);
        fmpz_mod_poly_compose_mod_brent_kung_vec_preinv (*(h + 1 +
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
    fmpz_mod_poly_set(H[0], h[l], ctx);
    fmpz_mat_init(HH, n_sqrt(v->length - 1) + 1, v->length - 1);
    fmpz_mod_poly_precompute_matrix(HH, H[0], v, vinv, ctx);
    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[i]=x^{p^(li)}mod v */
        if (j > 0)
            fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(H[j],
                                                   H[j - 1], HH, v, vinv, ctx);
        /* compute interval polynomials */
        fmpz_mod_poly_set_coeff_ui(I[j], 0, 1, ctx);
        for (i = l - 1; (i >= 0) && (2*d <= v->length - 1); i--, d++)
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
            result = 0;
            break;
        }
    }

    fmpz_mod_poly_clear(f, ctx);
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

    return result;
}
