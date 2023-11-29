/*
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Martin Lee
    Copyright (C) 2020, 2022 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "nmod_poly_factor.h"

#ifdef __GNUC__
# define ceil __builtin_ceil
# define log __builtin_log
# define pow __builtin_pow
#else
# include <math.h>
#endif

int nmod_poly_is_irreducible_ddf(const nmod_poly_t poly)
{

    nmod_poly_t f, v, vinv, tmp;
    nmod_poly_struct * h, * H, * I;
    nmod_mat_t HH;
    slong i, j, l, m, n, d;
    double beta;
    int result = 1;
    n = nmod_poly_degree(poly);

    if (n < 2)
        return 1;

    if (!nmod_poly_is_squarefree(poly))
        return 0;

    beta = 0.5 * (1. - (log(2)/log(n)));
    l = ceil(pow (n, beta));
    m = ceil(0.5*n/l);

    /* initialization */
    nmod_poly_init_mod(f, poly->mod);
    nmod_poly_init_mod(v, poly->mod);
    nmod_poly_init_mod(vinv, poly->mod);
    nmod_poly_init_mod(tmp, poly->mod);

    h =  flint_malloc((2 * m + l + 1) * sizeof(nmod_poly_struct));
    H = h + (l + 1);
    I = H + m;

    for (i = 0; i < 2*m + l + 1; i++)
        nmod_poly_init_mod(h + i, poly->mod);

    nmod_poly_make_monic(v, poly);

    nmod_poly_reverse(vinv, v, v->length);
    nmod_poly_inv_series(vinv, vinv, v->length);

    /* compute baby steps: h[i] = x^{p^i}mod v */
    nmod_poly_set_coeff_ui(h + 0, 1, 1);
    nmod_poly_powmod_x_ui_preinv(h + 1, poly->mod.n, v, vinv);

    if (FLINT_BIT_COUNT(poly->mod.n) > ((n_sqrt(v->length - 1) + 1)*3)/4)
    {
        for (i = 1; i < FLINT_BIT_COUNT (l); i++)
            nmod_poly_compose_mod_brent_kung_vec_preinv(h + 1 +
                            (1 << (i - 1)), h + 1, (1 << (i - 1)),
                            (1 << (i - 1)), h + (1 << (i - 1)), v, vinv);

        nmod_poly_compose_mod_brent_kung_vec_preinv(h + 1 + (1 << (i - 1)),
                            h + 1, (1 << (i - 1)), l - (1 << (i - 1)),
						    h + (1 << (i - 1)), v, vinv);
    }
    else
    {
        for (i = 2; i < l + 1; i++)
        {
            nmod_poly_init_mod(h + i, poly->mod);

            nmod_poly_powmod_ui_binexp_preinv(h + i, h + i - 1, poly->mod.n,
                                              v, vinv);
        }
    }

    /* compute coarse distinct-degree factorisation */
    nmod_poly_set(H + 0, h + l);
    nmod_mat_init(HH, n_sqrt(v->length - 1) + 1, v->length - 1, poly->mod.n);
    nmod_poly_precompute_matrix(HH, H + 0, v, vinv);

    d = 1;
    for (j = 0; j < m; j++)
    {
        /* compute giant steps: H[j] = x^{p^(lj)}mod s */
        if (j > 0)
            nmod_poly_compose_mod_brent_kung_precomp_preinv(H + j, H + j - 1, HH,
                                                            v, vinv);
        /* compute interval polynomials */
        nmod_poly_set_coeff_ui(I + j, 0, 1);

        for (i = l - 1; i >= 0 && 2*d <= v->length - 1; i--, d++)
        {
            nmod_poly_rem(tmp, h + i, v);
            nmod_poly_sub(tmp, H + j, tmp);
            nmod_poly_mulmod_preinv (I + j, tmp, I + j, v, vinv);
        }

        /* compute F_j=f^{[j*l+1]} * ... * f^{[j*l+l]} */
        /* F_j is stored on the place of I_j */
        nmod_poly_gcd(I + j, v, I + j);

        if (I[j].length > 1)
        {
            result = 0;
            break;
        }
    }

    nmod_poly_clear(f);
    nmod_poly_clear(v);
    nmod_poly_clear(vinv);
    nmod_poly_clear(tmp);

    nmod_mat_clear (HH);

    for (i = 0; i < l + 1; i++)
        nmod_poly_clear(h + i);

    for (i = 0; i < m; i++)
    {
        nmod_poly_clear(H + i);
        nmod_poly_clear(I + i);
    }

    flint_free (h);

    return result;
}
