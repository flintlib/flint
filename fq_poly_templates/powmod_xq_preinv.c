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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/

void
TEMPLATE(T, poly_powmod_xq_preinv)(TEMPLATE(T, poly_t) rop,
                                   const TEMPLATE(T, poly_t) f,
                                   const TEMPLATE(T, poly_t) finv,
                                   const TEMPLATE(T, ctx_t) ctx)
{
    int i;
    slong j, degree;
    
    TEMPLATE(T, t) ap, api, temp;
    TEMPLATE(T, poly_t) xp, xpi, new_xpi;
    TEMPLATE(T, ground_mat_t) ap_mat, api_mat;
    TEMPLATE(T, mat_t) xp_mat;

    degree = TEMPLATE(T, ctx_degree)(ctx);
    TEMPLATE(T, init)(temp, ctx);
    
    /* ap = a ^ p */
    TEMPLATE(T, init)(ap, ctx);
    TEMPLATE(T, init)(api, ctx);
    TEMPLATE(T, gen)(ap, ctx);
    TEMPLATE(T, pow)(ap, ap, TEMPLATE(T, ctx_prime)(ctx), ctx);
    TEMPLATE(T, set)(api, ap, ctx);
    
    /* xp = x ^ p mod f */
    TEMPLATE(T, poly_init)(xp, ctx);
    TEMPLATE(T, poly_init)(xpi, ctx);
    TEMPLATE(T, poly_gen)(xp, ctx);
    TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(xp, xp, TEMPLATE(T, ctx_prime)(ctx), 0,
                                                 f, finv, ctx);
    TEMPLATE(T, poly_set)(xpi, xp, ctx);

    /* Precompute compose_mod matrices */
    TEMPLATE(T, pow_pn_init_precomp_matrix)(ap_mat, ctx);
    TEMPLATE3(T, pow_pn_precompute_matrix, T)(ap_mat, ap, ctx);

    /* TODO: Just copy the matrix */
    TEMPLATE(T, pow_pn_init_precomp_matrix)(api_mat, ctx);
    TEMPLATE3(T, pow_pn_precompute_matrix, T)(api_mat, api, ctx);

    TEMPLATE(T, mat_init) (xp_mat, n_sqrt(f->length-1) + 1, f->length-1, ctx);
    TEMPLATE(T, poly_precompute_matrix)(xp_mat, xp, f, finv, ctx);
    
    TEMPLATE(T, poly_init)(new_xpi, ctx);
    TEMPLATE(T, poly_fit_length)(new_xpi, f->length - 1, ctx);

    for (i = ((int) FLINT_BIT_COUNT(degree) - 2); i >= 0; i--)
    {
        for (j = 0; j < xpi->length; j++)
        {
            TEMPLATE(T, pow_pn_precomp)(new_xpi->coeffs + j, xpi->coeffs + j, api_mat, ctx);
        }
        _TEMPLATE(T, poly_set_length)(new_xpi, xpi->length, ctx);
        TEMPLATE(T, poly_compose_mod_preinv)(xpi, new_xpi, xpi, f, finv, ctx);

        /* api = api \circ api */
        if (i)
            TEMPLATE(T, pow_pn_precomp)(api, api, api_mat, ctx);
        
        if (degree & (WORD(1) << i))
        {

            for (j = 0; j < xpi->length; j++)
            {
                TEMPLATE(T, pow_pn_precomp)(new_xpi->coeffs + j, xpi->coeffs + j, ap_mat, ctx);
            }
            _TEMPLATE(T, poly_set_length)(new_xpi, xpi->length, ctx);
            TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(xpi, new_xpi, xp_mat,
                                                                    f, finv, ctx);

            /* api = api \circ ap */
            if (i)
                TEMPLATE(T, pow_pn_precomp)(api, api, ap_mat, ctx);
        }

        if (i)
            TEMPLATE3(T, pow_pn_precompute_matrix, T)(api_mat, api, ctx);

    }

    TEMPLATE(T, poly_set)(rop, xpi, ctx);

/*    
    for (i = 1; i < degree; i++)
    {
        for (j = 0; j < rop->length; j++)
        {
            TEMPLATE(T, pow_pn_precomp)(temp, rop->coeffs + j, ap_mat, ctx);
            TEMPLATE(T, poly_set_coeff)(rop, j, temp, ctx);
        }
        TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(rop, rop, xp_mat,
                                                                f, finv, ctx);
    }
*/

    TEMPLATE(T, mat_clear)(xp_mat, ctx);
    TEMPLATE(T, pow_pn_clear_precomp_matrix)(ap_mat, ctx);
    TEMPLATE(T, pow_pn_clear_precomp_matrix)(api_mat, ctx);
    TEMPLATE(T, poly_clear)(xp, ctx);
    TEMPLATE(T, poly_clear)(xpi, ctx);
    TEMPLATE(T, clear)(ap, ctx);
    TEMPLATE(T, clear)(api, ctx);
    TEMPLATE(T, clear)(temp, ctx);
}
