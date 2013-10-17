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

#include "ulong_extras.h"
#include "fq_poly.h"

int
fq_poly_factor_equal_deg_prob(fq_poly_t factor, flint_rand_t state,
                              const fq_poly_t pol, slong d, const fq_ctx_t ctx)
{
    fq_poly_t a, b, c;
    fq_t t;
    fmpz_t exp, q;
    int res = 1;
    slong i, k;

    if (pol->length <= 1)
    {
        flint_printf("Exception (fq_poly_factor_equal_deg_prob): \n");
        flint_printf("Input polynomial is linear.\n");
        abort();
    }

    fmpz_init(q);
    fq_ctx_order(q, ctx);

    fq_poly_init(a, ctx);

    do
    {
        fq_poly_randtest(a, state, pol->length - 1, ctx);
    } while (a->length <= 1);

    fq_poly_gcd(factor, a, pol, ctx);

    if (factor->length != 1)
    {
        fq_poly_clear(a, ctx);
        return 1;
    }

    fq_poly_init(b, ctx);

    fmpz_init(exp);
    if (fmpz_cmp_ui(fq_ctx_prime(ctx), 2) > 0)
    {
        /* compute a^{(q^d-1)/2} rem pol */
        fmpz_pow_ui(exp, q, d);
        fmpz_sub_ui(exp, exp, 1);
        fmpz_fdiv_q_2exp(exp, exp, 1);

        fq_poly_powmod_fmpz_binexp(b, a, exp, pol, ctx);
    }
    else
    {
        /* compute b = (a^{2^{k*d-1}}+a^{2^{k*d-2}}+...+a^4+a^2+a) rem pol */
        k = d * fq_ctx_degree(ctx); /* TODO: Handle overflow? */
        fq_poly_rem(b, a, pol, ctx);
        fq_poly_init(c, ctx);
        fq_poly_set(c, b, ctx);
        for (i = 1; i < k; i++)
        {
            /* c = a^{2^i} = (a^{2^{i-1}})^2 */
            fq_poly_powmod_ui_binexp(c, c, 2, pol, ctx);
            fq_poly_add(b, b, c, ctx);
        }
        fq_poly_rem(b, b, pol, ctx);
        fq_poly_clear(c, ctx);
    }
    fmpz_clear(exp);

    fq_init(t, ctx);
    fq_sub_one(t, b->coeffs + 0, ctx);
    fq_poly_set_coeff(b, 0, t, ctx);
    fq_clear(t, ctx);

    fq_poly_gcd(factor, b, pol, ctx);

    if ((factor->length <= 1) || (factor->length == pol->length))
        res = 0;

    fq_poly_clear(a, ctx);
    fq_poly_clear(b, ctx);
    fmpz_clear(q);

    return res;
}
