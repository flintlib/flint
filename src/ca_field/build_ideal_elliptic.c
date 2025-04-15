/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "ca.h"
#include "ca_ext.h"
#include "ca_field.h"
#include "fmpz_mpoly.h"

void _ca_field_ideal_insert_clear_mpoly(ca_field_t K, fmpz_mpoly_t poly, fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);

/* This merely implements the Legendre relation between the generators i, j (elliptic K-functions), k, l (elliptic E-functions) and pi_index (pi), if it exists. */
void
ca_field_build_elliptic_relation(ca_field_t K, slong i, slong j, slong k, slong l, slong pi_index, ca_ctx_t ctx)
{
    ca_ptr xi, xj, xk, xl;
    ca_t t, u;
    truth_t e;
    
    ca_init(t, ctx);
    ca_init(u, ctx);
    
    xi = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, i));
    xj = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, j));
    xk = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, k));
    xl = CA_EXT_FUNC_ARGS(CA_FIELD_EXT_ELEM(K, l));
    
    ca_add(t, xi, xj, ctx);
    ca_add(u, xk, xl, ctx);
    
    /* swap indices if necessary */
    if (ca_check_equal(xi, xl, ctx) == T_TRUE){
        ca_swap(xk, xl, ctx);
        slong d = k;
        k = l;
        l = d;
    }
    
    e = truth_and(truth_and(ca_check_equal(xi, xk, ctx),
        ca_check_is_one(t, ctx)), ca_check_is_one(u, ctx));
    
    if (e == T_TRUE)
    {
        fmpz_mpoly_t poly;
        fmpz_mpoly_init(poly, CA_FIELD_MCTX(K, ctx));
        
        slong len = CA_FIELD_LENGTH(K);
        ulong exp[len];
        for (slong v = 0; v < len; v++){
            exp[v] = 0;
        }
        
        exp[i] = 1;
        exp[l] = 1;
        fmpz_mpoly_set_coeff_si_ui(poly, 2, exp, CA_FIELD_MCTX(K, ctx));
        exp[i] = 0;
        exp[l] = 0;
        exp[j] = 1;
        exp[k] = 1;
        fmpz_mpoly_set_coeff_si_ui(poly, 2, exp, CA_FIELD_MCTX(K, ctx));
        exp[k] = 0;
        exp[i] = 1;
        fmpz_mpoly_set_coeff_si_ui(poly, -2, exp, CA_FIELD_MCTX(K, ctx));
        exp[i] = 0;
        exp[j] = 0;
        exp[pi_index] = 1;
        fmpz_mpoly_set_coeff_si_ui(poly, -1, exp, CA_FIELD_MCTX(K, ctx));
        
        _ca_field_ideal_insert_clear_mpoly(K, poly, CA_FIELD_MCTX(K, ctx), ctx);
    }
    
    ca_clear(t, ctx);
    ca_clear(u, ctx);
}

void
ca_field_build_ideal_elliptic(ca_field_t K, ca_ctx_t ctx)
{
    calcium_func_code Fi, Fj, Fk, Fl;
    slong i, j, k, l, len, num_ellK, num_ellE, pi_index;

    len = CA_FIELD_LENGTH(K);

    if (len < 5)
    {
        return;
    }

    num_ellK = 0;
    num_ellE = 0;
    pi_index = -1;
    for (i = 0; i < len; i++)
    {
        Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

        if (Fi == CA_EllipticK)
        {
            num_ellK++;
        }
        else if (Fi == CA_EllipticE)
        {
            num_ellE++;
        }
        else if (Fi == CA_Pi)
        {
            pi_index = i;
        }
    }

    if (num_ellK < 2 || num_ellE < 2 || pi_index == -1)
        return;

    for (i = 0; i < len; i++)
    {
        Fi = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, i));

        if (Fi == CA_EllipticK)
        {
            for (j = i + 1; j < len; j++)
            {
                Fj = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, j));

                if (Fj == CA_EllipticK)
                {
                    for (k = 0; k < len; k++)
                    {
                        Fk = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, k));

                        if (Fk == CA_EllipticE)
                        {
                            for (l = k + 1; l < len; l++)
                            {
                                Fl = CA_EXT_HEAD(CA_FIELD_EXT_ELEM(K, l));

                                if (Fl == CA_EllipticE)
                                {
                                    ca_field_build_elliptic_relation(K, i, j, k, l, pi_index, ctx);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
