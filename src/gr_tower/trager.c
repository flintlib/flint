/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq_poly.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/*
    Trager's method: the norm N(M) of a polynomial M over the field
    F_{k-1} of a tower over QQ is the iterated resultant of M with the
    moduli of the steps below, a polynomial over QQ of degree
    deg(M) [F_{k-1} : QQ]. If the norm of M_s(x) = m(x - s theta) is
    squarefree for a shift s (theta a combination of the generators),
    then m is irreducible over F_{k-1} if and only if the norm is
    irreducible over QQ, and in general the factors of m are the gcds of
    M_s with the irreducible factors of the norm (shifted back).

    This is a complete but expensive proof (the norm has degree
    deg(m) [F_{k-1} : QQ]), used when the modular test fails, for
    moduli of moderate total degree. The norms and the factorization are
    those of factor.c; with transcendental generators, the norm is a
    polynomial over QQ(t_1, ..., t_r), factored as a multivariate
    polynomial.
*/

/*
    Decides the irreducibility of the modulus of step k over F_{k-1} by
    Trager's method, when the steps below are proven and
    deg(m_k) [F_{k-1} : F_0] <= degree_limit. Sets the status of the step
    to PROVEN if it is irreducible; otherwise refines the tower with the
    irreducible factor vanishing at the generator, which is then proven.
    Returns 1 if the step is proven (after refinement, if any), 0 if the
    method does not apply or fails (a shift making the norm squarefree
    not found within a few tries, a non-squarefree modulus).
*/
int
gr_tower_prove_step_trager(gr_tower_t T, slong k, slong degree_limit)
{
    gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
    gr_ctx_struct * F = gr_tower_field_at(T, k - 1);
    const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
    slong d = m->length - 1, D = 1, j;
    gr_ctx_t pctx;
    gr_vec_t fac;
    int result = 0, status;

    if (g->status == GR_TOWER_STATUS_PROVEN)
        return 1;
    if (d < 2)
    {
        _gr_tower_gen_set_status(g, GR_TOWER_STATUS_PROVEN);
        return 1;
    }
    if (!GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_FACTOR))
        return 0;
    for (j = 1; j < k; j++)
    {
        if (GR_TOWER_STEP(T, j - 1)->status != GR_TOWER_STATUS_PROVEN)
            return 0;
    }
    D = gr_tower_degree_at(T, k - 1);
    if (D > degree_limit / d)
        return 0;

    gr_ctx_init_gr_poly(pctx, F);
    gr_vec_init(fac, 0, pctx);

    if (k == 1)
    {
        fmpz_vec_t e;
        fmpz_vec_init(e, 0);
        status = _gr_tower_base_poly_factor(fac, e, m, T);
        for (j = 0; j < e->length; j++)
            if (!fmpz_is_one(e->entries + j))
                status = GR_UNABLE;
        fmpz_vec_clear(e);
    }
    else
        status = _gr_tower_poly_factor_squarefree_trager(fac, m, k - 1, T);

    if (status == GR_SUCCESS)
    {
        if (fac->length == 1)
        {
            _gr_tower_gen_set_status(g, GR_TOWER_STATUS_PROVEN);
            result = 1;
        }
        else
        {
            /* the modulus is divided by the factors not vanishing at the
               generator, one at a time, until a single (irreducible)
               factor is left: after refining with f_j the modulus is
               either f_j itself (irreducible) or the product of the
               factors not yet tried */
            slong remaining = fac->length;

            for (j = 0; j < fac->length; j++)
            {
                const gr_poly_struct * fj = gr_vec_entry_ptr(fac, j, pctx);
                const gr_poly_struct * cur;
                slong cur_len = gr_tower_step_minpoly(T, k)->length;

                if (gr_tower_refine_step(T, k, fj) != GR_SUCCESS)
                    break;

                cur = gr_tower_step_minpoly(T, k);
                remaining--;

                if (cur->length != fj->length)
                {
                    /* the modulus is the product of the remaining factors */
                    if (remaining == 1)
                    {
                        result = 1;
                        break;
                    }
                    continue;
                }

                if (cur_len - fj->length + 1 != fj->length || remaining == 1)
                {
                    /* f_j, or a single remaining factor of the same degree */
                    result = 1;
                    break;
                }
                else
                {
                    /* the same degree as the remaining product: decide by
                       divisibility (the modulus is f_j up to a unit, or
                       the product of at least two other factors) */
                    gr_poly_t q, r;
                    int divides;
                    gr_poly_init(q, F);
                    gr_poly_init(r, F);
                    divides = (gr_poly_divrem(q, r, fj, cur, F) == GR_SUCCESS &&
                               gr_poly_is_zero(r, F) == T_TRUE);
                    gr_poly_clear(q, F);
                    gr_poly_clear(r, F);
                    if (divides)
                    {
                        result = 1;
                        break;
                    }
                }
            }

            if (result)
                _gr_tower_gen_set_status(g, GR_TOWER_STATUS_PROVEN);
        }
    }

    gr_vec_clear(fac, pctx);
    gr_ctx_clear(pctx);
    return result;
}

POP_OPTIONS
