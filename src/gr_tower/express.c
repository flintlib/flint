/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "qqbar.h"
#include "gr_vec.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* Enclosures of the monomial basis elements of F_n over F_0, in the
   same order as gr_tower_get_coeffs_at. */
static int
_basis_get_acb(acb_ptr res, slong prec, gr_tower_t T)
{
    slong k, i, j, D = 1;
    acb_t z;
    int status = GR_SUCCESS;

    acb_init(z);
    acb_one(res);

    for (k = 1; k <= T->length && status == GR_SUCCESS; k++)
    {
        slong d = gr_tower_step_degree(T, k);

        status |= gr_tower_step_get_acb(z, T, k, prec);

        /* res[i + j D] = res[i] * z^j for 0 <= i < D, 1 <= j < d */
        for (j = 1; j < d; j++)
            for (i = 0; i < D; i++)
                acb_mul(res + i + j * D, res + i + (j - 1) * D, z, prec);

        D *= d;
    }

    acb_clear(z);
    return status;
}

/*
    Given that q(x) = 0 exactly and that z isolates one root of q, decides
    whether x is that root. The enclosure w of x is refined until either
    w is contained in z, or w is disjoint from z, or the disk z inflated
    to contain w still isolates a single root of q (verified by an
    interval Newton step). The last case handles enclosures of different
    shapes, e.g. an exactly real z and a complex w with a tiny imaginary
    radius.
*/
static int
_root_matches(gr_srcptr x, const gr_poly_t q, const acb_t z, gr_tower_t T)
{
    acb_t w, z2;
    acb_poly_t qz;
    mag_t r;
    slong prec;
    int res = -1;

    acb_init(w);
    acb_init(z2);
    acb_poly_init(qz);
    mag_init(r);

    for (prec = GR_TOWER_DEFAULT_PREC; prec <= 1000000 && res == -1; prec *= 2)
    {
        if (gr_tower_get_acb(w, x, prec, T) != GR_SUCCESS)
            break;

        if (acb_contains(z, w))
        {
            res = 1;
        }
        else if (!acb_overlaps(z, w))
        {
            res = 0;
        }
        else
        {
            /* inflate z by the diameter of w so that z2 contains w */
            acb_set(z2, z);
            mag_max(r, arb_radref(acb_realref(w)), arb_radref(acb_imagref(w)));
            mag_mul_2exp_si(r, r, 2);
            acb_add_error_mag(z2, r);

            if (acb_contains(z2, w) &&
                _gr_tower_poly_get_acb_poly(qz, q, T->length, prec, T) == GR_SUCCESS &&
                _gr_tower_newton_step(z2, qz, z2, prec))
            {
                res = 1;
            }
        }
    }

    acb_clear(w);
    acb_clear(z2);
    acb_poly_clear(qz);
    mag_clear(r);
    return res == 1;
}

int
gr_tower_express_limit(gr_ptr res, const gr_poly_t q, const acb_t z, slong prec_limit_in, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_ctx_struct * F = T->base;
    slong D = gr_tower_degree(T);
    slong i, prec, prec_limit;
    acb_ptr vec;
    fmpz * rel;
    gr_ptr coeffs, val;
    fmpq_t c;
    int status = GR_UNABLE;
    slong rnd, max_rounds;

    /* the search is for a Q-linear combination of the monomial basis;
       over a transcendental base only such combinations are found */
    if (F->which_ring != GR_CTX_FMPQ && F->which_ring != GR_CTX_FMPZ_MPOLY_Q)
        return GR_UNABLE;

    if (q->length < 2)
        return GR_DOMAIN;

    /* default precision limit: the coordinates can have about
       deg(q) * height(q) bits, times D for the common denominator */
    if (prec_limit_in <= 0)
        prec_limit = 64 * (D + 1) * 8;
    else
        prec_limit = prec_limit_in;
    if (prec_limit > 100000)
        prec_limit = 100000;

    vec = _acb_vec_init(D + 1);
    rel = _fmpz_vec_init(D + 1);
    coeffs = gr_heap_init_vec(D, F);
    fmpq_init(c);
    GR_TMP_INIT(val, top);

    /* (at most max_rounds precisions, GR_TOWER_OPT_EXPRESS_ROUNDS: a
       relation is nearly always found at the first or second, and a
       search which finds none runs up to the limit, the last rounds
       being the costliest) */
    max_rounds = GR_TOWER_OPTION(T, GR_TOWER_OPT_EXPRESS_ROUNDS);
    for (prec = FLINT_MAX(GR_TOWER_DEFAULT_PREC, 2 * (D + 1)), rnd = 0;
         prec <= prec_limit && (max_rounds <= 0 || rnd < max_rounds); prec *= 2, rnd++)
    {
        acb_t zz;
        int st;

        if (_basis_get_acb(vec, prec, T) != GR_SUCCESS)
            break;

        /* the target value at the working precision; refine z using q */
        acb_init(zz);
        {
            acb_poly_t qz;
            acb_poly_init(qz);
            acb_set(zz, z);
            if (_gr_tower_poly_get_acb_poly(qz, q, T->length, prec, T) == GR_SUCCESS)
            {
                slong iter;
                for (iter = 0; iter < 10 && acb_rel_accuracy_bits(zz) < prec; iter++)
                {
                    acb_t t;
                    acb_init(t);
                    if (_gr_tower_newton_step(t, qz, zz, prec))
                        acb_set(zz, t);
                    else
                    {
                        acb_clear(t);
                        break;
                    }
                    acb_clear(t);
                }
            }
            acb_poly_clear(qz);
        }
        acb_set(vec + D, zz);
        acb_clear(zz);

        if (!_qqbar_acb_lindep(rel, vec, D + 1, 1, prec) || fmpz_is_zero(rel + D))
            continue;

        /* candidate: res = -sum rel[i] b_i / rel[D] */
        st = GR_SUCCESS;
        for (i = 0; i < D; i++)
        {
            fmpz_neg(fmpq_numref(c), rel + i);
            fmpz_set(fmpq_denref(c), rel + D);
            fmpq_canonicalise(c);
            st |= gr_set_fmpq(GR_ENTRY(coeffs, i, F->sizeof_elem), c, F);
        }
        st |= gr_tower_set_coeffs_at(res, coeffs, T->length, T);

        if (st != GR_SUCCESS)
            continue;

        /* exact verification: q(res) = 0 and res lies in the isolating disk */
        st = gr_poly_evaluate(val, q, res, top);
        if (st == GR_SUCCESS && gr_tower_is_zero(val, T) == T_TRUE && _root_matches(res, q, z, T))
        {
            status = GR_SUCCESS;
            break;
        }
    }

    _acb_vec_clear(vec, D + 1);
    _fmpz_vec_clear(rel, D + 1);
    gr_heap_clear_vec(coeffs, D, F);
    fmpq_clear(c);
    GR_TMP_CLEAR(val, top);

    return status;
}

int
gr_tower_express(gr_ptr res, const gr_poly_t q, const acb_t z, gr_tower_t T)
{
    return gr_tower_express_limit(res, q, z, 0, T);
}

int
gr_tower_express_qqbar(gr_ptr res, const qqbar_t x, gr_tower_t T)
{
    gr_ctx_struct * top = gr_tower_field(T);
    gr_poly_t q;
    int status;

    if (!GR_TOWER_BASE_IS_CONSTS(T) || !GR_TOWER_HAS_CAP(T, GR_TOWER_CAP_RATIONAL))
        return GR_UNABLE;

    if (qqbar_is_rational(x))
    {
        fmpq_t c;
        fmpq_init(c);
        qqbar_get_fmpq(c, x);
        status = gr_set_fmpq(res, c, top);
        fmpq_clear(c);
        return status;
    }

    /* the degree of x over QQ must divide the degree of the field over
       QQ, which is gr_tower_degree(T) when the tower is algebraic with
       all steps proven irreducible, and otherwise only known to be at
       most that */
    if (gr_tower_degree(T) % qqbar_degree(x) != 0)
        return (T->num_trans == 0 && _gr_tower_all_proven(T)) ? GR_DOMAIN : GR_UNABLE;

    gr_poly_init(q, top);
    status = gr_poly_set_fmpz_poly(q, QQBAR_POLY(x), top);
    if (status == GR_SUCCESS)
        status = gr_tower_express(res, q, QQBAR_ENCLOSURE(x), T);
    gr_poly_clear(q, top);

    return status;
}

POP_OPTIONS
