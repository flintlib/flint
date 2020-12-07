/*
    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "ulong_extras.h"


typedef struct {
    fmpz_mod_poly_t f, xp, a, g;
} split_struct;

typedef struct {
    fmpz_mod_poly_t f, xp;
} queue_struct;

static void _queue_vec_fit_length(queue_struct ** Q, slong * Qalloc,
                                    slong new_length, const fmpz_mod_ctx_t ctx)
{
    if (new_length > *Qalloc)
    {
        slong i;
        slong old_alloc = *Qalloc;
        slong new_alloc = FLINT_MAX(2*old_alloc, new_length);
        *Q = (queue_struct *) flint_realloc(*Q, new_alloc*sizeof(queue_struct));
        *Qalloc = new_alloc;
        for (i = old_alloc; i < new_alloc; i++)
        {
            fmpz_mod_poly_init((*Q)[i].f, ctx);
            fmpz_mod_poly_init((*Q)[i].xp, ctx);
        }
    }
}

/* given g|f, update Q and res with g and f/g */
void _add_split(
    fmpz_mod_poly_factor_t res,
    queue_struct ** Q_, slong * Qlen_, slong * Qalloc_,
    fmpz_mod_poly_t f,
    fmpz_mod_poly_t g,
    slong d,
    const fmpz_mod_poly_t xp,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_t tmp)
{
    queue_struct * Q = *Q_;
    slong Qlen = *Qlen_;
    slong Qalloc = *Qalloc_;
    slong i, inc = 0;

    _queue_vec_fit_length(&Q, &Qalloc, Qlen + 2, ctx);

    fmpz_mod_poly_divrem(Q[Qlen + 0].f, tmp, f, g, ctx);
    fmpz_mod_poly_swap(Q[Qlen + 1].f, g, ctx);

    if (Q[Qlen + 0].f->length < Q[Qlen + 1].f->length)
        fmpz_mod_poly_swap(Q[Qlen + 0].f, Q[Qlen + 1].f, ctx);

    /* {Q[Qlen + 0], Q[Qlen + 1]} has been sorted by length */
    for (i = 0; i < 2; i++)
    {
        if (fmpz_mod_poly_degree(Q[Qlen + i].f, ctx) > d)
        {
            inc++;
            fmpz_mod_poly_divrem(tmp, Q[Qlen + i].xp, xp, Q[Qlen + i].f, ctx);
        }
        else if (fmpz_mod_poly_degree(Q[Qlen + i].f, ctx) == d)
        {
            fmpz_mod_poly_factor_fit_length(res, res->num + 1, ctx);
            res->exp[res->num] = 1;
            fmpz_mod_poly_set(res->poly + res->num, Q[Qlen + i].f, ctx);
            res->num++;
        }
    }

    *Q_ = Q;
    *Qlen_ = Qlen + inc;
    *Qalloc_ = Qalloc;
}

/* a = b^p^0 + ... + b^p^(d-1) in berlekamp subalgebra, i.e. a^p = a mod f */
static void _compute_trace(
    fmpz_mod_poly_t a,
    fmpz_mod_poly_t b,
    slong d,
    const fmpz_mod_poly_t xp,
    const fmpz_mod_poly_t f,
    const fmpz_mod_poly_t finv,
    const fmpz_mod_ctx_t ctx,
    fmpz_mod_poly_t xi,
    fmpz_mod_poly_t t)
{
    slong i;
    fmpz_mat_t H;
#if FLINT_WANT_ASSERT
    fmpz_mod_poly_t a_check;

    fmpz_mod_poly_init(a_check, ctx);
    fmpz_mod_poly_set(t, b, ctx);
    fmpz_mod_poly_set(a_check, b, ctx);
    for (i = 1; i < d; i++)
    {
        fmpz_mod_poly_powmod_fmpz_binexp(t, t, fmpz_mod_ctx_modulus(ctx), f, ctx);
        fmpz_mod_poly_add(a_check, a_check, t, ctx);
    }
#endif

    fmpz_mat_init(H, n_sqrt(f->length - 1) + 1, f->length - 1);

    if (d < 2)
    {
        fmpz_mod_poly_swap(a, b, ctx);
    }
    else if (d < 16)
    {
        fmpz_mod_poly_precompute_matrix(H, xp, f, finv, ctx);

        fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(xi, b, H, f, finv, ctx);
        fmpz_mod_poly_add(a, b, xi, ctx);
        for (i = 2; i < d; i++)
        {
            fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(t, xi, H, f, finv, ctx);
            fmpz_mod_poly_swap(xi, t, ctx);
            fmpz_mod_poly_add(a, a, xi, ctx);
        }
    }
    else
    {
        /*
            Optimized trace calculation from
                Computing Frobenius maps and factoring polynomials
                    Gathen and Shoup
            underwhelming for small d
        */
        fmpz_mod_poly_zero(a, ctx);
        fmpz_mod_poly_set(xi, xp, ctx); /* x^p^2^i*/

        while (1)
        {
            FLINT_ASSERT(d > 1);

            fmpz_mod_poly_precompute_matrix(H, xi, f, finv, ctx);

            if (d % 2 == 0)
            {
                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(t, b, H, f,
                                                                    finv, ctx);
                fmpz_mod_poly_add(b, b, t, ctx);
            }
            else if (fmpz_mod_poly_is_zero(a, ctx))
            {
                fmpz_mod_poly_swap(a, b, ctx);
                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(b, a, H, f,
                                                                    finv, ctx);
                fmpz_mod_poly_add(b, b, a, ctx);
            }
            else
            {
                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(t, a, H, f,
                                                                    finv, ctx);
                fmpz_mod_poly_add(a, b, t, ctx);
                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(t, b, H, f,
                                                                    finv, ctx);
                fmpz_mod_poly_add(b, b, t, ctx);
            }

            d = d/2;
            if (d <= 1)
            {
                if (fmpz_mod_poly_is_zero(a, ctx))
                {
                    fmpz_mod_poly_swap(a, b, ctx);
                    break;

                }

                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(t, xi, H, f,
                                                                    finv, ctx);
                fmpz_mod_poly_swap(xi, t, ctx);

                fmpz_mod_poly_precompute_matrix(H, xi, f, finv, ctx);

                fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(t, a, H, f,
                                                                    finv, ctx);
                fmpz_mod_poly_add(a, t, b, ctx);

                break;
            }

            fmpz_mod_poly_compose_mod(t, xi, xi, f, ctx);
            fmpz_mod_poly_swap(xi, t, ctx);
        }
    }

    fmpz_mat_clear(H);

    FLINT_ASSERT(fmpz_mod_poly_equal(a, a_check, ctx));
#if FLINT_WANT_ASSERT
    fmpz_mod_poly_clear(a_check, ctx);
#endif
}


/*
    Two (non-deterministic) algorithms Algorithm 3 and Theorem 2.3.1 from
        Deterministic Factorization of Polynomials over Finite Fields
            thesis by David Marquis
*/
static void _fmpz_mod_poly_factor_equal_deg_via_trace(
    fmpz_mod_poly_factor_t res,
    const fmpz_mod_poly_t ff,
    slong d,
    const fmpz_mod_poly_t frob,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j, k;
    const fmpz * p = fmpz_mod_ctx_modulus(ctx);
    fmpz_t halfp;   /* floor((p-1)/2) */
    slong n = fmpz_mod_poly_degree(ff, ctx);
    slong r = n/d;
    flint_bitcnt_t r_cutoff;
    fmpz_mod_poly_t t, tq, tr, finv;
    slong Qlen, Qalloc;
    queue_struct * Q;
    split_struct S[FLINT_BITS + 1];
    fmpz_mod_berlekamp_massey_t bma;
    flint_rand_t state;

    FLINT_ASSERT(r > 1 && d > 1);

    r_cutoff = 8*FLINT_BIT_COUNT(d + 1) + fmpz_bits(p);

    res->num = 0;

    flint_randinit(state);

    fmpz_init(halfp);
    fmpz_sub_ui(halfp, p, 1);
    fmpz_fdiv_q_2exp(halfp, halfp, 1);

    fmpz_mod_berlekamp_massey_init(bma, ctx);

    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_init(S[i].f, ctx);
        fmpz_mod_poly_init(S[i].xp, ctx);
        fmpz_mod_poly_init(S[i].a, ctx);
        fmpz_mod_poly_init(S[i].g, ctx);
    }

    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_poly_init(tq, ctx);
    fmpz_mod_poly_init(tr, ctx);
    fmpz_mod_poly_init(finv, ctx);

    Qalloc = 1;
    Q = FLINT_ARRAY_ALLOC(Qalloc, queue_struct);
    for (j = 0; j < Qalloc; j++)
    {
        fmpz_mod_poly_init(Q[j].f, ctx);
        fmpz_mod_poly_init(Q[j].xp, ctx);
    }

    fmpz_mod_poly_set(Q[0].f, ff, ctx);
    fmpz_mod_poly_set(Q[0].xp, frob, ctx);
    Qlen = 1;

    /* init done */

next_queued:

    /* try to factor Q[end].f */
    Qlen--;
    if (Qlen < 0)
        goto cleanup;

    n = fmpz_mod_poly_degree(Q[Qlen].f, ctx);
    r = n/d;

    FLINT_ASSERT(n == r*d);
    FLINT_ASSERT(r > 1);

    fmpz_mod_poly_swap(S->f, Q[Qlen].f, ctx);
    fmpz_mod_poly_swap(S->xp, Q[Qlen].xp, ctx);
    fmpz_mod_poly_reverse(finv, S->f, S->f->length, ctx);
    fmpz_mod_poly_inv_series_newton(finv, finv, S->f->length, ctx);

#if FLINT_WANT_ASSERT
    fmpz_mod_poly_powmod_x_fmpz_preinv(tq, p, S->f, finv, ctx);
    FLINT_ASSERT(fmpz_mod_poly_equal(S->xp, tq, ctx));
#endif

next_alpha:

    fmpz_mod_poly_randtest(t, state, S->f->length - 1, ctx);
    _compute_trace(S->a, t, d, S->xp, S->f, finv, ctx, tq, tr);
    if (fmpz_mod_poly_degree(S->a, ctx) < 1)
        goto next_alpha;

    if (r > r_cutoff)
    {
        if (fmpz_is_zero(halfp))
            fmpz_mod_poly_set(tr, S->a, ctx);
        else
            fmpz_mod_poly_powmod_fmpz_binexp_preinv(tr, S->a, halfp,
                                                              S->f, finv, ctx);
        fmpz_mod_poly_sub_si(tr, tr, 1, ctx);
        fmpz_mod_poly_gcd(t, tr, S->f, ctx);
        if (t->length <= 1 || t->length >= S->f->length)
            goto next_alpha;

        _add_split(res, &Q, &Qlen, &Qalloc, S->f, t, d, S->xp, ctx, tr);
        goto next_queued;
    }

    /*
        If r is not too big, try the "deterministic simplification" to root
        finding. Instead of constructing h(y) = Res_x(f(x), a(x) - y), find
        probablistically a non-trivial divisor g of h, which is good enough.
        If p is large, g^d = h is likely.
    */
    fmpz_mod_berlekamp_massey_start_over(bma, ctx);
    fmpz_mod_poly_one(t, ctx);
    fmpz_mod_berlekamp_massey_add_point_ui(bma, 1, ctx);
    for (i = 1; i < 2*r; i++)
    {
        fmpz_mod_poly_mulmod_preinv(tq, t, S->a, S->f, finv, ctx);
        fmpz_mod_poly_swap(t, tq, ctx);
        FLINT_ASSERT(t->length > 0);
        fmpz_mod_berlekamp_massey_add_point(bma, t->coeffs + 0, ctx);
    }

    /* something went wrong if V does not kill the whole sequence */
    fmpz_mod_berlekamp_massey_reduce(bma, ctx);
    FLINT_ASSERT(bma->R1->length < bma->V1->length);

    fmpz_mod_poly_make_monic(S->g, bma->V1, ctx);

    /* S->g should be squarefree and factor into linears */
    FLINT_ASSERT(1 <= fmpz_mod_poly_degree(S->g, ctx));
    FLINT_ASSERT(fmpz_mod_poly_degree(S->g, ctx) <= r);
#if FLINT_WANT_ASSERT
    fmpz_mod_poly_gen(t, ctx);
    fmpz_mod_poly_powmod_fmpz_binexp(tq, t, p, S->g, ctx);
    fmpz_mod_poly_sub(tq, tq, t, ctx);
    fmpz_mod_poly_gcd(t, tq, S->g, ctx);
    FLINT_ASSERT(fmpz_mod_poly_equal(t, S->g, ctx));
#endif

    /* split S->g into linears and split S->f along the way */
    k = 1;
    while (k > 0)
    {
        k--;
        FLINT_ASSERT(k + 1 < FLINT_BITS);
        FLINT_ASSERT(fmpz_mod_poly_degree(S[k].g, ctx) >= 1);

        if (fmpz_mod_poly_degree(S[k].g, ctx) < 2)
        {
            if (fmpz_mod_poly_degree(S[k].f, ctx) == d)
            {
                fmpz_mod_poly_factor_fit_length(res, res->num + 1, ctx);
                res->exp[res->num] = 1;
                fmpz_mod_poly_set(res->poly + res->num, S[k].f, ctx);
                res->num++;
            }
            else if (fmpz_mod_poly_degree(S[k].f, ctx) > d)
            {
                FLINT_ASSERT(S[k].g->length == 2);
                fmpz_mod_poly_scalar_mul_fmpz(tr, S[k].a, S[k].g->coeffs + 1, ctx);
                fmpz_mod_poly_add_fmpz(tr, S[k].a, S[k].g->coeffs + 0, ctx);
                fmpz_mod_poly_gcd(t, tr, S[k].f, ctx);
                _add_split(res, &Q, &Qlen, &Qalloc, S[k].f, t, d, S[k].xp, ctx, tr);
            }
            continue;
        }

        if (fmpz_mod_poly_is_zero(finv, ctx))
        {
            fmpz_mod_poly_reverse(finv, S[k].f, S[k].f->length, ctx);
            fmpz_mod_poly_inv_series_newton(finv, finv, S[k].f->length, ctx);
        }

        fmpz_mod_poly_swap(t, S[k].g, ctx);
        _fmpz_mod_poly_split_rabin(S[k+0].g, S[k+1].g, t, halfp,
                                                           tq, tr, state, ctx);
        fmpz_mod_poly_compose_mod_brent_kung_preinv(tr, S[k+1].g, S[k].a,
                                                            S[k].f, finv, ctx);
        fmpz_mod_poly_gcd(t, tr, S[k].f, ctx);
        fmpz_mod_poly_divrem(tq, tr, S[k].f, t, ctx);
        fmpz_mod_poly_swap(S[k+0].f, tq, ctx);
        fmpz_mod_poly_swap(S[k+1].f, t, ctx);
        fmpz_mod_poly_zero(finv, ctx); /* finv no longer inv(S[0].f) */

        fmpz_mod_poly_divrem(tq, S[k+1].a, S[k].a, S[k+1].f, ctx);
        fmpz_mod_poly_divrem(tq, tr, S[k].a, S[k+0].f, ctx);
        fmpz_mod_poly_swap(S[k+0].a, tr, ctx);

        fmpz_mod_poly_divrem(tq, S[k+1].xp, S[k+0].xp, S[k+1].f, ctx);
        fmpz_mod_poly_divrem(tq, tr, S[k+0].xp, S[k+0].f, ctx);
        fmpz_mod_poly_swap(S[k+0].xp, tr, ctx);

        k += 2;
    }

    goto next_queued;

cleanup:

    fmpz_mod_berlekamp_massey_clear(bma, ctx);

    for (i = 0; i < Qalloc; i++)
    {
        fmpz_mod_poly_clear(Q[i].f, ctx);
        fmpz_mod_poly_clear(Q[i].xp, ctx);
    }
    flint_free(Q);

    for (i = 0; i <= FLINT_BITS; i++)
    {
        fmpz_mod_poly_clear(S[i].f, ctx);
        fmpz_mod_poly_clear(S[i].xp, ctx);
        fmpz_mod_poly_clear(S[i].a, ctx);
        fmpz_mod_poly_clear(S[i].g, ctx);
    }

    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_poly_clear(tq, ctx);
    fmpz_mod_poly_clear(tr, ctx);
    fmpz_mod_poly_clear(finv, ctx);

    flint_randclear(state);
    fmpz_clear(halfp);

    return;
}

void fmpz_mod_poly_factor_equal_deg_with_frob(
    fmpz_mod_poly_factor_t factors,
    const fmpz_mod_poly_t f,
    slong d,
    const fmpz_mod_poly_t frob,
    const fmpz_mod_ctx_t ctx)
{
    slong r = fmpz_mod_poly_degree(f, ctx)/d;

    FLINT_ASSERT(fmpz_mod_poly_degree(f, ctx) == r*d);

    if (r == 1)
    {
        factors->num = 0;
        fmpz_mod_poly_factor_insert(factors, f, 1, ctx);
    }
    else if (d == 1)
    {
        fmpz_mod_poly_roots(factors, f, 0, ctx);
    }
    else
    {
        _fmpz_mod_poly_factor_equal_deg_via_trace(factors, f, d, frob, ctx);
    }
}

void fmpz_mod_poly_factor_equal_deg(
    fmpz_mod_poly_factor_t factors,
    const fmpz_mod_poly_t f,
    slong d,
    const fmpz_mod_ctx_t ctx)
{
    slong r = fmpz_mod_poly_degree(f, ctx)/d;

    FLINT_ASSERT(fmpz_mod_poly_degree(f, ctx) == r*d);

    if (r == 1)
    {
        factors->num = 0;
        fmpz_mod_poly_factor_insert(factors, f, 1, ctx);
    }
    else if (d == 1)
    {
        fmpz_mod_poly_roots(factors, f, 0, ctx);
    }
    else
    {
        fmpz_mod_poly_t xp, t;
        fmpz_mod_poly_init(xp, ctx);
        fmpz_mod_poly_init(t, ctx);
        fmpz_mod_poly_reverse(t, f, f->length, ctx);
        fmpz_mod_poly_inv_series_newton(t, t, f->length, ctx);
        fmpz_mod_poly_powmod_x_fmpz_preinv(xp, fmpz_mod_ctx_modulus(ctx), f, t, ctx);
        fmpz_mod_poly_clear(t, ctx);
        _fmpz_mod_poly_factor_equal_deg_via_trace(factors, f, d, xp, ctx);
        fmpz_mod_poly_clear(xp, ctx);
    }
}
