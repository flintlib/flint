/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"

/*
typedef struct {
    slong npoints;
    fmpz_mod_poly_t R0, R1;
    fmpz_mod_poly_t V0, V1;
    fmpz_mod_poly_t qt, rt; temporaries
    fmpz_mod_poly_t points;
} fmpz_mod_berlekamp_massey_struct;
typedef fmpz_mod_berlekamp_massey_struct nmod_berlekamp_massey_t[1];

    n = B->npoints is the number of points a_1, ..., a_n that have been added
    to the sequence. The polynomials A and S are then defined as

        A = x^n
        S = a_1*x^(n-1) + a_2*x^(n-2) + ... + a_n

    We maintain polynomials U0, V0, U1, V1 such that

        U0*A + V0*S = R0   deg(R0) >= n/2
        U1*A + V1*S = R1   deg(R1) < n/2

    where R0 and R1 are consecutive euclidean remainders and U0, V0, U1, V1 are
    the corresponding Bezout coefficients. Note that
        deg(U1) < deg(V1) = deg(A) - deg(R0) <= n/2

    The U0 and U1 are not stored explicitly. The points a_1, ..., a_n are stored
    in B->points, which is used merely as a resizable array.

    The main usage of this function is the rational reconstruction of a series

     a1    a2    a3         -U1
    --- + --- + --- + ... = ---- maybe
     x    x^2   x^3          V1

    It can be seen that

     a1    a2        an   -U1     R1
    --- + --- + ... --- = --- + -------
     x    x^2       x^n    V1   V1*x^n

    Thus the error is O(1/x^(n+1)) iff deg(R1) < deg(V1).
*/
void fmpz_mod_berlekamp_massey_init(
    fmpz_mod_berlekamp_massey_t B,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_init(B->V0, ctx);
    fmpz_mod_poly_init(B->R0, ctx);
    fmpz_mod_poly_set_ui(B->R0, 1, ctx);
    fmpz_mod_poly_init(B->V1, ctx);
    fmpz_mod_poly_set_ui(B->V1, 1, ctx);
    fmpz_mod_poly_init(B->R1, ctx);
    fmpz_mod_poly_init(B->rt, ctx);
    fmpz_mod_poly_init(B->qt, ctx);
    fmpz_mod_poly_init(B->points, ctx);
    B->npoints = 0;
    B->points->length = 0;
}

void fmpz_mod_berlekamp_massey_start_over(
    fmpz_mod_berlekamp_massey_t B,
    const fmpz_mod_ctx_t ctx)
{
    B->npoints = 0;
    B->points->length = 0;
    fmpz_mod_poly_zero(B->V0, ctx);
    fmpz_mod_poly_set_ui(B->R0, 1, ctx);
    fmpz_mod_poly_set_ui(B->V1, 1, ctx);
    fmpz_mod_poly_zero(B->R1, ctx);
}

void fmpz_mod_berlekamp_massey_clear(
    fmpz_mod_berlekamp_massey_t B,
    const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_clear(B->R0, ctx);
    fmpz_mod_poly_clear(B->R1, ctx);
    fmpz_mod_poly_clear(B->V0, ctx);
    fmpz_mod_poly_clear(B->V1, ctx);
    fmpz_mod_poly_clear(B->rt, ctx);
    fmpz_mod_poly_clear(B->qt, ctx);
    fmpz_mod_poly_clear(B->points, ctx);
}

void fmpz_mod_berlekamp_massey_print(
    const fmpz_mod_berlekamp_massey_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_mod_poly_print_pretty(B->V1, "#", ctx);
    flint_printf(",");
    for (i = 0; i < B->points->length; i++)
    {
        flint_printf(" ");
        fmpz_print(B->points->coeffs + i);
    }
}

void fmpz_mod_berlekamp_massey_add_points(
    fmpz_mod_berlekamp_massey_t B,
    const fmpz * a,
    slong count,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    slong old_length = B->points->length;
    fmpz_mod_poly_fit_length(B->points, old_length + count, ctx);
    for (i = 0; i < count; i++)
    {
        FLINT_ASSERT(fmpz_mod_is_canonical(a + i, ctx));
        fmpz_set(B->points->coeffs + old_length + i, a + i);
    }
    B->points->length = old_length + count;
}

void fmpz_mod_berlekamp_massey_add_zeros(
    fmpz_mod_berlekamp_massey_t B,
    slong count,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    slong old_length = B->points->length;
    fmpz_mod_poly_fit_length(B->points, old_length + count, ctx);
    for (i = 0; i < count; i++)
    {
        fmpz_zero(B->points->coeffs + old_length + i);
    }
    B->points->length = old_length + count;
}

void fmpz_mod_berlekamp_massey_add_point(
    fmpz_mod_berlekamp_massey_t B,
    const fmpz_t a,
    const fmpz_mod_ctx_t ctx)
{
    slong old_length = B->points->length;
    fmpz_mod_poly_fit_length(B->points, old_length + 1, ctx);
    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
    fmpz_set(B->points->coeffs + old_length, a);
    B->points->length = old_length + 1;
}

void fmpz_mod_berlekamp_massey_add_point_ui(
    fmpz_mod_berlekamp_massey_t B,
    ulong a,
    const fmpz_mod_ctx_t ctx)
{
    slong old_length = B->points->length;
    fmpz_mod_poly_fit_length(B->points, old_length + 1, ctx);
    FLINT_ASSERT(fmpz_cmp_ui(fmpz_mod_ctx_modulus(ctx), a) > 0);
    fmpz_set_ui(B->points->coeffs + old_length, a);
    B->points->length = old_length + 1;
}

/* return 1 if reduction changed the master poly, 0 otherwise */
int fmpz_mod_berlekamp_massey_reduce(
    fmpz_mod_berlekamp_massey_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i, l, k, queue_len, queue_lo, queue_hi;

    /*
        the points in B->points->coeffs[j] for queue_lo <= j < queue_hi need
        to be added to the internal polynomials.
        These are first reversed into rt. deg(rt) < queue_len.
    */
    queue_lo = B->npoints;
    queue_hi = B->points->length;
    queue_len = queue_hi - queue_lo;
    FLINT_ASSERT(queue_len >= 0);
    fmpz_mod_poly_zero(B->rt, ctx);
    for (i = 0; i < queue_len; i++)
    {
        fmpz_mod_poly_set_coeff_fmpz(B->rt, queue_len - i - 1,
                                        B->points->coeffs + queue_lo + i, ctx);
    }
    B->npoints = queue_hi;

    /* Ri = Ri * x^queue_len + Vi*rt */
    fmpz_mod_poly_mul(B->qt, B->V0, B->rt, ctx);
    fmpz_mod_poly_shift_left(B->R0, B->R0, queue_len, ctx);
    fmpz_mod_poly_add(B->R0, B->R0, B->qt, ctx);
    fmpz_mod_poly_mul(B->qt, B->V1, B->rt, ctx);
    fmpz_mod_poly_shift_left(B->R1, B->R1, queue_len, ctx);
    fmpz_mod_poly_add(B->R1, B->R1, B->qt, ctx);

    /* now start reducing R0, R1 */
    if (2*fmpz_mod_poly_degree(B->R1, ctx) < B->npoints)
    {
        /* already have deg(R1) < B->npoints/2 */
        return 0;
    }

    /* one iteration of euclid to get deg(R0) >= B->npoints/2 */
    fmpz_mod_poly_divrem(B->qt, B->rt, B->R0, B->R1, ctx);
    fmpz_mod_poly_swap(B->R0, B->R1, ctx);
    fmpz_mod_poly_swap(B->R1, B->rt, ctx);
    fmpz_mod_poly_mul(B->rt, B->qt, B->V1, ctx);
    fmpz_mod_poly_sub(B->qt, B->V0, B->rt, ctx);
    fmpz_mod_poly_swap(B->V0, B->V1, ctx);
    fmpz_mod_poly_swap(B->V1, B->qt, ctx);

    l = fmpz_mod_poly_degree(B->R0, ctx);
    FLINT_ASSERT(B->npoints <= 2*l && l < B->npoints);

    k = B->npoints - l;
    FLINT_ASSERT(0 <= k && k <= l);

    /*
        (l - k)/2 is the expected number of required euclidean iterations.
        Either branch is OK anytime. TODO: find cutoff
    */
    if (l - k < 10)
    {
        while (B->npoints <= 2*fmpz_mod_poly_degree(B->R1, ctx))
        {
            fmpz_mod_poly_divrem(B->qt, B->rt, B->R0, B->R1, ctx);
            fmpz_mod_poly_swap(B->R0, B->R1, ctx);
            fmpz_mod_poly_swap(B->R1, B->rt, ctx);
            fmpz_mod_poly_mul(B->rt, B->qt, B->V1, ctx);
            fmpz_mod_poly_sub(B->qt, B->V0, B->rt, ctx);
            fmpz_mod_poly_swap(B->V0, B->V1, ctx);
            fmpz_mod_poly_swap(B->V1, B->qt, ctx);
        }
    }
    else
    {
        /* TODO: get hgcd working in this branch */
        while (B->npoints <= 2*fmpz_mod_poly_degree(B->R1, ctx))
        {
            fmpz_mod_poly_divrem(B->qt, B->rt, B->R0, B->R1, ctx);
            fmpz_mod_poly_swap(B->R0, B->R1, ctx);
            fmpz_mod_poly_swap(B->R1, B->rt, ctx);
            fmpz_mod_poly_mul(B->rt, B->qt, B->V1, ctx);
            fmpz_mod_poly_sub(B->qt, B->V0, B->rt, ctx);
            fmpz_mod_poly_swap(B->V0, B->V1, ctx);
            fmpz_mod_poly_swap(B->V1, B->qt, ctx);
        }
    }

    FLINT_ASSERT(fmpz_mod_poly_degree(B->V1, ctx) >= 0);
    FLINT_ASSERT(2*fmpz_mod_poly_degree(B->V1, ctx) <= B->npoints);
    FLINT_ASSERT(2*fmpz_mod_poly_degree(B->R0, ctx) >= B->npoints);
    FLINT_ASSERT(2*fmpz_mod_poly_degree(B->R1, ctx) <  B->npoints);

    return 1;
}
