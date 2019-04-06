/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "nmod_poly.h"
#include "mpn_extras.h"

/*
    Assuming deg(a) > deg(b) (in particular b could be zero, but a must not be)
    compute a matrix M = [[m11, m12] [m21 m22]] and r1, r2 so that

    [ m11  m12 ] [ a ] = [ r1 ]
    [ m21  m22 ] [ b ]   [ r2 ]

    with

    (1) r1 and r2 are consecutive remainders in the euclidean remainder
        sequence for a, b satsifying 2*deg(r1) >= deg(a) > 2*deg(r2)

    (2) M is a product of [[0 1][1 -qi]] where the qi are the quotients
        obtained in (1)

    No aliasing: r1 is written to a, r2 is written to b.
    All of the moduli of the arguments should be the same and prime.
*/
void nmod_poly_hgcd(
    nmod_poly_t m11, nmod_poly_t m12,
    nmod_poly_t m21, nmod_poly_t m22,
    nmod_poly_t a, nmod_poly_t b)
{
    slong dega = nmod_poly_degree(a);
    nmod_poly_t q, r, t;

    if (nmod_poly_degree(a) <= nmod_poly_degree(b))
    {
        flint_printf("Exception in nmod_poly_hgcd: bad input degrees\n");
        flint_abort();
    }

    nmod_poly_init_mod(q, a->mod);
    nmod_poly_init_mod(r, a->mod);
    nmod_poly_init_mod(t, a->mod);

    nmod_poly_one(m11);
    nmod_poly_zero(m12);
    nmod_poly_zero(m21);
    nmod_poly_one(m22);

    while (dega <= 2*nmod_poly_degree(b))
    {
        nmod_poly_divrem(q, r, a, b);
        nmod_poly_swap(a, b);
        nmod_poly_swap(b, r);

        /* multipliy M by [[0 1] [1 -q]] on the left */
        nmod_poly_mul(t, q, m21);
        nmod_poly_sub(r, m11, t);
        nmod_poly_swap(m11, m21);
        nmod_poly_swap(m21, r);

        nmod_poly_mul(t, q, m22);
        nmod_poly_sub(r, m12, t);
        nmod_poly_swap(m12, m22);
        nmod_poly_swap(m22, r);
    }

    nmod_poly_clear(q);
    nmod_poly_clear(r);
    nmod_poly_clear(t);
}

/*
typedef struct {
    slong npoints;
    nmod_poly_t R0, R1;
    nmod_poly_t V0, V1;
    nmod_poly_t qt, rt;
    nmod_poly_t points;
} nmod_bma_struct;
typedef nmod_bma_struct nmod_bma_t[1];

    n = B->npoints is the number of points a_1, ..., a_n that have been added
    to the sequence. The polynomials A and S are then defined as

        A = x^n
        S = a_1*x^(n-1) + a_2*x^(n-2) + ... + a_n

    We maintain polynomials U0, V0, U1, V1 such that

        U0*A + V0*S = R0   2*deg(R0) >= n
        U1*A + V1*S = R1   2*deg(R1) < n

    where R0 and R1 are consecutive euclidean remainders and U0, V0, U1, V1 are
    the corresponding Bezout coefficients.

    The U0 and U1 are not stored explicitly. The points a_1, ..., a_n are stored
    in B->values, which is used merely as a resizable array.

    The main usage of this function is the rational reconstruction of a series

     a1    a2    a3         -U1
    --- + --- + --- + ... = ---- maybe
     x    x^2   x^3          V1
*/
void nmod_bma_init(nmod_bma_t B, mp_limb_t p)
{
    nmod_t fpctx;
    nmod_init(&fpctx, p);
    nmod_poly_init_mod(B->V0, fpctx);
    nmod_poly_init_mod(B->R0, fpctx);
    nmod_poly_one(B->R0);
    nmod_poly_init_mod(B->V1, fpctx);
    nmod_poly_one(B->V1);
    nmod_poly_init_mod(B->R1, fpctx);
    nmod_poly_init_mod(B->rt, fpctx);
    nmod_poly_init_mod(B->qt, fpctx);
    nmod_poly_init_mod(B->points, fpctx);
    B->npoints = 0;
    B->points->length = 0;
}


void nmod_bma_start_over(nmod_bma_t B)
{
    B->npoints = 0;
    B->points->length = 0;
    nmod_poly_zero(B->V0);
    nmod_poly_one(B->R0);
    nmod_poly_one(B->V1);
    nmod_poly_zero(B->R1);
}

void nmod_bma_clear(nmod_bma_t B)
{
    nmod_poly_clear(B->R0);
    nmod_poly_clear(B->R1);
    nmod_poly_clear(B->V0);
    nmod_poly_clear(B->V1);
    nmod_poly_clear(B->rt);
    nmod_poly_clear(B->qt);
    nmod_poly_clear(B->points);
}

/* setting the prime also starts over */
void nmod_bma_set_prime(nmod_bma_t B, mp_limb_t p)
{
    nmod_t fpctx;
    nmod_init(&fpctx, p);
    nmod_poly_set_mod(B->V0, fpctx);
    nmod_poly_set_mod(B->R0, fpctx);
    nmod_poly_set_mod(B->V1, fpctx);
    nmod_poly_set_mod(B->R1, fpctx);
    nmod_poly_set_mod(B->rt, fpctx);
    nmod_poly_set_mod(B->qt, fpctx);
    nmod_poly_set_mod(B->points, fpctx);
    nmod_bma_start_over(B);
}

void nmod_bma_print(const nmod_bma_t B)
{
    slong i;
    nmod_poly_print_pretty(B->V1, "#");
    for (i = 0; i < B->points->length; i++)
    {
        flint_printf(" %wu", B->points->coeffs[i]);
    }
}

void nmod_bma_add_points(nmod_bma_t B, const mp_limb_t * a, slong count)
{
    slong i;
    slong old_length = B->points->length;
    nmod_poly_fit_length(B->points, old_length + count);
    for (i = 0; i < count; i++)
    {
        B->points->coeffs[old_length + i] = a[i];
    }
    B->points->length = old_length + count;
}

void nmod_bma_add_zeros(nmod_bma_t B, slong count)
{
    slong i;
    slong old_length = B->points->length;
    nmod_poly_fit_length(B->points, old_length + count);
    for (i = 0; i < count; i++)
    {
        B->points->coeffs[old_length + i] = 0;
    }
    B->points->length = old_length + count;
}

void nmod_bma_add_point(nmod_bma_t B, mp_limb_t a)
{
    slong old_length = B->points->length;
    nmod_poly_fit_length(B->points, old_length + 1);
    B->points->coeffs[old_length] = a;
    B->points->length = old_length + 1;
}

/* return 1 if reduction changed the master poly, 0 otherwise */
int nmod_bma_reduce(nmod_bma_t B)
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
    nmod_poly_zero(B->rt);
    for (i = 0; i < queue_len; i++)
    {
        nmod_poly_set_coeff_ui(B->rt, queue_len - i - 1,
                                      B->points->coeffs[queue_lo + i]);
    }
    B->npoints = queue_hi;

    /* Ri = Ri * x^queue_len + Vi*rt */
    nmod_poly_mul(B->qt, B->V0, B->rt);
    nmod_poly_shift_left(B->R0, B->R0, queue_len);
    nmod_poly_add(B->R0, B->R0, B->qt);
    nmod_poly_mul(B->qt, B->V1, B->rt);
    nmod_poly_shift_left(B->R1, B->R1, queue_len);
    nmod_poly_add(B->R1, B->R1, B->qt);

    /* now start reducing R0, R1 */
    if (2*nmod_poly_degree(B->R1) < B->npoints)
    {
        /* already have deg(R1) < n/2 */
        return 0;
    }

    /* one iteration of euclid to get deg(R0) >= n/2 */
    nmod_poly_divrem(B->qt, B->rt, B->R0, B->R1);
    nmod_poly_swap(B->R0, B->R1);
    nmod_poly_swap(B->R1, B->rt);
    nmod_poly_mul(B->rt, B->qt, B->V1);
    nmod_poly_sub(B->qt, B->V0, B->rt);
    nmod_poly_swap(B->V0, B->V1);
    nmod_poly_swap(B->V1, B->qt);

    l = nmod_poly_degree(B->R0);
    FLINT_ASSERT(B->npoints <= 2*l && l < B->npoints);

    k = B->npoints - l;
    FLINT_ASSERT(0 <= k && k <= l);

    /*
        (l - k)/2 is the expected number of required euclidean iterations.
        Either branch is OK anytime. TODO: find cutoff
    */
    if (l - k < 4)
    {
        while (B->npoints <= 2*nmod_poly_degree(B->R1))
        {
            nmod_poly_divrem(B->qt, B->rt, B->R0, B->R1);
            nmod_poly_swap(B->R0, B->R1);
            nmod_poly_swap(B->R1, B->rt);
            nmod_poly_mul(B->rt, B->qt, B->V1);
            nmod_poly_sub(B->qt, B->V0, B->rt);
            nmod_poly_swap(B->V0, B->V1);
            nmod_poly_swap(B->V1, B->qt);
        }
    }
    else
    {
        nmod_poly_t m11, m12, m21, m22, r0, r1;
        nmod_poly_init_mod(m11, B->V1->mod);
        nmod_poly_init_mod(m12, B->V1->mod);
        nmod_poly_init_mod(m21, B->V1->mod);
        nmod_poly_init_mod(m22, B->V1->mod);
        nmod_poly_init_mod(r0, B->V1->mod);
        nmod_poly_init_mod(r1, B->V1->mod);

        nmod_poly_shift_right(r0, B->R0, k);
        nmod_poly_shift_right(r1, B->R1, k);
        nmod_poly_hgcd(m11, m12, m21, m22, r0, r1);

        /* multiply [[V0 R0] [V1 R1]] by M on the left */

        nmod_poly_mul(B->rt, m11, B->V0);
        nmod_poly_mul(B->qt, m12, B->V1);
        nmod_poly_add(r0, B->rt, B->qt);
        nmod_poly_mul(B->rt, m21, B->V0);
        nmod_poly_mul(B->qt, m22, B->V1);
        nmod_poly_add(r1, B->rt, B->qt);
        nmod_poly_swap(B->V0, r0);
        nmod_poly_swap(B->V1, r1);

        nmod_poly_mul(B->rt, m11, B->R0);
        nmod_poly_mul(B->qt, m12, B->R1);
        nmod_poly_add(r0, B->rt, B->qt);
        nmod_poly_mul(B->rt, m21, B->R0);
        nmod_poly_mul(B->qt, m22, B->R1);
        nmod_poly_add(r1, B->rt, B->qt);
        nmod_poly_swap(B->R0, r0);
        nmod_poly_swap(B->R1, r1);

        nmod_poly_clear(m11);
        nmod_poly_clear(m12);
        nmod_poly_clear(m21);
        nmod_poly_clear(m22);
        nmod_poly_clear(r0);
        nmod_poly_clear(r1);
    }

    FLINT_ASSERT(nmod_poly_degree(B->V1) >= 0);
    FLINT_ASSERT(2*nmod_poly_degree(B->V1) <= B->npoints);
    FLINT_ASSERT(2*nmod_poly_degree(B->R0) >= B->npoints);
    FLINT_ASSERT(2*nmod_poly_degree(B->R1) <  B->npoints);

    return 1;
}
