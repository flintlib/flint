/*
    Copyright (C) 2011, 2014 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_poly.h"
#include "gr_poly.h"

/*
    XXX: Currently supports aliasing between {A,a} and {B,b}.
 */

slong _nmod_poly_hgcd(mp_ptr *M, slong *lenM,
                     mp_ptr A, slong *lenA, mp_ptr B, slong *lenB,
                     mp_srcptr a, slong lena, mp_srcptr b, slong lenb,
                     nmod_t mod)
{
    slong sgnM;
    gr_ctx_t ctx;

    _gr_ctx_init_nmod(ctx, &mod);
    GR_MUST_SUCCEED(_gr_poly_hgcd(NULL, &sgnM, (gr_ptr *) M, lenM, A, lenA, B, lenB, a, lena, b, lenb, NMOD_POLY_HGCD_CUTOFF, ctx));

    return sgnM;
}

/*
    Assuming deg(a) > deg(b) (in particular b could be zero, but a must not be)
    compute a matrix M = [[m11, m12] [m21 m22]] and A, B so that

    [ a ] = [ m11  m12 ][ A ]
    [ b ]   [ m21  m22 ][ B ]

    with

    (1) A and B are consecutive remainders in the euclidean remainder
        sequence for a, b satsifying 2*deg(A) >= deg(a) > 2*deg(B)

    (2) M is a product of [[qi 1][1 0]] where the qi are the quotients
        obtained in (1)

    No aliasing. The return is det(M), which is +-1. All of the moduli of the
    arguments should be the same and prime. Here is a reference implementation
    in case something is wrong with _nmod_poly_hgcd (there doesn't seem to be).
*/
slong nmod_poly_hgcd_ref(
    nmod_poly_t m11, nmod_poly_t m12,
    nmod_poly_t m21, nmod_poly_t m22,
    nmod_poly_t A, nmod_poly_t B,
    const nmod_poly_t a, const nmod_poly_t b)
{
    slong sgnM;
    slong dega = nmod_poly_degree(a);
    nmod_poly_t q, r, t;

    if (nmod_poly_degree(a) <= nmod_poly_degree(b))
    {
        flint_throw(FLINT_ERROR, "Exception in nmod_poly_hgcd_ref:"
                                              " Input degrees are invalid.\n");
    }

    nmod_poly_init_mod(q, a->mod);
    nmod_poly_init_mod(r, a->mod);
    nmod_poly_init_mod(t, a->mod);

    nmod_poly_one(m11);
    nmod_poly_zero(m12);
    nmod_poly_zero(m21);
    nmod_poly_one(m22);
    nmod_poly_set(A, a);
    nmod_poly_set(B, b);

    sgnM = 1;
    while (dega <= 2*nmod_poly_degree(B))
    {
        nmod_poly_divrem(q, r, A, B);
        nmod_poly_swap(A, B);
        nmod_poly_swap(B, r);

        /* multipliy M by [[q 1] [1 0]] on the right */
        nmod_poly_mul(t, q, m11);
        nmod_poly_add(r, m12, t);
        nmod_poly_swap(m11, m12);
        nmod_poly_swap(m11, r);

        nmod_poly_mul(t, q, m21);
        nmod_poly_add(r, m22, t);
        nmod_poly_swap(m21, m22);
        nmod_poly_swap(m21, r);

        sgnM = -sgnM;
    }

    nmod_poly_clear(q);
    nmod_poly_clear(r);
    nmod_poly_clear(t);

    return sgnM;
}


slong nmod_poly_hgcd(
    nmod_poly_t m11, nmod_poly_t m12,
    nmod_poly_t m21, nmod_poly_t m22,
    nmod_poly_t A, nmod_poly_t B,
    const nmod_poly_t a, const nmod_poly_t b)
{
    mp_limb_t * M[4];
    slong lenM[4];
    slong sgnM;

    if (nmod_poly_degree(a) <= nmod_poly_degree(b))
    {
        flint_throw(FLINT_ERROR, "Exception in nmod_poly_hgcd:"
                                              " Input degrees are invalid.\n");
    }

    if (nmod_poly_length(b) == 0)
    {
        nmod_poly_one(m11);
        nmod_poly_zero(m12);
        nmod_poly_zero(m21);
        nmod_poly_one(m22);
        nmod_poly_set(A, a);
        nmod_poly_set(B, b);
        return 1;
    }

    /* now nmod_poly_length(a) > nmod_poly_length(b) > 0 */

    nmod_poly_fit_length(m11, nmod_poly_length(a));
    nmod_poly_fit_length(m12, nmod_poly_length(a));
    nmod_poly_fit_length(m21, nmod_poly_length(a));
    nmod_poly_fit_length(m22, nmod_poly_length(a));
    nmod_poly_fit_length(A, nmod_poly_length(a));
    nmod_poly_fit_length(B, nmod_poly_length(a));

    /*
        Looks like _nmod_poly_hgcd produces an M with det(M) = sgnM = +-1 and
            [ a ] = [ m11  m12 ] [ A ]
            [ b ]   [ m21  m22 ] [ B ]
    */

    M[0] = m11->coeffs;
    M[1] = m12->coeffs;
    M[2] = m21->coeffs;
    M[3] = m22->coeffs;
    sgnM = _nmod_poly_hgcd(M, lenM, A->coeffs, &A->length, B->coeffs, &B->length,
                          a->coeffs, a->length, b->coeffs, b->length, A->mod);
    m11->length = lenM[0];
    m12->length = lenM[1];
    m21->length = lenM[2];
    m22->length = lenM[3];

    FLINT_ASSERT(2*nmod_poly_degree(A) >= nmod_poly_degree(a));
    FLINT_ASSERT(nmod_poly_degree(a) > 2*nmod_poly_degree(B));

    return sgnM;
}
