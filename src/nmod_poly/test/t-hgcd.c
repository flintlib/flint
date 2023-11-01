/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#define __mul(C, lenC, A, lenA, B, lenB)                        \
do {                                                            \
    if ((lenA) != 0 && (lenB) != 0)                             \
    {                                                           \
        if ((lenA) >= (lenB))                                   \
            _nmod_poly_mul((C), (A), (lenA), (B), (lenB), mod); \
        else                                                    \
            _nmod_poly_mul((C), (B), (lenB), (A), (lenA), mod); \
        (lenC) = (lenA) + (lenB) - 1;                           \
    }                                                           \
    else                                                        \
    {                                                           \
        (lenC) = 0;                                             \
    }                                                           \
} while (0)

TEST_FUNCTION_START(nmod_poly_hgcd, state)
{
    slong i, j, result;

    /* check that [c1,d1] := M^{-1} [a,b] assuming that deg(M) = sgnM */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d, c1, d1, s, t;

        mp_ptr M[4];
        slong lenM[4];
        slong sgnM;

        mp_limb_t n = n_randprime(state, FLINT_BITS, 0);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_init(d, n);
        nmod_poly_init(c1, n);
        nmod_poly_init(d1, n);
        nmod_poly_init(s, n);
        nmod_poly_init(t, n);

        do {
            nmod_poly_randtest_not_zero(a, state, n_randint(state, 800) + 1);
            nmod_poly_randtest_not_zero(b, state, n_randint(state, 800) + 1);
        } while (a->length == b->length);

        if (a->length < b->length)
            nmod_poly_swap(a, b);

        M[0] = _nmod_vec_init(a->length);
        M[1] = _nmod_vec_init(a->length);
        M[2] = _nmod_vec_init(a->length);
        M[3] = _nmod_vec_init(a->length);

        nmod_poly_fit_length(c, a->length);
        nmod_poly_fit_length(d, b->length);

        sgnM = _nmod_poly_hgcd(M, lenM,
                        c->coeffs, &(c->length), d->coeffs, &(d->length),
                        a->coeffs, a->length, b->coeffs, b->length, a->mod);

        nmod_poly_fit_length(s, 2 * a->length);
        nmod_poly_fit_length(t, 2 * a->length);

        /* [c1, d1] := M^{-1} [a,b] */
        {
            const nmod_t mod = a->mod;

            MPN_SWAP(M[0], lenM[0], M[3], lenM[3]);
            _nmod_vec_neg(M[1], M[1], lenM[1], mod);
            _nmod_vec_neg(M[2], M[2], lenM[2], mod);

            __mul(s->coeffs, s->length, M[0], lenM[0], a->coeffs, a->length);
            __mul(t->coeffs, t->length, M[1], lenM[1], b->coeffs, b->length);
            nmod_poly_add(c1, s, t);
            __mul(s->coeffs, s->length, M[2], lenM[2], a->coeffs, a->length);
            __mul(t->coeffs, t->length, M[3], lenM[3], b->coeffs, b->length);
            nmod_poly_add(d1, s, t);
        }

        if (sgnM < 0)
        {
            nmod_poly_neg(c1, c1);
            nmod_poly_neg(d1, d1);
        }

        result = (nmod_poly_equal(c, c1) && nmod_poly_equal(d, d1));
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a  = "), nmod_poly_print(a), flint_printf("\n\n");
            flint_printf("b  = "), nmod_poly_print(b), flint_printf("\n\n");
            flint_printf("c  = "), nmod_poly_print(c), flint_printf("\n\n");
            flint_printf("d  = "), nmod_poly_print(d), flint_printf("\n\n");
            flint_printf("c1 = "), nmod_poly_print(c1), flint_printf("\n\n");
            flint_printf("d1 = "), nmod_poly_print(d1), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
        nmod_poly_clear(c1);
        nmod_poly_clear(d1);
        nmod_poly_clear(s);
        nmod_poly_clear(t);

        _nmod_vec_clear(M[0]);
        _nmod_vec_clear(M[1]);
        _nmod_vec_clear(M[2]);
        _nmod_vec_clear(M[3]);
    }

    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        mp_limb_t p;
        slong sgnM, sgnMr;
        nmod_poly_t m11, m12, m21, m22, A, B, a, b;
        nmod_poly_t m11r, m12r, m21r, m22r, Ar, Br;
        nmod_poly_t q, t;

        p = n_randtest_prime(state, 1);
        nmod_poly_init(m11, p);
        nmod_poly_init(m12, p);
        nmod_poly_init(m21, p);
        nmod_poly_init(m22, p);
        nmod_poly_init(A, p);
        nmod_poly_init(B, p);
        nmod_poly_init(m11r, p);
        nmod_poly_init(m12r, p);
        nmod_poly_init(m21r, p);
        nmod_poly_init(m22r, p);
        nmod_poly_init(Ar, p);
        nmod_poly_init(Br, p);
        nmod_poly_init(a, p);
        nmod_poly_init(b, p);
        nmod_poly_init(q, p);
        nmod_poly_init(t, p);

        /* prepare inputs a, b */
        nmod_poly_randtest_monic(a, state, 1 + n_randint(state, 10));
        nmod_poly_scalar_mul_nmod(a, a, 1 + n_randint(state, p - 1));
        nmod_poly_zero(b);
        for (j = 1 + n_randint(state, 100); j >= 0; j--)
        {
            nmod_poly_randtest_monic(q, state, 2 + n_randint(state, 20));
            nmod_poly_scalar_mul_nmod(q, q, 1 + n_randint(state, p - 1));
            nmod_poly_mul(t, a, q);
            nmod_poly_add(b, b, t);
            nmod_poly_swap(a, b);
        }

        sgnMr = nmod_poly_hgcd_ref(m11r, m12r, m21r, m22r, Ar, Br, a, b);
        sgnM  = nmod_poly_hgcd(m11, m12, m21, m22, A, B, a, b);

        if (sgnMr != sgnM
            || !nmod_poly_equal(m11r, m11)
            || !nmod_poly_equal(m12r, m12)
            || !nmod_poly_equal(m21r, m21)
            || !nmod_poly_equal(m22r, m22)
            || !nmod_poly_equal(Ar, A)
            || !nmod_poly_equal(Br, B))
        {
            flint_printf("check reference match\n i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        if (!(2*nmod_poly_degree(A) >= nmod_poly_degree(a)
            && nmod_poly_degree(a) > 2*nmod_poly_degree(B)))
        {
            flint_printf("check degrees\n i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_mul(q, m11, m22);
        nmod_poly_mul(t, m12, m21);
        nmod_poly_sub(t, q, t);
        if (sgnM < 0)
        {
            nmod_poly_neg(t, t);
        }
        if (!nmod_poly_is_one(t))
        {
            flint_printf("check sign of determinant\n i = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(m11);
        nmod_poly_clear(m12);
        nmod_poly_clear(m21);
        nmod_poly_clear(m22);
        nmod_poly_clear(A);
        nmod_poly_clear(B);
        nmod_poly_clear(m11r);
        nmod_poly_clear(m12r);
        nmod_poly_clear(m21r);
        nmod_poly_clear(m22r);
        nmod_poly_clear(Ar);
        nmod_poly_clear(Br);
        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(q);
        nmod_poly_clear(t);
    }

    TEST_FUNCTION_END(state);
}

#undef __mul
