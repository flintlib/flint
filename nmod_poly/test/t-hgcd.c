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

    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "mpn_extras.h"

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

int
main(void)
{
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("hgcd....");
    fflush(stdout);

    /* 
       Find coprime polys, multiply by another poly 
       and check the GCD is that poly 
    */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c, d, c1, d1, s, t;

        mp_ptr M[4];
        len_t lenM[4];
        len_t sgnM;

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

        /* [c1,d1] := sgnM * M^{-1} [a,b] */
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
            printf("FAIL:\n");
            printf("a  = "), nmod_poly_print(a), printf("\n\n");
            printf("b  = "), nmod_poly_print(b), printf("\n\n");
            printf("c  = "), nmod_poly_print(c), printf("\n\n");
            printf("d  = "), nmod_poly_print(d), printf("\n\n");
            printf("c1 = "), nmod_poly_print(c1), printf("\n\n");
            printf("d1 = "), nmod_poly_print(d1), printf("\n\n");
            abort();
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

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}

#undef __mul

