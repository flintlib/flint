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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include "fmpz_poly.h"

static void 
__fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R, ulong * d, 
                                     const fmpz * A, len_t lenA, 
                                     const fmpz * B, len_t lenB)
{
    if (lenB <= 16 || (lenA > 2 * lenB - 1 && lenA < 128))
    {
        _fmpz_poly_pseudo_divrem_basecase(Q, R, d, A, lenA, B, lenB);
    }
    else
    {
        const len_t n2 = lenB / 2;
        const len_t n1 = lenB - n2;

        const fmpz * d1 = B + n2;
        const fmpz * d2 = B;
        const fmpz * d3 = B + n1;
        const fmpz * d4 = B;
    
        if (lenA <= lenB + n2 - 1)
        {
            fmpz *p1, *r1, *d2q1;
            fmpz *f;

            /*
               Shift A right by n1, zero the bottom n2 - 1 coeffs; call this p1
             */

            p1 = (fmpz *) flint_malloc((lenA - n1) * sizeof(fmpz));
            {
                len_t i;
                flint_mpn_zero((mp_ptr) p1, n2 - 1);
                for (i = n2 - 1; i < lenA - n1; i++)
                    p1[i] = (A + n1)[i];
            }

            /*
               Compute p1 div d3, at most a 2 n2 - 1 by n2 division, leaving 
               lenA - lenB + 1 <= n2 terms in the quotient
             */

            r1 = R + n1;
            _fmpz_poly_pseudo_divrem_divconquer(Q, r1, d, p1, lenA - n1, d3, n2);

            flint_free(p1);

            /*
               Push the relevant {n2 - 1} terms of the remainder to the 
               top of {R, lenA}
             */

            {
                len_t i;
                for (i = n2 - 2; i >= 0; i--)
                    fmpz_swap(R + lenA - (n2 - 1) + i, r1 + i);
                r1 = R + lenA - (n2 - 1);
            }

            /*
               Compute d2q1 = Q d4 of length lenA - n2, which is 
               at most n1 + n2 - 1 terms
             */

            d2q1 = R;
            _fmpz_poly_mul(d2q1, d4, n1, Q, lenA - lenB + 1);

            /*
               Compute R = L^d R', where R' is the terms of A we have not dealt, 
               of which there are at most n1 + n2 - 1; that is, 

               Set R to {A, n1 + n2 - 1} * f + r1 x^n1 - d2q1
             */

            _fmpz_vec_neg(R, R, lenA - n2);
            _fmpz_vec_add(R + n1, R + n1, R + lenA - n2 + 1, lenA - lenB);
            _fmpz_vec_swap(R + lenA - n2, R + 2 * lenA - lenB + 1 - n2, n2 - (lenA - lenB + 1));

            f = R + lenB - 1;
            fmpz_pow_ui(f, B + (lenB - 1), *d);
            _fmpz_vec_scalar_addmul_fmpz(R, A, n1 + n2 - 1, f);
        }
        else if (lenA > 2 * lenB - 1)
        {
            /*
               XXX:  In this case, we expect A to be modifiable
             */

            ulong s1, s2;
            const len_t shift = lenA - 2 * lenB + 1;

            fmpz * q1 = Q + shift;
            fmpz * q2 = Q;
            fmpz * r1 = R;

            fmpz *p1, *t;
            fmpz_t f;

            fmpz_init(f);

            /*
               Shift A right until it is of length 2 lenB - 1, call this p1; 
               zero the bottom lenB - 1 coeffs
             */

            p1 = (fmpz *) flint_malloc((2 * lenB - 1) * sizeof(fmpz));
            {
                len_t i;
                flint_mpn_zero((mp_ptr) p1, lenB - 1);
                for (i = lenB - 1; i < 2*lenB - 1; i++)
                    p1[i] = (A + shift)[i];
            }

            /*
               Set q1 to p1 div B, a 2 lenB - 1 by lenB division, so q1 ends up 
               being at most length lenB; r1 is of length at most lenB - 1
             */

            _fmpz_poly_pseudo_divrem_divconquer(q1, r1, &s1, p1, 2 * lenB - 1, B, lenB);

            flint_free(p1);

            /*
               Compute t = L^s1 a2 + r1 x^shift, of length at most lenA - lenB 
               since r1 is of length at most lenB - 1.  Here a2 is what remains
               of A after the first lenR coefficients are removed
             */

            t = (fmpz *) A;

            fmpz_pow_ui(f, B + (lenB - 1), s1);

            _fmpz_vec_scalar_mul_fmpz(t, A, lenA - lenB, f);
            _fmpz_vec_add(t + shift, t + shift, r1, lenB - 1);

            /*
               Compute q2 = t div B; it is a smaller division than the original 
               since len(t) <= lenA - lenB, and r2 has length at most lenB - 1
             */

            _fmpz_poly_pseudo_divrem_divconquer(q2, R, &s2, t, lenA - lenB, B, lenB);

            /*
               Write out Q = L^s2 q1 x^shift + q2, of length at most 
               lenB + shift.  Note q2 has length at most shift since it is at 
               most an lenA - lenB by lenB division; q1 cannot have length zero
               since we are doing pseudo division
             */

            fmpz_pow_ui(f, B + (lenB - 1), s2);

            _fmpz_vec_scalar_mul_fmpz(q1, q1, lenB, f);

            *d = s1 + s2;

            fmpz_clear(f);
        }
        else  /* n1 + 2 n2 - 1 < lenA <= 2 lenB - 1 */
        {
            fmpz * q1   = Q + n2;
            fmpz * q2   = Q;
            fmpz * r1   = R;
            fmpz * d2q1 = R + (n1 - 1);
            fmpz *p1, *t;
            fmpz_t f;
            ulong s1, s2;

            fmpz_init(f);

            /*
               Set p1 to the top lenA - 2 n2 coeffs of A, clearing the bottom 
               n1 - 1 coeffs
             */

            p1 = (fmpz *) flint_malloc((lenA - 2 * n2) * sizeof(fmpz));
            {
                len_t i;
                flint_mpn_zero((mp_ptr) p1, n1 - 1);
                for (i = n1 - 1; i < lenA - 2 * n2; i++)
                    p1[i] = (A + 2 * n2)[i];
            }

            /*
               Set q1 to p1 div d1, at most a 2 n1 - 1 by n1 division, so q1 ends 
               up being of length at most n1; r1 is of length n1 - 1
             */

            _fmpz_poly_pseudo_divrem_divconquer(q1, r1, &s1, p1, lenA - 2 * n2, d1, n1);

            flint_free(p1);

            /*
               Compute d2q1 = d2q1, of length lenA - lenB

               Note lenA - lenB <= lenB - 1 <= 2 n2 and lenA - (n1 - 1) > 2 n2, 
               so we can store d2q1 in the top 2 n2 coeffs of R
             */

            if (n2 >= lenA - n1 - 2 * n2 + 1)
                _fmpz_poly_mul(d2q1, d2, n2, q1, lenA - (n1 + 2 * n2 - 1));
            else
                _fmpz_poly_mul(d2q1, q1, lenA - (n1 + 2 * n2 - 1), d2, n2);

            /*
               Compute
                   t = L^s1 * (a2 x^{n1 + n2 - 1} + a3) 
                       + r1 x^{2 n2} - d2q1 x^n2
               of length at most lenB + n2 - 1, since r1 is of length at most 
               n1 - 1 and d2q1 is of length at most n1 + n2 - 1
             */

            t = _fmpz_vec_init(n1 + 2 * n2 - 1);

            fmpz_pow_ui(f, B + (lenB - 1), s1);

            _fmpz_vec_scalar_mul_fmpz(t, A, n1 + 2 * n2 - 1, f);
            _fmpz_vec_add(t + 2 * n2, t + 2 * n2, r1, n1 - 1);
            _fmpz_vec_sub(t + n2, t + n2, d2q1, lenA - lenB);

            /*
               Compute q2 = t div B and set R to the remainder, at most a 
               lenB + n2 - 1 by lenB division, so q2 is of length at most n2
             */

            _fmpz_poly_pseudo_divrem_divconquer(q2, R, &s2, t, lenB + n2 - 1, B, lenB);

            _fmpz_vec_clear(t, n1 + 2 * n2 - 1);

            /*
               Write Q = L^s2 q1 x^n2 + q2; note len(q1) is non-zero since 
               we are performing pseudo division
             */

            fmpz_pow_ui(f, B + (lenB - 1), s2);

            _fmpz_vec_scalar_mul_fmpz(q1, q1, lenA - lenB + 1 - n2, f);

            *d = s1 + s2;

            fmpz_clear(f);
        }
    }
}

void 
_fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R, ulong * d, 
                                    const fmpz * A, len_t lenA, 
                                    const fmpz * B, len_t lenB)
{
    if (lenA <= 2 * lenB - 1)
    {
        __fmpz_poly_pseudo_divrem_divconquer(Q, R, d, A, lenA, B, lenB);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        fmpz *S = _fmpz_vec_init(lenA);

        _fmpz_vec_set(S, A, lenA);

        __fmpz_poly_pseudo_divrem_divconquer(Q, R, d, S, lenA, B, lenB);

        _fmpz_vec_clear(S, lenA);
    }
}

void
fmpz_poly_pseudo_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R,
                                   ulong * d, const fmpz_poly_t A,
                                   const fmpz_poly_t B)
{
    len_t lenq, lenr;
    fmpz *q, *r;

    if (B->length == 0)
    {
        printf("Exception (fmpz_poly_pseudo_divrem_divconquer). Division by zero.\n");
        abort();
    }
    if (Q == R)
    {
        printf("Exception (fmpz_poly_pseudo_divrem_divconquer). \n"
               "Output arguments Q and R may not be aliased.\n");
        abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        fmpz_poly_set(R, A);
        *d = 0;
        return;
    }

    lenq = A->length - B->length + 1;
    lenr = A->length;
    if (Q == A || Q == B)
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }
    if (R == A || R == B)
        r = _fmpz_vec_init(lenr);
    else
    {
        fmpz_poly_fit_length(R, lenr);
        r = R->coeffs;
    }
    
    _fmpz_poly_pseudo_divrem_divconquer(q, r, d, A->coeffs, A->length, 
                                                 B->coeffs, B->length);

    lenr = B->length - 1;
    FMPZ_VEC_NORM(r, lenr);

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);
    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc = A->length;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);
}

