/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPinterpE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Authored 2016 by Daniel S. Roche; US Government work in the public domain

******************************************************************************/

#include <stdio.h>
#include "flint.h"
#include "fmpz_mod_poly.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

/* checks that poly actually annihilates the given sequence. */
int check(const fmpz_mod_poly_t poly, const fmpz* seq, slong len)
{
    fmpz_t sum, temp;
    slong d = fmpz_mod_poly_degree(poly);
    int i, j;

    if (d < 0) return 0;

    fmpz_init(sum);
    fmpz_init(temp);

    for (i=0; i < len-d; ++i)
    {
        fmpz_zero(sum);
        for (j=0; j<d; ++j)
        {
            fmpz_mod_poly_get_coeff_fmpz(temp, poly, j);
            fmpz_addmul(sum, temp, seq+(i+j));
        }
        fmpz_add(sum, sum, seq+(i+d));
        fmpz_mod(sum, sum, &poly->p);
        if (!fmpz_is_zero(sum)) 
        {
            fmpz_clear(sum);
            fmpz_clear(temp);
            return 0;
        }
    }

    fmpz_clear(sum);
    fmpz_clear(temp);
    return 1;
}

int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("minpoly....");
    fflush(stdout);

    /* test random sequences */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz* seq;
        slong len;
        fmpz_mod_poly_t poly1, poly2;
        fmpz_t p;
        int j;

        fmpz_init(p);
        fmpz_randprime(p, state, 100, 0);

        len = n_randtest(state) % UWORD(100);
        seq = _fmpz_vec_init(len);
        for (j=0; j<len; ++j) fmpz_randtest_mod(seq+j, state, p);

        fmpz_mod_poly_init(poly1, p);
        fmpz_mod_poly_init(poly2, p);

        fmpz_mod_poly_minpoly_bm(poly1, seq, len);
        fmpz_mod_poly_minpoly_hgcd(poly2, seq, len);

        if (!check(poly1, seq, len) 
            || fmpz_mod_poly_degree(poly1) > fmpz_mod_poly_degree(poly2))
        {
            flint_printf("FAIL 1:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly1); flint_printf("\n\n");
            abort();
        }

        if (!check(poly2, seq, len) 
            || fmpz_mod_poly_degree(poly2) > fmpz_mod_poly_degree(poly1))
        {
            flint_printf("FAIL 2:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly2); flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_clear(poly1);
        fmpz_mod_poly_clear(poly2);
        fmpz_clear(p);
        _fmpz_vec_clear(seq, len);
    }

    /* test sequences with a known generator */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz* seq;
        slong len, d;
        fmpz_mod_poly_t poly1, poly2, gen, rem;
        fmpz_t p, temp;
        int j, k;

        fmpz_init(p);
        fmpz_init(temp);
        fmpz_randprime(p, state, 100, 0);

        len = n_randtest(state) % UWORD(200) + 2;
        seq = _fmpz_vec_init(len);

        fmpz_mod_poly_init(poly1, p);
        fmpz_mod_poly_init(poly2, p);
        fmpz_mod_poly_init(gen, p);
        fmpz_mod_poly_init(rem, p);
        fmpz_mod_poly_randtest_monic(gen, state, n_randint(state, len/2)+2);
        d = fmpz_mod_poly_degree(gen);
        FLINT_ASSERT (d > 0);

        for (j=0; j<d; ++j) fmpz_randtest_mod(seq+j, state, p);

        for (; j<len; ++j)
        {
            fmpz_zero(seq+j);
            for (k=0; k<d; ++k)
            {
                fmpz_mod_poly_get_coeff_fmpz(temp, gen, k);
                fmpz_submul(seq+j, temp, seq + (j-d+k));
            }
            fmpz_mod(seq+j, seq+j, p);
        }
        FLINT_ASSERT(check(gen, seq, len));

        fmpz_mod_poly_minpoly_bm(poly1, seq, len);

        if (!check(poly1, seq, len))
        {
            flint_printf("FAIL 3:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly1); flint_printf("\n\n");
            fmpz_mod_poly_print(gen);  flint_printf("\n\n");
            abort();
        }

        result = fmpz_mod_poly_degree(poly1) <= fmpz_mod_poly_degree(gen);
        if (result && fmpz_mod_poly_degree(gen) <= len/2)
        {
            fmpz_mod_poly_rem(rem, gen, poly1);
            result = fmpz_mod_poly_is_zero(rem);
        }

        if (!result)
        {
            flint_printf("FAIL 4:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly1); flint_printf("\n\n");
            fmpz_mod_poly_print(gen);  flint_printf("\n\n");
            abort();
        }

        fmpz_mod_poly_minpoly_hgcd(poly2, seq, len);

        if (!fmpz_mod_poly_equal(poly1, poly2))
        {
            flint_printf("FAIL 5:\n");
            _fmpz_vec_print(seq, len); flint_printf("\n\n");
            fmpz_mod_poly_print(poly2); flint_printf("\n\n");
            fmpz_mod_poly_print(gen);  flint_printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(seq, len);
        fmpz_clear(p);
        fmpz_clear(temp);
        fmpz_mod_poly_clear(poly1);
        fmpz_mod_poly_clear(poly2);
        fmpz_mod_poly_clear(gen);
        fmpz_mod_poly_clear(rem);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
