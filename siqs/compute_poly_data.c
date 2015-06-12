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

    Copyright (C) 2015 Nitin Kumar

******************************************************************************/

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#include <stdio.h>
#undef ulong
#define ulong mp_limb_t

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "qsieve.h"
#include "fmpz.h"


/*
    try to compute A0, which will remain fixed, for A = q0 * A0
    goal is to make sure that,
    (i) A is close to it's ideal value which is, sqrt(2 * kn) / M,
        where kn is the number to be factored and M is the size of sieve
   (ii) A should not contain too many small primes factor,
        it will decrease number of relation
  (iii) A should not contain too large primes factor for reason (i)
   (iv) A should have as many factor as possible, to obtain large number of
        polynomial
*/


void qsieve_compute_A0(qs_t qs_inf)
{
    slong i, s = 0, j = 0, k;
    mp_limb_t p, q;
    prime_t * factor_base = qs_inf->factor_base;
    slong small_primes = qs_inf->small_primes;
    mp_limb_t * A_ind = qs_inf->A_ind;
    fmpz_t prod;
    fmpz_t temp;
    fmpz_t temp2;
    fmpz_init(temp);
    fmpz_init(temp2);
    fmpz_init_set_ui(prod, qs_inf->q0);

    for (i = small_primes; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        fmpz_mul_ui(temp, prod, p);

        if (fmpz_cmp(temp, qs_inf->target_A) <= 0)
        {
            fmpz_set(prod, temp);
            s++;
            A_ind[j++] = i;
        }
        else
        {
            for (k = 0; k < j; k++)
            {
                q = factor_base[A_ind[k]].p;
                fmpz_divexact_ui(temp2, temp, q);
                if (fmpz_cmp(temp2, qs_inf->target_A) <= 0)
                {
                    fmpz_set(prod, temp2);
                    A_ind[k] = i;
                    break;
                }
            }
        }
    }

    fmpz_init_set(qs_inf->A, prod);
    fmpz_divexact_ui(qs_inf->A0, prod, qs_inf->q0);
    qs_inf->s = s;

    fmpz_clear(prod);
    fmpz_clear(temp);
    fmpz_clear(temp2);
}

/*
    calculate A = q0 * A0, where q0 are primes immediately following
    prime bound
*/


mp_limb_t qsieve_compute_A(qs_t qs_inf)
{
    mp_limb_t p, kron, pmod, pinv;
    slong i = 0;

    fmpz_t upper_bound;
    fmpz_t temp;
    n_primes_t iter;
    n_primes_init(iter);
    n_primes_jump_after(iter, qs_inf->factor_base[qs_inf->num_primes - 1].p);

    fmpz_init(upper_bound);
    fmpz_init(temp);
    fmpz_fdiv_q_ui(temp, qs_inf->target_A, UWORD(2));
    fmpz_add(upper_bound, qs_inf->target_A, temp);

    /* find first 'q0', such that 'kn' is quadratic residue modulo 'q0' */

    while(1)
    {
        p = n_primes_next(iter);
        pinv = n_preinvert_limb(p);
        pmod = fmpz_fdiv_ui(qs_inf->kn, p);

        if (pmod == 0)
            return p;

        kron = 1;

        while (pmod % 2 == 0)
        {
            if ((pmod % 8) == 3 || (pmod % 8) == 5) kron *= -1;
            pmod /= 2;
        }

        kron *= n_jacobi(pmod, p);

        if (kron == 1)
        {
            qs_inf->q0 = p;
            break;
        }
    }

    qsieve_compute_A0(qs_inf);

    /* collect all primes 'q0' for which 'kn' is quadratic residue modulo 'q0'
       and value of 'A' is at most 1.5 times the optimal value of 'A'
    */

    while (fmpz_cmp(qs_inf->A, upper_bound) <= 0)
    {
        qs_inf->q0_values = flint_realloc(qs_inf->q0_values, (i + 1) * sizeof(mp_limb_t));
        qs_inf->q0_values[i++] = qs_inf->q0;

        while (1)
        {
            p = n_primes_next(iter);
            pinv = n_preinvert_limb(p);
            pmod = fmpz_fdiv_ui(qs_inf->kn, p);

            if (pmod == 0)
                return p;

            kron = 1;

            while (pmod % 2 == 0)
            {
                if ((pmod % 8) == 3 || (pmod % 8) == 5) kron *= -1;
                pmod /= 2;
            }

            kron *= n_jacobi(pmod, p);

            if (kron == 1)
            {
                qs_inf->q0 = p;
                break;
            }

        }

        fmpz_mul_ui(qs_inf->A, qs_inf->A0, qs_inf->q0);
    }

    qs_inf->num_q0 = i;

    n_primes_clear(iter);
    fmpz_clear(upper_bound);
    fmpz_clear(temp);
    return 0;
}

/*
  precompute data associated with A0 which will remain fixed

*/


void qsieve_compute_pre_data(qs_t qs_inf)
{
    slong i, j;
    slong s = qs_inf->s;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t A0 = qs_inf->A0;
    fmpz_t * A0_divp = qs_inf->A0_divp;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    mp_limb_t p, pinv, temp, temp2;

    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;
        fmpz_init(A0_divp[i]);
        fmpz_divexact_ui(A0_divp[i], qs_inf->A0, p);
        temp2 = fmpz_fdiv_ui(A0_divp[i], p);
        temp2 = n_invmod(temp2, p);
        temp2 = n_mulmod2_preinv(temp2, sqrts[A_ind[i]], p, pinv);
        B0_terms[i] = temp2;
    }
}

/*
   initialized the data for first polynomial after, A = q0 * A0
   (q0 are immediate primes following prime bound) are calculated
   and also precompute data to be used for initialization of
   subsequent polynomial
*/

void qsieve_init_poly_first(qs_t qs_inf)
{
    slong i, j;
    slong s = qs_inf->s;
    mp_limb_t q0 = qs_inf->q0;
    mp_limb_t * A_ind = qs_inf->A_ind;
    fmpz_t * B_terms = qs_inf->B_terms;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    fmpz_t * A0_divp = qs_inf->A0_divp;
    fmpz_t * A_divp = qs_inf->A_divp;
    fmpz_t * B = qs_inf->B;
    prime_t * factor_base = qs_inf->factor_base;
    mp_limb_t p, pinv, temp, temp2;
    fmpz_t temp3, temp4;
    fmpz_init(temp3);
    B = flint_malloc((1 << s) * sizeof(fmpz_t));
    fmpz_init(B[0]);
    fmpz_zero(B[0]);
    fmpz_init(temp4);

    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;
        fmpz_set(temp3, A0_divp[i]);
        fmpz_init(A_divp[i]);
        fmpz_mul_ui(A_divp[i], temp3, q0);
        temp2 = n_invmod(q0, p);
        temp2 = n_mulmod2_preinv(temp2, B0_terms[i], p, pinv);

        if (temp2 >= p / 2) temp2 = p - temp2;

        fmpz_init(B_terms[i]);
        fmpz_mul_ui(B_terms[i], A_divp[i], temp2);
        fmpz_add(temp4, B[0], B_terms[i]);
        fmpz_set(B[0], temp4);
    }

    qs_inf->B = B;

    fmpz_clear(temp3);
    fmpz_clear(temp4);
}

/*
    generate all possible values of coefficient 'B', using
    gray-code formula, for current value of 'A'
*/

void qsieve_init_poly_next(qs_t qs_inf)
{
    slong i, j, v;
    slong s = qs_inf->s;
    prime_t * factor_base = qs_inf->factor_base;
    mp_limb_t sign, p, pinv, pmod, mod_inv;
    fmpz_t temp, temp2, temp3;
    fmpz_init(temp);
    fmpz_init(temp2);
    fmpz_init(temp3);

    for (i = 1; i < (1 << s); i++)
    {
        for (v = 0; v < s; v++)
        {
            if (((i >> v) & 1)) break;
        }

        sign = (i >> v) & 2;

        fmpz_mul_ui(temp, qs_inf->B_terms[v], UWORD(2));
        fmpz_init(qs_inf->B[i]);

        if (sign) fmpz_add(qs_inf->B[i], qs_inf->B[i - 1], temp);
        else fmpz_sub(qs_inf->B[i], qs_inf->B[i - 1], temp);

    }

/* add component of prime 'q0' to 'B' values */
/*
    p = qs_inf->q0;
    pinv = n_preinvert_limb(p);
    pmod = fmpz_fdiv_ui(qs_inf->A0, p);
    mod_inv = n_invmod(pmod, p);
    pmod = fmpz_fdiv_ui(qs_inf->kn, p);
    while (pmod % 2 == 0) pmod /= 2;
    pmod = n_sqrtmod(pmod, p);
    pmod = n_mulmod2_preinv(pmod, mod_inv, p, pinv);
    fmpz_mul_ui(temp2, qs_inf->A0, pmod);
    fmpz_neg(temp3, temp2);

    for (i = 0; i < (1 << s); i++)
    {
        fmpz_set(temp, qs_inf->B[i]);
        fmpz_add(qs_inf->B[i], temp3, temp);
    }
*/
    fmpz_clear(temp);
    fmpz_clear(temp2);
    fmpz_clear(temp3);
}

