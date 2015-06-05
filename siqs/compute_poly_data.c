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


/***
   New field added to qs_t struct,

   qs_inf->q0 : factor of A, immediately following largest factor-base prime
   qs_inf->A0 : value of coeff A excluding non-factor-base prime
   qs_inf->A_divp : store (A0 / p), for factor p of A0
   qs_inf->B0_terms : B0_terms[i] = (kn^(1 / 2) * (A0 / p)^(-1)) mod p

***/

void qsieve_compute_A(qs_t qs_inf)
{
    n_primes_t iter;
    n_primes_init(iter);
    n_primes_jump_after(iter, qs_inf->factor_base[qs_inf->num_primes - 1]);
    qs_inf->A0 = qsieve_compute_A0();     /* calculate A0, to add */

    qsieve_compute_pre_data(qs_inf);     /* pre-compute data for A0 which will be fixed */

    do
    {
        qs_inf->q0 = n_primes_next(iter);         /* iterate through primes immediately following largest factor base primes */
        qsieve_init_poly_first(qs_inf);
        qsieve_init_poly_next(qs_inf);

    }while(); /* we have required number of relation, condition to put */

    n_primes_clear(iter);
}

void qsieve_compute_pre_data(qs_t qs_inf)   /* pre-computes all the data associated with A0 */
{
    slong i, j;
    slong s = qs_inf->s;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t A0 = qs_inf->A0;
    mp_limb_t * A0_divp = qs_inf->A0_divp;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    mp_limb_t p, pinv, temp, temp2;

    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;
        A0_divp[i] = (temp2 = A0 / p);
        temp2 = n_invmod(temp2, p);
        temp2 = n_mulmod2_preinv(temp2, sqrts[A_ind[i]], p, pinv);
        B0_terms[i] = temp2;
    }

}

void qsieve_init_poly_first(qs_t qs_inf)        /* initialize value for first polynomial and */
{                                               /* precompute data for subsequent polynomial */
    slong i, j;
    slong s = qs_inf->s;
    mp_limb_t A0 = qs_inf->A0;
    mp_limb_t A = qs_inf->A;
    mp_limb_t q0 = qs_inf->q0;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * A_inv = qs_inf->A_inv;
    mp_limb_t * B_terms = qs_inf->B_terms;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    mp_limb_t * A0_modp = qs_inf->A0_modp;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    mp_limb_t p, pinv, temp, temp2, B = 0;

    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;
        temp = A0_divp[i] * q0;
        A_modp[i] = n_mod2_preinv(temp, p, pinv);
        temp2 = n_invmod(q0, p);
        temp2 = n_mulmod2_preinv(temp2, B0_terms[i], p, pinv);

        if (temp2 >= p / 2) temp2 = p - temp2;

        B_terms[i] = temp * temp2;
        B += B_terms[i];
    }

    qs_inf->B = B;

    for (i = 0; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        pinv = factor_base[i].pinv;
        A_inv[i] = n_invmod(n_mod2_preinv(A, p, pinv), p);

        for (j = 0; j < s; j++)
        {
            temp = n_mod2_preinv(B_terms[j], p, pinv);
            temp = n_mulmod2_preinv(temp, A_inv[i], p, pinv);
            temp *= 2;
            if (temp >= p) temp -= p;
            A_inv2B[j][i] = temp;
        }

        temp = n_mod2_preinv(B, p, pinv);
        temp = sqrts[i] + p - temp;
        temp *= A_inv[i];
        soln1[i] = n_mod2_preinv(temp, p, pinv);
        temp = p - sqrts[i];
        if (temp == p) temp -= p;
        temp = n_mulmod2_preinv(temp, A_inv[i], p, pinv);
        temp *= 2;
        if (temp >= p) temp -= p;
        soln2[i] = temp + soln1[i];
        if (soln2[i] >= p) soln2[i] -= p;
    }

    /* call for sieving */
}

void qsieve_init_poly_next(qs_t qs_inf)
{
    slong i, v;
    slong s = qs_inf->s;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    prime_t * factor_base = qs_inf->factor_base;
    mp_limb_t sign, p, pinv;

    for (i = 1; i < (1 << (s - 1)); i++)              /* iterate through all the possible polynomial */
    {
        for (v = 0; v < s; v++)
        {
            if (((i >> v) & 1) break;
        }

        sign = (i >> v) & 2;

        if (sign) qs_inf->B += (2 * qs_inf->B_terms[v]);
        else qs_inf->B -= (2 * qs_inf->B_terms[v]);

        for (i = 0; i < qs_inf->num_primes; i++)
        {
            p = factor_base[i].p;
            pinv = factor_base[i].pinv;

            if (sign)
            {
                soln1[i] += A_inv2B[v][i];
                soln2[i] += A_inv2B[v][i];
            }
            else
            {
                soln1[i] += (p - A_inv2B[v][i]);
                soln2[i] += (p - A_inv2B[v][i]);
            }

            if (soln1[i] >= p) soln1[i] -= p;
            if (soln2[i] >= p) soln2[i] -= p;
        }

        qsieve_compute_C(qs_inf);

        /* call for sieving */

    }
}
