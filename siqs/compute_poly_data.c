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

    New field added to qs_t struct,

    mp_limb_t q0         : factor of A, immediately following largest factor-base prime
    mp_limb_t A0         : value of coefficient A excluding non-factor-base prime
    mp_limb_t * A_divp   : store (A0 / p), for factor p of A0
    mp_limb_t * B0_terms : B0_terms[i] = (kn^(1 / 2) * (A0 / p)^(-1)) mod p
                           where root and inverse are taken modulo p

*/


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
    mp_limb_t prod = 1, p, q, temp, temp2;
    prime_t * factor_base = qs_inf->factor_base;
    slong small_primes = qs_inf->small_primes;
    mp_limb_t * A_ind = qs_inf->A_ind;
    prod = qs_inf->q0;

    for (i = small_primes; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        temp = p * prod;

        if (temp <= qs_inf->target_A / 2 || temp <= qs_inf->target_A * 2)
        {
            prod = temp;
            s++;
            A_ind[j++] = i;
        }
        else
        {
            for (k = j - 1; k >= 0; k--)
            {
                q = factor_base[A_ind[k]].p;
                temp2 = temp / q;
                if (temp2 <= qs_inf->target_A * 2)
                {
                    A_ind[k] = i;
                    prod = temp2;
                    break;
                }
            }
        }
    }

    qs_inf->A0 = prod;
    qs_inf->s = s;
}

/*
    calculate A = q0 * A0, where q0 are primes immediately following
    prime bound
*/
void qsieve_compute_A(qs_t qs_inf)
{
    n_primes_t iter;
    n_primes_init(iter);
    n_primes_jump_after(iter, qs_inf->factor_base[qs_inf->num_primes - 1].p);
    qs_inf->q0 = n_primes_next(iter);
    qsieve_compute_A0(qs_inf);
    qsieve_compute_pre_data(qs_inf);
    qs_inf->A = qs_inf->q0 * qs_inf->A0;

    while (qs_inf->A >= qs_inf->target_A / 2 && qs_inf <= 2 * qs_inf->target_A)
    {
       qsieve_init_poly_first(qs_inf);
       qsieve_init_poly_next(qs_inf);
       qs_inf->q0 = n_primes_next(iter);
       qs_inf->A = qs_inf->q0 * qs_inf->A0;
    }

    n_primes_clear(iter);
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
    mp_limb_t A0 = qs_inf->A0;
    mp_limb_t A = qs_inf->A;
    mp_limb_t q0 = qs_inf->q0;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * A_inv = qs_inf->A_inv;
    mp_limb_t * B_terms = qs_inf->B_terms;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    mp_limb_t * A0_divp = qs_inf->A0_divp;
    mp_limb_t * A_modp = qs_inf->A_modp;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    prime_t * factor_base = qs_inf->factor_base;
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
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

/*
    iterate through possible values of coefficient 'B', using
    gray-code formula
*/
void qsieve_init_poly_next(qs_t qs_inf)
{
    slong i, j, v;
    slong s = qs_inf->s;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    prime_t * factor_base = qs_inf->factor_base;
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
    mp_limb_t sign, p, pinv;

    for (i = 1; i < (1 << (s - 1)); i++)
    {
        for (v = 0; v < s; v++)
        {
            if (((i >> v) & 1)) break;
        }

        sign = (i >> v) & 2;

        if (sign) qs_inf->B += (2 * qs_inf->B_terms[v]);
        else qs_inf->B -= (2 * qs_inf->B_terms[v]);

        for (j = 0; j < qs_inf->num_primes; j++)
        {
            p = factor_base[j].p;
            pinv = factor_base[j].pinv;

            if (sign)
            {
                soln1[j] += A_inv2B[v][j];
                soln2[j] += A_inv2B[v][j];
            }
            else
            {
                soln1[j] += (p - A_inv2B[v][j]);
                soln2[j] += (p - A_inv2B[v][j]);
            }

            if (soln1[j] >= p) soln1[j] -= p;
            if (soln2[j] >= p) soln2[j] -= p;
        }

        qsieve_compute_C(qs_inf);

        /* call for sieving */

    }
}

/*
     calculate coefficient 'C', once we have coefficient
     'A' and 'B'
*/

void qsieve_compute_C(qs_t qs_inf)
{
   mp_limb_t A = qs_inf->A;
   mp_limb_t B = qs_inf->B;

   if ((mp_limb_signed_t) B < WORD(0)) B = -B;
   fmpz_set_ui(qs_inf->C, B);
   fmpz_mul_ui(qs_inf->C, qs_inf->C, B);
   fmpz_sub(qs_inf->C, qs_inf->C, qs_inf->kn);
   fmpz_divexact_ui(qs_inf->C, qs_inf->C, A);
}
