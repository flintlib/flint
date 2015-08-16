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


mp_limb_t qsieve_next_A0(qs_t qs_inf)
{
    slong i, j, mid;
    slong s = qs_inf->s;
    slong low = qs_inf->low;
    slong high = qs_inf->high;
    slong span = qs_inf->span;
    slong h = qs_inf->h;
    slong m = qs_inf->m;
    mp_limb_t ret = 1;
    mp_limb_t * curr_subset = qs_inf->curr_subset;
    mp_limb_t * A_ind = qs_inf->A_ind;
    prime_t * factor_base = qs_inf->factor_base;
    fmpz_t prod, temp, lower_bound, upper_bound;
    fmpz_init(prod);
    fmpz_init(temp);
    fmpz_init_set(lower_bound, qs_inf->low_bound);
    fmpz_init_set(upper_bound, qs_inf->upp_bound);

    /* generate next value of 'A0' */

    if (s <= 3)
    {
        if (curr_subset[0] != span - s + 1)
        {
            h = (m >= span - h) ? h + 1 : 1;
            m = curr_subset[s - h];

            for (j = 1; j <= h; j++)
                curr_subset[s + j - h - 1] = m + j;

            fmpz_set_ui(prod, UWORD(1));

            for (j = 0; j < s; j++)
                fmpz_mul_ui(prod, prod, factor_base[curr_subset[j] + low - 1].p);

            for (j = 0; j < s; j++)
                A_ind[j] = curr_subset[j] + low - 1;

            for (j = 0; j < s; j++)
            {
                if (A_ind[j] >= qs_inf->high)
                {
                    flint_printf("index out of bound !!");
                    abort();
                }
            }
        }
        else { ret = 0; }
    }
    else
    {
        if (curr_subset[0] != span - s + 2)
        {
            h = (m >= span - h) ? h + 1 : 1;
            m = curr_subset[s - 1 - h];

            for (j = 1; j <= h; j++)
                curr_subset[s + j - h - 2] = m + j;

            fmpz_set_ui(prod, UWORD(1));

            for (j = 0; j < s - 1; j++)
                fmpz_mul_ui(prod, prod, factor_base[2 * curr_subset[j] - 1 + low - 1].p);

            i = 1;
            j = qs_inf->num_primes - 1;

            while (i < j)
            {
                mid = i + (j - i) / 2;

                fmpz_mul_ui(temp, prod, factor_base[2 * mid - 1 + low].p);

                if (fmpz_cmp(lower_bound, temp) > 0)
                {
                    i = mid;
                }
                else if (fmpz_cmp(temp, upper_bound) > 0)
                {
                    j = mid;
                }
                else
                {
                    j = 2 * mid - 1 + low;
                    break;
                }
            }

            if (j >= qs_inf->num_primes)
            {
                flint_printf("A0 doesn't fit in bounds\n");
                abort();
            }

            A_ind[s - 1] = j;

            for (j = 0; j < s - 1; j++)
                A_ind[j] = 2 * curr_subset[j] - 1 + low - 1;
            fmpz_set(prod, temp);

         }
         else { ret = 0; }
    }

    qs_inf->h = h;
    qs_inf->m = m;
    fmpz_set(qs_inf->A0, prod);

    fmpz_clear(prod);
    fmpz_clear(temp);

    return ret;
}

/* re-initialize parameter for 'A0' after factor base increment */

void qsieve_re_init_A0(qs_t qs_inf)
{
    slong m, h, span, low, high, s, j;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * curr_subset = qs_inf->curr_subset;
    prime_t * factor_base = qs_inf->factor_base;

    fmpz_t prod;
    fmpz_init(prod);

    low =  qs_inf->low;
    high = qs_inf->num_primes;
    span = high - low;
    s = qs_inf->s;
    m = 0;
    h = s;

    fmpz_set_ui(prod, UWORD(1));

    for (j = 1; j <= h; j++)
        curr_subset[s + j - h - 1] = m + j;

    for (j = 0; j < s; j++)
        fmpz_mul_ui(prod, prod, factor_base[curr_subset[j] + low - 1].p);

    for (j = 0; j < s; j++)
        A_ind[j] = curr_subset[j] + low - 1;

    qs_inf->m = m;
    qs_inf->h = h;
    qs_inf->low = low;
    qs_inf->high = high;
    qs_inf->span = span;
    fmpz_set(qs_inf->A0, prod);

    fmpz_clear(prod);
}

void qsieve_init_A0(qs_t qs_inf)
{
    slong i, j = 0;
    slong s, low, high, span, m, h;
    mp_limb_t bits, num_factor, rem, mid;
    mp_limb_t * factor_bound = flint_malloc(32 * sizeof(mp_limb_t));
    mp_limb_t * A_ind;
    mp_limb_t * curr_subset;
    prime_t * factor_base = qs_inf->factor_base;
    fmpz_t target, prod, temp, temp2, upper_bound, lower_bound;
    fmpz_init(target);
    fmpz_init(temp);
    fmpz_init(temp2);
    fmpz_init(upper_bound);
    fmpz_init(lower_bound);
    fmpz_init_set_ui(prod, UWORD(1));
    fmpz_fdiv_q_ui(target, qs_inf->target_A, UWORD(2) * factor_base[qs_inf->q_idx].p);
    fmpz_fdiv_q_ui(lower_bound, target, UWORD(5));
    fmpz_mul_ui(upper_bound, target, UWORD(5));
    bits = fmpz_bits(target);

    for (i = qs_inf->small_primes; i < qs_inf->num_primes; i++)
    {
        if (qs_inf->factor_base[i].size != j)
        {
            factor_bound[j] = i;
            j = qs_inf->factor_base[i].size;
        }
    }

    /*
       following 'for' loop is taken from 'msieve' implementation,
       try to determine number of factors, such that we have enough
       primes to choose from
    */

    if (bits > 210) i = 15;
    else if (bits > 190) i = 13;
    else if (bits > 180) i = 12;
    else i = 11;

    for (num_factor = rem = 0; i > 7; i--)
    {
        num_factor = bits / i;
        rem = bits % i;

        if (factor_bound[i] == 0 || num_factor == 1) continue;
        if (rem == 0 && num_factor > 2 && factor_bound[i + 1] > 0)
        {
            low = factor_bound[i - 1];
            high = factor_bound[i + 1];
            break;
        }
        else if (rem <= num_factor)
        {
            if (num_factor > 2 && factor_bound[i + 1] > 0 && factor_bound[i + 2] > 0)
            {
                low = factor_bound[i - 1];
                high = factor_bound[i + 2];
                break;
            }
        }
        else if (i - rem <= num_factor)
        {
            if (factor_bound[i + 1] > 0 && factor_bound[i - 1] > 0)
            {
                num_factor += 1;
                low = factor_bound[i - 1];
                high = factor_bound[i + 1];
                break;
            }
        }
    }


    if (i == 7)
    {
        num_factor = (bits >= 15) ? 3 : 2;
        low = qs_inf->small_primes;
        high = qs_inf->num_primes;
    }

    s = num_factor;
    qs_inf->s = s;

    qsieve_poly_init(qs_inf);
    A_ind = qs_inf->A_ind;
    curr_subset = qs_inf->curr_subset;

    span = high - low;


    /*
       factor of 'A0' :
       if number of factor is less than 3 just generate 3-subset of
       possible range from factor base lexicographically
       else generate (s - 1) - subset of odd indices from possible
       range from factor base lexicographically and choose last factor
       from even index so that product of values is close to the
       approximated value of hyper-cube
    */

    if (s <= 3)
    {
        m = 0;
        h = s;

        for (j = 1; j <= h; j++)
            curr_subset[s + j - h - 1] = m + j;

        fmpz_set_ui(prod, UWORD(1));

        for (j = 0; j < s; j++)
            fmpz_mul_ui(prod, prod, factor_base[curr_subset[j] + low - 1].p);

        for (j = 0; j < s; j++)
            A_ind[j] = curr_subset[j] + low - 1;
    }
    else
    {
        m = 0;
        h = s - 1;

        for (j = 1; j <= h; j++)
            curr_subset[s - 1 + j - h - 1] = m + j;

        while (curr_subset[0] != span - s + 2)
        {
            fmpz_set_ui(prod, UWORD(1));

            for (j = 0; j < s - 1; j++)
            {
                fmpz_mul_ui(prod, prod, factor_base[2 * curr_subset[j] - 1 + low - 1].p);
            }

            i = 1;
            j = qs_inf->num_primes - 1;

            while (i < j)
            {
                mid = i + (j - i) / 2;

                fmpz_mul_ui(temp, prod, factor_base[2 * mid - 1 + low].p);

                if (fmpz_cmp(lower_bound, temp) > 0)
                {
                    i = mid;
                }
                else if (fmpz_cmp(temp, upper_bound) > 0)
                {
                    j = mid;
                }
                else
                {
                    j = 2 * mid - 1 + low;
                    break;
                }
            }

            if (j < qs_inf->num_primes) break;

            h = (m >= span - h) ? h + 1 : 1;
            m = curr_subset[s - 1 - h];

            for (j = 1; j <= h; j++)
                curr_subset[s + j - h - 2] = m + j;
        }

        high = j;

        A_ind[s - 1] = (j == qs_inf->num_primes) ? j - 1 : j;

        for (j = 0; j < s - 1; j++)
            A_ind[j] = 2 * curr_subset[j] - 1 + low - 1;

        fmpz_set(prod, temp);
    }

    qs_inf->h = h;
    qs_inf->m = m;
    qs_inf->low = low;
    qs_inf->high = high;
    qs_inf->span = span;
    fmpz_set(qs_inf->A0, prod);
    fmpz_set(qs_inf->low_bound, lower_bound);
    fmpz_set(qs_inf->upp_bound, upper_bound);

#if QS_DEBUG
    flint_printf("number of factor in hypercube = %wd\n", qs_inf->s);
    flint_printf("factor base size = %wd max prime = %wu\n", qs_inf->num_primes, qs_inf->factor_base[qs_inf->num_primes - 1].p);
    flint_printf("possible candidate in range [ %wd, %wd ]\n", qs_inf->low, qs_inf->high);
    flint_printf("optimal value of hypercube = ");
    fmpz_print(target);
    flint_printf("\n");
    flint_printf("lower bound       = ");
    fmpz_print(lower_bound);
    flint_printf("\n");
    flint_printf("upper bound       = ");
    fmpz_print(upper_bound);
    flint_printf("\n");
    flint_printf("initial hypercube = ");
    fmpz_print(qs_inf->A0);
    flint_printf("\n");
#endif

    fmpz_clear(target);
    fmpz_clear(prod);
    fmpz_clear(temp);
    fmpz_clear(temp2);
    fmpz_clear(upper_bound);
    fmpz_clear(lower_bound);
    flint_free(factor_bound);
}

/* precompute data associated with A0 which will remain fixed */

void qsieve_compute_pre_data(qs_t qs_inf)
{
    slong i, j;
    slong s = qs_inf->s;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * A0_inv = qs_inf->A0_inv;
    fmpz_t * A0_divp = qs_inf->A0_divp;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    mp_limb_t p, pinv, temp, temp2;

    /*
       calculate $A0 / p$ and $B0_i = \sqrt{kn} * (A0 / p)^{-1} modulo p$
       where 'p' is prime factor of 'A0'
    */
    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;

        fmpz_divexact_ui(A0_divp[i], qs_inf->A0, p);
        temp2 = fmpz_fdiv_ui(A0_divp[i], p);
        temp2 = n_invmod(temp2, p);
        temp2 = n_mulmod2_preinv(temp2, sqrts[A_ind[i]], p, pinv);
        B0_terms[i] = temp2;
    }

    /* calculate $A0^{-1} modulo p$ where $p \in factor_base$ */
    for (i = 2; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        temp = fmpz_fdiv_ui(qs_inf->A0, p);
        A0_inv[i] = n_invmod(temp, p);
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
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * A0_inv = qs_inf->A0_inv;
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    fmpz_t * B_terms = qs_inf->B_terms;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    fmpz_t * A0_divp = qs_inf->A0_divp;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    mp_limb_t q, qinv, qsqrt, p, pinv, temp, temp2, pmod, mod_inv, Ainv;
    fmpz_t temp3;
    fmpz_init(temp3);

    q = factor_base[qs_inf->q_idx].p;
    fmpz_zero(qs_inf->B);
    fmpz_mul_ui(qs_inf->A, qs_inf->A0, q);

    /*
       compute $(A/p) = (A0 / p) * q0$, and
       $B_i = (A / p) * \sqrt{kn} (A / p)^(-1) modulo p$
       where 'p' are prime factor of 'A0';
    */

    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;

        fmpz_mul_ui(temp3, A0_divp[i], q);
        temp2 = n_invmod(n_mod2_preinv(q, p, pinv), p);
        temp2 = n_mulmod2_preinv(temp2, B0_terms[i], p, pinv);

        if (temp2 > p / 2) temp2 = p - temp2;

        fmpz_mul_ui(B_terms[i], temp3, temp2);
        fmpz_add(qs_inf->B, qs_inf->B, B_terms[i]);
    }

    /*
       adding $A0 * \sqrt{kn} * A0^{-1} modulo q0$
       to $b_1 = \sum _{i=1} ^{s} B_i$
    */

    qinv = factor_base[qs_inf->q_idx].pinv;
    pmod = fmpz_fdiv_ui(qs_inf->A0, q);
    mod_inv = n_invmod(pmod, q);

    qsqrt = sqrts[qs_inf->q_idx];
    pmod = n_mulmod2_preinv(mod_inv, qsqrt, q, qinv);
    fmpz_mul_ui(temp3, qs_inf->A0, pmod);
    fmpz_add(qs_inf->B, qs_inf->B, temp3);

    /*
        calculating $A^(-1) modulo p$, $p \in factor_base$,
        calculate A_inv2B[j][i] = $2 * B_j * A^(-1) modulo p$
        for factor base prime 'p',
        compute roots of base polynomial modulo factor base prime
    */

    for (i = 2; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;

        for (j = 0; j < s; j++)
        {
            temp = fmpz_fdiv_ui(B_terms[j], p);
            temp *= 2;
            if (temp >= p) temp -= p;
            A_inv2B[j][i] = temp;
        }
    }

    for (i = 2; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        pinv = factor_base[i].pinv;
        temp = n_invmod(n_mod2_preinv(q, p, pinv), p);
        Ainv = n_mulmod2_preinv(A0_inv[i], temp, p, pinv);

        temp = fmpz_fdiv_ui(qs_inf->B, p);
        temp2 = sqrts[i] + p - temp;
        temp2 = n_mulmod2_preinv(temp2, Ainv, p, pinv);
        temp2 += qs_inf->sieve_size / 2;
        temp2 = n_mod2_preinv(temp2, p, pinv);
        soln1[i] = temp2;

        temp = n_mulmod2_preinv(sqrts[i], Ainv, p, pinv);
        temp *= 2;
        if (temp >= p) temp -= p;
        temp = soln1[i] + p - temp;
        if (temp >= p) temp -= p;
        soln2[i] = temp;

        for (j = 0; j < s; j++)
            A_inv2B[j][i] = n_mulmod2_preinv(A_inv2B[j][i], Ainv, p, pinv);
    }

    for (i = 2; i < qs_inf->num_primes; i++)
    {
        if (soln1[i] > soln2[i])
            soln1[i] = (soln1[i] + soln2[i]) - (soln2[i] = soln1[i]);
    }

    for (i = 0; i < s; i++)
    {
        soln1[A_ind[i]] = soln2[A_ind[i]] = 0;
    }

    qs_inf->curr_poly = 1;

    fmpz_clear(temp3);
}

/*
    generate next possible value of coefficient 'B', using
    gray-code formula, for current value of 'A'
*/

void qsieve_init_poly_next(qs_t qs_inf)
{
    slong i, j, v;
    slong s = qs_inf->s;
    prime_t * factor_base = qs_inf->factor_base;
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    mp_limb_t sign, p;
    fmpz_t temp;
    fmpz_init(temp);

    /* we have $b_i$, calculating $b_{i+1}$ using gray code formula */
    i = qs_inf->curr_poly;

    for (v = 0; v < s; v++)
    {
        if (((i >> v) & 1)) break;
    }

    sign = (i >> v) & 2;

    fmpz_mul_ui(temp, qs_inf->B_terms[v], UWORD(2));

    if (sign) fmpz_add(qs_inf->B, qs_inf->B, temp);
    else fmpz_sub(qs_inf->B, qs_inf->B, temp);

    fmpz_mul(temp, qs_inf->B, qs_inf->B);
    fmpz_mod(temp, temp, qs_inf->A);

    fmpz_t temp2;
    fmpz_init_set(temp2, qs_inf->kn);
    fmpz_mod(temp2, temp2, qs_inf->A);

    /* updating roots for $b_{i+1}$ */
    for (j = 2; j < qs_inf->num_primes; j++)
    {
        p = factor_base[j].p;

        if (sign)
        {
            soln1[j] = soln1[j] + p - A_inv2B[v][j];
            soln2[j] = soln2[j] + p - A_inv2B[v][j];
        }
        else
        {
            soln1[j] = soln1[j] + A_inv2B[v][j];
            soln2[j] = soln2[j] + A_inv2B[v][j];
        }

        if (soln1[j] >= p) soln1[j] -= p;
        if (soln2[j] >= p) soln2[j] -= p;

        if (soln1[j] > soln2[j])
            soln1[j] = (soln1[j] + soln2[j]) - (soln2[j] = soln1[j]);
    }

    i++;
    qs_inf->curr_poly = i;

    fmpz_clear(temp);
}

/* calculate coefficient 'C' for current 'A' and 'B' */
void qsieve_compute_C(qs_t qs_inf)
{
    fmpz_mul(qs_inf->C, qs_inf->B, qs_inf->B);
    fmpz_sub(qs_inf->C, qs_inf->C, qs_inf->kn);
    fmpz_divexact(qs_inf->C, qs_inf->C, qs_inf->A);
}

