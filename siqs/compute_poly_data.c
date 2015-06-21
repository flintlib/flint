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


void qsieve_next_A0(qs_t qs_inf)
{
    slong j;
    slong s = qs_inf->s;
    slong low = qs_inf->low;
    slong span = qs_inf->span;
    slong h = qs_inf->h;
    slong m = qs_inf->m;
    mp_limb_t * current_subset = qs_inf->current_subset;
    mp_limb_t * A_ind = qs_inf->A_ind;
    prime_t * factor_base = qs_inf->factor_base;
    fmpz_t prod, temp;
    fmpz_init(prod);
    fmpz_init(temp);

    /*
       generate lexicographically next s-subset from odd indices
       of possible range of factor's of 'A0'
    */
    if (current_subset[0] != span - s + 1)
    {
        if (m >= span - h) ++h;
        else h = 1;

        m = current_subset[s - h];

        for (j = 1; j <= h; j++)
            current_subset[s + j - h - 1] = m + j;

        fmpz_set_ui(prod, UWORD(1));

        for (j = 1; j <= s; j++)
        {
            fmpz_mul_ui(temp, prod, factor_base[2 * current_subset[j - 1] - 1 + low - 1].p);
            fmpz_set(prod, temp);
        }

        for (j = 0; j < s; j++)
            A_ind[j] = 2 * current_subset[j] - 1 + low - 1;

    }
    else
    {
        /* exception to add as we ran out of 'A0' coefficients */
    }

    fmpz_set(qs_inf->A0, prod);
    qs_inf->s = s;
    qs_inf->low = low;
    qs_inf->h = h;
    qs_inf->m = m;
    qs_inf->current_subset = current_subset;

    fmpz_clear(prod);
    fmpz_clear(temp);
}

void qsieve_init_A0(qs_t qs_inf)
{
    slong i, j = 0;
    slong s = qs_inf->s;
    slong low = qs_inf->low;
    slong high = qs_inf->high;
    slong span = qs_inf->span;
    slong m = qs_inf->m;
    slong h = qs_inf->h;
    mp_limb_t bits, num_factor, rem, mid;
    mp_limb_t * factor_bound = flint_malloc(32 * sizeof(mp_limb_t));
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * current_subset = qs_inf->current_subset;
    prime_t * factor_base = qs_inf->factor_base;
    fmpz_t target, prod, temp, temp2, upper_bound, lower_bound;
    fmpz_init(target);
    fmpz_init(temp);
    fmpz_init(temp2);
    fmpz_init(upper_bound);
    fmpz_init(lower_bound);
    fmpz_init_set_ui(prod, UWORD(1));
    fmpz_tdiv_q_ui(target, qs_inf->target_A, UWORD(2) * qs_inf->q0);
    fmpz_mul_ui(upper_bound, target, UWORD(2));
    fmpz_tdiv_q_ui(lower_bound, target, UWORD(2));
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
                low = factor_bound[i];
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
        if (bits >= 15) num_factor = 3;
        else num_factor = 2;

        low = qs_inf->small_primes;
        high = qs_inf->num_primes;
    }

    s = num_factor;
    span = high - low;
    if (span > 100) span = 100;

    /* generate first s-subset from range of possible factors */

    current_subset = flint_malloc( s * sizeof(mp_limb_t));
    m = 0;
    h = s;

    for (j = 1; j <= h; j++)
        current_subset[s + j - h - 1] = m + j;

    fmpz_set_ui(prod, UWORD(1));

    for (j = 1; j <= s; j++)
    {
        fmpz_mul_ui(temp, prod, factor_base[2 * current_subset[j - 1] - 1 + low - 1].p);
        fmpz_set(prod, temp);
    }

    qs_inf->current_subset = current_subset;
    qs_inf->h = h;
    qs_inf->m = m;
    qs_inf->low = low;
    qs_inf->span = span;
    qs_inf->s = s;
    fmpz_set(qs_inf->A0, prod);

    for (j = 0; j < s; j++)
        A_ind[j] = 2 * current_subset[j] - 1 + low - 1;

    fmpz_clear(target);
    fmpz_clear(prod);
    fmpz_clear(temp);
    fmpz_clear(temp2);
    fmpz_clear(upper_bound);
    fmpz_clear(lower_bound);
    flint_free(factor_bound);
}


/*
    generate q0, q0 are primes immediately following
    prime bound such that 'kn' is quadratic residue modulo 'q0'
*/


mp_limb_t qsieve_compute_q0(qs_t qs_inf)
{
    mp_limb_t p, kron, pmod, pinv, q0, qr;
    ulong i = 0;

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
            if ((p % 8) == 3 || (p % 8) == 5) kron *= -1;
            pmod /= 2;
        }

        kron *= n_jacobi(pmod, p);

        if (kron == 1)
        {
            q0 = qr = p;
            break;
        }
    }

    /*
       collect all primes 'q0' for which 'kn' is quadratic residue modulo 'q0'
    */

    while (qr <= 3 * q0)
    {
        qs_inf->q0_values = flint_realloc(qs_inf->q0_values, (i + 1) * sizeof(mp_limb_t));
        qs_inf->q0_values[i++] = qr;

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
                if ((p % 8) == 3 || (p % 8) == 5) kron *= -1;
                pmod /= 2;
            }

            kron *= n_jacobi(pmod, p);

            if (kron == 1)
            {
                qr = p;
                break;
            }

        }
    }

    qs_inf->q0 = q0;
    qs_inf->num_q0 = i;

    n_primes_clear(iter);
    fmpz_clear(upper_bound);
    fmpz_clear(temp);
    return 0;
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
        fmpz_init(A0_divp[i]);
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
    mp_limb_t q0 = qs_inf->q0;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * A0_inv = qs_inf->A0_inv;
    mp_limb_t * A_inv = qs_inf->A_inv;
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    fmpz_t * B_terms = qs_inf->B_terms;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    fmpz_t * A0_divp = qs_inf->A0_divp;
    fmpz_t * A_divp = qs_inf->A_divp;
    fmpz_t * B = qs_inf->B;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    mp_limb_t p, pinv, temp, temp2, pmod, mod_inv;
    fmpz_t temp3, temp4;
    fmpz_init(temp3);
    B = flint_malloc((1 << s) * sizeof(fmpz_t));
    fmpz_init(B[0]);
    fmpz_zero(B[0]);
    fmpz_init(temp4);
    fmpz_init(qs_inf->A);

    fmpz_mul_ui(qs_inf->A, qs_inf->A0, q0);

    /*
       compute $(A/p) = (A0 / p) * q0$, and
       $b_i = (A / p) * \sqrt{kn} (A / p)^(-1) modulo p$
       where 'p' are prime factor of 'A0';
    */
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

    /*
       adding $A0 * \sqrt{kn} * A0^{-1} modulo q0$
       to $b_1 = \sum _{i=1} ^{s} B_i$
    */
    p = qs_inf->q0;
    pinv = n_preinvert_limb(p);
    pmod = fmpz_fdiv_ui(qs_inf->A0, p);
    mod_inv = n_invmod(pmod, p);
    pmod = fmpz_fdiv_ui(qs_inf->kn, p);
    pmod = n_sqrtmod(pmod, p);
    pmod = n_mulmod2_preinv(mod_inv, pmod, p, pinv);
    fmpz_mul_ui(temp3, qs_inf->A0, pmod);
    fmpz_set(temp4, B[0]);
    fmpz_add(B[0], temp3, temp4);

    qs_inf->B = B;

    /* calculating $A^(-1) modulo p$, $p \in factor_base$*/
    for (j = 2; j < qs_inf->num_primes; j++)
    {
        p = factor_base[j].p;
        pinv = n_preinvert_limb(p);
        temp = n_invmod(q0, p);
        A_inv[j] = n_mulmod2_preinv(A0_inv[j], temp, p, pinv);
    }

    /*
      calculate A_inv2B[j][i] = $2 * b_j * A^(-1) modulo p$
      for factor base prime 'p'
    */
    for (j = 0; j < s; j++)
    {
        for (i = 2; i < qs_inf->num_primes; i++)
        {
            p = factor_base[i].p;
            pinv = n_preinvert_limb(p);
            temp = fmpz_fdiv_ui(B_terms[j], p);
            temp *= 2;
            if (temp >= p) temp -= p;
            temp = n_mulmod2_preinv(temp, A_inv[i], p, pinv);
            A_inv2B[j][i] = temp;
        }
    }

    /* compute roots of base polynomial modulo factor base prime */
    for (i = 2; i < qs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        pinv = n_preinvert_limb(p);
        temp = fmpz_fdiv_ui(B[0], p);
        temp2 = sqrts[i];
        temp2 = temp2 + p - temp;
        temp2 = n_mulmod2_preinv(temp2, A_inv[i], p, pinv);
        soln1[i] = temp2;
        temp2 = sqrts[i];
        temp2 = p - temp2;
        temp2 = temp2 + p - temp;
        temp2 = n_mulmod2_preinv(temp2, A_inv[i], p, pinv);
        soln2[i] = temp2;
    }

    qs_inf->curr_poly = 1;

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
    mp_limb_t * soln1 = qs_inf->soln1;
    mp_limb_t * soln2 = qs_inf->soln2;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    mp_limb_t sign, p;
    fmpz_t temp;
    fmpz_init(temp);

    /* we have b_i, calculating b_{i+1} using gray code formula */
    i = qs_inf->curr_poly;

    for (v = 0; v < s; v++)
    {
        if (((i >> v) & 1)) break;
    }

    sign = (i >> v) & 2;

    fmpz_mul_ui(temp, qs_inf->B_terms[v], UWORD(2));
    fmpz_init(qs_inf->B[i]);

    if (sign) fmpz_add(qs_inf->B[i], qs_inf->B[i - 1], temp);
    else fmpz_sub(qs_inf->B[i], qs_inf->B[i - 1], temp);

    /* updating roots for , b_{i+1} */
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
    }

    i++;
    qs_inf->curr_poly = i;

    fmpz_clear(temp);
}

/* calculate coefficient 'C' for current 'A' and 'B' */
void qsieve_compute_C(qs_t qs_inf)
{
    fmpz_t temp;
    fmpz_init(temp);
    fmpz_init(qs_inf->C);
    fmpz_pow_ui(qs_inf->C, qs_inf->B[qs_inf->curr_poly - 1], UWORD(2));
    fmpz_sub(temp, qs_inf->C, qs_inf->kn);
    fmpz_divexact(qs_inf->C, temp, qs_inf->A);
    fmpz_clear(temp);
}

