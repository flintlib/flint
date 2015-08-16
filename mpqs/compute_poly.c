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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "mpqs.h"

void
mpqs_sqrtmod_psq(fmpz_t root, fmpz_t n, mp_limb_t prime)
{
    fmpz_t p, psq, a, b, k;
    
    fmpz_init(psq);
    fmpz_init(b);
    fmpz_init(k);
    fmpz_init_set(a, n);
    fmpz_init_set_ui(p, prime);

    fmpz_mod(b, a, p);
    fmpz_sqrtmod(root, a, p);
    fmpz_mul(psq, p, p);

    /* b = root^2 mod p^2 */
    fmpz_mul(b, root, root);
    fmpz_mod(b, b, psq);

    if (fmpz_cmp(a, b) > 0)
        fmpz_sub(k, a, b);
    else
        fmpz_sub(k, b, a);

    fmpz_mod(k, k, psq);
    fmpz_divexact_ui(k, k, prime);

    if (fmpz_cmp(a, b) < 0)
        fmpz_negmod(k, k, p);

    fmpz_mul_2exp(b, root, 1); /* b = root * 2 */
    fmpz_invmod(b, b, p);
    fmpz_mul(k, k, b);
    fmpz_mod(k, k, p);
    fmpz_mul(k, k, p);
    fmpz_add(root, root, k);

    /* if first root is even, second root is odd */
    
    if (fmpz_is_even(root) == 1)
        fmpz_sub(root, psq, root);

    fmpz_clear(psq);
    fmpz_clear(b);
    fmpz_clear(k);
    fmpz_clear(p);
}

void
mpqs_compute_poly_roots(mpqs_t mpqs_inf)
{
    mp_limb_t * A_inv = mpqs_inf->A_inv;
    mp_limb_t * soln1 = mpqs_inf->soln1;
    mp_limb_t * soln2 = mpqs_inf->soln2;
    prime_t * factor_base = mpqs_inf->factor_base;
    int * sqrts = mpqs_inf->sqrts, i;
    mp_limb_t p, pinv, temp;

    for (i = 2; i < mpqs_inf->num_primes; i++)
    {
        p = factor_base[i].p;
        pinv = factor_base[i].pinv;

        /* A_inv[i] = A^-1 % p */
        A_inv[i] = n_invmod(fmpz_fdiv_ui(mpqs_inf->A, p), p);

        /* soln1 & soln2 are roots of Q(x) = N mod p for each p in F.B. */

        /* root 1 = (-b + n^1/2)/a mod p */

        temp = fmpz_fdiv_ui(mpqs_inf->B, p);  /* temp = B % p */
        temp = sqrts[i] + p - temp;
        temp = n_mulmod2_preinv(temp, A_inv[i], p, pinv);
        temp += mpqs_inf->sieve_size / 2;
        temp = n_mod2_preinv(temp, p, pinv);
        soln1[i] = temp;

        /* root 2 = (-b - n^1/2)/a mod p */

        temp = n_mulmod2_preinv(sqrts[i], A_inv[i], p, pinv);
        temp *= 2;
        if (temp >= p) temp -= p;
        temp = soln1[i] + p - temp;
        if (temp >= p) temp -= p;
        soln2[i] = temp;

    }

    for (i = 2; i < mpqs_inf->num_primes; i++)
    {
        if (soln1[i] > soln2[i])
            soln1[i] = (soln1[i] + soln2[i]) - (soln2[i] = soln1[i]);
    }
}

void
mpqs_find_A_targetprime(mpqs_t mpqs_inf)
{
    /* prime near d = (N/(2 * M^2))^1/4 */ 

    slong sieve_size;
    fmpz_t target;

    sieve_size = mpqs_inf->sieve_size;
    fmpz_init_set(target, mpqs_inf->kn);
    fmpz_fdiv_q_ui(target, target, 2 * sieve_size * sieve_size);
    fmpz_root(target, target, UWORD(4));
    mpqs_inf->A_targetprime = fmpz_get_ui(target);

    mpqs_inf->A_targetprime = n_nextprime(mpqs_inf->A_targetprime, 0);
    
    fmpz_set_ui(mpqs_inf->A, mpqs_inf->A_targetprime);
    fmpz_mul(mpqs_inf->A, mpqs_inf->A, mpqs_inf->A);

}

void
mpqs_compute_A(mpqs_t mpqs_inf)
{
    mpqs_inf->A_targetprime = n_nextprime(mpqs_inf->A_targetprime, 0);

    fmpz_set_ui(mpqs_inf->A, mpqs_inf->A_targetprime);
    fmpz_mul(mpqs_inf->A, mpqs_inf->A, mpqs_inf->A);
}

void
mpqs_compute_B(mpqs_t mpqs_inf)
{
    mp_limb_t prime = mpqs_inf->A_targetprime;

    /* odd root of x^2 = n mod A */
    /* A is A_targetprime ^ 2 */

    mpqs_sqrtmod_psq(mpqs_inf->B, mpqs_inf->kn, prime);
}

void
mpqs_compute_C(mpqs_t mpqs_inf)
{
    /* C = (B^2 - kn)/A */

    fmpz_set(mpqs_inf->C, mpqs_inf->B);
    fmpz_mul(mpqs_inf->C, mpqs_inf->C, mpqs_inf->C);      /* B^2 */
    fmpz_sub(mpqs_inf->C, mpqs_inf->C, mpqs_inf->kn);     /* B^2 - kn */
    fmpz_fdiv_q(mpqs_inf->C, mpqs_inf->C, mpqs_inf->A);   /* (B^2 - kn)/A */
}

void
mpqs_poly_init(mpqs_t mpqs_inf)
{
    mp_limb_t num_primes = mpqs_inf->num_primes;

    fmpz_init(mpqs_inf->A);
    fmpz_init(mpqs_inf->C);
    fmpz_init(mpqs_inf->B);

    mpqs_inf->A_inv = flint_malloc(num_primes * sizeof(mp_limb_t));
    mpqs_inf->soln1 = flint_malloc(num_primes * sizeof(mp_limb_t));
    mpqs_inf->soln2 = flint_malloc(num_primes * sizeof(mp_limb_t));

    mpqs_find_A_targetprime(mpqs_inf);
    mpqs_compute_poly_roots(mpqs_inf);
    mpqs_compute_B(mpqs_inf);
    mpqs_compute_C(mpqs_inf);    
}

void
mpqs_compute_next_poly(mpqs_t mpqs_inf)
{
    mpqs_compute_A(mpqs_inf);
    mpqs_compute_poly_roots(mpqs_inf);
    mpqs_compute_B(mpqs_inf);
    mpqs_compute_C(mpqs_inf);
}
