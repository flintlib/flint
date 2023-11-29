/*
    Copyright (C) 2016, 2020 William Hart
    Copyright (C) 2015 Nitin Kumar

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "qsieve.h"

/*
    * compute bounds for optimal A (coeff of poly f(x) = Ax^2 + 2Bx + C, where
      C = (B^2 - kn)/A; note (Ax + B)^2 - kn = A*f(x))
    * compute optimal number of factors of A
    * allocate space for polynomial data
    * compute first A value (a product of factor base primes)
    * set up data for lexicographic ordering for subsequent A coeffs

    We try to compute A such that:
    (i) A is close to it's ideal value which is, sqrt(2 * kn) / M,
        where kn is the number to be factored and M is the size of sieve
   (ii) A should not contain "small" prime factors,
        it will decrease number of relations
  (iii) A should not contain primes factor that are too large, or (i) will be
        difficult
   (iv) A should have as many factors as possible, to obtain a large number of
        polynomials

    If the optimal number of factors of A is <= 3 we simply take every possible
    product of three factor base primes in the allowed range, as we will need as
    many possible A's as we can get, as the number of polynomials we get from
    any given A is just 2^(s-1) where s is the number of prime factors of A.

    Otherwise, we pick s-1 factors using a lexicographic ordering so that we
    don't get any duplicate A's. The factors are chosen from the array of
    factor base primes such that their indices (starting with index 0
    corresponding to the smallest allowed factor) aren't divisible by 4. Then
    we choose the final factor of A to optimise the value of A as above, and
    such that it's index (as just described) is divisible by 4.

    Return value indicates success. Should always succeed unless tuning
    parameters are really screwed up!
*/
int qsieve_init_A(qs_t qs_inf)
{
    slong i, j;
    slong s, low, high, span, m, h;
    mp_limb_t bits, num_factors, rem, mid;
    mp_limb_t factor_bound[40];
    mp_limb_t * A_ind;
    mp_limb_t * curr_subset, * first_subset;
    prime_t * factor_base = qs_inf->factor_base;
    fmpz_t prod, temp, upper_bound, lower_bound;
    int ret = 1, found_j;

    fmpz_init(temp);
    fmpz_init(upper_bound);
    fmpz_init(lower_bound);
    fmpz_init_set_ui(prod, 1);

    fmpz_fdiv_q_ui(lower_bound, qs_inf->target_A, 2);
    fmpz_mul_ui(upper_bound, qs_inf->target_A, 2);

    bits = fmpz_bits(qs_inf->target_A);

    /*
       compute indices of points in factor base where number of bits of primes
       increases, starting with the first non-small prime, i.e. factor_bound[0]
       will be the index of the first non-small prime, and then factor_bound[j]
       will be the first prime that has more than j bits
    */
    for (i = 0; i < 40; i++)
       factor_bound[i] = 0;

    for (j = 0, i = qs_inf->small_primes; i < qs_inf->num_primes; i++)
    {
        if (qs_inf->factor_base[i].size != j)
        {
            factor_bound[j] = i;
            j = qs_inf->factor_base[i].size;
        }
    }

    /*
       the following loop is taken from msieve,
       it tries to determine number of factors, such that we have enough
       primes to choose factors of A from
    */

    if (bits > 210) i = 15;
    else if (bits > 190) i = 13;
    else if (bits > 180) i = 12;
    else i = 11;

    high = low = 0;

    for ( ; i > 7; i--)
    {
        num_factors = bits / i;
        rem = bits % i;

        if (factor_bound[i] == 0 || num_factors == 1)
            continue; /* factor too large or not enough factors */

        /*
            let n = bits, num_factors = floor(n/i), rem = n - i*num_factors
            we can only guarantee that bounds cover [2^(n - 1), 2^n]
            when factor base is very small, the algorithm is forced to use
            very small primes which grow too rapidly; to compensate we
            allow primes that are half the size on the lower end
        */
        if (rem == 0 && num_factors > 2 && factor_bound[i + 1] > 0 && factor_bound[i - 1 - (i <= 10)] > 0)
        {
            /*
                optimal case, primes in range [2^(i-1), 2^(i+1)]
                products cover [2^(n - n/i), 2^(n + n/i)]
            */
            low = factor_bound[i - 1 - (i <= 10)];
            high = factor_bound[i + 1];
            break;
        }
        else if (rem <= num_factors)
        {
            /*
                some left over bits, let primes be in [2^i, 2^(i+2)]
                products in [2^(n - rem), 2^(n - rem + 2*num_factors)]
            */
            if (num_factors > 2 && factor_bound[i + 2] > 0 && factor_bound[i + - (i <= 9)] > 0)
            {
                low = factor_bound[i - (i <= 9)];
                high = factor_bound[i + 2];
                break;
            }
        }
        else if (i - rem <= num_factors)
        {
            /*
                nearly enough bits for extra factor
                primes will be in range [2^(i-1), 2^(i+1)]
                product will be in [2^(n + (i - rem) - num_factors - 1),
                                    2^(n + (i - rem) + num_factors + 1)]
            */
            if (factor_bound[i + 1] > 0 && factor_bound[i - 1 - (i <= 10)] > 0)
            {
                num_factors++;
                low = factor_bound[i - 1 - (i <= 10)];
                high = factor_bound[i + 1];
                break;
            }
        }
    }

    /* not successful, must be a small factorisation */
    if (i == 7)
    {
        num_factors = (bits >= 15) ? 3 : 2;
        low = qs_inf->small_primes;
        high = qs_inf->num_primes;
    }

    s = num_factors;
    qs_inf->s = s;

#if QS_DEBUG
    flint_printf("s = %wd\n", s);
    flint_printf("high = %wd\n", high);
    flint_printf("low = %wd\n", low);
    flint_printf("span = %wd\n", high - low);
#endif

    qsieve_poly_init(qs_inf); /* allocate space for poly data */

    A_ind = qs_inf->A_ind; /* indices of primes dividing A */
    curr_subset = qs_inf->curr_subset;
    first_subset = qs_inf->first_subset;

    span = high - low;

    /*
       factors of A :
       if number of factors is at most 3 just generate 3-subset of
       possible range from factor base lexicographically
       else generate (s - 1) - subset of indices 1, 2, 3 mod 4 from possible
       range from factor base lexicographically and choose last factor
       from index 0 mod 4 so that product of values is close to the
       approximate optimal value for hypercube
    */

    if (s <= 3) /* small factorisation, will need all possible triples */
    {
        m = 0;
        h = s;

        for (j = 0; j < s; j++) /* start with first s allowed primes */
        {
           curr_subset[j] = j;
           first_subset[j] = j; /* save first tuple in case of restart */
        }

        fmpz_set_ui(prod, 1);

        for (j = 0; j < s; j++)
        {
           fmpz_mul_ui(prod, prod, factor_base[curr_subset[j] + low].p);
           A_ind[j] = curr_subset[j] + low;
        }
    }
    else
    {
        m = 0;
        h = s - 1;

        for (j = 0; j < s - 1; j++) /* first s - 1 allowed primes, indices not 0 mod 4 */
            curr_subset[j] = j;

        /* search until we find A of the right size, or until we run out of allowed primes */
        while (1)
        {
            if (4*(curr_subset[0] + s - 2)/3 >= span)
            {
                ret = 0;
                goto init_A_cleanup;
            }

            fmpz_set_ui(prod, 1);

            for (j = 0; j < s - 1; j++)
            {
                /* only pick primes whose indices are not 0 mod 4 */
                fmpz_mul_ui(prod, prod, factor_base[4*curr_subset[j]/3 + 1 + low].p);
            }

            /* binary search for final prime */
            i = 0;
            j = span/4 - 1;

            found_j = 0;

            while (i < j)
            {
                mid = i + (j - i) / 2;

                fmpz_mul_ui(temp, prod, factor_base[4*mid + low].p);

                if (fmpz_cmp(lower_bound, temp) > 0)
                {
                    i = mid + (i == mid);
                }
                else if (fmpz_cmp(temp, upper_bound) > 0)
                {
                    j = mid - (j == mid);
                }
                else
                {
                    j = 4*mid + low;
                    found_j = 1;
                    break;
                }
            }

            if (found_j) break; /* success */

            /* (s - 1)-tuple failed, step to next (s - 1)-tuple */
            h = (4*(m + h + 1)/3 >= span) ? h + 1 : 1;
            m = curr_subset[s - h - 1] + 1;

            for (j = 0; j < h; j++)
                curr_subset[s + j - h - 1] = m + j;
        }

        A_ind[s - 1] = j;

        fmpz_mul_ui(prod, prod, qs_inf->factor_base[A_ind[s - 1]].p);

        for (j = 0; j < s - 1; j++)
            A_ind[j] = 4*curr_subset[j]/3 + 1 + low;

        for (j = 0; j < s; j++)
            first_subset[j] = curr_subset[j]; /* save 1st tuple for restarts */

#if QS_DEBUG
        flint_printf("First A_ind = (");
        for (i = 0; i < s - 1; i++)
            flint_printf("%ld, ", A_ind[i]);
        flint_printf("%ld", A_ind[s - 1]);
        flint_printf(")\n");
#endif
    }

    if (s > 3)
    {
        qs_inf->j = A_ind[s - 1]; /* save s-th factor of A if s > 3 (restart) */
        qs_inf->A_ind_diff = 1;
    }

    qs_inf->h = h;
    qs_inf->m = m;
    qs_inf->low = low;
    qs_inf->high = high;
    qs_inf->span = span;

    fmpz_set(qs_inf->A, prod);

    fmpz_set(qs_inf->low_bound, lower_bound);

    fmpz_set(qs_inf->upp_bound, upper_bound);

#if QS_DEBUG
    flint_printf("number of factors in hypercube = %wd\n", qs_inf->s);
    flint_printf("factor base size = %wd max prime = %wu\n", qs_inf->num_primes, qs_inf->factor_base[qs_inf->num_primes - 1].p);
    flint_printf("possible candidate in range [ %wd, %wd ]\n", qs_inf->low, qs_inf->high);
    flint_printf("optimal value of hypercube = ");
    fmpz_print(qs_inf->target_A);
    flint_printf("\n");
    flint_printf("lower bound       = ");
    fmpz_print(lower_bound);
    flint_printf("\n");
    flint_printf("upper bound       = ");
    fmpz_print(upper_bound);
    flint_printf("\n");
    flint_printf("initial hypercube = ");
    fmpz_print(qs_inf->A);
    flint_printf("\n");
#endif

init_A_cleanup:

    fmpz_clear(prod);
    fmpz_clear(temp);
    fmpz_clear(upper_bound);
    fmpz_clear(lower_bound);

    return ret; /* success */
}

/*
   If the factorisation fails due to the factor base being too small we restart
   the computation with a new factor base. As the size of factors and the
   number of factors of A will not change, we simply recompute A from the
   data that was already computed for the original factor base, starting again
   from the first A.
*/
void qsieve_reinit_A(qs_t qs_inf)
{
    slong low, s, j;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * curr_subset = qs_inf->curr_subset;
    mp_limb_t * first_subset = qs_inf->first_subset;
    prime_t * factor_base = qs_inf->factor_base;

    low = qs_inf->low;
    s = qs_inf->s;

    fmpz_set_ui(qs_inf->A, UWORD(1));

    if (s <= 3)
    {
        for (j = 0; j < s; j++)
        {
           curr_subset[j] = first_subset[j]; /* restore first tuple */
           A_ind[j] = curr_subset[j] + low;
        }
    } else
    {
        for (j = 0; j < s - 1; j++)
        {
            curr_subset[j] = first_subset[j]; /* restore first tuple */
            A_ind[j] = 4*curr_subset[j]/3 + low;
        }

        A_ind[s - 1] = qs_inf->j;
    }

    for (j = 0; j < s; j++)
        fmpz_mul_ui(qs_inf->A, qs_inf->A, factor_base[A_ind[j]].p);

    qs_inf->m = 0;
    qs_inf->h = s;
}

/*
    Compute next A coefficient (differing from previous tuple in at least two
    positions). This is called every time we run out of B coefficients
    for a given A coefficient and need to compute a new A coefficient. If the
    return value is 0, we have run out of possible coefficients (should not
    ever happen unless tuning parameters are way off).
*/
int qsieve_next_A(qs_t qs_inf)
{
    slong i, j, mid, diff;
    slong s = qs_inf->s;
    slong low = qs_inf->low;
    slong span = qs_inf->span;
    slong h = qs_inf->h;
    slong m = qs_inf->m;
    mp_limb_t ret = 1;
    mp_limb_t * curr_subset = qs_inf->curr_subset;
    mp_limb_t * A_ind = qs_inf->A_ind;
    prime_t * factor_base = qs_inf->factor_base;
    fmpz_t prod, temp;
    int found_j, inc_diff;

    fmpz_init(prod);
    fmpz_init(temp);

    /* generate next value of A */

    if (s <= 3)
    {
        if (curr_subset[0] != span - s + 1)
        {
            h = (m >= span - h) ? h + 1 : 1;
            m = curr_subset[s - h] + 1;

            for (j = 0; j < h; j++)
                curr_subset[s + j - h] = m + j;

            fmpz_set_ui(prod, UWORD(1));

            for (j = 0; j < s; j++)
                fmpz_mul_ui(prod, prod, factor_base[curr_subset[j] + low].p);

            for (j = 0; j < s; j++)
                A_ind[j] = curr_subset[j] + low;
        }
        else
           ret = 0;
    } else
    {
        diff = qs_inf->A_ind_diff;
        inc_diff = 0;

        while (1)
        {
            if (4*(curr_subset[0] + s + diff)/3 + 1 >= span) /* have run out of A's */
            {
                ret = 0;
                goto next_A_cleanup;
            }

            h = (4*(m + diff + h + 1)/3 >= span) ? h + 1 : 1;
            m = curr_subset[s - 2 - h] + 1 + ((m%diff) == 0);
            if (h == 2)
               inc_diff = 1;
            else if (h > 2)
               diff = 1;

            for (j = 0; j < h; j++)
                curr_subset[s + j - h - 2] = m + j;

            curr_subset[s - 2] = curr_subset[s - 3] + diff;

            fmpz_set_ui(prod, 1);

            for (j = 0; j < s - 1; j++)
                fmpz_mul_ui(prod, prod, factor_base[4*curr_subset[j]/3 + 1 + low].p);

            /* binary search for final prime */
            i = 0;
            j = span/4 - 1;

            found_j = 0;

            while (i < j)
            {
                mid = i + (j - i) / 2;

                fmpz_mul_ui(temp, prod, factor_base[4*mid + low].p);

                if (fmpz_cmp(qs_inf->low_bound, temp) > 0)
                {
                    i = mid + (i == mid);
                }
                else if (fmpz_cmp(temp, qs_inf->upp_bound) > 0)
                {
                    j = mid - (j == mid);
                }
                else
                {
                    j = 4*mid + low;
                    found_j = 1;
                    if (inc_diff)
                    {
                        diff += 1;
                        qs_inf->A_ind_diff = diff;
                    }
                    break;
                }
            }

            if (found_j)
            {
                A_ind[s - 1] = j;

                fmpz_mul_ui(prod, prod, qs_inf->factor_base[j].p);

                for (j = 0; j < s - 1; j++)
                    A_ind[j] = 4*curr_subset[j]/3 + 1 + low;

                break;
            }
        }
    }

#if QS_DEBUG
    flint_printf("A_ind = (");
    for (i = 0; i < s - 1; i++)
        flint_printf("%ld, ", A_ind[i]);
    flint_printf("%ld", A_ind[s - 1]);
    flint_printf(")\n");
#endif

    qs_inf->h = h;
    qs_inf->m = m;
    fmpz_set(qs_inf->A, prod);

next_A_cleanup:

    fmpz_clear(prod);
    fmpz_clear(temp);

    return ret;
}

/*
   * Precompute data for polynomials which depends only on current A and the
     factor base primes and not on the current B coeff:
     A_inv, A_divp, B0_terms, B_terms, A_inv2B
   * Compute the first B coeff
*/
void qsieve_init_poly_first(qs_t qs_inf)
{
    slong i, k;
    slong s = qs_inf->s;
    mp_limb_t * A_ind = qs_inf->A_ind;
    mp_limb_t * A_inv = qs_inf->A_inv;
    mp_limb_t * B0_terms = qs_inf->B0_terms;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    fmpz_t * B_terms = qs_inf->B_terms;
    fmpz_t * A_divp = qs_inf->A_divp;
    prime_t * factor_base = qs_inf->factor_base;
    int * sqrts = qs_inf->sqrts;
    int * soln1 = qs_inf->soln1;
    int * soln2 = qs_inf->soln2;
    mp_limb_t p, pinv, temp, temp2;

#if QS_DEBUG
    qs_inf->poly_count += 1;
#endif

    fmpz_zero(qs_inf->B);

    /*
       calculate:
       * A_divp[i] = A/p_i and
       * B0_terms[i] = min(gamma_i, p - gamma_i) where
         gamma_i = (sqrt(kn)*(A_divp[i])^(-1)) mod p_i,
         where the p_i are the prime factors of A
       * B_terms[i] = A_divp[i]*B0_terms[i] (multi prec.)
       * first B coeff for current A
    */
    for (i = 0; i < s; i++)
    {
        p = factor_base[A_ind[i]].p;
        pinv = factor_base[A_ind[i]].pinv;

        /* compute A_divp[i] */
        fmpz_divexact_ui(A_divp[i], qs_inf->A, p);

        /* compute B0_terms[i] */
        temp = fmpz_fdiv_ui(A_divp[i], p);
        temp = n_invmod(temp, p);
        B0_terms[i] = n_mulmod2_preinv(temp, sqrts[A_ind[i]], p, pinv);

        if (B0_terms[i] > p/2)
           B0_terms[i] = p - B0_terms[i];

        /* compute B_terms[i] */
        fmpz_mul_ui(B_terms[i], A_divp[i], B0_terms[i]);

        /* compute first B */
        fmpz_add(qs_inf->B, qs_inf->B, B_terms[i]);
    }

    /* calculate A_inv[k] = A^-1 modulo p_k for p_k in the factor base */
    for (k = 3; k < qs_inf->num_primes; k++)
    {
        p = factor_base[k].p;
        temp = fmpz_fdiv_ui(qs_inf->A, p);
        A_inv[k] = temp == 0 ? 0 : n_invmod(temp, p);
    }

    /*
        calculate A_inv2B[i][k] = 2 * B_terms[i] * A^-1 modulo p_k
        for factor base primes p_k, will be multiplied by A_inv[k] below
    */
    for (k = 3; k < qs_inf->num_primes; k++)
    {
        p = factor_base[k].p;
        pinv = factor_base[k].pinv;

        for (i = 0; i < s; i++)
        {
            temp = fmpz_fdiv_ui(B_terms[i], p);
            temp *= 2;
            if (temp >= p)
               temp -= p;
            A_inv2B[i][k] = n_mulmod2_preinv(temp, A_inv[k], p, pinv);
        }
    }

    /*
        compute roots of first polynomial modulo factor base primes
    */
    for (k = 3; k < qs_inf->num_primes; k++)
    {
        p = factor_base[k].p;
        pinv = factor_base[k].pinv;

        /* compute first root of poly */
        temp = fmpz_fdiv_ui(qs_inf->B, p);
        temp = sqrts[k] + p - temp;
        temp = n_mulmod2_preinv(temp, A_inv[k], p, pinv);
        temp += qs_inf->sieve_size / 2;
        temp = n_mod2_preinv(temp, p, pinv);

        /* compute second root of poly */
        temp2 = n_mulmod2_preinv(sqrts[k], A_inv[k], p, pinv);
        temp2 *= 2;
        if (temp2 >= p)
           temp2 -= p;
        temp2 = temp + p - temp2;
        if (temp2 >= p)
           temp2 -= p;

        if (temp2 > temp)
        {
            soln1[k] = temp;
            soln2[k] = temp2;
        } else
        {
            soln1[k] = temp2;
            soln2[k] = temp;
        }
    }

    /* zero out roots corresponding to factors of A */
    for (i = 0; i < s; i++)
    {
        soln1[A_ind[i]] = soln2[A_ind[i]] = 0;
    }
}

/*
    Generate (i + 1)-th B coefficient from i-th B coefficient using Gray-code
*/
void qsieve_init_poly_next(qs_t qs_inf, slong i)
{
    slong j, v;
    slong s = qs_inf->s;
    prime_t * factor_base = qs_inf->factor_base;
    int * soln1 = qs_inf->soln1;
    int * soln2 = qs_inf->soln2;
    mp_limb_t ** A_inv2B = qs_inf->A_inv2B;
    mp_limb_t sign, p, r1, r2;
    fmpz_t temp;

    fmpz_init(temp);

#if QS_DEBUG
    qs_inf->poly_count += 1;
#endif

    /* we have B_i, calculating B_{i+1} using Gray code */
    for (v = 0; v < s; v++)
    {
        if (((i >> v) & 1)) break;
    }

    sign = (i >> v) & 2;

    fmpz_mul_ui(temp, qs_inf->B_terms[v], UWORD(2));

    if (sign)
       fmpz_add(qs_inf->B, qs_inf->B, temp);
    else
       fmpz_sub(qs_inf->B, qs_inf->B, temp);

    /* updating roots of poly for B_{i+1} */
    for (j = 3; j < qs_inf->num_primes; j++)
    {
        p = factor_base[j].p;

        if (sign)
        {
            r1 = soln1[j] + p - A_inv2B[v][j];
            r2 = soln2[j] + p - A_inv2B[v][j];
        } else
        {
            r1 = soln1[j] + A_inv2B[v][j];
            r2 = soln2[j] + A_inv2B[v][j];
        }

        if (r1 >= p)
            r1 -= p;
        if (r2 >= p)
            r2 -= p;

        if (r1 < r2)
        {
           soln1[j] = r1;
           soln2[j] = r2;
        } else
        {
            soln1[j] = r2;
            soln2[j] = r1;
        }
    }

    fmpz_clear(temp);
}

/*
    Calculate coefficient C for current A and B where
    C = (B^2 - kn)/A
*/
void qsieve_compute_C(fmpz_t C, qs_t qs_inf, qs_poly_t poly)
{
    fmpz_t r;

    fmpz_init(r);

    fmpz_mul(C, poly->B, poly->B);
    fmpz_sub(C, C, qs_inf->kn);

#if QS_DEBUG
    fmpz_mod(r, C, qs_inf->A);
    if (!fmpz_is_zero(r))
    {
       flint_throw(FLINT_ERROR, "B^2 - kn not divisible by A\n");
    }
#endif

    fmpz_divexact(C, C, qs_inf->A);

    fmpz_clear(r);
}

/*
    Make a copy of a poly structure. Used by threads so they have their own
    poly structure to modify.
*/
void qsieve_poly_copy(qs_poly_t poly, qs_t qs_inf)
{
   slong i;

   fmpz_set(poly->B, qs_inf->B);

   for (i = 0; i < qs_inf->num_primes; i++)
   {
      poly->soln1[i] = qs_inf->soln1[i];
      poly->soln2[i] = qs_inf->soln2[i];
   }
}
