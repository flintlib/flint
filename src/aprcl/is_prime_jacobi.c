/*
    Copyright (C) 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "thread_support.h"
#if FLINT_USES_PTHREAD
# include <pthread.h>
#endif
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "aprcl.h"

/*
    Below is the implementation of primality test using Jacobi sums.
    The steps are described in:
        [1] "A Course in Computational Algebraic Number Theory" by H. Cohen

    For a different version of the algorithm, also see:
        [2] "Implementation of a New Primality Test" by H. Cohen and A.K. Lenstra

    The algorithm consist of 4 steps:
        (1.) Precomutation;
        (2.) Pseudoprime tests with Jacobi sums;
        (3.) Additional tests;
        (4.) Final trial division and primality proving.

    This file contains the implementation of steps (2.) and (3.).
    It also contains the Jacobi sum primality test.
    A small part of implementation of step (1.) is here and most are in
    config_jacobi.c file.
    (4.) implemented in function aprcl_is_prime_final_division().

    Standard variables:
        n - number to check for primality;

        R - configuration parameter;
        s - configuration parameter depends on the R; s^2 > n;
        a^R = 1 mod s for all integer a coprime to s;

        q - prime number such that (q | s);
        p - prime number such that (p | q - 1);
        r - prime power p^k such that (r | q - 1) and not (r * p | q - 1);
        k - power of p in r;

        v - n % r;
        u - n / r;

    Standard notation:
        \zeta_p - p-th root of unity;
        \sigma_x - automorphism of Z[\zeta_p] for which
                \sigma_x(\zeta_p) = \zeta_p^x;
        \sigma_x^{-1} - inverse of \sigma_x;

        \chi_{p, q} - character defined by
                \chi_{p, q}(g^x) = \zeta_{p^k}^x,
                there g is a primitive root modulo q;

        J(p, q) - jacobi sum j(\chi_{p, q}, \chi_{p, q});
        J_2(q)  - jacobi sum
                ( j(\chi_{2, q}^{2^{k - 3}}, \chi_{2, q}^{3 * 2^{k - 3}}) )^2;
        J_3(q)  - jacobi sum
                j(\chi_{2, q}, \chi_{2, q}, \chi_{2, q}) =
                J(2, q) * j(\chi_{2, q}^2, \chi_{2, q});
*/

/*----------------------------------------------------------------------------*/

/*
    Checks the case p != 2.
    Computes j0 = j_{0, p, q}, jv = j_{v, p, q} and checks that
    j0^u * jv is root of unity.

    Parameters:
        j = J(p, q);
        u, v from standard variables;

    Returns:
        If there exist h such that j0^u * jv = \zeta_{p^k}^h returns h;
        otherwise returns -1.

    For details about j0 and jv see (i1a) in [2] or
    algorithm (9.1.28) step 4.a in [1].
*/
slong
_aprcl_is_prime_jacobi_check_pk(const unity_zp j, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, r, e;
    unity_zp j0, jv, temp, aut, jr, ji, je;

    /* use mpn_mod arithmetic whenever the modulus is in range */
    if (_aprcl_is_prime_jacobi_check_pk_mpn(&h, j, u, v))
        return h;

    /* initialization */
    r = n_pow(j->p, j->exp);            /* r = p^k */
    unity_zp_init(j0, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(jv, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(temp, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(aut, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(jr, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(ji, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(je, j->p, j->exp, fmpz_mod_ctx_modulus(j->ctx));

    unity_zp_coeff_set_ui(j0, 0, 1);    /* j0 = 1 */
    unity_zp_coeff_set_ui(jv, 0, 1);    /* jv = 1 */

    /*
        The powers j^i and j^{(v * i) / r} are maintained incrementally:
        ji = j^i and je = j^e with e = (v * i) / r, which grows by at most
        one per step since v < r.
    */
    unity_zp_reduce_cyclotomic(jr, j);  /* jr = j reduced mod Phi_{p^k} */
    unity_zp_coeff_set_ui(ji, 0, 1);
    unity_zp_coeff_set_ui(je, 0, 1);
    e = 0;

    /* for i in 1..p^k */
    for (i = 1; i <= r; i++)
    {
        /* ji = j^i */
        unity_zp_mul(temp, ji, jr);
        unity_zp_swap(temp, ji);

        /* je = j^{(v * i) / r}; note: overflow in v * i is not a concern */
        while (e < (v * i) / r)
        {
            unity_zp_mul(temp, je, jr);
            unity_zp_swap(temp, je);
            e++;
        }

        /* only for i that coprime to p */
        if (i % j->p == 0) continue;

        /* j0 *= \sigma_i^{-1}(j^i) */
        unity_zp_aut_inv(aut, ji, i);
        unity_zp_mul(temp, j0, aut);
        unity_zp_swap(temp, j0);

        /* jv *= \sigma_i^{-1}(j^{(v * i) / r}) */
        unity_zp_aut_inv(aut, je, i);
        unity_zp_mul(temp, jv, aut);
        unity_zp_swap(temp, jv);
    }

    /* temp = j0^u */
    unity_zp_pow_sliding_fmpz(temp, j0, u);
    /* j0 = j0^u * jv */
    unity_zp_mul(j0, jv, temp);

    /* try to find h */
    h = unity_zp_is_unity(j0);

    /* clear */
    unity_zp_clear(aut);
    unity_zp_clear(j0);
    unity_zp_clear(jv);
    unity_zp_clear(temp);
    unity_zp_clear(jr);
    unity_zp_clear(ji);
    unity_zp_clear(je);

    return h;
}

/*
    Check the case p = 2 and k = 1.

    Computes j0 = j_{0, 2, q}, jv = j_{v, 2, q} and checks that
    j0^u * jv is root of unity. j^0^u * jv = (-q)^{(n - 1) / 2}.

    Parameters:
        q, n from standard variables;

    Returns:
        if j0^u * jv = 1 returns 0;
        if j0^u * jv = -1 returns 1;
        otherwise returns -1.

    For details see algorithm (9.1.28) step 4.d in [1].
*/
slong
_aprcl_is_prime_jacobi_check_21(ulong q, const fmpz_t n)
{
    slong h;
    fmpz_t qpow, ndec, temp;

    /* initialization */
    fmpz_init(temp);
    fmpz_init_set_ui(qpow, q);
    fmpz_init_set(ndec, n);

    /* qpow = -q mod n */
    fmpz_sub(qpow, n, qpow);
    fmpz_sub_ui(ndec, ndec, 1);
    /* temp = (n - 1) / 2 */
    fmpz_fdiv_q_2exp(temp, ndec, 1);
    /* qpow = (-q)^{(n - 1) / 2} */
    fmpz_powm(qpow, qpow, temp, n);

    h = -1;
    /* check if qpow == +-1 */
    if (fmpz_equal_ui(qpow, 1))
        h = 0;
    if (fmpz_equal(qpow, ndec))
        h = 1;

    /* clear */
    fmpz_clear(temp);
    fmpz_clear(qpow);
    fmpz_clear(ndec);

    return h;
}

/*
    Check the case p = 2 and k = 2.

    Computes j0 = j_{0, 2, q}, jv = j_{v, 2, q} and checks that
    j0^u * jv is root of unity.

    Parameters:
        j = J(2, q);
        u, v, q from standard variables;

    Returns:
        If there exist h such that j0^u * jv = \zeta_{4}^h returns h
        \zeta_4^h \in (1, i, -1, -i);
        otherwise returns -1.

    For details see algorithm (9.1.28) step 4.c in [1].
*/
slong
_aprcl_is_prime_jacobi_check_22(const unity_zp j, const fmpz_t u, ulong v, ulong q)
{
    slong h;
    unity_zp j0, jv, j_pow;

    /* use mpn_mod arithmetic whenever the modulus is in range */
    if (_aprcl_is_prime_jacobi_check_22_mpn(&h, j, u, v, q))
        return h;

    /* initialization */
    unity_zp_init(j_pow, 2, 2, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(j0, 2, 2, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(jv, 2, 2, fmpz_mod_ctx_modulus(j->ctx));

    /* set j0 = q * j^2 */
    unity_zp_mul(j_pow, j, j);
    unity_zp_mul_scalar_ui(j0, j_pow, q);

    /*
        if v == 1 jv = 1
        if v == 3 jv = j^2
    */
    if (v == 1)
        unity_zp_coeff_set_ui(jv, 0, 1);
    else if (v == 3)
        unity_zp_swap(jv, j_pow);

    /* j0 = j0^u * jv */
    unity_zp_pow_sliding_fmpz(j_pow, j0, u);
    unity_zp_mul(j0, jv, j_pow);

    /* try to find h */
    h = unity_zp_is_unity(j0);

    /* clear */
    unity_zp_clear(j_pow);
    unity_zp_clear(j0);
    unity_zp_clear(jv);

    return h;
}


/*
    Check the case p = 2 and k >= 3.

    Computes j0 = j_{0, 2, q}, jv = j_{v, 2, q} and checks that
    j0^u * jv is root of unity.

    Parameters:
        j = J(2, q);
        j2_1 = J_3(q);
        j2_2 = j_2(q);

        u, v from standard variables;

    Returns:
        If there exist h such that j0^u * jv = \zeta_{2^k}^h returns h
        otherwise returns -1.

    For details see algorithm (9.1.28) step 4.b in [1].
*/
slong
_aprcl_is_prime_jacobi_check_2k(const unity_zp j, const unity_zp j2_1,
        const unity_zp j2_2, const fmpz_t u, ulong v)
{
    slong h;
    ulong i, r, e;
    unity_zp j_j0, j0, jv, temp, aut, ji, je, jj2, jj6;

    /* use mpn_mod arithmetic whenever the modulus is in range */
    if (_aprcl_is_prime_jacobi_check_2k_mpn(&h, j, j2_1, j2_2, u, v))
        return h;

    /* initialization */
    r = n_pow(j->p, j->exp);
    unity_zp_init(temp, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(j_j0, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(aut, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(j0, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(jv, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(ji, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(je, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(jj2, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));
    unity_zp_init(jj6, 2, j->exp, fmpz_mod_ctx_modulus(j->ctx));

    unity_zp_coeff_set_ui(j0, 0, 1);    /* j0 = 1 */
    unity_zp_coeff_set_ui(jv, 0, 1);    /* jv = 1 */

    /* j_j0 = J(2, q) * J_3(q) */
    unity_zp_mul(j_j0, j, j2_1);

    /*
        The loop runs over i = 1, 3, 9, 11, 17, 19, ... (i = 1 or 3 mod 8).
        The powers j_j0^i and j_j0^{(v * i) / r} are maintained
        incrementally: ji = j_j0^i is multiplied by jj2 = j_j0^2 or
        jj6 = j_j0^6 as i steps by 2 or 6, and je = j_j0^e picks up
        one factor of j_j0 each time e = (v * i) / r increases.
    */
    unity_zp_mul(jj2, j_j0, j_j0);
    unity_zp_mul(temp, jj2, jj2);
    unity_zp_mul(jj6, temp, jj2);
    unity_zp_copy(ji, j_j0);
    unity_zp_coeff_set_ui(je, 0, 1);
    e = 0;

    for (i = 1; i < r; i += (i % 8 == 1) ? 2 : 6)
    {
        /* je = j_j0^{(v * i) / r} */
        while (e < (v * i) / r)
        {
            unity_zp_mul(temp, je, j_j0);
            unity_zp_swap(temp, je);
            e++;
        }

        /* j0 *= \sigma_i^{-1}(j_j0^i) */
        unity_zp_aut_inv(aut, ji, i);
        unity_zp_mul(temp, j0, aut);
        unity_zp_swap(temp, j0);

        /* jv *= \sigma_i^{-1}(j_j0^{(v * i) / r}) */
        unity_zp_aut_inv(aut, je, i);
        unity_zp_mul(temp, jv, aut);
        unity_zp_swap(temp, jv);

        /* ji = j_j0^{next i} */
        unity_zp_mul(temp, ji, (i % 8 == 1) ? jj2 : jj6);
        unity_zp_swap(temp, ji);
    }

    /* if v % 8 not congruent 1 or 3 modulo 8 then jv *= J_2(q)^2 */
    if (v % 8 != 1 && v % 8 != 3)
    {
        unity_zp_mul(temp, j2_2, j2_2);
        unity_zp_mul(j_j0, jv, temp);
        unity_zp_swap(j_j0, jv);
    }

    /* set j0 = j0^u * jv */
    unity_zp_pow_sliding_fmpz(temp, j0, u);
    unity_zp_mul(j0, jv, temp);

    /* try to find h */
    h = unity_zp_is_unity(j0);

    /* clear */
    unity_zp_clear(aut);
    unity_zp_clear(j0);
    unity_zp_clear(jv);
    unity_zp_clear(j_j0);
    unity_zp_clear(temp);
    unity_zp_clear(ji);
    unity_zp_clear(je);
    unity_zp_clear(jj2);
    unity_zp_clear(jj6);

    return h;
}

/*
    Try to find prime number q such that:
        q == mod 2p;
        n^{(q - 1) / p} != 1 mod q;
        if p == 2 then 4 | q - 1 and 8 not | q - 1.

    If the q is found we have two cases:
        - if p == 2 verify
            1) h from above is a primitive root;
            2) q^{(n - 1) / 2} = -1 mod n;
        - if p != 2 verify
            1) h from above is a primitive root;
    if this is not the case for (p, q) then n is composite.

    If we can't find q check p | n or if n is a perfect power; if so
    then n is composite.
    Otherwise n can be prime but we can't prove its primality.

    Parameters:
        n, p from standard variables;

    Returns:
        if n composite : returns 2;
        if n can be prime : returns 1;
        if we do not find q and n can be prime : returns 0.

    For details see algorithm (9.1.28) step 5 in [1].
*/
int
_aprcl_is_prime_jacobi_additional_test(const fmpz_t n, ulong p)
{
    int result, p_counter, m;
    ulong q;
    fmpz_t npow, qmod;

    /* initialization */
    result = 0;
    p_counter = 50;     /* check only first 50 primes q = 2mp + 1 */
    m = 9;              /* begin from q = 2*9*p + 1 */

    fmpz_init(npow);
    fmpz_init(qmod);

    /* check first 50 primes */
    while (p_counter > 0)
    {
        /*
            q = 2mp + 1 and m is odd, so if p == 2
            then q - 1 | 4 and not q - 1 | 8.
        */
        q = 2 * m * p + 1;
        /* if q prime */
        if (n_is_prime(q) && fmpz_fdiv_ui(n, q) != 0)
        {
            fmpz_set_ui(qmod, q);
            /* npow = n^{(q - 1) / p} */
            fmpz_powm_ui(npow, n, (q - 1) / p, qmod);
            /* if n^{(q - 1) / p} != 1 mod n then we find q */
            if (!fmpz_equal_ui(npow, 1))
                break;
            /* else decrease the prime counter */
            p_counter--;
        }
        /* m is still odd */
        m += 2;
    }

    /* if we find q */
    if (p_counter != 0)
    {
        if (fmpz_fdiv_ui(n, q) == 0 && !fmpz_equal_ui(n, q))
           result = 2;
        else
        {
           ulong v, k;
           slong h;
           fmpz_t u;
           unity_zp jacobi_sum;

           fmpz_init(u);

           /* find max k such that p^k | q - 1; if p = 2 => k = 2 */
           k = aprcl_p_power_in_q(q - 1, p);

           /* compute J(p, q) */
           unity_zp_init(jacobi_sum, p, k, n);
           unity_zp_jacobi_sum_pq(jacobi_sum, q, p);

           /* compute u and v */
           fmpz_tdiv_q_ui(u, n, n_pow(p, k));
           v = fmpz_tdiv_ui(n, n_pow(p, k));

           /* if p == 2 */
           if (p == 2)
           {
               /* find h using p = 2, k = 2 */
               h = _aprcl_is_prime_jacobi_check_22(jacobi_sum, u, v, q);
               /* if h not found or h not primitive root then n is composite */
               if (h < 0 || h % 2 == 0)
                  result = 2;
               else  /* else verify condition for p = 2 */
               {
                  fmpz_t ndec, ndecdiv, qpow;

                   fmpz_init_set(ndec, n);
                   fmpz_init(ndecdiv);
                   fmpz_init_set_ui(qpow, q);

                   /* ndec = n - 1 */
                   fmpz_sub_ui(ndec, ndec, 1);
                   /* ndecdiv = (n - 1) / 2 */
                   fmpz_fdiv_q_2exp(ndecdiv, ndec, 1);
                   /* qpow = q^{(n - 1) / 2} */
                   fmpz_powm(qpow, qpow, ndecdiv, n);

                   /* if q^{(n - 1) / 2} = -1 mod n then n can be prime */
                   if (fmpz_equal(qpow, ndec))
                      result = 1;
                   else /* else n is composite */
                      result = 2;

                   fmpz_clear(ndec);
                   fmpz_clear(ndecdiv);
                   fmpz_clear(qpow);
               }
           }
           else  /* if p != 2 */
           {
              /* find h using (2.a) */
              h = _aprcl_is_prime_jacobi_check_pk(jacobi_sum, u, v);
              /* if h not found or h not primitive root then n is composite */
              if (h < 0 || h % p == 0)
                 result = 2;
              else /* else n can be prime */
                 result = 1;
           }

           fmpz_clear(u);
           unity_zp_clear(jacobi_sum);
       }
    }

    /* if we do not find a q then check (4.c) */
    if (p_counter == 0)
    {
       fmpz_t root;

       if (fmpz_tdiv_ui(n, p) == 0) /* if p | n then n is composite */
            result = 2;

       fmpz_init(root);

       if (fmpz_is_perfect_power(root, n)) /* if n is perfect power, composite */
          result = 2;

       fmpz_clear(root);
    }
    /* otherwise we can't prove composite or prime */

    /* clear */
    fmpz_clear(npow);
    fmpz_clear(qmod);

    return result;
}

/*
    One (q, p) pair of the Jacobi sum test. The pairs are independent, so
    they are run as a task list, possibly in parallel. Each task computes h
    (or -1 if n is proven composite) and, for p = 2, k >= 2, whether
    q^{(n - 1) / 2} = -1 mod n, which is needed for the (L_p) check.
*/
typedef struct
{
    ulong q;
    ulong p;
    ulong k;
    int pind;           /* index of p in config->rs */
    slong h;            /* result: h, or -1 if n is composite */
    int qpow_is_minus1; /* result: q^{(n - 1) / 2} = -1 mod n (p = 2, k >= 2) */
} _aprcl_jacobi_task;

typedef struct
{
    const fmpz * n;
    const fmpz * ndec;      /* n - 1 */
    const fmpz * ndecdiv;   /* (n - 1) / 2 */
    _aprcl_jacobi_task * tasks;
    slong num_tasks;
    slong next_task;
    int composite;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;
#endif
} _aprcl_jacobi_args;

static void
_aprcl_jacobi_run_task(_aprcl_jacobi_task * t, const _aprcl_jacobi_args * args)
{
    ulong p = t->p, k = t->k, q = t->q;
    ulong r = n_pow(p, k);
    ulong v;
    fmpz_t u;

    t->qpow_is_minus1 = 0;

    if (p == 2 && k == 1)
    {
        t->h = _aprcl_is_prime_jacobi_check_21(q, args->n);
        return;
    }

    fmpz_init(u);
    fmpz_tdiv_q_ui(u, args->n, r);      /* u = n / r */
    v = fmpz_tdiv_ui(args->n, r);       /* v = n % r */

    if (p == 2)
    {
        unity_zp jacobi_sum, jacobi_sum2_1, jacobi_sum2_2;
        fmpz_t q_pow;

        /* q_pow = q^{(n - 1) / 2} mod n, needed for the (L_2) check */
        fmpz_init_set_ui(q_pow, q);
        fmpz_powm(q_pow, q_pow, args->ndecdiv, args->n);
        t->qpow_is_minus1 = fmpz_equal(q_pow, args->ndec);
        fmpz_clear(q_pow);

        unity_zp_init(jacobi_sum, p, k, args->n);
        unity_zp_jacobi_sum_pq(jacobi_sum, q, p);

        if (k == 2)
        {
            t->h = _aprcl_is_prime_jacobi_check_22(jacobi_sum, u, v, q);
        }
        else
        {
            /* for k >= 3 we also need J_2(q) and J_3(q) */
            unity_zp_init(jacobi_sum2_1, p, k, args->n);
            unity_zp_init(jacobi_sum2_2, p, k, args->n);
            unity_zp_jacobi_sum_2q_one(jacobi_sum2_1, q);
            unity_zp_jacobi_sum_2q_two(jacobi_sum2_2, q);
            t->h = _aprcl_is_prime_jacobi_check_2k(jacobi_sum,
                                jacobi_sum2_1, jacobi_sum2_2, u, v);
            unity_zp_clear(jacobi_sum2_1);
            unity_zp_clear(jacobi_sum2_2);
        }

        unity_zp_clear(jacobi_sum);
    }
    else
    {
        unity_zp jacobi_sum;

        unity_zp_init(jacobi_sum, p, k, args->n);
        unity_zp_jacobi_sum_pq(jacobi_sum, q, p);
        t->h = _aprcl_is_prime_jacobi_check_pk(jacobi_sum, u, v);
        unity_zp_clear(jacobi_sum);
    }

    fmpz_clear(u);
}

/* worker: grab the next unclaimed task until none are left or n is composite */
static void
_aprcl_jacobi_worker(slong FLINT_UNUSED(i), void * varg)
{
    _aprcl_jacobi_args * args = (_aprcl_jacobi_args *) varg;

    while (1)
    {
        slong t;

#if FLINT_USES_PTHREAD
        pthread_mutex_lock(&args->mutex);
#endif
        t = (args->composite) ? args->num_tasks : args->next_task;
        if (t < args->num_tasks)
            args->next_task = t + 1;
#if FLINT_USES_PTHREAD
        pthread_mutex_unlock(&args->mutex);
#endif

        if (t >= args->num_tasks)
            break;

        _aprcl_jacobi_run_task(args->tasks + t, args);

        if (args->tasks[t].h < 0)
        {
#if FLINT_USES_PTHREAD
            pthread_mutex_lock(&args->mutex);
#endif
            args->composite = 1;
#if FLINT_USES_PTHREAD
            pthread_mutex_unlock(&args->mutex);
#endif
        }
    }
}

/* order tasks by decreasing phi(p^k), i.e. decreasing cost, for load balance */
static int
_aprcl_jacobi_task_cmp(const void * a, const void * b)
{
    const _aprcl_jacobi_task * x = (const _aprcl_jacobi_task *) a;
    const _aprcl_jacobi_task * y = (const _aprcl_jacobi_task *) b;
    ulong px = n_pow(x->p, x->k - 1) * (x->p - 1);
    ulong py = n_pow(y->p, y->k - 1) * (y->p - 1);

    if (px != py)
        return (px < py) ? 1 : -1;
    if (x->q != y->q)
        return (x->q < y->q) ? 1 : -1;
    return 0;
}

/*
    Steps (2.), (3.) and (4.) of the test. The pairs (q, p) with q | s and
    p | q - 1 are independent and are processed as a task list using all
    available threads; the (L_p) bookkeeping is done afterwards.
*/
primality_test_status
_aprcl_is_prime_jacobi(const fmpz_t n, const aprcl_config config)
{
    int * lambdas;
    ulong i, j, nmod4;
    slong num_tasks, t;
    primality_test_status result;
    fmpz_t temp, p2, ndec, ndecdiv;
    _aprcl_jacobi_task * tasks;
    _aprcl_jacobi_args args;

    /* deal with primes that can divide R */
    if (fmpz_cmp_ui(n, 2) == 0)
       return PRIME;
    if (fmpz_cmp_ui(n, 3) == 0)
       return PRIME;

    /* if n is one of the primes q, then n is prime */
    for (i = 0; i < config->qs->num; i++)
        if (config->qs_used[i] && fmpz_equal(n, config->qs->p + i))
            return PRIME;

    /* check that s*R and n are coprime; if not then n is composite */
    if (aprcl_is_mul_coprime_ui_fmpz(config->R, config->s, n) == 0)
        return COMPOSITE;

    /* initialization */
    fmpz_init(temp);
    fmpz_init(p2);
    fmpz_init(ndecdiv);
    fmpz_init_set(ndec, n);
    fmpz_sub_ui(ndec, ndec, 1);
    fmpz_fdiv_q_2exp(ndecdiv, ndec, 1);
    result = PROBABPRIME;

    /*
        Condition (Lp) is satisfied iff:
        For each p | R we must show that for all prime r | n and all
        positive integers a there exists l such that:
            r^(p-1) congruent n^(l(p-1)) mod p^a
        If (Lp), then lambdas_p = 1.
    */
    lambdas = (int*) flint_malloc(sizeof(int) * config->rs.num);

    /* nmod4 = n % 4 */
    nmod4 = fmpz_tdiv_ui(n, 4);

    /*
        For every prime p | R, set lambdas_p:
            to 1 if p >= 3 and n^{p - 1} != 1 mod p^2;
            to 0 otherwise.
    */
    for (i = 0; i < config->rs.num; i++)
    {
        ulong p = config->rs.p[i];
        lambdas[i] = 0;
        if (p > 2)
        {
            fmpz_set_ui(p2, p * p);
            fmpz_powm_ui(temp, n, p - 1, p2);
            if (fmpz_equal_ui(temp, 1) == 0)
                lambdas[i] = 1;
        }
    }

    /* Build the list of tasks: one for every prime q | s and prime p | q - 1 */
    num_tasks = 0;
    for (i = 0; i < config->qs->num; i++)
    {
        n_factor_t q_factors;

        if (config->qs_used[i] == 0)
            continue;

        n_factor_init(&q_factors);
        n_factor(&q_factors, fmpz_get_ui(config->qs->p + i) - 1, 1);
        num_tasks += q_factors.num;
    }

    tasks = (_aprcl_jacobi_task *) flint_malloc(sizeof(_aprcl_jacobi_task) * num_tasks);
    t = 0;
    for (i = 0; i < config->qs->num; i++)
    {
        n_factor_t q_factors;
        ulong q;

        if (config->qs_used[i] == 0)
            continue;

        q = fmpz_get_ui(config->qs->p + i); /* q must fit into an ulong */
        n_factor_init(&q_factors);
        n_factor(&q_factors, q - 1, 1);

        for (j = 0; j < q_factors.num; j++)
        {
            tasks[t].q = q;
            tasks[t].p = q_factors.p[j];
            tasks[t].k = q_factors.exp[j];
            tasks[t].pind = _aprcl_p_ind(config, q_factors.p[j]);
            tasks[t].h = 0;
            tasks[t].qpow_is_minus1 = 0;
            t++;
        }
    }

    qsort(tasks, num_tasks, sizeof(_aprcl_jacobi_task), _aprcl_jacobi_task_cmp);

    /* Begin pseudoprime tests with Jacobi sums step: run the tasks */
    args.n = n;
    args.ndec = ndec;
    args.ndecdiv = ndecdiv;
    args.tasks = tasks;
    args.num_tasks = num_tasks;
    args.next_task = 0;
    args.composite = 0;
#if FLINT_USES_PTHREAD
    pthread_mutex_init(&args.mutex, NULL);
#endif

    flint_parallel_do(_aprcl_jacobi_worker, &args,
            FLINT_MIN(num_tasks, flint_get_num_threads()), 0, FLINT_PARALLEL_STRIDED);

#if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&args.mutex);
#endif

    if (args.composite)
        result = COMPOSITE;

    /* Collect the (Lp) information from the tasks */
    if (result == PROBABPRIME)
    {
        for (t = 0; t < num_tasks; t++)
        {
            ulong p = tasks[t].p, k = tasks[t].k;
            slong h = tasks[t].h;
            int pind = tasks[t].pind;

            if (lambdas[pind] == 1)
                continue;

            if (p == 2 && k == 1)
            {
                /* if h == 1 (unity root = -1) and n % 4 == 1 then lambdas_2 = 1 */
                if (h == 1 && nmod4 == 1)
                    lambdas[pind] = 1;
            }
            else if (p == 2)
            {
                /*
                    if h is odd (primitive unity root)
                    and q^{(n - 1) / 2} = -1 mod n then lambdas_2 = 1
                */
                if (h % 2 != 0 && tasks[t].qpow_is_minus1)
                    lambdas[pind] = 1;
            }
            else
            {
                /* if h % p != 0 (primitive unity root) then lambdas_p = 1 */
                if (h % p != 0)
                    lambdas[pind] = 1;
            }
        }
    }

    flint_free(tasks);

    /* Begin L_p tests */
    /* if n can be prime */
    if (result == PROBABPRIME)
    {
        /* for every lambdas_p */
        for (i = 0; i < config->rs.num; i++)
        {
            /* if lambdas_p == 0 need run additional test for p */
            if (lambdas[i] == 0)
            {
                int r = _aprcl_is_prime_jacobi_additional_test(n, config->rs.p[i]);

                /* if r == 2 then we prove that n is composite */
                if (r == 2)
                {
                    result = COMPOSITE;
                }
                else if (r == 1)
                {
                    /* if r == 1 then we check that n still can be prime */
                    lambdas[i] = 1;
                }
                else
                {
                    /*
                        if r == 0 then we can't find q and we can't
                        prove n prime or composite
                    */
                    result = UNKNOWN;
                }
            }
        }
    }

    /* Trial division and primality proving step */
    /*
        If result = PROBAPRIME then all lambdas_p = 1.
        Using this information we know that:
            for every divisor p of f of n there exists i from 0..R-1 such that
            f == n^i mod s.
        This is done in aprcl_is_prime_final_division() function.
    */
    if (result == PROBABPRIME)
    {
        /* if we don't find divisors of n then n is composite */
        if (aprcl_is_prime_final_division(n, config->s, config->R) == 1)
            result = PRIME;
        else
            result = COMPOSITE;
    }

    /* clear */
    flint_free(lambdas);
    fmpz_clear(p2);
    fmpz_clear(ndec);
    fmpz_clear(ndecdiv);
    fmpz_clear(temp);

    return result;
}

int
aprcl_is_prime_jacobi(const fmpz_t n)
{
    primality_test_status result;
    aprcl_config config;

    /*
        Choose R and s values for n and store its factorisation.
        See definition in config_jacobi functions documentation. s^2 > n.
    */
    aprcl_config_jacobi_init(config, n);

    result = _aprcl_is_prime_jacobi(n, config);

    aprcl_config_jacobi_clear(config);

    if (result == PROBABPRIME || result == UNKNOWN)
    {
        flint_throw(FLINT_ERROR, "aprcl_is_prime_jacobi: failed to prove n prime for n = %s\n", fmpz_get_str(NULL, 10, n));
    }

    /* if we prove primality, returns 1 */
    if (result == PRIME)
        return 1;
    return 0;
}
