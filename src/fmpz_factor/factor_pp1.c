/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_factor.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

#define DEBUG 0 /* turn on some trace information */

#define pp1_mulmod(rxx, axx, bxx, nnn, nxx, ninv, norm)             \
   flint_mpn_mulmod_preinvn(rxx, axx, bxx, nnn, nxx, ninv, norm)

#ifdef FLINT64
static
ulong pp1_primorial[15] =
{
    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030), UWORD(510510), UWORD(9699690),
    UWORD(223092870), UWORD(6469693230), UWORD(200560490130), UWORD(7420738134810),
    UWORD(304250263527210), UWORD(13082761331670030), UWORD(614889782588491410)
};
#define num_primorials 15
#else
static
ulong pp1_primorial[9] =
{
    UWORD(2), UWORD(6), UWORD(30), UWORD(210), UWORD(2310), UWORD(30030), UWORD(510510), UWORD(9699690),
    UWORD(223092870)
};
#define num_primorials 9
#endif

void pp1_set(mp_ptr x1, mp_ptr y1,
              mp_srcptr x2, mp_srcptr y2, mp_size_t nn)
{
    flint_mpn_copyi(x1, x2, nn);
    flint_mpn_copyi(y1, y2, nn);
}

void pp1_set_ui(mp_ptr x, mp_size_t nn, ulong norm, ulong c)
{
    mpn_zero(x, nn);
    x[0] = (c << norm);
    if (nn > 1 && norm)
        x[1] = (c >> (FLINT_BITS - norm));
}

void pp1_2k(mp_ptr x, mp_ptr y, mp_size_t nn, mp_srcptr n,
            mp_srcptr ninv, mp_srcptr x0, ulong norm)
{
    pp1_mulmod(y, y, x, nn, n, ninv, norm);
    if (mpn_sub_n(y, y, x0, nn))
        mpn_add_n(y, y, n, nn);

    pp1_mulmod(x, x, x, nn, n, ninv, norm);
    if (mpn_sub_1(x, x, nn, UWORD(2) << norm))
        mpn_add_n(x, x, n, nn);
}

void pp1_2kp1(mp_ptr x, mp_ptr y, mp_size_t nn, mp_srcptr n,
              mp_srcptr ninv, mp_srcptr x0, ulong norm)
{
    pp1_mulmod(x, x, y, nn, n, ninv, norm);
    if (mpn_sub_n(x, x, x0, nn))
        mpn_add_n(x, x, n, nn);

    pp1_mulmod(y, y, y, nn, n, ninv, norm);
    if (mpn_sub_1(y, y, nn, UWORD(2) << norm))
        mpn_add_n(y, y, n, nn);
}

void pp1_pow_ui(mp_ptr x, mp_ptr y, mp_size_t nn,
                ulong exp, mp_srcptr n, mp_srcptr ninv, ulong norm)
{
    mp_limb_t t[30];
    mp_ptr x0 = t;
    ulong bit = ((UWORD(1) << FLINT_BIT_COUNT(exp)) >> 2);

    if (nn > 30)
        x0 = flint_malloc(nn*sizeof(mp_limb_t));
    flint_mpn_copyi(x0, x, nn);

    pp1_mulmod(y, x, x, nn, n, ninv, norm);
    if (mpn_sub_1(y, y, nn, UWORD(2) << norm))
        mpn_add_n(y, y, n, nn);

    while (bit)
    {
        if (exp & bit)
            pp1_2kp1(x, y, nn, n, ninv, x0, norm);
        else
            pp1_2k(x, y, nn, n, ninv, x0, norm);

        bit >>= 1;
    }

    if (nn > 30)
        flint_free(x0);
}

mp_size_t pp1_factor(mp_ptr factor, mp_srcptr n,
                     mp_srcptr x, mp_size_t nn, ulong norm)
{
    mp_size_t ret = 0, xn = nn;

    mp_ptr n2 = flint_malloc(nn*sizeof(mp_limb_t));
    mp_ptr x2 = flint_malloc(nn*sizeof(mp_limb_t));

    if (norm)
        mpn_rshift(n2, n, nn, norm);
    else
        flint_mpn_copyi(n2, n, nn);

    if (norm)
        mpn_rshift(x2, x, nn, norm);
    else
        flint_mpn_copyi(x2, x, nn);

    if (mpn_sub_1(x2, x2, nn, 2))
        mpn_add_n(x2, x2, n2, nn);

    MPN_NORM(x2, xn);

    if (xn == 0)
        goto cleanup;

    ret = flint_mpn_gcd_full(factor, n2, nn, x2, xn);

cleanup:

    flint_free(n2);
    flint_free(x2);

    return ret;
}

mp_size_t pp1_find_power(mp_ptr factor, mp_ptr x, mp_ptr y, mp_size_t nn,
                          ulong p, mp_srcptr n, mp_srcptr ninv, ulong norm)
{
    mp_size_t ret;

    do
    {
        pp1_pow_ui(x, y, nn, p, n, ninv, norm);
        ret = pp1_factor(factor, n, x, nn, norm);
    } while (ret == 1 && factor[0] == 1);

    return ret;
}

int fmpz_factor_pp1(fmpz_t fac, const fmpz_t n_in, ulong B1, ulong B2sqrt, ulong c)
{
    slong i, j;
    int ret = 0;
    mp_size_t nn = fmpz_size(n_in), r;
    mp_ptr x, y, oldx, oldy, n, ninv, factor, ptr_0, ptr_1, ptr_2, ptr_k;
    ulong pr, oldpr, sqrt, bits0, norm;
    n_primes_t iter;

    if (fmpz_is_even(n_in))
    {
        fmpz_set_ui(fac, 2);
        return 1;
    }

#if DEBUG
    flint_printf("starting stage 1\n");
#endif

    n_primes_init(iter);

    sqrt = n_sqrt(B1);
    bits0 = FLINT_BIT_COUNT(B1);

    x      = flint_malloc(nn*sizeof(mp_limb_t));
    y      = flint_malloc(nn*sizeof(mp_limb_t));
    oldx   = flint_malloc(nn*sizeof(mp_limb_t));
    oldy   = flint_malloc(nn*sizeof(mp_limb_t));
    n      = flint_malloc(nn*sizeof(mp_limb_t));
    ninv   = flint_malloc(nn*sizeof(mp_limb_t));
    factor = flint_malloc(nn*sizeof(mp_limb_t));

    if (nn == 1)
    {
        n[0] = fmpz_get_ui(n_in);
        norm = flint_clz(n[0]);
        n[0] <<= norm;
    } else
    {
        mp_ptr np = COEFF_TO_PTR(*n_in)->_mp_d;
        norm = flint_clz(np[nn - 1]);
        if (norm)
            mpn_lshift(n, np, nn, norm);
        else
            flint_mpn_copyi(n, np, nn);
    }

    flint_mpn_preinvn(ninv, n, nn);

    pp1_set_ui(x, nn, norm, c);

    /* mul by various prime powers */
    pr = 0;
    oldpr = 0;

    for (i = 0; pr < B1; )
    {
        j = i + 1024;
        oldpr = pr;
        pp1_set(oldx, oldy, x, y, nn);
        for ( ; i < j; i++)
        {
            pr = n_primes_next(iter);
            if (pr < sqrt)
            {
                ulong bits = FLINT_BIT_COUNT(pr);
                ulong exp = bits0 / bits;
                pp1_pow_ui(x, y, nn, n_pow(pr, exp), n, ninv, norm);
            } else
                pp1_pow_ui(x, y, nn, pr, n, ninv, norm);
        }

        r = pp1_factor(factor, n, x, nn, norm);
        if (r == 0)
            break;
        if (r != 1 || factor[0] != 1)
        {
            ret = 1;
            goto cleanup;
        }
    }

    if (pr < B1) /* factor = 0 */
    {
        n_primes_jump_after(iter, oldpr);
        pp1_set(x, y, oldx, oldy, nn);

        do
        {
            pr = n_primes_next(iter);
            pp1_set(oldx, oldy, x, y, nn);
            if (pr < sqrt)
            {
                ulong bits = FLINT_BIT_COUNT(pr);
                ulong exp = bits0 / bits;
                pp1_pow_ui(x, y, nn, n_pow(pr, exp), n, ninv, norm);
            } else
                pp1_pow_ui(x, y, nn, pr, n, ninv, norm);

            r = pp1_factor(factor, n, x, nn, norm);
            if (r == 0)
                break;
            if (r != 1 || factor[0] != 1)
            {
                ret = 1;
                goto cleanup;
            }
        } while (1);

        /* factor is still 0 */
        ret = pp1_find_power(factor, oldx, oldy, nn, pr, n, ninv, norm);
    } else /* stage 2 */
    {
        double quot;
        int num;
        char * sieve = flint_malloc(32768);
        slong * sieve_index = flint_malloc(32768*sizeof(slong));
        mp_ptr diff = flint_malloc(16384*nn*sizeof(mp_limb_t));
        ulong offset[15], num_roots;
        slong k, index = 0, s;
        fmpz * roots, * roots2, * evals;
        fmpz_poly_struct ** tree2;

#if DEBUG
        ulong primorial;
        flint_printf("starting stage 2\n");
#endif

        /* find primorial <= B2sqrt ... */
        for (num = 1; num < num_primorials; num++)
        {
            if (pp1_primorial[num] > B2sqrt)
                break;
        }
        num--;

        /* ... but not too big */
        quot = (double) B2sqrt / (double) pp1_primorial[num];
        if (quot < 1.1 && num > 0)
            num--;

#if DEBUG
        primorial = pp1_primorial[num];
        flint_printf("found primorial %wu\n", primorial);
#endif

        /* adjust B2sqrt to multiple of primorial */
        B2sqrt = (((B2sqrt - 1)/ pp1_primorial[num]) + 1) * pp1_primorial[num];

#if DEBUG
        flint_printf("adjusted B2sqrt %wu\n", B2sqrt);
#endif

        /* compute num roots */
        num++; /* number of primes is 1 more than primorial index */
        pr = 2;
        num_roots = B2sqrt;
        for (i = 0; i < num; i++)
        {
            num_roots = (num_roots*(pr - 1))/pr;
            pr = n_nextprime(pr, 0);
        }

#if DEBUG
        flint_printf("computed num_roots %wu\n", num_roots);
        flint_printf("B2 = %wu\n", num_roots * B2sqrt);
#endif

        /* construct roots */
        roots = _fmpz_vec_init(num_roots);
        for (i = 0; i < num_roots; i++)
        {
            __mpz_struct * m = _fmpz_promote(roots + i);
            mpz_realloc(m, nn);
        }

        roots2 = _fmpz_vec_init(num_roots);
        for (i = 0; i < num_roots; i++)
        {
            __mpz_struct * m = _fmpz_promote(roots2 + i);
            mpz_realloc(m, nn);
        }

        evals = _fmpz_vec_init(num_roots);

#if DEBUG
        flint_printf("constructed roots\n");
#endif

        /* compute differences table v0, ... */
        mpn_zero(diff, nn);
        diff[0] = (UWORD(2) << norm);

        /* ... v2, ... */
        pp1_mulmod(diff + nn, x, x, nn, n, ninv, norm);
        if (mpn_sub_1(diff + nn, diff + nn, nn, UWORD(2) << norm))
            mpn_add_n(diff + nn, diff + nn, n, nn);

        /* ... the rest ... v_{k+2} = v_k v_2 - v_{k-2} */
        k = 2*nn;
        for (i = 2; i < 16384; i++, k += nn)
        {
            pp1_mulmod(diff + k, diff + k - nn, diff + nn, nn, n, ninv, norm);
            if (mpn_sub_n(diff + k, diff + k, diff + k - 2*nn, nn))
                mpn_add_n(diff + k, diff + k, n, nn);
        }

#if DEBUG
        flint_printf("conputed differences table\n");
#endif

        /* initial positions */
        pr = 2;
        for (i = 0; i < num; i++)
        {
            offset[i] = pr/2;
            pr = n_nextprime(pr, 0);
        }

        s = 0;
        while (2*s + 1 < B2sqrt)
        {
            /* sieve */
            memset(sieve, 1, 32768);
            pr = 3;
            for (i = 1; i < num; i++)
            {
                j = offset[i];
                while (j < 32768)
                    sieve[j] = 0, j += pr;

                /* store offset for start of next sieve run */
                offset[i] = j - 32768;
                pr = n_nextprime(pr, 0);
            }

            /* compute roots */
            for (i = 0; i < 32768 && 2*(s + i) + 1 < B2sqrt; i++)
            {
                if (sieve[i])
                {
                    ptr_2 = COEFF_TO_PTR(roots[index])->_mp_d;
                    k = (i + 1)/2;
                    for (j = i - 1; j >= k; j--)
                    {
                        if (sieve[j] && sieve[2*j - i])
                        {
                            /* V_{n+k} = V_n V_k - V_{n-k} */
                            ptr_0 = COEFF_TO_PTR(roots[sieve_index[2*j - i]])->_mp_d;
                            ptr_1 = COEFF_TO_PTR(roots[sieve_index[j]])->_mp_d;
                            ptr_k = diff + (i - j)*nn;
                            pp1_mulmod(ptr_2, ptr_1, ptr_k, nn, n, ninv, norm);
                            if (mpn_sub_n(ptr_2, ptr_2, ptr_0, nn))
                                mpn_add_n(ptr_2, ptr_2, n, nn);
                            break;
                        }
                    }

                    if (j < k) /* pair not found, compute using pow_ui */
                    {
                        flint_mpn_copyi(ptr_2, x, nn);
                        pp1_pow_ui(ptr_2, y, nn, 2*(s + i) + 1, n, ninv, norm);
                    }

                    sieve_index[i] = index;
                    index++;
                }
            }

            s += 32768;
        }

#if DEBUG
        flint_printf("roots computed %wd\n", index);
#endif

        /* v_1 */
        flint_mpn_copyi(oldx, x, nn);
        pp1_pow_ui(oldx, y, nn, B2sqrt, n, ninv, norm);
        ptr_0 = COEFF_TO_PTR(roots2[0])->_mp_d;
        flint_mpn_copyi(ptr_0, oldx, nn);

        /* v_2 */
        ptr_1 = COEFF_TO_PTR(roots2[1])->_mp_d;
        pp1_mulmod(ptr_1, ptr_0, ptr_0, nn, n, ninv, norm);
        if (mpn_sub_1(ptr_1, ptr_1, nn, UWORD(2) << norm))
            mpn_add_n(ptr_1, ptr_1, n, nn);

        for (i = 2; i < num_roots; i++)
        {
            /* V_{k+n} = V_k V_n - V_{k-n} */
            ptr_2 = COEFF_TO_PTR(roots2[i])->_mp_d;

            pp1_mulmod(ptr_2, ptr_1, oldx, nn, n, ninv, norm);
            if (mpn_sub_n(ptr_2, ptr_2, ptr_0, nn))
                mpn_add_n(ptr_2, ptr_2, n, nn);

            ptr_0 = ptr_1;
            ptr_1 = ptr_2;
        }

#if DEBUG
        flint_printf("roots2 computed %wu\n", num_roots);
#endif

        for (i = 0; i < num_roots; i++)
        {
            mp_size_t sn;
            __mpz_struct * m1 = COEFF_TO_PTR(roots[i]);
            __mpz_struct * m2 = COEFF_TO_PTR(roots2[i]);

            ptr_1 = m1->_mp_d;
            ptr_2 = m2->_mp_d;

            if (norm)
            {
                mpn_rshift(ptr_1, ptr_1, nn, norm);
                mpn_rshift(ptr_2, ptr_2, nn, norm);
            }

            sn = nn;
            MPN_NORM(ptr_1, sn);
            m1->_mp_size = sn;

            sn = nn;
            MPN_NORM(ptr_2, sn);
            m2->_mp_size = sn;

            _fmpz_demote_val(roots + i);
            _fmpz_demote_val(roots2 + i);
        }

#if DEBUG
        flint_printf("normalised roots\n");
#endif

        fmpz_mod_ctx_t ctx;
        fmpz_mod_ctx_init(ctx, n_in);

        tree2 = _fmpz_mod_poly_tree_alloc(num_roots);
        _fmpz_mod_poly_tree_build(tree2, roots2, num_roots, ctx);

        /* todo: use fmpz_mod_poly_mul */
        fmpz_poly_mul(tree2[FLINT_CLOG2(num_roots)], tree2[FLINT_CLOG2(num_roots)-1], tree2[FLINT_CLOG2(num_roots)-1]+1);
        fmpz_poly_scalar_mod_fmpz(tree2[FLINT_CLOG2(num_roots)], tree2[FLINT_CLOG2(num_roots)], n_in);

#if DEBUG
        flint_printf("built trees\n");
#endif

        _fmpz_mod_poly_evaluate_fmpz_vec_fast(evals, tree2[FLINT_CLOG2(num_roots)]->coeffs, tree2[FLINT_CLOG2(num_roots)]->length, roots, num_roots, ctx);
        _fmpz_mod_poly_tree_free(tree2, num_roots);

        fmpz_mod_ctx_clear(ctx);

#if DEBUG
        flint_printf("evaluated at roots\n");
#endif

        for (i = 0; i < num_roots; i++)
        {
            fmpz_gcd(fac, n_in, evals + i);
            if (!fmpz_is_zero(fac) && !fmpz_is_one(fac))
            {
                ret = 1;
                break;
            }
        }

        _fmpz_vec_clear(evals, num_roots);
        _fmpz_vec_clear(roots, num_roots);
        _fmpz_vec_clear(roots2, num_roots);
        flint_free(sieve);
        flint_free(sieve_index);
        flint_free(diff);

        if (i < num_roots)
            goto cleanup2;
    }

#if DEBUG
    flint_printf("done stage2\n");
#endif

cleanup:

    if (ret)
    {
        __mpz_struct * fm = _fmpz_promote(fac);
        mpz_realloc(fm, r);
        flint_mpn_copyi(fm->_mp_d, factor, r);
        fm->_mp_size = r;
        _fmpz_demote_val(fac);
    }

cleanup2:

    flint_free(x);
    flint_free(y);
    flint_free(oldx);
    flint_free(oldy);
    flint_free(n);
    flint_free(ninv);
    flint_free(factor);

    n_primes_clear(iter);

    return ret;
}
