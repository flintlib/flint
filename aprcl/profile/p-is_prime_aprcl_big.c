/*
    Copyright 2015 Vladimir Glazachev

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "profiler.h"
#include "aprcl.h"

void p_is_prime_repeat(fmpz_t n)
{
    TIMEIT_START
    if (aprcl_is_prime(n) == 0)
    {
        flint_printf("Primality test failed\n");
        flint_abort();
    }
    TIMEIT_STOP
}

int main(void)
{
    FLINT_TEST_INIT(state);

    /*
        Using the primes from
        mpz_aprcl implementation readme.txt file
        link : https://sourceforge.net/projects/mpzaprcl/
    */

    flint_printf("Primality test profiling for numbers from 350 to 2000 digits\n");
    flint_printf("All timings given for one number\n");

    /* 350 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 1160);
        fmpz_pow_ui(n, n, 114);
        fmpz_add_ui(n, n, 7);

        flint_printf("350 digit prime : 1160^114 + 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 400 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 291);
        fmpz_pow_ui(n, n, 163);
        fmpz_sub_ui(n, n, 1);
        fmpz_fdiv_q_ui(n, n, 290);

        flint_printf("400 digit prime : (291^163 - 1) / 290\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 450 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 232);
        fmpz_pow_ui(n, n, 190);
        fmpz_add_ui(n, n, 7);

        flint_printf("450 digit prime : 232^190 + 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 500 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 1014);
        fmpz_pow_ui(n, n, 166);
        fmpz_add_ui(n, n, 7);

        flint_printf("500 digit prime : 1014^166 + 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 550 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 549);
        fmpz_mul_ui(n, n, 9);
        fmpz_sub_ui(n, n, 7);

        flint_printf("550 digit prime : 10^549 * 9 - 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 600 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 1432);
        fmpz_pow_ui(n, n, 190);
        fmpz_add_ui(n, n, 7);

        flint_printf("600 digit prime : 1432^190 + 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 650 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 2159);
        fmpz_add_ui(n, n, 375);

        flint_printf("650 digit prime : 2^2159 + 375\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 700 digits prime */
    {
        fmpz_t n1, n2;

        fmpz_init(n1);
        fmpz_init(n2);

        fmpz_set_ui(n1, 157);
        fmpz_pow_ui(n1, n1, 319);

        fmpz_set_ui(n2, 319);
        fmpz_pow_ui(n2, n2, 157);

        fmpz_add(n1, n1, n2);
        fmpz_fdiv_q_ui(n1, n1, 28);

        flint_printf("700 digit prime : (157^319 + 319^157) / 28\n");
        p_is_prime_repeat(n1);

        fmpz_clear(n1);
        fmpz_clear(n2);
    }

    /* 750 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 749);
        fmpz_mul_ui(n, n, 2);
        fmpz_add_ui(n, n, 89);

        flint_printf("750 digit prime : 10^749 * 2 + 89\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 800 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 799);
        fmpz_mul_ui(n, n, 61);
        fmpz_sub_ui(n, n, 7);
        fmpz_fdiv_q_ui(n, n, 9);

        flint_printf("800 digit prime : (10^799 * 61 - 7) / 9\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 850 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 2821);
        fmpz_sub_ui(n, n, 183);

        flint_printf("850 digit prime : 2^2821 - 183\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 900 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 24);
        fmpz_pow_ui(n, n, 653);
        fmpz_sub_ui(n, n, 1);
        fmpz_fdiv_q_ui(n, n, 23);

        flint_printf("900 digit prime : (24^653 - 1) / 23\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 950 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 949);
        fmpz_mul_ui(n, n, 4);
        fmpz_sub_ui(n, n, 9);

        flint_printf("950 digit prime : 10^949 * 4 - 9\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1000 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 999);
        fmpz_add_ui(n, n, 7);

        flint_printf("1000 digit prime : 10^999 + 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1100 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 3653);
        fmpz_add_ui(n, n, 41);

        flint_printf("1100 digit prime : 2^3653 + 41\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1200 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 1199);
        fmpz_mul_ui(n, n, 5);
        fmpz_add_ui(n, n, 9);

        flint_printf("1200 digit prime : 10^1199 * 5 + 9\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1300 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 4318);
        fmpz_add_ui(n, n, 165);

        flint_printf("1300 digit prime : 2^4318 + 165\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1400 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 187);
        fmpz_pow_ui(n, n, 617);
        fmpz_sub_ui(n, n, 1);
        fmpz_fdiv_q_ui(n, n, 186);

        flint_printf("1400 digit prime : (187^617 - 1) / 186\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1500 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 4972);
        fmpz_mul_ui(n, n, 1779);
        fmpz_sub_ui(n, n, 1);

        flint_printf("1500 digit prime : 2^4972 * 1779 - 1\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1600 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 12);
        fmpz_pow_ui(n, n, 1483);
        fmpz_add_ui(n, n, 1);
        fmpz_fdiv_q_ui(n, n, 13);

        flint_printf("1600 digit prime : (12^1483 + 1) / 13\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1700 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 5644);
        fmpz_sub_ui(n, n, 227);

        flint_printf("1700 digit prime : 2^5644 - 227\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1800 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 1800);
        fmpz_sub_ui(n, n, 87);

        flint_printf("1800 digit prime : 10^1800 - 87\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 1900 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 10);
        fmpz_pow_ui(n, n, 1900);
        fmpz_add_ui(n, n, 3);
        fmpz_fdiv_q_ui(n, n, 7);

        flint_printf("1900 digit prime : (10^1900 + 3) / 7\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    /* 2000 digits prime */
    {
        fmpz_t n;

        fmpz_init(n);
        fmpz_set_ui(n, 2);
        fmpz_pow_ui(n, n, 6643);
        fmpz_mul_ui(n, n, 113);
        fmpz_add_ui(n, n, 1);
        fmpz_fdiv_q_ui(n, n, 115);

        flint_printf("2000 digit prime : (2^6643 * 113 + 1) / 115\n");
        p_is_prime_repeat(n);

        fmpz_clear(n);
    }

    return 0;
}

