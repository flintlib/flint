/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "gr_special.h"
#include "thread_support.h"

#ifdef __GNUC__
# define sqrt __builtin_sqrt
#else
# include <math.h>
#endif

typedef struct
{
    arb_t P;
    arb_t R;
    arb_t Q;
    slong a;
    slong b;
}
bsplit_struct;

typedef bsplit_struct bsplit_t[1];

typedef struct
{
    slong N;
    slong prec;
}
bsplit_args_struct;

static void
bsplit_init(bsplit_t x, void * args)
{
    arb_init(x->P);
    arb_init(x->R);
    arb_init(x->Q);
    x->a = x->b = 0;
}

static void
bsplit_clear(bsplit_t x, void * args)
{
    arb_clear(x->P);
    arb_clear(x->R);
    arb_clear(x->Q);
}

static void
bsplit_basecase(bsplit_t res, slong n, slong n1, void * args)
{
    fmpz_t f1, f2, f1f2, f22, f4, t;

    FLINT_ASSERT(n1 == n + 1);

    fmpz_init(f1); fmpz_init(f2); fmpz_init(f1f2);
    fmpz_init(f22); fmpz_init(f4); fmpz_init(t);

    /* We compute the Fibonacci numbers needed in each basecase node from
       scratch. They could be recycled from each node to the next, but
       the speedup would be marginal. */

    /* f1, f2 = fib(2*n), fib(2*n+1) */
    {
        fmpz nn = 2 * n + 1;
        gr_ctx_t ctx;
        gr_ctx_init_fmpz(ctx);
        GR_MUST_SUCCEED(gr_generic_fib2_fmpz(f2, f1, &nn, ctx));
    }

    fmpz_mul(f1f2, f1, f2);
    fmpz_mul(f22, f2, f2);

    /* f4 = fib(4n+3) */
    fmpz_add(f4, f1f2, f22);
    fmpz_mul_2exp(f4, f4, 1);
    fmpz_addmul(f4, f1, f1);

    /* P = 2*f1 + f2 */
    fmpz_mul_2exp(t, f1, 1);
    fmpz_add(t, t, f2);
    arb_set_fmpz(res->P, t);

    /* Q = (f1f2 + f22) * P */
    fmpz_add(f1f2, f1f2, f22);
    fmpz_mul(t, t, f1f2);
    arb_set_fmpz(res->Q, t);

    /* R = (-1)^((n & 3)>>1) * (f4 + (-1)^n * (f1 + f2)) */
    fmpz_add(t, f1, f2);
    if (n & 1)
        fmpz_sub(t, f4, t);
    else
        fmpz_add(t, f4, t);
    if ((n & 3) >> 1)
        fmpz_neg(t, t);
    arb_set_fmpz(res->R, t);

    res->a = n;
    res->b = n1;

    fmpz_clear(f1); fmpz_clear(f2); fmpz_clear(f1f2);
    fmpz_clear(f22); fmpz_clear(f4); fmpz_clear(t);
}

static void
bsplit_merge(bsplit_t res, bsplit_t X1, bsplit_t X2, bsplit_args_struct * args)
{
    slong prec = args->prec;
    slong a = X1->a;
    slong b = X2->b;

    FLINT_ASSERT(res == X1);

    /*
        P1, R1, Q1 = bsplit(a, m)
        P2, R2, Q2 = bsplit(m, b)
        Q2 *= P1; P1 *= P2; R1 *= Q2; R2 *= Q1; R1 += R2; Q1 *= Q2
    */

    arb_mul(X2->Q, X2->Q, X1->P, prec);

    if (b != args->N)  /* this product is only needed internally */
        arb_mul(X1->P, X1->P, X2->P, prec);
    else
        arb_one(X1->P);

    arb_mul(X1->R, X1->R, X2->Q, prec);
    arb_mul(X2->R, X2->R, X1->Q, prec);
    arb_add(X1->R, X1->R, X2->R, prec);
    arb_mul(X1->Q, X1->Q, X2->Q, prec);

    X1->a = a;  /* actually a no-op because of aliasing */
    X1->b = b;
}

void
arb_const_reciprocal_fibonacci(arb_t res, slong prec)
{
    slong N, max_threads, wp;
    bsplit_t s;
    bsplit_args_struct args;

    max_threads = (prec < 30000) ? 1 : 0;
    args.prec = wp = prec + 10 + FLINT_BIT_COUNT(prec);
    /* Truncation error < (1/phi)^(N^2) < 2^(-wp) */
    args.N = N = (1.200175024907849 + 1e-12) * sqrt(wp) + 2;

    bsplit_init(s, NULL);

    flint_parallel_binary_splitting(s,
        (bsplit_basecase_func_t) bsplit_basecase,
        (bsplit_merge_func_t) bsplit_merge,
        sizeof(bsplit_struct),
        (bsplit_init_func_t) bsplit_init,
        (bsplit_clear_func_t) bsplit_clear,
        &args, 0, N, 1, max_threads, FLINT_PARALLEL_BSPLIT_LEFT_INPLACE);

    arb_div(res, s->R, s->Q, prec);
    arb_add_error_2exp_si(res, -wp);

    bsplit_clear(s, NULL);
}
