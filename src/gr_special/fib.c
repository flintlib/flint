/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "gr_special.h"

/* TODO: Binet formula + exp/log should be faster for huge
   numerical Fibonacci values */

/* TODO: single iteration applied directly to table; see fmpz_fib_ui */

/* TODO: for exact infinite rings, it is worth calling
   mpz_fib_ui instead of using generic code for n < 10000 or so.
   fmpz currently overrides the generic code, but would also
   help over QQ, etc. */

#if FLINT_BITS == 64
#define FIBTAB_NUM 94
#define FIBTAB_BITS 6
#else
#define FIBTAB_NUM 48
#define FIBTAB_BITS 5
#endif

static const ulong fibtab[] = {
    0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987,
    1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393,
    196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887,
    9227465, 14930352, 24157817, 39088169, 63245986, 102334155,
    165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903,
    UWORD(2971215073),
#if FLINT_BITS == 64
    UWORD(4807526976), UWORD(7778742049), UWORD(12586269025), UWORD(20365011074), UWORD(32951280099),
    UWORD(53316291173), UWORD(86267571272), UWORD(139583862445), UWORD(225851433717),
    UWORD(365435296162), UWORD(591286729879), UWORD(956722026041), UWORD(1548008755920),
    UWORD(2504730781961), UWORD(4052739537881), UWORD(6557470319842), UWORD(10610209857723),
    UWORD(17167680177565), UWORD(27777890035288), UWORD(44945570212853), UWORD(72723460248141),
    UWORD(117669030460994), UWORD(190392490709135), UWORD(308061521170129), UWORD(498454011879264),
    UWORD(806515533049393), UWORD(1304969544928657), UWORD(2111485077978050),
    UWORD(3416454622906707), UWORD(5527939700884757), UWORD(8944394323791464),
    UWORD(14472334024676221), UWORD(23416728348467685), UWORD(37889062373143906),
    UWORD(61305790721611591), UWORD(99194853094755497), UWORD(160500643816367088),
    UWORD(259695496911122585), UWORD(420196140727489673), UWORD(679891637638612258),
    UWORD(1100087778366101931), UWORD(1779979416004714189), UWORD(2880067194370816120),
    UWORD(4660046610375530309), UWORD(7540113804746346429), UWORD(12200160415121876738)
#endif
};

/* Essentially the same algorithm as
    https://gmplib.org/manual/Fibonacci-Numbers-Algorithm */

int
gr_generic_fib2_fmpz(gr_ptr v, gr_ptr u, const fmpz_t n, gr_ctx_t ctx)
{
    gr_method_unary_op sqr = GR_UNARY_OP(ctx, SQR);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    gr_method_binary_op sub = GR_BINARY_OP(ctx, SUB);
    gr_method_binary_op_ui add_ui = GR_BINARY_OP_UI(ctx, ADD_UI);
    gr_method_binary_op_ui sub_ui = GR_BINARY_OP_UI(ctx, SUB_UI);
    gr_method_binary_op_si mul_2exp_si = GR_BINARY_OP_SI(ctx, MUL_2EXP_SI);
    slong bit, nbits, i;
    int status = GR_SUCCESS;
    gr_ptr t;
    slong real_prec = 0;
    mp_srcptr np;
    mp_limb_t ntmp;

    if (fmpz_sgn(n) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, n);
        fmpz_add_ui(t, t, 1);

        status |= gr_generic_fib2_fmpz(u, v, t, ctx);
        if (fmpz_is_even(t))
            status |= gr_neg(u, u, ctx);
        else
            status |= gr_neg(v, v, ctx);

        fmpz_clear(t);
        return status;
    }

    nbits = fmpz_bits(n);

    if (nbits <= FIBTAB_BITS + 1 && *n < FIBTAB_NUM)
    {
        i = *n;
        status |= gr_set_ui(u, i ? fibtab[i - 1] : 1, ctx);
        status |= gr_set_ui(v, fibtab[i], ctx);
        return status;
    }

    bit = nbits - FIBTAB_BITS;

    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_tdiv_q_2exp(t, n, bit);
        i = *t;
        fmpz_clear(t);
    }

    if (!COEFF_IS_MPZ(*n))
    {
        ntmp = *n;
        np = &ntmp;
    }
    else
        np = COEFF_TO_PTR(*n)->_mp_d;

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        GR_IGNORE(gr_ctx_get_real_prec(&real_prec, ctx));
        GR_IGNORE(gr_ctx_set_real_prec(ctx, real_prec + 2 * nbits));
    }

    status |= gr_set_ui(u, i ? fibtab[i - 1] : 1, ctx);
    status |= gr_set_ui(v, fibtab[i], ctx);

    GR_TMP_INIT(t, ctx);

    while (1)
    {
        /* v = f(k), u = f(k-1) */
        status |= sqr(t, v, ctx);
        status |= sqr(v, u, ctx);

        /* f(2k-1) = f(k)^2 + f(k-1)^2 */
        status |= add(u, t, v, ctx);

        /* f(2k+1) = 4f(k)^2 - f(k-1)^2 + 2(-1)^k */
        if ((np[bit / FLINT_BITS] >> (bit % FLINT_BITS)) & 1)
            status |= add_ui(v, v, 2, ctx);
        else
            status |= sub_ui(v, v, 2, ctx);

        status |= mul_2exp_si(t, t, 2, ctx);
        status |= sub(v, t, v, ctx);

        bit--;

        /* restore precision before the final operation */
        if (bit == 0 && real_prec != 0)
            GR_IGNORE(gr_ctx_set_real_prec(ctx, real_prec));

        /* f(2k) = f(2k+1) - f(2k-1) */
        if ((np[bit / FLINT_BITS] >> (bit % FLINT_BITS)) & 1)
            status |= sub(u, v, u, ctx);
        else
            status |= sub(v, v, u, ctx);

        if (bit == 0)
            break;
    }

    GR_TMP_CLEAR(t, ctx);

    return status;
}

int
gr_generic_fib_fmpz(gr_ptr res, const fmpz_t n, gr_ctx_t ctx)
{
    fmpz_t n2;
    gr_ptr u, v;
    int status = GR_SUCCESS;
    slong real_prec = 0;

    if (fmpz_sgn(n) < 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_neg(t, n);

        status |= gr_generic_fib_fmpz(res, t, ctx);
        if (fmpz_is_even(t))
            status |= gr_neg(res, res, ctx);

        fmpz_clear(t);
        return status;
    }

    if (fmpz_cmp_ui(n, FIBTAB_NUM) < 0)
        return gr_set_ui(res, fibtab[*n], ctx);

    fmpz_init(n2);
    GR_TMP_INIT2(u, v, ctx);

    fmpz_tdiv_q_2exp(n2, n, 1);

    if (gr_ctx_has_real_prec(ctx) == T_TRUE)
    {
        GR_IGNORE(gr_ctx_get_real_prec(&real_prec, ctx));
        GR_IGNORE(gr_ctx_set_real_prec(ctx, real_prec + 10));
    }

    status = gr_generic_fib2_fmpz(v, u, n2, ctx);

    if (fmpz_is_odd(n))
    {
        /* f(2k+1) = (2f(k)+f(k-1))*(2f(k)-f(k-1)) + 2*(-1)^k */
        status |= gr_mul_two(v, v, ctx);
        status |= gr_add(res, v, u, ctx);
        status |= gr_sub(u, v, u, ctx);
        status |= gr_mul(v, res, u, ctx);

        /* restore precision before last operation */
        if (real_prec != 0)
            GR_IGNORE(gr_ctx_set_real_prec(ctx, real_prec));

        if (fmpz_is_odd(n2))
            status |= gr_sub_ui(res, v, 2, ctx);
        else
            status |= gr_add_ui(res, v, 2, ctx);
    }
    else
    {
        /* f(2k) = f(k)*(f(k)+2f(k-1)) */
        status |= gr_mul_two(u, u, ctx);
        status |= gr_add(u, u, v, ctx);

        /* restore precision before the final operation */
        if (real_prec != 0)
            GR_IGNORE(gr_ctx_set_real_prec(ctx, real_prec));

        status |= gr_mul(res, u, v, ctx);
    }

    fmpz_clear(n2);
    GR_TMP_CLEAR2(u, v, ctx);

    return status;
}

int
gr_generic_fib_ui(gr_ptr res, ulong n, gr_ctx_t ctx)
{
    fmpz_t t;
    int status;

    if (n < FIBTAB_NUM)
        return gr_set_ui(res, fibtab[n], ctx);

    fmpz_init_set_ui(t, n);
    status = gr_generic_fib_fmpz(res, t, ctx);
    fmpz_clear(t);
    return status;
}

int
gr_generic_fib_vec(gr_ptr res, slong len, gr_ctx_t ctx)
{
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    gr_method_unary_op_ui set_ui = GR_UNARY_OP_UI(ctx, SET_UI);
    int status = GR_SUCCESS;
    slong i;
    slong sz = ctx->sizeof_elem;

    for (i = 0; i < FLINT_MIN(FIBTAB_NUM, len); i++)
        status |= set_ui(GR_ENTRY(res, i, sz), fibtab[i], ctx);

    /* todo: consider increasing numerical precision */
    for (; i < len; i++)
        status |= add(GR_ENTRY(res, i, sz), GR_ENTRY(res, i - 1, sz), GR_ENTRY(res, i - 2, sz), ctx);

    return status;
}
