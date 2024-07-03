/*
    Copyright 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <ulong_extras.h>
#include <stdlib.h>  // for atoi

#include "profiler.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_mat.h"
#include "nmod_poly.h"
#include "gr_poly.h"

// utility (nmod vec uniform random)
static inline
void _nmod_vec_rand(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
{
    for (slong i = 0; i < len; i++)
        vec[i] = n_randint(state, mod.n);
}

// uniform (nmod mat uniform random)
static inline
void nmod_mat_rand(nmod_mat_t mat, flint_rand_t state)
{
    _nmod_vec_rand(mat->entries, state, mat->r * mat->c, mat->mod);
}

/*------------------------------------*/
/* direct: dot / dot_rev / dot expr   */
/*------------------------------------*/

// timings excluding dot_params
void time_dot(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile ulong FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod_vec_dot(v1, v2, len, mod, params);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

void time_dot_rev(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile ulong FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod_vec_dot_rev(v1, v2, len, mod, params);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

void time_dot_ptr(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    const dot_params_t params = _nmod_vec_dot_params(len, mod);

    const ulong offset = UWORD(7);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2tmp = _nmod_vec_init(len);
    _nmod_vec_rand(v2tmp, state, len, mod);
    nn_ptr * v2 = flint_malloc(sizeof(nn_ptr) * len);
    for (ulong i = 0; i < len; i++)
        v2[i] = &v2tmp[i] + offset;

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile ulong FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    res = _nmod_vec_dot_ptr(v1, v2, offset, len, mod, params);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2tmp);
    flint_free(v2);
}

// timings including dot_params
void time_dot_incparams(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile ulong FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    const dot_params_t params = _nmod_vec_dot_params(len, mod);
    res = _nmod_vec_dot(v1, v2, len, mod, params);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

void time_dot_rev_incparams(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile ulong FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    const dot_params_t params = _nmod_vec_dot_params(len, mod);
    res = _nmod_vec_dot_rev(v1, v2, len, mod, params);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

void time_dot_ptr_incparams(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);

    const ulong offset = UWORD(7);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2tmp = _nmod_vec_init(len);
    _nmod_vec_rand(v2tmp, state, len, mod);
    nn_ptr * v2 = flint_malloc(sizeof(nn_ptr) * len);
    for (ulong i = 0; i < len; i++)
        v2[i] = &v2tmp[i] + offset;

    // store results in volatile variable to avoid that they
    // are "optimized away" (especially for inlined part)
    volatile ulong FLINT_SET_BUT_UNUSED(res);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    const dot_params_t params = _nmod_vec_dot_params(len, mod);
    res = _nmod_vec_dot_ptr(v1, v2, offset, len, mod, params);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2tmp);
    flint_free(v2);
}

/*-------------------------*/
/* indirect: poly          */
/*-------------------------*/

void time_dot_poly_mul(ulong len, ulong n, flint_rand_t state)
{
    if (len > 10000)
    {
        printf("        ");
        return;
    }

    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr p1 = _nmod_vec_init(len);
    _nmod_vec_rand(p1, state, len, mod);
    nn_ptr p2 = _nmod_vec_init(len);
    _nmod_vec_rand(p2, state, len, mod);
    nn_ptr res = _nmod_vec_init(2*len);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    _nmod_poly_mul_classical(res, p1, len, p2, len, mod);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(p1);
    _nmod_vec_clear(p2);
    _nmod_vec_clear(res);
}

void time_dot_poly_inv_series(ulong len, ulong n, flint_rand_t state)
{
    if (len > 10000 || n % 2 == 0)
    {
        printf("        ");
        return;
    }

    nmod_t mod;
    nmod_init(&mod, n);

    nn_ptr p = _nmod_vec_init(len);
    _nmod_vec_rand(p, state, len, mod);
    nn_ptr res = _nmod_vec_init(len);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    _nmod_poly_inv_series_basecase(res, p, len, len, mod);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(p);
    _nmod_vec_clear(res);
}

void time_dot_poly_exp_series(ulong len, ulong n, flint_rand_t state)
{
    if (len > 10000 || n % 2 == 0)
    {
        printf("        ");
        return;
    }

    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, n);
    int FLINT_SET_BUT_UNUSED(status);

    gr_poly_t p;
    gr_poly_init(p, ctx);
    gr_poly_randtest(p, state, len, ctx);
    status = gr_poly_set_coeff_si(p, 0, WORD(0), ctx);
    gr_poly_t res;
    gr_poly_init(res, ctx);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    status |= gr_poly_exp_series_basecase(res, p, len, ctx);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    gr_poly_clear(p, ctx);
    gr_poly_clear(res, ctx);
}

/*-------------------------*/
/* indirect: mat           */
/*-------------------------*/

void time_dot_mat_mul(ulong len, ulong n, flint_rand_t state)
{
    if (len > 2000)
    {
        printf("        ");
        return;
    }

    nmod_mat_t mat1; nmod_mat_init(mat1, len, len, n); nmod_mat_rand(mat1, state);
    nmod_mat_t mat2; nmod_mat_init(mat2, len, len, n); nmod_mat_rand(mat2, state);
    nmod_mat_t mat; nmod_mat_init(mat, len, len, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    _nmod_mat_mul_classical_op(mat, mat, mat1, mat2, 1);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    nmod_mat_clear(mat);
    nmod_mat_clear(mat2);
    nmod_mat_clear(mat1);
}

void time_dot_mat_mul_vec(ulong len, ulong n, flint_rand_t state)
{
    if (len > 10000)
    {
        printf("        ");
        return;
    }

    nmod_t mod;
    nmod_init(&mod, n);

    nmod_mat_t mat; nmod_mat_init(mat, len, len, n); nmod_mat_rand(mat, state);
    nn_ptr v; v = _nmod_vec_init(len); _nmod_vec_rand(v, state, len, mod);
    nn_ptr u; u = _nmod_vec_init(len);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    nmod_mat_mul_nmod_vec(u, mat, v, len);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    nmod_mat_clear(mat);
    _nmod_vec_clear(u);
    _nmod_vec_clear(v);
}

void time_dot_mat_solve_tril(ulong len, ulong n, flint_rand_t state)
{
    if (len > 4000 || n % 2 == 0)
    {
        printf("        ");
        return;
    }

    nmod_mat_t mat1; nmod_mat_init(mat1, len, len, n); nmod_mat_randtril(mat1, state, 0);
    nmod_mat_t mat2; nmod_mat_init(mat2, len, len, n); nmod_mat_rand(mat2, state);
    nmod_mat_t mat; nmod_mat_init(mat, len, len, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    nmod_mat_solve_tril_classical(mat, mat1, mat2, 0);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    nmod_mat_clear(mat);
    nmod_mat_clear(mat2);
    nmod_mat_clear(mat1);
}

void time_dot_mat_solve_tril_vec(ulong len, ulong n, flint_rand_t state)
{
    if (len > 10000 || n % 2 == 0)
    {
        printf("        ");
        return;
    }

    nmod_mat_t mat1; nmod_mat_init(mat1, len, len, n); nmod_mat_randtril(mat1, state, 0);
    nmod_mat_t mat2; nmod_mat_init(mat2, len, 1, n); nmod_mat_rand(mat2, state);
    nmod_mat_t mat; nmod_mat_init(mat, len, len, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    nmod_mat_solve_tril_classical(mat, mat1, mat2, 0);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    nmod_mat_clear(mat);
    nmod_mat_clear(mat2);
    nmod_mat_clear(mat1);
}

void time_dot_mat_solve_triu(ulong len, ulong n, flint_rand_t state)
{
    if (len > 4000 || n % 2 == 0)
    {
        printf("        ");
        return;
    }

    nmod_mat_t mat1; nmod_mat_init(mat1, len, len, n); nmod_mat_randtriu(mat1, state, 0);
    nmod_mat_t mat2; nmod_mat_init(mat2, len, len, n); nmod_mat_rand(mat2, state);
    nmod_mat_t mat; nmod_mat_init(mat, len, len, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    nmod_mat_solve_triu_classical(mat, mat1, mat2, 0);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    nmod_mat_clear(mat);
    nmod_mat_clear(mat2);
    nmod_mat_clear(mat1);
}

void time_dot_mat_solve_triu_vec(ulong len, ulong n, flint_rand_t state)
{
    if (len > 10000 || n % 2 == 0)
    {
        printf("        ");
        return;
    }

    nmod_mat_t mat1; nmod_mat_init(mat1, len, len, n); nmod_mat_randtriu(mat1, state, 0);
    nmod_mat_t mat2; nmod_mat_init(mat2, len, 1, n); nmod_mat_rand(mat2, state);
    nmod_mat_t mat; nmod_mat_init(mat, len, len, n);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    nmod_mat_solve_triu_classical(mat, mat1, mat2, 0);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    nmod_mat_clear(mat);
    nmod_mat_clear(mat2);
    nmod_mat_clear(mat1);
}


/*-------------------------*/
/*  main                   */
/*-------------------------*/

int main(int argc, char ** argv)
{
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, time(NULL), time(NULL)+129384125L);

    // modulus bitsize
    const slong nbits = 14;
    const ulong bits[] = {12, 28, 30, 31, 32, 40, 50, 60, 61, 62, 63, 64, 232, 263};

    // vector lengths
    const slong nlens = 20;
    const ulong lens[] = {1, 2, 3, 4, 5, 7, 10, 15, 25, 35, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 100000, 1000000};

    // bench functions
    const slong nfuns = 15;
    typedef void (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_dot,                    // 0
        time_dot_incparams,          // 1
        time_dot_rev,                // 2
        time_dot_rev_incparams,      // 3
        time_dot_ptr,                // 4
        time_dot_ptr_incparams,      // 5
        time_dot_poly_mul,           // 6
        time_dot_poly_inv_series,    // 7
        time_dot_poly_exp_series,    // 8
        time_dot_mat_mul,            // 9
        time_dot_mat_solve_tril,     // 10
        time_dot_mat_solve_triu,     // 11
        time_dot_mat_mul_vec,        // 12
        time_dot_mat_solve_tril_vec, // 13
        time_dot_mat_solve_triu_vec, // 14
    };

    const char * description[] = {
        "#0  --> vec dot               ",
        "#1  --> vec dot inc params    ",
        "#2  --> vec dot rev           ",
        "#3  --> vec dot rev inc params",
        "#4  --> vec dot ptr           ",
        "#5  --> vec dot ptr inc params",
        "#6  --> poly_mul              ",
        "#7  --> poly_inv_series       ",
        "#8  --> poly_exp_series       ",
        "#9  --> mat_mul               ",
        "#10 --> mat_solve_tril        ",
        "#11 --> mat_solve_triu        ",
        "#12 --> mat_mul_vec           ",
        "#13 --> mat_solve_tril_vec    ",
        "#14 --> mat_solve_triu_vec    "
    };

    if (argc == 1)  // show usage
    {
        printf("Usage: `%s [fun] [nbits] [len]`\n", argv[0]);
        printf("   Each argument is optional; no argument shows this help.\n");
        printf("   - fun: id number of the timed function (see below),\n");
        printf("          exception: fun == -1 times all available functions successively\n");
        printf("   - nbits: number of bits for the modulus, chosen as nextprime(2**(nbits-1))\n");
        printf("          exception: nbits == 232 and 263 (moduli 2**32, 2**63)\n");
        printf("   - len: length for the vector, row and column dimension for the matrices\n");
        printf("\nAvailable functions:\n");
        for (slong j = 0; j < nfuns; j++)
            printf("   %s\n", description[j]);

        return 0;
    }

    printf("#warmup... ");
    for (slong i = 0; i < 10; i++)
    {
        time_dot(1000000, (UWORD(1)<<40)+5, state);
        printf(" ");
    }
    printf("\n");

    if (argc == 2 && atoi(argv[1]) == -1)  // launching full suite
    {
        for (slong ifun = 0; ifun < nfuns; ifun++)
        {
            const timefun tfun = funs[ifun];

            printf("\n%s\n", description[ifun]);
            printf("#bits\\len");
            for (slong i = 0; i < nlens; i++)
                printf("%9ld", lens[i]);
            printf("\n");

            for (slong j = 0; j < nbits; j++)
            {
                const slong b = bits[j];

                printf("%-10ld", b);
                ulong n;
                if (b == 232)
                    n = UWORD(1) << 32;
                else if (b == 263)
                    n = UWORD(1) << 63;
                else
                    n = n_nextprime(UWORD(1) << (b-1), 0);
                for (slong i = 0; i < nlens; i++)
                {
                    tfun(lens[i], n, state);
                    printf(" ");
                }
                printf("\n");
            }
        }
    }
    else if (argc == 2)  // function is given
    {
        const timefun tfun = funs[atoi(argv[1])];

        printf("\n%s\n", description[atoi(argv[1])]);
        printf("#bits\\len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");

        for (slong j = 0; j < nbits; j++)
        {
            const slong b = bits[j];

            printf("%-10ld", b);
            ulong n;
            if (b == 232)
                n = UWORD(1) << 32;
            else if (b == 263)
                n = UWORD(1) << 63;
            else
                n = n_nextprime(UWORD(1) << (b-1), 0);
            for (slong i = 0; i < nlens; i++)
            {
                tfun(lens[i], n, state);
                printf(" ");
            }
            printf("\n");
        }
    }
    else if (argc == 3)  // function + nbits given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);

        printf("\n%s\n", description[atoi(argv[1])]);
        printf("#bits\\len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");

        printf("%-10ld", b);
        ulong n;
        if (b == 232)
            n = UWORD(1) << 32;
        else if (b == 263)
            n = UWORD(1) << 63;
        else
            n = n_nextprime(UWORD(1) << (b-1), 0);
        for (slong i = 0; i < nlens; i++)
        {
            tfun(lens[i], n, state);
            printf(" ");
        }
        printf("\n");
    }
    else if (argc == 4)  // function + nbits + len given
    {
        const timefun tfun = funs[atoi(argv[1])];
        const slong b = atoi(argv[2]);
        const slong len = atoi(argv[3]);

        printf("\n%s\n", description[atoi(argv[1])]);
        printf("#bits\\len");
        for (slong i = 0; i < nlens; i++)
            printf("%9ld", lens[i]);
        printf("\n");

        printf("%-10ld", b);
        ulong n;
        if (b == 232)
            n = UWORD(1) << 32;
        else if (b == 263)
            n = UWORD(1) << 63;
        else
            n = n_nextprime(UWORD(1) << (b-1), 0);

        tfun(len, n, state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
