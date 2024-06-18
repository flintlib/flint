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

void time_dot(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    _nmod_vec_dot(v1, v2, len, mod, n_limbs);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

void time_dot_rev(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

    nn_ptr v1 = _nmod_vec_init(len);
    _nmod_vec_rand(v1, state, len, mod);
    nn_ptr v2 = _nmod_vec_init(len);
    _nmod_vec_rand(v2, state, len, mod);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    _nmod_vec_dot_rev(v1, v2, len, mod, n_limbs);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

void time_dot_expr(ulong len, ulong n, flint_rand_t state)
{
    nmod_t mod;
    nmod_init(&mod, n);
    const int n_limbs = _nmod_vec_dot_bound_limbs(len, mod);

    nn_ptr v1 = _nmod_vec_init(9*len);
    _nmod_vec_rand(v1, state, 9*len, mod);
    nn_ptr v2 = _nmod_vec_init(9*len);
    _nmod_vec_rand(v2, state, 9*len, mod);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    nn_srcptr v1i = v1;
    nn_srcptr v2i = v2;
    ulong i, FLINT_SET_BUT_UNUSED(res);

    TIMEIT_START
    NMOD_VEC_DOT(res, i, len, v1i[9*len - 1 - 9*i], v2i[9*len - 1 - 9*i], mod, n_limbs);
    TIMEIT_STOP_VALUES(tcpu, twall)

    printf("%.2e", twall);

    _nmod_vec_clear(v1);
    _nmod_vec_clear(v2);
}

/*-------------------------*/
/* indirect: poly          */
/*-------------------------*/

// void _nmod_poly_inv_series_basecase_preinv1(nn_ptr Qinv, nn_srcptr Q, slong Qlen, slong n, ulong q, nmod_t mod)
// void _nmod_poly_exp_series(nn_ptr f, nn_srcptr h, slong hlen, slong n, nmod_t mod)

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
    if (len > 10000 || n == (UWORD(1) << 63))
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
    if (len > 10000 || n == (UWORD(1) << 63))
    {
        printf("        ");
        return;
    }

    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, n);
    int status;

    gr_poly_t p;
    gr_poly_init(p, ctx);
    gr_poly_randtest(p, state, len, ctx);
    status = gr_poly_set_coeff_si(p, 0, WORD(0), ctx);
    gr_poly_t res;
    gr_poly_init(res, ctx);

    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START
    status |= gr_poly_exp_series_basecase(res, p, len, ctx);

// int gr_poly_exp_series_basecase(gr_poly_t f, const gr_poly_t h, slong n, gr_ctx_t ctx)

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
    if (len > 4000 || n == (UWORD(1) << 63))
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
    if (len > 10000 || n == (UWORD(1) << 63))
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
    if (len > 4000 || n == (UWORD(1) << 63))
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
    if (len > 10000 || n == (UWORD(1)<<63))
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
    const slong nbits = 12;
    const ulong bits[] = {0, 12, 28, 30, 31, 32, 40, 50, 60, 61, 62, 63, 64};

    // vector lengths
    const slong nlens = 14;
    const ulong lens[] = {1, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 100000, 1000000};

    // bench functions
    const slong nfuns = 12;
    typedef void (*timefun) (ulong, ulong, flint_rand_t);
    const timefun funs[] = {
        time_dot,                    // 0
        time_dot_rev,                // 1
        time_dot_expr,               // 2
        time_dot_poly_mul,           // 3
        time_dot_poly_inv_series,    // 4
        time_dot_poly_exp_series,    // 5
        time_dot_mat_mul,            // 6
        time_dot_mat_solve_tril,     // 7
        time_dot_mat_solve_triu,     // 8
        time_dot_mat_mul_vec,        // 9
        time_dot_mat_solve_tril_vec, // 10
        time_dot_mat_solve_triu_vec, // 11
    };

    const char * description[] = {
        "#0  --> vec dot           ",
        "#1  --> vec dot rev       ",
        "#2  --> vec dot expr      ",
        "#3  --> poly_mul          ",
        "#4  --> poly_inv_series   ",
        "#5  --> poly_exp_series   ",
        "#6  --> mat_mul           ",
        "#7  --> mat_solve_tril    ",
        "#8  --> mat_solve_triu    ",
        "#9  --> mat_mul_vec       ",
        "#10 --> mat_solve_tril_vec",
        "#11 --> mat_solve_triu_vec"
    };


    printf("#warmup... ");
    for (slong i = 0; i < 10; i++)
    {
        time_dot(1000000, (UWORD(1)<<40)+5, state);
        printf(" ");
    }
    printf("\n");

    if (argc == 1)  // launching full suite
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
                const ulong n = (b==0) ? (UWORD(1) << 63) : n_nextprime(UWORD(1) << (b-1), 0);
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
            const ulong n = (b==0) ? (UWORD(1) << 63) : n_nextprime(UWORD(1) << (b-1), 0);
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
        const ulong n = (b==0) ? (UWORD(1) << 63) : n_nextprime(UWORD(1) << (b-1), 0);
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
        const ulong n = (b==0) ? (UWORD(1) << 63) : n_nextprime(UWORD(1) << (b-1), 0);

        tfun(len, n, state);
        printf("\n");
    }

    flint_rand_clear(state);
    return 0;
}
