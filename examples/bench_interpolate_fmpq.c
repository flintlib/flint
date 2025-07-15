#include <stdlib.h>
#include <stdio.h>

#include <flint/profiler.h>
#include <flint/fmpz.h>
#include "flint/fmpq.h"
#include <flint/fmpz_vec.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/ulong_extras.h>
#include "flint/gr.h"
#include "flint/gr_poly.h"
#include "flint/gr_vec.h"

// Génère des données aléatoires xs != xs[j], et ys quelconques
void generate_unique_fmpq_vec(fmpq *xs, fmpq *ys, slong n, slong nbits, flint_rand_t state)
{
    slong i = 0;
    fmpq_t tmp;
    fmpq_init(tmp);

    while (i < n)
    {
        fmpq_randtest(tmp, state, nbits);
        int duplicate = 0;
        for (slong j = 0; j < i; j++)
        {
            if (fmpq_equal(tmp, xs + j))
            {
                duplicate = 1;
                break;
            }
        }
        if (!duplicate)
        {
            fmpq_set(xs + i, tmp);
            fmpq_randtest(ys + i, state, nbits);
            i++;
        }
    }

    fmpq_clear(tmp);
}

void generate_unique_fmpz_fmpq_vec(fmpz *xs, fmpq *ys, slong n, slong nbits, flint_rand_t state)
{
    slong i = 0;
    fmpz_t tmp;
    fmpz_init(tmp);

    while (i < n)
    {
        fmpz_randtest(tmp, state, nbits);
        int duplicate = 0;
        for (slong j = 0; j < i; j++)
        {
            if (fmpz_equal(tmp, xs + j))
            {
                duplicate = 1;
                break;
            }
        }
        if (!duplicate)
        {
            fmpz_set(xs + i, tmp);
            fmpq_randtest(ys + i, state, nbits);
            i++;
        }
    }

    fmpz_clear(tmp);
}

void generate_unique_fmpz_vec(fmpz *xs, fmpz *ys, slong n, slong nbits, flint_rand_t state)
{
    slong i = 0;
    fmpz_t tmp;
    fmpz_init(tmp);

    while (i < n)
    {
        fmpz_randtest(tmp, state, nbits);
        int duplicate = 0;
        for (slong j = 0; j < i; j++)
        {
            if (fmpz_equal(tmp, xs + j))
            {
                duplicate = 1;
                break;
            }
        }
        if (!duplicate)
        {
            fmpz_set(xs + i, tmp);
            fmpz_randtest(ys + i, state, nbits);
            i++;
        }
    }

    fmpz_clear(tmp);
}

// génération d'un polynôme P(x) aléatoire entier de degré n-1
void generate_random_polynomial(fmpz_poly_t poly, slong n, flint_rand_t state)
{
    fmpz_poly_zero(poly);
    fmpz_poly_fit_length(poly, n);
    for (slong i = 0; i < n; i++)
        fmpz_randtest(poly->coeffs + i, state, 100);
    _fmpz_poly_set_length(poly, n);
    _fmpz_poly_normalise(poly);
}

// évaluation entière : ys[i] = P(xs[i])
void evaluate_poly_at_points(fmpz *ys, const fmpz_poly_t poly, const fmpz *xs, slong n)
{
    fmpz_t tmp;
    fmpz_init(tmp);
    for (slong i = 0; i < n; i++)
        fmpz_poly_evaluate_fmpz(ys + i, poly, xs + i);
    fmpz_clear(tmp);
}

// génère xs uniques
void generate_unique_xs(fmpz *xs, slong n, flint_rand_t state)
{
    slong i = 0;
    fmpz_t x;
    fmpz_init(x);

    while (i < n)
    {
        fmpz_randtest(x, state, 100);
        int duplicate = 0;
        for (slong j = 0; j < i; j++)
        {
            if (fmpz_equal(x, xs + j))
            {
                duplicate = 1;
                break;
            }
        }
        if (!duplicate)
        {
            fmpz_set(xs + i, x);
            i++;
        }
    }

    fmpz_clear(x);
}

int main(int argc, char* argv[])
{
    slong j;

    flint_rand_t state;
    flint_rand_init(state);

    slong n = 300;  // taille du système
    slong reps = 1; // number of repetition
    slong nbits = 100; // bit size of input points

    fmpq *xs = _fmpq_vec_init(n);
    fmpq *ys = _fmpq_vec_init(n);

    fmpq_poly_t poly1;
    fmpq_poly_t poly2;
    fmpq_poly_init(poly1);
    fmpq_poly_init(poly2);


    // Étapes
    /*fmpz_poly_t P;
    fmpz_poly_init(P);
    generate_random_polynomial(P, n, state);
    generate_unique_xs(xs, n, state);
    evaluate_poly_at_points(ys, P, xs, n);*/

    //_fmpq_vec_print(xs, n); printf("\n");
    //_fmpq_vec_print(ys, n); printf("\n\n");

    /*gr_ctx_t ctx;
    gr_ctx_init_fmpq(ctx);

    gr_poly_t poly3;
    gr_vec_t x3, y3;
    gr_vec_init(x3, n, ctx);
    gr_vec_init(y3, n, ctx);
    gr_poly_init(poly3, ctx);
    for (j = 0; j < n; j++) {
            gr_set_fmpq(gr_vec_entry_ptr(x3, j, ctx), xs + j, ctx);
            gr_set_fmpq(gr_vec_entry_ptr(y3, j, ctx), ys + j, ctx);
    }
    */

    printf("Benchmark: interpolation rationnelle (%ld points; %ld repetitions;  %ld bitsize)\n",
        n, reps, nbits);

    slong tnaive, tfast;
    timeit_t t0;

    tnaive = 0;
    tfast = 0;

    for (j = 0; j < reps; j++) {
        generate_unique_fmpq_vec(xs, ys, n, nbits, state);

        // Classic
        timeit_start(t0);
        //fmpq_poly_interpolate_fmpq_vec(poly1, xs, ys, n);
        timeit_stop(t0);
        //printf("Classic fmpq : cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);
        tnaive += t0->cpu;
        /*
        fmpq *xs1, *ys1;
        xs1 = _fmpq_vec_init(n);
        ys1 = _fmpq_vec_init(n);
        for (j = 0; j < n; j++) {
            fmpq_set_fmpz(xs1 + j, xs + j);
            fmpq_set_fmpz(ys1 + j, ys + j);
        }*/
        // Classic bis
        timeit_start(t0);
        fmpq_poly_interpolate_fast(poly2, xs, ys, n);
        timeit_stop(t0);
        //printf("Fast fmpq: cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);
        tfast += t0->cpu;

        int match = fmpq_poly_equal(poly1, poly2);
        if (!match)
            printf("match ? : %s\n", match ? "OUI" : "NON");
        //match = fmpz_equal(poly1->den, poly2->den);
        //printf("match den ? : %s\n", match ? "OUI" : "NON");
    }

    printf("Classic fmpq : %ld ms   Fast: %ld ms\n", tnaive, tfast);

    // Gr based
    /*timeit_start(t0);
    ret = gr_poly_interpolate_fast(poly3, x3, y3, ctx);
    timeit_stop(t0);
    printf("Gr  : cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);

    // Vérification (optionnelle)
    gr_poly_t test;
    gr_poly_init(test, ctx);
    gr_poly_set_fmpq_poly(test, poly1, ctx);
    truth_t matchgr = gr_poly_equal(poly3, test, ctx);
    printf("GR worked ? : %s\n", ret == GR_SUCCESS ? "OUI" : "NON");
    printf("match ? : %s\n", matchgr == T_TRUE ? "OUI" : "NON");*/


    // Clean
    //_fmpq_vec_clear(xs, n);
    //_fmpq_vec_clear(ys, n);
    fmpq_poly_clear(poly1);
    fmpq_poly_clear(poly2);
    flint_rand_clear(state);

    return 0;
}
