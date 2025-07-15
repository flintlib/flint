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

int main()
{
    slong d = 500;     // degré du polynôme
    slong n = 10000;    // nombre de points d’évaluation

    FLINT_TEST_INIT(state);
    flint_rand_t rstate;
    flint_rand_init(rstate);

    fmpq_poly_t P;
    fmpq_poly_init(P);
    fmpq_poly_randtest(P, state, d + 1, 64);  // gros coefficients (512 bits)

    fmpq *points = _fmpq_vec_init(n);
    fmpq *values_vec = _fmpq_vec_init(n);
    fmpq *values_naive = _fmpq_vec_init(n);
    fmpq_t tmp;

    fmpq_init(tmp);

    // Génère n points aléatoires x_i ∈ Q
    for (slong i = 0; i < n; i++)
    {
        fmpz_t a, b;
        fmpz_init(a);
        fmpz_init(b);
        fmpz_randtest_not_zero(a, state, 64);
        fmpz_randtest_not_zero(b, state, 64);
        fmpq_set_fmpz_frac(points + i, a, b);
        fmpz_clear(a);
        fmpz_clear(b);
    }

    printf("⏱ Évaluation naïve (itérative)...\n");
    timeit_t t2;
    timeit_start(t2);
    for (slong i = 0; i < n; i++)
        fmpq_poly_evaluate_fmpq(values_naive + i, P, points + i);
    timeit_stop(t2);
    flint_printf("Temps (naïf): %ld ms\n", t2->cpu);

    printf("⏱ Évaluation rapide (vectorisée)...\n");
    timeit_t t1;
    timeit_start(t1);
    fmpq_poly_evaluate_vec_fast(values_vec, P, points, n);
    timeit_stop(t1);
    flint_printf("Temps (fast): %ld ms\n", t1->cpu);


    // Vérification
    slong errors = 0;
    for (slong i = 0; i < n; i++)
    {
        if (!fmpq_equal(values_vec + i, values_naive + i))
        {
            errors++;
            if (errors < 5)
            {
                flint_printf("Erreur à l’indice %wd\n", i);
                //flint_printf("  fast : "); fmpq_print(values_vec + i); flint_printf("\n");
                //flint_printf("  naïf : "); fmpq_print(values_naive + i); flint_printf("\n");
            }
        }
    }

    if (errors == 0)
        flint_printf("✅ Toutes les valeurs sont correctes.\n");
    else
        flint_printf("❌ %wd erreurs détectées.\n", errors);

    // Nettoyage
    fmpq_poly_clear(P);
    _fmpq_vec_clear(points, n);
    _fmpq_vec_clear(values_vec, n);
    _fmpq_vec_clear(values_naive, n);
    fmpq_clear(tmp);
    flint_rand_clear(rstate);
    FLINT_TEST_CLEAR(state);

    return 0;
}
