/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "arb.h"
#include "acb.h"
#include "arb_poly.h"
#include "flint/profiler.h"

void
keiper_li_series(arb_ptr z, slong len, slong prec)
{
    arb_ptr t, u, v;

    t = _arb_vec_init(len);
    u = _arb_vec_init(len);
    v = _arb_vec_init(len);

    /* -zeta(s) */
    flint_printf("zeta: ");
    TIMEIT_ONCE_START
    arb_zero(t + 0);
    arb_one(t + 1);
    arb_one(u);
    _arb_poly_zeta_series(v, t, 2, u, 0, len, prec);
    _arb_vec_neg(v, v, len);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    /* logarithm */
    flint_printf("log: ");
    TIMEIT_ONCE_START
    _arb_poly_log_series(t, v, len, len, prec);
    TIMEIT_ONCE_STOP

    /* add log(gamma(1+s/2)) */
    flint_printf("gamma: ");
    TIMEIT_ONCE_START
    arb_one(u);
    arb_one(u + 1);
    arb_mul_2exp_si(u + 1, u + 1, -1);
    _arb_poly_lgamma_series(v, u, 2, len, prec);
    _arb_vec_add(t, t, v, len, prec);
    TIMEIT_ONCE_STOP

    /* subtract 0.5 s log(pi) */
    arb_const_pi(u, prec);
    arb_log(u, u, prec);
    arb_mul_2exp_si(u, u, -1);
    arb_sub(t + 1, t + 1, u, prec);

    /* add log(1-s) */
    arb_one(u);
    arb_set_si(u + 1, -1);
    _arb_poly_log_series(v, u, 2, len, prec);
    _arb_vec_add(t, t, v, len, prec);

    /* binomial transform */
    flint_printf("binomial transform: ");
    TIMEIT_ONCE_START
    arb_set(z, t);
    _arb_vec_neg(t + 1, t + 1, len - 1);
    _arb_poly_binomial_transform(z + 1, t + 1, len - 1, len - 1, prec);
    TIMEIT_ONCE_STOP

    _arb_vec_clear(t, len);
    _arb_vec_clear(u, len);
    _arb_vec_clear(v, len);
}

int main(int argc, char *argv[])
{
    slong i, len, prec, num_threads;
    char * out_file;
    arb_ptr z;

    if (argc < 2)
    {
        flint_printf("keiper_li n [-prec prec] [-threads num_threads] [-out out_file]\n");
        return 1;
    }

    len = atol(argv[1]) + 1;
    prec = 1.1 * len + 50;
    num_threads = 1;
    out_file = NULL;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-prec"))
            prec = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-out"))
            out_file = argv[i+1];
    }

    flint_set_num_threads(num_threads);

    z = _arb_vec_init(len);

    keiper_li_series(z, len, prec);

    for (i = 0; i < len; i++)
    {
        if (i <= 10 || len - i <= 10)
        {
            flint_printf("%wd: ", i); arb_printd(z + i, 50); flint_printf("\n");
        }
    }

    SHOW_MEMORY_USAGE

    if (out_file != NULL)
    {
        fmpz_t man, exp;
        arf_t t;

        FILE * fp = fopen(out_file, "w");

        fmpz_init(man);
        fmpz_init(exp);
        arf_init(t);

        for (i = 0; i < len; i++)
        {
            arf_get_fmpz_2exp(man, exp, arb_midref(z + i));

            flint_fprintf(fp, "%wd ", i);
            fmpz_fprint(fp, man);
            flint_fprintf(fp, " ");
            fmpz_fprint(fp, exp);
            flint_fprintf(fp, " ");

            arf_set_mag(t, arb_radref(z + i));
            arf_get_fmpz_2exp(man, exp, t);

            fmpz_fprint(fp, man);
            flint_fprintf(fp, " ");
            fmpz_fprint(fp, exp);
            flint_fprintf(fp, "\n");
        }

        fclose(fp);

        fmpz_clear(man);
        fmpz_clear(exp);
        arf_clear(t);
    }

    _arb_vec_clear(z, len);
    flint_cleanup();
    return 0;
}

