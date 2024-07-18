/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_poly.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arith.h>
#include <flint/profiler.h>

int main(int argc, char *argv[])
{
    fmpz_poly_t f, g;
    fmpz_poly_factor_t fac;
    fmpz_t t;
    acb_ptr roots;
    slong compd, printd, i, j, deg;
    int flags;

    if (argc < 2)
    {
        flint_printf("poly_roots [-refine d] [-print d] <poly>\n\n");

        flint_printf("Isolates all the complex roots of a polynomial with integer coefficients.\n\n");

        flint_printf("If -refine d is passed, the roots are refined to a relative tolerance\n");
        flint_printf("better than 10^(-d). By default, the roots are only computed to sufficient\n");
        flint_printf("accuracy to isolate them. The refinement is not currently done efficiently.\n\n");

        flint_printf("If -print d is passed, the computed roots are printed to d decimals.\n");
        flint_printf("By default, the roots are not printed.\n\n");

        flint_printf("The polynomial can be specified by passing the following as <poly>:\n\n");

        flint_printf("a <n>          Easy polynomial 1 + 2x + ... + (n+1)x^n\n");
        flint_printf("t <n>          Chebyshev polynomial T_n\n");
        flint_printf("u <n>          Chebyshev polynomial U_n\n");
        flint_printf("p <n>          Legendre polynomial P_n\n");
        flint_printf("c <n>          Cyclotomic polynomial Phi_n\n");
        flint_printf("s <n>          Swinnerton-Dyer polynomial S_n\n");
        flint_printf("b <n>          Bernoulli polynomial B_n\n");
        flint_printf("w <n>          Wilkinson polynomial W_n\n");
        flint_printf("e <n>          Taylor series of exp(x) truncated to degree n\n");
        flint_printf("m <n> <m>      The Mignotte-like polynomial x^n + (100x+1)^m, n > m\n");
        flint_printf("coeffs <c0 c1 ... cn>        c0 + c1 x + ... + cn x^n\n\n");

        flint_printf("Concatenate to multiply polynomials, e.g.: p 5 t 6 coeffs 1 2 3\n");
        flint_printf("for P_5(x)*T_6(x)*(1+2x+3x^2)\n\n");

        return 1;
    }

    compd = 0;
    printd = 0;
    flags = ARB_FMPZ_POLY_ROOTS_VERBOSE;

    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_init(t);
    fmpz_poly_one(f);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-refine"))
        {
            compd = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "-print"))
        {
            printd = atol(argv[i+1]);
            i++;
        }
        else if (!strcmp(argv[i], "a"))
        {
            slong n = atol(argv[i+1]);
            fmpz_poly_zero(g);
            for (j = 0; j <= n; j++)
                fmpz_poly_set_coeff_ui(g, j, j+1);
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "t"))
        {
            fmpz_poly_chebyshev_t(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "u"))
        {
            fmpz_poly_chebyshev_u(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "p"))
        {
            fmpq_poly_t h;
            fmpq_poly_init(h);
            fmpq_poly_legendre_p(h, atol(argv[i+1]));
            fmpq_poly_get_numerator(g, h);
            fmpz_poly_mul(f, f, g);
            fmpq_poly_clear(h);
            i++;
        }
        else if (!strcmp(argv[i], "c"))
        {
            fmpz_poly_cyclotomic(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "s"))
        {
            fmpz_poly_swinnerton_dyer(g, atol(argv[i+1]));
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "b"))
        {
            fmpq_poly_t h;
            fmpq_poly_init(h);
            arith_bernoulli_polynomial(h, atol(argv[i+1]));
            fmpq_poly_get_numerator(g, h);
            fmpz_poly_mul(f, f, g);
            fmpq_poly_clear(h);
            i++;
        }
        else if (!strcmp(argv[i], "w"))
        {
            slong n = atol(argv[i+1]);
            fmpz_poly_zero(g);
            fmpz_poly_fit_length(g, n+2);
            arith_stirling_number_1_vec(g->coeffs, n+1, n+2);
            _fmpz_poly_set_length(g, n+2);
            fmpz_poly_shift_right(g, g, 1);
            fmpz_poly_mul(f, f, g);
            i++;
        }
        else if (!strcmp(argv[i], "e"))
        {
            fmpq_poly_t h;
            fmpq_poly_init(h);
            fmpq_poly_set_coeff_si(h, 0, 0);
            fmpq_poly_set_coeff_si(h, 1, 1);
            fmpq_poly_exp_series(h, h, atol(argv[i+1]) + 1);
            fmpq_poly_get_numerator(g, h);
            fmpz_poly_mul(f, f, g);
            fmpq_poly_clear(h);
            i++;
        }
        else if (!strcmp(argv[i], "m"))
        {
            fmpz_poly_zero(g);
            fmpz_poly_set_coeff_ui(g, 0, 1);
            fmpz_poly_set_coeff_ui(g, 1, 100);
            fmpz_poly_pow(g, g,  atol(argv[i+2]));
            fmpz_poly_set_coeff_ui(g, atol(argv[i+1]), 1);
            fmpz_poly_mul(f, f, g);
            i += 2;
        }
        else if (!strcmp(argv[i], "coeffs"))
        {
            fmpz_poly_zero(g);
            i++;
            j = 0;
            while (i < argc)
            {
                if (fmpz_set_str(t, argv[i], 10) != 0)
                {
                    i--;
                    break;
                }

                fmpz_poly_set_coeff_fmpz(g, j, t);
                i++;
                j++;
            }
            fmpz_poly_mul(f, f, g);
        }
    }

    fmpz_poly_factor_init(fac);

    flint_printf("computing squarefree factorization...\n");
    TIMEIT_ONCE_START
    fmpz_poly_factor_squarefree(fac, f);
    TIMEIT_ONCE_STOP

    TIMEIT_ONCE_START
    for (i = 0; i < fac->num; i++)
    {
        deg = fmpz_poly_degree(fac->p + i);

        flint_printf("%wd roots with multiplicity %wd\n", deg, fac->exp[i]);
        roots = _acb_vec_init(deg);

        arb_fmpz_poly_complex_roots(roots, fac->p + i, flags, compd * 3.32193 + 2);

        if (printd)
        {
            for (j = 0; j < deg; j++)
            {
                acb_printn(roots + j, printd, 0);
                flint_printf("\n");
            }
        }

        _acb_vec_clear(roots, deg);
    }
    TIMEIT_ONCE_STOP

    fmpz_poly_factor_clear(fac);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_clear(t);

    flint_cleanup();
    return 0;
}

