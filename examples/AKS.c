/*
   Naive AKS primality test [1,2].

   This is only an educational program to illustrate the AKS test
   and FLINT polynomial arithmetic. AKS is not a practical primality test
   (fmpz_is_prime is millions of times faster), and this program
   has not been optimized for speed. An exercise for the reader is to
   improve this program.

   [1] Agrawal, Kayal and Saxena, "PRIMES is in P",
       Annals of Mathematics. 160 (2): 781–793
   [2] https://en.wikipedia.org/wiki/AKS_primality_test

   This file is public domain. Author: Fredrik Johansson.
*/

#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "gr.h"
#include "profiler.h"

/* Check (X+a)^n = X^n+a  (mod X^r-1, n). */
int AKS_main(const fmpz_t n, ulong a, ulong r)
{
    fmpz_mod_poly_t Xr1, Xan, Xna, t;
    fmpz_mod_ctx_t modn_ctx;
    int result;

    fmpz_mod_ctx_init(modn_ctx, n);

    fmpz_mod_poly_init(Xr1, modn_ctx);
    fmpz_mod_poly_init(Xan, modn_ctx);
    fmpz_mod_poly_init(Xna, modn_ctx);
    fmpz_mod_poly_init(t, modn_ctx);

    /* X^r-1 */
    fmpz_mod_poly_set_coeff_si(Xr1, 0, -1, modn_ctx);
    fmpz_mod_poly_set_coeff_ui(Xr1, r, 1, modn_ctx);

    /* (X+a)^n */
    fmpz_mod_poly_set_coeff_ui(Xan, 0, a, modn_ctx);
    fmpz_mod_poly_set_coeff_ui(Xan, 1, 1, modn_ctx);
    fmpz_mod_poly_powmod_fmpz_binexp(Xan, Xan, n, Xr1, modn_ctx);

    /* X^n+a */
    fmpz_mod_poly_set_coeff_ui(Xna, 1, 1, modn_ctx);
    fmpz_mod_poly_powmod_fmpz_binexp(Xna, Xna, n, Xr1, modn_ctx);
    fmpz_mod_poly_set_ui(t, a, modn_ctx);
    fmpz_mod_poly_add(Xna, Xna, t, modn_ctx);

    result = fmpz_mod_poly_equal(Xan, Xna, modn_ctx);

    fmpz_mod_poly_clear(Xr1, modn_ctx);
    fmpz_mod_poly_clear(Xan, modn_ctx);
    fmpz_mod_poly_clear(Xna, modn_ctx);
    fmpz_mod_poly_clear(t, modn_ctx);

    fmpz_mod_ctx_clear(modn_ctx);

    return result;
}

/* Multiplicative order of a mod m. There is currently no
   n_multiplicative_order in FLINT, so here is a completely naive
   implementation. */
ulong multiplicative_order(ulong a, ulong m)
{
    ulong e;

    for (e = 1; ; e++)
    {
        if (n_powmod2(a, e, m) == 1)
            return e;
    }
}

/* Wrapper because we don't want the root computed by fmpz_is_perfect_power. */
int
is_perfect_power(const fmpz_t n)
{
    fmpz_t t;
    int result;
    fmpz_init(t);
    result = fmpz_is_perfect_power(t, n);
    fmpz_clear(t);
    return result;
}

int fmpz_is_prime_AKS(const fmpz_t n, int verbose)
{
    ulong a, r, n_mod_r, ord, log2n, bound;

    if (fmpz_cmp_ui(n, 1) <= 0)
        return 0;

    if (is_perfect_power(n))
        return 0;

    if (verbose)
        flint_printf("n = %{fmpz}\n", n);

    /* Upper bound for log2(n) */
    log2n = fmpz_clog_ui(n, 2);

    for (r = 2; ; r++)
    {
        if (verbose && (r & (r - 1)) == 0)
            flint_printf("  Testing candidate r = %wu\n", r);

        n_mod_r = fmpz_fdiv_ui(n, r);
        if (n_gcd(n_mod_r, r) != 1)
            continue;

        ord = multiplicative_order(n_mod_r, r);
        if (ord > log2n * log2n)
            break;
    }

    if (verbose)
        flint_printf("  Found r = %wu\n", r);

    for (a = 2; a <= r && fmpz_cmp_ui(n, a) > 0; a++)
    {
        if (verbose && (a & (a - 1)) == 0)
            flint_printf("  Trial division %wu / %wu\n", a, r);

        if (fmpz_fdiv_ui(n, a) == 0)
            return 0;
    }

    if (fmpz_cmp_ui(n, r) <= 0)
        return 1;

    bound = (n_sqrt(n_euler_phi(r)) + 1) * log2n;

    for (a = 1; a <= bound; a++)
    {
        if (verbose && (a & (a - 1)) == 0)
            flint_printf("  AKS test %wu / %wu\n", a, bound);

        if (!AKS_main(n, a, r))
            return 0;
    }

    return 1;
}

int main(int argc, char * argv[])
{
    fmpz_t n;
    fmpz_init(n);

    if (argc == 2)
    {
        /* allow expression input like "2^64+1" */
        gr_ctx_t ZZ;
        gr_ctx_init_fmpz(ZZ);

        if (gr_set_str(n, argv[1], ZZ) != GR_SUCCESS)
        {
            flint_printf("unable to parse integer\n");
            return 1;
        }

        gr_ctx_clear(ZZ);

        TIMEIT_ONCE_START
        flint_printf("%{fmpz} : %s\n", n, fmpz_is_prime_AKS(n, 1) ? "PRIME" : "COMPOSITE");
        TIMEIT_ONCE_STOP

        flint_printf("\nfmpz_is_prime for comparison:\n");
        int is_prime;
        TIMEIT_START
        is_prime = fmpz_is_prime(n);
        TIMEIT_STOP
        flint_printf("%{fmpz} : %s\n", n, is_prime ? "PRIME" : "COMPOSITE");
    }
    else
    {
        flint_printf("No input given, running examples:\n");

        flint_printf("\nVerifying the algorithm up to n = 1000...\n");
        TIMEIT_ONCE_START
        for (fmpz_one(n); fmpz_cmp_ui(n, 1000) <= 0; fmpz_add_ui(n, n, 1))
        {
            if (fmpz_fdiv_ui(n, 100) == 0)
                flint_printf("%{fmpz}...\n", n);

            if (fmpz_is_prime_AKS(n, 0) != fmpz_is_prime(n))
                flint_printf("FAIL: wrong result for n = %{fmpz}\n", n);
        }
        flint_printf("OK\n");
        TIMEIT_ONCE_STOP

        flint_printf("\nTesting a 10-digit prime\n");
        fmpz_set_str(n, "6772298791", 10);
        TIMEIT_ONCE_START
        flint_printf("%{fmpz} : %s\n", n, fmpz_is_prime_AKS(n, 1) ? "PRIME" : "COMPOSITE");
        TIMEIT_ONCE_STOP

        flint_printf("\nTesting a 50-digit semiprime\n");
        fmpz_set_str(n, "6792336996936278623737332130584104553341324999517", 10);
        TIMEIT_ONCE_START
        flint_printf("%{fmpz} : %s\n", n, fmpz_is_prime_AKS(n, 1) ? "PRIME" : "COMPOSITE");
        TIMEIT_ONCE_STOP
    }

    fmpz_clear(n);

    flint_cleanup_master();
    return 0;
}
