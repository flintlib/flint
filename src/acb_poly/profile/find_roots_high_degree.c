#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "acb_poly.h"

double normal_rand()
{
    double u1, u2, z;
    double pi=acos(-1);
    do {
        u1 = (double)rand() / RAND_MAX;
    } while (u1 == 0.0);
    u2 = (double) rand() / RAND_MAX;
    z = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);
    return z;
}

void set_coeffs_easy(acb_ptr coeffs, slong n)
{
    for(slong i=0; i<n; i++) {
        acb_set_si(coeffs + i, i+1);
    }
}

void set_coeffs_flat(acb_ptr coeffs, slong n, slong prec)
{
    slong i;
    arb_t inv_sqrt_fac;
    arb_init(inv_sqrt_fac);

    srand(42);
    for(i=0; i<n; i++) {
        acb_set_d_d(coeffs+i, normal_rand(), normal_rand());
        arb_fac_ui(inv_sqrt_fac, n, prec);
        arb_inv(inv_sqrt_fac, inv_sqrt_fac, prec);
        arb_sqrt(inv_sqrt_fac, inv_sqrt_fac, prec);
        acb_mul_arb(coeffs+i, coeffs+i, inv_sqrt_fac, prec);
    }
    arb_clear(inv_sqrt_fac);
}

void set_coeffs_hyperbolic(acb_ptr coeffs, slong n, slong prec)
{
    slong i;
    srand(42);
    for(i=0; i<n; i++) {
        acb_set_d_d(coeffs+i, normal_rand(), normal_rand());
    }
}

void set_coeffs_elliptic(acb_ptr coeffs, slong n, slong prec)
{
    slong i;
    arb_t sqrt_binom;
    arb_init(sqrt_binom);

    srand(42);
    for(i=0; i<n; i++) {
        acb_set_d_d(coeffs+i, normal_rand(), normal_rand());
        arb_bin_uiui(sqrt_binom, n-1, i, prec);
        arb_sqrt(sqrt_binom, sqrt_binom, prec);
        acb_mul_arb(coeffs+i, coeffs+i, sqrt_binom, prec);
    }
    arb_clear(sqrt_binom);
}

ulong timed_find_roots(acb_ptr roots, acb_srcptr coeffs, slong n, slong prec, char* name)
{
    ulong elapsed;
    struct timespec start, stop;
    flint_printf("%-10s (scaled) ... ", name);
    fflush(stdout);
    clock_gettime(CLOCK_MONOTONIC, &start);
    _acb_poly_find_roots_double(roots, coeffs, NULL, n, 500, prec);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    elapsed = (stop.tv_sec - start.tv_sec)*1000 + (stop.tv_nsec - start.tv_nsec) / 1000000;
    flint_printf("%lu ms | ", elapsed);
    return elapsed;
}

ulong timed_find_roots_double(acb_ptr roots, acb_srcptr coeffs, slong n, slong prec, char* name)
{
    ulong elapsed;
    struct timespec start, stop;
    double *cdz, *cdp;
    arb_t scale;
    slong i;

    arb_init(scale);
    cdp = flint_malloc(2*n*sizeof(double));
    cdz = flint_malloc(2*(n - 1)*sizeof(double));

    for(i=0; i<n; i++) {
        cdp[2*i]   = arf_get_d(arb_midref(acb_realref(coeffs + i)), ARF_RND_NEAR);
        cdp[2*i+1] = arf_get_d(arb_midref(acb_imagref(coeffs + i)), ARF_RND_NEAR);
    }

    flint_printf("%-10s (direct) ... ", name);
    fflush(stdout);
    clock_gettime(CLOCK_MONOTONIC, &start);
    cd_poly_find_roots(cdz, cdp, NULL, n, 500, 1e-53);
    clock_gettime(CLOCK_MONOTONIC, &stop);
    elapsed = (stop.tv_sec - start.tv_sec)*1000 + (stop.tv_nsec - start.tv_nsec) / 1000000;
    flint_printf("%lu ms | ", elapsed);
    for(i=0; i<n-1; i++) {
        acb_set_d_d(roots + i, cdz[2*i], cdz[2*i+1]);
    }

    flint_free(cdp);
    flint_free(cdz);
    return elapsed;
}


slong rel_accuracy(acb_srcptr roots, acb_srcptr coeffs, slong n, slong prec)
{
    slong i, ra=LONG_MAX;
    acb_ptr dcoeffs;
    acb_t r;
    flint_printf("  minimal accuracy bits ... ");
    fflush(stdout);
    dcoeffs = _acb_vec_init(n-1);
    acb_init(r);
    for(i=0; i<n-1; i++) {
        acb_mul_ui(dcoeffs+i, coeffs+i+1, i+1, 2*prec);
    }
    for(i=0; i<n-1; i++) {
        _acb_poly_root_inclusion(r, roots+i, coeffs, dcoeffs, n, 2*prec);
        ra = FLINT_MIN(ra, acb_rel_accuracy_bits(r));
    }
    flint_printf("%ld\n", ra);

    _acb_vec_clear(dcoeffs, n-1);
    acb_clear(r);
    return ra;
}

int main(int argc, char *argv[])
{
    int n;
    acb_ptr coeffs, roots;
    slong prec = 64;

    if(argc<2) {
        flint_printf("find_roots_high_degree <n>\n\n"); //n=2000 should lead to 38 or more accuracy bits 
        flint_printf("Find the roots of polynomials of degree n. For random coefficients, we use\n");
        flint_printf("random variables N_k with independent complex normal distribution.\n");
        flint_printf("The polynomials tested are:\n");
        flint_printf("- Easy polynomials 1 + 2x + ... + (d+1)x^d\n");
        flint_printf("- Hyperbolic polynomials: N_0 + ... + N_n x^d\n");
        flint_printf("- Flat polynomials:       N_0 + ... + 1/sqrt(d!) N_d x^d\n");
        flint_printf("- Elliptic polynomials:   N_0 + ... + sqrt(binom(d,k)) x^k + ... + N_d x^d\n\n");
        return 1;
    }

    n = atoi(argv[1]);
    coeffs = _acb_vec_init(n);
    roots = _acb_vec_init(n-1);

    set_coeffs_easy(coeffs, n);
    timed_find_roots(roots, coeffs, n, prec, "Easy");
    rel_accuracy(roots, coeffs, n, prec);
    timed_find_roots_double(roots, coeffs, n, prec, "Easy");
    rel_accuracy(roots, coeffs, n, prec);

    set_coeffs_hyperbolic(coeffs, n, prec);
    timed_find_roots(roots, coeffs, n, prec, "Hyperbolic");
    rel_accuracy(roots, coeffs, n, prec);
    timed_find_roots_double(roots, coeffs, n, prec, "Hyperbolic");
    rel_accuracy(roots, coeffs, n, prec);

    set_coeffs_flat(coeffs, n, prec);
    timed_find_roots(roots, coeffs, n, prec, "Flat");
    rel_accuracy(roots, coeffs, n, prec);
    timed_find_roots_double(roots, coeffs, n, prec, "Flat");
    rel_accuracy(roots, coeffs, n, prec);

    set_coeffs_elliptic(coeffs, n, prec);
    timed_find_roots(roots, coeffs, n, prec, "Elliptic");
    rel_accuracy(roots, coeffs, n, prec);
    timed_find_roots_double(roots, coeffs, n, prec, "Elliptic");
    rel_accuracy(roots, coeffs, n, prec);

    _acb_vec_clear(coeffs, n);
    _acb_vec_clear(roots, n-1);
    
    return 0;
}



