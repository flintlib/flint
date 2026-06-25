/* This file is public domain. Author: Hartmut Monien. */

#include <stdlib.h>

#include "flint.h"

#include "ulong_extras.h"

#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"

#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

#include "fq.h"
#include "fq_vec.h"
#include "fq_poly.h"
#include "fq_poly_factor.h"

#include "padic.h"

typedef struct
{
    fq_struct *x0;
    slong *multiplicity;
    slong num;
    slong alloc;
} fmpz_poly_roots_fq_struct;

typedef fmpz_poly_roots_fq_struct fmpz_poly_roots_fq_t[1];

void fmpz_poly_roots_fq_init(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
void fmpz_poly_roots_fq_clear(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
int fmpz_poly_roots_fq_print(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
int fmpz_poly_roots_fq_fprint_pretty(FILE * file, fmpz_poly_roots_fq_t roots,
                                     fq_ctx_t fctx);
int fmpz_poly_roots_fq_print_pretty(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx);
void fmpz_poly_roots_fq(fmpz_poly_roots_fq_t roots, fmpz_poly_t poly,
                        fq_ctx_t fctx);


void
fmpz_poly_roots_fq_init(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
    roots->x0 = NULL;
    roots->multiplicity = NULL;
    roots->num = 0;
    roots->alloc = 0;
}

void
fmpz_poly_roots_fq_clear(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
    _fq_vec_clear(roots->x0, roots->alloc, fctx);
    flint_free(roots->multiplicity);
}


int
fmpz_poly_roots_fq_print(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
    _fq_vec_print(roots->x0, roots->num, fctx);
    return 1;
}

int
fmpz_poly_roots_fq_fprint_pretty(FILE *file, fmpz_poly_roots_fq_t roots,
                                 fq_ctx_t fctx)
{
    slong j;

    for (j = 0; j < roots->num; j++)
    {
        fq_fprint_pretty(file, roots->x0 + j, fctx);
        flint_fprintf(file, " %wd\n", roots->multiplicity[j]);
    }

    return 1;
}

int
fmpz_poly_roots_fq_print_pretty(fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
{
    fmpz_poly_roots_fq_fprint_pretty(stdout, roots, fctx);
    return 1;
}

void
fmpz_poly_roots_fq(fmpz_poly_roots_fq_t roots, fmpz_poly_t poly, fq_ctx_t fctx)
{
    slong j, k, num;

    fmpz_mod_ctx_t fmctx;
    fmpz_mod_poly_t mpoly;
    fq_poly_factor_t f;
    fq_poly_t fpoly;

    fmpz_mod_ctx_init(fmctx, fq_ctx_prime(fctx));
    fmpz_mod_poly_init(mpoly, fmctx);
    fmpz_mod_poly_set_fmpz_poly(mpoly, poly, fmctx);

    fq_poly_factor_init(f, fctx);
    fq_poly_init(fpoly, fctx);
    fq_poly_set_fmpz_mod_poly(fpoly, mpoly, fctx);

    fq_poly_roots(f, fpoly, 1, fctx);

    fmpz_mod_poly_clear(mpoly, fmctx);
    fmpz_mod_ctx_clear(fmctx);
    fq_poly_clear(fpoly, fctx);

    num = 0;

    for (j = 0; j < f->num; j++)
    {
        if (fq_poly_degree(f->poly + j, fctx) == 1)
            num++;
    }

    roots->x0 = _fq_vec_init(num, fctx);
    roots->multiplicity = flint_malloc(sizeof(slong) * num);
    roots->num = num;
    roots->alloc = num;

    k = 0;

    for (j = 0; j < f->num; j++)
    {
        if (fq_poly_degree(f->poly + j, fctx) == 1)
        {
            fq_poly_get_coeff(roots->x0 + k, f->poly + j, 0, fctx);
            fq_neg(roots->x0 + k, roots->x0 + k, fctx);
            roots->multiplicity[k] = f->exp[j];
            k++;
        }
    }
}

static void
padic_hensel_iteration(fmpz_poly_t poly,
                       padic_t x, padic_ctx_t ctx, slong prec)
{
    padic_t tmp, y0, y1;
    padic_init2(tmp, prec);
    padic_init2(y0, prec);
    padic_init2(y1, prec);
    do
    {
        /* Horner evaluation of poly and poly' at x */
        padic_set_fmpz(y0, poly->coeffs + poly->length - 1, ctx);
        padic_zero(y1);
        for (slong j = poly->length - 2; j >= 0; j--)
        {
            padic_mul(y1, y1, x, ctx);
            padic_add(y1, y1, y0, ctx);
            padic_mul(y0, y0, x, ctx);
            padic_set_fmpz(tmp, poly->coeffs + j, ctx);
            padic_add(y0, y0, tmp, ctx);
        }
        /* Newton step: x -> x - poly / poly' */
        padic_inv(y1, y1, ctx);
        padic_mul(y1, y1, y0, ctx);
        padic_sub(x, x, y1, ctx);
    }
    while (padic_val(y0));
    padic_clear(tmp);
    padic_clear(y0);
    padic_clear(y1);
}

static padic_struct *
_padic_vec_init2(slong len, slong prec)
{
    slong j;
    padic_struct *vec = flint_malloc(len * sizeof(padic_t));
    for (j = 0; j < len; j++)
    {
        padic_init2(vec + j, prec);
    }
    return vec;
}

static void
_padic_vec_clear(padic_struct *vec, slong len)
{
    slong j;
    for (j = 0; j < len; j++)
    {
        padic_clear(vec + j);
    }
    flint_free(vec);
}

static fmpz_poly_struct *
_fmpz_poly_vec_init(slong len)
{
    slong j;
    fmpz_poly_struct *vec = flint_malloc(len * sizeof(fmpz_poly_t));
    for (j = 0; j < len; j++)
    {
        fmpz_poly_init(vec + j);
    }
    return vec;
}

static void
_fmpz_poly_vec_clear(fmpz_poly_struct *vec, slong len)
{
    slong j;
    for (j = 0; j < len; j++)
    {
        fmpz_poly_clear(vec + j);
    }
    flint_free(vec);
}

void
_padic_roots(fmpz_poly_t poly, fq_ctx_t fctx, padic_ctx_t pctx, slong prec)
{
    slong j, level, js, ns = 1, nr = 0, nz = 0, n = fmpz_poly_degree(poly);
    fmpz_poly_struct *s = _fmpz_poly_vec_init(n);
    padic_struct
        * x0 = _padic_vec_init2(n, prec), *xs = _padic_vec_init2(n, prec);
    fmpz_poly_roots_fq_t froots;
    fmpz_poly_t xi;
    fmpz_poly_init(xi);
    fmpz_poly_set_coeff_fmpz(xi, 1, pctx->p);
    fmpz_t zero;
    fmpz_init(zero);
    fmpz_poly_set(s, poly);
    padic_zero(xs);
    for (level = 0; level < PADIC_DEFAULT_PREC; level++)
    {
        for (js = 0; js < ns; js++)
        {
            if (!fmpz_poly_is_zero(s + js))
            {
                fmpz_poly_roots_fq_init(froots, fctx);
                fmpz_poly_roots_fq(froots, s + js, fctx);
                for (j = 0; j < froots->num; j++)
                {
                    fmpz_poly_get_coeff_fmpz(zero, froots->x0 + j, 0);
                    if (*(froots->multiplicity + j) > 1)
                    {
                        padic_set_fmpz(xs + n - 1 - nr, zero, pctx);
                        padic_shift(xs + n - 1 - nr, xs + n - 1 - nr, level,
                                    pctx);
                        padic_add(xs + n - 1 - nr, xs + n - 1 - nr, xs + js,
                                  pctx);
                        if (level + 1 < PADIC_DEFAULT_PREC)
                        {
                            fmpz_poly_set_coeff_fmpz(xi, 0, zero);
                            fmpz_poly_compose(s + n - 1 - nr, s + js, xi);
                            fmpz_pow_ui(zero, pctx->p,
                                        *(froots->multiplicity + j));
                            if (fmpz_divisible
                                (fmpz_poly_get_coeff_ptr(s + n - 1 - nr, 0),
                                 zero))
                            {
                                fmpz_poly_scalar_divexact_fmpz(s + n - 1 - nr,
                                                               s + n - 1 - nr,
                                                               zero);
                                nr++;
                            }
                        }
                        else
                        {
                            padic_print(xs + n - 1 - nr, pctx);
                            flint_printf(" (%wd)\n",
                                         *(froots->multiplicity + j));
                        }
                    }
                    else
                    {
                        padic_set_fmpz(x0 + nz, zero, pctx);
                        if (!fmpz_is_zero(zero))
                        {
                            padic_shift(x0 + nz, x0 + nz, level, pctx);
                        }
                        padic_add(x0 + nz, x0 + nz, xs + js, pctx);
                        padic_hensel_iteration(poly, x0 + nz, pctx, prec);
                        padic_print(x0 + nz, pctx);
                        flint_printf(" (1)\n");
                        nz++;
                    }
                }
                fmpz_poly_roots_fq_clear(froots, fctx);
            }
            else
            {
                padic_print(xs + n - 1 - nr, pctx);
                flint_printf(" (%wd)\n", fmpz_poly_degree(s + js));
            }
        }
        for (j = 0; j < nr; j++)
        {
            fmpz_poly_set(s + j, s + n - 1 - j);
            padic_set(xs + j, xs + n - 1 - j, pctx);
        }
        ns = nr;
        nr = 0;
    }
    fmpz_clear(zero);
    fmpz_poly_clear(xi);
    _padic_vec_clear(x0, n);
    _padic_vec_clear(xs, n);
    _fmpz_poly_vec_clear(s, n);
}

void
padic_roots(fmpz_poly_t poly, padic_ctx_t pctx, slong prec)
{
    slong j, deg;
    fmpz_t tmp;
    fmpz_poly_factor_t f;
    padic_t x;
    fq_ctx_t fctx;
    fq_ctx_init(fctx, pctx->p, 1, "a");
    fmpz_init(tmp);
    padic_init(x);
    fmpz_poly_factor_init(f);
    fmpz_poly_factor(f, poly);
    for (j = 0; j < f->num; j++)
    {
        deg = fmpz_poly_degree(f->p + j);
        if (deg == 1)
        {
            fmpz_set(tmp, fmpz_poly_get_coeff_ptr(f->p + j, 0));
            fmpz_neg(tmp, tmp);
            padic_set_fmpz(x, tmp, pctx);
            padic_print(x, pctx);
            flint_printf(" (%wd)\n", *(f->exp + j));
        }
        else
        {
            _padic_roots(f->p + j, fctx, pctx, prec);
        }
    }
    fq_ctx_clear(fctx);
}

/* example polynomials */

char *polys[] = {
    "3  -2774119 -2468 1",
    "3  6 -7 1",
    "4  -156 188 -33 1",
    "3  -11 0 1",
    "3  -30 1 1",
    "4  -17576 2028 -78 1",
    "10  -362880 1026576 -1172700 723680 -269325 63273 -9450 870 -45 1",
    "12  -39916800 120543840 -150917976 105258076 -45995730 13339535 -2637558 357423 -32670 1925 -66 1",
    "9  44100 -103740 103429 -57034 19019 -3928 491 -34 1",
    "5  83521 -19652 1734 -68 1",
    "12  22370117 15978655 10666271 5010005 1846306 575366 142702 28538 4585 523 35 1",
    "4  0 -11 0 1",
    "7  3 1 0 0 1 0 1",
    "11  23 -74 89 -68 35 0 -14 8 -2 -1 1",
    "4  -12 0 0 1",
};

int
main(int argc, char *argv[])
{

    const slong np = sizeof(polys) / sizeof(polys[0]);
    flint_printf("# examples: %d\n", np);

    ulong n = 7;

    if (argc != 2)
    {
        flint_printf("usage: %s p (prime)\n", argv[0]);
        flint_printf("find roots of polynomials over p-adic field Z[p].\n");
        exit(1);
    }
    else
    {
        n = atoi(argv[1]);
        if (!n_is_prime(n))
        {
            flint_printf("%d is not a prime as required for p-adic field.\n",
                         n);
            exit(1);
        }
    };

    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, n);

    padic_ctx_t pctx;
    padic_ctx_init(pctx, p, 128, 128, PADIC_SERIES);

    fmpz_poly_t poly;
    fmpz_poly_init(poly);

    /* use this to read from "poly_data" string. */

    flint_printf("reading polynomials:\n");

    for (slong j = 0; j < sizeof(polys) / sizeof(polys[0]); j++)
    {
        flint_printf("polynomial:\n\n");
        fmpz_poly_set_str(poly, polys[j]);
        fmpz_poly_print_pretty(poly, "x");
        flint_printf("\n\nroots with multiplicity:\n\n");
        padic_roots(poly, pctx, 64);
        flint_printf("\n");
    }

    fmpz_poly_clear(poly);
    fmpz_clear(p);

}
