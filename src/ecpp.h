/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef ECPP_H
#define ECPP_H

#include <stdio.h>
#include "fmpz_types.h"
#include "fmpz_mod_types.h"
#include "fmpz_mod_poly.h"
#include "gr_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    Elliptic curve primality proving (Atkin-Morain).

    A certificate is a chain of steps. Step i proves that n_i is prime
    assuming that q_i is prime: the curve y^2 = x^3 + a x + b over Z/n_i
    has a point P with m P = O and (m / q) P != O, where q | m and
    q > (n^{1/4} + 1)^2. The next step has n_{i+1} = q_i. The last q
    is small enough to be proved prime by other means (q < 2^64).
*/

typedef struct
{
    fmpz_t n;
    slong D;    /* discriminant used (informational) */
    fmpz_t a;
    fmpz_t b;
    fmpz_t m;
    fmpz_t q;
    fmpz_t x;   /* point P = (x, y) */
    fmpz_t y;
}
ecpp_step_struct;

typedef struct
{
    ecpp_step_struct * steps;
    slong num;
    slong alloc;
}
ecpp_cert_struct;

typedef ecpp_cert_struct ecpp_cert_t[1];

/* affine or Jacobian point on y^2 = x^3 + a x + b modulo n */
typedef struct
{
    fmpz_t X;
    fmpz_t Y;
    fmpz_t Z;   /* Z = 0 is the point at infinity */
}
ecpp_point_struct;

typedef ecpp_point_struct ecpp_point_t[1];

/* certificate */

void ecpp_cert_init(ecpp_cert_t cert);
void ecpp_cert_clear(ecpp_cert_t cert);
ecpp_step_struct * ecpp_cert_push(ecpp_cert_t cert);
void ecpp_cert_pop(ecpp_cert_t cert);
#define ECPP_CERT_FORMAT_FLINT 0
#define ECPP_CERT_FORMAT_PARI 1

void ecpp_cert_print(const ecpp_cert_t cert, int format);
void ecpp_cert_fprint(FILE * file, const ecpp_cert_t cert, int format);
char * ecpp_cert_get_str(const ecpp_cert_t cert, int format);
int ecpp_cert_set_str(ecpp_cert_t cert, const char * str);

/* elliptic curves modulo n */

void ecpp_point_init(ecpp_point_t P);
void ecpp_point_clear(ecpp_point_t P);
void ecpp_point_set_affine(ecpp_point_t P, const fmpz_t x, const fmpz_t y);
int ecpp_point_is_zero(const ecpp_point_t P);

void ecpp_gr_ctx_init(gr_ctx_t gctx, const fmpz_mod_ctx_t ctx);

int ecpp_point_mul_gr(gr_ptr R, gr_srcptr P, const fmpz_t k, gr_srcptr a,
                                                gr_ptr acc, gr_ctx_t ctx);

int ecpp_point_mul(ecpp_point_t R, const ecpp_point_t P, const fmpz_t k,
                        const fmpz_t a, fmpz_t acc, const fmpz_mod_ctx_t ctx);

/* discriminants */

#define ECPP_DISC_MAXFAC 8

/* Tuning parameters of the prover (see prove.c) */

/* expected prime cofactors per round of the discriminant search */
#define ECPP_MIN_PRIME 4.0
/* from this size, the pool admits all discriminants with a cheap class
   field tower rather than only those with a small genus factor */
#define ECPP_BIGPOOL_BITS 2500
/* class number and odd part bounds of the discriminant pool */
#define ECPP_POOL_HMAX 400
#define ECPP_POOL_OMAX 16
/* bound on the size of the discriminants, as a power of two */
#define ECPP_DMAX_BITS(bits) ((bits) < 1500 ? 20 : (bits) < 3000 ? 21 : (bits) < 4000 ? 22 : 23)
/* the primorial for the batch factoring is reduced in this many chunks
   when threads are used */
#define ECPP_PRIMORIAL_CHUNKS 8

/* from this size the class field tower carries Kummer data: levels of
   prime degree p >= 5 are descended with a p-th root when n = +-1 mod p */
#define ECPP_KUMMER_BITS 1500

typedef struct
{
    slong D;
    slong h;
    slong g;                /* number of genus characters */
    slong o;                /* odd part h / 2^(g-1) of the class number */
    slong q0;               /* 1, -4, 8 or -8: D = q0 prod p^* */
    slong nfac;             /* odd primes dividing D, as indices into the prime list */
    unsigned short fac[ECPP_DISC_MAXFAC];
    double cost;            /* estimated realisation cost, see ecpp_disc_cost */
}
ecpp_disc_struct;

#ifndef ECPP_DISC_COST_MAX
#define ECPP_DISC_COST_MAX 20.0
#endif

double ecpp_disc_cost(const ecpp_disc_struct * d);
double ecpp_disc_cost_v(const ecpp_disc_struct * d, int veven);
double ecpp_disc_cost_n(const ecpp_disc_struct * d, int veven, const ulong * nmodp, slong bits);
int ecpp_disc_use_tower(const ecpp_disc_struct * d, int veven, slong * Dt);

slong ecpp_disc_table(ecpp_disc_struct ** table, const ulong * primes,
            slong nprimes, slong Dmax, slong hmax, slong omax, double costmax);

/* number theory helpers */

int ecpp_cornacchia(fmpz_t t, fmpz_t v, const fmpz_t n, slong D, const fmpz_t sqrtD);

int ecpp_root_radicals(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx);

int ecpp_poly_root(fmpz_t x, const fmpz_mod_poly_t f, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx);

/* options of _ecpp_class_poly_tower: compute the Kummer data (done by
   ecpp_class_poly_tower from ECPP_KUMMER_BITS), use j rather than the Weber
   invariant */
#define ECPP_TOWER_KUMMER 1
#define ECPP_TOWER_J 2

int _ecpp_class_poly_tower(fmpz_t j, slong D, int flags, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx);

int ecpp_class_poly_tower(fmpz_t j, slong D, flint_rand_t state,
                                                        const fmpz_mod_ctx_t ctx);

int ecpp_class_poly_genus(fmpz_mod_poly_t F, slong D, const slong * pstar, slong g,
                                const fmpz * sqrts, const fmpz_mod_ctx_t ctx);

/* proving and verifying */

int ecpp_prove(ecpp_cert_t cert, const fmpz_t n);

/* prints the steps of the proof as they are found (0 to disable) */
void ecpp_set_verbose(int verbose);
int ecpp_verify(const ecpp_cert_t cert, const fmpz_t n);
int ecpp_verify_step(const ecpp_step_struct * s);
int ecpp_is_prime(const fmpz_t n);

#ifdef __cplusplus
}
#endif

#endif
