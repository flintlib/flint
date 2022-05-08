/*
    Copyright (C) 2022 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* NOTE: Assumes that flint.h has been called. */

#ifndef FLINT_TYPES_H
#define FLINT_TYPES_H

#ifdef __cplusplus
extern "C" {
#endif

/* elementary data types ******************************************************/

#ifdef LONGSLONG
typedef long unsigned int       ulong;
typedef long int                slong;
#else
typedef long long unsigned int  ulong;
typedef long long int           slong;
#endif

typedef ulong * ulong_ptr;
typedef const ulong * ulong_srcptr;

typedef slong * slong_ptr;
typedef const slong * slong_srcptr;

typedef ulong flint_bitcnt_t;

/* mock GMP types  *****************************************************/

#ifdef LONGSIZE
typedef long int    mp_mock_size_t;
#else
typedef int         mp_mock_size_t;
#endif

typedef struct
{
    int _mp_alloc;
    int _mp_size;
    ulong *_mp_d;
} __mpz_mock_struct;
typedef __mpz_mock_struct mpz_mock_t[1];
typedef __mpz_mock_struct * mpz_mock_ptr;
typedef const __mpz_mock_struct * mpz_mock_srcptr;

typedef enum
{
    GMP_MOCK_RAND_ALG_DEFAULT = 0,
    GMP_MOCK_RAND_ALG_LC = GMP_MOCK_RAND_ALG_DEFAULT
} gmp_mock_randalg_t;

typedef struct
{
    mpz_mock_t _mp_seed;
    gmp_mock_randalg_t _mp_alg;
    union {
        void *_mp_lc;
    } _mp_algdata;
} __gmp_mock_randstate_struct;
typedef __gmp_mock_randstate_struct gmp_mock_randstate_t[1];

/* pseudo-random numbers ******************************************************/

typedef struct
{
    gmp_mock_randstate_t gmp_state;
    int gmp_init;
    ulong __randval;
    ulong __randval2;
} flint_rand_s;
typedef flint_rand_s flint_rand_t[1];

/* matrices of doubles ********************************************************/

typedef struct
{
    double * entries;
    slong r;
    slong c;
    double ** rows;
} d_mat_struct;
typedef d_mat_struct d_mat_t[1];

/* modular arithmetic (word-size) *********************************************/

typedef struct
{
   ulong n;
   ulong ninv;
   flint_bitcnt_t norm;
} nmod_t;

typedef struct
{
    ulong * coeffs;
    slong alloc;
    slong length;
    nmod_t mod;
} nmod_poly_struct;
typedef nmod_poly_struct nmod_poly_t[1];

typedef struct
{
    ulong * entries;
    slong r;
    slong c;
    ulong ** rows;
    nmod_t mod;
}
nmod_mat_struct;
typedef nmod_mat_struct nmod_mat_t[1];

/* integers *******************************************************************/

typedef slong fmpz;
typedef fmpz fmpz_t[1];

typedef struct
{
   ulong * dinv;
   slong n;
   flint_bitcnt_t norm;
} fmpz_preinvn_struct;
typedef fmpz_preinvn_struct fmpz_preinvn_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
} fmpz_poly_struct;
typedef fmpz_poly_struct fmpz_poly_t[1];

typedef struct
{
    fmpz * entries;
    slong r;
    slong c;
    fmpz ** rows;
} fmpz_mat_struct;
typedef fmpz_mat_struct fmpz_mat_t[1];

/* modular arithmetic (multiprecision) ****************************************/

typedef struct fmpz_mod_ctx {
    fmpz_t n;
    void (* add_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    void (* sub_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    void (* mul_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const struct fmpz_mod_ctx *);
    nmod_t mod;
    ulong n_limbs[3];
    ulong ninv_limbs[3];
} fmpz_mod_ctx_struct;
typedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
} fmpz_mod_poly_struct;
typedef fmpz_mod_poly_struct fmpz_mod_poly_t[1];

typedef struct
{
    fmpz_mat_t mat;
    fmpz_t mod;
}
fmpz_mod_mat_struct;
typedef fmpz_mod_mat_struct fmpz_mod_mat_t[1];

/* rationals ******************************************************************/

typedef struct
{
    fmpz num;
    fmpz den;
}
fmpq;
typedef fmpq fmpq_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    fmpz_t den;
} fmpq_poly_struct;
typedef fmpq_poly_struct fmpq_poly_t[1];

typedef struct
{
    fmpq * entries;
    slong r;
    slong c;
    fmpq ** rows;
} fmpq_mat_struct;
typedef fmpq_mat_struct fmpq_mat_t[1];

/* rational functions over the integers ***************************************/

typedef struct
{
    fmpz_poly_struct * num;
    fmpz_poly_struct * den;
} fmpz_poly_q_struct;
typedef fmpz_poly_q_struct fmpz_poly_q_t[1];

/* p-adic numbers *************************************************************/

typedef struct
{
    fmpz u;
    slong v;
    slong N;
} padic_struct;
typedef padic_struct padic_t[1];

enum padic_print_mode
{
    PADIC_TERSE, 
    PADIC_SERIES, 
    PADIC_VAL_UNIT
};
typedef struct
{
    fmpz_t p;
    double pinv;
    fmpz * pow;
    slong min;
    slong max;
    enum padic_print_mode mode;
} padic_ctx_struct;
typedef padic_ctx_struct padic_ctx_t[1];

typedef struct
{
    fmpz * coeffs;
    slong alloc;
    slong length;
    slong val;
    slong N;
} padic_poly_struct;
typedef padic_poly_struct padic_poly_t[1];

typedef struct
{
    fmpz_mat_struct mat;
    slong val;
    slong N;
} padic_mat_struct;
typedef padic_mat_struct padic_mat_t[1];

/* q-adic numbers *************************************************************/

typedef padic_poly_struct qadic_struct;
typedef padic_poly_t qadic_t;

typedef struct
{
    padic_ctx_struct pctx;
    fmpz * a;
    slong * j;
    slong len;
    char * var;
} qadic_ctx_struct;
typedef qadic_ctx_struct qadic_ctx_t[1];

/* finite fields *************************************************************/

typedef fmpz_poly_struct fq_struct;
typedef fmpz_poly_t fq_t;

typedef struct
{
    fmpz_mod_ctx_t ctxp;
    int sparse_modulus;
    int is_conway;
    fmpz * a;
    slong * j;
    slong len;
    fmpz_mod_poly_t modulus;
    fmpz_mod_poly_t inv;
    char *var;
} fq_ctx_struct;
typedef fq_ctx_struct fq_ctx_t[1];

typedef struct
{
    fq_struct * coeffs;
    slong alloc;
    slong length;
} fq_poly_struct;
typedef fq_poly_struct fq_poly_t[1];

typedef struct
{
    fq_struct * entries;
    slong r;
    slong c;
    fq_struct ** rows;
} fq_mat_struct;
typedef fq_mat_struct fq_mat_t[1];

typedef nmod_poly_t fq_nmod_t;
typedef nmod_poly_struct fq_nmod_struct;

typedef struct
{
    fmpz p;
    nmod_t mod;
    int sparse_modulus;
    int is_conway;
    ulong * a;
    slong * j;
    slong len;
    nmod_poly_t modulus;
    nmod_poly_t inv;
    char * var;
} fq_nmod_ctx_struct;
typedef fq_nmod_ctx_struct fq_nmod_ctx_t[1];

typedef struct
{
    fq_nmod_struct * coeffs;
    slong alloc;
    slong length;
} fq_nmod_poly_struct;
typedef fq_nmod_poly_struct fq_nmod_poly_t[1];

typedef struct
{
    fq_nmod_struct * entries;
    slong r;
    slong c;
    fq_nmod_struct ** rows;
} fq_nmod_mat_struct;
typedef fq_nmod_mat_struct fq_nmod_mat_t[1];

typedef struct
{
    ulong value;
} fq_zech_struct;
typedef fq_zech_struct fq_zech_t[1];

typedef struct
{
    ulong qm1;              /* q - 1 */
    ulong qm1o2;            /* (q - 1) / 2 or 1 when p == 2 */
    ulong qm1opm1;          /* (q - 1) / (p - 1) */
    ulong p;
    double ppre;
    ulong prime_root;       /* primitive root for prime subfield */
    ulong * zech_log_table;
    ulong * prime_field_table;
    ulong * eval_table;
    fq_nmod_ctx_struct * fq_nmod_ctx;
    int owns_fq_nmod_ctx;
    int is_conway;
} fq_zech_ctx_struct;
typedef fq_zech_ctx_struct fq_zech_ctx_t[1];

typedef struct
{
    fq_zech_struct * coeffs;
    slong alloc;
    slong length;
} fq_zech_poly_struct;
typedef fq_zech_poly_struct fq_zech_poly_t[1];

typedef struct
{
    fq_zech_struct * entries;
    slong r;
    slong c;
    fq_zech_struct ** rows;
} fq_zech_mat_struct;
typedef fq_zech_mat_struct fq_zech_mat_t[1];

typedef union fq_default_struct
{
    fq_t fq;
    fq_nmod_t fq_nmod;
    fq_zech_t fq_zech;
    ulong nmod;
    fmpz_t fmpz_mod;
} fq_default_struct;
typedef fq_default_struct fq_default_t[1];

typedef struct
{
    int type;
    union ctx
    {
        fq_ctx_t fq;
        fq_nmod_ctx_t fq_nmod;
        fq_zech_ctx_t fq_zech;
        struct {
            nmod_t mod;
            ulong a;    /* minpoly is x - a */
        } nmod;
        struct {
            fmpz_mod_ctx_t mod;
            fmpz_t a;       /* minpoly is x - a */
        } fmpz_mod;
    } ctx;
} fq_default_ctx_struct;
typedef fq_default_ctx_struct fq_default_ctx_t[1];

typedef union fq_default_poly_struct
{
    fq_poly_t fq;
    fq_nmod_poly_t fq_nmod;
    fq_zech_poly_t fq_zech;
    nmod_poly_t nmod;
    fmpz_mod_poly_t fmpz_mod;
}
fq_default_poly_struct;
typedef fq_default_poly_struct fq_default_poly_t[1];

typedef union fq_default_mat_struct
{
    fq_mat_t fq;
    fq_nmod_mat_t fq_nmod;
    fq_zech_mat_t fq_zech;
    nmod_mat_t nmod;
    fmpz_mod_mat_t fmpz_mod;
} fq_default_mat_struct;
typedef fq_default_mat_struct fq_default_mat_t[1];

/* matrices of polynomials ****************************************************/

typedef struct
{
    nmod_poly_struct * entries;
    slong r;
    slong c;
    nmod_poly_struct ** rows;
    ulong modulus;
} nmod_poly_mat_struct;
typedef nmod_poly_mat_struct nmod_poly_mat_t[1];

typedef struct
{
    fmpz_poly_struct * entries;
    slong r;
    slong c;
    fmpz_poly_struct ** rows;
} fmpz_poly_mat_struct;
typedef fmpz_poly_mat_struct fmpz_poly_mat_t[1];

/* factorization **************************************************************/

typedef struct
{
    int sign;
    fmpz * p;
    ulong * exp;
    slong alloc;
    slong num;
} fmpz_factor_struct;
typedef fmpz_factor_struct fmpz_factor_t[1];

typedef struct
{
    nmod_poly_struct * p;
    slong * exp;
    slong num;
    slong alloc;
} nmod_poly_factor_struct;
typedef nmod_poly_factor_struct nmod_poly_factor_t[1];

typedef struct {
    fmpz c;
    fmpz_poly_struct *p;
    slong *exp;
    slong num;
    slong alloc;
} fmpz_poly_factor_struct;
typedef fmpz_poly_factor_struct fmpz_poly_factor_t[1];

typedef struct
{
    fmpz_mod_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fmpz_mod_poly_factor_struct;
typedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1];

typedef struct
{
    fq_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fq_poly_factor_struct;
typedef fq_poly_factor_struct fq_poly_factor_t[1];

typedef struct
{
    fq_nmod_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fq_nmod_poly_factor_struct;
typedef fq_nmod_poly_factor_struct fq_nmod_poly_factor_t[1];

typedef struct
{
    fq_zech_poly_struct * poly;
    slong * exp;
    slong num;
    slong alloc;
} fq_zech_poly_factor_struct;
typedef fq_zech_poly_factor_struct fq_zech_poly_factor_t[1];

typedef union fq_default_poly_factor_struct
{
    fq_poly_factor_t fq;
    fq_nmod_poly_factor_t fq_nmod;
    fq_zech_poly_factor_t fq_zech;
    nmod_poly_factor_t nmod;
    fmpz_mod_poly_factor_t fmpz_mod;
} fq_default_poly_factor_struct;
typedef fq_default_poly_factor_struct fq_default_poly_factor_t[1];

/* multi-variate polynomials **************************************************/

typedef enum {
    ORD_LEX,
    ORD_DEGLEX,
    ORD_DEGREVLEX
} ordering_t;
typedef struct
{
    slong nvars;    /* number of variables */
    slong nfields;  /* number of fields in exponent vector */
    ordering_t ord; /* monomial ordering */
    int deg;        /* is ord a degree ordering? */
    int rev;        /* is ord a reversed ordering? */
    slong lut_words_per_exp[FLINT_BITS];
    unsigned char lut_fix_bits[FLINT_BITS]; /* FLINT_BITS < 256 */
} mpoly_ctx_struct;
typedef mpoly_ctx_struct mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    nmod_t mod;
} nmod_mpoly_ctx_struct;
typedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1];

typedef struct
{
    ulong * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in ulong units */
    slong exps_alloc;       /* abs size in ulong units */
} nmod_mpoly_struct;
typedef nmod_mpoly_struct nmod_mpoly_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
} fmpz_mpoly_ctx_struct;
typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

typedef struct
{
   fmpz * coeffs; /* alloc fmpzs */
   ulong * exps;
   slong alloc;
   slong length;
   flint_bitcnt_t bits;     /* number of bits per exponent */
} fmpz_mpoly_struct;
typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fmpz_mod_ctx_t ffinfo;
} fmpz_mod_mpoly_ctx_struct;
typedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1];

typedef struct
{
    fmpz * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in ulong units */
    slong exps_alloc;       /* abs size in ulong units */
} fmpz_mod_mpoly_struct;
typedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1];

typedef struct
{
    fmpz_mpoly_ctx_t zctx;
} fmpq_mpoly_ctx_struct;
typedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fq_nmod_ctx_t fqctx;
} fq_nmod_mpoly_ctx_struct;
typedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1];

typedef struct {
    ulong * coeffs;
    ulong * exps;
    slong length;
    flint_bitcnt_t bits;    /* number of bits per exponent */
    slong coeffs_alloc;     /* abs size in ulong units */
    slong exps_alloc;       /* abs size in ulong units */
} fq_nmod_mpoly_struct;
typedef fq_nmod_mpoly_struct fq_nmod_mpoly_t[1];

typedef struct
{
    mpoly_ctx_t minfo;
    fq_zech_ctx_t fqctx;
} fq_zech_mpoly_ctx_struct;
typedef fq_zech_mpoly_ctx_struct fq_zech_mpoly_ctx_t[1];

typedef struct
{
    fq_zech_struct * coeffs;
    ulong * exps;
    slong alloc;
    slong length;
    flint_bitcnt_t bits;     /* number of bits per exponent */
} fq_zech_mpoly_struct;
typedef fq_zech_mpoly_struct fq_zech_mpoly_t[1];

/* multi-variate polynomial factorization *************************************/

typedef struct
{
    ulong constant;
    nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} nmod_mpoly_factor_struct;
typedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1];

typedef struct
{
    fmpz_t constant;
    fmpz_t constant_den;        /* should be one after normal operations */
    fmpz_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fmpz_mpoly_factor_struct;
typedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1];

typedef struct
{
    fq_nmod_t constant;
    fq_nmod_mpoly_struct * poly;
    fmpz * exp;
    slong num;
    slong alloc;
} fq_nmod_mpoly_factor_struct;
typedef fq_nmod_mpoly_factor_struct fq_nmod_mpoly_factor_t[1];

/* APRCL primality testing ****************************************************/

typedef struct
{
    fmpz_mod_poly_t *polys;
    ulong p;
    ulong q;
    fmpz_mod_ctx_t ctx;
} _unity_zpq;
typedef _unity_zpq unity_zpq[1];

typedef struct
{
    fmpz_mod_poly_t poly;
    ulong p;
    ulong exp;
    fmpz_mod_ctx_t ctx;
} _unity_zp;
typedef _unity_zp unity_zp[1];

/* LLL reduction **************************************************************/

typedef enum
{
    GRAM,
    Z_BASIS
} rep_type;
typedef enum
{
    APPROX,
    EXACT
} gram_type;
typedef struct
{
    double delta;
    double eta;
    rep_type rt;
    gram_type gt;
} fmpz_lll_struct;
typedef fmpz_lll_struct fmpz_lll_t[1];

/* miscellaneous **************************************************************/

typedef struct
{
    ulong * coeffs;
    slong alloc;
    slong length;
} n_poly_struct;
typedef n_poly_struct n_poly_t[1];

#ifdef __cplusplus
}
#endif

#endif
