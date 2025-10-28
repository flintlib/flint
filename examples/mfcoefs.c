/*
    Copyright (C) 2025 Pascal Molin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    FLINT program demonstrating computation of modular form coefficients
*/

#include <stdlib.h>
#include <flint.h>

#if FLINT_HAVE_FFT_SMALL

#include <stdio.h>
#include <string.h>
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fft_small.h"
#include "fmpz.h"
#include "arith.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "fmpz_mat.h"
#include "dirichlet.h"
#include "profiler.h"

/*
 Fast computation of modular forms coefficients via
 representation as products of Eisenstein series.

 A modular form f in S_k(N) is given as

 f = c_0 E2(N) + sum_{i=1}^{n} c_i E1(chi_i)E1(chi_i^-1)

 with E2(N) the usual Eisenstein series of weight 2 and level N
 and E1(chi) the weight one Eisenstein series of character chi.


 The coefficients of Eisenstein series are computed as lines in
 a nmod matrix (the modulus having enough roots of unity).

 For big lengths this matrix may be compressed to keep only
 prime indices.

 E1 a2 a3 a5 ...
 E2 a2 a3 a5 ...

 A form is a combination of Ei whose coefficients belong
 to some number ring.

 It is given as a matrix whose rows correspond to integral
 basis (and columns to generators).

 When several forms are required the matrices may be
 concatenated.

 A form is then given by a matrix whose lines correspond to an
 integral basis of the Hecke field.

 We choose to transpose the output so that each line corresponds
 to a Fourier coefficient. By default only coefficients of
 prime index are output.
*/

/* data format, one character */
struct mf_eis_space {
    const slong N;       // level
    const slong k;       // weight
                         // TODO: space character
    const slong nchi;    // number of Dirichlet character used
    const slong * chi;   // characters by Conrey index mod N

    const slong ord;     // order of root of unity
    const ulong modp;    // fft prime used for expression
    const ulong z;       // root of unity for character

    const slong num;     // number of Eisenstein generators
    const slong * l;     // weight
    const slong * c;     // character index (-1 for Ek)
    const slong * d;     // Bd operator

    const slong rank;    // rank of output basis
    const ulong * basis; // conversion matrix from generators to basis
                         // rank * num, could be nmod_mat
                         // each *line* encodes one coefficient on hecke basis

                         // FIXME: may skip this, depends on form
    const slong dim;     // number of actual forms
    const slong * deg;   // degree of form
    const char ** poly;  // Hecke polynomial
    const char ** zk;    // integral basis used as string
};
typedef struct mf_eis_space mf_space_t[1];

// [11, [ [2, Mod(1,11)], -3/2; [1, Mod(-1,11), 1, -1], 5/2 ], [-1,2] ]
const struct mf_eis_space mf11 = {
    11, 2,

    2, (const slong[]){ 1, 10 },

    2, 928284166586369, 928284166586368,

    2,
    (const slong[]){ 2, 1 },
    (const slong[]){ 0, 1 },
    (const slong[]){ 1, 1 },

    1,
    (const ulong[]){ 464142083293183,464142083293187 },

    1,
    (const slong[]){ 1 },
    (const char*[]){ "y-1" },
    (const char*[]){ "[1]" }
};
const struct mf_eis_space mf23 = {
    23, 2,
    3, (const slong[]){ 1,5,22 },
    22, 744953487556609, 327890428415624,
    3,
    (const slong[]){ 2,1,1 }, // weight
    (const slong[]){ 0,2,1 }, // index
    (const slong[]){ 1,1,1 }, // d
    2, // rank
    (const ulong[]){ 716793017771430, 239561039144421, 699801604575949, 456958153133828, 398747113901656, 135455648941980 },
    1, // forms
    (const slong[]){ 2 }, // hecke degree
    (const char *[]){ "y^2 - y - 1" }, // hecke pol
    (const char *[]){ "[1, y]" } // basis
};
const struct mf_eis_space mf41 = {
    41, 2,
    4, (const slong[]){ 1,3,6,11 },
    40, 587207928709121, 227071490884881,
    4,
    (const slong[]){ 2,1,1,1 }, // weight
    (const slong[]){ 0,1,2,3 }, // index
    (const slong[]){ 1,1,1,1 }, // d
    3, // rank
    (const ulong[]){ 69233414194803,56495677524383,395984113667253,583670866952449,275432891021865,53973456700211,237459359376836,179534281490698,24071217223589,301634451880304,383565831106330,246106924681630 },
    1, // forms
    (const slong[]){ 3 }, // hecke degree
    (const char *[]){ "y^3 - y^2 - 3*y + 1" }, // hecke pol
    (const char *[]){ "[1, y, y^2 - y - 2]" } // basis
};
const struct mf_eis_space mf131 = {
    131, 2,
    11, (const slong[]){ 1,2,8,10,14,17,22,23,26,29,130 },
    130, 708959514132481, 586822695664532,
    11,
    (const slong[]){ 2,1,1,1,1,1,1,1,1,1,1 }, // weight
    (const slong[]){ 0,10,1,2,3,4,5,6,7,8,9 }, // index
    (const slong[]){ 1,1,1,1,1,1,1,1,1,1,1 }, // d
    10, // rank
    (const ulong[]){ 419689470091886,255574599084099,699803187026726,18954044708734,429280161050971,651038748314273,36012914137727,615032388621403,455482254176853,440054925650024,453324867859397,707300356200848,289775764014946,126106408262275,567757971889846,172729569381415,178457093228947,610454617435762,354337675902414,627244157401223,200118908326692,351243427992810,504517361880407,532681765365538,700793904701229,675002619139711,250817945653665,573230180729341,304022589480280,605743041656159,389199218609257,117990933144826,65583587489595,148277197435256,584617780236622,359091558741812,503610707037325,44996173133172,597019225625155,464875015435713,262216128236979,327810432402918,31571984391542,489521021552563,326818751124043,252971092489530,574752599399524,404758443217259,499159734767711,707216834848433,537073232091478,21707873187690,501314817530805,387593518576946,382371590806484,358499269571204,566247500049862,33492497523770,164694472266478,30352692409355,23058362049849,173681518747944,568962409034955,268330386457800,697323618141921,460137928554939,592029366340193,445464762732680,158353244673956,410738896840342,101682818633393,290875952427265,240079916583990,565712846920067,106148979826246,217795480112622,668526374487082,182225090104431,518982715478679,262043890461584,588594365589034,615083302178597,7048890429295,403880696424916,409785585213979,270807854777441,224948373031878,104855449448709,137014408822430,450391725186438,431316661114234,153676391761881,468255031512777,608535734629677,490132552672039,389020265991828,627503543965891,191160835527493,383466618974635,477424209326765,371231961617645,399333872263146,108653058455188,702583732943232,251852418190874,248298225578617,53018682677347,422080036481016,376211299668354,142509208003682 },
    1, // forms
    (const slong[]){ 10 }, // hecke degree
    (const char *[]){ "y^10 - 18*y^8 + 2*y^7 + 111*y^6 - 18*y^5 - 270*y^4 + 28*y^3 + 232*y^2 + 16*y - 32" }, // hecke pol
    (const char *[]){ "[1, y, 1/8*y^8 - 2*y^6 + 81/8*y^4 - 67/4*y^2 + 5, 1/16*y^9 - 9/8*y^7 + 1/8*y^6 + 111/16*y^5 - 9/8*y^4 - 135/8*y^3 + 7/4*y^2 + 27/2*y + 1, 1/16*y^9 + 1/8*y^8 - 7/8*y^7 - 15/8*y^6 + 55/16*y^5 + 9*y^4 - 17/8*y^3 - 14*y^2 - 4*y + 2, -1/8*y^9 - 1/8*y^8 + 2*y^7 + 3/2*y^6 - 85/8*y^5 - 41/8*y^4 + 81/4*y^3 + 17/4*y^2 - 9*y + 1, -1/8*y^9 - 1/8*y^8 + 2*y^7 + 3/2*y^6 - 81/8*y^5 - 41/8*y^4 + 63/4*y^3 + 17/4*y^2 - 2*y + 3, 1/16*y^9 + 1/8*y^8 - 7/8*y^7 - 13/8*y^6 + 59/16*y^5 + 27/4*y^4 - 27/8*y^3 - 19/2*y^2 - 9/2*y + 1, 1/4*y^6 + 1/4*y^5 - 9/4*y^4 - 9/4*y^3 + 7/2*y^2 + 9/2*y + 2, 1/4*y^6 + 1/4*y^5 - 11/4*y^4 - 9/4*y^3 + 8*y^2 + 9/2*y - 4]" } // basis
};
const struct mf_eis_space mf13_4 = {
    13, 4,
    4, (const slong[]){ 2,3,4,12 },
    12, 802974200758273, 102569179906582,
    4,
    (const slong[]){ 1,2,2,2 }, // weight
    (const slong[]){ 0,1,2,3 }, // index
    (const slong[]){ 1,1,1,1 }, // d
    2, // rank
    (const ulong[]){ 290398498851043,710805682556811,622389199807714,493494834874998,0,154418115530437,669145167298561,350014395202324 },
    2, // forms
    (const slong[]){ 1, 2 }, // hecke degree
    (const char *[]){ "y", "y^2 - y - 4" }, // hecke pol
    (const char *[]){ "[1]", "[1, y]" } // basis
};

/* character values mod N */
struct mf_char_ctx {
    ulong q;
    ulong a;
    ulong * chivec;
};
typedef struct mf_char_ctx mf_char_ctx_t[1];
typedef struct mf_char_ctx * mf_char_ctx_ptr;

/* Compute the list of character values as powers of z mod p,
 * assume z has exact order ord */
void
dirichlet_chi_vec_nmod(ulong *v, slong nv, const dirichlet_group_t G,
        const dirichlet_char_t chi, ulong ord, ulong z, nmod_t mod)
{
    slong k;
    dirichlet_chi_vec_order(v, G, chi, ord, nv);
    for (k = 0; k < nv; k++)
    {
        if (v[k] == DIRICHLET_CHI_NULL)
            v[k] = 0;
        else
            v[k] = nmod_pow_ui(z, v[k], mod);
    }
}

/* store values chi(k) */
void
mf_char_ctx_init(mf_char_ctx_t ctx, const dirichlet_group_t G, slong a, ulong ord, ulong z, nmod_t mod)
{
    dirichlet_char_t chi;
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, a);
    ctx->q = G->q;
    ctx->a = a;
    ctx->chivec = (ulong *)flint_malloc(G->q*sizeof(ulong));
    dirichlet_chi_vec_nmod(ctx->chivec, G->q, G, chi, ord, z, mod);
    dirichlet_char_clear(chi);
}
void
mf_char_ctx_clear(mf_char_ctx_t ctx)
{
    flint_free(ctx->chivec);
}

/* change chi -> chi^(-1) */
void
mf_char_ctx_dual(mf_char_ctx_t ctx, nmod_t mod)
{
    slong k;
    for (k = 0; k < ctx->q; k++)
        if (ctx->chivec[k])
            ctx->chivec[k] = nmod_inv(ctx->chivec[k], mod);
}

/* precompute coprime decomposition k = p^e*m, p smallest prime */
typedef struct {
    slong pe;
    slong m;
} pem_struct;
typedef pem_struct * pem_ptr;
typedef const pem_ptr pem_srcptr;
struct rough {
    ulong m;
    struct rough * prev;
    struct rough * next;
};
typedef struct rough * rough_ptr;

void
pem_init_rough(pem_ptr tab, slong lim, slong len)
{
    ulong m, pe, pem, len1 = (len-1) / 2;
    rough_ptr rough, p1, m1;

    /* p=2 can be done separately */
    for (m = 3; m <= len1; m += 2)
        for (pe = 2, pem = 2*m; pem < len; pem <<= 1, pe <<= 1)
            tab[pem].pe = pe, tab[pem].m = m;

    /* now need 2-rough (ie odd) numbers up to len / 3 */
    len1 = (len-1) / 3;
    rough = flint_malloc((2 + len1/2) * sizeof(struct rough));
    rough->m = 1;
    rough->prev = NULL;
    rough->next = rough + 1;
    for (m1 = rough + 1, m = 3; m <= len1; m1++, m += 2)
    {
        m1->m = m;
        m1->prev = m1 - 1;
        m1->next = m1 + 1;
    }
    /* terminate */
    m1->m = len;
    m1->prev = m1 - 1;
    m1->next = NULL;

    for (p1 = rough + 1; p1->m <= lim; p1 = p1->next)
    {
        slong p = p1->m, pe;
        /* skip prime powers p^e up to len1 */
        for (pe = p; pe <= len1; pe *= p)
        {
            rough_ptr pe1 = rough + (pe>>1);
            pe1->next->prev = pe1->prev;
            pe1->prev->next = pe1->next;
        }
        /* then loop on p-rough numbers (p>=3) */
        for (pe = p; pe <= len1; pe *= p)
        {
            ulong lim = (len-1) / pe;
            /* loop m in p-rough numbers */
            for (m1 = p1->next; m1->m <= lim; m1 = m1->next)
            {
                slong pem = pe * m1->m;
                tab[pem].pe = pe;
                tab[pem].m = m1->m;
                /* update links to skip pem */
                if (pem >= len1) continue;
                rough_ptr pem1 = rough + (pem>>1);
                pem1->next->prev = pem1->prev;
                pem1->prev->next = pem1->next;
            }
        }
    }
    /* if lim = sqrt(len), rough now links primes > lim */
    flint_free(rough);
}

typedef struct {
    ulong n;
    ulong a;
    ulong b;
} coprime_t;
typedef coprime_t * coprime_ptr;

coprime_ptr
coprime_table_init(slong * size, slong len)
{
    slong k, n;
    coprime_ptr fac;
    pem_ptr tab = flint_malloc(len * sizeof(pem_struct));
    for (k = 0; k < len; k++)
        tab[k].pe = tab[k].m = 0;
    pem_init_rough(tab, n_sqrt(len), len);
    fac = flint_malloc(len * sizeof(coprime_t));
    for (n = 0, k = 1; k < len; k++)
        if(tab[k].m)
            fac[n++] = (coprime_t){ .n = k, .a = tab[k].pe, .b = tab[k].m };
    flint_free(tab);
    *size = n;
    fac = flint_realloc(fac, n * sizeof(coprime_t));
    return fac;
}

/* stupid: doing slong -> fmpq -> ulong internally, no matter */
void
nmod_poly_bernoulli(nmod_poly_t pol, slong n)
{
    fmpq_poly_t polq;
    fmpq_poly_init(polq);
    arith_bernoulli_polynomial(polq, n);
    fmpq_poly_get_nmod_poly(pol, polq);
    fmpq_poly_clear(polq);
}

ulong
eisenstein_constant(slong k, mf_char_ctx_t psi, nmod_t mod)
{
    slong a;
    ulong m, e0 = 0, N = psi->q;
    nmod_poly_t pol;

    if (psi->a == 1)
        return 0;

    nmod_poly_init_mod(pol, mod);
    nmod_poly_bernoulli(pol, k);

    /* - 1/2 * N^(k-1)/k * sum psi(a) pol(a/N) */
    for (a = 0; a < psi->q; a++)
    {
        ulong paN = nmod_poly_evaluate_nmod(pol, nmod_div(a, N, mod));
        e0 = nmod_add(e0, nmod_mul(psi->chivec[a], paN, mod), mod);
    }
    m = nmod_div(nmod_pow_ui(N, k-1, mod), mod.n - 2 * k, mod);
    e0 = nmod_mul(e0, m, mod);
    nmod_poly_clear(pol);
    return e0;
}

/* Eisenstein series: direct algorithm */
void
_nmod_poly_eisenstein_series(nn_ptr z, slong len, slong k, mf_char_ctx_t psi, const coprime_ptr tab, slong size, nmod_t mod)
{
    slong p;
    n_primes_t iter;
    n_primes_init(iter);

    z[0] = eisenstein_constant(k, psi, mod);
    z[1] = 1;
    /* first expand Euler factors */
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        slong pe, pe1 = 1;
        ulong cp = psi ? psi->chivec[p % psi->q] : 1;
        if (k > 1)
            cp = nmod_mul(cp, nmod_pow_ui(p, k-1, mod), mod);
        z[p] = nmod_add(cp, 1, mod);
        for (pe = p, pe1 = 1; pe < len; pe1 = pe, pe *= p)
            z[pe] = nmod_add(nmod_mul(z[pe1], cp, mod), 1, mod);
    }
    /* then fill composite */
    for (k = 0; k < size; k++)
        z[tab[k].n] = nmod_mul(z[tab[k].a], z[tab[k].b], mod);
}

typedef struct {
    slong count_euler;
    slong count_prod;
    slong cpu_euler;
    slong wall_euler;
    slong cpu_prod;
    slong wall_prod;
    slong cpu_total;
    slong wall_total;
} mf_timer;

/* Modular form */
void
nmod_vec_get_nmod_vec_primes(nn_ptr a, nn_srcptr g, slong len)
{
    n_primes_t iter;
    slong j, p;

    /* assume a and g large enough */
    n_primes_init(iter);
    for (j = 0, p = n_primes_next(iter); p < len; j++, p = n_primes_next(iter))
        a[j] = g[p];
    n_primes_clear(iter);
}

void
nmod_mat_modular_form_series(nmod_mat_t a, const mf_space_t mf, slong len, mf_timer * timer)
{
    nmod_t mod;
    coprime_ptr tab = NULL;
    dirichlet_group_t G;
    mf_char_ctx_ptr char_ctx;
    nn_ptr g1, g2, g12;
    nmod_mat_t eis, basis;
    mpn_ctx_t fft_ctx;
    slong size, cols, i, j;
    timeit_t t;
    /* init */

    nmod_mat_set_mod(a, mf->modp);
    nmod_init(&mod, mf->modp);

    /* sanity checks */
    cols = n_prime_pi(len);
    FLINT_ASSERT(n_prime_pi(len) == nmod_mat_ncols(a));
    FLINT_ASSERT(mf->rank == nmod_mat_nrows(a));
    FLINT_ASSERT(n_trailing_zeros(mf->modp-1) > n_clog2(len));


    /* precompute chars */
    char_ctx = flint_malloc(2 * mf->nchi * sizeof(struct mf_char_ctx));
    dirichlet_group_init(G, mf->N);
    for (i = 0; i < mf->nchi; i++)
    {
        mf_char_ctx_init(char_ctx + 2*i, G, mf->chi[i], mf->ord, mf->z, mod);
        mf_char_ctx_init(char_ctx + 2*i+1, G, n_invmod(mf->chi[i],mf->N), mf->ord, mf->z, mod);
    }
    dirichlet_group_clear(G);

    /* store eisenstein expansions */
    nmod_mat_init(eis, mf->num, cols, mf->modp);

    /* tabulate composite */
    tab = coprime_table_init(&size, len);
    /*
     Critical: force mpn_mul (_nmod_poly_mul_mid_mpn_ctx)
     to use custom 50 bits prime p = mod.n
     Initialize context fft_ctx accordingly.
     */
    mpn_ctx_init(fft_ctx, mod.n);
    g1 = _nmod_vec_init(len);
    g2 = _nmod_vec_init(len);
    g12 = _nmod_vec_init(len);
    for (i = 0; i < mf->num; i++)
    {
        slong k = mf->l[i], c = mf->c[i];
        nn_ptr row = nmod_mat_entry_ptr(eis, i, 0);
        mf_char_ctx_ptr psi1 = (c == -1) ? NULL : char_ctx + 2 * mf->c[i];
        mf_char_ctx_ptr psi2 = (c == -1) ? NULL : char_ctx + 2 * mf->c[i] + 1;
        if (k < mf->k)
        {
             if (timer)
                timeit_start(t);

            _nmod_poly_eisenstein_series(g1, len, k, psi1, tab, size, mod);
            _nmod_poly_eisenstein_series(g2, len, mf->k - k, psi2, tab, size, mod);

            if (timer)
            {
                timeit_stop(t);
                timer->count_euler += 2;
                timer->cpu_euler += t->cpu;
                timer->wall_euler += t->wall;
                timeit_start(t);
            }

            _nmod_poly_mul_mid_mpn_ctx(g12, 0, len, g1, len, g2, len, mod, fft_ctx);

            if (timer)
            {
                timeit_stop(t);
                timer->count_prod += 1;
                timer->cpu_prod += t->cpu;
                timer->wall_prod += t->wall;
            }

            nmod_vec_get_nmod_vec_primes(row, g12, len);
        }
        else
        {
            if (timer)
                timeit_start(t);

            /* TODO: need only prime indices */
            _nmod_poly_eisenstein_series(g12, len, k, psi1, tab, size, mod);

            if (timer)
            {
                timeit_stop(t);
                timer->count_euler += 1;
                timer->cpu_euler += t->cpu;
                timer->wall_euler += t->wall;
                timeit_start(t);
            }

            nmod_vec_get_nmod_vec_primes(row, g12, len);
        }
    }
    _nmod_vec_clear(g1);
    _nmod_vec_clear(g2);
    _nmod_vec_clear(g12);
    mpn_ctx_clear(fft_ctx);
    flint_free(tab);
    for (i = 0; i < 2 * mf->nchi; i++)
        mf_char_ctx_clear(char_ctx + i);
    flint_free(char_ctx);

    /* convert to basis */
    nmod_mat_init(basis, mf->rank, mf->num, mf->modp);
    for (i = 0; i < mf->rank; i++)
        for (j = 0; j < mf->num; j++)
            nmod_mat_entry(basis, i, j) = mf->basis[i*mf->num+j];

    nmod_mat_mul(a, basis, eis);

    nmod_mat_clear(basis);
    nmod_mat_clear(eis);
}
void
fmpz_mat_set_transpose_nmod_mat(fmpz_mat_t b, const nmod_mat_t a)
{
    slong i, j;

    for (i = 0; i < a->r; i++)
        for (j = 0; j < a->c; j++)
            fmpz_set_ui_smod(fmpz_mat_entry(b, j, i),
                             nmod_mat_entry(a, i, j), a->mod.n);
}

int usage(int count, const char * fname[])
{
    int i;
    flint_printf("mfcoefs [options] <level> <length>\n");
    flint_printf("where <level> is in");
    for (i = 0; i < count; i++) flint_printf(" %s,", fname[i]);
    flint_printf("\n and <length> is the number of terms to compute\n");
    flint_printf("output coefficients as a matrix, one row per coefficient a_p,");
    flint_printf(" columns indexed by an integral basis of the value field\n");
    flint_printf("options:\n");
    flint_printf(" --raw: raw flint output (matrix size followed by space separated values)\n");
    flint_printf(" --tail <n>: output only last <n> coefficients\n");
    flint_printf(" --time: time each step (implies --tail 0)\n");
    flint_printf(" --bench <nmin> <nmul> <nmax> benchmark all forms.\n");
    return EXIT_FAILURE;
}

int main(int argc, char * argv[])
{
    slong i, j;
    slong len0 = 1000, lmax = 1L<<30, lmul = 100;

    nmod_mat_t a;
    fmpz_mat_t m;
    slong rows, cols, count = 5;
    const char * mf_name[] = { "11", "23", "41" , "131", "13_4" };
    const struct mf_eis_space *f = NULL, mf[] = { mf11 , mf23 , mf41, mf131, mf13_4 };
    int opt_raw = 0, opt_time = 0, opt_bench = 0;
    slong opt_tail = -1;
    mf_timer timer_struct, * timer = NULL;
    timeit_t total_time;

    /* options */
    for (i = 1; i < argc;)
    {
        if (strcmp(argv[i], "--raw") == 0)
            opt_raw = 1, i++;
        else if (strcmp(argv[i], "--tail") == 0 && i + 1 < argc)
            opt_tail = atol(argv[i+1]), i +=2 ;
        else if (strcmp(argv[i], "--time") == 0)
            opt_time = 1, i++;
        else if (strcmp(argv[i], "--bench") == 0)
            opt_bench = 1, i++;
        else break;
    }

    if (opt_bench) opt_time = 1, j = 0;
    if (opt_time) timer = &timer_struct;

    /* form and length */
    if (opt_bench && argc == i + 3)
    {
        len0 = atol(argv[i++]);
        lmul = atol(argv[i++]);
        lmax = atol(argv[i++]);
    }
    else if (!opt_bench && argc == i + 2)
    {
        for (j = 0; j < count; j++)
           if (strcmp(argv[i], mf_name[j]) == 0)
               break;
        i++;
        len0 = lmax = atol(argv[i++]);
    }

    if (argc != i || len0 < 1 || j >= count)
        return usage(count, mf_name);

    for (; j < count; j = opt_bench ? j+1 : count)
    {

        slong len;
        f = mf + j;

        for (len = len0; len <= lmax; len = opt_bench ? len*lmul : lmax+1)
        {

            if (opt_time) {
                timer->count_euler = 0;
                timer->count_prod = 0;
                timer->cpu_euler = 0;
                timer->wall_euler = 0;
                timer->cpu_prod = 0;
                timer->wall_prod = 0;
                timeit_start(total_time);
            }

            cols = n_prime_pi(len);
            rows = f->rank;
            nmod_mat_init(a, rows, cols, f->modp);

            nmod_mat_modular_form_series(a, f, len, timer);

            if (opt_time)
            {
                timeit_stop(total_time);

                if (opt_bench && len == len0)
                {
                    flint_printf("%s & length & euler (x%ld) & prod (x%ld) & total \n",
                            mf_name[j], timer->count_euler, timer->count_prod);
                    flint_printf("   & %ld & %ld & %ld & %ld\n", len,
                            timer->wall_euler, timer->wall_prod, total_time->wall);
                }
                else if (opt_bench)
                {
                    flint_printf("   & %ld & %ld & %ld & %ld\n", len,
                            timer->wall_euler, timer->wall_prod, total_time->wall);
                }
                else
                {
                    flint_printf("euler (x%ld) cpu = %wd ms  wall = %wd ms\n", timer->count_euler, timer->cpu_euler, timer->wall_euler);
                    flint_printf("prod  (x%ld) cpu = %wd ms  wall = %wd ms\n", timer->count_prod, timer->cpu_prod, timer->wall_prod);
                    flint_printf("total cpu = %wd ms  wall = %wd ms\n", total_time->cpu, total_time->wall);
                }
            }
            else
            {
                nmod_mat_t a2;
                nmod_mat_struct * pa = a;
                slong size = cols;

                if (opt_tail >= 0)
                {
                    size = (opt_tail < cols) ? opt_tail : cols;
                    nmod_mat_window_init(a2, a, 0, cols - size, rows, cols);
                    pa = a2;
                }

                fmpz_mat_init(m, size, rows);
                fmpz_mat_set_transpose_nmod_mat(m, pa);

                if (opt_tail >= 0)
                    nmod_mat_window_clear(a2);

                if (opt_raw)
                    fmpz_mat_print(m);
                else
                    fmpz_mat_print_pretty(m);
                flint_printf("\n");

                fmpz_mat_clear(m);
            }
            nmod_mat_clear(a);
        }
    }
}

#else

int main(int argc, char * argv[])
{
    flint_printf("mfcoefs requires the fft_small module.\n");
    return EXIT_FAILURE;
}

#endif
