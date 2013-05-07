/* Pi computation using Chudnovsky's algortithm.

 * Copyright 2002, 2005 Hanhong Xue (macroxue at yahoo dot com)

 * Modified 2005 by Torbjorn Granlund (tege at swox dot com) to allow more than
   2G digits to be computed.  Modified 2009 by Torbjorn Granlund for GMPbench.

 * Modified 2011 by Fredrik Johansson to make reentrant and adapt for
   use in FLINT.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "arith.h"

#define A   13591409
#define B   545140134
#define C   640320
#define D   12

#define BITS_PER_DIGIT   3.32192809488736234787
#define DIGITS_PER_ITER  14.1816474627254776555
#define DOUBLE_PREC      53

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

typedef struct {
    unsigned long max_facs;
    unsigned long num_facs;
    unsigned long *fac;
    unsigned long *pow;
} fac_t[1];

typedef struct {
    long int fac;
    long int pow;
    long int nxt;
} sieve_t;

typedef struct 
{
    sieve_t *sieve;
    long int sieve_size;
    fac_t ftmp, fmul;
    mpz_t gcd;
    mpz_t *pstack, *qstack, *gstack;
    fac_t *fpstack, *fgstack;
    long int top;
}
pi_state_struct;

typedef pi_state_struct pi_state[1];

#define INIT_FACS 32

#define p1 (state->pstack[state->top])
#define q1 (state->qstack[state->top])
#define g1 (state->gstack[state->top])
#define fp1 (state->fpstack[state->top])
#define fg1 (state->fgstack[state->top])

#define p2 (state->pstack[state->top+1])
#define q2 (state->qstack[state->top+1])
#define g2 (state->gstack[state->top+1])
#define fp2 (state->fpstack[state->top+1])
#define fg2 (state->fgstack[state->top+1])

/* r = sqrt(x) */
void
my_sqrt_ui(pi_state state, mpf_t t1, mpf_t t2, mpf_t r, unsigned long x)
{
    unsigned long prec, bits, prec0;

    prec0 = mpf_get_prec(r);

    if (prec0 <= DOUBLE_PREC)
    {
        mpf_set_d(r, sqrt(x));
        return;
    }

    bits = 0;
    for (prec = prec0; prec > DOUBLE_PREC; )
    {
        int bit = prec&1;
        prec = (prec+bit)/2;
        bits = bits*2+bit;
    }

    mpf_set_prec_raw(t1, DOUBLE_PREC);
    mpf_set_d(t1, 1/sqrt(x));

    while (prec < prec0)
    {
        prec *=2;
        if (prec < prec0)
        {
            /* t1 = t1+t1*(1-x*t1*t1)/2; */
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, t1, t1);         /* half x half -> full */
            mpf_mul_ui(t2, t2, x);
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_div_2exp(t2, t2, 1);
            mpf_mul(t2, t2, t1);         /* half x half -> half */
            mpf_set_prec_raw(t1, prec);
            mpf_add(t1, t1, t2);
        }
        else
        {
            break;
        }
        prec -= (bits&1);
        bits /=2;
    }

    /* t2=x*t1, t1 = t2+t1*(x-t2*t2)/2; */
    mpf_set_prec_raw(t2, prec0/2);
    mpf_mul_ui(t2, t1, x);
    mpf_mul(r, t2, t2);          /* half x half -> full */
    mpf_ui_sub(r, x, r);
    mpf_mul(t1, t1, r);          /* half x half -> half */
    mpf_div_2exp(t1, t1, 1);
    mpf_add(r, t1, t2);
}

/* r = y/x   WARNING: r cannot be the same as y. */
void
my_div(pi_state state, mpf_t t1, mpf_t t2, mpf_t r, mpf_t y, mpf_t x)
{
    unsigned long prec, bits, prec0;

    prec0 = mpf_get_prec(r);

    if (prec0 <= DOUBLE_PREC)
    {
        mpf_set_d(r, mpf_get_d(y) / mpf_get_d(x));
        return;
    }

    bits = 0;
    for (prec = prec0; prec > DOUBLE_PREC; )
    {
        int bit = prec & 1;
        prec = (prec + bit) / 2;
        bits = bits*2 + bit;
    }

    mpf_set_prec_raw(t1, DOUBLE_PREC);
    mpf_ui_div(t1, 1, x);

    while (prec < prec0)
    {
        prec *= 2;
        if (prec < prec0)
        {
            /* t1 = t1+t1*(1-x*t1); */
            mpf_set_prec_raw(t2, prec);
            mpf_mul(t2, x, t1);          /* full x half -> full */
            mpf_ui_sub(t2, 1, t2);
            mpf_set_prec_raw(t2, prec/2);
            mpf_mul(t2, t2, t1);         /* half x half -> half */
            mpf_set_prec_raw(t1, prec);
            mpf_add(t1, t1, t2);
        }
        else
        {
            prec = prec0;
            /* t2=y*t1, t1 = t2+t1*(y-x*t2); */
            mpf_set_prec_raw(t2, prec / 2);
            mpf_mul(t2, t1, y);          /* half x half -> half */
            mpf_mul(r, x, t2);           /* full x half -> full */
            mpf_sub(r, y, r);
            mpf_mul(t1, t1, r);          /* half x half -> half */
            mpf_add(r, t1, t2);
            break;
        }
        prec -= (bits & 1);
        bits /= 2;
    }
}

/*///////////////////////////////////////////////////////////////////////////*/



static __inline__ void
fac_reset(fac_t f)
{
    f[0].num_facs = 0;
}

static __inline__ void
fac_init_size(fac_t f, long int s)
{
    if (s < INIT_FACS)
        s = INIT_FACS;

    f[0].fac  = flint_malloc(s*sizeof(unsigned long)*2);
    f[0].pow  = f[0].fac + s;
    f[0].max_facs = s;

    fac_reset(f);
}

static __inline__ void
fac_init(fac_t f)
{
    fac_init_size(f, INIT_FACS);
}

static __inline__ void
fac_clear(fac_t f)
{
    flint_free(f[0].fac);
}

static __inline__ void
fac_resize(fac_t f, long int s)
{
    if (f[0].max_facs < s)
    {
        fac_clear(f);
        fac_init_size(f, s);
    }
}

/* f = base^pow */
static __inline__ void
fac_set_bp(pi_state state, fac_t f, unsigned long base, long int pow)
{
    long int i;
    assert(base<state->sieve_size);
    for (i=0, base/=2; base>0; i++, base = state->sieve[base].nxt)
    {
        f[0].fac[i] = state->sieve[base].fac;
        f[0].pow[i] = state->sieve[base].pow*pow;
    }
    f[0].num_facs = i;
    assert(i<=f[0].max_facs);
}

/* r = f*g */
static __inline__ void
fac_mul2(pi_state state, fac_t r, fac_t f, fac_t g)
{
    long int i, j, k;

    for (i=j=k=0; i<f[0].num_facs && j<g[0].num_facs; k++)
    {
        if (f[0].fac[i] == g[0].fac[j])
        {
            r[0].fac[k] = f[0].fac[i];
            r[0].pow[k] = f[0].pow[i] + g[0].pow[j];
            i++; j++;
        }
        else if (f[0].fac[i] < g[0].fac[j])
        {
            r[0].fac[k] = f[0].fac[i];
            r[0].pow[k] = f[0].pow[i];
            i++;
        }
        else
        {
            r[0].fac[k] = g[0].fac[j];
            r[0].pow[k] = g[0].pow[j];
            j++;
        }
    }

    for (; i<f[0].num_facs; i++, k++)
    {
        r[0].fac[k] = f[0].fac[i];
        r[0].pow[k] = f[0].pow[i];
    }

    for (; j<g[0].num_facs; j++, k++)
    {
        r[0].fac[k] = g[0].fac[j];
        r[0].pow[k] = g[0].pow[j];
    }

    r[0].num_facs = k;
    assert(k<=r[0].max_facs);
}

/* f *= g */
static __inline__ void
fac_mul(pi_state state, fac_t f, fac_t g)
{
    fac_t tmp;
    fac_resize(state->fmul, f[0].num_facs + g[0].num_facs);
    fac_mul2(state, state->fmul, f, g);
    tmp[0]  = f[0];
    f[0]    = state->fmul[0];
    state->fmul[0] = tmp[0];
}

/* f *= base^pow */
static __inline__ void
fac_mul_bp(pi_state state, fac_t f, unsigned long base, unsigned long pow)
{
    fac_set_bp(state, state->ftmp, base, pow);
    fac_mul(state, f, state->ftmp);
}

/* remove factors of power 0 */
static __inline__ void
fac_compact(fac_t f)
{
    long int i, j;
    for (i=0, j=0; i<f[0].num_facs; i++)
    {
        if (f[0].pow[i]>0)
        {
            if (j < i)
            {
                f[0].fac[j] = f[0].fac[i];
                f[0].pow[j] = f[0].pow[i];
            }
            j++;
        }
    }
    f[0].num_facs = j;
}

/* convert factorized form to number */
void
bs_mul(pi_state state, mpz_t r, long int a, long int b)
{
    long int i, j;
    if (b-a<=32)
    {
        mpz_set_ui(r, 1);
        for (i=a; i<b; i++)
            for (j=0; j<state->fmul[0].pow[i]; j++)
        mpz_mul_ui(r, r, state->fmul[0].fac[i]);
    }
    else
    {
        mpz_t r2;
        mpz_init(r2);
        bs_mul(state, r2, a, (a+b)/2);
        bs_mul(state, r, (a+b)/2, b);
        mpz_mul(r, r, r2);
        mpz_clear(r2);
    }
}


/* f /= gcd(f,g), g /= gcd(f,g) */
void
fac_remove_gcd(pi_state state, mpz_t p, fac_t fp, mpz_t g, fac_t fg)
{
    long int i, j, k, c;

    fac_resize(state->fmul, min(fp->num_facs, fg->num_facs));

    for (i=j=k=0; i < fp->num_facs && j < fg->num_facs; )
    {
        if (fp->fac[i] == fg->fac[j])
        {
            c = min(fp->pow[i], fg->pow[j]);
            fp->pow[i] -= c;
            fg->pow[j] -= c;
            state->fmul->fac[k] = fp->fac[i];
            state->fmul->pow[k] = c;
            i++; j++; k++;
        }
        else if (fp->fac[i] < fg->fac[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }

    state->fmul->num_facs = k;
    assert(k <= state->fmul->max_facs);

    if (state->fmul->num_facs)
    {
        bs_mul(state, state->gcd, 0, state->fmul->num_facs);

        mpz_divexact(p, p, state->gcd);
        mpz_divexact(g, g, state->gcd);

        fac_compact(fp);
        fac_compact(fg);
    }
}

/*///////////////////////////////////////////////////////////////////////////*/





/* binary splitting */
void
bs(pi_state state, unsigned long a, unsigned long b,
    unsigned gflag, long int level)
{
    unsigned long i, mid;

    if (b - a == 1)
    {
        /*
          g(b-1,b) = (6b-5)(2b-1)(6b-1)
          p(b-1,b) = b^3 * C^3 / 24
          q(b-1,b) = (-1)^b*g(b-1,b)*(A+Bb).
        */
        mpz_set_ui(p1, b);
        mpz_mul_ui(p1, p1, b);
        mpz_mul_ui(p1, p1, b);
        mpz_mul_ui(p1, p1, (C/24)*(C/24));
        mpz_mul_ui(p1, p1, C*24);

        mpz_set_ui(g1, 2*b-1);
        mpz_mul_ui(g1, g1, 6*b-1);
        mpz_mul_ui(g1, g1, 6*b-5);

        mpz_set_ui(q1, b);
        mpz_mul_ui(q1, q1, B);
        mpz_add_ui(q1, q1, A);
        mpz_mul   (q1, q1, g1);

        if (b%2)
            mpz_neg(q1, q1);

        i=b;
        while ((i&1)==0) i>>=1;

        fac_set_bp(state, fp1, i, 3);	/*  b^3 */
        fac_mul_bp(state, fp1, 3*5*23*29, 3);
        fp1[0].pow[0]--;

        fac_set_bp(state, fg1, 2*b-1, 1);	/* 2b-1 */
        fac_mul_bp(state, fg1, 6*b-1, 1);	/* 6b-1 */
        fac_mul_bp(state, fg1, 6*b-5, 1);	/* 6b-5 */
    }
    else
    {
        /*
          p(a,b) = p(a,m) * p(m,b)
          g(a,b) = g(a,m) * g(m,b)
          q(a,b) = q(a,m) * p(m,b) + q(m,b) * g(a,m)
        */
        mid = a+(b-a)*0.5224;     /* tuning parameter */
        bs(state, a, mid, 1, level+1);

        state->top++;
        bs(state, mid, b, gflag, level+1);
        state->top--;

        if (level>=4) {           /* tuning parameter */
            fac_remove_gcd(state, p2, fp2, g1, fg1);
        }

        mpz_mul(p1, p1, p2);
        mpz_mul(q1, q1, p2);
        mpz_mul(q2, q2, g1);
        mpz_add(q1, q1, q2);
        fac_mul(state, fp1, fp2);

        if (gflag)
        {
            mpz_mul(g1, g1, g2);
            fac_mul(state, fg1, fg2);
        }
    }
}

void
build_sieve(pi_state state, long int n, sieve_t *s)
{
    long int m, i, j, k;

    state->sieve_size = n;
    m = (long int)sqrt(n);
    memset(s, 0, sizeof(sieve_t)*n/2);

    s[1/2].fac = 1;
    s[1/2].pow = 1;

    for (i=3; i<=n; i+=2)
    {
        if (s[i/2].fac == 0)
        {
            s[i/2].fac = i;
            s[i/2].pow = 1;
            if (i <= m)
            {
                for (j=i*i, k=i/2; j<=n; j+=i+i, k++)
                {
                    if (s[j/2].fac==0)
                    {
                        s[j/2].fac = i;
                        if (s[k].fac == i)
                        {
                            s[j/2].pow = s[k].pow + 1;
                            s[j/2].nxt = s[k].nxt;
                        }
                        else
                        {
                            s[j/2].pow = 1;
                            s[j/2].nxt = k;
                        }
                    }
                }
            }
        }
    }
}

void
mpfr_pi_chudnovsky(mpfr_t res, mpfr_rnd_t rnd)
{
    mpf_t  pi, qi, t1, t2;
    mpfr_prec_t prec;
    long int i, depth=1, terms;
    pi_state state;

    prec = mpfr_get_prec(res) + 64;
    terms = prec / (BITS_PER_DIGIT * DIGITS_PER_ITER);
    while ((1L<<depth)<terms)
        depth++;
    depth++;

    state->top = 0;
    state->sieve_size = max(3*5*23*29+1, terms*6);
    state->sieve = (sieve_t *)flint_malloc(sizeof(sieve_t)*(state->sieve_size)/2);
    build_sieve(state, state->sieve_size, state->sieve);

    /* allocate stacks */
    state->pstack = flint_malloc(sizeof(mpz_t)*depth);
    state->qstack = flint_malloc(sizeof(mpz_t)*depth);
    state->gstack = flint_malloc(sizeof(mpz_t)*depth);
    state->fpstack = flint_malloc(sizeof(fac_t)*depth);
    state->fgstack = flint_malloc(sizeof(fac_t)*depth);
    for (i=0; i<depth; i++)
    {
        mpz_init(state->pstack[i]);
        mpz_init(state->qstack[i]);
        mpz_init(state->gstack[i]);
        fac_init(state->fpstack[i]);
        fac_init(state->fgstack[i]);
    }

    mpz_init(state->gcd);
    fac_init(state->ftmp);
    fac_init(state->fmul);

    /* begin binary splitting process */
    if (terms<=0)
    {
        mpz_set_ui(p2,1);
        mpz_set_ui(q2,0);
        mpz_set_ui(g2,1);
    }
    else
    {
        bs(state, 0,terms,0,0);
    }

    /* free some resources */
    flint_free(state->sieve);

    mpz_clear(state->gcd);
    fac_clear(state->ftmp);
    fac_clear(state->fmul);

    for (i=1; i<depth; i++)
    {
        mpz_clear(state->pstack[i]);
        mpz_clear(state->qstack[i]);
        mpz_clear(state->gstack[i]);
        fac_clear(state->fpstack[i]);
        fac_clear(state->fgstack[i]);
    }

    mpz_clear(state->gstack[0]);
    fac_clear(state->fpstack[0]);
    fac_clear(state->fgstack[0]);

    flint_free(state->gstack);
    flint_free(state->fpstack);
    flint_free(state->fgstack);

      /*
        p*(C/D)*sqrt(C)
        pi = -----------------
        (q+A*p)
      */

    mpz_addmul_ui(q1, p1, A);
    mpz_mul_ui(p1, p1, C/D);

    mpf_init2(pi, prec);
    mpf_set_z(pi, p1);
    mpz_clear(p1);

    mpf_init2(qi, prec);
    mpf_set_z(qi, q1);
    mpz_clear(q1);

    flint_free(state->pstack);
    flint_free(state->qstack);

    /* initialize temp float variables for sqrt & div */
    mpf_init2(t1, prec);
    mpf_init2(t2, prec);

    /* final step */
    my_div(state, t1, t2, qi, pi, qi);
    my_sqrt_ui(state, t1, t2, pi, C);
    mpf_mul(qi, qi, pi);

    mpfr_set_f(res, qi, rnd);

    /* free float resources */
    mpf_clear(pi);
    mpf_clear(qi);

    mpf_clear(t1);
    mpf_clear(t2);
}
