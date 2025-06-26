#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fft_small.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "dirichlet.h"
#include "profiler.h"

/* precomputed smoothed table */
typedef struct {
    slong pe;
    slong m;
} pem_struct;
typedef pem_struct * pem_ptr;

/* factor table */
typedef struct {
    slong a;
    slong b;
    slong prev;
    slong next;
} coprime_struct;
typedef coprime_struct * coprime_ptr;
typedef const coprime_ptr coprime_srcptr;

struct rough {
    ulong m;
    struct rough * prev;
    struct rough * next;
};
//typedef struct {
//    ulong m;
//    ulong prev;
//    ulong next;
//} rough_struct;
typedef struct rough * rough_ptr;

void
pem_init(pem_ptr tab, slong len)
{
    ulong k;
    for (k = 0; k < len; k++)
        tab[k].pe = 1, tab[k].m = k;
}

/* assume tab has been initialized */
void
pem_init_rough_lim(pem_ptr tab, slong lim, slong len)
{
    ulong m, pe, pem, len1;
    rough_ptr rough, p1, m1;

    /* p=2 can be done separately */
    for (m = 3; m < len; m += 2)
        for (pe = 2, pem = 2*m; pem < len; pem <<= 1, pe <<= 1)
            tab[pem].pe = pe, tab[pem].m = m;

    /* now need 2-rough (ie odd) numbers up to len / 3 */
    len1 = len / 3;
    rough = flint_malloc((1 + len1/2) * sizeof(struct rough));
    rough->m = 1;
    rough->prev = NULL;
    rough->next = rough + 1;
    for (m1 = rough + 1, m = 3; m < len1; m1++, m += 2)
    {
        m1->m = m;
        m1->prev = m1 - 1;
        m1->next = m1 + 1;
    }
    /* terminate */
    m1->m = len;
    m1->prev = m1 - 1;
    m1->next = NULL;

    for (p1 = rough + 1; p1->m < lim; p1 = p1->next)
    {
        slong p = p1->m, pe;
        /* skip prime powers p^e up to len1 */
        for (pe = p; pe < len1; pe *= p)
        {
            rough_ptr pe1 = rough + (pe>>1);
            pe1->next->prev = pe1->prev;
            pe1->prev->next = pe1->next;
        }
        /* then loop on p-rough numbers */
        for (pe = p; pe < len; pe *= p)
        {
            ulong lim = len / pe;
            /* loop m in p-rough numbers */
            for (m1 = p1->next; m1->m < lim; m1 = m1->next)
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

void
pem_init_rough(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_rough_lim(tab, n_sqrt(len), len);
}

void
pem_init_rough_all(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_rough_lim(tab, len-1, len);
}
/* can ignore even numbers */
/* identify pmax-pem numbers in complete table of numbers less than len */
void
pem_init_sieve_range(pem_ptr tab, slong pmin, slong pmax, slong len)
{
    slong p, m;
    n_primes_t iter;
    n_primes_init(iter);
    if (pmin > 2) n_primes_jump_after(iter, pmin-1);
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                    if (tab[pem].pe == 1) tab[pem].m = m, tab[pem].pe = pe;
    }
    n_primes_clear(iter);
}
void
pem_init_sieve_all(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_sieve_range(tab, 2, len, len);
}
void
pem_init_sieve_sqrt(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_sieve_range(tab, 2, n_sqrt(len), len);
}
/* visit p^e m with m p-smooth */
void
pem_init_smooth_lim(pem_ptr tab, slong pmax, slong len)
{
    slong p, m;
    char * smooth;
    n_primes_t iter;
    n_primes_init(iter);

    smooth = flint_malloc(len * sizeof(char));
    memset(smooth, 0, len);

    smooth[1] = 1;
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                    if (smooth[m])
                    {
                        tab[pem].m = m, tab[pem].pe = pe;
                        smooth[pem] = 1;
                    }
    }
    n_primes_clear(iter);
    flint_free(smooth);
}
/* sieve up to sqrt(n) */
void
pem_init_smooth_sqrt(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_smooth_lim(tab, n_sqrt(len), len);
}
/* smooth then usual primes */
void
pem_init_smooth_all(pem_ptr tab, slong len)
{
    n_primes_t iter;
    slong p, pmax = n_sqrt(len);
    pem_init(tab, len);
    pem_init_smooth_lim(tab, pmax, len);
    n_primes_init(iter);
    n_primes_jump_after(iter, pmax);
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        ulong m, pm;
        for (m = 1, pm = p; pm < len; m++, pm += p)
        {
            tab[pm].pe = p;
            tab[pm].m = m;
        }
    }
    n_primes_clear(iter);
}

/* try to revert loop */
void
pem_init_sieve_smooth2(pem_ptr tab, slong len)
{
    slong p, m, pmax = n_sqrt(len);
    char * smooth;
    n_primes_t iter;

    smooth = flint_malloc(len * sizeof(char));
    memset(smooth, 0, len);

    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;
    tab[1].m = tab[1].pe = 1;
    smooth[1] = 1;

    n_primes_init(iter);
    p = n_primes_next(iter);
    for (ulong pe = 2; pe < len; pe *= 2)
        tab[pe].pe = pe, tab[pe].m = 1, smooth[pe] = 1;
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        ulong m, pm;
        /* We must take care: m must be strictly p-smooth in
         * this loop but we declare p-smooth numbers at the
         * same time.
         * This is not a big issue when we loop on p^e then on
         * m, since the last decomposition written is the good
         * one (p^e maximal).
         * When reversing the loops we must loop backwards on m
         * to solve the problem. */
        ulong e, pe[32], lim[32];
        pe[0] = p;
        /* precompute valid exponents */
        for (e = 0; pe[e] < len; e++)
        {
            lim[e] = len / pe[e];
            pe[e+1] = pe[e] * p;
        }
        lim[e] = 0;
        /* loop */
        for (m = lim[0], pm = p*m; m; m--, pm-=p)
        {
            for (; pm && !smooth[m]; m--, pm-=p)
            /* by definition of smooth m is not divisible by p */
            smooth[pm] = 1;
            tab[pm].m = m;
            tab[pm].pe = p;
            for (e = 1; m < lim[e]; e++)
            {
               ulong pem = pe[e]*m;
               tab[pem].m = m;
               tab[pem].pe = pe[e];
            }
        }
        /* set smooth p^e last */
        for (e = 0; lim[e]; e++)
        {
            smooth[pe[e]] = 1;
        }
    }
    n_primes_clear(iter);
    flint_free(smooth);
}

typedef struct {
    char * name;
    void (*func)(pem_ptr, slong);
} smooth_func;

int main(int argc, char * argv[])
{
    slong i;
    int e, opt_min = 15, opt_max = 28;
    slong opt_print = 0;

#define NUM 6
    const smooth_func func[NUM] = {
        (const smooth_func){ "all", &pem_init_sieve_all },
        (const smooth_func){ "sqrt", &pem_init_sieve_sqrt },
        (const smooth_func){ "smoothtab", &pem_init_smooth_sqrt },
        (const smooth_func){ "tab all", &pem_init_smooth_all },
        (const smooth_func){ "rough", &pem_init_rough },
        (const smooth_func){ "rough all", &pem_init_rough_all }
    };

    /* options */
    for (i = 1; i < argc;)
    {
        if (strcmp(argv[i], "--min") == 0 && i + 1 < argc)
            opt_min = atol(argv[i+1]), i +=2 ;
        else if (strcmp(argv[i], "--max") == 0 && i + 1 < argc)
            opt_max = atol(argv[i+1]), i +=2 ;
        else if (strcmp(argv[i], "--print") == 0 && i + 1 < argc)
            opt_print = atol(argv[i+1]), i +=2 ;
        else break;
    }

    if (opt_print)
    {
        slong i, j, len;
        pem_ptr tab[NUM];
        len = opt_print;

        for (i = 0; i < NUM; i++)
        {
            tab[i] = flint_malloc(len * sizeof(pem_struct));
            (func[i].func)(tab[i], len);
        }
        flint_printf("index");
        for (i = 0; i < NUM; i++)
            flint_printf("  # %9s  ", func[i].name);
        flint_printf("\n");
        for (j = 1; j < len; j++)
        {
            flint_printf("%5ld", j);
            for (i = 0; i < NUM; i++)
            {
                pem_struct pem = tab[i][j];
                if (pem.pe > 1)
                    flint_printf("  #  %3ld * %4ld", pem.pe, pem.m);
                else 
                    flint_printf("  #    . *    .");
            }
            flint_printf("\n");
        }
        for (i = 0; i < NUM; i++)
            flint_free(tab[i]);

        return 0;
    }

    flint_printf("2^k# none");
    for (i = 0; i < NUM; i++)
        flint_printf(" #  %8s ", func[i].name);
    flint_printf("\n");
    for (e = opt_min; e <= opt_max; e++)
    {
        timeit_t t0, t1;
        slong i, len = 1UL<<e;
        double ref;
        pem_ptr tab;

        timeit_start(t0);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_init(tab, len);
        flint_free(tab);
        timeit_stop(t0);
        flint_printf("%2d # %4wd", e, t0->wall);
        ref = 1. / (1 + t0->wall);
        for (i = 0; i < NUM; i++)
        {
            const smooth_func f = func[i];
            timeit_start(t1);
            tab = flint_malloc(len * sizeof(pem_struct));
            (f.func)(tab, len);
            flint_free(tab);
            timeit_stop(t1);
            flint_printf(" # %4wd [%02.1f]", t1->wall, t1->wall * ref);
        }
        flint_printf("\n");
    }
}
