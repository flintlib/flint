/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "flint.h"
#include "templates.h"

#include <sys/stat.h>
#include <math.h>
#include "profiler.h"
#include "fmpz_mat.h"

#define nalgs 2
#define cpumin 1
#define ncases 2

int
get_timings(double* s, slong degree, flint_bitcnt_t bits, slong length)
{
    TEMPLATE(T, ctx_t) ctx;
    TEMPLATE(T, poly_t) f, *h, finv;
    TEMPLATE(T, mat_t) HH;
    fmpz_t p, q;
    double beta;
    slong i, l;
    int n, c, reps = 0;
    FLINT_TEST_INIT(state);
    
    fmpz_init(p);
    fmpz_init(q);

    beta = 0.5 * (1. - (log(2) / log(length)));
    l = ceil(pow(length, beta));

    if (!(h = flint_malloc((l + 1) * sizeof(TEMPLATE(T, poly_struct)))))
    {
        flint_printf("Exception (p-iterated_frobenius):\n");
        flint_printf("Not enough memory.\n");
        flint_abort();
    }

    flint_printf("Trying %d %d %d\n", degree, bits, length);
    
    for (c = 0; c < nalgs; c++)
        s[c] = 0.0;

    reps = 0;
    /* Compute the timings */
    for (n = 0; n < ncases; n++)
    {
        double t[nalgs];
        int lo, loops = 1;

#ifdef FQ_ZECH_VEC_NORM
        do
        {
            fmpz_set_ui(p, n_randprime(state, bits, 1));
            fmpz_pow_ui(q, p, degree);
        } while (fmpz_cmp_ui(q, 1048576) > 0);
#else
        fmpz_set_ui(p, n_randprime(state, bits, 1));
        fmpz_pow_ui(q, p, degree);
#endif
        
        TEMPLATE(T, ctx_init)(ctx, p, degree, "a");
        TEMPLATE(T, poly_init)(f, ctx);
        TEMPLATE(T, poly_init)(finv, ctx);

        TEMPLATE(T, ctx_order)(q, ctx);

#ifdef FQ_ZECH_VEC_NORM
        if (fmpz_cmp_ui(q, 1048576) > 0)
        {
            flint_printf("Order too big for zech representation: ");
            fmpz_print(q);
            flint_printf("\n");
            flint_abort();            
        }
#endif        
        
        for (i = 0; i < l + 1; i++)
            TEMPLATE(T, poly_init)(h[i], ctx);

        /*
          Construct random elements of fq
        */
        {
            TEMPLATE(T, poly_randtest_monic)(f, state, length, ctx);
            TEMPLATE(T, poly_reverse)(finv, f, f->length, ctx);
            TEMPLATE(T, poly_inv_series_newton)(finv, finv, f->length, ctx);
        }
                
    loop:

        t[0] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gen)(h[0], ctx);
            TEMPLATE(T, mat_init)(HH, n_sqrt(f->length - 1) + 1, f->length - 1, ctx);
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[1], h[0], q, 0, f, finv, ctx);
            TEMPLATE(T, poly_precompute_matrix)(HH, h[1], f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                TEMPLATE(T, poly_compose_mod_brent_kung_precomp_preinv)(h[i], h[i - 1],
                                                                        HH, f, finv, ctx);
            TEMPLATE(T, mat_clear)(HH, ctx);
        }
        prof_stop();
        t[0] += get_clock(0);

        
        t[1] = 0.0;
        init_clock(0);
        prof_start();
        for (lo = 0; lo < loops; lo++)
        {
            TEMPLATE(T, poly_gen)(h[0], ctx);
            TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[1], h[0], q, 0, f, finv, ctx);
            for (i = 2; i < l + 1; i++)
                TEMPLATE(T, poly_powmod_fmpz_sliding_preinv)(h[i], h[i-1], q, 0, f, finv, ctx);
        }
        prof_stop();
        t[1] += get_clock(0);

        for (c = 0; c < nalgs; c++)
            if (t[c] * FLINT_CLOCK_SCALE_FACTOR <= cpumin)
            {
                loops *= 2;
                goto loop;
            }
                
        for (c = 0; c < nalgs; c++)
            s[c] += t[c];
        reps += loops;

        TEMPLATE(T, poly_clear)(f, ctx);
        TEMPLATE(T, poly_clear)(finv, ctx);
        for (i = 0; i < l + 1; i++)
            TEMPLATE(T, poly_clear)(h[i], ctx);
        TEMPLATE(T, ctx_clear)(ctx);

    }

    for (c = 0; c < nalgs; c++)
    {
        s[c] = s[c] / (double) reps;
    }

    fmpz_clear(p);
    fmpz_clear(q);


    flint_free(h);

    FLINT_TEST_CLEANUP(state);
    
    return s[0] > s[1];
}

long
a(fmpz_mat_t array, slong i, slong j)
{
    return fmpz_get_si(fmpz_mat_entry(array, i, j));
}

int
file_exists(char *filename)
{
  struct stat buffer;   
  return (stat (filename, &buffer) == 0);
}

int
init_array(fmpz_mat_t array, slong max_degree, slong max_bits, slong max_length, char* filename)
{
    int bigger_length = 0;
    fmpz_mat_t old_array;
    slong i, j;
    FILE * old_file;
    
    if( file_exists(filename) )
    {
        flint_printf(filename);
        /* old file exists */
        fmpz_mat_init(old_array, max_degree, max_bits);
        old_file = fopen(filename, "r");
        fmpz_mat_fread(old_file, old_array);
        fclose(old_file);
        if (fmpz_get_ui(fmpz_mat_entry(old_array, 0, 2)) < max_length)
            bigger_length = 1;
    }
    fmpz_mat_init(array, max_degree, max_bits);
    max_bits = FLINT_MAX(max_bits, 3);
    fmpz_set_si(fmpz_mat_entry(array, 0, 0), max_degree);
    fmpz_set_si(fmpz_mat_entry(array, 0, 1), max_bits);
    fmpz_set_si(fmpz_mat_entry(array, 0, 2), max_length);

    if( file_exists(filename) )
    {
        for (i = 2; i < max_degree; i++)
        {
            for (j = 0; j < max_bits; j++)
            {
                fmpz_set(fmpz_mat_entry(array, i, j),
                         fmpz_mat_entry(old_array, i, j));
                         
            }
        }
        fmpz_mat_clear(old_array);
    }
    fmpz_mat_print_pretty(array);
    return bigger_length;
}


void write_array(fmpz_mat_t array, char * filename)
{
    FILE * tmp;
    tmp = fopen(filename, "w");
    fmpz_mat_fprint(tmp, array);
    fclose(tmp);
}



int
main(int argc, char** argv)
{
    flint_bitcnt_t bits, max_bits, max_bits_used, max_bits_e;
    int is_hit, bigger_length;
    slong degree, length, max_degree, max_length, imin, imax, imid, diff;
    fmpz_mat_t array;
    char* filename;

    double s[nalgs];

    max_degree = atol(argv[1]);
    max_bits = atol(argv[2]);
    max_length = atol(argv[3]);
    filename = argv[4];

    bigger_length = init_array(array, max_degree, max_bits, max_length, filename);

    max_bits_used = 0;
    max_bits_e = 0;
    for (degree = 2; degree < max_degree; degree++)
    {
        flint_printf("Degree %d\n", degree);
        fflush(stdout);
        bits = 2;
        length = 3;
        while (bits < max_bits && (max_bits_e == 0 || bits < max_bits_e) )
        {
            if (a(array, degree, bits) != 0 ||
                (!bigger_length && bits >= a(array, degree, 0)))
            {
                bits += 1;
                continue;
            }

#ifdef FQ_ZECH_VEC_NORM
            /* Don't make zech fields too big */
            if (degree * bits >= 20)
            {
                bits += 1;
                continue;
            }
#endif
            
            /* Set the initial state */
            if (bits == 2 || bits == 3)
            {
                if (degree == 2)
                {
                    length = 3;
                }
                else
                {
                    if (a(array, degree - 1, bits) > length)
                        length = a(array, degree - 1, bits);
                }
                diff = length;
            }
            else
            {
                length = a(array, degree, bits - 1);
                diff = length - a(array, degree, bits - 2);
                if (diff < degree)
                {
                    diff = degree;
                }
            }

            /* Set the min */
            imax = 0;
            imin = length;
            is_hit = get_timings(s, degree, bits, imin);
            while (is_hit != 0)
            {
                imax = imin;
                imin -= 1;
                if (imin < 3)
                {
                    break;
                }
                is_hit = get_timings(s, degree, bits, imin);
            }

            if (imin < 3)
            {
                bits += 1;
                continue;
            }

            /* Set the max */
            if (imax == 0)
            {
                imax = FLINT_MIN(imin + 2 * diff, max_length);
                is_hit = get_timings(s, degree, bits, imax);
                while (is_hit != 1)
                {
                    flint_printf("Finding max, %d\n", diff);
                    if (imax == max_length)
                    {
                        imax = max_length + 1;
                        break;
                    }
                    imin = imax;
                    imax += 2 * diff;
                    imax = FLINT_MIN(imax, max_length);
                    is_hit = get_timings(s, degree, bits, imax);
                }
            }
            if (imax > max_length)
            {
                max_bits_e = bits;
                fmpz_set_si(fmpz_mat_entry(array, degree, 0), bits);
                write_array(array, filename);
                break;
            }

            flint_printf("Min - Max: %d - %d\n", imin, imax);
            
            while (imin < imax)
            {
                imid = imin + ((imax - imin) / 2);
                if (imid >= imax)
                {
                    flint_printf("Error in computing midpoint\n");
                    flint_abort();
                }
                
                is_hit = get_timings(s, degree, bits, imid);

                if (is_hit)
                {
                    imax = imid;
                }
                else
                {
                    imin = imid + 1;
                }
            }

            length = imin;

            /* Set the array */
            flint_printf("%d - %d %d\n", degree, bits, length);
            fmpz_set_si(fmpz_mat_entry(array, degree, bits), length);
            fmpz_set_si(fmpz_mat_entry(array, degree, 0), bits);
            if (degree == 2 && bits > max_bits_used)
                max_bits_used = bits + 1;
            fflush(stdout);
            write_array(array, filename);

            bits += 1;
        }


    }


    fmpz_mat_print_pretty(array);
    fmpz_mat_clear(array);
    return 0;
}

#endif
