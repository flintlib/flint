/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <sys/types.h>
#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER) 
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)

/*
    The function fdopen is declared in stdio.h.  It is POSIX.1 compliant, 
    but not ANSI compliant.  The following line enables compilation with 
    the "-ansi" flag.
 */
extern FILE * fdopen(int fildes, const char *mode);

int main(void)
{
    int i, j, n, result;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;
    fmpz_t two;

    FLINT_TEST_INIT(state);

    n = 100 * flint_test_multiplier();

    fmpz_init(two);
    fmpz_set_ui(two,2);

    flint_printf("print/ read....");
    fflush(stdout);

    /* Randomise n polynomials, write to and read from a pipe */
    {
        fmpz_mod_poly_t *a;

        a = flint_malloc(n * sizeof(fmpz_mod_poly_t));
        for (i = 0; i < n; i++)
        {
            fmpz_mod_poly_init(a[i], two);
            fmpz_mod_poly_randtest(a[i], state, n_randint(state, 100));
        }

        if (pipe(fd))
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to set-up the pipe.\n");
            abort();
        }

        if((childpid = fork()) == -1)
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to fork the process.\n");
            abort();
        }

        if(childpid == 0)  /* Child process */
        {
            int r;

            close(fd[0]);
            out = fdopen(fd[1], "w");
            if (out == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open output file at the pipe.\n");
                abort();
            }

            for (j = 0; j < n; j++)
            {
                r = fmpz_mod_poly_fprint(out, a[j]);
                if ((j < n - 1) && (r > 0))
                    r = flint_fprintf(out, "\n");

                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Write error.\n");
                    abort();
                }
            }

            for (j = 0; j < n; j++)
                fmpz_mod_poly_clear(a[j]);
            flint_free(a);
            fclose(out);
            exit(0);
        }
        else  /* Parent process */
        {
            int r;
            fmpz_mod_poly_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_mod_poly_init(t,two);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_mod_poly_fread(in, t);
                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Read error.\n");
                    abort();
                }

                result = fmpz_mod_poly_equal(t, a[i]);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("a[i] = "), fmpz_mod_poly_print(a[i]), flint_printf("\n");
                    flint_printf("t    = "), fmpz_mod_poly_print(t), flint_printf("\n");
                    abort();
                }

                ++i;
            }

            fmpz_mod_poly_clear(t);
            fclose(in);
        }

        if (i != n)
        {
            flint_printf("FAIL:\n");
            flint_printf("Only %d out of %d objects were processed.\n", i, n);
            abort();
        }

        for (i = 0; i < n; i++)
            fmpz_mod_poly_clear(a[i]);
        flint_free(a);
    }

    /* Write bad data to a pipe and read it */
    {
        if (pipe(fd))
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to set-up the pipe.\n");
            abort();
        }

        if((childpid = fork()) == -1)
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to fork the process.\n");
            abort();
        }

        if(childpid == 0)  /* Child process */
        {
            int r;

            close(fd[0]);
            out = fdopen(fd[1], "w");
            if (out == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open output file at the pipe.\n");
                abort();
            }

            r = flint_fprintf(out, "blah");
            if (r <= 0)
            {
                flint_printf("FAIL:\n");
                flint_printf("Write error.\n");
                abort();
            }

            fclose(out);
            exit(0);
        }
        else  /* Parent process */
        {
            int r;
            fmpz_mod_poly_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_mod_poly_init(t,two);

            i = 0;
            /* Only four junk bytes are sent and our read
               doesn't consume invalid bytes, so eof is never reached */
            for(i = 0; i < 500; i++)
            {
                r = fmpz_mod_poly_fread(in, t);
                if (r > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("r = %d\n", r);
                    abort();
                }
            }

            fmpz_mod_poly_clear(t);
            fclose(in);
        }
    }

    fmpz_clear(two);
    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

#else

int main(void)
{
    flint_printf("print/ read....");
    fflush(stdout);
    flint_printf("SKIPPED\n");
    return EXIT_SUCCESS;
}

#endif
