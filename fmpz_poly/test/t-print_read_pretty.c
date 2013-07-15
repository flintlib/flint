/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#if !defined (__WIN32) || defined(__CYGWIN__)

/*
    The function fdopen is declared in stdio.h.  It is POSIX.1 compliant, 
    but not ANSI compliant.  The following line enables compilation with 
    the "-ansi" flag.
 */
extern FILE * fdopen(int fildes, const char *mode);

int main(void)
{
    int i, j, n = 1000, result;
    flint_rand_t state;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;

    printf("print/ read_pretty....");
    fflush(stdout);

    flint_randinit(state);

    /* Randomise n polynomials, write to and read from a pipe */
    {
        fmpz_poly_t *a;
        char *var = "x";

        a = flint_malloc(n * sizeof(fmpz_poly_t));
        for (i = 0; i < n; i++)
        {
            fmpz_poly_init(a[i]);
            fmpz_poly_randtest(a[i], state, 100, 100);
        }

        if (pipe(fd))
        {
            printf("FAIL:\n");
            printf("Failed to set-up the pipe.\n");
            abort();
        }

        if((childpid = fork()) == -1)
        {
            printf("FAIL:\n");
            printf("Failed to fork the process.\n");
            abort();
        }

        if(childpid == 0)  /* Child process */
        {
            int r;

            close(fd[0]);
            out = fdopen(fd[1], "w");
            if (out == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open output file at the pipe.\n");
                abort();
            }

            for (j = 0; j < n; j++)
            {
                r = fmpz_poly_fprint_pretty(out, a[j], var);
                if ((j < n - 1) && (r > 0))
                    r = fprintf(out, "\n");

                if (r <= 0)
                {
                    printf("FAIL:\n");
                    printf("Write error.\n");
                    abort();
                }
            }

            fclose(out);
            exit(0);
        }
        else  /* Parent process */
        {
            int r;
            fmpz_poly_t t;
            char *rvar;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_poly_init(t);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_poly_fread_pretty(in, t, &rvar);
                if (r <= 0)
                {
                    printf("FAIL:\n");
                    printf("Read error.\n");
                    abort();
                }

                result = fmpz_poly_equal(t, a[i]) &&
                    (t->length <= 1 || (strcmp(var, rvar) == 0));
                if (!result)
                {
                    printf("FAIL:\n");
                    printf("a[i] = "), fmpz_poly_print_pretty(a[i], var), printf("\n");
                    printf("t    = "), fmpz_poly_print_pretty(t, rvar), printf("\n");
                    printf("rvar = %s\n", rvar);
                    abort();
                }
                flint_free(rvar);

                ++i;
            }

            fmpz_poly_clear(t);
            fclose(in);
        }

        if (i != n)
        {
            printf("FAIL:\n");
            printf("Only %d out of %d objects were processed.\n", i, n);
            abort();
        }

        for (i = 0; i < n; i++)
            fmpz_poly_clear(a[i]);
        flint_free(a);
    }

    /* Write "blah" to the pipe and see it read as a variable */
    {
        char str[5] = {'b', 'l', 'a', 'h', '\0'};

        if (pipe(fd))
        {
            printf("FAIL:\n");
            printf("Failed to set-up the pipe.\n");
            abort();
        }

        if((childpid = fork()) == -1)
        {
            printf("FAIL:\n");
            printf("Failed to fork the process.\n");
            abort();
        }

        if(childpid == 0)  /* Child process */
        {
            int r;

            close(fd[0]);
            out = fdopen(fd[1], "w");
            if (out == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open output file at the pipe.\n");
                abort();
            }

            r = fputs(str, out);
            if (r == EOF)
            {
                printf("FAIL:\n");
                printf("Write error.\n");
                abort();
            }

            fclose(out);
            exit(0);
        }
        else  /* Parent process */
        {
            char *rvar = NULL;
            int r;
            fmpz_poly_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_poly_init(t);

            while (!feof(in))
            {
                r = fmpz_poly_fread_pretty(in, t, &rvar);
                result = (r > 0) && rvar && (strcmp(str, rvar) == 0) && 
                         (t->length == 2) && (t->coeffs[0] == 0L) && 
                         (t->coeffs[1] == 1L);
                if (!result)
                {
                    printf("FAIL:\n");
                    printf("r    = %d\n", r);
                    printf("str  = {%s}\n", str);
                    printf("rvar = {%s}\n", rvar);
                    printf("t    = "), fmpz_poly_print(t), printf("\n");
                    abort();
                }
                if (rvar)
                    flint_free(rvar);
            }

            fmpz_poly_clear(t);
            fclose(in);
        }
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

#else

int main(void)
{
    printf("print/ read....");
    fflush(stdout);
    printf("SKIPPED\n");
    return EXIT_SUCCESS;
}

#endif
