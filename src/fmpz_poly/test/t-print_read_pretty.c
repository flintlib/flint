/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* try to get fdopen declared */
#if defined __STRICT_ANSI__
# undef __STRICT_ANSI__
#endif

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)
# include <stdlib.h>
# include <sys/wait.h>
# include <unistd.h>
#endif
#include <string.h>

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)

TEST_FUNCTION_START(fmpz_poly_print_read_pretty, state)
{
    int i, j, n = 1000, result;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;

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
            flint_printf("FAIL:\n");
            flint_printf("Failed to set-up the pipe.\n");
            fflush(stdout);
            flint_abort();
        }

        if((childpid = fork()) == -1)
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to fork the process.\n");
            fflush(stdout);
            flint_abort();
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
                fflush(stdout);
                flint_abort();
            }

            for (j = 0; j < n; j++)
            {
                r = fmpz_poly_fprint_pretty(out, a[j], var);
                if ((j < n - 1) && (r > 0))
                    r = flint_fprintf(out, "\n");

                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Write error.\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            fclose(out);
            for (i = 0; i < n; ++i)
                fmpz_poly_clear(a[i]);
            flint_free(a);
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
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_init(t);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_poly_fread_pretty(in, t, &rvar);
                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Read error.\n");
                    fflush(stdout);
                    flint_abort();
                }

                result = fmpz_poly_equal(t, a[i]) &&
                    (t->length <= 1 || (strcmp(var, rvar) == 0));
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("a[i] = "), fmpz_poly_print_pretty(a[i], var), flint_printf("\n");
                    flint_printf("t    = "), fmpz_poly_print_pretty(t, rvar), flint_printf("\n");
                    flint_printf("rvar = %s\n", rvar);
                    fflush(stdout);
                    flint_abort();
                }
                flint_free(rvar);

                ++i;
            }

            fmpz_poly_clear(t);
            fclose(in);
            waitpid(childpid, NULL, 0);
        }

        if (i != n)
        {
            flint_printf("FAIL:\n");
            flint_printf("Only %d out of %d objects were processed.\n", i, n);
            fflush(stdout);
            flint_abort();
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
            flint_printf("FAIL:\n");
            flint_printf("Failed to set-up the pipe.\n");
            fflush(stdout);
            flint_abort();
        }

        if((childpid = fork()) == -1)
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to fork the process.\n");
            fflush(stdout);
            flint_abort();
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
                fflush(stdout);
                flint_abort();
            }

            r = fputs(str, out);
            if (r == EOF)
            {
                flint_printf("FAIL:\n");
                flint_printf("Write error.\n");
                fflush(stdout);
                flint_abort();
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
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_poly_init(t);

            while (!feof(in))
            {
                r = fmpz_poly_fread_pretty(in, t, &rvar);
                result = (r > 0) && rvar && (strcmp(str, rvar) == 0) &&
                         (t->length == 2) && (t->coeffs[0] == WORD(0)) &&
                         (t->coeffs[1] == WORD(1));
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("r    = %d\n", r);
                    flint_printf("str  = {%s}\n", str);
                    flint_printf("rvar = {%s}\n", rvar);
                    flint_printf("t    = "), fmpz_poly_print(t), flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
                if (rvar)
                    flint_free(rvar);
            }

            fmpz_poly_clear(t);
            fclose(in);
            waitpid(childpid, NULL, 0);
        }
    }

    TEST_FUNCTION_END(state);
}

#else

TEST_FUNCTION_START(fmpz_poly_print_read_pretty, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}

#endif
