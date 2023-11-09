/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2023 Albin Ahlbäck

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

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)

TEST_FUNCTION_START(fmpz_poly_print_read, state)
{
    int i, j, n = 1000, result;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;

    /* Randomise n polynomials, write to and read from a pipe */
    {
        fmpz_poly_t *a;

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

        if ((childpid = fork()) == -1)
        {
            flint_printf("FAIL:\n");
            flint_printf("Failed to fork the process.\n");
            fflush(stdout);
            flint_abort();
        }

        if (childpid == 0)  /* Child process */
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
                r = fmpz_poly_fprint(out, a[j]);
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
                r = fmpz_poly_fread(in, t);
                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Read error.\n");
                    fflush(stdout);
                    flint_abort();
                }

                result = fmpz_poly_equal(t, a[i]);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("a[i] = "), fmpz_poly_print(a[i]), flint_printf("\n");
                    flint_printf("t    = "), fmpz_poly_print(t), flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }

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

    /* Write bad data to a pipe and read it */
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

            r = flint_fprintf(out, "blah");
            if (r <= 0)
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

            i = 0;
            while (!feof(in))
            {
                r = fmpz_poly_fread(in, t);
                if (r > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("r = %d\n", r);
                    fflush(stdout);
                    flint_abort();
                }
                ++i;
            }

            fmpz_poly_clear(t);
            fclose(in);
            waitpid(childpid, NULL, 0);
        }

        /* For {'b','l','a','h','\0'} we expect 5 reads */
        if (i != 5)
        {
            flint_printf("FAIL:\n");
            flint_printf("Carried out %d reads, but \"%s\" has only 4 characters.\n", i, str);
            fflush(stdout);
            flint_abort();
        }
    }

    TEST_FUNCTION_END(state);
}

#else

TEST_FUNCTION_START(fmpz_poly_print_read, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}

#endif
