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
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"

TEST_FUNCTION_START(fmpz_mod_poly_print_read, state)
{
#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)
    int i, j, n, result;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;
    fmpz_mod_ctx_t ctx;

    n = 100 * flint_test_multiplier();

    fmpz_mod_ctx_init_ui(ctx, 2);

    /* Randomise n polynomials, write to and read from a pipe */
    {
        fmpz_mod_poly_t *a;

        a = flint_malloc(n * sizeof(fmpz_mod_poly_t));
        for (i = 0; i < n; i++)
        {
            fmpz_mod_poly_init(a[i], ctx);
            fmpz_mod_poly_randtest(a[i], state, n_randint(state, 100), ctx);
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
                r = fmpz_mod_poly_fprint(out, a[j], ctx);
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

            for (j = 0; j < n; j++)
                fmpz_mod_poly_clear(a[j], ctx);
            flint_free(a);
            fclose(out);
            exit(0);
        }
        else  /* Parent process */
        {
            int r;
            fmpz_mod_poly_t t;
            fmpz_mod_ctx_t newctx;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_ctx_init_ui(newctx, 3);
            fmpz_mod_poly_init(t, newctx);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_mod_poly_fread(in, t, newctx);
                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Read error.\n");
                    fflush(stdout);
                    flint_abort();
                }

                result = fmpz_equal(fmpz_mod_ctx_modulus(ctx),
                                    fmpz_mod_ctx_modulus(newctx));
                result = result && fmpz_mod_poly_equal(t, a[i], newctx);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("a[i] = "), fmpz_mod_poly_print(a[i], ctx), flint_printf("\n");
                    flint_printf("t    = "), fmpz_mod_poly_print(t, newctx), flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }

                ++i;
            }

            fmpz_mod_poly_clear(t, newctx);
            fmpz_mod_ctx_clear(newctx);
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
            fmpz_mod_poly_clear(a[i], ctx);
        flint_free(a);
    }

    /* Write bad data to a pipe and read it */
    {
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
            fmpz_mod_poly_t t;
            fmpz_mod_ctx_t newctx;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_poly_init(t, ctx);
            fmpz_mod_ctx_init_ui(newctx, 3);

            i = 0;
            /* Only four junk bytes are sent and our read
               doesn't consume invalid bytes, so eof is never reached */
            for(i = 0; i < 500; i++)
            {
                r = fmpz_mod_poly_fread(in, t, newctx);
                if (r > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("r = %d\n", r);
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_mod_poly_clear(t, newctx);
            fmpz_mod_ctx_clear(newctx);
            fclose(in);
            waitpid(childpid, NULL, 0);
        }
    }

    fmpz_mod_ctx_clear(ctx);

    TEST_FUNCTION_END(state);
#else
    TEST_FUNCTION_END_SKIPPED(state);
#endif
}
