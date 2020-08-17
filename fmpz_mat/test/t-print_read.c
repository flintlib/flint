/*
    Copyright (C) 2011 Andy Novocin

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
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
#include "fmpz_mat.h"

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)

/*
    The function fdopen is declared in stdio.h.  It is POSIX.1 compliant, 
    but not ANSI compliant.  The following line enables compilation with 
    the "-ansi" flag.
 */
extern FILE * fdopen(int fildes, const char *mode);

int main(void)
{
    int i, j, m, n, k = 1000, result;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;

    FLINT_TEST_INIT(state);

    flint_printf("print/ read....");
    fflush(stdout);    

    /* Randomise k mats, write to and read from a pipe */
    {
        fmpz_mat_t *M;

        M = flint_malloc(k * sizeof(fmpz_mat_t));
        for (i = 0; i < k; i++)
        {
            m = n_randint(state, 10);
            n = n_randint(state, 10);
            fmpz_mat_init(M[i], m, n);
            fmpz_mat_randtest(M[i], state, 100);
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

            for (j = 0; j < k; j++)
            {
                r = fmpz_mat_fprint(out, M[j]);
                if ((j < k - 1) && (r > 0))
                    r = flint_fprintf(out, "\n");

                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Write error.\n");
                    abort();
                }
            }

            fclose(out);
            for (i = 0; i < k; ++i)
                fmpz_mat_clear(M[i]);
            flint_free(M);
            exit(0);
        }
        else  /* Parent process */
        {
            int r;
            fmpz_mat_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                abort();
            }

            for (i = 0; i < k && !feof(in); i++)
            {
                fmpz_mat_init(t, 0, 0);

                r = fmpz_mat_fread(in, t);
                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Read error.\n");
                    abort();
                }

                result = fmpz_mat_equal(t, M[i]);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("M[i] = "), fmpz_mat_print(M[i]), flint_printf("\n");
                    flint_printf("t    = "), fmpz_mat_print(t), flint_printf("\n");
                    abort();
                }

                fmpz_mat_clear(t);
            }

            fclose(in);
        }

        if (i != k)
        {
            flint_printf("FAIL:\n");
            flint_printf("Only %d out of %d objects were processed.\n", i, n);
            abort();
        }

        for (i = 0; i < k; i++)
            fmpz_mat_clear(M[i]);
        flint_free(M);
    }

    /* Write bad data to a pipe and read it */
    {
        char str[5] = {'b', 'l', 'a', 'h', '\0'};

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
            fmpz_mat_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_mat_init(t,0,0);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_mat_fread(in, t);
                if (r > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("r = %d\n", r);
                    abort();
                }
                ++i;
            }

            fmpz_mat_clear(t);
            fclose(in);
        }

        /* For {'b','l','a','h','\0'} we expect 5 reads */
        if (i != 5)
        {
            flint_printf("FAIL:\n");
            flint_printf("Carried out %d reads, but \"%s\" has only 4 characters.\n", i, str);
            abort();
        }
    }

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
