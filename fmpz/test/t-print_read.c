/*
    Copyright (C) 2010 Sebastian Pancratz

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

#if (!defined (__WIN32) || defined(__CYGWIN__)) && !defined(_MSC_VER)

/*
    The function fdopen is declared in stdio.h.  It is POSIX.1 compliant, 
    but not ANSI compliant.  The following line enables compilation with 
    the "-ansi" flag.
 */
extern FILE * fdopen(int fildes, const char *mode);

int main(void)
{
    int i, j, n = 10000, result;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;

    FLINT_TEST_INIT(state);

    flint_printf("print/ read....");
    fflush(stdout);   

    /* Randomise n integers, write to and read from a pipe */
    {
        fmpz *a;

        a = flint_calloc(n, sizeof(fmpz));
        for (i = 0; i < n; i++)
            fmpz_randtest(a + i, state, 200);

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
                r = fmpz_fprint(out, a + j);
                if ((j < n - 1) && (r > 0))
                    r = flint_fprintf(out, "\n");

                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Write error.\n");
                    abort();
                }
            }

            fclose(out);
            for (i = 0; i < n; ++i)
                fmpz_clear(a + i);
            flint_free(a);
            exit(0);
        }
        else  /* Parent process */
        {
            int r;
            fmpz_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_init(t);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_fread(in, t);
                if (r <= 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("Read error.\n");
                    abort();
                }

                result = fmpz_equal(t, a + i);
                if (!result)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("a[i] = "), fmpz_print(a + i), flint_printf("\n");
                    flint_printf("t    = "), fmpz_print(t), flint_printf("\n");
                    abort();
                }

                ++i;
            }

            fmpz_clear(t);
            fclose(in);
        }

        if (i != n)
        {
            flint_printf("FAIL:\n");
            flint_printf("Only %d out of %d objects were processed.\n", i, n);
            abort();
        }

        for (i = 0; i < n; i++)
            fmpz_clear(a + i);
        flint_free(a);
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
            fmpz_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                flint_printf("FAIL:\n");
                flint_printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_init(t);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_fread(in, t);
                if (r > 0)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("r = %d\n", r);
                    abort();
                }
                ++i;
            }

            fmpz_clear(t);
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
