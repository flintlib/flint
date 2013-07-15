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

    Copyright (C) 2011 Andy Novocin

******************************************************************************/


#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"

#if !defined (__WIN32) || defined(__CYGWIN__)

/*
    The function fdopen is declared in stdio.h.  It is POSIX.1 compliant, 
    but not ANSI compliant.  The following line enables compilation with 
    the "-ansi" flag.
 */
extern FILE * fdopen(int fildes, const char *mode);

int main(void)
{
    int i, j, m, n, k = 1000, result;
    flint_rand_t state;

    FILE *in, *out;
    int fd[2];
    pid_t childpid;

    printf("print/ read....");
    fflush(stdout);

    flint_randinit(state);

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

            for (j = 0; j < k; j++)
            {
                r = fmpz_mat_fprint(out, M[j]);
                if ((j < k - 1) && (r > 0))
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
            fmpz_mat_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open input file at the pipe.\n");
                abort();
            }

            for (i = 0; i < k && !feof(in); i++)
            {
                fmpz_mat_init(t, 0, 0);

                r = fmpz_mat_fread(in, t);
                if (r <= 0)
                {
                    printf("FAIL:\n");
                    printf("Read error.\n");
                    abort();
                }

                result = fmpz_mat_equal(t, M[i]);
                if (!result)
                {
                    printf("FAIL:\n");
                    printf("M[i] = "), fmpz_mat_print(M[i]), printf("\n");
                    printf("t    = "), fmpz_mat_print(t), printf("\n");
                    abort();
                }

                fmpz_mat_clear(t);
            }

            fclose(in);
        }

        if (i != k)
        {
            printf("FAIL:\n");
            printf("Only %d out of %d objects were processed.\n", i, n);
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

            r = fprintf(out, "blah");
            if (r <= 0)
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
            int r;
            fmpz_mat_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_mat_init(t,0,0);

            i = 0;
            while (!feof(in))
            {
                r = fmpz_mat_fread(in, t);
                if (r > 0)
                {
                    printf("FAIL:\n");
                    printf("r = %d\n", r);
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
            printf("FAIL:\n");
            printf("Carried out %d reads, but \"%s\" has only 4 characters.\n", i, str);
            abort();
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
