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
    Copyright (C) 2013 Qingwen GUAN

******************************************************************************/

#undef ulong
#define ulong ulongxx /* interferes with standard libraries */
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"

#if !defined (__WIN32) || defined(__CYGWIN__)

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

    printf("out_raw/inp_raw....");
    fflush(stdout);

    /* Randomise n integers, write to and read from a pipe */
    {
        fmpz *a;

        a = flint_calloc(n, sizeof(fmpz));
        for (i = 0; i < n; i++)
            fmpz_randtest(a + i, state, 200);

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
            size_t r;

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
                r = fmpz_out_raw(out, a + j);

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
            size_t r;
            fmpz_t t;

            close(fd[1]);
            in = fdopen(fd[0], "r");
            if (in == NULL)
            {
                printf("FAIL:\n");
                printf("Could not open input file at the pipe.\n");
                abort();
            }

            fmpz_init(t);

            i = 0;
            while ( (r = fmpz_inp_raw(t, in)) != 0 )
            {
                result = fmpz_equal(t, a + i);
                if (!result)
                {
                    printf("FAIL:\n");
                    printf("a[i] = "), fmpz_print(a + i), printf("\n");
                    printf("t    = "), fmpz_print(t), printf("\n");
                    abort();
                }

                ++i;
            }

            fmpz_clear(t);
            fclose(in);
        }

        if (i != n)
        {
            printf("FAIL:\n");
            printf("Only %d out of %d objects were processed.\n", i, n);
            abort();
        }

        for (i = 0; i < n; i++)
            fmpz_clear(a + i);
        flint_free(a);
    }

    /* Write bad data to a pipe and read it */
    /* Not necessary */

    FLINT_TEST_CLEANUP(state);
    
    printf("PASS\n");
    return EXIT_SUCCESS;
}

#else

int main(void)
{
    printf("out_raw/ inp_raw....");
    fflush(stdout);
    printf("SKIPPED\n");
    return EXIT_SUCCESS;
}

#endif
