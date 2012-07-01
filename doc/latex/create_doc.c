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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

static char * docsin[] = {
    "../../fmpz/doc/fmpz.txt", 
    "../../fmpz_vec/doc/fmpz_vec.txt", 
    "../../fmpz_factor/doc/fmpz_factor.txt", 
    "../../fmpz_mat/doc/fmpz_mat.txt", 
    "../../fmpz_poly/doc/fmpz_poly.txt", 
    "../../fmpz_poly_factor/doc/fmpz_poly_factor.txt", 
    "../../fmpq/doc/fmpq.txt", 
    "../../fmpq_mat/doc/fmpq_mat.txt", 
    "../../fmpq_poly/doc/fmpq_poly.txt", 
    "../../fmpz_poly_q/doc/fmpz_poly_q.txt", 
    "../../fmpz_poly_mat/doc/fmpz_poly_mat.txt", 
    "../../nmod_vec/doc/nmod_vec.txt",
    "../../nmod_mat/doc/nmod_mat.txt",
    "../../nmod_poly/doc/nmod_poly.txt",
    "../../nmod_poly_mat/doc/nmod_poly_mat.txt",
    "../../fmpz_mod_poly/doc/fmpz_mod_poly.txt",
    "../../padic/doc/padic.txt", 
    "../../arith/doc/arith.txt", 
    "../../ulong_extras/doc/ulong_extras.txt",
    "../../long_extras/doc/long_extras.txt",
    "../../doc/longlong.txt",
    "../../mpn_extras/doc/mpn_extras.txt",
    "../../doc/profiler.txt", 
    "../../interfaces/doc/interfaces.txt",
    "../../fft/doc/fft.txt",
    "../../qsieve/doc/qsieve.txt",
    "../../perm/doc/perm.txt",
};

static char * docsout[] = {
    "input/fmpz.tex", 
    "input/fmpz_vec.tex", 
    "input/fmpz_factor.tex", 
    "input/fmpz_mat.tex",
    "input/fmpz_poly.tex", 
    "input/fmpz_poly_factor.tex", 
    "input/fmpq.tex", 
    "input/fmpq_mat.tex", 
    "input/fmpq_poly.tex", 
    "input/fmpz_poly_q.tex", 
    "input/fmpz_poly_mat.tex", 
    "input/nmod_vec.tex",
    "input/nmod_mat.tex",
    "input/nmod_poly.tex",
    "input/nmod_poly_mat.tex",
    "input/fmpz_mod_poly.tex",
    "input/padic.tex", 
    "input/arith.tex", 
    "input/ulong_extras.tex",
    "input/long_extras.tex",
    "input/longlong.tex", 
    "input/mpn_extras.tex",
    "input/profiler.tex", 
    "input/interfaces.tex",
    "input/fft.tex",
    "input/qsieve.tex",
    "input/perm.tex",
};

static const int ndocs = sizeof(docsin) / sizeof(char *);

#define DOCS_WIDTH     79

#define DOCS_EOF      (-1)
#define DOCS_IOE      1
#define DOCS_SUCCESS  0

static FILE *in, *out;              /* Input and output handles           */

static char buf[DOCS_WIDTH + 9];    /* Buffer for one line                */

static char *name;                  /* Current file name                  */
static int line;                    /* Current line number                */

static char grp[DOCS_WIDTH + 1];    /* Group title                        */

struct fn_t {
    char mods[DOCS_WIDTH + 1];
    char name[DOCS_WIDTH + 1];
    char args[3 * DOCS_WIDTH + 1];
};

struct fn_t fnc;          /* Function data                      */

static int grp_open = 0;  /* Whether a group section is open    */
static int fnc_open = 0;  /* Whether a function section is open */
static int dsc_open = 0;  /* Whether a description is open      */

#define FSM
#define STATE(x)      s_ ## x :
#define NEXTSTATE(x)  goto s_ ## x

#define next_event()                                                    \
do {                                                                    \
    r = readline(in, buf, &n);                                          \
    ++line;                                                             \
    if (r == DOCS_IOE)                                                  \
        NEXTSTATE(ioe);                                                 \
    if (r == DOCS_SUCCESS)                                              \
    {                                                                   \
        if (n > DOCS_WIDTH)                                             \
        {                                                               \
            printf("\n");                                               \
            printf("Warning:\n");                                       \
            printf("The parser encountered a line of length %d\n", n);  \
            printf("in line %d in file %s.\n\n", line, name);           \
        }                                                               \
    }                                                                   \
} while (0)

/*
    Reads one line from the file into the buffer c (of length at 
    least DOCS_WIDTH + 2).  The number of characters, excluding any 
    newline characters or the terminating '\0' character, written to 
    the buffer c is written to n.

    Returns zero in case of success.  Otherwise, returns one of 
    DOCS_EOF and DOCS_IOE.
 */

static int readline(FILE *file, char *c, int *n)
{
    if (fgets(c, DOCS_WIDTH + 2, file))
    {
        int i;

        for (i = 0; i < DOCS_WIDTH && c[i] != '\n'; i++) ;
        c[i] = '\0';
        *n = i;
        return DOCS_SUCCESS;
    }
    else
    {
        return feof(file) ? DOCS_EOF : DOCS_IOE;
    }
}

/*
    Returns DOCS_SUCCESS if successful and DOCS_IOE otherwise.
 */

static int printline(FILE *file, char *c)
{
    int r = 0;

    r = r || (fputs(c, file) < 0);        /* fputs >= 0 if successful     */
    r = r || (fputc('\n', file) == EOF);  /* fputc == EOF if unsuccessful */

    return r ? DOCS_IOE : DOCS_SUCCESS;
}

/*
    A really, really simple quadratic time implementation that removes 
    leading, trailing and duplicate whitespace from a string.
 */

static void _str_cleanup(char * str)
{
    int i, j, k, len = strlen(str);

    /* Remove trailing whitespace */
    for ( ; len >= 0 && str[len] == ' '; len--)
        str[len] = '\0';

    /* Remove leading whitespace */
    for (j = 0; j < len && str[j] == ' '; j++) ;
    if (j > 0)
        for (k = j; k <= len; k++)
            str[k - j] = str[k];
    len = len - j;

    /* Remove intermediate whitespace */
    for (i = 1; i < len; i++)
    {
        if (str[i] == ' ' && i < len)
            i++;
        for (j = i; j < len && str[j] == ' '; j++) ;
        if (j - i > 0)
            for (k = j; k <= len; k++)
                str[k - (j - i)] = str[k];
        len = len - (j - i);
    }
}

/*****************************************************************************/

static int open_group()
{
    grp_open = 1;

    fprintf(out, "\n");
    fprintf(out, "\\section{%s}\n\n", grp);

    return DOCS_SUCCESS;
}

static int close_group()
{
    grp_open = 0;

    return DOCS_SUCCESS;
}

static int open_function()
{
    fnc_open = 1;

    _str_cleanup(fnc.mods);
    _str_cleanup(fnc.args);

    fprintf(out, "\n");
    fprintf(out, "\\vspace*{0.5em}\n");
    fprintf(out, "\\begin{lstlisting}\n");
    if (fnc.mods[0] != '\0')
        fprintf(out, "%s %s(%s)\n", fnc.mods, fnc.name, fnc.args);
    else
        fprintf(out, "%s(%s)\n", fnc.name, fnc.args);
    fprintf(out, "\\end{lstlisting}\n");
    fprintf(out, "\\vspace*{-0.5em}\n");

    return DOCS_SUCCESS;
}

static int close_function()
{
    fnc_open = 0;

    return DOCS_SUCCESS;
}

static int open_description()
{
    dsc_open = 1;

    return DOCS_SUCCESS;
}

static int close_description()
{
    dsc_open = 0;

    return DOCS_SUCCESS;
}

/*****************************************************************************/

/*
    Proceeds until the first line with 79 characters equal to '*' has 
    been read.  The return values are the same as for readline().
 */

static void processfile(void)
{
    int i, j, k, m, n, r;

    printf("Open file %s\n", name);
    fflush(stdout);
    line = 0;

    STATE(0)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == DOCS_WIDTH)
            {
                for (i = 0; i < n && buf[i] == '*'; i++) ;
                if (i == n)
                    NEXTSTATE(1);
            }
            NEXTSTATE(0);
        }
        NEXTSTATE(pe);
    }

    STATE(1)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == 0)
                NEXTSTATE(1);
            if (n > 4 && (buf[0] == ' ' && buf[1] == ' ' && buf[2] == ' ' 
                                        && buf[3] == ' ' && (isalpha(buf[4]) || buf[4] == '_')))
            {
                strncpy(grp, buf + 4, n - 4);
                grp[n - 4] = '\0';
                open_group();

                NEXTSTATE(2);
            }
        }
        NEXTSTATE(pe);
    }

    STATE(2)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == 0 || (n == 4 && buf[0] == ' ' && buf[1] == ' ' 
                                  && buf[2] == ' ' && buf[3] == ' '))
                NEXTSTATE(3);
        }
        NEXTSTATE(pe);
    }

    STATE(3)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == DOCS_WIDTH)
            {
                for (i = 0; i < n && buf[i] == '*'; i++) ;
                if (i == n)
                {
                    close_group();
                    NEXTSTATE(4);
                }
            }
            if (n == 0 || (n == 4 && buf[0] == ' ' && buf[1] == ' ' 
                                  && buf[2] == ' ' && buf[3] == ' '))
            {
                printline(out, "");
                NEXTSTATE(3);
            }
            if (n > 4 && (buf[0] == ' ' && buf[1] == ' ' && buf[2] == ' ' 
                                        && buf[3] == ' '))
            {
                printline(out, buf + 4);
                NEXTSTATE(3);
            }
        }
        NEXTSTATE(pe);
    }

    STATE(4)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == 0)
                NEXTSTATE(5);
        }
        NEXTSTATE(pe);
    }

    STATE(5)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == 0)
                NEXTSTATE(5);
            if (n == DOCS_WIDTH)
            {
                for (i = 0; i < n && buf[i] == '*'; i++) ;
                if (i == n)
                {
                    /* Open group */
                    NEXTSTATE(1);
                }
            }
            if (isalpha(buf[0]) || buf[0] == '_')
                NEXTSTATE(5a);
        }
        NEXTSTATE(eof);
    }

    /* isalpha(buf[0]) ---> Open a new function */
    STATE(5a)
    {
        close_description();
        close_function();

        for (j = 0; j < n && buf[j] != '('; j++) ;
        for (k = j; k < n && buf[k] != ')'; k++) ;

        /* No opening bracket. */
        if (j == n)
        {
            strncpy(fnc.mods, buf, n + 1);
            NEXTSTATE(6);
        }

        for (m = j; m > 0 && buf[m - 1] != ' '; m--) ;
        if (m == 0)
        {
            /* No modifiers */
            fnc.mods[0] = '\0';
        }
        else
        {
            strncpy(fnc.mods, buf, m - 1);
            fnc.mods[m - 1] = '\0';
        }
        strncpy(fnc.name, buf + m, j - m);
        fnc.name[j - m] = '\0';

        if (k - (j + 1) > 0)
        {
            strncpy(fnc.args, buf + j + 1, k - (j + 1));
            fnc.args[k - (j + 1)] = '\0';
        }

        if (k == n)
            NEXTSTATE(7);

        NEXTSTATE(5z);
    }

    STATE(5z)
    {
        open_function();
        NEXTSTATE(8);
    }

    STATE(6)
    {
        next_event();

        if (r == DOCS_SUCCESS && n > 0)
        {
            for (j = 0; j < n && buf[j] != '('; j++) ;
            for (k = j; k < n && buf[k] != ')'; k++) ;
            if (j == n)
                NEXTSTATE(pe);
            strncpy(fnc.name, buf, j);
            fnc.name[j] = '\0';
            if (j < n)
            {
                m = (k < n) ? k - 1 : n;
                strncpy(fnc.args, buf + (j + 1), m - (j + 1));
                fnc.args[n - (j + 1)] = '\0';
                if (k < n)
                    NEXTSTATE(5z);
                else
                    NEXTSTATE(7);
            }
        }
        NEXTSTATE(pe);
    }

    STATE(7)
    {
        next_event();

        if (r == DOCS_SUCCESS && n > 0)
        {
            size_t len = strlen(fnc.args);

            for (k = 0; k < n && buf[k] != ')'; k++) ;
            strncpy(fnc.args + len, buf, k);
            fnc.args[len + k] = '\0';
            if (k < n)
                NEXTSTATE(5z);
            else
                NEXTSTATE(7);
        }
        NEXTSTATE(pe);
    }

    STATE(8)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == 0 || (n == 4 && buf[0] == ' ' && buf[1] == ' ' 
                                  && buf[2] == ' ' && buf[3] == ' '))
            {
                open_description();
                NEXTSTATE(9);
            }
            if (isalpha(buf[0]) || buf[0] == '_')
                NEXTSTATE(5a);
            NEXTSTATE(pe);
        }
        NEXTSTATE(eof);
    }

    STATE(9)
    {
        next_event();

        if (r == DOCS_SUCCESS)
        {
            if (n == 0)
            {
                /* Print empty line to LaTeX */
                printline(out, "");
                NEXTSTATE(9);
            }
            if (n >= 4 && buf[0] == ' ' && buf[1] == ' ' 
                       && buf[2] == ' ' && buf[3] == ' ')
            {
                /* Print {buf + 4, n - 4} to LaTeX */
                printline(out, buf + 4);
                NEXTSTATE(9);
            }
            if (n == DOCS_WIDTH)
            {
                for (i = 0; i < n && buf[i] == '*'; i++) ;
                if (i == n)
                    NEXTSTATE(1);
            }
            if (isalpha(buf[0]) || buf[0] == '_')
                NEXTSTATE(5a);
            NEXTSTATE(pe);
        }
        NEXTSTATE(eof);
    }

    STATE(ioe)
    {
        printf("\n");
        printf("Input/ output exception:\n");
        printf("Raised when reading line %d in file %s.\n\n", line, name);

        NEXTSTATE(end);
    }

    STATE(pe)
    {
        printf("\n");
        printf("Parse exception:\n");
        printf("Encountered malformed input on line %d \n", line);
        printf("in file %s.\n\n", name);

        NEXTSTATE(end);
    }

    STATE(eof)
    {
        NEXTSTATE(end);
    }

    STATE(end)
    {
        close_description();
        close_function();
        close_group();

        printf("Close file %s\n", name);
        fflush(stdout);
    }
}

int main(void)
{
    int i;

    for (i = 0; i < ndocs; i++)
    {
        name = docsin[i];
        line = 0;
        in   = fopen(docsin[i], "r");
        out  = fopen(docsout[i], "w");

        processfile();

        fclose(in);
        fclose(out);
    }

    return EXIT_SUCCESS;
}

