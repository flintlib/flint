---
title: Handling FLINT exceptions
functions: flint_throw flint_abort flint_set_abort flint_set_throw
---

When FLINT meets an error it cannot report through a return value, such
as an impossible inverse, a division by zero or a failed allocation, it
calls `flint_throw`. By default this prints a message to standard output
and aborts the program through `flint_abort`. A program that must survive
such errors, for example an interpreter or a library binding, can install
its own handler with `flint_set_throw`, or replace the final `abort` with
`flint_set_abort`.

## Writing a handler {#handler}

A handler installed with `flint_set_throw` receives the kind of error,
the format string of the message and the argument list. It must call
`va_end` on the argument list and must not return to its caller, since
the function that raised the exception cannot continue. The usual
choices are to terminate the process with a custom message, or to jump
back to a recovery point with `longjmp`, as below.

```c example name=flint_throw_longjmp
#include <setjmp.h>
#include <stdarg.h>
#include <flint/flint.h>
#include <flint/nmod.h>

static jmp_buf env;
static flint_err_t last_error;

static void my_throw(flint_err_t exc, const char * msg, va_list ap)
{
    char buf[200];
    flint_vsnprintf(buf, sizeof(buf), msg, ap);
    va_end(ap);
    last_error = exc;
    flint_printf("caught: %s", buf);
    longjmp(env, 1);
}

int main(void)
{
    nmod_t mod;
    ulong a;

    nmod_init(&mod, 15);
    flint_set_throw(my_throw);

    if (setjmp(env) == 0)
    {
        a = nmod_inv(7, mod);
        flint_printf("1/7 = %wu mod 15\n", a);
        a = nmod_inv(6, mod);
        flint_printf("1/6 = %wu mod 15\n", a);   /* not reached */
    }
    else
    {
        flint_printf("error kind is FLINT_IMPINV: %d\n", last_error == FLINT_IMPINV);
    }

    return 0;
}
```

```output
1/7 = 13 mod 15
caught: Cannot invert modulo 3*5
error kind is FLINT_IMPINV: 1
```

## Caveats {#caveats}

Jumping out of FLINT is safe only in a limited sense. The function that
threw does not clean up: memory it allocated for temporaries is leaked,
locks it holds (for example in the thread pool) are not released, and the
objects it was writing to may be left in an inconsistent state. After a
`longjmp`, such objects should only be cleared, not used. An exception
raised in a worker thread runs the handler in that thread, where a
`longjmp` to a buffer set up in another thread is undefined behaviour.
The handler is global, not per thread.

Exceptions signal bugs or resource exhaustion more often than
recoverable conditions. Where failure is expected, FLINT provides
functions that return a status instead: for example `nmod_divides` instead
of `nmod_div`, `fmpz_invmod` instead of throwing inverses, and the
`GR_DOMAIN` and `GR_UNABLE` status codes of the [generic rings](gr).
