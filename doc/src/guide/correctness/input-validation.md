---
title: Invalid input
summary: Why most FLINT functions do not check their input, which classes of functions do check (and how they report errors), and typical preconditions that are not checked.
order: 3
---

Most FLINT functions assume that their input is valid: that objects are
initialised and in canonical form, that residues are reduced, that
arguments respect the documented aliasing rules and preconditions. They
do not check this, and the behaviour for invalid input is undefined in
the sense of C: the function may return a wrong result, corrupt memory,
crash, or appear to work. Checking every precondition would cost time
in the innermost loops, where FLINT spends most of its time, and many
properties (such as primality of a modulus) are expensive to check at
all.

## What is checked {#checked}

Some errors are checked in all builds, because they are cheap to detect
and easy to make, or because the function cannot know in advance
whether the input is valid:

- *Division by zero and impossible inverses* in the basic arithmetic of
  the exact types raise an exception with `flint_throw`: for example
  `fmpz_fdiv_q`, `fmpz_divexact` and `fmpq_div` with a zero divisor, `nmod_inv` and `fmpz_mod_inv` of a non-invertible residue,
  `fmpz_poly_div` and `nmod_poly_divrem` by the zero polynomial, and
  `nmod_poly_divrem` when the leading coefficient of the divisor is not
  invertible modulo a composite modulus. Context constructors check
  their parameters in the same way, for instance `fmpz_mod_ctx_init`
  with a nonpositive modulus. By default an exception prints a message
  and aborts the program; [](flint/exceptions) explains how to install a
  handler instead.
- *The generic rings* report errors in their return value instead:
  every `gr` function returns `GR_SUCCESS`, `GR_DOMAIN` (the result is
  not defined, such as $1/0$ or $1/2$ in $\Z$) or `GR_UNABLE` (the
  implementation cannot decide or compute it). Functions of the other
  modules that are wrappers of `gr` functions often turn such a status
  into an exception, as `nmod_poly_divrem` does.
- *Ball arithmetic* never fails on numerical grounds: dividing by a
  ball that contains zero gives an indeterminate ball (`arb_div` returns
  NaN with infinite radius), and functions outside their domain give
  indeterminate results; see [](guide/types/ball-arithmetic).
- *Exact numbers.* The `ca_t` type has special values for undefined
  and infinite results: `ca_div` by zero gives `UnsignedInfinity`, and
  $0/0$ gives `Undefined`. `qqbar_div` and `qqbar_inv` throw on
  division by zero.
- *Parsers* return an error code for strings they cannot read, but the
  conventions differ: `fmpz_set_str`, `fmpq_set_str` and
  `fmpz_poly_set_str` return $-1$ on failure and 0 on success,
  `arb_set_str` returns nonzero on failure, `nmod_poly_set_str` returns
  0 on failure and 1 on success, and the generic `gr_set_str` returns a
  status other than `GR_SUCCESS`. Checking the syntax is not the same as
  checking the value: `nmod_poly_set_str` does not reduce the
  coefficients it reads, and parsing is not hardened against malicious
  input; see [](guide/correctness/security#untrusted).
- *Allocation failure* is detected by `flint_malloc` and its relatives,
  which throw instead of returning `NULL`; see
  [](guide/correctness/limits#memory).

The documentation of each function says when it checks its input and
how it reports errors.

## What is not checked {#unchecked}

The following preconditions are typical of those that are assumed and
not checked (except possibly by assertions, see
[](guide/correctness/testing#assertions)):

- *Reduced residues.* The functions of [](nmod) and [](fmpz_mod)
  require operands in $[0, n)$. `nmod_add(10, 3, mod)` with $n = 7$
  returns 13, which is not a residue at all, and `fmpz_mod_add` behaves
  in the same way. Use `nmod_set_ui`, `NMOD_RED` or `fmpz_mod_set_fmpz`
  to reduce values first.
- *Canonical form.* Polynomials must be normalised (no zero leading
  coefficient), rational numbers must be in lowest terms with a
  positive denominator, and so on; see
  [](guide/conventions/types-and-naming#normalisation). Code that writes
  to the internal fields of an object must restore the canonical form.
- *Nonzero arguments where the function does not check.* `fmpq_inv`
  of zero silently produces the invalid value $1/0$, although `fmpq_div`
  by zero throws. The documentation of each function says which
  arguments it requires to be nonzero, prime, positive and so on.
- *Mathematical preconditions.* Functions for prime moduli or finite
  fields assume that the modulus is prime and do not test it; functions
  computing in a quotient ring assume that the modulus polynomial
  satisfies the documented conditions (monic, irreducible, ...).
  Violating these gives wrong results or failures in unrelated places.
- *Dimensions and lengths.* `fmpz_mat_mul` and `nmod_mat_mul` assume
  compatible dimensions; the underscore functions assume that arrays
  have the stated lengths, that lengths are positive where required and
  that outputs have room for the result.
- *Aliasing.* Arguments must respect the aliasing rules of
  [](guide/conventions/types-and-naming#aliasing). An underscore
  function given overlapping input and output arrays where this is not
  allowed overwrites its input while still reading it.
- *Initialisation.* Every object must have been initialised and not yet
  cleared, and must be used with the context it belongs to.

```c example name=guide_correctness_validation
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/fmpz.h>
#include <flint/gr.h>

int main(void)
{
    nmod_t mod;
    fmpz_t a;
    gr_ctx_t ZZ;
    gr_ptr x, y;
    int status;

    /* not checked: nmod functions assume reduced operands */
    nmod_init(&mod, 7);
    flint_printf("nmod_add(10, 3) mod 7 = %wu\n", nmod_add(10, 3, mod));
    flint_printf("reduced first:          %wu\n",
        nmod_add(nmod_set_ui(10, mod), 3, mod));

    /* checked: parsers report invalid strings */
    fmpz_init(a);
    flint_printf("fmpz_set_str(\"12a\") returns %d\n", fmpz_set_str(a, "12a", 10));
    fmpz_clear(a);

    /* checked: generic rings return a status code */
    gr_ctx_init_fmpz(ZZ);
    x = gr_heap_init(ZZ);
    y = gr_heap_init(ZZ);
    status = gr_set_si(x, 1, ZZ);
    status |= gr_set_si(y, 2, ZZ);
    status |= gr_div(x, x, y, ZZ);
    flint_printf("1/2 in ZZ: status %d (GR_DOMAIN = %d)\n", status, GR_DOMAIN);
    gr_heap_clear(x, ZZ);
    gr_heap_clear(y, ZZ);
    gr_ctx_clear(ZZ);

    return 0;
}
```

```output
nmod_add(10, 3) mod 7 = 13
reduced first:          6
fmpz_set_str("12a") returns -1
1/2 in ZZ: status 1 (GR_DOMAIN = 1)
```

## Finding invalid input {#finding}

An assertion build (`--enable-assert`) checks many of the preconditions
above, and is the first thing to try when a program gives strange
results; see [](guide/correctness/testing#assertions). Valgrind and the
address sanitizer find uses of uninitialised or cleared objects and
out-of-bounds accesses. Programs that must cope with arbitrary input,
for example interpreters and language bindings, should check the
preconditions themselves before calling FLINT, use the generic `gr`
interface where errors are reported through status codes, and install
an exception handler with `flint_set_throw`.
