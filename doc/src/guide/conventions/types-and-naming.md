---
title: Types, naming and aliasing
summary: The word types, how FLINT objects are declared and passed, how functions are named, what the leading underscore means, the aliasing rules and the normalisation of objects.
order: 1
---

All FLINT modules follow the same conventions for types, function names,
argument order and aliasing. The module documentation only mentions
deviations. The code conventions for
contributors are in [](guide/development/contributing), and which type
to use for a given task is discussed in [](guide/types/choosing).

## Word types {#words}

FLINT uses its own word types `ulong` and `slong`, unsigned and signed
integers of `FLINT_BITS` bits (64 on 64-bit platforms, including 64-bit
Windows, where `long` has only 32 bits). They have the size of a GMP limb.
Lengths, degrees, indices, precisions and exponents are `slong`, even when
they cannot be negative; residues modulo word-size moduli, limbs and
values that need the full word range are `ulong`. Bit counts have type
`flint_bitcnt_t`. Arrays of limbs, as in the `mpn` layer of GMP, are
passed as `nn_ptr` (writable) or `nn_srcptr` (read-only) together with
their length. Word constants are written `UWORD(123)` and `WORD(-1)`, and
words are printed with `flint_printf` using `%wu` and `%wd`. See [](flint)
and [](guide/conventions/portability).

## Objects {#objects}

Most FLINT types come in two forms. The structure type `foo_struct`
holds the data, and `foo_t` is defined as an array of one `foo_struct`:

```
typedef struct { fmpz * coeffs; slong alloc; slong length; } fmpz_poly_struct;
typedef fmpz_poly_struct fmpz_poly_t[1];
```

A variable declared as `fmpz_poly_t f;` therefore reserves the structure
itself, and passing `f` to a function passes a pointer to it, without an
`&`. Functions take `foo_t` parameters, read-only ones as `const foo_t`.
Arrays of objects are arrays of structures: a vector of polynomials is a
`fmpz_poly_struct *`, and its entries are passed as `vec + i`. Some
modules name these pointer types: `arb_ptr` and `arb_srcptr` are
`arb_struct *` and `const arb_struct *`, and similarly for `acb`, `mag`,
`arf`, `gr` and others.

Two basic types are exceptions. An `fmpz` is a single word that either
holds a small integer directly or encodes a pointer to a GMP integer;
`fmpz_t` is an array of one `fmpz`, and a vector of integers is a plain
`fmpz *`. A `nmod_t` (a word-size modulus with precomputed data) is
passed by value.

Every object must be initialised before use, with `foo_init` or a
variant, and cleared with `foo_clear` when no longer needed; see
[](guide/conventions/memory). Objects should not be copied with `=` or
`memcpy`: the copy would share the allocated memory of the original. Use
`foo_set` to copy, and `foo_swap` to exchange two objects in constant
time.

### Contexts {#contexts}

Objects whose meaning depends on parameters that are not stored in the
object itself take a *context* as their last argument: the modulus of
`fmpz_mod` and `fmpz_mod_poly` (`fmpz_mod_ctx_t`), the finite field of
`fq` (`fq_ctx_t`), the number and ordering of variables of multivariate
polynomials (`fmpz_mpoly_ctx_t`), the ring of the generic types
(`gr_ctx_t`), and so on. A context is initialised and cleared like an
object, must outlive every object used with it, and the same context
must be passed to all functions operating on an object. Objects of the
word-size modular types (`nmod_poly`, `nmod_mat`) instead store their
modulus.

## Function names {#naming}

Function names have the form `module_operation`, optionally followed by
suffixes that name the types of some arguments or the algorithm:

- *Argument types.* `fmpz_add_ui(r, a, b)` adds a `ulong` $b$;
  likewise `_si` (`slong`), `_d` (`double`), `_fmpz`, `_fmpq`, `_str`
  (a string) and `_ui_ui` or `_ui_2exp` for several such arguments.
  `fmpz_get_si` and `fmpz_set_si` convert to and from `slong`.
- *Algorithms.* Most operations have a main function that selects an
  algorithm and versions for each algorithm, such as `fmpz_poly_mul`
  and `fmpz_poly_mul_karatsuba`, `fmpz_poly_mul_KS` (Kronecker
  substitution) or `nmod_mat_mul_classical`. The main function is
  usually the right choice; the others are for experiments and special
  situations.
- *Variants.* `_preinv` and `_precomp` take precomputed data (an inverse
  or a table) as an extra argument; `_series` functions work with power
  series truncated to a given length; `_threaded` functions use threads
  explicitly, and `_threaded_pool` functions take threads already
  obtained from the thread pool (see [](guide/conventions/threading));
  `_mod` functions reduce modulo a polynomial or integer; `_inplace`
  functions overwrite an argument.
- *Common operations.* `init`, `clear`, `set`, `swap`, `zero`, `one`,
  `equal`, `is_zero`, `is_one`, `print`, `fprint`, `get_str`,
  `set_str`, `randtest` (random objects for testing, biased towards
  special cases) and `randtest_not_zero` exist for most types.

Output arguments come first, then inputs, then parameters such as the
precision or the length, and the context last. A function that may fail
returns an `int`: 1 (true) on success and 0 on failure, or for
predicates whether the property holds. The generic ring functions return
a status code (`GR_SUCCESS`, `GR_DOMAIN` or `GR_UNABLE`) and the
predicates of uncertain outcome a `truth_t` (`T_TRUE`, `T_FALSE` or
`T_UNKNOWN`); see [](gr). Errors that cannot be reported through the
return value raise an exception; see [](flint/exceptions).

## Underscore functions {#underscore}

A function name with a leading underscore, such as `_fmpz_poly_mul`,
marks a lower-level version of the function without the leading
underscore. Underscore functions typically:

- operate on raw arrays (coefficient vectors and their lengths, limbs)
  instead of on `foo_t` objects, and do not allocate or resize the
  output, which must have enough space;
- have stricter preconditions, for example nonzero lengths, normalised or
  reduced input, or a particular order of the lengths;
- do not allow the output to overlap the inputs, unless stated otherwise;
- do not check their preconditions, except in builds with
  `--enable-assert`.

They avoid the overhead of the checks and of memory management in inner
loops and allow working on parts of arrays. The documentation of each
underscore function states its preconditions. For example,
`_fmpz_poly_mul(res, poly1, len1, poly2, len2)` requires
`len1, len2 >= 1`, writes `len1 + len2 - 1` coefficients to `res`, and
requires `res` not to overlap `poly1` or `poly2` (which may coincide).

Some names start with an underscore for another reason: the function is
an internal helper, or a faster version with an extra restriction (for
example `_nmod_add` requires $n \le 2^{\bits - 1}$).
Functions with a double underscore, such as
`__flint_set_memory_functions`, are low-level hooks.

FLINT does not separate a public from a private interface strictly:
nearly all functions declared in the headers are documented, tested and
may be used. Common sense indicates which ones concern implementation
details and are more likely to change between releases. The interface
of `arb_add` or `fmpz_poly_mul` is not going to change, but that of an
internal helper such as `_arb_get_mpn_fixed_mod_pi4` might. Incompatible
changes are listed in the release notes.

## Aliasing {#aliasing}

Unless documented otherwise, the functions without a leading underscore
allow any aliasing between their arguments: an output may be the same
object as one or more inputs, and the same object may be passed for
several inputs. `fmpz_poly_mul(f, f, f)` squares `f` in place, and
`fmpz_addmul(a, a, a)` sets $a$ to $a + a^2$. Distinct outputs must be
distinct objects: in `fmpz_poly_divrem(Q, R, A, B)`, `Q` and `R` must
differ. Underscore functions generally do not allow the output to
overlap an input; where they do, the documentation says so. When objects
are aliased, functions may allocate temporary copies, so avoiding
aliasing can save time in performance-critical code.

The rule concerns whole objects of the same type, ignoring `const`, in
the spirit of the strict aliasing rule of C. An output must not be a
*part* of an input of another type unless the documentation allows it:
in `fmpz_poly_evaluate_fmpz(res, f, a)`, `res` may be the same object as
`a`, but must not be a coefficient of `f`. With $f = 3x^2 + 2x + 1$ and
$a = 10$, passing `f->coeffs` as `res` gives 352 instead of 321, since
Horner's rule overwrites the constant coefficient before reading it.
Where input and output arrays may be aliased, they must be identical or
disjoint, never partially overlapping.

The library itself sometimes relies on aliasing that a particular
implementation happens to support without documenting it. Code outside
the library should rely only on documented aliasing, since an
implementation may change.

## Normalisation {#normalisation}

Objects are kept in a canonical form, which functions may assume on
input and must produce on output:

- an `fmpz` holding a value that fits in a small word is stored directly,
  never as a GMP integer;
- an `fmpq` has coprime numerator and denominator and a positive
  denominator (`fmpq_canonicalise` restores this after setting the parts
  directly);
- a polynomial has a nonzero leading coefficient, or length zero, and its
  coefficients beyond the length are zero (`_fmpz_poly_normalise` and
  its analogues restore this after writing coefficients directly);
- residues modulo $n$ are reduced to $[0, n)$;
- the dimensions of a matrix are fixed when it is initialised; functions
  require compatible dimensions of their arguments.

Writing into the internal fields of an object is possible, but the code
must then restore the canonical form before passing the object to
library functions.

## Example {#example}

```c example name=guide_conventions_aliasing
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>

int main(void)
{
    fmpz_t a;
    fmpz_poly_t f, g;
    fmpz * u, * v, * w;

    /* aliasing: a = a + a * a */
    fmpz_init_set_ui(a, 12);
    fmpz_addmul(a, a, a);
    flint_printf("a = %{fmpz}\n", a);

    /* the same product with the polynomial type and with raw arrays */
    fmpz_poly_init(f);
    fmpz_poly_init(g);
    fmpz_poly_set_str(f, "3  1 2 3");
    fmpz_poly_mul(g, f, f);              /* aliased inputs are fine */
    fmpz_poly_mul(f, f, f);              /* so is an aliased output */
    flint_printf("g = %{fmpz_poly}, equal: %d\n", g, fmpz_poly_equal(f, g));

    u = _fmpz_vec_init(3);
    w = _fmpz_vec_init(5);
    fmpz_set_ui(u + 0, 1);
    fmpz_set_ui(u + 1, 2);
    fmpz_set_ui(u + 2, 3);
    v = u;                               /* the inputs may coincide */
    _fmpz_poly_mul(w, u, 3, v, 3);       /* w: 3 + 3 - 1 entries, no overlap */
    flint_printf("w = %{fmpz*}\n", w, (slong) 5);

    _fmpz_vec_clear(u, 3);
    _fmpz_vec_clear(w, 5);
    fmpz_poly_clear(f);
    fmpz_poly_clear(g);
    fmpz_clear(a);
    return 0;
}
```

```output
a = 156
g = 9 * x^4 + 12 * x^3 + 10 * x^2 + 4 * x + 1, equal: 1
w = [1, 4, 10, 12, 9]
```
