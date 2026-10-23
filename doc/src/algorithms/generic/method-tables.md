---
title: Contexts, method tables and specialisation
summary: How rings and elements are represented, how a generic call is dispatched through the method table, the three levels of operations with their defaults, specialised parametric types, and why contexts are cheap to create.
order: 10
---

The generic rings framework [`gr`](gr) rests on three data structures,
the context object, the method table and the element layout, and on a
layering of generic defaults and specialised methods that lets generic
code run at the speed of type-specific code where it matters
[@Joh2025gr, §3, §5–7].

## Contexts and method tables {#tables}

A `gr_ctx_struct` consists of six words of ring parameters (`data`), the
implementation code `which_ring`, the element size `sizeof_elem`, a
pointer `methods` to the method table, and a soft size limit
`size_limit`. The method table is an array of `GR_METHOD_TAB_SIZE`
function pointers indexed by the enumeration `gr_method`, one entry for
each operation of the interface: currently 584, from
`GR_METHOD_INIT` and `GR_METHOD_ADD` to `GR_METHOD_MAT_CHARPOLY` and the
special functions.

A wrapper such as `gr_add` is an inline function that loads the entry
`GR_METHOD_ADD` from `ctx->methods`, casts it to the signature of the
operation and calls it:

```c
static inline int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return GR_BINARY_OP(ctx, ADD)(res, x, y, ctx);
}
```

All contexts of one implementation share a single static table, built on
the first call to the constructor by `gr_method_tab_init`. That function
fills every entry with `gr_not_implemented`, then installs the defaults
of [`gr_generic`](gr_generic) (345 entries), and finally the methods
listed by the implementation, so that an implementation lists only what
it provides. Building a table takes about half a microsecond
[@Joh2025gr] and happens once per type per process.

## Elements {#elements}

Elements carry no header: an element of a ring is exactly an object of
the ring's element type, `sizeof_elem` bytes, and all type information
is in the context. An array of elements is therefore a packed C array,
the same array that the type-specific modules use. This is what allows
type-specific functions to call generic algorithms on their own data
without conversion: `_nmod_poly_gcd_hgcd` builds an `nmod` context on
the stack from its `nmod_t` with `_gr_ctx_init_nmod` and calls
`_gr_poly_gcd_hgcd` on its coefficient arrays, and many functions of
`nmod_poly`, `fmpz_mod_poly`, `fmpz_poly`, `arb_poly`, `acb_poly` and
`fmpz_mat` are thin wrappers of this kind.

Element types are *shallow* when they own no heap memory (`ulong` for
$\Z/n\Z$, the limb arrays of `mpn_mod`, `nfloat`) and *deep* otherwise
(`fmpq`, `arb`, polynomials); `fmpz` is a hybrid, storing integers below
$2^{62}$ inline and larger ones on the heap. Shallow elements need no
initialisation in the type-specific code, can be allocated on the stack
in large numbers and copied with `memcpy`. All element types must be
relocatable (movable with `memcpy`), since `gr_vec_t` grows by
reallocation.

## The cost of dispatch {#dispatch}

A generic call costs a load from the method table and an indirect call
that the compiler cannot inline: a few nanoseconds. In the measurements
of [@Joh2025gr, §5] for the loop `z_i = x_i + y_i` over `fmpz`
integers, the overhead of an elementwise generic call is about 30% for
single-word integers and about 10% for multiword integers, and
disappears completely at the vector level, where `_gr_vec_add` is
overridden by `_fmpz_vec_add`. This motivates the layering below.

## Three levels of operations {#levels}

Operations are grouped in three levels, and each method has a default
built on the levels below it (or on other methods at its level):

| Level | Examples | Defaults |
|---|---|---|
| 1. elementwise | `gr_add`, `gr_mul`, `gr_addmul`, `gr_mul_ui`, `gr_pow_ui` | required (`add`, `mul`, ...) or built from required methods |
| 2. vectors | `_gr_vec_add`, `_gr_vec_mul_scalar`, `_gr_vec_dot`, `_gr_vec_dot_rev`, `_gr_vec_sum` | loops over level 1 through the method table |
| 3. higher | `_gr_poly_mullow`, `gr_mat_mul`, `_gr_poly_divrem`, `gr_mat_lu`, `_gr_poly_gcd` | algorithms built on levels 2 and 3 |

A new ring has to implement only sixteen methods (memory management,
`randtest`, `write`, `zero`, `one`, `equal`, four conversions, `neg`,
`add`, `sub`, `mul`; see [](gr/implementing)), and every operation of
levels 2 and 3 then works. The chain for a polynomial product shows how
the levels interact: `gr_poly_mul` calls the method `POLY_MULLOW`, whose
default is classical multiplication (`_gr_poly_mullow_classical`); that
computes each coefficient with `_gr_vec_dot_rev`, which calls the
method `VEC_DOT_REV`, whose default multiplies and adds through the
`MUL` and `ADD` methods.

Overriding a method at any level speeds up everything above it:

* Overriding level 1 (for example a fused `ADDMUL`) saves temporaries.
* Overriding level 2 removes the indirect call per element, allows SIMD
  code, and allows delayed reduction or normalisation: a dot product
  modulo $n$ accumulates in two or three words and reduces once
  ([](algorithms/linear-algebra/dot-products)), and a dot product of
  balls sums exactly and bounds the error once [@Joh2019dot]. Since
  classical polynomial multiplication, matrix multiplication and
  triangular solving consist of dot products, this single override
  makes them as fast as type-specific code.
* Overriding level 3 replaces an algorithm by an asymptotically faster
  one where the ring allows it: `nmod` contexts install
  `nmod_poly_mullow` (Kronecker substitution and FFT) for `POLY_MULLOW`
  and `nmod_mat_mul` for `MAT_MUL`. Generic algorithms of higher
  complexity (division by Newton iteration, half-GCD, block LU
  decomposition, characteristic polynomials) are written in terms of
  polynomial and matrix multiplication and inherit the speed-up.

## Parametric types {#parametric}

Because the element size is fixed only when a context is created, a
single implementation can represent elements of a size chosen at run
time, store them shallowly, and switch internally to kernels for that
size. Two families use this [@Joh2025gr, §6]:

* [`mpn_mod`](mpn_mod) represents $\Z/n\Z$ for
  $2^{\bits} \le n < 2^{16 \bits}$ with elements of exactly as many limbs
  as $n$ (2 to 16, `MPN_MOD_MIN_LIMBS` to `MPN_MOD_MAX_LIMBS`).
  Compared with `fmpz_mod`, whose elements are `fmpz` values that may
  live on the heap, it avoids memory management and uses fixed-length
  multiplication and delayed reduction
  ([](algorithms/arithmetic/multiword-overview)). It gains from faster
  basecase arithmetic, less temporary allocation in divide-and-conquer
  algorithms, and cheaper additions that move the cutoffs of Karatsuba
  and Waksman multiplication down.
* [`nfloat`](nfloat) represents floating-point numbers with 1 to 66 limbs
  of precision (`NFLOAT_MAX_LIMBS` on 64-bit machines, 64 to 4224 bits)
  as a fixed-length mantissa with an inline exponent and sign, and
  `nfloat_complex` complex numbers.

All sizes share one method table, and the methods dispatch on the
number of limbs stored in the context. Templates or just-in-time
compilation would instead generate a copy of every generic algorithm for
every size; here only the critical kernels are specialised and the rest
of the code runs through run-time dispatch, without a compile-time
explosion. Elements should stay within a few kilobytes, since generic
code allocates temporaries on the stack; both limits above respect this.

## Cheap contexts {#cheap}

Multimodular algorithms may create a very large number of rings: a
determinant or a resultant over $\Z$ reduces modulo thousands of primes,
each a ring $\Z/p\Z$, possibly with a polynomial ring $(\Z/p\Z)[x]$ on
top. Most constructors therefore allocate nothing and store only a few
words: `gr_ctx_init_nmod` sets the code, the size and the method table
pointer and calls `nmod_init`, which takes nanoseconds once the shared
table exists, and `_gr_ctx_init_nmod` copies an existing `nmod_t`
instead of recomputing it.
Type-specific code exploits this by creating contexts on the stack
inside its functions, as in `fmpz_poly_interpolate_multi_mod`, which
creates an `nmod` context for each prime and calls
`_gr_poly_interpolate_fast`. Exceptions are rings that need
precomputation: number fields, finite fields (which store the defining
polynomial and reduction data), Calcium fields, and quotient rings
(which precompute an inverse of the modulus).
