.. _ca-ext:

**ca_ext.h** -- real and complex extension numbers
===============================================================================

A :type:`ca_ext_t` represents a fixed real or complex number *a*.
The content of a :type:`ca_ext_t` can be one of the following:

* An algebraic constant represented
  in canonical form by a :type:`qqbar_t` instance (example: `i`, represented
  as the root of `x^2+1` with positive imaginary part).
* A constant of the form `f(x_1, \ldots, x_n)` where *f* is
  a builtin symbolic function and `x_1, \ldots, x_n` are given :type:`ca_t`
  instances.
* A builtin symbolic constant such as `\pi`. (This is just a special
  case of the above with a zero-length argument list.)
* (Not implemented): a user-defined constant or function defined by suppling
  a function pointer for Arb numerical evaluation to specified precision.

The :type:`ca_ext_t` structure is heavy-weight object, not just meant to act
as a node in a symbolic expression. It will cache various data to support
repeated computation with this particular number, including its numerical
enclosure and number field data in the case of algebraic numbers.

Extension numbers are used internally by the :type:`ca_t` type to
define the embeddings `\mathbb{Q}(a) \to \mathbb{C}` of formal fields.
The user does not normally need to create :type:`ca_ext_t` instances
directly; the intended way for the user to work with the extension number *a*
is to create a :type:`ca_t` representing the field element `1 \cdot a`.
The underlying :type:`ca_ext_t` may be accessed to
determine symbolic and numerical properties of this number.

Since extension numbers may depend recursively on nontrivial fields for
function arguments, :type:`ca_ext_t` operations require a :type:`ca_ctx_t`
context object.

Type and macros
-------------------------------------------------------------------------------

For all types, a *type_t* is defined as an array of length one of type
*type_struct*, permitting a *type_t* to be passed by reference.

.. type:: ca_ext_struct

.. type:: ca_ext_t

    An extension number object contains a header, a hash value,
    data (a :type:`qqbar_t`
    instance and an Antic :type:`nf_t` in the case of algebraic numbers, and
    a pointer to arguments
    in the case of a symbolic function), and a cached :type:`acb_t` enclosure
    (in the case of a :type:`qqbar_t`, the enclosure internal to that
    structure is used).

.. type:: ca_ext_ptr

   Alias for ``ca_ext_struct *``.

.. type:: ca_ext_srcptr

   Alias for ``const ca_ext_struct *``.

.. macro:: CA_EXT_HEAD(x)

    Accesses the head (a :type:`calcium_func_code`) of *x*.
    This is *CA_QQBar* if *x* represents an algebraic constant in
    canonical form, and *CA_Exp*, *CA_Pi*, etc. for symbolic functions
    and constants.

.. macro:: CA_EXT_HASH(x)

    Accesses the hash value of *x*.

.. macro:: CA_EXT_QQBAR(x)

    Assuming that *x* represents an algebraic constant in canonical form,
    accesses this :type:`qqbar_t` object.

.. macro:: CA_EXT_QQBAR_NF(x)

    Assuming that *x* represents an algebraic constant in canonical form,
    accesses the corresponding Antic number field :type:`nf_t` object.

.. macro:: CA_EXT_FUNC_ARGS(x)

    Assuming that *x* represents a symbolic constant or function,
    accesses the argument list (as a :type:`ca_ptr`).

.. macro:: CA_EXT_FUNC_NARGS(x)

    Assuming that *x* represents a symbolic constant or function,
    accesses the number of function arguments.

.. macro:: CA_EXT_FUNC_ENCLOSURE(x)

    Assuming that *x* represents a symbolic constant or function,
    accesses the cached :type:`acb_t` numerical enclosure.

.. macro:: CA_EXT_FUNC_PREC(x)

    Assuming that *x* represents a symbolic constant or function,
    accesses the working precision of the cached numerical enclosure.

Memory management
-------------------------------------------------------------------------------

.. function:: void ca_ext_init_qqbar(ca_ext_t res, const qqbar_t x, ca_ctx_t ctx)

    Initializes *res* and sets it to the algebraic constant *x*.

.. function:: void ca_ext_init_const(ca_ext_t res, calcium_func_code func, ca_ctx_t ctx)

    Initializes *res* and sets it to the constant defined by *func*
    (example: *func* = *CA_Pi* for `x = \pi`).

.. function:: void ca_ext_init_fx(ca_ext_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)

    Initializes *res* and sets it to the univariate function value `f(x)`
    where *f* is defined by *func*  (example: *func* = *CA_Exp* for `e^x`).

.. function:: void ca_ext_init_fxy(ca_ext_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)

    Initializes *res* and sets it to the bivariate function value `f(x, y)`
    where *f* is defined by *func*  (example: *func* = *CA_Pow* for `x^y`).

.. function:: void ca_ext_init_fxn(ca_ext_t res, calcium_func_code func, ca_srcptr x, slong nargs, ca_ctx_t ctx)

    Initializes *res* and sets it to the multivariate function value
    `f(x_1, \ldots, x_n)` where *f* is defined by *func* and *n* is
    given by *nargs*.

.. function:: void ca_ext_init_set(ca_ext_t res, const ca_ext_t x, ca_ctx_t ctx)

    Initializes *res* and sets it to a copy of *x*.

.. function:: void ca_ext_clear(ca_ext_t res, ca_ctx_t ctx)

    Clears *res*.

Structure
-------------------------------------------------------------------------------

.. function:: slong ca_ext_nargs(const ca_ext_t x, ca_ctx_t ctx)

    Returns the number of function arguments of *x*.
    The return value is 0 for any algebraic constant and for any built-in
    symbolic constant such as `\pi`.

.. function:: void ca_ext_get_arg(ca_t res, const ca_ext_t x, slong i, ca_ctx_t ctx)

    Sets *res* to argument *i* (indexed from zero) of *x*.
    This calls *flint_abort* if *i* is out of range.

.. function:: ulong ca_ext_hash(const ca_ext_t x, ca_ctx_t ctx)

    Returns a hash of the structural representation of *x*.

.. function:: int ca_ext_equal_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)

    Tests *x* and *y* for structural equality, returning 0 (false) or 1 (true).

.. function:: int ca_ext_cmp_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)

    Compares the representations of *x* and *y* in a canonical sort order,
    returning -1, 0 or 1. This only performs a structural comparison
    of the symbolic representations; the return value does not say
    anything meaningful about the numbers represented by *x* and *y*.

Input and output
-------------------------------------------------------------------------------

.. function:: void ca_ext_print(const ca_ext_t x, const ca_ctx_t ctx)

    Prints a description of *x* to standard output.

Numerical evaluation
-------------------------------------------------------------------------------

.. function:: void ca_ext_get_acb_raw(acb_t res, ca_ext_t x, slong prec, ca_ctx_t ctx)

    Sets *res* to an enclosure of the numerical value of *x*.
    A working precision of *prec* bits is used for the evaluation,
    without adaptive refinement.

Cache
-------------------------------------------------------------------------------

.. type:: ca_ext_cache_struct

.. type:: ca_ext_cache_t

    Represents a set of structurally distinct :type:`ca_ext_t` instances.
    This object contains an array of pointers to individual heap-allocated
    :type:`ca_ext_struct` objects as well as a hash table for quick
    lookup.

.. function:: void ca_ext_cache_init(ca_ext_cache_t cache, ca_ctx_t ctx)

    Initializes *cache* for use.

.. function:: void ca_ext_cache_clear(ca_ext_cache_t cache, ca_ctx_t ctx)

    Clears *cache*, freeing the memory allocated internally.

.. function:: ca_ext_ptr ca_ext_cache_insert(ca_ext_cache_t cache, const ca_ext_t x, ca_ctx_t ctx)

    Adds *x* to *cache* without duplication. If a structurally identical
    instance already exists in *cache*, a pointer to that instance is returned.
    Otherwise, a copy of *x* is inserted into *cache* and a pointer to that new
    instance is returned.


.. raw:: latex

    \newpage

