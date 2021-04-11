.. _fmpz-poly-roots:

**fmpz_poly_roots.h** -- roots of univariate polynomials over the integers
==========================================================================

Find roots of univariate polynomials over the integers over finite
fields, p-adic or q-adic extensions.

-----------------------------

Roots over finite fields
------------------------

Types, macros and constants
___________________________

Type containing the information of roots over a finite field.

.. type:: fmpz_polynomial_roots_fq_t
	  

Memory management
_________________

.. function:: void fmpz_poly_roots_fq_init2 (fmpz_poly_roots_fq_t roots, slong n, fq_ctx_ fctx)

	Initializes ``roots`` for use, with context ``fctx``
	to contain at most ``n`` roots. A corresponding call to
	:func:`fmpz_poly_roots_fq_clear` must be made to free the memory
	used bye the roots.

.. function:: void fmpz_poly_roots_fq_clear (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)

	Clears the given roots, releasing any memory used. It must be
	reinitialized in order to be used again.

Calculation of roots
____________________

.. function:: void fmpz_poly_roots_fq (fmpz_poly_roots_fq_t roots, fmpz_poly_t poly,  fq_ctx_t fctx)
	      
	 Takes as input a polynomial ``poly`` and a context of of a
	 finite field ``fctx`` and initialized ``roots``.

Input and output
________________

.. function:: char * fmpz_poly_roots_fq_get_str (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)
	      
	 Returns  a representation of ``roots``
	 as a zero terminated string.

.. function:: int fmpz_poly_roots_fq_print (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)

	 Prints a representation of ``roots`` to ``stdout``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.

.. function:: char * fmpz_poly_roots_fq_get_str_pretty (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)

	 Returns  a pretty representation of ``roots``
	 as a zero terminated string.

.. function:: int fmpz_poly_roots_fq_fprint_pretty (FILE * file, fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)

	 Prints a pretty representation of ``roots`` to ``file``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.
	 
.. function:: int fmpz_poly_roots_fq_print_pretty (fmpz_poly_roots_fq_t roots, fq_ctx_t fctx)

	 Prints a pretty representation of ``roots`` to ``stdout``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.


Roots over padic fields
-----------------------

Types, macros and constants
___________________________

Type containing the information of roots over a p-adic field.

.. type:: fmpz_polynomial_roots_padic_t
	  
Memory management
_________________

.. function:: void fmpz_poly_roots_padic_init2 (fmpz_poly_roots_padic_t roots, slong n, fq_ctx_ fctx)

	Initializes ``roots`` for use, with context ``fctx``
	to contain at most ``n`` roots. A corresponding call to :func:`fmpz_poly_roots_padic_clear` must be made to free the memory
	used bye the roots.

.. function:: void fmpz_poly_roots_padic_clear (fmpz_poly_roots_padic_t roots, fq_ctx_t fctx)

	Clears the given roots, releasing any memory used. It must be
	reinitialized in order to be used again.

Calculation of roots
____________________

.. function:: void fmpz_poly_roots_padic_fmpz_poly (fmpz_poly_roots_padic_t roots, fmpz_poly_t poly,  fq_ctx_t fctx)
	      
	 Takes as input a polynomial ``poly`` and a context of of a
	 finite field ``fctx`` and initialized ``roots``.

Input and output
________________


.. function:: char * fmpz_poly_roots_padic_get_str (fmpz_poly_roots_padic_t roots, fq_ctx_t fctx)

	 Returns  a representation of ``roots``
	 as a zero terminated string.
	 
.. function:: int fmpz_poly_roots_padic_fprint_pretty (FILE * file, fmpz_poly_roots_padic_t roots, fq_ctx_t fctx)

	 Prints a representation of ``roots`` to ``file``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.
	 
.. function:: int fmpz_poly_roots_padic_print_pretty (fmpz_poly_roots_padic_t roots, fq_ctx_t fctx)

	 Prints a representation of ``roots`` to ``stdout``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.

Roots over q-adic fields
------------------------

Types, macros and constants
___________________________

Type containing the information of roots over an unramified
extension over p-adic field.

.. type:: fmpz_polynomial_roots_qadic_t
 
Memory management
_________________

.. function:: void fmpz_poly_roots_qadic_init2 (fmpz_poly_roots_qadic_t roots, slong n, fq_ctx_ fctx)

	Initializes ``roots`` for use, with context ``fctx``
	to contain at most ``n`` roots. A corresponding call to
	:func:`fmpz_poly_roots_qadic_clear` must be made to free the memory
	used bye the roots.

.. function:: void fmpz_poly_roots_qadic_clear (fmpz_poly_roots_qadic_t roots, fq_ctx_t fctx)

	Clears the given roots, releasing any memory used. It must be
	reinitialized in order to be used again.

Calculation of roots
________________________

.. function:: void fmpz_poly_roots_qadic (fmpz_poly_roots_qadic_t roots, fmpz_poly_t poly,  fq_ctx_t fctx)
	      
	 Takes as input a polynomial ``poly`` and a context of of a
	 finite field ``fctx`` and initialized ``roots``.

Input and output
________________

.. function:: char * fmpz_poly_roots_qadic_get_str_pretty (fmpz_poly_roots_qadic_t roots, fq_ctx_t fctx)

	 Returns a pretty representation of ``roots`` as zero
	 terminated string.

.. function:: int fmpz_poly_roots_qadic_fprint_pretty (FILE *file, fmpz_poly_roots_qadic_t roots, fq_ctx_t fctx)

	 Prints a pretty representation of ``roots`` to ``file``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.

.. function:: int fmpz_poly_roots_qadic_print_pretty (fmpz_poly_roots_qadic_t roots, fq_ctx_t fctx)

	 Prints a pretty representation of ``roots`` to ``stdout``.
	 In case of success, returns a positive value.  In case of
	 failure, returns a non-positive value.
