.. _portability:

**Portability**
===============================================================================

Portable FLINT types
-------------------------------------------------------------------------------

For platform independence, FLINT provides two types ``ulong``
and ``slong`` to replace ``unsigned long`` and ``long``
respectively. These are guaranteed to be the same size as GMP's
``mp_limb_t`` and ``mp_limb_signed_t`` types, respectively.

A full list of types provided by FLINT is available in
``code_conventions.txt`` in the top-level source tree.

As FLINT supports Windows 64 on which the FLINT ``ulong`` and
``slong`` types are 64 bits, whilst ``unsigned long`` and
``long`` are only 32 bits, it is necessary to have a special
format specifier which is 64 bits on Windows 64 instead of the usual
``"%lu"`` and ``"%ld"``.

For this purpose FLINT provides its own I/O functions, ``flint_printf``,
``flint_fprintf``, ``flint_sprintf``, ``flint_scanf``,
``flint_fscanf`` and ``flint_sscanf``, which work exactly as the
usual system versions, but which take the ``"%wu"`` and ``"%wd"``
format specifiers, which support FLINT ``ulong`` and ``slong``
types respectively.

Also, instead of using constants ``123UL`` and ``123L``, FLINT
provides the macros ``UWORD(123)`` and ``WORD(123)`` respectively
for constants of type ``ulong`` and ``slong`` respectively.

The maximum and minimum values that can be represented by these types
are given by ``UWORD_MAX`` and ``WORD_MAX`` respectively.

