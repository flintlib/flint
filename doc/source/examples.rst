.. _examples:

Example programs
===============================================================================

.. highlight:: text

The *examples* directory
(https://github.com/fredrik-johansson/calcium/tree/master/examples)
contains complete C programs illustrating use of Calcium.
Running::

    make examples

will compile the programs and place the binaries in ``build/examples``.

machin.c
-------------------------------------------------------------------------------

This program checks several variations of Machin's formula

.. math ::

    \frac{\pi}{4} = 4 \operatorname{atan}\left(\frac{1}{5}\right) - \operatorname{atan}\left(\frac{1}{239}\right)

expressing `\pi` or logarithms of small integers in terms of
arctangents or hyperbolic arctangents of rational numbers.
The program actually evaluates 
`4 \operatorname{atan}\left(\tfrac{1}{5}\right) - \operatorname{atan}\left(\tfrac{1}{239}\right) - \tfrac{\pi}{4}`
(etc.) and prints the result, which should be precisely 0
(``0  in  QQ``), proving the identity.
Inverse trigonometric functions are not yet implemented in Calcium,
so the example program evaluates them using logarithms.

Sample output::

    > build/examples/machin 
    [(1)*atan(1/1) - pi/4]   =   0  in  QQ
    [(1)*atan(1/2) + (1)*atan(1/3) - pi/4]   =   0  in  QQ
    [(2)*atan(1/2) + (-1)*atan(1/7) - pi/4]   =   0  in  QQ
    [(2)*atan(1/3) + (1)*atan(1/7) - pi/4]   =   0  in  QQ
    [(4)*atan(1/5) + (-1)*atan(1/239) - pi/4]   =   0  in  QQ
    [(1)*atan(1/2) + (1)*atan(1/5) + (1)*atan(1/8) - pi/4]   =   0  in  QQ
    [(1)*atan(1/3) + (1)*atan(1/4) + (1)*atan(1/7) + (1)*atan(1/13) - pi/4]   =   0  in  QQ
    [(12)*atan(1/49) + (32)*atan(1/57) + (-5)*atan(1/239) + (12)*atan(1/110443) - pi/4]   =   0  in  QQ

    [(14)*atanh(1/31) + (10)*atanh(1/49) + (6)*atanh(1/161) - log(2)]   =   0  in  QQ
    [(22)*atanh(1/31) + (16)*atanh(1/49) + (10)*atanh(1/161) - log(3)]   =   0  in  QQ
    [(32)*atanh(1/31) + (24)*atanh(1/49) + (14)*atanh(1/161) - log(5)]   =   0  in  QQ
    [(144)*atanh(1/251) + (54)*atanh(1/449) + (-38)*atanh(1/4801) + (62)*atanh(1/8749) - log(2)]   =   0  in  QQ
    [(228)*atanh(1/251) + (86)*atanh(1/449) + (-60)*atanh(1/4801) + (98)*atanh(1/8749) - log(3)]   =   0  in  QQ
    [(334)*atanh(1/251) + (126)*atanh(1/449) + (-88)*atanh(1/4801) + (144)*atanh(1/8749) - log(5)]   =   0  in  QQ
    [(404)*atanh(1/251) + (152)*atanh(1/449) + (-106)*atanh(1/4801) + (174)*atanh(1/8749) - log(7)]   =   0  in  QQ

    cpu/wall(s): 0.03 0.029
    virt/peak/res/peak(MB): 35.57 35.57 8.80 8.80


.. raw:: latex

    \newpage

