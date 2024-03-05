dnl
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

dnl TODO
dnl  * Can we get a better asymptotical algorithm? Currently only ~8% faster
dnl    asymptotically on M1 processors compared to mul_1 + addmul_1.
dnl  * Would be nice to avoid pushing r13 to stack.

define(`rp',  `x0')
define(`ap',  `x1')
define(`n',   `x2')
define(`bp',  `x3')

define(`b0',  `x3')
define(`b1',  `x4')

define(`r0',  `x5')
define(`r1',  `x6')
define(`r2',  `x7')
define(`r3',  `x8')
define(`r4',  `x9')
define(`r5',  `x10')
define(`r6',  `x11')
define(`r7',  `x12')
define(`r8',  `x13')
define(`r9',  `x14')
define(`r10', `x15')

define(`cy0', `x16')
define(`cy1', `x17')

define(`r11', `x19')
define(`r12', `x20')
define(`r13', `x21')

define(`zr', `xzr')

PROLOGUE(flint_mpn_mul_2)
	tst	n, #1		C Also clears C, needed if n % 4 == 0
	ldp	b0, b1, [bp]
	mov	cy0, zr		C Set carries to zero, needed if n % 4 == 0
	mov	cy1, zr
	b.eq	L(m2)

	C First shave off n % 4 to end with a nice unrolled 4-loop (no pun intended)

L(m1):	ldr	r0, [ap], #1*8
	mul	r1, r0, b0	C a0 b0
	umulh	r2, r0, b0
	mul	cy0, r0, b1	C a0 b1
	umulh	cy1, r0, b1
	str	r1, [rp], #1*8
	adds	cy0, cy0, r2
	cinc	cy1, cy1, cs

L(m2):	tst	n, #2
	ldp	r0, r1, [ap], #2*8
	lsr	n, n, #2
	b.eq	L(pre)

	mul	r2, r0, b0	C a0 b0
	umulh	r3, r0, b0
	mul	r4, r1, b0	C a1 b0
	umulh	r5, r1, b0

	mul	r6, r0, b1	C a0 b1
	umulh	r0, r0, b1
	mul	r7, r1, b1	C a1 b1
	umulh	r8, r1, b1
	C (cy0, 2), (cy1, 3, 4, 6), (5, 0, 7), 8

	adds	r2, r2, cy0
	adcs	r3, r3, cy1
	adcs	r5, r5, r0
	cinc	r8, r8, cs
	C 2, (3, 4, 6), (5, 7), 8

	adds	r4, r4, r3
	adcs	cy0, r5, r7
	cinc	cy1, r8, cs
	C 2, (4, 6), cy0, cy1

	adds	r4, r4, r6
	adcs	cy0, cy0, zr	C Needed if n = 2
	adcs	cy1, cy1, zr	C Clear C, needed by loop
	stp	r2, r4, [rp], #2*8

	cbz	n, L(end)	C Could be placed in L(pre), but we assume n >= 2.

	ldp	r0, r1, [ap], #2*8
L(pre):	ldp	r2, r3, [ap], #2*8
	stp	r11, r12, [sp,#-3*8]
	str	r13, [sp,#-1*8]
	b	L(ent)

	ALIGN(16)
L(top):	ldp	r0, r1, [ap], #2*8
	ldp	r2, r3, [ap], #2*8

	stp	r4, r5, [rp], #2*8
	stp	r7, r10, [rp], #2*8

L(ent):	sub	n, n, #1

	mul	r4, r0, b0	C a0 b0
	umulh	r5, r0, b0

	mul	r6, r1, b0	C a1 b0
	umulh	r7, r1, b0
	mul	r8, r0, b1	C a0 b1
	umulh	r0, r0, b1

	mul	r9, r2, b0	C a2 b0
	umulh	r10, r2, b0

	adcs	r4, r4, cy0	C Add carry from last loop
	adcs	r5, r5, cy1
	adcs	r7, r7, r0

	mul	r11, r1, b1	C a1 b1
	umulh	r1, r1, b1

	mul	r0, r2, b1	C a2 b1
	umulh	r2, r2, b1

	mul	cy0, r3, b1	C a3 b1
	umulh	cy1, r3, b1
	mul	r12, r3, b0	C a3 b0
	umulh	r13, r3, b0

	adcs	r10, r10, r1
	adcs	cy0, cy0, r2
	adc	cy1, cy1, zr
	C 4, (5, 6, 8), (7, [0], 9, 11), (10, [1], 12, 0), (13, [2], cy0), cy1

	adds	r5, r5, r8
	adcs	r7, r7, r11
	adcs	r10, r10, r0
	adcs	cy0, cy0, r13
	adc	cy1, cy1, zr
	C 4, (5, 6), (7, 9), (10, 12), cy0, cy1

	adds	r5, r5, r6
	adcs	r7, r7, r9
	adcs	r10, r10, r12
	C 4, 5, 7, 10, cy0 (+ carry), cy1

	cbnz	n, L(top)

	adcs	cy0, cy0, zr
	adc	cy1, cy1, zr

	stp	r4, r5, [rp], #2*8
	stp	r7, r10, [rp], #2*8

	ldp	r11, r12, [sp,#-3*8]
	ldr	r13, [sp,#-1*8]

L(end):	stp	cy0, cy1, [rp]
	mov	x0, cy1
	ret
EPILOGUE()
