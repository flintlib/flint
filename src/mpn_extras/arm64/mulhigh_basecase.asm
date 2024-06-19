dnl
dnl Copyright (C) 2024 Albin Ahlbäck
dnl
dnl Uses code adopted from the GNU MP Library.
dnl
dnl     Copyright 2020 Free Software Foundation, Inc.
dnl     Contributed to the GNU project by Torbjörn Granlund.
dnl
dnl This file is part of FLINT.
dnl
dnl FLINT is free software: you can redistribute it and/or modify it under
dnl the terms of the GNU Lesser General Public License (LGPL) as published
dnl by the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.  See <https://www.gnu.org/licenses/>.
dnl

include(`config.m4')

define(`rp',  `x0')
define(`ap',  `x1')
define(`bp',  `x2')
define(`n',   `x3')

define(`r0',  `x4')
define(`r1',  `x5')
define(`r2',  `x6')
define(`r3',  `x7')
define(`r4',  `x8')
define(`r5',  `x9')
define(`r6',  `x10')
define(`r7',  `x11')
define(`r8',  `x12')
define(`r9',  `x13')
define(`r10', `x14')
define(`r11', `x15')
define(`r12', `x16')
define(`r13', `x17')

dnl NOTE: The following has to be pushed and popped
dnl NOTE: At least on Apple, we need to preserve both x18 and x29 at all times!
define(`r14', `x19')
define(`r15', `x20')
define(`r16', `x21')
define(`r17', `x22')
define(`r18', `x23')
define(`r19', `x24')
define(`r20', `x25')
define(`r21', `x26')
define(`r22', `x27')
define(`r23', `x28')
define(`r24', `x30')

define(`res', `rp') dnl NOTE: Synonymous with rp!

dnl Idea is to compute the first 9 by 8 block, then continue via addmul chains.

dnl n = 10:
dnl
dnl   0 1 2 3 4 5 6 7 8 9
dnl 0  _ _ _ _ _ _ _  h x
dnl 1               h¦x x
dnl 2             h x¦x x
dnl 3           h x x¦x x
dnl 4         h x x x¦x x
dnl 5       h x x x x¦x x
dnl 6     h x x x x x¦x x
dnl 7   h x x x x x x¦x x
dnl 8 h x x x x x x x¦x x
dnl 9 x x x x x x x x¦x x

PROLOGUE(_flint_mpn_mulhigh_basecase)
	add	ap, ap, n, lsl #3
	sub	ap, ap, #9*8
	stp	r14, r15, [sp, #-10*8]!
	ldp	r9, r10, [bp]
	ldp	r6, r7, [ap, #6*8]
	ldr	r8, [ap, #8*8]

	ldp	r11, r12, [bp, #2*8]
	stp	r16, r17, [sp, #2*8]
	ldp	r13, r14, [bp, #4*8]
	stp	r18, r19, [sp, #4*8]

	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]

	stp	r20, r21, [sp, #6*8]
	ldp	r15, r16, [bp, #6*8]
	stp	r22, r23, [sp, #8*8]

	C ap: 0, 1, ..., 8
	C bp: 9, ..., 16
	C Free: 17, ..., 23

	umulh	r17, r7, r9
	mul	r18, r8, r9
	umulh	r9, r8, r9
	umulh	r19, r7, r10
	umulh	r20, r8, r10
	umulh	r21, r6, r10
	mul	r22, r7, r10
	mul	r10, r8, r10
	C (17, 18, 21, 22), (9, 19, 10), 20
	C Free: 23

	umulh	r23, r5, r11
	adds	r17, r17, r18
	umulh	r18, r6, r11
	adcs	r9, r9, r19
	umulh	r19, r7, r11
	cinc	r20, r20, cs
	adds	r17, r17, r21
	umulh	r21, r8, r11
	adcs	r9, r9, r10
	mul	r10, r6, r11
	cinc	r20, r20, cs
	adds	r17, r17, r22
	mul	r22, r7, r11
	mul	r11, r8, r11

	adcs	r9, r9, r18
	umulh	r18, r4, r12
	adcs	r20, r20, r19
	umulh	r19, r5, r12
	cinc	r21, r21, cs
	adds	r17, r17, r23
	umulh	r23, r6, r12
	adcs	r9, r9, r22
	umulh	r22, r7, r12
	adcs	r20, r20, r11
	umulh	r11, r8, r12
	cinc	r21, r21, cs
	adds	r17, r17, r10
	mul	r10, r5, r12
	adcs	r9, r9, r19
	mul	r19, r6, r12
	adcs	r20, r20, r23
	mul	r23, r7, r12
	mul	r12, r8, r12

	adcs	r21, r21, r22
	umulh	r22, r3, r13
	cinc	r11, r11, cs
	adds	r17, r17, r18
	umulh	r18, r4, r13
	adcs	r9, r9, r19
	umulh	r19, r5, r13
	adcs	r20, r20, r23
	umulh	r23, r6, r13
	adcs	r21, r21, r12
	umulh	r12, r7, r13
	cinc	r11, r11, cs
	adds	r17, r17, r10
	umulh	r10, r8, r13
	adcs	r9, r9, r18
	mul	r18, r4, r13
	adcs	r20, r20, r19
	mul	r19, r5, r13
	adcs	r21, r21, r23
	mul	r23, r6, r13
	adcs	r11, r11, r12
	mul	r12, r7, r13
	mul	r13, r8, r13

	cinc	r10, r10, cs
	adds	r17, r17, r22
	umulh	r22, r2, r14
	adcs	r9, r9, r19
	umulh	r19, r3, r14
	adcs	r20, r20, r23
	umulh	r23, r4, r14
	adcs	r21, r21, r12
	umulh	r12, r5, r14
	adcs	r11, r11, r13
	umulh	r13, r6, r14
	cinc	r10, r10, cs
	adds	r17, r17, r18
	umulh	r18, r7, r14
	adcs	r9, r9, r19
	umulh	r19, r8, r14
	adcs	r20, r20, r23
	mul	r23, r3, r14
	adcs	r21, r21, r12
	mul	r12, r4, r14
	adcs	r11, r11, r13
	mul	r13, r5, r14
	adcs	r10, r10, r18
	mul	r18, r6, r14
	cinc	r19, r19, cs
	adds	r17, r17, r22
	mul	r22, r7, r14
	mul	r14, r8, r14

	adcs	r9, r9, r12
	umulh	r12, r1, r15
	adcs	r20, r20, r13
	umulh	r13, r2, r15
	adcs	r21, r21, r18
	umulh	r18, r3, r15
	adcs	r11, r11, r22
	umulh	r22, r4, r15
	adcs	r10, r10, r14
	umulh	r14, r5, r15
	cinc	r19, r19, cs
	adds	r17, r17, r23
	umulh	r23, r6, r15
	adcs	r9, r9, r13
	umulh	r13, r7, r15
	adcs	r20, r20, r18
	umulh	r18, r8, r15
	adcs	r21, r21, r22
	mul	r22, r2, r15
	adcs	r11, r11, r14
	mul	r14, r3, r15
	adcs	r10, r10, r23
	mul	r23, r4, r15
	adcs	r19, r19, r13
	mul	r13, r5, r15
	cinc	r18, r18, cs
	adds	r17, r17, r12
	mul	r12, r6, r15
	adcs	r9, r9, r14
	mul	r14, r7, r15
	mul	r15, r8, r15

	adcs	r20, r20, r23
	umulh	r23, r1, r16
	adcs	r21, r21, r13
	umulh	r13, r2, r16
	adcs	r11, r11, r12
	umulh	r12, r3, r16
	adcs	r10, r10, r14
	umulh	r14, r4, r16
	adcs	r19, r19, r15
	umulh	r15, r5, r16
	cinc	r18, r18, cs
	adds	r17, r17, r22
	umulh	r22, r6, r16
	adcs	r9, r9, r23
	umulh	r23, r7, r16
	adcs	r20, r20, r13
	umulh	r13, r8, r16
	mul	r1, r1, r16
	adcs	r21, r21, r12
	mul	r2, r2, r16
	adcs	r11, r11, r14
	mul	r3, r3, r16
	adcs	r10, r10, r15
	mul	r4, r4, r16
	adcs	r19, r19, r22
	mul	r5, r5, r16
	adcs	r18, r18, r23
	mul	r6, r6, r16
	cinc	r13, r13, cs
	mul	r7, r7, r16
	adds	r1, r17, r1	C Store lowest limb (return value) in r1
	mul	r8, r8, r16
	adcs	r9, r9, r2
	umulh	r2, r0, r16	C Save r2 for next addmul
	adcs	r20, r20, r3
	stp	r9, r20, [rp]
	adcs	r21, r21, r4
	adcs	r11, r11, r5
	stp	r21, r11, [rp, #2*8]
	adcs	r10, r10, r6
	adcs	r19, r19, r7
	stp	r10, r19, [rp, #4*8]
	adcs	r18, r18, r8
	cinc	r13, r13, cs
	stp	r18, r13, [rp, #6*8]

define(`b0',	`r0')	C loaded bp value
define(`rx',	`r1')	C return value
define(`cy',	`r2')	C carry for next iteration in loop

define(`ix',	`r3')	C iterator
define(`ixsave',`r4')	C saved iterator
define(`ns',    `r5')	C resetter for ap and rp

define(`jmp',	`r6')	C jump values
dnl 7 and 20, ..., 23 will remain unused for the rest of the program

L(hi):	LBL_HI(	jmp, L(p1h))
	mov	ixsave, #2
	mov	ix, #2
	ldp	r20, r21, [sp, #6*8]	C No longer in use
	ldp	r22, r23, [sp, #8*8]
L(lo):	LBL_LO(	jmp, L(p1h))
	LOH_ADRPADD(L(hi), L(lo))
	ldr	b0, [bp, #8*8]!
	mov	ns, #8*8
	sub	n, n, #8

L(p1h):	ldr	r8, [ap], #1*8
	mul	r12, r8, b0
	umulh	r8, r8, b0
	add	jmp, jmp, #9*4		C jmp goes from p1h to p2h
	adds	rx, rx, cy
	cinc	r8, r8, cs
	adds	rx, rx, r12
	cinc	cy, r8, cs
	b	L(top)

L(p2h):	ldp	r8, r9, [ap], #2*8
	ldr	r17, [rp]
	mul	r12, r8, b0
	umulh	r8, r8, b0
	mul	r13, r9, b0
	umulh	r9, r9, b0
	add	jmp, jmp, #15*4		C jmp goes from p2h to p3h
	adds	rx, rx, cy
	adcs	r8, r17, r8
	cinc	r9, r9, cs
	adds	rx, rx, r12
	adcs	r8, r8, r13
	cinc	cy, r9, cs
	str	r8, [rp], #1*8
	b	L(top)

L(p3h):	ldp	r8, r9, [ap], #2*8
	ldr	r10, [ap], #1*8
	ldp	r17, r18, [rp]
	mul	r12, r8, b0
	umulh	r8, r8, b0
	mul	r13, r9, b0
	umulh	r9, r9, b0
	mul	r14, r10, b0
	umulh	r10, r10, b0
	add	jmp, jmp, #20*4		C jmp goes from p3h to p0h
	adds	rx, rx, cy
	adcs	r8, r17, r8
	adcs	r9, r18, r9
	cinc	r10, r10, cs
	adds	rx, rx, r12
	adcs	r8, r8, r13
	adcs	r9, r9, r14
	cinc	cy, r10, cs
	stp	r8, r9, [rp], #2*8
	b	L(top)

L(p0h):	ldp	r8, r9, [ap], #2*8
	ldp	r10, r11, [ap], #2*8
	ldp	r17, r18, [rp]
	ldr	r19, [rp, #2*8]
	add	ixsave, ixsave, #1
	sub	jmp, jmp, #44*4		C jmp goes from p0h to p1h
	mul	r12, r8, b0
	umulh	r8, r8, b0
	mul	r13, r9, b0
	umulh	r9, r9, b0
	mul	r14, r10, b0
	umulh	r10, r10, b0
	mul	r15, r11, b0
	umulh	r11, r11, b0
	adds	rx, rx, cy
	adcs	r8, r17, r8
	adcs	r9, r18, r9
	adcs	r10, r19, r10
	cinc	r11, r11, cs
	adds	rx, rx, r12
	adcs	r8, r8, r13
	adcs	r9, r9, r14
	adcs	r10, r10, r15
	cinc	cy, r11, cs
	stp	r8, r9, [rp], #2*8
	str	r10, [rp], #1*8

L(top):	ldp	r8, r9, [ap], #2*8
	ldp	r10, r11, [ap], #2*8
	ldp	r16, r17, [rp]
	ldp	r18, r19, [rp, #2*8]
	mul	r12, r8, b0
	umulh	r8, r8, b0
	mul	r13, r9, b0
	umulh	r9, r9, b0
	mul	r14, r10, b0
	umulh	r10, r10, b0
	mul	r15, r11, b0
	umulh	r11, r11, b0
	adds	r12, r16, r12
	adcs	r8, r17, r8
	adcs	r9, r18, r9
	adcs	r10, r19, r10
	cinc	r11, r11, cs
	adds	r12, r12, cy
	adcs	r8, r8, r13
	adcs	r9, r9, r14
	adcs	r10, r10, r15
	cinc	cy, r11, cs
	stp	r12, r8, [rp], #2*8
	stp	r9, r10, [rp], #2*8
	sub	ix, ix, #1
	cbnz	ix, L(top)

L(end):	str	cy, [rp]

	sub	n, n, #1
	cbz	n, L(fin)

	sub	rp, rp, ns		C Reset rp
	sub	ap, ap, ns		C Reset ap, part 1
	add	ns, ns, #1*8		C Increase reseter

	ldr	cy, [ap, #-2*8]!	C Reset ap, part 2
	umulh	cy, cy, b0

	mov	ix, ixsave
	ldr	b0, [bp, #1*8]!
	br	jmp

L(fin):	ldp	r16, r17, [sp, #2*8]
	ldp	r18, r19, [sp, #4*8]
	ldp	r14, r15, [sp], #10*8

	mov	res, rx

	ret
EPILOGUE()
