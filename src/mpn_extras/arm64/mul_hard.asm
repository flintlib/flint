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

define(`zr',  `xzr')

define(`res', `rp') dnl NOTE: Synonymous with rp!

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

PROLOGUE(flint_mpn_mul_1n)
	ldr	r0, [ap]
	ldr	r1, [bp]

	mul	r3, r0, r1
	umulh	r4, r0, r1

	stp	r3, r4, [rp]
	mov	res, r4

	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_2n)
	cmp	n, #2
	ldp	r0, r1, [ap]
	b.eq	L(m22)

L(m21):	ldr	r2, [bp]

	mul	r4, r0, r2	C a0 b0
	umulh	r5, r0, r2
	mul	r6, r1, r2	C a1 b0
	umulh	r7, r1, r2

	str	r4, [rp]
	adds	r5, r5, r6
	cinc	r7, r7, cs
	stp	r5, r7, [rp,#1*8]
	mov	res, r7
	ret

L(m22):	ldp	r2, r3, [bp]

	mul	r4, r0, r2	C a0 b0
	umulh	r5, r0, r2
	mul	r6, r1, r2	C a1 b0
	umulh	r7, r1, r2

	mul	r8, r0, r3	C a0 b1
	umulh	r9, r0, r3
	mul	r10, r1, r3	C a1 b1
	umulh	r11, r1, r3

	adds	r5, r5, r6
	cinc	r7, r7, cs
	C 5, 7

	adds	r9, r9, r10
	cinc	r11, r11, cs
	C 8, 9, 11

	adds	r5, r5, r8
	adcs	r7, r7, r9
	cinc	r11, r11, cs
	C 5, 7, 11

	stp	r4, r5, [rp]
	stp	r7, r11, [rp,#2*8]
	mov	res, r11
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_3n)
	ldr	r3, [bp], #1*8
	ldp	r0, r1, [ap]
	ldr	r2, [ap,#2*8]

L(m3):	mul	r4, r0, r3	C a0 b0
	umulh	r5, r0, r3
	mul	r6, r1, r3	C a1 b0
	umulh	r7, r1, r3
	mul	r8, r2, r3	C a2 b0
	umulh	r9, r2, r3

	str	r4, [rp], #1*8
	sub	n, n, #1
	adds	r6, r6, r5
	adcs	r8, r8, r7
	cinc	r9, r9, cs
	C 6, 8, 9

	cbz	n, L(f3)

	tst	n, #1
	b.eq	L(p3)

L(a3):	ldr	r3, [bp], #1*8
	mul	r4, r0, r3	C a0 b0
	umulh	r5, r0, r3
	mul	r7, r1, r3	C a1 b0
	umulh	r10, r1, r3
	mul	r11, r2, r3	C a2 b0
	umulh	r12, r2, r3

	sub	n, n, #1
	adds	r7, r7, r5
	adcs	r11, r11, r10
	cinc	r12, r12, cs
	C 4, 7, 11, 12

	C 4,  (6, 8,  9)
	C  <- (6, 8,  9)
	C   + (4, 7, 11, 12)
	adds	r4, r6, r4
	adcs	r6, r8, r7
	adcs	r8, r9, r11
	cinc	r9, r12, cs
	str	r4, [rp], #1*8

	cbz	n, L(f3)

L(p3):	lsr	n, n, #1
L(b3):	ldp	r3, r4, [bp], #2*8
	sub	n, n, #1

	mul	r5, r0, r3	C a0 b0
	umulh	r7, r0, r3
	mul	r10, r1, r3	C a1 b0
	umulh	r11, r1, r3
	mul	r12, r2, r3	C a2 b0
	umulh	r3, r2, r3

	mul	r13, r0, r4	C a0 b1
	umulh	ap, r0, r4

	adds	r5, r5, r6
	adcs	r7, r7, r8
	adcs	r11, r11, r9

	mul	r6, r1, r4	C a1 b1
	umulh	r8, r1, r4
	umulh	r9, r2, r4	C a2 b1
	mul	r4, r2, r4

	cinc	r3, r3, cs
	adds	r10, r10, r7
	adcs	r12, r12, r11
	cinc	r3, r3, cs	C Carry cannot happen here!
	C 5, 10, 12, 3

	adds	r6, r6, ap
	adcs	r4, r4, r8
	cinc	r9, r9, cs
	C 13, 6, 4, 9

	C   (5, 10, 12, 3)
	C    + (13,  6, 4, 9)
	C -> 5, 13, (6, 8, 9)
	adds	r13, r13, r10
	adcs	r6, r6, r12
	adcs	r8, r4, r3
	cinc	r9, r9, cs

	stp	r5, r13, [rp], #2*8

	cbnz	n, L(b3)

L(f3):	stp	r6, r8, [rp]
	str	r9, [rp,#2*8]
	mov	res, r9
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_4n)
	ldr	r4, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]

L(m4):	mul	r5, r0, r4	C a0 b0
	umulh	r6, r0, r4
	mul	r7, r1, r4	C a1 b0
	umulh	r8, r1, r4
	mul	r9, r2, r4	C a2 b0
	umulh	r10, r2, r4
	mul	r11, r3, r4	C a3 b0
	umulh	r12, r3, r4
	C 5, (6, 7), (8, 9), (10, 11), 12

	str	r5, [rp], #1*8
	adds	r6, r6, r7
	sub	n, n, #1
	adcs	r8, r8, r9
	adcs	r10, r10, r11
	cinc	r12, r12, cs
	C 6, 8, 10, 12

	cbz	n, L(f4)

	stp	r14, r15, [sp, #-2*8]!
L(a4):	ldr	r4, [bp], #1*8
	mul	r5, r0, r4	C a0 b0
	umulh	r7, r0, r4
	mul	r9, r1, r4	C a1 b0
	umulh	r11, r1, r4
	mul	r13, r2, r4	C a2 b0
	umulh	r14, r2, r4
	mul	r15, r3, r4	C a3 b0
	umulh	ap, r3, r4
	C 5, (7, 9), (11, 13), (14, 15), ap

	adds	r7, r7, r9
	adcs	r11, r11, r13
	adcs	r14, r14, r15
	cinc	ap, ap, cs
	C 5, 7, 11, 14, 4

	C (5, (6, 8, 10, 12)) <- (6, 8, 10, 12) + (5, 7, 11, 14, ap)
	adds	r5, r6, r5
	adcs	r6, r8, r7
	sub	n, n, #1
	adcs	r8, r10, r11
	adcs	r10, r12, r14
	cinc	r12, ap, cs

	str	r5, [rp], #1*8

	cbnz	n, L(a4)
	ldp	r14, r15, [sp], #2*8

L(f4):	stp	r6, r8, [rp]
	stp	r10, r12, [rp,#2*8]
	mov	res, r12
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_5n)
	ldr	r5, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldr	r4, [ap,#4*8]

L(m5):	mul	r6, r0, r5	C a0 b0
	umulh	r7, r0, r5
	mul	ap, r1, r5	C a1 b0
	umulh	r8, r1, r5
	mul	r9, r2, r5	C a2 b0
	umulh	r10, r2, r5
	mul	r11, r3, r5	C a3 b0
	str	r6, [rp], #1*8
	umulh	r12, r3, r5
	mul	r6, r4, r5	C a4 b0
	umulh	r5, r4, r5
	C (7, ap), (8, 9), (10, 11), (12, 6), 5

	adds	r7, r7, ap
	adcs	r8, r8, r9
	adcs	r10, r10, r11
	sub	n, n, #1
	adcs	r12, r12, r6
	cinc	r9, r5, cs
	C 7, 8, 10, 12, 9

	cbz	n, L(f5)

	stp	r14, r15, [sp, #-4*8]!
	stp	r16, r17, [sp, #2*8]
L(a5):	ldr	r5, [bp], #1*8
	mul	r6, r0, r5	C a0 b0
	umulh	r11, r0, r5
	mul	ap, r1, r5	C a1 b0
	umulh	r13, r1, r5
	mul	r14, r2, r5	C a2 b0
	umulh	r15, r2, r5
	mul	r16, r3, r5	C a3 b0
	umulh	r17, r3, r5
	adds	r11, r11, ap
	adcs	r13, r13, r14
	mul	ap, r4, r5	C a4 b0
	umulh	r5, r4, r5
	sub	n, n, #1
	adcs	r15, r15, r16
	adcs	r17, r17, ap
	cinc	r5, r5, cs
	C 6, 11, 13, 15, 17, 5

	C (6, (7, 8, 10, 12, 9)) <- (7, 8, 10, 12, 9) + (6, 11, 13, 15, 17, 5)
	adds	r6, r7, r6
	adcs	r7, r8, r11
	adcs	r8, r10, r13
	adcs	r10, r12, r15
	adcs	r12, r9, r17
	cinc	r9, r5, cs

	str	r6, [rp], #1*8

	cbnz	n, L(a5)
	ldp	r16, r17, [sp, #2*8]
	ldp	r14, r15, [sp], #4*8

L(f5):	stp	r7, r8, [rp]
	stp	r10, r12, [rp,#2*8]
	str	r9, [rp,#4*8]
	mov	res, r9
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_6n)
	ldr	r6, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]

L(m6):	mul	r7, r0, r6	C a0 b0
	umulh	r8, r0, r6
	mul	ap, r1, r6	C a1 b0
	umulh	r9, r1, r6
	mul	r10, r2, r6	C a2 b0
	umulh	r11, r2, r6
	mul	r12, r3, r6	C a3 b0
	umulh	r13, r3, r6
	str	r7, [rp], #1*8
	sub	n, n, #1
	adds	r8, r8, ap
	adcs	r9, r9, r10
	mul	r7, r4, r6	C a4 b0
	umulh	ap, r4, r6
	mul	r10, r5, r6	C a5 b0
	umulh	r6, r5, r6
	adcs	r11, r11, r12
	adcs	r13, r13, r7
	adcs	r10, ap, r10
	cinc	r7, r6, cs
	C 8, 9, 11, 13, 10, 7

	cbz	n, L(f6)

	stp	r14, r15, [sp, #-6*8]!
	stp	r16, r17, [sp, #2*8]
	stp	r18, r19, [sp, #4*8]
L(a6):	ldr	r6, [bp], #1*8
	mul	r12, r0, r6	C a0 b0
	umulh	r14, r0, r6
	mul	r15, r1, r6	C a1 b0
	umulh	r16, r1, r6
	mul	r17, r2, r6	C a2 b0
	umulh	r18, r2, r6
	mul	r19, r3, r6	C a3 b0
	umulh	ap, r3, r6
	adds	r14, r14, r15
	adcs	r16, r16, r17
	adcs	r18, r18, r19
	mul	r15, r4, r6	C a4 b0
	umulh	r17, r4, r6
	mul	r19, r5, r6	C a5 b0
	umulh	r6, r5, r6
	sub	n, n, #1
	adcs	ap, ap, r15
	adcs	r17, r17, r19
	cinc	r6, r6, cs
	C 12, 14, 16, 18, ap, 17, 6

	C (12, (8, 9, 11, 13, 10, 7))
	C <- (8, 9, 11, 13, 10, 7) + (12, 14, 16, 18, ap, 17, 6)
	adds	r12, r8, r12
	adcs	r8, r9, r14
	adcs	r9, r11, r16
	adcs	r11, r13, r18
	adcs	r13, r10, ap
	adcs	r10, r7, r17
	cinc	r7, r6, cs

	str	r12, [rp], #1*8

	cbnz	n, L(a6)
	ldp	r16, r17, [sp, #2*8]
	ldp	r18, r19, [sp, #4*8]
	ldp	r14, r15, [sp], #6*8

L(f6):	stp	r8, r9, [rp]
	stp	r11, r13, [rp,#2*8]
	stp	r10, r7, [rp,#4*8]
	mov	res, r7
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_7n)
	ldr	r7, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]
	ldr	r6, [ap,#6*8]

	stp	r14, r15, [sp, #-2*8]!

L(m7):	mul	r8, r0, r7	C a0 b0
	umulh	r9, r0, r7
	mul	ap, r1, r7	C a1 b0
	umulh	r14, r1, r7
	mul	r10, r2, r7	C a2 b0
	umulh	r15, r2, r7
	mul	r11, r3, r7	C a3 b0
	umulh	r12, r3, r7
	mul	r13, r4, r7	C a4 b0
	str	r8, [rp], #1*8
	sub	n, n, #1
	adds	r9, r9, ap
	adcs	r14, r14, r10
	adcs	r15, r15, r11
	umulh	r8, r4, r7
	mul	ap, r5, r7	C a5 b0
	umulh	r10, r5, r7
	mul	r11, r6, r7	C a6 b0
	umulh	r7, r6, r7
	adcs	r12, r12, r13
	adcs	r8, r8, ap
	adcs	r10, r10, r11
	cinc	r13, r7, cs
	C 9, 14, 15, 12, 8, 10, 13

	cbz	n, L(f7)

	stp	r16, r17, [sp, #-6*8]!
	stp	r18, r19, [sp, #2*8]
	stp	r20, r21, [sp, #4*8]
L(a7):	ldr	r7, [bp], #1*8
	mul	r11, r0, r7	C a0 b0
	umulh	ap, r0, r7
	mul	r16, r1, r7	C a1 b0
	umulh	r17, r1, r7
	mul	r18, r2, r7	C a2 b0
	umulh	r19, r2, r7
	mul	r20, r3, r7	C a3 b0
	umulh	r21, r3, r7
	adds	ap, ap, r16
	adcs	r17, r17, r18
	adcs	r19, r19, r20
	mul	r16, r4, r7	C a4 b0
	umulh	r18, r4, r7
	mul	r20, r5, r7	C a5 b0
	adcs	r21, r21, r16
	adcs	r18, r18, r20
	umulh	r16, r5, r7
	mul	r20, r6, r7	C a6 b0
	umulh	r7, r6, r7
	sub	n, n, #1
	adcs	r16, r16, r20
	cinc	r7, r7, cs
	C 11, ap, 17, 19, 21, 18, 16, 7

	C (11, (9, 14, 15, 12, 8, 10, 13))
	C <- (9, 14, 15, 12, 8, 10, 13) + (11, ap, 17, 19, 21, 18, 16, 7)
	adds	r11, r9, r11
	adcs	r9, r14, ap
	adcs	r14, r15, r17
	adcs	r15, r12, r19
	adcs	r12, r8, r21
	adcs	r8, r10, r18
	adcs	r10, r13, r16
	cinc	r13, r7, cs

	str	r11, [rp], #1*8

	cbnz	n, L(a7)
	ldp	r18, r19, [sp, #2*8]
	ldp	r20, r21, [sp, #4*8]
	ldp	r16, r17, [sp], #6*8

L(f7):	stp	r9, r14, [rp]
	stp	r15, r12, [rp,#2*8]
	stp	r8, r10, [rp,#4*8]
	str	r13, [rp,#6*8]
	ldp	r14, r15, [sp], #2*8
	mov	res, r13
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_8n)
	ldr	r8, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	stp	r14, r15, [sp, #-6*8]!
	ldp	r4, r5, [ap, #4*8]
	ldp	r6, r7, [ap, #6*8]
	stp	r16, r17, [sp, #2*8]
	stp	r18, r19, [sp, #4*8]

L(m8):	mul	r9, r0, r8
	umulh	r10, r0, r8
	mul	r11, r1, r8
	umulh	r12, r1, r8
	mul	r13, r2, r8
	umulh	r14, r2, r8
	mul	r15, r3, r8
	umulh	r16, r3, r8
	mul	r17, r4, r8
	umulh	r18, r4, r8
	mul	r19, r5, r8
	umulh	ap, r5, r8
	C 9, (10, 11), (12, 13), (14, 15), (16, 17), (18, 19), ap

	str	r9, [rp], #1*8
	sub	n, n, #1
	adds	r10, r10, r11

	mul	r9, r6, r8
	umulh	r11, r6, r8

	adcs	r12, r12, r13
	adcs	r14, r14, r15

	mul	r13, r7, r8
	umulh	r15, r7, r8
	C 10, 12, 14, (16, 17), (18, 19), (ap, 9), (11, 13), 15

	cbz	n, L(f8)

	stp	r20, r21, [sp, #-6*8]!
	stp	r22, r23, [sp, #2*8]
	str	r24, [sp, #4*8]

L(a8):	ldr	r8, [bp], #1*8
	adcs	r16, r16, r17
	adcs	r18, r18, r19
	adcs	ap, ap, r9
	adcs	r11, r11, r13
	cinc	r15, r15, cs
	C 10, 12, 14, 16, 18, ap, 11, 15
	C Free: 20, 21, 22, 23, 24, 17, 19, 9, 13

	mul	r20, r0, r8
	mul	r21, r1, r8
	mul	r22, r2, r8
	mul	r23, r3, r8
	mul	r17, r4, r8
	mul	r19, r5, r8
	mul	r9, r6, r8
	mul	r13, r7, r8

	adds	r20, r20, r10
	umulh	r24, r7, r8

	C (10, 20), (12, 21), (14, 22), (16, 23), (18, 17), (ap, 19), (11, 9), (15, 13), 24

	adcs	r21, r21, r12
	umulh	r10, r0, r8

	adcs	r22, r22, r14
	umulh	r12, r1, r8

	adcs	r23, r23, r16
	umulh	r14, r2, r8

	adcs	r17, r17, r18
	umulh	r16, r3, r8

	adcs	r19, r19, ap
	umulh	r18, r4, r8

	adcs	r9, r9, r11
	umulh	ap, r5, r8

	adcs	r13, r13, r15
	umulh	r11, r6, r8

	str	r20, [rp], #1*8

	cinc	r15, r24, cs

	adds	r10, r10, r21
	adcs	r12, r12, r22
	adcs	r14, r14, r23
	C 10, 12, 14, (16, 17), (18, 19), (ap, 9), (11, 13), 15

	sub	n, n, #1
	cbnz	n, L(a8)

	ldp	r22, r23, [sp, #2*8]
	ldr	r24, [sp, #4*8]
	ldp	r20, r21, [sp], #6*8
L(f8):	adcs	r16, r16, r17
	adcs	r18, r18, r19
	stp	r10, r12, [rp]
	adcs	ap, ap, r9
	adcs	r11, r11, r13
	stp	r14, r16, [rp, #2*8]
	cinc	r15, r15, cs
	ldp	r16, r17, [sp, #2*8]
	stp	r18, ap, [rp, #4*8]
	ldp	r18, r19, [sp, #4*8]
	stp	r11, r15, [rp, #6*8]
	mov	res, r15
	ldp	r14, r15, [sp], #6*8

	ret
EPILOGUE()

dnl Before we loaded everything from ap into r0, ..., r{n - 1}. However, we are
dnl now running out of registers. Instead, we will now load parts of ap
dnl continuously, and at the same time use those registers. Note that this will
dnl give worse performance for larger n.

PROLOGUE(flint_mpn_mul_9n)
	ldr	r4, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]

L(m9):	mul	r5, r0, r4	C a0 b0
	umulh	r6, r0, r4
	mul	r7, r1, r4	C a1 b0
	umulh	r8, r1, r4
	mul	r9, r2, r4	C a2 b0
	umulh	r10, r2, r4
	mul	r11, r3, r4	C a3 b0
	umulh	r12, r3, r4
	str	r5, [rp], #1*8
	ldp	r0, r1, [ap,#4*8]
	ldp	r2, r3, [ap,#6*8]
	adds	r6, r6, r7
	adcs	r8, r8, r9
	adcs	r10, r10, r11
	mul	r13, r0, r4	C a4 b0
	umulh	r5, r0, r4
	mul	r7, r1, r4	C a5 b0
	umulh	r1, r1, r4
	mul	r9, r2, r4	C a6 b0
	umulh	r2, r2, r4
	mul	r11, r3, r4	C a7 b0
	umulh	r3, r3, r4
	ldr	r0, [ap,#8*8]
	adcs	r12, r12, r13
	adcs	r5, r5, r7
	adcs	r9, r9, r1
	adcs	r11, r11, r2
	mul	r13, r0, r4	C a8 b0
	umulh	r0, r0, r4
	sub	n, n, #1
	adcs	r13, r13, r3
	cinc	r7, r0, cs
	C 6, 8, 10, 12, 5, 9, 11, 13, 7

	cbz	n, L(f9)

	stp	r14, r15, [sp, #-8*8]!
	stp	r16, r17, [sp, #2*8]
	stp	r18, r19, [sp, #4*8]
	stp	r20, r21, [sp, #6*8]
	ldp	r0, r1, [ap]
L(a9):	ldr	r4, [bp], #1*8
	ldp	r2, r3, [ap,#2*8]
	mul	r14, r0, r4	C a0 b0
	umulh	r15, r0, r4
	mul	r16, r1, r4	C a1 b0
	umulh	r17, r1, r4
	mul	r18, r2, r4	C a2 b0
	umulh	r2, r2, r4
	mul	r19, r3, r4	C a3 b0

	ldp	r0, r1, [ap,#4*8]

	umulh	r3, r3, r4
	mul	r20, r0, r4	C a4 b0
	umulh	r0, r0, r4

	adds	r15, r15, r16
	adcs	r17, r17, r18
	adcs	r19, r19, r2
	adcs	r20, r20, r3
	ldp	r2, r3, [ap,#6*8]

	mul	r16, r1, r4	C a5 b0
	umulh	r1, r1, r4
	mul	r18, r2, r4	C a6 b0
	umulh	r2, r2, r4

	adcs	r16, r16, r0
	adcs	r18, r18, r1
	ldr	r0, [ap,#8*8]

	mul	r21, r3, r4	C a7 b0
	umulh	r3, r3, r4
	mul	r1, r0, r4	C a8 b0
	umulh	r4, r0, r4

	sub	n, n, #1
	adcs	r21, r21, r2
	adcs	r3, r3, r1
	ldp	r0, r1, [ap]
	cinc	r4, r4, cs
	C 14, 15, 17, 19, 20, 16, 18, 21, 3, 4

	C 14, ( 6,  8, 10, 12,  5,  9, 11, 13, 7)
	C  <- ( 6,  8, 10, 12,  5,  9, 11, 13, 7)
	C   + (14, 15, 17, 19, 20, 16, 18, 21, 3, 4)
	adds	r14, r6, r14
	adcs	r6, r8, r15
	adcs	r8, r10, r17
	adcs	r10, r12, r19
	adcs	r12, r5, r20
	adcs	r5, r9, r16
	adcs	r9, r11, r18
	adcs	r11, r13, r21
	adcs	r13, r7, r3
	cinc	r7, r4, cs

	str	r14, [rp], #1*8

	cbnz	n, L(a9)
	ldp	r16, r17, [sp, #2*8]
	ldp	r18, r19, [sp, #4*8]
	ldp	r20, r21, [sp, #6*8]
	ldp	r14, r15, [sp], #8*8

L(f9):	stp	r6, r8, [rp]
	stp	r10, r12, [rp,#2*8]
	stp	r5, r9, [rp,#4*8]
	stp	r11, r13, [rp,#6*8]
	str	r7, [rp,#8*8]
	mov	res, r7
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_10n)
	ldr	r6, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]

L(m10):	mul	r7, r0, r6	C a0 b0
	umulh	r0, r0, r6
	mul	r8, r1, r6	C a1 b0
	umulh	r1, r1, r6
	mul	r9, r2, r6	C a2 b0
	umulh	r2, r2, r6
	mul	r10, r3, r6	C a3 b0
	umulh	r3, r3, r6
	mul	r11, r4, r6	C a4 b0
	umulh	r4, r4, r6
	mul	r12, r5, r6	C a5 b0
	umulh	r5, r5, r6

	str	r7, [rp], #1*8
	adds	r8, r8, r0
	adcs	r9, r9, r1
	ldp	r0, r1, [ap,#6*8]
	adcs	r10, r10, r2
	adcs	r11, r11, r3
	adcs	r12, r12, r4
	ldp	r2, r3, [ap,#8*8]

	mul	r13, r0, r6	C a6 b0
	umulh	r0, r0, r6
	mul	r7, r1, r6	C a7 b0
	umulh	r1, r1, r6
	mul	r4, r2, r6	C a8 b0

	adcs	r13, r13, r5

	umulh	r2, r2, r6
	mul	r5, r3, r6	C a9 b0
	umulh	r3, r3, r6

	sub	n, n, #1
	adcs	r7, r7, r0
	adcs	r4, r4, r1
	adcs	r5, r5, r2
	cinc	r3, r3, cs
	C 8, 9, 10, 11, 12, 13, 7, 4, 5, 3
	C 0, 1, 2, 6, 14, ...

	cbz	n, L(f10)

	stp	r14, r15, [sp, #-10*8]!
	stp	r16, r17, [sp, #2*8]
	stp	r18, r19, [sp, #4*8]
	stp	r20, r21, [sp, #6*8]
	stp	r22, r23, [sp, #8*8]
L(a10):	ldp	r0, r1, [ap]
	ldr	r14, [bp], #1*8
	ldp	r2, r6, [ap,#2*8]
	ldp	r15, r16, [ap,#4*8]

	mul	r17, r0, r14	C a0 b0
	umulh	r0, r0, r14
	mul	r18, r1, r14	C a1 b0
	umulh	r1, r1, r14
	mul	r19, r2, r14	C a2 b0
	umulh	r2, r2, r14
	mul	r20, r6, r14	C a3 b0
	umulh	r6, r6, r14
	mul	r21, r15, r14	C a4 b0
	umulh	r15, r15, r14
	mul	r22, r16, r14	C a5 b0
	umulh	r16, r16, r14

	adds	r18, r18, r0
	adcs	r19, r19, r1
	ldp	r0, r1, [ap,#6*8]
	adcs	r20, r20, r2
	adcs	r21, r21, r6
	adcs	r22, r22, r15
	ldp	r2, r6, [ap,#8*8]

	mul	r23, r0, r14	C a6 b0
	umulh	r0, r0, r14
	mul	r15, r1, r14	C a7 b0
	umulh	r1, r1, r14

	adcs	r23, r23, r16
	adcs	r15, r15, r0

	mul	r16, r2, r14	C a8 b0
	umulh	r2, r2, r14
	mul	r0, r6, r14	C a9 b0
	umulh	r6, r6, r14

	sub	n, n, #1
	adcs	r16, r16, r1
	adcs	r0, r0, r2
	cinc	r6, r6, cs
	C 17, 18, 19, 20, 21, 22, 23, 15, 16, 0, 6

	C 17, ( 8,  9, 10, 11, 12, 13,  7,  4,  5, 3)
	C  <- ( 8,  9, 10, 11, 12, 13,  7,  4,  5, 3)
	C   + (17, 18, 19, 20, 21, 22, 23, 15, 16, 0, 6)
	adds	r17, r8, r17
	adcs	r8, r9, r18
	adcs	r9, r10, r19
	adcs	r10, r11, r20
	adcs	r11, r12, r21
	adcs	r12, r13, r22
	adcs	r13, r7, r23
	adcs	r7, r4, r15
	adcs	r4, r5, r16
	adcs	r5, r3, r0
	cinc	r3, r6, cs

	str	r17, [rp], #1*8

	cbnz	n, L(a10)

	ldp	r16, r17, [sp, #2*8]
	ldp	r18, r19, [sp, #4*8]
	ldp	r20, r21, [sp, #6*8]
	ldp	r22, r23, [sp, #8*8]
	ldp	r14, r15, [sp], #10*8

L(f10):	stp	r8, r9, [rp]
	stp	r10, r11, [rp,#2*8]
	stp	r12, r13, [rp,#4*8]
	stp	r7, r4, [rp,#6*8]
	stp	r5, r3, [rp,#8*8]
	mov	res, r3
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_11n)
	ldr	r6, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]

	stp	r14, r15, [sp, #-2*8]!
L(m11):	mul	r7, r0, r6	C a0 b0
	umulh	r0, r0, r6
	mul	r8, r1, r6	C a1 b0
	umulh	r1, r1, r6
	mul	r9, r2, r6	C a2 b0
	umulh	r2, r2, r6
	mul	r10, r3, r6	C a3 b0
	umulh	r3, r3, r6
	mul	r11, r4, r6	C a4 b0
	umulh	r4, r4, r6
	mul	r12, r5, r6	C a5 b0
	umulh	r5, r5, r6

	str	r7, [rp], #1*8
	adds	r8, r8, r0
	adcs	r9, r9, r1
	ldp	r0, r1, [ap,#6*8]
	adcs	r10, r10, r2
	adcs	r11, r11, r3
	adcs	r12, r12, r4
	ldp	r2, r3, [ap,#8*8]
	ldr	r4, [ap,#10*8]

	mul	r13, r0, r6	C a6 b0
	umulh	r0, r0, r6
	mul	r14, r1, r6	C a7 b0
	umulh	r1, r1, r6
	mul	r15, r2, r6	C a8 b0
	umulh	r2, r2, r6
	mul	r7, r3, r6	C a9 b0
	umulh	r3, r3, r6

	sub	n, n, #1
	adcs	r13, r13, r5
	adcs	r14, r14, r0

	mul	r5, r4, r6	C a10 b0
	umulh	r4, r4, r6

	adcs	r15, r15, r1
	adcs	r7, r7, r2
	adcs	r5, r5, r3
	cinc	r4, r4, cs
	C 8, 9, 10, 11, 12, 13, 14, 15, 7, 5, 4
	C 0, 1, 2, 3, 6, 16, 17, ...

	cbz	n, L(f11)

	stp	r16, r17, [sp, #-10*8]!
	stp	r18, r19, [sp, #2*8]
	stp	r20, r21, [sp, #4*8]
	stp	r22, r23, [sp, #6*8]
	str	r24, [sp, #8*8]
L(a11):	ldp	r0, r1, [ap]
	ldr	r16, [bp], #1*8
	ldp	r2, r3, [ap,#2*8]
	ldp	r6, r17, [ap,#4*8]

	mul	r18, r0, r16	C a0 b0
	umulh	r0, r0, r16
	mul	r19, r1, r16	C a1 b0
	umulh	r1, r1, r16
	mul	r20, r2, r16	C a2 b0
	umulh	r2, r2, r16
	mul	r21, r3, r16	C a3 b0
	umulh	r3, r3, r16
	mul	r22, r6, r16	C a4 b0
	umulh	r6, r6, r16
	mul	r23, r17, r16	C a5 b0
	umulh	r17, r17, r16

	adds	r19, r19, r0
	adcs	r20, r20, r1
	ldp	r0, r1, [ap,#6*8]
	adcs	r21, r21, r2
	adcs	r22, r22, r3
	adcs	r23, r23, r6
	ldp	r2, r3, [ap,#8*8]

	mul	r24, r0, r16	C a6 b0
	umulh	r0, r0, r16
	mul	r6, r1, r16	C a7 b0
	umulh	r1, r1, r16

	sub	n, n, #1
	adcs	r24, r24, r17
	adcs	r0, r0, r6

	mul	r17, r2, r16	C a8 b0
	umulh	r2, r2, r16
	mul	r6, r3, r16	C a9 b0
	umulh	r3, r3, r16

	adcs	r17, r17, r1
	adcs	r2, r2, r6
	ldr	r6, [ap,#10*8]

	mul	r1, r6, r16	C a10 b0
	umulh	r6, r6, r16

	adcs	r1, r1, r3
	cinc	r6, r6, cs
	C 18, 19, 20, 21, 22, 23, 24, 0, 17, 2, 1, 6

	C 18, ( 8,  9, 10, 11, 12, 13, 14, 15,  7, 5, 4)
	C  <- ( 8,  9, 10, 11, 12, 13, 14, 15,  7, 5, 4)
	C   + (18, 19, 20, 21, 22, 23, 24,  0, 17, 2, 1, 6)
	adds	r18, r8, r18
	adcs	r8, r9, r19
	adcs	r9, r10, r20
	adcs	r10, r11, r21
	adcs	r11, r12, r22
	adcs	r12, r13, r23
	adcs	r13, r14, r24
	adcs	r14, r15, r0
	adcs	r15, r7, r17
	adcs	r7, r5, r2
	adcs	r5, r4, r1
	cinc	r4, r6, cs

	str	r18, [rp], #1*8

	cbnz	n, L(a11)

	ldp	r18, r19, [sp, #2*8]
	ldp	r20, r21, [sp, #4*8]
	ldp	r22, r23, [sp, #6*8]
	ldr	r24, [sp, #8*8]
	ldp	r16, r17, [sp], #10*8

L(f11):	stp	r8, r9, [rp]
	stp	r10, r11, [rp,#2*8]
	stp	r12, r13, [rp,#4*8]
	stp	r14, r15, [rp,#6*8]
	stp	r7, r5, [rp,#8*8]
	str	r4, [rp,#10*8]
	ldp	r14, r15, [sp], #2*8
	mov	res, r4
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_12n)
	ldr	r8, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]
	ldp	r6, r7, [ap,#6*8]

	stp	r14, r15, [sp, #-4*8]!
	stp	r16, r17, [sp, #2*8]
L(m12):	mul	r9, r0, r8	C a0 b0
	umulh	r0, r0, r8
	mul	r10, r1, r8	C a1 b0
	umulh	r1, r1, r8
	mul	r11, r2, r8	C a2 b0
	umulh	r2, r2, r8
	mul	r12, r3, r8	C a3 b0
	umulh	r3, r3, r8
	mul	r13, r4, r8	C a4 b0
	umulh	r4, r4, r8
	mul	r14, r5, r8	C a5 b0
	umulh	r5, r5, r8
	mul	r15, r6, r8	C a6 b0
	umulh	r6, r6, r8
	mul	r16, r7, r8	C a7 b0
	umulh	r7, r7, r8

	str	r9, [rp], #1*8
	adds	r10, r10, r0
	adcs	r11, r11, r1
	ldp	r0, r1, [ap,#8*8]
	adcs	r12, r12, r2
	adcs	r13, r13, r3
	adcs	r14, r14, r4
	adcs	r15, r15, r5
	ldp	r2, r3, [ap,#10*8]

	mul	r17, r0, r8	C a8 b0
	umulh	r0, r0, r8
	mul	r9, r1, r8	C a9 b0
	umulh	r1, r1, r8
	mul	r4, r2, r8	C a10 b0
	umulh	r2, r2, r8
	mul	r5, r3, r8	C a11 b0
	umulh	r3, r3, r8

	adcs	r16, r16, r6
	adcs	r17, r17, r7
	sub	n, n, #1
	adcs	r9, r9, r0
	adcs	r4, r4, r1
	adcs	r5, r5, r2
	cinc	r3, r3, cs
	C 10, 11, 12, 13, 14, 15, 16, 17, 9, 4, 5, 3
	C 0, 1, 2, 6, 7, 8, 18, 19, ...

	cbz	n, L(f12)

	stp	r18, r19, [sp, #-12*8]!
	stp	r20, r21, [sp, #2*8]
	stp	r22, r23, [sp, #4*8]
	str	r24, [sp, #6*8]
L(a12):	ldp	r0, r1, [ap]
	ldr	r18, [bp], #1*8
	ldp	r2, r6, [ap,#2*8]
	ldp	r7, r8, [ap,#4*8]

	stp	r9, r4, [sp, #8*8]
	stp	r5, r3, [sp, #10*8]

	mul	r19, r0, r18	C a0 b0
	umulh	r0, r0, r18
	mul	r20, r1, r18	C a1 b0
	umulh	r1, r1, r18
	mul	r21, r2, r18	C a2 b0
	umulh	r2, r2, r18
	mul	r22, r6, r18	C a3 b0
	umulh	r6, r6, r18
	mul	r23, r7, r18	C a4 b0
	umulh	r7, r7, r18
	mul	r24, r8, r18	C a5 b0
	umulh	r8, r8, r18

	adds	r20, r20, r0
	adcs	r21, r21, r1
	adcs	r22, r22, r2
	ldp	r0, r1, [ap,#6*8]
	adcs	r23, r23, r6
	adcs	r24, r24, r7
	ldp	r2, r6, [ap,#8*8]

	mul	r9, r0, r18	C a6 b0
	umulh	r0, r0, r18
	mul	r4, r1, r18	C a7 b0
	umulh	r1, r1, r18
	mul	r5, r2, r18	C a8 b0
	umulh	r2, r2, r18
	mul	r3, r6, r18	C a9 b0
	umulh	r6, r6, r18

	adcs	r9, r9, r8
	adcs	r4, r4, r0
	ldp	r8, r0, [ap,#10*8]
	adcs	r5, r5, r1
	adcs	r3, r3, r2

	mul	r1, r8, r18	C a10 b0
	umulh	r8, r8, r18
	mul	r2, r0, r18	C a11 b0
	umulh	r0, r0, r18

	sub	n, n, #1
	adcs	r1, r1, r6
	adcs	r2, r2, r8
	cinc	r0, r0, cs
	C 19, 20, 21, 22, 23, 24, 9, 4, 5, 3, 1, 2, 0

	ldp	r6, r7, [sp, #8*8]
	ldp	r8, r18, [sp, #10*8]

	C 19, (10, 11, 12, 13, 14, 15, 16, 17, 9, 4, 5,  3)
	C  <- (10, 11, 12, 13, 14, 15, 16, 17, 6, 7, 8, 18)
	C   + (19, 20, 21, 22, 23, 24,  9,  4, 5, 3, 1,  2, 0)
	adds	r19, r10, r19
	adcs	r10, r11, r20
	adcs	r11, r12, r21
	adcs	r12, r13, r22
	adcs	r13, r14, r23
	adcs	r14, r15, r24
	adcs	r15, r16, r9
	adcs	r16, r17, r4
	adcs	r17, r6, r5
	adcs	r9, r7, r3
	adcs	r4, r8, r1
	adcs	r5, r18, r2
	cinc	r3, r0, cs

	str	r19, [rp], #1*8

	cbnz	n, L(a12)

	ldp	r20, r21, [sp, #2*8]
	ldp	r22, r23, [sp, #4*8]
	ldr	r24, [sp, #6*8]
	ldp	r18, r19, [sp], #12*8

L(f12):	stp	r10, r11, [rp]
	stp	r12, r13, [rp,#2*8]
	stp	r14, r15, [rp,#4*8]
	stp	r16, r17, [rp,#6*8]
	stp	r9, r4, [rp,#8*8]
	stp	r5, r3, [rp,#10*8]
	ldp	r16, r17, [sp, #2*8]
	ldp	r14, r15, [sp], #4*8
	mov	res, r3
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_13n)
	ldr	r8, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]
	ldp	r6, r7, [ap,#6*8]

	stp	r14, r15, [sp, #-4*8]!
	stp	r16, r17, [sp, #2*8]
L(m13):	mul	r9, r0, r8	C a0 b0
	umulh	r0, r0, r8
	mul	r10, r1, r8	C a1 b0
	umulh	r1, r1, r8
	mul	r11, r2, r8	C a2 b0
	umulh	r2, r2, r8
	mul	r12, r3, r8	C a3 b0
	umulh	r3, r3, r8
	mul	r13, r4, r8	C a4 b0
	umulh	r4, r4, r8
	mul	r14, r5, r8	C a5 b0
	umulh	r5, r5, r8
	mul	r15, r6, r8	C a6 b0
	umulh	r6, r6, r8
	mul	r16, r7, r8	C a7 b0
	umulh	r7, r7, r8

	str	r9, [rp], #1*8
	adds	r10, r10, r0
	adcs	r11, r11, r1
	ldp	r0, r1, [ap,#8*8]
	adcs	r12, r12, r2
	adcs	r13, r13, r3
	adcs	r14, r14, r4
	adcs	r15, r15, r5
	adcs	r16, r16, r6
	ldp	r2, r3, [ap,#10*8]
	ldr	r4, [ap,#12*8]

	mul	r17, r0, r8	C a8 b0
	umulh	r0, r0, r8
	mul	r9, r1, r8	C a9 b0
	umulh	r1, r1, r8
	mul	r5, r2, r8	C a10 b0
	umulh	r2, r2, r8
	mul	r6, r3, r8	C a11 b0
	umulh	r3, r3, r8

	sub	n, n, #1
	adcs	r17, r17, r7
	adcs	r9, r9, r0

	mul	r7, r4, r8	C a12 b0
	umulh	r4, r4, r8

	adcs	r5, r5, r1
	adcs	r6, r6, r2
	adcs	r7, r7, r3
	cinc	r4, r4, cs
	C 10, 11, 12, 13, 14, 15, 16, 17, 9, 5, 6, 7, 4
	C 0, 1, 2, 3, 8, 18, 19, ...

	cbz	n, L(f13)

	stp	r18, r19, [sp, #-14*8]!
	stp	r20, r21, [sp, #2*8]
	stp	r22, r23, [sp, #4*8]
	str	r24, [sp, #6*8]
L(a13):	ldp	r0, r1, [ap]
	ldr	r8, [bp], #1*8
	ldp	r2, r3, [ap,#2*8]
	ldp	r18, r19, [ap,#4*8]
	ldp	r20, r21, [ap,#6*8]

	stp	r17, r9, [sp, #8*8]
	stp	r5, r6, [sp, #10*8]
	stp	r7, r4, [sp, #12*8]

	mul	r22, r0, r8	C a0 b0
	umulh	r0, r0, r8
	mul	r23, r1, r8	C a1 b0
	umulh	r1, r1, r8
	mul	r24, r2, r8	C a2 b0
	umulh	r2, r2, r8
	mul	r17, r3, r8	C a3 b0
	umulh	r3, r3, r8
	mul	r9, r18, r8	C a4 b0
	umulh	r18, r18, r8
	mul	r5, r19, r8	C a5 b0
	umulh	r19, r19, r8
	mul	r6, r20, r8	C a6 b0
	umulh	r20, r20, r8
	mul	r7, r21, r8	C a7 b0
	umulh	r21, r21, r8

	adds	r23, r23, r0
	adcs	r24, r24, r1
	adcs	r17, r17, r2
	ldp	r0, r1, [ap,#8*8]
	adcs	r9, r9, r3
	adcs	r5, r5, r18
	ldp	r2, r3, [ap,#10*8]
	ldr	r18, [ap,#12*8]
	adcs	r6, r6, r19
	adcs	r7, r7, r20

	mul	r4, r0, r8	C a8 b0
	umulh	r0, r0, r8
	mul	r19, r1, r8	C a9 b0
	umulh	r1, r1, r8
	mul	r20, r2, r8	C a10 b0
	umulh	r2, r2, r8

	sub	n, n, #1
	adcs	r4, r4, r21
	adcs	r19, r19, r0

	mul	r21, r3, r8	C a11 b0
	umulh	r3, r3, r8
	mul	r0, r18, r8	C a12 b0
	umulh	r18, r18, r8

	adcs	r20, r20, r1
	adcs	r21, r21, r2
	adcs	r0, r0, r3
	cinc	r18, r18, cs
	C 22, 23, 24, 17, 9, 5, 6, 7, 4, 19, 20, 21, 0, 18

	ldp	r1, r2, [sp, #8*8]
	ldp	r3, r8, [sp, #10*8]

	C 22, (10, 11, 12, 13, 14, 15, 16, 17, 9,  5,  6,  7, 4)
	C  <- (10, 11, 12, 13, 14, 15, 16,  1, 2,  3,  8,  X, X)
	C   + (22, 23, 24, 17,  9,  5,  6,  7, 4, 19, 20, 21, 0, 18)
	adds	r22, r10, r22
	adcs	r10, r11, r23
	adcs	r11, r12, r24
	ldp	r23, r24, [sp, #12*8]
	adcs	r12, r13, r17
	adcs	r13, r14, r9
	adcs	r14, r15, r5
	adcs	r15, r16, r6
	adcs	r16, r1, r7
	adcs	r17, r2, r4
	adcs	r9, r3, r19
	adcs	r5, r8, r20
	adcs	r6, r23, r21
	adcs	r7, r24, r0
	cinc	r4, r18, cs

	str	r22, [rp], #1*8

	cbnz	n, L(a13)

	ldp	r20, r21, [sp, #2*8]
	ldp	r22, r23, [sp, #4*8]
	ldr	r24, [sp, #6*8]
	ldp	r18, r19, [sp], #14*8

L(f13):	stp	r10, r11, [rp]
	stp	r12, r13, [rp,#2*8]
	stp	r14, r15, [rp,#4*8]
	stp	r16, r17, [rp,#6*8]
	stp	r9, r5, [rp,#8*8]
	stp	r6, r7, [rp,#10*8]
	str	r4, [rp,#12*8]
	ldp	r16, r17, [sp, #2*8]
	ldp	r14, r15, [sp], #4*8
	mov	res, r4
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_14n)
	ldr	r8, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]
	ldp	r6, r7, [ap,#6*8]

	stp	r14, r15, [sp, #-6*8]!
	stp	r16, r17, [sp, #2*8]
	stp	r18, r19, [sp, #4*8]
L(m14):	mul	r9, r0, r8	C a0 b0
	umulh	r0, r0, r8
	mul	r10, r1, r8	C a1 b0
	umulh	r1, r1, r8
	mul	r11, r2, r8	C a2 b0
	umulh	r2, r2, r8
	mul	r12, r3, r8	C a3 b0
	umulh	r3, r3, r8
	mul	r13, r4, r8	C a4 b0
	umulh	r4, r4, r8
	mul	r14, r5, r8	C a5 b0
	umulh	r5, r5, r8
	mul	r15, r6, r8	C a6 b0
	umulh	r6, r6, r8
	mul	r16, r7, r8	C a7 b0
	umulh	r7, r7, r8

	str	r9, [rp], #1*8
	adds	r10, r10, r0
	adcs	r11, r11, r1
	ldp	r0, r1, [ap,#8*8]
	adcs	r12, r12, r2
	adcs	r13, r13, r3
	adcs	r14, r14, r4
	adcs	r15, r15, r5
	adcs	r16, r16, r6
	ldp	r2, r3, [ap,#10*8]
	ldp	r4, r5, [ap,#12*8]

	mul	r17, r0, r8	C a8 b0
	umulh	r0, r0, r8
	mul	r18, r1, r8	C a9 b0
	umulh	r1, r1, r8
	mul	r19, r2, r8	C a10 b0
	umulh	r2, r2, r8
	mul	r9, r3, r8	C a11 b0
	umulh	r3, r3, r8
	mul	r6, r4, r8	C a12 b0
	umulh	r4, r4, r8

	sub	n, n, #1
	adcs	r17, r17, r7
	adcs	r18, r18, r0

	mul	r7, r5, r8	C a13 b0
	umulh	r5, r5, r8

	adcs	r19, r19, r1
	adcs	r9, r9, r2
	adcs	r6, r6, r3
	adcs	r7, r7, r4
	cinc	r5, r5, cs
	C 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 9, 6, 7, 5
	C 0, 1, 2, 3, 4, 8, 20, 21, ...

	cbz	n, L(f14)

	stp	r20, r21, [sp, #-14*8]!
	stp	r22, r23, [sp, #2*8]
	str	r24, [sp, #4*8]
L(a14):	ldp	r0, r1, [ap]
	ldr	r8, [bp], #1*8
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r20, [ap,#4*8]
	ldp	r21, r22, [ap,#6*8]

	stp	r16, r17, [sp, #6*8]
	stp	r18, r19, [sp, #8*8]
	stp	r9, r6, [sp, #10*8]
	stp	r7, r5, [sp, #12*8]

	mul	r23, r0, r8	C a0 b0
	umulh	r0, r0, r8
	mul	r24, r1, r8	C a1 b0
	umulh	r1, r1, r8
	mul	r16, r2, r8	C a2 b0
	umulh	r2, r2, r8
	mul	r17, r3, r8	C a3 b0
	umulh	r3, r3, r8
	mul	r18, r4, r8	C a4 b0
	umulh	r4, r4, r8
	mul	r19, r20, r8	C a5 b0
	umulh	r20, r20, r8
	mul	r9, r21, r8	C a6 b0
	umulh	r21, r21, r8
	mul	r6, r22, r8	C a7 b0
	umulh	r22, r22, r8

	adds	r24, r24, r0
	adcs	r16, r16, r1
	adcs	r17, r17, r2
	ldp	r0, r1, [ap,#8*8]
	adcs	r18, r18, r3
	adcs	r19, r19, r4
	ldp	r2, r3, [ap,#10*8]
	adcs	r9, r9, r20
	adcs	r6, r6, r21

	mul	r7, r0, r8	C a8 b0
	umulh	r0, r0, r8
	mul	r5, r1, r8	C a9 b0
	umulh	r1, r1, r8
	mul	r4, r2, r8	C a10 b0
	umulh	r2, r2, r8
	mul	r20, r3, r8	C a11 b0
	umulh	r3, r3, r8

	adcs	r7, r7, r22
	adcs	r5, r5, r0
	ldp	r22, r0, [ap,#12*8]
	adcs	r4, r4, r1
	adcs	r20, r20, r2

	mul	r21, r22, r8	C a12 b0
	umulh	r22, r22, r8
	mul	r1, r0, r8	C a13 b0
	umulh	r0, r0, r8

	sub	n, n, #1
	adcs	r21, r21, r3
	adcs	r1, r1, r22
	cinc	r0, r0, cs
	C 23, 24, 16, 17, 18, 19, 9, 6, 7, 5, 4, 20, 21, 1, 0

	C 2, 3, 8, 22
	ldp	r2, r3, [sp, #6*8]
	ldp	r8, r22, [sp, #8*8]

	C 23, (10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 9,  6,  7, 5)
	C  <- (10, 11, 12, 13, 14, 15,  2,  3,  8, 22, X,  X,  X, X)
	C   + (23, 24, 16, 17, 18, 19,  9,  6,  7,  5, 4, 20, 21, 1, 0)
	adds	r23, r10, r23
	adcs	r10, r11, r24
	adcs	r11, r12, r16
	str	r23, [rp], #1*8
	adcs	r12, r13, r17
	adcs	r13, r14, r18
	adcs	r14, r15, r19
	ldp	r23, r24, [sp,#10*8]
	adcs	r15, r2, r9
	adcs	r16, r3, r6
	adcs	r17, r8, r7
	ldp	r6, r7, [sp,#12*8]
	adcs	r18, r22, r5
	adcs	r19, r23, r4
	adcs	r9, r24, r20
	adcs	r6, r6, r21
	adcs	r7, r7, r1
	cinc	r5, r0, cs

	cbnz	n, L(a14)

	ldp	r22, r23, [sp, #2*8]
	ldr	r24, [sp, #4*8]
	ldp	r20, r21, [sp], #14*8

L(f14):	stp	r10, r11, [rp]
	stp	r12, r13, [rp,#2*8]
	stp	r14, r15, [rp,#4*8]
	stp	r16, r17, [rp,#6*8]
	stp	r18, r19, [rp,#8*8]
	stp	r9, r6, [rp,#10*8]
	stp	r7, r5, [rp,#12*8]
	ldp	r16, r17, [sp, #2*8]
	ldp	r18, r19, [sp, #4*8]
	ldp	r14, r15, [sp], #6*8
	mov	res, r5
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mul_15n)
	ldr	r8, [bp], #1*8
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap,#2*8]
	ldp	r4, r5, [ap,#4*8]
	ldp	r6, r7, [ap,#6*8]

	stp	r14, r15, [sp, #-6*8]!
	stp	r16, r17, [sp, #2*8]
	stp	r18, r19, [sp, #4*8]
L(m15):	mul	r9, r0, r8	C a0 b0
	umulh	r0, r0, r8
	mul	r10, r1, r8	C a1 b0
	umulh	r1, r1, r8
	mul	r11, r2, r8	C a2 b0
	umulh	r2, r2, r8
	mul	r12, r3, r8	C a3 b0
	umulh	r3, r3, r8
	mul	r13, r4, r8	C a4 b0
	umulh	r4, r4, r8
	mul	r14, r5, r8	C a5 b0
	umulh	r5, r5, r8
	mul	r15, r6, r8	C a6 b0
	umulh	r6, r6, r8
	mul	r16, r7, r8	C a7 b0
	umulh	r7, r7, r8

	str	r9, [rp], #1*8
	adds	r10, r10, r0
	adcs	r11, r11, r1
	ldp	r0, r1, [ap,#8*8]
	adcs	r12, r12, r2
	adcs	r13, r13, r3
	adcs	r14, r14, r4
	adcs	r15, r15, r5
	adcs	r16, r16, r6
	ldp	r2, r3, [ap,#10*8]
	ldp	r4, r5, [ap,#12*8]
	ldr	r6, [ap,#14*8]

	mul	r17, r0, r8	C a8 b0
	umulh	r0, r0, r8
	mul	r18, r1, r8	C a9 b0
	umulh	r1, r1, r8
	mul	r19, r2, r8	C a10 b0
	umulh	r2, r2, r8
	mul	r9, r3, r8	C a11 b0
	umulh	r3, r3, r8

	sub	n, n, #1
	adcs	r17, r17, r7
	adcs	r18, r18, r0
	adcs	r19, r19, r1

	mul	r7, r4, r8	C a12 b0
	umulh	r4, r4, r8
	mul	r0, r5, r8	C a13 b0
	umulh	r5, r5, r8
	mul	r1, r6, r8	C a14 b0
	umulh	r6, r6, r8

	adcs	r9, r9, r2
	adcs	r7, r7, r3
	adcs	r0, r0, r4
	adcs	r1, r1, r5
	cinc	r6, r6, cs
	C 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 9, 7, 0, 1, 6
	C 2, 3, 4, 5, 8

	cbz	n, L(f15)

	stp	r20, r21, [sp, #-14*8]!
	stp	r22, r23, [sp, #2*8]
	str	r24, [sp, #4*8]
L(a15):	ldp	r2, r3, [ap]
	ldr	r8, [bp], #1*8
	ldp	r4, r5, [ap,#2*8]
	ldp	r20, r21, [ap,#4*8]
	ldp	r22, r23, [ap,#6*8]

	stp	r17, r18, [sp, #6*8]
	stp	r19, r9, [sp, #8*8]
	stp	r7, r0, [sp, #10*8]
	stp	r1, r6, [sp, #12*8]

	mul	r17, r2, r8	C a0 b0
	umulh	r2, r2, r8
	mul	r18, r3, r8	C a1 b0
	umulh	r3, r3, r8
	mul	r19, r4, r8	C a2 b0
	umulh	r4, r4, r8
	mul	r9, r5, r8	C a3 b0
	umulh	r5, r5, r8
	mul	r7, r20, r8	C a4 b0
	umulh	r20, r20, r8
	mul	r0, r21, r8	C a5 b0
	umulh	r21, r21, r8
	mul	r1, r22, r8	C a6 b0
	umulh	r22, r22, r8
	mul	r6, r23, r8	C a7 b0
	umulh	r23, r23, r8

	adds	r18, r18, r2
	adcs	r19, r19, r3
	adcs	r9, r9, r4
	ldp	r2, r3, [ap,#8*8]
	adcs	r7, r7, r5
	ldp	r4, r5, [ap,#10*8]
	adcs	r0, r0, r20
	adcs	r1, r1, r21
	adcs	r6, r6, r22

	mul	r24, r2, r8	C a8 b0
	umulh	r2, r2, r8
	mul	r20, r3, r8	C a9 b0
	umulh	r3, r3, r8
	mul	r21, r4, r8	C a10 b0
	umulh	r4, r4, r8
	mul	r22, r5, r8	C a11 b0
	umulh	r5, r5, r8

	adcs	r24, r24, r23
	adcs	r20, r20, r2
	ldp	r23, r2, [ap,#12*8]
	adcs	r21, r21, r3
	adcs	r22, r22, r4

	mul	r3, r23, r8	C a12 b0
	umulh	r23, r23, r8
	mul	r4, r2, r8	C a13 b0
	umulh	r2, r2, r8

	adcs	r3, r3, r5
	ldr	r5, [ap,#14*8]
	adcs	r4, r4, r23
	sub	n, n, #1

	mul	r23, r5, r8	C a14 b0
	umulh	r5, r5, r8

	adcs	r23, r23, r2
	cinc	r5, r5, cs
	C 17, 18, 19, 9, 7, 0, 1, 6, 24, 20, 21, 22, 3, 4, 23, 5

	C 2, 8
	ldp	r2, r8, [sp, #6*8]

	C 17, (10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  9,  7, 0, 1,  6)
	C  <- (10, 11, 12, 13, 14, 15, 16,  X,  X,  X,  X,  X, X, X,  X)
	C   + (17, 18, 19,  9,  7,  0,  1,  6, 24, 20, 21, 22, 3, 4, 23, 5)
	adds	r17, r10, r17
	adcs	r10, r11, r18
	adcs	r11, r12, r19
	str	r17, [rp], #1*8
	ldp	r18, r19, [sp, #8*8]
	adcs	r12, r13, r9
	adcs	r13, r14, r7
	ldp	r9, r7, [sp, #10*8]
	adcs	r14, r15, r0
	adcs	r15, r16, r1
	ldp	r0, r1, [sp, #12*8]
	adcs	r16, r2, r6	C X
	adcs	r17, r8, r24	C X
	adcs	r18, r18, r20	C X
	adcs	r19, r19, r21	C X
	adcs	r9, r9, r22	C X
	adcs	r7, r7, r3	C X
	adcs	r0, r0, r4	C X
	adcs	r1, r1, r23	C X
	cinc	r6, r5, cs

	cbnz	n, L(a15)

	ldp	r22, r23, [sp, #2*8]
	ldr	r24, [sp, #4*8]
	ldp	r20, r21, [sp], #14*8

L(f15):	stp	r10, r11, [rp]
	stp	r12, r13, [rp,#2*8]
	stp	r14, r15, [rp,#4*8]
	stp	r16, r17, [rp,#6*8]
	stp	r18, r19, [rp,#8*8]
	stp	r9, r7, [rp,#10*8]
	stp	r0, r1, [rp,#12*8]
	str	r6, [rp,#14*8]
	ldp	r16, r17, [sp, #2*8]
	ldp	r18, r19, [sp, #4*8]
	ldp	r14, r15, [sp], #6*8
	mov	res, r6
	ret
EPILOGUE()
