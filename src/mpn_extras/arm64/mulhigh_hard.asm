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

define(`r0',  `x3')
define(`r1',  `x4')
define(`r2',  `x5')
define(`r3',  `x6')
define(`r4',  `x7')
define(`r5',  `x8')
define(`r6',  `x9')
define(`r7',  `x10')
define(`r8',  `x11')
define(`r9',  `x12')
define(`r10', `x13')
define(`r11', `x14')
define(`r12', `x15')
define(`r13', `x16')
define(`r14', `x17')

define(`r15', `ap')	dnl Beware!
define(`r16', `bp')

dnl NOTE: The following has to be pushed and popped
dnl NOTE: At least on Apple, we need to preserve both x18 and x29 at all times!
define(`r17', `x19')
define(`r18', `x20')
define(`r19', `x21')
define(`r20', `x22')
define(`r21', `x23')
define(`r22', `x24')
define(`r23', `x25')
define(`r24', `x26')
define(`r25', `x27')
define(`r26', `x28')
define(`r27', `x30')

define(`res', `rp') dnl NOTE: Synonymous with rp!

PROLOGUE(flint_mpn_mulhigh_1)
	ldr	r0, [ap]
	ldr	r1, [bp]

	umulh	r3, r0, r1
	mul	r2, r0, r1

	str	r3, [rp]
	mov	res, r2

	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_2)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [bp]

	umulh	r4, r0, r2

	mul	r5, r0, r3
	umulh	r6, r0, r3

	mul	r7, r1, r2
	umulh	r8, r1, r2

	mul	r9, r1, r3
	umulh	r10, r1, r3
	C (4, 5, 7), (6, 8, 9), 10

	adds	r4, r4, r5
	adcs	r6, r6, r8
	cinc	r10, r10, cs
	adds	r4, r4, r7
	adcs	r6, r6, r9
	cinc	r10, r10, cs

	stp	r6, r10, [rp]
	mov	res, r4

	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_3)
	ldp	r0, r1, [ap]
	ldp	r3, r4, [bp]
	ldr	r2, [ap, #2*8]
	ldr	r5, [bp, #2*8]

	umulh	r6, r1, r3
	mul	r7, r2, r3
	umulh	r3, r2, r3
	umulh	r8, r1, r4
	umulh	r9, r2, r4
	umulh	r10, r0, r4
	mul	r11, r2, r4
	mul	r12, r0, r5

	adds	r6, r6, r7

	umulh	r0, r0, r5
	umulh	r13, r1, r5

	adcs	r3, r3, r8

	umulh	r14, r2, r5
	mul	r4, r1, r4

	cinc	r9, r9, cs
	adds	r6, r6, r10

	mul	r1, r1, r5
	mul	r2, r2, r5

	adcs	r3, r3, r11
	cinc	r9, r9, cs

	adds	r6, r6, r12
	adcs	r3, r3, r0
	adcs	r9, r9, r13
	cinc	r14, r14, cs

	adds	r6, r6, r4
	adcs	r3, r3, r1
	adcs	r9, r9, r2
	cinc	r14, r14, cs

	stp	r3, r9, [rp]
	str	r14, [rp, #2*8]
	mov	res, r6

	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_4)
	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [bp]
	ldp	r6, r7, [bp, #2*8]

	umulh	r8, r2, r4
	mul	r9, r3, r4
	umulh	r4, r3, r4
	umulh	r12, r2, r5
	umulh	r14, r3, r5
	umulh	r10, r1, r5
	mul	r11, r2, r5
	mul	r13, r3, r5
	C (8, 9, 10, 11), (4, 12, 13), 14
	C Free: 5, 15, 16

	adds	r8, r8, r9

	umulh	r9, r0, r6
	umulh	r5, r1, r6

	adcs	r4, r4, r12

	umulh	r12, r2, r6
	umulh	r15, r3, r6

	cinc	r14, r14, cs
	adds	r8, r8, r10

	mul	r16, r1, r6
	mul	r10, r2, r6

	adcs	r4, r4, r13
	cinc	r14, r14, cs
	adds	r8, r8, r11

	mul	r13, r3, r6
	umulh	r11, r0, r7

	adcs	r4, r4, r5
	adcs	r14, r14, r12

	umulh	r5, r1, r7
	umulh	r12, r2, r7

	cinc	r15, r15, cs
	adds	r8, r8, r9
	adcs	r4, r4, r10

	umulh	r9, r3, r7
	mul	r10, r0, r7

	adcs	r14, r14, r13
	cinc	r15, r15, cs
	adds	r8, r8, r16

	mul	r13, r1, r7
	mul	r16, r2, r7

	adcs	r4, r4, r11
	adcs	r14, r14, r5

	mul	r11, r3, r7

	adcs	r15, r15, r12
	cinc	r9, r9, cs
	C (8, 10), (4, 13), (14, 16), (15, 11), 9

	adds	r8, r8, r10
	adcs	r4, r4, r13
	adcs	r14, r14, r16
	adcs	r15, r15, r11
	cinc	r9, r9, cs

	stp	r4, r14, [rp]
	stp	r15, r9, [rp, #2*8]
	mov	res, r8

	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_5)
	ldp	r2, r3, [ap, #2*8]
	ldr	r4, [ap, #4*8]
	ldp	r5, r6, [bp]
	ldp	r0, r1, [ap]
	ldp	r7, r8, [bp, #2*8]
	ldr	r9, [bp, #4*8]

	umulh	r10, r3, r5
	mul	r11, r4, r5
	umulh	r5, r4, r5
	umulh	r12, r3, r6
	umulh	r13, r4, r6
	umulh	r14, r2, r6
	mul	r15, r3, r6
	mul	r6, r4, r6
	C (10, 11, 14, 15), (5, 12, 6), 13
	C Free: 16

	adds	r10, r10, r11
	umulh	r16, r1, r7
	umulh	r11, r2, r7
	adcs	r5, r5, r12
	umulh	r12, r3, r7
	cinc	r13, r13, cs
	adds	r10, r10, r14
	umulh	r14, r4, r7
	adcs	r5, r5, r6
	mul	r6, r2, r7
	cinc	r13, r13, cs
	adds	r10, r10, r15
	mul	r15, r3, r7
	adcs	r5, r5, r11
	mul	r7, r4, r7
	umulh	r11, r0, r8
	adcs	r13, r13, r12
	umulh	r12, r1, r8
	cinc	r14, r14, cs
	adds	r10, r10, r16
	umulh	r16, r2, r8
	adcs	r5, r5, r15
	umulh	r15, r3, r8
	adcs	r13, r13, r7
	umulh	r7, r4, r8
	cinc	r14, r14, cs
	adds	r10, r10, r6
	mul	r6, r1, r8
	adcs	r5, r5, r12
	mul	r12, r2, r8
	adcs	r13, r13, r16
	mul	r16, r3, r8
	adcs	r14, r14, r15
	mul	r8, r4, r8
	umulh	r15, r0, r9
	cinc	r7, r7, cs
	adds	r10, r10, r11
	umulh	r11, r1, r9
	adcs	r5, r5, r12
	umulh	r12, r2, r9
	adcs	r13, r13, r16
	umulh	r16, r3, r9
	adcs	r14, r14, r8
	umulh	r8, r4, r9
	mul	r0, r0, r9
	cinc	r7, r7, cs
	mul	r1, r1, r9
	mul	r2, r2, r9
	mul	r3, r3, r9
	mul	r4, r4, r9

	adds	r10, r10, r6
	adcs	r5, r5, r15
	adcs	r13, r13, r11
	adcs	r14, r14, r12
	adcs	r7, r7, r16
	cinc	r8, r8, cs
	adds	r10, r10, r0
	adcs	r5, r5, r1
	adcs	r13, r13, r2
	adcs	r14, r14, r3
	stp	r5, r13, [rp]
	adcs	r7, r7, r4
	cinc	r8, r8, cs
	stp	r14, r7, [rp, #2*8]
	str	r8, [rp, #4*8]
	mov	res, r10
	
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_6)
	ldp	r2, r3, [ap, #2*8]
	ldp	r4, r5, [ap, #4*8]
	ldp	r6, r7, [bp]

	ldp	r0, r1, [ap]
	ldp	r8, r9, [bp, #2*8]
	ldp	r10, r11, [bp, #4*8]

	stp	r17, r18, [sp, #-2*8]!

	umulh	r12, r4, r6
	mul	r13, r5, r6
	umulh	r6, r5, r6
	umulh	r14, r4, r7
	umulh	r15, r5, r7
	umulh	r16, r3, r7
	mul	r17, r4, r7
	mul	r7, r5, r7
	C (12, 13, 16, 17), (6, 14, 7), 15
	C Free: 18

	umulh	r18, r2, r8
	adds	r12, r12, r13
	umulh	r13, r3, r8
	adcs	r6, r6, r14
	umulh	r14, r4, r8
	cinc	r15, r15, cs
	adds	r12, r12, r16
	umulh	r16, r5, r8
	adcs	r6, r6, r7
	mul	r7, r3, r8
	cinc	r15, r15, cs
	adds	r12, r12, r17
	mul	r17, r4, r8
	mul	r8, r5, r8
	adcs	r6, r6, r13
	umulh	r13, r1, r9
	adcs	r15, r15, r14
	umulh	r14, r2, r9
	cinc	r16, r16, cs
	adds	r12, r12, r18
	umulh	r18, r3, r9
	adcs	r6, r6, r17
	umulh	r17, r4, r9
	adcs	r15, r15, r8
	umulh	r8, r5, r9
	cinc	r16, r16, cs
	adds	r12, r12, r7
	mul	r7, r2, r9
	adcs	r6, r6, r14
	mul	r14, r3, r9
	adcs	r15, r15, r18
	mul	r18, r4, r9
	mul	r9, r5, r9
	adcs	r16, r16, r17
	umulh	r17, r0, r10
	cinc	r8, r8, cs
	adds	r12, r12, r13
	umulh	r13, r1, r10
	adcs	r6, r6, r14
	umulh	r14, r2, r10
	adcs	r15, r15, r18
	umulh	r18, r3, r10
	adcs	r16, r16, r9
	umulh	r9, r4, r10
	cinc	r8, r8, cs
	adds	r12, r12, r7
	umulh	r7, r5, r10
	adcs	r6, r6, r13
	mul	r13, r1, r10
	adcs	r15, r15, r14
	mul	r14, r2, r10
	adcs	r16, r16, r18
	mul	r18, r3, r10
	adcs	r8, r8, r9
	mul	r9, r4, r10
	mul	r10, r5, r10
	cinc	r7, r7, cs
	adds	r12, r12, r17
	umulh	r17, r0, r11
	adcs	r6, r6, r14
	umulh	r14, r1, r11
	adcs	r15, r15, r18
	umulh	r18, r2, r11
	adcs	r16, r16, r9
	umulh	r9, r3, r11
	adcs	r8, r8, r10
	umulh	r10, r4, r11
	cinc	r7, r7, cs
	adds	r12, r12, r13
	umulh	r13, r5, r11
	mul	r0, r0, r11
	adcs	r6, r6, r17
	mul	r1, r1, r11
	mul	r2, r2, r11
	adcs	r15, r15, r14
	mul	r3, r3, r11
	mul	r4, r4, r11
	adcs	r16, r16, r18
	mul	r5, r5, r11
	adcs	r8, r8, r9
	adcs	r7, r7, r10
	cinc	r13, r13, cs
	ldp	r17, r18, [sp], #2*8
	adds	r12, r12, r0
	adcs	r6, r6, r1
	adcs	r15, r15, r2
	adcs	r16, r16, r3
	stp	r6, r15, [rp]
	adcs	r8, r8, r4
	adcs	r7, r7, r5
	stp	r16, r8, [rp, #2*8]
	cinc	r13, r13, cs
	stp	r7, r13, [rp, #4*8]

	mov	res, r12

	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_7)
	ldp	r4, r5, [ap, #4*8]
	ldr	r6, [ap, #6*8]
	ldp	r7, r8, [bp]
	stp	r17, r18, [sp, #-4*8]!

	ldp	r9, r10, [bp, #2*8]
	ldp	r11, r12, [bp, #4*8]
	ldr	r13, [bp, #6*8]
	stp	r19, r20, [sp, #2*8]

	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]

	umulh	r14, r5, r7
	mul	r15, r6, r7
	umulh	r7, r6, r7
	umulh	r16, r5, r8
	umulh	r17, r6, r8
	umulh	r18, r4, r8
	mul	r19, r5, r8
	mul	r8, r6, r8
	C (14, 15, 18, 19), (7, 16, 8), 17
	C Free: 20

	umulh	r20, r3, r9
	adds	r14, r14, r15
	umulh	r15, r4, r9
	adcs	r7, r7, r16
	umulh	r16, r5, r9
	cinc	r17, r17, cs
	adds	r14, r14, r18
	umulh	r18, r6, r9
	adcs	r7, r7, r8
	mul	r8, r4, r9
	cinc	r17, r17, cs
	adds	r14, r14, r19
	mul	r19, r5, r9
	adcs	r7, r7, r15
	mul	r9, r6, r9

	umulh	r15, r2, r10
	adcs	r17, r17, r16
	umulh	r16, r3, r10
	cinc	r18, r18, cs
	adds	r14, r14, r20
	umulh	r20, r4, r10
	adcs	r7, r7, r19
	umulh	r19, r5, r10
	adcs	r17, r17, r9
	umulh	r9, r6, r10
	cinc	r18, r18, cs
	adds	r14, r14, r8
	mul	r8, r3, r10
	adcs	r7, r7, r16
	mul	r16, r4, r10
	adcs	r17, r17, r20
	mul	r20, r5, r10
	adcs	r18, r18, r19
	mul	r10, r6, r10

	umulh	r19, r1, r11
	cinc	r9, r9, cs
	adds	r14, r14, r15
	umulh	r15, r2, r11
	adcs	r7, r7, r16
	umulh	r16, r3, r11
	adcs	r17, r17, r20
	umulh	r20, r4, r11
	adcs	r18, r18, r10
	umulh	r10, r5, r11
	cinc	r9, r9, cs
	adds	r14, r14, r8
	umulh	r8, r6, r11
	adcs	r7, r7, r15
	mul	r15, r2, r11
	adcs	r17, r17, r16
	mul	r16, r3, r11
	adcs	r18, r18, r20
	mul	r20, r4, r11
	adcs	r9, r9, r10
	mul	r10, r5, r11
	mul	r11, r6, r11

	cinc	r8, r8, cs
	adds	r14, r14, r19
	umulh	r19, r0, r12
	adcs	r7, r7, r16
	umulh	r16, r1, r12
	adcs	r17, r17, r20
	umulh	r20, r2, r12
	adcs	r18, r18, r10
	umulh	r10, r3, r12
	adcs	r9, r9, r11
	umulh	r11, r4, r12
	cinc	r8, r8, cs
	adds	r14, r14, r15
	umulh	r15, r5, r12
	adcs	r7, r7, r16
	umulh	r16, r6, r12
	adcs	r17, r17, r20
	mul	r20, r1, r12
	adcs	r18, r18, r10
	mul	r10, r2, r12
	adcs	r9, r9, r11
	mul	r11, r3, r12
	adcs	r8, r8, r15
	mul	r15, r4, r12
	cinc	r16, r16, cs
	adds	r14, r14, r19
	mul	r19, r5, r12
	mul	r12, r6, r12

	adcs	r7, r7, r10
	umulh	r10, r0, r13
	adcs	r17, r17, r11
	umulh	r11, r1, r13
	adcs	r18, r18, r15
	umulh	r15, r2, r13
	adcs	r9, r9, r19
	umulh	r19, r3, r13
	adcs	r8, r8, r12
	umulh	r12, r4, r13
	cinc	r16, r16, cs
	adds	r14, r14, r20
	umulh	r20, r5, r13
	adcs	r7, r7, r10
	umulh	r10, r6, r13
	mul	r0, r0, r13
	adcs	r17, r17, r11
	mul	r1, r1, r13
	mul	r2, r2, r13
	adcs	r18, r18, r15
	mul	r3, r3, r13
	mul	r4, r4, r13
	adcs	r9, r9, r19
	mul	r5, r5, r13
	mul	r6, r6, r13
	adcs	r8, r8, r12
	adcs	r16, r16, r20
	ldp	r19, r20, [sp, #2*8]
	cinc	r10, r10, cs
	adds	r14, r14, r0
	adcs	r7, r7, r1
	adcs	r11, r17, r2
	adcs	r15, r18, r3
	adcs	r9, r9, r4
	stp	r7, r11, [rp]
	ldp	r17, r18, [sp], #4*8
	adcs	r8, r8, r5
	stp	r15, r9, [rp, #2*8]
	adcs	r16, r16, r6
	cinc	r10, r10, cs
	stp	r8, r16, [rp, #4*8]
	str	r10, [rp, #6*8]
	mov	res, r14
	ret
EPILOGUE()

PROLOGUE(flint_mpn_mulhigh_8)
	ldp	r4, r5, [ap, #4*8]
	ldp	r6, r7, [ap, #6*8]
	ldp	r8, r9, [bp]
	stp	r17, r18, [sp, #-6*8]!

	ldp	r0, r1, [ap]
	ldp	r2, r3, [ap, #2*8]
	ldp	r10, r11, [bp, #2*8]
	stp	r19, r20, [sp, #2*8]

	ldp	r12, r13, [bp, #4*8]
	ldp	r14, r15, [bp, #6*8]
	stp	r21, r22, [sp, #4*8]

	umulh	r16, r6, r8
	mul	r17, r7, r8
	umulh	r8, r7, r8
	umulh	r18, r6, r9
	umulh	r19, r7, r9
	umulh	r20, r5, r9
	mul	r21, r6, r9
	mul	r9, r7, r9
	C (16, 17, 20, 21), (8, 18, 9), 19
	C Free: 22

	umulh	r22, r4, r10
	adds	r16, r16, r17
	umulh	r17, r5, r10
	adcs	r8, r8, r18
	umulh	r18, r6, r10
	cinc	r19, r19, cs
	adds	r16, r16, r20
	umulh	r20, r7, r10
	adcs	r8, r8, r9
	mul	r9, r5, r10
	cinc	r19, r19, cs
	adds	r16, r16, r21
	mul	r21, r6, r10
	mul	r10, r7, r10

	adcs	r8, r8, r17
	umulh	r17, r3, r11
	adcs	r19, r19, r18
	umulh	r18, r4, r11
	cinc	r20, r20, cs
	adds	r16, r16, r22
	umulh	r22, r5, r11
	adcs	r8, r8, r21
	umulh	r21, r6, r11
	adcs	r19, r19, r10
	umulh	r10, r7, r11
	cinc	r20, r20, cs
	adds	r16, r16, r9
	mul	r9, r4, r11
	adcs	r8, r8, r18
	mul	r18, r5, r11
	adcs	r19, r19, r22
	mul	r22, r6, r11
	mul	r11, r7, r11

	adcs	r20, r20, r21
	umulh	r21, r2, r12
	cinc	r10, r10, cs
	adds	r16, r16, r17
	umulh	r17, r3, r12
	adcs	r8, r8, r18
	umulh	r18, r4, r12
	adcs	r19, r19, r22
	umulh	r22, r5, r12
	adcs	r20, r20, r11
	umulh	r11, r6, r12
	cinc	r10, r10, cs
	adds	r16, r16, r9
	umulh	r9, r7, r12
	adcs	r8, r8, r17
	mul	r17, r3, r12
	adcs	r19, r19, r18
	mul	r18, r4, r12
	adcs	r20, r20, r22
	mul	r22, r5, r12
	adcs	r10, r10, r11
	mul	r11, r6, r12
	mul	r12, r7, r12

	cinc	r9, r9, cs
	adds	r16, r16, r21
	umulh	r21, r1, r13
	adcs	r8, r8, r18
	umulh	r18, r2, r13
	adcs	r19, r19, r22
	umulh	r22, r3, r13
	adcs	r20, r20, r11
	umulh	r11, r4, r13
	adcs	r10, r10, r12
	umulh	r12, r5, r13
	cinc	r9, r9, cs
	adds	r16, r16, r17
	umulh	r17, r6, r13
	adcs	r8, r8, r18
	umulh	r18, r7, r13
	adcs	r19, r19, r22
	mul	r22, r2, r13
	adcs	r20, r20, r11
	mul	r11, r3, r13
	adcs	r10, r10, r12
	mul	r12, r4, r13
	adcs	r9, r9, r17
	mul	r17, r5, r13
	cinc	r18, r18, cs
	adds	r16, r16, r21
	mul	r21, r6, r13
	mul	r13, r7, r13

	adcs	r8, r8, r11
	umulh	r11, r0, r14
	adcs	r19, r19, r12
	umulh	r12, r1, r14
	adcs	r20, r20, r17
	umulh	r17, r2, r14
	adcs	r10, r10, r21
	umulh	r21, r3, r14
	adcs	r9, r9, r13
	umulh	r13, r4, r14
	cinc	r18, r18, cs
	adds	r16, r16, r22
	umulh	r22, r5, r14
	adcs	r8, r8, r12
	umulh	r12, r6, r14
	adcs	r19, r19, r17
	umulh	r17, r7, r14
	adcs	r20, r20, r21
	mul	r21, r1, r14
	adcs	r10, r10, r13
	mul	r13, r2, r14
	adcs	r9, r9, r22
	mul	r22, r3, r14
	adcs	r18, r18, r12
	mul	r12, r4, r14
	cinc	r17, r17, cs
	adds	r16, r16, r11
	mul	r11, r5, r14
	adcs	r8, r8, r13
	mul	r13, r6, r14
	mul	r14, r7, r14

	adcs	r19, r19, r22
	umulh	r22, r0, r15
	adcs	r20, r20, r12
	umulh	r12, r1, r15
	adcs	r10, r10, r11
	umulh	r11, r2, r15
	adcs	r9, r9, r13
	umulh	r13, r3, r15
	adcs	r18, r18, r14
	umulh	r14, r4, r15
	cinc	r17, r17, cs
	adds	r16, r16, r21
	umulh	r21, r5, r15
	adcs	r8, r8, r22
	umulh	r22, r6, r15
	adcs	r19, r19, r12
	umulh	r12, r7, r15
	mul	r0, r0, r15
	adcs	r20, r20, r11
	mul	r1, r1, r15
	adcs	r10, r10, r13
	mul	r2, r2, r15
	adcs	r9, r9, r14
	mul	r3, r3, r15
	adcs	r18, r18, r21
	mul	r4, r4, r15
	adcs	r17, r17, r22
	mul	r5, r5, r15
	cinc	r12, r12, cs
	ldp	r21, r22, [sp, #4*8]
	mul	r6, r6, r15
	adds	r16, r16, r0
	mul	r7, r7, r15
	adcs	r8, r8, r1
	adcs	r0, r19, r2
	adcs	r1, r20, r3
	stp	r8, r0, [rp]
	ldp	r19, r20, [sp, #2*8]
	adcs	r10, r10, r4
	adcs	r9, r9, r5
	stp	r1, r10, [rp, #2*8]
	adcs	r6, r18, r6
	adcs	r7, r17, r7
	ldp	r17, r18, [sp], #6*8
	stp	r9, r6, [rp, #4*8]
	cinc	r12, r12, cs
	stp	r7, r12, [rp, #6*8]
	mov	res, r16
	ret
EPILOGUE()
