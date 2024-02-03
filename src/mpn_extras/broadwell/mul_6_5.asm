#
#   Copyright (C) 2023 Albin Ahlbäck
#
#   This file is part of FLINT.
#
#   FLINT is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Lesser General Public License (LGPL) as published
#   by the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.  See <https://www.gnu.org/licenses/>.
#

include(`config.m4')dnl
dnl
.text

.macro	m6 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	mulx	(4+\ap_offset)*8(\ap), \scr1, \r4
	mulx	(5+\ap_offset)*8(\ap), \scr2, \r5
	adcx	\scr1, \r3
	adcx	\scr2, \r4
	adcx	\zero, \r5
.endm

.macro	am6 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, r5, r6, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r6, \scr
	adcx	\r6, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r6
	adcx	\scr, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r6, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \r6
	adcx	\scr, \r3
	adox	\r0, \r3
	mulx	(4+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r6, \r4
	adox	\r0, \r4
	mulx	(5+\ap_offset)*8(\ap), \r0, \r6
	adcx	\scr, \r5
	adox	\r0, \r5
	adcx	\zero, \r6
	adox	\zero, \r6
.endm

.global	FUNC(flint_mpn_mul_6_5)
.p2align	4, 0x90
TYPE(flint_mpn_mul_6_5)

FUNC(flint_mpn_mul_6_5):
	.cfi_startproc
	mov	%rdx, %rcx
	mov	0*8(%rcx), %rdx
	push	%rbx
	push	%rbp
	push	%r12
	push	%r13

	xor	%r13d, %r13d

	m6	%rdi, 0, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13

	mov	1*8(%rcx), %rdx
	am6	%rdi, 1, %rsi, 0, %r9, %r8, %rax, %r10, %r11, %rbx, %rbp, %r12, %r13
	mov	2*8(%rcx), %rdx
	am6	%rdi, 2, %rsi, 0, %r8, %rax, %r10, %r11, %rbx, %rbp, %r9, %r12, %r13
	mov	3*8(%rcx), %rdx
	am6	%rdi, 3, %rsi, 0, %rax, %r10, %r11, %rbx, %rbp, %r9, %r8, %r12, %r13
	mov	4*8(%rcx), %rdx
	am6	%rdi, 4, %rsi, 0, %r10, %r11, %rbx, %rbp, %r9, %r8, %rax, %r12, %r13

	mov	%r11, 5*8(%rdi)
	mov	%rbx, 6*8(%rdi)
	mov	%rbp, 7*8(%rdi)
	mov	%r9, 8*8(%rdi)
	mov	%r8, 9*8(%rdi)
	mov	%rax, 10*8(%rdi)

	pop	%r13
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mul_6_5_end:
SIZE(flint_mpn_mul_6_5, .flint_mpn_mul_6_5_end)
.cfi_endproc
