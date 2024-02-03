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

.macro	m4 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	mov	\scr1, \res_offset*8(\res)
	adcx	\scr2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	adcx	\zero, \r3
.endm

.macro	am4 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, r3, r4, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r4, \scr
	adcx	\r4, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r4
	adcx	\scr, \r1
	adox	\r0, \r1
	mulx	(2+\ap_offset)*8(\ap), \r0, \scr
	adcx	\r4, \r2
	adox	\r0, \r2
	mulx	(3+\ap_offset)*8(\ap), \r0, \r4
	adcx	\scr, \r3
	adox	\r0, \r3
	adcx	\zero, \r4
	adox	\zero, \r4
.endm

.macro	m4_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, r2, r3, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	mulx	(2+\ap_offset)*8(\ap), \scr1, \r2
	mulx	(3+\ap_offset)*8(\ap), \scr2, \r3
	adcx	\scr1, \r1
	adcx	\scr2, \r2
	adcx	\zero, \r3
.endm

.macro	m2_chain res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, ip1, ip2, r0, r1, scr1, scr2, zero
	mulx	(0+\ap_offset)*8(\ap), \scr1, \r0
	adcx	\ip1, \scr1
	mov	\scr1, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \scr2, \r1
	adcx	\scr2, \r0
	adox	\ip2, \r0
	adcx	\zero, \r1
.endm

.macro	am2 res=%rdi, res_offset=0, ap=%rsi, ap_offset=0, r0, r1, r2, scr, zero
	mulx	(0+\ap_offset)*8(\ap), \r2, \scr
	adcx	\r2, \r0
	mov	\r0, \res_offset*8(\res)
	mulx	(1+\ap_offset)*8(\ap), \r0, \r2
	adcx	\scr, \r1
	adox	\r0, \r1
	adcx	\zero, \r2
	adox	\zero, \r2
.endm

.global	FUNC(flint_mpn_mul_10_2)
.p2align	4, 0x90
TYPE(flint_mpn_mul_10_2)

FUNC(flint_mpn_mul_10_2):
	.cfi_startproc
	mov	0*8(%rdx), %rcx
	mov	1*8(%rdx), %r8
	push	%rbx
	push	%rbp
	push	%r12

	mov	%rcx, %rdx
	xor	%r12d, %r12d
	m4	%rdi, 0, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r8, %rdx
	am4	%rdi, 1, %rsi, 0, %r9, %r10, %r11, %rbx, %rax, %rbp, %r12
	mov	%r10, 2*8(%rdi)
	mov	%r11, 3*8(%rdi)

	mov	%rcx, %rdx
	m4_chain	%rdi, 4, %rsi, 4, %rbx, %rax, %r9, %r10, %r11, %rax, %rbp, %rbx, %r12
	mov	%r8, %rdx
	am4	%rdi, 5, %rsi, 4, %r9, %r10, %r11, %rax, %rbx, %rbp, %r12
	mov	%r10, 6*8(%rdi)
	mov	%r11, 7*8(%rdi)

	mov	%rcx, %rdx
	m2_chain	%rdi, 8, %rsi, 8, %rax, %rbx, %r9, %r10, %rbp, %rax, %r12
	mov	%r8, %rdx
	am2	%rdi, 9, %rsi, 8, %r9, %r10, %rax, %rbp, %r12
	mov	%r10, 10*8(%rdi)
	mov	%rax, 11*8(%rdi)
	pop	%r12
	pop	%rbp
	pop	%rbx

	ret
.flint_mpn_mul_10_2_end:
SIZE(flint_mpn_mul_10_2, .flint_mpn_mul_10_2_end)
.cfi_endproc
