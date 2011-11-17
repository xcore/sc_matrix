/*
 * Pointer wrangling to cheat in XC
 *
 * Copyright (C) 2011 Steve Kerrison <github@stevekerrison.com>
 *
 * This software is freely distributable under a derivative of the
 * University of Illinois/NCSA Open Source License posted in
 * LICENSE.txt and at <http://github.xcore.com/>
 */

#ifndef XC_POINTER_H_
#define XC_POINTER_H_

#ifdef __XC__
#define NULLABLE ?
#else
#define NULLABLE
#endif

#define POINTER(dst,src) asm("mov %0,%1":"=r"(dst):"r"(src))

unsigned int pointer_int(int NULLABLE p[]);
unsigned int pointer_short(short NULLABLE p[]);
unsigned int pointer_char(char NULLABLE p[]);

#endif //XC_POINTER_H_
