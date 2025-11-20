#ifndef _generatormt_h_
#define _generatormt_h_
#include <cmath>
/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/
//extern "C"{
/// initialize with a long seed
void init_genrand(unsigned long s);

/// nitialize with an array
void init_by_array(unsigned long init_key[], int key_length);

unsigned long genrand_int32(void);

/// generates a random number on [0,0x7fffffff]-interval 
long genrand_int31(void);

/// generates a random number on [0,1]-real-interval 
double genrand_real1(void);

/// generates a random number on [0,1)-real-interval 
double genrand_real2(void);

/// generates a random number on (0,1)-real-interval 
double genrand_real3(void);

/// generates a random number on [0,1) with 53-bit resolution
double genrand_res53(void); 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/// 53-bit resolution number from [0,1)
inline double frandom () { return genrand_res53();}

/// 32-bit resolution number from [0,1]
inline double frandom11 () {  return genrand_real1();}

/// 32-bit resolution number from [0,1)
inline double frandom10 (){ return genrand_real2();}

/// 32-bit resolution number from (0,1]
inline double frandom01 (){ return 1.0-genrand_real2();}

/// 32-bit resolution number from (0,1)
inline double frandom00 (){   return genrand_real3(); }

/// x from [0,1) with probability porportional to  x*x
inline double frandom_sqr () { return pow(frandom(),1.0/3);}

/// x from [0,1) with probability porportional to  x^n
inline double frandom_pow (int n) { return pow(frandom(),1.0/(n+1));}

/// initialize the random seed_generator
void frandom_init(int option);

/// read genrand state from file
void genrand_read_state(const char* file="random_seed");

/// save genrand state in file
void genrand_write_state(const char* file="random_seed");

#endif
