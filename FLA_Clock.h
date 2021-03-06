/*

libflame license
----------------

libflame is available as free software under the following "3-clause"
(also known as "new" or "modified") BSD license.

------
Copyright (C) 2014, The University of Texas at Austin

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
 - Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
 - Neither the name of The University of Texas nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
------

Select parts of libflame's f2c implementation were taken from:

  https://github.com/juanjosegarciaripoll/f2c

which uses the following license:

------
Copyright (C) 1990 - 1997 by AT&T, Lucent Technologies and Bellcore.

Permission to use, copy, modify, and distribute this software
and its documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appear in all
copies and that both that the copyright notice and this
permission notice and warranty disclaimer appear in supporting
documentation, and that the names of AT&T, Bell Laboratories,
Lucent or Bellcore or any of their entities not be used in
advertising or publicity pertaining to distribution of the
software without specific, written prior permission.

AT&T, Lucent and Bellcore disclaim all warranties with regard to
this software, including all implied warranties of
merchantability and fitness.  In no event shall AT&T, Lucent or
Bellcore be liable for any special, indirect or consequential
damages or any damages whatsoever resulting from loss of use,
data or profits, whether in an action of contract, negligence or
other tortious action, arising out of or in connection with the
use or performance of this software.
------

*/

#ifndef _FLA_CLOCK_H_
#define _FLA_CLOCK_H_

#if defined(__APPLE__) || defined(__MACH__)
#include <AvailabilityMacros.h>
  #include <mach/mach_time.h>
#else
  #include <time.h>
#endif

double FLA_Clock_helper( void );

// A global variable used when FLA_Clock_helper() is defined in terms of
// clock_gettime()/gettimeofday().
double gtod_ref_time_sec = 0.0;

double FLA_Clock( void )
{
	return FLA_Clock_helper();
}

#if defined(__APPLE__) || defined(__MACH__)

double FLA_Clock_helper()
{
    mach_timebase_info_data_t timebase;
    mach_timebase_info( &timebase );

    uint64_t nsec = mach_absolute_time();

    double the_time = (double) nsec * 1.0e-9 * timebase.numer / timebase.denom;

    if ( gtod_ref_time_sec == 0.0 )
        gtod_ref_time_sec = the_time;

    return the_time - gtod_ref_time_sec;
}

#else

double FLA_Clock_helper()
{
    double the_time, norm_sec;
    struct timespec ts;

    clock_gettime( CLOCK_MONOTONIC, &ts );

    if ( gtod_ref_time_sec == 0.0 )
        gtod_ref_time_sec = ( double ) ts.tv_sec;

    norm_sec = ( double ) ts.tv_sec - gtod_ref_time_sec;

    the_time = norm_sec + ts.tv_nsec * 1.0e-9;

    return the_time;
}

#endif

#endif // _FLA_CLOCK_H_
