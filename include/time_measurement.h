/*
 * Copyright (c) 2018 Maxime Schmitt <max.schmitt@unistra.fr>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __TIME_MEASURING_H
#define __TIME_MEASURING_H

#include <time.h>

#ifdef CLOCK_MONOTONIC_RAW
#define MEASURING_CLOCK CLOCK_MONOTONIC_RAW
#else
#define MEASURING_CLOCK CLOCK_MONOTONIC
#endif

typedef struct timespec time_measure;

__attribute__((always_inline)) static inline void
get_current_time(time_measure *time) {
  clock_gettime(MEASURING_CLOCK, time);
}

__attribute__((const)) static inline double
measuring_difftime(time_measure t0, time_measure t1) {
  double secdiff = difftime(t1.tv_sec, t0.tv_sec);
  if (t1.tv_nsec < t0.tv_nsec) {
    long val = 1000000000l - t0.tv_nsec + t1.tv_nsec;
    secdiff += (double)val / 1e9 - 1.;
  } else {
    long val = t1.tv_nsec - t0.tv_nsec;
    secdiff += (double)val / 1e9;
  }
  return secdiff;
}

#endif // __TIME_MEASURING_H
