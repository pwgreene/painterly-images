#ifndef _TIMING_H_
#define _TIMING_H_

#include <Halide.h>
using Halide::Func;

#define N_TIMES 5

unsigned long millisecond_timer(void);


/**
 * Print the runtime and throughout and return the runtime.
 */
float profile(Func myFunc, int w, int h);


/**
 * Print the runtime and throughout and return the runtime.
 */
float profile(Func myFunc, int w, int h, int c);

#endif // _TIMING_H_
