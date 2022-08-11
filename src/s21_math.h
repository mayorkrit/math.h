#ifndef SRC_S21_MATH_H_
#define SRC_S21_MATH_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <stdbool.h>
#include <stdint.h>

// #define S21_HUGE_VAL (__builtin_inff())
#define S21_PI 3.14159265358979323846264338327950288
#define S21_ME 2.71828182845904523536028747135266250
#define S21_LN2 0.693147180559945309417232121458176568
#define S21_PI2 1.5707963267948966192313216916397514
#define s21_SQRT2 1.41421356237309504880

#define EPS 1e-6
// #define S21_NAN (__builtin_nanf(""))
// #define S21_INF (__builtin_inff())
// #define is_nan(x) __builtin_isnan(x)
// #define is_inf(x) __builtin_isinf(x)



#define S21_HUGE_VAL (__builtin_inff())
/**
  * returns float NaN value
  */
#define S21_NAN (__builtin_nanf(""))
/**
  * check for infinity: returns 1 if infinite, -1 if -infinite and 0 if finite
  */
#define is_inf(x) __builtin_isinf_sign(x)
/**
 * check for nan value: returns 1 if NaN and 0 if not 
 */
#define is_nan(x) __builtin_isnan(x)
#define is_finite(x) __builtin_isfinite(x)


int s21_abs(int x);
long double s21_fabs(double x);
long double s21_exp(double x);
long double s21_log(double x);
long double s21_pow(double base, double exp);
long double s21_sqrt(double x);
long double s21_sin(double x);
long double s21_cos(double x);
long double s21_tan(double x);
long double s21_fmod(double x, double y);
long double s21_floor(double x);
long double s21_ceil(double x);
long double s21_asin(double x);
long double s21_acos(double x);
long double s21_atan(double x);


#endif  //  SRC_S21_MATH_H_
