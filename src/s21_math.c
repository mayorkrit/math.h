#include "./s21_math.h"

int s21_abs(int x) {
    x = (x > 0) ? x : -x;  //  если x > 0 то х, иначе -х
    return x;
}


long double s21_acos(double x) {
    long double result = S21_NAN;

    if (s21_fabs(x) == 0) {
        result = 1.5707963;
    } else if (s21_fabs(x) < 1) {
        if (s21_fabs(x) == (double)(s21_SQRT2 / 2)) {
            result = (x > 0) ? 0.7853981633974483 : 2.3561944901923448;
        } else {
            result = s21_atan(s21_sqrt(1 - x * x) / x);
            if (-1 < x && x < 0) {
                result += S21_PI;
            }
        }
    } else if (s21_fabs(x) == 1) {
        result = (x > 0) ? 0 : S21_PI;
    }
    return result;
}

long double s21_asin(double x) {
    long double result = 0;
    if (s21_fabs(x) < 1) {
        if (s21_fabs(x) == (double)(s21_SQRT2 / 2)) {
            result = 0.7853981634;
            if (x < 0) {
                result = -result;
            }
        } else {
            result = s21_atan(x / (s21_sqrt(1 - x * x)));
        }
    } else if (s21_fabs(x) == 1) {
            result = 1.570796;
        if (x < 0) {
            result = -result;
        }
    } else {
        result = S21_NAN;
    }

    return result;
}



long double s21_atan(double x) {
    long double result = 0.0;

    if (s21_fabs(x) == S21_HUGE_VAL) {
        result = S21_PI / 2;
        if (x < 0) {
            result = -result;
        }
    } else if (s21_fabs(x) == 1) {
        result = 0.7853981634;
        if (x < 0) {
            result = -result;
        }
    } else if (s21_fabs(x) < 1) {
        for (int i = 0; i < 500; i++) {
            result += s21_pow(-1, i) * s21_pow(x, (1 + (2 * i))) / (1 + (2 * i));
        }
    } else {
        for (int i = 0; i < 500; i++) {
            result += s21_pow(-1, i) * s21_pow(x, (-1 - (2 * i))) / (1 + (2 * i));
        }
        result = (S21_PI * s21_fabs(x)) / (2 * x) - result;
    }

    return result;
}


long double s21_ceil(double x) {
    long double res = x;

    res = (long long)x;
    if (x - res >= 1.0 || x - res <= -1.0) {
        res = x;
    } else {
        if (res < x) res += 1.0;
    }
    return (is_nan(x) || is_inf(x)) ? x : res;
}

long double s21_cos(double x) {
    // можно сделать через сдвиг и деление пи пополам
    int znak = 1;
    x = s21_fmod(x, 2 * S21_PI);
    if (x > (S21_PI/2.0) && x <= S21_PI) {
        x = S21_PI - x;
        znak = -znak;
    } else if (x > S21_PI && x <= S21_PI * 3.0 / 2.0) {
        x -= S21_PI;
        znak = -znak;
    } else if (x > (S21_PI * 3.0) / 2.0 && x <= 2.0 * S21_PI) {
        x = 2 * S21_PI - x;
    }

    long double sum = 1.0;
    long double tailor = 1.0;
    int i;

    for (i = 1; s21_fabs(tailor/sum) > 1e-100; i++) {
        tailor = (-tailor * x * x) / ((2.0 * i - 1) * (2.0 * i));
        sum += tailor;
    }
    return sum * znak;
}

long double s21_exp(double x) {
    long double r = 1, y = 1;
    long double i = 1;
    int z = 0;
    if (x < 0) {
        x *= -1;
        z = 1;
    }
    while (s21_fabs(r) > 1e-17) {
    r *= x / i;
        i++;
        y += r;
        if (y > 1.7976931348623157e308) {
            y = 1.0 / 0.0;
            break;
        }
    }
    if (z == 1) {
        if (y > 1.7976931348623157e308)
            y = 0;
        else
            y = 1. / y;
    }
    if (y > 1.7976931348623157e308)
        y = 1.0 / 0.0;
    return y;
}

long double s21_fabs(double x) {
    x = (x > 0) ? x : -x;  //  если x > 0 то х, иначе -х
    return x;
}

long double s21_floor(double x) {
    long long res = x;
    if (is_nan(x) || is_inf(x) || x == 0) return x;
    if (x < 0.) {
        if (res > x) {
            res--;
        }
    }
    return x == 0 ? 0: res;
}

long double s21_log(double x) {
    if (!x || s21_fabs(x) < EPS) return -S21_HUGE_VAL;
    if (is_inf(x) == 1) return S21_HUGE_VAL;
    if (x == 2) return S21_LN2;
    if (s21_fabs(x - 1.0) < EPS / 10) return 0.0L;
    if (s21_fabs(x - S21_ME) <= 1e-15) return 1.0L;
    if (__builtin_signbit(x) || is_nan(x)) return S21_NAN;
    long double res = 0.0L;
    if (x > 0 && x <= 2) {
        long double num = (long double)x - 1.0L;
        long double sign = 1.0L;
        long double delim = 1.0L;
        long double ch = num;
        while (s21_fabs(ch) > 1E-50) {  // 0 < x < 2
            res += ch / delim * sign;
            ch *= num;
            sign *= -1.0L;
            delim += 1.0L;
        }
    } else {
        double how_much = 0;
        while (x >= 2.0) {
            x /= 2.0L;
            how_much++;
        }
        res = how_much * S21_LN2 + s21_log(x);
    }
    return res;
}


long double s21_sin(double x) {
    int znak = 1;
    x = s21_fmod(x, 2 * S21_PI);

    if (x > (S21_PI/2.0) && x <= S21_PI) {
        x = S21_PI - x;
    } else if (x > S21_PI && x <= S21_PI * 3.0 / 2.0) {
        x = (x - S21_PI);
        znak = -znak;
    } else if (x > (S21_PI * 3.0) / 2.0 && x <= 2.0 * S21_PI) {
        x = 2 * S21_PI - x;
        znak = -znak;
    }

    long double sum = (long double)x;
    long double tailor = (long double)x;
    int i;

    for (i = 1; s21_fabs(tailor/sum) > 1e-100; i++) {
        tailor = (-tailor * x * x) / ((2.0 * i + 1) * (2.0 * i));
        sum += tailor;
    }
    return sum * znak;
}

long double s21_sqrt(double x) {
    if (__builtin_signbit(x)) {
        return S21_NAN;
    }
    if (is_inf(x)) {
        return S21_HUGE_VAL;
    }
    return s21_pow(x, 0.5);
}

long double s21_tan(double x) {
    if (is_nan(x) || !is_finite(x)) {
        return S21_NAN;
    }

    return (s21_sin(x) / s21_cos(x));
}

long double s21_pow(double base, double exp) {
    long double res;

    if (base == 1.0 || exp == 0.0) {
        res = 1.0;
    } else if (base < 0.0 && s21_fmod(exp, 1.0) != 0) {
        res = -S21_NAN;
    } else {
        res = s21_exp(exp * s21_log(s21_fabs(base)));
    }
    if (base < 0.0 && s21_fmod(exp, 2.0) != 0.0) {
        res *= -1;
    }
    if (base == 0.0 && exp < 0.0) res = S21_HUGE_VAL;
    if (base == -1 && exp == -S21_HUGE_VAL) res = 1.0;
    return res;
}

long double s21_fmod(double x, double y) {
    if (is_nan(x) || is_nan(y) || is_inf(x) || s21_fabs(y) < EPS || (is_inf(x) && is_inf(y))) {
        return S21_NAN;
    }
    if (s21_fabs(x) < EPS) return 0;
    if (is_inf(y)) return x;
    double a, b, result;
    a = x/y;
    if (x < 0 && y < 0) {
        b = s21_abs(a);
    } else if (x < 0 || y < 0) {
        b = (-1) * s21_abs(a);
    } else {
        b = s21_abs(a);
    }
    result = x - b * y;
    if (y == 0) {
        result = 0 / 0.0;
    }
    return result;
}
