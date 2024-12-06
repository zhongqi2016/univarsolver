"""This module implements correctly rounded interval functions
"""
import decimal as dec
import math

import interval_arithmetics as ia


# -------------- Factorial function --------

def factorial(n):
    """
    Computes interval extension for a factorial function
    
    Parameters
    ----------
    n : Interval 
        
    Returns:
    --------
    The interval that enclose the range of the factorial function
    """
    a = ia.c_one
    ia._set_rounding_mode_floor()
    for i in range(2, int(n.a) + 1):
        a = a * dec.Decimal(i)
    b = ia.c_one
    ia._set_rounding_mode_ceil()
    for i in range(2, int(n.b) + 1):
        b = b * dec.Decimal(i)
    return ia.Interval(a, b)


# --------------- Exponential ---------------

#     Number of Taylor's series terms for exponential
#     Notice, that this number will be multiplied by two to bound the exp(x) between the consecutive sums
#     1 + x + ... + x^(2n - 1) < exp(x) < 1 + x + ... + x^(2n)
exp_taylor_terms_number = 5


def get_exp_taylor_terms_number():
    """
    Retrieves the number of Taylor's series terms
    
    Notice, that this number will be multiplied by two to bound the exp(x) between the consecutive sums
    1 + x + ... + x^(2n - 1) < exp(x) < 1 + x + ... + x^(2n)
    
    Returns:
    ----------
    the  value of exp_taylor_terms_number    
    """
    global exp_taylor_terms_number
    return exp_taylor_terms_number


def set_exp_taylor_terms_number(number):
    """
    Sets new number of Taylor's series members to approximate the exponential. 
    
    Notice, that this number will be multiplied by two to bound the exp(x) between the consecutive sums
    1 + x + ... + x^(2n - 1) < exp(x) < 1 + x + ... + x^(2n)
    
    Parameters:
    ----------
    number : new value of exp_taylor_terms_number
    
    Returns:
    ----------
    the previous value of exp_taylor_terms_number
    """
    global exp_taylor_terms_number
    old = exp_taylor_terms_number
    exp_taylor_terms_number = number
    return old


# Exponential for x in a range [-1, 0], lower bound, x is a dot interval
def _exp_f_m1_to_0_lb(x):
    global exp_taylor_terms_number
    n = 2 * exp_taylor_terms_number
    s = ia.Interval(ia.c_one, ia.c_one)
    for i in range(1, n):
        ii = ia.Interval(dec.Decimal(i), dec.Decimal(i))
        term = (x ** i) / factorial(ii)
        s += term
    return s.a


# Exponential for x in a range [-1, 0], upper bound, x is a dot interval
def _exp_f_m1_to_0_ub(x):
    global exp_taylor_terms_number
    n = 2 * exp_taylor_terms_number + 1
    s = ia.Interval(ia.c_one, ia.c_one)
    for i in range(1, n):
        ii = ia.Interval(dec.Decimal(i), dec.Decimal(i))
        term = (x ** i) / factorial(ii)
        s += term
    return s.b


# Exponential for x in a range [-inf, -1), lower bound, x is a dot interval
def _exp_f_minf_to_m1_lb(x):
    z = -(x.a.to_integral(rounding=dec.ROUND_FLOOR))
    ivz = ia.Interval(z, z)
    ivf = x / ivz
    etv = _exp_f_m1_to_0_lb(ivf)
    ietv = ia.Interval(etv, etv)
    iev = ietv ** int(z)
    return iev.a


# Exponential for x in a range [-inf, -1), upper bound, x is a dot interval
def _exp_f_minf_to_m1_ub(x):
    z = -(x.a.to_integral(rounding=dec.ROUND_FLOOR))
    ivz = ia.Interval(z, z)
    ivf = x / ivz
    etv = _exp_f_m1_to_0_ub(ivf)
    ietv = ia.Interval(etv, etv)
    iev = ietv ** int(z)
    return iev.b


# Exponential for x in a range (0, +inf], lower bound, x is a dot interval
def _exp_f_0_to_inf_lb(x):
    if x.a > ia.c_one:
        irev = _exp_f_minf_to_m1_ub(-x)
    else:
        irev = _exp_f_m1_to_0_ub(-x)
    ione = ia.Interval(dec.Decimal('1'), dec.Decimal('1'))
    iev = ione / irev
    return iev.a


# Exponential for x in a range (0, +inf], upper bound, x is a dot interval
def _exp_f_0_to_inf_ub(x):
    if x.a > ia.c_one:
        irev = _exp_f_minf_to_m1_lb(-x)
    else:
        irev = _exp_f_m1_to_0_lb(-x)
    ione = ia.Interval(ia.c_one, ia.c_one)
    iev = ione / irev
    return iev.b


# Exponential lower bound for x, x is a dot interval
def _exp_lb(x):
    if x.a == ia.c_minf:
        return ia.c_zero
    elif x.a < ia.c_mone:
        return _exp_f_minf_to_m1_lb(x)
    elif ia.c_mone <= x.a < ia.c_zero:
        return _exp_f_m1_to_0_lb(x)
    elif x.a == ia.c_zero:
        return ia.c_one
    elif x.a == ia.c_inf:
        return ia.c_inf
    else:
        return _exp_f_0_to_inf_lb(x)


# Exponential upper bound for x, x is a dot interval
def _exp_ub(x):
    if x.b == ia.c_minf:
        return ia.c_zero
    elif x.b < ia.c_mone:
        return _exp_f_minf_to_m1_ub(x)
    elif ia.c_mone <= x.b < ia.c_zero:
        return _exp_f_m1_to_0_ub(x)
    elif x.b == ia.c_zero:
        return ia.c_one
    elif x.b == ia.c_inf:
        return ia.c_inf
    else:
        return _exp_f_0_to_inf_ub(x)


def exp(x):
    """
    Computes reliable bounds for exponential.
    
    Parameters:
    -----------
    x : an interval
    
    Returns:
    --------
    Enclosing interval for exp(x)
    """
    a = _exp_lb(x)
    b = _exp_ub(x)
    return ia.Interval(a, b)


# ------------ Natural logarithm ------------------------
#  Natural logarithm approximation based on Taylor series of ln ((1+x)/(1-x))
#

# Number of Taylor expansion terms

log_taylor_terms_number = 5


def get_log_taylor_terms_number():
    """
    Retrieves number of Taylor's series members to approximate the natural logarithm.
    
    Notice, that this number will be multiplied by two to bound the ln(1 + x) between the consecutive sums
    x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n < ln(1 + x) <  x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n + (x^(2n+1))/(2n+1)
    
    Returns:
    ----------
    the  value of log_taylor_terms_number    
    """
    global log_taylor_terms_number
    return log_taylor_terms_number


def set_log_taylor_terms_number(number):
    """
    Sets new number of Taylor's series members to approximate the natural logarithm. 
    
    Notice, that this number will be multiplied by two to bound the ln(1 + x/1 - x) 
    
    Parameters:
    ----------
    number : new value of ln_taylor_terms_number
    
    Returns:
    ----------
    the previous value of ln_taylor_terms_number
    """
    global log_taylor_terms_number
    old = log_taylor_terms_number
    log_taylor_terms_number = number
    return old


# Point x from 1 to infinity
def _log_f_1_to_inf(xp):
    x = ia.Interval(xp, xp)
    inul = ia.Interval(ia.c_zero, ia.c_zero)
    ione = ia.Interval(ia.c_one, ia.c_one)
    itwo = ia.Interval(ia.c_two, ia.c_two)
    z = (x - ione) / (x + ione)
    s = inul
    zn = z
    zq = pow(z, 2)
    global log_taylor_terms_number
    k = log_taylor_terms_number
    for i in range(0, k + 1):
        ii = ia.Interval(dec.Decimal(2 * i + 1), dec.Decimal(2 * i + 1))
        s += zn / ii
        zn *= zq
    t = ia.Interval(dec.Decimal(2 * k + 3), dec.Decimal(2 * k + 3)) * (ione - zq)
    su = s + zn / t
    s *= itwo
    su *= itwo
    iret = ia.Interval(s.a, su.b)
    return iret


# Point x from 0 to 1
def _log_f_0_to_1(xp):
    x = ia.Interval(xp, xp)
    inul = ia.Interval(ia.c_zero, ia.c_zero)
    ione = ia.Interval(ia.c_one, ia.c_one)
    itwo = ia.Interval(ia.c_two, ia.c_two)
    z = (x - ione) / (x + ione)
    s = inul
    zn = z
    zq = pow(z, 2)
    global log_taylor_terms_number
    k = log_taylor_terms_number
    for i in range(0, k + 1):
        ii = ia.Interval(dec.Decimal(2 * i + 1), dec.Decimal(2 * i + 1))
        s += zn / ii
        zn *= zq
    t = ia.Interval(dec.Decimal(2 * k + 3), dec.Decimal(2 * k + 3)) * (ione - zq)
    su = s + zn / t
    s *= itwo
    su *= itwo
    iret = ia.Interval(su.a, s.b)
    return iret


# logarithm for a point
def _log_point(xp):
    if xp == ia.c_zero:
        ival = ia.Interval(ia.c_minf, ia.c_minf)
    elif xp == ia.c_one:
        ival = ia.Interval(ia.c_zero, ia.c_zero)
    elif xp < ia.c_one:
        ival = _log_f_0_to_1(xp)
    elif xp == ia.c_inf:
        ival = ia.Interval(ia.c_inf, ia.c_inf)
    else:
        ival = _log_f_1_to_inf(xp)
    return ival


def log(x: ia.Interval):
    """
    Computes reliable bounds for natural logarithm.
    
    Parameters:
    -----------
    x : an interval
    
    Returns:
    --------
    Enclosing interval for ln(x)
    """
    loga = _log_point(x.a)
    logb = _log_point(x.b)
    iret = ia.Interval(loga.a, logb.b)
    return iret


# ------------ Natural logarithm ------------------------

#     Number of Taylor's series terms for ln(1 + x)
#     Notice, that this number will be multiplied by two to bound the ln(1 + x) between the consecutive sums
#     x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n < ln(1 + x) <  x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n + (x^(2n+1))/(2n+1)


# # Lower bound for ln(x), x in (1, 2]
# def _log_f_1_to_2_lb(x):
#     global log_taylor_terms_number
#     n = 2 * log_taylor_terms_number
#     y = x - ia.Interval(dec.Decimal(1), dec.Decimal(1))
#     s = y
#     for i in range(2, n + 1):
#         ii = ia.Interval(dec.Decimal(i), dec.Decimal(i))
#         term = (y ** i) / ii
#         if i % 2 == 0:
#             s -= term
#         else: 
#             s += term
#     return s.a

# # Upper bound for ln(x), x in (1, 2]
# def _log_f_1_to_2_ub(x):
#     global log_taylor_terms_number
#     n = 2 * log_taylor_terms_number + 1
#     y = x - ia.Interval(dec.Decimal(1), dec.Decimal(1))
#     s = y
#     for i in range(2, n + 1):
#         ii = ia.Interval(dec.Decimal(i), dec.Decimal(i))
#         term = (y ** i) / ii
#         if i % 2 == 0:
#             s -= term
#         else: 
#             s += term
#     return s.b

# # Lower for ln(x)
# def _log_lb(x):
#     return _log_f_1_to_2_lb(x)

# # Upper for ln(x)
# def _log_ub(x):
#     return _log_f_1_to_2_ub(x)


# def log(x):
#     """
#     Computes reliable bounds for logarithm.

#     Parameters:
#     -----------
#     x : an interval

#     Returns:
#     --------
#     Enclosing interval for ln(x)
#     """
#     a = _log_lb(ia.Interval(x.a, x.a))
#     b = _log_ub(ia.Interval(x.b, x.b))
#     return ia.Interval(a, b)


# ------------ sin ------------------------
#  Natural sin approximation based on Taylor series
#
pi = dec.Decimal(math.pi)
pi2 = dec.Decimal(2 * math.pi)
pi05 = dec.Decimal(math.pi / 2)
# Number of Taylor expansion terms

sin_taylor_terms_number = 9


def get_sin_taylor_terms_number():
    """
    Retrieves number of Taylor's series members to approximate the natural logarithm.

    Notice, that this number will be multiplied by two to bound the ln(1 + x) between the consecutive sums
    x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n < ln(1 + x) <  x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n + (x^(2n+1))/(2n+1)

    Returns:
    ----------
    the  value of log_taylor_terms_number
    """
    global sin_taylor_terms_number
    return sin_taylor_terms_number


def set_sin_taylor_terms_number(number):
    """
    Sets new number of Taylor's series members to approximate the natural logarithm.

    Notice, that this number will be multiplied by two to bound the ln(1 + x/1 - x)

    Parameters:
    ----------
    number : new value of ln_taylor_terms_number

    Returns:
    ----------
    the previous value of ln_taylor_terms_number
    """
    global sin_taylor_terms_number
    old = sin_taylor_terms_number
    sin_taylor_terms_number = number
    return old


def _sin_lb(x):
    ia._set_rounding_mode_floor()
    global sin_taylor_terms_number
    n = sin_taylor_terms_number
    s = ia.Interval(x.a, x.b)
    flag = False
    if n % 2 == 1:
        n += 1
    for i in range(1, n):
        ii = ia.Interval(dec.Decimal(2 * i + 1), dec.Decimal(2 * i + 1))
        term = (x ** (2 * i + 1)) / factorial(ii)
        if flag:
            s += term
            flag = False
        else:
            s -= term
            flag = True
    return s.a


def _sin_ub(x):
    ia._set_rounding_mode_ceil()
    global sin_taylor_terms_number
    n = sin_taylor_terms_number
    s = ia.Interval(x.a, x.b)
    flag = False
    if n % 2 == 1:
        n += 1
    for i in range(1, n + 1):
        ii = ia.Interval(dec.Decimal(2 * i + 1), dec.Decimal(2 * i + 1))
        term = (x ** (2 * i + 1)) / factorial(ii)
        if flag:
            s += term
            flag = False
        else:
            s -= term
            flag = True
    return s.b


def _sin(xp):
    xp = xp % pi2
    if xp < 0:
        xp2 = xp + pi2
        if xp2 < -xp:
            xp = xp2
    else:
        xp2 = xp - pi2
        if -xp2 < xp:
            xp = xp2
    x = ia.Interval(xp, xp)
    return ia.Interval(_sin_lb(x), _sin_ub(x))


def sin(x_input: ia.Interval):
    a = _sin(x_input.a)
    b = _sin(x_input.b)
    if ((x_input.a - pi05) / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_CEILING) <= (
            (x_input.b - pi05) / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_FLOOR):
        max_val = dec.Decimal(1)
    else:
        max_val = max(a.b, b.b)
    if ((x_input.a + pi05) / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_CEILING) <= (
            (x_input.b + pi05) / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_FLOOR):
        min_val = dec.Decimal(-1)
    else:
        min_val = min(a.a, b.a)
    return ia.Interval(min_val, max_val)


# ------------ cos ------------------------
#  Natural logarithm approximation based on Taylor series
#

# Number of Taylor expansion terms

cos_taylor_terms_number = 9


def get_cos_taylor_terms_number():
    """
    Retrieves number of Taylor's series members to approximate the cos

    Notice, that this number will be multiplied by two to bound the ln(1 + x) between the consecutive sums
    x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n < ln(1 + x) <  x - (x^2)/2 + (x^3)/3 - ... - (x^2n)/2n + (x^(2n+1))/(2n+1)

    Returns:
    ----------
    the  value of cos_taylor_terms_number
    """
    global cos_taylor_terms_number
    return cos_taylor_terms_number


def set_cos_taylor_terms_number(number):
    """
    Sets new number of Taylor's series members to approximate the natural logarithm.

    Notice, that this number will be multiplied by two to bound the ln(1 + x/1 - x)

    Parameters:
    ----------
    number : new value of ln_taylor_terms_number

    Returns:
    ----------
    the previous value of ln_taylor_terms_number
    """
    global cos_taylor_terms_number
    old = cos_taylor_terms_number
    cos_taylor_terms_number = number
    return old


def _cos_lb(x):
    ia._set_rounding_mode_floor()
    global cos_taylor_terms_number
    n = cos_taylor_terms_number
    s = ia.Interval(ia.c_one, ia.c_one)
    flag = False
    if n % 2 == 1:
        n += 1
    for i in range(1, n):
        ii = ia.Interval(dec.Decimal(2 * i), dec.Decimal(2 * i))
        term = (x ** (2 * i)) / factorial(ii)
        if flag:
            s += term
            flag = False
        else:
            s -= term
            flag = True
    return s.a


def _cos_ub(x):
    ia._set_rounding_mode_ceil()
    global cos_taylor_terms_number
    n = cos_taylor_terms_number
    s = ia.Interval(ia.c_one, ia.c_one)
    flag = False
    if n % 2 == 1:
        n -= 1
    for i in range(1, n + 1):
        ii = ia.Interval(dec.Decimal(2 * i), dec.Decimal(2 * i))
        term = (x ** (2 * i)) / factorial(ii)
        if flag:
            s += term
            flag = False
        else:
            s -= term
            flag = True
    return s.b


def _cos(xp):
    xp = xp % pi2
    if xp < 0:
        xp2 = xp + pi2
        if xp2 < -xp:
            xp = xp2
    else:
        xp2 = xp - pi2
        if -xp2 < xp:
            xp = xp2
    x = ia.Interval(xp, xp)
    return ia.Interval(_cos_lb(x), _cos_ub(x))


def cos(x_input: ia.Interval):
    a = _cos(x_input.a)
    b = _cos(x_input.b)
    if (x_input.a / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_CEILING) <= (
            x_input.b / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_FLOOR):
        max_val = dec.Decimal(1)
    else:
        max_val = max(a.b, b.b)
    if ((x_input.a - pi) / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_CEILING) <= (
            (x_input.b - pi) / pi2).quantize(dec.Decimal('1.'), rounding=dec.ROUND_FLOOR):
        min_val = dec.Decimal(-1)
    else:
        min_val = min(a.a, b.a)
    return ia.Interval(min_val, max_val)


def sqrt(x_input):
    if isinstance(x_input, ia.Interval):
        return x_input.sqrt()
    else:
        return math.sqrt(x_input)
