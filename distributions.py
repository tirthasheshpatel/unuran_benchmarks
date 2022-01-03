import numpy as np
from math import pi, sqrt, exp, erf, gamma
from scipy import special

class StdNorm:
    def pdf(self, x):
        return 1 / sqrt(2 * pi) * exp(-0.5 * x * x)

    def dpdf(self, x):
        return -1 / sqrt(2 * pi) * x * exp(-0.5 * x * x)

    def cdf(self, x):
        return (1 + erf(x / sqrt(2))) / 2

    def mode(self):
        return 0

    def __repr__(self):
        return "stdnormal"


class Gamma:
    def __init__(self, a):
        self.a = a

    def pdf(self, x):
        try:
            return x**(self.a - 1) * exp(-x) / gamma(self.a)
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def dpdf(self, x):
        try:
            first = (self.a - 1) * x ** (self.a - 2) * exp(-x) / gamma(self.a)
            second = -x**(self.a - 1) * exp(-x) / gamma(self.a)
            return first + second
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def cdf(self, x):
        return special.gammainc(self.a, x)

    def support(self):
        return 0, np.inf

    def mode(self):
        return self.a - 1 if self.a >= 1 else 0.0

    def __repr__(self):
        return "gamma"


class Beta:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def pdf(self, x):
        try:
            return (
                gamma(self.a + self.b)
                * x ** (self.a - 1)
                * (1 - x) ** (self.b - 1)
                / gamma(self.a)
                / gamma(self.b)
            )
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def dpdf(self, x):
        try:
            z = gamma(self.a) * gamma(self.b) / gamma(self.a + self.b)
            return (
                (self.a - 1) * x ** (self.a - 2) * (1 - x) ** (self.b - 1)
                - x ** (self.a - 1) * (self.b - 1) * (1 - x) ** (self.b - 2)
            ) / z
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def cdf(self, x):
        return special.btdtr(self.a, self.b, x)

    def support(self):
        return 0, 1

    def mode(self):
        if self.a > 1 and self.b > 1:
            return (self.a - 1) / (self.a + self.b - 2)
        if self.a == 1 and self.b == 1:
            return 0.5
        if self.a <= 1 and self.b > 1:
            return 0
        return 1

    def __repr__(self):
        return "beta"


class GenNorm:
    def __init__(self, beta):
        self.beta = beta

    def pdf(self, x):
        try:
            z = 2 * gamma(1 / self.beta) / self.beta
            return exp(-abs(x) ** self.beta) / z
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def dpdf(self, x):
        try:
            z = 2 * gamma(1 / self.beta) / self.beta
            if x > 0:
                return (exp(-(x ** self.beta)) / z) * (-self.beta * x ** (self.beta - 1))
            return (exp(-((-x) ** self.beta)) / z) * (self.beta * (-x) ** (self.beta - 1))
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def cdf(self, x):
        c = 0.5 * (abs(x) / x)
        return (0.5 + c) - c * special.gammaincc(1 / self.beta, abs(x) ** self.beta)

    def mode(self):
        return 0

    def __repr__(self):
        return "gennorm"


class Nakagami:
    def __init__(self, nu):
        self.nu = nu

    def pdf(self, x):
        try:
            z = gamma(self.nu) / 2 / self.nu ** self.nu
            return (x ** (2 * self.nu - 1) * exp(-self.nu * x * x)) / z
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def dpdf(self, x):
        try:
            z = gamma(self.nu) / 2 / self.nu ** self.nu
            first = (2 * self.nu - 1) * x ** (2 * self.nu - 2) * exp(-self.nu * x * x)
            second = x ** (2 * self.nu - 1) * exp(-self.nu * x * x) * (-self.nu * 2 * x)
            return (first + second) / z
        except (ZeroDivisionError, OverflowError):
            return np.inf

    def cdf(self, x):
        return special.gammainc(self.nu, self.nu * x * x)

    def support(self):
        return 0, np.inf

    def mode(self):
        return sqrt(2) / 2 * sqrt((2 * self.nu - 1) / self.nu)

    def __repr__(self):
        return "nakagami"
