from math import pi, sin, cos, asin, atan, sqrt

import numpy as np
from matplotlib.pyplot import imshow, show

n1 = 1
n2 = 1.5
e = 1e-3  # m
lam = 400e-9  # m
# theta_i = Angle(1)  # degrees
f = 5  # focal lens


def generate_wavelength_dict(wavelengths):
    return {wavelength: wavelength_to_rgb(wavelength) for wavelength in wavelengths}


wavelength_rgb_dict = generate_wavelength_dict(range(380, 780, 50))


def get_coefficients(n1, n2):
    r = (n1 - n2) / (n1 + n2)
    t = 2 * n1 / (n1 + n2)
    rr = -r
    tt = 2 * n2 / (n1 + n2)

    return r, t, rr, tt


def phi(radius, n1, n2, f, lam, thickness):
    return 2 * pi * 2 * n2 * thickness * cos(asin(n1 * sin(atan(radius / f)) / n2)) / lam


def phi0(radius, n1, n2, f, thickness):
    return 2 * pi * 2 * n2 * thickness * cos(asin(n1 * sin(atan(radius / f)) / n2))


def cumulative_intensity(A, B, radius, n1, n2, f, thickness):
    res = np.zeros(3)
    phase0 = phi0(radius, n1, n2, f, thickness)
    for wavelength in wavelength_rgb_dict:
        I = A * (1 - B * cos(phase0 * 1e9 / wavelength)) / (A * (1 + B))
        res += I * wavelength_rgb_dict[wavelength]
    return res


def intensity(A, B, radius, n1, n2, f, lam, thickness):
    return A * (1 - B * cos(phi(radius, n1, n2, f, lam, thickness))) / (A * (1 + B))


def wavelength_to_rgb(wavelength):
    wavelength = float(wavelength)
    if 380 <= wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation)
        G = 0.0
        B = (1.0 * attenuation)
    elif 440 <= wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440))
        B = 1.0
    elif 490 <= wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490))
    elif 510 <= wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510))
        G = 1.0
        B = 0.0
    elif 580 <= wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580))
        B = 0.0
    elif 645 <= wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation)
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return np.array([R, G, B], dtype=float)


if __name__ == "__main__":
    r, t, rr, tt = get_coefficients(n1, n2)
    alpha = abs(t * tt)
    A = r ** 2 * (1 + alpha ** 2)
    B = 2 * alpha / (1 + alpha ** 2)

    side = 1.0
    points = 1000

    spacing = side / points
    xi = np.zeros([points, points, 3], float)

    temp = {}

    for i in range(points):
        y = spacing * i
        for j in range(points):
            x = spacing * j
            r = sqrt((side / 2 - x) ** 2 + (side / 2 - y) ** 2)
            r = round(r, 3)
            temp[r] = temp.get(r, cumulative_intensity(A, B, r, n1, n2, f, e))
            xi[i, j, :] = temp[r]

    imshow(xi, extent=[0, side, 0, side])
    show()
