from math import pi, sin, cos, asin, atan, sqrt
from numpy import array, zeros
from matplotlib.pyplot import imshow, show


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
    return array([R, G, B], dtype=float)


class InterferenceGenerator:
    def __init__(self, **kwargs):
        self.n1 = kwargs.get("n1", 1)
        self.n2 = kwargs.get("n2", 1.5)
        self.thickness = kwargs.get("thickness", 1e-3)  # m
        self.f = kwargs.get("f", 5)  # focal lens

        self.wavelengths = kwargs.get("wavelengths", [500])
        self.wavelength_rgb_dict = self.__generate_wavelength_dict()

        self.r, self.t, self.rr, self.tt = self.__get_coefficients()
        self.alpha = abs(self.t * self.tt)
        self.A = self.r ** 2 * (1 + self.alpha ** 2)
        self.B = 2 * self.alpha / (1 + self.alpha ** 2)
        self.side = kwargs.get("side", 1.0)
        self.points = kwargs.get("points", 1000)
        self.normalization = self.A * (1 + self.B)

    def __generate_wavelength_dict(self):
        return {wavelength: wavelength_to_rgb(wavelength) for wavelength in self.wavelengths}

    def __get_coefficients(self):
        r = (self.n1 - self.n2) / (self.n1 + self.n2)
        t = 2 * self.n1 / (self.n1 + self.n2)
        rr = -r
        tt = 2 * self.n2 / (self.n1 + self.n2)

        return r, t, rr, tt

    def __phi0(self, radius):
        return 2 * pi * 2 * self.n2 * self.thickness * cos(asin(self.n1 * sin(atan(radius / self.f)) / self.n2))

    def __cumulative_intensity(self, radius):
        res = zeros(3)
        phase0 = self.__phi0(radius)
        for wavelength in self.wavelength_rgb_dict:
            I = self.A * (1 - self.B * cos(phase0 * 1e9 / wavelength)) / self.normalization
            res += I * self.wavelength_rgb_dict[wavelength]
        return res

    def run(self):

        spacing = self.side / self.points
        xi = zeros([self.points, self.points, 3], float)

        temp = {}

        for i in range(self.points):
            y = spacing * i
            for j in range(self.points):
                x = spacing * j
                radius = sqrt((self.side / 2 - x) ** 2 + (self.side / 2 - y) ** 2)
                radius = round(radius, 3)
                temp[radius] = temp.get(radius, self.__cumulative_intensity(radius))
                xi[i, j, :] = temp[radius]

        imshow(xi, extent=[0, self.side, 0, self.side])
        show()