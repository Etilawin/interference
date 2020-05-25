from InterferenceGenerator import *


if __name__ == "__main__":
    options = {
        "wavelengths": range(380, 780, 200)
    }

    system = InterferenceGenerator(wavelengths=range(380, 780, 200))
    system.run()


