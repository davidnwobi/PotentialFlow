import numpy as np
from math import ceil
import matplotlib.pyplot as plt
from scipy import stats


def mean_camber_1(m: float, p: float, x: float) -> float:

    return  m / p ** 2 * (2 * p * x - x ** 2) if m != 0 else 0


def mean_camber_2(m: float, p: float, x: float) -> float:
    return m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2) if m != 0 else 0


def mean_camber(m: float, p: float, x: float) -> float:
    return np.where(x < p, mean_camber_1(m, p, x), mean_camber_2(m, p, x)) if m != 0 else 0


def mean_camber_derivative_1(m: float, p: float, x: float) -> float:
    return 2 * m / p ** 2 * (p - x) if m != 0 else 0


def mean_camber_derivative_2(m: float, p: float, x: float) -> float:
    return 2 * m / (1 - p) ** 2 * (p - x) if m != 0 else 0


def mean_camber_derivative(m: float, p: float, x: float) -> float:
    return np.where(x < p, mean_camber_derivative_1(m, p, x), mean_camber_derivative_2(m, p, x)) if m != 0 else 0


def thickness(x: float, t: float) -> float:
    return 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1015 * x ** 4)


def generate_four_digit_NACA(num_NACA: str, num_points: int, chord_length: float) -> np.ndarray:
    if len(num_NACA) != 4:
        raise ValueError("NACA number must be 4 digits long")

    num_points = int(ceil(num_points/2))
    scale = chord_length / 1
    num1 = int(num_NACA[0])
    num2 = int(num_NACA[1])
    num34 = int(num_NACA[2:])

    x = 1 / 2 * (1 - np.cos(np.linspace(0, np.pi, num_points))) # circle transformation
    y = mean_camber(num1 / 100, num2 / 10, x)
    y_t = thickness(x, num34 / 100 * chord_length)
    theta = np.arctan2(mean_camber_derivative(num1 / 100, num2 / 10, x), 1)

    x = x * scale
    y = y * scale
    y_t = y_t * scale

    x_upper = x - y_t * np.sin(theta)
    y_upper = y + y_t * np.cos(theta)
    x_upper = np.flip(x_upper)
    y_upper = np.flip(y_upper)

    x_lower = x + y_t * np.sin(theta)
    y_lower = y - y_t * np.cos(theta)

    x, y = np.concatenate((x_upper, x_lower[1:])), np.concatenate((y_upper, y_lower[1:]))

    return x, y
