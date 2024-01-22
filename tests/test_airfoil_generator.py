from src.airfoil_generator import generate_four_digit_NACA
from scipy import stats
import numpy as np


def test_generate_four_digit_NACA():
    data = np.loadtxt('airfoils/NACA2412.txt')
    XB = data[:, 0]
    YB = data[:, 1]
    x, y = generate_four_digit_NACA(2412, len(XB), 1)

    assert stats.spearmanr(YB[:-1], y).statistic > 0.99


def test_generate_four_digit_NACA_2():
    data = np.loadtxt('airfoils/NACA4412.txt')
    XB = data[:, 0]
    YB = data[:, 1]
    x, y = generate_four_digit_NACA(2412, len(XB), 1)

    assert stats.spearmanr(XB[:-1], x).statistic > 0.99
