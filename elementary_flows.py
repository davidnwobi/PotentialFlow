import abc
import numpy as np
import typing as tp
from dataclasses import dataclass
from abc import abstractmethod


class ElementaryFlow:
    pass


@dataclass
class UniformFlow(ElementaryFlow):
    horizontal_vel: float
    vertical_vel: float

    def stream_function(self, x: np.ndarray, y: np.ndarray):
        return self.horizontal_vel * y - self.vertical_vel * x

    def velocity(self, x: np.ndarray, y: np.ndarray):
        horizontal_vel = np.full_like(x, self.horizontal_vel)
        vertical_vel = np.full_like(y, self.vertical_vel)
        return horizontal_vel, vertical_vel


@dataclass
class NonUniformFlow(ElementaryFlow):
    x_pos: float
    y_pos: float
    mask_tol: float = 1e-6

    @abstractmethod
    def stream_function(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        raise NotImplementedError("Not Implemented")

    @abstractmethod
    def velocity(self, x: np.ndarray, y: np.ndarray) -> tp.Tuple[np.ndarray, np.ndarray]:
        raise NotImplementedError("Not Implemented")

    def r_squared(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return (x - self.x_pos) ** 2 + (y - self.y_pos) ** 2

    def theta(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.arctan((y - self.y_pos) / (x - self.x_pos))

    def theta_mask(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.abs(self.theta(x, y)) < self.mask_tol

    def r_squared_mask(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return self.r_squared(x, y) < self.mask_tol


@dataclass
class Source(NonUniformFlow):
    strength: float = 0.0

    def stream_function(self, x: np.ndarray, y: np.ndarray):
        mask = self.theta_mask(x, y)
        x, y = (np.ma.masked_where(mask, x), np.ma.masked_where(mask, y))
        return self.strength / (2 * np.pi) * self.theta(x, y)

    def velocity(self, x: np.ndarray, y: np.ndarray):
        mask = self.r_squared_mask(x, y)
        x, y = (np.ma.masked_where(mask, x), np.ma.masked_where(mask, y))

        r_squared = self.r_squared(x, y)
        u = self.strength / (2 * np.pi) * (1 / r_squared) * (x - self.x_pos)
        v = self.strength / (2 * np.pi) * (1 / r_squared) * (y - self.y_pos)

        return u, v


@dataclass
class Vortex(NonUniformFlow):
    circulation: float = 0.0

    def stream_function(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        r_squared_mask = (x - self.x_pos) ** 2 + (y - self.y_pos) ** 2

        def mask():
            return r_squared_mask < self.mask_tol

        x, y = (np.ma.masked_where(mask(), x), np.ma.masked_where(mask(), y))

        r = np.log(np.sqrt((x - self.x_pos) ** 2 + (y - self.y_pos) ** 2))
        return self.circulation / (2 * np.pi) * r

    def velocity(self, x: np.ndarray, y: np.ndarray) -> tp.Tuple[np.ndarray, np.ndarray]:
        r_squared_mask = (x - self.x_pos) ** 2 + (y - self.y_pos) ** 2

        def mask():
            return r_squared_mask < self.mask_tol

        x, y = (np.ma.masked_where(mask(), x), np.ma.masked_where(mask(), y))

        r = np.sqrt((x - self.x_pos) ** 2 + (y - self.y_pos) ** 2)
        u = self.circulation / (2 * np.pi) * (1 / r) * ((y - self.y_pos) / r)
        v = -self.circulation / (2 * np.pi) * (1 / r) * ((x - self.x_pos) / r)

        return u, v


@dataclass
class Doublet(NonUniformFlow):
    kappa: float = 0.0

    def stream_function(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        mask = self.r_squared_mask(x, y)
        x, y = (np.ma.masked_where(mask, x), np.ma.masked_where(mask, y))

        r_squared = self.r_squared(x, y)
        return -self.kappa / (2 * np.pi) * ((y-self.y_pos) / r_squared)

    def velocity(self, x: np.ndarray, y: np.ndarray) -> tp.Tuple[np.ndarray, np.ndarray]:

        mask = self.r_squared_mask(x, y)
        x, y = (np.ma.masked_where(mask, x), np.ma.masked_where(mask, y))

        r_squared_squared = self.r_squared(x, y) ** 2
        U = -self.kappa / (2 * np.pi) * (x-self.x_pos+y-self.y_pos) * (x-self.x_pos-y+self.y_pos)/ (r_squared_squared)
        V = -self.kappa / (2 * np.pi) * 2 * (x-self.x_pos)*(y-self.y_pos)/ (r_squared_squared)
        return U, V
