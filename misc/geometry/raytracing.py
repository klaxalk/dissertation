import numpy as np
import math as m
import sympy as sp

class Plane:

    def __init__(self, planePoint, planeNormal):

        self.planePoint = planePoint
        self.planeNormal = planeNormal

    def intersectionRay(self, ray, epsilon=1e-6):

        ndotu = self.planeNormal.dot(ray.rayDirection)
        if abs(ndotu) < epsilon:
            raise RuntimeError("no intersection or line is within plane")

        w = ray.rayPoint - self.planePoint
        si = -self.planeNormal.dot(w) / ndotu
        Psi = w + si * ray.rayDirection + self.planePoint

        return Psi

class Ray:

    def __init__(self, rayPoint, ray2point, energy):

        self.rayPoint = rayPoint
        self.rayDirection = ray2point - rayPoint
        self.energy = energy
