import os

import time
import copy

import random
# random.seed(2)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import math as m
import numpy as np
# import sympy as sp
import shapely as sh
from shapely.geometry import Polygon
from shapely import speedups
speedups.enable()
# from geometry.convex_polygon import ConvexPolygon
from pyquaternion import Quaternion

# my stuff
import geometry.solid_angle
import photon_attenuation.materials as materials
import photon_attenuation.physics as physics
from geometry.raytracing import Plane, Ray

# #{ class Timer

class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart)

# #} end of class Timer

# #{ class Polygon3

class Polygon3D:

    # #{ __init__()

    def __init__(self, points):

        if len(points) < 3:

            print("Polygon3D: cannot create a polygon from < 3 vertices!")

            return False

        # basis of the linear subspace (the plane moved to the origin)
        self.basis = []

        # projector to the linear subspace (the plane moved to the origin)
        self.projector = []

        # projection of the origin
        self.zero_point = []

        # all the 2d points on the polygon
        self.points_2d = []

        # sympy polygon
        self.polygon_2d = []

        # find the normal vector of the plane
        v1 = points[1] - points[0]
        v2 = points[2] - points[0]
        self.normal_vector = np.cross(v1, v2)
        self.normal_vector = self.normal_vector / np.linalg.norm(self.normal_vector)

        # the defining point of the plane
        self.zero_point = points[0]

        self.plane = Plane(self.zero_point, self.normal_vector)

        # # find the orthogonal basis of the plane

        # first vector of the orthonormal basis
        if np.abs(self.normal_vector[0]) > 1e-3 or np.abs(self.normal_vector[1]) > 1e-3:
            basis1 = np.array([-self.normal_vector[1], self.normal_vector[0], self.normal_vector[2]])
        else:
            basis1 = np.array([self.normal_vector[0], -self.normal_vector[2], self.normal_vector[1]])

        # second vector of the orthonormal basis
        basis2 = np.cross(self.normal_vector, basis1)
        basis2 = basis2 / np.linalg.norm(basis2)

        self.basis = np.matrix([basis1, basis2, self.normal_vector]).transpose()

        # projector to orthogonal adjucent of the planes normal
        self.projector = self.basis*self.basis.transpose()

        for i,point in enumerate(points):
            self.points_2d.append(self.projectOn2D(point))

        self.polygon_2d = Polygon(self.points_2d)
        # self.polygon_2d = ConvexPolygon(self.points_2d)

    # #}

    # #{ projectOn2D()

    def projectOn2D(self, point_in):

        point = np.matrix([point_in - self.zero_point]).transpose()

        # project the point on the affine subspace of the plane
        projection = self.projector*point

        # express in the orthonormal basis of the subspace
        projection = np.dot(np.linalg.inv(self.basis), projection)

        return np.array([projection[0, 0], projection[1, 0]])

    # #} end of

    # #{ intersection()

    def intersection(self, intersector):

        intersectee_3d = self.plane.intersectionRay(intersector)

        if not isinstance(intersectee_3d, np.ndarray):
            return False

        intersectee_copy = copy.deepcopy(intersectee_3d)

        intersectee_2d = sh.geometry.Point(self.projectOn2D(intersectee_copy))

        if self.polygon_2d.contains(intersectee_2d):
            return intersectee_3d
        else:
            return False

    # #} end of rayCollision()

# #} end of Polygon3D

# #{ class Source

class Source:

    def __init__(self, energy, activity, position):

        self.energy = energy
        self.activity = activity
        self.position = position

# #} end of class Source

# #{ class Detector

class Detector:

    def __init__(self, material, thickness, position):

        self.size = 0.01408 # [mm]
        self.vertices = [] # in the frame of the sensor

        self.front = []
        self.back = []
        self.sides = []

        self.material = material
        self.thickness = thickness
        self.position = position

        # FRONT

        # give me the 4 front_vertices of the detector
        a = np.array([self.thickness/2, -self.size/2.0, self.size/2.0]) + position
        b = np.array([self.thickness/2, self.size/2.0, self.size/2.0]) + position
        c = np.array([self.thickness/2, self.size/2.0, -self.size/2.0]) + position
        d = np.array([self.thickness/2, -self.size/2.0, -self.size/2.0]) + position

        # create the sympy front plane
        self.front = Polygon3D([a, b, c, d])

        # BACK

        # give me the 4 back_vertices of the detector
        e = np.array([-self.thickness/2, -self.size/2.0, self.size/2.0]) + position
        f = np.array([-self.thickness/2, self.size/2.0, self.size/2.0]) + position
        g = np.array([-self.thickness/2, self.size/2.0, -self.size/2.0]) + position
        h = np.array([-self.thickness/2, -self.size/2.0, -self.size/2.0]) + position

        # create the sympy back plane
        self.back = Polygon3D([e, f, g, h])

        # SIDES

        # create the sympy side planes
        self.sides.append(Polygon3D([a, e, h, d]))
        self.sides.append(Polygon3D([b, f, e, a]))
        self.sides.append(Polygon3D([c, g, f, b]))
        self.sides.append(Polygon3D([d, h, c, g]))

        self.sides[0].polygon_2d = self.sides[0].polygon_2d.buffer(0.002)
        self.sides[1].polygon_2d = self.sides[1].polygon_2d.buffer(0.002)
        self.sides[2].polygon_2d = self.sides[2].polygon_2d.buffer(0.002)
        self.sides[3].polygon_2d = self.sides[3].polygon_2d.buffer(0.002)

        self.vertices.append(a)
        self.vertices.append(b)
        self.vertices.append(c)
        self.vertices.append(d)
        self.vertices.append(e)
        self.vertices.append(f)
        self.vertices.append(g)
        self.vertices.append(h)

    def getVertices(self):

        # detector front_vertices in the world frame
        a = self.vertices[0]
        b = self.vertices[1]
        c = self.vertices[2]
        d = self.vertices[3]
        e = self.vertices[4]
        f = self.vertices[5]
        g = self.vertices[6]
        h = self.vertices[7]

        return [a, b, c, d, e, f, g, h]

# #} end of class Detector

# define the source and the detector
source = Source(661000.0, 1e9, np.array([0.1, -0.0, 0.0]))
source_point = source.position
detector_1 = Detector(materials.Si, 0.001, np.array([0, 0, 0]))
detector_2 = Detector(materials.CdTe, 0.002, np.array([-0.003, 0, 0]))

[a1, b1, c1, d1, e1, f1, g1, h1] = detector_1.getVertices()
[a2, b2, c2, d2, e2, f2, g2, h2] = detector_2.getVertices()

# calculate the relative activity to the detector
detector_solid_angle = geometry.solid_angle.quadrilateral_solid_angle(a1, b1, c1, d1, source.position)

# apparent activity
aparent_activity = source.activity*(detector_solid_angle/(4*m.pi))

# #{ comptonScattering()

cs_cross_section = physics.comptonCrossSection(physics.conversions.energy_ev_to_J(source.energy))
cs_density_indeces, cs_density = physics.cs_distribution_function(detector_1.material, source.energy)

def comptonScattering(from_point, to_point, energy, material):

    distance = np.linalg.norm(to_point - from_point)

    # calculate the probability of the scattering

    # draw and deside wether it should happend
    prob_cs = 1.0 - np.exp(-detector_1.material.electron_density * cs_cross_section * distance)

    if random.uniform(0.0, 1.0) > prob_cs:
        return False, 0

    # calculate the point of scattering in the detector
    # scattering_point = (from_point + to_point)/2.0
    position_weight = random.uniform(0.0, 1.0)
    scattering_point = position_weight*from_point + (1 - position_weight)*to_point

    # calculate the azimuthal and radial angle
    phi = random.uniform(0.0, 2.0*m.pi)
    theta = cs_density[int(m.floor(random.uniform(0.0, len(cs_density))))]

    # calculate the point on a unit sphere
    x1 = m.cos(theta)
    y1 = m.sin(theta)*m.sin(phi)
    z1 = m.sin(theta)*m.cos(phi)

    # pitch rotational matrix
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([x1, y1, z1])
    my_axis = np.cross(v1, v2)

    original_direction = to_point - from_point
    original_direction = original_direction/np.linalg.norm(original_direction)

    try:
        my_quaternion = Quaternion(axis=my_axis, angle=geometry.solid_angle.vector_angle(v1, v2))
        new_ray_point = scattering_point + my_quaternion.rotate(original_direction)
    except:
        # no rotation should be applied
        new_ray_point = scattering_point + original_direction
        theta = 0.0

    # calculate the energy of the new photon
    new_photon_energy = source.energy * physics.comptonRatio(physics.conversions.energy_ev_to_J(source.energy), theta)
    electron_energy = source.energy - new_photon_energy

    new_ray = Ray(scattering_point, new_ray_point, new_photon_energy)

    return new_ray, electron_energy

# #} end of comptonScattering()

# #{ sample_detector()

def sample_detector(detector):

    [a, b, c, d, e, f, g, h] = detector.getVertices()

    ab = b-a
    ad = d-a

    k1 = random.uniform(1e-2, 1.0-1e-2)
    k2 = random.uniform(1e-2, 1.0-1e-2)

    return a + ab*k1 + ad*k2

# #} end of sample_detector()

# sample the 1st detector
hit_points = []
rays = []
n_particles = 50000
for i in range(0, n_particles):

   hit_point = sample_detector(detector_1)
   hit_points.append(hit_point)

   ray = Ray(source_point, hit_point, source.energy)
   rays.append(ray)

# plotting
def plot_everything(*args):

    fig = plt.figure(1)
    ax  = fig.add_subplot(111, projection = '3d')

    # #{ plot detector 1

    xs = [a1[0], b1[0], c1[0], d1[0], a1[0], e1[0], f1[0], b1[0], f1[0], g1[0], c1[0], g1[0], h1[0], d1[0], h1[0], e1[0]]
    ys = [a1[1], b1[1], c1[1], d1[1], a1[1], e1[1], f1[1], b1[1], f1[1], g1[1], c1[1], g1[1], h1[1], d1[1], h1[1], e1[1]]
    zs = [a1[2], b1[2], c1[2], d1[2], a1[2], e1[2], f1[2], b1[2], f1[2], g1[2], c1[2], g1[2], h1[2], d1[2], h1[2], e1[2]]
    plt.plot(xs, ys, zs)

    # #} end of plot detector 1

    # #{ plot detector 2

    xs = [a2[0], b2[0], c2[0], d2[0], a2[0], e2[0], f2[0], b2[0], f2[0], g2[0], c2[0], g2[0], h2[0], d2[0], h2[0], e2[0]]
    ys = [a2[1], b2[1], c2[1], d2[1], a2[1], e2[1], f2[1], b2[1], f2[1], g2[1], c2[1], g2[1], h2[1], d2[1], h2[1], e2[1]]
    zs = [a2[2], b2[2], c2[2], d2[2], a2[2], e2[2], f2[2], b2[2], f2[2], g2[2], c2[2], g2[2], h2[2], d2[2], h2[2], e2[2]]
    plt.plot(xs, ys, zs)

    # #} end of plot detector 2

    with Timer("test"):

        for index,point in enumerate(hit_points):

            # intersection with the back side of the 1st detector
            intersect1_second = detector_1.back.intersection(rays[index])

            # no collision with the back face
            if not isinstance(intersect1_second, np.ndarray):

                # check intersection with the sides
                for i,side in enumerate(detector_1.sides):

                    intersect1_second = side.intersection(rays[index])

                    if isinstance(intersect1_second, np.ndarray):
                        break

            # if the ray came from the oposite direction, discard it
            if np.linalg.norm(intersect1_second - source.position) < np.linalg.norm(point - source.position):
                continue

            # if there is a collision
            if not isinstance(intersect1_second, np.ndarray):
                print("! no intersection with the back/side face of the first detector")
                continue

            intersection_len = np.linalg.norm(point - intersect1_second)

            scattered_ray, electron_energy = comptonScattering(point, intersect1_second, source.energy, detector_1.material)
            # if not scattering happened, just leave
            if not isinstance(scattered_ray, Ray):
                continue

            # # plot every scattered ray
            # ax.plot(
            #     [scattered_ray.rayPoint[0], scattered_ray.rayPoint[0]+0.01*scattered_ray.rayDirection[0]],
            #     [scattered_ray.rayPoint[1], scattered_ray.rayPoint[1]+0.01*scattered_ray.rayDirection[1]],
            #     [scattered_ray.rayPoint[2], scattered_ray.rayPoint[2]+0.01*scattered_ray.rayDirection[2]],
            #     color='b'
            # )

            # check the collision with the other detector's front side
            intersect2_first = detector_2.front.intersection(scattered_ray)

            # if there is no intersection with the front face
            if not isinstance(intersect2_first, np.ndarray):
                continue

            # check the collision with the other detector's back side
            intersect2_second = detector_2.back.intersection(scattered_ray)

            # no collision with the back face
            if not isinstance(intersect2_second, np.ndarray):

                ## intersection with the sides
                for i,side in enumerate(detector_2.sides):

                    intersect2_second = side.intersection(scattered_ray)

                    if isinstance(intersect2_second, np.ndarray):
                        break

            if not isinstance(intersect2_second, np.ndarray):
                print("!! no intersection with the back/side face of the second detector")

                ax.scatter(intersect2_first[0], intersect2_first[1], intersect2_first[2], color='blue')
                ray_len = 0.005
                ax.plot([scattered_ray.rayPoint[0], scattered_ray.rayPoint[0]+ray_len*scattered_ray.rayDirection[0]], [scattered_ray.rayPoint[1], scattered_ray.rayPoint[1]+ray_len*scattered_ray.rayDirection[1]], [scattered_ray.rayPoint[2], scattered_ray.rayPoint[2]+ray_len*scattered_ray.rayDirection[2]], color='red')

                continue

            pe_cross_section = [physics.pe_cs_gavrila_pratt_simplified(mat, scattered_ray.energy) for mat in detector_2.material.elements]
            pe_thickness = np.linalg.norm(intersect2_second - intersect2_first)

            prob_pe = 1.0
            for index,cross_section in enumerate(pe_cross_section):
                prob_pe *= np.exp(-detector_2.material.element_quantities[index] * detector_2.material.molecular_density * cross_section * pe_thickness)
            prob_pe = 1.0 - prob_pe

            if random.uniform(0.0, 1.0) > prob_pe:
                continue

            # sample the interraction point in the 2nd detector
            position_weight = random.uniform(0.0, 1.0)
            absorption_point = position_weight*intersect2_first + (1 - position_weight)*intersect2_second

            print("energy: {}, thickness: {}, prob_pe: {}".format(scattered_ray.energy, pe_thickness, prob_pe))

            # ax.scatter(intersect2_first[0], intersect2_first[1], intersect2_first[2], color='red')
            ax.plot([source.position[0], point[0]], [source.position[1], point[1]], [source.position[2], point[2]], color='grey')
            ax.plot([point[0], intersect1_second[0]], [point[1], intersect1_second[1]], [point[2], intersect1_second[2]], color='r')
            ax.plot([scattered_ray.rayPoint[0], absorption_point[0]], [scattered_ray.rayPoint[1], absorption_point[1]], [scattered_ray.rayPoint[2], absorption_point[2]], color='b')

    # #{ set axis

    axis_range = 0.01

    ax.set_xlim(-axis_range + detector_1.position[0], axis_range + detector_1.position[0])
    ax.set_ylim(-axis_range + detector_1.position[1], axis_range + detector_1.position[1])
    ax.set_zlim(-axis_range + detector_1.position[2], axis_range + detector_1.position[2])

    plt.gca().set_aspect('equal', adjustable='box')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # #} end of set axis

    ax.grid(True)

    plt.show()

pid = os.fork()
if pid == 0:
    plot_everything()
