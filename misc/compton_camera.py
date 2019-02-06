import os
import random
import time

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import math as m
import numpy as np
# import sympy as sp
import shapely as sh
from shapely.geometry import Polygon
# from geometry.convex_polygon import ConvexPolygon

# my stuff
import geometry.solid_angle
import photon_attenuation.materials as materials
from  geometry.raytracing import Plane, Ray

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

# #{ class Polygon3D

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

    # #} end of 

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

        try:
            intersectee_3d = self.plane.intersectionRay(intersector)
        except:
            return False

        # intersectee_2d = self.projectOn2D(intersectee_3d)
        intersectee_2d = sh.geometry.Point(self.projectOn2D(intersectee_3d))


        if self.polygon_2d.contains(intersectee_2d):
        # if self.polygon_2d.isPointIn(intersectee_2d):
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

        a = np.array([self.thickness/2, -self.size/2.0, self.size/2.0])
        b = np.array([self.thickness/2, self.size/2.0, self.size/2.0])
        c = np.array([self.thickness/2, self.size/2.0, -self.size/2.0])
        d = np.array([self.thickness/2, -self.size/2.0, -self.size/2.0])

        # create the sympy front plane
        self.front = Polygon3D([a, b, c, d])

        # BACK

        # give me the 4 back_vertices of the detector
        e = np.array([-self.thickness/2, -self.size/2.0, self.size/2.0])
        f = np.array([-self.thickness/2, self.size/2.0, self.size/2.0])
        g = np.array([-self.thickness/2, self.size/2.0, -self.size/2.0])
        h = np.array([-self.thickness/2, -self.size/2.0, -self.size/2.0])

        # create the sympy back plane
        self.back = Polygon3D([e, f, g, h])

        # SIDES

        # create the sympy side planes
        self.sides.append(Polygon3D([a, e, h, d]))
        self.sides.append(Polygon3D([b, f, e, a]))
        self.sides.append(Polygon3D([c, g, f, b]))
        self.sides.append(Polygon3D([d, h, c, g]))
        
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
        a = self.position + self.vertices[0]
        b = self.position + self.vertices[1]
        c = self.position + self.vertices[2]
        d = self.position + self.vertices[3]
        e = self.position + self.vertices[4]
        f = self.position + self.vertices[5]
        g = self.position + self.vertices[6]
        h = self.position + self.vertices[7]

        return [a, b, c, d, e, f, g, h]

# #} end of class Detector

# #{ comptonScattering()

def comptonScattering(from_point, to_point, energy, material):

    distance = np.linalg.norm(to_point - from_point)

    # calculate the probability of the scattering
    prob = 1.0 # TODO

    # draw and deside wether it should happend
    # it happends # TODO

    # calculate the point of scattering in the detector
    scattering_point = (from_point + to_point)/2.0

    # calculate the azimuthal and radial angle
    # phi = 1.57 # TODO, uniform distribution
    # theta = m.pi/6.0 # TODO, draw from Klein-Nishina
    phi = 0.0 # TODO, uniform distribution
    theta = 0.0 # TODO, draw from Klein-Nishina

    # calculate the point on a unit sphere
    r = 0.01
    x1 = r*m.cos(theta)
    y1 = r*m.cos(theta)*m.sin(phi)
    z1 = r*m.sin(theta)*m.sin(phi)

    through_point = scattering_point + np.array([x1, y1, z1])

    new_ray = Ray(scattering_point, through_point, source.energy)

    return new_ray

# #} end of comptonScattering()

# define the source and the detector
source = Source(662000.0, 1e9, np.array([1.0, -1.0, 0.0]))
source_point = source.position
detector_1 = Detector(materials.Si, 0.001, np.array([0, 0, 0]))
detector_2 = Detector(materials.CdTe, 0.001, np.array([-0.005, 0, 0]))

[a1, b1, c1, d1, e1, f1, g1, h1] = detector_1.getVertices()
[a2, b2, c2, d2, e2, f2, g2, h2] = detector_2.getVertices()

# calculate the relative activity to the detector
detector_solid_angle = geometry.solid_angle.quadrilateral_solid_angle(a1, b1, c1, d1, source.position)

# apparent activity
aparent_activity = source.activity*(detector_solid_angle/(4*m.pi))

# #{ sample_detector()

def sample_detector(detector):

    [a, b, c, d, e, f, g, h] = detector.getVertices()

    ab = b-a
    ad = d-a

    k1 = random.uniform(0.0, 1.0)
    k2 = random.uniform(0.0, 1.0)

    return a + ab*k1 + ad*k2

# #} end of sample_detector()

# sample the 1st detector
hit_points = []
rays = []
for i in range(0, 500):

   hit_point = sample_detector(detector_1)
   hit_points.append(hit_point) 

   ray = Ray(source_point, hit_point, source.energy)
   rays.append(ray)

# plotting
def plot_everything(*args):

    fig = plt.figure(1)
    ax  = fig.add_subplot(111, projection = '3d')

    # plot the detector
    xs = [a1[0], b1[0], c1[0], d1[0], a1[0], e1[0], f1[0], b1[0], f1[0], g1[0], c1[0], g1[0], h1[0], d1[0], h1[0], e1[0]]
    ys = [a1[1], b1[1], c1[1], d1[1], a1[1], e1[1], f1[1], b1[1], f1[1], g1[1], c1[1], g1[1], h1[1], d1[1], h1[1], e1[1]]
    zs = [a1[2], b1[2], c1[2], d1[2], a1[2], e1[2], f1[2], b1[2], f1[2], g1[2], c1[2], g1[2], h1[2], d1[2], h1[2], e1[2]]
    plt.plot(xs, ys, zs)

    xs = [a2[0], b2[0], c2[0], d2[0], a2[0], e2[0], f2[0], b2[0], f2[0], g2[0], c2[0], g2[0], h2[0], d2[0], h2[0], e2[0]]
    ys = [a2[1], b2[1], c2[1], d2[1], a2[1], e2[1], f2[1], b2[1], f2[1], g2[1], c2[1], g2[1], h2[1], d2[1], h2[1], e2[1]]
    zs = [a2[2], b2[2], c2[2], d2[2], a2[2], e2[2], f2[2], b2[2], f2[2], g2[2], c2[2], g2[2], h2[2], d2[2], h2[2], e2[2]]
    plt.plot(xs, ys, zs)

    with Timer("test"):

        for index,point in enumerate(hit_points):

            ## intersection with the back side

            # plot the sampled points
            intersect = detector_1.back.intersection(rays[index])
            if np.linalg.norm(intersect - source.position) < np.linalg.norm(point - source.position):
                intersect = []
            my_color = 'g'

            if not isinstance(intersect, np.ndarray):

                ## intersection with the sides
                for i,side in enumerate(detector_1.sides):

                    intersect = side.intersection(rays[index])

                    if np.linalg.norm(intersect - source.position) < np.linalg.norm(point - source.position):
                        intersect = []

                    if isinstance(intersect, np.ndarray):

                        my_color = 'b'
                        break

            if isinstance(intersect, np.ndarray):

                ax.plot([point[0], intersect[0]], [point[1], intersect[1]], [point[2], intersect[2]], color=my_color)
                intersection_len = np.linalg.norm(point - intersect)

                scattered_ray = comptonScattering(point, intersect, source.energy, detector_1.material)

                test_point = scattered_ray.rayPoint + scattered_ray.rayDirection

                ax.scatter(test_point[0], test_point[1], test_point[2], color='b', marker='o')


    ax.set_xlim(-0.01 + detector_1.position[0], 0.01 + detector_1.position[0])
    ax.set_ylim(-0.01 + detector_1.position[1], 0.01 + detector_1.position[1])
    ax.set_zlim(-0.01 + detector_1.position[2], 0.01 + detector_1.position[2])

    ax.grid(True)
    plt.show()

pid = os.fork()
if pid == 0:
    plot_everything()
