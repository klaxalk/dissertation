import numpy as np
import math as m
import photon_attenuation.materials as materials
import geometry.solid_angle
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

from sympy import Point3D, Point2D, Polygon, Ray3D, Plane, linsolve, Matrix, Line3D
import time

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

# #{ Polygon3D

class Polygon3D:

    plane = []

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

        self.plane = Plane(points[0], points[1], points[2])

        normal = self.plane.normal_vector

        # convert the normal vector to numpy array and normalize it
        normal_vector = np.array([float(normal[0]), float(normal[1]), float(normal[2])])
        normal_vector = normal_vector / np.linalg.norm(normal_vector)

        # find the point closest to the origin
        zero_line = self.plane.perpendicular_line(Point3D(0, 0, 0, evaluate=False))

        self.zero_point = self.plane.intersection(zero_line)[0]

        # # find the orthogonal basis of the plane

        # first vector of the orthonormal basis
        if np.abs(normal_vector[0]) > 1e-3 or np.abs(normal_vector[1]) > 1e-3:
            basis1 = np.array([-normal_vector[1], normal_vector[0], normal_vector[2]])
        else:
            basis1 = np.array([normal_vector[0], -normal_vector[2], normal_vector[1]])

        # second vector of the orthonormal basis
        basis2 = np.cross(normal_vector, basis1)
        basis2 = basis2 / np.linalg.norm(basis2)

        self.basis = np.matrix([basis1, basis2, normal_vector]).transpose()

        # projector to orthogonal adjucent of the planes normal
        self.projector = self.basis*self.basis.transpose()
        # projector_normal = np.eye(3, dtype=float) - self.basis*self.basis.transpose()

        for i,point in enumerate(points):
            self.points_2d.append(self.projectOn2D(point))

        print("len(self.points_2d): {}".format(len(self.points_2d)))

        self.polygon_2d = Polygon(*self.points_2d)


    # #} end of 

    # #{ projectOn2D()

    def projectOn2D(self, point_in):

        point = np.matrix([point_in - self.zero_point]).transpose()

        # project the point on the affine subspace of the plane
        projection = self.projector*point

        # express in the orthonormal basis of the subspace 
        projection = np.dot(np.linalg.inv(self.basis), projection)

        return Point2D(projection[0, 0], projection[1, 0], evaluate=False)

    # #} end of 

    # #{ intersection()
    
    def intersection(self, intersector): 

        with Timer("intersection"):

            with Timer("subintersection"):

                print("intersector: {}".format(intersector))

                print("self.plane: {}".format(self.plane))

                # create the sympy intersection with the plane
                intersectee_3d = self.plane.intersection(intersector)

                print("intersectee_3d: {}".format(intersectee_3d))

            if not isinstance(intersectee_3d, Point3D):

                print("Polygon3D: the result of an intersection is not of Point3D type")

                return False

            with Timer("project2d"):
                intersectee_2d = self.projectOn2D(intersectee_3d)

            if self.polygon_2d.encloses_point(intersectee_2d):
                return intersectee_2d
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
        self.sp_vertices = []

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

        # create the sympy representation of the front vertices
        sp_a = Point3D(a, evaluate=False)
        sp_b = Point3D(b, evaluate=False)
        sp_c = Point3D(c, evaluate=False)
        sp_d = Point3D(d, evaluate=False)

        # create the sympy front plane
        self.front = Polygon3D([sp_a, sp_b, sp_c, sp_d])

        # BACK

        # give me the 4 back_vertices of the detector
        e = np.array([-self.thickness/2, -self.size/2.0, self.size/2.0])
        f = np.array([-self.thickness/2, self.size/2.0, self.size/2.0])
        g = np.array([-self.thickness/2, self.size/2.0, -self.size/2.0])
        h = np.array([-self.thickness/2, -self.size/2.0, -self.size/2.0])

        # create the sympy representation of the front vertices
        sp_e = Point3D(e, evaluate=False)
        sp_f = Point3D(f, evaluate=False)
        sp_g = Point3D(g, evaluate=False)
        sp_h = Point3D(h, evaluate=False)

        # create the sympy back plane
        self.back = Polygon3D([sp_e, sp_f, sp_g, sp_h])

        # SIDES

        # create the sympy side planes
        self.sides.append(Polygon3D([sp_a, sp_e, sp_h, sp_d]))
        self.sides.append(Polygon3D([sp_b, sp_f, sp_e, sp_a]))
        self.sides.append(Polygon3D([sp_c, sp_g, sp_f, sp_b]))
        self.sides.append(Polygon3D([sp_d, sp_h, sp_c, sp_g]))

        # #{ push to sp_vertices
        
        self.sp_vertices.append(a)
        self.sp_vertices.append(b)
        self.sp_vertices.append(c)
        self.sp_vertices.append(d)
        self.sp_vertices.append(e)
        self.sp_vertices.append(f)
        self.sp_vertices.append(g)
        self.sp_vertices.append(h)

        # #} end of push to sp_vertices

        # #{ push to vertices
        
        self.vertices.append(a)
        self.vertices.append(b)
        self.vertices.append(c)
        self.vertices.append(d)
        self.vertices.append(e)
        self.vertices.append(f)
        self.vertices.append(g)
        self.vertices.append(h)
        
        # #} end of push to vertices

    def getWorldVertices(self):

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

# define the source and the detector
source = Source(662000.0, 1e9, np.array([10.0, 0.0, 0.0]))
source_point = Point3D(source.position, evaluate=False)
detector_1 = Detector(materials.Si, 0.001, np.array([0, 0, 0]))
detector_2 = Detector(materials.CdTe, 0.001, np.array([-0.002, 0, 0]))

[a1, b1, c1, d1, e1, f1, g1, h1] = detector_1.getWorldVertices()
[a2, b2, c2, d2, e2, f2, g2, h2] = detector_2.getWorldVertices()

# calculate the relative activity to the detector
detector_solid_angle = geometry.solid_angle.quadrilateral_solid_angle(a1, b1, c1, d1, source.position)

# apparent activity
aparent_activity = source.activity*(detector_solid_angle/(4*m.pi))

# #{ sample_detector()

def sample_detector(detector):

    # get the detector vertices
    # [a, b, c, d] = detector.getWorldVertices()

    # # generate the convex combination of the vertices
    # k = []
    # k.append(random.uniform(0.0, 1.0))
    # k.append(random.uniform(0.0, 1.0))
    # k.append(random.uniform(0.0, 1.0))

    # k = np.sort(k)

    # sampled_point = a*k[0] + b*(k[1] - k[0]) + c*(k[2] - k[1]) + d*(1.0 - k[2])

    [a, b, c, d, e, f, g, h] = detector.getWorldVertices()

    ab = b-a
    ad = d-a

    k1 = random.uniform(0.0, 1.0)
    k2 = random.uniform(0.0, 1.0)

    return a + ab*k1 + ad*k2

# #} end of sample_detector()

# sample the 1st detector
hit_points = []
rays = []
for i in range(0, 10):

   hit_point = sample_detector(detector_1)
   hit_points.append(hit_point) 

   ray = Line3D(source_point, Point3D(hit_point, evaluate=False))
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

    # plot the sampled points
    for index,point in enumerate(hit_points):

        # ax.scatter(point[0], point[1], point[2], color='g', marker='o')

        intersect = detector_1.front.intersection(rays[index])

        if isinstance(intersect, Point3D):

            ax.scatter(intersect[0], intersect[1], intersect[2], color='g', marker='o')

    # plot the source
    # plt.scatter(source.position[0], source.position[1], source.position[2], color='r', marker='o')

    ax.set_xlim(-0.01 + detector_1.position[0], 0.01 + detector_1.position[0])
    ax.set_ylim(-0.01 + detector_1.position[1], 0.01 + detector_1.position[1])
    ax.set_zlim(-0.01 + detector_1.position[2], 0.01 + detector_1.position[2])

    ax.grid(True)
    plt.show()

pid = os.fork()
if pid == 0:
    plot_everything()
