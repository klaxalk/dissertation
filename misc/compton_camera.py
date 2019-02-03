import numpy as np
import math as m
import photon_attenuation.materials as materials
import geometry.solid_angle
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

from sympy import Point3D, Point2D, Polygon, Ray3D, Plane, linsolve, Matrix
from sympy.physics.vector import ReferenceFrame, Vector
from sympy.abc import x, y, z

class MyPlane:

    # the sympy object
    plane = []

    # basis of the linear subspace (the plane moved to the origin)
    basis = []

    # projector to the linear subspace (the plane moved to the origin)
    projector = []

    # projection of the origin
    zero_point = []

    def __init__(self, point1, point2, point3):

        self.plane = Plane(point1, point2, point3)

        normal = self.plane.normal_vector

        # convert the normal vector to numpy array and normalize it
        normal_vector = np.array([float(normal[0]), float(normal[1]), float(normal[2])])
        normal_vector = normal_vector / np.linalg.norm(normal_vector)

        # find the point closest to the origin
        zero_line = self.plane.perpendicular_line(Point3D(0, 0, 0))
        self.zero_point = self.plane.intersection(zero_line)[0]

        # # find the orthogonal basis of the plane

        # first vector of the orthonormal basis
        if np.abs(normal_vector[0]) > 1e-3 or np.abs(normal_vector[1]) > 1e-3:
            basis1 = np.array([-normal_vector[1], normal_vector[0], normal_vector[2]])
        else:
            basis1 = np.array([normal_vector[0], -normal_vector[2], normal_vector[1]])

        print("basis1: {}".format(basis1))
        print("normal_vector: {}".format(normal_vector))
        # second vector of the orthonormal basis
        basis2 = np.cross(normal_vector, basis1)
        print("basis2: {}".format(basis2))
        basis2 = basis2 / np.linalg.norm(basis2)

        self.basis = np.matrix([basis1, basis2, normal_vector]).transpose()

        print("self.basis: {}".format(self.basis))

        # projector to orthogonal adjucent of the planes normal
        self.projector = self.basis*self.basis.transpose()
        # projector_normal = np.eye(3, dtype=float) - self.basis*self.basis.transpose()

    def projectOn2D(self, point_in):

        point = np.matrix([point_in - self.zero_point]).transpose()

        # project the point on the affine subspace of the plane
        projection = self.projector*point

        # express in the orthonormal basis of the subspace 
        projection = np.dot(np.linalg.inv(self.basis), projection)

        return Point2D(projection[0, 0], projection[1, 0])

# #{ class Source

class Source:

    position = [] # [m, m, m]
    activity = 0 # [Bc]
    energy = 0 # [eV]

    def __init__(self, energy, activity, position):

        self.energy = energy
        self.activity = activity
        self.position = position

# #} end of class Source

# #{ class Detector

class Detector:

    material = materials.CdTe
    thickness = 0 # [m]
    position = [] # [m, m, m]
    size = 0.01408 # [mm]
    vertices = [] # in the frame of the sensor
    sp_vertices = []

    # sympy
    sp_front_plane = []
    sp_front_polygon = []

    sp_back_plane = []
    sp_back_polygon = []

    sp_sides_planes = []
    sp_sides_polygones = []

    sp_all_planes = []
    sp_all_polygones = []

    def __init__(self, material, thickness, position):

        self.material = material
        self.thickness = thickness
        self.position = position

        # give me the 4 front_vertices of the detector
        a = np.array([self.thickness/2, -self.size/2.0, self.size/2.0])
        b = np.array([self.thickness/2, self.size/2.0, self.size/2.0])
        c = np.array([self.thickness/2, self.size/2.0, -self.size/2.0])
        d = np.array([self.thickness/2, -self.size/2.0, -self.size/2.0])

        # create the sympy representation of the front vertices
        sp_a = Point3D(a)
        sp_b = Point3D(b)
        sp_c = Point3D(c)
        sp_d = Point3D(d)

        # create the sympy front plane
        self.sp_front_plane = MyPlane(sp_a, sp_b, sp_c)

        # project the front vertices onto the front plane
        sp_a_2d = self.sp_front_plane.projectOn2D(sp_a)
        sp_b_2d = self.sp_front_plane.projectOn2D(sp_b)
        sp_c_2d = self.sp_front_plane.projectOn2D(sp_c)
        sp_d_2d = self.sp_front_plane.projectOn2D(sp_d)

        print("sp_a_2d: {}".format(sp_a_2d))
        print("sp_b_2d: {}".format(sp_b_2d))
        print("sp_c_2d: {}".format(sp_c_2d))
        print("sp_d_2d: {}".format(sp_d_2d))

        # create the front polygon in the front plane
        self.sp_front_polygon = Polygon(sp_a_2d, sp_b_2d, sp_c_2d, sp_d_2d)

        # give me the 4 back_vertices of the detector
        e = np.array([-self.thickness/2, -self.size/2.0, self.size/2.0])
        f = np.array([-self.thickness/2, self.size/2.0, self.size/2.0])
        g = np.array([-self.thickness/2, self.size/2.0, -self.size/2.0])
        h = np.array([-self.thickness/2, -self.size/2.0, -self.size/2.0])

        # create the sympy representation of the front vertices
        sp_e = Point3D(e)
        sp_f = Point3D(f)
        sp_g = Point3D(g)
        sp_h = Point3D(h)

        # create the sympy back plane
        self.sp_back_plane = MyPlane(sp_e, sp_f, sp_g)

        # create the sympy side planes
        self.sp_sides_planes.append(MyPlane(sp_a, sp_e, sp_h))
        self.sp_sides_planes.append(MyPlane(sp_b, sp_f, sp_e))
        self.sp_sides_planes.append(MyPlane(sp_c, sp_g, sp_f))
        self.sp_sides_planes.append(MyPlane(sp_d, sp_h, sp_c))

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
source_point = Point3D(source.position)
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
for i in range(0, 100):

   hit_point = sample_detector(detector_1)
   hit_points.append(hit_point) 

   ray = Ray3D(source_point, Point3D(hit_point))

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
        ax.scatter(point[0], point[1], point[2], color='g', marker='o')

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
