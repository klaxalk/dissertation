import numpy as np
import math as m
import photon_attenuation.materials as materials
import geometry.solid_angle
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random

class Source:

    position = [] # [m, m, m]
    activity = 0 # [Bc]
    energy = 0 # [eV]

    def __init__(self, energy, activity, position):

        self.energy = energy
        self.activity = activity
        self.position = position

class Detector:

    material = materials.CdTe
    thickness = 0 # [m]
    position = [] # [m, m, m]
    size = 0.01408 # [mm]
    vertices = [] # in the frame of the sensor

    def __init__(self, material, thickness, position):

        self.material = material
        self.thickness = thickness
        self.position = position

        # give me the 4 front_vertices of the detector
        a = np.array([self.thickness/2, -self.size/2.0, self.size/2.0])
        b = np.array([self.thickness/2, self.size/2.0, self.size/2.0])
        c = np.array([self.thickness/2, self.size/2.0, -self.size/2.0])
        d = np.array([self.thickness/2, -self.size/2.0, -self.size/2.0])

        self.vertices.append(a)
        self.vertices.append(b)
        self.vertices.append(c)
        self.vertices.append(d)

    def getWorldVertices(self):

        # detector front_vertices in the world frame
        a = self.position + self.vertices[0]
        b = self.position + self.vertices[1]
        c = self.position + self.vertices[2]
        d = self.position + self.vertices[3]

        return [a, b, c, d]

# define the source and the detector
source = Source(662000.0, 1e9, np.array([10.0, 0.0, 0.0]))
detector_1 = Detector(materials.Si, 0.001, np.array([0, 0, 0]))
detector_2 = Detector(materials.CdTe, 0.001, np.array([-0.002, 0, 0]))

[a, b, c, d] = detector_1.getWorldVertices()
[e, f, g, h] = detector_2.getWorldVertices()

# calculate the relative activity to the detector
detector_solid_angle = geometry.solid_angle.quadrilateral_solid_angle(a, b, c, d, source.position)

# apparent activity
aparent_activity = source.activity*(detector_solid_angle/(4*m.pi))

def sample_detector(detector):

    # get the detector vertices
    [a, b, c, d] = detector.getWorldVertices()

    # generate the convex combination of the vertices
    k = []
    k.append(random.uniform(0.0, 1.0))
    k.append(random.uniform(0.0, 1.0))
    k.append(random.uniform(0.0, 1.0))

    k = np.sort(k)

    sampled_point = a*k[0] + b*(k[1] - k[0]) + c*(k[2] - k[1]) + d*(1.0 - k[2])

    return sampled_point
    # return

hit_points = []
for i in range(0, 1000):

   hit_points.append(sample_detector(detector_1)) 

# plotting
def plot_everything(*args):

    fig = plt.figure(1)
    ax  = fig.add_subplot(111, projection = '3d')

    # plot the detector
    xs = [a[0], b[0], c[0], d[0], a[0]]
    ys = [a[1], b[1], c[1], d[1], a[1]]
    zs = [a[2], b[2], c[2], d[2], a[2]]
    plt.plot(xs, ys, zs)

    xs = [e[0], f[0], g[0], h[0], e[0]]
    ys = [e[1], f[1], g[1], h[1], e[1]]
    zs = [e[2], f[2], g[2], h[2], e[2]]
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
