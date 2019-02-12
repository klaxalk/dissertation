import os

import time
import copy

import random
random.seed(1005)

import math as m
import numpy as np
import shapely as sh
from shapely.geometry import Polygon
from shapely import speedups
speedups.enable()
from pyquaternion import Quaternion

# plotly
import plotly.offline as py
import plotly.graph_objs as go
from scipy.spatial import Delaunay # for triangulation of the cones

# my stuff
import geometry.solid_angle
import photon_attenuation.materials as materials
import photon_attenuation.physics as physics
from geometry.raytracing import Plane, Ray

simulate_energy_noise = True
simulate_pixel_uncertainty = True
rich_plot = True

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

        # orthogonal basis of the detector
        self.b1 = np.array([0, self.size/2.0 - (-self.size/2.0), self.size/2.0 - self.size/2.0])
        self.b2 = np.array([0, -self.size/2.0 - (-self.size/2.0), -self.size/2.0 - self.size/2.0])

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

    def plotVertices(self):

        [a, b, c, d, e, f, g, h] = self.getVertices()

        xx = [a[0], b[0], c[0], d[0], a[0], e[0], f[0], b[0], f[0], g[0], c[0], g[0], h[0], d[0], h[0], e[0]]
        yy = [a[1], b[1], c[1], d[1], a[1], e[1], f[1], b[1], f[1], g[1], c[1], g[1], h[1], d[1], h[1], e[1]]
        zz = [a[2], b[2], c[2], d[2], a[2], e[2], f[2], b[2], f[2], g[2], c[2], g[2], h[2], d[2], h[2], e[2]]

        return xx, yy, zz

    def plotPixelVertices(self, x, y):

        x = x/256.0
        y = y/256.0

        dpx = 1.0/256.0

        thickness_offset = np.array([self.thickness/2, 0, 0])

        a = self.vertices[0] + self.b1*x + self.b2*y
        b = self.vertices[0] + self.b1*(x+dpx) + self.b2*y
        c = self.vertices[0] + self.b1*(x+dpx) + self.b2*(y+dpx)
        d = self.vertices[0] + self.b1*x + self.b2*(y+dpx)

        e = a - 2*thickness_offset
        f = b - 2*thickness_offset
        g = c - 2*thickness_offset
        h = d - 2*thickness_offset

        xx = [a[0], b[0], c[0], d[0], a[0], e[0], f[0], b[0], f[0], g[0], c[0], g[0], h[0], d[0], h[0], e[0]]
        yy = [a[1], b[1], c[1], d[1], a[1], e[1], f[1], b[1], f[1], g[1], c[1], g[1], h[1], d[1], h[1], e[1]]
        zz = [a[2], b[2], c[2], d[2], a[2], e[2], f[2], b[2], f[2], g[2], c[2], g[2], h[2], d[2], h[2], e[2]]

        return xx, yy, zz

    def getPixelMidPoint(self, x, y):

        x = x/256.0
        y = y/256.0

        dpx = 0.5/256.0

        thickness_offset = np.array([self.thickness/2, 0, 0])

        a = self.vertices[0] + self.b1*(x+dpx) + self.b2*(y+dpx) - thickness_offset

        return a

    def point2pixel(self, point):

        x = (point[1] - self.vertices[0][1])
        y = (point[2] - self.vertices[0][2])

        x = int(m.floor((x/self.b1[1])*256))
        y = int(m.floor((y/self.b2[2])*256))

        return x, y

# #} end of class Detector

# define the source and the detector
source = Source(621000.0, 1e9, np.array([0.1, 0.1, -0.05]))
source_distance = np.linalg.norm(source.position)
source_point = source.position
detector_1 = Detector(materials.Si, 0.001, np.array([0, 0, 0]))
detector_2 = Detector(materials.CdTe, 0.002, np.array([-0.005, 0, 0]))

[a1, b1, c1, d1, e1, f1, g1, h1] = detector_1.getVertices()
[a2, b2, c2, d2, e2, f2, g2, h2] = detector_2.getVertices()

# calculate the relative activity to the detector
detector_solid_angle = geometry.solid_angle.quadrilateral_solid_angle(a1, b1, c1, d1, source.position)

# apparent activity
aparent_activity = source.activity*(detector_solid_angle/(4*m.pi))

# prepare the Compton scattering cross section and probability density for the source's energy
cs_cross_section = physics.comptonCrossSection(physics.conversions.energy_ev_to_J(source.energy))
cs_density_indeces, cs_density = physics.cs_distribution_function(detector_1.material, source.energy)

# #{ comptonScattering()

def comptonScattering(from_point, to_point, energy, material, cs_cross_section, cs_density):

    distance = np.linalg.norm(to_point - from_point)

    # calculate the probability of the scattering

    # draw and decide wether it should happend
    prob_cs = 1.0 - np.exp(-material.electron_density * cs_cross_section * distance)

    if random.uniform(0.0, 1.0) > prob_cs:
        return False, 0, 0

    # calculate the point of scattering in the detector
    # scattering_point = (from_point + to_point)/2.0
    position_weight = random.uniform(0.0, 1.0)
    scattering_point = position_weight*from_point + (1 - position_weight)*to_point

    # calculate the azimuthal and radial angle
    phi = random.uniform(0.0, 2.0*m.pi)
    theta = cs_density[int(m.floor(random.uniform(0.0, len(cs_density))))]
    # theta = 10.0*(m.pi/180.0)

    # calculate the point on a unit sphere
    x1 = m.cos(theta)
    y1 = m.sin(theta)*m.sin(phi)
    z1 = m.sin(theta)*m.cos(phi)

    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([x1, y1, z1])
    my_axis = np.cross(v1, v2)

    original_direction = to_point - from_point
    original_direction = original_direction/np.linalg.norm(original_direction)

    my_axis2 = np.cross(v1, original_direction)
    angle2 = geometry.solid_angle.vector_angle(v1, original_direction)

    my_quaternion_1 = Quaternion(axis=my_axis2, angle=angle2)
    my_axis = my_quaternion_1.rotate(my_axis)

    # print("pes: {}".format(geometry.solid_angle.vector_angle(v1, v2)))

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

    return new_ray, electron_energy, theta

# #} end of 

# #{ sampleDetector()

def sampleDetector(detector):

    [a, b, c, d, e, f, g, h] = detector.getVertices()

    ab = b-a
    ad = d-a

    k1 = random.uniform(1e-2, 1.0-1e-2)
    k2 = random.uniform(1e-2, 1.0-1e-2)

    return a + ab*k1 + ad*k2

# #} end of sampleDetector()

# sample the 1st detector
n_particles = 10000

py_traces = []

time_start = time.time()

for i in range(0, n_particles):

   point = sampleDetector(detector_1)

   ray = Ray(source_point, point, source.energy)

   # intersection with the back side of the 1st detector
   intersect1_second = detector_1.back.intersection(ray)

   # no collision with the back face
   if not isinstance(intersect1_second, np.ndarray):
       continue

       # check intersection with the sides
       for i,side in enumerate(detector_1.sides):

           intersect1_second = side.intersection(ray)

           if isinstance(intersect1_second, np.ndarray):
               break

   # if the ray came from the oposite direction, discard it
   if np.linalg.norm(intersect1_second - source.position) < np.linalg.norm(point - source.position):
       continue

   # if there is not a collission with any other facet of the detector
   if not isinstance(intersect1_second, np.ndarray):
       print("! no intersection with the back/side face of the first detector")
       continue

   # calculate the length of the intersection with the detector
   intersection_len = np.linalg.norm(point - intersect1_second)

   # scatter the ray
   scattered_ray, electron_energy, theta = comptonScattering(point, intersect1_second, source.energy, detector_1.material, cs_cross_section, cs_density)

   # if not scattering happened, just leave
   if not isinstance(scattered_ray, Ray):
       continue

   # check the collision with the other detector's front side
   intersect2_first = detector_2.front.intersection(scattered_ray)

   # if there is no intersection with the front face, just leave
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

   # if there is no intersection with the other facets of the detector
   if not isinstance(intersect2_second, np.ndarray):
       print("!! no intersection with the back/side face of the second detector")
       continue

   # calculate the photo-electric cross section for the scattered photon
   pe_cross_section = [physics.pe_cs_gavrila_pratt_simplified(mat, scattered_ray.energy) for mat in detector_2.material.elements]
   # and the effective thickness of the material for the PE effect
   pe_thickness = np.linalg.norm(intersect2_second - intersect2_first)

   prob_pe = 1.0
   for index,cross_section in enumerate(pe_cross_section):
       prob_pe *= np.exp(-detector_2.material.element_quantities[index] * detector_2.material.molecular_density * cross_section * pe_thickness)
   prob_pe = 1.0 - prob_pe

   # do a coin toss for the photo-electric effect
   if random.uniform(0.0, 1.0) > prob_pe:
       continue

   # sample the interraction point in the 2nd detector
   position_weight = random.uniform(0.0, 1.0)
   absorption_point = position_weight*intersect2_first + (1 - position_weight)*intersect2_second

   print("Complete compton: p_energy: {}, absorber_thickness: {}, prob_pe: {}".format(scattered_ray.energy, pe_thickness, prob_pe))

   # #{ plot the rays from the source to the scatterer
   
   if rich_plot:
     xs = [source.position[0], point[0]]
     ys = [source.position[1], point[1]]
     zs = [source.position[2], point[2]]
     
     py_traces.append(go.Scatter3d(
         x=xs, y=ys, z=zs,
         mode='lines',
         marker=dict(
             size = 0,
             color = 'rgba(0, 0, 0, 1.0)',
         ),
         line=dict(
             color='rgb(128, 128, 128)',
             width=3
         ),
         showlegend=False
     ))
   
   # #} end of plot the rays from the source to the scatterer

   # #{ plot the rays from the scatterer to the absorber
   
   if rich_plot:
     xs = [scattered_ray.rayPoint[0], absorption_point[0]]
     ys = [scattered_ray.rayPoint[1], absorption_point[1]]
     zs = [scattered_ray.rayPoint[2], absorption_point[2]]
     
     py_traces.append(go.Scatter3d(
         x=xs, y=ys, z=zs,
         mode='lines',
         marker=dict(
             size = 0,
             color = 'rgba(0, 0, 0, 1.0)',
         ),
         line=dict(
             color='rgb(0, 0, 255)',
             width=3
         ),
         showlegend=False
     ))
   
   # #} end of plot the rays from the scatterer to the absorber

   # #{ plot the scattering pixel
   
   if rich_plot:
     x, y = detector_1.point2pixel(scattered_ray.rayPoint)
     xs, ys, zs = detector_1.plotPixelVertices(x, y)
     scatterer_mid_point = detector_1.getPixelMidPoint(x, y)
     
     py_traces.append(go.Scatter3d(
         x=xs, y=ys, z=zs,
         mode='lines',
         marker=dict(
             size = 0,
             color = 'rgba(0, 0, 0, 1.0)',
         ),
         line=dict(
             color='rgb(255, 0, 0)',
             width=3
         ),
         showlegend=False
     ))
     
   # #} end of plot the absorption pixels

   # #{ plot the absorption pixel
   
   if rich_plot:
     x, y = detector_2.point2pixel(absorption_point)
     xs, ys, zs = detector_2.plotPixelVertices(x, y)
     absorber_mid_point = detector_2.getPixelMidPoint(x, y)
     
     py_traces.append(go.Scatter3d(
         x=xs, y=ys, z=zs,
         mode='lines',
         marker=dict(
             size = 0,
             color = 'rgba(0, 0, 0, 1.0)',
         ),
         line=dict(
             color='rgb(255, 0, 0)',
             width=3
         ),
         showlegend=False
     ))
   
   # #} end of plot the absorption pixels

   # #{ plot the cone axis
   
   if rich_plot:
     xs = [scatterer_mid_point[0], absorber_mid_point[0]] 
     ys = [scatterer_mid_point[1], absorber_mid_point[1]] 
     zs = [scatterer_mid_point[2], absorber_mid_point[2]] 
     
     py_traces.append(go.Scatter3d(
         x=xs, y=ys, z=zs,
         mode='lines',
         marker=dict(
             size = 0,
             color = 'rgba(0, 0, 0, 1.0)',
         ),
         line=dict(
             color='rgb(0, 255, 0)',
             width=3
         ),
         showlegend=False
     ))
   
   # #} end of plot the cone axis

   # #{ plot the cone
   
   if rich_plot:

     # add noise to the energies
     if simulate_energy_noise:
         e_noise = random.gauss(0, 7000)
         f_noise = random.gauss(0, 7000)
     else:
         e_noise = 0
         f_noise = 0

     # calculate the direction of the cone
     if simulate_pixel_uncertainty:
         cone_direction = scatterer_mid_point - absorber_mid_point
     else: 
         cone_direction = scattered_ray.rayPoint - absorption_point

     # normalize the direction vector
     cone_direction = cone_direction/np.linalg.norm(cone_direction)

     # estimate the scattering angle theta
     theta_estimate = physics.getComptonAngle(electron_energy+e_noise, scattered_ray.energy+f_noise)

     # swap the cone, if its angle is > 90 deg
     if theta_estimate > m.pi/2:
         theta_estimate = m.pi - theta_estimate
         cone_direction *= -1.0
     
     # sample the end points of the cone (the cone is pointing along the x axis)
     u = np.linspace(source_distance, 0, 2)
     
     # sweeping along the cone's axis (vector v1)
     v = np.linspace(-np.pi, np.pi, 180)
     
     v1 = np.array([1.0, 0.0, 0.0])
     
     # rotate the cone to its final orientation
     my_axis = np.cross(v1, cone_direction) # axis of rotation
     
     my_quaternion = Quaternion(axis=my_axis, angle=geometry.solid_angle.vector_angle(v1, cone_direction))
     
     # prepare the meshgrid to plot the cone
     u, v = np.meshgrid(u,v)
     
     u = u.flatten()
     v = v.flatten()
     
     # # define the 3D coorinates
     x = u*m.cos(theta_estimate)
     y = u*m.sin(theta_estimate)*np.sin(v)
     z = u*m.sin(theta_estimate)*np.cos(v)

     # x = u
     # y = u*np.cos(v)*m.tan(theta_estimate)
     # z = u*np.sin(v)*m.tan(theta_estimate)
     
     # triangulate the surface
     points2D = np.vstack([u,v]).T
     tri = Delaunay(points2D)
     simplices = tri.simplices
     I,J,K=([triplet[c] for triplet in simplices] for c in range(3))
     
     # rotate and move the cone to its final position and orientation
     for index,point in enumerate(x):
     
         new_point = my_quaternion.rotate(np.array([x[index], y[index], z[index]]))
         x[index] = new_point[0] + scatterer_mid_point[0]
         y[index] = new_point[1] + scatterer_mid_point[1]
         z[index] = new_point[2] + scatterer_mid_point[2]

         if z[index] < source.position[2]:
            z[index] = source.position[2]
     
     # create the cone's meshgrid
     py_traces.append(go.Mesh3d(x=x,
                         y=y,
                         z=z,
                         i=I,
                         j=J,
                         k=K,
                         name='',
                         opacity=0.1
                         ))

   # #} end of plot the cone

duration = (time.time() - time_start)
print("duration: {}".format(duration))

# plotting
def plot_everything(*args):

    # #{ plot detector 1

    xs, ys, zs = detector_1.plotVertices()

    py_traces.append(go.Scatter3d(
        x=xs, y=ys, z=zs,
        marker=dict(
            size=1,
            color='rgb(0, 0, 0)',
        ),
        line=dict(
            color='rgb(0, 0, 255)',
            width=3
        ),
        name='Scatterer'
    ))

    # #} end of plot detector 1

    # #{ plot detector 2

    xs, ys, zs = detector_2.plotVertices()

    py_traces.append(go.Scatter3d(
        x=xs, y=ys, z=zs,
        marker=dict(
            size=1,
            color='rgb(0, 0, 0)',
        ),
        line=dict(
            color='rgb(255, 0, 0)',
            width=2
        ),
        name='Absorber'
    ))

    # #} end of plot detector 2

    # #{ plot the ground
    
    a = np.array([-0.1, 2*source.position[0], 1.0*source.position[2]])
    b = np.array([-0.1, -2*source.position[0], 1.0*source.position[2]])
    c = np.array([2*source.position[0], -2*source.position[0], 1.0*source.position[2]])
    d = np.array([2*source.position[0], 2*source.position[0], 1.0*source.position[2]])

    xs = [a[0], b[0], c[0], d[0]]
    ys = [a[1], b[1], c[1], d[1]]
    zs = [a[2], b[2], c[2], d[2]]

    py_traces.append(go.Mesh3d(
         x = xs,
         y = ys,
         z = zs,
         i = [0, 0],
         j = [1, 2],
         k = [2, 3],
         color='rgb(0, 256, 0)',
         opacity=0.3,
         name = 'Ground plane',
         showscale = True
    ))

    # #} end of plot ground

    # #{ plotly layout
    
    layout = dict(
        width=1920,
        height=1080,
        margin=dict(
            r=20, l=10,
            b=10, t=10),
        # autosize=False,
        title='Compton camera',
        scene=dict(
           xaxis=dict(
                # gridcolor='rgb(255, 255, 255)',
                # zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            yaxis=dict(
                # gridcolor='rgb(255, 255, 255)',
                # zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            zaxis=dict(
                # gridcolor='rgb(255, 255, 255)',
                # zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(230, 230,230)'
            ),
            camera=dict(
                up=dict(
                    x=0,
                    y=0,
                    z=1
                ),
                center=dict(
                    x=-0.01,
                    y=0.0,
                    z=0.0,
                ),
                eye=dict(
                    x=-2.5,
                    y=0.1,
                    z=0.1
                )
            ),
            aspectratio = dict( x=1, y=1, z=1 ),
            aspectmode = 'data'
        ),
    )
    
    # #} end of plotly layout

    fig2 = dict(data=py_traces, layout=layout)

    filename = "compton{}{}.html".format("_energy_noise" if simulate_energy_noise else "", "_pixel_unc" if simulate_pixel_uncertainty else "")
    filename = "compton.html"

    py.plot(fig2, filename=filename, auto_open=False)

    print("finished with {}".format(filename))

pid = os.fork()
if pid == 0:
    plot_everything()
