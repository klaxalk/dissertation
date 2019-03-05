import os

import time
import copy

import random
random.seed(1005)

import math as m
import numpy as np
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
from geometry.cone import Cone
# from geometry.polygon3d import Polygon3D
from geometry.rectangle3d import Rectangle3D

simulate_energy_noise = True
simulate_pixel_uncertainty = True
rich_plot = True
plot_res_x = 1920
plot_res_y = 1000

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
        # self.front = Polygon3D([a, b, c, d])
        self.front = Rectangle3D([a, b, c, d])

        # BACK

        # give me the 4 back_vertices of the detector
        e = np.array([-self.thickness/2, -self.size/2.0, self.size/2.0]) + position
        f = np.array([-self.thickness/2, self.size/2.0, self.size/2.0]) + position
        g = np.array([-self.thickness/2, self.size/2.0, -self.size/2.0]) + position
        h = np.array([-self.thickness/2, -self.size/2.0, -self.size/2.0]) + position

        # create the sympy back plane
        # self.back = Polygon3D([e, f, g, h])
        self.back = Rectangle3D([e, f, g, h])

        # orthogonal basis of the detector
        self.b1 = np.array([0, self.size/2.0 - (-self.size/2.0), self.size/2.0 - self.size/2.0])
        self.b2 = np.array([0, -self.size/2.0 - (-self.size/2.0), -self.size/2.0 - self.size/2.0])

        # SIDES

        # create the sympy side planes
        # self.sides.append(Polygon3D([a, e, h, d]))
        # self.sides.append(Polygon3D([b, f, e, a]))
        # self.sides.append(Polygon3D([c, g, f, b]))
        # self.sides.append(Polygon3D([d, h, c, g]))
        self.sides.append(Rectangle3D([a, e, h, d]))
        self.sides.append(Rectangle3D([b, f, e, a]))
        self.sides.append(Rectangle3D([c, g, f, b]))
        self.sides.append(Rectangle3D([d, h, c, g]))

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
source = Source(611000.0, 50*1e9, np.array([0.1, 0.0, 0.0]))
source_distance = np.linalg.norm(source.position)
source_point = source.position
detector_1 = Detector(materials.Si, 0.001, np.array([0, 0, 0]))
detector_2 = Detector(materials.CdTe, 0.002, np.array([-0.005, 0, 0]))

a = np.array([-0.1, 2.0*source.position[0], 1.0*source.position[2]])
b = np.array([-0.1, -2.0*source.position[0], 1.0*source.position[2]])
c = np.array([2.0*source.position[0], -2.0*source.position[0], 1.0*source.position[2]])
d = np.array([2.0*source.position[0], 2.0*source.position[0], 1.0*source.position[2]])

# ground_polygon = Polygon3D([a, b, c, d])
ground_polygon = Rectangle3D([a, b, c, d])

[a1, b1, c1, d1, e1, f1, g1, h1] = detector_1.getVertices()
[a2, b2, c2, d2, e2, f2, g2, h2] = detector_2.getVertices()

# calculate the relative activity to the detector
detector_solid_angle = geometry.solid_angle.quadrilateral_solid_angle(a1, b1, c1, d1, source.position)

# apparent activity
aparent_activity = source.activity*(detector_solid_angle/(4*m.pi))

print("aparent_activity: {}".format(aparent_activity))

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
n_particles = 30000

py_traces = []

hypo = [np.array([150.0, 150.0, source.position[2]])]

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
         cone_origin = scatterer_mid_point
         cone_direction = cone_origin - absorber_mid_point
     else: 
         cone_origin = scattered_ray.rayPoint
         cone_direction = cone_origin - absorption_point

     # normalize cone direction
     cone_direction = cone_direction/np.linalg.norm(cone_direction)

     # estimate the scattering angle theta
     theta_estimate = physics.getComptonAngle(electron_energy+e_noise, scattered_ray.energy+f_noise)

     print("theta_estimate: {}".format(theta_estimate))

     # swap the cone, if its angle is > 90 deg
     if theta_estimate > m.pi/2:
         theta_estimate = m.pi - theta_estimate
         cone_direction *= -1.0

     print("theta_estimate: {}".format(theta_estimate))

     cone = Cone(cone_origin, cone_direction, theta_estimate)

     # plotting the cone
     # sample the end points of the cone (the cone is pointing along the x axis)
     u = np.linspace(source_distance, 0, 2)
     
     # sweeping along the cone's axis (vector v1)
     v = np.linspace(-np.pi, np.pi, 180)
     
     v1 = np.array([1.0, 0.0, 0.0])
     
     # rotate the cone to its final orientation
     my_axis = np.cross(v1, cone.direction) # axis of rotation
     
     my_quaternion = Quaternion(axis=my_axis, angle=geometry.solid_angle.vector_angle(v1, cone.direction))
     
     # prepare the meshgrid to plot the cone
     u, v = np.meshgrid(u,v)
     
     u = u.flatten()
     v = v.flatten()
     
     # # define the 3D coorinates
     x = u*m.cos(cone.angle)
     y = u*m.sin(cone.angle)*np.sin(v)
     z = u*m.sin(cone.angle)*np.cos(v)

     # triangulate the surface
     points2D = np.vstack([u,v]).T
     tri = Delaunay(points2D)
     simplices = tri.simplices
     I,J,K=([triplet[c] for triplet in simplices] for c in range(3))
     
     # rotate and move the cone to its final position and orientation
     for index,point in enumerate(x):
     
         new_point = my_quaternion.rotate(np.array([x[index], y[index], z[index]]))
         x[index] = new_point[0] + cone.origin[0]
         y[index] = new_point[1] + cone.origin[1]
         z[index] = new_point[2] + cone.origin[2]

         # if z[index] < source.position[2]:
         #    z[index] = source.position[2]
     
     # create the cone's meshgrid
     py_traces.append(go.Mesh3d(x=x,
                         y=y,
                         z=z,
                         i=I,
                         j=J,
                         k=K,
                         name='',
                         opacity=0.1,
                         ))

   # #} end of plot the cone

   # #{ estimation

   axis_ground_proj = cone.projectPointOnPlane(hypo[-1], ground_polygon.plane)

   if isinstance(axis_ground_proj, np.ndarray):

       print("axis_ground_proj: {}".format(axis_ground_proj))

       coef = 0.9
       hypo.append(axis_ground_proj*(1.0-coef) + coef*hypo[-1])
       print("hypo[-1]: {}".format(hypo[-1]))
       # hypo.append(axis_ground_proj)

       # py_traces.append(go.Scatter3d(
       #      x=[axis_ground_proj[0]], y=[axis_ground_proj[1]], z=[axis_ground_proj[2]],
       #      marker=dict(
       #          size=3,
       #          color='rgb(255, 0, 0)',
       #      ),
       #      name=''
       # ))

   # hypo_proj, color = cone.projectPoint(hypo[-1])
   # coef = 0.9
   # hypo.append(hypo_proj*(1.0-coef) + coef*hypo[-1])

   # #} end of estimation

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

    # # #{ plot the ground

    # [xs, ys, zs] = ground_polygon.plotVertices()

    # py_traces.append(go.Mesh3d(
    #      x = xs,
    #      y = ys,
    #      z = zs,
    #      i = [0, 0],
    #      j = [1, 2],
    #      k = [2, 3],
    #      color='rgb(0, 256, 0)',
    #      opacity=0.3,
    #      name = 'Ground plane',
    #      showscale = True
    # ))

    # # #} end of plot ground

    # # #{ plot the filtration

    # xs = [a[0] for a in hypo]
    # ys = [a[1] for a in hypo]
    # zs = [a[2] for a in hypo]
    
    # py_traces.append(go.Scatter3d(
    #      x=xs, y=ys, z=zs,
    #      marker=dict(
    #          size=3,
    #          color='rgb(0, 0, 0)',
    #      ),
    #      line=dict(
    #          color='rgb(0, 0, 0)',
    #          width=1
    #      ),
    #      name=''
    # ))
    
    # # py_traces.append(go.Scatter3d(
    # #      x=[hypo_proj[0]], y=[hypo_proj[1]], z=[hypo_proj[2]],
    # #      marker=dict(
    # #          size=3,
    # #          color=color,
    # #      ),
    # #      line=dict(
    # #          color='rgb(0, 0, 255)',
    # #          width=0
    # #      ),
    # #      name=''
    # # ))
    
    # # #} end of plot the filtration

    # #{ plotly layout
    
    layout = dict(
        width=plot_res_x,
        height=plot_res_y,
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
