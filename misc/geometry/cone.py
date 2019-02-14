import numpy as np
import math as m

import solid_angle
from pyquaternion import Quaternion
from raytracing import Ray
import geometry.solid_angle

class Cone:

    def __init__(self, origin=np.array([0, 0, 0]), direction=np.array([1.0, 0, 0]), angle=m.pi/6.0):

        self.origin = origin
        self.direction = direction/np.linalg.norm(direction)
        self.angle = angle

        self.cone_axis_projector = np.dot(np.matrix([direction]).transpose(), np.matrix([direction]))

        self.cone_ray = Ray(self.origin, self.direction, 0)

    # #{ projectPoint()
    
    def projectPoint(self, point):
    
        # shift the point by the origin
        # so we can use the cone_axis_projector as an orthogonal projector
        # on a subspace
        point_vec = point - self.origin
        point_axis_angle = solid_angle.vector_angle(point_vec, self.direction)
    
        # project the point to the axis, and shift it back
        axis_projection = np.dot(self.cone_axis_projector, point_vec) + self.origin
    
        axis_rot = np.cross(self.direction, point_vec)
        axis_rot = axis_rot/np.linalg.norm(axis_rot)
    
        my_quaternion = Quaternion(axis=axis_rot, angle=self.angle - point_axis_angle)
    
        point_on_cone = my_quaternion.rotate(point_vec) + self.origin
    
        vec_point_on_cone = point_on_cone - self.origin
        vec_point_on_cone = vec_point_on_cone/np.linalg.norm(vec_point_on_cone)
    
        beta = self.angle - point_axis_angle
    
        if point_axis_angle < self.angle:
    
            ortho_projection = self.origin + vec_point_on_cone*m.cos(beta)*np.linalg.norm(point_vec)
            color = 'rgb(255, 0, 0)'
    
        elif point_axis_angle >= self.angle and point_axis_angle <= m.pi - self.angle:
    
            ortho_projection = self.origin + vec_point_on_cone*m.cos(point_axis_angle-self.angle)*np.linalg.norm(point_vec)
            color = 'rgb(0, 255, 0)'
        else:
    
            ortho_projection = self.origin
            color = 'rgb(0, 0, 255)'
    
        return np.array([ortho_projection[0], ortho_projection[1], ortho_projection[2]]), color
    
    # #} end of projectPoint()

    # #{ projectPointOnPlane()
    
    def projectPointOnPlane(self, point, plane):
    
        # find the intersection of the plane and the cone's axis
        intersect = plane.intersectionRay(self.cone_ray)
    
        if not isinstance(intersect, np.ndarray):
            return False

        # point to polygon axis angle
        point_axis_angle = geometry.solid_angle.vector_angle(self.direction, point - self.origin) 

        if point_axis_angle >= self.angle:
    
            alpha = point_axis_angle - self.angle
            beta = geometry.solid_angle.vector_angle(intersect - point, self.origin - point) 
            gamma = m.pi - alpha - beta

            print("point_axis_angle: {}".format(point_axis_angle))
            print("self.angle: {}".format(self.angle))

            c = np.linalg.norm(point - self.origin) 
            a = m.sin(alpha)*(c/m.sin(gamma))

            print("a: {}".format(a))

            vec = intersect - point
            # vec = vec/np.linalg.norm(vec) 
    
            return point + 1.0*vec
    
        else:

            vec = intersect - point
    
            # the ratio of angle
            ratio = (point_axis_angle - self.angle)/(point_axis_angle)

            return False

            return point + ratio*vec
    
    # #} end of projectPointOnPlane()
