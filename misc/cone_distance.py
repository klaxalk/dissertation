import sympy as sy
import numpy as np
import math as m
import scipy as sc
import scipy.optimize as op
import timeit

sy_x = sy.symbols('x y z')
sy_center = sy.symbols('a b c')
sy_direction = sy.symbols('o p q')
sy_theta = sy.symbols('Theta')

centers = [[0, 0, 0], [0, 5, 0], [0, 2.5, 0], [5, 0, 0]]
directions = [[1, 0, 0], [1, 0, 0], [0.707, -0.707, 0], [-1, 0, 0]]
thetas = [m.pi/4.0, m.pi/8.0, m.pi/4.0, m.pi/4]

# centers = [[0, 0, 0], [-5, 0, 0]]
# directions = [[1, 0, 0], [-1, 0, 0]]
# thetas = [m.pi/4.0, m.pi/4.0]

norm = sy.Pow(sy.Pow(sy_x[0] - sy_center[0], 2) + sy.Pow(sy_x[1] - sy_center[1], 2) + sy.Pow(sy_x[2] - sy_center[2], 2), 0.5)
cone_dist = norm*sy.sin(sy.acos(((sy_x[0]-sy_center[0])*sy_direction[0] + (sy_x[1]-sy_center[1])*sy_direction[1] + (sy_x[2]-sy_center[2])*sy_direction[2])/(norm)) - sy_theta)

sy.pprint(sy.simplify(cone_dist), use_unicode=True)

def coneDist(center, direction, theta):

    out = cone_dist

    for idx,x in enumerate(sy_center):
        out = out.subs(x, center[idx]) 

    for idx,x in enumerate(sy_direction):
        out = out.subs(x, direction[idx]) 

    out = out.subs(sy_theta, theta) 

    return out

def multiConeDist(cc, dd, tt):
    
    out = 1
    for idx,i in enumerate(cc):

        out = out + sy.Pow(coneDist(cc[idx], dd[idx], tt[idx]), 2)

    return out

def constraints(cc, dd):

    cons = []

    for idx,i in enumerate(cc):

        center = cc[idx]
        direction = dd[idx]

        cons.append({'type': 'ineq', 'fun': lambda x: direction[0]*x[0] + direction[1]*x[1] + direction[2]*x[2] - (direction[0]*center[0] + direction[1]*center[1] + direction[2]*center[2])})

    return cons

# test = coneDist(centers[1], directions[1], thetas[1])
# test = test.subs(sy_x[0], 0)
# test = test.subs(sy_x[1], 0)
# test = test.subs(sy_x[2], 0)
# print("test: {}".format(test))

func = multiConeDist(centers, directions, thetas)
sy.pprint(sy.simplify(func), use_unicode=True)
bnds = ((-100,-100,-100), (100,100,100))
cons = constraints(centers, directions)
print("cons: {}".format(cons))

# sy.pprint(sy.simplify(lamb_diff), use_unicode=True)

# sy_x = sy.symbols('x y')
# # func = sy.sin(x[0]) + sy.cos(x[1])
# func = sy.Pow(sy_x[0], 2) + sy.Pow(sy_x[1], 2)

# bnds = ((1,2), (2,2))
# cons = ({'type': 'ineq', 'fun': lambda x: - x[0] + 5},
#        {'type': 'ineq', 'fun': lambda x: x[0] - 2})

J = [func.diff(var) for var in sy_x]
H = [[func.diff(var1).diff(var2) for var1 in sy_x] for var2 in sy_x]

def subs_scal(scal, x, values):

    for idx1,var in enumerate(x):

        scal = scal.subs(var, values[idx1]) 

    return sy.N(scal)

def subs_vec(vec, x, values):

    out = np.zeros((len(vec)))

    for idx1,cell in enumerate(vec):

        out[idx1] = subs_scal(cell, x, values)

    return out

def subs_mat(mat, x, values):

    out = np.zeros((size(vec)))

    for idx1,vec in enumerate(mat):

        out[idx1] = subs_vec(vec, x, values)

    return out

# print("J: {}".format(J))
# print("H: {}".format(H))

def f(point):
    return subs_scal(func, sy_x, point)

def jacobian(point):
    return np.array(subs_vec(J, sy_x, point))

def hessian(point):
    return np.matrix(subs_mat(H, sy_x, point))

# test = [2, -1]

# print("test: {}".format(test))
# print("jacobian([1, 1]): {}".format(jacobian(test)))
# print("test: {}".format(test))

# res = minimize(f, [20, 10], method="CG")
# res = minimize(f, [1, 1], method="Newton-CG", jac=jacobiang
start = timeit.timeit()
res = op.minimize(f, [50, 50, 50], method="SLSQP", jac=jacobian, options={'gtol': 1e-6, 'disp': True}, constraints=cons)
end = timeit.timeit()
# res = op.minimize(f, [50, 50, 50], method="SLSQP", jac=jacobian, options={'gtol': 1e-6, 'disp': True}, constraints=cons)
# res = op.minimize(f, [1, 1], method="SLSQP", jac=jacobian, bounds=bnds)
# res = minimize(f, [1, 1], method="SLSQP", jac=jacobian, hess=hessian)

print("res.x: {}".format(res))

print("end - start: {}".format(end - start))

# sy.pprint(sy.simplify(func), use_unicode=True)
# sy.pprint(sy.simplify(J), use_unicode=True)
# sy.pprint(sy.simplify(H), use_unicode=True)
