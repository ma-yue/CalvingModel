"""
    This version solves the Stokes equation:
    - div( nu( grad(u) + grad(u).T ) - p I ) = f
    ie - nu div ( grad(u) + grad(u).T ) + grad(p) = f
    """

from dolfin import *
import numpy as np
import math
#import timeit
import os


# Physical Constants
steps = 365               # steps
n = 6.0                   # a constant
thickness = 800.0         # thickness of mesh
length = n*thickness      # length of mesh
depth = 700.0             # m ice under water
gridsize = 16.0           # m resolution

B0 = 15.77                # Pa*day^1/3
temp = 253.0              # K temperature
rho_i = 910               # kg/m3 ice density
rho_w = 1020              # kg/m3 seawater density
g = 9.8                   # m/s2
noblowup = 1E-15          # a small constant preventing inf viscosity values

# Load mesh from file
mesh = Mesh("mesh730.xml")

# Obtain x,z coordinates of vertices
x = mesh.coordinates()[:,0]
z = mesh.coordinates()[:,1]

# Define function spaces
degree = 1
scalar = FunctionSpace(mesh, "CG", degree)
vector = VectorFunctionSpace(mesh, "CG", degree)
system = vector * scalar

# Finding grounding line
grounding = length - gridsize
for v in vertices(mesh):
    k = v.index()
    if z[k] < 1E-3:
        for f in facets(v):
            if f.exterior():
                if x[k] > grounding:
                    grounding = x[k]

# Create mesh function over cell facets
boundary_parts = FacetFunction('size_t', mesh, 0)
boundary_parts.set_all(0)

# DOLFIN_EPS = 3e-16
# Mark right below water boundary as subdomain 5
right_below = []
for f in facets(mesh):
    count = 0
    for c in cells(f):
        count = count + 1
    if count == 1 and f.midpoint().y()<depth:
        right_below.append(f)
for rb in right_below:
    boundary_parts[rb] = 5

# Mark bottom boundary facets as subdomain 1
bottom = [f for f in facets(mesh)
          if f.midpoint().y()<1E-3]
for b in bottom:
    boundary_parts[b] = 1

# Mark left boundary as subdomain 3
left = [f for f in facets(mesh)
        if f.midpoint().x()<1E-3]
for l in left:
    boundary_parts[l] = 3

# No-slip boundary condition at bottom
bcb1 = DirichletBC(system.sub(0), Constant((0.0,0.0)), boundary_parts, 1)

# Free-slip boundary condition at bottom
bcb2 = DirichletBC(system.sub(0).sub(1), Constant(0.0), boundary_parts, 1)

# Free slip boundary condition on left hand side
bcl = DirichletBC(system.sub(0).sub(0), Constant(0.0), boundary_parts, 3)

# Traction-free boundary condition on surface
# do nothing

# Collect Dirichlet boundary conditions
bcs = [bcb2, bcl]


# Define new measures associated with the interior domains and
# exterior boundaries
ds = Measure("ds")[boundary_parts]


# Define strain rate tensor and viscosity
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
# see glacier dynamics van der veen P33
def nu(u,temp):
    return B0*exp(3155.0/temp - 0.16612/(273.39-temp)**1.17)*\
        (0.5*(epsilon(u)[0,0])**2 + (epsilon(u)[0,1])**2 + \
         0.5*(epsilon(u)[1,1])**2 + noblowup)**(-1.0/3.0)

# Define water pressure
#waterpressure = Expression(("A*B*(D - x[1])*(x[1] < D)*(fabs(x[0]-C) < DOLFIN_EPS | fabs(c*(x[0]-a)-(b-a)*x[1]) < DOLFIN_EPS | fabs(c*(x[0]-d)-(b-c)*x[1]) < DOLFIN_EPS)","0.0"),A = rho_w, B = g, C = length, D = depth, a = tipx - width/2., b = tipx, c = tipz, d = tipx + width/2.)
waterpressure = Expression(("A*B*(D - x[1])","0.0"),A = rho_w, B = g, D = depth)

# Define closest_vertex for a point
def closest_vertex(point):
    d = 2.*gridsize
    for v in vertices(mesh):
        deltax = x[v.index()]-point[0]
        deltaz = z[v.index()]-point[1]
        if abs(deltax) <= d and abs(deltaz) <= d:
            distance = math.hypot(deltax,deltaz)
            if distance < d:
                d = distance
                index = v.index()
    return index

# basal
bindex = []
basal = np.zeros([2,2])
bverts_x = []
bverts_z = []
for k in range(len(basal[0])):
    bindex.append([])
    bverts_x.append([])
    bverts_z.append([])
# surface
sindex = []
surface = np.zeros([2,2])
sverts_x = []
sverts_z = []
for k in range(len(surface[0])):
    sindex.append([])
    sverts_x.append([])
    sverts_z.append([])

# Read in saved crevasse paths
xxx = list()
zzz = list()
fhand = open('/Users/yuema/Documents/free_slip_800-700/NoRemesh/crevs730.out','r')
for i, line in enumerate(fhand):
    if i == 0:
        xxx = line.strip().split()
        xxx = map(float,xxx)
    if i == 1:
        zzz = line.strip().split()
        zzz = map(float,zzz)
        for k in range(len(xxx)):
            bindex[0].append(closest_vertex(Point(xxx[k],zzz[k])))
            bverts_x[0].append(xxx[k])
            bverts_z[0].append(zzz[k])
    if i == 2:
        xxx = line.strip().split()
        xxx = map(float,xxx)
    if i == 3:
        zzz = line.strip().split()
        zzz = map(float,zzz)
        for k in range(len(xxx)):
            bindex[1].append(closest_vertex(Point(xxx[k],zzz[k])))
            bverts_x[1].append(xxx[k])
            bverts_z[1].append(zzz[k])
    if i == 4:
        xxx = line.strip().split()
        xxx = map(float,xxx)
    if i == 5:
        zzz = line.strip().split()
        zzz = map(float,zzz)
        for k in range(len(xxx)):
            sindex[0].append(closest_vertex(Point(xxx[k],zzz[k])))
            sverts_x[0].append(xxx[k])
            sverts_z[0].append(zzz[k])
    if i == 6:
        xxx = line.strip().split()
        xxx = map(float,xxx)
    if i == 7:
        zzz = line.strip().split()
        zzz = map(float,zzz)
        for k in range(len(xxx)):
            sindex[1].append(closest_vertex(Point(xxx[k],zzz[k])))
            sverts_x[1].append(xxx[k])
            sverts_z[1].append(zzz[k])
fhand.close()

#start = timeit.default_timer()

for day in range(steps):
    # Define variational problem
    w = TrialFunction(system)
    y = TestFunction(system)
    (u,p) = split(w)
    (v,q) = split(y)
    u_k = interpolate(Constant((0.0,0.0)),vector)
    ux_k = interpolate(Constant(0.0),scalar)
    uz_k = interpolate(Constant(0.0),scalar)
    p_k = interpolate(Constant(0.0),scalar)

    f = Constant((0, -rho_i*g))   # with gravity
    h = CellSize(mesh)
    beta  = 0.2
    delta = beta*h*h
    a = (inner(nabla_grad(v), nu(u_k,temp)*epsilon(u)) - div(v)*p +\
         1.0E10*q*div(u) + delta*inner(nabla_grad(q), nabla_grad(p)))*dx
    L = inner(v + delta*nabla_grad(q), f)*dx - inner(waterpressure, v + delta*nabla_grad(q))*ds(5)

    # Picard iterations
    w = Function(system)  # new unknown function
    eps = 1.0             # error measure ||u-u_k|| and ||p-p_k||
    tol = 1.0E-6          # tolerance
    count = 0             # iteration counter
    max = 30              # max number of iterations allowed

    while eps>tol and count<max:
        count += 1
        solve(a == L, w, bcs)
        u,p = w.split(deepcopy=True)
        ux,uz = u.split(deepcopy=True)
        diffx = ux.vector().array() - ux_k.vector().array()
        diffz = uz.vector().array() - uz_k.vector().array()
        diffp = p.vector().array() - p_k.vector().array()
        epsx = np.linalg.norm(diffx)/np.linalg.norm(ux.vector().array())
        epsz = np.linalg.norm(diffz)/np.linalg.norm(uz.vector().array())
        epsp = np.linalg.norm(diffp)/np.linalg.norm(p.vector().array())
        if epsx > epsz:
            eps = epsx
        else:
            eps = epsz
        if eps < epsp:
            eps = epsp
#        print "count = %d, error = %g" % (count,eps)
        assign(ux_k,ux)  # update for next iteration
        assign(uz_k,uz)
        assign(u_k,u)
        assign(p_k,p)

    convergence = "convergence after %d Picard iterations" % count
    if count >= max:
        convergence = "no" + convergence

    # Plot final solution
    u,p = w.split()
    ux,uz = u.split(deepcopy=True)
    tensor = TensorFunctionSpace(mesh, "Lagrange", degree)

    # Full stress
    sigma = project(nu(u,temp)*epsilon(u)-p*Identity(tensor.cell().topological_dimension()),tensor)
    sigmaxx = Function(scalar)
    sigmaxz = Function(scalar)
    sigmazz = Function(scalar)
    assign(sigmaxx,sigma.sub(0))
    assign(sigmaxz,sigma.sub(1))
    assign(sigmazz,sigma.sub(3))
    
    # Deviatoric stress
    #tau = project(nu(u,temp)*epsilon(u),tensor)
    #tauxx = Function(scalar)
    #tauzz = Function(scalar)
    #assign(tauxx,tau.sub(0))
    #assign(tauzz,tau.sub(3))

    # Find the path for each basal
    for j in range(len(basal[0])):
        while True:
            k = len(bindex[j])
            eigenvalue,eigenvector = np.linalg.eig(np.array([[sigmaxx.vector().array()[bindex[j][k-1]],sigmaxz.vector().array()[bindex[j][k-1]]],[sigmaxz.vector().array()[bindex[j][k-1]],sigmazz.vector().array()[bindex[j][k-1]]]]))
            if z[bindex[j][k-1]] < depth:
                if eigenvalue[0] > eigenvalue[1]:
                    if eigenvalue[0] + rho_w*g*(depth-z[bindex[j][k-1]]) > 0:
                        propagation = eigenvector[:,1]
                    else:
                        break
                else:
                    if eigenvalue[1] + rho_w*g*(depth-z[bindex[j][k-1]]) > 0:
                        propagation = eigenvector[:,0]
                    else:
                        break
            else:
                if eigenvalue[0] > eigenvalue[1]:
                    if eigenvalue[0] > 0:
                        propagation = eigenvector[:,1]
                    else:
                        break
                else:
                    if eigenvalue[1] > 0:
                        propagation = eigenvector[:,0]
                    else:
                        break

            if propagation[1] < 0:
                propagation *= -1
        
            propagation *= gridsize
            bindex[j].append(closest_vertex([x[bindex[j][k-1]],z[bindex[j][k-1]]]+propagation))

    # Find the path for each surface
    for j in range(len(surface[0])):
        while True:
            k = len(sindex[j])
            eigenvalue,eigenvector = np.linalg.eig(np.array([[sigmaxx.vector().array()[sindex[j][k-1]],sigmaxz.vector().array()[sindex[j][k-1]]],[sigmaxz.vector().array()[sindex[j][k-1]],sigmazz.vector().array()[sindex[j][k-1]]]]))
            if eigenvalue[0] > eigenvalue[1]:
                if eigenvalue[0] > 0:
                    propagation = eigenvector[:,1]
                else:
                    break
            else:
                if eigenvalue[1] > 0:
                    propagation = eigenvector[:,0]
                else:
                    break
        
            if propagation[1] > 0:
                propagation *= -1
            
            propagation *= gridsize
            sindex[j].append(closest_vertex([x[sindex[j][k-1]],z[sindex[j][k-1]]]+propagation))

    # Save crevasse paths
    for i in range(len(basal[0])):
        for l in range(len(bindex[i])):
            if l < len(bverts_x[i]):
                bverts_x[i][l] = x[bindex[i][l]]
                bverts_z[i][l] = z[bindex[i][l]]
            else:
                bverts_x[i].append(x[bindex[i][l]])
                bverts_z[i].append(z[bindex[i][l]])

    for i in range(len(surface[0])):
        for l in range(len(sindex[i])):
            if l < len(sverts_x[i]):
                sverts_x[i][l] = x[sindex[i][l]]
                sverts_z[i][l] = z[sindex[i][l]]
            else:
                sverts_x[i].append(x[sindex[i][l]])
                sverts_z[i].append(z[sindex[i][l]])

    b0_x = np.array(bverts_x[0])
    b0_z = np.array(bverts_z[0])
    b1_x = np.array(bverts_x[1])
    b1_z = np.array(bverts_z[1])
    s0_x = np.array(sverts_x[0])
    s0_z = np.array(sverts_z[0])
    s1_x = np.array(sverts_x[1])
    s1_z = np.array(sverts_z[1])
    filename = "crevs"+str(day+730)+".out"
    
    with open(filename,"w") as f:
        f.write("\n".join(" ".join(map(str, x)) for x in (b0_x,b0_z,b1_x,b1_z,s0_x,s0_z,s1_x,s1_z)))

    # Stresses
    ddelta = project(sqrt(pow((sigmaxx-sigmazz),2) + 4*pow(sigmaxz,2)), scalar)
    sigma1 = project(0.5*(sigmaxx + sigmazz + ddelta), scalar)
    tau_max = project(sqrt(pow(0.5*(sigmaxx-sigmazz),2) + pow(sigmaxz,2)), scalar)
    #open_x = project(sqrt((ddelta-sigmaxx+sigmazz)/2.0/ddelta), scalar)
    #open_vector_z = project(sqrt(0.5+0.5*(sigmaxx-sigmazz)/ddelta), scalar)
    #open_vector = Function(vector)
    #assign(open_vector.sub(0), open_vector_x)
    #assign(open_vector.sub(1), open_vector_z)

    # Save mesh and solution to file
    File("mesh"+str(day+730)+".xml") << mesh
    File("tensile"+str(day+730)+".pvd") << sigma1
    File("shear"+str(day+730)+".pvd") << tau_max

    dt = 0.25         # 0.5 day per time step
    u1 = ux.compute_vertex_values()*dt
    u2 = uz.compute_vertex_values()*dt
    du = Function(vector)
    dux = Function(scalar)
    duz = Function(scalar)
    dux.vector()[:] = u1
    duz.vector()[:] = u2
    assign(du.sub(0), dux)
    assign(du.sub(1), duz)
    mesh.move(du)

    scalar = FunctionSpace(mesh, "CG", degree)
    vector = VectorFunctionSpace(mesh, "CG", degree)
    system = vector * scalar

#    if (day+1)%50 == 0:
#        print day
#        mesh = refine(mesh)
#        plt.figure()
#        plot(mesh)
#        plt.savefig('mesh'+str(day+1)+'.png')
#
    x = mesh.coordinates()[:,0]
    z = mesh.coordinates()[:,1]

    # Create mesh function over cell facets
    boundary_parts = FacetFunction('size_t', mesh, 0)
    boundary_parts.set_all(0)
    
    # DOLFIN_EPS = 3e-16
    # Mark right below water boundary as subdomain 5
    right_below = [f for f in facets(mesh)
                   if (f.exterior() and (f.midpoint().y()<depth))]
    for rb in right_below:
        boundary_parts[rb] = 5

    # Mark bottom boundary facets as subdomain 1
    bottom = [f for f in facets(mesh)
              if (f.exterior() and (f.midpoint().y()<1E-3))]
    for b in bottom:
        boundary_parts[b] = 1

    # Mark left boundary as subdomain 3
    left = [f for f in facets(mesh)
            if (f.exterior() and (f.midpoint().x()<1E-3))]
    for l in left:
        boundary_parts[l] = 3

    # No-slip boundary condition at bottom
    bcb1 = DirichletBC(system.sub(0), Constant((0.0,0.0)), boundary_parts, 1)
    
    # Free-slip boundary condition at bottom
    bcb2 = DirichletBC(system.sub(0).sub(1), Constant(0.0), boundary_parts, 1)
    
    # Free slip boundary condition on left hand side
    bcl = DirichletBC(system.sub(0).sub(0), Constant(0.0), boundary_parts, 3)
    
    # Traction-free boundary condition on surface
    # do nothing
    
    # Collect Dirichlet boundary conditions
    bcs = [bcb2, bcl]

    # Define new measures associated with the interior domains and
    # exterior boundaries
    ds = Measure("ds")[boundary_parts]

