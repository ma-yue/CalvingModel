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
steps = 2*365             # steps, ~0.75 yr
n = 6.0                   # a constant
thickness = 800.0         # thickness of mesh
length = n*thickness      # length of mesh
depth = 700.0             # m ice under water
gridsize = 16.0           # m resolution
numofb = 2                # init. number of basal crevasses
numofs = 2                # init. number of surface crevasses

B0 = 15.77                # time-dependent viscosity constant, Pa ...
temp = 253.0              # K temperature
rho_i = 910               # kg/m3 ice density
rho_w = 1020              # kg/m3 seawater density
g = 9.8                   # m/s2
noblowup = 1E-14          # a small constant preventing inf viscosity values
meltrate = 0.3            # m/day avg, 0 at water line, 2 at the bed, linear gradient
grounding = length-2*gridsize     # init. grounding line position guess (needed for re-meshing)
cliff = length-2*gridsize         # init. cliff position guess (needed for surface crevasse placing)

# Create empty Mesh
mesh = Mesh()

# Create list of polygonal domain vertices
domain_vertices = [Point(0.0, 0.0),
                   Point(0.0, thickness),
                   Point(length, thickness),
                   Point(length, 0.0),
                   Point(0.0, 0.0)]

# Generate mesh
PolygonalMeshGenerator.generate(mesh, domain_vertices, gridsize)

# Obtain x,z coordinates of vertices
x = mesh.coordinates()[:,0]
z = mesh.coordinates()[:,1]

# Define function spaces
degree = 1
scalar = FunctionSpace(mesh, "CG", degree)
vector = VectorFunctionSpace(mesh, "CG", degree)
system = vector * scalar

# Create mesh function over cell facets
boundary_parts = FacetFunction('size_t', mesh, 0)
boundary_parts.set_all(0)

# DOLFIN_EPS = 3e-16
# Mark bottom boundary facets as subdomain 1
class bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[1]) < DOLFIN_EPS

gamma_1 = bottom()
gamma_1.mark(boundary_parts, 1)

# Mark top boundary facets as subdomain 2
class top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[1] - thickness) < DOLFIN_EPS

gamma_2 = top()
gamma_2.mark(boundary_parts, 2)

# Mark left boundary as subdomain 3
class left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0]) < DOLFIN_EPS

gamma_3 = left()
gamma_3.mark(boundary_parts, 3)

# Mark right above water boundary as subdomain 4
class right_above(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] - length) < DOLFIN_EPS and x[1] >= depth

gamma_4 = right_above()
gamma_4.mark(boundary_parts, 4)

# Mark right below water boundary as subdomain 5
class right_below(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] - length) < DOLFIN_EPS and x[1] <= depth

gamma_5 = right_below()
gamma_5.mark(boundary_parts, 5)


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
#dx = Measure("dx")[domains]
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
    d = 3.*gridsize
    for v in vertices(mesh):
        deltax = x[v.index()]-point[0]
        deltaz = z[v.index()]-point[1]
        if abs(deltax) <= d and abs(deltaz) <= d:
            distance = math.hypot(deltax,deltaz)
            if distance < d:
                d = distance
                id = v.index()
    return id

#start = timeit.default_timer()

# Create empty lists to store info for propagation paths
# basal
#bindex = []
#basal = np.zeros([2,numofb])
#bverts_x = []
#bverts_z = []
#for k in range(len(basal[0])):
#    bindex.append([])
#    bverts_x.append([])
#    bverts_z.append([])
# surface
#sindex = []
#surface = np.zeros([2,numofs])
#sverts_x = []
#sverts_z = []
#for k in range(len(surface[0])):
#    sindex.append([])
#    sverts_x.append([])
#    sverts_z.append([])

# Sampling the bottom and surface to initiate crevasses
#for j in range(len(basal[0])):
#    basal[:,j] = [length-2*gridsize-0.1*thickness*j,0]
#    bindex[j].append(closest_vertex(basal[:,j]))
#    bverts_x[j].append(x[bindex[j][0]])
#    bverts_z[j].append(z[bindex[j][0]])
#
#for j in range(len(surface[0])):
#    surface[:,j] = [length-2*gridsize-0.1*thickness*j,thickness]
#    sindex[j].append(closest_vertex(surface[:,j]))
#    sverts_x[j].append(x[sindex[j][0]])
#    sverts_z[j].append(z[sindex[j][0]])


for step in range(steps):
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
    max = 32              # max number of iterations allowed

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
#    for j in range(len(basal[0])):
#        while True:
#            k = len(bindex[j])
#            eigenvalue,eigenvector = np.linalg.eig(np.array([[sigmaxx.vector().array()[bindex[j][k-1]],sigmaxz.vector().array()[bindex[j][k-1]]],[sigmaxz.vector().array()[bindex[j][k-1]],sigmazz.vector().array()[bindex[j][k-1]]]]))
#            if z[bindex[j][k-1]] < depth:
#                if eigenvalue[0] > eigenvalue[1]:
#                    if eigenvalue[0] + rho_w*g*(depth-z[bindex[j][k-1]]) > 0:
#                        propagation = eigenvector[:,1]
#                    else:
#                        break
#                else:
#                    if eigenvalue[1] + rho_w*g*(depth-z[bindex[j][k-1]]) > 0:
#                        propagation = eigenvector[:,0]
#                    else:
#                        break
#            else:
#                if eigenvalue[0] > eigenvalue[1]:
#                    if eigenvalue[0] > 0:
#                        propagation = eigenvector[:,1]
#                    else:
#                        break
#                else:
#                    if eigenvalue[1] > 0:
#                        propagation = eigenvector[:,0]
#                    else:
#                        break
#
#            if propagation[1] < 0:
#                propagation *= -1
#
#            propagation *= gridsize
#            bindex[j].append(closest_vertex([x[bindex[j][k-1]],z[bindex[j][k-1]]]+propagation))

    # Find the path for each surface
#    for j in range(len(surface[0])):
#        while True:
#            k = len(sindex[j])
#            eigenvalue,eigenvector = np.linalg.eig(np.array([[sigmaxx.vector().array()[sindex[j][k-1]],sigmaxz.vector().array()[sindex[j][k-1]]],[sigmaxz.vector().array()[sindex[j][k-1]],sigmazz.vector().array()[sindex[j][k-1]]]]))
#            if z[sindex[j][k-1]] < depth:
#                if eigenvalue[0] > eigenvalue[1]:
#                    if eigenvalue[0] + rho_w*g*(depth-z[sindex[j][k-1]]) > 0:
#                        propagation = eigenvector[:,1]
#                    else:
#                        break
#                else:
#                    if eigenvalue[1] + rho_w*g*(depth-z[sindex[j][k-1]]) > 0:
#                        propagation = eigenvector[:,0]
#                    else:
#                        break
#            else:
#                if eigenvalue[0] > eigenvalue[1]:
#                    if eigenvalue[0] > 0:
#                        propagation = eigenvector[:,1]
#                    else:
#                        break
#                else:
#                    if eigenvalue[1] > 0:
#                        propagation = eigenvector[:,0]
#                    else:
#                        break
#            
#            if propagation[1] > 0:
#                propagation *= -1
#            
#            propagation *= gridsize
#            sindex[j].append(closest_vertex([x[sindex[j][k-1]],z[sindex[j][k-1]]]+propagation))

    # Save crevasse paths
#    for i in range(len(basal[0])):
#        for l in range(len(bindex[i])):
#            if l < len(bverts_x[i]):
#                bverts_x[i][l] = x[bindex[i][l]]
#                bverts_z[i][l] = z[bindex[i][l]]
#            else:
#                bverts_x[i].append(x[bindex[i][l]])
#                bverts_z[i].append(z[bindex[i][l]])
#
#    for i in range(len(surface[0])):
#        for l in range(len(sindex[i])):
#            if l < len(sverts_x[i]):
#                sverts_x[i][l] = x[sindex[i][l]]
#                sverts_z[i][l] = z[sindex[i][l]]
#            else:
#                sverts_x[i].append(x[sindex[i][l]])
#                sverts_z[i].append(z[sindex[i][l]])
#
#    b0_x = np.array(bverts_x[0])
#    b0_z = np.array(bverts_z[0])
#    b1_x = np.array(bverts_x[1])
#    b1_z = np.array(bverts_z[1])
#    s0_x = np.array(sverts_x[0])
#    s0_z = np.array(sverts_z[0])
#    s1_x = np.array(sverts_x[1])
#    s1_z = np.array(sverts_z[1])
#    if numofs == 3:
#        s2_x = np.array(sverts_x[2])
#        s2_z = np.array(sverts_z[2])
#    elif numofs == 4:
#        s2_x = np.array(sverts_x[2])
#        s2_z = np.array(sverts_z[2])
#        s3_x = np.array(sverts_x[3])
#        s3_z = np.array(sverts_z[3])
#    elif numofs == 5:
#        s2_x = np.array(sverts_x[2])
#        s2_z = np.array(sverts_z[2])
#        s3_x = np.array(sverts_x[3])
#        s3_z = np.array(sverts_z[3])
#        s4_x = np.array(sverts_x[4])
#        s4_z = np.array(sverts_z[4])
#    elif numofs == 6:
#        s2_x = np.array(sverts_x[2])
#        s2_z = np.array(sverts_z[2])
#        s3_x = np.array(sverts_x[3])
#        s3_z = np.array(sverts_z[3])
#        s4_x = np.array(sverts_x[4])
#        s4_z = np.array(sverts_z[4])
#        s5_x = np.array(sverts_x[5])
#        s5_z = np.array(sverts_z[5])
#
#    filename = "crevs"+str(step+1)+".out"
#
#    if numofs == 3:
#        with open(filename,"w") as f:
#            f.write("\n".join(" ".join(map(str, x)) for x in (b0_x,b0_z,b1_x,b1_z,s0_x,s0_z,s1_x,s1_z,s2_x,s2_z)))
#    elif numofs == 4:
#        with open(filename,"w") as f:
#            f.write("\n".join(" ".join(map(str, x)) for x in (b0_x,b0_z,b1_x,b1_z,s0_x,s0_z,s1_x,s1_z,s2_x,s2_z,s3_x,s3_z)))
#    elif numofs == 5:
#        with open(filename,"w") as f:
#            f.write("\n".join(" ".join(map(str, x)) for x in (b0_x,b0_z,b1_x,b1_z,s0_x,s0_z,s1_x,s1_z,s2_x,s2_z,s3_x,s3_z,s4_x,s4_z)))
#    elif numofs == 6:
#        with open(filename,"w") as f:
#            f.write("\n".join(" ".join(map(str, x)) for x in (b0_x,b0_z,b1_x,b1_z,s0_x,s0_z,s1_x,s1_z,s2_x,s2_z,s3_x,s3_z,s4_x,s4_z,s5_x,s5_z)))
#    else:
#        with open(filename,"w") as f:
#            f.write("\n".join(" ".join(map(str, x)) for x in (b0_x,b0_z,b1_x,b1_z,s0_x,s0_z,s1_x,s1_z)))

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
    File("mesh"+str(step+1)+".xml") << mesh
    File("tensile"+str(step+1)+".pvd") << sigma1
    File("shear"+str(step+1)+".pvd") << tau_max
    
    dt = 0.25         # 0.25 day per time step
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
    
    # Obtain x,z coordinates of vertices
    x = mesh.coordinates()[:,0]
    z = mesh.coordinates()[:,1]

    # Mark boundary vertices
    bdry_label = []
    bdry_v = []
    bdry_order = []   # bdry_order contains the boundary vertices, in order
    for v in vertices(mesh):
        for f in facets(v):
            if f.exterior():
                bdry_v.append(v)
                bdry_label.append(v.index())
                break

    # Ordering the boundary vertices
    bdry_no = len(bdry_v)
    head = bdry_v[0]
    bdry_order.append(head)                  # keep in mind that here head is a fenics object so head.index()
    count = 1                                # is different from the bdry_label.index(...) below
    while count < bdry_no:
        flag = 0
        id = bdry_label.index(head.index())   # returns the lowest index in list that head.index() appears
        del bdry_label[id]
        del bdry_v[id]
        for f in facets(head):
            if f.exterior():
                for v in vertices(f):
                    if v in bdry_v:
                        bdry_order.append(v)
                        head = v
                        flag = 1
                        count += 1
                        break
            if flag:
                break

    domain_vertices = []           # a list to store final boundary points, in order
    # identify the gounding line position and surface "cliff" position
    for v in bdry_order:
        k = v.index()
        if x[k] > grounding:
            if z[k] < 0.1:
                grounding = x[k]
                z[k] = 0
                pointer = v  # point along the calving front to be melted, initially at grounding line
            elif z[k] > depth:
                if x[k] > cliff:
                    cliff = x[k]
        if z[k] < 0:
            z[k] = 0

    # Deleting all the nodes on the bed except (0,0) and the grounding line
    for v in bdry_order:
        k = v.index()
        if x[k] > 0.1 and z[k] < 0.1 and x[k] < grounding:
            del v
            continue

    print grounding

    # Pick out the boundary vertices under water and melt according to z position
    x[pointer.index()] -= 2*meltrate*dt  # melt away the grounding line first
    end = 0   # mark if the end of all vertices on calving front is found
    while True:
        next = 0  # mark if the next vertex on calving front is found
        for f in facets(pointer):
            if f.exterior():
                for v in vertices(f):
                    k = v.index()
                    if z[k] > depth:
                        end = 1
                        break
                    elif z[k] > z[pointer.index()]:
                        pointer = v
                        x[k] -= (1-z[k]/depth)*2*meltrate*dt
                        next = 1
                        break
            if next or end:
                break
        if end:
            break

    # The new outline of the mesh, in order
    for v in bdry_order:
        k = v.index()
        # getting rid of vertices too close to the grounding line
        if z[k] < depth and x[k] > 0.5*grounding:
            d = math.hypot(x[k]+(1-z[k]/depth)*2*meltrate*dt-grounding,z[k])
            if d > 1E-3 and d < gridsize:
                continue
        domain_vertices.append(Point(x[k],z[k]))
    
    # Reset guess of grounding line position
    grounding -= 2*meltrate*dt

    # Generate mesh
    PolygonalMeshGenerator.generate(mesh, domain_vertices, gridsize)
    
    # Obtain x,z coordinates of vertices
    x = mesh.coordinates()[:,0]
    z = mesh.coordinates()[:,1]
    
    # Define function spaces
    scalar = FunctionSpace(mesh, "CG", degree)
    vector = VectorFunctionSpace(mesh, "CG", degree)
    system = vector * scalar

#    if (day+1)%50 == 0:
#        print day
#        mesh = refine(mesh)
#        plt.figure()
#        plot(mesh)
#        plt.savefig('mesh'+str(day+1)+'.png')

    # Create mesh function over cell facets
    boundary_parts = FacetFunction('size_t', mesh, 0)
    boundary_parts.set_all(0)
    
    # DOLFIN_EPS = 3e-16
    # Note that when no boundary was labeled, facet.exterior() does not return meaningful values
    # Mark right below water boundary as subdomain 5
    right_below = []
    for f in facets(mesh):
        count = 0
        for c in cells(f):
            count += 1
        if count == 1 and f.midpoint().y()<depth:
            right_below.append(f)
    for rb in right_below:
        boundary_parts[rb] = 5
    
    # Mark bottom boundary facets as subdomain 1
    bottom = [f for f in facets(mesh) if f.midpoint().y()<1E-3]
    for b in bottom:
        boundary_parts[b] = 1

    # Mark left boundary as subdomain 3
    left = [f for f in facets(mesh) if f.midpoint().x()<1E-3]
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

    # Determine number of surface crevasses based on cliff position
#    overhang = cliff - grounding
#    if overhang > 0 and overhang < gridsize:
#        numofs = 3
#    elif overhang >= gridsize and overhang < 2*gridsize:
#        numofs = 4
#    elif overhang >= 2*gridsize and overhang < 3*gridsize:
#        numofs = 5
#    else:
#        numofs = 6

    # Create empty lists to store info for propagation paths
    # basal
#    bindex = []
#    basal = np.zeros([2,numofb])
#    bverts_x = []
#    bverts_z = []
#    for k in range(len(basal[0])):
#        bindex.append([])
#        bverts_x.append([])
#        bverts_z.append([])
    # surface
#    sindex = []
#    surface = np.zeros([2,numofs])
#    sverts_x = []
#    sverts_z = []
#    for k in range(len(surface[0])):
#        sindex.append([])
#        sverts_x.append([])
#        sverts_z.append([])

    # Sampling the bottom and surface to initiate crevasses
    # Be careful with the use of grounding here, for both surface and basal crevasses
#    for j in range(len(basal[0])):
#        basal[:,j] = [grounding-2*gridsize-0.1*thickness*j,0]
#        bindex[j].append(closest_vertex(basal[:,j]))
#        bverts_x[j].append(x[bindex[j][0]])
#        bverts_z[j].append(z[bindex[j][0]])
#    
#    for j in range(len(surface[0])):
#        surface[:,j] = [grounding-2*gridsize-0.1*thickness*j,thickness/2]
#        for v in vertices(mesh):
#            if abs(v.x(0)-surface[:,j][0]) < 0.5*gridsize:
#                if v.x(1) > surface[:,j][1]:
#                    surface[:,j][1] = v.x(1)
#        sindex[j].append(closest_vertex(surface[:,j]))
#        sverts_x[j].append(x[sindex[j][0]])
#        sverts_z[j].append(z[sindex[j][0]])

    # Place the extra surface crevasse in the overhang
#    if numofs > 2:
#        for i in range(numofs-2):
#            surface[:,2+i] = [grounding+i*gridsize,depth+3*gridsize]
#            for v in vertices(mesh):
#                if abs(v.x(0)-surface[:,2+i][0]) < 0.5*gridsize:
#                    if v.x(1) > surface[:,2+i][1]:
#                        surface[:,2+i][1] = v.x(1)
#            sindex[2+i][0] = closest_vertex(surface[:,2+i])
#            sverts_x[2+i][0] = x[sindex[2+i][0]]
#            sverts_z[2+i][0] = z[sindex[2+i][0]]

#os.system('say "your program has finished"')
