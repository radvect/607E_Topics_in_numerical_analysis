import numpy as np
from dolfinx.fem import assemble_scalar, form
import ufl
from petsc4py import PETSc
from mpi4py import MPI
import gmsh
from dolfinx import fem, mesh, io
from dolfinx.fem.petsc import assemble_vector, assemble_matrix, create_vector, apply_lifting, set_bc
from scipy.special import jv, jn_zeros
gmsh.initialize() 
gmsh.option.setNumber("General.Verbosity", 5)

def interpolate(u1, u2):

    rmesh = u2.function_space.mesh
    rtopology = rmesh.topology
    cmap = rtopology.index_map(rtopology.dim)
    num_cells = cmap.size_local + cmap.num_ghosts
    all_cells = np.arange(num_cells, dtype=np.int32)
    nmmid_all = fem.create_interpolation_data(
        u2.function_space,
        u1.function_space,
        all_cells,
        padding=1e-14,
    )
    u2.interpolate_nonmatching(u1, cells=all_cells, interpolation_data=nmmid_all)

def create_area(element_size):
    
    gmsh.option.setNumber("General.Terminal", 1)
    model = gmsh.model
    model.add("circle")
    circle = model.occ.addCircle(0, 0, 0, 1)
    loop = model.occ.addCurveLoop([circle])
    surface = model.occ.addPlaneSurface([loop])
    model.occ.synchronize()
    
  
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", element_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", element_size)
    
    model.addPhysicalGroup(2, [surface], tag=1, name="Interior")
    model.addPhysicalGroup(1, [circle], tag=2, name="Boundary")
    model.mesh.generate(2)
    

    return model

def initial_condition_line(x):
    return np.cos(2 * np.pi * x[0])  
def create_area_line(element_size):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    model = gmsh.model()
    model.add("line")

    p1 = model.occ.addPoint(0, 0, 0)
    p2 = model.occ.addPoint(1, 0, 0)
    line = model.occ.addLine(p1, p2)

    model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", element_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", element_size)
    
    model.add_physical_group(dim=1, tags=[line], name="Boundary")
    
    gmsh.model.mesh.generate(dim=1)

    return model

def create_area_square(element_size):

    gmsh.option.setNumber("General.Terminal", 1)  

    model = gmsh.model
    model.add("square")  

  
    p1 = model.occ.addPoint(0, 0, 0, element_size)
    p2 = model.occ.addPoint(1, 0, 0, element_size)
    p3 = model.occ.addPoint(1, 1, 0, element_size)
    p4 = model.occ.addPoint(0, 1, 0, element_size)

 
    l1 = model.occ.addLine(p1, p2)
    l2 = model.occ.addLine(p2, p3)
    l3 = model.occ.addLine(p3, p4)
    l4 = model.occ.addLine(p4, p1)

 
    cl = model.occ.addCurveLoop([l1, l2, l3, l4])

    
    s = model.occ.addPlaneSurface([cl])


    model.occ.synchronize()


 
    boundary = model.addPhysicalGroup(1, [l1, l2, l3, l4], tag=1, name="Boundary")

  
    interior = model.addPhysicalGroup(2, [s], tag=2, name="Interior")


    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", element_size)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", element_size)

 
    model.mesh.generate(2)

    return model

def initial_condition_square(x):
    return np.cos(2 * np.pi * x[0])*np.cos(2 * np.pi * x[1])


def L2_error(uh, u_ref):
    print("shape uh beginning",uh.x.array.shape)
    Vref = u_ref.function_space
    uh_Vref = fem.Function(Vref)
    interpolate(uh, uh_Vref)

    e_Vref = fem.Function(Vref)

    print("shape_ref", uh_Vref.x.array.shape)
    print("shape_num", u_ref.x.array.shape)
    
    e_Vref.x.array[:] = uh_Vref.x.array - u_ref.x.array 
    
    error = form(e_Vref**2 * ufl.dx)
    E = np.sqrt(assemble_scalar(error))
    return E

def solve_diffusion(element_size, num_steps=10, T=0.1):
    t = 0  
    dt = T / num_steps  

    model = create_area_square(element_size)
    partitioner = mesh.create_cell_partitioner(mesh.GhostMode.shared_facet)
    msh, _, _ = io.gmshio.model_to_mesh(
        model, MPI.COMM_WORLD, 0, gdim=2, partitioner=partitioner)
    domain = msh

    V = fem.functionspace(domain, ("Lagrange", 1))
    
    def initial_condition(x):
        first_dirichlet_zero = jn_zeros(1, 1)[0]  # j_0
        r = np.sqrt(x[0]**2 + x[1]**2)
        theta_cos = np.zeros_like(r)
        mask = r > 1e-14
        theta_cos[mask] = x[0][mask] / r[mask]
        scaled_r = r * first_dirichlet_zero
        return 3000*jv(0, scaled_r) * theta_cos


    
    u_n = fem.Function(V)
    u_n.name = "u_n"
    u_n.interpolate(initial_condition_square)

    # fdim = domain.topology.dim - 1
    # boundary_facets = mesh.locate_entities_boundary(
    #     domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool))
    # bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)

    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    f = fem.Constant(domain, PETSc.ScalarType(0))
    
    C= 1
    k=1
    a = C*u * v * ufl.dx + k*dt * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = (u_n - dt * f) * v * ufl.dx 

    bilinear_form = fem.form(a)
    linear_form = fem.form(L)

    A = assemble_matrix(bilinear_form)  #  bcs=[bc]
    A.assemble()
    b = create_vector(linear_form)

    solver = PETSc.KSP().create(domain.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU)

    uh = fem.Function(V)
    uh.name = "uh"
    with io.XDMFFile(domain.comm, "test_mesh.xdmf", "w") as mesh_file:
        mesh_file.write_mesh(domain)
    xdmf = io.XDMFFile(domain.comm, f"diffusion_circle_1234.xdmf", "w")
    xdmf.write_mesh(domain)
    xdmf.write_function(u_n, t)

    for i in range(num_steps):
        t += dt

        with b.localForm() as loc_b:
            loc_b.set(0)
        assemble_vector(b, linear_form)


        # apply_lifting(b, [bilinear_form], [[bc]])
        b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
        # set_bc(b, [bc])

      
        solver.solve(b, uh.x.petsc_vec)
        uh.x.scatter_forward()

   
        u_n.x.array[:] = uh.x.array

     
        xdmf.write_function(uh, t)

    xdmf.close()
    return uh


reference_solution = solve_diffusion(element_size=0.001)
reference_values = reference_solution.x.array
print(reference_values)
print(len(reference_values))

errors = []


mesh_sizes = sorted([0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125])

for h in mesh_sizes:
    numerical_solution = solve_diffusion(element_size=h)
    numerical_values = numerical_solution.x.array

    error = L2_error(numerical_solution, reference_solution)
    errors.append(error)

h_values = np.array(mesh_sizes)

import matplotlib.pyplot as plt

plt.figure()
plt.loglog(h_values, errors, marker='o', label='L2 error')


log_h = np.log(h_values)
log_errors = np.log(errors)


from scipy.stats import linregress


slope, intercept, _, _, _ = linregress(log_h, log_errors)
approx_line = np.exp(intercept) * h_values**slope
plt.loglog(h_values, approx_line, linestyle='--', label=f'Approx: O(h^{slope:.2f})')

plt.xlabel('h')
plt.ylabel('L2 error')
plt.title('Convergence rate')
plt.legend()
plt.grid(True, which='both', ls='--')
plt.savefig('convergence_rate.png', dpi=300)
plt.show()

