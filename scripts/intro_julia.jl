using Ferrite, SparseArrays


###############################
######## HEAT EQUATION ########
###############################

# Generate a mesh grid with 20 x 20 elements
grid = generate_grid(Quadrilateral, (20,20))

dim = 2
ip = Lagrange{dim, RefCube, 1}()
qr = QuadratureRule{dim, RefCube}(2)
cellvalues = CellScalarValues(qr, ip);

# Dofhandler - Manages distribution and numbering of degrees of freedom
# We can use it to query information about dofs, i.e. what dofs belong to a certain cell
# This is something we need when we are assembling the global stiffness matrix
dh = DofHandler(grid);

# Before distributing dofs, we need to add our Fields
# Here, we add a displacement field, u
add!(dh, :u, 1) 

# We distribute dofs with dh (dofhandler based on grid), with the variable u
# Each dof has 2 degrees of freedom (2dimensional problem)
close!(dh)

# Now we create a stiffness matrix for assembly
K = create_sparsity_pattern(dh);

# Next, dirichlet boundary conditions are dealt with by a "ConstraintHandler"
ch = ConstraintHandler(dh);

# define boundary for boundary conditions
∂Ω = union(
    getfaceset(grid, "left"),
    getfaceset(grid, "right"),
    getfaceset(grid, "top"),
    getfaceset(grid, "bottom"),
)

# Now let us define the boundary condition
# Here, we first define the field the bc is applied to,
# then, we define the boundary ∂Ω at which the BC is applied
# and lastly, the type of function. We say that the function is dependant on
# x (coordinates) and t (time), and whatever these are, we return a value 0
# i.e. the displacement value at the boundary is always 0

dbc = Dirichlet(:u, ∂Ω, (x,t) -> 0)
add!(ch, dbc);
close!(ch)

# assembly
# we define a function which takes the following as input: Ke, fe (pre-allocated empty matrix/vector) and cellvalues
function assemble_element!(Ke::Matrix, fe::Vector, cellvalues::CellScalarValues)
    n_basefuncs = getnbasefunctions(cellvalues)
    # Reset to 0
    fill!(Ke, 0)
    fill!(fe, 0)
    # Loop over quadrature points
    for q_point in 1:getnquadpoints(cellvalues)
        # Get the quadrature weight
        dΩ = getdetJdV(cellvalues, q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            δu  = shape_value(cellvalues, q_point, i)
            ∇δu = shape_gradient(cellvalues, q_point, i)
            # Add contribution to fe
            fe[i] += δu * dΩ
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = shape_gradient(cellvalues, q_point, j)
                # Add contribution to Ke
                Ke[i, j] += (∇δu ⋅ ∇u) * dΩ
            end
        end
    end
    return Ke, fe
end

function assemble_global(cellvalues::CellScalarValues, K::SparseMatrixCSC, dh::DofHandler)
    # Allocate the element stiffness matrix and element force vector
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    fe = zeros(n_basefuncs)
    # Allocate global force vector f
    f = zeros(ndofs(dh))
    # Create an assembler
    assembler = start_assemble(K, f)
    # Loop over all cels
    for cell in CellIterator(dh)
        # Reinitialize cellvalues for this cell
        reinit!(cellvalues, cell)
        # Compute element contribution
        assemble_element!(Ke, fe, cellvalues)
        # Assemble Ke and fe into K and f
        assemble!(assembler, celldofs(cell), Ke, fe)
    end
    return K, f
end

# We can now assemble our stiffness matrix and force vector
K, f = assemble_global(cellvalues, K, dh);

# We account for boundary conditions by using the apply! function
apply!(K, f, ch)

# Back-slash operator to solve system of equations
u = K \ f;

# Export to vtk
vtk_grid("heat_equation", dh) do vtk
    vtk_point_data(vtk, dh, u)
end
