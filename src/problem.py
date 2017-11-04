import geometry
import solver

pitch = 1
spacing = 0.1       #mesh spacing
fwidth = 0.4        #fuel width/height
q_fuel = 1          #fuel source
q_mod = 0           #moderator source
sigma_fuel = 1      #fuel macro XS
sigma_mod = 0       #moderator macro XS

#Sn order
order = 4

#set material objects
fuel = geometry.Material('fuel', q_fuel, sigma_fuel)
moderator = geometry.Material('moderator', q_mod, sigma_mod)

#setup mesh cells
mesh = geometry.Geometry(pitch, spacing, fwidth, fuel, moderator)
mesh.setMesh()

#give order, mesh to solver
solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod)
solve.solveSN()
solve.getAvgScalarFlux()



