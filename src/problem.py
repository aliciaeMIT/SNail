import geometry
import solver
from math import pi

#########################################
############# PROBLEM SETUP #############
#########################################
pitch = 1.25
spacing = 0.05                              #mesh spacing
fwidth = 0.8                                #fuel width/height
q_fuel = 10 / 4 * pi                        #fuel source
q_mod = 0                                   #moderator source

#Sn order
order = 4

#convergence tolerance
tol = 1e-5
max_iters = 100

#########################################
########## MATERIAL PROPERTIES ##########
#########################################
density_uo2 = 10.4                          #g/cc
density_h2o = 0.7                           #g/cc
molwt_u238 = 238                            #g/mol
molwt_o = 16                                #g/mol
molwt_h2 = 2                                #g/mol

molwt_uo2 = molwt_u238 + 2 * molwt_o
molwt_h2o = molwt_h2 + molwt_o

N_A = 6.022e23                              #1/mol; Avogadro's number

num_density_uo2 = density_uo2 * (N_A / molwt_uo2) #at/cc
num_density_h2o = density_h2o * (N_A / molwt_h2o) #at/cc

#microscopic XS, in barns
xs_scatter_238 = 11.29
xs_scatter_o = 3.888
xs_scatter_h = 20.47
xs_absorption_h = 1.0

barns = 1e-24 #cm; conversion factor

sigma_fuel = (xs_scatter_238 + 2 * xs_scatter_o) * num_density_uo2 * barns                             #fuel macro XS
sigma_mod = (2 * xs_absorption_h + xs_scatter_o + 2 * xs_absorption_h) * num_density_h2o * barns       #moderator macro XS


#########################################
############## SETUP SOLVE ##############
#########################################
#set material objects
fuel = geometry.Material('fuel', q_fuel, sigma_fuel)
moderator = geometry.Material('moderator', q_mod, sigma_mod)

#setup mesh cells
mesh = geometry.Geometry(pitch, spacing, fwidth, fuel, moderator)
mesh.setMesh()

#give order, mesh to solver
solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod, tol, mesh.fuel_area, mesh.mod_area)
solve.solveSN(max_iters)
solve.getAvgScalarFlux()



