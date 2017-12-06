import geometry
import solver
import plotter
#import PIL
from math import pi
import time

#########################################
############# PROBLEM SETUP #############
#########################################
pitch = 1.26
fwidth = 0.7                                #fuel width/height
q_fuel = 10 / (4 * pi )                       #fuel source
q_mod = 0.0                                  #moderator source



#convergence tolerance
tol = 1e-3
max_iters = 50

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
sigma_mod = (2 * xs_absorption_h + xs_scatter_o + 2 * xs_scatter_h) * num_density_h2o * barns       #moderator macro XS
sigma_mod_scatter = (xs_scatter_o + 2 * xs_scatter_h) * num_density_h2o * barns


#########################################
############## SETUP SOLVE ##############
#########################################
#set material objects
fuel = geometry.Material('fuel', q_fuel, sigma_fuel, sigma_fuel)
moderator = geometry.Material('moderator', q_mod, sigma_mod, sigma_mod_scatter)


#Sn order
orders = [2, 4, 8]
spacings = [0.005]  #mesh spacing
results = []

manyspacings = True
manyorders = False


#create directory to store plots in
timestr = 'plots/' + time.strftime("%Y-%m-%d_%H-%M")
plotter.mkdir_p(timestr)
savepath = timestr

print "Created new directory for plots.."


if manyspacings:
    print "Iterating over spacings....\n\n"
    order = 4
    i=0
    sptol = 0.1
    rat_old = 0
    for spacing in spacings:
    #for order in orders:
        print "\nSolving SN, order %d, spacing %g" % (order, spacing)
        #setup mesh cells
        mesh = geometry.Geometry(pitch, spacing, fwidth, fuel, moderator)
        mesh.setMesh()

        cell_width = int(pitch / spacing)
        fuel_width = int(fwidth / spacing)

        plot_cells = [int((cell_width / 2.0) * cell_width + cell_width / 2.0),
                  int((cell_width / 2.0 + fuel_width / 2 - 1) * cell_width + cell_width / 2.0) + fuel_width / 2 - 1,
                  int(cell_width * cell_width - 1),
                  int((cell_width / 2.0) * cell_width + cell_width / 2.0) + fuel_width / 2 - 1,
                  int((cell_width / 2.0 + 1) * cell_width - 1)]
        plotter.plotMaterial(mesh, spacing, plot_cells, savepath)

        #give order, mesh to solver
        solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod, tol, mesh.fuel_area, mesh.mod_area)
        solve.solveSN(max_iters, plotter, mesh, savepath)

        solve.getAvgScalarFlux()
        results.append(solve.results)
        fluxchg = (solve.results[-1] - rat_old)/solve.results[-1]
        print "flux ratio change: %g" %(fluxchg)
        if fluxchg <= sptol:
            print "Spatial mesh is converged!"
            break
        rat_old = solve.results[-1]
        i+=1

if manyorders:
    print "Iterating over orders...\n\n"
    results = []
    for order in orders:
        print "\nSolving SN, order %d, spacing %g" % (order, spacing)
        # setup mesh cells
        mesh = geometry.Geometry(pitch, spacing, fwidth, fuel, moderator)
        mesh.setMesh()

        # give order, mesh to solver
        solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod, tol, mesh.fuel_area,
                          mesh.mod_area)
        solve.solveSN(max_iters)
        solve.getAvgScalarFlux()
        results.append(solve.results)
        i=0
        print "\n\n---------------------------------\nSOLVE RESULTS\n---------------------------------"
        for result in results:
            print "SN Order: %d\n" %(orders[i])
            it, ff, mf, af, rat = result

            print "Iterations to convergence: %d \nModerator flux: %g\nFuel flux: %g\nAvg flux: %g\nFlux ratio: %g"\
                  %(it, mf, ff, af, rat)
            print "---------------------------------\n"
            i+=1