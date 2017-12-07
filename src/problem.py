import geometry
import solver
import plotter
#import PIL
from math import pi
import time

#########################################
############ RESULTS STORAGE ############
#########################################
#create directory to store plots, results in
pathname = 'plots/' + time.strftime("%Y-%m-%d_%H-%M")
timestr = time.strftime("%Y-%m-%d_%H-%M")
plotter.mkdir_p(pathname)
savepath = pathname
resultsfile = pathname + '/' + timestr + '_results'

f = open('%s.txt' %resultsfile, 'w+')

print "Created new directory for plots.."

#########################################
############# PROBLEM SETUP #############
#########################################
pitch = 1.24
fwidth = 0.8                               #fuel width/height
q_fuel = 10 / (4 * pi )                       #fuel source
q_mod = 0.0                                  #moderator source


#convergence tolerance
tol = 1e-6
max_iters = 50

f.write("********PROBLEM SETUP********\n")
f.write("cell pitch \t %g\nfuel width \t %g\nfuel source \t %g\nmod source \t %g\n\n" %(pitch, fwidth, q_fuel, q_mod))
f.write("converge tol \t %g\n\n" %(tol))

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

f.write("fuel total xs \t %g\nmod total xs \t %g\nmod scatter xs \t %g\n\n*****************************\n\n" %(sigma_fuel, sigma_mod, sigma_mod_scatter))

#########################################
############## SETUP SOLVE ##############
#########################################
#set material objects
fuel = geometry.Material('fuel', q_fuel, sigma_fuel, sigma_fuel)
moderator = geometry.Material('moderator', q_mod, sigma_mod, sigma_mod_scatter)


#Sn order
orders = [2, 4, 8]
spacings = [0.02, 0.01, 0.005, 0.004]  #mesh spacing
results = []

manyspacings = True
manyorders = False
convergemesh = False

if manyspacings:
    print "Iterating over spacings....\n\n"
    f.write("Iterating over spacings...\n\n")
    order = 4
    i=0
    sptol = 0.1
    rat_old = 0
    for spacing in spacings:
    #for order in orders:
        print "\nSolving SN, order %d, spacing %g" % (order, spacing)
        f.write("\nSolving SN, order %d, spacing %g" % (order, spacing))

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
        solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod, tol, mesh.fuel_area, mesh.mod_area, timestr)
        solve.solveSN(max_iters, plotter, mesh, savepath)

        #solve.getAvgScalarFlux()
        results.append(solve.results)
        f.write("\nConverged in %d iterations! \nAvg fuel flux = %f \nAvg mod flux = %f \nAverage Flux  = %f \nFlux ratio = %f\n\n"
                %(solve.results[0], solve.results[1], solve.results[2], solve.results[3], solve.results[4]))
        fluxchg = (solve.results[-1] - rat_old) / solve.results[-1]
        print "flux ratio change: %g" % (fluxchg)
        f.write("flux ratio change: %g\n\n\n" % (fluxchg))

        if convergemesh:
            if fluxchg <= sptol:
                print "Spatial mesh is converged with spacing of %g" %(spacing)
                f.write("Spatial mesh is converged with spacing of %g" %(spacing))
                break
        rat_old = solve.results[-1]
        i+=1

if manyorders:
    print "Iterating over orders...\n\n"
    f.write("Iterating over orders...\n\n")
    results = []
    spacing = 0.02

    for order in orders:
        print "\nSolving SN, order %d, spacing %g" % (order, spacing)
        f.write("\nSolving SN, order %d, spacing %g" % (order, spacing))

        # setup mesh cells
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

        # give order, mesh to solver
        solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod, tol, mesh.fuel_area,
                          mesh.mod_area, timestr)
        solve.solveSN(max_iters, plotter, mesh, savepath)


        results.append(solve.results)
        f.write("\nConverged in %d iterations! \nAvg fuel flux = %f \nAvg mod flux = %f \nAverage Flux  = %f \nFlux ratio = %f\n\n"
            % (solve.results[0], solve.results[1], solve.results[2], solve.results[3], solve.results[4]))
