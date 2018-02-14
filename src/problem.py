import geometry
import solver
import plotter
from math import pi
import time

#########################################
############# PROBLEM SETUP #############
#########################################
pitch = 1.6
fwidth = 0.8                               #fuel width/height
q_fuel = 10 / (4 * pi )                    #fuel source
q_mod = 0.0                                #moderator source

#convergence tolerance
tol = 1e-6
max_iters = 200

#Sn order, mesh spacings
orders = [4]
spacings = [0.2, 0.1, 0.08, 0.05, 0.04, 0.02]

getnumunknowns = False
singlesolve = False
manyspacings = True
manyorders = False
convergemesh = False
update_source = True
tally_fuel_corner = True

if singlesolve:
    order = 4
    spacing = 0.2

if manyspacings:
    order = 4
    if convergemesh:
        sptol = 0.1  # percent change in flux ratio

elif manyorders:
    spacing = 0.005

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

#macroscopic XS, cm^-1
sigma_fuel = (xs_scatter_238 + 2 * xs_scatter_o) * num_density_uo2 * barns
sigma_mod = (2 * xs_absorption_h + xs_scatter_o + 2 * xs_scatter_h) * num_density_h2o * barns
sigma_mod_scatter = (xs_scatter_o + 2 * xs_scatter_h) * num_density_h2o * barns

#set material objects
fuel = geometry.Material('fuel', q_fuel, sigma_fuel, sigma_fuel)
moderator = geometry.Material('moderator', q_mod, sigma_mod, sigma_mod_scatter)


#########################################
############ RESULTS STORAGE ############
#########################################
#create directory to store plots, results in
timestr = time.strftime("%Y-%m-%d_%H-%M")
pathname = 'plots/' + timestr
plotter.mkdir_p(pathname)
savepath = pathname
resultsfile = pathname + '/' + timestr + '_results'

f = open('%s.txt' %resultsfile, 'w+')

f.write("********PROBLEM SETUP********\n")
f.write("cell pitch \t %g\nfuel width \t %g\nfuel source \t %g\nmod source \t %g\n\n" %(pitch, fwidth, q_fuel, q_mod))
f.write("converge tol \t %g\n\n" %(tol))
f.write("fuel total xs \t %g\nmod total xs \t %g\nmod scatter xs \t %g\n\n"
        "*****************************\n\n" %(sigma_fuel, sigma_mod, sigma_mod_scatter))
f.close()
#########################################
################  SOLVE  ################
#########################################

def solveSN(spacing, order, savepath):
    f = open('%s.txt' %resultsfile, 'a+')
    print "\nSolving SN, order %d, spacing %g" % (order, spacing)
    f.write("\nSolving SN, order %d, spacing %g" % (order, spacing))
    f.write("\nTotal number of angles: \t%g" % (order * (order + 2) / 2))
    # setup mesh cells
    mesh = geometry.Geometry(pitch, spacing, fwidth, fuel, moderator)
    mesh.setMesh(tally_fuel_corner)

    cell_width = mesh.getWidth(pitch, spacing)
    fuel_width = mesh.getWidth(fwidth, spacing)
    plot_cells = mesh.getPlotCells(cell_width, fuel_width)
    plotter.plotMaterial(mesh, spacing, plot_cells, savepath)
    f.write("\n#ang\tcells\t#unkn\n%g\t%g\t%g\n" % ( order * (order + 2) / 2 , mesh.n_cells ** 2, (mesh.n_cells ** 2) * (order*(order+2)/2)))
    print "\nTotal number of mesh cells: \t%g" % (mesh.n_cells ** 2)
    #f.write("\nTotal number of unknowns: \t%g" % ((mesh.n_cells ** 2) * (order*(order+2)/2)))
    print "\nTotal number of unknowns: \t%g" % ((mesh.n_cells ** 2) * (order*(order+2)/2))

    # give order, mesh to solver
    solve = solver.SN(order, mesh.cells, spacing, mesh.n_cells, mesh.n_fuel, mesh.n_mod, tol, mesh.fuel_area,
                      mesh.mod_area, timestr, mesh.top_right_corner_fuel)
    SNresults = solve.solveSN(max_iters, plotter, mesh, savepath, update_source)

    if SNresults:
        """
        f.write(
            "\nAvg fuel flux \t %f \nAvg mod flux \t %f \nAverage Flux  \t %f \nFlux ratio \t %f"
            "\nfuel corner flux\t %g\nCorner normalized\t %g\n\n"
            "Converged in %d iterations! \tL2 \t%g \n"
            % (  solve.results[1], solve.results[2], solve.results[3],
               solve.results[4], solve.results[5], solve.results[5]/q_fuel, solve.results[0], solve.l2))
        """
        f.write(
            "Avgfuelfl\t  Avgmodfl\t AvgFlux\t  "
            "Fluxrat\t  cornerfl\t Corner_norm\t \n"
            "%f\t%f\t%f\t%f\t%g\t%g\n\n\n"
            "Converged in %d iterations! \tL2 \t%g \n"
            % (solve.results[1], solve.results[2], solve.results[3],
               solve.results[4], solve.results[5], solve.results[5] / q_fuel, solve.results[0], solve.l2))
        #Corner flux over fuel source is normalized
        return solve.results[4]

    elif not SNresults:
        f.write("\n*********************************\n"
                "Not converged after %d iterations. Rerun this case with more iterations. \nL2 \t%g" %(max_iters, solve.l2))
        f.write(
            "\nAvg fuel flux \t %f \nAvg mod flux \t %f \nAverage Flux  \t %f \nFlux ratio \t %f"
            "\nfuel corner flux\t %g\nCorner normalized\t %g\n"
            "\n ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** * \n"
            % (solve.results[1], solve.results[2], solve.results[3],
               solve.results[4], solve.results[5], solve.results[5] / q_fuel))
        f.close()
        return solve.results[4]


def solveSpacings(spacings, order):
    f = open('%s.txt' % resultsfile, 'a+')
    ratio_old = 0
    print "Iterating over spacings....\n\n"
    f.write("Iterating over spacings...\n\n")

    for spacing in spacings:
        savepath = pathname + '/spacing_' + str(spacing) + 'order_' + str(order)
        plotter.mkdir_p(savepath)

        result = solveSN(spacing, order, savepath)

        fluxchg = ((result - ratio_old) / result) * 100
        print "flux ratio change: %g" % (fluxchg)
        f.write("flux ratio change: %g\n\n\n" % (fluxchg))

        if convergemesh:
            if fluxchg <= sptol:
                print "Spatial mesh is converged with spacing of %g" % (spacing)
                f.write("Spatial mesh is converged with spacing of %g" % (spacing))
                break
        ratio_old = result
    f.close()

def solveOrders(spacing, orders):
    f = open('%s.txt' % resultsfile, 'a+')
    ratio_old = 0
    print "Iterating over orders...\n\n"
    f.write("Iterating over orders...\n\n")

    for order in orders:
        savepath = pathname + '/order_' + str(order) + 'spacing_' + str(spacing)
        plotter.mkdir_p(savepath)

        result = solveSN(spacing, order, savepath)
        fluxchg = ((result - ratio_old) / result) * 100

        print "flux ratio change: %g" % (fluxchg)
        f.write("flux ratio change: %g\n\n\n" % (fluxchg))
        ratio_old = result
    f.close()

def getNumUnknowns(spacing, order):
    # setup mesh cells
    f = open('%s.txt' % resultsfile, 'a+')
    mesh = geometry.Geometry(pitch, spacing, fwidth, fuel, moderator)
    mesh.setMesh(tally_fuel_corner)
    f.write("SN order\t%g\tspacing\t%g\tnum cells\t%g\tnum unknowns\t%g\n"
            % (order, spacing, (mesh.n_cells ** 2), (mesh.n_cells ** 2) * order))
    print "\nTotal number of unknowns: \t%g" % ((mesh.n_cells ** 2) * order)
    f.close()

if singlesolve:
    solveSN(spacing, order, savepath)


if getnumunknowns and manyspacings and manyorders:
    for ord in orders:
        for space in spacings:
            getNumUnknowns(space, ord)
elif manyspacings and manyorders:
    for item in spacings:
        solveOrders(item, orders)
elif manyorders:
    solveOrders(spacing, orders)
elif manyspacings:
    solveSpacings(spacings, order)


f.close()

plotter.saveInputFile(savepath)
