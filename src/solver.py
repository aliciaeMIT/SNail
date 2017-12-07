from quadrature import LevelSymmetricQuadrature
from convergence import ConvergenceTest
from os import mkdir
import time


class SN(object):

    def __init__(self, order, cells, mesh_spacing, num_cells, num_fuel, num_mod, tol, fuel_area, mod_area, resultsfile):
        self.quadrature = LevelSymmetricQuadrature().getQuadrature(order)
        self.order = order
        self.n_angles = self.quadrature['n_angles']
        self.cells = cells
        self.d_cell = mesh_spacing
        self.n_cells = num_cells
        self.n_fuel = num_fuel
        self.n_mod = num_mod
        self.tol = tol
        self.fuel_area = fuel_area
        self.mod_area = mod_area
        self.results = 0
        self.resultsfile = resultsfile + '_results'


    def solveSN(self, max_iters, plotter, mesh, savepath):

        check = ConvergenceTest()
        converged = False
        scalar_flux_old = [0,0]
        iters = 0
        na_oct = self.n_angles // 4 #number of angles per octant
        #print "fuel boundaries: %g, %g" %(self.n_mod, self.n_mod + self.n_fuel)

        #initialize scalar flux guesses for source calculation: phi = q / sigma_absorption for mod, phi = q for fuel
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                cell = self.cells[i][j]
                if cell.region == 'fuel':
                    absorpt = 1
                else:
                    absorpt = cell.material.xs - cell.material.scatter
                cell.flux = cell.material.q / absorpt
                #print "Initial cell flux guess: %g" %(cell.flux)

        while not converged:

            iters +=1
            print "Iteration %d" %(iters)

            #update source for each cell
            self.updateAllCellSource()

            #reinitialize fluxes to zero
            for i in range(self.n_cells):
                for cell in self.cells[i]:
                    cell.flux = 0

            #sweeps across y cells for given x, over all x cells
            #sweep starting from bottom left

            for i in range(self.n_cells):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):

                        #angles from [0, n_ang/4]
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        cell.avg_angular[angles] = self.getCellAvgFlux(cell.source, eta, xi, cell.material.xs, cell.angular[0][angles], cell.angular[1][angles])

                        # set right flux out (i+i cell left flux in)
                        cell.angular[2][angles] = self.getAngFluxOut(cell.avg_angular[angles], cell.angular[0][angles])

                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0][angles] = cell.angular[2][angles]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3][angles] = self.getAngFluxOut(cell.avg_angular[angles], cell.angular[1][angles])
                        if j < self.n_cells - 1:
                            self.cells[i][j+1].angular[1][angles] = cell.angular[3][angles]
            
            #reflect on right, top boundaries
            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[2][angles + na_oct] = self.cells[i][j].angular[2][angles]

            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[3][angles + na_oct] = self.cells[i][j].angular[3][angles]


            # sweep starting from bottom right
            for i in reversed(range(self.n_cells)):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        # angles from [n_ang/4, n_ang/2]

                        cell.avg_angular[angles + na_oct] = self.getCellAvgFlux(cell.source, eta, xi,
                                                                                cell.material.xs,
                                                                                cell.angular[2][
                                                                                    angles + na_oct],
                                                                                cell.angular[1][
                                                                                    angles + na_oct])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0][angles + na_oct] = self.getAngFluxOut(cell.avg_angular[angles + na_oct],
                                                                              cell.angular[2][angles + na_oct])
                        if i >= 1:  # self.n_cells + 1:
                            self.cells[i - 1][j].angular[2][angles + na_oct] = cell.angular[0][angles + na_oct]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3][angles + na_oct] = self.getAngFluxOut(cell.avg_angular[angles + na_oct],
                                                                              cell.angular[1][angles + na_oct])
                        if j < self.n_cells - 1:
                            self.cells[i][j + 1].angular[1][angles + na_oct] = cell.angular[3][angles + na_oct]

            # reflect on left, top boundaries

            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[0][j].angular[0][angles] = self.cells[0][j].angular[0][angles + na_oct]

            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][self.n_cells - 1].angular[3][angles] = \
                    self.cells[i][self.n_cells - 1].angular[3][angles + na_oct]

            #sweep from top right
            for i in reversed(range(self.n_cells)):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]
                        # angles from [n/2 , 3/4n]


                        cell.avg_angular[angles+2*na_oct] = self.getCellAvgFlux(cell.source, eta, xi, cell.material.xs,
                                                                                cell.angular[2][angles], cell.angular[3][angles])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0][angles] = self.getAngFluxOut(cell.avg_angular[angles+2*na_oct], cell.angular[2][angles])
                        if i >= 1:
                            self.cells[i-1][j].angular[2][angles] = cell.angular[0][angles]

                        # set bottom flux out (j-1 cell top flux in)
                        cell.angular[1][angles] = self.getAngFluxOut(cell.avg_angular[angles+2*na_oct], cell.angular[3][angles])
                        if j >= 1:
                            self.cells[i][j-1].angular[3][angles] = cell.angular[1][angles]

            # reflect on left, bottom boundaries


            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[0][j].angular[0][angles + na_oct] = self.cells[0][j].angular[0][angles]


            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][0].angular[1][angles + na_oct] = self.cells[i][0].angular[1][angles]


            # sweep from top left

            for i in range(self.n_cells):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):
                        # angles from [n/2, 3/4n]

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]


                        cell.avg_angular[angles + 3 * na_oct] = self.getCellAvgFlux(cell.source, eta, xi,
                                                                                    cell.material.xs,
                                                                                    cell.angular[0][
                                                                                        angles + na_oct],
                                                                                    cell.angular[3][angles + na_oct])

                        # set right flux out (and pass to i+i cell as left flux in)
                        cell.angular[2][angles + na_oct] = self.getAngFluxOut(cell.avg_angular[angles + 3 * na_oct],
                                                                              cell.angular[0][angles + na_oct])
                        if i < self.n_cells - 1:
                            self.cells[i + 1][j].angular[0][angles + na_oct] = cell.angular[2][angles + na_oct]

                        # set bottom flux out (and pass to j-1 cell as top flux in)
                        cell.angular[1][angles + na_oct] = self.getAngFluxOut(cell.avg_angular[angles + 3 * na_oct],
                                                                     cell.angular[3][angles + na_oct])
                        if j >= 1:
                            self.cells[i][j - 1].angular[3][angles + na_oct] = cell.angular[1][angles + na_oct]

            # reflect on right, bottom boundaries
            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[self.n_cells - 1][j].angular[2][angles] = self.cells[self.n_cells - 1][j].angular[2][angles + na_oct]


            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][0].angular[1][angles] = self.cells[i][0].angular[1][angles + na_oct]

            #update scalar flux for each cell

            for i in range(self.n_cells):
                for cell in self.cells[i]:
                    for angles in range(na_oct):
                        cell.flux +=self.quadrature['weight'][angles] * (cell.avg_angular[angles] + cell.avg_angular[angles+na_oct]
                                                                         + cell.avg_angular[angles+2*na_oct] + cell.avg_angular[angles+3*na_oct])

            getfluxes = self.getAvgScalarFlux()
            scalar_flux = getfluxes[:2]

            #print "plotting scalar flux for iteration %d" % (iters)

            midpt = self.n_cells/2 - 1

            plotter.plotScalarFlux(mesh, self.order, self.d_cell, iters, savepath)
            plotter.plotCenterFlux(mesh, self.cells, midpt, iters, self.order, savepath)
            plotter.plotCenterFluxY(mesh, self.cells, midpt, iters, self.order, savepath)

            converged = check.isConverged(scalar_flux, scalar_flux_old, self.tol)

            if not converged:
                scalar_flux_old = scalar_flux[:]

                for i in range(self.n_cells):
                    for cell in self.cells[i]:
                        cell.flux = 0
                if iters == max_iters:
                    converged = True
                    print "Not converged after %d iterations." %(iters)
            else:
                print "Converged in %d iterations\n" %(iters)
                self.results = self.returnSolveResults(iters, getfluxes[0], getfluxes[1], getfluxes[2], getfluxes[3])

                # plot angular flux for center, center-edge fuel, corner fuel, and corner pincell cells
                # center
                ii = self.n_cells / 2
                #fuel right edge
                jj = self.n_cells/2 + self.n_fuel / 2 - 1

                cells = (self.cells[ii][ii], self.cells[jj][ii], self.cells[jj][self.n_mod], self.cells[self.n_cells-1][0])
                location = ['center', 'center-edge', 'fuel-corner', 'pincell-corner']
                for k, cell in enumerate(cells):
                    plotter.plotAngularFlux(cell, self.quadrature, location[k], savepath)


    def getCellAvgFlux(self, q, eta, xi, sigma, psi_h, psi_v):
        delta = self.d_cell
        tet = 2 * eta / delta
        txi = 2 * xi / delta
        num = q + tet * psi_h + txi * psi_v
        den = sigma + tet + txi
        return num/den

    def getAngFluxOut(self, avgflux, influx):
        return 2 * avgflux - influx

    def updateSource(self, q, flux, scatter):
        #calc = ((1/2) * (scatter * flux + q))
        return ((1./2.) * (scatter * flux + q))

    def updateAllCellSource(self):
        for i in range(self.n_cells):
            for j in range(self.n_cells):
                cell = self.cells[i][j]
               # print "cell [%d][%d] source update:\nOld source; %g" % (i, j, cell.source)

                cell.source = self.updateSource(cell.material.q, cell.flux, cell.material.scatter)
                #sourceval = cell.source
                #print "new source; %g" %(cell.source)

    def getAvgScalarFlux(self):
        fuelflux = 0.0
        modflux = 0.0
        scalarflux = 0.0

        for i in range(self.n_cells):
            for cell in self.cells[i]:
               scalarflux += cell.flux
        avgflux = scalarflux / (self.n_cells ** 2)
        maxflux = 0.0
        fuelcell = 0
        modcell = 0

        #normalize flux
        for i in range(self.n_cells):
            for cell in self.cells[i]:
                cell.flux /= avgflux
                if cell.flux > maxflux:
                    maxflux = cell.flux
                if cell.region == 'fuel':
                    # accumulate scalar flux avg for fuel
                    fuelflux += cell.flux
                    fuelcell+=1
                else:
                    # accumulate scalar flux avg for mod
                    modflux += cell.flux
                    modcell +=1

        fuelflux /= fuelcell
        #fuelflux /= maxflux
        modflux /= modcell

        #modflux /= maxflux
        ratio = fuelflux / modflux

        print "Avg fuel flux = %f \nAvg mod flux = %f \nAverage Flux  = %f \nFlux ratio = %f" %(fuelflux, modflux, avgflux, ratio)
        return fuelflux, modflux, avgflux, ratio

    def returnSolveResults(self, iters, fuelflux, modflux, avgflux, ratio):
       # print "Iterations to convergence: %d \nModerator flux: %g\nFuel flux: %g\n Avg flux: %g\n Flux ratio: %g" \
        #      %(iters, modflux, fuelflux, avgflux, ratio)
        return iters, fuelflux, modflux, avgflux, ratio