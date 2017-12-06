from quadrature import LevelSymmetricQuadrature
from convergence import ConvergenceTest
from os import mkdir
import time


class SN(object):

    def __init__(self, order, cells, mesh_spacing, num_cells, num_fuel, num_mod, tol, fuel_area, mod_area):
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


    def solveSN(self, max_iters, plotter, mesh):
        check = ConvergenceTest()
        converged = False
        scalar_flux_old = [0,0]
        iters = 0
        na_oct = self.n_angles // 4 #number of angles per octant
        print "fuel boundaries: %g, %g" %(self.n_mod, self.n_mod + self.n_fuel)

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
            #"""
            for i in range(self.n_cells):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):

                        #angles from [0, n_ang/4]
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        #avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[0][angles], cell.angular[1][angles])
                        avg = self.getCellAvgFlux(cell.source, eta, xi, cell.material.xs, cell.angular[0][angles], cell.angular[1][angles])

                        # set right flux out (i+i cell left flux in)
                        cell.angular[2][angles] = self.getAngFluxOut(avg, cell.angular[0][angles])

                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0][angles] = cell.angular[2][angles]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3][angles] = self.getAngFluxOut(avg, cell.angular[1][angles])
                        if j < self.n_cells - 1:
                            self.cells[i][j+1].angular[1][angles] = cell.angular[3][angles]

            #reflect on right, top boundaries
            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[2][angles + na_oct] = self.cells[i][j].angular[2][angles]

            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[3][angles + na_oct*3] = self.cells[i][j].angular[3][angles]
            
            #"""

            # sweep starting from bottom right
            for i in reversed(range(self.n_cells)):
                for j in (range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        # angles from [n_ang/4, n_ang/2]
                        angles += na_oct
                        avg = self.getCellAvgFlux(cell.source, eta, xi, cell.material.xs, cell.angular[2][angles-na_oct], cell.angular[1][angles])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0][angles] = self.getAngFluxOut(avg, cell.angular[2][angles])
                        if i > self.n_cells + 1:
                            self.cells[i - 1][j].angular[2][angles] = cell.angular[0][angles]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3][angles] = self.getAngFluxOut(avg, cell.angular[1][angles])
                        if j < self.n_cells - 1:
                            self.cells[i][j + 1].angular[1][angles] = cell.angular[3][angles]

            # reflect on left, top boundaries
            i=0
            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[0][angles] = self.cells[i][j].angular[0][angles + na_oct]

            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[3][angles + na_oct * 2] = self.cells[i][j].angular[3][angles + na_oct]

            #sweep from top right
            for i in reversed(range(self.n_cells)):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]
                        # angles from [3/4n , n]
                        angles += 2*na_oct

                        avg = self.getCellAvgFlux(cell.source, eta, xi, cell.material.xs, cell.angular[2][angles-na_oct], cell.angular[3][angles-na_oct])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0][angles] = self.getAngFluxOut(avg, cell.angular[2][angles])
                        if i >= 1:
                            self.cells[i-1][j].angular[2][angles] = cell.angular[0][angles]

                        # set bottom flux out (j-1 cell top flux in)
                        cell.angular[1][angles] = self.getAngFluxOut(avg, cell.angular[3][angles])
                        if j >= 1:
                            self.cells[i][j-1].angular[3][angles] = cell.angular[1][angles]

            # reflect on left, bottom boundaries
            i = 0
            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[0][angles + 3 * na_oct] = self.cells[i][j].angular[0][angles + 2* na_oct]

            j = 0
            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[1][angles + na_oct] = self.cells[i][j].angular[1][angles + 2 * na_oct]


            #sweep from top left

            for i in range(self.n_cells):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(na_oct):
                        #angles from [n/2, 3/4n]

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]
                        angles += 3 * na_oct

                        avg = self.getCellAvgFlux(cell.source, eta, xi, cell.material.xs, cell.angular[0][angles-na_oct/2], cell.angular[3][angles-3*na_oct])

                        # set right flux out (and pass to i+i cell as left flux in)
                        cell.angular[2][angles] = self.getAngFluxOut(avg, cell.angular[0][angles])
                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0][angles] = cell.angular[2][angles]

                        # set bottom flux out (and pass to j-1 cell as top flux in)
                        cell.angular[1][angles] = self.getAngFluxOut(avg, cell.angular[3][angles])
                        if j >= 1:
                            self.cells[i][j - 1].angular[3][angles] = cell.angular[1][angles]

            # reflect on right, bottom boundaries
            for j in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[2][angles + 2 * na_oct] = self.cells[i][j].angular[2][angles + 3 * na_oct]

            j = 0
            for i in range(self.n_cells):
                for angles in range(na_oct):
                    self.cells[i][j].angular[1][angles] = self.cells[i][j].angular[1][angles + 2 * na_oct]


            #update scalar flux for each cell

            for i in range(self.n_cells):
                for cell in self.cells[i]:
                    for angles in range(na_oct):
                        for allangles in range(na_oct*4):
                            cell.flux += self.quadrature['weight'][angles] * \
                                 (cell.angular[0][allangles] + cell.angular[2][allangles]) / 2

            getfluxes = self.getAvgScalarFlux()
            scalar_flux = getfluxes[:2]

            #timestr = time.strftime("%Y-%m-%d %H:%M")
            #mkdir(timestr)

            #print "plotting scalar flux for iteration %d" % (iters)
            plotter.plotScalarFlux(mesh, self.order, self.d_cell, iters)

            #midpt = self.n_cells/2 - 5
            midpt = self.n_cells/2 - 1

            plotter.plotCenterFlux(mesh, self.cells, midpt, iters, self.order)
            plotter.plotCenterFluxY(mesh, self.cells, midpt, iters, self.order)

            converged = check.isConverged(scalar_flux, scalar_flux_old, self.tol)

            if not converged:
                scalar_flux_old = scalar_flux[:]#self.getAvgScalarFlux()

                for i in range(self.n_cells):
                    for cell in self.cells[i]:
                        cell.flux = 0
                if iters == max_iters:
                    converged = True
                    print "Not converged after %d iterations." %(iters)
            else:
                print "Converged in %d iterations\n" %(iters)
                self.results = self.returnSolveResults(iters, getfluxes[0], getfluxes[1], getfluxes[2], getfluxes[3])


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
        #avgflux = (fuelflux + modflux) / self.n_cells
        avgflux = scalarflux / self.n_cells
        maxflux = 0.0

        #normalize flux
        for i in range(self.n_cells):
            for cell in self.cells[i]:
                #cell.flux /= avgflux
                if cell.flux > maxflux:
                    maxflux = cell.flux
                if cell.region == 'fuel':
                    # accumulate scalar flux avg for fuel
                    fuelflux += cell.flux
                else:
                    # accumulate scalar flux avg for mod
                    modflux += cell.flux

        """
        for i in range(self.n_cells):
            for cell in self.cells[i]:
                if cell.region == 'fuel':
                    # accumulate scalar flux avg for fuel
                    fuelflux += cell.flux
                else:
                    # accumulate scalar flux avg for mod
                    modflux += cell.flux
        """

        #fuelflux /= self.fuel_area
        fuelflux /= self.n_fuel
        #fuelflux /= maxflux
        #modflux /= self.mod_area
        modflux /= self.n_mod
        #modflux /= maxflux
        ratio = fuelflux / modflux

        print "Avg fuel flux = %f \nAvg mod flux = %f \nAverage Flux  = %f \nFlux ratio = %f" %(fuelflux, modflux, avgflux, ratio)
        return fuelflux, modflux, avgflux, ratio

    def returnSolveResults(self, iters, fuelflux, modflux, avgflux, ratio):
       # print "Iterations to convergence: %d \nModerator flux: %g\nFuel flux: %g\n Avg flux: %g\n Flux ratio: %g" \
        #      %(iters, modflux, fuelflux, avgflux, ratio)
        return iters, fuelflux, modflux, avgflux, ratio