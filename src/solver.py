from quadrature import LevelSymmetricQuadrature
from convergence import ConvergenceTest


class SN(object):

    def __init__(self, order, cells, mesh_spacing, num_cells, num_fuel, num_mod, tol, fuel_area, mod_area):
        self.quadrature = LevelSymmetricQuadrature().getQuadrature(order)
        self.n_angles = self.quadrature['n_angles']
        self.cells = cells
        self.d_cell = mesh_spacing
        self.n_cells = num_cells
        self.n_fuel = num_fuel
        self.n_mod = num_mod
        self.tol = tol
        self.fuel_area = fuel_area
        self.mod_area = mod_area

    def solveSN(self, max_iters):
        check = ConvergenceTest()
        converged = False
        scalar_flux_old = [0,0]
        iters = 0

        while not converged:
            iters +=1
            #sweeps across y cells for given x, over all x cells


            #sweep starting from bottom left

            for i in range(self.n_cells):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]

                    for angles in range(self.n_angles // 4):
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[0][angles], cell.angular[1][angles])

                        # set right flux out (i+i cell left flux in)
                        cell.angular[2][angles] = self.getAngFluxOut(avg, cell.angular[0][angles])
                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0][angles] = cell.angular[2][angles]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3][angles] = self.getAngFluxOut(avg, cell.angular[1][angles])
                        if j < self.n_cells - 1:
                            self.cells[i][j+1].angular[1] = cell.angular[3]
            """
            # reflect on right and top boundaries
            for angles in range(self.n_angles // 4):
                for j in range(self.n_cells):
                    cell = self.cells[-1][j]
                    cell.angular[2][angles] = cell.angular[0][angles]
                for i in range(self.n_cells):
                    cell = self.cells[i][-1]
                    cell.angular[3][angles] = cell.angular[1][angles]
            """
            #sweep starting from bottom right
            for i in reversed(range(self.n_cells)):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]

                    for angles in range(self.n_angles // 4):
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[2][angles], cell.angular[1][angles])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0][angles] = self.getAngFluxOut(avg, cell.angular[2][angles])
                        if i > self.n_cells+1:
                            self.cells[i-1][j].angular[2][angles] = cell.angular[0][angles]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3][angles] = self.getAngFluxOut(avg, cell.angular[1][angles])
                        if j < self.n_cells - 1:
                            self.cells[i][j + 1].angular[1][angles] = cell.angular[3][angles]
            """
            # reflect on left and top boundaries
            for angles in range(self.n_angles // 4):
                for j in range(self.n_cells):
                    cell = self.cells[0][j]
                    cell.angular[0][angles] = cell.angular[2][angles]
                for i in range(self.n_cells):
                    cell = self.cells[i][-1]
                    cell.angular[1][angles + self.n_angles//4] = cell.angular[3][angles + self.n_angles//4]
            """
            #sweep from top right
            for i in reversed(range(self.n_cells)):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(self.n_angles // 4):
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[2][angles], cell.angular[3][angles])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0][angles] = self.getAngFluxOut(avg, cell.angular[2][angles])
                        if i >= 1:
                            self.cells[i-1][j].angular[2][angles] = cell.angular[0][angles]

                        # set bottom flux out (j-1 cell top flux in)
                        cell.angular[1][angles] = self.getAngFluxOut(avg, cell.angular[3][angles])
                        if j >= 1:
                            self.cells[i][j-1].angular[3][angles] = cell.angular[1][angles]
            """
            # reflect on left and bottom boundaries
            for angles in range(self.n_angles // 4):
                for j in range(self.n_cells):
                    cell = self.cells[0][j]
                    cell.angular[0][angles] = cell.angular[2][angles]
                for i in range(self.n_cells):
                    cell = self.cells[i][0]
                    cell.angular[1][angles] = cell.angular[3][angles]
            """
            #sweep from top left

            for i in range(self.n_cells):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]

                    for angles in range(self.n_angles // 4):

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[0][angles], cell.angular[3][angles])

                        # set right flux out (i+i cell left flux in)
                        cell.angular[2][angles] = self.getAngFluxOut(avg, cell.angular[0][angles])
                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0][angles] = cell.angular[2][angles]

                        # set bottom flux out (j-1 cell top flux in)
                        cell.angular[1][angles] = self.getAngFluxOut(avg, cell.angular[3][angles])
                        if j >= 1:
                            self.cells[i][j - 1].angular[3][angles] = cell.angular[1][angles]

            #update scalar flux for each cell

            for i in range(self.n_cells):
                for cell in self.cells[i]:
                    for angles in range(self.n_angles // 4):
                        cell.flux += self.quadrature['weight'][angles] * \
                                 (cell.angular[0][angles] + cell.angular[2][angles]) / 2

            scalar_flux = self.getAvgScalarFlux()

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
                print "Converged in %d iterations" %(iters)


    def getCellAvgFlux(self, q, eta, xi, sigma, psi_h, psi_v):
        delta = self.d_cell
        tet = 2 * eta / delta
        txi = 2 * xi / delta
        num = q + tet * psi_h + txi * psi_v
        den = sigma + tet + txi
        return num/den

    def getAngFluxOut(self, avgflux, influx):
        return 2 * avgflux - influx

    def getAvgScalarFlux(self):
        fuelflux = 0
        modflux = 0

        for i in range(self.n_cells):
            for cell in self.cells[i]:
                if cell.region == 'fuel':
                    #accumulate scalar flux avg for fuel
                    fuelflux += cell.flux
                else:
                    #accumulate scalar flux avg for mod
                    modflux += cell.flux

        fuelflux /= self.fuel_area#self.n_fuel #* cell_area
        modflux /= self.mod_area #self.n_mod #* cell_area
        print "Avg fuel flux = %f \nAvg mod flux = %f" %(fuelflux, modflux)
        return fuelflux, modflux