from quadrature import LevelSymmetricQuadrature
from convergence import ConvergenceTest


class SN(object):
    """
    Need to:
    - add iteration to convergence for solveSN
    - import L2 convergence check class
    """
    def __init__(self, order, cells, mesh_spacing, num_cells, num_fuel, num_mod, tol):
        self.quadrature = LevelSymmetricQuadrature().getQuadrature(order)
        self.n_angles = self.quadrature['n_angles']
        self.cells = cells
        self.d_cell = mesh_spacing
        self.n_cells = num_cells
        self.n_fuel = num_fuel
        self.n_mod = num_mod
        self.tol = tol

    def solveSN(self):
        check = ConvergenceTest()
        converged = False
        scalar_flux_old = [0,0]

        while not converged:

            #sweeps across y cells for given x, over all x cells

            #sweep starting from bottom left
            print "Sweeping from bottom left..."
            for i in range(self.n_cells):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]
                    for angles in range(self.n_angles / 4):
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[0], cell.angular[1])

                        # set right flux out (i+i cell left flux in)
                        cell.angular[2] = self.getAngFluxOut(avg, cell.angular[0])
                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0] = cell.angular[2]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3] = self.getAngFluxOut(avg, cell.angular[1])
                        if j < self.n_cells - 1:
                            self.cells[i][j+1].angular[1] = cell.angular[3]

            #sweep starting from bottom right
            print "Sweeping from bottom right..."
            for i in reversed(range(self.n_cells)):
                for j in range(self.n_cells):
                    cell = self.cells[i][j]
                    for angles in range(self.n_angles / 4):
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[2], cell.angular[1])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0] = self.getAngFluxOut(avg, cell.angular[2])
                        if i > self.n_cells+1:
                            self.cells[i-1][j].angular[2] = cell.angular[0]

                        # set top flux out (j+1 cell bottom flux in)
                        cell.angular[3] = self.getAngFluxOut(avg, cell.angular[1])
                        if j < self.n_cells - 1:
                            self.cells[i][j + 1].angular[1] = cell.angular[3]

            #sweep from top right
            print "Sweeping from top right..."
            for i in reversed(range(self.n_cells)):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]
                    for angles in range(self.n_angles / 4):
                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[2], cell.angular[3])

                        # set left flux out (i-1 cell right flux in)
                        cell.angular[0] = self.getAngFluxOut(avg, cell.angular[2])
                        if i >= 1:
                            self.cells[i-1][j].angular[2] = cell.angular[0]

                        # set bottom flux out (j-1 cell top flux in)
                        cell.angular[1] = self.getAngFluxOut(avg, cell.angular[3])
                        if j >= 1:
                            self.cells[i][j-1].angular[3] = cell.angular[1]


            #sweep from top left
            print "Sweeping from top left..."
            for i in range(self.n_cells):
                for j in reversed(range(self.n_cells)):
                    cell = self.cells[i][j]
                    for angles in range(self.n_angles / 4):

                        eta = self.quadrature['eta'][angles]
                        xi = self.quadrature['xi'][angles]

                        avg = self.getCellAvgFlux(cell.material.q, eta, xi, cell.material.xs, cell.angular[0], cell.angular[3])

                        # set right flux out (i+i cell left flux in)
                        cell.angular[2] = self.getAngFluxOut(avg, cell.angular[0])
                        if i < self.n_cells-1:
                            self.cells[i+1][j].angular[0] = cell.angular[2]

                        # set bottom flux out (j-1 cell top flux in)
                        cell.angular[1] = self.getAngFluxOut(avg, cell.angular[3])
                        if j >= 1:
                            self.cells[i][j - 1].angular[3] = cell.angular[1]

            #update scalar flux for each cell, then zero angular fluxes

            for i in range(self.n_cells):
                for cell in self.cells[i]:
                    cell.flux += self.quadrature['weight'][angles] * \
                                 (cell.angular[0] + cell.angular[2]) / 2
                    print(cell.flux)
                    for m in range(4):
                        cell.angular[m] = 0
            scalar_flux = self.getAvgScalarFlux()

            converged = check.isConverged(scalar_flux, scalar_flux_old, self.tol)
            if not converged:
                scalar_flux_old = self.getAvgScalarFlux()


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
        fuelflux /= self.n_fuel
        modflux /= self.n_mod
        print "Avg fuel flux = %f \nAvg mod flux = %f" %(fuelflux, modflux)
        return fuelflux, modflux