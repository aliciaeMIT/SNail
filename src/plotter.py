
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from math import *
import matplotlib.pyplot as plt
from errno import EEXIST
from os import makedirs,path
import shutil
"""
Plotting tool for SN to plot mesh, scalar flux, angular flux

plotMaterial, plotScalarFlux, and plotAngularFlux methods are adapted from Samuel Shaner's Sn solver, DiscOrd:
https://github.com/samuelshaner/DiscOrd/

"""

def plotCenterFlux(mesh, solved, j, iter, order, savepath):
    fignum = 1
    plt.figure(fignum)
    ivals = []
    cellfluxes = []
    fuelbounds = [(mesh.n_mod-0.5), (mesh.n_fuel + mesh.n_mod - 0.5)]

    for xc in fuelbounds:
        plt.axvline(x=xc, color='k', linestyle = '--')

    for i in range(mesh.n_cells):
        ivals.append(i)
        cellfluxes.append(solved[i][j].flux)
    plt.plot(ivals, cellfluxes)
    plt.ylabel('Scalar Flux')
    plt.xlabel('X node # across centerline (Y node num ' + str(j) + ' )')
    plt.title('Horizontal centerline flux, iteration ' + str(iter) + ', order ' + str(order))
    plt.savefig(savepath + '/horiz_iter' + str(iter) + '_flux_center.png')
    #plt.close()

def plotCenterFluxY(mesh, solved, i, iter, order, savepath):
    fignum = 2
    plt.figure(fignum)
    jvals = []
    cellfluxes = []
    fuelbounds = [(mesh.n_mod - 0.5), (mesh.n_fuel + mesh.n_mod - 0.5)]
    #0.5 is subtracted to account for boundary being between cells for fuel and mod indices

    for xc in fuelbounds:
        plt.axvline(x=xc, color='k', linestyle = '--')

    for j in range(mesh.n_cells):
        jvals.append(j)
        cellfluxes.append(solved[i][j].flux)
    plt.plot(jvals, cellfluxes)
    plt.ylabel('Scalar Flux')
    plt.xlabel('Y node # across centerline (X node num ' + str(i) + ' )')
    plt.title('Vertical centerline flux, iteration ' + str(iter) + ', order ' + str(order))
    plt.savefig(savepath + '/vert_iter' + str(iter) + '_flux_center.png')
    #plt.close()

def plotMaterial(mesh, spacing, plot_cells, savepath):


    bit_size = 500.0

    # create image
    bit_length = round(bit_size / mesh.n_cells)
    size = int(bit_length * mesh.n_cells)
    img = Image.new('RGB', (size,size), 'white')
    draw = ImageDraw.Draw(img)

    # draw cell interior
    for i in range(mesh.n_cells):
        for j in range(mesh.n_cells):
            cell = mesh.cells[i][j]
            # fuel red; moderator blue
            if cell.region == 'fuel':
                draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (255,0,0))
            else:
                draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (0,0,255))

            # find cells where angular flux needs to be plotted; make them white
            #if mesh.cells[i][j].id in plot_cells:
            #    draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (255,255,255))

    # draw horizontal grid lines
    for j in range(1,mesh.n_cells):
        draw.line((0, j*bit_length, size,j*bit_length), fill=(0,0,0), width=1)

    # draw vertical grid lines
    for i in range(1,mesh.n_cells):
        draw.line((i*bit_length, 0, i*bit_length, size), fill=(0,0,0), width=1)

    # save image
    img.save(savepath + '/material_' + str(spacing)[2:] + '.png')

def plotScalarFlux(mesh, order, spacing, iteration, savepath):

    # create image
    bit_length = round(500.0 / mesh.n_cells)
    size = int(bit_length * mesh.n_cells)
    img = Image.new('RGB', (size,size), 'white')
    draw = ImageDraw.Draw(img)

    # find max and min flux
    max_flux = mesh.cells[0][0].flux
    min_flux = mesh.cells[0][0].flux
    for i in range(mesh.n_cells):
        for j in range(mesh.n_cells):
            max_flux = max(mesh.cells[i][j].flux, max_flux)
            min_flux = min(mesh.cells[i][j].flux, min_flux)

    flux_range = max_flux - min_flux

    # draw flux map
    for i in range(mesh.n_cells):
        for j in range(mesh.n_cells):
            cell = mesh.cells[i][j]


            # get color
            if ((cell.flux-min_flux) / flux_range <= 1.0/3.0):
                red = 0.0
                green = 3.0 * (cell.flux-min_flux) / flux_range
                blue = 1.0
            elif ((cell.flux-min_flux) / flux_range <= 2.0/3.0):
                red = 3.0 * (cell.flux-min_flux) / flux_range - 1.0
                green = 1.0
                blue = -3.0 * (cell.flux-min_flux) / flux_range + 2.0
            else:
                red = 1.0
                green = -3.0 * (cell.flux-min_flux) / flux_range + 3.0
                blue = 0.0

            # convert color to RGB triplet
            red = int(255*red)
            green = int(255*green)
            blue = int(255*blue)

            # draw pin and pin power
            draw.rectangle([i*bit_length, size - j*bit_length - bit_length, (i+1)*bit_length, size - j*bit_length], (red,green,blue))

    # save image
    img.save(savepath + '/flux_' + str(spacing)[2:] + '_' + str(int(floor(order/10))) + str(order % 10) + '_' + str(int(floor(iteration/10))) + str(iteration % 10) + '.png')


def plotAngularFlux(cell, quad, location, savepath):

    # create image
    img = Image.new('RGB', (500,500), 'white')
    draw = ImageDraw.Draw(img)

    # get xi, eta's, and index for smallest polar angle (mu)
    mu_min = min(quad['mu'])
    min_mu_indices = []
    for i,a in enumerate(quad['mu']):
        if abs(a - mu_min) < .0001:
            min_mu_indices.append(i)

    # get the angular fluxes for the selected angles and make theta array
    ang_fluxes = []
    theta = []

    for i in min_mu_indices:
        ang_fluxes.append(cell.avg_angular[i])
        theta.append(acos(quad['eta'][i]))

    for i in min_mu_indices:
        ang_fluxes.append(cell.avg_angular[quad['n_angles_per_octant'] + i])
        theta.append(pi - acos(quad['eta'][i]))

    for i in min_mu_indices:
        ang_fluxes.append(cell.avg_angular[2*quad['n_angles_per_octant'] + i])
        theta.append(pi + acos(quad['eta'][i]))

    for i in min_mu_indices:
        ang_fluxes.append(cell.avg_angular[3*quad['n_angles_per_octant'] + i])
        theta.append(2*pi - acos(quad['eta'][i]))


    # rarrange theta array from min to max
    theta_min = min(theta)
    theta_sort = []
    ang_fluxes_sort = []
    for j in range(len(theta)):
        for i,a in enumerate(theta):
            if a == theta_min:
                theta_sort.append(a)
                ang_fluxes_sort.append(ang_fluxes[i])
                theta.remove(theta_min)
                ang_fluxes.remove(ang_fluxes[i])
                if len(theta) != 0:
                    theta_min = min(theta)

    # append 0 and 2 pi values
    ang_flux_max = max(ang_fluxes_sort)
    ang_flux_0 = ang_fluxes_sort[0]
    theta_0 = - theta_sort[0]
    ang_fluxes_sort.insert(0, ang_flux_0)
    theta_sort.insert(0, theta_0)
    ang_flux_2pi = ang_fluxes_sort[-1]
    theta_2pi = 2*pi + (2*pi - theta_sort[-1])
    ang_fluxes_sort.append(ang_flux_2pi)
    theta_sort.append(theta_2pi)

    # draw slices for each theta
    for i, angle in enumerate(theta_sort[1:-1]):
        i = i+1

        radius = int(ang_fluxes_sort[i]/ang_flux_max*200)
        ang_start = int((theta_sort[i]+theta_sort[i-1])/2*180/pi)
        ang_end   = int((theta_sort[i+1]+theta_sort[i])/2*180/pi)

        #               print 'plotting slice radius ' + str(radius) + ' theta_min ' + str(ang_start) + ' theta_max ' + str(ang_end)

        ang_start = 360 - ang_start
        ang_end = 360 - ang_end

        draw.pieslice((250 - radius, 250 - radius, 250 + radius, 250 + radius), ang_end, ang_start, fill=(0,0,0), outline=(255,255,255))

    fig = plt.figure()
    plt.plot(theta_sort, ang_fluxes_sort)

    fig.savefig(savepath + '/ang_flux_' + location +'.png')
    img.save(savepath + '/ang_flux_circle_' + location + '.png')


def mkdir_p(mypath):
    '''Creates a directory. equivalent to using mkdir -p on the command line'''

    try:
        makedirs(mypath)
    except OSError as exc: # Python >2.5
        if exc.errno == EEXIST and path.isdir(mypath):
            pass
        else: raise

def saveInputFile(savepath):
    """saves copy of input file in same directory as plots"""
    filename = savepath + '/input.py'
    try:
        shutil.copy('problem.py', filename)
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)