import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyfftw
import random
import time
from astropy.io import fits


start_time = time.time()

def main():
    N = 101  # number of spatial lattice sites (keep odd)
    t = 0  # initial time
    delta = 10 ** -50  # small number to handle the limiting form of ~ sin(x)/x

    # Domain [0,L] x [0,L] x [0,L]
    L = 25
    dx = L / N
    xlin = np.linspace(0, L, num=N + 1)  # array in any spatial direction
    xlin = xlin[0:N]  # to erase the last point, as it is the same as the first point in our periodic box
    xx, yy, zz = np.meshgrid(xlin, xlin, xlin, sparse=True)
    eta = 24  # number of points to sample a full 2 pi rotation
    fullrot = 2 * np.pi * (dx) ** 2 * (1 / 3)  # one full rotation of the laplacian phase
    dt = (1 / eta) * fullrot  # timestep
    # print(dt)
    tEnd = 20  # time at which simulation ends
    Nt = int(np.ceil(tEnd / dt))  # array of counts

    plotrealtime = True  # switch on for plotting as the simulation goes along
    tdraw = 10 * dt  # draw frequency for plotting in real time
    plotcount = 1  # initial count for plotting in real time

    Lambda = -0.01  # self interaction strength; positive means attractive, negative means repulsive

    # initial field configuration; modify as desired
    sigma = 4.0 * dx
    psi1 = np.exp(1.j * random.uniform(3 * dx, L - 3 * dx) * 2 * np.pi) * np.exp(
        -((xx - random.uniform(3 * dx, L - 3 * dx)) ** 2 + (yy - random.uniform(3 * dx, L - 3 * dx)) ** 2 + (
                    zz - random.uniform(3 * dx, L - 3 * dx)) ** 2) / (
                2 * sigma ** 2))
    psi2 = np.exp(1.j * random.uniform(3 * dx, L - 3 * dx) * 2 * np.pi) * np.exp(
        -((xx - random.uniform(3 * dx, L - 3 * dx)) ** 2 + (yy - random.uniform(3 * dx, L - 3 * dx)) ** 2 + (
                    zz - random.uniform(3 * dx, L - 3 * dx)) ** 2) / (
                2 * sigma ** 2))
    psi3 = np.exp(1.j * random.uniform(3 * dx, L - 3 * dx) * 2 * np.pi) * np.exp(
        -((xx - random.uniform(3 * dx, L - 3 * dx)) ** 2 + (yy - random.uniform(3 * dx, L - 3 * dx)) ** 2 + (
                    zz - random.uniform(3 * dx, L - 3 * dx)) ** 2) / (
                2 * sigma ** 2))

    rho1 = np.abs(psi1) ** 2
    rho2 = np.abs(psi2) ** 2
    rho3 = np.abs(psi3) ** 2

    # creating k-space and associated nabla and laplacian matrix
    nvalues = np.arange(-(N - 1) / 2, (N + 1) / 2)
    klin = 2.0 * N / L * np.sin(np.pi * nvalues / N)
    kx, ky, kz = np.meshgrid(klin, klin, klin, sparse=True)
    kx = pyfftw.interfaces.numpy_fft.ifftshift(kx)
    ky = pyfftw.interfaces.numpy_fft.ifftshift(ky)
    kz = pyfftw.interfaces.numpy_fft.ifftshift(kz)
    ksq = kx ** 2 + ky ** 2 + kz ** 2

    # #  This is relevant for Schrodinger current/momentum
    # klin2 = N / L * np.sin(2 * np.pi * nvalues / N)
    # kx2, ky2, kz2 = np.meshgrid(klin2, klin2, klin2, sparse=True)
    # kx2 = pyfftw.interfaces.numpy_fft.ifftshift(kx2)
    # ky2 = pyfftw.interfaces.numpy_fft.ifftshift(ky2)
    # kz2 = pyfftw.interfaces.numpy_fft.ifftshift(kz2)

    # prep figure
    fig = plt.figure(figsize=(8, 6), dpi=80)
    grid = plt.GridSpec(2, 3, wspace=0.5, hspace=0.2)
    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1])
    ax3 = plt.subplot(grid[1, 0])
    ax4 = plt.subplot(grid[1, 1])
    # ax5 = plt.subplot(grid[0, 2])
    # ax6 = plt.subplot(grid[1, 2])

    # creating arrays for fractional change in conserved quantities
    num = np.sum(rho1 + rho2 + rho3) * dx ** 3
    numfrac = []

    spinx = (1.j * np.sum(psi2 * np.conjugate(psi3)) - 1.j * np.sum(psi3 * np.conjugate(psi2))) * dx ** 3
    spiny = (-1.j * np.sum(psi1 * np.conjugate(psi3)) + 1.j * np.sum(psi3 * np.conjugate(psi1))) * dx ** 3
    spinz = (1.j * np.sum(psi1 * np.conjugate(psi2)) - 1.j * np.sum(psi2 * np.conjugate(psi1))) * dx ** 3
    spinfrac = []

    rtime = []

    # energy = []
    # momentumx = []
    # momentumy = []
    # momentumz = []

    for i in range(Nt):

        # performing half drift
        psi1k = pyfftw.interfaces.numpy_fft.fftn(psi1)
        psi2k = pyfftw.interfaces.numpy_fft.fftn(psi2)
        psi3k = pyfftw.interfaces.numpy_fft.fftn(psi3)
        psi1k = np.exp(dt * (-1.j * ksq / 4.0)) * psi1k
        psi2k = np.exp(dt * (-1.j * ksq / 4.0)) * psi2k
        psi3k = np.exp(dt * (-1.j * ksq / 4.0)) * psi3k
        psi1 = pyfftw.interfaces.numpy_fft.ifftn(psi1k)
        psi2 = pyfftw.interfaces.numpy_fft.ifftn(psi2k)
        psi3 = pyfftw.interfaces.numpy_fft.ifftn(psi3k)

        # Calculating quantities for the full kick evolution
        rho = np.abs(psi1) ** 2 + np.abs(psi2) ** 2 + np.abs(psi3) ** 2
        phik = -pyfftw.interfaces.numpy_fft.fftn(0.5 * (rho - np.mean(rho))) / (
                ksq + (ksq == 0))
        phi = np.real(pyfftw.interfaces.numpy_fft.ifftn(phik))

        fsq = psi1 ** 2 + psi2 ** 2 + psi3 ** 2
        magfsq = np.abs(fsq)
        magfsq = np.where(magfsq > rho, rho, magfsq)  # getting rid of possible error due to machine precision
        phase = (np.angle(fsq) + 2.0 * np.pi) % (2.0 * np.pi)
        spin = np.sqrt(rho ** 2 - magfsq ** 2)

        U11 = np.cos(dt * Lambda * rho) * np.cos(dt * Lambda * spin) + dt * Lambda * np.sinc(
            dt * Lambda * spin / np.pi) * (rho * np.sin(dt * Lambda * rho) - magfsq * np.sin(phase + dt * Lambda * rho))
        U12 = -np.sin(dt * Lambda * rho) * np.cos(dt * Lambda * spin) + dt * Lambda * np.sinc(
            dt * Lambda * spin / np.pi) * (rho * np.cos(dt * Lambda * rho) + magfsq * np.cos(phase + dt * Lambda * rho))
        U21 = np.sin(dt * Lambda * rho) * np.cos(dt * Lambda * spin) + dt * Lambda * np.sinc(
            dt * Lambda * spin / np.pi) * (
                      -rho * np.cos(dt * Lambda * rho) + magfsq * np.cos(phase + dt * Lambda * rho))
        U22 = np.cos(dt * Lambda * rho) * np.cos(dt * Lambda * spin) + dt * Lambda * np.sinc(
            dt * Lambda * spin / np.pi) * (rho * np.sin(dt * Lambda * rho) + magfsq * np.sin(phase + dt * Lambda * rho))

        # Performing full kick
        psi1 = np.exp(-1.j * dt * (phi - 2.0 * Lambda * rho)) * (
                    (U11 + 1.j * U21) * (np.real(psi1)) + (U12 + 1.j * U22) * np.imag(
                psi1))
        psi2 = np.exp(-1.j * dt * (phi - 2.0 * Lambda * rho)) * (
                    (U11 + 1.j * U21) * (np.real(psi2)) + (U12 + 1.j * U22) * np.imag(
                psi2))
        psi3 = np.exp(-1.j * dt * (phi - 2.0 * Lambda * rho)) * (
                    (U11 + 1.j * U21) * (np.real(psi3)) + (U12 + 1.j * U22) * np.imag(
                psi3))


        #  Performing half drift
        psi1k = pyfftw.interfaces.numpy_fft.fftn(psi1)
        psi2k = pyfftw.interfaces.numpy_fft.fftn(psi2)
        psi3k = pyfftw.interfaces.numpy_fft.fftn(psi3)
        psi1k = np.exp(dt * (-1.j * ksq / 4.)) * psi1k
        psi2k = np.exp(dt * (-1.j * ksq / 4.)) * psi2k
        psi3k = np.exp(dt * (-1.j * ksq / 4.)) * psi3k
        psi1 = pyfftw.interfaces.numpy_fft.ifftn(psi1k)
        psi2 = pyfftw.interfaces.numpy_fft.ifftn(psi2k)
        psi3 = pyfftw.interfaces.numpy_fft.ifftn(psi3k)

        # #  Calculating energy
        # Egrav = 0.5 * np.sum(phi * rho) * dx ** 3
        # Eself = -0.5 * Lambda * np.sum(3 * rho ** 2 - spin ** 2) * dx ** 3
        # Ekin = 0.5 * np.sum(ksq * np.abs(psi1k) ** 2 + ksq * np.abs(psi2k) ** 2 + ksq * np.abs(psi3k) ** 2) * dx**3/N**3
        # Ekin = np.abs(Ekin)
        # Etot = Egrav + Eself + Ekin
        #
        # #  Calculating momentum
        # px = np.sum(kx2 * (np.abs(psi1k) ** 2 + np.abs(psi2k) ** 2 + np.abs(psi3k) ** 2)) * dx ** 3 / N ** 3
        # py = np.sum(ky2 * (np.abs(psi1k) ** 2 + np.abs(psi2k) ** 2 + np.abs(psi3k) ** 2)) * dx ** 3 / N ** 3
        # pz = np.sum(kz2 * (np.abs(psi1k) ** 2 + np.abs(psi2k) ** 2 + np.abs(psi3k) ** 2)) * dx ** 3 / N ** 3

        t += dt
        rtime.append(t)
        numfrac.append(np.abs(np.sum(rho) * dx ** 3 - num)/num)
        spinfracxinst = np.abs(
            (1.j * np.sum(psi2 * np.conjugate(psi3)) - 1.j * np.sum(psi3 * np.conjugate(psi2))) * dx ** 3 - spinx) / np.abs(spinx)
        spinfracyinst = np.abs(
            (-1.j * np.sum(psi1 * np.conjugate(psi3)) + 1.j * np.sum(psi3 * np.conjugate(psi1))) * dx ** 3 - spiny) / np.abs(spiny)
        spinfraczinst = np.abs(
            (1.j * np.sum(psi1 * np.conjugate(psi2)) - 1.j * np.sum(psi2 * np.conjugate(psi1))) * dx ** 3 - spinz) / np.abs(spinz)
        spinfracinst = (spinfracxinst+spinfracyinst+spinfraczinst)/3
        spinfrac.append(spinfracinst)
        # energy.append(Etot)
        # momentumx.append(px)
        # momentumy.append(py)
        # momentumz.append(pz)

        # # check CFL condition
        # if dt > 2*np.pi*eta*np.minimum(np.abs(1/phi)) or dt > 2*np.pi*eta*np.minimum(np.abs(1/(2*Lambda*rho))):
        #     print('Caution, high dense regions are appearing')

        # plot in real time
        plotthisturn = False
        if t + dt > plotcount * tdraw:
            plotthisturn = True
        if (plotrealtime and plotthisturn) or (i == Nt - 1):
            plt.sca(ax1)
            plt.cla()
            plt.imshow(np.log10(np.sum(rho, axis=2)/N), cmap='afmhot')
            plt.clim(-4, 1)
            plt.gca().invert_yaxis()
            ax1.get_xaxis().set_visible(False)
            ax1.get_yaxis().set_visible(False)
            ax1.set_aspect('equal')

            plt.sca(ax2)
            plt.cla()
            plt.imshow(np.log10(np.sum(spin, axis=2) / N), cmap='magma')
            plt.clim(-4, 1)
            plt.gca().invert_yaxis()
            ax2.get_xaxis().set_visible(False)
            ax2.get_yaxis().set_visible(False)
            ax2.set_aspect('equal')

            plt.sca(ax3)
            plt.cla()
            plt.plot(rtime, numfrac)
            plt.yscale("log")
            plt.xscale("linear")

            plt.sca(ax4)
            plt.cla()
            plt.plot(rtime, spinfrac)
            plt.yscale("log")
            plt.xscale("linear")

            # plt.sca(ax5)
            # plt.cla()
            # plt.plot(rtime, energy)
            #
            # plt.sca(ax6)
            # plt.cla()
            # plt.plot(rtime, momentumx)
            # plt.plot(rtime, momentumy)
            # plt.plot(rtime, momentumz)

            plt.pause(0.00001)
            plotcount += 1


    # Save figure
    plt.sca(ax1)
    plt.title(r'rho')
    plt.sca(ax2)
    plt.title(r'spin')
    plt.sca(ax3)
    plt.title(r'Delta N')
    plt.sca(ax4)
    plt.title(r'Delta S')
    # plt.sca(ax5)
    # plt.title(r'E_tot at 2999')
    # plt.sca(ax6)
    # plt.title(r'P at 2999')

    return 0

if __name__ == "__main__":
    main()