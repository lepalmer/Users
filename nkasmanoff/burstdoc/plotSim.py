#!/usr/bin/env python

try:
    import matplotlib.pyplot as plt

    # Set the default title font dict
    titleFormat = {'fontsize': 12,
                   'fontweight': plt.rcParams['axes.titleweight'],
                   'verticalalignment': 'baseline',
                   'horizontalalignment': 'center'}

except ImportError:
    print("\n**** Warning: matplotlib not found. " +
          "Do not try to make plots or bad things will happen! ****")
    exit()

    
def plotAeffvsEnergy(energy, aeff, aeff_eres, aeff_eres_modfrac, theta=0):

    plt.figure(figsize=(8, 6))
    plt.title(r'Effective Area vs. Energy ($\theta$ = {:,.0f}$^\circ$)'
              .format(theta))
    plt.scatter(energy, aeff, color='black')
    plt.plot(energy, aeff, color='black', alpha=0.5, linestyle='--',
             lw=2, label='BurstCube')
    plt.scatter(energy, aeff_eres, color='blue')
    plt.plot(energy, aeff_eres, color='blue', alpha=0.5, linestyle='--',
             lw=2, label='BurstCube with E$_{\mathrm{res}}$')
    plt.scatter(energy, aeff_eres_modfrac, color='red')
    plt.plot(energy, aeff_eres_modfrac, color='red', alpha=0.5, linestyle='--',
             lw=2, label='BurstCube with E$_{\mathrm{res}}$ + escape')

    plt.xscale('log')
    plt.xlabel('Energy (keV)', fontsize=16)
    plt.yscale('log')
    plt.ylabel('Effective Area (cm$^2$)', fontsize=16)
    plt.legend(loc='lower center', prop={'size': 16}, numpoints=1,
               frameon=False)

    
def plotAeffvsTheta(theta, aeff, aeff_eres, aeff_eres_modfrac, energy=100.):

    plt.figure(figsize=(8, 6))
    plt.title(r'Effective Area vs. Angle (E = {:,.0f} keV)'
              .format(energy))
    plt.scatter(theta, aeff, color='black')
    plt.plot(theta, aeff, color='black', alpha=0.5, linestyle='--',
             lw=2, label='BurstCube')
    plt.scatter(theta, aeff_eres, color='blue')
    plt.plot(theta, aeff_eres, color='blue', alpha=0.5, linestyle='--',
             lw=2, label='BurstCube with E$_{\mathrm{res}}$')
    plt.scatter(theta, aeff_eres_modfrac, color='red')
    plt.plot(theta, aeff_eres_modfrac, color='red', alpha=0.5, linestyle='--',
             lw=2, label='BurstCube with E$_{\mathrm{res}}$ + escape')

    plt.xlabel('Incident Angle (deg)', fontsize=16)
    plt.ylabel('Effective Area (cm$^2$)', fontsize=16)
    plt.legend(loc='lower center', scatterpoints=1, prop={'size': 16},
               frameon=False)
    plt.grid(True)

    
def plotAeff(simFiles):

    aeffs = simFiles.calculateAeff()

    for angle in set(aeffs['theta']):
        mask = aeffs['theta'] == angle
        plotAeffvsEnergy(aeffs['keV'][mask],
                         aeffs['aeff'][mask],
                         aeffs['aeff_eres'][mask],
                         aeffs['aeff_eres_modfrac'][mask],
                         angle)
    plt.show()

    for energy in set(aeffs['keV']):
        mask = aeffs['keV'] == energy
        plotAeffvsTheta(aeffs['theta'][mask],
                        aeffs['aeff'][mask],
                        aeffs['aeff_eres'][mask],
                        aeffs['aeff_eres_modfrac'][mask],
                        energy)
        plt.grid(True)

    plt.show()
