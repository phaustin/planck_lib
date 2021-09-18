# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     notebook_metadata_filter: all,-language_info,-toc,-latex_envs
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.3.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
#   latex_metadata:
#     chead: center head
#     lhead: left head
# ---

# %%
import numpy as np
from planck import planckwavelen, planckwavenum, clight
from matplotlib import pyplot as plt
import matplotlib
from itertools import cycle

c1 = 3.74e-16
c2 = 1.44e-2

__version__ = "$Id: plotplanck.py,v 1.1 2006/02/21 19:18:06 phil Exp phil $"


# %%
#this trick will run  the following script if
#the file planck.py is run as a program, but won't
#if  planck.py is imported from another  module

if __name__ == '__main__':

    #
    # fig1
    #
    fig1, ax1 = plt.subplots(1, 1)
    #check to see how good this approx. is:
    wavel = np.arange(1, 20., 0.1) * 1.e-6
    TempList = [
        260.,
        270,
        280,
        290,
        300,
    ]
    bList = []
    for Temp in TempList:
        bbr = planckwavelen(wavel, Temp)
        bList.append(bbr)
    for bbr in bList:
        plt.plot(wavel * 1.e6, bbr * np.pi * 1.e-6)
    ax1.grid(True)
    ax1.set_ylabel('$M_\lambda\ {(W\,m^{-2}\,\mu m^{-1})}$')
    ax1.legend(('260 K', '270 K', '280 K', '290 K', '300 K'))
    ax1.set_title(
        'Planck function $M_\lambda\ {(W\,m^{-2}\,\mu m)}$ for 5 blackbody temperatures'
    )
    ax1.set_xlabel('wavelength $\lambda\ {\mu m}$')
    fig1.savefig('planckII.png')


# %%
if __name__ == '__main__':
    #
    # fig 2
    #
    fig2, ax2 = plt.subplots(1,1)
    themarks = cycle(['-+', '-,', '-o', '-.', '->'])
    theLines = []
    for bbr in bList:
        theLines.append(
            ax2.plot(wavel * 1.e6, bbr * 1.e-6, next(themarks), markersize=7))
    ax2.plot(wavel * 1.e6, bbr * 1.e-6)
    ax2.set_ylabel('$L_{bb}\ \mathrm{(W\,m^{-2}\,\mu m^{-1}\,sr^{-1})}$')
    ax2.legend(('260 K', '270 K', '280 K', '290 K', '300 K'))
    ax2.set_title('Planck function $L_{bb}$ for 5 blackbody temperatures')
    ax2.set_xlabel('wavelength $\lambda\ {\mu m}$')
    fig2.savefig('planckLtrun.pdf')


# %%
if __name__ == '__main__':
#
# fig 3
#
    fig3, ax3 = plt.subplots(1,1)
    TempList = [
        260.,
        270,
        280,
        290,
        300,
    ]
    wavenum = np.arange(25, 2500, 20)  #in inverse cm
    wavenum_m = wavenum*100.  #in inverse meters
    bList = []
    for Temp in TempList:
        #
        # convert from cm-1 to m-1
        #
        bbr = planckwavenum(wavenum_m, Temp)*100.
        bList.append(bbr)
    for bbr in bList:
        ax3.plot(wavenum, bbr* np.pi)  #Watts/m^2/cm^-1
    ax3.set_ylabel('$B_n\ {(W\,(m^{-2}\,cm^{-1})}$')
    ax3.grid(b='off', linewidth=1, linestyle='-', which='both')
    ax3.grid(b=False)
    ax3.figure.canvas.draw()
    xminorLocator = matplotlib.ticker.MaxNLocator(nbins=5)
    ax3.xaxis.set_minor_locator(xminorLocator)
    yminorLocator = matplotlib.ticker.MaxNLocator(nbins=10)
    ax3.yaxis.set_minor_locator(yminorLocator)
    ax3.legend(('260 K', '270 K', '280 K', '290 K', '300 K'))
    ax3.set_title(
        'Planck function $B_n\ {(W\,/(m^{2}\,cm^{-1})}$ for 5 blackbody temperatures'
    )
    ax3.set_xlabel('wavenumber $n\ {cm^{-1}}$')
    fig3.savefig('q3_retrievalIIb.png', dpi=150)


# %%
if __name__ == '__main__':
    Temp = 300.
    #
    # fig 4
    #
    fig4, ax4 = plt.subplots(1,1,figsize=(12,10))
    bbr = planckwavenum(wavenum_m, 300.)
    first_fill = ax4.fill_between(wavenum, 0,
                                  bbr * np.pi* 100.)  #Watts/m^2/cm^-1
    xlab = ax4.set_xlabel('$wavenumber\ (cm^{-1})$')
    xlab.set_fontsize(16)
    theTextTitle = ax4.text(1800,
                            0.4,
                            '300 K blackbody flux',
                             ha='center')
    theTextTitle.set_fontsize(25)
    ylab = ax4.set_ylabel('$I\ (W/(m^{2}\,cm^{-1})$')
    ylab.set_fontsize(16)
    theText = ax4.text(1500, 0.2, '$\sigma T^4 = 460\ W\,m^{-2}$')
    theText.set_fontsize(20)
    fig4.savefig('fig41.png', dpi=300)


# %%
    sb = 5.67e-8 * 300**4.
    level = sb / 2500.
    ax4.plot([0, 2500], [level, level], 'r-')
    theFill = ax4.fill_between(wavenum, 0, level, facecolor='red')
    theFill.set_alpha(0.2)
    fig4.savefig('fig42.png', dpi=300)
    display(fig4)


# %%
    ax2 = fig4.add_axes(ax4.get_position(), frameon=False)
    ax2.set_ylim(ax4.get_ylim())
    ax2.set_xlim(ax4.get_xlim())
    topticks = np.asarray([1, 2, 3, 4, 5, 7, 10., 12, 16, 25])
    tickLabels = ["%d" % i for i in topticks]
    topticks = topticks * 1.e-4
    ticloc = 1. / topticks
    ax4.xaxis.set_ticks_position('bottom')
    ax2.xaxis.set_ticks(ticloc)
    ax2.xaxis.set_ticklabels(tickLabels)
    ax2.xaxis.tick_top()
    ax2.yaxis.tick_right()
    ax2.set_xlim(ax4.get_xlim())
    ax2.xaxis.set_label_position('top')
    xlab2 = ax2.set_xlabel('$wavelength\ (\mu m)$')
    xlab2.set_fontsize(16)
    display(fig4)
    fig4.savefig('fig43.png', dpi=300)


# %%
    ax4.fill_between(wavenum,
                     0,
                     bbr * np.pi*100.,
                     where=np.logical_and(wavenum > 600, wavenum < 800),
                     facecolor='green')
    display(fig4)
    fig4.savefig('fig44.png', dpi=300)


# %%
