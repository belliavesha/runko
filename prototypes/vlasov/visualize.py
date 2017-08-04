import numpy as np
from pylab import *
from matplotlib import cm
import os, sys


def plot_phasespace(ax, xxi, vxi, ffi, kk):

    #pick kth species
    xx = xxi
    vx = vxi[:,kk]
    ff = ffi[:,:,kk]

    # xx vx[:, kk], ff[:,:,kk] ff[vx:, x:, species]
    ax.set_xlim(xx[0], xx[-1])
    ax.set_ylim(vx[0], vx[-1])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$v_{x}$')
    
    X, Y = np.meshgrid(xx, vx)
    ax.pcolormesh(X, Y, ff, 
            cmap='Reds', 
            vmin=ff.min(),
            vmax=ff.max(),
            )



def plot_mean_velocity_pdf(ax, vx, ff, kk):
    ax.cla()

    ax.set_xlim(-10, 10)
    #ax.set_ylim( )

    ax.set_xlabel(r'$v_{x}$')
    #ax.set_ylabel(r'pdf')

    fv = np.mean(ff[:,3:-3, kk], 1)
    ax.plot(vx, fv, "k-")



def plot_field(ax, xx, f, quantity):
    ax.cla()

    ax.set_xlim(xx[0], xx[-1])
    #ax.set_ylim(vx[0], vx[-1])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(quantity)

    ax.plot(xx, f, "b-")



class visualize:

    initialized = False

    def __init__(self, path, xx, vx):

        self.path = path
        if not os.path.exists(path):
            os.makedirs(path)

        self.xx = xx
        self.vx = vx


        self.gs = GridSpec(3, 3)
        self.gs.update(hspace = 0.2)
        self.gs.update(wspace = 0.4)

        self.ax1a = subplot(self.gs[0, 0:2])
        self.ax1b = subplot(self.gs[0, 2])

        self.ax2a = subplot(self.gs[1, 0:2])
        self.ax2b = subplot(self.gs[1, 2])

        self.ax3a = subplot(self.gs[2, 0])
        self.ax3b = subplot(self.gs[2, 1])
        self.ax3c = subplot(self.gs[2, 2])



    def plot(self, step, ff, ex, ajx, rho):

        sstep = str(step).rjust(4, '0')
        fname = self.path+'/vlasov'+sstep+'.png'

        #phase spaces // species 0
        plot_phasespace(self.ax1a, self.xx, self.vx, ff, 0)
        plot_mean_velocity_pdf(self.ax1b, self.vx, ff, 0)
            
        #phase spaces // species 1
        plot_phasespace(self.ax2a, self.xx, self.vx, ff, 1)
        plot_mean_velocity_pdf(self.ax2b, self.vx, ff, 1)
            
        #fields
        plot_field(self.ax3a, self.xx, ex, r'$E_x$')
        plot_field(self.ax3b, self.xx, rho, r'$\rho_q$')
        plot_field(self.ax3c, self.xx, ajx, r'$J_x$')
        

        savefig(fname)

