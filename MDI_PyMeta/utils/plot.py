import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class AnimatedPlot():
    def __init__(self, kcalmol_to_atomic, angstrom_to_atomic):
        self.kb = 0.0019872041
        self.T = 1.0
        self.kcalmol_to_atomic = kcalmol_to_atomic
        self.angstrom_to_atomic = angstrom_to_atomic

        self.fig = plt.figure()

        #creating a subplot
        self.ax1 = self.fig.add_subplot(1,1,1)

        self.x_min = 1.0 * angstrom_to_atomic
        self.x_max = 20.0 * angstrom_to_atomic
        self.x = np.linspace(self.x_min, self.x_max, 200)

        # The most recent nsmooth Gaussians will be scaled in a running average
        self.nsmooth = 2500


    def show(self, s_of_t, width, height, bias, dgrid):

        #f_x = np.zeros(len(self.x))
        #n = len(s_of_t)

        #for i in range(len(s_of_t)):
        #    scaling = 1.0
        #    for ix in range(len(self.x)):
        #        x = self.x[ix]
        #        f_x[ix] -= scaling * height * np.exp(-(x - s_of_t[i])**2/(2*width**2))
        #f_x -= f_x.min()

        min_bias = min(bias)
        xbias = [ (i * dgrid ) / self.angstrom_to_atomic for i in range(len(bias)) ]
        ybias = [ ( bias[i] - min_bias ) / self.kcalmol_to_atomic for i in range(len(bias)) ]


        
        epsilon = 0.2
        sigma = 3.88
        f_ref = 4*epsilon * ((sigma/self.x)**12 - (sigma/self.x)**6) + 2*self.kb*self.T*np.log(self.x)
        #f_ref = 2*self.kb*self.T*np.log(self.x)
        f_ref -= f_ref.min()

        
        self.ax1.clear()
        #self.ax1.plot( self.x / self.angstrom_to_atomic, f_x / self.kcalmol_to_atomic,
        #               xbias, ybias )
        self.ax1.plot( xbias, ybias, self.x, f_ref )

        self.ax1.set_ylim([0.0, 0.5])

        plt.show(block = False)
        plt.pause(0.0001)

    def finalize(self):
        plt.show()
