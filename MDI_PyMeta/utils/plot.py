import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class AnimatedPlot():
    def __init__(self, kcalmol_to_atomic, angstrom_to_atomic):
        self.kb = 0.0019872041
        self.T = 130.0
        self.kcalmol_to_atomic = kcalmol_to_atomic
        self.angstrom_to_atomic = angstrom_to_atomic

        self.fig = plt.figure()

        #creating a subplot
        self.ax1 = self.fig.add_subplot(1,1,1)

        self.x_min = 1.0 * angstrom_to_atomic
        self.x_max = 14.0 * angstrom_to_atomic
        self.x = np.linspace(self.x_min, self.x_max, 100)

        # The most recent nsmooth Gaussians will be scaled in a running average
        self.nsmooth = 2500


    def show(self, s_of_t, width, height):


        f_x = np.zeros(len(self.x))
        n = len(s_of_t)

        #gauss = lambda x, n: np.sum(height * np.exp(-(x - s_of_t[:n])**2/(2*width**2)))
        #f_x = np.array([-gauss(xx,n) for xx in self.x])

        for i in range(len(s_of_t)):
            scaling = min( float(len(s_of_t) - i) / float(self.nsmooth), 1.0)
            for ix in range(len(self.x)):
                x = self.x[ix]
                f_x[ix] -= scaling * height * np.exp(-(x - s_of_t[i])**2/(2*width**2))
        f_x -= f_x.min()

        self.ax1.clear()
        self.ax1.plot( self.x / self.angstrom_to_atomic, f_x / self.kcalmol_to_atomic )

        plt.show(block = False)
        plt.pause(0.0001)

    def finalize(self):
        plt.show()
