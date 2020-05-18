import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class AnimatedPlot():
    def __init__(self):
        self.kb = 0.0019872041
        self.T = 0.0
        self.kcalmol_to_atomic = 0.00015936014383657207 / 0.1
        self.angstrom_to_atomic = 0.37794522509156564 / 0.2
        self.kcalmol_to_kjmol = 4.184

        self.fig = plt.figure()

        #creating a subplot
        self.ax1 = self.fig.add_subplot(1,1,1)

        self.x_min = 1.0 * self.angstrom_to_atomic
        self.x_max = 20.0 * self.angstrom_to_atomic
        self.x_points = 200
        self.x = np.linspace(self.x_min, self.x_max, self.x_points)

        # The most recent nsmooth Gaussians will be scaled in a running average
        self.nsmooth = 30000


    def show(self):

        hills = np.loadtxt('work/s_of_t.out')
        
        f_x = np.zeros(len(self.x))
        n = len(hills)
        print("N: " + str(n))

        for i in range(len(hills)):
            s_of_t = hills[i,1]
            width = hills[i,2]
            height = hills[i,3]
            scaling = min( float(len(hills) - i) / float(self.nsmooth), 1.0)

            # Try scaling the Gaussians based on the volume of the nearby shell
            #scaling *= 1.0 / s_of_t
            
            for ix in range(len(self.x)):
                x = self.x[ix]
                #print("s, w, h: " + str(s_of_t) + " " + str(width) + " " + str(height))
                f_x[ix] -= scaling * height * np.exp(-(x - s_of_t)**2/(2*width**2))
        max_index = round( 10.0 * self.angstrom_to_atomic / ( (self.x_max - self.x_min) / self.x_points ) )
        print("MMM: " + str(max_index))
        f_x -= f_x[:max_index].min()

        #min_bias = min(bias)
        #xbias = [ (i * dgrid ) / self.angstrom_to_atomic for i in range(len(bias)) ]
        #ybias = [ ( bias[i] - min_bias ) / self.kcalmol_to_atomic for i in range(len(bias)) ]

        epsilon = 0.2
        sigma = 3.88
        f_ref = 4*epsilon * ((sigma/self.x)**12 - (sigma/self.x)**6) + 2*self.kb*self.T*np.log(self.x)
        f_ref -= f_ref.min()

        self.ax1.clear()
        self.ax1.plot( self.x / self.angstrom_to_atomic, f_x / self.kcalmol_to_atomic,
                       self.x, f_ref)
        #self.ax1.plot( xbias, ybias )
        self.ax1.set_ylim([0.0, 0.5])

        #plt.show(block = False)
        #plt.pause(0.0001)
        plt.show()

    def finalize(self):
        plt.show()

my_plot = AnimatedPlot()
my_plot.show()

