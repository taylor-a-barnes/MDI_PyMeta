import math

def Crossp(a, b):
    cross = [ 0.0 for i in range(3) ]
    cross[0] = a[1]*b[2] - a[2]*b[1]
    cross[1] = - (a[0]*b[2] - a[2]*b[0])
    cross[2] = a[0]*b[1] - a[1]*b[0];
    return cross

def Minimum_Image(xyz, box_length):
    minimage = [ 0.0 for i in range(3) ]
    for i in range(3):
        minimage[i] = xyz[i] - float( round(xyz[i] / box_length[i]) ) * box_length[i]
    return minimage

def Norm(array):
    size = 0.0
    for element in array:
        size += element * element
    size = math.sqrt(size)
    return size

def Gaussian(ds, width, height):
    arg = ( ds * ds ) / (2.0 * width * width )
    return height * math.exp(-arg)

def Gaussian_derv(ds, width, height):
    arg = ( ds * ds ) / (2.0 * width * width )
    pre = - ds / ( width * width )
    return height * pre * math.exp(-arg)
