from . import collectivevariable as cv
from . import utils

class Distance(cv.CollectiveVariable):
    def __init__(self, atomi, atomj):
        super(cv.CollectiveVariable, self).__init__()

        self.natoms_ = 2
        self.atomi_ = atomi - 1
        self.atomj_ = atomj - 1
        self.gradi_ = [ 0.0 for i in range(3) ]
        self.gradj_ = [ 0.0 for i in range(3) ]
    
    def GetValue(self):
        return self.value_

    def GetGradient(self):
        return [ self.gradi_, self.gradj_ ]

    def GetAtoms(self):
        return [ self.atomi_, self.atomj_ ]

    def Evaluate(self, xyz, natoms, box_len):
        atomi_xyz = [ xyz[ 3 * self.atomi_ + i ] for i in range(3) ]
        atomj_xyz = [ xyz[ 3 * self.atomj_ + i ] for i in range(3) ]
        rij = [ atomi_xyz[i] - atomj_xyz[i] for i in range(3) ]

        rij = utils.Minimum_Image(rij, box_len)

        self.value_ = utils.Norm(rij)

        self.gradi_ = [ rij[i] / self.value_ for i in range(3) ]
        self.gradj_ = [ - rij[i] / self.value_ for i in range(3) ]

        
