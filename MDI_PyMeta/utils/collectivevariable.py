class CollectiveVariable(object):
    def __init__(self):
        self.value_ = 0.0
        self.natoms_ = 0

    def Evaluate(xyz, natoms, box_len):
        raise NotImplementedError()

    def GetValue():
        raise NotImplementedError()

    def GetGradient():
        raise NotImplementedError()

    def GetAtoms():
        raise NotImplementedError()
