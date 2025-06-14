import random

############################################################
############################################################
##################### DecisionTree Class ###################
############################################################
############################################################

class IsoformTree:
    def __init__(self, dons, accs, dict, mtry, max_depth=None):
        self.dons = dons
        self.accs = accs
        self.dict = dict
        self.mtry = mtry
        self.max_depth = max_depth
    
    def btstrapping(dons, accs):
        ''' boostrapping original donor and acceptor sites '''
    
def btstrapping(dons, accs, size=None):
    ''' boostrapping original donor and acceptor sites '''
    
    if size is None:    size = len(dons)
    dons = sorted(random.choices(dons, k=size))
    
    if size is None:    size = len(accs)
    accs = sorted(random.choices(accs, k=size))
    
    return dons, accs
        