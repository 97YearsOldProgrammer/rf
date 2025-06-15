import random
import numpy as np

def finding_outlier(doa, pos2info, z_threshold=None):
    '''
    doa     :   dons or accs
    pos2info:   dictionary for pos : val, typ
    '''
    assert(doa)

    mdoa = []

    for pos in doa:
        val, typ = pos2info[pos]
        mdoa.append((pos, typ, val))

    mdoa = sorted(mdoa, key=lambda x: x[2])
    vals  = np.asarray([val for _, _, val in mdoa], dtype=float)
    mean  = vals.mean()
    std   = vals.std(ddof=0)
    zscr  = (vals - mean)/std

    if not z_threshold: otest = np.abs(zscr) >= 3.0
    else:               otest = np.abs(zscr) >= z_threshold

    otl   = [info for info, tf in zip(mdoa, otest) if tf]
    doa   = [info for info, tf in zip(mdoa, otest) if not tf]

    return otl, doa

############################################################
############################################################
##################### DecisionTree Class ###################
############################################################
############################################################

class IsoformTree:
    
    def __init__(self, dons, accs, dict, mtry, max_depth=None, outlier=False):
        self.dons = dons
        self.accs = accs
        self.dict = dict
        self.mtry = mtry
        self.max_depth = max_depth
        self.subset = []
        self.output = []
        
        self._btstrapping()
        self.dataset = self.dons + self.accs
        
    def _btstrapping(self):
        ''' boostrapping original donor and acceptor sites '''

        self.dons = sorted(random.choices(self.dons, k=len(self.dons) ))
        self.accs = sorted(random.choices(self.accs, k=len(self.accs) ))
        
    def f_subset(self, node):
        ''' find random subset of current dataset '''

        subset = np.random.choice(node, size=self.mtry, replace=False)
        return subset
    
    def compute_mse(self, array):
        ''' compute the mean squared error '''

        arr  = np.asarray(array, dtype=float)
        return np.mean( (arr - arr.mean())**2 ) if arr.size else 0.0
    
    def split_mse(self, node):
        ''' compute the mean squared error gain for split criteria '''

        subset     = self.f_subset(node)
        parent_mse = self.compute_mse(subset)
        
        all_psplit = dict()
        
        for i in range(len(subset)):
            gain_mse   = parent_mse
            left_node  = subset[:i+1]
            right_node = subset[i+1:]
            left_mse   = self.compute_mse([val for _, val in left_node])
            right_mse  = self.compute_mse([val for _, val in right_node])
            gain_mse   = gain_mse - left_mse - right_mse
            all_psplit[(left_node, right_node)] = gain_mse
            
    def _recursion_tree(self, dataset, depth):
        ''' recursion until reach terminal node or depth'''
        if self.depth:
            if self.depth == self.max_depth: return
            
        pass
        
        