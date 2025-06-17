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
    
    def __init__(self, dons, accs, pos2info, mtry, max_depth=None, outlier=False, gini_threshold=False):
        self.dons           = dons
        self.accs           = accs
        self.pos2info       = pos2info
        self.mtry           = mtry
        self.max_depth      = max_depth
        self.gini_threshold = gini_threshold
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

        return random.sample(node, self.mtry)
    
    def compute_mse(self, vals):
        ''' compute the mean squared error '''

        vals  = np.asarray(vals, dtype=float)
        return np.mean( (vals - vals.mean())**2 ) if vals.size else 0.0
    
    def compute_gini(self, typs):
        ''' compute the gini impurity '''

        count = dict()
        total = len(typs)

        for typ in typs:
            if typ not in count: count[typ]  = 1
            else:                count[typ] += 1
        
        gini = 1.0
        for frq in count.values():
            p       = frq/total
            gini   -= p*p

        return gini

    def split_mse(self, node, max_gini_test=False):
        ''' 
            compute the mean squared error gain for split criteria
            additionally gini impurity is maintained for diversity of split quality
        '''

        subset     = self.f_subset(node)
        parent_mse = self.compute_mse([val for _, _, val in subset])
        subset     = sorted(subset, key=lambda x: x[2])

        all_psplit = dict()
        
        for i in range(len(subset)):
            gain_mse   = parent_mse
            left_node  = subset[:i+1]
            right_node = subset[i+1:]
            left_mse   = self.compute_mse([val for _, _, val in left_node])
            right_mse  = self.compute_mse([val for _, _, val in right_node])
            gain_mse   = gain_mse - left_mse - right_mse
            if gain_mse not in all_psplit:  all_psplit[gain_mse] = [tuple(left_node) + tuple(right_node)]
            else:                           all_psplit[gain_mse].append(tuple(left_node)+tuple(right_node))
        
        max_gain = sorted([mse for mse in all_psplit], reverse=True)

        iter = 0
        for mse in max_gain:

            if max_gini_test and iter > max_gini_test:  return False

            gini_threshold = self.gini_threshold if self.gini_threshold else 0.1

            if len(all_psplit[mse]) == 2:
                left_split, right_split = all_psplit[mse]
                left_gini  = self.compute_gini([typ for _, typ, _ in left_split])
                right_gini = self.compute_gini([typ for _, typ, _ in right_split])

                if left_gini <= gini_threshold and right_gini <= gini_threshold: continue
                return left_split, right_split

            else:
                split_ginis = []

                for left_split, right_split in all_psplit[mse]:
                    left_gini  = self.compute_gini([typ for _, typ, _ in left_split])
                    right_gini = self.compute_gini([typ for _, typ, _ in right_split])
                    split_ginis.append(left_gini + right_gini)

                if max(split_ginis) <= 0.1: return False

                max_idx    = split_ginis.index(max(split_ginis))

                if max_idx % 2 == 0:
                    left_split  = all_psplit[mse][max_idx]
                    right_split = all_psplit[mse][max_idx+1]
                    return left_split, right_split
                else:
                    left_split  = all_psplit[mse][max_idx-1]
                    right_split = all_psplit[mse][max_idx]
                    return left_split, right_split

            iter += 1

    def _recursion_tree(self, dataset, depth):
        ''' recursion until death '''
        if self.depth:
            if self.depth == self.max_depth: return
            
        pass
        
        