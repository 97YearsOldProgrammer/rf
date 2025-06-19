import random
import numpy as np

def base2_to_int(bits) -> int:

    idx = 0
    for bit in bits:
        idx = (idx << 1) | bit
    return idx

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
    
    def __init__(self, dons, accs, pos2info, mtry, min_samples_split, max_depth=None, outlier=False, gini_threshold=False):
        self.dons               = dons
        self.accs               = accs
        self.pos2info           = pos2info
        self.mtry               = mtry
        self.max_depth          = max_depth
        self.gini_threshold     = gini_threshold
        self.min_samples_split  = min_samples_split

        self.subset = []
        self.output = dict()
        self.rules  = []
        
        self._btstrapping()
        self.dataset = self.dons + self.accs
        
    def _btstrapping(self):
        ''' boostrapping original donor and acceptor sites '''

        self.dons = sorted(random.choices(self.dons, k=len(self.dons) ))
        self.accs = sorted(random.choices(self.accs, k=len(self.accs) ))
        
    def _fsubset(self, node):
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

    def _split(self, node, max_gini_test=False):
        ''' 
            compute the mean squared error gain for split criteria
            additionally gini impurity is maintained for diversity of split quality
        '''

        subset     = self._fsubset(node)
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
            if gain_mse not in all_psplit:  all_psplit[gain_mse] = [(left_node, right_node)]
            else:                           all_psplit[gain_mse].append((left_node, right_node))
        
        max_gain = sorted([mse for mse in all_psplit], reverse=True)

        iter = 0
        for mse in max_gain:

            if max_gini_test and iter > max_gini_test:  return False

            gini_threshold = self.gini_threshold if self.gini_threshold else 0.1

            if len(all_psplit[mse]) == 2:
                left_split, right_split = all_psplit[mse]
                if len(left_split) < self.min_samples_split or len(right_split) < self.min_samples_split:continue
                left_gini  = self.compute_gini([typ for _, typ, _ in left_split])
                right_gini = self.compute_gini([typ for _, typ, _ in right_split])

                if left_gini <= gini_threshold and right_gini <= gini_threshold: continue
                return left_split, right_split

            else:
                split_ginis = []

                for left_split, right_split in all_psplit[mse]:
                    if len(left_split) < self.min_samples_split or len(right_split) < self.min_samples_split:continue
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

    def _store_output(self, node, path):
        ''' store the leaf node of decision tree into base2 keys'''

        id = base2_to_int(path)
        self.output[id] = node
    
    def _recursion_tree(self, node, depth, path):
        ''' recursion until death '''
        '''
            path:   list with 0 or 1 s.t. store the leaf node as decisions
                    instead data structure approach
        '''

        ''' save leaf node ; end recursion '''
        if self.max_depth:
            if depth == self.max_depth:
                self._store_output(node, path)
                return
        
        if len(node) < self.min_samples_split:
            self._store_output(node, path)
            return

        ''' continue splitting '''
        split = self._split(node)

        if split is False:
            self._store_output(node, path)
            return
        
        left_node, right_node = split
        right_node      = sorted(right_node, key=lambda x: x[2])
        _, _, val       = right_node[0]
        split_threshold = val
        self.rules.append(split_threshold)
        
        self._recursion_tree(left_node,  depth+1, [0]+path)
        self._recursion_tree(right_node, depth+1, [1]+path)
        pass

############################################################
############################################################
##################### DevData Generator ####################
############################################################
############################################################

def generate_dev_data(
        total_samples   = 20,
        don_ratio       = 0.5,
        value_range     = (0.0, 1.0),
        pos_range       = (100, 1000),
        seed            = None):

    if seed is not None:
        random.seed(seed)

    positions = random.sample(range(pos_range[0], pos_range[1] + 1), total_samples)
    random.shuffle(positions)

    num_dons = int(total_samples * don_ratio)
    dons     = positions[:num_dons]
    accs     = positions[num_dons:]

    pos2info = {
        p: (random.uniform(value_range[0], value_range[1]),
            'don' if p in dons else 'acc')
        for p in positions
    }

    return dons, accs, pos2info