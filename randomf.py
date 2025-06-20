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
##################### DevData Generator ####################
############################################################

def generate_dev_data(
        total_samples   = 20,
        don_ratio       = 0.6,
        value_range     = (0.0, 1.0),
        pos_range       = (0, 1000),
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

############################################################
##################### DecisionTree Class ###################
############################################################


'''
random forest parameter:
    mtry    : p/3
    nodesize: 30% of root node
'''
class IsoformTree:

    def __init__(
        self,
        dons,
        accs,
        pos2info,
    ):
        
        self.dons               = dons
        self.accs               = accs
        self.pos2info           = pos2info
        self.min_samples_split  = int(len(self.dons+self.accs) * 0.3)

        self.output = dict()
        self.rules = []

        self._btstrapping()
        self._recursion_tree(self.dons+self.accs, [])

    def _btstrapping(self):
        ''' Bootstrap original donor and acceptor sites '''
        self.dons = sorted(random.choices(self.dons, k=len(self.dons)))
        self.accs = sorted(random.choices(self.accs, k=len(self.accs)))

    def _fsubset(self, node, mtry):
        ''' Find random subset of current dataset '''
        return random.sample(node, mtry)

    def compute_mse(self, vals):
        ''' Compute the mean squared error '''
        vals = np.asarray(vals, dtype=float)

        return np.mean((vals - vals.mean())**2) if vals.size else 0.0

    def compute_gini(self, typs):
        ''' Compute the gini impurity '''
        count = {}
        total = len(typs)

        for typ in typs:
            count[typ] = count.get(typ, 0) + 1

        gini = 1.0
        for frq in count.values():
            p = frq / total
            gini -= p * p

        return gini

    def _gini_test(self, typs):
        ''' Compute the gini impurity test '''
        threshold   = self.gini_threshold if self.gini_threshold else 0.1
        gini        = self.compute_gini(typs)

        return gini < threshold

    def _find_threshold(self, right_split):
        ''' Find the split threshold for split original dataset '''
        right_split     = sorted(right_split, key=lambda x: x[2])
        _, _, threshold = right_split[0]

        self.rules.append(threshold)

    def _split(self, node):
        ''' 
            Compute the mean squared error gain for split criteria
            Additionally gini impurity is maintained for diversity of split quality
        '''
        subset      = self._fsubset(node, int(len(node)/3))
        parent_mse  = self.compute_mse([val for _, _, val in subset])
        subset      = sorted(subset, key=lambda x: x[2])

        all_psplit = {}

        for i in range(len(subset)):
            gain_mse    = parent_mse
            left_node   = subset[:i+1]
            right_node  = subset[i+1:]
            left_mse    = self.compute_mse([val for _, _, val in left_node])
            right_mse   = self.compute_mse([val for _, _, val in right_node])
            gain_mse    = gain_mse - left_mse - right_mse

            all_psplit.setdefault(gain_mse, []).append((left_node, right_node))

        max_gain = sorted(all_psplit.keys(), reverse=True)

        for mse in max_gain:
            splits = all_psplit[mse]

            if len(splits) == 1:
                left_split, right_split = splits[0]

                ltest = len(left_split)  < self.min_samples_split
                rtest = len(right_split) < self.min_samples_split
                if ltest and rtest: continue

                ltest = self._gini_test([typ for _, typ, _ in left_split])
                rtest = self._gini_test([typ for _, typ, _ in right_split])
                if ltest and rtest:  continue

                self._find_threshold(right_split)
                return

            else:
                ''' if mse collision happen '''
                split_ginis = []

                for left_split, right_split in splits:
                    ltest = len(left_split)  < self.min_samples_split
                    rtest = len(right_split) < self.min_samples_split
                    if ltest and rtest: continue

                    left_gini  = self.compute_gini([typ for _, typ, _ in left_split])
                    right_gini = self.compute_gini([typ for _, typ, _ in right_split])
                    split_ginis.append(left_gini)
                    split_ginis.append(right_gini)

                if not split_ginis or max(split_ginis) <= 0.1:  return False

                max_idx = split_ginis.index(max(split_ginis))

                if max_idx % 2 == 0:
                    right_split = splits[max_idx+1]
                    self._find_threshold(right_split)
                    return
                
                else:
                    right_split = splits[max_idx]
                    self._find_threshold(right_split)
                    return

    def _store_output(self, node, path):
        ''' Store the leaf node of decision tree into base2 keys '''
        id = base2_to_int(path)
        self.output[id] = node

    def _recursion_tree(self, node, path):
        ''' Recursion until death '''

        if len(node) < self.min_samples_split:
            self._store_output(node, path)
            return

        split = self._split(node)
        if split is False:
            self._store_output(node, path)
            return

        split_threshold = self.rules[-1]
        node            = sorted(node, key=lambda x: x[2])

        idx = 0
        for i, (*_, val) in enumerate(node):
            if val > split_threshold:
                idx = i
                break
        
        left_node  = node[:idx]
        right_node = node[idx+1:]

        self._recursion_tree(left_node,  [0] + path)
        self._recursion_tree(right_node, [1] + path)