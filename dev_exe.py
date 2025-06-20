import randomf

from randomf import IsoformTree

dons, accs, pos2info = randomf.generate_dev_data(60)
tree = IsoformTree(dons, accs, pos2info)
op   = tree.output
print(op)