def next_one(array, pos, len):
    ''' find nearest element in an array '''
    
    for i, num in enumerate(array):
        if  num <= pos: continue
        if  pos - num < len: continue
        del array[:i]   
        return num
    
    return None
          
def fw_rnatree(dons, accs, min_intron, min_exon):
    ''' forward RNA tree for random forest '''

    isoform = []
    def tree(pos):
        ''' greedy decision tree '''
        
        doner = next_one(dons, pos, min_exon   + 1)
        accer = next_one(accs, pos, min_intron - 1)
        
        if doner == None or accer == None: return
        
        if  doner > accer: 
            isoform.append(doner)
            pos = doner
        else:
            isoform.append(accer)
            pos = accer
        
        tree(pos)
    
    tree(min_exon)
    return isoform
        
def bw_rnatree(dons, accs, min_intron, min_exon, pexon, pintron):
    pass