from sketch_for_max_whole import sketch

def remove_supersets(sets):
    unique_sets = list( set( frozenset(s) for s in sets ) )
    result = list()
    for s in unique_sets:
        if not any( s > other for other in unique_sets if s != other ):
            result.append(s)
    return [ set(s) for s in result ]

def minimalize(GC, L, U, k, I_family):
    
    # remove dominated constraints
    I_family = remove_supersets(I_family)
    
    # set priority weights used to order the vertices when inspecting them for deletion
    weight = { i : 0 for i in GC.nodes }
    gtv = { GC.nodes[i]['GEOID20'] : i for i in GC.nodes }
    for I in I_family:
        for g in I:
            weight[gtv[g]] += (1/len(I))

    # minimalize
    minimal_I_family = list()
    for I in I_family:

        #print("Trying to minimalize I =", I)
        if len(I)==1:
            minimal_I_family.append([ i for i in I])
            continue
            
        IV = [ gtv[g] for g in I ] # I vertices
        WV = IV.copy()     # W vertices
        items = [(i,weight[i]) for i in IV ]
        for (i,iweight) in sorted(items, key=lambda item: item[1], reverse=True):
            #print(f"attempting to drop {i} with weight {iweight}")
            WV.remove(i)
            soln = sketch(GC, L, U, k, whole_counties=WV, draw=False, time_limit=10, statewide=False, verbose=False)
            if soln != 'infeasible':
                WV.append(i)
        W = [ GC.nodes[i]['GEOID20'] for i in WV ] # W geoids
        if len(W) < len(I):
            print(f"Dropped {len(I)-len(W)} and reduced from {I} to {W}.")
        minimal_I_family.append(W)

    return minimal_I_family          
