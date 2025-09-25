import time
import networkx as nx
from math import floor, ceil
from sketch_for_max_whole import sketch
from minimalize import minimalize 
from wcd_finder import wcd
from cluster_for_max_whole import is_within_window

def get_initial_constraints(GC, L, U, k, cover_threshold=5, vi=True, gvi=True, ai=True, gai=True, analysis=False, verbose=False):

    start = time.perf_counter()
    assert GC.number_of_nodes() <= 254, "Expecting county-level graph in add_initial_constraints()."

    # directed version of GC with population-weighted edges
    HC = nx.DiGraph(GC)
    for i,j in HC.edges:
        HC.edges[i,j]['pweight'] = HC.nodes[i]['TOTPOP']

    wc_districts = list()   # list of whole-county districts
    num_districts = 0
    cover = { i : list() for i in GC.nodes } # for each county, which wcds cover it?
    constraints = list() 

    for c in GC.nodes:
        # draw some whole-county districts that cover county c?
        # if covered less than cover_threshold, then keep searching for more
        while len(cover[c]) < cover_threshold:
            district = wcd(HC, L, U, k, c, wc_districts)
            if district is None:
                break
            for i in district:
                cover[i].append( num_districts )
            wc_districts.append( district )
            num_districts += 1

    # vicinity inequalities
    overpopulated_counties = [ i for i in GC.nodes if GC.nodes[i]['TOTPOP'] > U ]
    if vi:
        for c in GC.nodes:
            if len(cover[c]) == 0:
                dist = nx.shortest_path_length( HC, source=c, weight='pweight' )
                vicinity = [ i for i in dist.keys() if dist[i] <= U ]
                if c in overpopulated_counties or [ i for i in vicinity if i in overpopulated_counties ] == []:
                    #m.addConstr( sum( m._w[i] for i in vicinity ) <= len(vicinity) - 1 )
                    cf = GC.nodes[c]['GEOID20']
                    vf = { GC.nodes[i]['GEOID20'] for i in vicinity }
                    #printif(verbose, f"Adding vicinity inequality for c = {cf} and V_c = {vf}.")
                    constraints.append( vf )
    
    # generalized vicinity inequalities
    if gvi:
        for c1 in GC.nodes:
            if len(cover[c1])==0 or len(cover[c1])>=cover_threshold:
                continue
            for c2 in GC.nodes:
                if c1>=c2 or len(cover[c2])==0 or len(cover[c2])>=cover_threshold:
                    continue

                if all( is_incompatible(GC, L, U, k, wc_districts[cover[c1][i]], wc_districts[cover[c2][j]]) for i in range(len(cover[c1])) for j in range(len(cover[c2])) ):
                    dist1 = nx.shortest_path_length( HC, source=c1, weight='pweight' )
                    dist2 = nx.shortest_path_length( HC, source=c2, weight='pweight' )
                    vicinity = [ i for i in HC.nodes if ( i in dist1.keys() and dist1[i] <= U ) or ( i in dist2.keys() and dist2[i] <= U ) ]
                    if [ i for i in vicinity if i in overpopulated_counties ] == []:
                        #m.addConstr( sum( m._w[i] for i in vicinity ) <= len(vicinity) - 1 )
                        Cpf = {GC.nodes[c1]['GEOID20'], GC.nodes[c2]['GEOID20']}
                        vf = { GC.nodes[i]['GEOID20'] for i in vicinity }
                        #printif(verbose, f"Adding generalized vicinity inequality for C' = {Cpf} and V_c = {vf}.")
                        constraints.append( vf )
                
    # articulation inequalities
    if ai:
        added = { c : False for c in GC.nodes }
        for c in nx.articulation_points(GC):
            Sc = [ i for i in GC.nodes if i != c ]
            for comp in nx.connected_components( GC.subgraph(Sc) ):
                pop = sum( GC.nodes[i]['TOTPOP'] for i in comp )
                q = ceil( (pop+1) / L ) # smallest q satisfying pop < L*q
                if U*q - GC.nodes[c]['TOTPOP'] < pop and not added[c]:
                    #m._w[c].UB = 0
                    added[c] = True
                    #printif(verbose, f"Added articulation ineq. for c = {GC.nodes[c]['GEOID20']}.")
                    constraints.append( { GC.nodes[c]['GEOID20'] } )

    # vertex cut inequalities (a.k.a. generalized articulation inequalities)
    if gai:
        conflict_graph = generate_conflict_graph(GC, L, U, k)
        Cp_complement = set( GC.nodes )
        for Cp in enumerate_small_cliques(conflict_graph):
            if not is_vertex_cut(GC, Cp):
                continue
            for i in Cp:
                Cp_complement.remove(i)
            for S in nx.connected_components( GC.subgraph( Cp_complement) ):
                B = { i for i in Cp for j in GC.neighbors(i) if j in S }
                CpUS = set( list(Cp)+list(S) )
                Bp = { i for i in B if set( GC.neighbors(i) ).issubset(CpUS) }
                pB = sum( GC.nodes[i]['TOTPOP'] for i in B )
                pS = sum( GC.nodes[i]['TOTPOP'] for i in S )
                pBp = sum( GC.nodes[i]['TOTPOP'] for i in Bp )
                if floor( (pS+pBp)/L ) - len(Bp) < max( 0, ceil( (pS+pB)/U ) - len(B) ):
                    Cpf = { GC.nodes[c]['GEOID20'] for c in Cp }
                    constraints.append( Cpf )
                    #printif(verbose, f"Adding vertex cut inequality for C' = {Cpf}.")
                    break
            for i in Cp:
                Cp_complement.add(i)
    
    if analysis:
        for c in GC.nodes:
            cov = len(cover[c])
            if 0 < cov and cov < cover_threshold:
                print("\nPossible whole-county districts for c =", GC.nodes[c]['GEOID20'], "are:")
                for j in range(len(cover[c])):
                    fps = [ GC.nodes[i]['GEOID20'] for i in wc_districts[cover[c][j]] ]
                    print("\nDistrict =",fps)
                    dist = nx.shortest_path_length( HC, source=c, weight='pweight' )
                    vicinity = [ i for i in dist.keys() if dist[i] <= U ]
                    print("Counties in vicinity:",[ GC.nodes[i]['GEOID20'] for i in vicinity ])
                    print("Overpopulated counties in vicinity:",[ GC.nodes[i]['GEOID20'] for i in vicinity if GC.nodes[i]['TOTPOP'] > U ])
                print("")

    # remove the dominated constraints and minimalize the remaning ones
    minimal_constraints = minimalize(GC, L, U, k, constraints)
    
    print("Adding non-dominated, minimal constraints:")
    for constraint in minimal_constraints:
        print(constraint)
    
    elapsed = time.perf_counter() - start
    print(f"Time to generate initial inequalities: {round(elapsed, 2)} seconds.")
    return minimal_constraints

def is_incompatible(GC, L, U, k, district1, district2):
    return not is_compatible(GC, L, U, k, district1, district2)

def is_compatible(GC, L, U, k, district1, district2):
    assert k >= 1
    if k == 1:
        assert set(district1) == set(district2)
        assert set(district1) == set(GC.nodes)
        return True

    # do districts overlap, while being distinct?
    if [ i for i in district1 if i in district2 ] != [] and set(district1) != set(district2):
        return False

    # does removing these two districts lead to components with "bad" population?
    nondistrict = [ i for i in GC.nodes if i not in district1+district2 ]
    for component in nx.connected_components( GC.subgraph( nondistrict ) ):
        population = sum( GC.nodes[i]['TOTPOP'] for i in component )
        if not is_within_window(population, L, U):
            return False
            
    return True

def is_vertex_cut(GC, S):
    S_complement = [ i for i in GC.nodes if i not in S ]
    return not nx.is_connected( GC.subgraph(S_complement) )

def generate_conflict_graph(GC, L, U, k):
    conflict_graph = nx.Graph()
    conflict_graph.add_nodes_from( GC.nodes )
    articulation = set( nx.articulation_points(GC) )
    under_nonart = { i for i in GC.nodes if GC.nodes[i]['TOTPOP'] <= U and i not in articulation }

    # add trivial conflicts
    trivial_conflicts = [ (i,j) for i in under_nonart for j in under_nonart if i<j and GC.nodes[i]['TOTPOP'] + GC.nodes[j]['TOTPOP'] > U ]
    conflict_graph.add_edges_from( trivial_conflicts )

    # add nontrival conflicts
    other_nodes = set( GC.nodes )
    for i in under_nonart:
        other_nodes.remove(i)

        for j in nx.articulation_points( GC.subgraph(other_nodes) ):
            if i < j and j in under_nonart:
                other_nodes.remove(j)
                for S in nx.connected_components( GC.subgraph(other_nodes) ):
                    pS = sum( GC.nodes[v]['TOTPOP'] for v in S )
                    pc = GC.nodes[i]['TOTPOP'] + GC.nodes[j]['TOTPOP']
                    if floor( pS/L ) < ceil( (pS+pc)/U ) - 1:
                        conflict_graph.add_edge(i,j)
                        #print("Nontrivial conflict!", GC.nodes[i]['NAME20'], GC.nodes[j]['NAME20'])
                        break
                other_nodes.add(j)
        other_nodes.add(i)
            
    return conflict_graph      

def enumerate_small_cliques(GC):
    
    for i in GC.nodes: # cliques of size 1
        yield [i]
    
    for i,j in GC.edges: # cliques of size 2
        yield [i,j]
    
    for i,j in GC.edges: # cliques of size 3
        for v in GC.neighbors(i):
            if v < min(i,j) and v in GC.neighbors(j):
                yield [v,i,j]

    # for i,j in GC.edges: # cliques of size 4
    #     CN = [ v for v in GC.neighbors(i) if v < min(i,j) and v in GC.neighbors(j) ]
    #     for u,v in GC.subgraph(CN).edges:
    #         yield [u,v,i,j]
