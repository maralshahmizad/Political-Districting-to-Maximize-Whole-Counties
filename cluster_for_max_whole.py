from itertools import combinations
import matplotlib.pyplot as plt
import time
import math
import random
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from sketch_for_max_whole import sketch

from separation_for_max_whole import labeling_contiguity_callback



def sort_by_second(val):
    return val[1]

def find_ordering(G):
    V_with_population = [ (i, G.nodes[i]['TOTPOP']) for i in G.nodes ]
    V_with_population.sort( key=sort_by_second, reverse=True )
    return [ v for (v,p) in V_with_population ]

def construct_position(ordering):
    position = { i : -1 for i in ordering }
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position

# key:    ( tuple(sorted(G.nodes)), k, num_clusters )
# value:  (clusters, sizes, cut_edge_LB, cut_edge_UB)
min_cut_cache = dict()
        
def min_cut_county_clustering_via_labeling(G, L, U, k, num_clusters, time_limit=3600, carve=False, initial_clusters=None, initial_sizes=None, whole_counties=list(), randomized_carving=False, cache=True, verbose=False, sketch_intensity=0):
    
    # check if we've already solved this instance before.
    # if so, return our cached solution!
    key = ( tuple(sorted(G.nodes)), k, num_clusters )
    if cache and key in min_cut_cache:
        #print("Pulling min_cut solution from cache.")
        return min_cut_cache[key]
    
    # when carving, we only seek 2 clusters!
    assert num_clusters == 2 or carve == False
    
    DG = nx.DiGraph(G) # bidirected version of G
    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0
    m.Params.TimeLimit = time_limit

    # VARIABLES
    # x[i,j]=1 if vertex i is assigned to cluster j
    x = m.addVars(DG.nodes, num_clusters, vtype=GRB.BINARY) 

    # y[j] = # of districts in cluster j
    y = m.addVars(num_clusters, vtype=GRB.INTEGER, lb=1)
    
    # z[u,v,j] = 1 if directed edge (u,v) is cut because u->j but v!->j
    z = m.addVars(DG.edges, num_clusters, vtype=GRB.BINARY)
    
    # is_cut[u,v] = 1 if undirected edge {u,v} is cut
    is_cut = m.addVars(G.edges, vtype=GRB.BINARY)
    
    # BASE CONSTRAINTS
    # Each county i assigned to one cluster
    m.addConstrs( gp.quicksum( x[i,j] for j in range(num_clusters) ) == 1 for i in DG.nodes )

    # Cluster "sizes" should sum to k
    m.addConstr( y.sum() == k )

    # Population balance: population of cluster j should be an integer multiple of [L,U] 
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) <= U * y[j] for j in range(num_clusters) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) >= L * y[j] for j in range(num_clusters) )
    
    # CUT EDGE CONSTRAINTS
    m.addConstrs( x[u,j] - x[v,j] <= z[u,v,j] for u,v in DG.edges for j in range(num_clusters) )
    m.addConstrs( is_cut[u,v] == gp.quicksum( z[u,v,j] for j in range(num_clusters) ) for u,v in G.edges )
    m.addConstrs( is_cut[u,v] == gp.quicksum( z[v,u,j] for j in range(num_clusters) ) for u,v in G.edges )
                
    # SYMMETRY BREAKING, via diagonal fixing
    # Find vertex ordering
    ordering = find_ordering(G)

    for p in range(G.number_of_nodes()):
        i = ordering[p]
        for j in range(p+1,num_clusters):
            x[i,j].UB = 0
            
    # MIP START
    if initial_clusters is not None:

        assert len(initial_clusters) == num_clusters
        assert len(initial_sizes) == num_clusters
        
        root_found = [ False for j in range(num_clusters) ]
        labeling = { i : j for j in range(num_clusters) for i in initial_clusters[j] } 
        p = 0
        for i in ordering:
            j = labeling[i]
            if not root_found[j]:
                root_found[j] = True
                y[p].start = initial_sizes[j]
                for v in initial_clusters[j]:
                    x[v,p].start = 1
                p += 1

    # CONTIGUITY CONSTRAINTS
    m.Params.LazyConstraints = 1
    m._DG = DG
    m._parts = num_clusters
    m._x = x
    m._y = y
    m._k = k
    m._L = L
    m._U = U
    m._G = G
    m._whole_counties = whole_counties
    m._sketch_intensity = sketch_intensity
    m._callback = labeling_contiguity_callback
    if carve:
        # *temporarily* use large MIP gap, just to determine feasibility
        m.Params.MIPGap = 1.00
    
    # MINIMIZE CUT EDGES
    m.setObjective( is_cut.sum(), GRB.MINIMIZE )
    m.Params.IntFeasTol = 1.e-7
    m.Params.FeasibilityTol = 1.e-7
    m.optimize( m._callback )

    # NO SOLUTION FOUND? EXIT.
    if m.solCount == 0:
        clusters = [ list(G.nodes) ]
        sizes = [ k ]
        cut_edge_LB = 0
        cut_edge_UB = G.number_of_edges() + 1
        
        # add solution to our cache and return
        min_cut_cache[key] = (clusters, sizes, cut_edge_LB, cut_edge_UB)
        return min_cut_cache[key]
    
    if carve:
        m.Params.MIPGap = 0.00
        # give the complement solution as a warm start
        for j in range(2):
            y[j].start = k - round( y[j].x )
            for i in G.nodes:
                m._x[i,j].start = 1 - round( m._x[i,j].x )
        
        # undo diagonal fixing
        for p in range(G.number_of_nodes()):
            i = ordering[p]
            for j in range(p+1,num_clusters):
                x[i,j].UB = 1
        
        # try each possible size for first cluster, and return
        # the smallest feasible one, with fewest cut edges
        if randomized_carving:
            for i,j in G.edges:
                is_cut[i,j].obj = random.random() # random number between 0 and 1
        else:
            m.setObjective( is_cut.sum(), GRB.MINIMIZE )
        for size in range(1,k-1):
            y[0].LB = size
            y[0].UB = size
            m._callback = labeling_contiguity_callback
            m.optimize( m._callback )
            if m.solCount > 0:
                break
        
    # Retrieve solution
    clusters = [ [ i for i in G.nodes if m._x[i,j].x > 0.5 ] for j in range(num_clusters) ]
    sizes = [ round( y[j].x ) for j in range(num_clusters) ]
    cut_edge_LB = round( m.objBound )
    cut_edge_UB = round( m.objVal )

    # add solution to our cache and return
    if cache:
        min_cut_cache[key] = (clusters, sizes, cut_edge_LB, cut_edge_UB)
    return (clusters, sizes, cut_edge_LB, cut_edge_UB)


def is_within_window(population, L, U):
    min_size = math.ceil( population / U )
    max_size = math.floor( population / L )
    return min_size <= max_size

# Is it possible to build a county cluster 
#  that contains all nodes from 'required_nodes'?
#
def is_cluster_possible(G, required_nodes, L, U, k, verify=True):
    (cluster, size) = find_cluster(G, required_nodes, L, U, k, verify) 
    return cluster is not None
    
def find_cluster(G, required_nodes, L, U, k, verify=True):
    
    # To make contiguity easy, this function assumes that
    #   required_nodes is a connected dominating set.
    if verify:
        assert nx.is_connected( G.subgraph( required_nodes ) )
        assert nx.is_dominating_set( G, required_nodes )
    
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars( G.nodes, vtype=GRB.BINARY )
    for j in required_nodes:
        x[j].LB = 1
    
    z = m.addVar( vtype=GRB.INTEGER )
    z.LB = 1
    z.UB = k
    
    m.addConstr( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i] for i in G.nodes ) >= L * z )
    m.addConstr( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i] for i in G.nodes ) <= U * z )
    
    m.optimize()
    if m.solCount > 0:
        cluster = [ i for i in G.nodes if x[i].x > 0.5 ]
        size = round(z.x)
        return ( cluster, size ) 
    else:
        return (None, None)


# key:    ( tuple(sorted(G.nodes)), k)
# value:  (clusters, sizes, cluster_UB) 
max_cluster_cache = dict()


def number_of_cut_edges( G, clusters ):
    labeling = { i : j for j in range(len(clusters)) for i in clusters[j] }
    return sum( 1 for i,j in G.edges if labeling[i] != labeling[j] )


def carve_cluster(G, L, U, k, whole_counties, randomized_carving=False, sketch_intensity=0):
    (clusters, sizes, cut_edge_LB, cut_edge_UB) = min_cut_county_clustering_via_labeling(G, L, U, k, num_clusters=2, carve=True, whole_counties=whole_counties, randomized_carving=randomized_carving, cache=not randomized_carving, sketch_intensity=sketch_intensity)
    return (clusters[0], sizes[0])
    

def carve_heuristic(G, DG, L, U, k, whole_counties, randomized_carving=False, sketch_intensity=0):
    
    clusters = list()
    sizes = list()
    R = set(G.nodes)   # remaining vertices
    kr = k              # remaining size
    
    while kr > 0:
        (cluster, size) = carve_cluster( G.subgraph(R), L, U, kr, whole_counties, randomized_carving, sketch_intensity=sketch_intensity )
        if len(clusters)==0:
            print("carved cluster sizes = ",end="")
            
        clusters.append(cluster)
        sizes.append(size)
        R -= set( cluster )
        kr -= size
        print(size,end=', ')
        
    print("\ncarved LB =", len(clusters))
    print("carved cut edges =",number_of_cut_edges(G, clusters))
    return (clusters, sizes)

# converts our double objective into a single value
#   1st objective: max # clusters
#   2nd objective: min # cut edges
def objective(G, clusters):
    M = G.number_of_edges() + 1
    return M * len(clusters) - number_of_cut_edges(G, clusters)


# our callback
def vicinity_callback(m, where):
    
    # check if LP relaxation at this BB node is integer
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve data
    xval = m.cbGetSolution(m._x)
    k = m._k
    L = m._L
    U = m._U
    DG = m._DG
    
    # does each component of G-D have a population 
    #  within [L,U] or an integer multiple thereof?
    unselected = [ i for i in DG.nodes if xval[i] < 0.5 ]
    for component in nx.strongly_connected_components( DG.subgraph(unselected) ):
        population = sum( DG.nodes[i]['TOTPOP'] for i in component )
        multiple = False
        for q in range(1,k):
            if q*L <= population and population <= q*U:
                multiple = True
        if not multiple:
            # if not, then either:
            #   1. something from component needs to be picked (by x), or 
            #   2. one of its neighbors shouldn't be picked.
            neighbors = { j for i in component for j in DG.neighbors(i) if j not in component }
            m.cbLazy( gp.quicksum( m._x[i] for i in component ) + gp.quicksum( (1-m._x[i]) for i in neighbors ) >= 1 )

            component_names = [ DG.nodes[i]['NAME20'] for i in component ]
            neighbor_names = [ DG.nodes[i]['NAME20'] for i in neighbors ]
            #if m._verbose:
            #    print("adding cut for component =",component,"and neighbors =",neighbors)

            
def max_whole_vicinity_relaxation(G, DG, k, L, U, verbose=True):
    
    vicinity = dict()

    # collection of all whole-county districts found during our search
    districts = list() 
    num_districts = 0
    
    # which whole-county districts cover each given node?
    cover = { i : list() for i in DG.nodes }

    vic_start = time.time()

    for c in DG.nodes:

        # if covered less than twice, then keep searching for more
        while len(cover[c]) < 3:
            
            # is it possible to build a county-level feasible district containing county c?
            m = gp.Model()
            m.Params.OutputFlag = 0 # stop Gurobi from printing logs

            # Variables
            x = m.addVars( DG.nodes, vtype=GRB.BINARY )
            f = m.addVars( DG.edges )

            # root the district at c
            x[c].LB = 1

            # population between L and U
            m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i] for i in DG.nodes ) >= L )
            m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i] for i in DG.nodes ) <= U )

            # contiguity constraints: send one unit of flow to selected nodes
            m.addConstrs( gp.quicksum( f[j,i] - f[i,j] for j in DG.neighbors(i) ) == x[i] for i in DG.nodes if i != c )

            # contiguity constraints: accept inflow only at selected nodes, but not at c.
            M = DG.number_of_nodes() - 1
            m.addConstrs( gp.quicksum( f[j,i] for j in DG.neighbors(i) ) <= M * x[i] for i in DG.nodes if i != c )
            m.addConstr( gp.quicksum( f[j,c] for j in DG.neighbors(c) ) == 0 )

            # must differ from previous districts
            #if len(cover[c]) > 0:
            for j in range(len(cover[c])):
                district_number = cover[c][j]
                D = districts[district_number]
                VD = [ i for i in DG.nodes if i not in D ]
                m.addConstr( gp.quicksum( x[i] for i in VD ) + gp.quicksum( (1-x[i]) for i in D ) >= 1 )
            
            # callback info
            m._x = x
            m._L = L
            m._U = U
            m._DG = DG
            m._k = k
            m._verbose = verbose
            m.Params.LazyConstraints = 1
            m._callback = vicinity_callback
            
            # solve
            m.optimize(m._callback)

            # if not possible, what is V_c?
            if m.status == GRB.INFEASIBLE:
                break # exit the while-loop
            else:
                D = [ i for i in DG.nodes if m._x[i].x > 0.5 ]
                for i in D:
                    cover[i].append(num_districts)
                districts.append(D)
                num_districts += 1
                
    vic_end = time.time()
    vic_time = '{0:.2f}'.format(vic_end - vic_start)
    if verbose:
        print("Time to generate vicinity sets (in seconds) =",vic_time)
        print("cover =")
        for c in DG.nodes:
            if len(cover[c])>0 and len(cover[c])<3:
                print(c,cover[c])
    
    # articulation cuts
    articulation_cuts = list()
    for c in nx.articulation_points(G):
        Sc = [ i for i in G.nodes if i != c ]
        for component in nx.connected_components(G.subgraph(Sc)):
            population = sum( DG.nodes[i]['TOTPOP'] for i in component )
            for q in range(1,k+1):
                if q*U - DG.nodes[c]['TOTPOP'] <= population and population <= q*L:
                    articulation_cuts.append(c)
                    break
            if len(articulation_cuts)>0 and articulation_cuts[-1] == c:
                break
    
    # SOLVE VICINITY RELAXATION MODEL
    m = gp.Model()
    m.Params.OutputFlag = 0
    w = m.addVars(DG.nodes, vtype=GRB.BINARY)
    m.setObjective( w.sum(), GRB.MAXIMIZE)
    
    # add vicinity inequalities
    for c in DG.nodes:
        if len(cover[c]) != 0:
            continue
        dist = nx.shortest_path_length( DG, source = c, weight = 'pweight' )
        vicinity = [ i for i in dist.keys() if dist[i] < U ]
        vicinity_names = [ DG.nodes[i]['NAME20'] for i in vicinity ]
        if verbose:
            print("Adding vicinity inequality for c =",DG.nodes[c]['NAME20'],"county, where vicinity =",vicinity_names)
        m.addConstr( gp.quicksum( w[i] for i in vicinity ) <= len(vicinity) - 1 )
    
    # add pairwise vicinity inequalities
    if U < 500000: # skip for CD instances
        for c1 in DG.nodes:
            if len(cover[c1]) != 1:
                continue
            for c2 in DG.nodes:
                if c1>=c2 or len(cover[c2]) != 1:
                    continue
                D1 = set( districts[cover[c1][0]] )
                D2 = set( districts[cover[c2][0]] )
                if D1 == D2 or D1.isdisjoint(D2):
                    continue
                dist1 = nx.shortest_path_length( DG, source = c1, weight = 'pweight' )
                dist2 = nx.shortest_path_length( DG, source = c2, weight = 'pweight' )
                vicinity = [ i for i in DG.nodes if ( i in dist1.keys() and dist1[i] < U ) or ( i in dist2.keys() and dist2[i] < U ) ]
                vicinity_names = [ DG.nodes[i]['NAME20'] for i in vicinity ]
                if verbose:
                    print("Adding pairwise vicinity inequality for c1, c2 =",DG.nodes[c1]['NAME20'],DG.nodes[c2]['NAME20'],", where vicinity =",vicinity_names)
                m.addConstr( gp.quicksum( w[i] for i in vicinity ) <= len(vicinity) - 1 )

    # add articulation cuts
    for c in articulation_cuts:
        if verbose:
            print("adding articulation cut for",DG.nodes[c]['NAME20'])
        w[c].UB = 0
    
    m.optimize()
    num_whole = round(m.objVal)
    if verbose:
        print("# whole =",num_whole)
    
    # second solve to split the most populous counties (a hack to get "good" w vectors, that can be extended to full plans)
    # i.e., keep whole the least populous counties
    m.addConstr( w.sum() == num_whole )
    m.setObjective( gp.quicksum( DG.nodes[i]['TOTPOP'] * w[i] for i in DG.nodes ), GRB.MINIMIZE)
    m.optimize()
    
    # report solution
    split_counties = [ i for i in DG.nodes if w[i].x < 0.5 ]
    split_names = [ DG.nodes[i]['NAME20'] for i in split_counties ]
    if verbose:
        print("split_counties =",split_names)
    
    whole_counties = [ i for i in DG.nodes if w[i].x > 0.5 ]
    whole_names = [ DG.nodes[i]['NAME20'] for i in whole_counties ]
    if verbose:
        print("\nwhole_counties =",whole_names)
    
    # Find the N best solutions
    if verbose:
        print("How many optimal w vectors are there???")
        m.Params.OutputFlag = 1
        N = 100
        m.Params.PoolSearchMode = 2 # finds N best solutions
        m.Params.PoolSolutions = N
        m.setObjective( 0, GRB.MAXIMIZE)
        m.optimize()
        
        # retrieve optimal w*
        w_stars = list()
        for sol in range(m.SolCount):
            m.Params.SolutionNumber = sol
            w_star = [ i for i in G.nodes if w[i].Xn > 0.5 ]
            w_stars.append(w_star)
        print("All w* optima =",w_stars)
        print("There are",len(w_stars),"of them.")
        print("They are:")
        for ws in w_stars:
            print([ DG.nodes[i]['NAME20'] for i in ws])
    
    if U < 500000: # skip detailed analysis for CD instances
        for c in DG.nodes:
            if len(cover[c]) > 0: # no vicinity ineq added
                dist = nx.shortest_path_length( DG, source = c, weight = 'pweight' )
                vicinity = [ i for i in dist.keys() if dist[i] < U ]
                split_vicinity = [ i for i in vicinity if i in split_counties ]
                if split_vicinity == []:
                    vicinity_names = [ DG.nodes[i]['NAME20'] for i in vicinity ]
                    if verbose:
                        print("FYI. Nothing split in the vicinity of c =",DG.nodes[c]['NAME20'],"county, where vicinity =",vicinity_names)
                        print("Its cover is cover[c] =",cover[c])
                        for district_number in cover[c]:
                            district = districts[district_number]
                            nondistrict = [ i for i in DG.nodes if i not in district ]
                            (whole_counties2, split_counties2) = max_whole_vicinity_relaxation(G.subgraph(nondistrict), DG.subgraph(nondistrict), k-1, L, U, verbose=False)
                            print("if we force district =",district,"in our solution, then we can have total # whole counties at most:",len(district)+len(whole_counties2))
                
    return (whole_counties, split_counties)


def construct_position(ordering):
    position = { i : -1 for i in ordering }
    for p in range(len(ordering)):
        v = ordering[p]
        position[v] = p
    return position


def is_within_window(population, L, U):
    min_size = math.ceil( population / U )
    max_size = math.floor( population / L )
    return min_size <= max_size


def is_cluster_possible(G, required_nodes, L, U, k, verify=True):
    (cluster, size) = find_cluster(G, required_nodes, L, U, k, verify) 
    return cluster is not None
    
def find_cluster(G, required_nodes, L, U, k, verify=True):
    
    # To make contiguity easy, this function assumes that
    #   required_nodes is a connected dominating set.
    if verify:
        assert nx.is_connected( G.subgraph( required_nodes ) )
        assert nx.is_dominating_set( G, required_nodes )
    
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars( G.nodes, vtype=GRB.BINARY )
    for j in required_nodes:
        x[j].LB = 1
    
    z = m.addVar( vtype=GRB.INTEGER )
    z.LB = 1
    z.UB = k
    
    m.addConstr( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i] for i in G.nodes ) >= L * z )
    m.addConstr( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i] for i in G.nodes ) <= U * z )
    
    m.optimize()
    if m.solCount > 0:
        cluster = [ i for i in G.nodes if x[i].x > 0.5 ]
        size = round(z.x)
        return ( cluster, size ) 
    else:
        return (None, None)
    
    
def sort_by_second(val):
    return val[1]


def vertices_sorted_by_population(DG, vertices=list()):
    if vertices==list():
        vertices = [ i for i in DG.nodes ]
    nodes_with_population = [ ( i, DG.nodes[i]['TOTPOP'] ) for i in vertices ]
    nodes_with_population.sort(key=sort_by_second,reverse=True)
    return [v for (v,p) in nodes_with_population]


# our callback
def cut_callback(m, where):
    
    # check if LP relaxation at this BB node is integer
    if where != GRB.Callback.MIPSOL: 
        return

    # retrieve the LP relaxation solution at this BB node
    xval = m.cbGetSolution(m._x)
    yval = m.cbGetSolution(m._y)
    k = m._k
    L = m._L
    U = m._U
    G = m._G
    DG = m._DG

      
    split_roots = [ j for j in DG.nodes if xval[j,j] > 0.5 and j in m._split_counties ]
    
    for j in split_roots:
        
        # which counties are assigned to this split root?
        S = [ i for i in DG.nodes if xval[i,j] > 0.5 ]
        
        # check that there are no overpopulated appendages
        S_whole = [ i for i in m._whole_counties if xval[i,j] > 0.5 ]
        for component in nx.strongly_connected_components(DG.subgraph(S_whole)):
            population = sum( DG.nodes[v]['TOTPOP'] for v in component )
            if population > U:
                split_neighbors = [ u for v in component for u in DG.neighbors(v) if u in S and u in m._split_counties ]
                assert len(split_neighbors) > 0
                if len(split_neighbors) == 1:
                    split_neighbor = split_neighbors[0]
                    # if all of component is assigned to j, then need another one of its neighbors too
                    neighbors_not_sn = { u for v in component for u in DG.neighbors(v) if u not in S }
                    neighbors_not_sn = list(neighbors_not_sn)
                    m.cbLazy( gp.quicksum( (1-m._x[v,j]) for v in component ) + gp.quicksum( m._x[v,j] for v in neighbors_not_sn ) >= 1 )
                    component_names = [ DG.nodes[v]['NAME20'] for v in component ]
                    neighbor_names = [ DG.nodes[v]['NAME20'] for v in neighbors_not_sn ]
                    
            
        # check that each whole county i is "close" to a split county
        for i in S:
            if i in m._split_counties:
                continue
                
            dist = nx.shortest_path_length( DG.subgraph(S), source = i, weight = 'pweight' )
            nearby_split_counties = [ v for v in dist.keys() if dist[v] <= U and v in m._split_counties ]
            
            # i is close to a split county. We might be good.
            if len(nearby_split_counties) > 0:
                continue
                
            # otherwise, soln is definitely bad. Add a (minimal) cut.
            B = [ v for v in DG.nodes if xval[v,j] < 0.5 ] # "bad" nodes. We can write x[i,j] <= sum( x[v,j] for v in B )
            VB = [ v for v in S ] # V\B
            for v in B: # move v from B to VB?
                VB.append(v)
                dist = nx.shortest_path_length( DG.subgraph(VB), source = i, weight = 'pweight' )
                nearby_split_counties = [ v for v in dist.keys() if dist[v] <= U and v in m._split_counties ]
                if len(nearby_split_counties) > 0:
                    VB.pop() # remove v
            B = [ v for v in DG.nodes if v not in VB ]
            m.cbLazy( m._x[i,j] <= gp.quicksum( m._x[v,j] for v in B )  )
            B_names = [ DG.nodes[v]['NAME20'] for v in B ]
            
            
        # check that every split county is close to another split county (or is only split county in cluster)
        for i in S:
            if i in m._whole_counties or i == j:
                continue
                
            dist = nx.shortest_path_length( DG.subgraph(S), source = i, weight = 'pweight' )
            nearby_split_counties = [ v for v in dist.keys() if dist[v] <= U + DG.nodes[i]['TOTPOP'] and v in m._split_counties and v != i ]
            
            # i is close to a split county. We might be good.
            if len(nearby_split_counties) > 0:
                continue
                
            # otherwise, soln is definitely bad. Add a (minimal) cut.
            B = [ v for v in DG.nodes if xval[v,j] < 0.5 ] # "bad" nodes. We can write x[i,j] <= sum( x[v,j] for v in B )
            VB = [ v for v in S ] # V\B
            for v in B: # move v from B to VB?
                VB.append(v)
                dist = nx.shortest_path_length( DG.subgraph(VB), source = i, weight = 'pweight' )
                nearby_split_counties = [ v for v in dist.keys() if dist[v] <= U + DG.nodes[i]['TOTPOP'] and v in m._split_counties and v != i ]
                if len(nearby_split_counties) > 0:
                    VB.pop() # remove v
            B = [ v for v in DG.nodes if v not in VB ]
            m.cbLazy( m._x[i,j] <= gp.quicksum( m._x[v,j] for v in B )  )
            B_names = [ DG.nodes[v]['NAME20'] for v in B ]
            
            
        # check that no cluster has a pair of split counties, with an overpopulated "bridge" between them
        for component in nx.strongly_connected_components(DG.subgraph(S_whole)):
            population = sum( DG.nodes[v]['TOTPOP'] for v in component )
            if population > U:
                split_neighbors = [ u for v in component for u in DG.neighbors(v) if u in S and u in m._split_counties ]
                assert len(split_neighbors) > 0
                if len(split_neighbors) == 2 and split_neighbors[0] != split_neighbors[1]:
                    # if all of component is assigned to j, then need another one of its neighbors too
                    neighbors_not_sn = { u for v in component for u in DG.neighbors(v) if u not in S }
                    neighbors_not_sn = list(neighbors_not_sn)
                    m.cbLazy( gp.quicksum( (1-m._x[v,j]) for v in component ) + gp.quicksum( m._x[v,j] for v in neighbors_not_sn ) >= 1 )
                    component_names = [ DG.nodes[v]['NAME20'] for v in component ]
                    neighbor_names = [ DG.nodes[v]['NAME20'] for v in neighbors_not_sn ]
                    split_neighbor_names = [ DG.nodes[v]['NAME20'] for v in split_neighbors ]
                  
                    
        # check if any cut vertex c of the cluster neighbors a component G[S'] that has a bad population
        for c in nx.articulation_points(G.subgraph(S)):
            if c in m._split_counties:
                continue
                
            Sc = [ i for i in S if i != c ]
            for component in nx.connected_components(G.subgraph(Sc)):
                population = sum( DG.nodes[i]['TOTPOP'] for i in component )
                for q in range(1,k+1):
                    if q*U - DG.nodes[c]['TOTPOP'] <= population and population <= q*L:
                        neighbors_not_S = { u for v in component for u in DG.neighbors(v) if u not in S } 
                        m.cbLazy( gp.quicksum( (1-m._x[v,j]) for v in component ) + gp.quicksum( m._x[v,j] for v in neighbors_not_S ) >= 1)
                        component_names = [ DG.nodes[v]['NAME20'] for v in component ]
                        neighbor_names = [ DG.nodes[v]['NAME20'] for v in neighbors_not_S ]
                    
                        
    # Check sketch feasibility 
    centers = [ j for j in DG.nodes if xval[j,j] > 0.5 ]
    clusters = [ [ i for i in G.nodes if xval[i,j] > 0.5 ] for j in centers ]
    sizes = [ round( yval[j]) for j in centers ]
    
    for p in range(len(clusters)):
        cluster = clusters[p]
        size = sizes[p]
        GS = G.subgraph(cluster) 
        whole = [ i for i in cluster if i in m._whole_counties ]
        (xsoln, ysoln, zsoln) = sketch(GS, L, U, size, whole_counties=whole, time_limit=60)

        if xsoln is None:
            noncluster = [ i for i in DG.nodes if i not in cluster ]
            c = centers[p]
            m.cbLazy( gp.quicksum( (1-m._x[v,c]) for v in cluster ) + gp.quicksum( m._x[v,c] for v in noncluster ) >= 1 )
    
    return


def feasibility_callback(m, where):
    
    # check if LP relaxation at this BB node is integer
    if where != GRB.Callback.MIPSOL: 
        return

    # retrieve the LP relaxation solution at this BB node
    xval = m.cbGetSolution(m._x)
    yval = m.cbGetSolution(m._y)
    k = m._k
    L = m._L
    U = m._U
    G = m._G
    DG = m._DG
    centers = [ j for j in DG.nodes if xval[j,j] > 0.5 ]
    clusters = [ [ i for i in G.nodes if xval[i,j] > 0.5 ] for j in centers ]
    sizes = [ round( yval[j]) for j in centers ]
    
    for p in range(len(clusters)):
        cluster = clusters[p]
        size = sizes[p]
        GS = G.subgraph(cluster) 
        whole = [ i for i in cluster if i in m._whole_counties ]
        (xsoln, ysoln, zsoln) = sketch(GS, L, U, size, whole_counties=whole, time_limit=60)

        if xsoln is None:
            noncluster = [ i for i in DG.nodes if i not in cluster ]
            c = centers[p]
            m.cbLazy( gp.quicksum( (1-m._x[v,c]) for v in cluster ) + gp.quicksum( m._x[v,c] for v in noncluster ) >= 1 )
            

    return


def max_cluster_main(G, DG, L, U, k, whole_counties, split_counties, max_t=4, restarts=3, time_limit=3600, sketch_intensity=0):   
    
    results = dict()
    
    # incumbents
    (clusters, sizes) = ( [ list(G.nodes) ], [k] )
    cluster_UB = k - sum( math.floor( G.nodes[i]['TOTPOP'] / (U+1) ) for i in G.nodes )
    results['initial_UB'] = cluster_UB
    print("Initially, cluster_UB =",cluster_UB)
    
    if max_t > cluster_UB - 1:
        print("No need for t-opt local search, with t =",max_t,"; reducing to",cluster_UB - 1)
        max_t = cluster_UB - 1
        
    # apply carving construction heuristic and t-opt recom local search
    start_time = time.time()
    for iteration in range(restarts):
        print("\n****************************")
        print("Heuristic iteration #",iteration)
        print("****************************")
        (iclusters, isizes) = carve_heuristic(G, DG, L, U, k, whole_counties, randomized_carving=True, sketch_intensity=sketch_intensity))
        if objective(G, iclusters) > objective(G, clusters):
            (clusters, sizes) = (iclusters, isizes)
            print("new incumbent!")

    results['heuristic_time'] = '{0:.2f}'.format( time.time() - start_time )
    results['heuristic_num_clusters'] = len(clusters)
    results['heuristic_num_cut_edges'] = number_of_cut_edges(G, clusters)
    results['heuristic_iterations'] = restarts
    
    print("\n********************************************************")
    print("After local search, # clusters, #cut edges =",results['heuristic_num_clusters'],results['heuristic_num_cut_edges'] )
    print("********************************************************\n")
    
    # warm start our max_cluster MIP 
    start_time = time.time()
    (max_whole_clusters, max_whole_sizes) = clustering_for_max_whole(G, DG, k, L, U, whole_counties, split_counties, initial_clusters=clusters, initial_sizes=sizes)
    results['MIP_time'] = '{0:.2f}'.format( time.time() - start_time )

    
    print("********************************************************")
    print("MIP gives #clusters, #cut edges =",len(max_whole_clusters),number_of_cut_edges(G, max_whole_clusters) )
    print("********************************************************\n")

    # last, make MIP-generated clusters more compact (only if they are different from before)
    start_time = time.time()
    if len(max_whole_clusters) > results['heuristic_num_clusters']:
        (clusters, sizes) = (max_whole_clusters, max_whole_sizes)
     
    results['cleanup_time'] = '{0:.2f}'.format( time.time() - start_time )
    results['clusters'] = clusters
    results['sizes'] = sizes
    results['num_clusters'] = len(clusters)
    results['num_cut_edges'] = number_of_cut_edges(G, clusters)
    
    return (clusters, sizes, results)

def clustering_for_max_whole(G, DG, k, L, U, whole_counties, split_counties, initial_clusters=None, initial_sizes=None):  

    # setup    
    m = gp.Model()
    m.Params.TimeLimit = 3600

    # x[i,j]=1 if vertex i is assigned to (cluster centered at) vertex j
    x = m.addVars(DG.nodes, DG.nodes, vtype=GRB.BINARY) 

    # y[j] = # of districts in cluster centered at vertex j
    y = m.addVars(DG.nodes, vtype=GRB.INTEGER)

    # BASE CONSTRAINTS
    # Each county i assigned to one cluster
    m.addConstrs( gp.quicksum( x[i,j] for j in DG.nodes ) == 1 for i in DG.nodes )

    # Cluster "sizes" should sum to k
    m.addConstr( y.sum() == k )

    # Population balance: population of cluster j should be an integer multiple of [L,U] 
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) <= U * y[j] for j in DG.nodes )
    m.addConstrs( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i,j] for i in DG.nodes ) >= L * y[j] for j in DG.nodes )

    # Enforce relationships between x[i,j] and x[j,j]
    m.addConstrs( x[i,j] <= x[j,j] for i in DG.nodes for j in DG.nodes if i != j )

    # CONTIGUITY CONSTRAINTS
    # f[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    f = m.addVars(DG.nodes, DG.edges)
    M = DG.number_of_nodes() - 1
    m.addConstrs( gp.quicksum( f[j,u,j] for u in DG.neighbors(j) ) == 0 for j in DG.nodes )
    m.addConstrs( gp.quicksum( f[j,u,i] - f[j,i,u] for u in DG.neighbors(i) ) == x[i,j] for i in DG.nodes for j in DG.nodes if i != j )
    m.addConstrs( gp.quicksum( f[j,u,i] for u in DG.neighbors(i) ) <= M * x[i,j] for i in DG.nodes for j in DG.nodes if i != j )

    # Find vertex ordering
    # ordering = find_ordering(G)
    whole_counties = [ i for i in DG.nodes if i not in split_counties ]
    ordering = vertices_sorted_by_population(DG, split_counties)
    ordering += vertices_sorted_by_population(DG, whole_counties)
    position = construct_position( ordering )

    # Diagonal fixing
    for i in DG.nodes:
        for j in DG.nodes:
            if position[i] < position[j]:
                x[i,j].UB = 0
    m.update()

    # Add Chvatal-Gomory-like valid inequalities (t=1)
    for j in DG.nodes:
        if y[j].UB > 0.5:
            m.addConstr( gp.quicksum( math.floor( DG.nodes[i]['TOTPOP'] / (U+1) ) * x[i,j] for i in DG.nodes ) <= y[j] - x[j,j] )

    # Add cluster-separator inequalities (special cases, anyway)
    # for each vertex b, try a=b and S=N(b)
    for b in DG.nodes:
        if y[b].UB > 0.5 and not is_within_window(DG.nodes[b]['TOTPOP'], L, U):
            m.addConstr( x[b,b] <= sum( x[i,b] for i in DG.neighbors(b) if position[i] > position[b] ) )

    # for each edge {a,b}, try S=N({a,b})
    for a,b in DG.edges:
        if position[a] < position[b] or y[b].UB < 0.5: # a cannot be assigned to b anyway..
            continue
        if not is_within_window(DG.nodes[a]['TOTPOP']+DG.nodes[b]['TOTPOP'], L, U):
            S = nx.node_boundary(G, [a,b])
            S = [ i for i in S if position[i] > position[b] ]
            m.addConstr( x[a,b] <= sum( x[i,b] for i in S ) )

    # for each vertex b, try a=b and all S=N(component) for all components "around" b. 
    for b in G.nodes:
        if y[b].UB < 0.5:
            continue
        for size in range(2,1+G.degree(b)):
            for comb in combinations(G.neighbors(b),size):
                component = list(comb) + [b]
                component = [ i for i in component if position[i] >= position[b] ]
                if not is_cluster_possible( G.subgraph(component), [b], L, U, k, verify=False ):
                    S = nx.node_boundary(G, component)
                    S = [ i for i in S if position[i] > position[b] ]
                    m.addConstr( x[b,b] <= sum( x[i,b] for i in S ) )
                    
    # Provide MIP start
    if initial_clusters is not None:
        assert len(initial_clusters) == len(initial_sizes)
        for p in range(len(initial_clusters)):
            cluster = initial_clusters[p]
            min_pos = min( position[i] for i in cluster )
            j = [ i for i in cluster if position[i] == min_pos ][0]
            y[j].start = initial_sizes[p]
            for i in cluster:
                x[i,j].start = 1 
                
    # Parameter settings
    m.Params.IntFeasTol = 1.e-7
    m.Params.FeasibilityTol = 1.e-7
    for j in DG.nodes:
        x[j,j].BranchPriority = 1
        y[j].BranchPriority = 1
        
    # if non-split county, then can only root size 1 clusters
    m.addConstrs( y[i] == x[i,i] for i in whole_counties )

    # if whole-county i is assigned to split-county j, 
    # then at least one split-county v "nearby" i is as well.
    #no_near = list()
    for i in whole_counties:
        dist = nx.shortest_path_length( DG, source = i, weight = 'pweight' )
        nearby_split_counties = [ v for v in split_counties if dist[v] <= U ]
        m.addConstrs( x[i,j] <= gp.quicksum( x[v,j] for v in nearby_split_counties ) for j in split_counties )
    
    # callback info
    m._G = G
    m._DG = DG
    m._x = x
    m._k = k
    m._L = L
    m._U = U
    m._whole_counties = whole_counties
    m._split_counties = split_counties
    m._y = y
    m.Params.LazyConstraints = 1
    
    print("\n***First, solve using MOI compactness objective to quickly find a feasible county clustering.")
    for i in DG.nodes:
        dist = nx.shortest_path_length(DG, source = i)
        for j in DG.nodes:
            m._x[i,j].obj = dist[j] * dist[j] * round( DG.nodes[i]['TOTPOP'] / 1000 )
    m._callback = cut_callback
    m.optimize(m._callback)
    
    if m.solCount == 0:
        print("Unable to find clustering. Exiting.")
        clusters = None
        sizes = None
        return (clusters, sizes)

    print("\n**Second, re-solve using max clustering objective.")
    m.setObjective( gp.quicksum( m._x[j,j] for j in G.nodes ), GRB.MAXIMIZE )
    m._callback = cut_callback
    m.optimize(m._callback)
    
    # Retrieve solution
    centers = [ j for j in DG.nodes if m._x[j,j].x > 0.5 ]
    clusters = [ [ i for i in G.nodes if m._x[i,j].x > 0.5 ] for j in centers ]
    sizes = [ round( y[j].x ) for j in centers ]
    # double check the populations and sizes
    assert sum( sizes ) == k
    for p in range(len(sizes)):
        pop = sum( G.nodes[i]['TOTPOP'] for i in clusters[p] )
        assert pop >= L * sizes[p] and pop <= U * sizes[p]

    # z[u,v,j] = 1 if arc (u,v) is cut because u->j but v!->j
    z = m.addVars(DG.edges, G.nodes, vtype=GRB.BINARY)
    is_cut = m.addVars(G.edges, vtype=GRB.BINARY)

    m.addConstrs( m._x[u,j] - m._x[v,j] <= z[u,v,j] for u,v in DG.edges for j in G.nodes )
    m.addConstrs( is_cut[u,v] == gp.quicksum( z[u,v,j] for j in G.nodes ) for u,v in G.edges )
    m.addConstrs( is_cut[u,v] == gp.quicksum( z[v,u,j] for j in G.nodes ) for u,v in G.edges )

    m.setObjective( is_cut.sum(), GRB.MINIMIZE)

    # add constraint on # clusters
    m.addConstr( gp.quicksum( m._x[i,i] for i in G.nodes ) == len(sizes) )

    print("\n***Third, solve using cut edge compactness objective to get a pretty county clustering.")
    m._callback = cut_callback
    m.optimize(m._callback)

    # Retrieve solution
    if m.solCount > 0:
        centers = [ j for j in DG.nodes if m._x[j,j].x > 0.5 ]
        clusters = [ [ i for i in G.nodes if m._x[i,j].x > 0.5 ] for j in centers ]
        sizes = [ round( y[j].x ) for j in centers ]
        cut_edge_LB = round( m.objBound )
        cut_edge_UB = round( m.objVal )
        # double check the populations and sizes
        assert sum( sizes ) == k
        for p in range(len(sizes)):
            pop = sum( G.nodes[i]['TOTPOP'] for i in clusters[p] )
            assert pop >= L * sizes[p] and pop <= U * sizes[p]
    else:
        clusters = [ list(G.nodes) ]
        sizes = [k]

    print("for min cut edge model (with max # clusters), we have:")
    print("centers =",centers)
    print("clusters =",clusters)
    print("sizes = ",sizes)

    for p in range(len(clusters)):
        cluster = clusters[p]
        size = sizes[p]
        cluster_names = {i : G.nodes[i]['NAME20'] for i in cluster }
        print("cluster has these nodes:",cluster)
        print("this size:",size,"where L,U =",L,U)
        print("cluster_names =",cluster_names)

        GS = G.subgraph(cluster)
        node_pos = { i : ( GS.nodes[i]['X'], GS.nodes[i]['Y'] ) for i in GS.nodes }
        node_colors = [ "red" if i in split_counties else "green" for i in GS.nodes ]

        my_labels = { i : G.nodes[i]['NAME20'] + "\n" + f"{GS.nodes[i]['TOTPOP']:,}" for i in GS.nodes }
        nx.draw(GS, pos=node_pos, with_labels=True, labels=my_labels, node_color=node_colors, node_size=2000)
        plt.show()
             
    
    return (clusters, sizes)
