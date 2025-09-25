import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from itertools import combinations
from math import floor, ceil
import matplotlib.pyplot as plt

from utils import printif, check_clustering
from sketch_for_max_whole import sketch

def cluster_for_max_whole(G, L, U, k, whole_counties, sketch_time_limit=None, time_limit=3600, clusters_start=None, sizes_start=None, verbose=True):  

    # setup    
    m = gp.Model()

    # x[i,j]=1 if vertex i is assigned to (cluster centered at) vertex j
    x = m.addVars(G.nodes, G.nodes, vtype=GRB.BINARY) 

    # y[j] = # of districts in cluster centered at vertex j
    y = m.addVars(G.nodes, vtype=GRB.INTEGER)

    # *binarize* the y vars
    # z[j,t] = 1  iff  y[j] = t
    z = m.addVars(G.nodes, range(1,k+1), vtype=GRB.BINARY)
    m.addConstrs( y[j] == sum( t * z[j,t] for t in range(1,k+1) ) for j in G.nodes )
    m.addConstrs( sum( z[j,t] for t in range(1,k+1) ) <= x[j,j] for j in G.nodes )

    # BASE CONSTRAINTS
    # Each county i assigned to one cluster
    m.addConstrs( sum( x[i,j] for j in G.nodes ) == 1 for i in G.nodes )

    # Cluster "sizes" should sum to k
    m.addConstr( y.sum() == k )

    # Population balance: population of cluster j should be an integer multiple of [L,U] 
    m.addConstrs( sum( G.nodes[i]['TOTPOP'] * x[i,j] for i in G.nodes ) <= U * y[j] for j in G.nodes )
    m.addConstrs( sum( G.nodes[i]['TOTPOP'] * x[i,j] for i in G.nodes ) >= L * y[j] for j in G.nodes )

    # Coupling constraints between x[i,j] and x[j,j]
    m.addConstrs( x[i,j] <= x[j,j] for i in G.nodes for j in G.nodes if i != j )

    # CONTIGUITY CONSTRAINTS
    # f[j,u,v] tells how much flow (from source j) is sent across arc (u,v)
    DG = nx.DiGraph(G)
    f = m.addVars(G.nodes, DG.edges)
    M = G.number_of_nodes() - 1
    m.addConstrs( gp.quicksum( f[j,u,j] for u in G.neighbors(j) ) == 0 for j in G.nodes )
    m.addConstrs( gp.quicksum( f[j,u,i] - f[j,i,u] for u in G.neighbors(i) ) == x[i,j] for i in G.nodes for j in G.nodes if i != j )
    m.addConstrs( gp.quicksum( f[j,u,i] for u in G.neighbors(i) ) <= M * x[i,j] for i in G.nodes for j in G.nodes if i != j )

    # Find vertex ordering
    split_counties = [ i for i in G.nodes if i not in whole_counties ]
    ordering = vertices_sorted_by_population(G, split_counties)
    ordering += vertices_sorted_by_population(G, whole_counties)
    position = construct_position( ordering )

    # Diagonal fixing
    for i in G.nodes:
        for j in G.nodes:
            if position[i] < position[j]:
                x[i,j].UB = 0
    m.update()

    # Add Chvatal-Gomory-like ``Rounding inequalities'' (the case t=1)
    for j in G.nodes:
        if y[j].UB > 0.5:
            m.addConstr( sum( floor( DG.nodes[i]['TOTPOP'] / (U+1) ) * x[i,j] for i in G.nodes ) <= y[j] - x[j,j] )

    # Add ``cluster-separator inequalities'' (special cases, anyway)
    # for each vertex b, try a=b and S=N(b)
    for b in G.nodes:
        if y[b].UB > 0.5 and not is_within_window( G.nodes[b]['TOTPOP'], L, U ):
            m.addConstr( x[b,b] <= sum( x[i,b] for i in G.neighbors(b) if position[i] > position[b] ) )

    # for each edge {a,b}, try S=N({a,b})
    for a,b in DG.edges:
        if position[a] < position[b] or y[b].UB < 0.5: # a cannot be assigned to b anyway..
            continue
        if not is_within_window(G.nodes[a]['TOTPOP']+G.nodes[b]['TOTPOP'], L, U):
            S = nx.node_boundary(G, [a,b])
            S = [ i for i in S if position[i] > position[b] ]
            m.addConstr( x[a,b] <= sum( x[i,b] for i in S ) )

    # for each vertex b, try a=b and all S=N(component) for all components "around" b. 
    for b in G.nodes:
        if y[b].UB < 0.5:
            continue
        for size in range(2,1+G.degree(b)):
            for comb in combinations( G.neighbors(b), size ):
                component = list(comb) + [b]
                component = [ i for i in component if position[i] >= position[b] ]
                if not is_cluster_possible( G.subgraph(component), [b], L, U, k, verify=False ):
                    S = nx.node_boundary(G, component)
                    S = [ i for i in S if position[i] > position[b] ]
                    m.addConstr( x[b,b] <= sum( x[i,b] for i in S ) )
                    
    # Parameter settings
    m.Params.MIPGap = 0.00
    m.Params.IntFeasTol = 1.e-7
    m.Params.FeasibilityTol = 1.e-7
    for j in G.nodes:
        x[j,j].BranchPriority = 1
        y[j].BranchPriority = 1
        
    # if non-split county, then can only root size one clusters
    m.addConstrs( y[j] == x[j,j] for j in whole_counties )

    # nearby (or vicinity-like) inequalities:
    #   if whole-county i is assigned to split-county j, 
    #   then at least one split-county v "nearby" i is as well.
    for i in whole_counties:
        dist = nx.shortest_path_length( DG, source = i, weight = 'pweight' )
        nearby_split_counties = [ v for v in split_counties if dist[v] <= U ]
        m.addConstrs( x[i,j] <= gp.quicksum( x[v,j] for v in nearby_split_counties ) for j in split_counties )

    #   if split-county i is assigned to split-county j (j!=i),
    #   then at least one split-county v (v!=i) "nearby" i is as well.
    for i in split_counties:
        dist = nx.shortest_path_length( DG, source = i, weight = 'pweight' )
        nearby_split_counties = [ v for v in split_counties if v!=i and dist[v] <= U+DG.nodes[i]['TOTPOP'] ]
        m.addConstrs( x[i,j] <= gp.quicksum( x[v,j] for v in nearby_split_counties ) for j in split_counties if i!=j )
    
    # Provide MIP start
    gtv = { G.nodes[i]['GEOID20'] : i for i in G.nodes }
    if clusters_start is not None:
        assert len(clusters_start) == len(sizes_start)
        for p in range(len(clusters_start)):
            cluster = clusters_start[p]
            min_pos = min( position[gtv[g]] for g in cluster )
            j = [ gtv[g] for g in cluster if position[gtv[g]] == min_pos ][0]
            y[j].start = sizes_start[p]
            for g in cluster:
                x[gtv[g],j].start = 1 

    # callback info
    m._G = G
    m._DG = DG
    m._x = x
    m._y = y
    m._z = z
    m._k = k
    m._L = L
    m._U = U
    m._whole_counties = whole_counties
    m._split_counties = split_counties
    m._verbose = verbose
    m._sketch_time_limit = sketch_time_limit
    m.Params.LazyConstraints = 1
    
    # set objective: maximize # clusters
    #m.Params.MIPGap = 0.05 # we don't need county clustering to be "optimal"
    #m.setObjective( sum( m._x[j,j] for j in G.nodes ), GRB.MAXIMIZE )

    # set objective: cluster compactness
    for i in G.nodes:
        dist = nx.single_source_shortest_path_length( G, source = i )
        for j in dist.keys():
            m._x[i,j].obj = dist[j] * ceil( G.nodes[i]['TOTPOP'] / 1000 )

    # singleton starts
    for j in G.nodes:
        if is_within_window(G.nodes[j]['TOTPOP'], L, U):
            for i in G.nodes:
                m._x[i,j].start = 1 if i==j else 0
    
    # optimize
    m.Params.TimeLimit = time_limit
    m._callback = cluster_callback
    m.optimize(m._callback)

    if m.status == GRB.INFEASIBLE:
        printif(verbose, "Infeasible. Unable to find clustering. Exiting.")
        return 'infeasible'
    elif m.solCount == 0:
        printif(verbose, "Inconclusive. Unable to find a suitable clustering but did not prove infeasibility. Exiting.")
        return 'inconclusive'

    # Retrieve solution
    centers = [ j for j in G.nodes if m._x[j,j].x > 0.5 ]
    clusters = [ [ i for i in G.nodes if m._x[i,j].x > 0.5 ] for j in centers ]
    sizes = [ round( m._y[j].x ) for j in centers ]

    # check solution
    check_clustering( G, L, U, k, clusters, sizes )
    return (clusters, sizes)

def sort_by_second(val):
    return val[1]

def find_ordering(G):
    V_with_population = [ (i, G.nodes[i]['TOTPOP']) for i in G.nodes ]
    V_with_population.sort( key=sort_by_second, reverse=True )
    return [ v for (v,p) in V_with_population ]

def construct_position(ordering):
    return { ordering[p] : p for p in range(len(ordering)) }

def is_within_window(population, L, U):
    min_size = ceil( population / U )
    max_size = floor( population / L )
    return min_size <= max_size

def is_cluster_possible(G, required_nodes, L, U, k, verify=True):
    (cluster, size) = find_cluster(G, required_nodes, L, U, k, verify) 
    return cluster is not None

def vertices_sorted_by_population(DG, vertices=list()):
    if vertices==list():
        vertices = [ i for i in DG.nodes ]
    nodes_with_population = [ ( i, DG.nodes[i]['TOTPOP'] ) for i in vertices ]
    nodes_with_population.sort(key=sort_by_second,reverse=True)
    return [ v for (v,p) in nodes_with_population ]
    
def find_cluster(G, required_nodes, L, U, k, verify=True):
    
    # To make contiguity checks easy, this function assumes that
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
    if m.solCount > 0:  # ( cluster, size )
        return ( [ i for i in G.nodes if x[i].x > 0.5 ], round(z.x) ) 
    return (None, None)
    
def cluster_callback(m, where):

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
    added_cut = False

    num_clusters = sum( round( xval[j,j] ) for j in G.nodes )
    printif(m._verbose, f"Entered callback with # clusters = { num_clusters }.")

    centers = [ j for j in G.nodes if xval[j,j] > 0.5 ]
    clusters = [ [ i for i in G.nodes if xval[i,j] > 0.5 ] for j in centers ]
    sizes = [ round( yval[j]) for j in centers ]

    for (cluster, size, j) in zip(clusters, sizes, centers):

        printif(m._verbose, f"Seeking {size} districts for cluster = {cluster}.")
    
        # get sketch
        whole = [ c for c in m._whole_counties if c in cluster ]
        GS = G.subgraph( cluster )
        sketch_soln = sketch( GS, L, U, size, whole_counties=whole, time_limit=m._sketch_time_limit)
        
        if sketch_soln == 'infeasible':
            printif(m._verbose, "For this cluster, sketch is infeasible. Adding (stronger) no-good cut.")
            neigh = { v for i in cluster for v in G.neighbors(i) if v not in cluster }
            m.cbLazy( (1-m._z[j,size]) + sum( (1-m._x[i,j]) for i in cluster ) + sum( m._x[i,j] for i in neigh ) >= 1 )
        
        printif(m._verbose and sketch_soln == 'inconclusive', "FYI. For this cluster, sketch is inconclusive. Allow it, for now.")
    
    return
    