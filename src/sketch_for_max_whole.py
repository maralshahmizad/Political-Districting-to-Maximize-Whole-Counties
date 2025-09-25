import networkx as nx
import gurobipy as gp
from gurobipy import GRB  
from math import ceil
from utils import printif
from coarsen import graph_coarsen
import matplotlib.pyplot as plt

# An improved version of sketch. It works as follows:
# 1. create a coarsened graph where each whole county from W remains by itself and 
#      the split counties are merged together if they belong to same component of G-W.
# 2. seek a sketch in this coarsened graph (using symmetry breaking tricks); if no sketch exists, 
#      then the state-wide (or cluster-wide) sketch problem is infeasible, so exit.
# 3. otherwise, if the coarsened graph admits a sketch, then use it to find a statewide (or cluster-wide) sketch. 
#      Basically it just re-solves the sketch model over the whole graph, but requires the support of each 
#      district to be the same (or subset thereof). This essentially fixes the districts of the whole 
#      counties and has to figure out the districts for the split counties.

def sketch(G, L, U, k, whole_counties=list(), draw=False, statewide=True, time_limit=None, verbose=False):
    assert k > 0
    
    S = [ i for i in G.nodes if i not in whole_counties ]
    
    partition = [ list(comp) for comp in nx.connected_components( G.subgraph(S) ) ] + [ [c] for c in whole_counties ]
    coarse_graph = graph_coarsen(G, partition)
    coarse_graph_W = [ p for p in range(len(partition)) if len(partition[p])==1 and partition[p][0] in whole_counties ]
    
    soln = sketch_subroutine(coarse_graph, L, U, k, whole_counties=coarse_graph_W, symmetry_breaking=True, draw=draw, time_limit=time_limit, verbose=verbose)
    if soln in ['infeasible', 'inconclusive']:
        return soln

    # return early if we don't care about getting a statewide sketch
    if not statewide: 
        return 'inconclusive'

    printif(verbose, "Sketch on coarse graph was successful. Unpacking and seeking sketch for original graph.")
    (xsoln, ysoln, zsoln) = soln
    support = { v : list( xsoln[i].keys() ) for i in coarse_graph.nodes for v in partition[i] }
    return sketch_subroutine(G, L, U, k, whole_counties=whole_counties, allowable_support=support, draw=draw, time_limit=time_limit, verbose=verbose)
    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the sketch routine (without coarsening) but with all other MIP tricks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
def sketch_subroutine(G, L, U, k, whole_counties=list(), symmetry_breaking=False, allowable_support=None, draw=False, time_limit=None, verbose=False):

    if symmetry_breaking and allowable_support:
        print("WARNING: Symmetry breaking and support restrictions may conflict. Tread lightly.")
    
    m = gp.Model()
    m.Params.OutputFlag = 1 if verbose else 0
    if time_limit is not None:
        m.Params.TimeLimit = time_limit

    # x[i,j]=1 if county i is assigned to district j (in whole or in part)
    x = m.addVars( G.nodes, k, vtype=GRB.BINARY )

    # y[e,j]=1 if edge e is assigned to district j (in whole or in part)
    y = m.addVars( G.edges, k, vtype=GRB.BINARY )

    # ysum[e]=1 if edge e is used at all
    ysum = m.addVars( G.edges, vtype=GRB.BINARY )
    m.addConstrs( ysum[u,v] == gp.quicksum( y[u,v,j] for j in range(k) ) for u,v in G.edges )
    
    # z[i,j] = the proportion of county i that is assigned to district j
    z = m.addVars( G.nodes, k )
    #min_proportion_constraints = m.addConstrs( z[i,j] >= 0.0125 * x[i,j] for i in G.nodes for j in range(k) )
    
    # s[i] = the number of splits permitted in county i
    s = m.addVars( G.nodes, vtype=GRB.INTEGER )
    for i in G.nodes:
        s[i].LB = ceil( G.nodes[i]['TOTPOP'] / U ) - 1  
        
    for i in whole_counties:
        s[i].UB = 0
    
    # indegree extended formulation
    m.addConstrs( gp.quicksum( x[i,j] for i in G.nodes ) - gp.quicksum( y[u,v,j] for u,v in G.edges ) <= 1 for j in range(k) )

    # *tentatively* pick edge at most once
    edge_constraints = m.addConstrs( gp.quicksum( y[u,v,j] for j in range(k) ) <= 1 for u,v in G.edges )
    
    # if pick edge, then must pick endpoints
    for u,v in G.edges:
        m.addConstrs( y[u,v,j] <= x[u,j] for j in range(k) )
        m.addConstrs( y[u,v,j] <= x[v,j] for j in range(k) )
        
    # if pick endpoints, then must pick edge
    for u,v in G.edges:
        m.addConstrs( y[u,v,j] >= x[u,j] + x[v,j] - 1 for j in range(k) )

    # assignment constraints
    m.addConstrs( gp.quicksum( z[i,j] for j in range(k) ) == 1 for i in G.nodes )
        
    # population balance
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * z[i,j] for i in G.nodes ) >= L for j in range(k) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * z[i,j] for i in G.nodes ) <= U for j in range(k) )

    # coupling constraints
    m.addConstrs( z[i,j] <= x[i,j] for i in G.nodes for j in range(k) )

    # splitting constraints
    m.addConstrs( gp.quicksum( x[i,j] for j in range(k) ) == s[i] + 1 for i in G.nodes )

    # *tentatively* add min-split constraint 
    min_split_constraint = m.addConstr( s.sum() == k-1 )

    # restrict the support?
    if allowable_support:
        for i in G.nodes:
            for j in range(k):
                x[i,j].UB = 0
                z[i,j].UB = 0
        for i in allowable_support.keys():
            for j in allowable_support[i]:
                x[i,j].UB = 1
                z[i,j].UB = 1
        
    # break some model symmetry?
    if symmetry_breaking:
        
        # stablize some districts within split counties
        ss = 0
        for i in G.nodes:
            if i in whole_counties or not all( j in whole_counties for j in G.neighbors(i) ):
                continue
            remaining_pop = max( 0, G.nodes[i]['TOTPOP'] - sum( (U-G.nodes[v]['TOTPOP']) for v in G.neighbors(i) ))
            printif(verbose,f"for i = {i}, remaining_pop = {remaining_pop}, requiring at least {ceil( remaining_pop / U )} districts.")  
            if remaining_pop <= 0:
                continue
            for j in range( ceil(remaining_pop / U) ):
                x[i,ss+j].LB = 1
            for j in range(ceil(remaining_pop / U)-1):
                m.addConstr( z[i,ss+j] >= z[i,ss+j+1] )
                for v in G.nodes:
                    if v != i:
                        x[v,ss+j].UB = 0
                        z[v,ss+j].UB = 0
            ss += ceil( remaining_pop / U )

        # stabilize some districts at whole counties
        DG = nx.DiGraph(G)
        for (i,j) in DG.edges:
            DG.edges[i,j]['whole_weight'] = G.nodes[i]['TOTPOP'] if i in whole_counties else 0
        compatibility_graph = nx.Graph()
        compatibility_graph.add_nodes_from( whole_counties )
        for i in whole_counties:
            dist = nx.shortest_path_length(DG, source=i, weight='whole_weight')
            for j in dist.keys():
                if i<j and j in whole_counties and dist[j] + G.nodes[j]['TOTPOP'] <= U:
                    compatibility_graph.add_edge(i,j)
    
        ind_set = maximum_independent_set(compatibility_graph)
        if ss + len(ind_set) > k:
            printif(verbose,f"Symmetry breaking proves sketch infeasibility: ss = {ss} and ind_set = {ind_set} and ss + len(ind_set) > k = {k}.")
            return 'infeasible'
        for p in range(len(ind_set)):
            i = ind_set[p]
            x[i,ss+p].LB = 1
            z[i,ss+p].LB = 1
            dist = nx.shortest_path_length(DG, source=i, weight='whole_weight')
            for j in dist.keys():
                whole_weight = G.nodes[j]['TOTPOP'] if j in whole_counties else 0
                if dist[j] + whole_weight > U:
                    x[j,ss+p].UB = 0
                    z[j,ss+p].UB = 0

        printif(verbose,f"Symmetry breaking stabilizes {len(ind_set)} districts using conflicting whole counties and {ss} districts in split counties.") 
    
    # contiguity constraints
    m.Params.LazyConstraints = 1
    m._DG = nx.DiGraph(G)
    m._parts = k
    m._x = x

    # first, maximize preserved edges
    m.setObjective( y.sum(), GRB.MAXIMIZE )
    m._callback = labeling_contiguity_callback
    m.Params.MIPGap = 0.34  # it's not critical to get true minimum here.
    m.optimize(m._callback)

    # if no sketch has been found, remove min split and edge constraints and try again
    if m.solCount == 0:
        m.remove( min_split_constraint )
        m.remove( edge_constraints )
        m._callback = labeling_contiguity_callback
        m.optimize(m._callback)

    if m.status == GRB.INFEASIBLE:
        return 'infeasible'
    elif m.solCount == 0:
        return 'inconclusive'

    # fix the number of preserved edges
    obj = round( m.objVal )
    m.addConstr( y.sum() >= obj )
    
    # c[j,t] = 1 if district j touches t counties
    c = m.addVars( k, range(1,G.number_of_nodes()+1), vtype=GRB.BINARY )
    m.addConstrs( gp.quicksum( c[j,t] for t in range(1,G.number_of_nodes()+1) ) == 1 for j in range(k) )
    m.addConstrs( gp.quicksum( t * c[j,t] for t in range(1,G.number_of_nodes()+1) ) == gp.quicksum( x[i,j] for i in G.nodes ) for j in range(k) )

    # minimize sum of squares of counties touched 
    m.setObjective( gp.quicksum( t * t * c[j,t] for j in range(k) for t in range(1,G.number_of_nodes()+1) ), GRB.MINIMIZE )
    m._callback = labeling_contiguity_callback
    m.optimize(m._callback)

    xsoln = { i : { j : x[i,j].x for j in range(k) if x[i,j].x > 0 } for i in G.nodes }
    ysoln = { (u,v) : { j : y[u,v,j].x for j in range(k) if y[u,v,j].x > 0 } for (u,v) in G.edges }
    zsoln = { i : { j : z[i,j].x for j in range(k) if z[i,j].x > 0 } for i in G.nodes }

    if draw:
        node_pos = { i : ( G.nodes[i]['X'], G.nodes[i]['Y'] ) for i in G.nodes }
        node_colors = [ list(xsoln[i])[0] for i in G.nodes ]
        edge_colors = [ list(ysoln[u,v])[0] if len(ysoln[u,v]) > 0 else -1 for u,v in G.edges ]

        my_labels = { i : list(xsoln[i].keys()) for i in G.nodes }
        nx.draw(G, pos=node_pos, node_color=node_colors, edge_color=edge_colors, with_labels=True, labels=my_labels, font_color="red")
        plt.show()

    return (xsoln, ysoln, zsoln)

def maximum_independent_set(G):
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars(G.nodes, vtype=GRB.BINARY)
    m.setObjective( x.sum(), GRB.MAXIMIZE )
    m.addConstrs( x[i] + x[j] <= 1 for i,j in G.edges )
    m.optimize()
    assert m.solCount > 0
    return [ i for i in G.nodes if x[i].x > 0.5 ]

def find_fischetti_separator(DG, component, b):
    neighbors_component = { i : False for i in DG.nodes }
    for i in nx.node_boundary(DG, component, None):
        neighbors_component[i] = True
    
    visited = { i : False for i in DG.nodes }
    child = [b]
    visited[b] = True
    
    while child:
        parent = child
        child = list()
        for i in parent:
            if not neighbors_component[i]:
                for j in DG.neighbors(i):
                    if not visited[j]:
                        child.append(j)
                        visited[j] = True
    
    return [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]

def labeling_contiguity_callback(m, where):

    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve the LP solution
    xval = m.cbGetSolution(m._x)
    DG = m._DG

    for j in range(m._parts):

        # which nodes are selected?
        S = [ i for i in DG.nodes if xval[i,j] > 0.5 ]
        if len(S)==0:
            print("WARNING: part j =",j,"is empty.")
            continue
        if len(S)==1 or nx.is_strongly_connected( DG.subgraph(S) ):
            continue

        # find a maximum population component
        b = None
        max_component_population = -1
        for component in nx.strongly_connected_components( DG.subgraph(S) ):
            component_population = sum( DG.nodes[i]['TOTPOP'] for i in component )
            if component_population > max_component_population:

                # find a maximum population vertex 'b' from this component
                max_component_population = component_population
                max_vertex_population = max( DG.nodes[i]['TOTPOP'] for i in component )
                max_population_vertices = [ i for i in component if DG.nodes[i]['TOTPOP'] == max_vertex_population ]
                b = max_population_vertices[0]

        # each other component (besides b's), find some vertex 'a' and add cut.
        for component in nx.strongly_connected_components( DG.subgraph(S) ):
            if b in component:
                continue

            # find a maximum population vertex 'a' from this component
            max_vertex_population = max( DG.nodes[i]['TOTPOP'] for i in component )
            max_population_vertices = [ i for i in component if DG.nodes[i]['TOTPOP'] == max_vertex_population ]
            a = max_population_vertices[0]

            # add a,b-separator inequality
            C = find_fischetti_separator(DG, component, b)
            for t in range(m._parts):
                m.cbLazy( m._x[a,t] + m._x[b,t] <= 1 + gp.quicksum( m._x[c,t] for c in C ) )
