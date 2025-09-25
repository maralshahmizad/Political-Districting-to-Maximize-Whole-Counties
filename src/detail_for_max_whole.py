import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from math import isclose

from sketch_for_max_whole import find_fischetti_separator
from cluster_for_max_whole import construct_position, sort_by_second
from coarsen import graph_coarsen_by_county, graph_coarsen
from utils import check_plan, printif

def detail(G, L, U, k, county_geoid_support=None, verbose=False, time_limit=600):
    
    # before detailing, should we coarsen the graph first? (yes--if support is nontrivial and at least one whole county)
    do_coarsen = ( county_geoid_support is not None ) and any( len(supp) == 1 for supp in county_geoid_support.values() )
    if do_coarsen: 
        split_fips = { g for g, supp in county_geoid_support.items() if len( supp ) > 1 }
        (coarse_graph, partition) = graph_coarsen_by_county(G, granular_fips=split_fips, verbose=False)
        coarse_plan = detail_subroutine(coarse_graph, L, U, k, county_geoid_support=county_geoid_support, verbose=verbose, time_limit=time_limit)
        if coarse_plan is None:
            return coarse_plan
        plan = [ [ i for p in coarse_district for i in partition[p] ] for coarse_district in coarse_plan ]
    else:
        plan = detail_subroutine(G, L, U, k, county_geoid_support=county_geoid_support, verbose=verbose, time_limit=time_limit)
        if plan is None:
            return plan
    
    check_plan(G, L, U, k, plan)
    return plan

def detail_subroutine(G, L, U, k, county_geoid_support=None, allow_coarsen=True, verbose=False, time_limit=600):

    assert k >= 1 and nx.is_connected( G )
    if k == 1:
        total_population = sum( G.nodes[i]['TOTPOP'] for i in G.nodes )
        assert L <= total_population and total_population <= U
        if county_geoid_support is not None:
            for i in G.nodes:
                g = G.nodes[i]['GEOID20'][0:5]
                assert 0 in county_geoid_support[g] 
        return [ list(G.nodes) ]
    
    # tract- or block-level?
    max_len = max( len( G.nodes[v]['GEOID20'] ) for v in G.nodes )
    if max_len == 15:
        level = 'block'
    elif max_len == 11:
        level = 'tract'
    else:
        level = 'vtd'

    print(f"Trying {level}-level instance with n = {G.number_of_nodes()}.")

    # create a (sparse) directed "assignment graph" AG 
    # that has vertices V = V(G) + {0, 1, ..., k-1}  (which may overlap)
    # and edge (i,j) if node i can be assigned to district j in {0, 1, ..., k-1}
    nodes = list( set(G.nodes) | set(range(k)) )
    AG = nx.DiGraph()
    AG.add_nodes_from( nodes )
    
    # add edges to AG
    if county_geoid_support is None:
        AG.add_edges_from( [(i,j) for i in G.nodes for j in range(k) ])
    else:
        for i in G.nodes:
            g = G.nodes[i]['GEOID20'][0:5]
            AG.add_edges_from( [(i,j) for j in county_geoid_support[g] ] )
    
    # build model
    m = gp.Model()
    m.Params.OutputFlag = 0
    x = m.addVars( AG.edges, vtype=GRB.CONTINUOUS )
    
    # assignment constraints
    m.addConstrs( gp.quicksum( x[i,j] for j in AG.neighbors(i) ) == 1 for i in G.nodes )
    
    # population balance
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in AG.predecessors(j) ) >= L for j in range(k) )
    m.addConstrs( gp.quicksum( G.nodes[i]['TOTPOP'] * x[i,j] for i in AG.predecessors(j) ) <= U for j in range(k) )
    
    max_iterations = 100
    printif(verbose, "Iteration \t Objective \t iter_time")
    
    for iteration in range(max_iterations):
    
        if iteration == 0:
            old_total_cost = -1               # dummy value
            m.setObjective( 0, GRB.MINIMIZE ) # dummy objective
        else:
            old_total_cost = total_cost
            m.setObjective( gp.quicksum( (1+G.nodes[i]['TOTPOP']) * sq_eucl_dist(G.nodes[i]['X'],G.nodes[i]['Y'],means_x[j],means_y[j]) * x[i,j] for (i,j) in AG.edges ), GRB.MINIMIZE )

        m.optimize()
        if verbose:
            print(iteration,"\t",'{0:.2f}'.format(m.objVal),"\t",'{0:.2f}'.format(m.runtime))
    
        total_cost = m.objVal
        clusters = [ [ i for i in AG.predecessors(j) if x[i,j].x > 0.5 ] for j in range(k) ]

        if isclose( total_cost, old_total_cost):
            break
        else:
            # get next means
            population = [ sum( G.nodes[i]['TOTPOP'] * x[i,j].x for i in AG.predecessors(j) ) for j in range(k) ]
            means_x = [ sum( G.nodes[i]['TOTPOP'] * G.nodes[i]['X'] * x[i,j].x for i in AG.predecessors(j) ) / population[j] for j in range(k) ]
            means_y = [ sum( G.nodes[i]['TOTPOP'] * G.nodes[i]['Y'] * x[i,j].x for i in AG.predecessors(j) ) / population[j] for j in range(k) ]
    
    # retrieve solution
    lp_obj = m.objVal
    soln = { (i,j) : x[i,j].x for (i,j) in AG.edges }

    if allow_coarsen:
        buckets = list()
        bucketed = list()
        for j in range(k):
            district = [ i for i in AG.predecessors(j) if x[i,j].x > 0.99 ]
            printif(verbose, f"District {j} has components with # vertices = { [ len(comp) for comp in nx.connected_components( G.subgraph(district) ) ]}.")
            printif(verbose, f"Their interiors have # vertices = { [ len(interior(G, set(comp))) for comp in nx.connected_components( G.subgraph(district) ) ]}.")
            for comp in nx.connected_components( G.subgraph(district) ):
                I = interior(G, set(comp) )
                I_fips = { G.nodes[i]['GEOID20'][0:5] for i in I }
                for f in I_fips:
                    I_f = [ i for i in I if G.nodes[i]['GEOID20'][0:5] == f ]
                    for I_f_comp in nx.connected_components( G.subgraph(I_f) ):
                        if len(I_f_comp) >= 10:
                            buckets.append( list( I_f_comp ) )
                            bucketed += list( I_f_comp )
        
        if len(bucketed) > 0:
            
            unbucketed = set( G.nodes ) - set( bucketed )
            partition = buckets + [ [i] for i in unbucketed ]
            coarse_graph = graph_coarsen(G, partition)
            printif(verbose, f"In detail(), coarsened graph based on constrained k-means solution from {G.number_of_nodes()} to {len(partition)} vertices.")
            for p in range(len(partition)):
                v = partition[p][0]
                if len(partition[p]) == 1:
                    coarse_graph.nodes[p]['GEOID20'] = G.nodes[v]['GEOID20']
                else:
                    coarse_graph.nodes[p]['GEOID20'] = G.nodes[v]['GEOID20'][0:5] 
            coarse_plan = detail_subroutine(coarse_graph, L, U, k, county_geoid_support=county_geoid_support, allow_coarsen=False, verbose=verbose, time_limit=time_limit)
            if coarse_plan is not None:
                return [ [ i for p in coarse_district for i in partition[p] ] for coarse_district in coarse_plan ]
            else:
                printif(verbose,"Unable to find plan in coarsened graph. Trying uncoarsened next.")
    
    # convert to binary vars
    for i,j in AG.edges:
        x[i,j].vtype = GRB.BINARY

    # parameters
    m.Params.TimeLimit = time_limit
    m.Params.MIPGap = 0.10 
    m.Params.FeasibilityTol = 1e-7
    m.Params.IntFeasTol = 1e-7
    m.Params.OutputFlag = 1 if verbose else 0
    m.Params.Cutoff = 2 * lp_obj
    
    ########################
    # Add tree-based or DAG constraints
    ########################
    
    # first, identify roots
    #   required: abide by sketch
    #   preferred: min-weight distance
    roots = list()
    for j in range(k):
        root = None
        for i in AG.predecessors(j):
            
            g = G.nodes[i]['GEOID20'][0:5]
            if j not in county_geoid_support[g]:
                continue

            # no root yet?
            if root is None:
                root = i
                continue

            # is it closer than a known root?
            dist_i = sq_eucl_dist( G.nodes[i]['X'], G.nodes[i]['Y'], means_x[j], means_y[j] )
            dist_r = sq_eucl_dist( G.nodes[root]['X'], G.nodes[root]['Y'], means_x[j], means_y[j] )
            
            if dist_i < dist_r:
                root = i
                
        roots.append(root)

    # now add DAG constraints (using a special county-weighted distance)
    DG = nx.DiGraph(G)
    c = [ None for j in range(k) ]
    for j in range(k):
        root = roots[j]
        order = find_ordering( DG.subgraph( AG.predecessors(j) ), root, level )
        pos = construct_position( order )
        c[j] = m.addConstrs( x[i,j] <= gp.quicksum( x[p,j] for p in G.neighbors(i) if p in AG.predecessors(j) and pos[p] < pos[i] ) for i in AG.predecessors(j) if i != root )
    
    print("Trying DAG model, with root geoids =",[ G.nodes[i]['GEOID20'] for i in roots ])
    m.optimize()
    if m.solCount > 0:
        return [ [ i for i in AG.predecessors(j) if x[i,j].x > 0.5 ] for j in range(k) ]
    m.remove(c)

    #####################################
    # Add a,b-separator constraints
    #####################################
    
    m.Params.LazyConstraints = 1
    m._callback = AG_contiguity_callback
    m._AG = AG
    m._DG = DG
    m._x = x
    m._parts = k

    print("Trying a,b-separator model.")
    m.optimize(m._callback)
    if m.solCount > 0:
        return [ [ i for i in AG.predecessors(j) if m._x[i,j].x > 0.5 ] for j in range(k) ]
    return None

def interior(G, S):
    return S - set( nx.node_boundary(G, set(G.nodes) - S) )

def sq_eucl_dist(x1,y1,x2,y2):
    return (x1-x2)**2 + (y1-y2)**2

def find_ordering(DG, source, level):
    if level in { 'tract', 'vtd' }:
        n = DG.number_of_nodes()
        for (i,j) in DG.edges:
            gi = DG.nodes[i]['GEOID20'][0:5]
            gj = DG.nodes[j]['GEOID20'][0:5]
            if gi == gj:  # endpoints belong to same county
                DG[i][j]['weight'] = 1
            else:         # endpoints belong to different counties
                DG[i][j]['weight'] = n
    elif level == 'block':
        n = DG.number_of_nodes()
        n_squared = n * n
        for (i,j) in DG.edges:
            gi = DG.nodes[i]['GEOID20']
            gj = DG.nodes[j]['GEOID20']
            if gi[0:11] == gj[0:11]: # endpoints belong to same tract
                DG[i][j]['weight'] = 1
            elif gi[0:5] == gj[0:5]:  # endpoints belong to different tracts, but same county
                DG[i][j]['weight'] = n 
            else: # endpoints belong to different counties
                DG[i][j]['weight'] = n_squared
    else:
        for (i,j) in DG.edges:
            DG[i][j]['weight'] = 1
            
    dist = nx.shortest_path_length(DG,source=source,weight='weight')
    vertex_dist_pairs = [ (i,dist[i]) for i in dist.keys() ]
    vertex_dist_pairs.sort(key=sort_by_second)
    return [ v for (v,dist) in vertex_dist_pairs ]

def AG_contiguity_callback(m, where):

    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve the LP solution
    xval = m.cbGetSolution(m._x)
    DG = m._DG
    AG = m._AG

    # check if each district is connected
    for j in range(m._parts):

        # which nodes are selected?
        S = list( { i for i in AG.predecessors(j) if xval[i,j] > 0.5 } )
        if len(S)==1 or nx.is_strongly_connected( DG.subgraph(S) ):
            continue

        # for smaller connected components, find some vertex 'a' and add cut.
        max_comp = None
        
        for comp in sorted(nx.strongly_connected_components( DG.subgraph(S) ), key=len, reverse=True):
            if max_comp == None:
                max_comp = comp
                b = nx.utils.arbitrary_element( comp )
                continue

            # pick some vertex from this component
            a = nx.utils.arbitrary_element( comp )

            # add a,b-separator inequality
            C = find_fischetti_separator( DG, comp, b)
            m.cbLazy( m._x[a,j] + m._x[b,j] <= 1 + gp.quicksum( m._x[c,j] for c in C if c in AG.predecessors(j)) )
    