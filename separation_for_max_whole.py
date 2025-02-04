import networkx as nx
import gurobipy as gp
from gurobipy import GRB 
from sketch_for_max_whole import sketch

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
    
    C = [ i for i in DG.nodes if neighbors_component[i] and visited[i] ]
    return C


def labeling_contiguity_callback(m, where):

    # check if LP relaxation at this branch-and-bound node has an integer solution
    if where != GRB.Callback.MIPSOL: 
        return 
    # retrieve the LP solution
    xval = m.cbGetSolution(m._x)
    yval = m.cbGetSolution(m._y)
    k = m._k
    L = m._L
    U = m._U
    G = m._G
    DG = m._DG
    sketch_intensity = m._sketch_intensity

    # first, check connectivity
    disconnected = False
    for j in range(m._parts):

        # which nodes are selected?
        S = [ i for i in DG.nodes if xval[i,j] > 0.5 ]
        if len(S)==1 or nx.is_strongly_connected( DG.subgraph(S) ):
            continue

        disconnected = True
        
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
   
    if disconnected:
        return
    
    # If connected, then check sketch feasibility
    if sketch_intensity==0:
        sketch_parts = 0
    elif sketch_intensity==1:
        sketch_parts = 1
    elif sketch_intensity==2:
        sketch_parts = m._parts
    else:
        assert False, "sketch_intensity must be 0, 1, or 2."
        
    for j in range(sketch_parts):
        # which nodes are selected?
        S = [ i for i in DG.nodes if xval[i,j] > 0.5 ]
        GS = G.subgraph(S) 
        #print("whole_counties =",m._whole_counties)
        whole = [ i for i in S if i in m._whole_counties ]
        size = round(yval[j])
        if size == 1:
            continue
        (xsoln, ysoln, zsoln) = sketch(GS, L, U, size, whole_counties=whole, time_limit=60)
        if xsoln is None:
            VS = [ i for i in DG.nodes if i not in S]
            # no-good cut
            for t in range(m._parts):
                m.cbLazy( gp.quicksum( (1-m._x[v,t]) for v in S ) + gp.quicksum( m._x[v,t] for v in VS ) >= 1 )