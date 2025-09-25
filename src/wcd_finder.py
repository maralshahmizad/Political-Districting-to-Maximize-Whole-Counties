import gurobipy as gp
from gurobipy import GRB
import networkx as nx
from math import ceil

# whole-county district finder
def wcd(DG, L, U, k, c, forbidden_districts=list(), verbose=False):

    m = gp.Model()
    if not verbose:
        m.Params.OutputFlag = 0 # stop Gurobi from printing logs

    # Variables
    x = m.addVars( DG.nodes, vtype=GRB.BINARY )
    f = m.addVars( DG.edges )

    # root the district at c
    x[c].LB = 1

    # use MOI compactness objective to promote contiguous solutions
    # (divide by 1,000 and round to integer for better numerics)
    dist = nx.single_source_shortest_path_length(DG, source=c)
    coef = { i : ceil( dist[i] * dist[i] * DG.nodes[i]['TOTPOP'] / 1000 ) for i in DG.nodes }
    m.setObjective( sum( coef[i] * x[i] for i in DG.nodes ), GRB.MINIMIZE )

    # population between L and U
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i] for i in DG.nodes ) >= L )
    m.addConstr( gp.quicksum( DG.nodes[i]['TOTPOP'] * x[i] for i in DG.nodes ) <= U )

    # contiguity constraints: send one unit of flow to selected nodes
    m.addConstrs( gp.quicksum( f[j,i] - f[i,j] for j in DG.neighbors(i) ) == x[i] for i in DG.nodes if i != c )

    # contiguity constraints: accept inflow only at selected nodes, but not at c.
    M = DG.number_of_nodes() - 1
    m.addConstrs( gp.quicksum( f[j,i] for j in DG.neighbors(i) ) <= M * x[i] for i in DG.nodes if i != c )
    m.addConstr( gp.quicksum( f[j,c] for j in DG.neighbors(c) ) == 0 )

    # must differ from previous/forbidden districts
    for D in forbidden_districts:
        VD = [ i for i in DG.nodes if i not in D ]
        m.addConstr( gp.quicksum( x[i] for i in VD ) + gp.quicksum( (1-x[i]) for i in D ) >= 1 )
    
    # callback info
    m._x = x
    m._DG = DG
    m._L = L
    m._U = U
    m._k = k
    m.Params.LazyConstraints = 1
    m._callback = wcd_callback
    
    # solve
    m.Params.MIPGap = 100    # allow large gap, because we care about feasibility, not optimality
    m.optimize(m._callback)

    if m.status == GRB.INFEASIBLE:
        return None
    else:
        return [ i for i in DG.nodes if m._x[i].x > 0.5 ]


def wcd_callback(m, where):
    
    # check if LP relaxation at this BB node is integer
    if where != GRB.Callback.MIPSOL: 
        return
        
    # retrieve data
    xval = m.cbGetSolution(m._x)
    DG = m._DG
    L = m._L
    U = m._U
    k = m._k
    
    # does each component of G-D have a population 
    #  within [L,U] or an integer multiple thereof?
    unselected = [ i for i in DG.nodes if xval[i] < 0.5 ]
    for component in nx.strongly_connected_components( DG.subgraph(unselected) ):
        population = sum( DG.nodes[i]['TOTPOP'] for i in component )
        assert population > 0, "wcd_callback() works with county-level graphs, so no component should have population <= 0."
        q = ceil( population / U ) # smallest q satisfying population <= q*U
        if q*L > population:
            # if no q with Lq <= population <= Uq, then either:
            #   1. something from component needs to be picked (by x), or 
            #   2. one of its neighbors shouldn't be picked.
            neighbors = { j for i in component for j in DG.neighbors(i) if j not in component }
            m.cbLazy( gp.quicksum( m._x[i] for i in component ) + gp.quicksum( (1-m._x[i]) for i in neighbors ) >= 1 )
