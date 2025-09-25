import gurobipy as gp
from gurobipy import GRB
import time

from coarsen import graph_coarsen_by_county
from utils import printif, number_of_whole_counties, check_plan
from initialize import get_initial_constraints
from csd import cluster_sketch_detail

def solve_max_whole(G, L, U, k, incumbent_plan, verbose=True, time_limit=7200):
    '''
    Solves max whole districting problem in graph G 
    :param G: granular graph (e.g., precinct-level, tract-level, or block-level)
    :param L: smallest allowable district population
    :param U: largest allowable district population
    :param k: number of districts to draw
    :param incumbent_plan: an initial districting plan (stored as list of lists)
    :return: a dictionary called `results'
    '''
    start_time = time.perf_counter()
    incomplete_plans = list() # used to store partial plans, which may have multi-districts
    
    # process incumbent_plan 
    check_plan(G, L, U, k, incumbent_plan)
    print(f"Successfully processed plan with { number_of_whole_counties(G, incumbent_plan) } whole counties.")
    
    # create county-level graph 
    (GC, partition) = graph_coarsen_by_county(G, separate_county_components=False, verbose=verbose)

    # get initial set family and fail set
    I_family = get_initial_constraints(GC, L, U, k, verbose=verbose)
    Fail = list()

    # used to cache the detailed districts obtained in cluster_sketch_detail (by cluster)
    csd_cache = dict()

    #   1. solve main problem relaxation, and then
    #   2. add violated combinatorial Benders cut (if any).
    while number_of_whole_counties(G, incumbent_plan) < obj(GC, I_family + Fail) and time.perf_counter() - start_time < time_limit:

        W = soln(GC, I_family + Fail)
        S_fips = [ GC.nodes[c]['GEOID20'] for c in GC.nodes if c not in W ]
        print(f"Seeking a plan for G_W, where |W| = {len(W)} and S = {S_fips}.")

        (plan, incomplete_plan) = cluster_sketch_detail(G, GC, W, L, U, k, incumbent_plan=incumbent_plan, cache=csd_cache, verbose=verbose)
        if incomplete_plan is not None:
            incomplete_plans.append( incomplete_plan )

        if plan == 'infeasible':
            I_family.append( [ GC.nodes[c]['GEOID20'] for c in W ] )
            printif(verbose and plan=='infeasible', "Infeasible. Added combinatorial Benders cut.")
        elif plan == 'inconclusive':
            Fail.append( [ GC.nodes[c]['GEOID20'] for c in W ] )
            printif(verbose and plan=='inconclusive', "Inconclusive. Added (possibly invalid) combinatorial Benders cut.")
        else:
            check_plan(G, L, U, k, plan)
            incumbent_plan = plan
            print(f"Feasible. Found a plan with |W| = { number_of_whole_counties(G, incumbent_plan) } whole counties!")
                                     
    results = dict()
    results['L'] = L
    results['U'] = U
    results['k'] = k
    results['|C|'] = GC.number_of_nodes()
    results['LB'] = number_of_whole_counties(G, incumbent_plan)
    results['UB'] = obj(GC, I_family)
    results['I_family'] = I_family
    results['Fail'] = Fail
    results['incomplete_plans'] = incomplete_plans
    results['incumbent_plan'] = incumbent_plan
    return results

def obj(GC, I_family, verbose=False):
    return len( soln(GC, I_family, verbose=verbose) )

def soln(GC, I_family, verbose=True):
    m = gp.Model()
    m.Params.OutputFlag = 1 if verbose else 0
    w = m.addVars(GC.nodes, vtype=GRB.BINARY)

    # Add constraints
    gtv = { GC.nodes[c]['GEOID20'] : c for c in GC.nodes }
    #print("I_family =",I_family)
    for I in I_family:
        m.addConstr( sum( w[gtv[g]] for g in I ) <= len(I)-1  )
    
    # Primary objective is to maximize number of whole counties.
    # The search is *guided by* secondary objective of maximizing the population of the split counties. 
    m.ModelSense = GRB.MAXIMIZE
    m.setObjectiveN( w.sum(), index=0, priority=1)
    m.setObjectiveN( sum( GC.nodes[c]['TOTPOP'] * (1-w[c]) for c in GC.nodes ), index=1, priority=0 )

    m.optimize()
    assert m.solCount > 0, 'Model is infeasible'
    W = [ c for c in GC.nodes if w[c].x > 0.5 ]
    return W

def get_max_whole_UBs(GC, L, U, k, verbose=True):
    '''
    Gets various upper bounds for the max-whole problem in graph G
    :param GC: county-level graph 
    :param L: smallest allowable district population
    :param U: largest allowable district population
    :param k: number of districts to draw
    :return: a dictionary called `results'
    '''
    assert GC.number_of_nodes() <= 254, "Must pass county-level graph as input."

    results = dict()
    results['k'] = k
    results['L'] = L
    results['U'] = U
    results['|C|'] = GC.number_of_nodes()
    results['m'] = GC.number_of_edges()
    
    # trivial UB
    results['trivial_UB'] = sum( 1 for i in GC.nodes if GC.nodes[i]['TOTPOP'] <= U )
    
    # vicinity UB
    results['vicinity_constraints'] = get_initial_constraints(GC, L, U, k, vi=True, gvi=False, ai=False, gai=False, verbose=verbose)
    results['vicinity_UB'] = obj( GC, results['vicinity_constraints'] )

    # vicinity + generalized vicinity UB
    results['generalized_vicinity_constraints'] = get_initial_constraints(GC, L, U, k, vi=True, gvi=True, ai=False, gai=False, verbose=verbose)
    results['generalized_vicinity_UB'] = obj( GC, results['generalized_vicinity_constraints'] )

    # vicinity + generalized vicinity + articulation UB
    results['articulation_constraints'] = get_initial_constraints(GC, L, U, k, vi=True, gvi=True, ai=True, gai=False, verbose=verbose)
    results['articulation_UB'] = obj( GC, results['articulation_constraints'] )

    # vicinity + generalized vicinity + articulation + vertex cut UB
    results['generalized_articulation_constraints'] = get_initial_constraints(GC, L, U, k, vi=True, gvi=True, ai=True, gai=True, verbose=verbose)
    results['generalized_articulation_UB'] = obj( GC, results['generalized_articulation_constraints'] )

    return results
