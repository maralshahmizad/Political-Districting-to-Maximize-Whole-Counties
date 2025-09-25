import networkx as nx
import matplotlib.pyplot as plt

from utils import check_plan, printif, induced_county_clustering, get_whole_fips
from coarsen import graph_coarsen_by_tract
from cluster_for_max_whole import cluster_for_max_whole
from sketch_for_max_whole import sketch
from detail_for_max_whole import detail

# returns `infeasible', 'inconclusive', or a plan (as list of lists)
def cluster_sketch_detail(G, GC, W, L, U, k, incumbent_plan=None, cache=dict(), verbose=False):

    # key: (cluster_fips, whole_fips, L, U, size)
    # value: districts as list of lists
    print("after entering cluster_sketch_detail, len(cache) = ",len(cache))
    
    # get (sketch-feasible) clustering
    (clusters_start, sizes_start) = (None, None) if incumbent_plan is None else induced_county_clustering(G, incumbent_plan)
    cs = cluster_for_max_whole(GC, L, U, k, whole_counties=W, sketch_time_limit=60, clusters_start=clusters_start, sizes_start=sizes_start, verbose=verbose)
    if cs in ['infeasible', 'inconclusive']:
        return (cs, None)

    # add incumbent plan to cache?
    if incumbent_plan is not None:
        add_to_cache(cache, G, GC, L, U, incumbent_plan, clusters_start, sizes_start)
        print("after adding incumbent_plan, len(cache) = ",len(cache))

    # sketch and detail each cluster
    plan = list()
    (clusters, sizes) = cs
    assert k == sum( sizes )
    for p in range(len(sizes)):   

        # get cluster
        cluster = clusters[p]
        size = sizes[p]
        cluster_fips = { GC.nodes[i]['GEOID20'] for i in cluster }
        whole = [ c for c in W if c in cluster ]
        whole_fips = { GC.nodes[c]['GEOID20'] for c in W if c in cluster }

        # check that cluster is connected and has appropriate population
        vertices = [ i for i in G.nodes if G.nodes[i]['GEOID20'][0:5] in cluster_fips ]
        assert nx.is_connected( G.subgraph(vertices) )
        population = sum( G.nodes[i]['TOTPOP'] for i in vertices )
        assert L * size <= population and population <= U * size 

        print(f"Seeking {size} districts for cluster #{p+1}/{len(sizes)} over n = {len(vertices)} vertices.")
        print(f"This cluster has counties {cluster_fips} of which {whole_fips} are whole.")
        print(f"Compare [L,U] = {[L,U]} to p(cluster)/size = {round(population/size,2)}.")

        # check if we've already handled this cluster before
        key = ( frozenset(cluster_fips), frozenset(whole_fips), L, U, size )
        if key in cache:
            plan += cache[key]
            print("Used cached districts.")
            continue

        # trivial case (don't bother caching this...)
        if size == 1:
            plan.append( vertices )
            continue
                    
        # get sketch
        GS = GC.subgraph( cluster )
        sketch_soln = sketch( GS, L, U, size, whole_counties=whole, time_limit=60, verbose=verbose )

        if sketch_soln == 'infeasible':
            assert False, "By design, the clusters should be sketch feasible."
        elif sketch_soln == 'inconclusive':
            return ('inconclusive', None)
        else:
            (xsoln, ysoln, zsoln) = sketch_soln

        # get detailed districts
        districts = None
        sketch_support = { GS.nodes[i]['GEOID20'] : list( xsoln[i].keys() ) for i in GS.nodes }
        
        # if G is block-level, then first try to district at tract level (60 sec time limit)
        v = vertices[0]
        if len( G.nodes[v]['GEOID20'] ) >= 15:
            (G_tract, tract_buckets) = graph_coarsen_by_tract( G.subgraph(vertices), verbose=False )
            assert nx.is_connected( G_tract ), "Tract-level subgraph is not connected. Consider adding edges to input graph."
            print(f"First, trying to draw districts at tract-level where n = {len(tract_buckets)}.")
            tract_districts = detail( G_tract, L, U, size, sketch_support, time_limit=60, verbose=verbose )
            if tract_districts is None:
                print("Tract-level attempt failed. Returning to original graph (block-level).")
            else:
                districts = [ [ v for t in tract_district for v in tract_buckets[t] ] for tract_district in tract_districts ]

        # second, try to district at original level (block or precinct)
        if districts is None:
            assert nx.is_connected( G.subgraph(vertices) ), "Original-level subgraph is not connected. Consider adding edges to input graph."
            districts = detail( G.subgraph(vertices), L, U, size, sketch_support, verbose=verbose )

        # successful?
        if districts is not None:
            assert len(districts) == size
            check_plan(G.subgraph(vertices), L, U, size, districts)
            plan += districts
            cache[key] = districts
        else:
            plan.append( vertices )
            ideal = sum( GS.nodes[i]['TOTPOP'] for i in GS.nodes ) / size
            print(f"Failed to detail cluster = {[GS.nodes[i]['GEOID20'] for i in GS.nodes]}. Keeping as multi-district of size = {size}.")
            print(f"Compare [L, U] = [{L},{U}] to p(cluster)/size = {ideal}.")
            node_pos = { i : ( GS.nodes[i]['X'], GS.nodes[i]['Y'] ) for i in GS.nodes }
            node_colors = [ 'green' if i in whole else 'red' for i in GS.nodes ]
            my_labels = { i : GS.nodes[i]['TOTPOP'] for i in GS.nodes }
            nx.draw(GS, pos=node_pos, node_color=node_colors, with_labels=True, labels=my_labels, font_color="black")
            plt.show()

    # return ordered pair of the form: (complete_plan, incomplete_plan)
    if len(plan) == k:
        return (plan, None)
    elif len(plan) < k:
        return ('inconclusive', plan)
    else:
        assert False, "We should not have len(plan) > k."

def add_to_cache(cache, G, GC, L, U, plan, clusters, sizes):
    all_whole_fips = get_whole_fips(G, plan)
    for p in range(len(sizes)):
        whole_fips = { f for f in all_whole_fips if f in clusters[p] }
        key = ( frozenset(clusters[p]), frozenset(whole_fips), L, U, sizes[p] )
        cache[key] = [ district for district in plan if G.nodes[district[0]]['GEOID20'][0:5] in clusters[p] ]
    return
