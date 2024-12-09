import math
from pyproj import Proj
from number_of_districts import *
import csv

def update_attributes(G, state):
    
    code = get_epsg(state)          
    p = Proj("EPSG:"+code, preserve_units=True)

    # update lat, long, and pop
    for i in G.nodes:
        G.nodes[i]['C_X'] = float(G.nodes[i]['INTPTLON20'])
        G.nodes[i]['C_Y'] = float(G.nodes[i]['INTPTLAT20'])
        G.nodes[i]['TOTPOP'] = int(G.nodes[i]['P0010001'])
        
        # projection to get (x,y) coordinates (units of KM)
        (x,y) = p( G.nodes[i]['C_X'], G.nodes[i]['C_Y'] )
        G.nodes[i]['X'] = x / 1000
        G.nodes[i]['Y'] = y / 1000  

    return
        
        
def get_k_L_U(GC, state, district_type):
    
    if district_type == 'CD':
        deviation = 0.01
        k = congressional_districts_2020[state]
    elif district_type == 'SS':
        deviation = 0.10
        k = state_senate_districts_2020[state]
    elif district_type == 'SH':
        deviation = 0.10
        k = state_house_districts_2020[state]
    else:
        print("ERROR: option district_type =",district_type,"not supported.")

    total_population = sum( GC.nodes[i]['TOTPOP'] for i in GC.nodes )
    if k==0:  # case: NE SH
        ideal_population = 0
    else:
        ideal_population = total_population / k
    L = math.ceil( ( 1 - deviation / 2 ) * ideal_population )
    U = math.floor( ( 1 + deviation / 2 ) * ideal_population )
    
    print("Starting",state,"with k =",k,"and deviation =",deviation)
    print("Thus, we have L =",L,"and U =",U)
    return (k, L, U)


states = {
'01': {'abbr': 'AL', 'epsg': '3465', 'name': 'Alabama'},
'02': {'abbr': 'AK', 'epsg': '3471', 'name': 'Alaska'},
'04': {'abbr': 'AZ', 'epsg': '3478', 'name': 'Arizona'},
'05': {'abbr': 'AR', 'epsg': '3484', 'name': 'Arkansas'},
'06': {'abbr': 'CA', 'epsg': '3493', 'name': 'California'},
'08': {'abbr': 'CO', 'epsg': '3501', 'name': 'Colorado'},
'09': {'abbr': 'CT', 'epsg': '3507', 'name': 'Connecticut'},
'10': {'abbr': 'DE', 'epsg': '3509', 'name': 'Delaware'},
'12': {'abbr': 'FL', 'epsg': '3514', 'name': 'Florida'},
'13': {'abbr': 'GA', 'epsg': '3518', 'name': 'Georgia'},
'15': {'abbr': 'HI', 'epsg': '2784', 'name': 'Hawaii'},
'16': {'abbr': 'ID', 'epsg': '3524', 'name': 'Idaho'},
'17': {'abbr': 'IL', 'epsg': '3528', 'name': 'Illinois'},
'18': {'abbr': 'IN', 'epsg': '3532', 'name': 'Indiana'},
'19': {'abbr': 'IA', 'epsg': '3536', 'name': 'Iowa'},
'20': {'abbr': 'KS', 'epsg': '3540', 'name': 'Kansas'},
'21': {'abbr': 'KY', 'epsg': '3544', 'name': 'Kentucky'},
'22': {'abbr': 'LA', 'epsg': '3550', 'name': 'Louisiana'},
'23': {'abbr': 'ME', 'epsg': '3557', 'name': 'Maine'},
'24': {'abbr': 'MD', 'epsg': '3559', 'name': 'Maryland'},
'25': {'abbr': 'MA', 'epsg': '3585', 'name': 'Massachusetts'},
'26': {'abbr': 'MI', 'epsg': '3587', 'name': 'Michigan'},
'27': {'abbr': 'MN', 'epsg': '3594', 'name': 'Minnesota'},
'28': {'abbr': 'MS', 'epsg': '3597', 'name': 'Mississippi'},
'29': {'abbr': 'MO', 'epsg': '3602', 'name': 'Missouri'},
'30': {'abbr': 'MT', 'epsg': '3604', 'name': 'Montana'},
'31': {'abbr': 'NE', 'epsg': '3606', 'name': 'Nebraska'},
'32': {'abbr': 'NV', 'epsg': '3607', 'name': 'Nevada'},
'33': {'abbr': 'NH', 'epsg': '3613', 'name': 'NewHampshire'},
'34': {'abbr': 'NJ', 'epsg': '3615', 'name': 'NewJersey'},
'35': {'abbr': 'NM', 'epsg': '3617', 'name': 'NewMexico'},
'36': {'abbr': 'NY', 'epsg': '3623', 'name': 'NewYork'},
'37': {'abbr': 'NC', 'epsg': '3631', 'name': 'NorthCarolina'},
'38': {'abbr': 'ND', 'epsg': '3633', 'name': 'NorthDakota'},
'39': {'abbr': 'OH', 'epsg': '3637', 'name': 'Ohio'},
'40': {'abbr': 'OK', 'epsg': '3639', 'name': 'Oklahoma'},
'41': {'abbr': 'OR', 'epsg': '3645', 'name': 'Oregon'},
'42': {'abbr': 'PA', 'epsg': '3649', 'name': 'Pennsylvania'},
'44': {'abbr': 'RI', 'epsg': '3653', 'name': 'RhodeIsland'},
'45': {'abbr': 'SC', 'epsg': '3655', 'name': 'SouthCarolina'},
'46': {'abbr': 'SD', 'epsg': '3657', 'name': 'SouthDakota'},
'47': {'abbr': 'TN', 'epsg': '3661', 'name': 'Tennessee'},
'48': {'abbr': 'TX', 'epsg': '3669', 'name': 'Texas'},
'49': {'abbr': 'UT', 'epsg': '3675', 'name': 'Utah'},
'50': {'abbr': 'VT', 'epsg': '3684', 'name': 'Vermont'},
'51': {'abbr': 'VA', 'epsg': '3685', 'name': 'Virginia'},
'53': {'abbr': 'WA', 'epsg': '3689', 'name': 'Washington'},
'54': {'abbr': 'WV', 'epsg': '3693', 'name': 'WestVirginia'},
'55': {'abbr': 'WI', 'epsg': '3695', 'name': 'Wisconsin'},
'56': {'abbr': 'WY', 'epsg': '3703', 'name': 'Wyoming'}
}


def get_epsg(state):
    for code in states.keys():
        if states[code]['abbr'] == state:
            return states[code]['epsg']
    assert False # this shouldn't happen
    return None

def get_state(fips):
    assert fips in states.keys()
    return states[fips]['abbr']
    
def get_fips(state):
    for code in states.keys():
        if states[code]['abbr'] == state:
            return code
    assert False # this shouldn't happen
    return None

def number_of_county_splits(districts_geoid):
    county_geoids = { geoid[0:5] for district in districts_geoid for geoid in district }
    county_assignment = { county_geoid : list() for county_geoid in county_geoids }
    for j in range(len(districts_geoid)):
        for geoid in districts_geoid[j]:
            county_geoid = geoid[0:5]
            if j not in county_assignment[county_geoid]:
                county_assignment[county_geoid].append(j)
    return sum( len(county_assignment[cg]) - 1 for cg in county_assignment.keys() )

def export_to_baf(GB, geoid_assignment, outfilename):
    
    with open(outfilename, 'w') as csvfile: 
        csvwriter = csv.writer(csvfile) 
        for i in GB.nodes:
            geoid = GB.nodes[i]['GEOID20']
            county = geoid[0:5]
            tract = geoid[0:11]
            block = geoid
            if county in geoid_assignment:
                row = [geoid, geoid_assignment[county]+1]
            elif tract in geoid_assignment:
                row = [geoid, geoid_assignment[tract]+1]
            else:
                row = [geoid, geoid_assignment[block]+1]
            csvwriter.writerow(row) 
    
    return
    