# File Formats
This document describes the file formats used by this code.

## JSON to store graphs
Encodes the geographic adjacency graph of a state at the county, tract, or block level.
Each vertex represents a geographic unit. An edge connects two units if they share a border of positive length. 
Each vertex has attributes like population ('TOTPOP') and [geographic identifier or GEOID](https://www.census.gov/programs-surveys/geography/guidance/geo-identifiers.html) ('GEOID20').

Our file naming convention is
```
{STATE}_{LEVEL}.json
```
and uses [2-letter postal abbreviations](https://en.wikipedia.org/wiki/List_of_U.S._state_and_territory_abbreviations#Postal_codes). For example, Iowa's tract-level graph is stored in IA_tract.json. 

JSON files are read using [NetworkX](https://networkx.org/documentation/stable/reference/readwrite/json_graph.html), which appears in our [utils.py code as read_graph_from_json()](https://github.com/maralshahmizad/Political-Districting-to-Maximize-Whole-Counties/blob/367c43005faa93110c4b9a2dc62b8fd5b8ddb986/src/utils.py#L10).

## BAF to store plans
Districting plans are stored as [block-assignment files](https://www.census.gov/geographies/reference-files/time-series/geo/block-assignment-files.html) (.baf), which are essentially just CSV files.
Each row gives a census block's GEOID and its district's number.
BAFs can be imported into [Dave's Redistricting App](https://davesredistricting.org/) for further analysis.

## GIS Shapefiles
Geographic data is stored as [shapefiles](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html). 

Our file naming convention is
```
{STATE}_{LEVEL}.shp
```

It is common to have multiple files with different extensions (.shp, .dbf, .prj, .cpg, .shx), which may 
all be needed even when the code only refers to the file with .shp file extension.
