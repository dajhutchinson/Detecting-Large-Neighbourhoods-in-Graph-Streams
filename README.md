# Detecting-Large-Neighbourhoods-in-Graph-Streams
Streaming Frequent Items with Timestamps and Detecting Large Neighbourhoods in Graph Streams

## Compiling
`g++ read.cpp -o read`
`g++ insertionDeletionStreams.cpp lib/MurmurHash3.cpp -o ids`

## Execution
`read`
`ids`

## Data
Graphs are stored as a stream of edges, in no particular order.
Each line represents an edge with a space separating each node id.

| Table Name     | Type                  | # Edges   | # Vertices | Max Degree |
| -------------- | --------------------- | --------: | ---------: | ---------: |
| facebook_small | Undirected, Insertion | 292       | 52         | 36         |
| facebook       | Undirected, Insertion | 60,050    | 747        | 586        |
| gplus          | Undirected, Insertion | 1,179,613 | 12,417     | 5,948      |

# TODO
Find/Create some insertion-deletion graphs.
