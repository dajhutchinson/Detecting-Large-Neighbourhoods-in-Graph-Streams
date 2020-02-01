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

| Table Name              | Type               | # Edges   | # Vertices | Max Degree |
| ----------------------- | ------------------ | --------: | ---------: | ---------: |
| facebook_small          | Insertion-Only     | 292       | 52         | 36         |
| facebook_small_deletion | Insertion-Deletion | 258       | 52         | 33         |
| facebook                | Insertion-Only     | 60,050    | 747        | 586        |
| facebook_deletion       | Insertion-Deletion | 53,457    | 747        | 522        |
| gplus                   | Insertion-Only     | 1,179,613 | 12,417     | 5,948      |
| gplus_deletion          | Insertion-Deletion | 1,049,309 | 12,417     | 4,998      |
**NOTE** - For *Insertion-Deletion Graphs*: '# Edges' refers to the final number of edges in the graph after a full stream; '# Vertices' refers to the number of different vertices encountered during the whole stream so some vertices may have 0-degree after the stream has finished.

# TODO
Find/Create some insertion-deletion graphs.
