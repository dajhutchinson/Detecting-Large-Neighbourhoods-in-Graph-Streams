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

| Table Name              | Type               | # Edges    | # Vertices | Max Degree | Size         |
| ----------------------- | ------------------ | ---------: | ---------: | ---------: | -----------: |
| facebook_small          | Insertion-Only     | 292        | 52         | 36         | 3 KB         |
| facebook_small_deletion | Insertion-Deletion | 258        | 52         | 33         | 5 KB         |
| facebook                | Insertion-Only     | 60,050     | 747        | 586        | 587 KB       |
| facebook_deletion       | Insertion-Deletion | 53,457     | 747        | 522        | 846 KB       |
| gplus                   | Insertion-Only     | 1,179,613  | 12,417     | 5,948      | 52,839 KB    |
| gplus_deletion          | Insertion-Deletion | 1,049,309  | 12,417     | 4,998      | 60,124 KB    |
| gplus_large             | Insertion-Only     | 30,238,035 | 102,100    | 104,947    | 1,328,820 KB |
**NOTE** - For *Insertion-Deletion Graphs*: '# Edges' refers to the final number of edges in the graph after a full stream; '# Vertices' refers to the number of different vertices encountered during the whole stream so some vertices may have 0-degree after the stream has finished.

# TODO
Find/Create some insertion-deletion graphs.
