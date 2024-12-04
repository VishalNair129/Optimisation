//VISHAL NAIR  NIU - 1740105
//NICO WUSTNER NIU - 1736863


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>  // Include this header for ULONG_MAX
#include <math.h>
#include <float.h>
#include <stdbool.h>

// Define the structure for nodes in the graph
typedef struct {
    unsigned long id; // Node identification
    double lat, lon; // Node position
    unsigned short nsucc; // Number of node successors; i.e., length of successors
    unsigned long *successors; // Successor nodes (dynamically allocated)
    char *name; // Node name (dynamically allocated)

} node;

// Define the structure for A* pathfinding data
typedef struct {
    double gscore;        // The cost from the start node
    double fscore;       // The estimated cost to the goal node (f = g + h)
    unsigned long parent; // Parent node in the path
    bool isOpen;       // boolean to check if a node is in the Queue or not
} AstarPath;

// Define the structure for vertices in the priority queue
typedef struct {
    unsigned id;
    double value;
} Vertex;

// Defines a node in the ordered queue (a linked list node) which contains:
// a pointer to a Vertex and
// a pointer to the next element.
typedef struct OQueueElement {
    Vertex *vertex;
    struct OQueueElement *next;
} OQueueElement;

// Defines the structure of an ordered queue:
// only keeps a pointer to the first element
typedef struct {
    OQueueElement *first;
} OQueue;

//The functions are defined in the end but we need to declare the functions used within main
int Oenqueue(Vertex *vert, OQueue *Q);
Vertex Odequeue(OQueue *Q);
int Orequeue(Vertex *vert, OQueue *Q);
double haversine_distance(double lat1, double lon1, double lat2, double lon2);
unsigned long searchNode(unsigned long id, node *nodes, unsigned long nnodes);
bool astar(node *nodes, unsigned long nnodes, unsigned long start_id, unsigned long goal_id,AstarPath *PathData);
void ExitError(const char *miss, int errcode);

int main(int argc, char *argv[]) {
    clock_t start_time;
    unsigned long nnodes;

    start_time = clock();

    char binmapname[80];
    strcpy(binmapname, "andorra.csv.bin"); // Default binary map file name

    if (argc < 3) {
        printf("Usage: %s <start_id> <goal_id>\n", argv[0]);
        return 1;
    }

    unsigned long start_id = strtoul(argv[1], NULL, 10); // Convert start_id from argument
    unsigned long goal_id = strtoul(argv[2], NULL, 10);  // Convert goal_id from argument
    if (argc > 3) strcpy(binmapname, argv[3]); // Use the provided binary map file name if given

    FILE *binmapfile;
    start_time = clock();

    binmapfile = fopen(binmapname, "rb");// Open the binary file for reading
    if (binmapfile == NULL) {
        printf("Error when opening the binary file\n");
        return 1;
    }

    fread(&nnodes, sizeof(unsigned long), 1, binmapfile); // Read the number of nodes from the binary file

    node *nodes;

    nodes = (node *)malloc(nnodes * sizeof(node));
    if (nodes == NULL) {
        printf("Error when allocating the memory for the nodes\n");
        return 2;
    }

    // Read the node data from the binary file
     for (unsigned long i = 0; i < nnodes; i++) {
       if (fread(&nodes[i], sizeof(node) - sizeof(unsigned long *) - sizeof(char *), 1, binmapfile) != 1) {
        ExitError("Error when reading node metadata from the binary file", 33);
        }
        // Allocate memory for successors
        nodes[i].successors = (unsigned long *)malloc(nodes[i].nsucc * sizeof(unsigned long));
        if (nodes[i].successors == NULL) {
            ExitError("Error when allocating memory for successors", 34);
        }

        // Read successors
        if (fread(nodes[i].successors, sizeof(unsigned long), nodes[i].nsucc, binmapfile) != nodes[i].nsucc) {
            ExitError("Error when reading node successors from the binary file", 35);
            }

        // Read the length of the name
        unsigned long name_length;
        if (fread(&name_length, sizeof(unsigned long), 1, binmapfile) != 1) {
            ExitError("Error when reading name length from the binary file", 36);
            }


    
        // Allocate memory for the name
        nodes[i].name = (char *)malloc(name_length * sizeof(char));
        if (nodes[i].name == NULL) {
            ExitError("Error when allocating memory for name", 38);
            }

         // Read the name itself
        if (fread(nodes[i].name, sizeof(char), name_length, binmapfile) != name_length) {
            ExitError("Error when reading name from the binary file", 39);
            }
    }

    fclose(binmapfile);// Close the binary file

    printf("Total number of nodes is %ld\n", nnodes);
    printf("Elapsed time: %f seconds\n", (float)(clock() - start_time) / CLOCKS_PER_SEC);

    // Optionally print a node with more than 4 successors to verify correct reading
    unsigned long index = 0;
    for (unsigned long i = 0; i < nnodes; i++) {
        if (nodes[i].nsucc > 4) {
            index = i;
            break;
        }
    }
    printf("Node %lu has id=%lu and %u successors:\n", index, nodes[index].id, nodes[index].nsucc);
    for (int i = 0; i < nodes[index].nsucc; i++) {
        printf("  Node %lu with id %lu.\n", nodes[index].successors[i], nodes[nodes[index].successors[i]].id);
    }

    
    // Run A* algorithm
    AstarPath *PathData = (AstarPath*)malloc(nnodes * sizeof(AstarPath));
    if (PathData == NULL) {
        printf("Error: cannot allocate memory for Astar Path.\n");
        return 1;
    }

    start_time = clock();

    bool result = astar(nodes, nnodes,start_id, goal_id,PathData);
    
    printf("Elapsed finding path time: %f seconds\n", (float)(clock() - start_time) / CLOCKS_PER_SEC);

    if (!result)
        ExitError("no solution found in AStar", 7);
        
    // Trace the path from goal to start and reverse it
    unsigned long start = searchNode(start_id, nodes, nnodes);
    unsigned long goal = searchNode(goal_id, nodes, nnodes);
    
    unsigned long v = goal;
    unsigned long pv = PathData[v].parent;
    unsigned long ppv;
    PathData[goal].parent = ULONG_MAX;

    while(v != start) {
        ppv = PathData[pv].parent;
        PathData[pv].parent = v;
        v = pv;
        pv = ppv;
    }
    
 // Print and write the optimal path and distance to a file named Astarpath.txt
    const char *filename = "Astarpath.txt";
    FILE *file; 
    file=fopen(filename, "w");
    if (file == NULL) {
        printf("Error: cannot open file for writing.\n");
    }

    printf("# Distance from %lu to %lu : %f metres.\n",start_id,goal_id,PathData[goal].gscore);
    fprintf(file, "# Distance from %lu to %lu : %f metres.\n", start_id, goal_id, PathData[goal].gscore);
    printf("# Optimal path:\n");
    fprintf(file, "# Optimal path:\n");
    printf("Id = %lu | %f | %f | Dist = %f\n", nodes[start].id, nodes[start].lat, nodes[start].lon, 0.0);
    fprintf(file, "Id = %lu | %f | %f | Dist = %f\n", nodes[start].id, nodes[start].lat, nodes[start].lon, 0.0);

    // Trace the path from the goal to the start
    for (v = PathData[start].parent; v != ULONG_MAX; v = PathData[v].parent) {
        printf("Id = %lu | %f | %f | Dist = %f\n", nodes[v].id, nodes[v].lat, nodes[v].lon, PathData[v].gscore);
        fprintf(file, "Id = %lu | %f | %f | Dist = %f\n", nodes[v].id, nodes[v].lat, nodes[v].lon, PathData[v].gscore);
        }
    
    
    fclose(file);

    // Free allocated memory
    free(PathData);

    for (unsigned long i = 0; i < nnodes; i++) {
        free(nodes[i].successors);
        free(nodes[i].name); 

    }

    free(nodes);

    return 0;
}


// Enqueue function for priority queue
int Oenqueue(Vertex *vert, OQueue *Q) {
    OQueueElement *element_aux, *element_iter;
    element_aux = (OQueueElement *) malloc(sizeof(OQueueElement));
    if (element_aux == NULL) {
        printf("Error: cannot allocate memory for queue element.");
        return 1;
    }
    element_aux->vertex = vert;
    element_aux->next = NULL;
    if (Q->first == NULL) {
        Q->first = element_aux;
        return 0;
    }
    if (Q->first->vertex->value > element_aux->vertex->value) {
        element_aux->next = Q->first;
        Q->first = element_aux;
        return 0;
    }
    element_iter = Q->first;
    while (element_iter->next != NULL && element_iter->next->vertex->value < vert->value) {
        element_iter = element_iter->next;
    }
    element_aux->next = element_iter->next;
    element_iter->next = element_aux;
    return 0;
}

// Dequeue function for priority queue
Vertex Odequeue(OQueue *Q) {
    Vertex v = {0, 0.0};  // Initialize a default Vertex
    if (Q->first == NULL) { // There is no queue!
        return v;
    }
    OQueueElement *element_aux = Q->first;
    v = *(element_aux->vertex);
    Q->first = Q->first->next;
    free(element_aux);
    return v;
}

// Requeue function for priority queue
int Orequeue(Vertex *vert, OQueue *Q) {
    OQueueElement *element_iter = Q->first;
    OQueueElement *prev_element = NULL;

    // Search for the element with the same id
    while (element_iter != NULL) {
        if (element_iter->vertex->id == vert->id) {
            // Remove the existing element from the queue
            if (prev_element == NULL) {
                Q->first = element_iter->next;
            } else {
                prev_element->next = element_iter->next;
            }
            free(element_iter);
            break;
        }
        prev_element = element_iter;
        element_iter = element_iter->next;
    }
    return Oenqueue(vert, Q);
}


// Function to calculate the Haversine distance between two points on Earth
double haversine_distance(double lat1, double lon1, double lat2, double lon2) {

    const double R = 6371e3;  // Earth radius in meters

    // Convert latitude and longitude from degrees to radians
    double phi1 = lat1 * M_PI / 180.0;
    double phi2 = lat2 * M_PI / 180.0;
    double delta_phi = (lat2 - lat1) * M_PI / 180.0;
    double delta_lambda = (lon2 - lon1) * M_PI / 180.0;

    // Calculate the Haversine formula components
    double a = sin(delta_phi / 2) * sin(delta_phi / 2) +
               cos(phi1) * cos(phi2) *
               sin(delta_lambda / 2) * sin(delta_lambda / 2);

    // Calculate the angular distance in radians    
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));

    // Return the distance in meters
    return R * c;
}

// Function to perform binary search for a node by its ID
unsigned long searchNode(unsigned long id, node *nodes, unsigned long nnodes) {
    unsigned long l = 0, r = nnodes - 1, m;
    while (l <= r) {
        m = l + (r - l) / 2;
        if (nodes[m].id == id) 
        return m;
        if (nodes[m].id < id)
            l = m + 1;
        else
            r = m - 1;
    }
    return nnodes + 1;  // id not found, return nnodes + 1
}

// A* Algorithm Implementation
bool astar(node *nodes, unsigned long nnodes, unsigned long start_id, unsigned long goal_id,AstarPath *PathData) {
    // Find the index of the start and goal nodes using their IDs
    unsigned long start = searchNode(start_id, nodes, nnodes);
    unsigned long goal = searchNode(goal_id, nodes, nnodes);
    
    // Check if the start or goal nodes are not found
    if (start == nnodes + 1 || goal == nnodes + 1) {
        printf("Start or goal node not found.\n");
        return false;
    }

    

    // Initialize PathData for all nodes
    for (unsigned long i = 0; i < nnodes; i++) {
         PathData[i].isOpen = false;
         PathData[i].gscore = DBL_MAX;
         PathData[i].fscore = DBL_MAX;
        
    }
    // Initialize the start node's gscore and fscore
    PathData[start].parent = ULONG_MAX;
    PathData[start].gscore = 0;
    PathData[start].fscore = haversine_distance(nodes[start].lat, nodes[start].lon, nodes[goal].lat, nodes[goal].lon);

    // Initialize the open set  
    OQueue *open_set = (OQueue*)malloc(sizeof(OQueue));
    if (open_set == NULL) {
        printf("Error: cannot allocate memory for open queue.\n");
        return false;
    }
    open_set->first = NULL;
    Vertex start_vertex = {start, PathData[start].fscore};
    // Add the start node to the open set
    Oenqueue(&start_vertex, open_set);
    PathData[start].isOpen = true;

    // Main loop of the A* algorithm
    while (open_set->first != NULL) {
        // Dequeue the node with the lowest fscore
        Vertex current_vertex = Odequeue(open_set);
        unsigned long current = current_vertex.id;
        PathData[current].isOpen = false;

// Check if we have reached the goal
        if (current == goal) {
            return true;
        }
        
        // Explore the neighbors of the current node
        for (unsigned short i = 0; i < nodes[current].nsucc; i++) {
            unsigned long neighbor = nodes[current].successors[i];

            // Calculate the tentative gscore for the neighbor
            double tentative_g_score = PathData[current].gscore + haversine_distance(nodes[current].lat, nodes[current].lon, nodes[neighbor].lat, nodes[neighbor].lon);
                        
            // Check if this path to the neighbor is better
            if (tentative_g_score < PathData[neighbor].gscore) {
                PathData[neighbor].parent = current;

                // Update the fscore
                if (PathData[neighbor].gscore == DBL_MAX) {
                    PathData[neighbor].fscore = tentative_g_score + haversine_distance(nodes[neighbor].lat, nodes[neighbor].lon, nodes[goal].lat, nodes[goal].lon);
                }
                else {
                    PathData[neighbor].fscore = tentative_g_score + PathData[neighbor].fscore - PathData[neighbor].gscore;
                }

                PathData[neighbor].gscore = tentative_g_score;

                // update the cost value or add new node to the queue
               Vertex *neighbor_vertex = (Vertex*)malloc(sizeof(Vertex));  // dynamically allocate memory
               if (neighbor_vertex == NULL) {
                   printf("Error: cannot allocate memory for neighbor vertex.\n");
                   return false;
                }
                neighbor_vertex->id = neighbor;
                neighbor_vertex->value = PathData[neighbor].fscore;

                if (!PathData[neighbor].isOpen) {
                    if (Oenqueue(neighbor_vertex, open_set))
                        return 0;
                    PathData[neighbor].isOpen = 1;
                }
                else {
                    if (Orequeue(neighbor_vertex, open_set))
                        return 0;
                }
            }
            
           }

    }
    free(open_set);
    // Return false if the goal was not reached
    return false;

}

// Function to print error and exit
void ExitError(const char *miss, int errcode) {
    fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
}