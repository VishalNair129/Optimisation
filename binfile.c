//VISHAL NAIR  NIU - 1740105
//NICO WUSTNER NIU - 1736863


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

// Define the node structure
typedef struct {
    unsigned long id; // Node identification
    double lat, lon; // Node position
    unsigned short nsucc; // Number of node successors; i.e., length of successors
    unsigned long *successors; // Successor nodes (dynamically allocated)
    char *name; // Node name (dynamically allocated)

} node;

unsigned long searchNode(unsigned long id, node *nodes, unsigned long nnodes);
void ExitError(const char *miss, int errcode);

int main(int argc, char *argv[]) {
    clock_t start_time; // Start the timer
    FILE *mapfile;
    unsigned long nnodes;
    char *line = NULL;
    size_t len;

    start_time = clock(); // Default map file name

    char mapname[80];
    strcpy(mapname, "andorra.csv");

    if (argc > 1) strcpy(mapname, argv[1]); // Use command line argument if provided

    mapfile = fopen(mapname, "r");
    if (mapfile == NULL) {
        printf("Error when opening the file\n");
        return 1;
    }
    // count the nodes
    nnodes = 0UL;
    while (getline(&line, &len, mapfile) != -1) {
        if (strncmp(line, "node", 4) == 0) {
            nnodes++;
        }
    }
    printf("Total number of nodes is %ld\n", nnodes);
    rewind(mapfile);
    
    node *nodes;
    char *tmpline, *field, *ptr;
    unsigned long index = 0;

    nodes = (node *) malloc(nnodes * sizeof(node));
    if (nodes == NULL) {
        printf("Error when allocating the memory for the nodes\n");
        return 2;
    }

    while (getline(&line, &len, mapfile) != -1) {
        if (strncmp(line, "#", 1) == 0) continue;
        tmpline = line;  // make a copy of line to tmpline to keep the pointer of line
        field = strsep(&tmpline, "|");
        if (strcmp(field, "node") == 0) {
            field = strsep(&tmpline, "|");
            nodes[index].id = strtoul(field, &ptr, 10);
            field = strsep(&tmpline, "|");

            // Allocate memory for node name
            nodes[index].name = (char *)malloc(200 * sizeof(char));
            if (nodes[index].name == NULL) {
                printf("Error when allocating the memory for node name\n");
                return 2;
            }
            strcpy(nodes[index].name, field);

            for (int i = 0; i < 7; i++)
                field = strsep(&tmpline, "|");
            nodes[index].lat = atof(field);
            field = strsep(&tmpline, "|");
            nodes[index].lon = atof(field);

            nodes[index].nsucc = 0; // start with 0 successors
            nodes[index].successors = NULL; // initialize to NULL

            index++;
        }
    }
    printf("Assigned data to %ld nodes\n", index);
    printf("Last node has:\n id=%lu\n GPS=(%lf,%lf)\n Name=%s\n", nodes[index - 1].id, nodes[index - 1].lat, nodes[index - 1].lon, nodes[index - 1].name);

    rewind(mapfile);

    int oneway;
    unsigned long nedges = 0, origin, dest, originId, destId;
    while (getline(&line, &len, mapfile) != -1) {
        if (strncmp(line, "#", 1) == 0) continue;
        tmpline = line;  // make a copy of line to tmpline to keep the pointer of line
        field = strsep(&tmpline, "|");
        if (strcmp(field, "way") == 0) {
            for (int i = 0; i < 7; i++) field = strsep(&tmpline, "|");  // skip 7 fields
            if (strcmp(field, "") == 0) oneway = 0;  // no oneway
            else if (strcmp(field, "oneway") == 0) oneway = 1;
            else continue;  // No correct information
            field = strsep(&tmpline, "|");  // skip 1 field
            field = strsep(&tmpline, "|");
            if (field == NULL) continue;
            originId = strtoul(field, &ptr, 10);
            origin = searchNode(originId, nodes, nnodes);
            while (1) {
                field = strsep(&tmpline, "|");
                if (field == NULL) break;
                destId = strtoul(field, &ptr, 10);
                dest = searchNode(destId, nodes, nnodes);
                if ((origin == nnodes + 1) || (dest == nnodes + 1)) {
                    originId = destId;
                    origin = dest;
                    continue;
                }
                if (origin == dest) continue;
                // Check if the edge did appear in a previous way
                int newdest = 1;
                for (int i = 0; i < nodes[origin].nsucc; i++)
                    if (nodes[origin].successors[i] == dest) {
                        newdest = 0;
                        break;
                    }
                if (newdest) {
                    nodes[origin].nsucc++;
                    nodes[origin].successors = realloc(nodes[origin].successors, nodes[origin].nsucc * sizeof(unsigned long));
                    nodes[origin].successors[nodes[origin].nsucc - 1] = dest;
                    nedges++;
                }
                if (!oneway) {
                    // Check if the edge did appear in a previous way
                    int newor = 1;
                    for (int i = 0; i < nodes[dest].nsucc; i++)
                        if (nodes[dest].successors[i] == origin) {
                            newor = 0;
                            break;
                        }
                    if (newor) {
                        nodes[dest].nsucc++;
                        nodes[dest].successors = realloc(nodes[dest].successors, nodes[dest].nsucc * sizeof(unsigned long));
                        nodes[dest].successors[nodes[dest].nsucc - 1] = origin;
                        nedges++;
                    }
                }
                originId = destId;
                origin = dest;
            }
        }
    }

    fclose(mapfile);
    printf("Assigned %ld edges\n", nedges);
    printf("Elapsed time: %f seconds\n", (float)(clock() - start_time) / CLOCKS_PER_SEC);

    // Look for a node with more than 4 successors
    for (unsigned long i = 0; i < nnodes; i++) {
        if (nodes[i].nsucc > 4) {
            index = i;
            break;
        }
    }
    // Optionally print a node with more than 4 successors to verify correct reading
    printf("Node %lu has id=%lu and %u successors:\n", index, nodes[index].id, nodes[index].nsucc);
    for (int i = 0; i < nodes[index].nsucc; i++) printf("  Node %lu with id %lu.\n", nodes[index].successors[i], nodes[nodes[index].successors[i]].id);

    // Open the binary file for writing
    FILE *binmapfile;
    char binmapname[80];
    strcpy(binmapname, mapname);
    strcat(binmapname, ".bin");

    if ((binmapfile = fopen (binmapname, "wb")) == NULL)
        ExitError("the output binary data file cannot be opened", 31);

    if(fwrite(&nnodes, sizeof(unsigned long), 1, binmapfile)!=1){
        ExitError("when initializing the output binary data file", 32);
        }

    // Write nodes and their successors to the binary file
    for (unsigned long i = 0; i < nnodes; i++) {
         if (fwrite(&nodes[i], sizeof(node) - sizeof(unsigned long *) - sizeof(char *), 1, binmapfile) != 1) {
        ExitError("Error when writing node metadata to the binary file", 33);
    }

    // Write successors separately
    if (fwrite(nodes[i].successors, sizeof(unsigned long), nodes[i].nsucc, binmapfile) != nodes[i].nsucc) {
        ExitError("Error when writing node successors to the binary file", 34);
    }

    // Write the length of the name
    unsigned long name_length = strlen(nodes[i].name) + 1;
    if (fwrite(&name_length, sizeof(unsigned long), 1, binmapfile) != 1) {
        ExitError("Error when writing name length to the binary file", 35);
    }

    // Write the name itself
    if (fwrite(nodes[i].name, sizeof(char), name_length, binmapfile) != name_length) {
        ExitError("Error when writing name to the binary file", 36);
    }
    }

    fclose(binmapfile);

    // Free dynamically allocated memory
    for (unsigned long i = 0; i < nnodes; i++) {
        free(nodes[i].successors);
        free(nodes[i].name); 
    }
    free(nodes);

    return 0;
}

// Function to perform binary search for a node by its ID
unsigned long searchNode(unsigned long id, node *nodes, unsigned long nnodes) {
    // we know that the nodes where numerically ordered by id, so we can do a binary search.
    unsigned long l = 0, r = nnodes - 1, m;
    while (l <= r) {
        m = l + (r - l) / 2;
        if (nodes[m].id == id) return m;
        if (nodes[m].id < id)
            l = m + 1;
        else
            r = m - 1;
    }
    // id not found, we return nnodes + 1
    return nnodes + 1;
}

// Function to print error and exit
void ExitError(const char *miss, int errcode) {
    fprintf (stderr, "\nERROR: %s.\nStopping...\n\n", miss); exit(errcode);
}