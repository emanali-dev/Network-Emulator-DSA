#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <climits>
#include <sstream>


#define MAX_MACHINES 15
#define MAX_ROUTERS 5
#define MAX_NODES (MAX_MACHINES + MAX_ROUTERS)
#define MAX_PRIORITY 10
#define MAX_MSG_ID 1000
#define MAX_STR 256
#define DELAY 1 // 1 sec delay
using namespace std ;

// Message structure
struct Message {
    int id; // msg identity
    int priority; // ;;;
    char src[MAX_STR];
    char dest[MAX_STR];
    char payload[MAX_STR];
    char trace[MAX_STR];
    Message* next; // For linked list in FIFO queue
    Message() : id(0), priority(0), next(nullptr) {
        src[0] = dest[0] = payload[0] = trace[0] = '\0';
    }
};

// node for the ququue fifo
struct QueueNode {
    Message* msg;
    QueueNode* next;
    QueueNode(Message* m) : msg(m), next(nullptr) {}
};

// queue fir outgoing fifo queues
struct FIFOQueue {
    QueueNode* head;
    QueueNode* tail;
    FIFOQueue() : head(nullptr), tail(nullptr) {}
    ~FIFOQueue() {
        while (head) {
            // basic destructor code
            QueueNode* temp = head;
            head = head->next;
            // deletinbg msg
            delete temp->msg; 
            delete temp;
        }
    }

    // baisc enque fun ction
    void enqueue(Message* msg) {
        if (!msg) return;
        QueueNode* node = new QueueNode(msg);
        if (!head) {
            head = tail = node;
        } else {
            tail->next = node;
            tail = node;
        }
    }
    Message* dequeue() {
        if (!head) return nullptr;
        QueueNode* temp = head;
        Message* msg = temp->msg;
        head = head->next;
        if (!head) tail = nullptr;
        delete temp;
        return msg;
    }
    bool isEmpty() const { return head == nullptr; }
};

// for incoming priority queue
struct BinaryHeap {
    Message** heap; // array of message pointers 
    int size; //no of msgs in the heap
    int capacity; // capacity of msgs that can be stored

    // constructor 
    BinaryHeap(int cap) : size(0), capacity(cap) {
        heap = new Message*[cap]; //dybamically adding space for msg* array
        for (int i = 0; i < cap; i++) heap[i] = nullptr; //iniyakizing all to null
    }
    // destructor
    ~BinaryHeap() {
        // prevents memory leaks....every new is paired with delete
        for (int i = 0; i < size; i++) delete heap[i]; // clean up msgs object
        delete[] heap; // clean up array itsekf
    }
    //basic swapping code 
    // swap 2 msgs in heap array
    void swap(int i, int j) {
        Message* temp = heap[i];
        heap[i] = heap[j];
        heap[j] = temp;
    }
    // helper func
    // calculated array index ofd node's parents, left child, right child
    int parent(int i) const { return (i - 1) / 2; }
    int left(int i) const { return 2 * i + 1; }
    int right(int i) const { return 2 * i + 2; }
    void heapify(int i) {
        int largest = i;
        int l = left(i), r = right(i);
        if (l < size && heap[l] && heap[l]->priority > heap[largest]->priority) largest = l;
        if (r < size && heap[r] && heap[r]->priority > heap[largest]->priority) largest = r;
        if (largest != i) {
            swap(i, largest);
            heapify(largest); // recursive call to fic the subtree recursively
        }
    }
    // adds msgs at end
    // after that it keeps swapping it upward if it has priorty highrer than its parent
    void insert(Message* msg) {
        if (!msg || size == capacity) {
            if (size == capacity) cout << "Warning: BinaryHeap overflow\n";
            return;
        }
        heap[size] = msg;
        int i = size++;
        while (i > 0 && heap[parent(i)] && heap[parent(i)]->priority < heap[i]->priority) {
            swap(i, parent(i));
            i = parent(i); // movce up the tree
        }
    }
// reomves and returns msg with hig prioroty
    Message* extractMax() {
        if (size == 0) return nullptr;
        Message* max = heap[0]; // saves msx element
        heap[0] = heap[--size]; // moves last element to root
        heap[size] = nullptr; // clears last slot
        heapify(0); // restore heap order
        return max;
    }
    // to check if heap has any element
    bool isEmpty() const { return size == 0; }
};

// Splay Tree Node for routing table
struct SplayNode {
    char dest[MAX_STR];
    char nextHop[MAX_STR];
    SplayNode* left;
    SplayNode* right;
    SplayNode(const char* d, const char* nh) : left(nullptr), right(nullptr) {
        strncpy(dest, d, MAX_STR);
        strncpy(nextHop, nh, MAX_STR);
    }
};

// Splay Tree for routing table
struct SplayTree {
    SplayNode* root;
    SplayTree() : root(nullptr) {}
    ~SplayTree() { destroy(root); } // destructore
    void destroy(SplayNode* node) {
        if (node) {
            destroy(node->left); //Recursively deletes all nodes (post-order traversal).
            //Ensures no memory leaks.
            destroy(node->right);
            delete node;
        }
    }
    //rotations are needed to restructure the tree during splaying.


    SplayNode* rightRotate(SplayNode* x) {
        if (!x || !x->left) return x;
        SplayNode* y = x->left;
        x->left = y->right;
        y->right = x;
        return y;
    }
    SplayNode* leftRotate(SplayNode* x) {
        if (!x || !x->right) return x;
        SplayNode* y = x->right;
        x->right = y->left;
        y->left = x;
        return y;
    }
    SplayNode* splay(SplayNode* root, const char* key) {
        if (!root || !key || strcmp(root->dest, key) == 0) return root;
        if (strcmp(key, root->dest) < 0) {
            if (!root->left) return root;
            if (strcmp(key, root->left->dest) < 0) {
                root->left->left = splay(root->left->left, key);
                root = rightRotate(root);
            } else if (strcmp(key, root->left->dest) > 0) {
                root->left->right = splay(root->left->right, key);
                if (root->left->right) root->left = leftRotate(root->left);
            }
            return root->left ? rightRotate(root) : root;
        } else {
            if (!root->right) return root;
            if (strcmp(key, root->right->dest) > 0) {
                root->right->right = splay(root->right->right, key);
                root = leftRotate(root);
            } else if (strcmp(key, root->right->dest) < 0) {
                root->right->left = splay(root->right->left, key);
                if (root->right->left) root->right = rightRotate(root->right);
            }
            return root->right ? leftRotate(root) : root;
        }
    }
    /*
    Rotations are needed to restructure the tree during splaying.
    Steps:

    If tree is empty, insert at root.
    Else, splay the tree with the key.
    If already present, just update nextHop.
    Else:
    Create new node.
    Rearrange existing tree:
    If new node < root, it becomes root with right = root.
    Else, left = root.


    */
    void insert(const char* dest, const char* nextHop) {
        if (!dest || !nextHop) return;
        if (!root) {
            root = new SplayNode(dest, nextHop);
            return;
        }
        root = splay(root, dest);
        if (strcmp(dest, root->dest) == 0) {
            strncpy(root->nextHop, nextHop, MAX_STR);
            return;
        }
        SplayNode* newNode = new SplayNode(dest, nextHop);
        if (strcmp(dest, root->dest) < 0) {
            newNode->right = root;
            newNode->left = root->left;
            root->left = nullptr;
        } else {
            newNode->left = root;
            newNode->right = root->right;
            root->right = nullptr;
        }
        root = newNode;
    }
    /*
    Splay with dest to bring it to the root.
    If root doesn’t match → key not found.
    Else:
    If no left child: just shift right child to root.
    Else:
    Find max node in left subtree.
    Splay it to root of left subtree.
    Attach right subtree.
    Delete old root.
    */
    void remove(const char* dest) {
        if (!root || !dest) return;
        root = splay(root, dest);
        if (strcmp(dest, root->dest) != 0) return;
        SplayNode* temp = root;
        if (!root->left) {
            root = root->right;
        } else {
            SplayNode* maxLeft = root->left;
            while (maxLeft->right) maxLeft = maxLeft->right;
            root->left = splay(root->left, maxLeft->dest);
            root->left->right = root->right;
            root = root->left;
        }
        delete temp;
    }
    char* find(const char* dest) {
        if (!root || !dest) return nullptr;
        root = splay(root, dest);
        if (strcmp(dest, root->dest) == 0) return root->nextHop;
        return nullptr;
    }
};

// Linear List for routing table
// Linear List for routing table
struct LinearList {
    struct Entry {
        char dest[MAX_STR];
        char nextHop[MAX_STR];
        Entry* next;
        Entry(const char* d, const char* nh) : next(nullptr) {
            strncpy(dest, d, MAX_STR);
            strncpy(nextHop, nh, MAX_STR);
        }
    };
    Entry* head;
    LinearList() : head(nullptr) {}
    ~LinearList() {
        while (head) {
            Entry* temp = head;
            head = head->next;
            delete temp;
        }
    }
    void insert(const char* dest, const char* nextHop) {
        if (!dest || !nextHop) return;
        Entry* curr = head;
        while (curr && strcmp(curr->dest, dest) != 0) curr = curr->next;
        if (curr) {
            strncpy(curr->nextHop, nextHop, MAX_STR);
            return;
        }
        Entry* newEntry = new Entry(dest, nextHop);
        newEntry->next = head;
        head = newEntry;
    }
    void remove(const char* dest) {
        if (!dest) return;
        Entry* curr = head;
        Entry* prev = nullptr;
        while (curr && strcmp(curr->dest, dest) != 0) {
            prev = curr;
            curr = curr->next;
        }
        if (!curr) return;
        if (prev) prev->next = curr->next;
        else head = curr->next;
        delete curr;
    }
    char* find(const char* dest) {
        if (!dest) return nullptr;
        Entry* curr = head;
        while (curr && strcmp(curr->dest, dest) != 0) curr = curr->next;
        return curr ? curr->nextHop : nullptr;
    }
};
// Router structure
struct Router {
    char id[MAX_STR];
    BinaryHeap* incomingQueue;
    FIFOQueue** outgoingQueues;
    int numNeighbors;
    char** neighbors;
    bool useSplayTree;
    SplayTree* splayTable;
    LinearList* linearTable;
        //router class construcror
    Router(const char* id_, int maxMsgs, bool useSplay) : incomingQueue(nullptr), numNeighbors(0), useSplayTree(useSplay), splayTable(nullptr), linearTable(nullptr) {
       // Copy router ID
        strncpy(id, id_, MAX_STR);
        //binary heap for incoming msgs
        incomingQueue = new BinaryHeap(maxMsgs);
       //allocation of arrays for outgoin queues & neighbor names
        outgoingQueues = new FIFOQueue*[MAX_NODES];
        neighbors = new char*[MAX_NODES];
       // initilizeing arrays
        for (int i = 0; i < MAX_NODES; i++) {
            outgoingQueues[i] = nullptr;
            neighbors[i] = nullptr;
        }
        //creation of splaytree/linearlist
        splayTable = useSplay ? new SplayTree() : nullptr;
        linearTable = useSplay ? nullptr : new LinearList();
    }
    //destructor
    ~Router() {
        delete incomingQueue;   //deleting binary heap
       
        //loop that deletes each individual outgoin queue & neighbor name
        for (int i = 0; i < MAX_NODES; i++) {
            delete outgoingQueues[i];
            delete[] neighbors[i];
        }

        //deleting arrays
        delete[] outgoingQueues;
        delete[] neighbors;

        //deleting routing table
        delete splayTable;
        delete linearTable;
    }

    //fn to add meibhor to router
    void addNeighbor(const char* neighbor) {
       
        if (!neighbor || numNeighbors >= MAX_NODES) return;
       
       //allocate memory and copy neighnor id
        neighbors[numNeighbors] = new char[MAX_STR];
        strncpy(neighbors[numNeighbors], neighbor, MAX_STR);
       
       //create corresponding outgoimg queue
        outgoingQueues[numNeighbors] = new FIFOQueue();
        numNeighbors++;
    }

    //finding index of given neighbor by id
    int findNeighborIndex(const char* neighbor) const {
        if (!neighbor) return -1;
       
       //loop to find match through neighbor array
        for (int i = 0; i < numNeighbors; i++)
            if (neighbors[i] && strcmp(neighbors[i], neighbor) == 0) return i;
        return -1;
    }
};

// Machine structure
struct Machine {
    char id[MAX_STR];   //id for machine
    FIFOQueue* incomingQueue;   //queue for recieving msges
    FIFOQueue* outgoingQueue;   //queue for sending msges
    char connectedRouter[MAX_STR];  //id of router 
    
    //initilaize machine with given id
    Machine(const char* id_) : incomingQueue(nullptr), outgoingQueue(nullptr) {
       
       //copy machine id
        strncpy(id, id_, MAX_STR);
       
        //create new queue for incomimg and outgoing msges
        incomingQueue = new FIFOQueue();
        outgoingQueue = new FIFOQueue();
       
       //initilize connected router id as empty
        connectedRouter[0] = '\0';
    }

    //destructor that cleans up dynamically allocated queues
    ~Machine() {
        delete incomingQueue;
        delete outgoingQueue;
    }
};

// Network Emulator
struct Network {
    Machine* machines[MAX_MACHINES];   // array of pointrs to machine obj 
    Router* routers[MAX_ROUTERS];   //array of poniters to router obj
    int machineCount;   //cuureent no. of machines
    int routerCount;    //cuurnnt no. of routers
    int graph[MAX_NODES][MAX_NODES];    //matrix representing link weights
    time_t lastInterrupt;   //timestamp for last event

    //constructor to initialize counts ...pointrs...matrix
    Network() : machineCount(0), routerCount(0), lastInterrupt(0) {
        for (int i = 0; i < MAX_MACHINES; i++) machines[i] = nullptr;
        for (int i = 0; i < MAX_ROUTERS; i++) routers[i] = nullptr;
        for (int i = 0; i < MAX_NODES; i++)
            for (int j = 0; j < MAX_NODES; j++)
                graph[i][j] = INT_MAX;
    }

    //destructor to clean all dynamically allocayed machine and routers
    ~Network() {
        for (int i = 0; i < MAX_MACHINES; i++) delete machines[i];
        for (int i = 0; i < MAX_ROUTERS; i++) delete routers[i];
    }

    // fn to return index of node (machine or roter) in graph
    int getNodeIndex(const char* id) const {
        if (!id) return -1;
        for (int i = 0; i < machineCount; i++)
            if (machines[i] && strcmp(machines[i]->id, id) == 0) return i;
        for (int i = 0; i < routerCount; i++)
            if (routers[i] && strcmp(routers[i]->id, id) == 0) return MAX_MACHINES + i;
        return -1;
    }

    //fn to check if node id is valid and it exists in network
    bool isValidNodeId(const char* id) const {
        if (!id || (id[0] != 'M' && id[0] != 'R')) return false;
        int idx = getNodeIndex(id);
        return idx != -1;
    }

    //fn to load network topoly from csv file
   bool loadNetwork(const char* filename, bool useSplayTree) {
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Cannot open " << filename << "\n";
        return false;
    }

    char line[1024];
    char nodeIds[MAX_NODES][MAX_STR]; // store node ids in header
    int nodeCount = 0;

    // Read header to determine node order
    if (!file.getline(line, 1024)) {
        cout << "Error: Empty file or unable to read header\n";
        file.close();
        return false;
    }

    stringstream headerStream(line);
    string cell;

    // Skip the first column
    getline(headerStream, cell, ',');  //skip 1st label cell

    // Read node IDs from header
    while (getline(headerStream, cell, ',') && nodeCount < MAX_NODES) {
        strncpy(nodeIds[nodeCount], cell.c_str(), MAX_STR);
        nodeCount++;
    }

    if (nodeCount == 0) {
        cout << "Error: Invalid header in " << filename << "\n";
        file.close();
        return false;
    }

    // Read matrix rows
    int row = 0;
    while (file.getline(line, 1024) && row < MAX_NODES) {
        stringstream rowStream(line);
        string value;

        // 1st value in row is node id
        if (!getline(rowStream, value, ',')) {
            cout << "Warning: Empty line in " << filename << "\n";
            continue;
        }

        char nodeId[MAX_STR];
        strncpy(nodeId, value.c_str(), MAX_STR);

        // Find node index in header
        int nodeIdx = -1;
        for (int i = 0; i < nodeCount; ++i) {
            if (strcmp(nodeIds[i], nodeId) == 0) {
                nodeIdx = i;
                break;
            }
        }

        if (nodeIdx == -1) {
            cout << "Warning: Node " << nodeId << " not in header, skipping\n";
            continue;
        }

        // Create machine or router 
        if (nodeId[0] == 'M' && machineCount < MAX_MACHINES) {
            machines[machineCount] = new Machine(nodeId);
            cout << "Created machine: " << nodeId << "\n";
            machineCount++;
        } else if (nodeId[0] == 'R' && routerCount < MAX_ROUTERS) {
            routers[routerCount] = new Router(nodeId, MAX_MSG_ID, useSplayTree);
            cout << "Created router: " << nodeId << "\n";
            routerCount++;
        } else {
            cout << "Warning: Invalid node " << nodeId << " or limit reached, skipping\n";
            continue;
        }

        // Read the rest of the row (weights)
        int col = 0;
        while (getline(rowStream, value, ',') && col < nodeCount) {
            if (value != "?" && !value.empty()) {
                int weight = atoi(value.c_str());
                if (weight > 0) {
                    graph[nodeIdx][col] = weight;

                    //for routers and register neignoring nodes
                    if (nodeId[0] == 'R') {
                        char neighbor[MAX_STR];
                        strncpy(neighbor, nodeIds[col], MAX_STR);
                        if (routers[routerCount - 1]) {
                            routers[routerCount - 1]->addNeighbor(neighbor);
                            cout << "Added neighbor " << neighbor << " to " << nodeId << "\n";
                        }

                        //for machine , store connected router id
                    } else if (nodeIds[col][0] == 'R') {
                        if (machines[machineCount - 1]) {
                            strncpy(machines[machineCount - 1]->connectedRouter, nodeIds[col], MAX_STR);
                            cout << "Connected " << nodeId << " to router " << nodeIds[col] << "\n";
                        }
                    }
                }
            }
            col++;
        }

        row++;
    }

    file.close();


    //final chexk whter at least one node loaded
    if (machineCount == 0 && routerCount == 0) {
        cout << "Error: No valid nodes loaded from " << filename << "\n";
        return false;
    }

    cout << "Loaded " << machineCount << " machines and " << routerCount << " routers\n";
    
    //building routing tables for routers
    buildRoutingTables();
    return true;
}
   // Runs Dijkstra's algorithm to compute the shortest path from a router to all machines
void dijkstra(int src, char routingTable[][2 * MAX_STR], int& rtSize) {
    int dist[MAX_NODES];         // Distance from source to each node
    int prev[MAX_NODES];         // Predecessor of each node in the shortest path
    bool visited[MAX_NODES] = {false};  // Tracks visited nodes

    // Initialize distances and predecessors
    for (int i = 0; i < MAX_NODES; i++) {
        dist[i] = INT_MAX;
        prev[i] = -1;
    }
    dist[src] = 0;  // Distance from source to itself is 0

    // Main loop of Dijkstra's algorithm
    for (int count = 0; count < MAX_NODES - 1; count++) {
        int minDist = INT_MAX, u = -1;

        // Find the unvisited node with the smallest distance
        for (int v = 0; v < MAX_NODES; v++)
            if (!visited[v] && dist[v] <= minDist) {
                minDist = dist[v];
                u = v;
            }

        if (u == -1) break; // All remaining nodes are unreachable

        visited[u] = true;

        // Update distances of neighbors
        for (int v = 0; v < MAX_NODES; v++)
            if (!visited[v] && graph[u][v] != INT_MAX &&
                dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
                prev[v] = u;
            }
    }

    // Build routing table entries for reachable machines
    rtSize = 0;
    for (int i = 0; i < MAX_MACHINES; i++) {
        if (dist[i] == INT_MAX) continue; // Skip unreachable machines

        int curr = i;
        int nextHop = -1;

        // Trace back from machine to source to find the next hop
        while (curr != src && prev[curr] != -1) {
            nextHop = prev[curr];
            curr = nextHop;
        }

        // Only add to routing table if next hop is valid and not the source itself
        if (nextHop != -1 && nextHop != src) {
            char dest[MAX_STR], next[MAX_STR];
            snprintf(dest, MAX_STR, "M%d", i + 1);  // Destination ID

            // Determine next hop's ID
            if (nextHop < MAX_MACHINES && machines[nextHop]) {
                snprintf(next, MAX_STR, "%s", machines[nextHop]->id);
            } else if (nextHop >= MAX_MACHINES && nextHop < MAX_NODES && routers[nextHop - MAX_MACHINES]) {
                snprintf(next, MAX_STR, "%s", routers[nextHop - MAX_MACHINES]->id);
            } else {
                continue; // Skip invalid next hop
            }

            // Store routing table entry: destination + next hop
            strncpy(routingTable[rtSize], dest, MAX_STR);
            strncpy(routingTable[rtSize] + MAX_STR, next, MAX_STR);
            rtSize++;

            // Debug output
            cout << "Router R" << (src - MAX_MACHINES + 1) << ": Added route " << dest << " -> " << next << "\n";
        }
    }
}

// Builds routing tables for all routers using Dijkstra
void buildRoutingTables() {
    for (int i = 0; i < routerCount; i++) {
        if (!routers[i]) continue;

        char routingTable[MAX_MACHINES][2 * MAX_STR]; // [dest | nextHop]
        int rtSize = 0;

        // Run Dijkstra from this router
        dijkstra(MAX_MACHINES + i, routingTable, rtSize);

        // Insert each route into the router's table (Splay or Linear)
        for (int j = 0; j < rtSize; j++) {
            if (routers[i]->useSplayTree)
                routers[i]->splayTable->insert(routingTable[j], routingTable[j] + MAX_STR);
            else
                routers[i]->linearTable->insert(routingTable[j], routingTable[j] + MAX_STR);
        }
    }
}

// Sends a message into the network from a machine
void sendMessage(Message* msg) {
    if (!msg || !isValidNodeId(msg->src) || !isValidNodeId(msg->dest)) {
        cout << "Error: Invalid message source or destination (src: "
             << (msg ? msg->src : "null") << ", dest: " << (msg ? msg->dest : "null") << ")\n";
        delete msg;
        return;
    }

    int srcIdx = getNodeIndex(msg->src);

    // Only machines can initiate sending
    if (srcIdx < MAX_MACHINES) {
        strncpy(msg->trace, msg->src, MAX_STR); // Initialize trace
        machines[srcIdx]->outgoingQueue->enqueue(msg);
        cout << "Message " << msg->id << " enqueued from " << msg->src << " to " << msg->dest << "\n";
    } else {
        delete msg;
        cout << "Error: Source is not a machine\n";
    }
}

// Processes a router's message queue, forwarding the highest-priority message
void processRouter(Router* router) {
    if (!router) return;

    // Extract the highest-priority message
    Message* msg = router->incomingQueue->extractMax();
    if (!msg) return;

    // Find the next hop using routing table
    char nextHop[MAX_STR];
    char* nh = router->useSplayTree ?
        router->splayTable->find(msg->dest) :
        router->linearTable->find(msg->dest);

    if (!nh) {
        cout << "Error: No route to destination " << msg->dest << " from router " << router->id << "\n";
        delete msg;
        return;
    }

    strncpy(nextHop, nh, MAX_STR);

    // Append router ID to trace
    char newTrace[MAX_STR];
    snprintf(newTrace, MAX_STR, "%s:%s", msg->trace, router->id);
    strncpy(msg->trace, newTrace, MAX_STR);

    // Find index of the next hop neighbor
    int neighborIdx = router->findNeighborIndex(nextHop);

    if (neighborIdx != -1 && router->outgoingQueues[neighborIdx]) {
        router->outgoingQueues[neighborIdx]->enqueue(msg);
        cout << "Router " << router->id << ": Forwarded message " << msg->id << " to " << nextHop << "\n";
    } else {
        cout << "Error: Invalid neighbor " << nextHop << " from router " << router->id << "\n";
        delete msg;
    }
}

// Processes a machine's incoming queue to check for message delivery
void processMachine(Machine* machine) {
    if (!machine) return;

    Message* msg = machine->incomingQueue->dequeue();
    if (!msg) return;

    if (strcmp(msg->dest, machine->id) == 0) {
        // Final destination reached
        printPath(msg);
        cout << "Machine " << machine->id << ": Delivered message " << msg->id << "\n";
        delete msg;
    } else {
        // Message not for this machine; requeue it
        machine->outgoingQueue->enqueue(msg);
        cout << "Machine " << machine->id << ": Re-enqueued message " << msg->id << "\n";
    }
}
void timerInterrupt() {
    // Skip this call if not enough time has passed since last run
    if (time(nullptr) - lastInterrupt < DELAY) return;

    // Update the timestamp of this interrupt
    lastInterrupt = time(nullptr);

    // Go through all machines
    for (int i = 0; i < machineCount; i++) {
        // If machine slot is empty, skip it
        if (!machines[i]) continue;

        // Try to get a message from this machine's outgoing queue
        Message* msg = machines[i]->outgoingQueue->dequeue();

        // If there's a message to send
        if (msg) {
            // Find the index of the router this machine is connected to
            int routerIdx = getNodeIndex(machines[i]->connectedRouter) - MAX_MACHINES;

            // Check if router index is valid and router exists
            if (routerIdx >= 0 && routerIdx < routerCount && routers[routerIdx]) {
                // Push the message into the router’s incoming queue
                routers[routerIdx]->incomingQueue->insert(msg);
                // Show confirmation in console
                cout << "Moved message " << msg->id << " from " << machines[i]->id << " to " << routers[routerIdx]->id << "\n";
            } else {
                // If router is invalid, show error and delete message
                cout << "Error: Invalid connected router " << machines[i]->connectedRouter << " for machine " << machines[i]->id << "\n";
                delete msg;
            }
        }
    }

    // Now go through each router
    for (int i = 0; i < routerCount; i++) {
        // Skip if router doesn’t exist
        if (!routers[i]) continue;

        // For each neighbor of the current router
        for (int j = 0; j < routers[i]->numNeighbors; j++) {
            // Skip if this outgoing queue is null
            if (!routers[i]->outgoingQueues[j]) continue;

            // Try to get a message from this outgoing queue
            Message* msg = routers[i]->outgoingQueues[j]->dequeue();

            // If we have a message
            if (msg) {
                // Find the index of the neighbor node
                int nextIdx = getNodeIndex(routers[i]->neighbors[j]);

                // If it's a machine, send message to its incoming queue
                if (nextIdx < MAX_MACHINES && machines[nextIdx]) {
                    machines[nextIdx]->incomingQueue->enqueue(msg);
                    cout << "Moved message " << msg->id << " from " << routers[i]->id << " to " << machines[nextIdx]->id << "\n";

                // If it's another router, send message to its incoming queue
                } else if (nextIdx >= MAX_MACHINES && nextIdx < MAX_NODES && routers[nextIdx - MAX_MACHINES]) {
                    routers[nextIdx - MAX_MACHINES]->incomingQueue->insert(msg);
                    cout << "Moved message " << msg->id << " from " << routers[i]->id << " to " << routers[nextIdx - MAX_MACHINES]->id << "\n";

                // If neighbor is invalid, delete the message
                } else {
                    cout << "Error: Invalid next node " << routers[i]->neighbors[j] << "\n";
                    delete msg;
                }
            }
        }
    }
}
void printPath(const Message* msg) {
    // If the message is null, just return
    if (!msg) return;

    // Open path.txt to add the message trace
    ofstream file("path.txt", ios::app);

    // If file couldn’t be opened, show error and return
    if (!file.is_open()) {
        cout << "Error: Cannot open path.txt\n";
        return;
    }

    // Write the message ID, trace, and destination to the file
    file << msg->id << ":" << msg->trace << ":" << msg->dest << "\n";

    // Close the file
    file.close();
}
void changeRoutingTable(const char* routerId, const char* filename, bool add) {
    // Get router index from ID
    int routerIdx = getNodeIndex(routerId) - MAX_MACHINES;

    // If router doesn't exist or index is out of bounds
    if (routerIdx < 0 || routerIdx >= routerCount || !routers[routerIdx]) {
        cout << "Error: Invalid router " << routerId << "\n";
        return;
    }

    // Open the routing table file
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Cannot open " << filename << "\n";
        return;
    }

    // Read file line by line
    char line[512];
    while (file.getline(line, 512)) {
        // Split line into destination and next hop
        char* dest = strtok(line, ":");
        char* nextHop = strtok(nullptr, ":");

        // If both values exist and next hop is valid
        if (dest && nextHop && isValidNodeId(nextHop)) {
            // Add or remove entry based on the flag
            if (add) {
                if (routers[routerIdx]->useSplayTree)
                    routers[routerIdx]->splayTable->insert(dest, nextHop);
                else
                    routers[routerIdx]->linearTable->insert(dest, nextHop);
            } else {
                if (routers[routerIdx]->useSplayTree)
                    routers[routerIdx]->splayTable->remove(dest);
                else
                    routers[routerIdx]->linearTable->remove(dest);
            }
        } else {
            // If format is wrong, show warning
            cout << "Warning: Invalid routing table entry in " << filename << "\n";
        }
    }

    // Close the file
    file.close();
}
void changeEdge(const char* filename) {
    // Open CSV file to read graph edges
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Cannot open " << filename << "\n";
        return;
    }

    char line[1024];
    file.getline(line, 1024); // skip header

    int row = 0;

    // Read each row for graph weights
    while (file.getline(line, 1024) && row < MAX_NODES) {
        strtok(line, ","); // skip first column (node ID)
        int col = 0;

        // For each column, update weight or set to infinity
        while (char* token = strtok(nullptr, ",")) {
            if (strcmp(token, "?") != 0 && token[0] != '\0') {
                int weight = atoi(token);
                graph[row][col] = (weight > 0) ? weight : INT_MAX;
            } else {
                graph[row][col] = INT_MAX;
            }
            col++;
        }
        row++;
    }

    // Rebuild all routing tables based on updated graph
    file.close();
    buildRoutingTables();
}
void changeEdge(const char* node1, const char* node2, int weight) {
    // If input values are invalid
    if (!node1 || !node2 || weight <= 0) {
        cout << "Error: Invalid edge parameters\n";
        return;
    }

    // Get indices for both nodes
    int idx1 = getNodeIndex(node1);
    int idx2 = getNodeIndex(node2);

    // If both indices are valid
    if (idx1 != -1 && idx2 != -1) {
        // Update edge weights both ways
        graph[idx1][idx2] = graph[idx2][idx1] = weight;

        // Rebuild routing tables
        buildRoutingTables();
    } else {
        // Invalid nodes
        cout << "Error: Invalid nodes " << node1 << " or " << node2 << "\n";
    }
}
void sendMessageFromFile(const char* filename) {
    // Open the file with message data
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Error: Cannot open " << filename << "\n";
        return;
    }

    char line[512];

    // Read messages line by line
    while (file.getline(line, 512)) {
        Message* msg = new Message();

        // Parse message ID
        char* token = strtok(line, ":");
        if (!token) { delete msg; continue; }
        msg->id = atoi(token);

        // Parse priority
        token = strtok(nullptr, ":");
        if (!token) { delete msg; continue; }
        msg->priority = atoi(token);

        // Validate priority range
        if (msg->priority < 0 || msg->priority > MAX_PRIORITY) {
            cout << "Warning: Invalid priority in " << filename << "\n";
            delete msg;
            continue;
        }

        // Parse source
        token = strtok(nullptr, ":");
        if (!token) { delete msg; continue; }
        strncpy(msg->src, token, MAX_STR);

        // Parse destination
        token = strtok(nullptr, ":");
        if (!token) { delete msg; continue; }
        strncpy(msg->dest, token, MAX_STR);

        // Parse payload
        token = strtok(nullptr, ":");
        if (!token) { delete msg; continue; }
        strncpy(msg->payload, token, MAX_STR);

        // Send the parsed message
        sendMessage(msg);
    }

    // Done reading
    file.close();
}
void printPathQuery(const char* src, const char* dest) {
    // Check if source or destination is null
    if (!src || !dest) {
        cout << "Error: Invalid source or destination\n";
        return;
    }

    // Open the path log file
    ifstream file("path.txt");

    // Show error if file couldn't open
    if (!file.is_open()) {
        cout << "Error: Cannot open path.txt\n";
        return;
    }

    // Temporary buffer to hold each line
    char line[512];

    // Flag to track if anything was found
    bool found = false;

    // Go through each line in path.txt
    while (file.getline(line, 512)) {
        // Extract message ID
        char* id = strtok(line, ":");

        // Extract source of the message
        char* pathSrc = strtok(nullptr, ":");

        // If ID or source is missing, skip this line
        if (!id || !pathSrc) continue;

        // Reconstruct the full path trace string
        char path[MAX_STR] = "";
        char* token = pathSrc;
        char* lastToken = pathSrc;

        // Append all tokens into the path string
        while (token) {
            strcat(path, token);         // Add token to path
            lastToken = token;           // Keep track of last token in case destination is missing
            token = strtok(nullptr, ":"); // Move to next part
            if (token) strcat(path, ":"); // Add separator
        }

        // Get the destination part (if any)
        char* pathDest = strtok(nullptr, ":");

        // If destination was missing, use the last seen token
        if (!pathDest) pathDest = lastToken;

        // Check if the path matches the query (wildcards allowed)
        if ((strcmp(src, "*") == 0 || strcmp(src, pathSrc) == 0) &&
            (strcmp(dest, "*") == 0 || strcmp(dest, pathDest) == 0)) {
            cout << id << ":" << path << ":" << pathDest << "\n";
            found = true;
        }
    }

    // If nothing matched the query
    if (!found) cout << "No matching paths found\n";

    // Close the file
    file.close();
}
};
int main() {
    Network network; // Create a network instance

    bool useSplayTree = false; // Flag to choose splay tree or linear list

    int input;
    cout << "Use Splay Tree for routing table? (1 for yes, 0 for no): ";
    cin >> input; // Ask user to pick the routing structure

    // Set the structure based on user input
    if (input == 1) useSplayTree = true;
    else if (input != 0) {
        cout << "Invalid input, using Linear List\n"; // Warn for invalid input
    }

    cin.ignore(); // Clear newline from input buffer

    // Load network topology from CSV
    if (!network.loadNetwork("network.csv", useSplayTree)) {
        cout << "Failed to load network. Exiting.\n";
        return 1; // Exit if loading fails
    }

    // Start main simulation loop
    while (true) {
        network.timerInterrupt(); // Move messages between nodes

        // Let each router process its messages
        for (int i = 0; i < network.routerCount; i++)
            network.processRouter(network.routers[i]);

        // Let each machine process its messages
        for (int i = 0; i < network.machineCount; i++)
            network.processMachine(network.machines[i]);

        // Display the interactive menu
        cout << "\nMenu:\n"
             << "1. Send message\n"
             << "2. Send message from file\n"
             << "3. Change routing table\n"
             << "4. Change edge\n"
             << "5. Print path\n"
             << "6. Exit\n"
             << "Choose: ";

        int choice;
        if (!(cin >> choice)) {
            // Handle bad input
            cout << "Error: Invalid input\n";
            cin.clear();
            cin.ignore(256, '\n');
            continue;
        }

        cin.ignore(); // Clear buffer

        // Option 1: Send single message manually
        if (choice == 1) {
            Message* msg = new Message(); // Create message object

            cout << "ID: ";
            if (!(cin >> msg->id)) {
                cout << "Error: Invalid ID\n";
                delete msg;
                cin.clear();
                cin.ignore(256, '\n');
                continue;
            }

            cout << "Priority (0-10): ";
            if (!(cin >> msg->priority) || msg->priority < 0 || msg->priority > MAX_PRIORITY) {
                cout << "Error: Invalid priority\n";
                delete msg;
                cin.clear();
                cin.ignore(256, '\n');
                continue;
            }

            cin.ignore(); // Clear newline

            // Get remaining fields
            cout << "Source: ";
            cin.getline(msg->src, MAX_STR);
            cout << "Destination: ";
            cin.getline(msg->dest, MAX_STR);
            cout << "Payload: ";
            cin.getline(msg->payload, MAX_STR);

            // Send the message into the network
            network.sendMessage(msg);
        }

        // Option 2: Send messages from a file
        else if (choice == 2) {
            char filename[256];
            cout << "Filename: ";
            cin.getline(filename, 256);
            network.sendMessageFromFile(filename);
        }

        // Option 3: Modify routing table
        else if (choice == 3) {
            char routerId[256], filename[256], op[10];
            cout << "Router ID: ";
            cin.getline(routerId, 256);
            cout << "Add or remove (add/remove): ";
            cin.getline(op, 10);
            cout << "Filename: ";
            cin.getline(filename, 256);

            // Validate operation type
            if (strcmp(op, "add") != 0 && strcmp(op, "remove") != 0) {
                cout << "Error: Invalid operation\n";
                continue;
            }

            // Apply routing table update
            network.changeRoutingTable(routerId, filename, strcmp(op, "add") == 0);
        }

        // Option 4: Modify a single edge weight
        else if (choice == 4) {
            char node1[256], node2[256];
            int weight;

            cout << "Node 1: ";
            cin.getline(node1, 256);
            cout << "Node 2: ";
            cin.getline(node2, 256);
            cout << "Weight: ";
            if (!(cin >> weight)) {
                cout << "Error: Invalid weight\n";
                cin.clear();
                cin.ignore(256, '\n');
                continue;
            }

            cin.ignore(); // Clear buffer
            network.changeEdge(node1, node2, weight); // Update the graph
        }

        // Option 5: Query paths based on source and destination
        else if (choice == 5) {
            char src[256], dest[256];
            cout << "Source (* for all): ";
            cin.getline(src, 256);
            cout << "Destination (* for all): ";
            cin.getline(dest, 256);
            network.printPathQuery(src, dest);
        }

        // Option 6: Exit program
        else if (choice == 6) {
            break;
        }

        // Handle invalid menu choices
        else {
            cout << "Error: Invalid choice\n";
        }
    }

    return 0; // End of program
}