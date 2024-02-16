//Main Code

#include <iostream>
#include <fstream>
#include <limits.h>
#include <queue>
#include <climits>
#include <Windows.h>
#include <string>
#include <thread>
#include <chrono>
#include <mmsystem.h>
using namespace std;
#pragma comment(lib,"winmm.lib")
// Constants for the size of the map

const int ROWS = 20;
const int COLS = 20;
const int INF = INT_MAX;


// Constants for the weights of the edges
const int DEFAULT_WEIGHT = 1;
const int OBSTACLE_WEIGHT = INT_MAX;
const int ENEMY_WEIGHT = INT_MAX;

// Struct to represent a cell in the map
struct Gmae {
    char data;
    int row;
    int col;
    bool visited;
    int weight;
    Gmae* up = nullptr;
    Gmae* down = nullptr;
    Gmae* left = nullptr;
    Gmae* right = nullptr;
    Gmae* parent;
    int rank;
    int x_axes;
    int y_axes;
    bool enemy;
    bool obstacle;
    int distance;
    int neighbors[4][2];
};

struct Edge {
    Gmae* starting_location;
    Gmae* destination;
    int edge_weight;
};


// Struct to represent the player's position
struct Player {
    int row = 0;
    int col = 0;
};


// Function to move the player up
void moveUp(Player& player, Gmae* map[ROWS][COLS]) {
    if (map[player.row][player.col]->up && map[player.row][player.col]->up->data != '#') {
        player.row--;
    }
}

// Function to move the player down
void moveDown(Player& player, Gmae* map[ROWS][COLS]) {
    if (map[player.row][player.col]->down && map[player.row][player.col]->down->data != '#') {
        player.row++;
    }
}

// Function to move the player left
void moveLeft(Player& player, Gmae* map[ROWS][COLS]) {
    if (map[player.row][player.col]->left && map[player.row][player.col]->left->data != '#') {
        player.col--;
    }
}

// Function to move the player right
void moveRight(Player& player, Gmae* map[ROWS][COLS]) {
    if (map[player.row][player.col]->right && map[player.row][player.col]->right->data != '#') {
        player.col++;
    }
}

// Function to interact with the current cell
void interact(Player& player, Gmae* map[ROWS][COLS])
{
    char cellType = map[player.row][player.col]->data;
    switch (cellType) {
    case 'J':
        cout << "You found a jewel!" << endl;
        break;
    case 'P':
        cout << "You found a potion!" << endl;
        break;
    case 'W':
        cout << "You found a weapon!" << endl;
        break;
    case '%':
        cout << "You died!" << endl;
        exit(0);
        break;
    case '&':
        cout << "You encountered a dragon!" << endl;
        break;
    case '$':
        cout << "You encountered a goblin!" << endl;
        break;
    case '@':
        // The player is already here
        break;
    case '*':
        // The player reached the crystal
        break;
    default:
        cout << "Nothing here" << endl;
        break;
    }
}



//kuruskal algo
//// Function to read the map from a file and create a 2D linked list of cells
void map_reading_for_kuruskal_algo(Gmae* map[ROWS][COLS]) {
    ifstream inputFile;
    inputFile.open("map.txt");
    if (!inputFile) {
        cerr << "Error: could not open file" << endl;
        exit(1);
    }
    char ch;
    int row = 0;
    int col = 0;
    while (inputFile.get(ch)) {
        if (ch == '\n') {
            row++;
            col = 0;
            continue;
        }
        map[row][col] = new Gmae;
        map[row][col]->data = ch;
        map[row][col]->row = row;
        map[row][col]->col = col;
        col++;
    }
    inputFile.close();

    // Set the pointers for the linked list
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            if (i > 0) {
                map[i][j]->up = map[i - 1][j];
            }
            if (i < ROWS - 1) {
                map[i][j]->down = map[i + 1][j];
            }
            if (j > 0) {
                map[i][j]->left = map[i][j - 1];
            }
            if (j < COLS - 1) {
                map[i][j]->right = map[i][j + 1];
            }
            map[i][j]->parent = map[i][j];
            map[i][j]->rank = 0;
        }
    }
}



// Function to find the parent of a cell
Gmae* finding_in_for_Kuruskal(Gmae* cell) {
    if (cell != cell->parent) {
        cell->parent = finding_in_for_Kuruskal(cell->parent);//recurssion
    }
    return cell->parent;
}
// Function to merge two sets of cells
void merge_for_KrurskalAlgo(Gmae* cell1, Gmae* cell2) {
    cell1 = finding_in_for_Kuruskal(cell1);
    cell2 = finding_in_for_Kuruskal(cell2);
    if (cell1 == cell2) {
        return;
    }
    if (cell1->rank < cell2->rank) {
        cell1->parent = cell2;
    }
    else if (cell1->rank > cell2->rank) {
        cell2->parent = cell1;
    }
    else {
        cell2->parent = cell1;
        cell1->rank++;
    }
}

// Function to get the edges of the forest
void getEdges_for_kuruskal_algo(Gmae* map[ROWS][COLS], Edge edges[]) {
    int numEdges = 0;
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            Gmae* cell = map[i][j];
            if (cell->down != nullptr) {
                edges[numEdges].starting_location = cell;
                edges[numEdges].destination = cell->down;
                edges[numEdges].edge_weight = (cell->data == '%' || cell->down->data == '%') ? INT_MAX : 1;
                numEdges++;
            }
            if (cell->right != nullptr) {
                edges[numEdges].starting_location = cell;
                edges[numEdges].destination = cell->right;
                edges[numEdges].edge_weight = (cell->data == '%' || cell->right->data == '%') ? INT_MAX : 1;
                numEdges++;
            }
        }
    }
}

// Function to sort the edges by weight
void sortEdges_for_kuruskal_algo(Edge edges[], int numEdges) {
    for (int i = 0; i < numEdges; i++) {
        for (int j = i + 1; j < numEdges; j++) {
            if (edges[i].edge_weight > edges[j].edge_weight) {
                Edge temp = edges[i];
                edges[i] = edges[j];
                edges[j] = temp;
            }
        }
    }
}

// Function to find the minimum spanning tree of the forest
void kruskal_Algorithm(Gmae* map[ROWS][COLS]) {
    Edge edges[(ROWS - 1) * COLS + (COLS - 1) * ROWS];
    getEdges_for_kuruskal_algo(map, edges);
    sortEdges_for_kuruskal_algo(edges, (ROWS - 1) * COLS + (COLS - 1) * ROWS);

    int numEdges = 0;
    Edge mstEdges[ROWS * COLS - 1];
    for (int i = 0; i < (ROWS - 1) * COLS + (COLS - 1) * ROWS; i++) {
        Gmae* cell1 = edges[i].starting_location;
        Gmae* cell2 = edges[i].destination;
        if (finding_in_for_Kuruskal(cell1) != finding_in_for_Kuruskal(cell2)) {
            merge_for_KrurskalAlgo(cell1, cell2);
            mstEdges[numEdges] = edges[i];
            numEdges++;
        }
        if (numEdges == ROWS * COLS - 1) {
            break;
        }
    }

    // Print minimum spanning tree edges
    for (int i = 0; i < ROWS * COLS - 1; i++) {
        cout << "(" << mstEdges[i].starting_location->row << "," << mstEdges[i].starting_location->col << ") -- "
            << "(" << mstEdges[i].destination->row << "," << mstEdges[i].destination->col << ") "
            /*<< mstEdges[i].weight */ << endl;
    }
}




///////////prims algo
void Map_reading_for_primsAlgo(char map[][COLS], const char* fileName) {
    ifstream inputFile;
    inputFile.open(fileName);
    if (!inputFile) {
        cerr << "Error: could not open file" << endl;
        exit(1);
    }
    char ch;
    int row = 0;
    int col = 0;
    while (inputFile.get(ch)) {
        if (ch == '\n') {
            row++;
            col = 0;
            continue;
        }
        map[row][col] = ch;
        col++;
    }
    inputFile.close();
}


void prim_Algo(Gmae graph[][COLS]) {
    // Initialize all cells to unvisited and infinite weight
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            graph[i][j].visited = false;
            graph[i][j].weight = INT_MAX;
        }
    }

    // Start at (0, 0) and set its weight to 0
    int row = 0;
    int col = 0;
    graph[row][col].weight = 0;

    // Main loop
    for (int i = 0; i < ROWS * COLS; i++) {
        // Find the unvisited cell with the smallest weight
        int minWeight = INT_MAX;
        int minRow = -1;
        int minCol = -1;
        for (int j = 0; j < ROWS; j++) {
            for (int k = 0; k < COLS; k++) {
                if (!graph[j][k].visited && graph[j][k].weight < minWeight) {
                    minWeight = graph[j][k].weight;
                    minRow = j;
                    minCol = k;
                }
            }
        }

        // Mark the cell as visited
        graph[minRow][minCol].visited = true;

        // Update the weights of the neighboring cells
        if (minRow > 0 && !graph[minRow - 1][minCol].visited) {
            graph[minRow - 1][minCol].weight = min(graph[minRow - 1][minCol].weight, graph[minRow][minCol].weight + 1);
        }
        if (minRow < ROWS - 1 && !graph[minRow + 1][minCol].visited) {
            graph[minRow + 1][minCol].weight = min(graph[minRow + 1][minCol].weight, graph[minRow][minCol].weight + 1);
        }
        if (minCol > 0 && !graph[minRow][minCol - 1].visited) {
            graph[minRow][minCol - 1].weight = min(graph[minRow][minCol - 1].weight, graph[minRow][minCol].weight + 1);
        }
        if (minCol < COLS - 1 && !graph[minRow][minCol + 1].visited) {
            graph[minRow][minCol + 1].weight = min(graph[minRow][minCol + 1].weight, graph[minRow][minCol].weight + 1);
        }
    }
}

void print_min_spanning_tree_for_prim_algo(Gmae graph[][COLS]) {
    // Print the edges of the MST and compute total weight
    int totalWeight = 0;
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            if (graph[i][j].visited) {
                if (i > 0 && graph[i - 1][j].visited) {
                    cout << "(" << i << "," << j << ") -- (" << i - 1 << "," << j << ") " << endl;
                    totalWeight += graph[i][j].weight;
                }
                if (i < ROWS - 1 && graph[i + 1][j].visited) {
                    cout << "(" << i << "," << j << ") -- (" << i + 1 << "," << j << ") " << endl;
                    totalWeight += graph[i][j].weight;
                }
                if (j > 0 && graph[i][j - 1].visited) {
                    cout << "(" << i << "," << j << ") -- (" << i << "," << j - 1 << ") " << endl;
                    totalWeight += graph[i][j].weight;
                }
                if (j < COLS - 1 && graph[i][j + 1].visited) {
                    cout << "(" << i << "," << j << ") -- (" << i << "," << j + 1 << ") " << endl;
                    totalWeight += graph[i][j].weight;
                }
            }
        }
    }

}


void printAdjacencyMatrix(Gmae* map[ROWS][COLS]) {

    int adjMatrix[ROWS * COLS][ROWS * COLS] = { 0 };

    // Set weights of edges towards obstacles to 100, and all other edges to 1
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            Gmae* current = map[i][j];
            int currentIdx = i * COLS + j;
            if (current->up) {
                int upIdx = (i - 1) * COLS + j;
                if (current->up->data == '#' || current->up->data == '&' || current->up->data == '$' || current->up->data == '@') {
                    adjMatrix[currentIdx][upIdx] = 100;
                }
                else {
                    adjMatrix[currentIdx][upIdx] = 1;
                }
            }
            if (current->down) {
                int downIdx = (i + 1) * COLS + j;
                if (current->down->data == '#' || current->down->data == '&' || current->down->data == '$' || current->down->data == '@') {
                    adjMatrix[currentIdx][downIdx] = 100;
                }
                else {
                    adjMatrix[currentIdx][downIdx] = 1;
                }
            }
            if (current->left) {
                int leftIdx = i * COLS + (j - 1);
                if (current->left->data == '#' || current->left->data == '&' || current->left->data == '$' || current->left->data == '@') {
                    adjMatrix[currentIdx][leftIdx] = 100;
                }
                else {
                    adjMatrix[currentIdx][leftIdx] = 1;
                }
            }
            if (current->right) {
                int rightIdx = i * COLS + (j + 1);
                if (current->right->data == '#' || current->right->data == '&' || current->right->data == '$' || current->right->data == '@') {
                    adjMatrix[currentIdx][rightIdx] = 100;
                }
                else {
                    adjMatrix[currentIdx][rightIdx] = 1;
                }
            }
        }
    }

    //// Print the adjacency matrix
    //cout << "Adjacency Matrix:" << endl;
    //for (int i = 0; i < ROWS * COLS; i++) {
    //    for (int j = 0; j < ROWS * COLS; j++) {
    //        if (adjMatrix[i][j] != 0)
    //        {
    //            cout << adjMatrix[i][j] << " ";
    //        }
    //    }
    //    cout << endl;
    //}
}
Gmae map[ROWS][COLS];

void map_initialization_for_dijkstra() {
    // Initialize map with cells, enemies, obstacles, and distances
    // ...

    // Add neighbors to each cell
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            int idx = 0;
            if (i > 0 && !map[i - 1][j].obstacle) {
                map[i][j].neighbors[idx][0] = i - 1;
                map[i][j].neighbors[idx][1] = j;
                idx++;
            }
            if (i < ROWS - 1 && !map[i + 1][j].obstacle) {
                map[i][j].neighbors[idx][0] = i + 1;
                map[i][j].neighbors[idx][1] = j;
                idx++;
            }
            if (j > 0 && !map[i][j - 1].obstacle) {
                map[i][j].neighbors[idx][0] = i;
                map[i][j].neighbors[idx][1] = j - 1;
                idx++;
            }
            if (j < COLS - 1 && !map[i][j + 1].obstacle) {
                map[i][j].neighbors[idx][0] = i;
                map[i][j].neighbors[idx][1] = j + 1;
                idx++;
            }
            // Mark end of neighbor array with -1
            map[i][j].neighbors[idx][0] = -1;
            map[i][j].neighbors[idx][1] = -1;
        }
    }
}

void Dijkstra_algorithm(int startX, int startY, int endX, int endY) {
    int jewl = 0, weapon = 0, potion = 0, totalScore = 0;
    //SocreAvlTree s;
    // Initialize distances and visited flags
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            map[i][j].distance = INF;
            map[i][j].visited = false;
        }
    }

    // Set distance of starting cell to 0 and add it to priority queue
    map[startX][startY].distance = 0;
    priority_queue<pair<int, pair<int, int>>> pq;
    pq.push(make_pair(0, make_pair(startX, startY)));

    // Run Dijkstra's algorithm
    while (!pq.empty()) {
        // Get cell with minimum distance from priority queue
        int currX = pq.top().second.first;
        int currY = pq.top().second.second;
        pq.pop();



        // Check if we've reached the destination cell
        if (currX == endX && currY == endY) {
            break;
        }

        // Check if cell has already been visited
        if (map[currX][currY].visited) {
            continue;
        }
        map[currX][currY].visited = true;

        // Update distances of neighboring cells
        int idx = 0;
        while (map[currX][currY].neighbors[idx][0] != -1) {
            int neighborX = map[currX][currY].neighbors[idx][0];
            int neighborY = map[currX][currY].neighbors[idx][1];
            int dist = map[currX][currY].distance + (map[neighborX][neighborY].enemy ? INF : (map[neighborX][neighborY].obstacle ? 3 : 1));
            if (dist < map[neighborX][neighborY].distance) {
                map[neighborX][neighborY].distance = dist;
                pq.push(make_pair(-dist, make_pair(neighborX, neighborY)));
            }
            idx++;
        }
    }

    // Print shortest distance to crystal
    cout << "Shortest distance to crystal: (" << startX << ", " << startY << ") to (" << endX << ", " << endY << ") crystal: " << map[endX][endY].distance << endl;
}

// Function to calculate the shortest path between any two areas using Floyd's algorithm
void Floyd_algo(int adjMatrix[ROWS * COLS][ROWS * COLS]) {
    for (int k = 0; k < ROWS * COLS; k++) {
        for (int i = 0; i < ROWS * COLS; i++) {
            for (int j = 0; j < ROWS * COLS; j++) {
                if (adjMatrix[i][k] != INT_MAX &&
                    adjMatrix[k][j] != INT_MAX &&
                    adjMatrix[i][j] > adjMatrix[i][k] + adjMatrix[k][j]) {
                    adjMatrix[i][j] = adjMatrix[i][k] + adjMatrix[k][j];
                }
            }
        }
    }

}


class Node
{
public:
    int id; //generated randomly
    int reward_score; //score of a particular reward you obtained
    int count; //to avoid duplicate id nodes (maintain count of the id)
    Node* left;
    Node* right;
    int height;
};

class SocreAvlTree
{
private:
    Node* rootNode;

    Node* Node_creation(int id, int reward_score)
    {
        Node* newNode = new Node();
        newNode->id = id;
        newNode->height = 1;
        newNode->count = 1;
        newNode->left = newNode->right = nullptr;
        newNode->reward_score = reward_score;
        return newNode;
    }

    int Tree_Height(Node* node)
    {
        if (node == nullptr)
        {
            return 0;
        }
        return node->height;
    }

    int Tree_balancing(Node* node)
    {
        if (node == nullptr)
        {
            return 0;
        }
        return Tree_Height(node->left) - Tree_Height(node->right);
    }

    Node* Rotating_tree_right(Node* y)
    {
        Node* x = y->left;
        Node* t2 = x->right;

        x->right = y;
        y->left = t2;
        y->height = max(Tree_Height(y->left), Tree_Height(y->right)) + 1;
        x->height = max(Tree_Height(x->left), Tree_Height(x->right)) + 1;
        return x;
    }

    Node* Rotating_tree_left(Node* x)
    {
        Node* y = x->right;
        Node* t2 = y->left;

        y->left = x;
        x->right = t2;
        x->height = max(Tree_Height(x->left), Tree_Height(x->right)) + 1;
        y->height = max(Tree_Height(y->left), Tree_Height(y->right)) + 1;

        return y;
    }

    Node* Inserting_Node_in_Tree(Node* node, int id, int reward_score)
    {
        if (node == nullptr)
            return Node_creation(id, reward_score);

        if (id < node->id)
            node->left = Inserting_Node_in_Tree(node->left, id, reward_score);
        else if (id > node->id)
            node->right = Inserting_Node_in_Tree(node->right, id, reward_score);
        else
            node->count++; // Already exists, increment count

        node->height = max(Tree_Height(node->left), Tree_Height(node->right)) + 1;

        int balance = Tree_balancing(node);

        if (balance > 1 && id < node->left->id)
            return Rotating_tree_right(node);

        if (balance < -1 && id > node->right->id)
            return Rotating_tree_left(node);

        if (balance > 1 && id > node->left->id)
        {
            node->left = Rotating_tree_left(node->left);
            return Rotating_tree_right(node);
        }

        if (balance < -1 && id < node->right->id)
        {
            node->right = Rotating_tree_right(node->right);
            return Rotating_tree_left(node);
        }
        return node;
    }

    Node* Finding_minimum_val_Node(Node* node)
    {
        Node* current = node;

        while (current->left != nullptr)
        {
            current = current->left;
        }
        return current;
    }

public:
    SocreAvlTree()
    {
        rootNode = nullptr;
        srand(time(nullptr)); // Seed the random number generator
    }

    void insertNode(int reward_score)
    {
        int id;

        if (rootNode == nullptr)
            id = 100; // First node, use id 100
        else
        {
            // Generate random id between 0 and 200, excluding existing ids
            do
            {
                id = rand() % 201;
            } while (searchNode(id) != nullptr);
        }

        rootNode = Inserting_Node_in_Tree(rootNode, id, reward_score);
    }

    Node* searchNode(int id)
    {
        Node* current = rootNode;

        while (current != nullptr)
        {
            if (id == current->id)
            {
                return current;
            }

            if (id < current->id)
            {
                current = current->left;
            }
            else
            {
                current = current->right;
            }
        }
        return nullptr;
    }

    void printInorder()
    {
        Printing_desired_node_in_tree(rootNode);
    }

    void Printing_desired_node_in_tree(Node* node)
    {
        if (node != nullptr)
        {
            Printing_desired_node_in_tree(node->left);
            cout << "ID: " << node->id << ", Score: " << node->reward_score << ", Count: " << node->count << endl;
            Printing_desired_node_in_tree(node->right);
        }
    }
};
void Play_Sound() {
    PlaySound(TEXT("song.wav"), NULL, SND_FILENAME | SND_ASYNC | SND_LOOP);
}
//-------------------------------------------------  MAIN FUNCTION  --------------------------------------------------------------

int main() {

    thread soundThread(Play_Sound);
    std::cout << "   ___                  _      __              ___                 _          _   _  __ _                 _             \n"
        "  / _ \\  _  _  ___  ___| |_   / _| ___  _ _   / __| _ _  _  _  ___| |_  __ _ | | | |/ /(_) _ _   __ _  __| | ___  _ __  \n"
        " | (_) || || |/ -_)(_-<|  _|  | _|/ _ \\| '_| | (__ | '_|| || |(_-<|  _|/ _` || | | ' < | || ' \\ / _` |/ _` |/ _ \\| '  \\ \n"
        "  \\__\\_\\ \\_,_|\___|/__/ \\__| |_|  \\___/|_|    \\___||_|   \\_, |/__/ \\__|\\__,_||| |_|\_\\|_||_||_|\__, |\\__,_|\___/|_|_|_|_|\n"
        "                                                         |__/                                   |___/                   \n";





   /* string name;
    cout << "\t\t\tEnter Your Name: ";
    cin >> name;*/

    cout << " \n\n\n\t\t\t\t\tWELCOME \n";
    cout << "\t\t\tPress any key to start your adventure.\n";

    cin.get();

    ////adding sound
    //cout << "playing sound" << endl; 
    //PlaySound(TEXT("song.wav"), NULL, SND_FILENAME | SND_SYNC );
    //
   /* string input;
    getline(cin, input);
    PlaySound(0, 0, 0);
    cout << "Music stopped: " << endl;*/


    Gmae grid[ROWS][COLS];
    // Create a 2D array of struct cells of size 20 x 20
    Gmae* map[ROWS][COLS];
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            map[i][j] = new Gmae;
        }
    }

    //printAdjacencyMatrix(map);

    // Read the map from the file
    ifstream inputFile;
    inputFile.open("map.txt");
    if (!inputFile) {
        cerr << "Error: could not open file" << endl;
        exit(1);
    }
    char ch;
    int row = 0;
    int col = 0;
    while (inputFile.get(ch)) {
        if (ch == '\n') {
            row++;
            col = 0;
            continue;
        }
        map[row][col]->data = ch;
        map[row][col]->row = row;
        map[row][col]->col = col;
        col++;
    }
    inputFile.close();

    // Set the pointers for the linked list
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            if (i > 0) {
                map[i][j]->up = map[i - 1][j];
            }
            if (i < ROWS - 1) {
                map[i][j]->down = map[i + 1][j];
            }
            if (j > 0) {
                map[i][j]->left = map[i][j - 1];
            }
            if (j < COLS - 1) {
                map[i][j]->right = map[i][j + 1];
            }
        }
    }

    //// Read the input grid from a file
    ifstream infile("map.txt");

    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            char c;
            infile >> c;
            grid[i][j] = { c, i, j };
        }
    }
    infile.close();

    //// Construct the adjacency matrix
    int adjMatrix[ROWS * COLS][ROWS * COLS];
    for (int i = 0; i < ROWS * COLS; i++) {
        for (int j = 0; j < ROWS * COLS; j++) {
            adjMatrix[i][j] = (i == j) ? 0 : INT_MAX;
        }
    }
    //adjency matrix
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++) {
            Gmae& curr = grid[i][j];
            if (i > 0) {
                Gmae& up = grid[i - 1][j];
                int weight = (up.data == 'X') ? OBSTACLE_WEIGHT :
                    (up.data == 'E') ? ENEMY_WEIGHT : DEFAULT_WEIGHT;
                adjMatrix[curr.row * COLS + curr.col][up.row * COLS + up.col] = weight;
            }
            if (i < ROWS - 1) {
                Gmae& down = grid[i + 1][j];
                int weight = (down.data == 'X') ? OBSTACLE_WEIGHT :
                    (down.data == 'E') ? ENEMY_WEIGHT : DEFAULT_WEIGHT;
                adjMatrix[curr.row * COLS + curr.col][down.row * COLS + down.col] = weight;
            }
            if (j > 0) {
                Gmae& left = grid[i][j - 1];
                int weight = (left.data == 'X') ? OBSTACLE_WEIGHT :
                    (left.data == 'E') ? ENEMY_WEIGHT : DEFAULT_WEIGHT;
                adjMatrix[curr.row * COLS + curr.col][left.row * COLS + left.col] = weight;
            }
            if (j < COLS - 1) {
                Gmae& right = grid[i][j + 1];
                int weight = (right.data == 'X') ? OBSTACLE_WEIGHT :
                    (right.data == 'E') ? ENEMY_WEIGHT : DEFAULT_WEIGHT;
                adjMatrix[curr.row * COLS + curr.col][right.row * COLS + right.col] = weight;
            }
        }
    }

    // Print the adjacency matrix
    //printAdjacencyMatrix(map);

    //// Calculate the shortest path between any two areas
    Floyd_algo(adjMatrix);

    // Example usage: find the shortest path from the player's current position to the bottom-right corner of the grid
    Player player;

    int star_x_axes, star_y_axes;
    for (int i = 0; i < ROWS; i++)
    {
        for (int j = 0; j < COLS; j++)
        {
            if (map[i][j]->data == '*')
            {
                star_x_axes = i;
                star_y_axes = j;
                cout << "crystal location found at -> (" << i << ", " << j << ")\n\n";
            }
        }
    }
    int jewl = 0, weapon = 0, potion = 0, totalScore = 0;
    SocreAvlTree s;

    int c;
    cout << "1. Shortest Path through Floyd algorithm\n";
    cout << "2. Shortest Path through Dijkstra. algorithm\n";
    cout << "3. Minimum Spanning Tree through Kruskal's algorithm\n";
    cout << "4. Minimum Spanning Tree through Prim's algorithm\n";
    cout << "5. Play Game\n";
    cout << "6. Score Track\n";
    cout << "7. Quit\n";
    cin >> c;


    if (c == 1)
    {
        //________________________________________________FLOYDS ALGORITHM__________________________________________________________

        cout << "\n\t\t\tFLOYD'S ALgorithm\n" << endl;

        player.row = 0;
        player.col = 0;
        int destRow = star_x_axes;
        int destCol = star_y_axes;
        int shortestPath = adjMatrix[player.row * COLS + player.col][destRow * COLS + destCol];
        cout << "Shortest distance from (" << player.row << ", " << player.col << ") to (" << destRow << ", " << destCol << ") crystal: " << shortestPath << endl << endl;
        int x, y;
        cout << "Enter the no of row: " << endl;
        cin >> x;
        while (x > 20 || x < 0)
        {
            cout << "Enter the value of row again and less than 21" << endl;
            cin >> x;
        }

        cout << "Enter the no of column: " << endl;
        cin >> y;
        while (y > 20 || y < 0)
        {
            cout << "Enter the value of column again and less than 21" << endl;
            cin >> y;
        }

        player.row = x;
        player.col = y;
        destRow = star_x_axes;
        destCol = star_y_axes;
        shortestPath = adjMatrix[player.row * COLS + player.col][destRow * COLS + destCol];
        cout << "Shortest path from (" << player.row << ", " << player.col << ") to (" << destRow << ", " << destCol << "): " << shortestPath << endl;


        for (int i = 0; i < ROWS; i++)
        {
            for (int j = 0; j < COLS; j++)
            {
                if (map[i][j]->data != '*')
                {
                    if (map[i][j]->data == 'J')
                    {
                        jewl += 50;
                    }
                    if (map[i][j]->data == 'W')
                    {
                        weapon += 30;
                    }
                    if (map[i][j]->data == 'P')
                    {
                        potion += 70;
                    }
                }
                else
                {
                    break;
                }

            }
        }

        for (int i = 0; i < ROWS; i++)
        {
            for (int j = 0; j < COLS; j++)
            {
                if (map[i][j]->data != '*')
                {
                    if (map[i][j]->data == '@')
                    {
                        weapon -= 30;
                    }
                    if (map[i][j]->data == '$')
                    {
                        potion -= 70;
                    }
                    if (map[i][j]->data == '&')
                    {
                        jewl -= 50;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        totalScore = jewl + weapon + potion;

        s.insertNode(totalScore);
        s.printInorder();
    }

    else if (c == 2)
    {
        //---------------------------------------------------Dijkstra's ALGORITHM---------------------------------------------------------
       //printing the path for the dijkstra's algorithm.

        cout << "\n\t\t\tDijkstra's Algo\n" << endl;
        map_initialization_for_dijkstra();


        int startX, startY, endX, endY;
        player.row = 0;
        player.col = 0;
        int destRow = star_x_axes;
        int destCol = star_y_axes;
        int shortestPath = adjMatrix[player.row * COLS + player.col][destRow * COLS + destCol];
        cout << "Shortest distance from (" << player.row << ", " << player.col << ") to (" << destRow << ", " << destCol << ") crystal: " << shortestPath << endl << endl;

        cout << "Enter starting coordinates (row, col): ";
        cin >> startX >> startY;

        while (startX > 20 || startY > 21 || startX < 0 || startY < 0)
        {
            cout << "Enter the values again less than 21 and graeter or equal to zero" << endl;
            cin >> startX >> startY;

        }

        Dijkstra_algorithm(startX, startY, star_x_axes, star_y_axes);

        for (int i = 0; i < ROWS; i++)
        {
            for (int j = 0; j < COLS; j++)
            {
                if (map[i][j]->data != '*')
                {
                    if (map[i][j]->data == 'J')
                    {
                        jewl += 50;
                    }
                    if (map[i][j]->data == 'W')
                    {
                        weapon += 30;
                    }
                    if (map[i][j]->data == 'P')
                    {
                        potion += 70;
                    }
                }
                else
                {
                    break;
                }

            }
        }

        for (int i = 0; i < ROWS; i++)
        {
            for (int j = 0; j < COLS; j++)
            {
                if (map[i][j]->data != '*')
                {
                    if (map[i][j]->data == '@')
                    {
                        weapon -= 30;
                    }
                    if (map[i][j]->data == '$')
                    {
                        potion -= 70;
                    }
                    if (map[i][j]->data == '&')
                    {
                        jewl -= 50;
                    }
                }
                else
                {
                    break;
                }
            }
        }
        totalScore = jewl + weapon + potion;

        s.insertNode(totalScore);
        s.printInorder();
    }
    else if (c == 3)
    {
        //kuruskal Algo
        Gmae* map2[ROWS][COLS];
        map_reading_for_kuruskal_algo(map2);
        kruskal_Algorithm(map2);
    }
    else if (c == 4)
    {
        //prims
        char map1[ROWS][COLS];
        Map_reading_for_primsAlgo(map1, "map.txt");
        //printMap_PrimsAlgo(map1); 
        Gmae graph[ROWS][COLS];
        prim_Algo(graph);
        print_min_spanning_tree_for_prim_algo(graph);

    }

    else if (c == 5)
    {
        //Create the player object
       //Player player;

       // Main game loop
        while (map[player.row][player.col]->data != '*') {
            // Print the map
            for (int i = 0; i < ROWS; i++) {
                for (int j = 0; j < COLS; j++) {
                    if (i == player.row && j == player.col) {
                        cout << "M ";
                    }
                    else {
                        cout << map[i][j]->data << " ";
                    }
                }
                cout << endl;
            }
            // Get the player's move
            int move;
            cout << "Enter your move : \n";
            cout << "1. Up\n";
            cout << "2. Down\n";
            cout << "3. Left\n";
            cout << "4. Right\n";
            cin >> move;
            system("CLS");
            //using the clear command to clear the terminal on evry command.
            // Move the player
            switch (move) {
            case 1:
                moveUp(player, map);
                break;
            case 2:
                moveDown(player, map);
                break;
            case 3:
                moveLeft(player, map);
                break;
            case 4:
                moveRight(player, map);
                break;
            default:
                cout << "Invalid move" << endl;
                break;
            }

            // Interact with the cell
            interact(player, map);
        }

        // Print the victory message
        cout << "Congratulations, you retrieved the crystal!" << endl;
    }
    bool c1 = false, c2 = false;
    if (c == 6)
    {
        int jewlS = 0, weaponS = 0, potionS = 0, totalS = 0;
        for (int i = 0; i < ROWS; i++)
        {
            for (int j = 0; j < COLS; j++)
            {
                if (c1 == true)
                {
                    break;
                }
                if (map[i][j]->data != '*')
                {
                    if (map[i][j]->data == 'J')
                    {
                        jewlS += 50;
                        cout << map[i][j]->data << " Found at [" << i << "][" << j << "], Jewl score = " << jewlS << endl;
                    }
                    if (map[i][j]->data == 'W')
                    {
                        weaponS += 30;
                        cout << map[i][j]->data << " Found at [" << i << "][" << j << "], Weapon score = " << weaponS << endl;
                    }
                    if (map[i][j]->data == 'P')
                    {
                        potionS += 70;
                        cout << map[i][j]->data << " Found at [" << i << "][" << j << "], Potion score = " << potionS << endl;
                    }

                }
                if (map[i][j]->data == '*')
                {
                    c1 = true;
                    totalS = jewlS + weaponS + potionS;
                    cout << "Total Score = " << jewlS << " + " << weaponS << " + " << potionS << " = " << totalS << endl;
                    //s.insertNode(totalS);
                    //s.printInorder();
                    //exit(0);
                }

            }
        }

        cout << "\n\nScore deduction\n";
        for (int i = 0; i < ROWS; i++)
        {
            for (int j = 0; j < COLS; j++)
            {
                if (map[i][j]->data != '*')
                {
                    if (map[i][j]->data == '@')
                    {
                        weaponS -= 30;
                        cout << map[i][j]->data << " Found at [" << i << "][" << j << "], Weapon score = " << weaponS << endl;
                    }
                    if (map[i][j]->data == '$')
                    {
                        potionS -= 70;
                        cout << map[i][j]->data << " Found at [" << i << "][" << j << "], Potion score = " << potionS << endl;
                    }
                    if (map[i][j]->data == '&')
                    {
                        jewlS -= 50;
                        cout << map[i][j]->data << " Found at [" << i << "][" << j << "], Jewl score = " << jewlS << endl;
                    }

                }
                if (map[i][j]->data == '*')
                {
                    totalS = jewlS + weaponS + potionS;
                    /*cout << "Total Score = " << totalS+10 << endl;*/
                    s.insertNode(totalS + 10);
                    s.printInorder();
                    exit(0);
                }

            }
        }

    }
    if (c == 7)
    {
        cout << "Have a Good day.\n";


    }
    PlaySound(NULL, NULL, SND_FILENAME);

    // Join the sound thread to ensure it has finished
    soundThread.join();

    return 0;
}
