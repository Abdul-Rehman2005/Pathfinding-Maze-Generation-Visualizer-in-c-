#include <raylib.h>
#include <vector>
#include <queue>
#include <stack>
#include <random>
#include <algorithm>
#include <chrono>
#include <string>
#include <tuple>
#include <thread>
#include <cstdio>
#include <iostream>
using namespace std;

const int ROWS = 40;
const int COLS = 60;
const int CELL_SIZE = 20;

const int MAZE_WIDTH = COLS * CELL_SIZE;
const int MAZE_HEIGHT = ROWS * CELL_SIZE;
const int SIDEBAR_WIDTH = 300;
const int SCREEN_WIDTH = MAZE_WIDTH + SIDEBAR_WIDTH;
const int SCREEN_HEIGHT = MAZE_HEIGHT;


const Color WALL_COLOR = {44,62,80,255}; 
const Color PATH_COLOR = {241,196,15,255};
const Color VISITED_COLOR = {52,152,219,255};
const Color SOLUTION_COLOR = {155,89,182,255};
const Color START_COLOR = {46,204,113,255};
const Color END_COLOR = {231,76,60,255};
const Color UI_BG_COLOR = {25,25,25,255};
const Color BUTTON_COLOR = {52,73,94,255};
const Color BUTTON_HOVER_COLOR = {41,128,185,255};
const Color TEXT_COLOR = {236,240,241,255};

bool isSolving = false;
bool isPaused = false;
float animationDelay = 3; 

enum CELL_TYPE {CELL_WALL=1, CELL_PATH=0};

struct STATISTICS{
    double time_taken = 0.0;
    int nodes_visited = 0;
    int comparisons = 0;
    int path_length = 0;
    double memory = 0.0;

    void reset(){
        time_taken = 0.0;
        nodes_visited = 0;
        comparisons = 0;
        path_length = 0;
        memory = 0.0;
    }
};

struct BUTTON {
    Rectangle Bound;
    string text;
    bool ishovered=false;
    bool enabled=true;

    BUTTON(int x,int y,int w,int h,string t){
        Bound={float(x),float(y),float(w),float(h)};
        text=t;
    }

    bool isCLicked(){
        // make sure the button is now active
        if(!enabled) {
            return false;
        }
        if(IsMouseButtonPressed(MOUSE_LEFT_BUTTON) && CheckCollisionPointRec(GetMousePosition(), Bound)){
            return true;
        }
        return false;
    }

    void update(){ 
        // update the hovering effect of mouse on the button
        if(enabled && CheckCollisionPointRec(GetMousePosition(), Bound)){
            ishovered=true;
        }
        else{
            ishovered=false;
        }
    }
    void draw(){
        // nested tenray if stament to first see if button is enable and then to see if the ishovered true if so
        // set the button to their correponding color and draw the test;
        Color c = enabled?(ishovered? BUTTON_HOVER_COLOR: BUTTON_COLOR) : Color{100,100,100,255};
        DrawRectangleRec(Bound,c);
        DrawRectangleLinesEx(Bound,2,ishovered? WHITE:GRAY);
        int w = MeasureText(text.c_str(),20);
        DrawText(text.c_str(), Bound.x + (Bound.width-w)/2, Bound.y+(Bound.height-20)/2, 20, TEXT_COLOR);
    }
    void setEnabled(bool e){
        // set the current button to true or false depnding on whether it is being used or not
        enabled=e;
    }
};

struct MAZE {
    vector<vector<int>> grid;

    MAZE() { 
        //create a 2d vector array of size given and inilise it entirely to wall or 1
        grid.resize(ROWS, vector<int>(COLS, CELL_WALL));
    }

    bool isInBound(int row, int col) {
        //see if the coordiante are in the ui 
        if(row >= 0 && row < ROWS && col >= 0 && col < COLS){
            return true;
        }
        return false;
    }
    bool isWall(int row, int col) {
        //check if the cell is wall or not 
        if(grid[row][col] == CELL_WALL){
            return true;
        }
        return false;
    }
    void setPath(int row, int col) { 
        //set the cell to path by setting it to 0
        grid[row][col] = CELL_PATH;
    }
    void setWall(int row, int col) {
        // set the cell to wall by setting it 1
        grid[row][col] = CELL_WALL;
    }
    void addFrontier(int row, int col, vector<tuple<int,int,int>>& frontier) {
        //                       up    right  down  left
        int directions[4][2] = {{-1,0},{0,1},{1,0},{0,-1}};
        // check each direction from that specific WALL ndex!!!!
        for(int dir=0; dir<4; dir++){
            int newRow = row + directions[dir][0]*2;
            int newCol = col + directions[dir][1]*2;
            int wallRow = row + directions[dir][0];
            int wallCol = col + directions[dir][1];
            // for each direction move 2 block in that to see if that cell is a wall or not
            // and move one block and ahead to see if that cell is wall or not
            if(isInBound(newRow,newCol) && isWall(newRow,newCol)){
                // if the 2 block ahead is wall we need to carve t so for now start by storing the 1 block ahead 
                frontier.push_back(make_tuple(wallRow,wallCol,dir));
                // this makes it that the the current wall next index one block ahead index will now be converted into path
                // kee;ing the format of wall path wall
            }
        }
    }
    void GENERATE_PRIMS() {
        grid = vector<vector<int>>(ROWS, vector<int>(COLS, CELL_WALL));
        //
        // this is supposed to create a random number included in the range the diffenrce is that we have 
        // make sure the number is a int and in the list of 
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> rowDist(0, ROWS-1);
        uniform_int_distribution<> colDist(0, COLS-1);
        // set the start row and column with random number form rowgen and colgen
        int startRow = rowDist(gen);
        int startCol = colDist(gen);
        // set it to 0 
        setPath(startRow, startCol);
        //create a vector to store the row column and in which direction up down left right
        vector<tuple<int,int,int>> frontier;
        // call frontier function to find possible path to carve
        addFrontier(startRow, startCol, frontier);

        while(!frontier.empty()){
            // pick random wall from frontier 
            uniform_int_distribution<> idxDist(0, frontier.size()-1);
            int randomIndex = idxDist(gen);
            auto [row, col, dir] = frontier[randomIndex];
            // tuple unpacking
            frontier.erase(frontier.begin() + randomIndex);
            // this is making sure that the random index you got is now delted from frontier
            // so it can be marked visitd and cannot be again added 
            int neighborRow=row; 
            int neighborCol=col;
            switch(dir){ 
                case 0:
                    neighborRow--;
                     break;
                case 1:
                    neighborCol++;
                    break;
                case 2:
                    neighborRow++;
                    break;
                case 3:
                    neighborCol--;
                    break; }
            //now we will chose next cell by seeing in which direction was it orignally be carved 
            //
            if(isInBound(neighborRow, neighborCol) && isWall(neighborRow, neighborCol)){
                // see if this valid and now make the the origanla and new cell both as path to make them connected
                setPath(row, col);
                setPath(neighborRow, neighborCol);
                // add the new cell to start from there as well;
                addFrontier(neighborRow, neighborCol, frontier);
            }
        }
            //How this avoids revisiting cells
            //Only walls that connect to an unvisited cell are added to the frontier in addFrontie.
            //The second check ensures we don’t carve into a cell that is already a path.
            //Therefore, every new cell is carved exactly once, and the maze grows like a tree no cycles.

        // mark these 2 as path regqrdless if they are wall or not
        setPath(0,0);
        setPath(ROWS-1,COLS-1);
    }
//
    int findParent(vector<int> &parent, int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]]; //path compression
            x = parent[x];
        }
        return x;
    }

    bool uniteCells(vector<int> &parent, int x, int y) {
        int rootX = findParent(parent, x);
        int rootY = findParent(parent, y);
        if (rootX != rootY) {
            parent[rootY] = rootX;
            return true;
        }
        return false;
    }

    int getCellIndex(vector<pair<int,int>> &cells, int r, int c) {
        for (int i = 0; i < cells.size(); i++) {
            if (cells[i].first == r && cells[i].second == c) {
                return i;
            }
        }
        return -1;
    }
// implmenetation of Union-Find (also called Disjoint Set Union, DSU) 
// is a data structure to keep track of which elements belong to which set.
// Each cell in the maze is considered a node.
// Initially, each cell is in its own separate set  no connections.
// When we remove a wall between two cells, we merge their sets → they are now connected.
// This helps prevent cycles (loops) in the maze. 
    void generateKruskals() {
        grid = vector<vector<int>>(ROWS, vector<int>(COLS, CELL_WALL));
        //create cells (odd rows/cols)
        vector<pair<int,int>> cells;
        for (int r = 1; r < ROWS; r += 2) {
            for (int c = 1; c < COLS; c += 2) {
                cells.push_back(make_pair(r, c));
            }
        }
        //create edges (horizontal and vertical)
        struct Edge { 
            int r1;
            int c1;
            int r2;
            int c2;
        };
        vector<Edge> edges;
        for (int r = 1; r < ROWS; r += 2) {
            for (int c = 1; c < COLS - 2; c += 2) {
                Edge e = {r, c, r, c + 2};
                edges.push_back(e);
            }
        }
        for (int r = 1; r < ROWS - 2; r += 2) {
            for (int c = 1; c < COLS; c += 2) {
                Edge e = {r, c, r + 2, c};
                edges.push_back(e);
            }
        }
        //shuffle edges randomly
        random_device rd;
        mt19937 gen(rd());
        shuffle(edges.begin(), edges.end(), gen);
        //initialize union-find parent
        vector<int> parent(cells.size());
        for (int i = 0; i < cells.size(); i++) parent[i] = i;
        //mark all cells as path
        for (int i = 0; i < cells.size(); i++) {
            setPath(cells[i].first, cells[i].second);
        }
        //process edges
        for (int i = 0; i < edges.size(); i++) {
            int idx1 = getCellIndex(cells, edges[i].r1, edges[i].c1);
            int idx2 = getCellIndex(cells, edges[i].r2, edges[i].c2);
            if (idx1 != -1 && idx2 != -1 && uniteCells(parent, idx1, idx2)) {
                int wallR = (edges[i].r1 + edges[i].r2) / 2;
                int wallC = (edges[i].c1 + edges[i].c2) / 2;
                setPath(wallR, wallC); // remove the wall
            }
        }
        //start and end are paths
        setPath(0, 0);
        setPath(ROWS - 1, COLS - 1);
    }

        void draw() {
            for(int r=0;r<ROWS;r++)
                for(int c=0;c<COLS;c++){
                    Color color = isWall(r,c)? PATH_COLOR: WALL_COLOR; 
                    DrawRectangle(c*CELL_SIZE, r*CELL_SIZE, CELL_SIZE, CELL_SIZE, color);
                    if(isWall(r,c)) {
                        DrawRectangleLines(c*CELL_SIZE, r*CELL_SIZE, CELL_SIZE, CELL_SIZE, WALL_COLOR);
                    }
                }
            DrawRectangle(0,0,CELL_SIZE,CELL_SIZE,START_COLOR);
            DrawRectangleLines(0,0,CELL_SIZE,CELL_SIZE,PATH_COLOR);

            DrawRectangle((COLS-1)*CELL_SIZE,(ROWS-1)*CELL_SIZE,CELL_SIZE,CELL_SIZE,END_COLOR);
            DrawRectangleLines((COLS-1)*CELL_SIZE,(ROWS-1)*CELL_SIZE,CELL_SIZE,CELL_SIZE,PATH_COLOR);
        }

    int reconstructPath(const vector<vector<pair<int,int>>>& parent, pair<int,int> end){
        vector<pair<int,int>> path;
        pair<int,int> current = end;
        while(current.first!=-1){
            path.push_back(current);
            current = parent[current.first][current.second];
        }
        
        for(int i=path.size()-1;i>=0;i--){
            auto [r,c]=path[i];
            if(!(r==0&&c==0) && !(r==ROWS-1&&c==COLS-1)){
                DrawRectangle(c*CELL_SIZE,r*CELL_SIZE,CELL_SIZE,CELL_SIZE,SOLUTION_COLOR);
                DrawRectangleLines(c*CELL_SIZE,r*CELL_SIZE,CELL_SIZE,CELL_SIZE,PATH_COLOR);
                BeginDrawing(); 
                EndDrawing();
            }
        }
        return path.size()-1;
    }

    bool bfs(STATISTICS& stats){
        auto starttime = chrono::high_resolution_clock::now();
        pair<int,int> start={1,1};
        vector<pair<int,int>> end= {{ROWS-2,COLS-2},{ROWS-1,COLS-2},{ROWS-2,COLS-1},{ROWS-1,COLS-1}};

        stats.reset(); stats.comparisons=0;

        queue<pair<int,int>> q;
        vector<vector<bool>> visited(ROWS, vector<bool>(COLS,false));
        vector<vector<pair<int,int>>> parent(ROWS, vector<pair<int,int>>(COLS,{-1,-1}));
        q.push(start); visited[start.first][start.second]=true; stats.nodes_visited=1;

        int dirs[4][2]={{-1,0},{1,0},{0,-1},{0,1}};

        while(!q.empty()){
            if(isPaused){ this_thread::sleep_for(chrono::milliseconds(1)); continue; }
            auto [row,col] = q.front(); q.pop();

            for(auto X:end){
                if(row==X.first && col==X.second){
                stats.path_length = reconstructPath(parent,X);
                stats.time_taken = chrono::duration<double>(chrono::high_resolution_clock::now()-starttime).count();
                return true;
                }
            }

            for(int i=0;i<4;i++){
                int nr=row+dirs[i][0];
                int nc=col+dirs[i][1];
                stats.comparisons++;
                if(isInBound(nr,nc) && !isWall(nr,nc) && !visited[nr][nc]){
                    q.push({nr,nc}); visited[nr][nc]=true;
                    parent[nr][nc]={row,col}; stats.nodes_visited++;
                    DrawRectangle(nc*CELL_SIZE,nr*CELL_SIZE,CELL_SIZE,CELL_SIZE,VISITED_COLOR);
                    DrawRectangleLines(nc*CELL_SIZE,nr*CELL_SIZE,CELL_SIZE,CELL_SIZE,PATH_COLOR);
                    BeginDrawing(); EndDrawing();
                    this_thread::sleep_for(chrono::milliseconds((int)animationDelay));
                }
            }
        }
        return false;
    }

    bool dfs(STATISTICS& stats){
        auto starttime = chrono::high_resolution_clock::now();
        pair<int,int> start ={1,1};
        vector<pair<int,int>> end= {{ROWS-2,COLS-2},{ROWS-1,COLS-2},{ROWS-2,COLS-1},{ROWS-1,COLS-1}};
        // pair<int,int> end1 = {ROWS-2,COLS-2};
        // pair<int,int> end2 = {ROWS-1,COLS-2};
        // pair<int,int> end3 = {ROWS-2,COLS-1};
        // pair<int,int> end4 = {ROWS-1,COLS-1};
        // end if rows -2 col-2 or rows-1 col-2 or row-2 col-1 or row-1 col-1

        stats.reset(); stats.comparisons=0;

        stack<pair<int,int>> st;
        vector<vector<bool>> visited(ROWS, vector<bool>(COLS,false));
        vector<vector<pair<int,int>>> parent(ROWS, vector<pair<int,int>>(COLS,{-1,-1}));
        st.push(start); visited[start.first][start.second]=true; stats.nodes_visited=1;

        int dirs[4][2]={{-1,0},{1,0},{0,-1},{0,1}};

        while(!st.empty()){
            if(isPaused){ this_thread::sleep_for(chrono::milliseconds(1)); continue; }
            auto [row,col] = st.top(); st.pop();

            for(auto X:end){
                if(row==X.first && col==X.second){
                stats.path_length = reconstructPath(parent,X);
                stats.time_taken = chrono::duration<double>(chrono::high_resolution_clock::now()-starttime).count();
                return true;
                }
            }


            for(int i=0;i<4;i++){
                int nr=row+dirs[i][0], nc=col+dirs[i][1];
                stats.comparisons++;
                if(isInBound(nr,nc) && !isWall(nr,nc) && !visited[nr][nc]){
                    st.push({nr,nc}); visited[nr][nc]=true;
                    parent[nr][nc]={row,col}; stats.nodes_visited++;

                    DrawRectangle(nc*CELL_SIZE,nr*CELL_SIZE,CELL_SIZE,CELL_SIZE,VISITED_COLOR);
                    DrawRectangleLines(nc*CELL_SIZE,nr*CELL_SIZE,CELL_SIZE,CELL_SIZE,PATH_COLOR);
                    BeginDrawing(); EndDrawing();
                    this_thread::sleep_for(chrono::milliseconds((int)animationDelay));
                }
            }
        }
        return false;
    }

    bool dijkstra(STATISTICS& stats){
        auto starttime = chrono::high_resolution_clock::now();
        pair<int,int> start={1,1};
        vector<pair<int,int>> end= {{ROWS-2,COLS-2},{ROWS-1,COLS-2},{ROWS-2,COLS-1},{ROWS-1,COLS-1}};

        stats.reset(); 
        stats.comparisons=0;

        priority_queue<tuple<int,int,int>,vector<tuple<int,int,int>>,greater<tuple<int,int,int>>> pq;
        vector<vector<int>> dist(ROWS,vector<int>(COLS,INT_MAX));
        vector<vector<bool>> visited(ROWS,vector<bool>(COLS,false));
        vector<vector<pair<int,int>>> parent(ROWS,vector<pair<int,int>>(COLS,{-1,-1}));

        dist[start.first][start.second]=0;
        pq.push({0,start.first,start.second}); 
        stats.nodes_visited=1;
        int dirs[4][2]={{-1,0},{1,0},{0,-1},{0,1}};

        while(!pq.empty()){
            if(isPaused){ this_thread::sleep_for(chrono::milliseconds(1)); continue; }
            auto [currentdist,row,col] = pq.top(); pq.pop();
            if(visited[row][col]) continue;
            visited[row][col]=true;

            for(auto X:end){
                if(row==X.first && col==X.second){
                stats.path_length = reconstructPath(parent,X);
                stats.time_taken = chrono::duration<double>(chrono::high_resolution_clock::now()-starttime).count();
                return true;
                }
            }

            for(int i=0;i<4;i++){
                int nr=row+dirs[i][0], nc=col+dirs[i][1];
                stats.comparisons++;
                if(isInBound(nr,nc) && !isWall(nr,nc)){
                    int newdist=currentdist+1;
                    if(newdist<dist[nr][nc]){
                        dist[nr][nc]=newdist;
                        parent[nr][nc]={row,col};
                        pq.push({newdist,nr,nc});
                        stats.nodes_visited++;

                        DrawRectangle(nc*CELL_SIZE,nr*CELL_SIZE,CELL_SIZE,CELL_SIZE,VISITED_COLOR);
                        DrawRectangleLines(nc*CELL_SIZE,nr*CELL_SIZE,CELL_SIZE,CELL_SIZE,PATH_COLOR);
                        BeginDrawing(); EndDrawing();
                        this_thread::sleep_for(chrono::milliseconds((int)animationDelay));
                    }
                }
            }
        }
        return false;
    }
};

void drawStats(const STATISTICS &stats){
    // use the globat statitc variable to draw the result of the stats;
    DrawRectangle(MAZE_WIDTH,0,SIDEBAR_WIDTH,SCREEN_HEIGHT,UI_BG_COLOR);
    DrawText(("Nodes visited: "+to_string(stats.nodes_visited)).c_str(), MAZE_WIDTH+20, 420, 20, TEXT_COLOR);
    DrawText(("Comparisons: "+to_string(stats.comparisons)).c_str(), MAZE_WIDTH+20, 450, 20, TEXT_COLOR);
    DrawText(("Path length: "+to_string(stats.path_length)).c_str(), MAZE_WIDTH+20, 480, 20, TEXT_COLOR);
    DrawText(("Time: "+to_string(stats.time_taken)).c_str(), MAZE_WIDTH+20, 510, 20, TEXT_COLOR);
}
int main(){
    InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Maze Solver with Raylib");
    SetTargetFPS(60);

    MAZE maze;
    maze.GENERATE_PRIMS();

    STATISTICS stats;

    enum MAZE_ALGORITHM{MAZE_PRIMS, MAZE_KRUSKALS};
    enum SEARCH_ALGORITHM{SEARCH_BFS, SEARCH_DFS, SEARCH_DIJKSTRA};

    MAZE_ALGORITHM currentMazeAlgorithm = MAZE_PRIMS;
    SEARCH_ALGORITHM currentSearchAlgorithm = SEARCH_BFS;

    // Buttons
    vector<BUTTON> buttons;
    buttons.push_back(BUTTON(MAZE_WIDTH + 20, 50, 250, 40, "Generate Prim's"));
    buttons.push_back(BUTTON(MAZE_WIDTH + 20, 100, 250, 40, "Generate Kruskal's"));
    buttons.push_back(BUTTON(MAZE_WIDTH + 20, 150, 250, 50, "Start Solving"));

    vector<BUTTON> algoButtons;
    algoButtons.push_back(BUTTON(MAZE_WIDTH + 20, 310, 80, 35, "BFS"));
    algoButtons.push_back(BUTTON(MAZE_WIDTH + 110, 310, 80, 35, "DFS"));
    algoButtons.push_back(BUTTON(MAZE_WIDTH + 200, 310, 80, 35, "Dijkstra"));

    while(!WindowShouldClose()){
        for(auto &btn: buttons) btn.update();
        for(auto &btn: algoButtons) btn.update();

        if(buttons[0].isCLicked()){ 
            maze.GENERATE_PRIMS(); 
            currentMazeAlgorithm = MAZE_PRIMS;  
            stats.reset(); 
            isSolving=false; 
        }
        if(buttons[1].isCLicked()){ 
            maze.generateKruskals(); 
            currentMazeAlgorithm = MAZE_KRUSKALS;  
            stats.reset(); 
            isSolving=false; 
        }
        if(buttons[2].isCLicked()){ 
            stats.reset(); 
            isSolving=true; 
        }

        if(algoButtons[0].isCLicked()){
            currentSearchAlgorithm = SEARCH_BFS;
        }
        if(algoButtons[1].isCLicked()){
            currentSearchAlgorithm = SEARCH_DFS;
        }
        if(algoButtons[2].isCLicked()) {
            currentSearchAlgorithm = SEARCH_DIJKSTRA;
        }

        BeginDrawing();
        ClearBackground(UI_BG_COLOR);

        // Draw maze area
        DrawRectangle(0, 0, MAZE_WIDTH, MAZE_HEIGHT, UI_BG_COLOR);
        maze.draw();

        // Draw sidebar background FIRST
        DrawRectangle(MAZE_WIDTH, 0, SIDEBAR_WIDTH, SCREEN_HEIGHT, UI_BG_COLOR);
        
        // Draw sidebar title
        DrawText("MAZE SOLVER", MAZE_WIDTH + 20, 20, 30, TEXT_COLOR);
        DrawRectangle(MAZE_WIDTH + 15, 55, SIDEBAR_WIDTH - 30, 2, BUTTON_COLOR);
        
        // Draw buttons
        for(auto &btn: buttons) btn.draw();
        for(auto &btn: algoButtons) btn.draw();
        
        // Draw algorithm selection title
        DrawText("Search Algorithm:", MAZE_WIDTH + 20, 280, 20, TEXT_COLOR);
        
        // Highlight selected algorithm
        Color highlightColor = {41, 128, 185, 100};
        if (currentSearchAlgorithm == SEARCH_BFS) {
            DrawRectangleRec(algoButtons[0].Bound, highlightColor);
        } else if (currentSearchAlgorithm == SEARCH_DFS) {
            DrawRectangleRec(algoButtons[1].Bound, highlightColor);
        } else if (currentSearchAlgorithm == SEARCH_DIJKSTRA) {
            DrawRectangleRec(algoButtons[2].Bound, highlightColor);
        }
        
        DrawText("Current Maze:", MAZE_WIDTH + 20, 360, 20, TEXT_COLOR);
        string mazeType;
        if(currentMazeAlgorithm == MAZE_PRIMS){
            mazeType = "Prim's Algorithm";
        }
        else if(currentMazeAlgorithm == MAZE_KRUSKALS){  
            mazeType = "Kruskal's Algorithm";
        }
        DrawText(mazeType.c_str(), MAZE_WIDTH + 20, 390, 20, GREEN);

        
        DrawText("Statistics:", MAZE_WIDTH + 20, 430, 20, TEXT_COLOR);
        
        char timeText[50];
        sprintf(timeText, "Time: %.3fs", stats.time_taken);
        DrawText(timeText, MAZE_WIDTH + 20, 460, 18, TEXT_COLOR);
        
        char nodesText[50];
        sprintf(nodesText, "Nodes Visited: %d", stats.nodes_visited);
        DrawText(nodesText, MAZE_WIDTH + 20, 485, 18, TEXT_COLOR);
        
        char pathText[50];
        sprintf(pathText, "Path Length: %d", stats.path_length);
        DrawText(pathText, MAZE_WIDTH + 20, 510, 18, TEXT_COLOR);
        
        char compText[50];
        sprintf(compText, "Comparisons: %d", stats.comparisons);
        DrawText(compText, MAZE_WIDTH + 20, 535, 18, TEXT_COLOR);

        EndDrawing();

        if(isSolving){
            switch(currentSearchAlgorithm){
                case SEARCH_BFS:
                    maze.bfs(stats);
                    break;
                case SEARCH_DFS: 
                    maze.dfs(stats); 
                    break;
                case SEARCH_DIJKSTRA: 
                    maze.dijkstra(stats); 
                    break;
            }
            isSolving=false; 
        }
    }

    CloseWindow();
    return 0;
}