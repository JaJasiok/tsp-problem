#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <unordered_set>

using namespace std;

vector<vector<int>> readKroaFile(const string &filename)
{
    ifstream file(filename);
    string line;
    vector<vector<int>> verticesCoords;

    if (file.is_open())
    {
        // Skip the header lines
        for (int i = 0; i < 6; i++)
        {
            getline(file, line);
        }

        // Read the coordinates and populate the cost matrix
        while (getline(file, line) && line != "EOF")
        {
            istringstream iss(line);
            int index, x, y;
            iss >> index >> x >> y;
            verticesCoords.push_back({x, y});
        }

        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }

    return verticesCoords;
}

vector<vector<int>> createDistanceMatrix(const vector<vector<int>> &verticesCoords)
{
    int numVertices = verticesCoords.size();
    vector<vector<int>> distanceMatrix(numVertices, vector<int>(numVertices));

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            double x1 = verticesCoords[i][0];
            double y1 = verticesCoords[i][1];
            double x2 = verticesCoords[j][0];
            double y2 = verticesCoords[j][1];

            double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
            distanceMatrix[i][j] = int(round(distance));
        }
    }

    return distanceMatrix;
}

class Vertex
{
public:
    int id;
    vector<Vertex *> neighbours;

    Vertex(int _id) : id(_id) {}
};

class Edge
{
public:
    Vertex *src;
    Vertex *dest;
    int distance;

    Edge(Vertex *_src, Vertex *_dest, int _distance) : src(_src), dest(_dest), distance(_distance) {}
};

class Graph
{
public:
    vector<Vertex *> vertices;
    vector<Edge *> edges;
    int distance = 0;

    void addVertex(Vertex *v)
    {
        if (v)
        {
            vertices.push_back(v);
        }
    }

    void addEdge(Vertex *src, Vertex *dest, int distance = 0)
    {
        if (src && dest)
        {
            Edge *e = new Edge(src, dest, distance);
            src->neighbours.push_back(dest);
            dest->neighbours.push_back(src);
            this->distance += distance;
            edges.push_back(e);
        }
    }

    void removeVertex(Vertex *v)
    {
        if (v)
        {
            vertices.erase(remove(vertices.begin(), vertices.end(), v), vertices.end());
            for (Vertex *neighbour : v->neighbours)
            {
                neighbour->neighbours.erase(remove(neighbour->neighbours.begin(), neighbour->neighbours.end(), v), neighbour->neighbours.end());
            }
            for (Edge *e : edges)
            {
                if (e->src == v || e->dest == v)
                {
                    removeEdge(e->src, e->dest);
                }
            }
        }
    }

    void removeEdge(Vertex *src, Vertex *dest)
    {
        if (src && dest)
        {
            src->neighbours.erase(remove(src->neighbours.begin(), src->neighbours.end(), dest), src->neighbours.end());
            dest->neighbours.erase(remove(dest->neighbours.begin(), dest->neighbours.end(), src), dest->neighbours.end());
            for (Edge *e : edges)
            {
                if (e->src == src && e->dest == dest)
                {
                    edges.erase(remove(edges.begin(), edges.end(), e), edges.end());
                    this->distance -= e->distance;
                }
            }
        }
    }

    Vertex *findVertex(int id)
    {
        for (Vertex *v : vertices)
        {
            if (v->id == id)
                return v;
        }
        return nullptr;
    }

    bool hasCycle()
    {
        unordered_set<Vertex *> visited;
        for (Vertex *v : vertices)
        {
            if (!visited.count(v) && hasCycleUtil(v, nullptr, visited))
                return true;
        }
        return false;
    }

    bool hasCycleUtil(Vertex *v, Vertex *parent, unordered_set<Vertex *> &visited)
    {
        visited.insert(v);
        for (Vertex *neighbour : v->neighbours)
        {
            if (!visited.count(neighbour))
            {
                if (hasCycleUtil(neighbour, v, visited))
                    return true;
            }
            else if (neighbour != parent)
            {
                return true;
            }
        }
        return false;
    }
};

void saveGraphs(const vector<Graph> &graphs, const string &filename)
{
    ofstream file(filename);
    if (file.is_open())
    {
        for (Graph g : graphs)
        {
            for (Edge *e : g.edges)
            {
                file << e->src->id << " " << e->dest->id << endl;
            }
            file << endl;
        }
        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }
}

vector<Graph> NNCycles(const vector<vector<int>> &distanceMatrix, int startId)
{
    int numVertices = distanceMatrix.size();
    vector<Graph> cycles(2);
    vector<bool> visited(numVertices, false);

    // Choose the first vertex randomly
    Vertex *startVertex1 = new Vertex(startId);
    cycles[0].addVertex(startVertex1);
    visited[startVertex1->id] = true;

    // Find the furthest vertex from startVertex1
    int furthestVertex = -1;
    int maxDistance = -1;
    for (int j = 0; j < numVertices; j++)
    {
        if (!visited[j] && distanceMatrix[startVertex1->id][j] > maxDistance)
        {
            maxDistance = distanceMatrix[startVertex1->id][j];
            furthestVertex = j;
        }
    }

    Vertex *startVertex2 = new Vertex(furthestVertex);
    cycles[1].addVertex(startVertex2);
    visited[startVertex2->id] = true;

    // Find the nearest neighbours
    for (int i = 0; i < numVertices - 2; i++)
    {
        int minDistance = INT_MAX;
        int nearestNeighbourId = -1;
        pair<Vertex *, Graph *> minPair;

        vector<pair<Vertex *, Graph *>> lastVertices;
        if (cycles[0].vertices.size() == numVertices / 2)
        {
            for (Vertex *v : cycles[1].vertices)
            {
                if (v->neighbours.size() < 2)
                {
                    lastVertices.push_back({v, &cycles[1]});
                }
            }
        }
        else if (cycles[1].vertices.size() == numVertices / 2)
        {
            for (Vertex *v : cycles[0].vertices)
            {
                if (v->neighbours.size() < 2)
                {
                    lastVertices.push_back({v, &cycles[0]});
                }
            }
        }
        else
        {
            for (Graph &cycle : cycles)
            {
                for (Vertex *v : cycle.vertices)
                {
                    if (v->neighbours.size() < 2)
                    {
                        lastVertices.push_back({v, &cycle});
                    }
                }
            }
        }

        for (pair<Vertex *, Graph *> singlePair : lastVertices)
        {
            for (int j = 0; j < numVertices; j++)
            {
                if (visited[j] == false)
                {
                    Vertex *v = singlePair.first;
                    if (distanceMatrix[v->id][j] < minDistance)
                    {
                        minDistance = distanceMatrix[v->id][j];
                        nearestNeighbourId = j;
                        minPair.first = v;
                        minPair.second = singlePair.second;
                    }
                }
            }
        }

        Graph *cyclePtr = minPair.second;
        Graph &cycle = *cyclePtr;

        Vertex *vertex = minPair.first;

        Vertex *nearestNeighbour = new Vertex(nearestNeighbourId);
        cycle.addVertex(nearestNeighbour);
        cycle.addEdge(vertex, nearestNeighbour, minDistance);

        visited[nearestNeighbourId] = true;
    }

    for (Graph &cycle : cycles)
    {
        vector<Vertex *> verticesWithOneNeighbour;
        for (Vertex *v : cycle.vertices)
        {
            if (v->neighbours.size() == 1)
            {
                verticesWithOneNeighbour.push_back(v);
            }
        }

        Vertex *vertex1 = verticesWithOneNeighbour[0];
        Vertex *vertex2 = verticesWithOneNeighbour[1];

        cycle.addEdge(vertex1, vertex2, distanceMatrix[vertex1->id][vertex2->id]);
    }

    return cycles;
}

vector<Graph> greedyCycles(const vector<vector<int>> &distanceMatrix, int startId)
{
    int numVertices = distanceMatrix.size();
    vector<Graph> cycles(2);
    vector<bool> visited(numVertices, false);

    // Choose the first vertex randomly
    Vertex *startVertex1 = new Vertex(startId);
    cycles[0].addVertex(startVertex1);
    visited[startVertex1->id] = true;

    // Find the furthest vertex from startVertex1
    int furthestVertex = -1;
    int maxDistance = -1;
    for (int j = 0; j < numVertices; j++)
    {
        if (!visited[j] && distanceMatrix[startVertex1->id][j] > maxDistance)
        {
            maxDistance = distanceMatrix[startVertex1->id][j];
            furthestVertex = j;
        }
    }

    Vertex *startVertex2 = new Vertex(furthestVertex);
    cycles[1].addVertex(startVertex2);
    visited[startVertex2->id] = true;

    // Find the nearest neighbour for each starting vertex
    for (Graph &cycle : cycles)
    {
        Vertex *vertex = cycle.vertices.front();
        int vertexId = vertex->id;
        int minDistance = INT_MAX;
        int nearestNeighbourId = -1;

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                if (distanceMatrix[vertexId][j] < minDistance)
                {
                    minDistance = distanceMatrix[vertexId][j];
                    nearestNeighbourId = j;
                }
            }
        }

        Vertex *nearestNeighbour = new Vertex(nearestNeighbourId);

        cycle.addVertex(nearestNeighbour);
        cycle.addEdge(vertex, nearestNeighbour, minDistance);
        visited[nearestNeighbourId] = true;
    }

    // Building the rest of the cycle
    for (int i = 0; i < numVertices - 4; i++)
    {
        int minDistance = INT_MAX;
        int vertexId = -1;
        pair<Edge *, Graph *> minPair;

        vector<pair<Edge *, Graph *>> edgesInGraphs;
        if (cycles[0].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[1].edges)
            {
                edgesInGraphs.push_back({e, &cycles[1]});
            }
        }
        else if (cycles[1].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[0].edges)
            {
                edgesInGraphs.push_back({e, &cycles[0]});
            }
        }
        else
        {
            for (Graph &cycle : cycles)
            {
                for (Edge *e : cycle.edges)
                {
                    edgesInGraphs.push_back({e, &cycle});
                }
            }
        }

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                for (pair<Edge *, Graph *> singlePair : edgesInGraphs)
                {
                    Edge *e = singlePair.first;
                    int distanceSum = cycles[0].distance + cycles[1].distance;
                    if (distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j] < minDistance)
                    {
                        minDistance = distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j];
                        vertexId = j;
                        minPair = singlePair;
                    }
                }
            }
        }

        Graph *cyclePtr = minPair.second;
        Graph &cycle = *cyclePtr;

        Edge *minEdge = minPair.first;

        Vertex *newVertex = new Vertex(vertexId);
        cycle.addVertex(newVertex);
        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);

        if (cycle.vertices.size() > 3)
        {
            cycle.removeEdge(minEdge->src, minEdge->dest);
        }

        visited[vertexId] = true;
    }
    return cycles;
}

vector<Graph> regretCycles(const vector<vector<int>> &distanceMatrix, int stardId)
{
    int numVertices = distanceMatrix.size();
    vector<Graph> cycles(2);
    vector<bool> visited(numVertices, false);

    // Choose the first vertex randomly
    Vertex *startVertex1 = new Vertex(stardId);
    cycles[0].addVertex(startVertex1);
    visited[startVertex1->id] = true;

    // Find the furthest vertex from startVertex1
    int furthestVertex = -1;
    int maxDistance = -1;
    for (int j = 0; j < numVertices; j++)
    {
        if (!visited[j] && distanceMatrix[startVertex1->id][j] > maxDistance)
        {
            maxDistance = distanceMatrix[startVertex1->id][j];
            furthestVertex = j;
        }
    }

    Vertex *startVertex2 = new Vertex(furthestVertex);
    cycles[1].addVertex(startVertex2);
    visited[startVertex2->id] = true;

    // Find the nearest neighbour for each starting vertex
    for (Graph &cycle : cycles)
    {
        Vertex *vertex = cycle.vertices.front();
        int vertexId = vertex->id;
        int minDistance = INT_MAX;
        int nearestNeighbourId = -1;

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                if (distanceMatrix[vertexId][j] < minDistance && distanceMatrix[vertexId][j] != 0)
                {
                    minDistance = distanceMatrix[vertexId][j];
                    nearestNeighbourId = j;
                }
            }
        }

        Vertex *nearestNeighbour = new Vertex(nearestNeighbourId);

        cycle.addVertex(nearestNeighbour);
        cycle.addEdge(vertex, nearestNeighbour, minDistance);
        visited[nearestNeighbourId] = true;
    }

    // Selecting third vertices
    for (Graph &cycle : cycles)
    {
        int minDistance = INT_MAX;
        int vertexId = -1;
        Edge *minEdge = nullptr;

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                for (Edge *e : cycle.edges)
                {
                    int distance = cycle.distance;
                    if (distance - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j] < minDistance)
                    {
                        minDistance = distance - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j];
                        vertexId = j;
                        minEdge = e;
                    }
                }
            }
        }

        Vertex *newVertex = new Vertex(vertexId);
        cycle.addVertex(newVertex);
        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);

        visited[vertexId] = true;
    }

    // Building the rest of the cycle
    for (int i = 0; i < numVertices - 6; i++)
    {
        int vertexId = -1;
        vector<pair<int, pair<int, pair<Edge *, Graph *>>>> regrets;
        pair<Edge *, Graph *> minPair;

        vector<pair<Edge *, Graph *>> edgesInGraphs;
        if (cycles[0].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[1].edges)
            {
                edgesInGraphs.push_back({e, &cycles[1]});
            }
        }
        else if (cycles[1].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[0].edges)
            {
                edgesInGraphs.push_back({e, &cycles[0]});
            }
        }
        else
        {
            for (Graph &cycle : cycles)
            {
                for (Edge *e : cycle.edges)
                {
                    edgesInGraphs.push_back({e, &cycle});
                }
            }
        }

        for (int j = 0; j < numVertices; j++)
        {
            vector<pair<int, pair<Edge *, Graph *>>> costs;
            if (!visited[j])
            {
                for (pair<Edge *, Graph *> singlePair : edgesInGraphs)
                {
                    Edge *e = singlePair.first;
                    costs.push_back(make_pair(distanceMatrix[e->src->id][j] + distanceMatrix[e->dest->id][j] - e->distance, singlePair));
                }
                sort(costs.begin(), costs.end(), less<pair<int, pair<Edge *, Graph *>>>());
                regrets.push_back({abs(costs[0].first - costs[1].first), {j, costs.front().second}});
            }
        }

        sort(regrets.begin(), regrets.end(), greater<pair<int, pair<int, pair<Edge *, Graph *>>>>());

        vertexId = regrets.front().second.first;
        minPair = regrets.front().second.second;

        Graph *cyclePtr = minPair.second;
        Graph &cycle = *cyclePtr;

        Edge *minEdge = minPair.first;

        Vertex *newVertex = new Vertex(vertexId);
        cycle.addVertex(newVertex);
        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);
        cycle.removeEdge(minEdge->src, minEdge->dest);

        visited[vertexId] = true;
    }

    return cycles;
}

int main()
{
    vector<vector<int>> verticesCoords = readKroaFile("kroA100.tsp");

    vector<vector<int>> distanceMatrix = createDistanceMatrix(verticesCoords);

    pair<int, int> bestNNValue = {INT_MAX, -1};
    pair<int, int> bestGreedyValue = {INT_MAX, -1};
    pair<int, int> bestRegretValue = {INT_MAX, -1};

    pair<int, int> worstNNValue = {-1, -1};
    pair<int, int> worstGreedyValue = {-1, -1};
    pair<int, int> worstRegretValue = {-1, -1};

    long averageNNValue = 0;
    long averageGreedyValue = 0;
    long averageRegretValue = 0;

    for (int i = 0; i < verticesCoords.size(); i++)
    {

        vector<Graph> NNCyclesResult = NNCycles(distanceMatrix, i);
        vector<Graph> greedyCyclesResult = greedyCycles(distanceMatrix, i);
        vector<Graph> regretCyclesResult = regretCycles(distanceMatrix, i);

        if (NNCyclesResult[0].distance + NNCyclesResult[1].distance < bestNNValue.first)
        {
            bestNNValue.first = NNCyclesResult[0].distance + NNCyclesResult[1].distance;
            bestNNValue.second = i;
        }

        if (greedyCyclesResult[0].distance + greedyCyclesResult[1].distance < bestGreedyValue.first)
        {
            bestGreedyValue.first = greedyCyclesResult[0].distance + greedyCyclesResult[1].distance;
            bestGreedyValue.second = i;
        }

        if (regretCyclesResult[0].distance + regretCyclesResult[1].distance < bestRegretValue.first)
        {
            bestRegretValue.first = regretCyclesResult[0].distance + regretCyclesResult[1].distance;
            bestRegretValue.second = i;
        }

        if (NNCyclesResult[0].distance + NNCyclesResult[1].distance > worstNNValue.first)
        {
            worstNNValue.first = NNCyclesResult[0].distance + NNCyclesResult[1].distance;
            worstNNValue.second = i;
        }

        if (greedyCyclesResult[0].distance + greedyCyclesResult[1].distance > worstGreedyValue.first)
        {
            worstGreedyValue.first = greedyCyclesResult[0].distance + greedyCyclesResult[1].distance;
            worstGreedyValue.second = i;
        }

        if (regretCyclesResult[0].distance + regretCyclesResult[1].distance > worstRegretValue.first)
        {
            worstRegretValue.first = regretCyclesResult[0].distance + regretCyclesResult[1].distance;
            worstRegretValue.second = i;
        }

        averageNNValue += NNCyclesResult[0].distance + NNCyclesResult[1].distance;
        averageGreedyValue += greedyCyclesResult[0].distance + greedyCyclesResult[1].distance;
        averageRegretValue += regretCyclesResult[0].distance + regretCyclesResult[1].distance;
    }

    averageNNValue /= verticesCoords.size();
    averageGreedyValue /= verticesCoords.size();
    averageRegretValue /= verticesCoords.size();

    vector<Graph> nNCyclesResult = NNCycles(distanceMatrix, bestNNValue.second);
    vector<Graph> greedyCyclesResult = greedyCycles(distanceMatrix, bestGreedyValue.second);
    vector<Graph> regretCyclesResult = regretCycles(distanceMatrix, bestRegretValue.second);

    cout << averageNNValue << " (" << bestNNValue.first << " – " << worstNNValue.first << ")" << endl;
    cout << averageGreedyValue << " (" << bestGreedyValue.first << " – " << worstGreedyValue.first << ")" << endl;
    cout << averageRegretValue << " (" << bestRegretValue.first << " – " << worstRegretValue.first << ")" << endl;

    saveGraphs(nNCyclesResult, "nNCyclesa.txt");
    saveGraphs(greedyCyclesResult, "greedyCyclesa.txt");
    saveGraphs(regretCyclesResult, "regretCyclesa.txt");

    return 0;
}