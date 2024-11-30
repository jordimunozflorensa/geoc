// BOTTOM-UP CODE STRUCTURE (to ease readability)

// ----------------------------------- General Functions -----------------------------------
function mod(a, n) {
    return ((a % n) + n) % n;
}

function det(p, q, r) {
    return (q.x - p.x) * (r.y - p.y) - (r.x - p.x) * (q.y - p.y);
}

function submatrix(M, row, column) {
    const ret = new Array(M.length);
    for (let i = 0; i < ret.length; i++) {
        ret[i] = [...M[i]];
    }

    ret.splice(row, 1);
    for (let i = 0; i < ret.length; i++) {
        ret[i].splice(column, 1);
    }

    return ret;
}


function matrixDeterminant(M) {
    if (M.length == 2) {
        return M[0][0] * M[1][1] - M[0][1] * M[1][0];
    } else {
        return M[0][0] * matrixDeterminant(submatrix(M, 0, 0)) - M[0][1] * matrixDeterminant(submatrix(M, 0, 1)) + M[0][2] * matrixDeterminant(submatrix(M, 0, 2));
    }
}

function detProjectedPoints(p, a, b, c, dcel) {
    const M = [
        [b.x - a.x, b.y - a.y, (b.x - a.x) * (b.x + a.x) + (b.y - a.y) * (b.y + a.y)],
        [c.x - a.x, c.y - a.y, (c.x - a.x) * (c.x + a.x) + (c.y - a.y) * (c.y + a.y)],
        [p.x - a.x, p.y - a.y, (p.x - a.x) * (p.x + a.x) + (p.y - a.y) * (p.y + a.y)]
    ];
    return matrixDeterminant(M);
}

// ----------------------------------- DCEL Structure -----------------------------------

class DCEL {
    constructor() {
        this.tableOfVertices = [];
        this.tableOfFaces = [];
        this.tableOfEdges = [];

        this.numberOfVertices = 0;
        this.numberOfFaces = 0;
        this.numberOfEdges = 0;

        this.enclosingTriangleAdded = false;
    }

    initDCEL(p0, p1, p2) {
        // reset just in case we execute the triangulation multiple times
        this.tableOfVertices = [];
        this.tableOfFaces = [];
        this.tableOfEdges = [];
        this.numberOfVertices = 0;
        this.numberOfFaces = 0;
        this.numberOfEdges = 0;

        // faces
        this.tableOfFaces[0] = { edge: 0 };
        this.numberOfFaces += 1;
        
        // vertices
        this.tableOfVertices[0] = { x: p0.x, y: p0.y, edge: 0 };
        this.tableOfVertices[1] = { x: p1.x, y: p1.y, edge: 2 };
        this.tableOfVertices[2] = { x: p2.x, y: p2.y, edge: 4 };
        this.numberOfVertices += 3;

        // edges1
        this.tableOfEdges[0] = {sourceVertex: 0, rightFace: 0, nextEdge: 2, twinEdge: 1};
        this.tableOfEdges[1] = {sourceVertex: 1, rightFace: -1, nextEdge: 5, twinEdge: 0};

        // edges2
        this.tableOfEdges[2] = {sourceVertex: 1, rightFace: 0, nextEdge: 4, twinEdge: 3};
        this.tableOfEdges[3] = {sourceVertex: 2, rightFace: -1, nextEdge: 1, twinEdge: 2};

        // edges3
        this.tableOfEdges[4] = {sourceVertex: 2, rightFace: 0, nextEdge: 0, twinEdge: 5};
        this.tableOfEdges[5] = {sourceVertex: 0, rightFace: -1, nextEdge: 3, twinEdge: 4};

        this.numberOfEdges += 6;
    }

    addPointDCEL(point) {
        const newVertex = this.numberOfVertices;
        this.tableOfVertices[newVertex] = { x: point.x, y: point.y };
        this.numberOfVertices += 1;

        // substitution of tree structure to find the face where the point stands
        const face = findFaceFromEdge(this.tableOfVertices[2].edge, point, this);

        const faceEdges = getEdgesFromFace(face, this);
        const newFaces = [face, this.numberOfFaces, this.numberOfFaces + 1];
        const newEdges = [this.numberOfEdges, this.numberOfEdges + 2, this.numberOfEdges + 4];
        const newTwinEdges = [this.numberOfEdges + 1, this.numberOfEdges + 3, this.numberOfEdges + 5];

        this.numberOfEdges += 6;
        this.numberOfFaces += 2;

        this.computeNewFaces(faceEdges, newVertex, newFaces, newEdges, newTwinEdges);
        
        // function to implement the loop for rotations
        transformToDelaunayDCEL(newVertex, this);
    }

    computeNewFaces(faceEdges, newVertex, newFaces, newEdges, newTwinEdges) {
        const N = newFaces.length;

        for (let i = 0; i < newFaces.length; i++) {
            this.tableOfVertices[newVertex].edge = newEdges[i];

            this.tableOfEdges[newEdges[i]] = { sourceVertex: newVertex, rightFace: newFaces[i], nextEdge: faceEdges[i], twinEdge: newTwinEdges[i] };
            this.tableOfEdges[newTwinEdges[i]] = { sourceVertex: this.tableOfEdges[faceEdges[i]].sourceVertex, rightFace: newFaces[mod(i - 1, N)], nextEdge: newEdges[mod(i - 1, N)], twinEdge: newEdges[i] };

            this.tableOfFaces[newFaces[i]] = { edge: newEdges[i] };

            this.tableOfEdges[faceEdges[i]].rightFace = newFaces[i];
            this.tableOfEdges[faceEdges[i]].nextEdge = newTwinEdges[mod(i + 1, N)];
        }
    }
}

// ----------------------------------- Delaunay Triangulation -----------------------------------

// recursive function that replaces the tree structure to finde the face
function findFaceFromEdge(auxEdge, point, dcel) {
    let currentEdge = auxEdge;

    for (let i = 0; i < dcel.numberOfEdges; i++) {
        const srcVertex = dcel.tableOfEdges[currentEdge].sourceVertex;
        const dstVertex = dcel.tableOfEdges[dcel.tableOfEdges[currentEdge].twinEdge].sourceVertex;

        // if point to the left (since we are in the right face)
        if (det(dcel.tableOfVertices[srcVertex], dcel.tableOfVertices[dstVertex], point) > 0) {
            const actualTwin = dcel.tableOfEdges[currentEdge].twinEdge;
            // recursive call
            return findFaceFromEdge(actualTwin, point, dcel);
        }

        currentEdge = dcel.tableOfEdges[currentEdge].nextEdge;
        if (currentEdge == auxEdge) break;
    }

    return dcel.tableOfEdges[currentEdge].rightFace;
}

function getEdgesFromFace(faceIndex, dcel) {
    const edges = new Array(3);

    edges[0] = dcel.tableOfFaces[faceIndex].edge;
    edges[1] = dcel.tableOfEdges[edges[0]].nextEdge;
    edges[2] = dcel.tableOfEdges[edges[1]].nextEdge;

    return edges;
}

function getFaceVertices(faceIndex, dcel) {
    const vertices = new Array(3);
    const edges = getEdgesFromFace(faceIndex, dcel);

    vertices[0] = dcel.tableOfEdges[edges[0]].sourceVertex;
    vertices[1] = dcel.tableOfEdges[edges[1]].sourceVertex;
    vertices[2] = dcel.tableOfEdges[edges[2]].sourceVertex;

    return vertices;
}

function isInsideTriangleCircle(vertex, face, dcel) {
    const faceVertices = getFaceVertices(face, dcel);
    return isPointInCircle(dcel.tableOfVertices[vertex], [dcel.tableOfVertices[faceVertices[0]], 
                                                        dcel.tableOfVertices[faceVertices[1]], 
                                                        dcel.tableOfVertices[faceVertices[2]]], dcel);
}

function isPointInCircle(p, circle_points, dcel) {
    const a = circle_points[0];
    const b = circle_points[1];
    const c = circle_points[2];

    const clockOrder = Math.sign(det(a, b, c));
    const planeSide = Math.sign(detProjectedPoints(p, a, b, c, dcel));

    // we have to check the order too
    return planeSide == 0 || planeSide != clockOrder;
}

function transformToDelaunayDCEL(newVertex, dcel) {
    const startingEdge = dcel.tableOfVertices[newVertex].edge;
    let currentEdge = startingEdge;
    let hasBeenRotation = false;

    for (let i = 0; i < dcel.numberOfEdges; i++) {
        hasBeenRotation = false;
        const oppositeEdge = dcel.tableOfEdges[currentEdge].nextEdge;
        const adjacentFace = dcel.tableOfEdges[dcel.tableOfEdges[oppositeEdge].twinEdge].rightFace;

        // case when we have to rotate the edge (point inside the circumcircle)
        // check also that the face is not the infinite one
        if (adjacentFace != -1 && isInsideTriangleCircle(newVertex, adjacentFace, dcel)) {
            edgeRotation(oppositeEdge, dcel);
            hasBeenRotation = true;
        } else {
            const lastEdge = dcel.tableOfEdges[oppositeEdge].nextEdge;
            currentEdge = dcel.tableOfEdges[lastEdge].twinEdge;
        }

        if (currentEdge == startingEdge && !hasBeenRotation) {
            break;
        }
    }
}

function edgeRotation(edge, dcel) {
    const twinEdge = dcel.tableOfEdges[edge].twinEdge;

    const edge0 = dcel.tableOfEdges[edge].nextEdge;
    const edge1 = dcel.tableOfEdges[edge0].nextEdge;
    const edge2 = dcel.tableOfEdges[twinEdge].nextEdge;
    const edge3 = dcel.tableOfEdges[edge2].nextEdge;

    const vertex0 = dcel.tableOfEdges[edge0].sourceVertex;
    const vertex1 = dcel.tableOfEdges[edge2].sourceVertex;

    const face0 = dcel.tableOfEdges[edge0].rightFace;
    const face1 = dcel.tableOfEdges[edge2].rightFace;

    // edge rotation in counter-clockwise direction
    const edgeAtributes = dcel.tableOfEdges[edge];
    const twinEdgeAtributes = dcel.tableOfEdges[twinEdge];

    edgeAtributes.sourceVertex = dcel.tableOfEdges[edge1].sourceVertex;
    edgeAtributes.nextEdge = edge3;

    twinEdgeAtributes.sourceVertex = dcel.tableOfEdges[edge3].sourceVertex;
    twinEdgeAtributes.nextEdge = edge1;

    // update the structures that are affected by the rotation
    dcel.tableOfEdges[edge0].nextEdge = edge;

    dcel.tableOfEdges[edge1].nextEdge = edge2;
    dcel.tableOfEdges[edge1].rightFace = face1;

    dcel.tableOfEdges[edge2].nextEdge = twinEdge;

    dcel.tableOfEdges[edge3].nextEdge = edge0;
    dcel.tableOfEdges[edge3].rightFace = face0;

    dcel.tableOfVertices[vertex0].edge = edge0;
    dcel.tableOfVertices[vertex1].edge = edge2;

    dcel.tableOfFaces[face0].edge = edge;
    dcel.tableOfFaces[face1].edge = twinEdge;
}

// ----------------------------------- Pruning Boundaries -----------------------------------

function removeEdge(dcel, edge) {
    const edgeTwin = dcel.tableOfEdges[edge].twinEdge;
    const faceRight = dcel.tableOfEdges[edge].rightFace;

    // we can use the twin face since we are in the right face
    const faceLeft = dcel.tableOfEdges[edgeTwin].rightFace;

    dcel.tableOfFaces[faceRight] = null;
    dcel.tableOfFaces[faceLeft] = null;
}

function getVertexEdges(dcel, vertex) {
    const edges = [];
    const startingEdge = dcel.tableOfVertices[vertex].edge;
    let currentEdge = startingEdge;

    for (let i = 0; i < dcel.numberOfEdges; i++) {
        edges.push(currentEdge);
        const triangleEdge2 = dcel.tableOfEdges[currentEdge].nextEdge;
        const triangleEdge3 = dcel.tableOfEdges[triangleEdge2].nextEdge;

        currentEdge = dcel.tableOfEdges[triangleEdge3].twinEdge;
        if (currentEdge == startingEdge) break;
    }

    return edges;
}

function pruneBoundaries(dcel, boundaries) {
    for (let b = 0; b < boundaries.length; b++) {
        const boundary = boundaries[b];

        if (boundary.length == 1) {
            // we have to remove all the edges that are connected to the vertex
            const vertex = boundary[0] + 3;
            const edges = getVertexEdges(dcel, vertex);
            for (let e = 0; e < edges.length; e++) {
                removeEdge(dcel, edges[e]);
            }
        } 
        // we have to remove the edges that are not convex
        else {
            for (let i = 0; i < boundary.length; i++) {
                const vertex = boundary[i] + 3;
                const previousVertex = boundary[mod(i - 1, boundary.length)] + 3;
                const nextVertex = boundary[mod(i + 1, boundary.length)] + 3;

                const edges = getVertexEdges(dcel, vertex);
                for (let e = 0; e < edges.length; e++) {
                    const edge = edges[e];
                    const dstVertex = dcel.tableOfEdges[dcel.tableOfEdges[edge].twinEdge].sourceVertex;

                    // we need to check if local boundary is convex or concave compared to the global boundary
                    const turn = det(dcel.tableOfVertices[previousVertex], dcel.tableOfVertices[vertex], dcel.tableOfVertices[nextVertex]);

                    const det1 = det(dcel.tableOfVertices[previousVertex], dcel.tableOfVertices[vertex], dcel.tableOfVertices[dstVertex]);
                    const det2 = det(dcel.tableOfVertices[vertex], dcel.tableOfVertices[nextVertex], dcel.tableOfVertices[dstVertex]);
                    
                    // if the edge is not convex compared to the boundary, we remove it
                    if ((turn < 0 && det1 < 0 && det2 < 0) || (turn >= 0 && (det1 < 0 || det2 < 0))) {
                        removeEdge(dcel, edge);
                    }
                }
            }
        }
    }
}

// ----------------------------------- Containing Triangle -----------------------------------

function getBoundingBox(points) {
    var xmin, xmax, ymin, ymax;
    xmin = xmax = points[0].x;
    ymin = ymax = points[0].y;

    for (var i = 1; i < points.length; i++) {
        var point = points[i];
        if (point.x < xmin) xmin = point.x;
        else if (point.x > xmax) xmax = point.x;

        if (point.y < ymin) ymin = point.y;
        else if (point.y > ymax) ymax = point.y;
    }

    return {xmin: xmin, xmax: xmax, ymin: ymin, ymax: ymax};
}

function computeContainingTriangle(dcel, points) {
    const boundingBox = getBoundingBox(points);

    const offsetSmallCases = 0.2;
    const offset = offsetSmallCases;

    const offsetX = (boundingBox.xmax - boundingBox.xmin) * offset;
    const offsetY = (boundingBox.ymax - boundingBox.ymin) * offset;

    const p0 = {
        x: boundingBox.xmin - (boundingBox.xmax - boundingBox.xmin) / 2 - offsetX,
        y: boundingBox.ymin - offsetY
    };
    const p1 = {
        x: (boundingBox.xmin + boundingBox.xmax) / 2,
        y: boundingBox.ymax + (boundingBox.ymax - boundingBox.ymin)
    };
    const p2 = {
        x: boundingBox.xmax + (boundingBox.xmax - boundingBox.xmin) / 2 + offsetX,
        y: boundingBox.ymin - offsetY
    };

    points.push(p0, p1, p2);
    dcel.enclosingTriangleAdded = true;
}

// -------------------------------------- Triangulation --------------------------------------

function computeTriangulation(points, boundaries) {
    const dcel = new DCEL();
    computeContainingTriangle(dcel, points);
    const N = points.length;

    dcel.initDCEL(points[N - 3], points[N - 2], points[N - 1]);
    for (let i = 0; i < N - 3; i++) {
        dcel.addPointDCEL(points[i]);
    }

    if (boundaries != null) {
        pruneBoundaries(dcel, boundaries);
    }

    const outputTriangles = [];
    for (let i = 0; i < dcel.numberOfFaces; i++) {
        if (dcel.tableOfFaces[i] != null) {
            const faceVertices = getFaceVertices(i, dcel);
            outputTriangles.push([mod(faceVertices[2] - 3, N), mod(faceVertices[1] - 3, N), mod(faceVertices[0] - 3, N)]); // Store INDICES, not points
        }
    }
    return outputTriangles;
}
