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

function transformToDelaunayDCEL(newVertex, dcel) {
    const startingEdge = dcel.tableOfVertices[newVertex].edge;
    let currentEdge = startingEdge;
    let hasBeenRotation = false;

    do {
        hasBeenRotation = false;
        const oppositeEdge = dcel.tableOfEdges[currentEdge].nextEdge;
        const adjacentFace = dcel.tableOfEdges[dcel.tableOfEdges[oppositeEdge].twinEdge].rightFace;

        // if no infinite face and triangle contains newVertex, then rotate edge
        if (adjacentFace != -1 && isInsideTriangleCircle(newVertex, adjacentFace, dcel)) {
            rotateEdge(oppositeEdge, dcel);
            hasBeenRotation = true;
        } else {
            const lastEdge = dcel.tableOfEdges[oppositeEdge].nextEdge;
            currentEdge = dcel.tableOfEdges[lastEdge].twinEdge;
        }
    } while (currentEdge != startingEdge || hasBeenRotation);
}

function isInsideTriangleCircle(vertex, face, dcel) {
    const faceVertices = getFaceVertices(face, dcel);
    return isPointInCircle(dcel.tableOfVertices[vertex], [dcel.tableOfVertices[faceVertices[0]], dcel.tableOfVertices[faceVertices[1]], dcel.tableOfVertices[faceVertices[2]]], dcel);
}

function rotateEdge(edge, dcel) {
    // identify data requiring update
    const twinEdge = dcel.tableOfEdges[edge].twinEdge;

    const e0 = dcel.tableOfEdges[edge].nextEdge;
    const e1 = dcel.tableOfEdges[e0].nextEdge;
    const e2 = dcel.tableOfEdges[twinEdge].nextEdge;
    const e3 = dcel.tableOfEdges[e2].nextEdge;

    const f0 = dcel.tableOfEdges[e0].rightFace;
    const f1 = dcel.tableOfEdges[e2].rightFace;

    const v0 = dcel.tableOfEdges[e0].sourceVertex;
    const v1 = dcel.tableOfEdges[e2].sourceVertex;

    // rotate edge counterclockwise
    const edgeData = dcel.tableOfEdges[edge];
    const twinEdgeData = dcel.tableOfEdges[twinEdge];

    edgeData.sourceVertex = dcel.tableOfEdges[e1].sourceVertex;
    edgeData.nextEdge = e3;

    twinEdgeData.sourceVertex = dcel.tableOfEdges[e3].sourceVertex;
    twinEdgeData.nextEdge = e1;

    // update affected edges
    dcel.tableOfEdges[e0].nextEdge = edge;

    dcel.tableOfEdges[e1].nextEdge = e2;
    dcel.tableOfEdges[e1].rightFace = f1;

    dcel.tableOfEdges[e2].nextEdge = twinEdge;

    dcel.tableOfEdges[e3].nextEdge = e0;
    dcel.tableOfEdges[e3].rightFace = f0;

    // update affected faces
    dcel.tableOfFaces[f0].edge = edge;
    dcel.tableOfFaces[f1].edge = twinEdge;

    // update affected vertices
    dcel.tableOfVertices[v0].edge = e0;
    dcel.tableOfVertices[v1].edge = e2;
}

function getVertexEdges(dcel, vertex) {
    const edges = [];

    const startingEdge = dcel.tableOfVertices[vertex].edge;
    let currentEdge = startingEdge;
    do {
        edges.push(currentEdge);
        const triangleEdge2 = dcel.tableOfEdges[currentEdge].nextEdge;
        const triangleEdge3 = dcel.tableOfEdges[triangleEdge2].nextEdge;

        currentEdge = dcel.tableOfEdges[triangleEdge3].twinEdge;
    } while (currentEdge != startingEdge);

    return edges;
}

function removeEdge(dcel, edge) {
    const edgeTwin = dcel.tableOfEdges[edge].twinEdge;
    const faceRight = dcel.tableOfEdges[edge].rightFace;
    const faceLeft = dcel.tableOfEdges[edgeTwin].rightFace;

    dcel.tableOfFaces[faceRight] = null;
    dcel.tableOfFaces[faceLeft] = null;

    // TODO: update other affected DCEL data 
    // if it needs to be used in the future after deletions
}

function pruneBoundaries(dcel, boundaries) {
    boundaries.forEach(boundary => {
        if (boundary.length == 1) {
            const vertex = boundary[0] + 3; // we need to consider the enclosing triangle points added first
            const edges = getVertexEdges(dcel, vertex);
            edges.forEach(edge => {
                removeEdge(dcel, edge);
            });
        } else {
            for (let i = 0; i < boundary.length; i++) {
                const vertex = boundary[i] + 3;
                const previousVertex = boundary[mod(i - 1, boundary.length)] + 3;
                const nextVertex = boundary[mod(i + 1, boundary.length)] + 3;

                const edges = getVertexEdges(dcel, vertex);
                edges.forEach(edge => {
                    const dstVertex = dcel.tableOfEdges[dcel.tableOfEdges[edge].twinEdge].sourceVertex;

                    // we need to check if local boundary is convex or concave
                    const turn = det(dcel.tableOfVertices[previousVertex], dcel.tableOfVertices[vertex], dcel.tableOfVertices[nextVertex]);

                    const det1 = det(dcel.tableOfVertices[previousVertex], dcel.tableOfVertices[vertex], dcel.tableOfVertices[dstVertex]);
                    const det2 = det(dcel.tableOfVertices[vertex], dcel.tableOfVertices[nextVertex], dcel.tableOfVertices[dstVertex]);

                    if ((turn < 0 && det1 < 0 && det2 < 0) || (turn >= 0 && (det1 < 0 || det2 < 0))) {
                        removeEdge(dcel, edge);
                    }
                });
            }
        }
    });
}

function addEnclosingTriangle(dcel, points) {
    if (dcel.enclosingTriangleAdded) {
        points.pop();
        points.pop();
        points.pop();
    }

    const box = getBoundingBox(points);

    const xOffset = (box.xmax - box.xmin) * 0.1;
    const yOffset = (box.ymax - box.ymin) * 0.1;

    const p0 = {
        x: (box.xmin + box.xmax) / 2,
        y: box.ymax + (box.ymax - box.ymin)
    };
    const p1 = {
        x: box.xmax + (box.xmax - box.xmin) / 2 + xOffset,
        y: box.ymin - yOffset
    };
    const p2 = {
        x: box.xmin - (box.xmax - box.xmin) / 2 - xOffset,
        y: box.ymin - yOffset
    };

    points.push(p0, p1, p2);
    dcel.enclosingTriangleAdded = true;
}

function getBoundingBox(points) {
    let xmin, xmax, ymin, ymax;
    xmin = xmax = points[0].x;
    ymin = ymax = points[0].y;

    points.forEach(point => {
        if (point.x < xmin) xmin = point.x;
        else if (point.x > xmax) xmax = point.x;

        if (point.y < ymin) ymin = point.y;
        else if (point.y > ymax) ymax = point.y;
    });

    return { xmin: xmin, xmax: xmax, ymin: ymin, ymax: ymax };
}

function isPointInTriangle(p, triangle, dcel) {
    // determinant signs of the triangles segments
    const d_s1 = Math.sign(det(triangle[0], triangle[1], p));
    const d_s2 = Math.sign(det(triangle[1], triangle[2], p));
    const d_s3 = Math.sign(det(triangle[2], triangle[0], p));

    // we check where the point lies for every segment (-1:'sideA, 0:'onLine', 1:sideB) to study its case
    const insideTriangle = (d_s1 == d_s2 && d_s1 == d_s3);
    const onEdge = (d_s1 == 0 && d_s2 == d_s3) || (d_s2 == 0 && d_s1 == d_s3) || (d_s3 == 0 && d_s1 == d_s2);
    const onCorner = Math.abs(d_s1) + Math.abs(d_s2) + Math.abs(d_s3) == 1;

    // TODO: change output when necessary
    return (insideTriangle || onEdge || onCorner);
}

function isPointInCircle(p, circle_points, dcel) {
    const a = circle_points[0];
    const b = circle_points[1];
    const c = circle_points[2];

    const clockOrder = Math.sign(det(a, b, c));
    const planeSide = Math.sign(det_projected_points(p, a, b, c, dcel));

    return planeSide == 0 || planeSide != clockOrder;
}

function mod(a, n) {
    return ((a % n) + n) % n;
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

function det(p, q, r) {
    return (q.x - p.x) * (r.y - p.y) - (r.x - p.x) * (q.y - p.y);
}

function matrixDeterminant(M) {
    if (M.length == 2) {
        return M[0][0] * M[1][1] - M[0][1] * M[1][0];
    } else {
        return M[0][0] * matrixDeterminant(submatrix(M, 0, 0)) - M[0][1] * matrixDeterminant(submatrix(M, 0, 1)) + M[0][2] * matrixDeterminant(submatrix(M, 0, 2));
    }
}

function det_projected_points(p, a, b, c, dcel) {
    const M = [
        [b.x - a.x, b.y - a.y, (b.x - a.x) * (b.x + a.x) + (b.y - a.y) * (b.y + a.y)],
        [c.x - a.x, c.y - a.y, (c.x - a.x) * (c.x + a.x) + (c.y - a.y) * (c.y + a.y)],
        [p.x - a.x, p.y - a.y, (p.x - a.x) * (p.x + a.x) + (p.y - a.y) * (p.y + a.y)]
    ];
    return matrixDeterminant(M);
}

function computeTriangulation(points, boundaries) {
    const dcel = new DCEL();
    addEnclosingTriangle(dcel, points);
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
        } else {
            console.log("Face ignored");
        }
    }
    return outputTriangles;
}
