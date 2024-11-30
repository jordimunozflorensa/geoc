class DCEL {
    constructor() {
        this.tableOfVertices = [];
        this.tableOfFaces = [];
        this.tableOfEdges = [];

        this.numberOfVertices = 0;
        this.numberOfFaces = 0;
        this.numberOfEdges = 0;

        this.Tree = {};
        // aux var to safe the index of a face to the leaf for updating the tree
        this.faceIndex = [];
    }

    det(p, q, r) {
        return (q.x - p.x) * (r.y - p.y) - (r.x - p.x) * (q.y - p.y);
    }

    mod(x, n) {
        return ((x % n) + n) % n;
    }

    orientationTestSegmentPoint(p, q, r) {
        const determinant = this.det(p, q, r);
        return determinant === 0 ? 0 : (determinant > 0 ? 1 : -1);
    }
    
    isInsideTriangle(D1, D2, D3){
      // all determinants have the same sign
      if (D1 == D2 && D2 == D3)
        return true;
      else
        return false;
    }
    
    isOnAVertice(D1, D2, D3){
      // the same as saying that both are 0 and the other is 1
      if (Math.abs(D1) + Math.abs(D2) + Math.abs(D3) == 1)
        return true;
      else
        return false;
    }
    
    isOnAnEdge(D1, D2, D3){
      // has to be 0 with one and the sign in the other two, since if we just check that is 0
      // in some one it could be the case in which it is collinear but not ON the edge
      if (( D1 == 0 && D2 == D3) || 
          ( D2 == 0 && D1 == D3) || 
          ( D3 == 0 && D1 == D2))
        return true;
      else
        return false;
    }

    initDCEL(p0, p1, p2) {
        this.tableOfFaces[0] = {edge: 0};
        this.numberOfFaces += 1;

        this.tableOfVertices[0] = {x: p0.x, y: p0.y, edge: 0};
        this.tableOfVertices[1] = {x: p1.x, y: p1.y, edge: 2};
        this.tableOfVertices[2] = {x: p2.x, y: p2.y, edge: 4};
        this.numberOfVertices += 3;

        // edges1
        this.tableOfEdges[0] = {sourceVertex: 0, rightFace: 0, nextEdge: 2, reverseEdge: 1};
        this.tableOfEdges[1] = {sourceVertex: 1, rightFace: -1, nextEdge: 5, reverseEdge: 0};

        // edges2
        this.tableOfEdges[2] = {sourceVertex: 1, rightFace: 0, nextEdge: 4, reverseEdge: 3};
        this.tableOfEdges[3] = {sourceVertex: 2, rightFace: -1, nextEdge: 1, reverseEdge: 2};

        // edges3
        this.tableOfEdges[4] = {sourceVertex: 2, rightFace: 0, nextEdge: 0, reverseEdge: 5};
        this.tableOfEdges[5] = {sourceVertex: 0, rightFace: -1, nextEdge: 3, reverseEdge: 4};

        this.numberOfEdges += 6;

        this.Tree = {edge: 0, childTriangles: []};
        this.faceIndex = [];
        this.faceIndex[0] = this.Tree;
    }

    findFace(point, triangle) {
        var childTriangles = triangle.childTriangles;

        // base case of the recursion (we have no more children)
        if (childTriangles.length == 0) {
            var faceEdges = this.getFaceFromEdges(this.tableOfEdges[triangle.edge].rightFace);
            var D = new Array(3);
            for (var i = 0; i < faceEdges.length; i++) {
                var sourceEdge = this.tableOfEdges[faceEdges[i]].sourceVertex;
                var destEdge = this.tableOfEdges[this.tableOfEdges[faceEdges[i]].reverseEdge].sourceVertex;

                D[i] = this.orientationTestSegmentPoint(this.tableOfVertices[sourceEdge], this.tableOfVertices[destEdge], point);
            }

            if (this.isOnAnEdge(D[0], D[1], D[2])) {
                for (var i = 0; i < faceEdges.length; i++) {
                    if (D[i] == 0) {
                        return {
                            edge: faceEdges[i],
                            isOnAnEdge: true,
                        };
                    }
                }
            } 
            else if (this.isOnAVertice(D[0], D[1], D[2])) {
                return null;
            }
            else if (this.isInsideTriangle(D[0], D[1], D[2])) {
                return {
                    face: this.tableOfEdges[triangle.edge].rightFace,
                    isOnAnEdge: false,
                };
            }
            else {
                return null;
            }
        }

        // case of
        //              O
        //             /|\
        //            / | \
        //           /  |  \
        //          /   |   \
        //         O----O----O
        if (childTriangles.length == 2) {
            var t0 = childTriangles[0];
            var t1 = childTriangles[1];

            var srcVertex = this.tableOfEdges[t1.edge].sourceVertex;
            var dstVertex = this.tableOfEdges[this.tableOfEdges[t1.edge].reverseEdge].sourceVertex;

            var D1 = this.det(this.tableOfVertices[srcVertex], this.tableOfVertices[dstVertex], point);

            if (D1 <= 0) {
                return this.findFace(point, t0);
            } 
            else {
                return this.findFace(point, t1);
            }
        } 

        // case of 
        //              O
        //             /|\
        //            / | \
        //           /  O  \
        //          /  / \  \
        //         |  /   \  |
        //         | /     \ |      
        //         O---------O
        else if (childTriangles.length == 3) {
            for (var i = 0; i < childTriangles.length; i++) {
                var currentTriangle = childTriangles[i];
                var srcVertex0 = this.tableOfEdges[currentTriangle.edge].sourceVertex;
                var dstVertex0 = this.tableOfEdges[this.tableOfEdges[currentTriangle.edge].reverseEdge].sourceVertex;

                var neighbourTriangle = childTriangles[(i + 1) % childTriangles.length];
                var srcVertex1 = this.tableOfEdges[neighbourTriangle.edge].sourceVertex;
                var dstVertex1 = this.tableOfEdges[this.tableOfEdges[neighbourTriangle.edge].reverseEdge].sourceVertex;

                var D1 = this.det(this.tableOfVertices[srcVertex0], this.tableOfVertices[dstVertex0], point);
                var D2 = this.det(this.tableOfVertices[srcVertex1], this.tableOfVertices[dstVertex1], point);

                if (D1 <= 0 && D2 > 0) {
                    return this.findFace(point, currentTriangle);
                }
            }
        }

        // point is duplicated, we ignore it
        return null;
    }

    getFaceFromEdges(faceIndex) {
        var edges = new Array(3);

        edges[0] = this.tableOfFaces[faceIndex].edge;
        edges[1] = this.tableOfEdges[edges[0]].nextEdge;
        edges[2] = this.tableOfEdges[edges[1]].nextEdge;

        return edges;
    }

    computeNewFaces(faceEdges, newVertex, newFaces, newEdges, newReversedEdges) {
        var N = newFaces.length;
        for (var i = 0; i < newFaces.length; i++) {
            this.tableOfEdges[newEdges[i]] = {
                sourceVertex: newVertex, 
                rightFace: newFaces[i], 
                nextEdge: faceEdges[i], 
                reverseEdge: newReversedEdges[i]
            };
            this.tableOfEdges[newReversedEdges[i]] = {
                sourceVertex: this.tableOfEdges[faceEdges[i]].sourceVertex,
                rightFace: newFaces[this.mod(i - 1, N)],
                nextEdge: newEdges[this.mod(i - 1, N)],
                reverseEdge: newEdges[i]
            };
            this.tableOfFaces[newFaces[i]] = {edge: newEdges[i]};

            this.tableOfEdges[faceEdges[i]].rightFace = newFaces[i];
            this.tableOfEdges[faceEdges[i]].nextEdge = newReversedEdges[this.mod(i + 1, N)];
        }
    }

    updateChagesToTree(isOnAnEdge, prevFaces, newFaces, newEdges) {
        if (isOnAnEdge == false) {
            var treeNode = this.faceIndex[prevFaces[0]];
            treeNode.childTriangles = [
                {edge: newEdges[0], childTriangles: []},
                {edge: newEdges[1], childTriangles: []},
                {edge: newEdges[2], childTriangles: []}
            ];
            this.faceIndex[newFaces[0]] = treeNode.childTriangles[0];
            this.faceIndex[newFaces[1]] = treeNode.childTriangles[1];
            this.faceIndex[newFaces[2]] = treeNode.childTriangles[2];
        } 
        
        else {
            var treeNode0 = this.faceIndex[prevFaces[0]];
            var treeNode1 = this.faceIndex[prevFaces[1]];

            treeNode0.childTriangles = [
                {edge: newEdges[0], childTriangles: []},
                {edge: newEdges[1], childTriangles: []}
            ];
            this.faceIndex[newFaces[0]] = treeNode0.childTriangles[0];
            this.faceIndex[newFaces[1]] = treeNode0.childTriangles[1];

            treeNode1.childTriangles = [
                {edge: newEdges[2], childTriangles: []},
                {edge: newEdges[3], childTriangles: []}
            ];
            this.faceIndex[newFaces[2]] = treeNode1.childTriangles[0];
            this.faceIndex[newFaces[3]] = treeNode1.childTriangles[1];
        }
    }

    addPointDCEL(point) {
        var newVertex = this.numberOfVertices;
        this.tableOfVertices[newVertex] = {x: point.x, y: point.y};
        this.numberOfVertices += 1;

        var searchData = this.findFace(point, this.Tree);
        if (searchData == null) return; // point is duplicate, we ignore it

        if (searchData.isOnAnEdge == false) {
            var face = searchData.face;
            var faceEdges = this.getFaceFromEdges(face);

            // from:                      to:
            //              O                       O
            //             / \                     /|\
            //            /   \    =========>     / | \  
            //           /  O  \                 /  O  \
            //          /       \               /  / \  \
            //         |         |             |  /   \  |
            //         |         |             | /     \ |
            //         O---------O             O---------O
            var newFaces = [
                face,
                this.numberOfFaces,
                this.numberOfFaces + 1
            ];
            var newEdges = [
                this.numberOfEdges,
                this.numberOfEdges + 2,
                this.numberOfEdges + 4
            ];
            var newReversedEdges = [
                this.numberOfEdges + 1,
                this.numberOfEdges + 3,
                this.numberOfEdges + 5
            ];

            this.numberOfFaces += 2;
            this.numberOfEdges += 6;

            this.computeNewFaces(faceEdges, newVertex, newFaces, newEdges, newReversedEdges);
            this.updateChagesToTree(false, [face], newFaces, newEdges);
        } 
        
        else {
            var edge = searchData.edge;
            var reverseEdge = this.tableOfEdges[edge].reverseEdge;

            // instead of one face we have two faces
            var prevFace0 = this.tableOfEdges[edge].rightFace;
            var prevFace1 = this.tableOfEdges[reverseEdge].rightFace;

            var faceEdges = new Array(4);
            faceEdges[0] = this.tableOfEdges[edge].nextEdge;
            faceEdges[1] = this.tableOfEdges[faceEdges[0]].nextEdge;
            faceEdges[2] = this.tableOfEdges[reverseEdge].nextEdge;
            faceEdges[3] = this.tableOfEdges[faceEdges[2]].nextEdge;

            // from:                      to:
            //              O                       O
            //             / \                     /|\
            //            /   \    =========>     / | \  
            //           /     \                 /  |  \
            //          /       \               /   |   \
            //         O----O----O             O----O----O
            var newFaces = [
                this.tableOfEdges[edge].rightFace,
                this.tableOfEdges[reverseEdge].rightFace,
                this.numberOfFaces, this.numberOfFaces + 1
            ];
            var newEdges = [
                edge,
                this.numberOfEdges,
                this.numberOfEdges + 2,
                this.numberOfEdges + 4
            ];
            var newReversedEdges = [
                reverseEdge,
                this.numberOfEdges + 1,
                this.numberOfEdges + 3,
                this.numberOfEdges + 5
            ];

            this.numberOfFaces += 2;
            this.numberOfEdges += 8;

            this.computeNewFaces(faceEdges, newVertex, newFaces, newEdges, newReversedEdges);
            this.updateChagesToTree(true, [prevFace0, prevFace1], newFaces, newEdges);
        }
    }
}

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

function computeContainingTriangle(points, dcel) {
    var boundingBox = getBoundingBox(points);

    var offsetSmallCases = 0.2;
    // var offsetBigCases = 0.05;
    var offset = offsetSmallCases;

    var offsetX = (boundingBox.xmax - boundingBox.xmin) * offset;
    var offsetY = (boundingBox.ymax - boundingBox.ymin) * offset;

    var p0 = {
        x: boundingBox.xmin - (boundingBox.xmax - boundingBox.xmin) / 2 - offsetX,
        y: boundingBox.ymin - offsetY
    };
    var p1 = {
        x: (boundingBox.xmin + boundingBox.xmax) / 2,
        y: boundingBox.ymax + (boundingBox.ymax - boundingBox.ymin)
    };
    var p2 = {
        x: boundingBox.xmax + (boundingBox.xmax - boundingBox.xmin) / 2 + offsetX,
        y: boundingBox.ymin - offsetY
    };

    points.push(p0, p1, p2);
}

function getFaceVertices(dcel, faceIndex) {
    var vertices = new Array(3);
    var edges = dcel.getFaceFromEdges(faceIndex);

    vertices[0] = dcel.tableOfEdges[edges[0]].sourceVertex;
    vertices[1] = dcel.tableOfEdges[edges[1]].sourceVertex;
    vertices[2] = dcel.tableOfEdges[edges[2]].sourceVertex;

    return vertices;
}

function shufflePoints(points) {
    for (let i = points.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [points[i], points[j]] = [points[j], points[i]];
    }
}

function computeTriangulation(points) {
    var dcel = new DCEL();

    // shufflePoints(points);

    // compute the Containing Triangle
    computeContainingTriangle(points, dcel);
    var N = points.length;

    // init DCEL
    dcel.initDCEL(points[N - 3], points[N - 2], points[N - 1]);

    // compute the triangulation
    for (var i = 0; i < N - 3; i++) {
        dcel.addPointDCEL(points[i]);
    }

    // get all the triangles of the DCEL
    var outputTriangles = new Array(dcel.tableOfFaces.length);

    for (var i = 0; i < outputTriangles.length; i++) {
        var faceVertices = getFaceVertices(dcel, i);
        // we take the indices of the vertices of the face
        outputTriangles[i] = [
            (faceVertices[0] + N - 3) % N,
            (faceVertices[1] + N - 3) % N,
            (faceVertices[2] + N - 3) % N
        ];
    }
    return outputTriangles;
}