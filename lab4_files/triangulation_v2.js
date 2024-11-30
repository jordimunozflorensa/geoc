function orientationTestSegmentPoint(p, q, r){
    const determinant = (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
    const orientation = Math.sign(determinant);

    return orientation;
}

function isInsideTriangle(D1, D2, D3){
  if (D1 == D2 && D2 == D3)
    return true;
  else
    return false;
}

function isOnAnEdge(D1, D2, D3){
    if (( D1 == 0 && D2 == D3) || 
        ( D2 == 0 && D1 == D3) || 
        ( D3 == 0 && D1 == D2))
      return true;
    else
      return false;
  }

// Triangle node class for the tree structure
class TriangleNode {
    constructor(p1, p2, p3) {
        this.points = [p1, p2, p3]; // store point indices
        this.children = []; // store child triangles
    }

    // Subdivide the triangle and add child nodes
    subdivide(newPointIndex, points) {
        const [p1, p2, p3] = this.points;
        this.children.push(new TriangleNode(p1, p2, newPointIndex));
        this.children.push(new TriangleNode(p2, p3, newPointIndex));
        this.children.push(new TriangleNode(p3, p1, newPointIndex));
    }

    // Check if a point is inside this triangle
    containsPoint(point, points) {
        const [p1, p2, p3] = this.points;
        var D1 = orientationTestSegmentPoint(points[p1], points[p2], point);
        var D2 = orientationTestSegmentPoint(points[p2], points[p3], point);
        var D3 = orientationTestSegmentPoint(points[p3], points[p1], point);
        return isInsideTriangle(D1, D2, D3);
    }

    // Check if a point is on the edge of this triangle
    isOnEdge(point, points) {
        const [p1, p2, p3] = this.points;
        var D1 = orientationTestSegmentPoint(points[p1], points[p2], point);
        var D2 = orientationTestSegmentPoint(points[p2], points[p3], point);
        var D3 = orientationTestSegmentPoint(points[p3], points[p1], point);
        return isOnAnEdge(D1, D2, D3);
    }
}

// Main function to compute triangulation
function computeTriangulation(points) {
    const outputTriangles = [[0, 1, 2]]; // initial triangle (using points 0, 1, 2)

    // Initialize tree with the root triangle
    const root = new TriangleNode(0, 1, 2);

    // Process each additional point for triangulation
    for (let i = 3; i < points.length; i++) {
        let current = root;
        const point = points[i];
        
        console.log("point: ", point, "i: ", i);
        console.log("current: ", current);
        // Traverse the tree to find the triangle containing the point
        let found = false;
        while (!found) {
            if (current.children.length === 0) {
                // If no children, check if the point is in this triangle or on its edge
                if (current.containsPoint(point, points)) {
                    current.subdivide(i, points); // Subdivide if point is inside
                    found = true;
                } else if (current.isOnEdge(point, points)) {
                    // Handle the case where the point is on an edge
                    // const [p1, p2, p3] = current.points;
                    // const edges = [
                    //     [p1, p2],
                    //     [p2, p3],
                    //     [p3, p1]
                    // ];

                    // for (const [e1, e2] of edges) {
                    //     if (orientationTestSegmentPoint(points[e1], points[e2], point) === 0) {
                    //         // Subdivide the triangle into two new triangles
                    //         current.children.push(new TriangleNode(e1, e2, i));
                    //         current.children.push(new TriangleNode(e2, p3, i));
                    //         current.children.push(new TriangleNode(p3, e1, i));
                    //         found = true;
                    //         break;
                    //     }
                    // }
                    found = true;
                    break;
                }
            } else {
                // Otherwise, continue traversing through children
                for (let child of current.children) {
                    if (child.containsPoint(point, points) || child.isOnEdge(point, points)) {
                        current = child;
                        break;
                    }
                }
            }
        }
    }

    // Helper function to extract triangles from the tree
    function extractTriangles(node) {
        if (node.children.length === 0) {
            outputTriangles.push(node.points);
        } else {
            for (const child of node.children) {
                extractTriangles(child);
            }
        }
    }

    extractTriangles(root);
    return outputTriangles;
}



