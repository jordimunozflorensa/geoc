/**
 TODO Replace this by your own, correct, triangulation function
 Triangles should be return as arrays of array of indexes
 e.g., [[1,2,3],[2,3,4]] encodes two triangles, where the indices are relative to the array points
**/

// Helper function to check if a point is inside a triangle
function pointInTriangle(px, py, ax, ay, bx, by, cx, cy) {
    let v0x = cx - ax, v0y = cy - ay;
    let v1x = bx - ax, v1y = by - ay;
    let v2x = px - ax, v2y = py - ay;

    let dot00 = v0x * v0x + v0y * v0y;
    let dot01 = v0x * v1x + v0y * v1y;
    let dot02 = v0x * v2x + v0y * v2y;
    let dot11 = v1x * v1x + v1y * v1y;
    let dot12 = v1x * v2x + v1y * v2y;

    let invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    let u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    let v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return u >= 0 && v >= 0 && (u + v) < 1;
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
        return pointInTriangle(
            point.x, point.y,
            points[p1].x, points[p1].y,
            points[p2].x, points[p2].y,
            points[p3].x, points[p3].y
        );
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

        // Traverse the tree to find the triangle containing the point
        let found = false;
        while (!found) {
            if (current.children.length === 0) {
                // If no children, check if the point is in this triangle
                if (current.containsPoint(point, points)) {
                    current.subdivide(i, points); // Subdivide if point is inside
                    found = true;
                }
            } else {
                // Otherwise, continue traversing through children
                for (let child of current.children) {
                    if (child.containsPoint(point, points)) {
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



