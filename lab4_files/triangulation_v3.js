//#region DCEL 

var vertexTable = [];
var vertexCount = 0;

var faceTable = [];
var faceCount = 0;

var edgeTable = [];
var edgeCount = 0;

//hierarchy of triangles used to find face of a point
var rootHierarchy = {};
var faceToLeaf = []; // map face indices to leaf nodes of the triangle hierarchy

// boolean to check cases
var containingTriangleAdded = false;

function det(p, q, r){
	return (q.x - p.x)*(r.y - p.y) - (r.x - p.x)*(q.y - p.y);
}

function mod(a, n) {
	return ((a % n ) + n ) % n;
}

// REMIND we are working in counterclockwise order
function initializeDCEL(p0, p1, p2){
	vertexTable = []; vertexCount = 0;
	faceTable = []; faceCount = 0;
	edgeTable = []; edgeCount = 0;

	vertexTable[0] = {x: p0.x, y: p0.y, e: 0};
	vertexTable[1] = {x: p1.x, y: p1.y, e: 2};
	vertexTable[2] = {x: p2.x, y: p2.y, e: 4};
	vertexCount += 3;

	edgeTable[0] = {vOrigin: 0, fRight: 0, eN: 2, eTwin: 1};
	edgeTable[1] = {vOrigin: 1, fRight: -1, eN: 5, eTwin: 0};

	edgeTable[2] = {vOrigin: 1, fRight: 0, eN: 4, eTwin: 3};
	edgeTable[3] = {vOrigin: 2, fRight: -1, eN: 1, eTwin: 2};

	edgeTable[4] = {vOrigin: 2, fRight: 0, eN: 0, eTwin: 5};
	edgeTable[5] = {vOrigin: 0, fRight: -1, eN: 3, eTwin: 4};

	edgeCount += 6;

	faceTable[0] = {e: 0}
	faceCount += 1;

	rootHierarchy = {e: 0, subtriangles: []};
	faceToLeaf = [];
	faceToLeaf[0] = rootHierarchy;
}

// Tree version with O(log n)
function searchFaceHierarchy(point, triangle){

	// we are in a leaf
	if (triangle.subtriangles.length == 0) {
		var faceEdges = getFaceEdges(edgeTable[triangle.e].fRight);
		for (var i = 0; i < faceEdges.length; i++){
			var eOrigin = edgeTable[faceEdges[i]].vOrigin;
			var eDest = edgeTable[edgeTable[faceEdges[i]].eTwin].vOrigin;
            
            // point on an edge
			if ( det(vertexTable[eOrigin], vertexTable[eDest], point) == 0) {
				return {
					degenerateCase: true,
					edge: faceEdges[i], 
				};
			}
		}
		return {
			degenerateCase: false,
			face: edgeTable[triangle.e].fRight, 
		}	
	}

	
	var subtriangles = triangle.subtriangles;

	// case of 
	if (subtriangles.length == 2){
		var triangle0 = subtriangles[0];
		var triangle1 = subtriangles[1];

		var e0_origin = edgeTable[triangle1.e].vOrigin;
		var e0_dest = edgeTable[edgeTable[triangle1.e].eTwin].vOrigin;

		if (det(vertexTable[e0_origin], vertexTable[e0_dest], point) >= 0) {
			return searchFaceHierarchy(point, triangle0);
		}
		else{
			return searchFaceHierarchy(point, triangle1);
		}
	}

	else if (subtriangles.length == 3) {
		for (var i = 0; i < subtriangles.length; i++) {

			var currentTriangle = subtriangles[i];
			var e0_origin = edgeTable[currentTriangle.e].vOrigin;
			var e0_dest = edgeTable[edgeTable[currentTriangle.e].eTwin].vOrigin;

			var nextAdjacentTriangle = subtriangles[(i+1)%subtriangles.length];
			var e1_origin = edgeTable[nextAdjacentTriangle.e].vOrigin;
			var e1_dest = edgeTable[edgeTable[nextAdjacentTriangle.e].eTwin].vOrigin;

			var det1 = det(vertexTable[e0_origin], vertexTable[e0_dest], point);
			var det2 = det(vertexTable[e1_origin], vertexTable[e1_dest], point);
			
			if ( det1 <= 0 && det2 > 0) {
				return searchFaceHierarchy(point, currentTriangle);
			}
		}
	}

    // if we reach this point, that means that the point is repeated, therefore we will ignore him
	return null;
}

function getFaceEdges(faceIndex){
	var edges = new Array(3);

	edges[0] = faceTable[faceIndex].e;
	edges[1] = edgeTable[edges[0]].eN;
	edges[2] = edgeTable[edges[1]].eN;

	return edges;
}

//faceEdges in clockwise order
function triangulateFaceDCEL(faceEdges, newVertex, newFaces, newEdges, newTwinEdges){

	// update DCEL
	var N = newFaces.length;
	for (var i = 0; i < newFaces.length; i++) {	
		// update new dual edge
		edgeTable[newEdges[i]] = {vOrigin: newVertex, fRight: newFaces[i], eN: faceEdges[i], eTwin: newTwinEdges[i]};
		edgeTable[newTwinEdges[i]] = {vOrigin: edgeTable[faceEdges[i]].vOrigin, fRight: newFaces[mod(i-1, N)], eN: newEdges[mod(i-1, N)], eTwin: newEdges[i]};

		// update new face
		faceTable[newFaces[i]] = {e: newEdges[i]};

		// update old edges
		edgeTable[faceEdges[i]].fRight = newFaces[i];
		edgeTable[faceEdges[i]].eN = newTwinEdges[mod(i+1, N)];
	}
}


function updateHierarchyDCEL(isDegenerateCase, oldFaces, newFaces, newEdges) {

	if (!isDegenerateCase){
		var nodeRef = faceToLeaf[oldFaces[0]];
		nodeRef.subtriangles = [
			{e: newEdges[0], subtriangles: []},
			{e: newEdges[1], subtriangles: []},
			{e: newEdges[2], subtriangles: []}
		]
		faceToLeaf[newFaces[0]] = nodeRef.subtriangles[0];
		faceToLeaf[newFaces[1]] = nodeRef.subtriangles[1];
		faceToLeaf[newFaces[2]] = nodeRef.subtriangles[2];
	}
	else{
		var node0Ref = faceToLeaf[oldFaces[0]];
		var node1Ref = faceToLeaf[oldFaces[1]];

		node0Ref.subtriangles = [
			{e: newEdges[0], subtriangles: []},
			{e: newEdges[1], subtriangles: []}
		]
		faceToLeaf[newFaces[0]] = node0Ref.subtriangles[0];
		faceToLeaf[newFaces[1]] = node0Ref.subtriangles[1];

		node1Ref.subtriangles = [
			{e: newEdges[2], subtriangles: []},
			{e: newEdges[3], subtriangles: []}
		]
		faceToLeaf[newFaces[2]] = node1Ref.subtriangles[0];
		faceToLeaf[newFaces[3]] = node1Ref.subtriangles[1];
	}
}

function addPointDCEL(point){

	// add vertex
	var newVertex = vertexCount;
	vertexTable[newVertex] = {x: point.x, y: point.y }; 
	vertexCount += 1;

	// find face belonging to the point
	var searchData = searchFaceHierarchy(point, rootHierarchy); 
	if (searchData == null) return; 	// repeated point, ignore triangulation
	var face = searchData.face; 
	var isDegenerateCase = searchData.degenerateCase;

	if(!isDegenerateCase) {
		// create identifiers
		var faceEdges = getFaceEdges(face);
		var newFaces = [face, faceCount, faceCount+1]; 
		var newEdges = [edgeCount, edgeCount+2, edgeCount+4];
		var newTwinEdges = [edgeCount+1, edgeCount+3, edgeCount+5];

		edgeCount += 6;
		faceCount += 2; 

		// triangulate face and update hierarchy
		triangulateFaceDCEL(faceEdges, newVertex, newFaces, newEdges, newTwinEdges);
		updateHierarchyDCEL(isDegenerateCase, [face], newFaces, newEdges);
	}
	else {	//degenerate case, point on edge
		var edge = searchData.edge;
		var twinEdge = edgeTable[edge].eTwin;
		var oldFace0 = edgeTable[edge].fRight;
		var oldFace1 = edgeTable[twinEdge].fRight;

		var faceEdges = new Array(4);
		faceEdges[0] = edgeTable[edge].eN;
		faceEdges[1] = edgeTable[faceEdges[0]].eN;
		faceEdges[2] = edgeTable[twinEdge].eN;
		faceEdges[3] = edgeTable[faceEdges[2]].eN;

		var newFaces = [edgeTable[edge].fRight, edgeTable[twinEdge].fRight, faceCount, faceCount+1];
		var newEdges = [edge, edgeCount, edgeCount+2, edgeCount+4];
		var newTwinEdges = [twinEdge, edgeCount+1, edgeCount+3, edgeCount+5];

		edgeCount += 8;
		faceCount += 2;

		triangulateFaceDCEL(faceEdges, newVertex, newFaces, newEdges, newTwinEdges);
		updateHierarchyDCEL(isDegenerateCase, [oldFace0, oldFace1], newFaces, newEdges);
	}
}

// function pointInTriangle(p, triangle) {

// 	//determinant signs of the triangles segments
// 	var d_s1 = Math.sign(det(triangle[0], triangle[1], p));
// 	var d_s2 = Math.sign(det(triangle[1], triangle[2], p));
// 	var d_s3 = Math.sign(det(triangle[2], triangle[0], p));
	
// 	//we check where the point lies for every segment (-1:'sideA, 0:'onLine', 1:sideB) to study its case
// 	var insideTriangle = (d_s1 == d_s2 && d_s1 == d_s3);
// 	var onEdge =  (d_s1 == 0 && d_s2 == d_s3) || (d_s2 == 0 && d_s1 == d_s3) || (d_s3 == 0 && d_s1 == d_s2);
// 	var onCorner = Math.abs(d_s1) + Math.abs(d_s2) + Math.abs(d_s3) == 1;

// 	return (insideTriangle || onEdge || onCorner);
//   }

// // linear time version
// function linearVersion(point){
// 	for (var i = 0; i < faceTable.length; i++) {
// 		var faceVertices = getFaceVertices(i);
// 		if (pointInTriangle(point, [ vertexTable[faceVertices[0]], vertexTable[faceVertices[1]], vertexTable[faceVertices[2]] ])) {
// 			return {face: i, degenerateCase: false};
// 		} 
// 	}
// 	return null;
// }

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

function computeContainingTriangle(points){
	if (containingTriangleAdded) {
		points.pop();
		points.pop();
		points.pop();
	} 
	
	var boundingBox = getBoundingBox(points);

	offsetSmallCases = 0.2;
	offsetBigCases = 0.05;
	offset = offsetSmallCases;

	var xOffset = (boundingBox.xmax - boundingBox.xmin) * offset;
	var yOffset = (boundingBox.ymax - boundingBox.ymin) * offset;

	var p0 = {
		x: boundingBox.xmin - (boundingBox.xmax - boundingBox.xmin)/2 - xOffset, 
		y: boundingBox.ymin - yOffset
	};
	var p1 = {
		x: (boundingBox.xmin + boundingBox.xmax)/2, 
		y: boundingBox.ymax + (boundingBox.ymax - boundingBox.ymin)
	};
	var p2 = {
		x: boundingBox.xmax + (boundingBox.xmax - boundingBox.xmin)/2 + xOffset, 
		y: boundingBox.ymin - yOffset
	};

	points.push(p0,p1,p2);
	containingTriangleAdded = true;
}

function getFaceVertices(faceIndex){
	var vertices = new Array(3);
	var edges = getFaceEdges(faceIndex);

	vertices[0] = edgeTable[edges[0]].vOrigin;
	vertices[1] = edgeTable[edges[1]].vOrigin;
	vertices[2] = edgeTable[edges[2]].vOrigin;

	return vertices;
}

function computeTriangulation(points) {

	// compute the Containing Triangle
	computeContainingTriangle(points);
	var N = points.length;

	// compute DCEL
	initializeDCEL(points[N-3], points[N-2], points[N-1]);
	
	for (var i = 0; i < N-3; i++){
		addPointDCEL(points[i]);
	}

	// get the all the triangles of the DCEL
	var outputTriangles = new Array(faceTable.length); 
	for (i=0; i<outputTriangles.length; i++) {
		var faceVertices = getFaceVertices(i);
		outputTriangles[i] = [mod(faceVertices[0] - 3, N), mod(faceVertices[1] -3, N), mod(faceVertices[2] - 3,N)]; // Store INDICES, not points
	}
	return outputTriangles;
}