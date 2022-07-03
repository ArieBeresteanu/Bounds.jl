# Todo list:

This directory contains notebooks that use the code in this repo. Here are some things to do:

[] The sort() function in the Polygon structure does not guarantee that after sorting the first vertex will have both the lowest y-coordinate and the lowest x-coordinate.

[] Minkowski summation doesn't work perfectly. There seem to be an extra vertex at the end.

[] Minkowski summation assumes the two polygons are sorted but doesn't check if this is the case. A possible solution is to add a boolean field that is equal to ```false``` at inception and turned to ```true``` after the sort() function is invoked on the polynomial. 
