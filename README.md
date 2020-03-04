WHAT IS THIS PROJECT ABOUT?
It is a small project completly written in C during my Computer Graphics's university course.
Surface Mesh Simplification is the process that reduces the number of faces of a mesh.

STRATEGIES IMPLEMENTED
In this work you can choose three possible strategies to simplify your mesh:
1. Midpoint solution. It is the simplest and fastest way (but be careful, degrades the quality of the mesh quickly. I suggest this strategy only if you need fast performances and have a huge mesh (500k triangles or more) and reduce to no more than 80/83% of triangles).
2. Lindstrom-Turk solution. This is based on "Fast and Memory Efficient Polygonal Simplification" written by Peter Lindstrom and Greg Turk. Overall, in my opinion, keeps quite a nice form of simplified mesh. Be careful not to reduce to a mesh with less than 200/300 triangles or it can show bizarre aspects.
3. Garland and Heckbert solution. This approach is based on "Surface Simplification Using Quadric Error Metrics", advantage respect to Lindstrom-Turk is that computes really faster. Sometimes I found this solution more accurate, other less. As strategy 2, don't reduce to less than 200/300 triangles. 


WHAT IS NEEDED TO RUN THIS CODE?
- C compiler, I used gcc  
- graphics libraries: OpenGL, glu, glut 

FIRST RUN GUIDE 
in the shell:
if you are under Linux, compile with BLABLA 
If you are under Mac, compile with BLABLA
Then run with "./a.out lib/bunny.off 1 0.1 0" [lib/bunny.off is the file you are simplifying, 1 indicates Lindstrom-Turk simplification alghoritm (0 for midpoint, 1 Lindstrom-Turk, 2 Garland), 0.1 indicates a check for distance from new solution to old one (a security guard to avoid drastic changes in shape), 0 indicates the number of simplifications per step you wish (0 means default which is approximatively 1/6 of all the triangles, you can choose any integer < n_tringles)]
If everything went smooth, you should have opened a new windows with your mesh. You can move with "w a s d" in every direction, zoom in or out "+ -", or can simplify clicking "k" on your keyboard. When you are satisifed with simplification, with the mesh in foreground, click "Esc". It will output and terminate your program.

WHICH MESHES CAN I INPUT?
Every mesh formed only of triangles can be used. It shall be in this format:
Number of vertices number of triangles 
coordinates for each vertex (one per row)
3 vertex1 vertex2 vertex3 (where a vertex is identified by its number, first vertex is 0)
(if you have any doubt, please read any mesh in input folder)

