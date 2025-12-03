# Rust-Maze-Generator
A 2d maze generation system written in Rust.

I wrote this maze generator as a fun little project to work with some different ways to generate mazes.

There's a few maze types that this project supports, which I'll explain below:

Basic Erosion Maze:
The basic erosion maze works by converting an certain amount of tiles (based on a given weight parameter in its generation) that can be turned into walls while safely avoiding the possibility of creating enclosed areas. Due to this, it lacks dense areas and is more clusters and spreads of single-few blocks. 

Tunnel Maze:
The tunnel maze works by creating pathways in a world that is originally set to fully be walls. It does this by creating "tunnels", which are somewhat randomly generated paths with a provided start and end point. It will create tunnels until a certain amount of tiles (depending on a weight parameter) are turned into paths. To help with promoting tunnel variety and nicer looking generations, while still not being extremely slow, parameters involving distances for tunnels are also provided. An option to fill 1x1 holes within tunnels is also provided.

Branch Maze:
Similar to the tunnel maze, the branch maze cuts pathways into a world initially set as all walls. However, unlike the tunnel maze, the branch maze uses "branches" to do this. Branches have a start point, but no definite end point--only a length. When a branch is generated, the tiles along it have a chance to "sway" (go left or right for vertical branches, or up and down for horizontal branches). Like the tunnel maze, once a certain portion of the world has become open tiles (based on the weight parameter), the process will stop. An option to fill 1x1 holes also exists for this maze generation.

Random Tunnel Maze:
The random tunnel maze also creates pathways in an originally all-wall world. However, where it differs from the branch and the tunnel mazes is that it doesn't necessarily have an exact general direction. This generation creates a random tunnel moves randomly (left, right, up, or down), which again creates pathways until a certain number depending on the weight parameter is reached. Weighting for the direction that the random tunnel makes can be specified, allowing for random tunnels to focus on moving more in a certain direction. Directional weighting is also an option, which allows for random tunnels to move more in the opposite direction of their original position (hence allowing for overall neater map coverage). Because it's rather random, the generations on this maze will not be consistent, but it's a fun one.

Branched Erosion Maze:
Another erosion maze, this creates walls in an open world. However, unlike the basic erosion maze, the branched erosion maze adds a twist that allows for the more dense and compact areas some may be looking for: it generates branches along a set span, allowing for more clusters.

In addition, there is a method called "square_quantize_maze", which can return the wall bodies at coordinates, which are grouped into squares. That is, instead of being many individual tiles, each wall in the result of this method will be given a magnitude to it, yielding a square of varying size that still equates to what the original grid was. This allows for the simplification and reduction of walls generated, while still yielding the same structure for every maze. If you're looking to use this in an environment that demands less walls (for performance reasons), this will do the trick while still making sure all walls maintain a square shape (and are not rectangular). 
