# Rust-Maze-Generator
A 2d maze generation system written in Rust.

I wrote this maze generator as a fun little project to work with some different ways to generate mazes.

There's a few maze types that this project supports, which I'll explain below:

Basic Erosion Maze:
The basic erosion maze works by converting an certain amount of tiles (based on a given weight parameter in its generation) that can be turned into walls while safely avoiding the possibility of creating enclosed areas. Due to this, it lacks dense areas and is more clusters and spreads of single-few blocks. 
<img width="762" height="757" alt="image" src="https://github.com/user-attachments/assets/5ff530ff-7e09-4d6e-8c94-a7681ab64148" />

Tunnel Maze:
The tunnel maze works by creating pathways in a world that is originally set to fully be walls. It does this by creating "tunnels", which are somewhat randomly generated paths with a provided start and end point. It will create tunnels until a certain amount of tiles (depending on a weight parameter) are turned into paths. To help with promoting tunnel variety and nicer looking generations, while still not being extremely slow, parameters involving distances for tunnels are also provided. An option to fill 1x1 holes within tunnels is also provided.
<img width="756" height="759" alt="image" src="https://github.com/user-attachments/assets/e473c6b7-a7a1-4429-b8bb-a4f91d8101d7" />

Branch Maze:
Similar to the tunnel maze, the branch maze cuts pathways into a world initially set as all walls. However, unlike the tunnel maze, the branch maze uses "branches" to do this. Branches have a start point, but no definite end point--only a length. When a branch is generated, the tiles along it have a chance to "sway" (go left or right for vertical branches, or up and down for horizontal branches). Like the tunnel maze, once a certain portion of the world has become open tiles (based on the weight parameter), the process will stop.
<img width="756" height="758" alt="image" src="https://github.com/user-attachments/assets/f49a9c3e-e0d7-400e-b5bd-af95d6cf3782" />

Random Tunnel Maze:
The random tunnel maze also creates pathways in an originally all-wall world. However, where it differs from the branch and the tunnel mazes is that it doesn't necessarily have an exact general direction. This generation creates a random tunnel moves randomly (left, right, up, or down), which again creates pathways until a certain number depending on the weight parameter is reached. Weighting for the direction that the random tunnel makes can be specified, allowing for random tunnels to focus on moving more in a certain direction. Directional weighting is also an option, which allows for random tunnels to move more in the opposite direction of their original position (hence allowing for overall neater map coverage). Because it's rather random, the generations on this maze will not be consistent, but it's a fun one.
<img width="759" height="760" alt="image" src="https://github.com/user-attachments/assets/c2796f28-4204-4b3a-9e70-852ba943196d" />

Branched Erosion Maze:
Another erosion maze, this creates walls in an open world. However, unlike the basic erosion maze, the branched erosion maze adds a twist that allows for the more dense and compact areas some may be looking for: it generates branches along a set span, allowing for more clusters.
<img width="747" height="758" alt="image" src="https://github.com/user-attachments/assets/f546fc10-a7a7-461c-9e45-1949099c806e" />





