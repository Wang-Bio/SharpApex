# 2. Geometrics and Physics of Cells in PlantVertexModel
## 2.1 Geometrics of Cells
Cells in 2D vertex model are defined as polygons, or the vertices and the connections of vertices. There are four different biological levels in vertex model: 
1. vertex, or point, which constitutes lines, cells, and organ
2. lines, which are constituted by vertices, and constitutes cells and organ
3. cells, which are constituted by vertices and line, and constitutes the organ
4. organ, which is made by vertex, line, and cells 

For examples in 61cells.txt, the primordia organ is made by 120 vertices, 180 lines, and 61 cells. By default setting, the simulation ends where cell number in the organ reaches 1,000, where vertex number reaches around 2,000 and line number reaches around 3,000.

## 2.2 Physics of Cells
There are several definition of cell physics in different vertex model (citation). In this PlantVertexModel, we choosed the vertex model from (citation Hayakawa et al., 2016, also seen in Kinoshita&Naito&Wang et al., 2022):


## 2.3 Geometric Attributes of Organ, Cells, Lines, and Vertices
### 2.3.1 Basic Geometric Attributes of Cells, Lines, and Vertices
Area of a cell/organ, Perimeter of a cell/organ, Length of a line, position of a vertex

### 2.3.2 Advance Geometric Attributes of Organ, Cells
Regularity of a cell/organ, Area rank of a cell/organ, Circularity of cell/organ, 
### 2.3.3 Epidermal Identity of Cells, Lines, and Vertices
Epidermal cell, Outermost line, Surface vertex
### 2.3.4 Layers of Cells
### 2.3.5 Overlap of Organ
### 2.3.6 Length, Width and length/width ratio of Organ
### 2.3.7 Cell arrangement index of organ
