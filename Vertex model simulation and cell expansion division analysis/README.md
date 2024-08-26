# Welcome to the PlantVertexModel !

# 0. What is PlantVertexModel

PlantVertexModel is a program used for simulate and analyze the plant organ morphogenesis based on cells behavior. The plant cells are abstracted as polygons (idealy, regular polygons), and the plant organs are abstracted as networks of polygons. Cell behavior including cell expansion, cell division in different spatial and temporal pattern can be performed to simulated the growing of plant organs (dynamics of polygon network).

<img src="/Guide/model_introduction.jpg" alt="PlantVertexModel introduction" title="PlantVertexModel introduction" width="500"/>

This Guide is divided into three parts:

1. [PlantVertexModel User's Guide ](Guide/User's%20Guide): covers the basic rules and usage of PlantVertexModel. The related mathematical, physical and coding details could be found in the reference mannual. 

2. [PlantVertexModel Reference Manual](Guide/Reference%20Manual): dives into the mathematical, physical details and coding algorithm of PlantVertexModel, corresponding to the User's Guide part.

3. [PlantVertexModel Tutorials](Guide/Tutorials): step-by-step tutorials to help you do the PlantVertexModel simulation and analysis by using the examples. If you are in hurry, it is also possible to directly start the tutorials and feel how the vertex model works. If you feel confused during tutorial, please refer to the User's Guide or Reference mannual.

# 1. Basics of PlantVertexModel
To do a PlantVertexModel simulation, we need the inputs, including mode.txt, which defines the mode we want to use; parameter.txt, gives the parameter values; InitialOrgan.txt, defines the initial shape of organ primordia. These inputs will initialize the simulation. The program will start to calculate cell expansion, cell division based on defined rules and input parameters. Finally, the program will output the organ geometric information on cell and organ, cell division records, system log and records. These outputs could be used to visualize and analyze the simulated organ morphogenesis. PlantVertexModel also provides analysis codes for the simulation results and also experiments results. 

<img src="/Guide/flowchart_simulation.jpg" alt="The flowchart for simulation" title="Flowchart for PlantVertexModel Simulation" width="500"/>

# 2. PlantVertexModel User's Guide
* 2.1 Introduction to Vertex Model (What is Vertex Model and a Review of Previous Research)
* 2.2 Physics of plant cells in Vertex Model 
* 2.3 Division rules of plant cells in Vertex Model 
* 2.4 Geometric analysis of cells and organs in Vertex Model
* 2.5 Creating your own simulated plant organ 

# 3. PlantVertexModel Reference Mannual
* 3.1 Basics of Vertex Model (C++ class and FlowChart)
* 3.2 Physics of plant cells
* 3.3 Division rules of plant cells in Vertex Model
* 3.4 Geometric analysis of cells and organs in Vertex Model
* 3.5 How to design rule to create your own simulated plant organ

# 4. PlantVertexModel Tutorials: 
* [4.0 Basics for C++ and CMake](Guide/Tutorials/0%20C++%20and%20CMake.md): teach you the basics of compiling C++ files by CMake and do a Hello World.
* [4.1 A circular organ generation](Guide/Tutorials/1%20A%20Circular%20Organ.md): teach you how to generate a circular organ shape by setting cell division direction and cell division frequency control. 
* [4.2 An elliptical organ generation](Guide/Tutorials/2%20An%20Elliptical%20Organ.md): teach you how to generate an elliptical organ shape by setting biased cell division direction.
* [4.3 A petal-like organ generation](Guide/Tutorials/3%20A%20Petal-like%20Shape.md): teach you how to generate a petal-like organ shape by setting apical meristem position. Reference: [https://doi.org/10.1242/dev.199773](https://doi.org/10.1242/dev.199773)
