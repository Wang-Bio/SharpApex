# SM1 Curvature 

In this part, I will briefly introduce the concept of curvature and the methods to calculate it, by differential geometrics and also by circle fitting. Please notice that circle fitting method is what I used in the model. 

## SM1.0 Prblem Statement

Assume that we have an organ, either with image files (eg., TIFF file, PNG file) or we have simulated organs (Class Organ) by this code, how can we calculate the curvature of the organ contour?

## SM1.1 Definition and Understanding of Curvature (in 2D space)
1. Intuition Understanding: Curvature is the amount by which a curve deviates from being a straight line, or a surface deviates from being a plane. This could be directly observed and intuitively percieved by human eyes and brain. 


2. Canonical definition of curvature of a circle: the curvature of a cirlce equal to the reciprocal of its radius. Smaller circles bend more sharply, and hence have higher curvature. 

3. Definition of the curvature at a point of a differentiable curve (circle fitting): the curvature of its osculating circle, that is the circle that best approximates the curve near this point. 

## SM1.2 Calculation of curvature in a continuous curve (curve defined by explicit function; calculated by differential geometrics)
1. arc-length parameterization

2. general parameterization

## SM1.3 Calculation of curvature for a discrete curve
1. Preparation: Sampling a series of equal-distance and clockwise points from the organ contour 

2. Circle Fitting: Circle fitting of three points based on Kasa Method

3. Analysis of curvature: the Spectrum of Curvature