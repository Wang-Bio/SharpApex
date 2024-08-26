# Define the circle equation as a parametric function
set parametric
set urange [-pi:pi]
set vrange [-1:1]

# Define the radius of the circle
radius = 0.5
center_x = 0
center_y = 0

# Parametric functions for the circle
fx(t) = center_x + radius*cos(t)
fy(t) = center_y + radius*sin(t)

# Define the Gaussian function
gauss(x, a, b, c) = a * exp(-(x-b)**2 / (2*c**2))

# Set the range for the plot to include the circle and the Gaussian
set xrange [-1:1]
set yrange [-1:1]

# Set an appropriate size for the plot window to ensure aspect ratio
set size ratio -1

# Plot the circle and add the Gaussian to the y-values
# Here 'a' is the amplitude, 'b' is the center, and 'c' controls the width of the Gaussian.
fy_mod(t) = fy(t) > 0 ? fy(t) + gauss(fx(t), 0.4, 0, 0.05) : fy(t)
plot fx(t), fy_mod(t) with points lc rgb "dark-green" pt 7 notitle