# 2D Hough transform for line detection

## Slope-intercept representation

The slope-intercept representation of the Hough transform transforms points in the $xy$ space to straight lines in the $ab$ space using the equation $b = -x_ia + y_i$. An accumulator array is generated to capture the local maxima of intersections between lines (ab plane) to find the pixels in the source image located along the line (xy plane) with gradient $a$ and intercept $b$.
