# 2D Hough transform for line detection

## Slope-intercept representation

The slope-intercept representation of the Hough transform transforms points in the $xy$ space to straight lines in the $ab$ space using the equation $b = -x_ia + y_i$. An accumulator array is generated to capture the local maxima of intersections between lines (ab plane) to find the pixels in the source image located along the line (xy plane) with gradient $a$ and intercept $b$.

As the size of the discretized ab accumulator array depends on the gradient $a$ and intercept $b$, it becomes incomputably large as the gradient of a line in the $xy$ plane approaches vertical (infinity). The runtime for the ab transform hence takes considerably longer than the second representation.

## Normal representation

The normal representation of the Hough transform uses the $\rho\theta$ space from the equation $\rho = x_i\cos(\theta)+y_i\sin(\theta)$ which solves the above issue as $\theta$ can be restricted to the range $-90\degree < \theta \leq 90\degree$. 

The `ht_plot_aligned` function in the [rt representation](hough_transform_rt_run.m) aligns the boundary of the original image to the lines found from the rt hough transform. 

