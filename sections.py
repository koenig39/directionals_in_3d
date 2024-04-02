import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Define points and normals for planes and line segment points
P1, N1 = np.array([0, 0, 1]), np.array([0, 0, 1])  # Plane 1: z = 1
P2, N2 = np.array([0, 0, -1]), np.array([0, 0, -1])  # Plane 2: z = -1
L1, L2 = np.array([-1, -1, -2]), np.array([1, 1, 2])  # Line segment


def line_plane_intersection(plane_point, plane_normal, line_point1, line_point2):
    line_direction = line_point2 - line_point1
    dot_product = np.dot(plane_normal, line_direction)
    if dot_product == 0:
        return None  # Line and plane are parallel
    t = np.dot(plane_normal, (plane_point - line_point1)) / dot_product
    if 0 <= t <= 1:  # Intersection point is within the line segment
        return line_point1 + t * line_direction
    return None


# Calculate intersection points again (reuse previous function)
intersection1 = line_plane_intersection(P1, N1, L1, L2)
intersection2 = line_plane_intersection(P2, N2, L1, L2)

# Setup for plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot original line
ax.plot([L1[0], L2[0]], [L1[1], L2[1]], [L1[2], L2[2]], color='gray', linewidth=1, linestyle='--')

# Plot segments if intersections exist
if intersection1 is not None and intersection2 is not None:
    # Plot from L1 to intersection2
    ax.plot([L1[0], intersection2[0]], [L1[1], intersection2[1]], [L1[2], intersection2[2]], color='blue', linewidth=3)
    # Plot from intersection2 to intersection1
    ax.plot([intersection2[0], intersection1[0]], [intersection2[1], intersection1[1]], [intersection2[2], intersection1[2]], color='red', linewidth=3)
    # Plot from intersection1 to L2
    ax.plot([intersection1[0], L2[0]], [intersection1[1], L2[1]], [intersection1[2], L2[2]], color='green', linewidth=3)

# Plot planes (as large flat surfaces for visualization)
# For Plane 1 (z=1)
X, Y = np.meshgrid(np.arange(-2, 3), np.arange(-2, 3))
Z = np.ones(X.shape)
ax.plot_surface(X, Y, Z, alpha=0.5, color='orange')

# For Plane 2 (z=-1)
Z = -np.ones(X.shape)
ax.plot_surface(X, Y, Z, alpha=0.5, color='purple')

# Setting labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Show plot
plt.show()
