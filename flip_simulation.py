import ctypes
import os
import platform

import matplotlib.pyplot as plt
import numpy as np


class FlipSimulation:
    def __init__(self):
        # Load the C library
        if platform.system() == "Windows":
            lib_name = "flip.dll"
        elif platform.system() == "Darwin":  # macOS
            lib_name = "libflip.so"  # or .dylib
        else:  # Linux and others
            lib_name = "libflip.so"

        lib_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), lib_name)
        self.lib = ctypes.CDLL(lib_path)

        # Define FLIP fluid function signatures
        self.lib.createFlipFluid.argtypes = [
            ctypes.c_float,  # density
            ctypes.c_float,  # width
            ctypes.c_float,  # height
            ctypes.c_float,  # spacing
            ctypes.c_float,  # particleRadius
            ctypes.c_int,  # maxParticles
        ]
        self.lib.createFlipFluid.restype = ctypes.c_void_p  # Returns FlipFluid*

        self.lib.destroyFlipFluid.argtypes = [ctypes.c_void_p]  # FlipFluid*
        self.lib.destroyFlipFluid.restype = None

        self.lib.simulate.argtypes = [
            ctypes.c_void_p,  # FlipFluid*
            ctypes.c_float,  # dt
            ctypes.c_float,  # gravity
            ctypes.c_float,  # flipRatio
            ctypes.c_int,  # numPressureIters
            ctypes.c_int,  # numParticleIters
            ctypes.c_float,  # overRelaxation
            ctypes.c_int,  # compensateDrift
            ctypes.c_int,  # separateParticles
            ctypes.c_float,  # obstacleX
            ctypes.c_float,  # obstacleY
            ctypes.c_float,  # obstacleRadius
            ctypes.c_float,  # obstacleVelX
            ctypes.c_float,  # obstacleVelY
        ]
        self.lib.simulate.restype = None

        # Add particle management functions
        self.lib.add_particle.argtypes = [
            ctypes.c_void_p,  # FlipFluid*
            ctypes.c_float,  # x
            ctypes.c_float,  # y
        ]
        self.lib.add_particle.restype = (
            ctypes.c_int
        )  # Returns 1 if successful, 0 otherwise

        self.lib.get_particle_count.argtypes = [ctypes.c_void_p]  # FlipFluid*
        self.lib.get_particle_count.restype = ctypes.c_int

        self.lib.get_particles.argtypes = [ctypes.c_void_p]  # FlipFluid*
        self.lib.get_particles.restype = ctypes.POINTER(ctypes.c_float)

        # Cell data functions
        self.lib.get_cell_type.argtypes = [
            ctypes.c_void_p,  # FlipFluid*
            ctypes.c_int,  # x
            ctypes.c_int,  # y
        ]
        self.lib.get_cell_type.restype = (
            ctypes.c_int
        )  # Returns cell type (FLUID_CELL, AIR_CELL, SOLID_CELL)

        self.lib.get_cell_density.argtypes = [
            ctypes.c_void_p,  # FlipFluid*
            ctypes.c_int,  # x
            ctypes.c_int,  # y
        ]
        self.lib.get_cell_density.restype = ctypes.c_float

        # Boundary setting function
        self.lib.set_solid_boundaries.argtypes = [ctypes.c_void_p]  # FlipFluid*
        self.lib.set_solid_boundaries.restype = None

        # Initialize instance variables
        self.width = 0
        self.height = 0
        self.fluid = None  # Will hold the FlipFluid pointer
        self.spacing = 1.0  # Default grid spacing

    def create_grid(self, width, height):
        """Initialize the simulation grid with the given dimensions."""
        self.width = width
        self.height = height
        return True

    def set_fluid_parameters(
        self, density=1000.0, particle_radius=0.1, max_particles=100000
    ):
        """Set up the FLIP fluid simulation with given parameters."""
        if self.fluid:
            self.lib.destroyFlipFluid(self.fluid)

        # Create a new fluid simulator
        self.fluid = self.lib.createFlipFluid(
            density,
            float(self.width),
            float(self.height),
            self.spacing,
            particle_radius,
            max_particles,
        )

        return self.fluid is not None

    def set_solid_boundaries(self):
        """Set up solid boundaries for the tank walls."""
        if not self.fluid:
            return False

        self.lib.set_solid_boundaries(self.fluid)
        return True

    def add_particle(self, x, y):
        """Add a fluid particle at the given position."""
        if not self.fluid:
            return False
        return bool(self.lib.add_particle(self.fluid, float(x), float(y)))

    def simulate(
        self,
        dt,
        gravity,
        flip_ratio,
        num_pressure_iters,
        num_particle_iters,
        over_relaxation,
        compensate_drift,
        separate_particles,
        obstacle_x,
        obstacle_y,
        obstacle_radius,
        obstacle_vel_x,
        obstacle_vel_y,
    ):
        """Run one step of the FLIP fluid simulation."""
        if not self.fluid:
            return False

        self.lib.simulate(
            self.fluid,
            dt,
            gravity,
            flip_ratio,
            num_pressure_iters,
            num_particle_iters,
            over_relaxation,
            1 if compensate_drift else 0,
            1 if separate_particles else 0,
            obstacle_x,
            obstacle_y,
            obstacle_radius,
            obstacle_vel_x,
            obstacle_vel_y,
        )
        return True

    def get_cell_type(self, x, y):
        """Get the type of a cell (FLUID_CELL, AIR_CELL, or SOLID_CELL)."""
        if not self.fluid:
            return -1
        return self.lib.get_cell_type(self.fluid, int(x), int(y))

    def get_cell_density(self, x, y):
        """Get the fluid density at a cell."""
        if not self.fluid:
            return 0.0
        return self.lib.get_cell_density(self.fluid, int(x), int(y))

    def get_particles(self):
        """Get all particle data (positions and colors)."""
        if not self.fluid:
            return np.array([])

        count = self.get_particle_count()
        if count == 0:
            return np.array([])

        # Get pointer to particle data (x, y, r, g, b for each particle)
        particles_ptr = self.lib.get_particles(self.fluid)
        if not particles_ptr:
            return np.array([])

        # Create numpy array from particles data (5 values per particle)
        particles = np.ctypeslib.as_array(particles_ptr, shape=(count * 5,))

        # Make a copy to avoid memory issues
        return np.copy(particles)

    def get_particle_count(self):
        """Get the current number of fluid particles."""
        if not self.fluid:
            return 0
        return self.lib.get_particle_count(self.fluid)

    def set_cell(self, x, y, value):
        """Placeholder for backward compatibility - does nothing in current implementation."""
        return True

    def get_cell(self, x, y):
        """Placeholder for backward compatibility - returns density if fluid exists, otherwise 0."""
        if self.fluid:
            return self.get_cell_density(x, y)
        return 0.0

    def get_grid_data(self):
        """Returns a grid of cell densities as a numpy array."""
        if not self.fluid or self.width == 0 or self.height == 0:
            return None

        # Create a new grid and fill it with cell densities
        grid = np.zeros((self.height, self.width), dtype=np.float32)

        for y in range(self.height):
            for x in range(self.width):
                grid[y, x] = self.get_cell_density(x, y)

        return grid

    def visualize_grid(self, cmap="viridis", vmin=None, vmax=None):
        """Visualize the density grid using matplotlib."""
        grid_data = self.get_grid_data()
        if grid_data is None:
            print("No grid data available for visualization")
            return

        plt.figure(figsize=(10, 8))
        plt.imshow(grid_data, cmap=cmap, origin="lower", vmin=vmin, vmax=vmax)
        plt.colorbar(label="Density")
        plt.title(f"Density Visualization ({self.width}x{self.height})")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    def visualize_particles(self):
        """Visualize particles using matplotlib."""
        particles = self.get_particles()
        if len(particles) == 0:
            print("No particles to visualize")
            return

        # Reshape to get positions and colors separately
        positions = particles.reshape(-1, 5)[:, :2]
        colors = particles.reshape(-1, 5)[:, 2:5]

        plt.figure(figsize=(10, 8))
        plt.scatter(positions[:, 0], positions[:, 1], c=colors, s=20, alpha=0.7)
        plt.title(f"Particle Visualization (Count: {len(positions)})")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis("equal")
        plt.show()

    def __del__(self):
        """Clean up resources when the object is destroyed."""
        if hasattr(self, "fluid") and self.fluid:
            self.lib.destroyFlipFluid(self.fluid)
