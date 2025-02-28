import ctypes
import os
import platform
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


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

        # Define function signatures
        self.lib.create_grid.argtypes = [ctypes.c_int, ctypes.c_int]
        self.lib.create_grid.restype = ctypes.c_int

        self.lib.set_cell.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_float]
        self.lib.set_cell.restype = ctypes.c_int

        self.lib.get_cell.argtypes = [ctypes.c_int, ctypes.c_int]
        self.lib.get_cell.restype = ctypes.c_float

        self.lib.get_grid_data.restype = ctypes.POINTER(ctypes.c_float)

        self.lib.get_grid_dimensions.argtypes = [
            ctypes.POINTER(ctypes.c_int),
            ctypes.POINTER(ctypes.c_int),
        ]

        self.lib.destroy_grid.argtypes = []
        self.lib.destroy_grid.restype = None

        # Initialize instance variables
        self.width = 0
        self.height = 0

    def create_grid(self, width, height):
        """Create a grid with the given dimensions."""
        result = self.lib.create_grid(width, height)
        if result:
            self.width = width
            self.height = height
            return True
        return False

    def set_cell(self, x, y, value):
        """Set the value of a cell in the grid."""
        return bool(self.lib.set_cell(x, y, value))

    def get_cell(self, x, y):
        """Get the value of a cell in the grid."""
        return self.lib.get_cell(x, y)

    def get_grid_data(self):
        """Get the entire grid as a numpy array."""
        width = ctypes.c_int()
        height = ctypes.c_int()
        self.lib.get_grid_dimensions(ctypes.byref(width), ctypes.byref(height))

        if width.value == 0 or height.value == 0:
            return None

        # Get pointer to the raw grid data in C
        data_ptr = self.lib.get_grid_data()
        if not data_ptr:
            return None

        # Create a numpy array that references the C array's memory
        # This is a view, not a copy
        buffer_size = width.value * height.value
        grid_data = np.ctypeslib.as_array(data_ptr, shape=(buffer_size,))

        # Reshape to 2D array and make a copy to avoid memory issues
        return np.copy(grid_data.reshape((height.value, width.value)))

    def visualize_grid(self, cmap="viridis", vmin=None, vmax=None):
        """Visualize the current grid state using matplotlib."""
        grid_data = self.get_grid_data()
        if grid_data is None:
            print("No grid data available for visualization")
            return

        plt.figure(figsize=(10, 8))
        plt.imshow(grid_data, cmap=cmap, origin="lower", vmin=vmin, vmax=vmax)
        plt.colorbar(label="Value")
        plt.title(f"Grid Visualization ({self.width}x{self.height})")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    def __del__(self):
        """Clean up resources when the object is destroyed."""
        self.lib.destroy_grid()
