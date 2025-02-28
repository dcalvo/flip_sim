import numpy as np

from flip_simulation import FlipSimulation


def main():
    # Create a simulation instance
    sim = FlipSimulation()

    # Create a 50x30 grid
    width, height = 50, 30
    if not sim.create_grid(width, height):
        print("Failed to create grid")
        return

    # Set some test values
    # Create a pattern for testing
    for x in range(width):
        for y in range(height):
            # Create a simple pattern (e.g., circular wave)
            distance = np.sqrt((x - width / 2) ** 2 + (y - height / 2) ** 2)
            value = np.sin(distance / 3) * 0.5 + 0.5
            sim.set_cell(x, y, value)

    # Read back a specific cell
    test_x, test_y = 25, 15
    cell_value = sim.get_cell(test_x, test_y)
    print(f"Value at ({test_x}, {test_y}): {cell_value}")

    # Get and display the entire grid
    grid_data = sim.get_grid_data()
    if grid_data is not None:
        print("Grid shape:", grid_data.shape)
        print("Grid min value:", grid_data.min())
        print("Grid max value:", grid_data.max())

    # Visualize the grid
    sim.visualize_grid()


if __name__ == "__main__":
    main()
