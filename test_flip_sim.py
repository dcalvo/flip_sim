import time

import numpy as np
import pygame

from flip_simulation import FlipSimulation


class FlipSimulationApp:
    def __init__(self, width=50, height=30, cell_size=15):
        # Simulation parameters
        self.width = width
        self.height = height
        self.cell_size = cell_size
        self.screen_width = width * cell_size
        self.screen_height = height * cell_size

        # Initialize simulation
        self.sim = FlipSimulation()
        self.sim.create_grid(width, height)

        # Initialize display
        pygame.init()
        self.screen = pygame.display.set_mode((self.screen_width, self.screen_height))
        pygame.display.set_caption("FLIP Simulation")
        self.clock = pygame.time.Clock()

        # FPS tracking
        self.frame_count = 0
        self.start_time = time.time()

        # Mouse interaction state
        self.dragging = False
        self.ball_radius = 5  # Radius in cells
        self.ball_value = 1.0  # Value to set for the ball cells

        # Initialize grid with pattern
        self.initialize_grid()

    def initialize_grid(self):
        """Initialize the grid with a pattern"""
        for x in range(self.width):
            for y in range(self.height):
                # Create a simple pattern (circular wave)
                distance = np.sqrt(
                    (x - self.width / 2) ** 2 + (y - self.height / 2) ** 2
                )
                value = np.sin(distance / 3) * 0.5 + 0.5
                self.sim.set_cell(x, y, value)

    def get_color(self, value):
        """Map a value (0-1) to a color gradient"""
        r = min(255, int(value * 255))
        g = min(255, int((1 - value) * 128))
        b = min(255, int((1 - value) * 255))
        return (r, g, b)

    def update_simulation(self):
        """Update the simulation state"""
        # For demonstration, we'll animate the pattern by shifting the phase
        # In a real simulation, you would call sim.update() here instead
        phase = time.time() * 2  # Time-based animation
        for x in range(self.width):
            for y in range(self.height):
                distance = np.sqrt(
                    (x - self.width / 2) ** 2 + (y - self.height / 2) ** 2
                )
                value = np.sin(distance / 3 + phase) * 0.5 + 0.5
                self.sim.set_cell(x, y, value)

    def create_ball(self, center_x, center_y):
        """Create a ball of cells with the specified value around the center point"""
        for x in range(
            max(0, int(center_x - self.ball_radius)),
            min(self.width, int(center_x + self.ball_radius + 1)),
        ):
            for y in range(
                max(0, int(center_y - self.ball_radius)),
                min(self.height, int(center_y + self.ball_radius + 1)),
            ):
                distance = np.sqrt((x - center_x) ** 2 + (y - center_y) ** 2)
                if distance <= self.ball_radius:
                    # You can use a falloff function if you want a gradient ball
                    # value = self.ball_value * (1 - distance / self.ball_radius)
                    self.sim.set_cell(x, y, self.ball_value)

    def handle_events(self):
        """Process pygame events"""
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    return False

            # Mouse events for ball creation/dragging
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:  # Left click
                    self.dragging = True
                    mouse_x, mouse_y = event.pos
                    cell_x = mouse_x / self.cell_size
                    cell_y = mouse_y / self.cell_size
                    self.create_ball(cell_x, cell_y)

            elif event.type == pygame.MOUSEMOTION:
                if self.dragging:
                    mouse_x, mouse_y = event.pos
                    cell_x = mouse_x / self.cell_size
                    cell_y = mouse_y / self.cell_size
                    self.create_ball(cell_x, cell_y)

            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:  # Left click
                    self.dragging = False

        return True

    def render(self):
        """Render the current state to the screen"""
        grid_data = self.sim.get_grid_data()

        # Clear the screen
        self.screen.fill((0, 0, 0))

        # Render each cell
        for x in range(self.width):
            for y in range(self.height):
                rect = pygame.Rect(
                    x * self.cell_size,
                    y * self.cell_size,
                    self.cell_size,
                    self.cell_size,
                )
                color = self.get_color(grid_data[y, x])
                pygame.draw.rect(self.screen, color, rect)

        # Update the display
        pygame.display.flip()

    def update_fps_display(self):
        """Update and display the FPS counter"""
        self.frame_count += 1
        if self.frame_count % 30 == 0:
            current_time = time.time()
            fps = self.frame_count / (current_time - self.start_time)
            pygame.display.set_caption(f"FLIP Simulation - FPS: {fps:.2f}")
            self.frame_count = 0
            self.start_time = current_time

    def print_stats(self):
        """Print statistics about the simulation"""
        grid_data = self.sim.get_grid_data()
        test_x, test_y = self.width // 2, self.height // 2
        cell_value = self.sim.get_cell(test_x, test_y)
        print(f"Value at ({test_x}, {test_y}): {cell_value}")
        print("Grid shape:", grid_data.shape)
        print("Grid min value:", grid_data.min())
        print("Grid max value:", grid_data.max())

    def run(self):
        """Main simulation loop"""
        running = True

        while running:
            # Handle user input
            running = self.handle_events()

            # Update simulation state
            self.update_simulation()

            # Render current state
            self.render()

            # Update FPS display
            self.update_fps_display()

            # Cap the frame rate
            self.clock.tick(60)  # Target 60 FPS

        # Clean up and print stats
        pygame.quit()
        self.print_stats()


def main():
    app = FlipSimulationApp()
    app.run()


if __name__ == "__main__":
    main()
