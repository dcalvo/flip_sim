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

        # Scene parameters matching original JS implementation
        self.gravity = -9.81  # Negative for downward gravity
        self.dt = 1.0 / 60.0  # Time step
        self.flip_ratio = 0.95  # Blend between PIC (0) and FLIP (1)
        self.num_pressure_iters = 50  # Match JS (was 20)
        self.num_particle_iters = 2  # Particle separation iterations
        self.frame_nr = 0
        self.over_relaxation = 1.9  # SOR parameter
        self.compensate_drift = True  # Flag to compensate for density drift
        self.separate_particles = True  # Flag to enable particle separation

        # Obstacle properties
        h = self.height / 100.0  # Cell size (scaled to match JS)
        self.obstacle_x = int(self.width * 0.8)  # Match JS
        self.obstacle_y = int(self.height * 0.8)  # Match JS
        self.obstacle_radius = 0.15 * min(self.width, self.height)  # Match JS scaling
        self.obstacle_vel_x = 0.0
        self.obstacle_vel_y = 0.0

        # Visualization flags matching JS
        self.paused = False
        self.show_obstacle = True
        self.show_particles = True
        self.show_grid = False  # Start with grid hidden like JS

        # Initialize simulation
        self.sim = FlipSimulation()
        self.setup_simulation()

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
        self.is_moving_obstacle = False
        self.ball_radius = 5  # Radius in cells
        self.ball_value = 1.0  # Value to set for the ball cells

    def setup_simulation(self):
        """Initialize the simulation with current parameters"""
        # Create fluid grid with appropriate dimensions and density
        self.sim.create_grid(self.width, self.height)

        # Calculate particle radius as a proportion of cell size (matching JS)
        h = 1.0  # Cell size in simulation units
        r = 0.3 * h  # Particle radius (matches JS setup)

        self.sim.set_fluid_parameters(
            density=1000.0, particle_radius=r, max_particles=100000  # Water density
        )

        # Set up tank boundaries
        self.setup_tank_boundaries()

        # Initialize particles
        self.initialize_particles()

    def setup_tank_boundaries(self):
        """Mark tank walls as solid cells like in the JS implementation"""
        # Set solid boundaries (tank walls)
        self.sim.set_solid_boundaries()

    def initialize_particles(self):
        """Initialize fluid particles in a hexagonal pattern matching the original JS"""
        # Calculate proportions for dam break setup
        rel_water_width = 0.4
        rel_water_height = 0.8

        # Calculate particle spacing based on cell size
        h = 1.0  # Cell size in simulation units
        r = 0.3 * h  # Particle radius relative to cell size
        dx = 2.0 * r  # Horizontal spacing
        dy = np.sqrt(3.0) / 2.0 * dx  # Vertical spacing for hexagonal packing

        # Calculate grid dimensions in simulation units
        tank_width = self.width
        tank_height = self.height
        water_width = rel_water_width * tank_width
        water_height = rel_water_height * tank_height

        # Calculate number of particles in each dimension
        num_x = int((water_width - 2.0 * h - 2.0 * r) / dx)
        num_y = int((water_height - 2.0 * h - 2.0 * r) / dy)

        print(
            f"Creating {num_x}×{num_y} = {num_x * num_y} particles in hexagonal pattern"
        )

        # Create particles
        for i in range(num_x):
            for j in range(num_y):
                # Calculate particle position with hexagonal offset
                x = h + r + dx * i + (0.0 if j % 2 == 0 else r)
                y = h + r + dy * j
                self.sim.add_particle(x, y)

    def update_simulation(self):
        """Update the simulation state using the current scene parameters"""
        if not self.paused:
            self.frame_nr += 1

            # Pass all scene parameters to the simulation
            self.sim.simulate(
                dt=self.dt,
                gravity=self.gravity,
                flip_ratio=self.flip_ratio,
                num_pressure_iters=self.num_pressure_iters,
                num_particle_iters=self.num_particle_iters,
                over_relaxation=self.over_relaxation,
                compensate_drift=self.compensate_drift,
                separate_particles=self.separate_particles,
                obstacle_x=self.obstacle_x,
                obstacle_y=self.obstacle_y,
                obstacle_radius=self.obstacle_radius,
                obstacle_vel_x=self.obstacle_vel_x,
                obstacle_vel_y=self.obstacle_vel_y,
            )

    def handle_events(self):
        """Process pygame events with expanded controls for simulation parameters"""
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                return False

            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    return False
                elif event.key == pygame.K_SPACE:
                    # Toggle pause
                    self.paused = not self.paused
                elif event.key == pygame.K_o:
                    # Toggle obstacle visibility
                    self.show_obstacle = not self.show_obstacle
                elif event.key == pygame.K_p:
                    # Toggle particle visibility
                    self.show_particles = not self.show_particles
                elif event.key == pygame.K_g:
                    # Toggle grid visibility
                    self.show_grid = not self.show_grid
                elif event.key == pygame.K_c:
                    # Toggle drift compensation
                    self.compensate_drift = not self.compensate_drift
                elif event.key == pygame.K_s:
                    # Toggle particle separation
                    self.separate_particles = not self.separate_particles
                elif event.key == pygame.K_r:
                    # Reset simulation
                    self.setup_simulation()
                    self.frame_nr = 0
                elif event.key == pygame.K_UP:
                    # Increase pressure iterations
                    self.num_pressure_iters += 10
                    print(f"Pressure iterations: {self.num_pressure_iters}")
                elif event.key == pygame.K_DOWN:
                    # Decrease pressure iterations
                    self.num_pressure_iters = max(10, self.num_pressure_iters - 10)
                    print(f"Pressure iterations: {self.num_pressure_iters}")

            # Mouse events for obstacle movement and particle creation
            elif event.type == pygame.MOUSEBUTTONDOWN:
                mouse_x, mouse_y = event.pos
                # Flip the y-coordinate to match the display transformation
                mouse_y = self.screen_height - mouse_y
                cell_x = mouse_x / self.cell_size
                cell_y = mouse_y / self.cell_size

                # Check if click is on obstacle
                dx = cell_x - self.obstacle_x
                dy = cell_y - self.obstacle_y
                if dx * dx + dy * dy <= self.obstacle_radius * self.obstacle_radius:
                    if event.button == 1:  # Left click
                        self.is_moving_obstacle = True
                else:
                    if event.button == 1:  # Left click on fluid
                        self.dragging = True
                        self.create_particles(cell_x, cell_y)
                    elif event.button == 3:  # Right click to adjust obstacle size
                        # Increase/decrease obstacle radius
                        if pygame.key.get_mods() & pygame.KMOD_SHIFT:
                            self.obstacle_radius = max(1.0, self.obstacle_radius - 0.5)
                        else:
                            self.obstacle_radius += 0.5

            elif event.type == pygame.MOUSEMOTION:
                mouse_x, mouse_y = event.pos
                # Flip the y-coordinate to match the display transformation
                mouse_y = self.screen_height - mouse_y
                cell_x = mouse_x / self.cell_size
                cell_y = mouse_y / self.cell_size

                if self.is_moving_obstacle:
                    # Calculate velocity based on movement
                    self.obstacle_vel_x = cell_x - self.obstacle_x
                    self.obstacle_vel_y = cell_y - self.obstacle_y
                    # Update position
                    self.obstacle_x = cell_x
                    self.obstacle_y = cell_y
                elif self.dragging:
                    self.create_particles(cell_x, cell_y)

            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:  # Left click
                    self.dragging = False
                    self.is_moving_obstacle = False
                    # Reset obstacle velocity when released
                    if self.is_moving_obstacle:
                        self.obstacle_vel_x *= 0.5  # Apply some damping
                        self.obstacle_vel_y *= 0.5  # Apply some damping

        return True

    def create_particles(self, center_x, center_y, radius=0.5):
        """Create fluid particles around the center point with denser packing

        Parameters:
            center_x, center_y: Center position of the particle circle
            radius: Radius of the particle circle, defaults to self.ball_radius if None
        """
        # Spacing between particles
        particle_spacing = 0.25

        # Calculate bounds of the box containing the circle
        x_min = max(0, int(center_x - radius))
        x_max = min(self.width, int(center_x + radius + 1))
        y_min = max(0, int(center_y - radius))
        y_max = min(self.height, int(center_y + radius + 1))

        # Create particles within the circle
        for i in range(int((x_max - x_min) / particle_spacing) + 1):
            x = x_min + i * particle_spacing

            for j in range(int((y_max - y_min) / particle_spacing) + 1):
                y = y_min + j * particle_spacing

                # Check if point is within the circle
                distance = ((x - center_x) ** 2 + (y - center_y) ** 2) ** 0.5

                if distance <= radius:
                    self.sim.add_particle(x, y)

    def render(self):
        """Render the current state to the screen"""
        # Clear the screen
        self.screen.fill((0, 0, 0))

        # Create a temporary surface to draw everything on
        temp_surface = pygame.Surface((self.screen_width, self.screen_height))
        temp_surface.fill((0, 0, 0))

        # Render grid if enabled
        if self.show_grid:
            grid_data = self.sim.get_grid_data()
            if grid_data is not None:
                for x in range(self.width):
                    for y in range(self.height):
                        rect = pygame.Rect(
                            x * self.cell_size,
                            y * self.cell_size,
                            self.cell_size,
                            self.cell_size,
                        )
                        # Get cell type and color accordingly
                        cell_type = self.sim.get_cell_type(x, y)
                        if cell_type == 2:  # SOLID_CELL
                            color = (128, 128, 128)  # Gray
                        elif cell_type == 0:  # FLUID_CELL
                            density = self.sim.get_cell_density(x, y)
                            color = self.get_density_color(density)
                        else:  # AIR_CELL
                            color = (0, 0, 0)  # Black

                        pygame.draw.rect(
                            temp_surface, color, rect, 0 if cell_type != 1 else 1
                        )

        # Render particles if enabled
        if self.show_particles:
            particles = self.sim.get_particles()
            for i in range(len(particles) // 5):  # Assuming format [x, y, r, g, b, ...]
                idx = i * 5
                x, y = particles[idx], particles[idx + 1]
                color = (
                    int(particles[idx + 2] * 255),
                    int(particles[idx + 3] * 255),
                    int(particles[idx + 4] * 255),
                )
                pygame.draw.circle(
                    temp_surface,
                    color,
                    (int(x * self.cell_size), int(y * self.cell_size)),
                    max(1, int(self.cell_size * 0.3)),
                )

        # Render obstacle if enabled
        if self.show_obstacle:
            pygame.draw.circle(
                temp_surface,
                (255, 0, 0),  # Red
                (
                    int(self.obstacle_x * self.cell_size),
                    int(self.obstacle_y * self.cell_size),
                ),
                int(self.obstacle_radius * self.cell_size),
                2,  # Line width
            )

        # Flip the temporary surface vertically and blit to the screen
        flipped_surface = pygame.transform.flip(temp_surface, False, True)
        self.screen.blit(flipped_surface, (0, 0))

        # Update the display
        pygame.display.flip()

    def get_density_color(self, density):
        """Map a density value to a color using scientific colormap (blue to cyan to green to yellow to red)"""
        if density <= 0:
            return (0, 0, 255)  # Blue for lowest

        norm_density = min(density / 2.0, 0.999)  # Normalize to [0, 1]

        # Map to color segments
        m = 0.25
        segment = int(norm_density / m)
        t = (norm_density - segment * m) / m  # Position within segment [0,1]

        if segment == 0:  # Blue to Cyan
            return (0, int(255 * t), 255)
        elif segment == 1:  # Cyan to Green
            return (0, 255, int(255 * (1 - t)))
        elif segment == 2:  # Green to Yellow
            return (int(255 * t), 255, 0)
        elif segment == 3:  # Yellow to Red
            return (255, int(255 * (1 - t)), 0)
        else:
            return (255, 0, 0)  # Full red for highest values

    def update_fps_display(self):
        """Update and display the FPS counter"""
        self.frame_count += 1
        if self.frame_count % 30 == 0:
            current_time = time.time()
            fps = self.frame_count / (current_time - self.start_time)
            paused_indicator = " [PAUSED]" if self.paused else ""
            pygame.display.set_caption(
                f"FLIP Simulation - FPS: {fps:.2f}{paused_indicator}"
            )
            self.frame_count = 0
            self.start_time = current_time

    def print_stats(self):
        """Print statistics about the simulation"""
        print(f"Frame number: {self.frame_nr}")
        print(f"Num particles: {self.sim.get_particle_count()}")
        print(
            f"Obstacle: pos=({self.obstacle_x:.1f}, {self.obstacle_y:.1f}), r={self.obstacle_radius:.1f}"
        )
        print(
            f"Settings: dt={self.dt}, flip_ratio={self.flip_ratio}, pressure_iters={self.num_pressure_iters}"
        )
        print(
            f"Features: compensate_drift={self.compensate_drift}, separate_particles={self.separate_particles}"
        )

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
    app = FlipSimulationApp(width=60, height=35, cell_size=35)
    app.run()


if __name__ == "__main__":
    main()
