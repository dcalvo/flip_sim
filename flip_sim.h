#ifndef FLIP_SIM_H
#define FLIP_SIM_H

// Constants
#define U_FIELD 0
#define V_FIELD 1

#define FLUID_CELL 0
#define AIR_CELL 1
#define SOLID_CELL 2

// Data structures
typedef struct FlipFluid
{
    // Fluid grid properties
    float density;
    int fNumX;
    int fNumY;
    float h;           // Cell size
    float fInvSpacing; // 1.0 / h
    int fNumCells;     // fNumX * fNumY

    // Fluid grid arrays
    float *u;         // X velocity component
    float *v;         // Y velocity component
    float *du;        // Delta X velocity
    float *dv;        // Delta Y velocity
    float *prevU;     // Previous X velocity
    float *prevV;     // Previous Y velocity
    float *p;         // Pressure
    float *s;         // Solid boundary indicator
    int *cellType;    // Type of each cell (FLUID_CELL, AIR_CELL, SOLID_CELL)
    float *cellColor; // RGB color values for visualization (3 * fNumCells)

    // Particle properties
    int maxParticles;
    int numParticles;
    float *particlePos;     // Particle positions (2 * maxParticles)
    float *particleColor;   // Particle colors (3 * maxParticles)
    float *particleVel;     // Particle velocities (2 * maxParticles)
    float *particleDensity; // Density at grid cells
    float particleRestDensity;
    float particleRadius;

    // Particle grid properties (for neighbor search)
    float pInvSpacing; // 1.0 / (2.2 * particleRadius)
    int pNumX;
    int pNumY;
    int pNumCells;
    int *numCellParticles;  // Number of particles in each cell
    int *firstCellParticle; // Index of first particle in each cell
    int *cellParticleIds;   // Particle indices sorted by cell
} FlipFluid;

// Utility function prototype
float clamp(float x, float min, float max);

// Core simulation function prototypes
FlipFluid *createFlipFluid(float density, float width, float height, float spacing,
                           float particleRadius, int maxParticles);
void destroyFlipFluid(FlipFluid *fluid);
void simulate(FlipFluid *fluid, float dt, float gravity, float flipRatio,
              int numPressureIters, int numParticleIters, float overRelaxation,
              int compensateDrift, int separateParticles,
              float obstacleX, float obstacleY, float obstacleRadius,
              float obstacleVelX, float obstacleVelY);

// Individual simulation step prototypes
void integrateParticles(FlipFluid *fluid, float dt, float gravity);
void pushParticlesApart(FlipFluid *fluid, int numIters);
void handleParticleCollisions(FlipFluid *fluid, float obstacleX, float obstacleY,
                              float obstacleRadius, float obstacleVelX, float obstacleVelY);
void updateParticleDensity(FlipFluid *fluid);
void transferVelocities(FlipFluid *fluid, int toGrid, float flipRatio);
void solveIncompressibility(FlipFluid *fluid, int numIters, float dt, float overRelaxation, int compensateDrift);
void updateParticleColors(FlipFluid *fluid);
void updateCellColors(FlipFluid *fluid);
void setSciColor(FlipFluid *fluid, int cellNr, float val, float minVal, float maxVal);

// Particle management functions
int add_particle(FlipFluid *fluid, float x, float y);
int get_particle_count(FlipFluid *fluid);
float *get_particles(FlipFluid *fluid);

// Cell information functions
int get_cell_type(FlipFluid *fluid, int x, int y);
float get_cell_density(FlipFluid *fluid, int x, int y);
void set_solid_boundaries(FlipFluid *fluid);

#endif /* FLIP_SIM_H */