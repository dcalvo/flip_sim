#include "flip_sim.h"
#include <stdlib.h> // For malloc, free
#include <math.h>   // For sqrt, floor
#include <string.h> // For memset, memcpy

// Utility function implementation
float clamp(float x, float min, float max)
{
    if (x < min)
        return min;
    else if (x > max)
        return max;
    else
        return x;
}

FlipFluid *createFlipFluid(float density, float width, float height, float spacing,
                           float particleRadius, int maxParticles)
{
    // Allocate the fluid structure
    FlipFluid *fluid = (FlipFluid *)malloc(sizeof(FlipFluid));
    if (!fluid)
        return NULL;

    // Initialize grid dimensions
    fluid->density = density;
    fluid->fNumX = (int)(width / spacing) + 1;
    fluid->fNumY = (int)(height / spacing) + 1;
    fluid->h = fmaxf(width / fluid->fNumX, height / fluid->fNumY);
    fluid->fInvSpacing = 1.0f / fluid->h;
    fluid->fNumCells = fluid->fNumX * fluid->fNumY;

    // Allocate grid arrays
    fluid->u = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->v = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->du = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->dv = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->prevU = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->prevV = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->p = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->s = (float *)calloc(fluid->fNumCells, sizeof(float));
    fluid->cellType = (int *)calloc(fluid->fNumCells, sizeof(int));
    fluid->cellColor = (float *)calloc(fluid->fNumCells * 3, sizeof(float));
    fluid->particleDensity = (float *)calloc(fluid->fNumCells, sizeof(float));

    // Check allocations
    if (!fluid->u || !fluid->v || !fluid->du || !fluid->dv ||
        !fluid->prevU || !fluid->prevV || !fluid->p || !fluid->s ||
        !fluid->cellType || !fluid->cellColor || !fluid->particleDensity)
    {
        destroyFlipFluid(fluid);
        return NULL;
    }

    // Initialize particle properties
    fluid->maxParticles = maxParticles;
    fluid->numParticles = 0;
    fluid->particleRadius = particleRadius;
    fluid->particleRestDensity = 0.0f;

    // Allocate particle arrays
    fluid->particlePos = (float *)calloc(maxParticles * 2, sizeof(float));
    fluid->particleColor = (float *)calloc(maxParticles * 3, sizeof(float));
    fluid->particleVel = (float *)calloc(maxParticles * 2, sizeof(float));

    // Check allocations
    if (!fluid->particlePos || !fluid->particleColor || !fluid->particleVel)
    {
        destroyFlipFluid(fluid);
        return NULL;
    }

    // Initialize particle color (setting blue component to 1.0)
    for (int i = 0; i < maxParticles; i++)
    {
        fluid->particleColor[3 * i + 2] = 1.0f;
    }

    // Set up particle neighborhood search grid
    fluid->pInvSpacing = 1.0f / (2.2f * particleRadius);
    fluid->pNumX = (int)(width * fluid->pInvSpacing) + 1;
    fluid->pNumY = (int)(height * fluid->pInvSpacing) + 1;
    fluid->pNumCells = fluid->pNumX * fluid->pNumY;

    // Allocate particle grid arrays
    fluid->numCellParticles = (int *)calloc(fluid->pNumCells, sizeof(int));
    fluid->firstCellParticle = (int *)calloc(fluid->pNumCells + 1, sizeof(int));
    fluid->cellParticleIds = (int *)calloc(maxParticles, sizeof(int));

    // Check allocations
    if (!fluid->numCellParticles || !fluid->firstCellParticle || !fluid->cellParticleIds)
    {
        destroyFlipFluid(fluid);
        return NULL;
    }

    return fluid;
}

void destroyFlipFluid(FlipFluid *fluid)
{
    if (!fluid)
        return;

    // Free grid arrays
    free(fluid->u);
    free(fluid->v);
    free(fluid->du);
    free(fluid->dv);
    free(fluid->prevU);
    free(fluid->prevV);
    free(fluid->p);
    free(fluid->s);
    free(fluid->cellType);
    free(fluid->cellColor);
    free(fluid->particleDensity);

    // Free particle arrays
    free(fluid->particlePos);
    free(fluid->particleColor);
    free(fluid->particleVel);

    // Free particle grid arrays
    free(fluid->numCellParticles);
    free(fluid->firstCellParticle);
    free(fluid->cellParticleIds);

    // Free the fluid structure itself
    free(fluid);
}

void simulate(FlipFluid *fluid, float dt, float gravity, float flipRatio,
              int numPressureIters, int numParticleIters, float overRelaxation,
              int compensateDrift, int separateParticles,
              float obstacleX, float obstacleY, float obstacleRadius,
              float obstacleVelX, float obstacleVelY) // Added parameters
{
    int numSubSteps = 1;
    float sdt = dt / numSubSteps;

    // Perform simulation sub-steps
    for (int step = 0; step < numSubSteps; step++)
    {
        // Step 1: Integrate particles with gravity
        integrateParticles(fluid, sdt, gravity);

        // Step 2: Optional particle separation for better distribution
        if (separateParticles)
        {
            pushParticlesApart(fluid, numParticleIters);
        }

        // Step 3: Handle collisions with boundaries and obstacle
        handleParticleCollisions(fluid, obstacleX, obstacleY, obstacleRadius,
                                 obstacleVelX, obstacleVelY); // Pass obstacle velocities

        // Step 4: Transfer particle velocities to grid
        transferVelocities(fluid, 1, 0.0f); // toGrid = true

        // Step 5: Update particle density field
        updateParticleDensity(fluid);

        // Step 6: Solve pressure and enforce incompressibility
        solveIncompressibility(fluid, numPressureIters, sdt, overRelaxation, compensateDrift);

        // Step 7: Transfer grid velocities back to particles
        transferVelocities(fluid, 0, flipRatio); // toGrid = false
    }

    // Update visual representation
    updateParticleColors(fluid);
    updateCellColors(fluid);
}

void integrateParticles(FlipFluid *fluid, float dt, float gravity)
{
    // Safety check
    if (!fluid)
        return;

    for (int i = 0; i < fluid->numParticles; i++)
    {
        // Update vertical velocity component with gravity
        fluid->particleVel[2 * i + 1] += dt * gravity;

        // Update position based on velocity (explicit Euler integration)
        fluid->particlePos[2 * i] += fluid->particleVel[2 * i] * dt;
        fluid->particlePos[2 * i + 1] += fluid->particleVel[2 * i + 1] * dt;
    }
}

void pushParticlesApart(FlipFluid *fluid, int numIters)
{
    if (!fluid)
        return;

    float colorDiffusionCoeff = 0.001f;

    // Count particles per cell
    memset(fluid->numCellParticles, 0, fluid->pNumCells * sizeof(int));

    for (int i = 0; i < fluid->numParticles; i++)
    {
        float x = fluid->particlePos[2 * i];
        float y = fluid->particlePos[2 * i + 1];

        int xi = (int)clamp(floorf(x * fluid->pInvSpacing), 0, fluid->pNumX - 1);
        int yi = (int)clamp(floorf(y * fluid->pInvSpacing), 0, fluid->pNumY - 1);
        int cellNr = xi * fluid->pNumY + yi;
        fluid->numCellParticles[cellNr]++;
    }

    // Calculate partial sums
    int first = 0;
    for (int i = 0; i < fluid->pNumCells; i++)
    {
        first += fluid->numCellParticles[i];
        fluid->firstCellParticle[i] = first;
    }
    fluid->firstCellParticle[fluid->pNumCells] = first; // Guard

    // Fill particles into cells
    for (int i = 0; i < fluid->numParticles; i++)
    {
        float x = fluid->particlePos[2 * i];
        float y = fluid->particlePos[2 * i + 1];

        int xi = (int)clamp(floorf(x * fluid->pInvSpacing), 0, fluid->pNumX - 1);
        int yi = (int)clamp(floorf(y * fluid->pInvSpacing), 0, fluid->pNumY - 1);
        int cellNr = xi * fluid->pNumY + yi;
        fluid->firstCellParticle[cellNr]--;
        fluid->cellParticleIds[fluid->firstCellParticle[cellNr]] = i;
    }

    // Push particles apart
    float minDist = 2.0f * fluid->particleRadius;
    float minDist2 = minDist * minDist;

    for (int iter = 0; iter < numIters; iter++)
    {
        for (int i = 0; i < fluid->numParticles; i++)
        {
            float px = fluid->particlePos[2 * i];
            float py = fluid->particlePos[2 * i + 1];

            int pxi = (int)floorf(px * fluid->pInvSpacing);
            int pyi = (int)floorf(py * fluid->pInvSpacing);
            int x0 = fmaxf(pxi - 1, 0);
            int y0 = fmaxf(pyi - 1, 0);
            int x1 = fminf(pxi + 1, fluid->pNumX - 1);
            int y1 = fminf(pyi + 1, fluid->pNumY - 1);

            for (int xi = x0; xi <= x1; xi++)
            {
                for (int yi = y0; yi <= y1; yi++)
                {
                    int cellNr = xi * fluid->pNumY + yi;
                    int first = fluid->firstCellParticle[cellNr];
                    int last = fluid->firstCellParticle[cellNr + 1];

                    for (int j = first; j < last; j++)
                    {
                        int id = fluid->cellParticleIds[j];
                        if (id == i)
                            continue;

                        float qx = fluid->particlePos[2 * id];
                        float qy = fluid->particlePos[2 * id + 1];

                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;

                        if (d2 > minDist2 || d2 == 0.0f)
                            continue;

                        float d = sqrtf(d2);
                        float s = 0.5f * (minDist - d) / d;
                        dx *= s;
                        dy *= s;

                        fluid->particlePos[2 * i] -= dx;
                        fluid->particlePos[2 * i + 1] -= dy;
                        fluid->particlePos[2 * id] += dx;
                        fluid->particlePos[2 * id + 1] += dy;

                        // Diffuse colors
                        for (int k = 0; k < 3; k++)
                        {
                            float color0 = fluid->particleColor[3 * i + k];
                            float color1 = fluid->particleColor[3 * id + k];
                            float color = (color0 + color1) * 0.5f;
                            fluid->particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
                            fluid->particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
                        }
                    }
                }
            }
        }
    }
}

void handleParticleCollisions(FlipFluid *fluid, float obstacleX, float obstacleY,
                              float obstacleRadius, float obstacleVelX, float obstacleVelY)
{
    if (!fluid)
        return;

    float h = 1.0f / fluid->fInvSpacing;
    float r = fluid->particleRadius;
    float or = obstacleRadius;
    float or2 = or * or ;
    float minDist = obstacleRadius + r;
    float minDist2 = minDist * minDist;

    float minX = h + r;
    float maxX = (fluid->fNumX - 1) * h - r;
    float minY = h + r;
    float maxY = (fluid->fNumY - 1) * h - r;

    for (int i = 0; i < fluid->numParticles; i++)
    {
        float x = fluid->particlePos[2 * i];
        float y = fluid->particlePos[2 * i + 1];

        // Obstacle collision
        float dx = x - obstacleX;
        float dy = y - obstacleY;
        float d2 = dx * dx + dy * dy;

        if (d2 < minDist2)
        {
            // Note: The original code has a commented section for position adjustment
            // that could be uncommented if needed:
            //
            // float d = sqrtf(d2);
            // float s = (minDist - d) / d;
            // x += dx * s;
            // y += dy * s;

            // Set particle velocity to match obstacle velocity
            fluid->particleVel[2 * i] = obstacleVelX;
            fluid->particleVel[2 * i + 1] = obstacleVelY;
        }

        // Wall collisions
        if (x < minX)
        {
            x = minX;
            fluid->particleVel[2 * i] = 0.0f;
        }
        if (x > maxX)
        {
            x = maxX;
            fluid->particleVel[2 * i] = 0.0f;
        }
        if (y < minY)
        {
            y = minY;
            fluid->particleVel[2 * i + 1] = 0.0f;
        }
        if (y > maxY)
        {
            y = maxY;
            fluid->particleVel[2 * i + 1] = 0.0f;
        }

        // Update particle position after collision handling
        fluid->particlePos[2 * i] = x;
        fluid->particlePos[2 * i + 1] = y;
    }
}

void updateParticleDensity(FlipFluid *fluid)
{
    if (!fluid)
        return;

    int n = fluid->fNumY;
    float h = fluid->h;
    float h1 = fluid->fInvSpacing;
    float h2 = 0.5f * h;

    // Clear the density array
    memset(fluid->particleDensity, 0, fluid->fNumCells * sizeof(float));

    for (int i = 0; i < fluid->numParticles; i++)
    {
        float x = fluid->particlePos[2 * i];
        float y = fluid->particlePos[2 * i + 1];

        // Clamp particle position to be within grid bounds
        x = clamp(x, h, (fluid->fNumX - 1) * h);
        y = clamp(y, h, (fluid->fNumY - 1) * h);

        // Calculate bilinear interpolation indices and weights
        int x0 = (int)floorf((x - h2) * h1);
        float tx = ((x - h2) - x0 * h) * h1;
        int x1 = (int)fminf(x0 + 1, fluid->fNumX - 2);

        int y0 = (int)floorf((y - h2) * h1);
        float ty = ((y - h2) - y0 * h) * h1;
        int y1 = (int)fminf(y0 + 1, fluid->fNumY - 2);

        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        // Distribute density to surrounding cells using bilinear weights
        if (x0 < fluid->fNumX && y0 < fluid->fNumY)
            fluid->particleDensity[x0 * n + y0] += sx * sy;
        if (x1 < fluid->fNumX && y0 < fluid->fNumY)
            fluid->particleDensity[x1 * n + y0] += tx * sy;
        if (x1 < fluid->fNumX && y1 < fluid->fNumY)
            fluid->particleDensity[x1 * n + y1] += tx * ty;
        if (x0 < fluid->fNumX && y1 < fluid->fNumY)
            fluid->particleDensity[x0 * n + y1] += sx * ty;
    }

    // If rest density hasn't been calculated yet, calculate it
    if (fluid->particleRestDensity == 0.0f)
    {
        float sum = 0.0f;
        int numFluidCells = 0;

        for (int i = 0; i < fluid->fNumCells; i++)
        {
            if (fluid->cellType[i] == FLUID_CELL)
            {
                sum += fluid->particleDensity[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            fluid->particleRestDensity = sum / numFluidCells;
    }

    /* Commented out as in the original code
    for (int xi = 1; xi < fluid->fNumX; xi++) {
        for (int yi = 1; yi < fluid->fNumY; yi++) {
            int cellNr = xi * n + yi;
            if (fluid->cellType[cellNr] != FLUID_CELL)
                continue;
            float hx = fluid->h;
            float hy = fluid->h;

            if (fluid->cellType[(xi - 1) * n + yi] == SOLID_CELL ||
                fluid->cellType[(xi + 1) * n + yi] == SOLID_CELL)
                hx -= fluid->particleRadius;
            if (fluid->cellType[xi * n + yi - 1] == SOLID_CELL ||
                fluid->cellType[xi * n + yi + 1] == SOLID_CELL)
                hy -= fluid->particleRadius;

            float scale = fluid->h * fluid->h / (hx * hy);
            fluid->particleDensity[cellNr] *= scale;
        }
    }
    */
}

void transferVelocities(FlipFluid *fluid, int toGrid, float flipRatio)
{
    if (!fluid)
        return;

    int n = fluid->fNumY;
    float h = fluid->h;
    float h1 = fluid->fInvSpacing;
    float h2 = 0.5f * h;

    if (toGrid)
    {
        // Save previous velocities
        memcpy(fluid->prevU, fluid->u, fluid->fNumCells * sizeof(float));
        memcpy(fluid->prevV, fluid->v, fluid->fNumCells * sizeof(float));

        // Clear velocity fields and weights
        memset(fluid->du, 0, fluid->fNumCells * sizeof(float));
        memset(fluid->dv, 0, fluid->fNumCells * sizeof(float));
        memset(fluid->u, 0, fluid->fNumCells * sizeof(float));
        memset(fluid->v, 0, fluid->fNumCells * sizeof(float));

        // Mark cells as solid or air based on s values
        for (int i = 0; i < fluid->fNumCells; i++)
            fluid->cellType[i] = fluid->s[i] == 0.0f ? SOLID_CELL : AIR_CELL;

        // Mark cells containing particles as fluid
        for (int i = 0; i < fluid->numParticles; i++)
        {
            float x = fluid->particlePos[2 * i];
            float y = fluid->particlePos[2 * i + 1];
            int xi = (int)clamp(floorf(x * h1), 0, fluid->fNumX - 1);
            int yi = (int)clamp(floorf(y * h1), 0, fluid->fNumY - 1);
            int cellNr = xi * n + yi;
            if (fluid->cellType[cellNr] == AIR_CELL)
                fluid->cellType[cellNr] = FLUID_CELL;
        }
    }

    // Process both velocity components (u and v)
    for (int component = 0; component < 2; component++)
    {
        float dx = component == 0 ? 0.0f : h2;
        float dy = component == 0 ? h2 : 0.0f;

        float *f = component == 0 ? fluid->u : fluid->v;
        float *prevF = component == 0 ? fluid->prevU : fluid->prevV;
        float *d = component == 0 ? fluid->du : fluid->dv;

        // Process all particles
        for (int i = 0; i < fluid->numParticles; i++)
        {
            float x = fluid->particlePos[2 * i];
            float y = fluid->particlePos[2 * i + 1];

            // Clamp particle position to grid bounds
            x = clamp(x, h, (fluid->fNumX - 1) * h);
            y = clamp(y, h, (fluid->fNumY - 1) * h);

            // Calculate bilinear interpolation weights
            int x0 = (int)fminf(floorf((x - dx) * h1), fluid->fNumX - 2);
            float tx = ((x - dx) - x0 * h) * h1;
            int x1 = (int)fminf(x0 + 1, fluid->fNumX - 2);

            int y0 = (int)fminf(floorf((y - dy) * h1), fluid->fNumY - 2);
            float ty = ((y - dy) - y0 * h) * h1;
            int y1 = (int)fminf(y0 + 1, fluid->fNumY - 2);

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            float d0 = sx * sy;
            float d1 = tx * sy;
            float d2 = tx * ty;
            float d3 = sx * ty;

            int nr0 = x0 * n + y0;
            int nr1 = x1 * n + y0;
            int nr2 = x1 * n + y1;
            int nr3 = x0 * n + y1;

            if (toGrid)
            {
                // Transfer particle velocity to grid
                float pv = fluid->particleVel[2 * i + component];
                f[nr0] += pv * d0;
                d[nr0] += d0;
                f[nr1] += pv * d1;
                d[nr1] += d1;
                f[nr2] += pv * d2;
                d[nr2] += d2;
                f[nr3] += pv * d3;
                d[nr3] += d3;
            }
            else
            {
                // Transfer grid velocity to particles using FLIP/PIC blend
                int offset = component == 0 ? n : 1;
                float valid0 = (fluid->cellType[nr0] == FLUID_CELL ||
                                fluid->cellType[nr0 - offset] == FLUID_CELL)
                                   ? 1.0f
                                   : 0.0f;
                float valid1 = (fluid->cellType[nr1] == FLUID_CELL ||
                                fluid->cellType[nr1 - offset] == FLUID_CELL)
                                   ? 1.0f
                                   : 0.0f;
                float valid2 = (fluid->cellType[nr2] == FLUID_CELL ||
                                fluid->cellType[nr2 - offset] == FLUID_CELL)
                                   ? 1.0f
                                   : 0.0f;
                float valid3 = (fluid->cellType[nr3] == FLUID_CELL ||
                                fluid->cellType[nr3 - offset] == FLUID_CELL)
                                   ? 1.0f
                                   : 0.0f;

                float v = fluid->particleVel[2 * i + component];
                float totalWeight = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                if (totalWeight > 0.0f)
                {
                    // PIC component (interpolated grid velocity)
                    float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] +
                                  valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) /
                                 totalWeight;

                    // FLIP component (current velocity plus grid velocity change)
                    float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) +
                                  valid1 * d1 * (f[nr1] - prevF[nr1]) +
                                  valid2 * d2 * (f[nr2] - prevF[nr2]) +
                                  valid3 * d3 * (f[nr3] - prevF[nr3])) /
                                 totalWeight;
                    float flipV = v + corr;

                    // Blend PIC and FLIP
                    fluid->particleVel[2 * i + component] = (1.0f - flipRatio) * picV + flipRatio * flipV;
                }
            }
        }

        if (toGrid)
        {
            // Normalize grid velocities
            for (int i = 0; i < fluid->fNumCells; i++)
            {
                if (d[i] > 0.0f)
                    f[i] /= d[i];
            }

            // Restore velocities at solid boundaries
            for (int i = 0; i < fluid->fNumX; i++)
            {
                for (int j = 0; j < fluid->fNumY; j++)
                {
                    int solid = fluid->cellType[i * n + j] == SOLID_CELL;
                    if (solid || (i > 0 && fluid->cellType[(i - 1) * n + j] == SOLID_CELL))
                        fluid->u[i * n + j] = fluid->prevU[i * n + j];
                    if (solid || (j > 0 && fluid->cellType[i * n + j - 1] == SOLID_CELL))
                        fluid->v[i * n + j] = fluid->prevV[i * n + j];
                }
            }
        }
    }
}

void solveIncompressibility(FlipFluid *fluid, int numIters, float dt, float overRelaxation, int compensateDrift)
{
    // Initialize pressure to 0
    memset(fluid->p, 0, fluid->fNumCells * sizeof(float));

    // Copy current velocities to previous velocities
    memcpy(fluid->prevU, fluid->u, fluid->fNumCells * sizeof(float));
    memcpy(fluid->prevV, fluid->v, fluid->fNumCells * sizeof(float));

    int n = fluid->fNumY;
    float cp = fluid->density * fluid->h / dt;

    // Note: The loop below is in the original JavaScript code but doesn't do anything useful
    // It's included here for completeness but could be removed
    for (int i = 0; i < fluid->fNumCells; i++)
    {
        float u = fluid->u[i];
        float v = fluid->v[i];
    }

    // Main solver loop
    for (int iter = 0; iter < numIters; iter++)
    {
        for (int i = 1; i < fluid->fNumX - 1; i++)
        {
            for (int j = 1; j < fluid->fNumY - 1; j++)
            {
                int center = i * n + j;

                // Skip non-fluid cells
                if (fluid->cellType[center] != FLUID_CELL)
                    continue;

                int left = (i - 1) * n + j;
                int right = (i + 1) * n + j;
                int bottom = i * n + j - 1;
                int top = i * n + j + 1;

                // In the original code, this line exists but the value is not used
                // float s_center = fluid->s[center];

                // Calculate weights from solid fraction values
                float sx0 = fluid->s[left];
                float sx1 = fluid->s[right];
                float sy0 = fluid->s[bottom];
                float sy1 = fluid->s[top];
                float s = sx0 + sx1 + sy0 + sy1;

                if (s == 0.0f)
                    continue;

                // Calculate divergence
                float div = fluid->u[right] - fluid->u[center] +
                            fluid->v[top] - fluid->v[center];

                // Apply density drift compensation if enabled
                if (fluid->particleRestDensity > 0.0f && compensateDrift)
                {
                    float k = 1.0f;
                    float compression = fluid->particleDensity[center] - fluid->particleRestDensity;
                    if (compression > 0.0f)
                        div = div - k * compression;
                }

                // Calculate pressure adjustment
                float p = -div / s;
                p *= overRelaxation;
                fluid->p[center] += cp * p;

                // Apply pressure to velocities
                fluid->u[center] -= sx0 * p;
                fluid->u[right] += sx1 * p;
                fluid->v[center] -= sy0 * p;
                fluid->v[top] += sy1 * p;
            }
        }
    }
}

void updateParticleColors(FlipFluid *fluid)
{
    float h1 = fluid->fInvSpacing;

    for (int i = 0; i < fluid->numParticles; i++)
    {
        // Color change rate
        float s = 0.01f;

        // Gradually shift colors (decrease red and green, increase blue)
        fluid->particleColor[3 * i] = clamp(fluid->particleColor[3 * i] - s, 0.0f, 1.0f);
        fluid->particleColor[3 * i + 1] = clamp(fluid->particleColor[3 * i + 1] - s, 0.0f, 1.0f);
        fluid->particleColor[3 * i + 2] = clamp(fluid->particleColor[3 * i + 2] + s, 0.0f, 1.0f);

        // Get particle position
        float x = fluid->particlePos[2 * i];
        float y = fluid->particlePos[2 * i + 1];

        // Determine cell that contains this particle
        int xi = (int)clamp(floorf(x * h1), 1, fluid->fNumX - 1);
        int yi = (int)clamp(floorf(y * h1), 1, fluid->fNumY - 1);
        int cellNr = xi * fluid->fNumY + yi;

        // Check for low density regions
        float d0 = fluid->particleRestDensity;

        if (d0 > 0.0f)
        {
            float relDensity = fluid->particleDensity[cellNr] / d0;
            if (relDensity < 0.7f)
            {
                // Special coloring for low-density regions (light blue)
                float s = 0.8f;
                fluid->particleColor[3 * i] = s;
                fluid->particleColor[3 * i + 1] = s;
                fluid->particleColor[3 * i + 2] = 1.0f;
            }
        }
    }
}

void setSciColor(FlipFluid *fluid, int cellNr, float val, float minVal, float maxVal)
{
    // Clamp the value to the given range
    val = fminf(fmaxf(val, minVal), maxVal - 0.0001f);

    // Normalize to [0,1]
    float d = maxVal - minVal;
    val = (d == 0.0f) ? 0.5f : (val - minVal) / d;

    // Find color segment and interpolation factor
    float m = 0.25f;
    int num = (int)(val / m);
    float s = (val - num * m) / m;

    float r, g, b;

    // Set RGB based on segment
    switch (num)
    {
    case 0:
        r = 0.0f;
        g = s;
        b = 1.0f;
        break; // Blue to Cyan
    case 1:
        r = 0.0f;
        g = 1.0f;
        b = 1.0f - s;
        break; // Cyan to Green
    case 2:
        r = s;
        g = 1.0f;
        b = 0.0f;
        break; // Green to Yellow
    case 3:
        r = 1.0f;
        g = 1.0f - s;
        b = 0.0f;
        break; // Yellow to Red
    default:
        r = 1.0f;
        g = 0.0f;
        b = 0.0f;
        break; // Default Red
    }

    // Assign color to cell
    fluid->cellColor[3 * cellNr] = r;
    fluid->cellColor[3 * cellNr + 1] = g;
    fluid->cellColor[3 * cellNr + 2] = b;
}

void updateCellColors(FlipFluid *fluid)
{
    // Clear all colors to black
    memset(fluid->cellColor, 0, 3 * fluid->fNumCells * sizeof(float));

    // Set colors based on cell type
    for (int i = 0; i < fluid->fNumCells; i++)
    {
        if (fluid->cellType[i] == SOLID_CELL)
        {
            // Solid cells are gray
            fluid->cellColor[3 * i] = 0.5f;
            fluid->cellColor[3 * i + 1] = 0.5f;
            fluid->cellColor[3 * i + 2] = 0.5f;
        }
        else if (fluid->cellType[i] == FLUID_CELL)
        {
            // Fluid cells based on density
            float d = fluid->particleDensity[i];
            if (fluid->particleRestDensity > 0.0f)
            {
                d /= fluid->particleRestDensity;
            }
            setSciColor(fluid, i, d, 0.0f, 2.0f);
        }
    }
}

// Create a particle at the specified position
int add_particle(FlipFluid *fluid, float x, float y)
{
    // Safety checks
    if (!fluid || fluid->numParticles >= fluid->maxParticles)
        return 0;

    // Set position and default velocity for the new particle
    int id = fluid->numParticles;
    fluid->particlePos[2 * id] = x;
    fluid->particlePos[2 * id + 1] = y;
    fluid->particleVel[2 * id] = 0.0f;
    fluid->particleVel[2 * id + 1] = 0.0f;

    // Default particle coloring (blue)
    fluid->particleColor[3 * id] = 0.0f;     // Red
    fluid->particleColor[3 * id + 1] = 0.0f; // Green
    fluid->particleColor[3 * id + 2] = 1.0f; // Blue

    // Increment the particle count
    fluid->numParticles++;
    return 1;
}

// Get the number of particles
int get_particle_count(FlipFluid *fluid)
{
    if (!fluid)
        return 0;

    return fluid->numParticles;
}

// Get particle data (returns a pointer to an array of floats: x, y, r, g, b for each particle)
float *get_particles(FlipFluid *fluid)
{
    if (!fluid || fluid->numParticles == 0)
        return NULL;

    // This function returns a pointer to a static buffer to avoid memory management issues
    // The buffer will contain [x, y, r, g, b] for each particle
    static float *particleData = NULL;
    static int lastSize = 0;

    // Allocate or reallocate the buffer if needed
    int requiredSize = fluid->numParticles * 5;
    if (!particleData || lastSize < requiredSize)
    {
        if (particleData)
            free(particleData);

        particleData = (float *)malloc(requiredSize * sizeof(float));
        if (!particleData)
            return NULL;

        lastSize = requiredSize;
    }

    // Fill the buffer with particle data
    for (int i = 0; i < fluid->numParticles; i++)
    {
        // Position
        particleData[5 * i] = fluid->particlePos[2 * i];         // x
        particleData[5 * i + 1] = fluid->particlePos[2 * i + 1]; // y

        // Color
        particleData[5 * i + 2] = fluid->particleColor[3 * i];     // r
        particleData[5 * i + 3] = fluid->particleColor[3 * i + 1]; // g
        particleData[5 * i + 4] = fluid->particleColor[3 * i + 2]; // b
    }

    return particleData;
}

// Get the type of cell at (x,y)
int get_cell_type(FlipFluid *fluid, int x, int y)
{
    if (!fluid || x < 0 || x >= fluid->fNumX || y < 0 || y >= fluid->fNumY)
        return -1; // Invalid coordinates

    int cellNr = x * fluid->fNumY + y;
    return fluid->cellType[cellNr];
}

// Get the fluid density at cell (x,y)
float get_cell_density(FlipFluid *fluid, int x, int y)
{
    if (!fluid || x < 0 || x >= fluid->fNumX || y < 0 || y >= fluid->fNumY)
        return 0.0f; // Invalid coordinates

    int cellNr = x * fluid->fNumY + y;
    return fluid->particleDensity[cellNr];
}

// Set solid boundaries (tank walls)
void set_solid_boundaries(FlipFluid *fluid)
{
    if (!fluid)
        return;

    int n = fluid->fNumY;

    // Set up tank walls as solid
    for (int i = 0; i < fluid->fNumX; i++)
    {
        for (int j = 0; j < fluid->fNumY; j++)
        {
            float s = 1.0f; // Default to fluid

            // Mark left, right and bottom walls as solid
            if (i == 0 || i == fluid->fNumX - 1 || j == 0)
                s = 0.0f; // Solid

            fluid->s[i * n + j] = s;
        }
    }
}