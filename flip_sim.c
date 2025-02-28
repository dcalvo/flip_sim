// flip_sim.c
#include <stdio.h>
#include <stdlib.h>

// Grid structure
typedef struct
{
    int width;
    int height;
    float *cells; // 1D array representing 2D grid for simplicity
} Grid;

// Global grid pointer (for this simple example)
static Grid *g_grid = NULL;

// Create a grid with specified dimensions
int create_grid(int width, int height)
{
    // Free existing grid if any
    if (g_grid != NULL)
    {
        free(g_grid->cells);
        free(g_grid);
    }

    // Allocate new grid
    g_grid = (Grid *)malloc(sizeof(Grid));
    if (!g_grid)
        return 0; // Allocation failed

    g_grid->width = width;
    g_grid->height = height;
    g_grid->cells = (float *)calloc(width * height, sizeof(float));

    if (!g_grid->cells)
    {
        free(g_grid);
        g_grid = NULL;
        return 0; // Allocation failed
    }

    printf("Created grid: %d x %d\n", width, height);
    return 1; // Success
}

// Set a value in the grid
int set_cell(int x, int y, float value)
{
    if (!g_grid)
        return 0; // Grid not initialized
    if (x < 0 || x >= g_grid->width || y < 0 || y >= g_grid->height)
        return 0; // Out of bounds

    g_grid->cells[y * g_grid->width + x] = value;
    return 1; // Success
}

// Get a value from the grid
float get_cell(int x, int y)
{
    if (!g_grid)
        return 0.0f; // Grid not initialized
    if (x < 0 || x >= g_grid->width || y < 0 || y >= g_grid->height)
        return 0.0f; // Out of bounds

    return g_grid->cells[y * g_grid->width + x];
}

// Get the entire grid data
float *get_grid_data()
{
    if (!g_grid)
        return NULL;
    return g_grid->cells;
}

// Get grid dimensions
void get_grid_dimensions(int *width, int *height)
{
    if (!g_grid)
    {
        *width = 0;
        *height = 0;
        return;
    }
    *width = g_grid->width;
    *height = g_grid->height;
}

// Clean up resources
void destroy_grid()
{
    if (g_grid)
    {
        free(g_grid->cells);
        free(g_grid);
        g_grid = NULL;
        printf("Grid destroyed\n");
    }
}