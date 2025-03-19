#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>

/* AUTHOR : Charles Bouillaguet <charles.bouillaguet@lip6.fr>
   USAGE  : compile with -lm (and why not -O3)
            redirect the standard output to a text file
            gcc heatsink.c -O3 -lm -o heatsink
            ./heatsink > steady_state.txt
            then run the indicated python script for graphical rendering

   DISCLAIMER : this code does not claim to an absolute realism.
                this code could be obviously improved, but has been written so as
				to make as clear as possible the physics principle of the simulation.
*/

/* EDITOR : Okba kharef */

/* one can change the matter of the heatsink, its size, the power of the CPU, etc. */
#define ALUMINIUM
#define FAST         /* MEDIUM is faster, and FAST is even faster (for debugging) */
#define DUMP_STEADY_STATE

const double L = 0.15;      /* length (x) of the heatsink (m) */
const double l = 0.12;      /* height (y) of the heatsink (m) */
const double E = 0.008;     /* width (z) of the heatsink (m) */
const double watercooling_T = 20;   /* temperature of the fluid for water-cooling, (°C) */
const double CPU_TDP = 280; /* power dissipated by the CPU (W) */

/* dl: "spatial step" for simulation (m) */
/* dt: "time step" for simulation (s) */
#ifdef FAST
double dl = 0.004;
double dt = 0.004;
#endif

#ifdef MEDIUM
double dl = 0.002;
double dt = 0.002;
#endif

#ifdef NORMAL
double dl = 0.001;
double dt = 0.001;
#endif

#ifdef CHALLENGE
double dl = 0.0001;
double dt = 0.00001;
#endif

/* sink_heat_capacity: specific heat capacity of the heatsink (J / kg / K) */
/* sink_density: density of the heatsink (kg / m^3) */
/* sink_conductivity: thermal conductivity of the heatsink (W / m / K) */
/* euros_per_kg: price of the matter by kilogram */
#ifdef ALUMINIUM
double sink_heat_capacity = 897;
double sink_density = 2710;
double sink_conductivity = 237;
double euros_per_kg = 1.594;
#endif

#ifdef COPPER
double sink_heat_capacity = 385;
double sink_density = 8960;
double sink_conductivity = 390;
double euros_per_kg = 5.469;
#endif

#ifdef GOLD
double sink_heat_capacity = 128;
double sink_density = 19300;
double sink_conductivity = 317;
double euros_per_kg = 47000;
#endif

#ifdef IRON
double sink_heat_capacity = 444;
double sink_density = 7860;
double sink_conductivity = 80;
double euros_per_kg = 0.083;
#endif

const double Stefan_Boltzmann = 5.6703e-8;  /* (W / m^2 / K^4), radiation of black body */
const double heat_transfer_coefficient = 10;    /* coefficient of thermal convection (W / m^2 / K) */
double CPU_surface;

/*
 * Return True if the CPU is in contact with the heatsink at the point (x,y).
 * This describes an AMD EPYC "Rome".
 */
static inline bool CPU_shape(double x, double y)
{
    x -= (L - 0.0754) / 2;
    y -= (l - 0.0585) / 2;
    bool small_y_ok = (y > 0.015 && y < 0.025) || (y > 0.0337 && y < 0.0437);
    bool small_x_ok = (x > 0.0113 && x < 0.0186) || (x > 0.0193 && x < 0.0266)
        || (x > 0.0485 && x < 0.0558) || (x > 0.0566 && x < 0.0639);
    bool big_ok = (x > 0.03 && x < 0.045 && y > 0.0155 && y < 0.0435);
    return big_ok || (small_x_ok && small_y_ok);
}

/* returns the total area of the surface of contact between CPU and heatsink (in m^2) */
double CPU_contact_surface()
{
    double S = 0;
    for (double x = dl / 2; x < L; x += dl)
        for (double y = dl / 2; y < l; y += dl)
            if (CPU_shape(x, y))
                S += dl * dl;
    return S;
}

/* Returns the new temperature of the cell (i, j, k). For this, there is an access to neighbor
 * cells (left, right, top, bottom, front, back), except if (i, j, k) is on the external surface. */
/*static inline double update_temperature(const double *T, int u, int n, int m, int o, int i, int j, int k)
{
		/* quantity of thermal energy that must be brought to a cell to make it heat up by 1°C
    const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl; /* J.K
    const double dl2 = dl * dl;
    double thermal_flux = 0;

    if (i > 0)
        thermal_flux += (T[u - 1] - T[u]) * sink_conductivity * dl; /* neighbor x-1
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (i < n - 1)
        thermal_flux += (T[u + 1] - T[u]) * sink_conductivity * dl; /* neighbor x+1
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (j > 0)
        thermal_flux += (T[u - n] - T[u]) * sink_conductivity * dl; /* neighbor y-1
    else {
        /* Bottom cell: does it receive it from the CPU ?
        if (CPU_shape(i * dl, k * dl))
            thermal_flux += CPU_TDP / CPU_surface * dl2;
        else {
            thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
            thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
        }
    }

    if (j < m - 1)
        thermal_flux += (T[u + n] - T[u]) * sink_conductivity * dl; /* neighbor y+1
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k > 0)
        thermal_flux += (T[u - n * m] - T[u]) * sink_conductivity * dl; /* neighbor z-1
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    if (k < o - 1)
        thermal_flux += (T[u + n * m] - T[u]) * sink_conductivity * dl; /* neighbor z+1
    else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
    }

    /* adjust temperature depending on the heat flux
    return T[u] + thermal_flux * dt / cell_heat_capacity;
}*/

/* Returns the new temperature of the cell (i, j, k). For this, there is an access to neighbor
 * cells (left, right, top, bottom, front, back), except if (i, j, k) is on the external surface. */
static inline double update_temperature_parallel(int cube_size_x, int cube_size_y, int cube_size_z, const double T[cube_size_x + 2][cube_size_y + 2][cube_size_z + 2],
    const int *blocks, int cube_index_x, int cube_index_y, int cube_index_z, int x, int y, int z, int n, int m, int o) {

    /* quantity of thermal energy that must be brought to a cell to make it heat up by 1°C */
    const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl; /* J.K */
    const double dl2 = dl * dl;
    double thermal_flux = 0;

    if (x > 1 || cube_index_x > 0) {
        thermal_flux += (T[x - 1][y][z] - T[x][y][z]) * sink_conductivity * dl; /* neighbor x-1 */
    } else {
        if (isnan(T[x][y][z])) fprintf(stderr, "err");
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[x][y][z], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[x][y][z] - (watercooling_T + 273.15));
    }

    if (x < cube_size_x || cube_index_x < blocks[0] - 1) {
        thermal_flux += (T[x + 1][y][z] - T[x][y][z]) * sink_conductivity * dl; /* neighbor x+1 */
    } else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[x][y][z], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[x][y][z] - (watercooling_T + 273.15));
    }

    if (y > 1 || cube_index_y > 0) {
        thermal_flux += (T[x][y - 1][z] - T[x][y][z]) *
                        sink_conductivity * dl; /* neighbor y-1 */
    } else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[x][y][z], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[x][y][z] - (watercooling_T + 273.15));
    }

    if (y < cube_size_y || cube_index_y < blocks[1] - 1) {
        thermal_flux += (T[x][y + 1][z] - T[x][y][z]) * sink_conductivity * dl; /* neighbor y+1 */
    } else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[x][y][z], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[x][y][z] - (watercooling_T + 273.15));
    }

    if (z > 1 || cube_index_z > 0) {
        thermal_flux += (T[x][y][z - 1] - T[x][y][z]) * sink_conductivity * dl; /* neighbor z-1 */
    } else {
        /* Bottom cell: does it receive it from the CPU ? */
        if (CPU_shape(((cube_index_x * n) / blocks[0] + x - 1) * dl, ((cube_index_y * m) / blocks[1] + y - 1) * dl)) {
            thermal_flux += CPU_TDP / CPU_surface * dl2;
        } else {
            thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[x][y][z], 4);
            thermal_flux -= heat_transfer_coefficient * dl2 * (T[x][y][z] - (watercooling_T + 273.15));
        }
    }

    if (z < cube_size_z || cube_index_z < blocks[2] - 1) {
        thermal_flux += (T[x][y][z + 1] - T[x][y][z]) * sink_conductivity * dl; /* neighbor z+1 */
    } else {
        thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[x][y][z], 4);
        thermal_flux -= heat_transfer_coefficient * dl2 * (T[x][y][z] - (watercooling_T + 273.15));
    }

    /* adjust temperature depending on the heat flux */
    return T[x][y][z] + thermal_flux * dt / cell_heat_capacity;
}

/* Run the simulation on the k-th xy plane.
 * v is the index of the start of the k-th xy plane in the arrays T and R. */
/*static inline void do_xy_plane(const double *T, double *R, int v, int n, int m, int o, int k)
{
    if (k == 0)
				// we do not modify the z = 0 plane: it is maintained at constant temperature via water-cooling
        return;

    for (int j = 0; j < m; j++) {   // y
        for (int i = 0; i < n; i++) {   // x
            int u = v + j * n + i;
            R[u] = update_temperature(T, u, n, m, o, i, j, k);
        }
    }
}*/

int main()
{
    MPI_Init(NULL, NULL);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    CPU_surface = CPU_contact_surface();
    double V = L * l * E;
    int n = ceil(L / dl);
    int m = ceil(l / dl);
    int o = ceil(E / dl);

    fprintf(stderr, "HEATSINK\n");
    fprintf(stderr, "\tDimension (cm) [x,y,z] = %.1f x %.1f x %.1f\n", 100 * L, 100 * E, 100 * l);
    fprintf(stderr, "\tVolume = %.1f cm^3\n", V * 1e6);
    fprintf(stderr, "\tWeight = %.2f kg\n", V * sink_density);
    fprintf(stderr, "\tPrice = %.2f €\n", V * sink_density * euros_per_kg);
    fprintf(stderr, "SIMULATION\n");
    fprintf(stderr, "\tGrid (x,y,z) = %d x %d x %d (%.1fMo)\n", n, m, o, 7.6293e-06 * n * m * o);
    fprintf(stderr, "\tdt = %gs\n", dt);
    fprintf(stderr, "CPU\n");
    fprintf(stderr, "\tPower = %.0fW\n", CPU_TDP);
    fprintf(stderr, "\tArea = %.1f cm^2\n", CPU_surface * 10000);

    int *blocks = malloc(sizeof(int) * 3);
    int *distances = malloc(sizeof(double) * 3);
    int *maxProc = malloc(sizeof(int) * 3);
    blocks[0] = blocks[1] = blocks[2] = 1;
    distances[0] = L;
    distances[1] = l;
    distances[2] = E;
    maxProc[0] = n;
    maxProc[1] = m;
    maxProc[2] = o;

    // determine blocks per dimension
    while (1) {
        int *order = malloc(sizeof (int) * 3);
        order[0] = 0;
        order[1] = 1;
        order[2] = 2;
        if (distances[order[0]] / blocks[order[0]] < distances[order[1]] / blocks[order[1]]) {
            int tmp = order[0];
            order[0] = order[1];
            order[1] = tmp;
        }
        if (distances[order[1]] / blocks[order[1]] < distances[order[2]] / blocks[order[2]]) {
            int tmp = order[1];
            order[1] = order[2];
            order[2] = tmp;
        }
        if (distances[order[0]] / blocks[order[0]] < distances[order[1]] / blocks[order[1]]) {
            int tmp = order[0];
            order[0] = order[1];
            order[1] = tmp;
        }

        if (blocks[order[0]] < maxProc[order[0]] && (blocks[order[0]] + 1) * blocks[order[1]] * blocks[order[2]] <= world_size) {
            blocks[order[0]]++;
        } else if (blocks[order[1]] < maxProc[order[1]] && blocks[order[0]] * (blocks[order[1]] + 1) * blocks[order[2]] <= world_size) {
            blocks[order[1]]++;
        } else if (blocks[order[2]] < maxProc[order[2]] && blocks[order[0]] * blocks[order[1]] * (blocks[order[2]] + 1) <= world_size) {
            blocks[order[2]]++;
        } else break;
    }

    free(maxProc);

    int used_processors = blocks[0] * blocks[1] * blocks[2];
    if (world_rank >= used_processors) {
        MPI_Finalize();
        return 0;
    }

    int cube_index_x = world_rank % blocks[0];
    int cube_index_y = (world_rank / blocks[0]) % blocks[1];
    int cube_index_z = world_rank / (blocks[0] * blocks[1]);

    int cube_size_x = ((cube_index_x + 1) * n) / blocks[0] - (cube_index_x * n) / blocks[0];
    int cube_size_y = ((cube_index_y + 1) * m) / blocks[1] - (cube_index_y * m) / blocks[1];
    int cube_size_z = ((cube_index_z + 1) * o) / blocks[2] - (cube_index_z * o) / blocks[2];

    double local_T[cube_size_x + 2][cube_size_y + 2][cube_size_z + 2];
    double local_R[cube_size_x][cube_size_y][cube_size_z];


    if (local_T == NULL || local_R == NULL) {
        perror("T or R could not be allocated");
        exit(1);
    }

    // initialise local_T
    for (int x = 0; x < cube_size_x + 2; x++) {
        for (int y = 1; y < cube_size_y + 2; y++) {
            for (int z = 1; z < cube_size_z + 2; z++) {
                if (x == 0 || y == 0 || z == 0 || x == cube_size_x + 1 || y == cube_size_y + 1 || z == cube_size_z + 1) {
                    local_T[x][y][z] = NAN;
                } else
                    local_T[x][y][z] = watercooling_T + 273.15;
            }
        }
    }



    double *x_left_send, *x_right_send, *y_left_send, *y_right_send, *z_left_send, *z_right_send,
        *x_left_recv, *x_right_recv, *y_left_recv, *y_right_recv, *z_left_recv, *z_right_recv;
    x_left_send = malloc(sizeof (double) * cube_size_y * cube_size_z);
    x_left_recv = malloc(sizeof (double) * cube_size_y * cube_size_z);
    x_right_send = malloc(sizeof (double) * cube_size_y * cube_size_z);
    x_right_recv = malloc(sizeof (double) * cube_size_y * cube_size_z);
    y_left_send = malloc(sizeof (double) * cube_size_x * cube_size_z);
    y_left_recv = malloc(sizeof (double) * cube_size_x * cube_size_z);
    y_right_send = malloc(sizeof (double) * cube_size_x * cube_size_z);
    y_right_recv = malloc(sizeof (double) * cube_size_x * cube_size_z);
    z_left_send = malloc(sizeof (double) * cube_size_x * cube_size_y);
    z_left_recv = malloc(sizeof (double) * cube_size_x * cube_size_y);
    z_right_send = malloc(sizeof (double) * cube_size_x * cube_size_y);
    z_right_recv = malloc(sizeof (double) * cube_size_x * cube_size_y);

    int has_x_left = cube_index_x != 0;
    int has_x_right = cube_index_x != blocks[0] - 1;
    int has_y_left = cube_index_y != 0;
    int has_y_right = cube_index_y != blocks[1] - 1;
    int has_z_left = cube_index_z != 0;
    int has_z_right = cube_index_z != blocks[2] - 1;
    int neighbours = has_x_left + has_x_right + has_y_left + has_y_right + has_z_left + has_z_right;
    MPI_Request *recv_reqs = malloc(sizeof (MPI_Request) * neighbours);

    double T_max, epsilon, local_T_max, local_epsilon, t = 0;

    int convergence = 0;
    while (!convergence) {
        MPI_Request x_left_send_req, x_right_send_req, y_left_send_req, y_right_send_req, z_left_send_req, z_right_send_req,
            x_left_recv_req, x_right_recv_req, y_left_recv_req, y_right_recv_req, z_left_recv_req, z_right_recv_req;

        if (has_x_left) {
            // send information to (cube_index_x - 1, cube_index_y, cube_index_z)
            for (int i = 0; i < cube_size_y * cube_size_z; i++) {
                int y = i % cube_size_y;
                int z = i / cube_size_y;
                x_left_send[i] = local_T[1][y+1][z+1];
            }
            MPI_Isend(x_left_send, cube_size_y * cube_size_z, MPI_DOUBLE, cube_index_x - 1 + cube_index_y * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &x_left_send_req);

            // receive
            MPI_Irecv(x_left_recv, cube_size_y * cube_size_z, MPI_DOUBLE, cube_index_x - 1 + cube_index_y * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &x_left_recv_req);
        }
        if (has_x_right) {
            // send information to (cube_index_x + 1, cube_index_y, cube_index_z)
            for (int i = 0; i < cube_size_y * cube_size_z; i++) {
                int y = i % cube_size_y;
                int z = i / cube_size_y;
                x_right_send[i] = local_T[cube_size_x][y+1][z+1];
            }
            MPI_Isend(x_right_send, cube_size_y * cube_size_z, MPI_DOUBLE, cube_index_x + 1 + cube_index_y * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &x_right_send_req);

            // receive
            MPI_Irecv(x_right_recv, cube_size_y * cube_size_z, MPI_DOUBLE, cube_index_x + 1 + cube_index_y * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &x_right_recv_req);
        }
        if (has_y_left) {
            // send information to (cube_index_x, cube_index_y - 1, cube_index_z)
            for (int i = 0; i < cube_size_x * cube_size_z; i++) {
                int x = i % cube_size_x;
                int z = i / cube_size_x;
                y_left_send[i] = local_T[x+1][1][z+1];
            }
            MPI_Isend(y_left_send, cube_size_x * cube_size_z, MPI_DOUBLE, cube_index_x + (cube_size_y - 1) * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &y_left_send_req);

            // receive
            MPI_Irecv(y_left_recv, cube_size_x * cube_size_z, MPI_DOUBLE, cube_index_x + (cube_size_y - 1) * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &y_left_recv_req);
        }
        if (has_y_right) {
            // send information to (cube_index_x, cube_index_y + 1, cube_index_z)
            for (int i = 0; i < cube_size_x * cube_size_z; i++) {
                int x = i % cube_size_x;
                int z = i / cube_size_x;
                y_right_send[i] = local_T[x+1][cube_size_y][z+1];
            }
            MPI_Isend(y_right_send, cube_size_x * cube_size_z, MPI_DOUBLE, cube_index_x + (cube_size_y + 1) * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &y_right_send_req);

            // receive
            MPI_Irecv(y_right_recv, cube_size_x * cube_size_z, MPI_DOUBLE, cube_index_x + (cube_size_y + 1) * blocks[0] + cube_index_z * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &y_right_recv_req);
        }
        if (has_z_left) {
            // send information to (cube_index_x, cube_index_y, cube_index_z - 1)
            for (int i = 0; i < cube_size_x * cube_size_y; i++) {
                int x = i % cube_size_x;
                int y = i / cube_size_x;
                z_left_send[i] = local_T[x+1][y+1][1];
            }
            MPI_Isend(z_left_send, cube_size_x * cube_size_y, MPI_DOUBLE, cube_index_x + cube_size_y * blocks[0] + (cube_index_z - 1) * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &z_left_send_req);

            // receive
            MPI_Irecv(z_left_recv, cube_size_x * cube_size_y, MPI_DOUBLE, cube_index_x + cube_size_y * blocks[0] + (cube_index_z - 1) * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &z_left_recv_req);
        }
        if (has_z_right) {
            // send information to (cube_index_x, cube_index_y, cube_index_z + 1)
            for (int i = 0; i < cube_size_x * cube_size_y; i++) {
                int x = i % cube_size_x;
                int y = i / cube_size_x;
                z_right_send[i] = local_T[x+1][y+1][cube_size_z];
            }
            MPI_Isend(z_right_send, cube_size_x * cube_size_y, MPI_DOUBLE, cube_index_x + cube_size_y * blocks[0] + (cube_index_z + 1) * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &z_right_send_req);

            // receive
            MPI_Irecv(z_right_recv, cube_size_x * cube_size_y, MPI_DOUBLE, cube_index_x + cube_size_y * blocks[0] + (cube_index_z + 1) * blocks[0] * blocks[1], 0, MPI_COMM_WORLD, &z_right_recv_req);
        }
        // calculate R interior
        for (int x = 1; x < cube_size_x - 1; x++) {
            for (int y = 1; y < cube_size_y - 1; y++) {
                for (int z = 1; z < cube_size_z - 1; z++) {
                    // local_R[x][y][z] = calc(T[x+1][y+1][z+1])
                    local_R[x][y][z] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, x+1, y+1, z+1, n, m, o);
                }
            }
        }

        // wait for receives
        for (int i = 0, j = 0; i < 6; ++i) {
            MPI_Request *req;
            switch (i) {
                case 0:
                    req = &x_left_recv_req;
                    if (!has_x_left)
                        continue;
                    break;
                case 1:
                    req = &x_right_recv_req;
                    if (!has_x_right)
                        continue;
                    break;
                case 2:
                    req = &y_left_recv_req;
                    if (!has_y_left)
                        continue;
                    break;
                case 3:
                    req = &y_right_recv_req;
                    if (!has_y_right)
                        continue;
                    break;
                case 4:
                    req = &z_left_recv_req;
                    if (!has_z_left)
                        continue;
                    break;
                case 5:
                     req = &z_right_recv_req;
                    if (!has_z_right)
                        continue;
                    break;
            }
            recv_reqs[j++] = *req;
        }
        MPI_Waitall(neighbours, recv_reqs, MPI_STATUSES_IGNORE);

        // store halo in T boundary
        if (has_x_left) {
            for (int i = 0; i < cube_size_y * cube_size_z; i++) {
                int y = i % cube_size_y;
                int z = i / cube_size_y;
                local_T[0][y+1][z+1] = x_left_recv[i];
            }
        }
        if (has_x_right) {
            for (int i = 0; i < cube_size_y * cube_size_z; i++) {
                int y = i % cube_size_y;
                int z = i / cube_size_y;
                local_T[cube_size_x+1][y+1][z+1] = x_right_recv[i];
            }
        }
        if (has_y_left) {
            for (int i = 0; i < cube_size_x * cube_size_z; i++) {
                int x = i % cube_size_x;
                int z = i / cube_size_x;
                local_T[x+1][0][z+1] = y_left_recv[i];
            }
        }
        if (has_y_right) {
            for (int i = 0; i < cube_size_x * cube_size_z; i++) {
                int x = i % cube_size_x;
                int z = i / cube_size_x;
                local_T[x+1][cube_size_y+1][z+1] = y_right_recv[i];
            }
        }
        if (has_z_left) {
            for (int i = 0; i < cube_size_x * cube_size_y; i++) {
                int x = i % cube_size_x;
                int y = i / cube_size_x;
                local_T[x+1][y+1][0] = z_left_recv[i];
            }
        }
        if (has_z_right) {
            for (int i = 0; i < cube_size_x * cube_size_y; i++) {
                int x = i % cube_size_x;
                int y = i / cube_size_x;
                local_T[x+1][y+1][cube_size_z+1] = z_right_recv[i];
            }
        }

        // calculate R boundary
        for (int y = 0; y < cube_size_y; y++) {
            for (int z = 0; z < cube_size_z; z++) {
                local_R[0][y][z] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, 1, y+1, z+1, n, m, o);
                local_R[cube_size_x-1][y][z] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, cube_size_x, y+1, z+1, n, m, o);
            }
        }
        for (int x = 1; x < cube_size_x - 1; x++) {
            for (int z = 0; z < cube_size_z; z++) {
                local_R[x][0][z] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, x+1, 1, z+1, n, m, o);
                local_R[x][cube_size_y-1][z] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, x+1, cube_size_y, z+1, n, m, o);
            }
        }
        for (int x = 1; x < cube_size_x - 1; x++) {
            for (int y = 1; y < cube_size_y - 1; y++) {
                local_R[x][y][0] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, x+1, y+1, 1, n, m, o);
                local_R[x][y][cube_size_z-1] = update_temperature_parallel(cube_size_x, cube_size_y, cube_size_z, local_T, blocks, cube_index_x, cube_index_y, cube_index_z, x+1, y+1, cube_size_z, n, m, o);
            }
        }

        // local T_max and epsilon
        local_T_max = local_R[0][0][0];
        double tmp;
        local_epsilon = 0;
        for (int x = 0; x < cube_size_x; ++x) {
            for (int y = 0; y < cube_size_y; ++y) {
                for (int z = 0; z < cube_size_z; ++z) {
                    if (local_R[x][y][z] > local_T_max) {
                        local_T_max = local_R[x][y][z];
                    }
                    tmp = local_R[x][y][z] - local_T[x+1][y+1][z+1];
                    local_epsilon += tmp * tmp;
                }
            }
        }

        // copy R to T
        for (int x = 0; x < cube_size_x; ++x) {
            for (int y = 0; y < cube_size_y; ++y) {
                for (int z = 0; z < cube_size_z; ++z) {
                    local_T[x+1][y+1][z+1] = local_R[x][y][z];
                }
            }
        }

        // renvoyer T_max et epsilon
        MPI_Allreduce(&local_T_max, &T_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // reduce T_max
        MPI_Allreduce(&local_epsilon, &epsilon, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // reduce epsilon

        // convergence = sqrt(epsilon)/dt < 0.1
        convergence = sqrt(epsilon)/dt < 0.1;
        if (world_rank == 0) {
            fprintf(stderr, "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n", t, T_max - 273.15, epsilon);
        }
        t += dt;
    }
    fprintf(stderr,"rank %d we good\n",world_rank);
    free(x_left_send);
    free(x_left_recv);
    free(x_right_send);
    free(x_right_recv);
    free(y_left_send);
    free(y_left_recv);
    free(y_right_send);
    free(y_right_recv);
    free(z_left_send);
    free(z_left_recv);
    free(z_right_send);
    free(z_right_recv);
    fprintf(stderr,"rank %d we still good\n",world_rank);

    /* temperature of each cell, in degree Kelvin. */
    /*double *T = malloc(n * m * o * sizeof(*T));
    double *R = malloc(n * m * o * sizeof(*R));
    if (T == NULL || R == NULL) {
        perror("T or R could not be allocated");
        exit(1);
    }*/

    /* initially the heatsink is at the temperature of the water-cooling fluid */
    /*for (int u = 0; u < n * m * o; u++)
        R[u] = T[u] = watercooling_T + 273.15;*/

    /* let's go! we switch the CPU on and launch the simulation until it reaches a stationary state. */
    // double t = 0;
    // int n_steps = 0;
    // int convergence = 0;
    //
    // /* simulating time steps */
    // while (convergence == 0) {
    //     /* Update all cells. xy planes are processed, for increasing values of z. */
    //     for (int k = 0; k < o; k++) {   // z
    //         int v = k * n * m;
    //         do_xy_plane(T, R, v, n, m, o, k);
    //     }
    //
    //     /* each second, we test the convergence, and print a short progress report */
    //     if (n_steps % ((int)(1 / dt)) == 0) {
    //         double delta_T = 0;
    //         double max = -INFINITY;
    //         for (int u = 0; u < n * m * o; u++) {
    //             delta_T += (R[u] - T[u]) * (R[u] - T[u]);
    //             if (R[u] > max)
    //                 max = R[u];
    //         }
    //         delta_T = sqrt(delta_T) / dt;
    //         fprintf(stderr, "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n", t, max - 273.15, delta_T);
    //         if (delta_T < 0.1)
    //             convergence = 1;
    //     }
    //
    //     /* the new temperatures are in R */
    //     double *tmp = R;
    //     R = T;
    //     T = tmp;
    //
    //     t += dt;
    //     n_steps += 1;

#ifdef DUMP_STEADY_STATE
    if (world_rank) {
        double *send = malloc(sizeof(double) * 3);
        send[0] = cube_size_x;
        send[1] = cube_size_y;
        send[2] = cube_size_z;
        MPI_Send(send, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_R, cube_size_x * cube_size_y * cube_size_z, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
        double global_T[n][m][o];
        int cube_size[3];
        double *recv;
        for (int z_i = 0, z_offset = 0; z_i < blocks[2]; z_i++) {
            for (int y_i = 0, y_offset = 0; y_i < blocks[1]; y_i++) {
                for (int x_i = 0, x_offset = 0; x_i < blocks[0]; x_i++) {
                    int rank = x_i + y_i*blocks[0] + z_i*blocks[0]*blocks[1];
                    if (rank == world_rank) {
                        cube_size[0] = cube_size_x;
                        cube_size[1] = cube_size_y;
                        cube_size[2] = cube_size_z;
                        for (int x = 0; x < cube_size[0]; x++) {
                            for (int y = 0; y < cube_size[1]; y++) {
                                for (int z = 0; z < cube_size[2]; z++) {
                                    int x_pos = x + x_offset;
                                    int y_pos = y + y_offset;
                                    int z_pos = z + z_offset;
                                    global_T[x_pos][y_pos][z_pos] = local_R[x][y][z];
                                }
                            }
                        }
                    } else {
                        MPI_Recv(cube_size, 3, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        recv = malloc(sizeof(double) * cube_size[0]*cube_size[1]*cube_size[2]);
                        MPI_Recv(recv, cube_size[0]*cube_size[1]*cube_size[2], MPI_DOUBLE, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        for (int x = 0; x < cube_size[0]; x++) {
                            for (int y = 0; y < cube_size[1]; y++) {
                                for (int z = 0; z < cube_size[2]; z++) {
                                    int x_pos = x + x_offset;
                                    int y_pos = y + y_offset;
                                    int z_pos = z + z_offset;
                                    global_T[x_pos][y_pos][z_pos] = recv[x + y*cube_size[0] + z*cube_size[0]*cube_size[1]];
                                }
                            }
                        }
                    }
                    x_offset += cube_size[0];
                    if (x_i == blocks[0] - 1) {
                        y_offset += cube_size[1];
                        if (y_i == blocks[1] - 1)
                            z_offset += cube_size[2];
                    }
                }
            }
        }
        printf("###### STEADY STATE; t = %.1f\n", t);
        for (int k = 0; k < o; k++) {   // z
            printf("# z = %g\n", k * dl);
            for (int j = 0; j < m; j++) {   // y
                for (int i = 0; i < n; i++) {   // x
                    printf("%.1f ", global_T[i][j][k] - 273.15);
                }
                printf("\n");
            }
        }
        printf("\n");
        fprintf(stderr, "For graphical rendering: python3 rendu_picture_steady.py [filename.txt] %d %d %d\n", n, m, o);
    }
#endif
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}