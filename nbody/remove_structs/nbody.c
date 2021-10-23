#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/*
 * pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
 */
#define MODULUS 2147483647
#define MULTIPLIER 48271
#define DEFAULT 123456789

static long seed = DEFAULT;

double Random(void)
/* ----------------------------------------------------------------
 * Random returns a pseudo-random real number uniformly distributed 
 * between 0.0 and 1.0. 
 * ----------------------------------------------------------------
 */
{
  const long Q = MODULUS / MULTIPLIER;
  const long R = MODULUS % MULTIPLIER;
  long t;

  t = MULTIPLIER * (seed % Q) - R * (seed / Q);
  if (t > 0)
    seed = t;
  else
    seed = t + MODULUS;
  return ((double)seed / MODULUS);
}

/*
 * End of the pRNG algorithm
 */

typedef struct
{
  double x, y, z;
  double mass;
} Particle;
typedef struct
{
  double xold, yold, zold;
  double fx, fy, fz;
} ParticleV;

void InitParticles(double particles_x[], 
                      double particles_y[], 
                      double particles_z[], 
                      double particles_mass[], 
                      double particles_v_fx[], 
                      double particles_v_fy[], 
                      double particles_v_fz[], 
                      double particles_v_old_x[], 
                      double particles_v_old_y[],
                      double particles_v_old_z[],
                      int particle_count);

double ComputeForces(double myparticles_x[], 
                      double myparticles_y[], 
                      double others_x[], 
                      double others_y[], 
                      double others_mass[], 
                      double particles_velocities_fx[], 
                      double particles_velocities_fy[], 
                      int particle_count);

double ComputeNewPos(double particles_x[], 
                      double particles_y[], 
                      double particles_velocities_fx[], 
                      double particles_velocities_fy[], 
                      double particles_velocities_old_x[], 
                      double particles_velocities_old_y[], 
                      int particle_count, 
                      double max_f);

int main()
{
  omp_set_dynamic(0);

  double time;
  int particle_count, i, j;
  int loop_count;      /* number of times in loop */
  double sim_t; /* Simulation time */
  int tmp;
  tmp = fscanf(stdin, "%d\n", &particle_count);
  tmp = fscanf(stdin, "%d\n", &loop_count);

  /* Allocate memory for particles */
  double *particles_x = (double *)malloc(sizeof(double) * particle_count);
  double *particles_y = (double *)malloc(sizeof(double) * particle_count);
  double *particles_z = (double *)malloc(sizeof(double) * particle_count);
  double *particles_mass = (double *)malloc(sizeof(double) * particle_count);

  double *particles_v_fx = (double *)malloc(sizeof(double) * particle_count);
  double *particles_v_fy = (double *)malloc(sizeof(double) * particle_count);
  double *particles_v_fz = (double *)malloc(sizeof(double) * particle_count);
  double *particles_v_old_x = (double *)malloc(sizeof(double) * particle_count);
  double *particles_v_old_y = (double *)malloc(sizeof(double) * particle_count);
  double *particles_v_old_z = (double *)malloc(sizeof(double) * particle_count);


  /* Generate the initial values */
  InitParticles(particles_x,
                particles_y,
                particles_z,
                particles_mass,
                particles_v_fx,
                particles_v_fy,
                particles_v_fz,
                particles_v_old_x,
                particles_v_old_y,
                particles_v_old_z, 
                particle_count);
  sim_t = 0.0;

  while (loop_count--)
  {
    double max_f;
    /* Compute forces (2D only) */
    max_f = ComputeForces(particles_x,
                          particles_y,
                          particles_x,
                          particles_y,
                          particles_mass,
                          particles_v_fx,
                          particles_v_fy, 
                          particle_count);
    /* Once we have the forces, we compute the changes in position */
    sim_t += ComputeNewPos(particles_x, 
                            particles_y, 
                            particles_v_fx, 
                            particles_v_fy, 
                            particles_v_old_x, 
                            particles_v_old_y,
                            particle_count, max_f);
  }
  for (i = 0; i < particle_count; i++)
    fprintf(stdout, "%.5lf %.5lf %.5lf\n", particles_x[i], particles_y[i], particles_z[i]);
  return 0;
}

void InitParticles(double particles_x[], 
                      double particles_y[], 
                      double particles_z[], 
                      double particles_mass[], 
                      double particles_v_fx[], 
                      double particles_v_fy[], 
                      double particles_v_fz[], 
                      double particles_v_old_x[], 
                      double particles_v_old_y[],
                      double particles_v_old_z[],
                      int particle_count)
{
  int i;
  // #pragma omp parallel for private(i)
  for (i = 0; i < particle_count; i++)
  {
    particles_x[i] = Random();
    particles_y[i] = Random();
    particles_z[i] = Random();
    particles_mass[i] = 1.0;
    particles_v_old_x[i] = particles_x[i];
    particles_v_old_y[i] = particles_y[i];
    particles_v_old_z[i] = particles_z[i];
    particles_v_fx[i] = 0;
    particles_v_fy[i] = 0;
    particles_v_fz[i] = 0;
  }
}

double ComputeForces(double myparticles_x[], 
                      double myparticles_y[], 
                      double others_x[], 
                      double others_y[], 
                      double others_mass[], 
                      double particles_v_fx[], 
                      double particles_v_fy[], 
                      int particle_count)
{
  double max_f, new_max_f;
  int i;
  max_f = 0.0;
  new_max_f = 0.0;
  int j;
  double xi, yi, mi, rx, ry, mj, r, fx, fy, rmin;
  #pragma omp parallel for simd schedule(guided) reduction(max:new_max_f) private(i, j, xi, yi, mi, rx, ry, mj, r, fx, fy, rmin) 
  for (i = 0; i < particle_count; i++)
  {
    rmin = 100.0;
    xi = myparticles_x[i];
    yi = myparticles_y[i];
    fx = 0.0;
    fy = 0.0;
    for (j = 0; j < particle_count; j++)
    {
      /* ignore overlap and same particle */
      if (i == j)
        continue;

      rx = xi - others_x[j];
      ry = yi - others_y[j];
      mj = others_mass[j];
      r = rx * rx + ry * ry;
      
      /* second check for overlap. might not be necessary */
      if (r == 0.0) 
        continue;
      
      if (r < rmin)
        rmin = r;

      r = r * sqrt(r);
      fx -= mj * rx / r;
      fy -= mj * ry / r;
    }
    particles_v_fx[i] += fx;
    particles_v_fy[i] += fy;
    new_max_f = sqrt(fx * fx + fy * fy) / rmin;
  }
  return new_max_f;
}

double ComputeNewPos(double particles_x[], 
                      double particles_y[], 
                      double particles_v_fx[], 
                      double particles_v_fy[], 
                      double particles_v_old_x[], 
                      double particles_v_old_y[], 
                      int particle_count, 
                      double max_f)
{
  int i;
  double a0, a1, a2;
  static double dt_old = 0.001, dt = 0.001;
  double dt_new;
  a0 = 2.0 / (dt * (dt + dt_old));
  a2 = 2.0 / (dt_old * (dt + dt_old));
  a1 = -(a0 + a2);

  double xi, yi;
  #pragma omp parallel for simd private(xi, yi) 
  for (i = 0; i < particle_count; i++)
  {
    xi = particles_x[i];
    yi = particles_y[i];
    particles_x[i] = (particles_v_fx[i] - a1 * xi - a2 * particles_v_old_x[i]) / a0;
    particles_y[i] = (particles_v_fy[i] - a1 * yi - a2 * particles_v_old_y[i]) / a0;
    particles_v_old_x[i] = xi;
    particles_v_old_y[i] = yi;
    particles_v_fx[i] = 0;
    particles_v_fy[i] = 0;
  }

  dt_new = 1.0 / sqrt(max_f);

  /* Set a minimum: */
  if (dt_new < 1.0e-6)
    dt_new = 1.0e-6;
  /* Modify time step */

  if (dt_new < dt)
  {
    dt_old = dt;
    dt = dt_new;
  }
  else if (dt_new > 4.0 * dt)
  {
    dt_old = dt;
    dt *= 2.0;
  }

  return dt_old;
}
