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

void InitParticles(Particle[], ParticleV[], int);
double ComputeForces(Particle[], Particle[], ParticleV[], int);
double ComputeNewPos(Particle[], ParticleV[], int, double);

int main()
{
  omp_set_dynamic(0);

  double time;
  Particle *particles; /* Particles */
  ParticleV *particles_velocities;       /* Particle velocity */
  int particle_count, i, j;
  int loop_count;      /* number of times in loop */
  double sim_t; /* Simulation time */
  int tmp;
  tmp = fscanf(stdin, "%d\n", &particle_count);
  tmp = fscanf(stdin, "%d\n", &loop_count);

  /* Allocate memory for particles */
  particles = (Particle *)malloc(sizeof(Particle) * particle_count);
  particles_velocities = (ParticleV *)malloc(sizeof(ParticleV) * particle_count);

  /* Generate the initial values */
  InitParticles(particles, particles_velocities, particle_count);
  sim_t = 0.0;

  while (loop_count--)
  {
    double max_f;
    /* Compute forces (2D only) */
    max_f = ComputeForces(particles, particles, particles_velocities, particle_count);
    /* Once we have the forces, we compute the changes in position */
    sim_t += ComputeNewPos(particles, particles_velocities, particle_count, max_f);
  }
  for (i = 0; i < particle_count; i++)
    fprintf(stdout, "%.5lf %.5lf %.5lf\n", particles[i].x, particles[i].y, particles[i].z);
  return 0;
}

void InitParticles(Particle particles[], ParticleV particles_velocities[], int particle_count)
{
  int i;
  // #pragma omp parallel for private(i)
  for (i = 0; i < particle_count; i++)
  {
    particles[i].x = Random();
    particles[i].y = Random();
    particles[i].z = Random();
    particles[i].mass = 1.0;
    particles_velocities[i].xold = particles[i].x;
    particles_velocities[i].yold = particles[i].y;
    particles_velocities[i].zold = particles[i].z;
    particles_velocities[i].fx = 0;
    particles_velocities[i].fy = 0;
    particles_velocities[i].fz = 0;
  }
}

// Note: fast but has precision issues. Use only as appendix or get proper terminology
double sqrt(double number)
{
	long i;
	double x2, y;
	const double threehalfs = 1.5F;

	x2 = number * 0.5F;
	y  = number;
	i  = * ( long * ) &y;                       // evil floating point bit level hacking
	i  = 0x5fe6eb50c7b537a9 - ( i >> 1 );               // what the fuck? 
	y  = * ( double * ) &i;
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
	y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed

	return 1 / y;
}

double ComputeForces(Particle myparticles[], Particle others[], ParticleV particles_velocities[], int particle_count)
{

  double xis[particle_count];
  double yis[particle_count];

  double rmins[particle_count];
  double fxs [particle_count];
  double fys[particle_count];
  
  omp_lock_t rmins_locks[particle_count];
  omp_lock_t fxs_locks[particle_count];
  omp_lock_t fys_locks[particle_count];

  int i;
  int j;
  int k;

  #pragma omp parallel for
  for (i = 0; i < particle_count; i++)
  {
    xis[i] = myparticles[i].x;
    yis[i] = myparticles[i].y;
    
    // locked ones
    rmins[i] = 100.0;
    fxs[i] = 0.0;
    fys[i] = 0.0;

    omp_init_lock(&rmins_locks[i]);
  }

  double mi, rx, ry, mj, r, fx_delta, fy_delta;
  #pragma omp parallel for schedule(guided) private(i, j, mi, rx, ry, mj, r, fx_delta, fy_delta) 
  for (k = 0; k < particle_count * particle_count; k++)
  {
    int i = k / particle_count;
    int j = k % particle_count;

    /* ignore overlap and same particle */
    if (i == j)
      continue;

    rx = xis[i] - others[j].x;
    ry = yis[i] - others[j].y;
    mj = others[j].mass;
    r = rx * rx + ry * ry;
    
    /* second check for overlap. might not be necessary */
    if (r == 0.0) 
      continue;
    
    // Need to guarantee that: 
    // - rmin is the minimum r
    // - no race condition on the sum of fxs[i] -= mj * rx / r;
    // - no race condition on the sum of fys[i] -= mj * ry / r;

    omp_set_lock(&rmins_locks[i]);
    if (r < rmins[i])
       rmins[i] = r;
    omp_unset_lock(&rmins_locks[i]);

    r = r * sqrt(r);
    fx_delta = mj * rx / r;
    fy_delta = mj * ry / r;

    #pragma omp atomic
    fxs[i] -= fx_delta;

    #pragma omp atomic
    fys[i] -= fy_delta;
  }

  double new_max_f = 0.0;
  #pragma omp parallel for reduction(max:new_max_f) 
  for (i = 0; i < particle_count; i++)
  {
    particles_velocities[i].fx += fxs[i];
    particles_velocities[i].fy += fys[i];
    new_max_f = sqrt(fxs[i] * fxs[i] +  fys[i] *  fys[i]) / rmins[i];

    // Also destroy locks since we're only reading
    omp_destroy_lock(&rmins_locks[i]);
  }
  
  return new_max_f;
}

double ComputeNewPos(Particle particles[], ParticleV particles_velocities[], int particle_count, double max_f)
{
  int i;
  double a0, a1, a2;
  static double dt_old = 0.001, dt = 0.001;
  double dt_new;
  a0 = 2.0 / (dt * (dt + dt_old));
  a2 = 2.0 / (dt_old * (dt + dt_old));
  a1 = -(a0 + a2);

  double xi, yi;
  #pragma omp parallel for private(xi, yi)
  for (i = 0; i < particle_count; i++)
  {
    xi = particles[i].x;
    yi = particles[i].y;
    particles[i].x = (particles_velocities[i].fx - a1 * xi - a2 * particles_velocities[i].xold) / a0;
    particles[i].y = (particles_velocities[i].fy - a1 * yi - a2 * particles_velocities[i].yold) / a0;
    particles_velocities[i].xold = xi;
    particles_velocities[i].yold = yi;
    particles_velocities[i].fx = 0;
    particles_velocities[i].fy = 0;
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
