#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define N M
#define M_LEN M + 1
#define N_LEN N + 1
#define ITMAX 4000
int main(int argc, char **argv) {
  // solution arrays
  double u[M_LEN][N_LEN],v[M_LEN][N_LEN],p[M_LEN][N_LEN];
  double unew[M_LEN][N_LEN],vnew[M_LEN][N_LEN],pnew[M_LEN][N_LEN];
  double uold[M_LEN][N_LEN],vold[M_LEN][N_LEN],pold[M_LEN][N_LEN];
  double cu[M_LEN][N_LEN],cv[M_LEN][N_LEN],z[M_LEN][N_LEN],h[M_LEN][N_LEN],psi[M_LEN][N_LEN];

  double dt,tdt,dx,dy,a,alpha,el,pi,tpi,di,dj,pcf;
  double tdts8,tdtsdx,tdtsdy,fsdx,fsdy;
  int mnmin,ncycle,i,j;
 
  // Note below that two delta t (tdt) is set to dt on the first
  // cycle after which it is reset to dt+dt.
  dt = 90.;
  tdt = dt;
  dx = 100000.;
  dy = 100000.;
  fsdx = 4. / dx;
  fsdy = 4. / dy;

  a = 1000000.;
  alpha = .001;
  el = N * dx;
  pi = 4. * atanf(1.);
  tpi = pi + pi;
  di = tpi / M;
  dj = tpi / N;
  pcf = pi * pi * a * a / (el * el);

  // Initial values of the stream function and p
  for (i=0;i<M_LEN;++i) {
    for (j=0;j<N_LEN;++j) {
      psi[i][j] = a * sin((i + .5) * di) * sin((j + .5) * dj);
      p[i][j] = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.;
    }
  }
    
  // Initialize velocities
  for (i=0;i<M;++i) {
    for (j=0;j<N;++j) {
      u[i + 1][j] = -(psi[i + 1][j + 1] - psi[i + 1][j]) / dy;
      v[i][j + 1] = (psi[i + 1][j + 1] - psi[i][j + 1]) / dx;
    }
  }
     
  // Periodic continuation
  for (j=0;j<N;++j) {
    u[0][j] = u[M][j];
    v[M][j + 1] = v[0][j + 1];
  }
  for (i=0;i<M;++i) {
    u[i + 1][N] = u[i + 1][0];
    v[i][0] = v[i][N];
  }
  u[0][N] = u[M][0];
  v[M][0] = v[0][N];
  for (i=0;i<M_LEN;++i) {
    for (j=0;j<N_LEN;++j) {
      uold[i][j] = u[i][j];
      vold[i][j] = v[i][j];
      pold[i][j] = p[i][j];
    }
  }
     
  // ** Start of time loop ** 
  for (ncycle=1;ncycle<=ITMAX;ncycle++) {
    // Compute capital u, capital v, z and h
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        cu[i + 1][j] = .5 * (p[i + 1][j] + p[i][j]) * u[i + 1][j];
      }
    }
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        cv[i][j + 1] = .5 * (p[i][j + 1] + p[i][j]) * v[i][j + 1];
      }
    }
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        z[i + 1][j + 1] = (fsdx * (v[i + 1][j + 1] - v[i][j + 1]) - fsdy * (u[i + 1][j + 1] - u[i + 1][j])) / (p[i][j] + p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);
      }
    }
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        h[i][j] = p[i][j] + .25 * (u[i + 1][j] * u[i + 1][j] + u[i][j] * u[i][j] + v[i][j + 1] * v[i][j + 1] + v[i][j] * v[i][j]);
      }
    }

    // Periodic continuation
    for (j=0;j<N;++j) {
      cu[0][j] = cu[M][j];
      cv[M][j + 1] = cv[0][j + 1];
      z[0][j + 1] = z[M][j + 1];
      h[M][j] = h[0][j];
    }
    for (i=0;i<M;++i) {
      cu[i + 1][N] = cu[i + 1][0];
      cv[i][0] = cv[i][N];
      z[i + 1][0] = z[i + 1][N];
      h[i][N] = h[i][0];
    }
    cu[0][N] = cu[M][0];
    cv[M][0] = cv[0][N];
    z[0][0] = z[M][N];
    h[M][N] = h[0][0];
     
    // Compute new values u,v and p
    tdts8 = tdt / 8.;
    tdtsdx = tdt / dx;
    tdtsdy = tdt / dy;
    
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        unew[i + 1][j] = uold[i + 1][j] + tdts8 * (z[i + 1][j + 1] + z[i + 1][j]) * (cv[i + 1][j + 1] + cv[i][j + 1] + cv[i][j] + cv[i + 1][j]) - tdtsdx * (h[i + 1][j] - h[i][j]);
      }
    }
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        vnew[i][j + 1] = vold[i][j + 1] - tdts8 * (z[i + 1][j + 1] + z[i][j + 1]) * (cu[i + 1][j + 1] + cu[i][j + 1] + cu[i][j] + cu[i + 1][j]) - tdtsdy * (h[i][j + 1] - h[i][j]);
      }
    }
    for (i=0;i<M;++i) {
      for (j=0;j<N;++j) {
        pnew[i][j] = pold[i][j] - tdtsdx * (cu[i + 1][j] - cu[i][j]) - tdtsdy * (cv[i][j + 1] - cv[i][j]); 
      }
    }

    // Periodic continuation
    for (j=0;j<N;++j) {
      unew[0][j] = unew[M][j];
      vnew[M][j + 1] = vnew[0][j + 1];
      pnew[M][j] = pnew[0][j];
    }
    for (i=0;i<M;++i) {
      unew[i + 1][N] = unew[i + 1][0];
      vnew[i][0] = vnew[i][N];
      pnew[i][N] = pnew[i][0];
    }
    unew[0][N] = unew[M][0];
    vnew[M][0] = vnew[0][N];
    pnew[M][N] = pnew[0][0];

    // Time smoothing and update for next cycle
    if ( ncycle > 1 ) {
      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          uold[i][j] = u[i][j] + alpha * (unew[i][j] - 2. * u[i][j] + uold[i][j]);
        }
      }
      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          vold[i][j] = v[i][j] + alpha * (vnew[i][j] - 2. * v[i][j] + vold[i][j]);
        }
      }
      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          pold[i][j] = p[i][j] + alpha * (pnew[i][j] - 2. * p[i][j] + pold[i][j]);
        }
      }
      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          u[i][j] = unew[i][j];
        }
      }
      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          v[i][j] = vnew[i][j];
        }
      }
      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          p[i][j] = pnew[i][j];
        }
      }
    } else {
      tdt = tdt + tdt;

      for (i=0;i<M_LEN;++i) {
        for (j=0;j<N_LEN;++j) {
          uold[i][j] = u[i][j];
          vold[i][j] = v[i][j];
          pold[i][j] = p[i][j];
          u[i][j] = unew[i][j];
          v[i][j] = vnew[i][j];
          p[i][j] = pnew[i][j];
        }
      }
    }
  }
}