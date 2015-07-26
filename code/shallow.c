#define N M
#define M_LEN (M+1)
#define N_LEN (N+1)
#define ITMAX   4000
#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

int main(int argc, char **argv) {
    // solution arrays
    double u[(M_LEN+1)*(N_LEN+1)];
    double tdt,a,alpha,el,pi,tpi,di,dj,pcf;
    double tdts8,tdtsdx,tdtsdy,fsdx,fsdy;
    int ncycle,i,j;

    // opencl variables
    cl_kernel kernel_init1, kernel_init2, kernel_init_pc, kernel_l100, kernel_l100_pc, kernel_l200, kernel_l200_pc, kernel_l300, kernel_l300_pc;

    cl_mem buf_u, buf_v, buf_p;
    cl_mem buf_uold, buf_vold, buf_pold, buf_unew, buf_vnew, buf_pnew;
    cl_mem buf_cu, buf_cv, buf_z, buf_h, buf_psi;

    size_t global_worksize[2] = {M_LEN, N_LEN};
#ifdef M_BLOCK_LEN
    size_t local_worksize[2] = {M_BLOCK_LEN, N_BLOCK_LEN};
#else
    size_t* local_worksize = NULL;
#endif

    int elements = (M_LEN+1) * (N_LEN+1);
    int datasize = elements * sizeof(double);

    // ** Initialise vars ** 
    // Note below that two delta t (tdt) is set to dt on the first
    // cycle after which it is reset to dt+dt.
    const double dt = 90.;
    tdt = dt;
    const double dx = 100000.;
    const double dy = 100000.;
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
    // Buffer and kernel initialisation omitted...

    // *** Initialise grids ***
    clEnqueueNDRangeKernel(cmd_queue, kernel_init1, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);
    clFinish(cmd_queue);
    clEnqueueNDRangeKernel(cmd_queue, kernel_init2, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);
    clEnqueueNDRangeKernel(cmd_queue, kernel_init_pc, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);

    // ** Start of time loop ** 
    for (ncycle=1; ncycle <= ITMAX; ncycle++) {
        // l100: Compute capital u, capital v, z and h
        clEnqueueNDRangeKernel(cmd_queue, kernel_l100, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);

        // Periodic continuation
        clEnqueueNDRangeKernel(cmd_queue, kernel_l100_pc, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);

        // Compute new values u,v and p
        tdts8 = tdt / 8.0;
        tdtsdx = tdt / dx;
        tdtsdy = tdt / dy;

        // l200
        clSetKernelArg(kernel_l200,  0, sizeof(double), (void *)&tdts8);
        clSetKernelArg(kernel_l200,  1, sizeof(double), (void *)&tdtsdx);
        clSetKernelArg(kernel_l200,  2, sizeof(double), (void *)&tdtsdy);
        clEnqueueNDRangeKernel(cmd_queue, kernel_l200, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);

        // Periodic continuation
        clEnqueueNDRangeKernel(cmd_queue, kernel_l200_pc, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);
 
        // Time smoothing and update for next cycle
        if ( ncycle > 1 ) {
            // l300
            clEnqueueNDRangeKernel(cmd_queue, kernel_l300, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);
        } else {
            tdt = tdt + tdt;
            clEnqueueNDRangeKernel(cmd_queue, kernel_l300_pc, 2, NULL, global_worksize, local_worksize, 0, NULL, NULL);
        }
    }
    // ** End of time loop ** 
    // Read output buffer back to host. Again, blocking reads will force sync to host
    clEnqueueReadBuffer(cmd_queue, buf_p, CL_TRUE, 0, datasize, p, 0, NULL, NULL);
    clEnqueueReadBuffer(cmd_queue, buf_u, CL_TRUE, 0, datasize, u, 0, NULL, NULL);
    clEnqueueReadBuffer(cmd_queue, buf_v, CL_TRUE, 0, datasize, v, 0, NULL, NULL);
}
