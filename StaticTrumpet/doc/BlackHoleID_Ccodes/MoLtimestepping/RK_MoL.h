// C code implementation of Euler Method of Lines timestepping.
// ***Euler timestepping only requires one RHS evaluation***

LOOP_ALL_GFS_GPS(i) {
  y_n_gfs[i] = y_n_gfs[i] + y_nplus1_running_total_gfs[i]*dt;
}


