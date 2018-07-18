#ifndef blackhole_lab_h
#define blackhole_lab_h

#ifdef KERR_SPACETIME
    #include "KBLspacetime.h"
#else
    #include "SBLspacetime.h"
#endif

/* Parameters for lab setting         */
#define C_PI            acos(-1.)
#define C_L1            10.                     /* distance between observer and black hole */
#define C_L2            10.                     /* distance between black hole and picture  */
#define C_l             30.                     /* length of side of picture                */
#define C_VISUAL_ANGLE  120.                    /* visual angle of vision                   */
#define C_d             (1. / (sqrt(2.) * tan(C_VISUAL_ANGLE / 2. / 180. * C_PI)))
                                                /* distance between observer and vision     */
/* Parameters for resolution of vision */
#define C_OUTPUT_FILENAME         "data"
#ifndef C_RESOLUTION_WIDTH
    #define C_RESOLUTION_WIDTH        100
#endif
#ifndef C_RESOLUTION_HEIGHT
    #define C_RESOLUTION_HEIGHT       100
#endif
#define C_RESOLUTION_TOTAL_PIXELS (C_RESOLUTION_WIDTH * C_RESOLUTION_HEIGHT)

/* Parameters for time integration     */
#define C_DELTA_TIME                -0.05
#define C_TOTAL_STEP                20000

/* For classification of rays          */
enum _RAY_STATUS {HIT, FALL, OUTSIDE, YET};

/* For initialization of quantities    */
#define INITIAL_X_r(_w,_h,v)          C_L1
#define INITIAL_X_theta(_w,_h,v)      (C_PI / 2.)
#define INITIAL_X_phi(_w,_h,v)        0.
#define INITIAL_U_r(_w,_h,v)          (C_d * U_unit(r,v))
#define INITIAL_U_theta(_w,_h,v)      (-(-0.5 + (double)(_h) / ((double)(C_RESOLUTION_HEIGHT - 1))) * U_unit(h,v))
#define INITIAL_U_phi(_w,_h,v)        (-(-0.5 + (double)(_w) / ((double)(C_RESOLUTION_WIDTH  - 1))) * U_unit(p,v))
#define INITIAL_U_t(_w,_h,v)          U_t(v)

/* Criteria for rays */
#define BOOL_OUTSIDE_BOUNDARY(v)    (X(v) > C_L1 || X(v) < -C_L2 || fabs(Y(v)) > C_l || fabs(Z(v)) > C_l)
#define BOOL_CROSS_PICTURE(v)       (X(v) < - C_L2 && fabs(Y(v)) < C_l / 2. && fabs(Z(v)) < C_l / 2.)
#define BOOL_STOP_CONDITION(v)      (BOOL_OUTSIDE_BOUNDARY(v) || BOOL_CROSS_PICTURE(v) || BOOL_NEAR_HORIZON(v))

#endif /* blackhole_lab_h */
