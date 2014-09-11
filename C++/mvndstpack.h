extern "C" void mvndst_(int *n, double lower[], double upper[], int infin[], double correl[], 
			int *maxpts, double *abseps, double *releps,
                        double *error, double *value, int *inform);

//      SUBROUTINE MVNDST( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
//     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )

// *     INFIN  INTEGER, array of integration limits flags:
// *            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
// *            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
// *            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
// *            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
