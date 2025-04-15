/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File run_erfcx.c:
 *   Interactive evaluation of erfcx(x).
 *
 * Copyright:
 *   (C) 2025 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Public domain.
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2025
 *
 * Website:
 *   http://apps.jcns.fz-juelich.de/libcerf
 */

#include <stdio.h>
#include <stdlib.h>
#include "cerf.h"

int main( int argc, char **argv )
{
    if( argc!=2 ){
        fprintf( stderr,  "usage:\n" );
        fprintf( stderr,  "   run_erfcx x\n" );
        exit(-1);
    }

    double x = atof( argv[1] );

    double y = erfcx(x);
    printf( "%21.16g %21.16g\n", x, y);
    return 0;
}
