/* Library libcerf:
 *   Compute complex error functions, based on a new implementation of
 *   Faddeeva's w_of_z. Also provide Dawson and Voigt functions.
 *
 * File run_imwofx.c:
 *   Interactive evaluation of Im w(x), which is proportional to dawson(x).
 *
 * Copyright:
 *   (C) 2013 Forschungszentrum Jülich GmbH
 *
 * Licence:
 *   Public domain.
 *
 * Author:
 *   Joachim Wuttke, Forschungszentrum Jülich, 2013
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
        fprintf( stderr,  "   run_imwofx x\n" );
        exit(-1);
    }

    double x = atof( argv[1] );

    double y = im_w_of_x(x);
    printf( "%21.16g %21.16g\n", x, y);
    return 0;
}
