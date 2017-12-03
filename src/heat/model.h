#pragma once

namespace model
{

    struct parameters
    {
        // cylinder params
        double R, L, d, lambda_m;

        // ring params
        double d_c, L_c, z_c;

        // heat source params
        double R_h, L_h, z_h, P_h;

        // other params
        double dt, dr, dz;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
            // cylinder params
            2, 1, 0.05, 10,

            // ring params
            0.05, 0.2, 0.5,

            // heat source params
            0.5, 0.1, 0.2, 1,

            // other params
            0.01
        };
    }
}
