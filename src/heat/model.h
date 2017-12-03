#pragma once

#include "plot.h"

namespace model
{

    using points_t = std::vector < plot::point < double > > ;

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

    struct plot_data
    {
        util::ptr_t < std::vector < points_t > > data;
        plot::multilist_drawable < points_t > :: ptr_t plot;
        plot::world_t world;
        plot::world_mapper_t world_mapper;
    };

    inline static plot_data make_plot_data()
    {
        plot_data pd;
        pd.data = util::create < std::vector < points_t > > ();
        pd.plot = plot::multilist_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr // no point painter
        );
        pd.world_mapper = plot::make_world_mapper(pd.world);
        return pd;
    }

    inline static void adjust(parameters & params,
                              plot_data  & data)
    {
        data.world.xmin = data.world.ymin = 0;
        data.world.xmax = params.L + 2 * params.d;
        data.world.ymax = params.R + params.d;
        if (params.z_c > 0)
        {
            data.world.ymax += params.d_c;
        }
        params.dr = params.dz = min(params.d_c, params.d) / 2;
    }
}
