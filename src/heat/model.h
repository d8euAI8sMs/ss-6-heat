#pragma once

#include <util/common/plot/plot.h>

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
            1, 2.5, 0.05, 10,

            // ring params
            0.05, 0.2, 0.5,

            // heat source params
            0.5, 0.1, 0.5, 1,

            // other params
            0.01
        };
    }

    struct plot_data
    {
        util::ptr_t < std::vector < points_t > > data;
        plot::multilist_drawable < points_t > :: ptr_t plot;
        plot::world_t::ptr_t world;
        plot::world_mapper_t world_mapper;
    };

    inline static plot_data make_plot_data()
    {
        plot_data pd;
        pd.data = util::create < std::vector < points_t > > ();
        pd.plot = plot::multilist_drawable < points_t > :: create
        (
            plot::make_data_source(pd.data),
            nullptr, // no point painter
            plot::palette::pen(0xffffff, 2)
        );
        pd.world = plot::world_t::create();
        pd.world_mapper = plot::make_world_mapper(pd.world);
        return pd;
    }

    inline static void adjust(parameters & params,
                              plot_data  & data)
    {
        data.world->xmin = data.world->ymin = 0;
        data.world->xmax = params.L + 2 * params.d;
        data.world->ymax = params.R + params.d;
        if (params.z_c > 0)
        {
            data.world->ymax += params.d_c;
        }
        params.dr = params.dz = min(params.d_c, params.d) / 2;
    }

    inline static plot::drawable::ptr_t make_root_drawable
    (
        const plot_data & data,
        std::vector < plot::drawable::ptr_t > layers
    )
    {
        using namespace plot;

        layers.push_back(data.plot);

        return viewporter::create(
            tick_drawable::create(
                layer_drawable::create(layers),
                const_n_tick_factory<axe::x>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                const_n_tick_factory<axe::y>::create(
                    make_simple_tick_formatter(2, 5),
                    0,
                    5
                ),
                palette::pen(RGB(80, 80, 80)),
                RGB(200, 200, 200)
            ),
            make_viewport_mapper(data.world_mapper)
        );
    }

    inline static RECT xyxy(const plot::viewport & vp, const plot::rect < double > & r)
    {
        return
        {
            vp.world_to_screen().x(r.xmin),
            vp.world_to_screen().y(r.ymin),
            vp.world_to_screen().x(r.xmax),
            vp.world_to_screen().y(r.ymax)
        };
    }

    inline static plot::painter_t make_system_painter(parameters & params)
    {
        return [&] (CDC & dc, const plot::viewport & vp)
        {
            auto metal_brush       = plot::palette::brush(RGB(100, 100, 100));
            auto heat_source_brush = plot::palette::brush(RGB(100, 100, 100));

            RECT r;

            // heat source
            r = xyxy(vp,
            {
                params.z_h - params.L_h / 2,
                params.z_h + params.L_h / 2,
                0,
                params.R_h
            });
            dc.FillRect(&r, heat_source_brush.get());

            // ring
            if (params.z_c > 0)
            {
                r = xyxy(vp,
                {
                    params.z_c - params.L_c / 2,
                    params.z_c + params.L_c / 2,
                    params.R + params.d,
                    params.R + params.d + params.d_c
                });
                dc.FillRect(&r, metal_brush.get());
            }

            // metal

            r = xyxy(vp,
            {
                0,
                params.d,
                0,
                params.R
            });
            dc.FillRect(&r, metal_brush.get());

            r = xyxy(vp,
            {
                params.L + params.d,
                params.L + 2 * params.d,
                0,
                params.R
            });
            dc.FillRect(&r, metal_brush.get());

            r = xyxy(vp,
            {
                0,
                params.L + 2 * params.d,
                params.R,
                params.R + params.d
            });
            dc.FillRect(&r, metal_brush.get());
        };
    }
}
