#pragma once

#include <array>

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

        // material params
        double c_l, k_l, c_m, k_m, c_h, k_h;

        // ext params
        double sigma;
    };

    inline static parameters make_default_parameters()
    {
        return
        {
            // cylinder params
            50, 50, 15, 1,

            // ring params
            5, 20, 50,

            // heat source params
            30, 15, 50, 1000,

            // other params
            0.3, 2.5, 2.5,

            // material params
            1000, 0.6, 450, 92, 450, 12,

            // ext params
            1e-8
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
        if ((params.z_c != 0) && (params.L_c != 0) && (params.d_c != 0))
        {
            data.world->ymax += params.d_c;
        }
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

    struct chasing_data;

    enum class chasing_dir
    {
        i, j
    };

    struct chasing_coefs
    {
        double A, B, C, D;
    };

    using chasing_fn = std::function < chasing_coefs (size_t i, size_t j, chasing_dir dir,
                                                      const chasing_data & d,
                                                      const parameters & p,
                                                      std::vector < std::vector < double > > & T) > ;

    using material_t = size_t;

    namespace material
    {
        const material_t border_i = 0x1 << 0;
        const material_t border_j = 0x1 << 1;
        const material_t ext      = 0x1 << 4;
        const material_t metal    = 0x1 << 5;
        const material_t heater   = 0x1 << 6;
        const material_t liquid   = 0x1 << 7;

        const material_t border   = border_i | border_j;
        const material_t flags    = border;
        const material_t no_flags = ~flags;
    };

    struct chasing_data
    {
        std::vector < std::vector < material_t > > area_map;
        std::vector < std::vector < double > > heat_src;
        chasing_fn coefs;
        size_t n, m;
    };

    inline material_t get_material_at(const chasing_data & d,
                                      const plot::point < int > & p,
                                      bool preserve_flags = true)
    {
        return ((p.x < 0) || (p.x >= d.n) || (p.y < 0) || (p.y >= d.m))
            ? material::ext
            : (preserve_flags ? d.area_map[p.x][p.y]
                              : (d.area_map[p.x][p.y] & material::no_flags));
    }

    inline std::array < material_t, 2 > get_nearest_materials(const chasing_data & d,
                                                              const plot::point < int > & p,
                                                              chasing_dir dir)
    {
        material_t m1, m2;

        if (dir == chasing_dir::i)
        {
            return {{ get_material_at(d, { p.x - 1, p.y }), get_material_at(d, { p.x + 1, p.y }) }};
        }
        else
        {
            return {{ get_material_at(d, { p.x, p.y - 1 }), get_material_at(d, { p.x, p.y + 1 }) }};
        }
    }

    inline std::pair < double, double > get_material_props(const parameters & p,
                                                           material_t m)
    {
        switch (m & material::no_flags)
        {
        case model::material::metal:
            return { p.c_m, p.k_m };
        case model::material::heater:
            return { p.c_h, p.k_h };
        case model::material::liquid:
            return { p.c_l, p.k_l };
        default:
            return { 0, 0 };
        }
    }

    inline bool is_in_rect(const plot::point < size_t > & p,
                           const plot::rect < size_t > & r)
    {
        return (r.xmin <= p.x) && (p.x <= r.xmax)
            && (r.ymin <= p.y) && (p.y <= r.ymax);
    }

    inline void make_chasing_data(chasing_data & d, const parameters & p)
    {
        bool has_ring = (p.z_c != 0) && (p.d_c != 0) && (p.L_c != 0);

        d.n = (size_t) std::ceil((p.R + p.d + (has_ring ? p.d_c : 0)) / p.dr) + 1;
        d.m = (size_t) std::ceil((p.L + 2 * p.d) / p.dz) + 1;

        d.area_map.clear();
        d.heat_src.clear();

        d.area_map.resize(d.n, std::vector < material_t > (d.m));
        d.heat_src.resize(d.n, std::vector < double > (d.m));

        size_t R_n  = (size_t) std::ceil(p.R / p.dr);
        size_t d_n  = (size_t) std::ceil(p.d / p.dr);
        size_t d_m  = (size_t) std::ceil(p.d / p.dz);
        size_t h_j1 = (size_t) std::ceil((p.z_h - p.L_h / 2) / p.dz);
        size_t h_j2 = (size_t) std::ceil((p.z_h + p.L_h / 2) / p.dz);
        size_t h_n  = (size_t) std::ceil(p.R_h / p.dr);
        size_t c_j1 = (!has_ring) ? 0 : (size_t) std::ceil((p.z_c - p.L_c / 2) / p.dz);
        size_t c_j2 = (!has_ring) ? 0 : (size_t) std::ceil((p.z_c + p.L_c / 2) / p.dz);
        size_t c_n  = (!has_ring) ? 0 : (size_t) std::ceil(p.d_c / p.dr);

        // set up materials and sketch out borders

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                bool is_border        = ((i == 0) || (i == (R_n + d_n)) && !((c_j1 < j) && (j < c_j2))
                                      || ((j == 0) || ((j + 1) == d.m))) && (i <= R_n + d_n)
                                      || (i == (d.n - 1)) && (((c_j1 < j) && (j < c_j2)))
                                      || ((j == c_j1) || (j == c_j2)) && (((R_n + d_n) <= i) && (i <= (d.n - 1)));
                bool is_heater        = is_in_rect({ i, j }, { 0, h_n, h_j1, h_j2 });
                bool is_heater_border = is_heater && !is_in_rect({ i, j }, { 1, h_n - 1, h_j1 + 1, h_j2 - 1 });
                bool is_liquid        = is_in_rect({ i, j }, { 0, R_n, d_m, d.m - d_m - 1 });
                bool is_liquid_border = is_liquid && !is_in_rect({ i, j }, { 1, R_n - 1, d_m + 1, d.m - d_m - 2 });
                bool is_ring          = is_in_rect({ i, j }, { R_n + d_n, d.n - 1, c_j1, c_j2 });
                bool is_metal         = is_in_rect({ i, j }, { 0, R_n + d_n, 0, d.m - 1 });

                if (is_border || is_heater_border || is_liquid_border)
                                    d.area_map[i][j]  = material::border;
                     if (is_heater) d.area_map[i][j] |= material::heater;
                else if (is_liquid) d.area_map[i][j] |= material::liquid;
                else if (is_metal)  d.area_map[i][j] |= material::metal;
                else if (is_ring)   d.area_map[i][j] |= material::metal;
                else                d.area_map[i][j]  = material::ext;

                if (is_heater) d.heat_src[i][j] = p.P_h;
            }
        }

        // detect border orientations

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                if (get_material_at(d, { (int) i, (int) j }) & material::border)
                {
                    d.area_map[i][j] &= ~material::border;
                    if ((get_material_at(d, { (int) i - 1, (int) j })
                        | get_material_at(d, { (int) i + 1, (int) j })) & material::border)
                    {
                        d.area_map[i][j] |= material::border_j;
                    }
                    if ((get_material_at(d, { (int) i, (int) j - 1 })
                        | get_material_at(d, { (int) i, (int) j + 1 })) & material::border)
                    {
                        d.area_map[i][j] |= material::border_i;
                    }
                }
            }
        }

        d.coefs = [&] (size_t i, size_t j, chasing_dir dir,
                       const chasing_data & d,
                       const parameters & p,
                       std::vector < std::vector < double > > & T) -> chasing_coefs
        {
            // make ordinary A B C D coefficients
            auto make_normal_coefs = [&] (
                size_t i, size_t j, chasing_dir s,
                double c, double k) -> chasing_coefs
            {
                double lambda = std::sqrt(k * p.c_l / p.k_l / c) * p.lambda_m;

                bool is_ext_i = (get_material_at(d, { (int) i - 1, ( int) j })
                                 | get_material_at(d, { (int) i + 1, ( int) j })) & material::ext;
                bool is_ext_j = (get_material_at(d, { (int) i, ( int) j - 1 })
                                 | get_material_at(d, { (int) i, ( int) j + 1 })) & material::ext;

                if (s == chasing_dir::i)
                {
                    return
                    {
                        - lambda * lambda * p.dt / 2. * (1 - 1. / 2. / i) / p.dr / p.dr,
                        - lambda * lambda * p.dt / 2. * (1 + 1. / 2. / i) / p.dr / p.dr,
                        1 + lambda * lambda * p.dt / 1. / p.dr / p.dr,
                        T[i][j] + lambda * lambda * p.dt / 2. *
                        (
                            (is_ext_j ? 0 : ((T[i][j + 1] + T[i][j - 1] - 2 * T[i][j]) / p.dz / p.dz)) +
                            // for now, examine it later
                            (is_ext_i ? 0 : (0 * (
                                (T[i + 1][j] + T[i - 1][j] - 2 * T[i][j]) / p.dr / p.dr +
                                (T[i + 1][j] - T[i - 1][j]) / 2. / i / p.dr / p.dr
                            )))
                        ) + p.dt / 2. * d.heat_src[i][j]
                    };
                }
                return
                {
                    - lambda * lambda * p.dt / 2. / p.dz / p.dz,
                    - lambda * lambda * p.dt / 2. / p.dz / p.dz,
                    1 + lambda * lambda * p.dt / 1. / p.dz / p.dz,
                    T[i][j] + lambda * lambda * p.dt / 2. *
                    (
                        (is_ext_i ? 0 : ((T[i + 1][j] + T[i - 1][j] - 2 * T[i][j]) / p.dr / p.dr +
                        (T[i + 1][j] - T[i - 1][j]) / 2. / i / p.dr / p.dr)) +
                        // for now, examine it later
                        (is_ext_j ? 0 : (0 * ((T[i][j + 1] + T[i][j - 1] - 2 * T[i][j]) / p.dz / p.dz)))
                    ) + p.dt / 2. * d.heat_src[i][j]
                };
            };

            // deal with boundary conditions
            if (d.area_map[i][j] & material::border)
            {
                // determine border orientations
                bool vertical   = d.area_map[i][j] & material::border_j;
                bool horizontal = d.area_map[i][j] & material::border_i;
                bool corner     = vertical && horizontal;

                // detect materials the border divides
                auto nearest = get_nearest_materials(d, { (int) i, (int) j }, dir);

                material_t m;

                if (dir == chasing_dir::i)
                {
                    // chasing direction is orthogonal to the border orientation
                    if (horizontal)
                    {
                        // handle external boundary conditions
                        if (i == 0)
                        {
                            // radial symmetry of the system
                            return { 0, -1, 1, 0 };
                        }
                        else if (((m = get_material_at(d, { (int) i - 1, (int) j })) & material::ext)
                            || corner && (std::get<0>(nearest) & material::ext))
                        {
                            // radiation
                            return { 0, -1, 1, - p.sigma * T[i][j] * T[i][j] * T[i][j] * T[i][j] };
                        }
                        else if (((m = get_material_at(d, { (int) i + 1, (int) j })) & material::ext)
                                 || corner && (std::get<1>(nearest) & material::ext))
                        {
                            // radiation
                            return { -1, 0, 1, - p.sigma * T[i][j] * T[i][j] * T[i][j] * T[i][j] };
                        }
                        // handle internal boundary conditions
                        else
                        {
                            auto p1 = get_material_props(p, std::get<0>(nearest));
                            auto p2 = get_material_props(p, std::get<1>(nearest));
                            double k = p1.second / p2.second;
                            return { -1, -k, (k + 1), 0 };
                        }
                    }
                    // chasing direction coincides with the border orientation,
                    // apply the ordinary equation instead of boundary conditions
                    else
                    {
                        material_t m1 = get_material_at(d, { (int) i, (int) j - 1 });
                        material_t m2 = get_material_at(d, { (int) i, (int) j + 1 });
                        auto p1 = get_material_props(p, ((m1 != material::ext) ? m1 : m2));
                        return make_normal_coefs(i, j, dir, p1.first, p1.second);
                    }
                }
                else
                {
                    // chasing direction is orthogonal to the border orientation
                    if (vertical)
                    {
                        // handle external boundary conditions
                        if (((m = get_material_at(d, { (int) i, (int) j - 1 })) & material::ext)
                            || corner && (std::get<0>(nearest) & material::ext))
                        {
                            // radiation
                            return { 0, -1, 1, - p.sigma * T[i][j] * T[i][j] * T[i][j] * T[i][j] };
                        }
                        else if (((m = get_material_at(d, { (int) i, (int) j + 1 })) & material::ext)
                                 || corner && (std::get<1>(nearest) & material::ext))
                        {
                            // radiation
                            return { -1, 0, 1, - p.sigma * T[i][j] * T[i][j] * T[i][j] * T[i][j] };
                        }
                        // handle internal boundary conditions
                        else
                        {
                            auto p1 = get_material_props(p, std::get<0>(nearest));
                            auto p2 = get_material_props(p, std::get<1>(nearest));
                            double k = p1.second / p2.second;
                            return { -1, -k, (k + 1), 0 };
                        }
                    }
                    // chasing direction coincides with the border orientation,
                    // apply the ordinary equation instead of boundary conditions
                    else
                    {
                        material_t m1 = get_material_at(d, { (int) i - 1, (int) j });
                        material_t m2 = get_material_at(d, { (int) i + 1, (int) j });
                        auto p1 = get_material_props(p, (!(m1 & material::ext) ? m1 : m2));
                        return make_normal_coefs(i, j, dir, p1.first, p1.second);
                    }
                }
            }
            // deal with ordinary space
            else if (!(get_material_at(d, { (int) i, (int) j }) & material::ext))
            {
                auto p1 = get_material_props(p, get_material_at(d, { (int) i, (int) j }));
                return make_normal_coefs(i, j, dir, p1.first, p1.second);
            }
            // will not anyhow be used
            else
            {
                return { 0, 0, 1, 0 };
            }
        };
    }

    // post-apply boundary conditions again
    inline void chasing_bc
    (
        const chasing_data & d,
        const parameters & p,
        std::vector < std::vector < double > > & T
    )
    {
        chasing_coefs c;

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                if (d.area_map[i][j] & material::border_i)
                {
                    if ((get_material_at(d, { (int) i - 1, (int) j })
                        | get_material_at(d, { (int) i + 1, (int) j })) & material::ext)
                    {
                        c = d.coefs(i, j, chasing_dir::i, d, p, T);
                        T[i][j] = c.D;
                        if (!(get_material_at(d, { (int) i - 1, (int) j }) & material::ext))
                            T[i][j] -= c.A * T[i - 1][j];
                        if (!(get_material_at(d, { (int) i + 1, (int) j }) & material::ext))
                            T[i][j] -= c.B * T[i + 1][j];
                        T[i][j] /= c.C;
                    }
                    else if ((get_material_at(d, { (int) i, (int) j - 1 })
                        | get_material_at(d, { (int) i, (int) j + 1 })) & material::ext)
                    {
                        c = d.coefs(i, j, chasing_dir::j, d, p, T);
                        T[i][j] = c.D;
                        if (!(get_material_at(d, { (int) i, (int) j - 1 }) & material::ext))
                            T[i][j] -= c.A * T[i][j - 1];
                        if (!(get_material_at(d, { (int) i, (int) j + 1 }) & material::ext))
                            T[i][j] -= c.B * T[i][j + 1];
                        T[i][j] /= c.C;
                    }
                }
            }
        }
    }

    inline void chasing_solve
    (
        const chasing_data & d,
        const parameters & p,
        std::vector < std::vector < double > > & T
    )
    {
        std::vector < std::vector < double > > a(d.n + 1), b(d.n + 1);

        for (size_t i = 0; i < d.n + 1; ++i)
        {
            a[i].resize(d.m + 1); b[i].resize(d.m + 1);
        }

        chasing_coefs c;

        // i-chasing

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                c = d.coefs(i, j, chasing_dir::i, d, p, T);
                a[i + 1][j] = - c.B / (c.C + c.A * a[i][j]);
                b[i + 1][j] = (c.D - c.A * b[i][j]) / (c.C + c.A * a[i][j]);
            }
        }

        for (size_t j = 0; j < d.m; ++j)
        {
            size_t s = d.n - 1, e = 0;
            while (get_material_at(d, { (int) s, (int) j }) & material::ext) --s;
            while (get_material_at(d, { (int) e, (int) j }) & material::ext) ++e;
            for (size_t i = s + 2; i-- > e + 1;)
            {
                if (get_material_at(d, { (int) i, (int) j }) & material::ext)
                {
                    T[i - 1][j] = b[i][j];
                }
                else
                {
                    T[i - 1][j] = a[i][j] * T[i][j] + b[i][j];
                }
            }
        }

        // boundary condition chasing

        chasing_bc(d, p, T);

        // j-chasing

        for (size_t i = 0; i < d.n; ++i)
        {
            for (size_t j = 0; j < d.m; ++j)
            {
                c = d.coefs(i, j, chasing_dir::j, d, p, T);
                a[i][j + 1] = - c.B / (c.C + c.A * a[i][j]);
                b[i][j + 1] = (c.D - c.A * b[i][j]) / (c.C + c.A * a[i][j]);
            }
        }

        for (size_t i = 0; i < d.n; ++i)
        {
            size_t s = d.m - 1, e = 0;
            while (get_material_at(d, { (int) i, (int) s }) & material::ext) --s;
            while (get_material_at(d, { (int) i, (int) e }) & material::ext) ++e;
            for (size_t j = s + 2; j-- > e + 1;)
            {
                if (get_material_at(d, { (int) i, (int) j }) & material::ext)
                {
                    T[i][j - 1] = b[i][j];
                }
                else
                {
                    T[i][j - 1] = a[i][j] * T[i][j] + b[i][j];
                }
            }
        }

        // boundary condition chasing

        chasing_bc(d, p, T);
    }

    using stencil_fn = std::function < bool (int i, int j) > ;

    inline static stencil_fn make_simple_stencil(size_t n, size_t m)
    {
        return [=] (int i, int j) { return (i >= 0) && (j >= 0) && (i < n) && (j < m); };
    }

    inline static stencil_fn make_material_based_stencil(const chasing_data & d)
    {
        return [&] (int i, int j) { return !(get_material_at(d, { i, j }) & material::ext); };
    }

    void find_isolines
    (
        const std::vector < std::vector < double > > & T,
        double dT,
        std::vector < std::vector < plot::point < double > > > & out,
        size_t n, size_t m,
        const parameters & p,
        stencil_fn stencil,
        size_t max_isolines = 100,
        size_t max_points_in_stack = 100000
    );

    inline static plot::painter_t make_system_painter(const parameters & params,
                                                      const chasing_data & d,
                                                      const std::vector < std::vector < double > > & T,
                                                      double & max_heat_value,
                                                      bool & draw_heat_map)
    {
        return [&] (CDC & dc, const plot::viewport & vp)
        {
            if (d.area_map.empty()) return;

            auto metal_brush       = plot::palette::brush(RGB(100, 100, 100));
            auto heat_source_brush = plot::palette::brush(RGB(100, 100, 100));
            auto border_brush      = plot::palette::brush(RGB( 50,  50,  50));
            auto heat_brush        = plot::palette::brush();

            RECT r;

            material_t m;

            for (size_t i = 0; i < d.n; ++i)
            {
                for (size_t j = 0; j < d.m; ++j)
                {
                    m = get_material_at(d, { (int) i, (int) j });

                    if (m & material::ext) continue;

                    if (m & (material::border | material::metal | material::heater))
                    {
                        r = xyxy(vp,
                        {
                            ((double)j - 0.5) * params.dz,
                            ((double)j + 0.5) * params.dz,
                            ((double)i - 0.5) * params.dr,
                            ((double)i + 0.5) * params.dr
                        });

                        if (m & material::border)
                        {
                            dc.FillRect(&r, border_brush.get());
                        }
                        else if (m & material::metal)
                        {
                            dc.FillRect(&r, metal_brush.get());
                        }
                        else if (m & material::heater)
                        {
                            dc.FillRect(&r, heat_source_brush.get());
                        }
                    }

                    if (draw_heat_map)
                    {
                        r = xyxy(vp,
                        {
                            ((double)j - 0.25) * params.dz,
                            ((double)j + 0.25) * params.dz,
                            ((double)i - 0.25) * params.dr,
                            ((double)i + 0.25) * params.dr
                        });

                        heat_brush->DeleteObject();
                        heat_brush->CreateSolidBrush(RGB(155 * max(0, min(1, T[i][j] / max_heat_value)), 0, 0));

                        dc.FillRect(&r, heat_brush.get());
                    }
                }
            }
        };
    }
}
