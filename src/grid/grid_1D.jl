
## 1D grid

function generate_grid1D(;parameter_file, demo_run)

    include(parameter_file)

    t          = 1.0e-13
    ######################

    h_active     = h_activePL
    heightLayers = heightLayersPL
    h_total      = h_totalPL

    if demo_run
        ## coarser grid (1113 nodes)
        h1 = 4.0e-3
        h2 = 3.0e-3
        h3 = 3.0e-3
       h4 = 3.0e-3
    else
        ## fine solution was calculated with (3834 nodes)
        h1 = 1.0e-3
        h2 = 1.0e-3
        h3 = 9.0e-4
        h4 = 1.0e-3
    end

    coord_n1_1 = geomspace(0.0, h_ETL1/2, h2 * h_ETL1, h3 * h_ETL1, tol=t)
    coord_n1_2 = geomspace(h_ETL1/2, h_ETL1, h3 * h_ETL1, h1 * h_ETL1, tol=t)

    coord_Ac_1 = geomspace(h_ETL1, h_ETL1 + 0.15 * h_active,
                            h4 * h_ETL1, h1 * h_active, tol=t)
    coord_Ac_2 = geomspace(h_ETL1 + 0.15 * h_active, h_ETL1 + 0.50 * h_active,
                            h1 * h_active, h3 * h_active, tol=t)
    coord_Ac_3 = geomspace(h_ETL1 + 0.50 * h_active, h_ETL1 + 0.85 * h_active,
                             h3 * h_active, h1 * h_active, tol=t)
    coord_Ac_4 = geomspace(h_ETL1 + 0.85 * h_active, h_ETL1 + 1.0 * h_active,
                            h1 * h_active, h4 * h_HTL, tol=t)

    coord_p_1  = geomspace(h_ETL1 + h_active, h_ETL1 + h_active + h_HTL/2,
                            h4 * h_HTL, h3 * h_HTL, tol=t)
    coord_p_2  = geomspace(h_ETL1 + h_active + h_HTL/2, h_total,
                            h3 * h_HTL, h4 * h_HTL, tol=t)


    coordn1    = glue(coord_n1_1, coord_n1_2, tol=t)
    coordAc    = glue(coord_Ac_1, coord_Ac_2, tol=t)
    coordAc    = glue(coordAc,    coord_Ac_3, tol=t)
    coordAc    = glue(coordAc,    coord_Ac_4, tol=t)
    coordp     = glue(coord_p_1,  coord_p_2,  tol=t)

    coord      = glue(coordn1,    coordAc,    tol=t)
    coord      = glue(coord,      coordp,     tol=t)
    grid       = simplexgrid(coord)

    ## cellmask! for setting different inner regions
    cellmask!(grid,  [0.0],             [heightLayers[1]], regionETL1,   tol=t) # n-doped region = 1
    cellmask!(grid,  [heightLayers[1]], [heightLayers[2]], regionPero,   tol=t) # pero region    = 2
    cellmask!(grid,  [heightLayers[2]], [heightLayers[3]], regionHTL,    tol=t) # p-doped region = 3

    ## bfacemask! for setting different boundary regions
    bfacemask!(grid, [0.0],             [0.0],             bregionLeft,  tol=t) # outer left boundary
    bfacemask!(grid, [h_total],         [h_total],         bregionRight, tol=t) # outer right boundary

    bfacemask!(grid, [heightLayers[1]], [heightLayers[1]], bregionJ1,    tol=t) # left pero boundary
    bfacemask!(grid, [heightLayers[2]], [heightLayers[2]], bregionJ2,    tol=t) # right pero boundary

    return grid

end