


## 2D grid

function generate_grid2D_planar(;parameter_file, demo_run)

    include(parameter_file)

    h_active     = h_activePL
    heightLayers = heightLayersPL
    h_total      = h_totalPL

    maxvolume    = 5.0e-17
    eps          = 2.0e-10

    ## ETL1  C60
    heightETL1Bottom     = 0.0             + 2.0e1 * eps
    heightETL1Top        = heightLayers[1] - 1.0  * eps
    heightETL1Top2       = heightLayers[1] - 2.0  * eps
    heightETL1Top3       = heightLayers[1] - 4.0  * eps
    heightETL1Top4       = heightLayers[1] - 6.0  * eps
    heightETL1Top5       = heightLayers[1] - 10.0 * eps
    heightETL1Top6       = heightLayers[1] - 14.0 * eps

    ## pero
    heightPeroBottom     = heightLayers[1] + 1.0  * eps
    heightPeroBottom2    = heightLayers[1] + 2.0  * eps
    heightPeroBottom3    = heightLayers[1] + 4.0  * eps
    heightPeroBottom4    = heightLayers[1] + 6.0  * eps
    heightPeroBottom5    = heightLayers[1] + 10.0 * eps
    heightPeroBottom6    = heightLayers[1] + 14.0 * eps

    heightPeroTop        = heightLayers[2] - 1.0  * eps
    heightPeroTop2       = heightLayers[2] - 2.0  * eps
    heightPeroTop3       = heightLayers[2] - 4.0  * eps
    heightPeroTop4       = heightLayers[2] - 6.0  * eps
    heightPeroTop5       = heightLayers[2] - 10.0 * eps
    heightPeroTop6       = heightLayers[2] - 14.0 * eps

    heightPero1MaxBottom = 60.0  * nm - 5.0 * nm
    heightPero1MaxTop    = 60.0  * nm + 5.0 * nm
    heightPero2MaxBottom = 200.0 * nm - 5.0 * nm
    heightPero2MaxTop    = 200.0 * nm + 5.0 * nm
    heightPero3MaxBottom = 350.0 * nm - 5.0 * nm
    heightPero3MaxTop    = 350.0 * nm + 5.0 * nm

    # HTL
    heightHTLBottom      = heightLayers[2] + 1.0  * eps
    heightHTLBottom2     = heightLayers[2] + 2.0  * eps
    heightHTLBottom3     = heightLayers[2] + 4.0  * eps
    heightHTLBottom4     = heightLayers[2] + 6.0  * eps
    heightHTLBottom5     = heightLayers[2] + 10.0 * eps
    heightHTLBottom6     = heightLayers[2] + 14.0 * eps

    heightHTLTop         = heightLayers[3] - 2.0e1 * eps

    #############################################
    fPlanar(x, heightLayer) = heightLayer

    XCoarse = collect(range(0.0, heightDev, length = 10))
    XFine   = collect(range(0.0, heightDev, length = 50))
    XFine2  = collect(range(0.0, heightDev, length = 50))

    b       = SimplexGridBuilder(Generator=Triangulate)

    ###############################################################
    ## first region (SnO2)
    ###############################################################

    ## nodes
    height_0 = point!(b, heightDev, 0.0)
    height_n = point!(b, heightDev, heightLayers[1])

    length_0 = point!(b, 0.0, 0.0)
    length_n = point!(b, 0.0, heightLayers[1])

    ## facets
    facetregion!(b, bregionLeft); facet!(b, height_0,  length_0)
    facetregion!(b, 0);           facet!(b, height_n,  length_n)
    facetregion!(b, 0);           facet!(b, length_0,  length_n)
    facetregion!(b, 0);           facet!(b, height_0,  height_n)

    # ## refinement
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Bottom))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Bottom))
        facetregion!(b, 0); facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Top))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Top))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Top2))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Top2))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Top3))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Top3))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Top4))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Top4))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Top5))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Top5))
        facet!(b, A, B)
    end

    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightETL1Top6))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightETL1Top6))
        facet!(b, A, B)
    end

    ## regions
    cellregion!(b, regionETL1)
    regionpoint!(b, 0.0, 0.6*heightLayers[1])
    regionpoint!(b, 0.0, 0.99999*heightETL1Bottom)
    regionpoint!(b, 0.0, 1.00001*heightETL1Top)
    regionpoint!(b, 0.0, 1.00001*heightETL1Top2)
    regionpoint!(b, 0.0, 1.00001*heightETL1Top3)
    regionpoint!(b, 0.0, 1.00001*heightETL1Top4)
    regionpoint!(b, 0.0, 1.00001*heightETL1Top5)
    regionpoint!(b, 0.0, 1.00001*heightETL1Top6)




    ## refinement

    ###############################################################
    ## second region (perovskite)
    ###############################################################

    ## nodes
    length_ni = point!(b, 0.0, heightLayers[2])
    height_ni = point!(b, heightDev, heightLayers[2])

    ## facets
    facet!(b, length_n,  length_ni)
    facet!(b, height_n,  height_ni)
    facetregion!(b, bregionJ2); facet!(b, length_ni, height_ni)

    # last facet
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightLayers[2]))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightLayers[2]))
        facetregion!(b, bregionJ2); facet!(b, A, B)
    end

    ## refinement

    ## bottom
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroBottom))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroBottom))
        facetregion!(b, 0); facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroBottom2))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroBottom2))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroBottom3))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroBottom3))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroBottom4))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroBottom4))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroBottom5))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroBottom5))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroBottom6))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroBottom6))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroTop))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroTop))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroTop2))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroTop2))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroTop3))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroTop3))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroTop4))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroTop4))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroTop5))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroTop5))
        facet!(b, A, B)
    end

    ## top
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPeroTop6))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPeroTop6))
        facet!(b, A, B)
    end

    ## maxima

    # for ix in 2:length(XFine)
    #     A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPero1MaxBottom))
    #     B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPero1MaxBottom))
    #     facet!(b, A, B)
    # end

    # for ix in 2:length(XFine)
    #     A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPero1MaxTop))
    #     B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPero1MaxTop))
    #     facet!(b, A, B)
    # end

    # for ix in 2:length(XFine)
    #     A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPero2MaxBottom))
    #     B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPero2MaxBottom))
    #     facet!(b, A, B)
    # end

    # for ix in 2:length(XFine)
    #     A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPero2MaxTop))
    #     B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPero2MaxTop))
    #     facet!(b, A, B)
    # end

    # for ix in 2:length(XFine)
    #     A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPero3MaxBottom))
    #     B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPero3MaxBottom))
    #     facet!(b, A, B)
    # end

    # for ix in 2:length(XFine)
    #     A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightPero3MaxTop))
    #     B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightPero3MaxTop))
    #     facet!(b, A, B)
    # end

    ## regions
    cellregion!(b, regionPero)
    regionpoint!(b, 0.0, 0.9999*heightPeroBottom)
    regionpoint!(b, 0.0, 0.9999*heightPeroBottom2)
    regionpoint!(b, 0.0, 0.9999*heightPeroBottom3)
    regionpoint!(b, 0.0, 0.9999*heightPeroBottom4)
    regionpoint!(b, 0.0, 0.9999*heightPeroBottom5)
    regionpoint!(b, 0.0, 0.9999*heightPeroBottom6)
    regionpoint!(b, 0.0, 1.5*heightLayers[1])
    regionpoint!(b, 0.0, 1.0001*heightPeroTop)
    regionpoint!(b, 0.0, 1.0001*heightPeroTop)
    regionpoint!(b, 0.0, 1.0001*heightPeroTop2)
    regionpoint!(b, 0.0, 1.0001*heightPeroTop3)
    regionpoint!(b, 0.0, 1.0001*heightPeroTop4)
    regionpoint!(b, 0.0, 1.0001*heightPeroTop5)
    regionpoint!(b, 0.0, 1.0001*heightPeroTop6)

    regionpoint!(b, 0.0, 0.9999*heightPero1MaxBottom)
    regionpoint!(b, 0.0, 0.9999*heightPero1MaxTop)
    regionpoint!(b, 0.0, 0.9999*heightPero2MaxBottom)
    regionpoint!(b, 0.0, 0.9999*heightPero2MaxTop)
    regionpoint!(b, 0.0, 0.9999*heightPero3MaxBottom)
    regionpoint!(b, 0.0, 1.0001*heightPero3MaxBottom)
    regionpoint!(b, 0.0, 1.0001*heightPero3MaxTop)

    ###############################################################
    ## third region (PTAA)
    ###############################################################

    ## nodes
    length_nip  = point!(b, 0.0, heightLayers[3])
    height_nip  = point!(b, heightDev, heightLayers[3])

    ## facets
    facet!(b, length_ni, length_nip)
    facet!(b, height_ni, height_nip)

    ## boundary
    for ix in 2:length(XFine)
        A = point!(b, XFine[ix-1], fPlanar(XFine[ix-1], heightLayers[3]))
        B = point!(b, XFine[ix],   fPlanar(XFine[ix],   heightLayers[3]))
        facetregion!(b, bregionRight); facet!(b, A, B)
    end

    # refinement

    ## bottom
    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLBottom))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLBottom))
        facetregion!(b, 0); facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLBottom2))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLBottom2))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLBottom3))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLBottom3))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLBottom4))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLBottom4))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLBottom5))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLBottom5))
        facet!(b, A, B)
    end

    ## bottom
    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLBottom6))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLBottom6))
        facet!(b, A, B)
    end


    for ix in 2:length(XFine2)
        A = point!(b, XFine2[ix-1], fPlanar(XFine2[ix-1], heightHTLTop))
        B = point!(b, XFine2[ix],   fPlanar(XFine2[ix],   heightHTLTop))
        facet!(b, A, B)
    end

    ## regions
    cellregion!(b, regionHTL)
    regionpoint!(b, 0.0, 0.9999*heightHTLBottom)
    regionpoint!(b, 0.0, 0.9999*heightHTLBottom2)
    regionpoint!(b, 0.0, 0.9999*heightHTLBottom3)
    regionpoint!(b, 0.0, 0.9999*heightHTLBottom4)
    regionpoint!(b, 0.0, 0.9999*heightHTLBottom5)
    regionpoint!(b, 0.0, 0.9999*heightHTLBottom6)

    regionpoint!(b, 0.0, 1.00001*heightHTLTop)

    regionpoint!(b, 0.0, heightLayers[2] + 0.5 * h_HTL)

    ###############################################################
    ## final
    ###############################################################

    options!(b, maxvolume=maxvolume)

    grid = simplexgrid(b)

    # delete unnecessary masks with setting region to zero
    bfacemask!(grid, [0.0, 0.0],       [0.0, heightLayers[3]],       0, tol = 1.0e-20)
    bfacemask!(grid, [heightDev, 0.0], [heightDev, heightLayers[3]], 0, tol = 1.0e-20)

    # this one is needed for calculation of photogeneration
    bfacemask!(grid, [0.0, heightLayers[1]], [heightDev, heightLayers[1]], bregionJ1, tol = 1.0e-20)

    #builderplot(Plotter = PyPlot, b; resolution = (500, 500))

    return grid

end
