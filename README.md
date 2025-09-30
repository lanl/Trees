The Trees program is a pre-processing tool used to virtually build the forest fuel arrays used by FIRETEC and QUIC-Fire. This is meant to be a living document to be updated as the program evolves and is modified.

## Language and Compilers
In its current form, the Trees program is written in fortran and has been tested/verified using the following compilers:
- gfortran (9.2.0) with flags -O2 -ffixed-line-length=none
- intel’s fortran compiler (19.0.4) with flags -g -extend_source
- cray fortran compiler

## Getting Started
The trees program is a git remote repository available both on yellow (git.lanl.gov) and open (github.com) networks. Included in the directory called trees is the source code for the program, both a default makefile and cmake files to build (chose one or the other), and an Inputs directory which contains example input files for a basic forest build.

To run the program one can either make adjustments to the makefile for the desired compiler and flags, then run ‘make’, or use cmake ('cmake .' at the top directory, followed by 'make'). This will create an executable file called ‘trees.exe’. Move this file to the location of your input files and execute to create the four fuel files -- treesrhof.dat, treesss.dat, treesmoist.dat, and treesfueldepth.dat -- needed by FIRETEC and/or QUIC-Fire.

## Inputs
There is one main input file with an example found in the Inputs directory.

The main input file which the program will look for is called ‘fuellist’ and is parsed with variables listed in the following table. NOTE: for each variable to be correctly parsed, they need to be found in 'fuellist' as VARIABLE = VALUE (ex: nx = 200) with spaces before and after the equality and only one variable per line.
### Fuellist Main Variables
    !-------------------------------------------------------------------------------------------------
    Variable -- Default Value -- Dependency --Description
    !-------------------------------------------------------------------------------------------------
    nx          -- 200      -- None     -- Number of cells in the x-dimension
    ny          -- 200      -- None     -- Number of cells in the y-dimension
    nz          -- 41       -- None     -- Number of cells in the z-dimension
    dx          -- 2.0      -- None     -- Spacing of cells in the x-dimension (m)
    dy          -- 2.0      -- None     -- Spacing of cells in the y-dimension (m)
    dz          -- 15.0     -- None     -- Spacing of cells in the z-dimension before stretching (m)
    aa1         -- 0.1      -- None     -- Vertical stretching component defined in zcart function of metryc.f, 0 if no stretching or no     topo to vertical size
    singlefuel  -- 0        -- None     -- Flag for averaging cells' fuels into a single fuel type
    lreduced    -- 0        -- None     -- Flag for writing files only up to lfuel rather than entire domain
    topofile    -- ‘flat’   -- None     -- Name of topo file to be read in if using  topography
    !-------------------------------------------------------------------------------------------------
    ifuelin     -- 0                          -- None       --Flag for reading existing fuel data files
    inx         -- nx                         -- ifuelin=1  --Number of cells in the existing fuel data x-dimension
    iny         -- ny                         -- ifuelin=1  --Number of cells in the existing fuel data y-dimension
    inz         -- nz                         -- ifuelin=1  --Number of cells in the existing fuel data z-dimension
    idx         -- dx                         -- ifuelin=1  --Spacing of cells in the existing fuel data x-dimension (m)
    idy         -- dy                         -- ifuelin=1  --Spacing of cells in the existing fuel data y-dimension (m)
    idz         -- dz                         -- ifuelin=1  --Spacing of cells in the existing fuel data z-dimension before stretching (m)
    iaa1        -- aa1                        -- ifuelin=1  --Vertical stretching component in existing fuel data
    infuel      -- 1                          -- ifuelin=1  --Number of fuel types in existing fuel data
    intopofile  -- ‘flat’                     -- ifuelin=1  --Name of existing topo file to be read in if using topography
    rhoffile    -- 'treesrhof.dat.orig'       -- ifuelin=1  --Existing rhof file for read-in fuel data
    ssfile      -- 'treesss.dat.orig'         -- ifuelin=1  --Existing sizescale file for read-in fuel data
    moistfile   -- 'treesmoist.dat.orig'      -- ifuelin=1  --Existing moisture file for read-in fuel data
    afdfile     -- 'treesfueldepth.dat.orig'  -- ifuelin=1  --Existing fueldepth file for read-in fuel data
    !-------------------------------------------------------------------------------------------------
    igrass        -- 0      -- None         -- Flag for grass, 0 if no grass, 1 if generalizd grassfile, 2 if ground fuel levels read directly
    ngrass        -- 1      -- igrass=1     -- Number of grass species (for multiple fuels)
    grassconstant -- 5      -- igrass=1     -- Exponential constant used to determine the decay of grass mass with tree shading
    grho          -- 1.18   -- ngrass>0     -- Array of ngrass size containing the bulk fuel density of grasses
    gmoisture     -- 0.06   -- ngrass>0     -- Array of ngrass size containing the moisture content of grasses
    gss           -- 0.0005 -- ngrass>0     -- Array of ngrass size containing the sizescale of grasses
    gdepth        -- 0.27   -- ngrass>0     -- Array of ngrass size containing the height of grasses
    !-------------------------------------------------------------------------------------------------
    itrees      -- 0                      -- None       -- Flag for trees, 0 if no trees, 2 if specific treefile w/ locations     randomized to fill the domain, 3 if specific treefile w/ random locations using base area to fill the domain, 7 if using FastFuels .csv file
    tfuelbins   -- 1                      -- itrees>0   -- Number of size bins to distribute canopy fuels (foliage, branches, etc)
    treefile    -- 'AJoseEglinTrees.txt'  -- itrees>0   -- Name of treefile with the additional information to be read
    ndatax      -- nx*dx                  -- itrees>0   -- Size of dataset domain in x direction (m); recalculated if itrees=7
    ndatay      -- ny*dy                  -- itrees>0   -- Size of dataset domain in y direction (m); recalculated if itrees=7
    datalocx    -- 0.                     -- itrees>0   -- Bottom left corner x coordinate for where dataset should be placed
    datalocy    -- 0.                     -- itrees>0   -- Bottom left corner y coordinate for where dataset should be placed
    trhofmax    -- 1.5*tBulkDenisty       -- itrees>0   -- Array of size ntspecies containing the maximum tree bulk density in any given cell
    !-------------------------------------------------------------------------------------------------
    ilitter         -- 0                    -- None             -- Flag for litter, 0 if no litter, 1 if litter distributed by vertical fuel-load, 2 if litter is filled with DUET
    litterconstant  -- 5                    -- ilitter=1        -- Exponential constant to increase of litter mass under trees
    lrho            -- 1.18                 -- ilitter=1        -- Array of size ntspecies*tfuelbins containing the normal bulk litter density
    lmoisture       -- 0.06                 -- ilitter=1        -- Array of size ntspecies*tfuelbins containing the litter moisture content
    lss             -- 0.0005               -- ilitter=1        -- Array of size ntspecies*tfuelbins containing the litter sizescale
    ldepth          -- 0.05                 -- ilitter=1        -- Array of size ntspecies*tfuelbins containing the litter depth
    windprofile     -- 2                    -- ilitter=2        -- Flag for wind profile type used by DUET, 1 if yearly averaged, 2 if imported windfield
    YearsSinceBurn  -- 1                    -- ilitter=2        -- Years since last fire
    StepsPerYear    -- 1                    -- ilitter=2        -- Number of DUET timesteps within a year
    relhum          -- 0.1                  -- ilitter=2        -- relative humidity of the area
    grassstep       -- 1                    -- ilitter=2        -- Step within a year in which grass grows
    iFIA            -- 1                    -- ilitter=2        -- Use FIA database
    iFIAspecies     -- 1                    -- ilitter=2        -- Number of species from FIA database
    litout          -- 0                    -- ilitter=2        -- flag to turn on extra output files for visualizations of DUET functionality
    gmoistoverride  -- 0.0                  -- ilitter=2        -- Overwrite grass moisture content
    uavg            -- 0.                   -- windprofile=0    -- Average u velocity throughout a year
    vavg            -- 0.                   -- windprofile=0    -- Average v velocity throughout a year
    ustd            -- 0.                   -- windprofile=0    -- Standard deviation u velocity throughout a year
    vstd            -- 0.                   -- windprofile=0    -- Standard deviation v velocity throughout a year
    randomwinds     -- 2                    -- windprofile=2    -- Flag for wind generation: 1 for up to 3m/s, 2 for up to 6m/s, 3 for up to 9m/s, etc.
    speciesfile     -- 'speciesfile.dat'    -- iFIA=0           -- Provide own species database instead of FIA database 
    FIA             -- 1                    -- iFIA=1           -- Array of size infuel+ntspecies*tfuelbins containing who knows what
    !-------------------------------------------------------------------------------------------------
    verbose -- 0    -- None     -- Option to output treelist and fuellist from trees run
    !-------------------------------------------------------------------------------------------------
    !---------------------  END  OF  FUELLIST  INPUTS  -----------------------------------------------
    !-------------------------------------------------------------------------------------------------

## Trees Input Information
If itrees is not 0, then the program will look for the treefile with a name specified in the main input file. However, different itrees values will read the file in different ways.

### Trees Specific
If itrees is 2, the treefile should specific tree information the program will use to populate the domain. For this treefile, include the data for each individual tree on a new line. Each line of data should contain 10 columns of data in the following order:
1. Tree species identity, a numerical value
2. x-coordinate of tree stem [m]
3. y-coordinate of tree stem [m]
4. Tree height [m]
5. Height to canopy base [m]
6. Canopy diameter [m]
7. Height to maximum canopy diameter [m]
8. Canopy bulk density [kg/m3]
9. Canopy moisture content
10. Canopy fuel radius (sizescale) [m]

If itrees is 3, you use the same file as for itrees equals 2, but the trees will randomly re-distributed instead of using the provided coordinates.

## Outputs
This program creates four binary files each containing a full fortran array which can be read directly by FIRETEC or QUIC-Fire:
1. treesrhof.dat contains bulk fuel density for entire HIGRAD array
2. treesfueldepth.dat contains actual fuel depths for entire HIGRAD array (3D array, but only bottom layer actually used in fire simulations)
3. treesss.dat contains sizescales for entire HIGRAD array
4. treesmoist.dat contains fuel moisture contents for entire HIGRAD array

## Program Outline
Main driver file for this program is located in main.f. From this file, all other subroutines are called which execute in this order:
1. fuellist_input (io.F90) is called to read input data and initialize variables used throughout the simulation
2. define_constant variables (define_variables.F90) defines all additional variables needed for the simulation and not specified by the user within input files.
3. define_grid_variables (define_variables.F90) initializes all data arrays needed throughout the simulation
4. grid_readin (fuel_read.F90) will readin existing trees.dat files if called
4. fuels_create (fuels_create.F90) first fills the grass arrays, then the trees arrays, and finally the litter arrays while modifying the other arrays as needed
5. output (io.F90) writes the .dat files and finalizes any simulation variables

## Building and Executing
A makefile is provided for building this software. Modify as needed for system and execute 'make' from the commandline in the location where the makefile is located. Alternatively, the program is set up for cmake as well. At the top level directory simply execute 'cmake .' followed by 'make'.
