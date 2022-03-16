The Trees program is a pre-processing tool used to virtually build forest fuel arrays to be read into FIRETEC or QUIC-Fire. This is meant to be a living document to be updated as the program evolves and is modified.


!---Language and Compilers---!
In its current form, the Trees program is written in fortran and has been tested/verified using the following compilers:
- gfortran (9.2.0) with flags -O2 -ffixed-line-length=none
- intel’s fortran compiler (19.0.4) with flags -g -extend_source

!---Getting Started---!
The trees program is attached to the FIRETEC repository held on the turquoise HPC system HIGRAD project space. Will be included on branch 1.4.1 and later. Included in the directory called trees is the source code for the program, a makefile to build, and a Inputs directory which contains example input files for a basic forest build.

To run the program make adjustments to the makefile for the desired compiler and flags, then simply run ‘make’. This will create an executable file called ‘trees’. Move this file to the location of your input files and execute to create the four fuel files: treesrhof.dat, treesss.dat, treesmoist.dat, and treesfueldepth.dat.


!---Inputs---!
There is one main input file and additional input files depending on the configurations of the main input file, and examples of each input file are found in the Inputs directory.

The main input file which the program will look for is called ‘fuellist’ and is a fortran namelist with variables listed in the following table.

!-------------------------------------------------------------------------------------------------
Variable -- Default Value -- Description
!-------------------------------------------------------------------------------------------------
nx            -- Required             -- Number of cells in the x-dimension
ny            -- Required             -- Number of cells in the y-dimension
nz            -- Required             -- Number of cells in the z-dimension
dx            -- Required             -- Spacing of cells in the x-dimension (m)
dy            -- Required             -- Spacing of cells in the y-dimension (m)
dz            -- Required             -- Spacing of cells in the z-dimension before stretching (m)
aa1           -- 0.1                  -- Vertical stretching component defined in zcart function of metryc.f, 0 if no stretching or no topo to vertical size
singlefuel    -- 0                    -- Flag for averaging cells' fuels into a single fuel type
topofile      -- ‘flat’               -- Name of topo file to be read in if using  topography
!-------------------------------------------------------------------------------------------------
ifuelin       -- 0                    -- Flag for reading existing fuel data files
inx           -- nx                   -- Number of cells in the existing fuel data x-dimension
iny           -- ny                   -- Number of cells in the existing fuel data y-dimension
inz           -- nz                   -- Number of cells in the existing fuel data z-dimension
idx           -- dx                   -- Spacing of cells in the existing fuel data x-dimension (m)
idy           -- dy                   -- Spacing of cells in the existing fuel data y-dimension (m)
idz           -- dz                   -- Spacing of cells in the existing fuel data z-dimension before stretching (m)
iaa1          -- aa1                  -- Vertical stretching component in existing fuel data
infuel        -- 1                    -- Number of fuel types in existing fuel data
rhoffile      -- Required if ifuelin=1-- Existing rhof file for read-in fuel data
ssfile        -- Required if ifuelin=1-- Existing sizescale file for read-in fuel data
moistfile     -- Required if ifuelin=1-- Existing moisture file for read-in fuel data
afdfile       -- Required if ifuelin=1-- Existing fueldepth file for read-in fuel data
!-------------------------------------------------------------------------------------------------
igrass        -- 0                    -- Flag for grass, 0 if no grass, 1 if generalizd grassfile, 2 if ground fuel levels read directly
grassconstant -- 5                    -- Exponential constant used to determine the decay of grass mass with tree shading
ngrass        -- Required if igrass=1 -- Number of grass species (for multiple fuels)
grassfile     -- Required if igrass=2 -- Name of grassfile with the additional information to be read
!-------------------------------------------------------------------------------------------------
itrees        -- 0                    -- Flag for trees, 0 if no trees, 1 if generalized treefile, 2 if specific treefile w/ locations randomized to fill the domain, 3 if specific treefile w/ random locations using base area to fill the domain
ntspecies     -- 1                    -- Number of tree species
tfuelbins     -- 1                    -- Number of size bins to distribute canopy fuels (foliage, branches, etc)
treefile      -- ‘ ‘                  -- Name of treefile with the additional information to be read
istem         -- 0                    -- Adds tree stems and bark to fuel arrays, DO NOT USE!!!!
ndatax        -- nx*dx                -- Size of dataset domain in x direction (m)
ndatay        -- ny*dx                -- Size of dataset domain in y direction (m)
datalocx      -- 0                    -- Bottom left corner x coordinate for where dataset should be placed
datalocy      -- 0                    -- Bottom left corner y coordinate for where dataset should be placed
tdnx          -- [0,nx*dx]            -- x range of the domain to be populated with trees (m)
tdny          -- [0,ny*dy]            -- y range of the domain to be populated with trees (m)
!-------------------------------------------------------------------------------------------------
ilitter       -- 0                    -- Flag for litter, 0 if no litter, 1 if litter distributed by vertical fuel-load
litterconstant-- 5                    -- Exponential constant to increase of litter mass under trees
litterfile    -- ‘ ‘                  -- Name of litterfile with additional information to be read
!-------------------------------------------------------------------------------------------------
itreatment    -- 0                    -- Flag for fuel treatments, 0 if no treatment, 1 if slashlayer, 2 if slashpiles, 3 if fuel removal
sdnx          -- [0,nx*dx]            -- x range of the domain over which a treatment is applied (m)
sdny          -- [0,ny*dy]            -- y range of the domain over which a treatment is applied (m)
sdepth        -- 0                    -- Depth of slashlayer
smoist        -- 0                    -- Moisture content of slashlayers or slashpiles
sdiameter     -- 0                    -- Diameter of slashpiles
sheight       -- 0                    -- Height of slashpiles
sprho         -- 0                    -- Bulk density of slashpiles

!---Grasses---
If igrass is 1, then the program will look for the grassfile with a name specified in the main input file. The first four rows of this file will be read by the program and should contain the following:
1. Fuel bulk density within the fuel layer (kg/m3)
2. Fuel moisture fraction
3. Radius of fuel cylinders (m)
4. Depth of the fuel layer (m)
Each line should contain ngrass columns with the above values for each grass species tracked in the simulation.

If igrass is 2, then the program will read in an existing data array to populate the ground fuels.

!---Trees---
If itrees is not 0, then the program will look for the treefile with a name specified in the main input file. However, different itrees values will read the file in different ways.

If itrees is 1, the treefile should be a ‘generalized’ file with general species information that the program will use to compute the total number of trees for each tracked species and recreate these trees based on a normal distribution. For this treefile, the first 9+3*tfuelbins lines will be read by the program and should contain the following:
1. frequency of trees [stem/ha]
2. average tree height [m]
3. standard deviation of tree height [m]
4. average height to live crown [m]
5. standard deviation height to live crown [m]
6. average crown diameter [m]
7. standard deviation crown diameter [m]
8. average height to maximum crown diameter [m]
9. standard deviation to maximum crown diameter [m]
10.  tfuel1 bulk density [kg/m3]
11.  tfuel1 moisture content [fraction]
12.  tfuel1 sizescale [m]
13.  tfuel2 bulk density [kg/m3]
14.  tfuel2 bin moisture content [unitless]
15.  tfuel2 bin sizescale [m]
16.  … and so on

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

If itrees is 7, use a FastFuels (FF) dataset. This code will take the FF data and calculate the 10 parameters listed in itrees equals 2. When using a FF dataset, do not alter the csv column order after downloading from FF website. In fuellist ndatax, ndatay, datalocx, datalocy, ntspecies and tfuelbins will be overwritten to match the FF data. 
(JO - As of 10/19/2021) The FF bounding box is not completely filled by data you expect (ex: a 400 x 400 m box will actually return a 350 x 450 rectangle domain). The ndatax, ndatay, datalocx, datalocy values are computed to be the FF dataset bounds and (0,0) respectively. Currently no moisture specified in FF data, so set all trees to 100% moisture content. FF data has 19 columns currently: 
1. sp - not used ; species number, Please see Jenkins et al. 2003 for a table of species codes with latin and common names. https://www.fs.fed.us/ne/newtown_square/publications/other_publishers/OCR/ne_2003jenkins01.pdf
2. dia - not used ; Diameter at Breast Height. Diameter of the tree in inches measured at 4.5 feet.
3. ht - tree height [ft]
4. crown_ratio - not used ; ratio of crown length to total tree height
5. sp_grp - species group number; see sp paper link for discription
6. mu - not used ; Mean of radial canopy distribution (as a proportion of crow radius)
7. sigma - not used ; Standard Deviation of radial canopy distribution (as proportion of crown radius) 
8. sav - used to find size scale [1/m]
9. crown_len - crown length [ft]
10. crown_base_ht - height to live crown [ft]
11. beta_a - beta params, used for finding height to max crown radius
12. beta_b - beta params, used for finding height to max crown radius
13. beta_c - beta params, used for finding height to max crown radius
14. beta_norm - beta params, used for finding height to max crown radius
15. weight - crown weight [kg]
16. volume - crown volume [m^3]
17. crad - maximum crown radius [m]
18. x - location coordinates from bounding box center [m]
19. y - location coordinates from bounding box center [m]

!---Outputs---
This program creates four binary files each containing a full fortran array which can be read directly by FIRETEC or QUIC-Fire:
1. treesrhof.dat contains bulk fuel density for entire HIGRAD array
2. treesfueldepth.dat contains actual fuel depths for entire HIGRAD array (3D array, but only bottom layer actually used in fire simulations)
3. treesss.dat contains sizescales for entire HIGRAD array
4. treesmoist.dat contains fuel moisture contents for entire HIGRAD array

!---Program Outline---
Main driver file for this program is located in main.f. From this file, all other subroutines are called which execute in this order:
1. namelist_input (io.f) is called to read input data and initialize variables used throughout the simulation
2. define_constant variables (define_variables.f) defines all additional variables needed for the simulation and not specified by the user within input files.
3. define_grid_variables (define_variables.f) initializes all data arrays needed throughout the simulation
4. baseline (baseline.f) first fills the grass arrays, then the trees arrays, and finally the litter arrays while modifying the other arrays as needed
5. treatment (treatments.f) executes any specified treatments to the base forest estabilished in the baseline function
6. output (io.f) writes the .dat files and finalizes any simulation variables

!---Modifications---
Under-development so modifications happening all the time and not recorded here for now.

!---Building and Executing---
A makefile is provided for building this software. Simply execute 'make' from the commandline in the location where the makefile is located. Specific configurations unique to a users build should be specified there.
There are a number of possible input configurations available to the user but in basic when executed, the program will look for a file called 'fuellist' in the location from which the executable is called. This fuellist contains input parameters are specified. Any additional input files (ie grassfile, treefile, or litterfile) are specified within the fuellist. The four output tree files will be written to the location from which the exectuable is called.
