! This program was created to take a FF input file and create a DUET input file

subroutine FFmakeDuetInputs

integer:: nx,ny,nz
real:: dx,dy,dz
integer:: seedchange
real:: winddirection,windvary

integer :: inputprogram=2, iFIA=1, iFIAspecies=1, windprofile=2
integer :: randomwinds=3,YearsSinceBurn=5,StepsPerYear=1,zmax
real :: PI=3.14159265,grassconstant=5.,litterconstant=5.
integer :: fueltotal=1,infuel=0,singlefuel=1,ntspecies=1
integer :: tfuelbins=1,ngrass=1,grassstep=1
real :: gmoistoverride=0,relhum=0.1
integer :: litout=0,controlseed=1

real :: grho=1.1766,gmoisture=0.06,gss=0.0005,gdepth=0.27

open(unit=1,file='duet.in',form='formatted',status='old')
read(1,*) nx
read(1,*) ny
read(1,*) nz
read(1,*) dx
read(1,*) dy
read(1,*) dz
read(1,*) seedchange
read(1,*) winddirection
read(1,*) windvary
close(1)

zmax=nz-1

open(unit=10,file='DUETInputs',form='formatted',status='unknown')
write(10,*) '&duetlist'
write(10,*) '                              '
write(10,*) '      inputprogram = ',inputprogram
write(10,*) '                              '
write(10,*) '!----------------------------!'
write(10,*) '! Grid variables'
write(10,*) '!----------------------------!'
write(10,*) '      nx = ',nx
write(10,*) '      ny = ',ny
write(10,*) '      nz = ',nz
write(10,*) '      dx = ',dx
write(10,*) '      dy = ',dy
write(10,*) '      dz = ',dz
write(10,*) '      zmax = ',zmax
write(10,*) '      PI = ',PI
write(10,*) '                               '
write(10,*) '!----------------------------!'
write(10,*) '! Fuel variables'
write(10,*) '!----------------------------!'
write(10,*) '      fueltotal = ',fueltotal
write(10,*) '      infuel = ',infuel
write(10,*) '      singlefuel = ',singlefuel
write(10,*) '      ntspecies = ',ntspecies
write(10,*) '      tfuelbins = ',tfuelbins
write(10,*) '      grassconstant = ',grassconstant
write(10,*) '      litterconstant = ',litterconstant
write(10,*) '      ngrass = ',ngrass
write(10,*) '      grassstep = ',grassstep
write(10,*) '      gmoistoverride = ',gmoistoverride
write(10,*) '                                '
write(10,*) '      iFIA = ',iFIA
write(10,*) '      iFIAspecies = ',iFIAspecies
write(10,*) '                                '
write(10,*) '!----------------------------!'
write(10,*) '! Wind and Time variables'
write(10,*) '!----------------------------!'
write(10,*) '      windprofile = ',windprofile
write(10,*) '      randomwinds = ',randomwinds
write(10,*) '      winddirection = ',winddirection
write(10,*) '      windvary = ',windvary
write(10,*) '                              '
write(10,*) '      YearsSinceBurn = ',YearsSinceBurn
write(10,*) '      StepsPerYear = ',StepsPerYear
write(10,*) '                              '
write(10,*) '      relhum = ',relhum
write(10,*) '                                '
write(10,*) '!----------------------------!'
write(10,*) '! Extra options'
write(10,*) '!----------------------------!'
write(10,*) '      litout = ',litout
write(10,*) '      controlseed = ',controlseed
write(10,*) '      seedchange = ',seedchange
write(10,*) '                              '
write(10,*) '/'
write(10,*) '                                '
write(10,*) '&grassdata'
write(10,*) '                                 '
write(10,*) '      grho = ',grho
write(10,*) '      gmoisture = ',gmoisture
write(10,*) '      gss = ',gss
write(10,*) '      gdepth = ',gdepth
write(10,*) '                                 '
write(10,*) '                                 '
write(10,*) '/'
close(10)




end subroutine FFmakeDuetInputs