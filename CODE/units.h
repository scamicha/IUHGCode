!    The following units: Mstar, Rdiskau, Tstar are model
!    dependent.  Must be set by the user.   
!    The rest are astronomical units.

      REAL*8 Mstar,Rstar,Rdiskau,Tstar,gridlim,phylim
      REAL*8 Msuncgs,Rsuncgs,AUcgs,sigmacgs,Gcgs,bkmpcgs
      PARAMETER (Mstar=1.0,Rstar=2.0,Rdiskau=40.0,Tstar=4000.0)
      PARAMETER (Msuncgs=1.989e33, Rsuncgs=6.96e10, AUcgs=1.496e13)
      PARAMETER (sigmacgs=5.670e-5, Gcgs=6.672e-8, phylim=1.e-8) 
      PARAMETER (bkmpcgs=8.254d7,gridlim=1.e-12)

!    Most of the following variables are used in the source routine only. 
!    torp is one ORP in code time units, defined as 2*pi/omega(200,2,1) 
!    of the starting model.
!    cq = 3 for AV on, 0 for AV off. 
!    irtype = 0 for no irr, 1 for irradiation.
!    ictype is the cooling type: 1 for const. Tcool, 2 for 
!    Edd. grey atm. only, 3 for diff. aprox., 4 for diff. 
!    aprox. with shining atm.
!    cct is only used when using type 1.  It is the cooling
!    time in ORPS.      
!    tcoollmt is the lower limit of cooling time in ORP.
!    tirrlmt  is the lower limit of irradiation time in ORP.
!    Tbgrnd is the lower limit of the temperature
!    irop is used to spread the irradiation to a few more cells by dividing the
!    opacities by some number. 
!    jirr is the first first j zone to be irradiated.
!    jcool is the j zone at which cooling starts.  Disabled in most ictypes.
!    fk = metallicity factor: fraction of standard metallicity to keep.
!    tenvk = constant envelop temperature
!    SETUNITS = 0 for polytropic units (G=K=M=1) and 1 for Aaron's units
!       (G=1, 1 R = 1 AU, 1 M = Msun)
!    Use H2STAT to select what type of mixture you want.
!    0 = pure para
!    1 = pure ortho -- this is not really physical for astro applications
!    2 = set mixture
!    3 = equilibrium
!    4 = single gamma law
!    Be sure to set ac and bc for mixture (ac: para component, bc: ortho component)
!    also pick a metallicity by adjusting xabun, yabun, and zabun.
!  
!    star_state can take 1 of 4 values they are:
!    fixed -- The traditional method of fixing the star to the center of 
!             the grid.
!    indirect -- Transform to a non-inertial frame with the star at the 
!                grid center and add the potential due to an accelerated
!                frame. If a lot of mass is leaving the grid you shouldn't
!                use this since it won't account for it
!    wiggle -- Has the same effect as star_wiggle=.true. and restart_wiggle
!              =.false.
!    wiggle_restart -- Same effect as star_wiggle=.true. and restart_wiggle
!                      =.true.
!
!    ISTOR2 sets how often saved and coolheat_full files are output.
!    make sure that ISTOR2 is a multiple of IDIAG in fort.5 or the output
!    won't work correctly.	
!    outpath is the path to the ouput files, use ./ if you want files in
!    the same directory as the binary you are executing.  

      REAL*8 torp,cct,Tbgrnd,tcoollmt,theatlmt,cq,irop,fk,tenvk,amp0
      real*8, parameter :: ac = 1.d0, bc = 3.d0
      real*8, parameter :: xabun = .73d0,yabun=.25d0,zabun=.02d0
      INTEGER ictype,irtype,jcool,jirr,SETUNITS,H2STAT
      PARAMETER(torp=1605.63,cct=2.,cq=3.0,tcoollmt=50.,theatlmt=50.)
      PARAMETER(ictype=69,irtype=0,jcool=12,jirr=28,irop=1.0,Tbgrnd=3.)
      PARAMETER(tenvk=3.d0,fk=1.d0)
      PARAMETER(SETUNITS=0,H2STAT=2) 
      PARAMETER(amp0=0.0001) ! amplitude of initial, random perturbation

      logical, parameter :: restart_fluid  = .false. 
      logical, parameter :: fluid_elements = .false.
      logical, parameter :: planet         = .true.
      REAL*8, parameter  :: mass_ratio     = 0.02d0
      character(len=14)  :: star_state     = "indirect"
      integer,parameter  :: istor2         = 1000
      character(len=80)  :: outpath        = "./DC_DATA/"
