      subroutine rthpest
      
!!     ~ ~ ~ PURPOSE ~ ~ ~
!!     this subroutine computes the hourly stream pesticide balance 
!!     (soluble and sorbed) 

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name          |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    ch_l2(:)      |km            |length of main channel
!!    ch_w(2,:)     |m             |average width of main channel
!!    chpst_conc(:) |mg/(m**3)     |initial pesticide concentration in reach
!!    chpst_koc(:)  |m**3/g        |pesticide partition coefficient between
!!                                 |water and sediment in reach
!!    chpst_mix(:)  |m/day         |mixing velocity (diffusion/dispersion) for
!!                                 |pesticide in reach
!!    chpst_rea(:)  |1/day         |pesticide reaction coefficient in reach
!!    chpst_rsp(:)  |m/day         |resuspension velocity in reach for pesticide
!!                                 |sorbed to sediment
!!    chpst_stl(:)  |m/day         |settling velocity in reach for pesticide
!!                                 |sorbed to sediment
!!    chpst_vol(:)  |m/day         |pesticide volatilization coefficient in 
!!                                 |reach
!!    drift(:)      |kg            |amount of pesticide drifting onto main
!!                                 |channel in subbasin
!!    hdepth(:)     |m             |depth of flow in hour
!!    hru_sub(:)    |none          |subbasin number where reach is located
!!    inum1         |none          |reach number
!!    inum2         |none          |inflow hydrograph storage location number
!!    rchwtr        |m^3 H2O       |water stored in reach at beginning of day
!!    rnum1         |none          |fraction of overland flow
!!    rtwtr         |m^3 H2O       |water leaving reach on day
!!    sedpst_act(:) |m             |depth of active sediment layer in reach for
!!                                 |pesticide
!!    sedpst_bry(:) |m/day         |pesticide burial velocity in river bed
!!                                 |sediment
!!    sedpst_conc(:)|mg/(m**3)     |inital pesticide concentration in river bed
!!                                 |sediment
!!    sedpst_rea(:) |1/day         |pesticide reaction coefficient in river bed
!!                                 |sediment
!!    varoute(11,:) |mg pst        |pesticide in solution
!!    varoute(12,:) |mg pst        |pesticide sorbed to sediment
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bury        |mg pst        |loss of pesticide from active sediment layer
!!                               |by burial
!!    difus       |mg pst        |diffusion of pesticide from sediment to reach
!!    reactb      |mg pst        |amount of pesticide in sediment that is lost
!!                               |through reactions
!!    reactw      |mg pst        |amount of pesticide in reach that is lost
!!                               |through reactions
!!    resuspst    |mg pst        |amount of pesticide moving from sediment to
!!                               |reach due to resuspension
!!    setlpst     |mg pst        |amount of pesticide moving from water to
!!                               |sediment due to settling
!!    hsolpst(:)  |mg pst/m^3    |soluble pesticide concentration in outflow
!!                               |on day
!!    hsorpst(:)  |mg pst/m^3    |sorbed pesticide concentration in outflow
!!                               |on day
!!    volatpst    |mg pst        |amount of pesticide in reach lost by
!!                               |volatilization
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    bedvol      |m^3           |volume of river bed sediment
!!    chpstmass   |mg pst        |mass of pesticide in reach
!!    depth       |m             |depth of water in reach
!!    fd2         |
!!    frsol       |none          |fraction of pesticide in reach that is soluble
!!    frsrb       |none          |fraction of pesticide in reach that is sorbed
!!    ii          |none          |counter
!!    jrch        |none          |reach number
!!    pstin       |mg pst        |total pesticide transported into reach
!!                               |during time step
!!    sedcon      |g/m^3         |sediment concentration
!!    sedpstmass  |mg pst        |mass of pesticide in bed sediment
!!    solpstin    |mg pst        |soluble pesticide entering reach during 
!!                               |time step
!!    sorpstin    |mg pst        |sorbed pesticide entering reach during
!!                               |time step
!!    thour       |hour          |flow duration
!!    wtrin       |m^3 H2O       |volume of water entering reach during time
!!                               |step
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: abs

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm

      implicit none             
      integer :: jrch, ii
      real*8 :: solpstin, sorpstin, pstin, depth, chpstmass, frsol, frsrb
      real*8 :: sedpstmass, bedvol, fd2, wtrin, solmax, sedcon, thour
      real*8 :: fr_stored, fr_routed, reactw1, reactw2

      jrch = 0
      jrch = inum1

!! calculate volume of active river bed sediment layer
      bedvol = 0.
      bedvol = ch_w(2,jrch) * ch_l2(jrch) * 1000. * sedpst_act(jrch)

      do ii = 1, nstep
!! initialize depth of water for pesticide calculations
      depth = 0.
      if (hdepth(ii) < 0.01) then  !!HR_2019: Wonder why 0.1 m is set to the lowest depth? Set to 0.01.
        depth = .01
      else
        depth = hdepth(ii)
      endif

!! calculate volume of water entering reach
      wtrin = 0.
      wtrin = hhvaroute(2,inum2,ii) * (1. - rnum1)
         
!! pesticide transported into reach during day
      solpstin = 0.
      sorpstin = 0.
      pstin = 0.
      solpstin = hhvaroute(11,inum2,ii) * (1. - rnum1)
      sorpstin = hhvaroute(12,inum2,ii) * (1. - rnum1)
      pstin = solpstin + sorpstin
!! get fraction of pesticde routed and stored based on incoming, outgoing, and stored water 
      !rtwtr         |m^3 H2O       |water leaving reach on day
      !rchwtr        |m^3 H2O       |water stored in reach at beginning of day
      !!MW 2019: I think fraction stored should be average of initial storage and ending storage for time step, divided by initial storage + inflow ... See below
      fr_routed = (hrtwtr(ii)) / (hrchwtr(ii) + wtrin) 
      fr_stored = 1. - fr_routed      
      !write(*,*) fr_routed

!! add pesticide drifting from HRUs in subbasin to reach
!      if (rtwtr > 0.) then
!        pstin = pstin + (drift(jrch) * 1.e6)
!      else
!        sedpst_conc(jrch) = sedpst_conc(jrch) + drift(jrch) * 1.e6 /    &
!     &                                                            bedvol
!      endif
 
      !! calculate mass of pesticide in reach
      chpstmass = 0.
      chpstmass = pstin + chpst_conc(jrch) * hrchwtr(ii)
      
      !! calculate mass of pesticide in bed sediment
      sedpstmass = 0.
      sedpstmass = sedpst_conc(jrch) * bedvol

      if (chpstmass + sedpstmass < 1.e-9) then
        chpst_conc(jrch) = 0.
        sedpst_conc(jrch) = 0.
      end if
      if (chpstmass + sedpstmass < 1.e-9) return

!!in-stream processes
      if (hrtwtr(ii) / (idt*60.) > 1.e-9) then !HR Aug 2019 - Change threshold from 0.01 to avoid all pesticide going to benthic if low flow
        !! calculated sediment concentration
        sedcon = 0.
        !sedcon = hsedyld(ii) / hrtwtr(ii) * 1.e6
        !!HR Sep 2019; The sed concentration is based on current initial sed (sedinorg). Resuspended (rchdy(56,jrch) from previous day is added in routing routine
        sedcon = (hsedinorg(ii)) / (hrchwtr(ii) + wtrin) * 1.e6  

        !! calculate fraction of soluble and sorbed pesticide
        frsol = 0.
        frsrb = 0.
        !if (solpstin + sorpstin > 1.e-6) then
        !! HR_MW_19: This calculation of Frsol/Frsrb needs to happen every timestep based on the amount of sediment 
        if (chpstmass > 1.e-9) then                           
          if (chpst_koc(jrch) > 0.) then
            frsol = 1. / (1. + chpst_koc(jrch)* sedcon) !HR updated / fixed equation for sedcon
          else
            frsol = solpstin / (solpstin + sorpstin)
          end if
          frsrb = 1. - frsol
        else
          !!drifting pesticide is only pesticide entering
          !!and none is sorbed
          frsol = 1.
          frsrb = 0.
        end if

        !! ASSUME POR=0.5; DENSITY=2.6E6; KD2=KD1
        !fd2 = 1. / (.5 + chpst_koc(jrch))
        !!HR_MW_2019 corrected equation based on manual
        fd2 = 1. / (.5 + ((1. - 0.5)* 2600000. * chpst_koc(jrch)))                                                      

        !! calculate flow duration
         thour = 0.
         thour = hhtime(ii)
         !write(*,*) thour                 
         if (thour > 1.0) thour = 1.0
         !!thour = 1.0 !!HR_2019 use actual routing time for routed elements

        !! calculate amount of pesticide that undergoes chemical or
        !! biological degradation on day in reach
        !! HR 16 Aug 2019 according to MFW, 3/12/12: modify decay to be 1st order
        !! HR Aug 2019 split storage (pestprev) and routed (pstin) time
        !! reactw = chpst_rea(jrch) * chpstmass * thour / 24.                       
        !!reactw = chpstmass - (chpstmass * EXP(-1. * chpst_rea(jrch) * (thour/24.)))      
        reactw1 = fr_stored*chpstmass  - (fr_stored*chpstmass * EXP(-1. * chpst_rea(jrch)     !!MW_2019: NOt sure why we aren't using fr_stored*chpstmass here
     &           * (1./24.)))
        reactw2 = fr_routed*chpstmass - (fr_routed*chpstmass * EXP(-1. * chpst_rea(jrch)          !!MW_2019: NOt sure why we aren't using fr_routed*chpstmass here
     &           * (thour/24.)))
        reactw = reactw1 + reactw2
        chpstmass = chpstmass - reactw

        !! calculate amount of pesticide that volatilizes from reach
		!!HR change to match daily routing (add frsol)
        ! HR August 2019; use full hour for fraction in storage and actual travel time for routed elements
        !volatpst = chpst_vol(jrch) * frsol * chpstmass * thour / (depth * 24.)
        volatpst = chpst_vol(jrch) * frsol * chpstmass  / depth * (((thour /24.) * fr_routed) + ((1./24)* fr_stored)) 
        if (volatpst > frsol * chpstmass) then
          volatpst = frsol * chpstmass !HR add frsol
          chpstmass = chpstmass - volatpst !hr update equation
        else
          chpstmass = chpstmass - volatpst
        end if

        !! calculate amount of pesticide removed from reach by
        !! settling
        ! HR Aug 2019; full hour for in storage actual time routed                                                         
        !setlpst = chpst_stl(jrch) * frsrb * chpstmass * thour / (depth * 24.)
        setlpst = chpst_stl(jrch) * frsrb * chpstmass / depth * (((thour /24.) * fr_routed) + ((1. /24.)* fr_stored))   !!MW_2019: This makes sense to me             
        if (setlpst > frsrb * chpstmass) then !HR modify to match daily routing; add frsb
          setlpst = frsrb * chpstmass
          chpstmass = chpstmass - setlpst
        else
          chpstmass = chpstmass - setlpst
        end if
        sedpstmass = sedpstmass + setlpst

        !! calculate resuspension of pesticide in reach
        ! HR Aug 2019; full timestep 
        resuspst = chpst_rsp(jrch) * sedpstmass / depth * ((1./24.))   !!MW_2019: Since we are not routing sediment in benthic I think this process is 100% for the full 1-day timestep
        if (resuspst > sedpstmass) then
          resuspst = sedpstmass
          sedpstmass = 0.
        else
          sedpstmass = sedpstmass - resuspst
        end if
        chpstmass = chpstmass + resuspst

        !! calculate diffusion of pesticide between reach and sediment
        ! HR Aug 2019; full hour for in storage, actual time for routed                                                              
        !difus = chpst_mix(jrch) * (fd2 * sedpstmass - frsol * chpstmass) * thour / (depth * 24.)
        difus = chpst_mix(jrch) * (fd2 * sedpstmass - frsol * chpstmass) / depth * (((thour /24.) * fr_routed) + ((1. /24.)* 
     &      fr_stored))   !!MW_2019: This makes sense to me        
        if (difus > 0.) then
          if (difus > sedpstmass) then
            difus = sedpstmass
            sedpstmass = 0.
          else
            sedpstmass = sedpstmass - abs(difus)
          end if
          chpstmass = chpstmass + abs(difus)
        else
          if (abs(difus) > chpstmass) then
            difus = -chpstmass
            chpstmass = 0.
          else
            chpstmass = chpstmass - abs(difus)
          end if
          sedpstmass = sedpstmass + abs(difus)
        end if

        !! calculate removal of pesticide from active sediment layer
        !! by burial
        bury = sedpst_bry(jrch) * sedpstmass / (sedpst_act(jrch) * 24.)
        if (bury > sedpstmass) then
          bury = sedpstmass
          sedpstmass = 0.
        else
          sedpstmass = sedpstmass - bury
        end if

        !! verify that water concentration is at or below solubility
        solmax = 0.
        solmax = pest_sol * (rchwtr + wtrin)
        if (solmax < chpstmass * frsol) then
         sedpstmass = sedpstmass + (chpstmass * frsol - solmax)
         chpstmass = chpstmass - (chpstmass * frsol - solmax)
        end if
        
      else   
!!insignificant flow
        sedpstmass = sedpstmass + chpstmass
        chpstmass = 0.
      end if

!! sediment processes
      !! calculate loss of pesticide from bed sediments by reaction
      !reactb = sedpst_rea(jrch) * sedpstmass / 24.
      !!HR_2019: Original equation replaced it with 1st order decay      
      reactb = sedpstmass - (sedpstmass * EXP(-1. * sedpst_rea(jrch) * (1./24.))) !! assume full one hour for bentich sediment
      if (reactb > sedpstmass) then
        reactb = sedpstmass
        sedpstmass = 0.
      else
        sedpstmass = sedpstmass - reactb
      end if

!! calculate pesticide concentrations at end of hour
      chpst_conc(jrch) = 0.
      sedpst_conc(jrch) = 0.
      if (hrchwtr(ii) + wtrin > 1.e-9) then 
        chpst_conc(jrch) = chpstmass / (hrchwtr(ii) + wtrin)
      else
        sedpstmass = sedpstmass + chpstmass !!MW_2019: Don't we also need to set chpst_conc(jrch) = 0 in this else?
        chpst_conc(jrch) = 0.                     
      end if
      sedpst_conc(jrch) = sedpstmass / bedvol


!! calculate amount of pesticide transported out of reach
      if (hrtwtr(ii) > 1.e-9) then !HR Aug 2019 - Change threshold from 0.01 to avoid all pesticide going to benthic if low flow
        hsolpst(ii) = chpst_conc(jrch) * frsol
        hsorpst(ii) = chpst_conc(jrch) * frsrb
      else
        hsolpst(ii) = 0.
        hsorpst(ii) = 0.
      end if
      end do

      return
      end