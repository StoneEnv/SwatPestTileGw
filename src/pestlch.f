      subroutine pestlch
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine calculates pesticides leached through each layer,
!!    pesticide transported with lateral subsurface flow, tile flow 
!!    and pesticide transported with surface runoff

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    flat(:,:)    |mm H2O        |lateral flow in soil layer on current day

!!    hrupest(:)   |none          |pesticide use flag:
!!                                | 0: no pesticides used in HRU
!!                                | 1: pesticides used in HRU
!!    ihru         |none          |HRU number
!!    npmx         |none          |number of different pesticides used in
!!                                |the simulation
!!    npno(:)      |none          |array of unique pesticides used in watershed
!!    percop       |none          |pesticide percolation coefficient (0-1)
!!                                |0: concentration of pesticide in surface
!!                                |   runoff is zero
!!                                |1: percolate has same concentration of
!!                                |   pesticide as surface runoff
!!    pst_wsol(:)  |mg/L (ppm)    |solubility of chemical in water
!!    sol_bd(:,:)  |Mg/m**3       |bulk density of the soil
!!    sol_kp(:,:,:)|(mg/kg)/(mg/L)|pesticide sorption coefficient, Kp; the
!!                 |  or m^3/ton  |ratio of the concentration in the solid
!!                                |phase to the concentration in solution
!!    sol_ul(:,:)  |mm H2O        |amount of water held in the soil layer at
!!                                |saturation (sat - wp water)
!!    sol_nly(:)   |none          |number of layers in soil profile
!!    sol_por(:,:) |none          |total porosity of soil layer expressed as
!!                                |a fraction of the total volume
!!    sol_prk(:,:) |mm H2O        |percolation from soil layer on current day
!!    sol_pst(:,:,:)|kg/ha        |amount of pesticide in layer
!!    sol_wpmm(:,:)|mm H20        |water content of soil at -1.5 MPa (wilting
!!                                |point)
!!    sol_z(:,:)   |mm            |depth to bottom of soil layer
!!    surfq(:)     |mm H2O        |surface runoff generated on day in HRU
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name         |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    lat_pst(:)   |kg pst/ha     |amount of pesticide in lateral flow in HRU
!!                                |for the day
!!    pst_surq(:,:)|kg/ha         |amount of pesticide type lost in surface
!!                                |runoff on current day in HRU
!!    pst_lat(:)   |kg pst/ha     |amount of pesticide in lateral flow in HRU
!!                                |for the day
!!    pst_sol(:,:) |kg/ha         |amount of pesticide type leached from soil
!!                                |profile on current day in HRU
!!    pstsol(:)    |kg/ha         |amount of pesticide type leached from soil
!!                                |profile on current day
!!    tile_pst(:)  |kg pst/ha     |amount of pesticide in tile flow in HRU for
!!                                |the day
!!    zdb(:,:)     |mm            |division term from net pesticide equation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    co          |kg/mm-ha      |concentration of pesticide in water
!!    cocalc      |kg/mm-ha      |calc concentration of pesticide in water
!!    csurf       |kg/mm-ha      |concentration of pesticide in surq and latq
!!    dg          |mm            |depth of soil layer
!!    j           |none          |HRU number
!!    k           |none          |counter
!!    kk          |none          |pesticide number from pest.dat
!!    ly          |none          |counter (soil layers)
!!    qsurf       |mm H2O        |surface runoff for layer
!!    vf          |
!!    xx          |kg/ha         |amount of pesticide removed from soil layer
!!    yy          |
!!    zdb1        |mm            |division term from net pesticide equation
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: Exp, Min

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~

      use parm
      implicit none

      integer :: j, ly, k, kk
      integer :: iprintpsthrulayer
      real*8 :: dg, yy, qsurf, vf, zdb1, xx, co, csurf, cocalc
      real*8 :: mobPest_frac
      real*8 :: tile_pst_conc
 
      j = 0
      j = ihru

      if (hrupest(j) /= 0) then

        do ly = 1, sol_nly(j)
          if (ly == 1) then
            yy = 0.
          else
            yy = 0.
            yy = sol_z(ly-1,j)
          end if
          dg = 0.
          dg = sol_z(ly,j) - yy

          do k = 1, npmx
            kk = 0
            kk = npno(k)

            if (kk > 0) then
              qsurf = 0.
              if (ly == 1) then
                qsurf = surfq(j)
              else
                qsurf = 0.
              endif

              zdb1 = 0.
              !zdb1 = sol_ul(ly,j) + sol_kp(k,j,ly) * sol_bd(1,j) * dg !HR always uses bulk density of the first layer
              zdb1 = sol_ul(ly,j) + sol_kp(k,j,ly) * sol_bd(ly,j) * dg

              !! units: mm + (m^3/ton)*(ton/m^3)*mm = mm
              if (ly == 1) zdb(k,j) = zdb1

              vf = 0.
              !vf = qsurf + sol_prk(ly,j) + flat(ly,j)
              !hendrik 1/2017 !If loop to distinguish between surface, lateral, and tile flow
              !!MW: This calculated the mobile water in the layer (vf)
              if (ly == 1) then
                 vf = qsurf + sol_prk(ly,j) + flat(ly,j)
              else if (ldrain(j) == ly) then
                vf = qtile + sol_prk(ly,j) + flat(ly,j)
              else
                vf = sol_prk(ly,j) + flat(ly,j)
              end if
              !hendrik 1/2017


              if (sol_pst(k,j,ly) >= 1.e-9 .and. vf > 0.) then
                xx = 0.
                mobPest_frac = 1. - Exp(-vf / (zdb1 + 1.e-6))
                mobPest_frac = mobPest_frac * 1.
                !mobPest_frac = Min(1.,mobPest_frac )
                xx = sol_pst(k,j,ly) * mobPest_frac   !Amount of pesticide leaving the soil layer kg/ha
                !if (kk==236) then
                !  write (*,*) (1. - Exp(-vf / (zdb1 + 1.e-6)))
                !end if
                cocalc = 0.
                co = 0.
                if (ly == 1) then
                  cocalc = xx /                                         
     &           (sol_prk(ly,j) + percop * (qsurf + flat(ly,j)) + 1.e-6) !Concentration of pesticde leaving soi layer
                !hendrik 1/2017 !Add tile-drain case
                !!MW: This qtile could be from more than only the layer with he tile; thus concentration underdone
                else if (ldrain(j) == ly) then
                 cocalc = xx / 
     &           (qtile + sol_prk(ly,j) + flat(ly,j) + 1.e-6)
                !hendrik 1/2017
                else
                  cocalc = xx / (sol_prk(ly,j) + flat(ly,j) + 1.e-6)
                end if
                co = Min(pst_wsol(kk) / 100., cocalc)
               
                !! calculate concentration of pesticide in surface runoff and lateral flow
                csurf = 0.
                if (ly == 1) then
                  csurf = percop * co             !Concentration of pesticde in surface, lateral, and tile runoff mgl/l
                else
                  csurf = co
                end if

                !! calculate pesticide leaching
                xx = 0.
                xx = co * sol_prk(ly,j)
                if (xx > sol_pst(k,j,ly)) xx = sol_pst(k,j,ly)
                sol_pst(k,j,ly) = sol_pst(k,j,ly) - xx

                if (ly < sol_nly(j)) then
                  sol_pst(k,j,ly+1) = sol_pst(k,j,ly+1) + xx
                else
                  pstsol(k) = xx
                end if

                !! calculate pesticide lost in surface runoff
                if (ly == 1) then
                  yy = 0.
                  yy = csurf * surfq(j)
                  if (yy > sol_pst(k,j,ly)) yy = sol_pst(k,j,ly)
                  sol_pst(k,j,ly) = sol_pst(k,j,ly) - yy
                  pst_surq(k,j) = yy 
                endif

                !HR 
                !! calculate pestcide in tile flow 
                if (ldrain(j) == ly) then				
					yy = 0.
					yy = csurf * qtile
					if (yy > sol_pst(k,j,ly)) yy = sol_pst(k,j,ly)
					sol_pst(k,j,ly) = sol_pst(k,j,ly) - yy
					tile_pst(k) = yy
                end if
                !end HR

                !! calculate pesticide lost in lateral flow
                yy = 0.
                yy = csurf * flat(ly,j)
                if (yy > sol_pst(k,j,ly)) yy = sol_pst(k,j,ly)
                sol_pst(k,j,ly) = sol_pst(k,j,ly) - yy
                lat_pst(k) = lat_pst(k) + yy 

              end if
            end if
            !hendrik 09/2017 custom subsurface pest output
            iprintpsthrulayer = 0
            if (curyr > 0 .and. iprintpsthrulayer==1) then
            !if (curyr > nyskip .and. iprintpsthrulayer==1) then
              !if (j==2 .and. kk==236) then
              !write(*,*) subnum(j)
              !write(*,*) hruno(j)
              if (subnum(j)=='00004' .and. hruno(j)=='0010' .and. kk==236) then
               if (qtile > 1.e-06) then
                 tile_pst_conc = tile_pst(k)*100000/qtile
               else
                 tile_pst_conc = -1
               endif
               write(32137,"(2i11, 2a11, 3i11, 19e11.3)") iyr, iida,      
     &             subnum(j), hruno(j), kk, ly, ldrain(j),              
     &             sol_pst(k,j,ly),lat_pst(k), tile_pst(k),         
     &             pstsol(k),                                        
     &             qtile, flat(ly,j), sol_prk(ly,j), gw_q(j),sw_excess,
     &             revapday,gwseep,
     &             (1. - Exp(-vf / (zdb1 + 1.e-6))),
     &             vf, sol_ul(ly,j), sol_kp(k,j,ly), sol_bd(ly,j), dg,
     &             mobPest_frac,
     &             tile_pst_conc
               if (ly==sol_nly(j)) then !print groundwater
                write(32137,"(2i11, 2a11, 3i11, 19e11.3)") iyr, iida,      
     &             subnum(j), hruno(j), kk, ly+1, 0,              
     &             pst_shallst(k,j),0, 0,         
     &             pstsol(k),                                        
     &             0, 0, sol_prk(ly,j), gw_q(j),sw_excess,
     &             revapday,gwseep,
     &             0,
     &             0, 0, 0, 0, 0,
     &             mobPest_frac,
     &             -1.
               endif          
             end if
            end if
          end do
        end do
      end if

      return
      end