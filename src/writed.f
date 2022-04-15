      subroutine writed
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine contains the daily output writes

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    da_ha       |ha            |area of watershed in hectares
!!    hrupest(:)  |none          |pesticide use flag:
!!                               | 0: no pesticides used in HRU
!!                               | 1: pesticides used in HRU
!!    hrupstd(:,1,:)|mg pst      |amount of pesticide type in surface runoff
!!                               |contribution to stream from HRU on day
!!                               |(in solution)
!!    hrupstd(:,2,:)|mg pst      |amount of pesticide type in surface runoff
!!                               |contribution to stream from HRU on day
!!                               |(sorbed to sediment)
!!    iida        |julian date   |current day of simulation 
!!    iprint      |none          |print code:
!!                               |0 monthly
!!                               |1 daily
!!                               |2 annually
!!    iprp        |none          |print code for output.pst file
!!                               |0 do not print pesticide output
!!                               |1 print pesticide output
!!    isproj      |none          |special project code:
!!                               |1 test rewind (run simulation twice)
!!    iyr         |year          |year being simulated (eg 1980)
!!    mstdo       |none          |watershed output array size
!!    nhru        |none          |number of HRUs in watershed
!!    npmx        |none          |number of different pesticides used in
!!                               |the simulation
!!    subtot      |none          |number of subbasins in watershed
!!    wshddayo(1) |mm H2O        |average amountof precipitation in watershed
!!                               |for the day
!!    wshddayo(3) |mm H2O        |surface runoff in watershed for day
!!    wshddayo(4) |mm H2O        |lateral flow contribution to streamflow in
!!                               |watershed for day
!!    wshddayo(5) |mm H2O        |water percolation past bottom of soil profile
!!                               |in watershed for day
!!    wshddayo(6) |mm H2O        |water yield to streamflow from HRUs in
!!                               |watershed for day
!!    wshddayo(7) |mm H2O        |actual evapotranspiration in watershed
!!                               |for day
!!    wshddayo(12)|metric tons   |sediment yield from HRUs in watershed 
!!                               |for day
!!    wshddayo(35)|mm H2O        |amount of water stored in soil profile in
!!                               |watershed for day
!!    wshddayo(40)|kg N/ha       |organic N loading to stream in watershed for
!!                               |day
!!    wshddayo(41)|kg P/ha       |organic P loading to stream in watershed for
!!                               |day
!!    wshddayo(42)|kg N/ha       |nitrate loading to stream in surface runoff
!!                               |in watershed for day
!!    wshddayo(43)|kg P/ha       |soluble P loading to stream in watershed for
!!                               |day
!!    wshddayo(44)|kg N/ha       |plant uptake of N in watershed for day
!!    wshddayo(45)|kg N/ha       |nitrate loading to stream in lateral flow
!!                               |in watershed for day
!!    wshddayo(46)|kg N/ha       |nitrate percolation past bottom of soil
!!                               |profile in watershed for day
!!    wshddayo(104)|mm H2O        |groundwater contribution to stream in
!!                               |watershed on day
!!    wshddayo(108)|mm H2O        |potential evapotranspiration in watershed
!!                               |on day
!!    wshddayo(109)|mm H2O        |drainage tile flow contribution to stream
!!                               |in watershed on day
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    hrupstm(:,1,:)|mg pst      |amount of pesticide type in surface runoff
!!                               |contribution to stream from HRU during month
!!                               |(in solution)
!!    hrupstm(:,2,:)|mg pst      |amount of pesticide type in surface runoff
!!                               |contribution to stream from HRU during month
!!                               |(sorbed to sediment)
!!    hrupstm(:,3,:)|mg pst/ha   |total pesticide loading to stream in surface
!!                               |runoff from HRU during month
!!    wshddayo(12)|metric tons/ha|sediment yield from HRUs in watershed 
!!                               |for day
!!    wshdmono(:) |varies        |watershed monthly output array
!!                               |(see definitions for wshddayo array elements)
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name        |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j           |none          |counter
!!    k           |none          |counter
!!    pstsum      |mg pst        |pesticide loading in watershed
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    SWAT: rchday

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
      use parm

      integer :: j, k
      real*8 :: pstsum


!!    write statement to new output file (output.swr)
!!    writes out the amount of water stored in the soil layer
      if (isto > 0) then 
        do j = 1, nhru
          write (129,5000) iida, j, subnum(j), hruno(j), (sol_st(j1,j), 
     &         j1 = 1, sol_nly(j))
!          write (129,5000) iida, subnum(j), hruno(j),                   
!     &             (sol_no3(j1,j), j1 = 1, sol_nly(j))
        enddo
      end if

      if (iprint == 1.or.iprint==3) then
        if (da_ha < 1.e-9) then
	    call rchday
	    call rseday
	    return
	  end if

        !! daily write to output.std
        if (iscen == 1) then
        write (26,6200) iida, wshddayo(1), wshddayo(3), wshddayo(4),    
     &                 wshddayo(104), wshddayo(5), wshddayo(109),       
     &                 wshddayo(35), wshddayo(7), wshddayo(108),        
     &                 wshddayo(6), wshddayo(12) / da_ha, wshddayo(42), 
     &                 wshddayo(45), wshddayo(46), wshddayo(44),        
     &                 wshddayo(40), wshddayo(43), wshddayo(41),        
     &                 wshddayo(111)
        else if (isproj == 1) then
        write (19,6200) iida, wshddayo(1), wshddayo(3), wshddayo(4),    
     &                 wshddayo(104), wshddayo(5), wshddayo(109),       
     &                 wshddayo(35), wshddayo(7), wshddayo(108),        
     &                 wshddayo(6), wshddayo(12) / da_ha, wshddayo(42), 
     &                 wshddayo(45), wshddayo(46), wshddayo(44),        
     &                 wshddayo(40), wshddayo(43), wshddayo(41)
        endif

        !! daily write to pesticide output file (output.pst) for HRUs
        do j = 1, nhru
          if (hrupest(j) == 1) then
          pstsum = 0.
          do k = 1, npmx
              pstsum = pstsum + hrupstd(k,1,j) + hrupstd(k,2,j) +                       
     &                          hrupstd(k,4,j) + hrupstd(k,5,j) +
     &                          hrupstd(k,6,j) + hrupstd(k,7,j)
          end do
          if (pstsum > 0.) then
            if (iprp /= 0) then
                write (30,5100) subnum(j), hruno(j), iyr, iida,         
     &                     (hrupstd(k,1,j), hrupstd(k,2,j), 
     &                      hrupstd(k,4,j), hrupstd(k,5,j),
     &                      hrupstd(k,6,j), hrupstd(k,7,j), k = 1, npmx) !HR added lateral, tile, and gw
            end if
            if (iprp == 2) then !HR write concentrations in ug/L (convert from mg/mm)
                do k=1, npmx
                if (surfq(j)>0) then
                    conc1(k) = hrupstd(k,1,j)/(surfq(j)*hru_ha(j)*10)
                    conc2(k) = hrupstd(k,2,j)/(surfq(j)*hru_ha(j)*10)
                else
                    conc1(k) = 0
                    conc2(k) = 0
                end if
                if (latq(j)>0) then
                    conc3(k) = hrupstd(k,4,j)/(latq(j)*hru_ha(j)*10)
                else
                    conc3(k) = 0
                end if
                if (tileq(j)>0) then
                    conc4(k) = hrupstd(k,5,j)/(tileq(j)*hru_ha(j)*10)
                else
                    conc4(k) = 0
                end if
                if (gw_q(j)>0) then
                    conc5(k) = hrupstd(k,6,j)/(gw_q(j)*hru_ha(j)*10)
                else
                    conc5(k) = 0
                end if
                if (gw_qdeep(j)>0) then
                    conc6(k) = hrupstd(k,7,j)/(gw_qdeep(j)*hru_ha(j)*10)
                else
                    conc6(k) = 0
                end if                
                end do
                write (32168,5100) subnum(j), hruno(j), iyr, iida,         
     &                     (conc1(k), conc2(k), conc3(k), conc4(k), conc5(k), conc6(k), k = 1, npmx)
            end if
          end if
          end if
        end do

        !! write daily reach output
        call rchday

        !! write daily sediment routing output (.sed)
        call rseday

      end if

      !! write velocities for steve/woody in temp file (Balaji)
      if (itemp == 1 .and. nrch > 0) then 
         write (141,5001) iida,iyr,(vel_chan(k),k= 1,nrch)
         write (142,5001) iida,iyr,(dep_chan(k),k= 1,nrch)
      end if 

!! monthly watershed output
      wshddayo(12) = wshddayo(12) / (da_ha + 1.e-6)

      wshdmono = wshdmono + wshddayo
      wpstmono = wpstmono + wpstdayo
      hrupstm = hrupstm + hrupstd


      return
!5000  format(i5,1x,a5,a4,1x,500e12.4)

!5000  format (i5,1x,i5,1x,500e12.4)
5000  format (i5,1x,i5,1x,a5,a4,1x,500e12.4)
5001  format(2i5,500f12.4)
5100  format(1x,a5,a4,1x,i4,1x,i3,1x,250(e17.4,1x))
5200  format(i7,i9,i6,i5,1x,e9.4,f12.3,f7.1,f14.3)
!!6200  format(i5,13f7.2,2f5.2,1x,5f8.2)
6200  format(i5,15f8.2,1x,4f8.2)
      end