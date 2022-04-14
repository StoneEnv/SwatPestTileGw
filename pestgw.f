      subroutine pestgw
      !hendrik 5/2017
      
!!    ~ ~ ~ PURPOSE ~ ~ ~
!!    this subroutine estimates groundwater pesticide contribution to
!!    streamflow

!!    ~ ~ ~ INCOMING VARIABLES ~ ~ ~
!!    name            |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    deepst(:)       |mm H2O        |depth of water in deep aquifer
!!    ihru            |none          |HRU number
!!    gw_delaye(:)    |none          |Exp(-1./(delay(:)) where delay(:) is the 
!!                                   |groundwater delay (time required for water
!!                                   |leaving the bottom of the root zone to reach
!!                                   |the shallow aquifer; units-days)
!!    gwseep          |mm H2O        |amount of water recharging deep aquifer on
!!                                   |current day in HRU
!!    gw_q(:)         |mm H2O        |groundwater contribution to streamflow from
!!                                   |HRU on current day
!!    gw_qdeep(:)     |mm H2O        |Deep aquifer groundwater contribution to streamflow 
!!                                   |from HRU on current day
!!    hrupest(:)      |none          |Flag whether HRU has pesticide.
!!    npno(:)         |none          |Array of unique pesticides used in watershed
!!    npmx            |none          |Total number of peticides in simulation.
!!    pst_rchrg(:,:)  |mg/ha         |Amount of pesticide entering the shallow aquifer
!!    pstsol(:)       |mg/ha         |Amount of pesticide type leached from soil
!!                                   |profile on current day
!!    revapday        |mm H2O        |Amount of water moving from the shallow 
!!                                   |aquifer into the soil profile or being taken
!!                                   |up by plant roots in the shallow aquifer
!!    rootzpest(:,:)  |kg/ha         |Amount of pesticide in vadose zone on previos day. 
!!    shallst(:)      |mm H2O        |depth of water in shallow aquifer

!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ OUTGOING VARIABLES ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    pst_deepst(:,:)  |kg/ha         |Amount of peticide stored in deep aquifer
!!    pst_gw(:,:)      |kg/ha         |Amount of pesticide entering the channel via  
!!                                    |shallow aquifer groundwater flow
!!    pst_gwdeep(:,:)  |kg/ha         |Amount of pesticide entering the channel via deep 
!!                                    |aquifer groundwater flow
!!    pst_shallst(:,:) |kg/ha         |Amount of peticide stored in shallow aquifer
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ LOCAL DEFINITIONS ~ ~ ~
!!    name             |units         |definition
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!    j                |none          |HRU number
!!    k                |none          |Pesticide number 
!!    kk               |none          |Pesticide number in pesticide database
!!    rchrgpest1       |kg/ha         |amount of pesticide entering shallow aquifer on
!!                                    |previous day
!!    revappst         |kg/ha         |Amount of pesticide moving from the shallow
!!                                    |aquifer int the soi profile or being taken 
!!                                    |up by plant roots in the shallow aquifer
!!                                    |currently not active / set to 0.
!!    gwseeppst        |kg/ha         |Amount of pesticide recharging deep aquifer
!!    rootzpest1       |kg/ha         |Amount of pesticide in vadose zone on previos day. 
!!    xx               |mm h20        |Water volume in shallow aquifer diluting with 
!!                                    |pesticide
!!    ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!!    ~ ~ ~ SUBROUTINES/FUNCTIONS CALLED ~ ~ ~
!!    Intrinsic: amax1

!!    ~ ~ ~ ~ ~ ~ END SPECIFICATIONS ~ ~ ~ ~ ~ ~
!!    revap is subtracted and rchrg is delayed (johnson, 1977)

      use parm
      implicit none

      
      integer :: j, k, kk
      real*8  :: xx
      real*8  :: rchrgpest1, revappst, gwseeppst, rootzpest1

      j = ihru

      if (hrupest(j) /= 0) then
        do k = 1, npmx
            kk = npno(k)

            rchrgpest1 = pst_rchrg(k,j)
            rootzpest1 = rootzpest(k,j)
            if (rchrgpest1< 1.e-6) rchrgpest1 = 0.0

            !! compute pesticide aquifer loading from recharge for current day
            pst_rchrg(k,j) = (1.- gw_delaye(j)) * pstsol(k) 
     &                    + gw_delaye(j) *  rchrgpest1
            pst_shallst(k,j) = pst_shallst(k,j) + pst_rchrg(k,j)
            rootzpest(k,j) = rootzpest1 + pstsol(k) - pst_rchrg(k,j) 

            if (pst_shallst(k,j) < 1.e-10) pst_shallst(k,j) = 0.0
            if (shallst(j) < 1.e-10) shallst(j) = 0.0
            if (gw_q(j) < 1.e-10) gw_q(j) = 0.0
            if (revapday < 1.e-10) revapday = 0.0
            if (gwseep < 1.e-10) gwseep = 0.0

            !! compute pesticide groundwater contribution to streamflow for day
            !xx = shallst(j) + gw_q(j) + revapday + gwseep !water volume in shallow aquifer (mm H2O)
            !xx = shallst(j)*0.02 + gw_q(j) + revapday + gwseep !water volume in shallow aquifer (mm H2O) 0.02
            xx = shallst(j)*pestgwfactor + gw_q(j) + revapday + gwseep !water volume in shallow aquifer (mm H2O) 0.02
            if (xx > 0.) then
                xx = pst_shallst(k,j) / xx                !pst concentration mg/ha/mm
            else
                xx = 0.
            end if
            if (xx < 1.e-9) xx = 0.0
            pst_gw(k,j) = xx * gw_q(j)                   !pst mass mg/ha

            revappst = xx * revapday                     !pst mass mg/ha
            gwseeppst = xx * gwseep                      !pst mass mg/ha

            revappst = amax1(0.,revappst)
            gwseeppst = amax1(0.,gwseeppst)
            revappst = 0 !!HR set pst revap to 0

            !! subtract pesticide losses from the shallow aquifer
            pst_shallst(k,j) = pst_shallst(k,j) - pst_gw(k,j)           !pst mass mg/ha
     &                                 - revappst - gwseeppst
            pst_shallst(k,j) = amax1 (0., pst_shallst(k,j))


            !!Deep aquifer
            !! compute pesticide deep aquifer loading from seepage for current day
            pst_deepst(k,j) =  pst_deepst(k,j) + gwseeppst             !pst mass mg/ha
            !! compute pesticide deep aquifer groundwater contribution to streamflow for day
            xx = deepst(j)*pestgwdeepfactor + gw_qdeep(j)   
            if (xx > 0.) then
                xx = pst_deepst(k,j) / xx                            !pst concentration mg/ha/mm
            else
                xx = 0.
            end if
            if (xx < 1.e-10) xx = 0.0
            pst_gwdeep(k,j) = xx * gw_qdeep(j)                       !pst mass mg/ha

            !! subtract pesticide losses from the deep 
            pst_deepst(k,j) = pst_deepst(k,j) - pst_gwdeep(k,j)      !pst mass mg/ha
            pst_deepst(k,j) = amax1 (0., pst_deepst(k,j))

            !! compute pesticide losses / decay in the shallow aquifer
            !! turned off because most chemicals are stable to hydrolysis
            !x1 = pst_shallst(k,j)
            !xx = x1 * decay_s(kk)
            !wshd_pstdg(k) = wshd_pstdg(k) + (x1 - xx) * hru_dafr(j)
            !pst_shallst(k,j) = amax1(0., xx)

            !! compute pesticide losses / decay in the deep aquifer
            !! ...

        end do
      end if
 
      return
      end