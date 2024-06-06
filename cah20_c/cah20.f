!***********************************************************************
!*****************************  KEK, High Energy Accelerator Research  *
!*****************************  Organization                           *
!*** u c n a i c g v *********                                         *
!*****************************     EGS5.0 USER CODE - 03 May 2006/1600 *
!***********************************************************************
!* This is a general User Code based on the cg geometry scheme.        *
!***********************************************************************
!                                                                      *
!  PROGRAMMERS:  H. Hirayama                                           *
!                Applied Research Laboratory                           *
!                KEK, High Energy Accelerator Research Organization    *
!                1-1, Oho, Tsukuba, Ibaraki, 305-0801                  *
!                Japan                                                 *
!                                                                      *
!                E-mail:      hideo.hirayama@kek.jp                    *
!                Telephone:   +81-29-864-5451                          *
!                Fax:         +81-29-864-4051                          *
!                                                                      *
!***********************************************************************
!***********************************************************************
! The ucnaicgv.f User Code requires a cg-input file only               *
! (e.g., ucnaicgv.data).                                               *
! The following shows the geometry for ucnaicg.data.                   *
! Input data for CG geometry must be written at the top of data-input  *
! file together with material assignment to each region.  Cg-data can  *
! be checked by CGview.                                                *
! This user code corresponds to ucnai3cgp.mor for egs4.                *
! Use Ranlux random number generator.                                  *
!***********************************************************************
!                                                                      *
!             -----------------------                                  *
!             cg Geometry (ucnaicgv)                                   *
!             -----------------------                                  *
!                                                                      *
!                                                                      *
!                R                                                     *
!                ^                                                     *
!                |                                                     *
!           +----+----+----+--------+------+---                        *
!           |                                                          *
!           |     Outer vacuum region                                  *
!           +    +----+----+--------+------+-----+ R=9.41              *
!           |    |         Air                   |                     *
!           |    |    +--------------------+     + R=4.41              *
!           |    |    |    Al cover        |     |                     *
!           |    +    +    +-----------+---+     + R=4.31              *
!           |    |    |    | Gap       |   |     |                     *
!           +    +    +    +   +-------+   +     + R=3.81              *
!           |    |    |    |   |       |Quartz   |                     *
!           |    |    |    |   |  NaI  |   |     |                     *
! 1.253 MeV |    |    |    |   |       |   |     |                     *
!      ============>--+----+---+-------+---+-----+--------> Z          *
!     photons  -5.6 -0.6 -0.5 0.0    7.62 8.12  13.12                  *
!                -5.0                                                  *
!                                                                      *
!                                                                      *
!                                                                      *
!***********************************************************************
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

!-----------------------------------------------------------------------
!------------------------------- main code -----------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Step 1: Initialization
!-----------------------------------------------------------------------

      implicit none

!     ------------
!     EGS5 COMMONs
!     ------------
      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_bounds.f'
      include 'include/egs5_brempr.f'
      include 'include/egs5_edge.f'
      include 'include/egs5_media.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_thresh.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_usersc.f'
      include 'include/egs5_userxt.f'
      include 'include/randomm.f'

!     ----------------------
!     Auxiliary-code COMMONs
!     ----------------------
      include 'auxcommons/aux_h.f'   ! Auxiliary-code "header" file

      include 'auxcommons/edata.f'
      include 'auxcommons/etaly1.f'
      include 'auxcommons/instuf.f'
      include 'auxcommons/lines.f'
      include 'auxcommons/nfac.f'
      include 'auxcommons/watch.f'

      include 'auxcommons/etaly2.f'      !  Added SJW for energy balance

!     ------------------
!     cg related COMMONs
!     ------------------
      include 'auxcommons/geom_common.f' ! geom-common file
      integer irinn

      common/totals/                          ! Variables to score
     * depe,deltae,spg(1,50),spe(1,50),spp(1,50),imode
      real*8 depe,deltae,spg,spe,spp
      integer imode

!**** real*8                                                 ! Arguments
      real*8 totke
      real*8 rn_disc1,rn_disc2,etot,rn1,rn2
      real*8 esumt
      
      real*8                                           ! Local variables
     * availke,avpe,avph,avspe,avspg,avspp,avte,desci2,ekin,pef,
     * rr0,sigpe,sigte,sigph,sigspg,sigspe,sigspp,tef,wtin,wtsum,
     * xi0,yi0,zi0

      real*8
     * ph(50),phpb(50,50),spgpb(1,50,50),spepb(1,50,50),
     * spppb(1,50,50),pefpb(50),tefpb(50)
     
      real                                             ! Local variables
     * elow,eup,rdet,rtcov,rtgap,tcov,tdet,tgap

      real
     * tarray(2),tt,tt0,tt1,cputime

      integer
     * i,icases,idin,ie,ifti,ifto,ii,iiz,imed,ireg,isam,
     * isot,izn,nlist,j,k,n,ndet,nd,nbatch,ncaspb,
     * ner,nofbat

      character*24 medarr(3)

!     ----------
!     Open files
!     ----------
!----------------------------------------------------------------
!     Units 7-26 are used in pegs and closed.  It is better not
!     to use as output file. If they are used, they must be opened
!     after getcg etc. Unit for pict must be 39.
!----------------------------------------------------------------
      open(2,FILE='positions.out',STATUS='unknown')
      open(3,FILE='trajectories.out',STATUS='unknown')
      open(1,FILE='egs5job.out',STATUS='unknown')
      open(UNIT= 4,FILE='egs5job.inp',STATUS='old')
      open(39,FILE='egs5job.pic',STATUS='unknown')

!     ====================
      call counters_out(0)
!     ====================

!-----------------------------------------------------------------------
! Step 2: pegs5-call
!-----------------------------------------------------------------------
!     ==============
      call block_set                 ! Initialize some general variables
!     ==============

!     ---------------------------------
!     Define media before calling PEGS5
!     ---------------------------------

      nmed=3 !# number of materials is 2
      if (nmed.gt.MXMED) then
        write(1,100) nreg,mxreg 
100     FORMAT(' nmed(=',I12,') must be less than MXMED(=',I12,')' /
     *  ' You must change MXMED in include/egs5_h.f.')
        stop
      end if

      medarr(1)='CA                      ' !# new materials
      medarr(2)='H2O                     '
      medarr(3)='AIR-AT-NTP              '

      do j=1,nmed
        do i=1,24
          media(i,j)=medarr(j)(i:i)
        end do
      end do  

      chard(1) = 7.00  ! # automatic step-size control
      chard(2) = 2.0d0
      chard(3) = 0.1d0

      write(1,*) 'chard =',(chard(j),j=1,3)

!     -----------------------------------
!     Run KEK PEGS5 before calling HATCH
!     -----------------------------------
      !write(1,109)
!109  !FORMAT('testing custom outputs'/) !# custom output
      write(1,110)
110   FORMAT(' PEGS5-call comes next'/)

!     ==========
      call pegs5
!     ==========

!-----------------------------------------------------------------------
! Step 3: Pre-hatch-call-initialization
!-----------------------------------------------------------------------
!-----------------------------------------------
!     Initialize cg related parameter 
!-----------------------------------------------
      npreci=2     ! PICT data mode for CGView

      itbody=0
      irppin=0
      isphin=0
      irccin=0
      itorin=0
      itrcin=0
      izonin=0
      itverr=0
      igmmax=0
      ifti = 4     ! Input unit number for cg-data
      ifto = 39    ! Output unit number for PICT

      write(39,120)
120   FORMAT('CSTA')
      call geomgt(ifti,ifto)
      write(39,130)
130   FORMAT('CEND')

!--------------------------------
!     Get nreg from cg input data
!--------------------------------
      nreg=izonin
      if (nreg.gt.mxreg) then
        write(1,140) nreg,mxreg 
140     FORMAT(' NREG(=',I12,') must be less than MXREG(=',I12,')' /
     *  ' You must change MXREG in include/egs5_h.f.')
        stop
      end if

!   Set medium index for each region
      med(1)=1      ! # medium one is material 1 (Ca)					
      iedgfl(1)=1   ! 1:Produce fluorescent X-rays
                    ! 0:Fluorescent X-ray is not produced

      med(2)=2      ! # medium two is material 2 (H2O)					
      !iedgfl(2)=1   ! 1:Produce fluorescent X-rays	
      !              ! 0:Fluorescent X-ray is not produced

      med(3)=3      ! # medium two is material 3 (Air)				
      !iedgfl(3)=1   ! 1:Produce fluorescent X-rays	
      !              ! 0:Fluorescent X-ray is not produced
      		
      med(nreg)=0   ! Vacuum region		

!     do i=1,nreg
!       if(i.eq.1) ecut(i)=0.561
!     end do

!     --------------------------------------------------------
!     Random number seeds.  Must be defined before call hatch
!     or defaults will be used.  inseed (1- 2^31)
!     --------------------------------------------------------
      luxlev = 1
      inseed=1
      write(1,160) inseed
160   FORMAT(/,' inseed=',I12,5X,
     *         ' (seed for generating unique sequences of Ranlux)')

!     =============
      call rluxinit  ! Initialize the Ranlux random-number generator
!     =============

!-----------------------------------------------------------------------
! Step 4:  Determination-of-incident-particle-parameters
!-----------------------------------------------------------------------
! Define initial variables for incident particle normally incident
! on the slab

      iqin=0             ! Incident particle charge - photons
      !ekein=1.253       ! Incident particle kinetic energy
	  ekein = 1.000      ! # New incident particle energy (1 MeV)
      xin=0.0            ! Source position
      yin=0.0
	  zin = -6.0         !# z-position of new incident particle
      uin=0.0            ! Moving along z axis
      vin=0.0
      win=1.0
      irin=0             ! Starting region (0: Automatic search in CG)
      wtin=1.0           ! Weight = 1 since no variance reduction used

!     pdf data for many source
      deltae=0.05        ! Energy bin of response

!-----------------------------------------------------------------------
! Step 5:   hatch-call
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Maximum total energy of an electron for this problem must be
! defined before hatch call
      emaxe = ekein + RM         

      write(1,170)
170   format(/' Start ucnaicgv '/
     *' Call hatch to get cross-section data')

!     ------------------------------
!     Open files (before HATCH call)
!     ------------------------------
      open(UNIT=KMPI,FILE='pgs5job.pegs5dat',STATUS='old')
      open(UNIT=KMPO,FILE='egs5job.dummy',STATUS='unknown')

      write(1,180)
180   FORMAT(/,' HATCH-call comes next',/)
   

!     ==========
!     write(1,185)
!185   FORMAT(/,' test test before hatch',/)
!     ==========


!     ==========
      call hatch
!     ==========

!     ------------------------------
!     Close files (after HATCH call)
!     ------------------------------
      close(UNIT=KMPI)
      close(UNIT=KMPO)

! ----------------------------------------------------------
! Print various data associated with each media (not region)
! ----------------------------------------------------------
      write(1,190)
190   FORMAT(/,' Quantities associated with each MEDIA:')
      do j=1,nmed
        write(1,200) (media(i,j),i=1,24)
200     FORMAT(/,1X,24A1)
        write(1,210) rhom(j),rlcm(j)
210     FORMAT(5X,' rho=',G15.7,' g/cu.cm     rlc=',G15.7,' cm')
        write(1,220) ae(j),ue(j)
220     FORMAT(5X,' ae=',G15.7,' MeV    ue=',G15.7,' MeV')
        write(1,230) ap(j),up(j)
230     FORMAT(5X,' ap=',G15.7,' MeV    up=',G15.7,' MeV',/)
      end do

! -------------------------------------------------------
! Print media and cutoff energies assigned to each region
! -------------------------------------------------------
      do i=1,nreg
        if (med(i) .eq. 0) then
          write(1,240) i
240       FORMAT(' medium(',I3,')=vacuum')
        else
          write(1,250) i,(media(ii,med(i)),ii=1,24),ecut(i),pcut(i)
250       FORMAT(' medium(',I3,')=',24A1,
     *           'ecut=',G10.5,' MeV, pcut=',G10.5,' MeV') !#check here
!        -----------------------------------------------
!        Print out energy information of K- and L-X-rays
!        -----------------------------------------------
          if (iedgfl(i) .ne. 0) then             ! Output X-ray energy
            ner = nne(med(i))
            do iiz=1,ner
              izn = zelem(med(i),iiz)  ! Atomic number of this element
              write(1,260) izn
260           FORMAT('   X-ray information for Z=',I3)
              write(1,270) (ekx(ii,izn),ii=1,10)
270           FORMAT('   K-X-ray energy in keV',/,
     *               4G15.5,/,4G15.5,/,2G15.5)
              write(1,280) (elx1(ii,izn),ii=1,8)
280           FORMAT('   L-1 X-ray in keV',/,4G15.5,/,4G15.5)
              write(1,290) (elx2(ii,izn),ii=1,5)
290           FORMAT('   L-2 X-ray in keV',/,5G15.5)
              write(1,300) (elx3(ii,izn),ii=1,7)
300           FORMAT('   L-3 X-ray in keV',/,4G15.5,/,3G15.5)
            end do
          end if
        end if
      end do

      write(39,310)
310   FORMAT('MSTA')
      write(39,320) nreg
320   FORMAT(I4)
      write(39,330) (med(i),i=1,nreg)
330   FORMAT(15I4)
      write(39,340)
340   FORMAT('MEND')

!-----------------------------------
!     Selection mode from Keyboard.
!-----------------------------------
      write(6,350)
350   FORMAT(' Key in mode. 0:trajectory display, 1:dose calculation')
      read(5,*) imode

!-----------------------------------------------------------------------
! Step 6:  Initialization-for-howfar
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Step 7:  Initialization-for-ausgab
!-----------------------------------------------------------------------

      ncount = 0
      ilines = 0
      nwrite = 10
      nlines = 10
      idin = -1
      totke = 0.
      wtsum = 0.

!     =========================
      call ecnsv1(0,nreg,totke)
      call ntally(0,nreg)
!     =========================

      write(1,360)
360   format(/,' Energy/coordinates/direction cosines/etc.',/,
     *        6X,'e',16X,'x',14X,'y',14X,'z'/
     *        1X,'u',14X,'v',14X,'w',9X,'iq',4X,'ir',3X,'iarg',/)

      ndet=1

!     Energy bin width 
      deltae=ekein / 50
      
!     Zero the variables
      depe=0.D0
      pef=0.D0
      tef=0.D0
      do j=1,50
        ph(j)=0.D0
        do nd=1,ndet
          spg(nd,j)=0.D0
          spe(nd,j)=0.D0
          spp(nd,j)=0.D0
        end do
      end do

!     Set histories, number of batch and histories per batch
      write(6,*) (' Key in number of cases.')
      read(5,*) ncases
      if(imode.eq.0) then
         write(6,*) (' Key in number of batch (=<50).')
         read(5,*)  nbatch
      else
        nbatch = 50
      end if
      ncaspb = ncases / nbatch !# HAVETO DO > 50 CASES
	  
      !write(6,*) (ncaspb) !# custom print to console

      tt=etime(tarray)
      tt0=tarray(1)
	  
      !write(6,*) (' Starting shower call.') !# custom print to console

!-----------------------------------------------------------------------
! Step 8:  Shower-call
!-----------------------------------------------------------------------
                                             ! -------------------------
      do nofbat=1,nbatch                     ! Start of batch -loop
                                             ! -------------------------
         write(39,370) nofbat
370      format('0',I5)                     

         do icases=1,ncaspb                  ! Start of CALL SHOWER loop

!       ----------------------
!       Select incident energy
!       ----------------------
          eparte = 0.d0                   ! Initialize some energy-balance
          epartd = 0.d0                   !      tallying parameters (SJW)

          ekin = ekein
          wtin = 1.0
			
          wtsum = wtsum + wtin               ! Keep running sum of weights
          etot = ekin + iabs(iqin)*RM        ! Incident total energy (MeV)
          availke = etot + iqin*RM        ! Available K.E. (MeV) in system 
          totke = totke + availke                 ! Keep running sum of KE
		  
      ! #Additional print statements for debugging
      !write(6, *) 'availke:', availke
      !write(6, *) 'icases:', icases
      !write(6, *) 'nofbat:', nofbat
      !write(6, *) 'wtsum:', wtsum
      !write(6, *) 'totke:', totke
! ----------------------
! Select incident angle 
! ---------------------- 
          
  
          call random_number(rn1)  ! Generate a random number between 0 and 1
          call random_number(rn2)  ! Generate another random number between 0 and 1
          uin = SIN(ACOS(2.0*rn1-1.0))*COS(2.0*3.14*rn2)
          vin = SIN(ACOS(2.0*rn1-1.0))*SIN(2.0*3.14*rn2)
          win = COS(ACOS(2.0*rn1-1.0))

          !# random sampling of x,y position in a disc 4 cm diameter
          call random_number(rn_disc1)
          call random_number(rn_disc2)
          xin = 2*SIN(ACOS(2.0*rn_disc1-1.0))*COS(2.0*3.14*rn_disc2)
          yin = 2*SIN(ACOS(2.0*rn_disc1-1.0))*SIN(2.0*3.14*rn_disc2)

!       ------------------------------------
!       Get source region from cg input data
!       ------------------------------------

          if(irin.le.0.or.irin.gt.nreg) then
            call srzone(xin,yin,zin,iqin+2,0,irinn)
            call rstnxt(iqin+2,0,irinn)
          else
            irinn=irin
          end if    
	  
!       ---------------------------------------------------
!       Print first NWRITE or NLINES, whichever comes first
!       ---------------------------------------------------
          if (ncount .le. nwrite .and. ilines .le. nlines) then
            ilines = ilines + 1
            write(1,390) etot,xin,yin,zin,uin,vin,win,iqin,irinn,idin
390         FORMAT(4G15.7/3G15.7,3I5)
          end if
!         ==========================================================
          ! #Additional print statements for debugging the segfault
          !write(6, *) 'iqin:', iqin 
          !write(6, *) 'etot:', etot 
          !write(6, *) 'xin:',  xin  
          !write(6, *) 'yin:',  yin  
          !write(6, *) 'zin:',  zin  
          !write(6, *) 'uin:',  uin  
          !write(6, *) 'vin:',  vin  
          !write(6, *) 'win:',  win  
          !write(6, *) 'irinn:',irinn
          !write(6, *) 'irin:',irin
          !write(6, *) 'wtin:', wtin 
          write(3, *) 'uin,vin,win=',uin,vin,win
          write(2, *) 'xin,yin,zin=',xin,yin,zin
          call shower (iqin,etot,xin,yin,zin,uin,vin,win,irinn,wtin)
          
          !         ==========================================================

!       Added for energy balance tests (SJW)
          if(DABS(eparte + epartd - ekin)/ekin .gt. 1.d-10) then
            write(1,400) icases, eparte, epartd
400         FORMAT('Error on # ',I6,' Escape = ',F9.5,
     *             ' Deposit = ',F9.5)
          end if

       !write(6, *) 'Checking for Energy, depe:', depe 

!      If some energy is deposited inside detector add pulse-height
!      and efficiency.

          if (depe .gt. 0.D0) then
            !write(6, *) 'Adding pulse height, depe:', depe
            ie=depe/deltae + 1
            if (ie .gt. 50)  ie = 50
            ph(ie)=ph(ie)+wtin
            tef=tef + wtin
            if(depe .ge. ekein*0.999) pef=pef +wtin
            depe = 0.D0
          end if

          ncount = ncount + 1         ! Count total number of actual cases

                                               ! -----------------------
        end do                                 ! End of CALL SHOWER loop
                        
        close(3,status='keep')               
        close(2,status='keep')        

!  Calculate average value for this BATCH
        do ie=1,50
          phpb(ie,nofbat) = ph(ie) /ncaspb
          ph(ie)=0.D0
        end do
        pefpb(nofbat)=pef / ncaspb
        tefpb(nofbat)=tef /ncaspb
        pef=0.D0
        tef=0.D0
        do nd=1,ndet
          do ie=1,50
            spgpb(nd,ie,nofbat)=spg(nd,ie)/ncaspb !photon spectrum
            spepb(nd,ie,nofbat)=spe(nd,ie)/ncaspb !electron spectrum
            spppb(nd,ie,nofbat)=spp(nd,ie)/ncaspb !positron spectrum
            spg(nd,ie)=0.D0
            spe(nd,ie)=0.D0
            spp(nd,ie)=0.D0
          end do
        end do

        call plotxyz(99,0,0,0.D0,0.D0,0.D0,0.D0,0,0.D0) 

        write(39,410)               ! Set end of batch for CG View
410     FORMAT('9')

                                               ! ------------------
      end do                                   ! End of batch loop
                                               ! -------------------
      tt=etime(tarray)
      tt1=tarray(1)
      cputime=tt1-tt0
      write(1,420) cputime
420   format(' Elapsed Time (sec)=',G15.5)

!-----------------------------------------------------------------------
! Step 9:  Output-of-results
!-----------------------------------------------------------------------
      write(1,430) ncount,ncases,totke
430   FORMAT(/,' Ncount=',I10,' (actual cases run)',/,
     *       ' Ncases=',I10,' (number of cases requested)',/,
     *       ' TotKE =',G15.5,' (total KE (MeV) in run)')

      if (totke .le. 0.D0) then
        write(1,440) totke,availke,ncount
440     FORMAT(//,' Stopped in MAIN with TotKE=',G15.5,/,
     *         ' AvailKE=',G15.5, /,' Ncount=',I10)
        stop
      end if
      !# might need to change
      tdet=7.62
      rdet=3.81
      tcov=0.1
      rtcov=0.1
      tgap=0.5
      rtgap=0.5
      write(1,450) tdet,rdet,tcov,rtcov,tgap,rtgap
450   FORMAT(/' Detector length=',G15.5,' cm'/
     *       ' Detector radius=',G15.5,' cm'/
     *       ' Al cover thickness=',G10.2,' cm'/
     *       ' Al cover side thickness=',G10.2,' cm'/
     *       ' Front gap =',G10.2,' cm'/' Side gap =',G10.2,' cm'/)
     
      write(1,460) ekin
460   FORMAT(' Results for ',G15.5,'MeV photon'/)
      
!     -----------------------------------
!     Calculate average and its deviation
!     -----------------------------------

!     ---------------
!     Peak efficiency
!     ---------------
      avpe = 0.D0
      desci2 = 0.D0
      do j = 1, nbatch
        avpe = avpe + pefpb(j)/nbatch
        desci2 = desci2 + pefpb(j)*pefpb(j)/nbatch
      end do
      sigpe = sqrt((desci2 - avpe*avpe)/(nbatch-1))
      avpe = avpe*100.0
      sigpe = sigpe*100.0
      write(1,470) avpe,sigpe
!470   FORMAT(' Peak efficiency =',G15.5,'+-',G15.5,' %')
470   FORMAT(' Peak efficiency =',G11.4,'+-',G9.2,' %')

!     ----------------
!     Total efficiency
!     ----------------
      avte = 0.D0
      desci2 = 0.D0
      do j = 1, nbatch
        avte = avte + tefpb(j)/nbatch
        desci2 = desci2 + tefpb(j)*tefpb(j)/nbatch
      end do
      sigte = sqrt((desci2 - avte*avte)/(nbatch-1))
      avte = avte*100.0
      sigte = sigte*100.0
      write(1,480) avte,sigte
!80   FORMAT(' Total efficiency =',G15.5,'+-',G15.5,' %')
480   FORMAT(' Total efficiency =',G11.4,'+-',G9.2,' %')

!     --------------------------
!     Pulse height distribution
!     --------------------------
      write(1,490)
490   FORMAT(/' Pulse height distribution ')
      do ie=1,50
        elow=deltae*(ie-1)
        eup=deltae*ie
        if (elow .gt. ekein ) go to 510
        
        avph = 0.D0
        desci2 = 0.D0
        do j = 1, nbatch
          avph = avph + phpb(ie,j)/nbatch
          desci2 = desci2 + phpb(ie,j)*phpb(ie,j)/nbatch
        end do
        sigph = sqrt((desci2 - avph*avph)/(nbatch-1))
        avph = avph/deltae
        sigph= sigph/deltae
        write(1,500) eup,avph,sigph
500     FORMAT(' E (upper-edge --',G10.4,' MeV )=',G15.5,'+-',G15.5,
     *         ' counts/MeV/incident');
       end do

510    continue

!     ----------------------------------------------------------
!     Particle spectrum.  Incident particle spectrum to detector.
!     ----------------------------------------------------------
      write(1,520)
520   FORMAT(/' Particle spectrum crossing the detector plane'/
     *       30X,'particles/MeV/source photon'/
     *       ' Upper energy',11X,'  Gamma',18X,' Electron',
     *       14X,' Positron')
     
      do nd=1,ndet
        do ie=1,50
          elow=deltae*(ie-1)
          eup=deltae*ie
          if (elow .gt. ekein ) go to 540
          
!     ----------------------------------
!     Gamma spectrum per MeV per source
!     ----------------------------------
         
          avspg = 0.D0
          desci2 = 0.D0
          do j = 1, nbatch
            avspg = avspg + spgpb(nd,ie,j)/nbatch
            desci2 = desci2 + spgpb(nd,ie,j)*spgpb(nd,ie,j)/nbatch
          end do
          sigspg = sqrt((desci2 - avspg*avspg)/(nbatch-1))
          avspg = avspg/deltae
          sigspg= sigspg/deltae

!     -------------------------------------
!     Electron spectrum per MeV per source
!     -------------------------------------
         
          avspe = 0.D0
          desci2 = 0.D0
          do j = 1, nbatch
            avspe = avspe + spepb(nd,ie,j)/nbatch
            desci2 = desci2 + spepb(nd,ie,j)*spepb(nd,ie,j)/nbatch
          end do
          sigspe = sqrt((desci2 - avspe*avspe)/(nbatch-1))
          avspe = avspe/deltae
          sigspe= sigspe/deltae
       
!     ------------------------------------
!     Positron spectrum per MeV per source
!     ------------------------------------
         
          avspp = 0.D0
          desci2 = 0.D0
          do j = 1, nbatch
            avspp = avspp + spppb(nd,ie,j)/nbatch
            desci2 = desci2 + spppb(nd,ie,j)*spppb(nd,ie,j)/nbatch
          end do
          sigspp = sqrt((desci2 - avspp*avspp)/(nbatch-1))
          avspp = avspp/deltae
          sigspp= sigspp/deltae

          write(1,530) eup,avspg,sigspg,avspe,sigspe,avspp,sigspp
530       FORMAT(G10.5,' MeV--',3(G12.5,'+-',G12.5))
        end do
      end do

540   continue

      nlist=1
      
      if(imode.ne.0) then
        close(1,status='keep')
        open(6,file='egs5job.out',access='append')

!       =============================
        call ecnsv1(nlist,nreg,totke)
        call ntally(nlist,nreg)
!       =============================
      end if

600   continue
!     ====================
      call counters_out(1)
!     ====================

      stop

      end

!-------------------------last line of main code------------------------

!-------------------------------ausgab.f--------------------------------
! Version:   030831-1300
! Reference: SLAC-265 (p.19-20, Appendix 2)
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

! ----------------------------------------------------------------------
! Required subroutine for use with the EGS5 Code System
! ----------------------------------------------------------------------
! A AUSGAB to:
!
!   1) Score energy deposition
!   2) Score particle information enter to detector from outside
!   3) Print out particle transport information 
!   4) call plotxyz if imode=0

! ----------------------------------------------------------------------

      subroutine ausgab(iarg)

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_useful.f'

      include 'auxcommons/aux_h.f'   ! Auxiliary-code "header" file

      include 'auxcommons/etaly1.f'        ! Auxiliary-code COMMONs
      include 'auxcommons/lines.f'
      include 'auxcommons/ntaly1.f'
      include 'auxcommons/watch.f'

      include 'auxcommons/etaly2.f'      !  Added SJW for energy balance

      common/totals/                          ! Variables to score
     * depe,deltae,spg(1,50),spe(1,50),spp(1,50),imode
      real*8 depe,deltae,spg,spe,spp
      integer imode

      integer                                                ! Arguments
     * iarg

      real*8                                           ! Local variables
     * edepwt
     
      integer 
     * ie,iql,irl

!     ------------------------
!     Set some local variables
!     ------------------------
      irl = ir(np)
      iql = iq(np)
      edepwt = edep*wt(np)

!     -----------------------------------------------------------
!     Keep track of energy deposition (for conservation purposes)
!     -----------------------------------------------------------
      if (iarg .lt. 5) then !# custom change
        esum(iql+2,irl,iarg+1) = esum(iql+2,irl,iarg+1) + edepwt
        nsum(iql+2,irl,iarg+1) = nsum(iql+2,irl,iarg+1) + 1

        
!  added SJW for particle by particle energy balance
        if(irl.eq.nreg) then
          eparte = eparte + edepwt
        else 
          epartd = epartd + edepwt
        endif
      end if

!     ----------------------------------------------
!     Score energy deposition inside NaI detector
!     ----------------------------------------------
      
      !write(6,*) 'med(irl):',med(irl)!#
      if (med(irl). eq. 1) then
        depe = depe + edepwt
        !write(6,*) 'Energy Deposited:',depe
        
!      ----------------------------------------------------
!      Score particle information if it enters from outside
!      ----------------------------------------------------
        if (irl .ne. irold .and. iarg .eq. 0) then
          if (iql .eq. 0) then             ! photon
            ie = e(np)/deltae +1
            if(ie .gt. 50) ie = 50
            spg(1,ie) = spg(1,ie) + wt(np)
          elseif (iql .eq. -1) then        ! electron
            ie = (e(np) - RM)/deltae +1
            if(ie .gt. 50) ie = 50
            spe(1,ie) = spe(1,ie) + wt(np)
          else                             ! positron
            ie = (e(np) - RM)/deltae +1
            if(ie .gt. 50) ie = 50
            spp(1,ie) = spp(1,ie) + wt(np)
          end if
        end if         
      end if


!     ----------------------------------------------------------------
!     Print out stack information (for limited number cases and lines)
!     ----------------------------------------------------------------
      if (ncount .le. nwrite .and. ilines .le. nlines) then
        ilines = ilines + 1
        write(1,100) e(np),x(np),y(np),z(np),u(np),v(np),w(np),
     *               iql,irl,iarg
 100    FORMAT(4G15.7/3G15.7,3I5)
      end if

!     -----------------------------------------------------------------
!     Print out particle transport information (if switch is turned on)
!     -----------------------------------------------------------------
!                        ========================
      if (iwatch .gt. 0) call swatch(iarg,iwatch)
!                        ========================

!     ------------------------------------
!     Output particle information for plot
!     ------------------------------------
      if (imode.eq.0) then
        call plotxyz(iarg,np,iq(np),x(np),y(np),z(np),e(np),ir(np),
     *       wt(np))
      end if


      return

      end

!--------------------------last line of ausgab.f------------------------
!-------------------------------howfar.f--------------------------------
! Version:   050716-1300
! Reference: T. Torii and T. Sugita, "Development of PRESTA-CG 
! Incorporating Combinatorial Geometry in EGS4/PRESTA", JNC TN1410 2002-201,
! Japan Nuclear Cycle Development Institute (2002).
! Improved version is provided by T. Sugita. 7/27/2004
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

! ----------------------------------------------------------------------
! Required (geometry) subroutine for use with the EGS5 Code System
! ----------------------------------------------------------------------
! This is a CG-HOWFAR. 
! ----------------------------------------------------------------------

      subroutine howfar
      implicit none
c
      include 'include/egs5_h.f'       ! Main EGS "header" file
      include 'include/egs5_epcont.f'  ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'
      include 'auxcommons/geom_common.f' ! geom-common file
c
c
      integer i,j,jjj,ir_np,nozone,jty,kno
      integer irnear,irnext,irlold,irlfg,itvlfg,ihitcg
      double precision xidd,yidd,zidd,x_np,y_np,z_np,u_np,v_np,w_np
      double precision tval,tval0,tval00,tval10,tvalmn,delhow
      double precision atvaltmp
      integer iq_np
c
      ir_np = ir(np)
      iq_np = iq(np) + 2
c
      if(ir_np.le.0) then
        write(6,*) 'Stopped in howfar with ir(np) <=0'
        stop
      end if
c
      if(ir_np.gt.izonin) then
        write(6,*) 'Stopped in howfar with ir(np) > izonin'
        stop
      end if
c
      if(ir_np.EQ.izonin) then
        idisc=1
        return
      end if
c
      tval=1.d+30
      itvalm=0
c
c     body check
      u_np=u(np)
      v_np=v(np)
      w_np=w(np)
      x_np=x(np)
      y_np=y(np)
      z_np=z(np)
c
      do i=1,nbbody(ir_np)
        nozone=ABS(nbzone(i,ir_np))
        jty=itblty(nozone)
        kno=itblno(nozone)
c     rpp check
        if(jty.eq.ityknd(1)) then
          if(kno.le.0.or.kno.gt.irppin) go to 190
          call rppcg1(kno,x_np,y_np,z_np,u_np,v_np,w_np)
c     sph check
        elseif(jty.eq.ityknd(2)) then
          if(kno.le.0.or.kno.gt.isphin) go to 190
          call sphcg1(kno,x_np,y_np,z_np,u_np,v_np,w_np)
c     rcc check
        elseif(jty.eq.ityknd(3)) then
          if(kno.le.0.or.kno.gt.irccin) go to 190
          call rcccg1(kno,x_np,y_np,z_np,u_np,v_np,w_np)
c     trc check
        elseif(jty.eq.ityknd(4)) then
          if(kno.le.0.or.kno.gt.itrcin) go to 190
          call trccg1(kno,x_np,y_np,z_np,u_np,v_np,w_np)
c     tor check
        elseif(jty.eq.ityknd(5)) then
          if(kno.le.0.or.kno.gt.itorin) go to 190
          call torcg1(kno,x_np,y_np,z_np,u_np,v_np,w_np)
c
c**** add new geometry in here
c
       end if
  190  continue
      end do
c
      irnear=ir_np
      if(itvalm.eq.0) then
        tval0=cgeps1
        xidd=x_np+tval0*u_np
        yidd=y_np+tval0*v_np
        zidd=z_np+tval0*w_np
  310   continue
          if(x_np.ne.xidd.or.y_np.ne.yidd.or.z_np.ne.zidd) goto 320
          tval0=tval0*10.d0
          xidd=x_np+tval0*u_np
          yidd=y_np+tval0*v_np
          zidd=z_np+tval0*w_np
          go to 310
  320   continue
c       write(*,*) 'srzone:1'
        call srzone(xidd,yidd,zidd,iq_np,ir_np,irnext)
c
        if(irnext.ne.ir_np) then
          tval=0.0d0
          irnear=irnext
        else
          tval00=0.0d0
          tval10=10.0d0*tval0
          irlold=ir_np
          irlfg=0
  330     continue
          if(irlfg.eq.1) go to 340
            tval00=tval00+tval10
            if(tval00.gt.1.0d+06) then
              write(6,9000) iq(np),ir(np),x(np),y(np),z(np),
     &                      u(np),v(np),w(np),tval00
 9000 format(' TVAL00 ERROR : iq,ir,x,y,z,u,v,w,tval=',
     &       2I3,1P7E12.5)
              stop
            end if
            xidd=x_np+tval00*u_np
            yidd=y_np+tval00*v_np
            zidd=z_np+tval00*w_np
            call srzold(xidd,yidd,zidd,irlold,irlfg)
            go to 330
  340     continue
c
          tval=tval00
          do j=1,10
            xidd=x_np+tval00*u_np
            yidd=y_np+tval00*v_np
            zidd=z_np+tval00*w_np
c           write(*,*) 'srzone:2'
            call srzone(xidd,yidd,zidd,iq_np,irlold,irnext)
            if(irnext.ne.irlold) then
              tval=tval00
              irnear=irnext
            end if
            tval00=tval00-tval0
          end do
          if(ir_np.eq.irnear) then
            write(0,*) 'ir(np),tval=',ir_np,tval
          end if
        end if
      else
        do j=1,itvalm-1
          do i=j+1,itvalm
            if(atval(i).lt.atval(j)) then
              atvaltmp=atval(i)
              atval(i)=atval(j)
              atval(j)=atvaltmp
            endif
          enddo
        enddo
        itvlfg=0
        tvalmn=tval
        do jjj=1,itvalm
          if(tvalmn.gt.atval(jjj)) then
            tvalmn=atval(jjj)
          end if
          delhow=cgeps2
          tval0=atval(jjj)+delhow
          xidd=x_np+tval0*u_np
          yidd=y_np+tval0*v_np
          zidd=z_np+tval0*w_np
  410     continue
          if(x_np.ne.xidd.or.y_np.ne.yidd.or.z_np.ne.zidd) go to 420
            delhow=delhow*10.d0
            tval0=atval(jjj)+delhow
            xidd=x_np+tval0*u_np
            yidd=y_np+tval0*v_np
            zidd=z_np+tval0*w_np
          go to 410
  420     continue
c         write(*,*) 'srzone:3'
          call srzone(xidd,yidd,zidd,iq_np,ir_np,irnext)
          if((irnext.ne.ir_np.or.atval(jjj).ge.1.).and.
     &        tval.gt.atval(jjj)) THEN
            tval=atval(jjj)
            irnear=irnext
            itvlfg=1
            goto 425
          end if
        end do
  425   continue
        if(itvlfg.eq.0) then
          tval0=cgmnst
          xidd=x_np+tval0*u_np
          yidd=y_np+tval0*v_np
          zidd=z_np+tval0*w_np
  430     continue
          if(x_np.ne.xidd.or.y_np.ne.yidd.or.z_np.ne.zidd) go to 440
            tval0=tval0*10.d0
            xidd=x_np+tval0*u_np
            yidd=y_np+tval0*v_np
            zidd=z_np+tval0*w_np
            go to 430
  440     continue
          if(tvalmn.gt.tval0) then
            tval=tvalmn
          else
            tval=tval0
          end if
        end if
      end if
      ihitcg=0
      if(tval.le.ustep) then
        ustep=tval
        ihitcg=1
      end if
      if(ihitcg.eq.1) THEN
        if(irnear.eq.0) THEN
          write(6,9200) iq(np),ir(np),x(np),y(np),z(np),
     &                  u(np),v(np),w(np),tval
 9200 format(' TVAL ERROR : iq,ir,x,y,z,u,v,w,tval=',2I3,1P7E12.5)
          idisc=1
          itverr=itverr+1
          if(itverr.ge.100) then
            stop
          end if
          return
        end if
        irnew=irnear
        if(irnew.ne.ir_np) then
          call rstnxt(iq_np,ir_np,irnew)
        endif
      end if
      return
      end
!--------------------last line of subroutine howfar---------------------

