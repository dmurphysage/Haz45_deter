      program haz45a

c     Probabilisitic Seismic Hazard Program (PSHA) 


      include 'pfrisk.h'
      include 'declare1.h'
      integer faultFlag(MAX_FLT,100,MAX_FLT), nDD(MAX_FLT)
      real segWt1(MAX_FLT)
      real fltGrid_X(MAXFLT_DD,MAXFLT_AS), fltGrid_y(MAXFLT_DD,MAXFLT_AS), 
     1     fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_fLen(MAXFLT_DD,MAXFLT_AS)
      real rupGrid_X(MAXFLT_DD,MAXFLT_AS), rupGrid_y(MAXFLT_DD,MAXFLT_AS), 
     1     rupGrid_z(MAXFLT_DD,MAXFLT_AS), hDD(MAX_FLT), hAS(MAX_FLT)
      integer nfltGrid(2), nRupGrid(2), hDDcell, hAScell
      real fltGrid_w(MAXFLT_DD,MAXFLT_AS),  fltGrid_a(MAXFLT_DD,MAXFLT_AS)
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_RJB(MAXFLT_DD,MAXFLT_AS)
      real testsum(1000), sum1(1000,10), dipaverage(1), Rx, Ry, Ry0
      real*8 p1_sum, wt, p1
      real*8 BR_haz(MAX_INTEN, MAX_PROB,MAX_BRANCH,MAX_NODE)
      integer BR_index(MAX_FLT,20,MAX_WIDTH,MAXPARAM), nNode(MAX_NODE)
      integer segModelFlag(MAX_FLT,100), nSegModel(MAX_FLT)
      integer icellRupStrike, icellRupDip, runflag
      real BR_wt(MAX_FLT,20,MAX_WIDTH,MAXPARAM), br_wt1(MAX_BRANCH,MAX_NODE)
      real segModelWt1(MAX_FLT,100), distDensity2(MAX_GRID), lnDir, lgIo, lgIntenscl
      real sigDirY, sigtemp, lg1, sig1, lat1, wt1
      integer n1AS(MAXFLT_AS), n2AS(MAXFLT_AS)
      real phi, tau, medadj, sigadj, phiSSS
      character*80 filebmode, siteFile, outfile, out3file
      integer bnum, bnumflag, coefcountRrup, coefcountRjb
      integer iMixture(4, MAX_PROB, MAX_ATTEN)
      real pLocY(MAXFLT_AS), sigmaTotal, sigma1, sigma2
      real*8 prock1, prock2
      real*8 sum0_Mo(MAXPARAM), sum1_Mo(MAXPARAM)
      real Pmag_all(MAXPARAM)
      integer rup1_flag
      
      real RefMag(MAX_FLT,MAXPARAM,MAX_WIDTH)
      real SA_med_SSC(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_PROB,MAX_FTYPE)
      real SA_sigma_SSC(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_PROB,MAX_FTYPE)
      real SA_med_GMC(MAX_FLT,4,MAX_ATTEN,MAX_PROB)
      real SA_sigma_GMC(MAX_FLT,4,MAX_ATTEN,MAX_PROB)
      real SA_SSC(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_PROB,MAX_FTYPE)
      real SA_GMC(MAX_FLT,4,MAX_ATTEN,MAX_PROB)
      real gm1, gm_max_SSC(MAX_FLT,MAX_PROB), gm_min_SSC(MAX_FLT,MAX_PROB)
      real gm_max_GMC(MAX_FLT,MAX_PROB), gm_min_GMC(MAX_FLT,MAX_PROB)
      real det_min_SSC(MAX_PROB), det_min_SSC1(MAX_PROB)
      real det_min_GMC(MAX_PROB), det_min_GMC1(MAX_PROB)
      real det_GMC(MAX_PROB,4,MAX_ATTEN), deter1(MAX_PROB)
      integer i100
      real wt100, det_mean(MAX_FLT,MAX_PROB)
      real deter_GMTYpe(4,MAX_PROB), UHS1(MAX_PROB), UHS2(MAX_PROB)
      integer iControlFlt(MAX_PROB)
      real SR0(4), SR1(4)
      data SR0 / 0., 0., 0.1, 100. /
      data SR1 / 0.1, 1.0, 9., 100. /
   
c     Write Program information to the screen.
      write (*,*) '*********************************'
      write (*,*) '*   Deterministic Hazard Code: Deter_45.2    *'
      write (*,*) '*         Mar, 2016          *'
      write (*,*) '*********************************'
      write (*,*)

      write (*,*) 'Enter number of cases to run in batch mode.'
      write (*,*) '         (For a single run enter 0)'
      read (*,*) bnum
      if (bnum .eq. 0) then
          bnumflag = 0 
          bnum = 1
      else
          bnumflag = 1
      endif

      if ( bnumflag .eq. 1) then
         write (*,*) 'Enter the batch mode filename.'
         read (*,'(a80)') filebmode
         open (77,file=filebmode,status='old')
      endif
      
c     Start loop over number of batch mode runs
      do 2000 ibnum=1,bnum
         if (bnumflag .eq. 1) then
            write (*,*) 'Looping over Number of Batch Mode Cases: ', ibnum, bnum
         endif

c     Read Run File
      call RdInput ( nProb, nAttenType, nAtten, jcalc, specT, sigTrunc,
     1               gmScale, dirFlag, nInten, testInten, lgTestInten, 
     2               psCorFlag, minlat, maxlat, minlong, maxlong, distmax,
     3               nMagBins, magBins, nDistBins, distBins, nepsBins, epsBins,
     4               nXcostBins, xcostBins, soilAmpFlag, gm_wt, runflag, sigvaradd,
     5               sCalc, sigfix, ssscalc, bnumflag, cfcoefRrup, cfcoefRjb, 
     6               coefcountRrup, coefcountRjb, iMixture )

c     read fault File
      call Rd_Fault_Data ( nFlt, fName, minMag, magStep, xStep,
     1     yStep, segModelWt, rateParam, rateParamWt, beta,
     2     magRecur, magRecurWt, faultWidth, faultWidthWt, 
     3     maxMag,  maxMagWt, fLong, fLat, fZ, dip, nfp, nMag, 
     4     ftype, sourceType, nRupArea, coeff_area, sigArea, nRupWidth, 
     5     coeff_width, sigWidth, nParamVar, iCoor,minDepth,
     6     fIndex, probAct, nWidth, mpdf_param, 
     7     al_segWt, attenType, sampleStep,
     8     grid_a,grid_dlong,grid_dlat,grid_n,grid_long, grid_lat,
     9     grid_top, minlat, maxlat, minlong, maxlong, scaleRate, fsys,
     1     mMagout, mMagoutWt, fltDirect, synchron, nsyn_Case, synjcalc,
     1     synmag, syndistRup, syndistJB, synDistSeismo, synHypo,
     2     synftype, synhwflag, synwt, RateType, iDepthModel, depthParam, 
     3     nMaxmag2, segWt1, faultFlag, nDD, nFtype, ftype_wt, 
     4     br_index, br_wt, segModelFlag, nSegModel, segModelWt1, runflag, 
     7     syn_dip, syn_zTOR, syn_RupWidth, syn_RX, syn_Ry0, refMag )

c     Read starting dam number
      read (13,*) iDam1
      
c     Open the file withe dam site info and class
      read (13,'( a80)') sitefile
      open (14,file=sitefile,status='old')

c     Open file containing the names of the out3 files
      read (13,'( a80)') out3file
      open (22,file=out3file,status='old')

c     Open output file for SSC Deter Tornado      
      read (13,'( a80)') outfile
      open (56,file=outfile,status='unknown')

c     Open output file for GMC Deter Tornado      
      read (13,'( a80)') outfile
      open (57,file=outfile,status='unknown')
                      
c     Open output file for the deterministic results      
      read (13,'( a80)') outfile
      open (59,file=outfile,status='unknown')

c     Open output file for the deterministic by Fault      
      read (13,'( a80)') outfile
      open (60,file=outfile,status='unknown')

c     Loop Over Number of Sites
      read (14,*) nSite   
      write (*,'( i5,2x,'' nsite'')') nSite
        
      do 1000 iSite = 1, nSite      
            
c      Read site coordinates and properties
       read (14,*) i1, SiteY, SiteX, iDamClass, vs, depthvs10, depthvs15, D25, vrup,
     1       forearc
      
C       Initialize temp GM arrays (all faults)
        call Init_SSC ( SA_med_SSC, SA_sigma_SSC, SA_SSC )
        call Init_GMC ( SA_med_GMC, SA_sigma_GMC, SA_GMC )

       do iPRob=1,nPRob     
         det_min_SSC(iPRob) = -1.E30
         det_min_GMC(iPRob) = -1.E30
       enddo
       
c      Read out3 file to get the hazard levels
       call RdOut3 ( UHS1, UHS2, sourceType ) 
       write (*,'( i5,2f10.4)') (k,UHS1(k), UHS2(k),k=1,nProb)

       if ( iSite .lt. iDam1 ) goto 1000

c      Loop Over Number of sources 
       do 900 iFlt = 1, nFlt


        i100 = 0
        wt100 = 0.

         do iProb=1,nProb
           det_mean(iFlt,iProb) = 0.
           gm_max_SSC(iFlt,iProb) = -1.E30
           gm_min_SSC(iFlt,iProb) = 1.E30
           gm_max_GMC(iFlt,iProb) = -1.E30
           gm_min_GMC(iFlt,iProb) = 1.E30
         enddo
       
        write (*,'(2x,''Site = '',i4,'', '',''iFlt ='',4i5)') iSite, iFlt, nFlt, sourceType(iFlt), ibnum
        wt_sum = 0.
          
c        skip areal zones
         if ( sourceType(iFlt) .eq. 2 .or. sourceType(iFlt) .eq. 3 ) then
           do iProb=1,nProb
             gm_min_SSC(iFlt,iProb) = -1.E30
             gm_max_SSC(iFlt,iProb) = -1.E30
             gm_min_GMC(iFlt,iProb) = -1.E30
             gm_max_GMC(iFlt,iProb) = -1.E30
             det_mean(iFlt,iProb) = -99.
           enddo
             goto 890        
         endif

c       Loop over alternative Fault Widths (epistemic)
c       (This changes the geometry of the source)
        do 860 iFltWidth=1,nWidth(iFlt)
        	
c        Set bottom of fault for standard faults (source type 1)
          if ( sourceType(iFlt) .eq. 1. ) then
            call SetFltBottom (iCoor, iFlt, nfp, dip(iFlt,iFltWidth,1), 
     1                         faultWidth(iFlt,iFltWidth), fZ, flat, flong, nDD)
          endif

c        Convert Long, Lat to x,y in km and put into new array (1-D)
         call ConvertCoordinates2 (nfp(iFlt), iFlt, iCoor, grid_n(iFlt), 
     1           sourceType(iFlt), nDD(iFlt), siteX, siteY, fLat, fLong, fZ, 
     2           grid_lat, grid_long, grid_dlat, grid_dlong, nPts, xFlt, yFlt, 
     3           zFlt, grid_x, grid_y, grid_dx, grid_dy, x0, y0, z0)      

c        Turn fault into a grid 
         if ( sourceType(iFlt) .eq. 1 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               minDist, xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )   
         elseif ( sourceType(iFlt) .eq. 5 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               minDist, xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )   
         elseif ( sourceType(iFlt) .eq. 6 ) then
           call calcFltGrid ( xFlt, yFlt, zFlt, nfp(iFlt), nDD(iFlt), fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, aveWidth, 
     3               minDist, xStep(iFlt), fltGrid_fLen, fltGrid_x1, fltGrid_y1, 
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3, 
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 ) 
         endif

c        Set Sampling of Rupture Area and Rupture Width Distributions
         call initRup ( sigArea, nRupArea, sigMaxArea, areaStep, iFlt)
         call initRup ( sigWidth, nRupWidth, sigMaxWidth, widthStep, iFlt)
                    
c         write (*,'( 2i5)') nRupArea(iFlt),nRupWidth(iFlt)
c         pause 'narea, nwid'  
          mag = 9.          
                    
c         Intergrate Over Rupture Area for this mag (aleatory)
          do 750 iArea = 1, nRupArea(iFlt)

c          Compute Rupture Area and Probability of Rupture Area
           call rupDimProb ( mag, coeff_area, sigArea, areaStep, 
     1             sigMaxArea, rupArea, pArea, iFlt, iArea )

c          Intergrate Over Rupture Width for this mag (aleatory)
           do 700 iWidth = 1, nRupWidth(iFlt)

c           Compute Rupture Width and Probability of Rupture Width
            call rupDimProb ( mag, coeff_width, sigWidth, 
     1         widthStep, sigMaxWidth, rupWidth, pWidth, iFlt, iWidth)

        call RupDims (sourcetype(iFlt), rupWidth, aveWidth, rupArea, faultLen,
     1                faultWidth(iFlt,iFltWidth), nLocYST1, yStep(iFlt), rupLen)

c      fix for deterministic
c      Just make a large dimension
       rupArea = 200000.
       rupWidth = 200.
       rupLen = 1000.

        call nLocXcells (sourceType(iFlt), nLocXAS, grid_n(iFlt), nfltgrid, fltgrid_w,
     1                   rupWidth, fltgrid_a, ruparea, nLocYST1, nLocX, n1AS, n2AS)

c           Integrate Over Rupture Location - along strike (aleatory)
c           This is along strike for faults and epicentral distance for source zones
            iDepthFlag = 0


            do 650 iLocX = 1, nLocX

               call nLocYcells (iLocX, n1AS, sourceType(iFlt), nLocX, distDensity, xStep(iFlt),
     1                          faultWidth(iFlt,iFltWidth), yStep(iFlt), distDensity2, grid_x,
     2                          grid_y, x0, y0, nLocY, pLocX, r_horiz)

c            Integrate Over Rupture Location - Down Dip (aleatory)
             do 600 iLocY = 1, nLocY

c             SourceType 1 fixed, assumes hypocenter is in middle of rupture
c             Set the hypocentral depth (is this really ztor??)
              if (sourceType(iFlt) .eq. 1 ) then
                hypoDepth = (iLocY-1.)*ystep(iFlt)*sin(abs(dip(iFlt,iWidth,1))*3.14159/180.0) + zFlt(1,1)
     1          + ((0.5*rupWidth)*sin(abs(dip(iFlt,iWidth,1))*3.14159/180.0))
              elseif (sourceType(iFlt) .eq. 5 ) then
                hypoDepth = fltgrid_Z(iLocY,iLocX)
              elseif ( sourceType(iFlt) .eq. 2 .or. sourceType(iFlt) .eq. 3 ) then
                hypoDepth = (iLocY-0.5)*ystep(iFlt) + grid_top(iFlt,1)
              elseif ( sourceType(iFlt) .eq. 4 ) then
                hypoDepth = (iLocY-0.5)*ystep(iFlt) + grid_top(iFlt,1)
              endif  

c            Find the Closest Distances for this rupture
c            Pass along fault grid locations for calculation of HW and Rx values within CalcDist subroutine.     
             call CalcDist (sourceType(iflt), pscorflag, hypoDepth, RupWidth, RupLen, 
     1             r_horiz, mindepth(iflt), nFltGrid, n1AS, iLocX, iLocY, n2AS,
     2             fltGrid_x, fltGrid_y, fltGrid_z, fltgrid_x1, fltgrid_y1, 
     3             fltgrid_z1, fltgrid_x2, fltgrid_y2, fltgrid_x3, fltgrid_y3,
     4             fltgrid_x4, fltgrid_y4, fltgrid_z4, fltGrid_Rrup, fltGrid_Rjb,
     5             distJB, distRup, ZTOR, distSeismo, distepi, disthypo, HWFlag,
     6             dipavg, n1, n2, Rx, Ry, Ry0, icellRupstrike, icellRupdip, dip, iFltWidth, iFlt) 

c             Loop over ftypes
              do 561 iFtype=1,nFtype(iFlt)

c             Loop Over Number of Problems (e.g. spectral periods)
              do 560 iProb=1,nProb
               jType = attenType(iFlt)
               do 550 iAtten = 1,nAtten(iProb,jType)

c               Check for negative jcalc values which will set the corresponding sigma to 
C               either a fixed value or sigma from another model.
                sigflag = 0
                if (jcalc(iProb,jType,iAtten) .lt. 0) then
                   jcalc1 = abs(jcalc(iProb,jType,iAtten) )
                   scalc1 = scalc(iProb,jtype,iAtten) 
                   sigfix1 = sigfix(iProb,jType,iAtten)
c                   ssscalc1 = ssscalc(iProb,jType,iAtten)
C                Check for either fixed sigma value (scalc1<0) or other sigma model
                   if (scalc1 .lt. 0) then
                      sigflag = 2
                   else
                      sigflag = 1
                   endif
                else
                   jcalc1 = jcalc(iProb,jType,iAtten) 
                endif

               dipaverage(1) = dipavg*180.0/3.14159  


c     temp fix
      nHypoX = 1
      nHypoZ = 1
                    
c              Loop over parameter variations (epistemic)
               do 500 iParam=1,nParamVar(iFlt,iFltWidth)
                      mag = refMag(iFlt, iParam, iFltWidth )

c                Compute the median and sigma of the ground motions
                 call meanInten ( distRup, distJB, distSeismo,
     1               HWFlag, mag, jcalc1, specT(iProb),  
     2               lgInten,sigmaY, ftype(iFlt,iFtype), attenName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipaverage(1),
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6               cfcoefrrup, cfcoefrjb, Ry0 )

c                Add epistemic uncertainty term (constant shift) to median
                 lgInten = lgInten + gmScale(iProb,jType,iAtten)

C                Second call goto GMPE for different sigma model 
                 if (sigflag .eq. 1) then
                  call meanInten ( distRup, distJB, distSeismo,
     1               hwflag, mag, scalc1, specT(iProb),  
     2               temp, sigmaY, ftype(iFlt,iFtype), sigmaName, period1, 
     3               iAtten, iProb, jType, vs, hypodepth, intflag, AR, dipaverage(1),
     4               disthypo, depthvs10, depthvs15, D25, tau,
     5               zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi, 
     6               cfcoefrrup, cfcoefrjb, Ry0 )

c                Check if a constant, user input sigma, is selected
                 elseif (sigflag .eq. 2) then
                  sigmaY = sigfix1
                 endif

C                Adjust sigma value if sigvaradd .ne. 0.0
                 if (sigvaradd(iProb,jType,iAtten) .ne. 0.0) then
                  sigmaY = sqrt(sigmaY*sigmaY + sigvaradd(iProb,jType,iAtten) )
                 endif

c                Check that sigma is not less than zero (0.0001)
                 if (sigmaY .lt. 0.0001 ) sigmaY = 0.0001

c                Reset SigmaTotal variable
                 sigmaTotal = sigmaY

c                Set the weight for this set of parameters (epistemic)
                     wt = RateParamWt(iFlt,iParam,iFltWidth) 
     1                 * magRecurWt(iFlt,iParam,iFltWidth) 
     2                 * faultWidthWt(iFlt,iFltWidth)
     2                 * maxMagWt(iFlt,iParam,iFltWidth) 
     2                 * ftype_wt(iFlt,iFtype) 

c      Is this supposed to be here?

                temp1 = sigmaTotal*wt

c               Set GM level based on dam class and slip rate
                if ( RateType(iFlt,iParam,iFltWidth) .eq. 1 ) then
                  SR = rateParam(iFlt,iParam,iFltWidth)
                  if ( SR .gt. SR1(iDamClass)) then
                    gm1 = lgInten + sigmaTotal
                  elseif ( SR .lt. SR0(iDamClass)) then
                    gm1 = lgInten
                  else
                    gm1 = alog(UHS1(iProb))
                    if ( gm1 .lt. lgInten ) gm1 = lgInten
                    if ( gm1 .gt. lgInten+sigmaTotal ) gm1 = lgInten+sigmaTotal
                  endif
                else
c                 This is just for cascadia for now
                  gm1 = lgInten + sigmaTotal
                endif


                i100 = i100 + 1
                if ( iProb .eq. 1. ) wt100 = wt100 + wt*gm_wt(iProb,jType,iAtten)
                det_mean(iFlt,iProb) = det_mean(iflt,iProb) + gm1*wt*gm_wt(iProb,jType,iAtten)
c                write (*,'( i5,2f8.1,3f10.4)') iFLt, mag, distRup, gm1, lgInten, sigmaTotal
                
                SA_SSC(iFlt,iFltWidth,iParam,iProb,iFtype) = 
     1                gm1 * gm_wt(iProb,jType,iAtten)
     2              + SA_SSC(iFlt,iFltWidth,iParam,iProb,iFtype)

                if ( SA_GMC(iFlt,jtype,iAtten,iProb) .eq. -999. ) SA_GMC(iFlt,jtype,iAtten,iProb) = 0.
                SA_GMC(iFlt,jtype,iAtten,iProb) = 
     1                gm1*wt + SA_GMC(iFlt,jtype,iAtten,iProb)
                                              
 500                continue

c 510               continue
c 530              continue
c 540             continue
c 549            continue
 550           continue
 560          continue
 561          continue
 600         continue
 650        continue
c 651        continue
 700       continue
 750      continue

c 800     continue
     

 860    continue

c       Set the minimum and maximum GM for this fault for SSC
        do iProb=1,nProb
         do iWidth=1,nWidth(iFlt)
          do iFtype=1,nFtype(iFlt)
           do iParam=1,nParamVar(iFlt,iWidth)
            if( SA_SSC(iFlt,iWidth,iParam,iProb,iFtype) .gt. gm_max_SSC(iFlt,iProb)) then
              gm_max_SSC(iFlt,iProb) = SA_SSC(iFlt,iWidth,iParam,iProb,iFtype)
            endif
            if( SA_SSC(iFlt,iWidth,iParam,iProb,iFtype) .lt. gm_min_SSC(iFlt,iProb)) then
              gm_min_SSC(iFlt,iProb) = SA_SSC(iFlt,iWidth,iParam,iProb,iFtype)
            endif
           enddo
          enddo
         enddo
        enddo
        goto 900

c       special coding for background
 890    do iProb=1,nProb
          det_mean(iFlt,iProb) = alog(UHS2(iProb))
        enddo

 900   continue
c       pause 'end flt loop'

c      Find the max of the min GMs for all faults
c      This is the smallest determinsitic value for SSC
       do iProb=1,nProb
         det_min_SSC(iProb) = -1.0E10
         do iFlt=1,nFlt
           if ( gm_min_SSC(iFlt,iProb) .gt. det_min_SSC(iProb) ) then
             det_min_SSC(iProb) = gm_min_SSC(iFlt,iProb)
c             write (*,'( 2i5,f10.4)') iProb, iFlt, det_min_SSC(iProb)
           endif
         enddo
         det_min_SSC1(iProb) = exp(det_min_SSC(iProb))
c         write (*,'( e12.4)') det_min_SSC(iPRob)
       enddo

c      Find the deterministic GM over all faults for each GM model
       do iProb=1,nProb
        do jType=1,nAttenType
         do iAtten=1,nAtten(iProb,jType)
          det_GMC(iProb,jType,iAtten) = -1.E30
          do iFlt=1,nFlt
           t1 = SA_GMC(iFlt,jtype,iAtten,iProb)
           if ( t1 .gt. det_GMC(iProb,jTYpe,iAtten) ) then
             det_GMC(iProb,jTYpe,iAtten) = t1
           endif
          enddo
         enddo
        enddo
       enddo   

       do iflt=1,nFlt
         write (60,'( 2i5,20f10.4)') iSite, iFlt, (det_mean(iFlt,iProb),iProb=1,nProb)
       enddo

c      Find the deterministic GM over all faults for mean determinsitic       
       do iProb=1,nProb
        deter1(iProb) = -1.E30
        do iFlt=1,nFlt
         if ( det_mean(iFlt,iProb) .gt. deter1(iProb) ) then
          deter1(iProb) = det_mean(iFlt,iProb)
          iControlFlt(iProb) = iFlt
         endif
        enddo
        write (59,'(2i5,2f10.4,i5,2x,a80 )') iSite, iProb, specT(iProb), exp(deter1(iProb)),
     1   iControlFlt(iProb), fname(iControlFlt(iProb))
       enddo 
c       pause 'deter1'

c      Find the mean deterministic by source type
       do iProb=1,nProb
        do jType=1,nAttenType
         deter_GMTYpe (jType,iProb) = 0.
         do iAtten=1,nAtten(iProb,jType)
           deter_GMTYpe(jType,iProb) = deter_GMTYpe(jType,iProb) 
     1          + det_GMC(iProb,jTYpe,iAtten)* gm_wt(iProb,jType,iAtten)
         enddo
        enddo 
       enddo

c      Write out SSC GM
       do iProb=1,nProb
        do iFlt=1,nFlt
         if ( gm_max_SSC(iflt,iProb) .gt. det_min_SSC(iProb) ) then

c         find the deterministic GM considering all other faults
          det_min1 = -1.E30
          do jFlt=1,nFlt
           if (jFlt .ne. iFlt) then
             if ( det_mean(jFlt,iProb) .ge. det_min1 ) then
               det_min1 = det_mean(jFlt,iProb)
             endif
           endif
          enddo

          do iWidth=1,nWidth(iFlt)
           do iFtype=1,nFtype(iFlt)
            do iParam=1,nParamVar(iFlt,iWidth)
             t1 = SA_SSC(iFlt,iWidth,iParam,iProb,iFtype)

c            check if this case would be controlling GM. 
c             If not, then use the deterministic GM from the other faults
c             if ( t1 .lt. det_min1) t1 = det_min1
             write (56,'( 6i5,4f10.4)') iSite, iProb, iFlt, iWidth, iFtype, iParam, 
     1         exp(t1), t1 - deter1(iProb), det_min1
            enddo
           enddo
          enddo
         endif
        enddo
       enddo

c      Write out GMC GM
       do iProb=1,nProb
        do jType=1,nAttenType
         if ( jTYpe .eq. 1 ) Jtype2 = 2
         if ( jTYpe .eq. 2 ) Jtype2 = 1
         
         do iAtten=1,nAtten(iProb,jType)
c         if ( det_GMC(iProb,jTYpe,iAtten) .gt. deter_GMType(jType2,iProb)) then
            t1 = det_GMC(iProb,jTYpe,iAtten)
c         else
c           t1 = deter_GMType(jType2,iProb)
c         endif
          write (57,'( 4i5,4f10.4)') iSite, iProb, jTYpe, iAtten, exp(t1),
     1      t1 - deter1(iProb), (exp(deter_GMType(k,iProb)),k=1,2)
         enddo
        enddo
      enddo
      flush (56)
      flush (57)
      flush (59)
      flush (60)

c      write the output
c       call SSC_out ( gm_max, det_min_all, SA_SSC, nFlt, nParamVar, nWidth, nFtype, nProb )
c       call SSC_out ( gm_max, det_min_all,SA_GMC, nFlt, nAtten, nattenType, nProb )
        
 1000 continue

 2000 continue
      close (77)

      stop
      end

c -------------

      subroutine RdOut3 ( UHS1, UHS2, sourceType ) 
c      implicit none
      include 'pfrisk.h'

      character*80 file1, dummy, fname(MAX_FLT)
      integer nRd, nFlt, nProb, nAmp(MAX_PROB)
      integer iProb, iAtten, iFlt, iAmp, k
      real lat, long
      real segModelwt(MAX_FLT),al_segwt(MAX_FLT),
     1     mindist(MAX_FLT),Haz(MAX_PROB,MAX_INTEN,MAX_FLT)
      real HazTotal(MAX_PROB,MAX_INTEN), amp(MAX_PROB,MAX_INTEN)
      real mBar(MAX_PROB,MAX_INTEN), dBar(MAX_PROB,MAX_INTEN),
     1  eBar(MAX_PROB,MAX_INTEN)
      real testHaz, x, UHS1(MAX_PROB), UHS2(MAX_PROB)
      real haz_zone(MAX_PROB,MAX_INTEN)
      integer sourceType(MAX_FLT)

c     open the out3 file
      read (22,'( a80)') file1  
      write (*,'( a80)') file1
c      pause
      nRd = 23
      open (nRd,file=file1,status='old')

c     Read the output 3 file, keeping the total hazard and the background 
c     zone hazard 

c     Read in the hazard curves for all periods in given Haz45 ouptut file.
      read (nRd,*) nFlt 
      read (nRd,'(a1)') dummy

      do j=1,4
        read (nRd,'( a1)') dummy
      enddo

      read (nRd,'( 21x,2f9.3)') long, lat
      read (nRd,*) nProb

      do iProb=1,nProb
          
        read (nRd,'( 15x,i5)')  iAtten
        read (nRd,*) nAmp(iProb)
        read (nRd,'( 60x,30f12.4)') (amp(iProb,k),k=1,nAmp(iProb))

        do iFlt=1,nFlt
          read (nRd,'( 2x,a38,2f6.3,f8.1,1x,30e12.4)') fname(iFlt),
     1              segModelwt(iFlt),al_segwt(iFlt),
     1              mindist(iFlt),(Haz(iProb,k,iFlt),k=1,nAmp(iProb))
        enddo

        read (nRd,'( 61x,50e12.4)') (HazTotal(iProb,k),k=1,nAmp(iProb))
        read (nRd,'(a1)') dummy

        read (nRd,'( 61x,50e12.3)') (mBar(iProb,k),k=1,nAmp(iProb))
        read (nRd,'( 61x,50e12.3)') (dBar(iProb,k),k=1,nAmp(iProb))
        read (nRd,'( 61x,50e12.3)') (eBar(iProb,k),k=1,nAmp(iProb))

        do idum=1,3
          read (nRd,'( a80)') dummy
        enddo
            
      enddo
      close (10)

c     Interpolate to desired return period for the total hazard
      testHaz = 1./2000.
      do iProb=1,nProb
	do iAmp=2,nAmp(iProb)

c         Check for zero values in hazard curve.
          if (hazTotal(iProb,iAmp) .eq. 0. ) then
            write (*,*) 'warning: Zero Values for hazard curve at desired haz level.'
            write (*,*) 'Setting UHS to last nonzero value'
            UHS1(iProb) = exp(x)
          endif
          
c         Interpolate the hazard curve.
          if ( hazTotal(iProb,iAmp) .lt. testHaz ) then
            x = ( alog(testHaz) - alog(hazTotal(iProb,iAmp-1)) )/
     1            ( alog(hazTotal(iProb,iAmp))-alog(hazTotal(iProb,iAmp-1))) 
     2          * (alog(amp(iProb,iAmp))-alog(amp(iProb,iAmp-1))) 
     3          + alog(amp(iProb,iAmp-1))
            UHS1(iProb) = exp(x)
            goto 10
          endif
        enddo
 10     continue
      enddo

c     Compute the hazard from background zones
      do iProb=1,nProb
	do iAmp=2,nAmp(iProb)
          haz_zone(iProb,iAmp) = 0.
        enddo
      enddo

      do iFlt=1,nFlt
        if ( sourceType(iFlt) .eq. 2 .or.  sourceType(iFlt) .eq. 3 ) then
          do iProb=1,nProb
  	    do iAmp=1,nAmp(iProb)
              haz_zone(iProb,iAmp)= haz_zone(iProb,iAmp) + Haz(iProb,iAmp,iFlt) 
            enddo
          enddo
        endif
      enddo

c      do iProb=1,nProb
c        write (*,'( i5,15e12.4 )') iProb, (haz_zone(iProb,iAmp), iAmp=1,nAmp(iProb))
c      enddo
      
c     Interpolate background to desired return period for the background zones
      testHaz = 1./3000.
      do iProb=1,nProb
	do iAmp=2,nAmp(iProb)

c         Check for zero values in hazard curve.
          if (haz_zone(iProb,iAmp) .eq. 0. ) then
            write (*,*) 'warning: Zero Values for hazard curve at desired haz level.'
            write (*,*) 'Setting UHS to last nonzero value'
            UHS2(iProb) = haz_zone(iProb,iAmp-1)
          endif
          
c           Interpolate the hazard curve.
          if ( haz_zone(iProb,iAmp) .lt. testHaz ) then
            x = ( alog(testHaz) - alog(haz_zone(iProb,iAmp-1)) )/
     1          ( alog(haz_zone(iProb,iAmp))-alog(haz_zone(iProb,iAmp-1))) 
     2        * (alog(amp(iProb,iAmp))-alog(amp(iProb,iAmp-1))) 
     3        + alog(amp(iProb,iAmp-1))
            UHS2(iProb) = exp(x)
            goto 20
          endif
        enddo
 20     continue
      enddo

      return
      end

    

      
      
      
