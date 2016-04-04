c fHaz_init.f                                                                  
                                                                                               
      subroutine InitFltMax ( amxmg, mindist, amxi )                            
      include 'pfrisk.h'                                                        
      real amxmg, mindist, amxi                                                 
                                                                                
      AMXMG = -99999.                                                           
                                                                                
      AMXI = -99999.                                                            
      minDist = 1.0e30                                                          
      return                                                                    
      end                                                                       
                                                                                
c -------------------------------------------------------------

      subroutine Init_SSC ( SA_med_SSC, SA_sigma_SSC, SA_SSC )
      include 'pfrisk.h'                                                        
      real SA_med_SSC(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_PROB,MAX_FTYPE)
      real SA_sigma_SSC(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_PROB,MAX_FTYPE)
      real SA_SSC(MAX_FLT,MAX_WIDTH,MAXPARAM,MAX_PROB,MAX_FTYPE)
      
      do iFLt=1,MAX_FLT
      do iParam=1,MAXPARAM
        do iFltWidth=1,MAX_WIDTH
            do iProb=1,MAX_PROB
              do iFtype=1,MAX_FTYPE
                SA_med_SSC(iFlt,iFltWidth,iParam,iProb,iFtype) = 0.0
                SA_sigma_SSC(iFlt,iFltWidth,iParam,iProb,iFtype) = 0.0
                SA_SSC(iFlt,iFltWidth,iParam,iProb,iFtype) = 0.0
              enddo
            enddo
        enddo
      enddo
      enddo
      return
      end                                                            

c -------------------------------------------------------------

      subroutine Init_GMC ( SA_med_GMC, SA_sigma_GMC, SA_GMC )
      include 'pfrisk.h'                                                        
      real SA_GMC(MAX_FLT,4,MAX_ATTEN,MAX_PROB)
      real SA_med_GMC(MAX_FLT,4,MAX_ATTEN,MAX_PROB)
      real SA_sigma_GMC(MAX_FLT,4,MAX_ATTEN,MAX_PROB)
      
      do iFLt=1,MAX_FLT
      do i=1,4  
         do iProb=1,MAX_PROB
           do iAtten=1,MAX_ATTEN
             SA_GMC(iFlt,i,iAtten,iProb) = -999.
             SA_med_GMC(iFlt,i,iAtten,iProb) = 0.0
             SA_sigma_GMC(iFlt,i,iAtten,iProb) = 0.0
           enddo
         enddo
      enddo
      enddo
      return
      end                                                            

                                                                                
c --------------------------------------------------------------                
                                                                                
      subroutine initRup ( sigma, nValue, sigmaMax, step, iFlt)                 
      real sigma(1), step                                                       
      integer nValue(1)                                                         
                                                                                
c     For Small Sigma, only one rupture length is needed                        
      IF ( sigma(iFlt) .LE. 0.02 ) nValue(iFlt) = 1                             
                                                                                
c     Set Rupture Length Integration parameters                                 
c     (maximum number of standard deviations and step size)                     
      sigmaMax = 2.0                                                            
      step = 2.*sigmaMax / nValue(iFlt)                                         
                                                                                
      return                                                                    
      end                                                                       
                                                                                
c --------------------------------------------------------------                
                                                                                
      subroutine SetBin ( nBins, bins, x, iBin )                                
                                                                                
      integer nBins, iBin                                                       
      real bins(1), x                                                           
                                                                                
      ibin = 0                                                                  
                                                                                
      do i=2,nBins                                                              
        if ( x .le. bins(i) .and. x .ge. bins(i-1)) then                               
          iBin = i-1                                                                   
          return                                                                       
        endif                                                                          
      enddo                                                                     

c     Reset Bin Range to include testing value.                                                                                 
c     First check for value less than lowest bin value.
      if (ibin .eq. 0 .and. x .lt. bins(1) ) then
         bins(1) = x
         iBin = 1
      endif

c     Now check for value greater than highest bin value.
      if (ibin .eq. 0 .and. x .gt. bins(nBins) ) then
         bins(nBins) = x
         ibin = nBins - 1
      endif

      return                                                                                
      end                                                                       
                                                                                
c --------------------------------------------------------------                
                                                                                
      subroutine SetBin_Eps ( nBins, bins, x, iBin )                            
                                                                                
      integer nBins, iBin                                                       
      real bins(1), x                                                                                                                                           
 
      ibin = 0                                                                  
                                                                         
      do i=2,nBins                                                              
        if ( x .le. bins(i) .and. x .ge. bins(i-1)) then                               
          iBin = i-1                                                                   
          return                                                                       
        endif                                                                          
      enddo                                                                     

c     Reset Bin Range to include testing value.                                                                                 
c     First check for value less than lowest bin value.
      if (ibin .eq. 0 .and. x .lt. bins(1) ) then
         bins(1) = x
         iBin =1
      endif

c     Now check for value greater than highest bin value.
      if (ibin .eq. 0 .and. x .gt. bins(nBins) ) then
         bins(nBins) = x
         iBin = nBins - 1
      endif

      return                                                                                
      end           

c -------------------------------------------------------------

      subroutine Init_tempHaz ( tempHaz )
      include 'pfrisk.h'                                                        
      real*8 tempHaz(MAXPARAM,MAX_WIDTH,MAX_INTEN, MAX_PROB, MAX_ATTEN,MAX_FTYPE)
      
      do iParam=1,MAXPARAM
        do iFltWidth=1,MAX_WIDTH
          do iInten=1,MAX_INTEN
            do iProb=1,MAX_PROB
              do iAtten=1,MAX_ATTEN
                do iFtype=1,MAX_FTYPE
                  tempHaz(iParam,iFltWidth,iInten,iProb,iAtten,iFtype) = 0.0
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      return
      end                                                            

c -------------------------------------------------------------

      subroutine Init_tempHaz1 ( tempHaz1 )
      include 'pfrisk.h'                                                        
      real*8 tempHaz1(MAXPARAM,MAX_WIDTH,MAX_INTEN, MAX_PROB,MAX_FTYPE)
      
      do iParam=1,MAXPARAM
        do iFltWidth=1,MAX_WIDTH
          do iInten=1,MAX_INTEN
            do iProb=1,MAX_PROB
                do iFtype=1,MAX_FTYPE
                  tempHaz1(iParam,iFltWidth,iInten,iProb,iFtype) = 0.0
                enddo
            enddo
          enddo
        enddo
      enddo
      return
      end                                                            

c -------------------------------------------------------------

      subroutine Init_tempHaz2 ( tempHaz2 )
      include 'pfrisk.h'                                                        
      real*8 tempHaz2(4, MAX_ATTEN, MAX_INTEN, MAX_PROB )
      
      do i=1,4  
       do iInten=1,MAX_INTEN
         do iProb=1,MAX_PROB
           do iAtten=1,MAX_ATTEN
             tempHaz2(i, iInten,iProb,iAtten) = 0.0
           enddo
         enddo
       enddo
      enddo
      return
      end                                                            

c -------------------------------------------------------------

      subroutine InitDeagg ( m_bar, d_bar, e_bar, Xcost_bar )
      include 'pfrisk.h'                                                        
      real*8 m_bar(MAX_PROB,MAX_INTEN), d_bar(MAX_PROB,MAX_INTEN)
      real*8 Xcost_bar(MAX_PROB,MAX_INTEN)
      real*8 e_bar(MAX_PROB,MAX_INTEN)

      do i=1,MAX_PROB
         do ii=1,MAX_INTEN
           m_bar(i,ii) = 0.
           d_bar(i,ii) = 0.
           e_bar(i,ii) = 0.
           Xcost_bar(i,ii) = 0.
         enddo
      enddo
      return
      end

c -------------------------------------------------------------

      subroutine InitBrhaz ( nProb, nInten, nNode, br_haz )
      include 'pfrisk.h'                                                        
      integer nInten(MAX_PROB), nProb, iBr, iNode, jInten, iProb 
      integer nNode(MAX_NODE)
      real*8 BR_haz(MAX_INTEN, MAX_PROB, MAX_BRANCH, MAX_NODE)

      do iProb=1,nProb
       do jInten=1,nInten(iProb)
        do iBr=1,MAX_BRANCH 
         nNode(iBR) = 0
         do iNode=1,MAX_Node
          br_haz(jInten,iProb,iBr,iNode) = 0.
         enddo
        enddo
       enddo
      enddo       

      return 
      end        

c -------------------------------------------------------------

      subroutine InitMinDist ( nWidth, FaultDist )
      include 'pfrisk.h'                                                        
      integer nInten(MAX_PROB), iWidth
      real faultDist(MAX_FLT,MAX_WIDTH,3), minDist

      minDist = 1.e10
      do iWidth=1,nWidth(iFlt)
        do j=1,3
           FaultDist(iFlt,iWidth,j) = 1.e10
        enddo
      enddo

      return 
      end        

c -------------------------------------------------------------


