!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also contains the following subrouines:
!          abq_UEL_3D_integrationpoints           - defines integration ponits for 3D continuum elements
!          abq_UEL_3D_shapefunctions              - defines shape functions for 3D continuum elements
!          abq_UEL_invert3D                       - computes the inverse and determinant of a 3x3 matrix
!          abq_facenodes_3D                       - returns list of nodes on the face of a 3D element
!
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,k,l,n_points,kint,ie
      integer      :: face_node_list(8)                       ! List of nodes on an element face
    !
      double precision  ::  xi(3,64)                          ! Volumetric Integration points
      double precision  ::  w(64)                             ! Integration weights
      double precision  ::  N(20)                             ! 3D Shape functions
      double precision  ::  dNdxi(20,3)                       ! 3D Shape function derivatives
      double precision  ::  dxdxi(3,3)                        ! Derivative of position wrt normalized coords
      double precision  ::  dNdx(20,3)                        ! Derivative of shape functions wrt reference coords
      double precision  ::  dNdy(20,3)                        ! Derivative of shape functions wrt deformed coords
    !
    !   Variables below are for computing integrals over element faces
      double precision  ::  xi2(2,9)                          ! Area integration points
      double precision  ::  N2(9)                             ! 2D shape functions
      double precision  ::  dNdxi2(9,2)                       ! 2D shape function derivatives
      double precision  ::  dxdxi2(3,2)                       ! Derivative of spatial coord wrt normalized areal coord
    !
      double precision  ::  sigma(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  S(3,3), Snew(3,3)
      double precision  ::  stress(6)
      double precision  ::  F(3,3)                            ! Deformation gradient
      double precision  ::  Finv(3,3)                         ! Inverse of deformation gradient
      double precision  ::  C(3,3)                            ! C-G deformation tensor
      double precision  ::  JJ                                ! det(F)
      double precision  ::  G(6,6)
      double precision  ::  D(6,6)                            ! Material tangent
      double precision  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
      double precision  ::  q(9) 
      double precision  ::  H(6,9)
      double precision  ::  Y(3*NNODE,3*NNODE) 
      double precision  ::  Stif(NNODE,NNODE)
    !
    !     Example ABAQUS UEL implementing 3D linear elastic elements
    !     El props are:

    !     PROPS(1)         Young's modulus
    !     PROPS(2)         Poisson's ratio

      if (NNODE == 4) n_points = 1               ! Linear tet
      if (NNODE == 10) n_points = 4              ! Quadratic tet
      if (NNODE == 8) n_points = 8               ! Linear Hex
      if (NNODE == 20) n_points = 27             ! Quadratic hex

      call abq_UEL_3D_integrationpoints(n_points, NNODE, xi, w)

      if (MLVARX<3*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 3*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0


      ENERGY(1:8) = 0.d0
   
! ============================ B* matrix ================================
      Bstar = 0.d0
      
      call abq_UEL_3D_shapefunctions(xi,NNODE,N,dNdxi)
      dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
      call abq_UEL_invert3d(dxdxi,dxidx,determinant)
      dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)
      
      Bstar(1,1:3*NNODE-2:3) = dNdx(1:NNODE,1)
      Bstar(2,2:3*NNODE-1:3) = dNdx(1:NNODE,2)
      Bstar(3,3:3*NNODE:3)   = dNdx(1:NNODE,3)
      Bstar(4,1:3*NNODE-2:3) = dNdx(1:NNODE,2)
      Bstar(5,2:3*NNODE-1:3) = dNdx(1:NNODE,1)
      Bstar(6,1:3*NNODE-2:3) = dNdx(1:NNODE,3)
      Bstar(7,3:3*NNODE:3)   = dNdx(1:NNODE,1)
      Bstar(8,2:3*NNODE-1:3) = dNdx(1:NNODE,3)
      Bstar(9,3:3*NNODE:3)   = dNdx(1:NNODE,2)
       
! ============================ q vector ================================
      q = 0.d0
       
      ! Loop over integration points
      do kint = 1, n_points
        call abq_UEL_3D_shapefunctions(xi(1:3,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:3,1:NNODE),dNdxi(1:NNODE,1:3))
        call abq_UEL_invert3d(dxdxi,dxidx,determinant)
        dNdx(1:NNODE,1:3) = matmul(dNdxi(1:NNODE,1:3),dxidx)

        ! Caculate the deformation gradient
        do k = 1,3
            ie = 3*(NNODE-1)+k
            F(k,1:3) = matmul(U(k:ie:3),dNdx(1:NNODE,1:3))
            F(k,k) = F(k,k) + 1.d0
        end do
        
        ! right Cauchy-Green deformation tensor C 
        C = matmul(transpose(F),F)
        
        ! JJ = det(F) 
        call abq_UEL_invert3d(F,Finv,JJ)
        
        ! the second Piola-Kirchhoff stress, sigma
        call secondPK(PROPS(1:NPROPS),NPROPS,F,JJ,sigma,D)
        S = 0.d0
        
        S(1,1) = sigma(1)
        S(1,2) = sigma(4)
        S(1,3) = sigma(5)
        S(2,1) = sigma(4)
        S(2,2) = sigma(2)
        S(2,3) = sigma(6)
        S(3,1) = sigma(5)
        S(3,2) = sigma(6)
        S(3,3) = sigma(3)
        
        ! q vector 
        do k=1,3
        	q(1) = q(1) + S(1,k)*F(1,k)
        	q(2) = q(2) + S(2,k)*F(2,k)
        	q(3) = q(2) + S(3,k)*F(3,k)
        	q(4) = q(2) + S(2,k)*F(1,k)
        	q(5) = q(2) + S(1,k)*F(2,k)
        	q(6) = q(2) + S(3,k)*F(1,k)
        	q(7) = q(2) + S(1,k)*F(3,k)
        	q(8) = q(2) + S(3,k)*F(2,k)
        	q(9) = q(2) + S(2,k)*F(3,k)
        end do 
        	
! ============================ H matrix ================================        
        H = 0.d0
        
        H(1,1) = F(1,1)
        H(1,5) = F(2,1)
        H(1,7) = F(3,1)
        H(2,2) = F(2,2)
        H(2,4) = F(1,2)
        H(2,9) = F(3,2)
        H(3,3) = F(3,3)
        H(3,6) = F(1,3)
        H(3,8) = F(2,3)
        H(4,1) = F(1,2)
        H(4,2) = F(2,1)
        H(4,4) = F(1,1)
        H(4,5) = F(2,2) 
        H(4,7) = F(3,2)
        H(4,9) = F(3,1)
        H(5,1) = F(1,3)
        H(5,3) = F(3,1)
        H(5,5) = F(2,3)
        H(5,6) = F(1,1)
        H(5,7) = F(3,3)
        H(5,8) = F(2,1)
        H(6,2) = F(2,3)
        H(6,3) = F(3,2)
        H(6,4) = F(1,3)
        H(6,6) = F(1,2)
        H(6,8) = F(2,2)
        H(6,9) = F(3,3)
        
 ! ============================ Y matrix ================================     
        ! geometric stiffness
        Stif = 0.d0
        do i = 1,NNODE
        	do j = 1,NNODE
        		do k = 1,3
        			do l = 1,3
        				Stif(i,j) = Stif(i,j) + dNdx(i,k)*S(k,l)*dNdx(j,l)
                     end do
                 end do
             end do
         end do
        
         ! Y matrix  
         Y = 0.d0
         do i=1,NNODE
             do j=1,NNODE
                 Y(3*i-2,3*j-2)=Stif(i,j)
                 Y(3*i-1,3*j-1)=Stif(i,j)
                 Y(3*i,3*j)=Stif(i,j)
             end do
         end do
        
 ! ============================ R and K ================================              
     
        	RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
	   1         - matmul(transpose(Bstar(1:9,1:3*NNODE)),q(1:9))*
  	   2                                          w(kint)*determinant

  
       		 AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
    	 1          + matmul(transpose(Bstar(1:9,1:3*NNODE)),
   	  2           matmul(H,matmul(D,matmul(H,Bstar(1:9,1:3*NNODE))))*
   	  3            w(kint)*determinant +
      4            Y(1:3*NNODE,1:3*NNODE)*w(kint)*determinant
                  
 ! ============================ Stress ================================
       		 Snew = 0.d0
        	 Snew = matmul(F,matmul(stress,transpose(F)))
        	
             stress = 0.d0
             stress(1) = Snew(1,1)
             stress(2) = Snew(2,2)
             stress(3) = Snew(3,3)
             stress(4) = Snew(1,2)
             stress(5) = Snew(1,3)
             stress(6) = Snew(2,3)
        
           if (NSVARS>=n_points*6) then   ! Store Cauchy stress at each integration point (if space was allocated to do so)
               SVARS(6*kint-5:6*kint) =  stress(1:6)/JJ
           endif
        end do
      
      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    
      return

      END SUBROUTINE UEL


      subroutine secondPK(element_properties,n_properties,F,JJ,sigma,D)

       implicit none

       integer, intent(in)           :: n_properties
       double precision, intent(in)  :: element_properties(n_properties)
       double precision, intent(in)  :: F(3,3)
       double precision, intent(in)  :: JJ
       double precision, intent(out) :: sigma(6)
       double precision, intent(out) :: D(6,6)
       
       double precision :: mu
       double precision :: K
       double precision :: ss
       double precision :: G11, G22, G33, G44
       double precision :: G(6,6)
       double precision :: C(3,3)
       double precision :: Cvec(6)
       double precision :: Cbar(6)
       double precision :: Cinv(3,3)
       double precision :: Cinvvec(6)
       double precision :: Cstar(6)
       double precision :: Cstarbar(6)
       double precision :: eyevec(6)
       double precision :: Q(1,1)
       double precision :: P(6,1)
       double precision :: omega(6,6)
       integer :: i,j,k
       
       mu = element_properties(1)
       K  = element_properties(2)
       G11 = element_properties(3)
       G22 = element_properties(4)
       G33 = element_properties(5)
       G44 = element_properties(6)
       
! ====================== G ==========================
       G = 0.d0
       G(1,1) = G11
       G(2,2) = G22
       G(3,3) = G33
       do k = 4,6
       	G(k,k) = G44
       end do
       
! ===================== Cvec ==========================
       C = matmul(F,transpose(F))
       Cvec(1) = C(1,1)
       Cvec(2) = C(2,2)
       Cvec(3) = C(3,3)
       Cvec(4) = C(1,2)
       Cvec(5) = C(1,3)
       Cvec(6) = C(2,3)

! ===================== Cbar ==========================       
       call abq_UEL_invert3d(C,Cinv,ss)      ! ss is just a dummy variable here
       ss = JJ**(-2.d0/3.d0)
       
       Cbar = 0.d0
       do i = 1,6
       	Cbar(i) = C(i)*ss
       end do

! ==================== Cinvvec =========================       
       Cinvvec(1) = Cinv(1,1)
       Cinvvec(2) = Cinv(2,2)
       Cinvvec(3) = Cinv(3,3)
       Cinvvec(4) = Cinv(1,2)
       Cinvvec(5) = Cinv(1,3)
       Cinvvec(6) = Cinv(2,3)

! ================= Cstar, Cstarbar ====================    
       Cstar = 0.d0
       Cstar(1) = C(1,1)
       Cstar(2) = C(2,2)
       Cstar(3) = C(3,3)
       Cstar(4) = 2.d0*C(1,2)
       Cstar(5) = 2.d0*C(1,3)
       Cstar(6) = 2.d0*C(2,3)
       
       Cstarbar = 0.d0
       do i = 1:6
       	Cstarbar(i) = Cstar(i)*ss
       end do
       	
! ================================ Q ===================================
       eyevec(1:3) = 1.d0
       eyevec(4:6) = 0.d0
       Q = 0.d0
       Q = matmul(transpose(Cstarbar-eyevec),matmul(G,(Cstarbar-eyevec)))/4.d0     

! ================================ P ===================================
       P = 0.d0
       P = (matmul(G,(Cstarbar-eyevec))
  1      -(matmul(transpose(Cstar),matmul(G,Cstarbar-eyevec))*Cinvvec/3.d0))
  2                  *ss/2.d0
 
! ============================== Sigma =================================       		       
       sigma = 0.d0
       
       do i = 1,6
       	sigma(i) = mu*exp(Q)*P(i)-K*JJ**(JJ-1.d0)*Cinvvec(i)
       end do
       
! ============================== Omega =================================       
       omega = 0.d0
       
       do j = 1:3
       	omega(j,j) = Cinvvec(j)*Cinvvec(j)
       end do 
       omega(4,4) = (Cinvvec(1)*Cinvvec(2)+Cinvvec(4)*Cinvvec(4))/2.d0
       omega(5,5) = (Cinvvec(1)*Cinvvec(3)+Cinvvec(5)*Cinvvec(5))/2.d0
       omega(6,6) = (Cinvvec(2)*Cinvvec(3)+Cinvvec(6)*Cinvvec(6))/2.d0
    
       omega(1,2) = Cinvvec(4)*Cinvvec(4)
       omega(2,1) = Cinvvec(4)*Cinvvec(4)
       omega(1,3) = Cinvvec(5)*Cinvvec(5)
       omega(3,1) = Cinvvec(5)*Cinvvec(5)
       omega(1,4) = Cinvvec(1)*Cinvvec(4)
       omega(4,1) = Cinvvec(1)*Cinvvec(4)
       omega(1,5) = Cinvvec(1)*Cinvvec(5)
       omega(5,1) = Cinvvec(1)*Cinvvec(5)
       omega(1,6) = Cinvvec(4)*Cinvvec(5)
       omega(6,1) = Cinvvec(4)*Cinvvec(5)
       omega(2,3) = Cinvvec(6)*Cinvvec(6)
       omega(3,2) = Cinvvec(6)*Cinvvec(6)
       omega(2,4) = Cinvvec(4)*Cinvvec(2)
       omega(4,2) = Cinvvec(4)*Cinvvec(2)
       omega(2,5) = Cinvvec(4)*Cinvvec(5)
       omega(5,2) = Cinvvec(4)*Cinvvec(5)
       omega(2,6) = Cinvvec(2)*Cinvvec(6)
       omega(6,2) = Cinvvec(2)*Cinvvec(6)
       omega(3,4) = Cinvvec(5)*Cinvvec(6)
       omega(4,3) = Cinvvec(5)*Cinvvec(6)
       omega(3,5) = Cinvvec(5)*Cinvvec(3)
       omega(5,3) = Cinvvec(5)*Cinvvec(3)
       omega(3,6) = Cinvvec(6)*Cinvvec(3)
       omega(6,3) = Cinvvec(6)*Cinvvec(3)
       omega(4,5) = (Cinvvec(1)*Cinvvec(6)+Cinvvec(5)*Cinvvec(4))/2.d0
       omega(5,4) = (Cinvvec(1)*Cinvvec(6)+Cinvvec(5)*Cinvvec(4))/2.d0
       omega(4,6) = (Cinvvec(4)*Cinvvec(6)+Cinvvec(5)*Cinvvec(2))/2.d0
       omega(6,4) = (Cinvvec(4)*Cinvvec(6)+Cinvvec(5)*Cinvvec(2))/2.d0
       omega(5,6) = (Cinvvec(4)*Cinvvec(3)+Cinvvec(5)*Cinvvec(6))/2.d0
       omega(6,5) = (Cinvvec(4)*Cinvvec(3)+Cinvvec(5)*Cinvvec(6))/2.d0
       
! ================================ D ====================================  
        D = 0.d0
        
        D = mu*exp(Q)*(ss**2)*(G-(1.d0/3.d0)*(matmul(G,
     1  (spread(Cstar,dim=2,ncopies=6)*spread(Cinvvec,dim=1,ncopies=6)))
     2  + spread(Cinvvec,dim=2,ncopies=6)*spread(matmul(G,Cstar),dim=1,ncopies=6))
     3  -(1.d0/(3.d0*ss))*dot_product(Cstar,matmul(G,Cstarbar-eyevec))*omega
     4  +(1.d0/9.d0)*dot_product(Cstar,matmul(G,Cstar))*spread(Cinvvec,dim=2,ncopies=6)
     5  *spread(Cinvvec,dim=1,ncopies=6))+ mu*exp(Q)*(2.d0*spread(P,dim=2,ncopies=6)
     6  *spread((P-(1.d0/3.d0)*Cinvvec),dim=1,ncopies=6)
     7  -(1.d0*ss/3.d0)*spread(Cinvvec,dim=2,ncopies=6)*
     8  spread((matmul(G,Cstarbar-eyevec)),dim=1,ncopies=6))+K*JJ*(
     9  (2.d0*JJ-1.d0)*spread(Cinvvec,dim=2,ncopies=6)*
     10  spread(Cinvvec,dim=1,ncopies=6)+2.d0*(JJ-1.d0)*omega)
       
       return

      end subroutine neohooke

