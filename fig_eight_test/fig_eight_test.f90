PROGRAM fig_eight
  IMPLICIT NONE
  INTEGER :: n, i, j, k, ndims
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IDs
  DOUBLE PRECISION :: t, G, pi, dt, twrite, tcount, tlim, E0, E, Eerr, Ek, Ep, vmagsqrd, dist, vdotr, eta, amagsqrd
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, dr, dv, dtOpt
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r, v, a, jrk, s, c, apred, jrkpred

  OPEN(7, file='xpositions.txt', status='replace')
  OPEN(8, file='ypositions.txt', status='replace')
  OPEN(9, file='energy_error.txt', status='replace')

  n = 3
  ndims = 3

  ALLOCATE(IDs(1:n))
  ALLOCATE(m(1:n))
  ALLOCATE(r(1:ndims, 1:n))
  ALLOCATE(v(1:ndims, 1:n))
  ALLOCATE(a(1:ndims, 1:n))
  ALLOCATE(jrk(1:ndims, 1:n))
  ALLOCATE(s(1:ndims, 1:n))
  ALLOCATE(c(1:ndims, 1:n))
  ALLOCATE(dr(1:ndims))
  ALLOCATE(dv(1:ndims))
  ALLOCATE(apred(1:ndims, 1:n))
  ALLOCATE(jrkpred(1:ndims, 1:n))
  ALLOCATE(dtOpt(1:n))

  t = 0.
  tcount = 0.
  twrite = 0.01
  tlim = 100.
  dt = 0. ! Timestep length will be made non-zero before the first time the end-of-timestep positions/velocities are predicted (i.e. before snap/crackle are calculated)
  dtOpt = 0.
  eta = 1e-10
  G = 1.
  pi = 4. * ATAN(1.)
  r = 0.
  v = 0.
  a = 0.
  jrk = 0. ! 'jrk' rather than 'j' to distinguish it from the iteration variable 'j'
  s = 0.
  c = 0.
  IDs = 0.
  m = 1.
  E0 = 0.
  E = 0.
  Eerr = 0.
  Ek = 0.
  Ep = 0.
  vmagsqrd = 0.
  amagsqrd = 0.
  dist = 0.
  dr = 0.
  dv = 0.
  vdotr = 0.
  apred = 0.
  jrkpred = 0.

  r(1,1) = 0.9700436
  r(2,1) = -0.24308753
  r(1,2) = -1. * r(1,1)
  r(2,2) = -1. * r(2,1)

  v(1,1) = 0.466203685
  v(2,1) = 0.43236573
  v(1,2) = v(1,1)
  v(2,2) = v(2,1)
  v(1,3) = -2. * v(1,1)
  v(2,3) = -2. * v(2,1)
  v(3,3) = -2. * v(3,1)

  DO i=1,n
     IDs(i) = i
  END DO

  DO i=1,n
     ! Kinetic Energy
     vmagsqrd = 0.
     DO k=1,ndims
        vmagsqrd = vmagsqrd + (v(k,i)**2)
     END DO
     Ek = Ek + (m(i) * vmagsqrd / 2.)
  END DO
  
  DO i=1,n-1
     ! Potential Energy
     DO j=i+1,n
        dist = 0.
        DO k=1,ndims
           dist = dist + ((r(k,i)-r(k,j))**2) ! Will take square root AFTER this DO loop
        END DO
        dist = dist**0.5 ! dist is now the magnitude of the distance between objects i and j
        Ep = Ep - (G * m(i) * m(j) / dist)
     END DO
  END DO

  E0 = Ek + Ep ! Initial total energy

  DO
     ! Main Loop

     ! Resetting values that are calculated through incrementation
     Ek = 0.
     Ep = 0.
     a = 0.
     jrk = 0.
     s = 0.
     c = 0.
     apred = 0.
     jrkpred = 0.

     DO i=1,n
        amagsqrd = 0.
        DO j=1,n
           IF (i==j) CYCLE
           dist = 0.
           vdotr = 0.
           DO k=1,ndims
              dr(k) = r(k,i) - r(k,j)
              dv(k) = v(k,i) - v(k,j)
              vdotr = vdotr + (dv(k) * dr(k))
              dist = dist + (dr(k)**2)
           END DO
           dist = dist**0.5
           DO k=1,ndims
              a(k,i) = a(k,i) - (G * m(j) * dr(k) / (dist**3))
              amagsqrd = amagsqrd + (a(k,i)**2)
              jrk(k,i) = jrk(k,i) + (G * m(j) * ((dv(k)/(dist**3)) - (3.*vdotr*dr(k)/(dist**5))))
           END DO
        END DO
        dtOpt(i) = ((eta/amagsqrd)**0.5)
     END DO

     dt = MINVAL(dtOpt(:)) ! Adaptive timestep length

     DO i=1,n
        ! Position/velocity predictions. Separate DO loop to the acceleration/jerk calculations because those calculations will not give the correct results if any positions/velocities
        ! change during them
        r(:,i) = r(:,i) + (v(:,i) * dt) + (a(:,i) * (dt**2) / 2.) + (jrk(:,i) * (dt**3) / 6.)
        v(:,i) = v(:,i) + (a(:,i) * dt) + (jrk(:,i) * (dt**2) / 2.)
     END DO

     DO i=1,n
        ! Predicted accelerations/jerks at predicted positions/velocities
        DO j=1,n
           IF (i==j) CYCLE
           dist = 0.
           vdotr = 0.
           DO k=1,ndims
              dr(k) = r(k,i) - r(k,j)
              dv(k) = v(k,i) - v(k,j)
              vdotr = vdotr + (dv(k) * dr(k))
              dist = dist + (dr(k)**2)
           END DO
           dist = dist**0.5
           DO k=1,ndims
              apred(k,i) = apred(k,i) - (G * m(j) * dr(k) / (dist**3))
              jrkpred(k,i) = jrkpred(k,i) + (G * m(j) * ((dv(k)/(dist**3)) - (3.*vdotr*dr(k)/(dist**5))))
           END DO
        END DO
     END DO

     DO i=1,n
        ! Finding snap and crackle. Separate DO loop because calculating the snap/crackle of a given object this way requires knowledge of the predicted accelerations/jerks of all other objects
        DO j=1,n
           IF (i==j) CYCLE
           DO k=1,ndims
              s(k,i) = s(k,i) - (((6.*(a(k,j)-apred(k,j))) + (((4.*jrk(k,j))+(2.*jrkpred(k,j)))*dt)) / (dt**2)) ! The subtraction is outside the brackets so the jerk contribution inside the brackets here
                                                                                                                ! is an addition rather than a subtraction (and similarly for the acceleration contribution)
              c(k,i) = c(k,i) + (((12.*(a(k,j)-apred(k,j))) + (6.*(jrk(k,j)+jrkpred(k,j))*dt)) / (dt**3))
           END DO
        END DO
     END DO

     DO i=1,n
        ! Position/velocity corrections. The snap/crackle calculations in the previous DO loop do not involve the current positions/velocities of the relevant objects, so this section could be included
        ! at the end of the previous DO loop instead (after cycling through all j's for a given i), but I have separated them to make the code easier to read

        ! Note that the end-of-timestep positions/velocities were calculated to second order for the predictions earlier, and so only the snap and crackle contributions need to be considered here
        r(:,i) = r(:,i) + ((s(:,i)*(dt**4))/24.) + ((c(:,i)*(dt**5))/120.)
        v(:,i) = v(:,i) + ((s(:,i)*(dt**3))/6.) + ((c(:,i)*(dt**4))/24.)
     END DO

     t = t + dt
     tcount = tcount + dt

     ! Energy error calculation
     DO i=1,n
        ! Kinetic Energy
        vmagsqrd = 0.
        DO k=1,ndims
           vmagsqrd = vmagsqrd + (v(k,i)**2)
        END DO
        Ek = Ek + (m(i) * vmagsqrd / 2.)
     END DO

     DO i=1,n-1
        ! Potential Energy
        DO j=i+1,n
           dist = 0.
           DO k=1,ndims
              dist = dist + ((r(k,i)-r(k,j))**2)
           END DO
           dist = dist**0.5
           Ep = Ep - (G * m(i) * m(j) / dist)
        END DO
     END DO

     E = Ek + Ep ! Total energy
     Eerr = (E-E0)/E0
     Eerr = ABS(Eerr)

     t = t + dt
     tcount = tcount + dt

     ! Write to files
     IF (t >= tlim) THEN
        WRITE(7,*) r(1,:)
        WRITE(8,*) r(2,:)
        WRITE(9,*) Eerr, t
        EXIT    
     ELSE IF (tcount >= twrite) THEN
        WRITE(7,*) r(1,:)
        WRITE(8,*) r(2,:)
        WRITE(9,*) Eerr, t
        tcount = tcount - twrite
     END IF

  END DO

END PROGRAM fig_eight

     
     
     
  
  
  
