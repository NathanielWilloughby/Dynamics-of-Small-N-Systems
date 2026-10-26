PROGRAM fig_eight_test
    IMPLICIT NONE
    INTEGER :: n, i, j, k, ndims
    DOUBLE PRECISION :: t, tcount, twrite, tlim, dt, dt2, dt3, dt4, dt5, eta, G, pi, mTot, E0, Ek, Ep, E, Eerr
    DOUBLE PRECISION :: vmag2, dist, dist2, dist3, dist5, amag2, vdotr, recip6, recip24, recip120, recipdt2, recipdt3
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, dtOpt, comPos, comVel, dr, dv
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pos, vel, acc, jrk, snp, crk, accPred, jrkPred

    OPEN(7, file='xpositions.txt', status='replace')
    OPEN(8, file='ypositions.txt', status='replace')
    OPEN(9, file='energy_error.txt', status='replace')

    n = 3               ! Number of objects
    ndims = 3           ! Number of dimensions

    ALLOCATE(m(1:n))
    ALLOCATE(pos(1:ndims, 1:n))
    ALLOCATE(vel(1:ndims, 1:n))
    ALLOCATE(acc(1:ndims, 1:n))
    ALLOCATE(jrk(1:ndims, 1:n))
    ALLOCATE(snp(1:ndims, 1:n))
    ALLOCATE(crk(1:ndims, 1:n))
    ALLOCATE(dtOpt(1:n))
    ALLOCATE(comPos(1:ndims))
    ALLOCATE(comVel(1:ndims))
    ALLOCATE(accPred(1:ndims, 1:n))
    ALLOCATE(jrkPred(1:ndims, 1:n))
    ALLOCATE(dr(1:ndims))
    ALLOCATE(dv(1:ndims))

    t = 0.              ! Time
    tcount = 0.         ! Tracks the time between file writes
    twrite = 0.         ! Time between file writes, to be given a non-zero value later
    tlim = 0.           ! Total simulated time, to be given a non-zero value later
    dt = 0.             ! Timestep length, will vary each timestep according to the acceleration of each object, and will be made non-zero before 
                        ! the first time it is used in calculations
    dt2 = 0.            ! dt^2
    dt3 = 0.            ! dt^3
    dt4 = 0.            ! dt^4
    dt5 = 0.            ! dt^5
    recipdt2 = 0.       ! dt^-2
    recipdt3 = 0.       ! dt^-3
    dtOpt = 0.          ! Optimal timestep lengths for each object, which depends on the current acceleration of that object
    eta = 1e-10         ! Used in calculating the optimal timestep lengths for each object, all of the timestep lengths are scaled by SQRT(eta)
    G = 1.              ! Gravitational constant, which is set to 1 in N-body units
    pi = 4. * ATAN(1.)  ! Pi to maximum computational precision
    pos = 0.            ! Positions
    vel = 0.            ! Velocities
    acc = 0.            ! Accelerations
    jrk = 0.            ! Jerks
    snp = 0.            ! Snaps
    crk = 0.            ! Crackles
    m = 0.              ! Masses
    mTot = 0.           ! Total mass, will be calculated when calculating the position and velocity of the origin of the COM frame
    comPos = 0.         ! COM origin position, will be calculated at the same time as mTot
    comVel = 0.         ! COM origin velocity, will be calculated at the same time as mTot
    E0 = 0.             ! Initial total energy, will be calculated after moving to COM frame
    Ek = 0.             ! Tracks the total kinetic energy
    Ep = 0.             ! Tracks the total potential energy
    E = 0.              ! Tracks the total energy
    Eerr = 0.           ! Fractional difference between E0 and E
    vmag2 = 0.          ! Velocity magnitude squared (vel.vel)
    dist = 0.           ! Temporarily stores the distance between two objects
    dist2 = 0.          ! Temporarily stores dist^2
    dist3 = 0.          ! Temporarily stores dist^3
    dist5 = 0.          ! Temporarily stores dist^5. Note that dist^4 is not used in the simulation at all so there is no dist4 variable
    accPred = 0.        ! Acceleration at initially predicted end-of-timestep positions and velocities (before accounting for snap and crackle)
    jrkPred = 0.        ! Same as accPred but for the jerk
    amag2 = 0.          ! Acceleration magnitude squared (acc.acc)
    dr = 0.             ! Disrance vector between two objects (pos_i - pos_j)
    dv = 0.             ! Difference in velocity between two objects (vel_i - vel_j)
    vdotr = 0.          ! Dot product of dv and dr
    recip6 = 1./6.      ! Multiply by this instead of dividing by 6
    recip24 = 1./24.    ! Multiply by this instead of dividing by 24
    recip120 = 1./120.  ! Multiply by this instead of dividing by 120

    ! Time limits
    twrite = 0.01   ! Write to files every 0.01 N-body time units
    tlim = 100.     ! Run for 100 N-body time units

    ! Initial conditions
    m = 1. ! All objects have the same mass of 1 N-body mass unit

    pos(1,1) = 0.9700436
    pos(2,1) = -0.24308753
    pos(1,2) = -1. * pos(1,1)
    pos(2,2) = -1. * pos(2,1)

    vel(1,1) = 0.466203685
    vel(2,1) = 0.43236573
    vel(1,2) = vel(1,1)
    vel(2,2) = vel(2,1)
    vel(1,3) = -2. * vel(1,1)
    vel(2,3) = -2. * vel(2,1)

    ! COM frame calculation
    DO i=1,n
        mTot = mTot + m(i)
        DO k=1,ndims
            comPos(k) = comPos(k) + (m(i) * pos(k,i))
            comVel(k) = comVel(k) + (m(i) * vel(k,i))
        END DO
    END DO
    comPos(:) = comPos(:) / mTot
    comVel(:) = comVel(:) / mTot

    ! Moving into COM frame
    DO i=1,n
        DO k=1,ndims
            pos(k,i) = pos(k,i) - comPos(k)
            vel(k,i) = vel(k,i) - comVel(k)
        END DO
    END DO

    ! Initial total energy. Note Ek and Ep do not need to be reset to zero here as they will not have been changed from zero as of this point in the simulation
    DO i=1,n
        vmag2 = 0. ! Needs to be reset to zero for each iteration
        DO k=1,ndims
            vmag2 = vmag2 + (vel(k,i)**2)
        END DO
        Ek = Ek + (0.5 * m(i) * vmag2)
        
        DO j=i+1,n ! These boundaries avoid double-counting pairs of objects in the potential energy calculations
            dist2 = 0. ! Needs to be reset to zero for each iteration
            DO k=1,ndims
                dist2 = dist2 + ((pos(k,i)-pos(k,j))**2)
            END DO
            dist = SQRT(dist2) ! Do not need to reset dist to zero each iteration here as it is not calculated directly through incrementation
            Ep = Ep - ((G * m(i) * m(j)) / dist)
        END DO
    END DO
    
    E0 = Ek + Ep

    ! Main Loop
    DO
        ! Resetting values that are calculated through incrementation over all objects
        Ek = 0.
        Ep = 0.
        acc = 0.
        jrk = 0.
        snp = 0.
        crk = 0.
        accPred = 0.
        jrkPred = 0.

        ! Accelerations, jerks and adaptive timestep length
        DO i=1,n
            amag2 = 0. ! Needs to be reset to zero for each iteration
            DO j=1,n
                IF (i==j) CYCLE ! Finding the contribution of particle i to the acceleration and jerk of particle i is redundant
                dist2 = 0.
                vdotr = 0.
                DO k=1,ndims
                    dr(k) = pos(k,i) - pos(k,j)
                    dv(k) = vel(k,i) - vel(k,j)
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                DO k=1,ndims
                    acc(k,i) = acc(k,i) - ((G * m(j) * dr(k)) / dist3)
                    amag2 = amag2 + (acc(k,i)**2)
                    jrk(k,i) = jrk(k,i) + (G * m(j) * ((dv(k) / dist3) - ((3. * vdotr * dr(k)) / dist5)))
                END DO
            END DO
            dtOpt(i) = SQRT(eta / amag2)
        END DO
        dt = MINVAL(dtOpt(:))
        dt2 = dt**2
        dt3 = dt**3
        dt4 = dt**4
        dt5 = dt**5
        recipdt2 = 1. / dt2
        recipdt3 = 1. / dt3

        ! Predicted end-of-timestep positions and velocities (before snap and crackle are accounted for)
        DO i=1,n
            ! Predicted positions/velocities do not to be in separate arrays to the pre-timestep positions/velocities.
            ! This is because the pre-timestep positions/velocities are not needed for the remainder of the timestep after the predicted
            ! positions/velocities have been calculated
            pos(:,i) = pos(:,i) + (vel(:,i) * dt) + (0.5 * acc(:,i) * dt2) + (recip6 * jrk(:,i) * dt3)
            vel(:,i) = vel(:,i) + (acc(:,i) * dt) + (0.5 * jrk(:,i) * dt2) ! Velocity MUST be updated AFTER position as updating the position requires
                                                                           ! the pre-timestep velocity
        END DO

        ! Predicted accelerations and jerks at the predicted positions and velocities
        DO i=1,n
            ! Predicted accelerations/jerks do need separate arrays to the accelerations/jerks used in predicting the end-of-timestep positions/velocities
            ! because both pairs are needed to estimate the snaps/crackles
            DO j=1,n
                IF(i==j) CYCLE
                dist2 = 0.
                vdotr = 0.
                DO k=1,ndims
                    dr(k) = pos(k,i) - pos(k,j)
                    dv(k) = vel(k,i) - vel(k,j)
                    vdotr = vdotr + (dv(k) * dr(k))
                    dist2 = dist2 + (dr(k)**2)
                END DO
                dist = SQRT(dist2)
                dist3 = dist**3
                dist5 = dist**5
                DO k=1,ndims
                    accPred(k,i) = accPred(k,i) - ((G * m(j) * dr(k)) / dist3)
                    jrkPred(k,i) = jrkPred(k,i) + (G * m(j) * ((dv(k) / dist3) - ((3. * vdotr * dr(k)) / dist5)))
                END DO
            END DO
        END DO

        ! Estimating the snap and crackle of each object. Separate DO loop to the calculation of accPred and jrkPred because the estimation of the snap and crackle
        ! for a given object requires knowledge of accPred and jrkPred for all other objects
        DO i=1,n
            DO j=1,n
                IF (i==j) CYCLE
                DO k=1,ndims
                    ! Note that the subtraction in the snap estimation is outside the parentheses so the jerk term is added rather than subtracted
                    snp(k,i) = snp(k,i) - (recipdt2 * ((6. * (acc(k,j) - accPred(k,j))) + (dt * ((4. * jrk(k,j)) + (2. * jrkPred(k,j))))))
                    crk(k,i) = crk(k,i) + (recipdt3 * ((12. * (acc(k,j) - accPred(k,j))) + (6. * dt * (jrk(k,j) + jrkPred(k,j)))))
                END DO
            END DO
        END DO

        ! Position and velocity corrections. This COULD be including in the previous DO loop after cycling through all j's for a given i as the snap/crackle
        ! calculations in the actual loop do not involve the predicted positions/velocities of the relevant objects (the contributions of the predicted positions/velocities
        ! are contained within accPred and jrkPred, which were calculated beforehand). Updating the positions and velocities during the same DO loop as
        ! the snaps and crackles are esimated shouldn't affect the snap and crackle estimations, and won't have any effect on the calculated positions and velocities
        ! as long as they are updated AFTER cycling through all j's for a given i. I have decided against putting them in the same DO loop
        ! in order to make the code easier to follow
        DO i=1,n
            ! The contributions from orders up to the jerk are already included within the current values of pos and vel for a given object before they are updated here
            ! so only the contributions from the snap and crackle need to be accounted for here
            pos(:,i) = pos(:,i) + (recip24 * snp(:,i) * dt4) + (recip120 * crk(:,i) * dt5)
            vel(:,i) = vel(:,i) + (recip6 * snp(:,i) * dt3) + (recip24 * crk(:,i) * dt4)
        END DO

        ! Time incrementation

        t = t + dt
        tcount = tcount + dt

        ! File writes
        IF (t >= tlim) THEN ! This takes priority over the condition (tcount >= twrite) in the event that both are satisfied
            
            ! Energy error calculation. This will be re-used in the ELSE IF segment after this, but I don't want it to be calculated on every single timestep
            ! so increasing the number of lines of code by re-writing this block of code again later will decrease the runtime
            Ek = 0.
            Ep = 0.
            DO i=1,n
                vmag2 = 0.
                DO k=1,ndims
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                Ek = Ek + (0.5 * m(i) * vmag2)

                DO j=i+1,n
                    dist2 = 0.
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    Ep = Ep - ((G * m(i) * m(j)) / dist)
                END DO
            END DO
            E = Ek + Ep
            Eerr = ABS((E-E0)/E0)


            ! File writes
            WRITE(7,*) pos(1,:)
            WRITE(8,*) pos(2,:)
            WRITE(9,*) Eerr, t
            EXIT

        ELSE IF (tcount >= twrite) THEN

            ! Energy error
            Ek = 0.
            Ep = 0.
            DO i=1,n
                vmag2 = 0.
                DO k=1,ndims
                    vmag2 = vmag2 + (vel(k,i)**2)
                END DO
                Ek = Ek + (0.5 * m(i) * vmag2)

                DO j=i+1,n
                    dist2 = 0.
                    DO k=1,ndims
                        dist2 = dist2 + ((pos(k,i) - pos(k,j))**2)
                    END DO
                    dist = SQRT(dist2)
                    Ep = Ep - ((G * m(i) * m(j)) / dist)
                END DO
            END DO
            E = Ek + Ep
            Eerr = ABS((E-E0)/E0)


            ! File writes
            WRITE(7,*) pos(1,:)
            WRITE(8,*) pos(2,:)
            WRITE(9,*) Eerr, t

            tcount = tcount - twrite ! Not just resetting tcount to zero because this could cause the file outputs to drift from the intended times over the course
                                     ! of the simulation. Subtracting twrite keeps the file outputs at times as close to multiples of twrite as possible
        END IF

    END DO

END PROGRAM fig_eight_test