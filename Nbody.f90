PROGRAM n_body_main
    IMPLICIT NONE
    INTEGER :: n, i, j, k, ndims
    DOUBLE PRECISION :: t, twrite, tlim, tcount, dt, dt2, dt3, dt4, dt5, recipdt2, recipdt3, eta, G, Greal, recip6, recip24, recip120
    DOUBLE PRECISION :: E, E0, Ek, Ep, Eerr, dist, dist2, dist3, dist5, vmag2, amag2, vdotr, mTot
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, dtOpt, comPos, comVel, dr, dv
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pos, vel, acc, jrk, snp, crk, accPred, jrkPred

    ! Opening the files that will be written to. This block of code can be modified as long as the file output section later on is modified correspondingly
    OPEN(7, file='xpositions.txt', status='replace')
    OPEN(8, file='ypositions.txt', status='replace')
    OPEN(9, file='zpositions.txt', status='replace')
    OPEN(10, file='energy_error.txt', status='replace')

    ! Variables needed for array size allocation - the values of these CAN be changed (although ndims probably shouldn't be)
    n = 1000        ! Number of objects
    ndims = 3       ! Number of dimensions

    ! Array size allocation
    ALLOCATE(m(1:n))
    ALLOCATE(pos(1:ndims, 1:n))
    ALLOCATE(vel(1:ndims, 1:n))
    ALLOCATE(acc(1:ndims, 1:n))
    ALLOCATE(jrk(1:ndims, 1:n))
    ALLOCATE(snp(1:ndims, 1:n))
    ALLOCATE(crk(1:ndims, 1:n))
    ALLOCATE(dtOpt(1:n))
    ALLOCATE(accPred(1:ndims, 1:n))
    ALLOCATE(jrkPred(1:ndims, 1:n))
    ALLOCATE(comPos(1:ndims))
    ALLOCATE(comVel(1:ndims))
    ALLOCATE(dr(1:ndims))
    ALLOCATE(dv(1:ndims))

    ! Constants
    G = 1.              ! Gravitational constant, which is set to 1 in N-body units
    recip6 = 1./6.      ! Multiply by this instead of dividing by 6
    recip24 = 1./24.    ! Multiply by this instead of dividing by 24
    recip120 = 1./120.  ! Multiply by this instead of dividing by 120

    ! Adjustable values - given in SI here as it is more intuitive to decide on SI values then adjust to N-body units than deciding on values that are
    ! already in N-body units
    twrite = 0.01       ! File write interval (in seconds)
    tlim = 100.         ! Total simulated time limit (in seconds)
    Greal = 6.67e-11    ! Gravitational constant in SI
    eta = 1e-5          ! Accuracy parameter used in calculating the optimal timestep lengths for each object. All timestep lengths are scaled by SQRT(eta)

    ! Initialising other variables
    t = 0.              ! Time
    tcount = 0.         ! Tracks the time between file writes
    dt = 0.             ! Timestep length, will change for each timestep
    dt2 = 0.            ! dt^2
    dt3 = 0.            ! dt^3
    dt4 = 0.            ! dt^4
    dt5 = 0.            ! dt^5
    recipdt2 = 0.       ! dt^-2
    recipdt3 = 0.       ! dt^-3
    dtOpt = 0.          ! Optimal timestep length for each object
    pos = 0.            ! Positions
    vel = 0.            ! Velocities
    acc = 0.            ! Accelerations
    jrk = 0.            ! Jerks
    snp = 0.            ! Snaps
    crk = 0.            ! Crackles
    m = 0.              ! Masses
    comPos = 0.         ! COM origin position
    comVel = 0.         ! COM origin velocity
    E0 = 0.             ! Initial total energy
    Ek = 0.             ! Instantaneous total kinetic energy
    Ep = 0.             ! Instantaneous total potential energy
    E = 0.              ! Instantaneous total energy
    Eerr = 0.           ! Instantaneous fractional difference between E and E0
    vmag2 = 0.          ! Velocity magnitude squared (vel . vel)
    amag2 = 0.          ! Acceleration magnitude squared (acc. acc)
    dr = 0.             ! Distance vector between two objects (pos_i - pos_j)
    dv = 0.             ! Difference in velocity between two objects (vel_i - vel_j)
    vdotr = 0.          ! Dot product of dv and dr
    dist = 0.           ! Absolute distance between two objects
    dist2 = 0.          ! dist^2
    dist3 = 0.          ! dist^3
    dist5 = 0.          ! dist^5. Note the absence of dist4 due to dist^4 not being used anywhere in the simulation
    mTot = 0.           ! Total mass


    ! Initial conditions
    m = 1.99e30         ! All objects have a mass of 1 solar mass for now
    CALL setup_plummer(n, pos, vel) ! Sets up a cluster of n stars distributed in a Plummer sphere

    ! Conversion to COM frame - Note that mTot, comPos and comVel do not need to be reset to zero here as they were initialised to be zero and will not have had their
    ! values modified yet
    DO i=1,n
        mTot = mTot + m(i)
        DO k=1,ndims
            comPos(k) = comPos(k) + (m(i) * pos(k,i))
            comVel(k) = comVel(k) + (m(i) * vel(k,i))
        END DO
    END DO
    
    DO k=1, ndims
        comPos(k) = comPos(k) / mTot
        comVel(k) = comVel(k) / mTot
    END DO
    
    DO i=1,n
        DO k=1,ndims
            pos(k,i) = pos(k,i) - comPos(k)
            vel(k,i) = vel(k,i) - comVel(k)
        END DO
    END DO


    ! N-body unit conversion

    ! Initial total energy. Note that Ek and Ep do not need to be reset to zero here as they will not have been changed from zero as of this point in the simulation
    ! because the Ek and Ep used in the n_body_conversion subroutine are variables that are local to that subroutine
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
            WRITE(9,*) pos(3,:)
            WRITE(10,*) Eerr, t
            EXIT

        ELSE IF (tcount >= twrite) THEN

            ! Energy error
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
            WRITE(9,*) pos(3,:)
            WRITE(10,*) Eerr, t

            tcount = tcount - twrite ! Not just resetting tcount to zero because this could cause the file outputs to drift from the intended times over the course
                                     ! of the simulation. Subtracting twrite keeps the file outputs at times as close to multiples of twrite as possible
        END IF

    END DO

END PROGRAM n_body_main

! Plummer sphere (3 dimensions only). Method from Aarseth, Hénon and Wielen (1974)
SUBROUTINE setup_plummer(n, pos, vel)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: i
    DOUBLE PRECISION, INTENT(OUT) :: pos(1:3, 1:n), vel(1:3, 1:n)
    DOUBLE PRECISION :: ranNum, ranNum2, twopi, twothirds, rad, sqrtr2z2, vesc, root2, q2, gq, vmag

    twopi = 8. * ATAN(1.)   ! 2 * pi
    ranNum = 0.             ! Random number in a uniform distribution between 0 and 1
    ranNum2 = 0.            ! Another random number in a uniform distribution between 0 and 1, used when a new random number needs to be generated while the value stored
                            ! in ranNum is also still needed
    twothirds = 2./3.       ! 2/3
    rad = 0.                ! Distance of a given object from the origin. Not an array as the generation of initial conditions for a given object only needs the value
                            ! of rad for the object itself
    sqrtr2z2 = 0.           ! SQRT(r^2 - z^2). Also stores SQRT(v^2 - v_z^2) as this is value is only needed for a given object once SQRT(r^2 - z^2) is no longer needed
    vesc = 0.               ! Escape velocity at a given object's position, which represents the maximum 'allowed' velocity at that position. Not an array for the same
                            ! reasons as rad
    vmag = 0.               ! Velocity magnitude of a given object. Not an array for the same reasons as rad
    q2 = 0.                 ! q^2, where q = vmag/vesc
    gq = 0.                 ! g(q) = q^2 - (1-q^2)^7/2
    root2 = SQRT(2.)        ! Square root of 2

    DO i=1,n
        CALL RANDOM_NUMBER(ranNum)
        rad = SQRT(1./((ranNum**(-1.*twothirds)) - 1.))
        CALL RANDOM_NUMBER(ranNum)
        pos(3,i) = rad * (1. - (2.*ranNum))
        CALL RANDOM_NUMBER(ranNum) ! Need a different random number for the x and y coordinates to the one used for the z coordinates
        sqrtr2z2 = SQRT((rad**2) - (pos(3,i)**2))
        pos(1,i) = COS(twopi * ranNum) * sqrtr2z2
        pos(2,i) = SIN(twopi * ranNum) * sqrtr2z2
        vesc = root2 * ((1. * (rad**2))**(-0.25))
        DO
            ! This will run continuously until a certain condition is met (this avoids having to use GOTO)
            CALL RANDOM_NUMBER(ranNum)
            q2 = ranNum**2
            gq = q2 * ((1. - q2)**3.5)
            CALL RANDOM_NUMBER(ranNum2)
            IF ((0.1 * ranNum2) < gq) THEN
                EXIT ! This ends the current DO loop (which would otherwise run infinitely)
            END IF
        END DO
        vmag = ranNum * vesc ! Note this uses ranNum, NOT ranNum2
        CALL RANDOM_NUMBER(ranNum)
        vel(3,i) = vmag * (1. - (2.*ranNum))
        CALL RANDOM_NUMBER(ranNum) ! Need a different random number for the x and y velocity components to the one used for the z velocity component
        sqrtr2z2 = SQRT((vmag**2) - (vel(3,i)**2)) ! Reusing this variable name as the value it was set to earlier in the loop is no longer needed at this point
        vel(1,i) = COS(twopi * ranNum) * sqrtr2z2
        vel(2,i) = SIN(twopi * ranNum) * sqrtr2z2
    END DO

    RETURN

END SUBROUTINE setup_plummer

SUBROUTINE n_body_conversion(n, ndims, mTot, pos, vel, m) ! CURRENTLY UNFINISHED AND NOT IMPLEMENTED
    IMPLICIT NONE
    INTEGER :: i, j, k
    INTEGER, INTENT(IN) :: n, ndims
    DOUBLE PRECISION :: E0, Ek, Ep, vmag2, dist2, G
    DOUBLE PRECISION, INTENT(INOUT) :: pos(1:ndims, 1:n), vel(1:ndims, 1:n), m(1:n), mTot

    ! Energies - these will not be passed out of this subroutine as the initial total energy will be recalculated in the main program as a sense check
    ! After converting to N-body units, the total energy should be -0.25 (Aarseth, 2003)
    E0 = 0.
    Ek = 0.
    Ep = 0.

    vmag2 = 0. ! Velocity magnitude squared (vel. vel)
    dist2 = 0. ! Distance between two objects squared

    G = 1. ! Gravitational constant in N-body units

    DO i=1,n
        m(i) = m(i) / mTot ! Scales masses so that the total mass is 1 N-body mass unit
    END DO

    mTot = 1. ! Scales the value of mTot to 1 N-body mass unit. Must be done AFTER scaling the individual masses

    ! Energy
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
            Ep = Ep - ((G * m(i) * m(j)) / SQRT(dist2))
        END DO
    END DO

    E0 = Ek + Ep

    RETURN

END SUBROUTINE n_body_conversion



