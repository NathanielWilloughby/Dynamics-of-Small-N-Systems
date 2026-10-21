PROGRAM plummer_setup_test
    IMPLICIT NONE
    INTEGER :: n, i, j, ndims
    DOUBLE PRECISION :: pi, ranNum, rad, mTot, twothirds, sqrtr2z2, vesc, root2, gq, q2, velmag
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: m, origDist, mwithin
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: pos, vel

    OPEN(7, file='xpositions.txt', status='replace')
    OPEN(8, file='ypositions.txt', status='replace')
    OPEN(9, file='zpositions.txt', status='replace')
    OPEN(10, file='xvelocities.txt', status='replace')
    OPEN(11, file='yvelocities.txt', status='replace')
    OPEN(12, file='zvelocities.txt', status='replace')
    OPEN(13, file='mass_within.txt', status='replace')
    OPEN(14, file='origin_distances.txt', status='replace')

    n = 1000
    ndims = 3

    ALLOCATE(m(1:n))
    ALLOCATE(pos(1:ndims, 1:n))
    ALLOCATE(vel(1:ndims, 1:n))
    ALLOCATE(origDist(1:n))
    ALLOCATE(mwithin(1:n))


    mTot = 1.           ! Total mass is 1 in N-body units
    m = mTot/n          ! For now, all masses are the same
    pos = 0.            ! Positions
    vel = 0.            ! Velocities
    pi = 4. * ATAN(1.)  ! Pi
    ranNum = 0.         ! Random number
    rad = 0.            ! Distance of a given object from the origin
    twothirds = 2./3.   ! 2/3
    sqrtr2z2 = 0.       ! SQRT((r**2) - (z**2)). ALSO used to represent SQRT((velmag**2) - (v_z**2)) despite the variable name not implying this
    vesc = 0.           ! Escape velocity
    root2 = SQRT(2.)    ! Square root of 2
    gq = 0.             ! g(q) = q^2 * (1-q^2)^7/2
    q2 = 0.             ! q squared (q is vel/vesc, and is a random value between 0 and 1 that is used to randomise initial velocities)
    velmag = 0.         ! Magnitude of velocity

    ! The following variables will not be used when this code is incorporated into a full N-body code
    origDist = 0.       ! Distance from origin
    mwithin = 0.        ! Mass enclosed within a sphere of radius origDist


    ! Plummer sphere setup method from Aarseth, Hénon and Wielen (1974)
    DO i=1,n
        CALL RANDOM_NUMBER(ranNum)
        rad = SQRT(1./((ranNum**(-1.*twothirds)) - 1.))
        CALL RANDOM_NUMBER(ranNum)
        pos(3,i) = rad * (1. - (2.*ranNum))
        CALL RANDOM_NUMBER(ranNum) ! Need a different random number for the x and y coordinates to the one used for the z coordinates
        sqrtr2z2 = SQRT((rad**2) - (pos(3,i)**2))
        pos(1,i) = COS(2. * pi * ranNum) * sqrtr2z2
        pos(2,i) = SIN(2. * pi * ranNum) * sqrtr2z2
        vesc = root2 * ((1. + (rad**2))**(-1.*0.25))
        DO
            ! This will run continuously until a certain condition is met (this avoids having to use GOTO)
            CALL RANDOM_NUMBER(ranNum)
            q2 = ranNum ** 2
            gq = q2 * ((1. - q2)**3.5)
            CALL RANDOM_NUMBER(ranNum)
            IF ((0.1 * ranNum) < gq) THEN
                EXIT ! Ends the current DO loop (which would otherwise run infinitely)
            END IF
        END DO
        velmag = ranNum * vesc
        CALL RANDOM_NUMBER(ranNum)
        vel(3,i) = velmag * (1. - (2. * ranNum))
        CALL RANDOM_NUMBER(ranNum)
        sqrtr2z2 = SQRT((velmag**2) - (vel(3,i)**2))
        vel(1,i) = COS(2. * pi * ranNum) * sqrtr2z2
        vel(2,i) = SIN(2. * pi * ranNum) * sqrtr2z2
        origDist(i) = SQRT((pos(1,i)**2) + (pos(2,i)**2) + (pos(3,i)**2))
    END DO

    ! Enclosed mass
    DO i=1,n
        ! Do not need to reset mwithin to zero as it should already be zero due to not having been used before this point in the execution of the code
        DO j=1,n
            IF (origDist(j) > origDist(i)) CYCLE
            mwithin(i) = mwithin(i) + m(j) ! This should only be carried out when either i=j or when object j is closer to the origin than object i
        END DO
    END DO


    WRITE(7,*) pos(1,:)
    WRITE(8,*) pos(2,:)
    WRITE(9,*) pos(3,:)
    WRITE(10,*) vel(1,:)
    WRITE(11,*) vel(2,:)
    WRITE(12,*) vel(3,:)
    WRITE(13,*) mwithin(:)
    WRITE(14,*) origDist(:)

END PROGRAM plummer_setup_test



