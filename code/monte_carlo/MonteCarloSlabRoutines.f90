module MonteCarloSlabRoutines
    implicit none
    type Input ! Data structure for holding all relevant problem information
        ! cell boundaries
        double precision, allocatable, dimension(:)   :: CellBounds 
        ! material id's in each cell
        integer, allocatable, dimension(:)            :: MatId
        ! cross-sections
        double precision, allocatable, dimension(:,:) :: Sigma
        integer :: NumCells, & ! number of cells
                   NumMats,  & ! number of different materials
                   NumHist     ! number of histories
    end type Input
    type Output ! Data structure for holding the relevant tally data
        double precision, allocatable, dimension(:) :: TrakEst, CollEst
    end type Output
contains

subroutine PlayGame(in,out)
!-------------------------------------------------------------------------------
! subroutine PlayGame(in,out)
! Play the game, or do the simulation, or let particles roam free.
!   Inputs:
!       in          -- data structure containing relevant input data
!   Outputs:
!       out         -- data structure of output indomation
!       Outp%trakest -- track length estimator tally
!       Outp%collest -- collision estimator tally
    type(Input), intent(in)         :: Inpt
    type(Output), intent(inout)     :: Outp
!       Local Variables
    double precision, allocatable, dimension(:) :: tmptrak,& ! temporary
                                                   tmpcoll   ! tallies
    integer ::        alive,        & ! Whether or not I'm alive
                      coll,         & ! Whether or not I'm exiting a collision
                      mycell,       & ! My current cell
                      neighbor,     & ! My neighbor
                      n,m=0

    double precision  mu,           & ! My direction
                      x,            & ! My location
                      SigT,         & ! SigmaT of mycell
                      mfps,         & ! Mean free paths to travel
                      dist,         & ! Distance along x axis to travel
                      d2neighbor      ! Distance (along x) to neighbors)
    double precision, parameter :: one = 1.0d0 ! for use in "sign" function
!-------------------------------------------------------------------------------
    print *, 'beginning...'
    ! Initialize built-in random number generator
    call random_seed()   
    ! Set the estimator tallies
    Outp%TrakEst = 0
    Outp%CollEst = 0
    ! Allocate the temporary tally arrays
    allocate( tmptrak(Inpt%NumCells), tmpcoll(Inpt%NumCells) )

    !--------------------------------------------------------------------------
    ! Let the histories begin!
    !--------------------------------------------------------------------------
    histories: do n = 1, Inpt%NumHist
        ! Get my starting location and direction
        call Source(in,x,mu)
        ! Reset temporary tallies
        tmptrak = 0.0
        tmpcoll = 0.0
        alive   = 1  
        coll    = 0
        call GetMyCell(x,in,mycell) ! Get my cell

        !----------------------------------------------------------------------
        ! Welcome to the life of a neutron!
        !----------------------------------------------------------------------
        do while ( alive .eq. 1 )
            if ( coll .eq. 1 ) then
                ! I need a new direction in life
                mu = 2*rand()-1
            end if
            ! Determine my neighbor and how far to her
            call GetMyNeighbor(mycell,x,mu,in,neighbor,d2neighbor)
            ! Sample how many mfp's I will go
            mfps    = -log(rand())
            ! Determine my total cross-section
            SigT = Inpt%Sigma(Inpt%MatId(mycell),1)
            ! Determine how far along x-axis I go
            dist = mfps*mu/SigT
            ! Determine whether I reach surface or collide
            if ( abs(dist) .gt. d2neighbor) then
                ! Then I pass to the next cell, so update track length tally
                tmptrak(mycell) = tmptrak(mycell) + abs(d2neighbor/mu)
                if  ( (neighbor.eq.0) .or. (neighbor.eq.(Inpt%NumCells+1)) ) then
                    ! I leaked out of left or right
                    alive = 0
                else
                    ! I get to the boundary
                    mycell  = neighbor
                    x       = x + sign(one,mu)*d2neighbor 
                    coll    = 0
                end if
            else
                ! Otherwise, I collide, so update both tallies
                tmpcoll(mycell) = tmpcoll(mycell) + 1
                tmptrak(mycell) = tmptrak(mycell) + mfps/SigT
                ! Update my location
                x = x + dist
                ! Sample the actual collision
                if ( rand() .lt. Inpt%Sigma(Inpt%MatId(mycell),2)/SigT ) then
                    ! I am absorbed
                    alive = 0
                else
                    ! I am scattered
                    coll = 1
                end if
            end if 
        end do ! while
        !----------------------------------------------------------------------
        ! R.I.P.
        !----------------------------------------------------------------------

        Outp%TrakEst = Outp%TrakEst + tmptrak ! Update tallies; these are "S1"
        Outp%CollEst = Outp%CollEst + tmpcoll ! counters; S2 would be here too!

    end do histories
    !--------------------------------------------------------------------------
    ! Game Over.
    !--------------------------------------------------------------------------
    print *, '...ending'

end subroutine PlayGame

subroutine Source(Inpt,x,mu)
!-------------------------------------------------------------------------------
! subroutine Source(in,x,mu)
! Defines the source distribution.
!   Inputs:
!       in  -- data structure containing relevant input data
!   Outputs:
!       x   -- my current position
!      mu   -- my current direction (i.e. cosine w/r to x axis).
    type(Input), intent(in)         :: Inpt
    double precision, intent(out)   :: x, mu
!-------------------------------------------------------------------------------
    ! Example: We are born in the middle 4 slabs unidomly and
    !          isotropically (in the LAB!)
    x  = rand()*4.0+4.0
    mu = 2.0*rand()-1.0
end subroutine Source

subroutine GetMyNeighbor(mycell,x,mu,Inpt,neighbor,d2neighbor)
!-------------------------------------------------------------------------------
! subroutine GetMyNeighbor(mycell,x,mu,Inpt,neighbor,d2neighbor)
! Finds out who my neighboring cell is.
!   Inputs:
!       mycell      -- my current home
!       x           -- my current x coordinate
!       mu          -- my current direction
!       in          -- data structure containing relevant input data
!   Outputs:
!       neighbor	-- the next cell over (in my direction)
!       d2neighbor  -- distance along x axis bedoe I enter neighbor
    integer, intent(in)             :: mycell
    double precision, intent(in)    :: x, mu
    type(Input), intent(in)         :: Inpt
    integer, intent(out)            :: neighbor
    double precision, intent(out)   :: d2neighbor
!-------------------------------------------------------------------------------
    if ( mu .gt. 0.0 ) then
        d2neighbor = Inpt%CellBounds(mycell+1) - x
        neighbor   = mycell + 1
    else
        d2neighbor = x - Inpt%CellBounds(mycell)
        neighbor   = mycell - 1;
    end if
end subroutine GetMyNeighbor

subroutine GetMyCell(x,Inpt,mycell)
!-------------------------------------------------------------------------------
! subroutine GetMyCell(x,in,mycell)
! Given my x coordinate, find out my cell.
!   Inputs:
!       x           -- my current x coordinate
!       in          -- data structure containing relevant input data
!   Outputs:
!       mycell      -- my current home
    double precision, intent(in)    :: x
    type(Input), intent(in)         :: Inpt
    integer, intent(out)            :: mycell
!-------------------------------------------------------------------------------
    do mycell = 1, Inpt%NumCells
        if ( x < Inpt%CellBounds(mycell+1) ) then
            exit
        end if
    end do
end subroutine GetMyCell

double precision function rand()
!-------------------------------------------------------------------------------
! Use a Matlab-like call for random numbers
!-------------------------------------------------------------------------------
    call random_number(rand)
end function rand

end module MonteCarloSlabRoutines
