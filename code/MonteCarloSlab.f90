program MonteCarloSlab
!-------------------------------------------------------------------------------
! MonteCarloSlab -- A Fortran Monte Carlo Code for Neutron Transport in Slabs
!-------------------------------------------------------------------------------
    use MonteCarloSlabRoutines
    implicit none
    type(Input)     :: Inpt  ! Problem input
    type(Output)    :: Outp ! Problem output
    integer         :: i
    !---------------------------------------------------------------------------
    ! SAMPLE PROBLEM:  (A crafty student will add input file handling!)
    !     vacuum | 4 cm, mat2 | 4 cm, mat1, source=1 | 4 cm, mat2 | vacuum
    ! We divide the domain into 24 cells of with 0.5 cm:
    Inpt%NumCells = 24
    ! Allocate the cell boundary and material id arrays
    allocate( Inpt%CellBounds(Inpt%NumCells+1), Inpt%MatId(Inpt%NumCells) )
    ! Set the boundaries and assign materials
    Inpt%CellBounds = 0.0
    do i = 2, Inpt%NumCells+1
        Inpt%CellBounds(i) = Inpt%CellBounds(i-1) + 0.5
    end do
    Inpt%MatId = (/2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2/)
    ! We have two materials, each with a SigmaT, SigmaA, and SigmaS
    Inpt%NumMats = 2
    allocate( Inpt%Sigma(Inpt%NumMats,3) ) 
    ! Sigma =     material 1  ( SigmaT | SigmaA | SigmaS )
    !             material 2  ( ...
    Inpt%Sigma(1,1) = 1.0d0; Inpt%Sigma(1,2) = 0.5d0; Inpt%Sigma(1,3) = 0.5d0
    Inpt%Sigma(2,1) = 1.5d0; Inpt%Sigma(2,2) = 1.2d0; Inpt%Sigma(2,3) = 0.3d0
    ! Let's do a million histories
    Inpt%NumHist = 1e6;
    ! Allocate the tallies
    allocate( Outp%TrakEst(Inpt%NumCells), Outp%CollEst(Inpt%NumCells) )
    call PlayGame(Inpt,Outp)
	!**** INSERT TALLY POST-PROCESSING HERE (i.e. flux computation, etc.)
    ! Deallocate input
    deallocate( Inpt%CellBounds, Inpt%MatId, Inpt%Sigma )
    ! Deallocate output
    deallocate( out%TrakEst, out%CollEst )
end program MonteCarloSlab






