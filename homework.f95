module homework
contains
subroutine FindMaxCoordinates(a, x1, y1, x2, y2)
implicit none
include "mpif.h"
real(8), dimension(:,:), intent(out) :: a
integer(4), intent(out) :: x1, y1, x2, y2
integer(4) :: mpiErr, mpiRank, mpisize
call mpi_comm_size(MPI_COMM_WORLD, mpisize, mpiErr)
if(mpisize == 1) then
    call helper(a, x1, y1, x2, y2)
else
  call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
  if(mpiRank == 0) then
      call master(a, x1, y1, x2, y2)
  else
      call helper(a, x1, y1, x2, y2)
  endif
endif
end subroutine

subroutine master(a, x1, y1, x2, y2)
implicit none
include "mpif.h"
real(8), dimension(:,:) :: a
integer(4), intent(out) :: x1, y1, x2, y2
real(8), dimension(5) :: res
integer(4) :: lines, columns, i, nextline, usedline, activeHelpers
integer(4) :: mpiErr, mpisize, mpiRank
integer(4), dimension(MPI_sTATUs_sIZE) :: status
real(8) :: now, maxsum
integer(4), dimension(4) :: answer
lines = size(a, dim = 1)
columns = size(a, dim = 2)
call mpi_comm_size(MPI_COMM_WORLD, mpisize, mpiErr)
call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)
activeHelpers = mpisize - 1
maxsum = a(1,1)
x1 = 1
y1 = 1
x2 = 1
y2 = 1
nextline = 1
usedline = 0
do while (usedline < lines .or. activeHelpers > 0)
    call mpi_recv(res, 5, MPI_REAL8, MPI_ANY_sOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpiErr)
    if (status(MPI_TAG) == 4) then
        usedline = usedline + 1        
        now = res(1)        
        if (now > maxsum) then
            maxsum = now
            x1 = int(res(2))
            y1 = int(res(3))
            x2 = int(res(4))
            y2 = int(res(5))
        end if           
    end if
    if (nextline <= lines) then
        call mpi_send(nextline, 1, MPI_INTEGER4, status(MPI_sOURCE), 5, MPI_COMM_WORLD, mpiErr)
        nextline = nextline + 1
    else
        call mpi_send(nextline, 1, MPI_INTEGER4, status(MPI_sOURCE), 6, MPI_COMM_WORLD, mpiErr)
        activeHelpers = activeHelpers - 1
    end if
end do
answer(1) = x1
answer(2) = y1
answer(3) = x2
answer(4) = y2
do i=1,mpisize-1
    call mpi_send(answer, 4, MPI_INTEGER4, i, 1, MPI_COMM_WORLD, mpiErr)
end do
end subroutine

subroutine helper(a, x1, y1, x2, y2)
implicit none
include "mpif.h"
integer(4) :: mpiErr, mpisize, mpiRank
integer(4), dimension(MPI_sTATUs_sIZE) :: status
real(8), dimension(:,:) :: a
real(8), dimension(:), allocatable :: p
real(8), dimension(5) :: res
integer(4), intent(out) :: x1, y1, x2, y2
integer(4) :: columns, lines, up, down, left, right, j
integer(4), dimension(4) :: answer
real(8) :: now, maxsum
lines = size(a, dim = 1)
columns = size(a, dim = 2)
allocate(p(columns))
call mpi_send(now, 1, MPI_REAL8, 0, 3, MPI_COMM_WORLD, mpiErr)
do
    call mpi_recv(up, 1, MPI_INTEGER4, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, mpiErr)
    if (status(MPI_TAG) == 6) then
        exit
    end if
    p = 0
    x1 = up
    do down = up,lines
        do j = 1,columns
            p(j) = p(j) + a(down, j)
        enddo
        call kadane(p, left, right, now)
        if (now > maxsum .or. up == down) then
            maxsum = now;
            y1 = left
            y2 = right
            x2 = down
        endif
    enddo
    res(1) = maxsum
    res(2) = x1
    res(3) = y1
    res(4) = x2
    res(5) = y2
    call mpi_send(res, 5, MPI_REAL8, 0, 4, MPI_COMM_WORLD, mpiErr)
end do
deallocate(p)
call mpi_recv(answer, 4, MPI_INTEGER4, 0, 1, MPI_COMM_WORLD, status, mpiErr)
x1 = answer(1)
y1 = answer(2)
x2 = answer(3)
y2 = answer(4)
end subroutine

subroutine kadane(a, x1, x2, maxsum)
implicit none
real(8), intent(in), dimension(:) :: a
integer(4), intent(out) :: x1, x2
real(8), intent(out) :: maxsum
integer(4) :: i, ind, n
real(8) :: maximum, sum1, sum2
n = size(a)
maxsum = a(1); x1 = 1; x2 = 1
maximum=a(1); ind=1
do i=2,n
  sum1 = a(i)
  sum2 = maximum + a(i)
  if (sum1 > sum2) then
    maximum = sum1
    ind = i
  else
    maximum = sum2
  endif
  if (maximum > maxsum) then
    maxsum = maximum
    x1 = ind
    x2 = i
  endif
enddo
end subroutine

end module
