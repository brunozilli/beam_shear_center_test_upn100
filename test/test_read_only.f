      program test_read_only
      
      implicit none
      
      character(len=256) :: filename
      integer :: nn, ne, unit, iargc
      integer :: i
      double precision, allocatable :: nodes(:,:)
      integer, allocatable :: elements(:,:)
      
      if (iargc() .lt. 1) then
         write(*,*) 'Usage: test_read_only <meshfile.unv>'
         stop
      end if
      
      call getarg(1, filename)
      
      write(*,*) 'Reading: ', trim(filename)
      
      open(newunit=unit, file=filename, status='old')
      call read_section_mesh_unv(unit, nn, ne, nodes, elements)
      close(unit)
      
      write(*,*) 'nn = ', nn
      write(*,*) 'ne = ', ne
      
      if (nn .gt. 0) then
         write(*,*) 'First 5 nodes:'
         do i = 1, min(5, nn)
            write(*,*) '  Node', i, ': y=', nodes(1,i), ' z=', nodes(2,i)
         end do
      end if
      
      if (ne .gt. 0) then
         write(*,*) 'First 5 elements:'
         do i = 1, min(5, ne)
            write(*,*) '  Elem', i, ': ', elements(1,i), elements(2,i),
     +                 elements(3,i)
         end do
      end if
      
      deallocate(nodes, elements)
      
      end
