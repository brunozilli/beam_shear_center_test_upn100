cat > test/test_all_meshes_safe.f << 'EOF'
c     ------------------------------------------------------------------
c     test_all_meshes_safe.f - Safe version with bounds checking
c     ------------------------------------------------------------------
c     MIT License
c     Copyright (c) 2024 Bruno Zilli & DeepSeek
c     ------------------------------------------------------------------

      program test_all_meshes_safe
      
      implicit none
      
c     .. Parameters ..
      integer, parameter :: max_filename_len = 256
      integer, parameter :: max_meshes = 100
      
c     .. Local variables ..
      character(len=max_filename_len) :: mesh_dir, filename
      character(len=max_filename_len) :: mesh_files(max_meshes)
      integer :: n_meshes, i, unit, ios
      integer :: nn, ne
      double precision, allocatable :: nodes(:,:)
      integer, allocatable :: elements(:,:)
      double precision :: area, y_c, z_c, Iy, Iz, Iyz
      double precision :: y_s, z_s
      
c     ------------------------------------------------------------------
      mesh_dir = 'meshes/'
      
      write(*,*) '========================================='
      write(*,*) '  BEAM SHEAR CENTRE - SAFE MULTI-MESH TEST'
      write(*,*) '========================================='
      write(*,*)
      
c     Get list of mesh files
      call get_unv_files_safe(mesh_dir, mesh_files, max_meshes, 
     +                         n_meshes)
      
      if (n_meshes .eq. 0) then
         write(*,*) 'ERROR: No .unv files found in ', mesh_dir
         stop
      end if
      
      write(*,*) 'Found ', n_meshes, ' mesh file(s):'
      do i = 1, n_meshes
         write(*,*) '  ', i, ': ', trim(mesh_files(i))
      end do
      write(*,*)
      
c     ------------------------------------------------------------------
c     Loop over all mesh files
c     ------------------------------------------------------------------
      do i = 1, n_meshes
      
         filename = trim(mesh_dir) // trim(mesh_files(i))
         
         write(*,*) '-----------------------------------------'
         write(*,*) 'TESTING: ', trim(mesh_files(i))
         write(*,*) '-----------------------------------------'
         
c        Read mesh with error checking
         call read_mesh_safe(filename, nn, ne, nodes, elements)
         
         if (nn .eq. 0 .or. ne .eq. 0) then
            write(*,*) '  SKIPPED: Invalid mesh'
            cycle
         end if
         
         write(*,*) '  Nodes: ', nn
         write(*,*) '  Elements: ', ne
         
c        Compute section properties
         call compute_section_properties(nn, ne, nodes, elements,
     +                                   area, y_c, z_c, Iy, Iz, Iyz)
         
         write(*,*) '  Area: ', area
         write(*,*) '  Centroid (y_c, z_c): ', y_c, ', ', z_c
         write(*,*) '  Iy: ', Iy, ' Iz: ', Iz, ' Iyz: ', Iyz
         
c        Compute shear centre
         call compute_shear_center(nn, ne, nodes, elements,
     +                             Iy, Iz, Iyz, y_s, z_s)
         
         write(*,*) '  Shear centre (y_s, z_s): ', y_s, ', ', z_s
         write(*,*)
         
c        Clean up
         if (allocated(nodes)) deallocate(nodes)
         if (allocated(elements)) deallocate(elements)
         
      end do
      
      write(*,*) '========================================='
      write(*,*) '  TEST COMPLETED'
      write(*,*) '========================================='
      
      stop
      end
      
c     ------------------------------------------------------------------
c     Subroutine: get_unv_files_safe
c     ------------------------------------------------------------------
      
      subroutine get_unv_files_safe(dir, file_list, max_files, n_files)
      
      implicit none
      
      character(len=*), intent(in) :: dir
      character(len=*), intent(out) :: file_list(*)
      integer, intent(in) :: max_files
      integer, intent(out) :: n_files
      
      character(len=256) :: command
      character(len=256) :: temp_file
      character(len=256) :: fullname
      integer :: unit, ios
      
      n_files = 0
      temp_file = 'temp_unv_list.txt'
      
c     Build command to list .unv files
      command = 'ls ' // trim(dir) // '*.unv 2>/dev/null > ' // 
     +          trim(temp_file)
      
c     Execute command
      call system(command)
      
c     Read the temporary file
      open(newunit=unit, file=temp_file, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,*) 'Warning: No .unv files found in ', dir
         return
      end if
      
      do while (n_files .lt. max_files)
         read(unit, '(A)', iostat=ios) fullname
         if (ios .ne. 0) exit
         
c        Extract just the filename
         call strip_path_safe(fullname)
         
         n_files = n_files + 1
         file_list(n_files) = trim(fullname)
      end do
      
      close(unit, status='delete')
      
      return
      end
      
c     ------------------------------------------------------------------
c     Subroutine: strip_path_safe
c     ------------------------------------------------------------------
      
      subroutine strip_path_safe(fullname)
      
      implicit none
      
      character(len=*), intent(inout) :: fullname
      
      integer :: i
      
      do i = len_trim(fullname), 1, -1
         if (fullname(i:i) .eq. '/') then
            fullname = fullname(i+1:)
            return
         end if
      end do
      
      return
      end
      
c     ------------------------------------------------------------------
c     Subroutine: read_mesh_safe
c     ------------------------------------------------------------------
      
      subroutine read_mesh_safe(filename, nn, ne, nodes, elements)
      
      implicit none
      
      character(len=*), intent(in) :: filename
      integer, intent(out) :: nn, ne
      double precision, allocatable, intent(out) :: nodes(:,:)
      integer, allocatable, intent(out) :: elements(:,:)
      
      integer :: unit, ios
      integer :: i, j, n1, n2, n3
      double precision :: x, y
      
      nn = 0
      ne = 0
      
c     Open file
      open(newunit=unit, file=filename, status='old', iostat=ios)
      if (ios .ne. 0) then
         write(*,*) '  ERROR: Cannot open file: ', trim(filename)
         return
      end if
      
c     Read UNV format (simplified)
c     Look for node section (marker 2411) and element section (2412)
      
      do
         read(unit, *, iostat=ios) i
         if (ios .ne. 0) exit
         
         if (i .eq. 2411) then
c           Node section found
            read(unit, *) nn
            allocate(nodes(2, nn))
            do i = 1, nn
               read(unit, *) j, x, y
               nodes(1, i) = x
               nodes(2, i) = y
            end do
         else if (i .eq. 2412) then
c           Element section found
            read(unit, *) ne
            allocate(elements(3, ne))
            do i = 1, ne
               read(unit, *) j, n1, n2, n3
               elements(1, i) = n1
               elements(2, i) = n2
               elements(3, i) = n3
            end do
         end if
      end do
      
      close(unit)
      
      write(*,*) '  Read mesh: nn=', nn, ' ne=', ne
      
      return
      end
EOF
