c=======================================================================
c read_section_mesh_unv.f - Fixed to support edge elements (type 11)
c=======================================================================
c Copyright (c) 2024 Bruno Zilli & DeepSeek
c MIT License
c=======================================================================
c C section with opening towards +X, CG centered in plane X-Y,
c X symmetry axis.
c=======================================================================

      subroutine read_section_mesh_unv(filename, max_nodes, 
     &                                 max_elements, n_nodes, 
     &                                 n_elements, conn, x, y)

      implicit none
      
c     INPUT
      character(len=*), intent(in) :: filename
      integer, intent(in) :: max_nodes, max_elements
      
c     OUTPUT
      integer, intent(out) :: n_nodes, n_elements
      integer, intent(out) :: conn(3, max_elements)
      double precision, intent(out) :: x(max_nodes), y(max_nodes)
      
c     LOCALS
      integer :: io, node_id, elem_id, elem_type
      integer :: n1, n2, n3, dummy, phys_prop
      integer :: n_edges, n_triangles
      character(len=256) :: line, coord_line
      logical :: in_nodes, in_elements, salome_format
      double precision :: xc, yc, zc
      
c     Initialise
      n_nodes = 0
      n_elements = 0
      n_edges = 0
      n_triangles = 0
      in_nodes = .false.
      in_elements = .false.
      salome_format = .false.
      
c     Open file
      open(10, file=filename, status='old', iostat=io)
      if (io .ne. 0) then
        print *, 'ERROR: Cannot open file ', filename
        return
      end if
      
c     Scan file for sections
      do
        read(10, '(A)', iostat=io) line
        if (io .ne. 0) exit
        
c       Detect node section (2411)
        if (index(line, '2411') .ne. 0) then
          in_nodes = .true.
          in_elements = .false.
c         Check format: read next line to see if it has 4 integers
          read(10, '(A)', iostat=io) line
          backspace(10)
          read(line, *, iostat=io) dummy, dummy, dummy, dummy
          if (io .eq. 0) then
            salome_format = .true.      ! Salome/ASTER format
          else
            salome_format = .false.     ! Simple format
          endif
          cycle
        endif
        
c       Detect element section (2412)
        if (index(line, '2412') .ne. 0) then
          in_nodes = .false.
          in_elements = .true.
          cycle
        endif
        
c       End of section marker
        if (line .eq. '    -1' .or. line .eq. '-1') then
          in_nodes = .false.
          in_elements = .false.
          cycle
        endif
        
c       Read nodes
        if (in_nodes) then
          if (salome_format) then
c           Salome format: "id 1 1 11" on first line, coordinates on next
            read(line, *, iostat=io) node_id, dummy, dummy, dummy
            if (io .eq. 0 .and. node_id .gt. 0) then
c             Read coordinate line as string
              read(10, '(A)', iostat=io) coord_line
              if (io .ne. 0) exit
c             Check for end of node section
              if (coord_line .eq. '    -1' .or. 
     &            coord_line .eq. '-1') then
                in_nodes = .false.
                cycle
              endif
c             Parse coordinates (UNV format: X, Y, Z)
              read(coord_line, *, iostat=io) xc, yc, zc
              if (io .eq. 0) then
                n_nodes = n_nodes + 1
                if (n_nodes .le. max_nodes) then
                  x(n_nodes) = xc
                  y(n_nodes) = yc
                endif
              endif
            endif
          else
c           Simple format: "id X Y Z" all on one line
            read(line, *, iostat=io) node_id, xc, yc, zc
            if (io .eq. 0 .and. node_id .gt. 0) then
              n_nodes = n_nodes + 1
              if (n_nodes .le. max_nodes) then
                x(n_nodes) = xc
                y(n_nodes) = yc
              endif
            endif
          endif
          cycle
        endif
        
c       Read elements (supports type 41 = triangle, type 11 = edge)
        if (in_elements) then
c         Element header: "id type phys_prop ..."
          read(line, *, iostat=io) elem_id, elem_type, phys_prop
          if (io .eq. 0) then
            if (elem_type .eq. 41) then
c             Triangle element (3 nodes)
              n_elements = n_elements + 1
              n_triangles = n_triangles + 1
              if (n_elements .le. max_elements) then
c               Read connectivity: n1 n2 n3
                read(10, *, iostat=io) n1, n2, n3
                if (io .eq. 0) then
                  conn(1, n_elements) = n1
                  conn(2, n_elements) = n2
                  conn(3, n_elements) = n3
                endif
              endif
            else if (elem_type .eq. 11) then
c             Edge element (2 nodes) - SALOME format
              n_elements = n_elements + 1
              n_edges = n_edges + 1
              if (n_elements .le. max_elements) then
c               Skip the "0 1 1" line (3 integers)
                read(10, *, iostat=io) dummy, dummy, dummy
c               Read connectivity on next line
                read(10, *, iostat=io) n1, n2
                if (io .eq. 0) then
                  conn(1, n_elements) = n1
                  conn(2, n_elements) = n2
                  conn(3, n_elements) = 0
                endif
              endif
            endif
          endif
        endif
      end do
      
      close(10)
      
      print *, 'Read mesh from: ', filename
      print *, '   Nodes:       ', n_nodes
      print *, '   Elements:    ', n_elements
      print *, '     Triangles: ', n_triangles
      print *, '     Edges:     ', n_edges
      
      end subroutine read_section_mesh_unv
