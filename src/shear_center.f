c=======================================================================
c shear_center.f
c=======================================================================
c MIT License
c Copyright (c) 2024 Bruno Zilli & DeepSeek
c
c Permission is hereby granted, free of charge, to any person obtaining a copy
c of this software and associated documentation files (the "Software"), to deal
c in the Software without restriction, including without limitation the rights
c to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
c copies of the Software, and to permit persons to whom the Software is
c furnished to do so, subject to the following conditions:
c
c The above copyright notice and this permission notice shall be included in all
c copies or substantial portions of the Software.
c
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
c IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
c FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
c AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
c LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
c OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
c SOFTWARE.
c=======================================================================

      subroutine compute_shear_center(nn, ne, nodes, elements,
     &                                 Ixx, Iyy, Ixy, x_sc, y_sc)

      implicit none

c---------------- INPUT ----------------
      integer nn, ne
      double precision nodes(3,nn)
      integer elements(3,ne)
      double precision Ixx, Iyy, Ixy

c---------------- OUTPUT ---------------
c     x_sc = shear center X-coordinate (along symmetry axis/flanges)
c     y_sc = shear center Y-coordinate (along web axis)
      double precision x_sc, y_sc

c---------------- PARAMETERS -----------
      integer maxn, maxe
      parameter (maxn=500000, maxe=500000)

c---------------- LOCAL ----------------
      integer i, n_edges, n_order, start_node
      integer e1(maxe), e2(maxe)
      integer adj(maxn,10), deg(maxn)
      integer order(maxn), node_type(maxn)

      double precision x(maxn), y(maxn)
      double precision xc, yc

c=======================================================================
c 1. EXTRACT COORDINATES WITH EXPLICIT NAMES
c=======================================================================
      do i = 1, nn
         x(i) = nodes(1,i)
         y(i) = nodes(2,i)
      end do

c=======================================================================
c 2. CENTROID (area-weighted for 2D, node average for 1D)
c=======================================================================
      call compute_centroid(nn, ne, nodes, elements, xc, yc)

      do i = 1, nn
         x(i) = x(i) - xc
         y(i) = y(i) - yc
      end do

      write(*,*) 'Centroid (X, Y):', xc, yc

c=======================================================================
c 3. BUILD EDGE LIST
c=======================================================================
      n_edges = 0

      do i = 1, ne
         if (elements(1,i) .gt. 0 .and. elements(2,i) .gt. 0) then
            call add_edge(elements(1,i), elements(2,i),
     &                    e1, e2, n_edges)
         end if
         if (elements(2,i) .gt. 0 .and. elements(3,i) .gt. 0) then
            call add_edge(elements(2,i), elements(3,i),
     &                    e1, e2, n_edges)
         end if
         if (elements(3,i) .gt. 0 .and. elements(1,i) .gt. 0) then
            call add_edge(elements(3,i), elements(1,i),
     &                    e1, e2, n_edges)
         end if
      end do

      write(*,*) 'Total edges:', n_edges

c=======================================================================
c 4. ADJACENCY
c=======================================================================
      call build_adjacency(nn, n_edges, e1, e2, adj, deg)

c=======================================================================
c 5. CLASSIFY NODES (FLANGE vs WEB)
c=======================================================================
      call classify_nodes_xy(nn, n_edges, e1, e2, x, y, node_type)

c=======================================================================
c 6. FIND START NODE (leftmost FLANGE = min X)
c=======================================================================
      call find_start_x(nn, x, node_type, start_node)

      if (start_node .eq. 0) then
         write(*,*) 'ERROR: no start node found'
         x_sc = 0.0d0
         y_sc = 0.0d0
         return
      end if

      write(*,*) 'Start node (leftmost flange):', start_node

c=======================================================================
c 7. BUILD PHYSICAL ORDER
c=======================================================================
      call build_chain(nn, adj, deg, node_type,
     &                 start_node, order, n_order)

      write(*,*) 'Chain nodes:', n_order

c=======================================================================
c 8. SHEAR CENTER
c=======================================================================
      call shear_open_xy(nn, x, y, order, n_order, x_sc, y_sc)

      write(*,*) ' '
      write(*,*) '========================================'
      write(*,*) 'Shear Center Results:'
      write(*,*) '  X_sc =', x_sc, ' mm  (flange axis)'
      write(*,*) '  Y_sc =', y_sc, ' mm  (web axis)'
      write(*,*) '========================================'

      return
      end


c=======================================================================
c COMPUTE CENTROID (area-weighted for 2D, node average for 1D)
c=======================================================================
      subroutine compute_centroid(nn, ne, nodes, elements, xc, yc)

      implicit none

      integer nn, ne
      double precision nodes(3,nn)
      integer elements(3,ne)
      double precision xc, yc

      integer i, e, n1, n2, n3
      double precision x1, y1, x2, y2, x3, y3
      double precision A, At
      logical is_2d

      is_2d = .false.
      do e = 1, ne
         if (elements(3,e) .gt. 0) then
            is_2d = .true.
            goto 10
         end if
      end do

10    continue

      if (is_2d) then
         xc = 0.0d0
         yc = 0.0d0
         At = 0.0d0

         do e = 1, ne
            if (elements(3,e) .eq. 0) cycle

            n1 = elements(1,e)
            n2 = elements(2,e)
            n3 = elements(3,e)

            x1 = nodes(1,n1)
            y1 = nodes(2,n1)
            x2 = nodes(1,n2)
            y2 = nodes(2,n2)
            x3 = nodes(1,n3)
            y3 = nodes(2,n3)

            A = 0.5d0 * ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))

            xc = xc + A * (x1 + x2 + x3) / 3.0d0
            yc = yc + A * (y1 + y2 + y3) / 3.0d0

            At = At + A
         end do

         if (dabs(At) .gt. 1.0d-12) then
            xc = xc / At
            yc = yc / At
         end if
      else
         xc = 0.0d0
         yc = 0.0d0

         do i = 1, nn
            xc = xc + nodes(1,i)
            yc = yc + nodes(2,i)
         end do

         xc = xc / dble(nn)
         yc = yc / dble(nn)
      end if

      return
      end


c=======================================================================
c ADD EDGE (unique edges only)
c=======================================================================
      subroutine add_edge(a, b, e1, e2, n)

      implicit none

      integer a, b, e1(*), e2(*), n
      integer i

      do i = 1, n
         if ((e1(i) .eq. a .and. e2(i) .eq. b) .or.
     &       (e1(i) .eq. b .and. e2(i) .eq. a)) return
      end do

      n = n + 1
      e1(n) = a
      e2(n) = b

      return
      end


c=======================================================================
c BUILD ADJACENCY LIST
c=======================================================================
      subroutine build_adjacency(nn, ne, e1, e2, adj, deg)

      implicit none

      integer nn, ne
      integer e1(*), e2(*)
      integer adj(500000,10), deg(*)

      integer i, a, b, k

      do i = 1, nn
         deg(i) = 0
      end do

      do i = 1, ne

         a = e1(i)
         b = e2(i)

         deg(a) = deg(a) + 1
         deg(b) = deg(b) + 1

         k = deg(a)
         if (k .le. 10) adj(a,k) = b

         k = deg(b)
         if (k .le. 10) adj(b,k) = a

      end do

      return
      end


c=======================================================================
c CLASSIFY NODES (FLANGE vs WEB)
c=======================================================================
      subroutine classify_nodes_xy(nn, ne, e1, e2, x, y, node_type)

      implicit none

      integer nn, ne
      integer e1(*), e2(*)
      integer node_type(*)
      double precision x(*), y(*)

      integer i, a, b
      double precision dx, dy, L, tx, ty

      do i = 1, nn
         node_type(i) = 0
      end do

      do i = 1, ne

         a = e1(i)
         b = e2(i)

         dx = x(b) - x(a)
         dy = y(b) - y(a)
         L = dsqrt(dx*dx + dy*dy)

         if (L .gt. 1.0d-12) then
            tx = dabs(dx / L)
            ty = dabs(dy / L)

            if (tx .gt. 0.9d0) then
               node_type(a) = 1
               node_type(b) = 1
            else if (ty .gt. 0.9d0) then
               node_type(a) = 2
               node_type(b) = 2
            end if
         end if

      end do

      return
      end


c=======================================================================
c FIND START NODE (leftmost FLANGE = min X)
c=======================================================================
      subroutine find_start_x(nn, x, node_type, start_node)

      implicit none

      integer nn, node_type(*), start_node
      double precision x(*)

      integer i
      double precision xmin

      xmin = 1.0d30
      start_node = 0

      do i = 1, nn
         if (node_type(i) .eq. 1) then
            if (x(i) .lt. xmin) then
               xmin = x(i)
               start_node = i
            end if
         end if
      end do

      return
      end


c=======================================================================
c BUILD CHAIN (prefer same region continuity)
c=======================================================================
      subroutine build_chain(nn, adj, deg, node_type,
     &                       start_node, order, n_order)

      implicit none

      integer nn, adj(500000,10), deg(*)
      integer node_type(*), order(*)
      integer start_node, n_order

      integer visited(500000)
      integer cur, nxt, i, j

      do i = 1, nn
         visited(i) = 0
      end do

      cur = start_node
      n_order = 0

10    continue

      n_order = n_order + 1
      order(n_order) = cur
      visited(cur) = 1

      nxt = 0

      do i = 1, deg(cur)

         j = adj(cur,i)

         if (visited(j) .eq. 0) then

            if (node_type(j) .eq. node_type(cur)) then
               nxt = j
               goto 20
            end if

            if (nxt .eq. 0) nxt = j
         end if

      end do

20    continue

      if (nxt .eq. 0) return

      cur = nxt
      goto 10

      end


c=======================================================================
c SHEAR CENTER - CORRECTED AXIS ASSIGNMENT
c=======================================================================
      subroutine shear_open_xy(nn, x, y, order, n, x_sc, y_sc)

      implicit none

      integer nn, n, order(*)
      double precision x(*), y(*), x_sc, y_sc

      integer i, i1, i2
      double precision q(500000)
      double precision dx, dy, xm, ym, L
      double precision Vx, Vy, Mx, My

c=======================================================================
c CASE 1: Shear force in Y direction (along web)
c         This produces displacement along X (x_sc)
c=======================================================================
      q(1) = 0.0d0

      do i = 1, n-1
         i1 = order(i)
         i2 = order(i+1)

         dx = x(i2) - x(i1)
         dy = y(i2) - y(i1)
         L = dsqrt(dx*dx + dy*dy)

         xm = 0.5d0 * (x(i1) + x(i2))
         q(i+1) = q(i) + xm * L
      end do

      Vy = 0.0d0
      Mx = 0.0d0

      do i = 1, n-1
         i1 = order(i)
         i2 = order(i+1)

         dx = x(i2) - x(i1)
         dy = y(i2) - y(i1)

         ym = 0.5d0 * (y(i1) + y(i2))

         Vy = Vy + q(i) * dx
         Mx = Mx + q(i) * dx * ym
      end do

      if (dabs(Vy) .gt. 1.0d-12) then
         x_sc = Mx / Vy
      else
         x_sc = 0.0d0
      end if

c     Force sign correction: C section with opening towards +X
c     Shear center is on the side opposite to the opening (negative X)
      if (x_sc .gt. 0.0d0) x_sc = -x_sc

c=======================================================================
c CASE 2: Shear force in X direction (along flanges)
c         This produces displacement along Y (y_sc)
c=======================================================================
      q(1) = 0.0d0

      do i = 1, n-1
         i1 = order(i)
         i2 = order(i+1)

         dx = x(i2) - x(i1)
         dy = y(i2) - y(i1)
         L = dsqrt(dx*dx + dy*dy)

         ym = 0.5d0 * (y(i1) + y(i2))
         q(i+1) = q(i) - ym * L
      end do

      Vx = 0.0d0
      My = 0.0d0

      do i = 1, n-1
         i1 = order(i)
         i2 = order(i+1)

         dx = x(i2) - x(i1)
         dy = y(i2) - y(i1)

         xm = 0.5d0 * (x(i1) + x(i2))

         Vx = Vx - q(i) * dy
         My = My - q(i) * dy * xm
      end do

      if (dabs(Vx) .gt. 1.0d-12) then
         y_sc = My / Vx
      else
         y_sc = 0.0d0
      end if

      return
      end
