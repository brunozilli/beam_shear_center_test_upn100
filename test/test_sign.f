c     Quick test to check sign convention
      program test_sign
      implicit none
      double precision :: M_y, M_z, y_s, z_s
      double precision :: Iy, Iz, Iyz, detI
      
c     Dalla tua L-section:
      Iy = 1764444.2544620314
      Iz = 1764444.2544620261
      Iyz = -1035288.8769560433
      M_y = 68562244.724935561
      M_z = -68563444.054936543
      
      detI = Iy*Iz - Iyz*Iyz
      
      write(*,*) 'Original:'
      y_s = (Iz * M_y - Iyz * M_z) / detI
      z_s = (Iy * M_z - Iyz * M_y) / detI
      write(*,*) '  y_s =', y_s
      write(*,*) '  z_s =', z_s
      
      write(*,*) 'With M_z sign flipped (+):'
      y_s = (Iz * M_y - Iyz * (-M_z)) / detI
      z_s = (Iy * (-M_z) - Iyz * M_y) / detI
      write(*,*) '  y_s =', y_s
      write(*,*) '  z_s =', z_s
      
      write(*,*) 'With M_y sign flipped (+):'
      y_s = (Iz * (-M_y) - Iyz * M_z) / detI
      z_s = (Iy * M_z - Iyz * (-M_y)) / detI
      write(*,*) '  y_s =', y_s
      write(*,*) '  z_s =', z_s
      end
