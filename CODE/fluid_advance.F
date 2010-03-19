      subroutine fluid_advance(i)
      use fluid_mod
      implicit  none
      integer :: i
      real(DOUBLE) :: k1r,k2r,k3r,k4r
      real(DOUBLE) :: k1phi,k2phi,k3phi,k4phi
      real(DOUBLE) :: k1z,k2z,k3z,k4z
      real(DOUBLE) :: dvr,domega,dvz,ro,zo,phio
      ro=r_e(TSTEP,i)
      phio=angle_e(TSTEP,i)
      zo=z_e(TSTEP,i)
      k1r=vr_e(TSTEP,i)*delta_t
      k1phi=omega_e(TSTEP,i)*delta_t
      k1z=vz_e(TSTEP,i)*delta_t
      r_e(TSTEP,i)=ro+k1r*0.5d0
      angle_e(TSTEP,i)=phio+k1phi*0.5d0
      z_e(TSTEP,i)=zo+k1z*0.5d0
      if (angle_e(TSTEP,i)>2.d0*pi) then
        angle_e(TSTEP,i)=angle_e(TSTEP,i)-2.d0*pi
      ENDIF
      call fluid_interp(i,dvr,domega,dvz,.false.)
      k2r=dvr*delta_t
      k2phi=domega*delta_t
      k2z=dvz*delta_t
      r_e(TSTEP,i)=ro+k2r*0.5d0
      angle_e(TSTEP,i)=phio+k2phi*0.5d0
      z_e(TSTEP,i)=zo+k2z*0.5d0
      if (angle_e(TSTEP,i)>2.d0*pi) then
        angle_e(TSTEP,i)=angle_e(TSTEP,i)-2.d0*pi
      ENDIF
      call fluid_interp(i,dvr,domega,dvz,.false.)
      k3r=dvr*delta_t
      k3phi=domega*delta_t
      k3z=dvz*delta_t
      r_e(TSTEP,i)=ro+k3r
      angle_e(TSTEP,i)=phio+k3phi
      z_e(TSTEP,i)=zo+k3z
      if (angle_e(TSTEP,i)>2.d0*pi) then
        angle_e(TSTEP,i)=angle_e(TSTEP,i)-2.d0*pi
      ENDIF
      call fluid_interp(i,dvr,domega,dvz,.false.)
      k4r=dvr*delta_t
      k4phi=domega*delta_t
      k4z=dvz*delta_t
      r_e(TSTEP,i)=ro
      angle_e(TSTEP,i)=phio
      z_e(TSTEP,i)=zo
      r_e(TSTEP+1,i)=ro+(k1r+k4r)/6.d0+(k2r+k3r)/3.d0
      angle_e(TSTEP+1,i)=phio+(k1phi+k4phi)/6.d0+(k2phi+k3phi)/3.d0
      z_e(TSTEP+1,i)=zo+(k1z+k4z)/6.d0+(k2z+k3z)/3.d0
      if (angle_e(TSTEP+1,i)>2.d0*pi) then
        angle_e(TSTEP+1,i)=angle_e(TSTEP+1,i)-2.d0*pi
      ENDIF
      return
      endsubroutine
