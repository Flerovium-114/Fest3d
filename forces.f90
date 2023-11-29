module forces

#include "error.h"
#include "debug.h"
use vartypes
use gradients
use geometry
use mpi
use viscosity, only : mu,mu_t
use face_interpolant

private

real(wp),public :: dudx,dudy,dudz
real(wp),public :: dvdx,dvdy,dvdz
real(wp),public :: dwdx,dwdy,dwdz
real(wp),public :: mu_l
real(wp),public :: cd,cl
real(wp),public :: cd_ovl,cl_ovl
real(wp),public :: cp,cpx,cpy,cpz
!real(wp),public :: Ifaces,Jfaces,Kfaces,f_qp_left,f_qp_right,flow
real(wp),public :: imin,imax,jmin,jmax,kmin,kmax,k
real(wp),dimension(:),allocatable,public :: buffer1

public :: main_forces

contains

subroutine  main_forces(qp,cells,dims,bc,control,flow,Ifaces,Jfaces,Kfaces,x_qp_left,x_qp_right,y_qp_left,y_qp_right,z_qp_left,z_qp_right,nodes)

implicit none

integer,dimension(3) :: flags
type(extent) :: dims
real(wp), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2,1:dims%n_var), intent(inout):: qp
type(celltype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: cells
type(boundarytype),intent(in) :: bc
type(controltype), intent(in) :: control
type(flowtype), intent(in) :: flow
type(facetype), dimension(-2:dims%imx+3,-2:dims%jmx+2,-2:dims%kmx+2), intent(in) :: Ifaces
!< Store face quantites for I faces 
            type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes
type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+3,-2:dims%kmx+2), intent(in) :: Jfaces
!< Store face quantites for J faces 
type(facetype), dimension(-2:dims%imx+2,-2:dims%jmx+2,-2:dims%kmx+3), intent(in) :: Kfaces
!< Store face quantites for K faces 
real(wp), dimension(0:dims%imx+1,1:dims%jmx-1,1:dims%kmx-1,1:dims%n_var), intent(inout) :: x_qp_left, x_qp_right
!< Store primitive state at the I-face 
real(wp), dimension(1:dims%imx-1,0:dims%jmx+1,1:dims%kmx-1,1:dims%n_var), intent(inout) :: y_qp_left, y_qp_right
!< Store primitive state at the J-face 
real(wp), dimension(1:dims%imx-1,1:dims%jmx-1,0:dims%kmx+1,1:dims%n_var), intent(inout) :: z_qp_left, z_qp_right

cd = 0.0
cl = 0.0
cl_ovl=0.
cd_ovl=0.

flags = (/1,0,0/) ! I faces
! I MIN 
if (bc%imin_id .eq. -005) then
   call calculate_forces(1,1,1,dims%jmx-1,1,dims%kmx-1,Ifaces,x_qp_left,x_qp_right,dims,flow,'x',flags,nodes,control)  
end if
! I MAX
if (bc%imax_id .eq. -005) then
   call calculate_forces(dims%imx,dims%imx,1,dims%jmx-1,1,dims%kmx-1,Ifaces,x_qp_left,x_qp_right,dims,flow,'x',flags,nodes,control)  
end if

! Jfaces
flags = (/0,1,0/)
!JMIN
if (bc%jmin_id .eq. -005) then
    call calculate_forces(1,dims%imx-1,1,1,1,dims%kmx-1,Jfaces,y_qp_left,y_qp_right,dims,flow,'y',flags,nodes,control)  
end if
! JMAX
if (bc%jmax_id .eq. -005) then
   call calculate_forces(1,dims%imx-1,dims%jmx,dims%jmx,1,dims%kmx-1,Jfaces,y_qp_left,y_qp_right,dims,flow,'y',flags,nodes,control)  
end if

if (dims%kmx .gt. 2) then
   flags = (/0,0,1/)
   ! KMIN
   if (bc%kmin_id .eq. -005) then
           call calculate_forces(1,dims%imx-1,1,dims%jmx-1,1,1,Kfaces,z_qp_left,z_qp_right,dims,flow,'z',flags,nodes,control)  
   end if
   ! KMAX
   if (bc%kmax_id .eq. -005) then
       call calculate_forces(1,dims%imx-1,1,dims%jmx-1,dims%kmx,dims%kmx,Kfaces,z_qp_left,z_qp_right,dims,flow,'z',flags,nodes,control)  
   end if
end if

call gather_and_sum_coefficients_from_all_blocks(control)

end subroutine main_forces


subroutine calculate_forces(imin,imax,jmin,jmax,kmin,kmax,faces,f_qp_left,f_qp_right,dims,flow,ftype,flags,nodes,control)        

implicit none

integer :: i,j,k
integer,intent(in) :: imin,imax,jmin,jmax,kmin,kmax
character,intent(in) :: ftype
integer,dimension(3),intent(in) :: flags
type(extent) :: dims
type(flowtype), intent(in) :: flow
type(facetype), dimension(-2:dims%imx+2+flags(1),-2:dims%jmx+2+flags(2),-2:dims%kmx+2+flags(3)), intent(in) :: faces
real(wp), dimension(1-flags(1):dims%imx-1+2*flags(1), 1-flags(2):dims%jmx-1+2*flags(2), 1-flags(3):dims%kmx-1+2*flags(3), 1:dims%n_var), intent(in) :: f_qp_left, f_qp_right
type(nodetype), dimension(-2:imx+3,-2:jmx+3,-2:kmx+3), intent(in)  :: nodes
type(controltype), intent(in) :: control
real(wp) :: chord
real(wp),dimension(3) :: nl,nd,ns
real(wp) :: fx,fy,fz,Fn
real(wp) :: Fwall_x,Fwall_y,Fwall_z,Fwall
real(wp) :: cfx,cfy,cfz,cf
real(wp) :: pface,pdiff,dyn_pressure,delv
real(wp) :: txx,tyy,tzz,txy,tyz,txz
real(wp) :: xfc,yfc,zfc
character(len=FILE_NAME_LENGTH):: cp_outfile
character(len=FILE_NAME_LENGTH):: cfy_outfile


chord = 0.2
nd = (/1,0,0/)
nl = (/0,1,0/)
ns = (/0,0,1/)

dyn_pressure = 0.5 * flow%density_inf * flow%vel_mag**2 
write(cp_outfile,'(A,I2.2)') 'cp_process_',process_id
write(cfy_outfile,'(A,I2.2)') 'cfy_process_',process_id
open(71,file=trim(cp_outfile)//trim('.dat'))
open(72,file=trim(cfy_outfile)//trim('.dat'))

do k = kmin,kmax
   do j = jmin,jmax
      do i =imin,imax
         select case(ftype)
             case('x')
                xfc = 0.5*(nodes(i,j,k)%x + nodes(i,j+1,k)%x) ! 2D
                yfc = 0.5*(nodes(i,j,k)%y + nodes(i,j+1,k)%y) ! 2D
                zfc = 0.5*(nodes(i,j,k)%z + nodes(i,j+1,k)%z)
                call approx_facedata('x',i,j,k)
             case('y')
                xfc = 0.5*(nodes(i,j,k)%x + nodes(i+1,j,k)%x)
                yfc = 0.5*(nodes(i,j,k)%y + nodes(i+1,j,k)%y)
                zfc = 0.5*(nodes(i,j,k)%z + nodes(i+1,j,k)%z)
                call approx_facedata('y',i,j,k)
             case('z')
                xfc = 0.25*(nodes(i,j,k)%x + nodes(i+1,j,k)%x + nodes(i,j+1,k)%x + nodes(i+1,j+1,k)%x)
                yfc = 0.25*(nodes(i,j,k)%y + nodes(i+1,j,k)%y + nodes(i,j+1,k)%y + nodes(i+1,j+1,k)%y)
                zfc = 0.25*(nodes(i,j,k)%z + nodes(i+1,j,k)%z + nodes(i,j+1,k)%z + nodes(i+1,j+1,k)%z)
                call approx_facedata('z',i,j,k)
         end select            

         ! Shear Stress Calculation
         delv = dudx + dvdy + dwdz
         txx = mu_l * ((dudx + dudx) - (2/3)*delv) 
         tyy = mu_l * ((dvdy + dvdy) - (2/3)*delv) 
         tzz = mu_l * ((dwdz + dwdz) - (2/3)*delv) 
         txy = mu_l * (dudy + dvdx)
         tyz = mu_l * (dvdz + dwdy)
         txz = mu_l * (dudz + dwdx)

         ! Forces along i,j,k
         fx =  txx * faces(i,j,k)%nx + txy*faces(i,j,k)%ny + txz*faces(i,j,k)%nz 
         fy =  txy * faces(i,j,k)%nx + tyy*faces(i,j,k)%ny + tyz*faces(i,j,k)%nz 
         fz =  txz * faces(i,j,k)%nx + tyz*faces(i,j,k)%ny + tzz*faces(i,j,k)%nz
         Fn = fx * faces(i,j,k)%nx + fy * faces(i,j,k)%ny + fz * faces(i,j,k)%nz 
         Fwall_x = fx - Fn * faces(i,j,k)%nx
         Fwall_y = fy - Fn * faces(i,j,k)%ny
         Fwall_z = fz - Fn * faces(i,j,k)%nz
         Fwall = sqrt(Fwall_x**2 + Fwall_y**2 + Fwall_z**2)

         ! Skin friction coefficient
         cf  = Fwall/dyn_pressure
         cfx = Fwall_x/dyn_pressure
         cfy = Fwall_y/dyn_pressure
         cfz = Fwall_z/dyn_pressure

         pface = 0.5 * (f_qp_left(i,j,k,5) + f_qp_right(i,j,k,5))
         pdiff = pface - flow%pressure_inf  
         cp = pdiff/dyn_pressure
         cpx = -cp * faces(i,j,k)%nx * faces(i,j,k)%A
         cpy = -cp * faces(i,j,k)%ny * faces(i,j,k)%A
         cpz = -cp * faces(i,j,k)%nz * faces(i,j,k)%A

         cd = cd +((fx * nd(1) + fy * nd(2) + fz *nd(3)) * faces(i,j,k)%A)/(dyn_pressure) 
         cd = cd + cpx * nd(1) + cpy * nd(2) + cpz * nd(3)
         cl = cl +((fx * nl(1) + fy * nl(2) + fz *nl(3)) * faces(i,j,k)%A)/(dyn_pressure) 
         cl = cl + cpx * nl(1) + cpy * nl(2) + cpz * nl(3)
     
         write(71,*) xfc,yfc,cp
         write(72,*) xfc,yfc,cf
     !    if (control%current_iter .lt. control%max_iters) then 
     !        write(71,*) xfc,yfc,cp
     !        write(72,*) xfc,yfc,cfy
     !    end if

      end do
   end do
end do
close(71)
close(72)

cd = cd/chord        
cl = cl/chord        
        
end subroutine calculate_forces

subroutine approx_facedata(ftype,i,j,k)

implicit none
integer,intent(in) :: i,j,k
character,intent(in) :: ftype
!real(wp) :: mu_l

select case(ftype)
   case ('x')
        dudx = 0.5*(gradu_x(i-1,j,k) + gradu_x(i,j,k)) 
        dudy = 0.5*(gradu_y(i-1,j,k) + gradu_y(i,j,k)) 
        dudz = 0.5*(gradu_z(i-1,j,k) + gradu_z(i,j,k)) 
        dvdx = 0.5*(gradv_x(i-1,j,k) + gradv_x(i,j,k)) 
        dvdy = 0.5*(gradv_y(i-1,j,k) + gradv_y(i,j,k)) 
        dvdz = 0.5*(gradv_z(i-1,j,k) + gradv_z(i,j,k)) 
        dwdx = 0.5*(gradw_x(i-1,j,k) + gradw_x(i,j,k)) 
        dwdy = 0.5*(gradw_y(i-1,j,k) + gradw_y(i,j,k)) 
        dwdz = 0.5*(gradw_z(i-1,j,k) + gradw_z(i,j,k)) 
        mu_l = 0.5*(mu(i-1,j,k) + mu(i,j,k))
   case('y')
        dudx = 0.5*(gradu_x(i,j-1,k) + gradu_x(i,j,k)) 
        dudy = 0.5*(gradu_y(i,j-1,k) + gradu_y(i,j,k)) 
        dudz = 0.5*(gradu_z(i,j-1,k) + gradu_z(i,j,k)) 
        dvdx = 0.5*(gradv_x(i,j-1,k) + gradv_x(i,j,k)) 
        dvdy = 0.5*(gradv_y(i,j-1,k) + gradv_y(i,j,k)) 
        dvdz = 0.5*(gradv_z(i,j-1,k) + gradv_z(i,j,k)) 
        dwdx = 0.5*(gradw_x(i,j-1,k) + gradw_x(i,j,k)) 
        dwdy = 0.5*(gradw_y(i,j-1,k) + gradw_y(i,j,k)) 
        dwdz = 0.5*(gradw_z(i,j-1,k) + gradw_z(i,j,k)) 
        mu_l = 0.5*(mu(i,j-1,k) + mu(i,j,k))
   case('z')
        dudx = 0.5*(gradu_x(i,j,k-1) + gradu_x(i,j,k)) 
        dudy = 0.5*(gradu_y(i,j,k-1) + gradu_y(i,j,k)) 
        dudz = 0.5*(gradu_z(i,j,k-1) + gradu_z(i,j,k)) 
        dvdx = 0.5*(gradv_x(i,j,k-1) + gradv_x(i,j,k)) 
        dvdy = 0.5*(gradv_y(i,j,k-1) + gradv_y(i,j,k)) 
        dvdz = 0.5*(gradv_z(i,j,k-1) + gradv_z(i,j,k)) 
        dwdx = 0.5*(gradw_x(i,j,k-1) + gradw_x(i,j,k)) 
        dwdy = 0.5*(gradw_y(i,j,k-1) + gradw_y(i,j,k)) 
        dwdz = 0.5*(gradw_z(i,j,k-1) + gradw_z(i,j,k)) 
        mu_l = 0.5*(mu(i,j,k-1) + mu(i,j,k))
end select           

end subroutine approx_facedata        
        
subroutine gather_and_sum_coefficients_from_all_blocks(control)
  !< MPI Communication to gather residual from all processes

  implicit none

  type(controltype), intent(in) :: control
  !< Control parameters: number of variables
  integer :: ierr,i,j

  call MPI_ALLGATHER(cl,1, MPI_DOUBLE_PRECISION, &
  buffer1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
 

  do i=1,control%total_process
      cl_ovl =  cl_ovl +buffer1(i)
  end do

  buffer1(:) = 0.0

  call MPI_ALLGATHER(cd,1, MPI_DOUBLE_PRECISION, &
  buffer1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
 

  do i=1,control%total_process
      cd_ovl =  cd_ovl +buffer1(i)
  end do

end subroutine gather_and_sum_coefficients_from_all_blocks

end module forces       
