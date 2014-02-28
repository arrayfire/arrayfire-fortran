program basic
  use arrayfire
  implicit none

  real, dimension(3,3) :: a
  real, dimension(:,:), allocatable :: b
  type(array) M1, M2, tmp

  a(1,:) = (/ 1, 0, 0 /)
  a(2,:) = (/ 0, 3, 0 /)
  a(3,:) = (/ 0, 0, 2 /)

  ! Copy data from host to device
  M1 = a
  write(*,*) "Showing the matrix after mem copy from host (M1)"
  call print(M1)

  ! Generate a uniformly random, single precision matrix
  M2 = randu(3,3, ty=f32)
  write(*,*) "Showing a randomly generated matrix (M2)"
  call print(M2)

  ! Transpose of matrix
  tmp = transpose(M2)
  call print(tmp, "Transpose of M2") ! Displays array after printing message

  ! Element wise addition
  tmp = M1 + M2
  call print(tmp, "M1 + M2")

  ! element wise subtraction
  call print(-M2, "-M2")

  ! Trignometric functions
  write(*,*) "Displaying sin(M2)**2 + cos(M2)**2"
  call print(sin(M2)**2.0 + cos(M2)**2.0)

  ! Multiplication of matrices
  ! Matrix multiply
  write(*, *) "Matrix multiply: matmul(M1, M2)"
  call print(matmul(M1, M2))

  ! Element wise multiplication
  write(*, *) "Element wise multiplication: M1 * M2"
  call print(M1 * M2)

  ! minimum value
  tmp = min(M2)
  call print(tmp, "min(M2)")

  ! Get back to host
  b = tmp
  write(*,*) "Showing min(M2) data back on host "
  write(*,*) b(:,:)

end program basic
