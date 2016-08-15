program template
  use arrayfire
  implicit none

  type(array) A1, A2, tmp
  integer :: I1(2), I2(2)

  ! 1D indexing
  A1 = randu(5, 1)
  A2 = constant(0, 5, 1)
  tmp = get(A1, [3,5])     ! Get elements 3 through 5
  call set(A2, tmp, [1,3]) ! Set elements 1 through 3 with values from tmp
  call print(A1,"A1")
  call print(tmp, "tmp")
  call print(A2,"A2")

  ! 2D indexing
  A1 = randu(3,3)
  A2 = constant(1,3,3)
  I1 = [1, 3]
  I2 = [2, 3]
  tmp = get(A1, idx(I1),  idx(I2))    ! Get rows 1 and 3 for columns 2 and 3
  call set(A2, tmp, idx(I2), idx(I1)) ! Set rows 2 and 3 for columns 1 and 3 with values from tmp
  call print(A1, "A1")
  call print(tmp, "tmp")
  call print(A2, "A2")

  ! 3D indexing
  A1 = randu(3,3,2)
  A2 = constant(1,3,3,2)
  I1 = [1, 3]
  I2 = [2, 3]
  tmp = get(A1, idx(I1),  [1,3,2], [1])  ! Get rows 1 and 3 for columns 1 and 3, tile 1
  call set(A2, tmp, idx(I2), [1,2], [2]) ! Set rows 2 and 3 for columns 1 and 2, tile 2 with tmp
  call print(A1, "A1")
  call print(tmp, "tmp")
  call print(A2, "A2")

end program template
