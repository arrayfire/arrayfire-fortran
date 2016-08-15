program dla
  use arrayfire
  implicit none

  type(array) A, B, C, X0, X1, X2
  type(array) l, u, p, r, diff
  integer :: N

  ! Set the size of the matrix
  N = 5

  ! Randomly generate a system
  A = randu(N, N)
  X0 = randu(N, 1)

  ! LU decomposition
  call lu(l, u, p, A)

  ! Construct a positive definite matrix
  C = A + transpose(A) + identity(N, N) * 100.0
  call cholesky(r, C)

  ! Solve a general system of equations
  B = matmul(A, X0)
  X1 = solve(A, B)
  call print(max(abs(X0 - X1)), "absolute error: Solving general system")

  ! Solve a positive definite system of equations
  B = matmul(C, X0)
  X2 = solve(C, B)
  call print(max(abs(X0 - X2)), "absolute error: Solving positive definite system")

  ! Invert a matrix
  r = inverse(A)
  diff = abs(identity(N, N) - matmul(r, A))
  call print(max(moddims(diff, N * N, 1)), "absolute error: Inverting a matrix")

end program dla
