program timer
  use arrayfire
  implicit none
  double precision elapsed
  type(array) A, B

  ! Generate a random matrix
  A = randu(1024, 2048)

  ! Start the timer
  call timer_start()

  ! Perform operations here.
  B = matmul(A, transpose(A))

  ! Stop the timer
  elapsed = timer_stop()

  ! Print time taken
  write (*,"(a15, d8.2)") "Time taken: ", elapsed
end program timer
