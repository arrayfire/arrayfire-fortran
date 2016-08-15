module arrayfire
  use, intrinsic :: ISO_C_Binding, only: C_ptr, C_NULL_ptr
  implicit none

  !> Contains the last known error in the arrayfire module
  integer :: err

  !> Single precision, real  type
  integer :: f32 = 1
  !> Single precision, complex type
  integer :: c32 = 2
  !> Double precision, real type
  integer :: f64 = 3
  !> Double precision, complex type
  integer :: c64 = 4
  !> Boolean type
  integer :: b8 = 5

  !> type(array) containing information about device
  type array
     !> Dimensions of array
     integer :: shape(4)
     !> Rank of the array
     integer :: rank
     !> Device pointer
     type(C_ptr) :: ptr = C_NULL_ptr
     ! not supported in gfortran 4.3
     ! contains
     !   procedure :: get => array_get
     !   procedure :: set => array_set
  end type array

  !> @defgroup basic Basics
  !! @{

  !> @defgroup mem Create arrays from host data
  !> @{
  !> Memory transfer from host to device, device to host.
  !> @param[in] rhs -- Can be type array, real, double precision, complex, double complex
  !> @returns lhs after assigning rhs to lhs
  !> @code
  !! type(array) arr
  !! real, dimension(2, 2) :: a
  !! a(1,1) = 1
  !! a(2,2) = 1
  !! ! a is now identity matrix of width 2
  !! arr = a         ! Data is now on device
  !! arr = arr + 1.0 ! Increment arr by 1
  !! a = log(arr)    ! Return log(arr) to host
  !! @endcode
  interface assignment (=)
     module procedure device1_s, device1_d, device1_c, device1_z
     module procedure device2_s, device2_d, device2_c, device2_z
     module procedure device3_s, device3_d, device3_c, device3_z
     module procedure device4_s, device4_d, device4_c, device4_z
     module procedure assign
     module procedure host1_s, host1_d, host1_c, host1_z
     module procedure host2_s, host2_d, host2_c, host2_z
     module procedure host3_s, host3_d, host3_c, host3_z
     module procedure host4_s, host4_d, host4_c, host4_z
  end interface assignment (=)

  interface getptr
     module procedure hostp1_s, hostp1_d, hostp1_c, hostp1_z
     module procedure hostp2_s, hostp2_d, hostp2_c, hostp2_z
     module procedure hostp3_s, hostp3_d, hostp3_c, hostp3_z
     module procedure hostp4_s, hostp4_d, hostp4_c, hostp4_z
  end interface getptr

  !> @}


  !> @defgroup gen Generate random or constant matrices
  !> Matrix generation
  !> @code
  !! type(array) res
  !! res = randu(3, 3)          ! Uniformly distributed random matrix, single precision by default
  !! res = randn(3, 3)          ! Normally distributed  random matrix, single precision by default
  !! res = constant(1,5, 5, ty = f32) ! Single precision matrix of all ones
  !! res = constant(0,4, 4,ty = f64) ! Double precision matrix of all zeros
  !! res = identity(2, 2)       ! Identity Matrix, single precision by default
  !! @endcode
  !> @{
  !> Generate matrices on the devices

  !> @{
  !> Random, uniformly distributed
  !> @param[in] x1 -- The 1st dimension in the array (should be integer)
  !> @param[in] x2 -- The 2nd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x3 -- The 3rd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x4 -- The 4th dimension in the array (should be integer, optional, default: 1)
  !> @param[in] ty -- Should be one of (f32, f64, c32, c64), (optional, default: f32)
  !> @returns output of size (x1, x2, x3, x4) filled with the required data type
  interface randu
     module procedure array_randu
  end interface randu
  !> @}

  !> @{
  !> Random, normally distributed
  !> @param[in] x1 -- The 1st dimension in the array (should be integer)
  !> @param[in] x2 -- The 2nd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x3 -- The 3rd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x4 -- The 4th dimension in the array (should be integer, optional, default: 1)
  !> @param[in] ty -- Should be one of (f32, f64, c32, c64), (optional, default: f32)
  !> @returns output of size (x1, x2, x3, x4) filled with the required data type
  interface randn
     module procedure array_randn
  end interface randn
  !> @}

  !> @{
  !> Constant value
  !> @param[in] val -- constant value to constant with
  !> @param[in] x1 -- The 1st dimension in the array (should be integer)
  !> @param[in] x2 -- The 2nd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x3 -- The 3rd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x4 -- The 4th dimension in the array (should be integer, optional, default: 1)
  !> @param[in] ty -- Should be one of (f32, f64, c32, c64), (optional, default: f32)
  !> @returns output of size (x1, x2, x3, x4) constanted with the required data type
  interface constant
     module procedure array_constant
  end interface constant
  !> @}

  !> @{
  !> Constant, identity
  !> @param[in] x1 -- The 1st dimension in the array (should be integer)
  !> @param[in] x2 -- The 2nd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] ty -- Should be one of (f32, f64, c32, c64), (optional, default: f32)
  !> @returns output of size (x1, x2, x3, x4) filled with the required data type
  interface identity
     module procedure array_identity
  end interface identity
  !> @}
  !> @}

  !> @defgroup indexing Array Indexing
  !> @{
  !> Setting and getting values from an array
  !> @code
  !> type(array) r, s, z
  !> r = randu(5,5)
  !> z = constant(0,5,5)
  !> s = get(a, idx(2), seq(1,5,2))   !Get the second row and every other column between 1 and 5
  !> call set(z, s, idx(1), seq(1,3)) !Set s as the first row and columsn 1 through 3
  !> @endcode

  !> @{
  !> @param[in] in Input array
  !> @param[in] d1 type(array) denoting indices along 1st dimension
  !> @param[in] d2 type(array) denoting indices along 2nd dimension. Optional.
  !> @param[in] d3 integer denoting the index of the 3rd dimension. Optional.
  !> @param[in] d4 integer denoting the index of the 4th dimension. Optional.
  !> @returns subarry of in referenced by d1,d2,d3,d4
  interface get
     module procedure array_get, array_get2, array_get_seq
  end interface get
  !> @}

  !> @{
  !> @param[in] lhs Array whos values are being set by rhs
  !> @param[in] rhs The value being set
  !> @param[in] d1 type(array) denoting indices along 1st dimension
  !> @param[in] d2 type(array) denoting indices along 2nd dimension. Optional.
  !> @param[in] d3 integer denoting the index of the 3rd dimension. Optional.
  !> @param[in] d4 integer denoting the index of the 4th dimension. Optional.
  interface set
     module procedure array_set, array_set2, array_set_seq
  end interface set
  !> @}

  !> @{
  !> @param[in] index Can be an integer scalar or array
  !> @returns type(array) holding the value of index
  interface idx
     module procedure idx_scalar, idx_vector
  end interface idx
  !> @}

  !> @{
  !> @param[in] first The first element of the sequence. Optional. Default: 0.
  !> @param[in] last  The last  element of the sequence. Optional. Default: 0.
  !> @param[in] setp  The step size. Optional. Default: 1.
  interface seq
     module procedure arr_seq
  end interface seq
  !> @}
  !> @}


  !> @defgroup tile Tiling and reshaping
  !> @{
  !> moddims, flat, tile, etc

  !> @{
  !> Moddims an array
  !> @param[in] A an array
  !> @param[in] x1 -- The 1st dimension in the array (should be integer)
  !> @param[in] x2 -- The 2nd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x3 -- The 3rd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x4 -- The 4th dimension in the array (should be integer, optional, default: 1)
  !> @returns output of size (x1, x2, x3, x4) with the same data as input
  !> @code
  !! type(array) A, B
  !! A = randu(2, 4)          ! 2 dimensional array
  !! B = moddims(A, 2, 2, 2)  ! 3 dimensional array
  !! @endcode
  interface moddims
     module procedure array_moddims
  end interface moddims
  !> @}

  !> @{
  !> Flatten an array
  !> @param[in] A an array
  !> @returns output of size with the same data as input, but as a column vector
  !> @code
  !! type(array) A, B
  !! A = randu(2, 4) ! 2 dimensional array
  !! B = flat(A)     ! 1 dimensional array
  !! @endcode
  interface flat
     module procedure array_flat
  end interface flat
  !> @}

  !> @{
  !> Tile an array
  !> @param[in] A an array
  !> @param[in] x1 -- The 1st dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x2 -- The 2nd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x3 -- The 3rd dimension in the array (should be integer, optional, default: 1)
  !> @param[in] x4 -- The 4th dimension in the array (should be integer, optional, default: 1)
  !> @returns output of size (x1, x2, x3, x4) with the same data as input
  !> @code
  !! type(array) A, B
  !! A = randu(2, 4)       ! 2 dimensional array
  !! B = tile(A, 1, 1, 2)  ! 3 dimensional array
  !! @endcode
  interface tile
     module procedure array_tile
  end interface tile
  !> @}

  !> @{
  !> Join two arrays
  !> @param[in] d -- The dimension along which two arrays are to be joined
  !> @param[in] first -- The first array
  !> @param[in] second -- The second array
  !> @returns output which is the combined array of first and second
  !> @code
  !! type(array) A, B, C
  !! A = randu(2, 4)       ! 2 dimensional array
  !! B = randu(3, 4)       ! 2 dimensional array
  !! C = join(1, A, B)     ! 3 dimensional array. Size: 5, 4
  !! @endcode
  interface join
     module procedure array_join
  end interface join
  !> @}
  !> @}


  !> @defgroup manip Reorder and sorting: (transpose, reorder, sort)
  !> @{
  !> Functions that perform data re-ordering and sorting

  !> @{
  !> Matrix transpose
  !> @param[in] A -- type array of size M x N
  !> @returns B of size N x M which is the matrix transpose of A
  !> @code
  !! type(array) A, B
  !! A = randu(20, 25) ! Random matrix
  !! B = transpose(A)  ! Matrix transpose
  !! @endcode
  interface transpose
     module procedure array_transpose
  end interface transpose
  !> @}

  !> @{
  !> Matrix hermitian transpose
  !> @param[in] A -- type array of size M x N
  !> @returns B of size N x M which is the matrix transpose of A
  !> @code
  !! type(array) Re, Im, Cplx, Cplx2, Im2
  !! Re = randu(5, 5)         ! Random matrix, real part
  !! Im = randu(5, 5)         ! Random matrix, imaginary part
  !! Cplx = complex(Re, Im)   ! Create complex matrix
  !! Cplx2 = htranspose(Cplx) ! Hermitian transpose
  !! Im2 = imag(Cplx2)        ! Extract imaginary part. same
  !! call print(sum(Im + Im2)) ! Im is equal to -Im2
  !! @endcode
  interface htranspose
     module procedure array_htranspose
  end interface htranspose
  !> @}

  !> @{
  !> Matrix reorder
  !> @param[in] A -- type array of size M x N
  !> @param[in] d1 -- integer, optional. 1st output dimension corresponding to an input dimension.
  !> @param[in] d2 -- integer, optional. 2nd output dimension corresponding to an input dimension.
  !> @param[in] d3 -- integer, optional. 3rd output dimension corresponding to an input dimension.
  !> @param[in] d4 -- integer, optional. 4th output dimension corresponding to an input dimension.
  !> @returns B of size N x M which is the matrix reorder of A
  !> @code
  !! type(array) A, B
  !! A = randu(5, 2, 3, 4)       ! Random matrix of size 5 x 2 x 3 x 4
  !! B = reorder(A, 4, 1, 3, 2)  ! reordered matrix of size 4 x 5 x 3 x 2
  !! @endcode
  interface reorder
     module procedure array_reorder
  end interface reorder
  !> @}

  !> @{
  !> Matrix sort
  !> @param[in] A - Matrix
  !> @returns B whose columns are in sorted order of A
  !> @code
  !! type(array) A, B
  !! A = randu(10, 10)  ! Random matrix
  !! B = sort(A)        ! sort each column
  !! @endcode
  interface sort
     module procedure array_sort
  end interface sort
  !> @}
  !> @}

  !> @defgroup extract Extraction: (lower, upper, real, imaginary)
  !> @{
  !> Functions that perform data re-ordering and sorting

  !> @{
  !> Matrix lower
  !> @param[in] A - Matrix
  !> @returns B the lower marix of A
  !> @code
  !! type(array) A, B
  !! A = randu(10, 10)  ! Random matrix
  !! B = lower(A)       ! Extract lower matrix of A
  !! @endcode
  interface lower
     module procedure array_lower
  end interface lower
  !> @}

  !> @{
  !> Matrix upper
  !> @param[in] A - Matrix
  !> @returns B the upper marix of A
  !> @code
  !! type(array) A, B
  !! A = randu(10, 10)  ! Random matrix
  !! B = upper(A)       ! Extract upper matrix of A
  !! @endcode
  interface upper
     module procedure array_upper
  end interface upper
  !> @}

  !> @{
  !> Matrix diag
  !> @param[in] A - Matrix
  !> @returns B the diag marix of A
  !> @code
  !! type(array) A, B
  !! A = randu(10, 10)  ! Random matrix
  !! B = diag(A)       ! Extract diag matrix of A
  !! @endcode
  interface diag
     module procedure array_diag
  end interface diag
  !> @}

  !> @{
  !> Matrix real
  !> @param[in] A - Matrix
  !> @returns B the real marix of A
  !> @code
  !! type(array) Re, Im, Cplx, Re2, Im2
  !! Re = randu(10, 10)     ! Random matrix, real part
  !! Im = randu(10, 10)     ! Random matrix, imaginary part
  !! Cplx = complex(Re, Im) ! Create complex matrix
  !! Re2 = real(Cplx)       ! Extract real part, Same as Re
  !! Im2 = imag(Cplx)       ! Extract imaginary part, Same as Im
  !! @endcode
  interface real
     module procedure array_real
  end interface real
  !> @}

  !> @{
  !> Matrix imag
  !> @param[in] A - Matrix
  !> @returns B the imag marix of A
  !> @code
  !! type(array) Re, Im, Cplx, Re2, Im2
  !! Re = randu(10, 10)     ! Random matrix, real part
  !! Im = randu(10, 10)     ! Random matrix, imaginary part
  !! Cplx = complex(Re, Im) ! Create complex matrix
  !! Re2 = real(Cplx)       ! Extract real part, Same as Re
  !! Im2 = imag(Cplx)       ! Extract imaginary part, Same as Im
  !! @endcode
  interface imag
     module procedure array_imag
  end interface imag
  !> @}

  !> @{
  !> Matrix complex
  !> @param[in] A - Matrix
  !> @returns B the complex marix of A
  !> @code
  !! type(array) Re, Im, Cplx, Re2, Im2
  !! Re = randu(10, 10)     ! Random matrix, real part
  !! Im = randu(10, 10)     ! Random matrix, imaginary part
  !! Cplx = complex(Re, Im) ! Create complex matrix
  !! Re2 = real(Cplx)       ! Extract real part, Same as Re
  !! Im2 = imag(Cplx)       ! Extract imaginary part, Same as Im
  !! @endcode
  interface complex
     module procedure array_complex, array_complex2
  end interface complex
  !> @}

  !> @{
  !> Matrix complex conjugate
  !> @param[in] A - Matrix (complex type)
  !> @returns B the complex conjugate marix of A
  !> @code
  !! type(array) Re, Im, Cplx, Conj
  !! Re = randu(10, 10)     ! Random matrix, real part
  !! Im = randu(10, 10)     ! Random matrix, imaginary part
  !! Cplx = complex(Re, Im) ! Create complex matrix
  !! Conj = conjg(Cplx)     ! Create complex conjugate
  !! @endcode
  interface conjg
     module procedure array_conjg
  end interface conjg
  !> @}
  !> @}

  !> @defgroup help Helper functions (show, info, eval, sync)
  !> @{
  !> Functions useful for debugging.

  !> @{
  !> Displays the contens of type array
  !> @param[in] in -- Conents of type array
  !> @param[in] STR -- String to be displayed before displaying in (Optional)
  !> @code
  !! call print(randu(3,3))
  !! call print(constant(1,3,3), "Matrix of all ones")
  !! @endcode
  interface print
     module procedure array_print
  end interface print
  !> @}

  !> @{
  !> Displays system and arayfire information
  !> @code
  !! call device_info()
  !! @endcode
  interface device_info
     module procedure device_info_
  end interface device_info
  !> @}

  !> @{
  !> Syncs all operations on the current GPU
  !> @code
  !! call device_sync()
  !! @endcode
  interface device_sync
     module procedure device_sync_
  end interface device_sync
  !> @}

  !> @{
  !> Forces evaluation on the device
  !> @code
  !! call device_eval(A + B)
  !! @endcode
  interface device_eval
     module procedure device_eval_
  end interface device_eval
  !> @}
  !> @}

  !> @defgroup gpu Multi-GPU and device handling
  !> @{
  !> Functions useful for Multi-GPU. Requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.

  !> @{
  !> Get the count of gpus available
  !> @code
  !! integer count = device_count
  !! @endcode
  interface device_count
     module procedure device_count_
  end interface device_count
  !> @}

  !> @{
  !> Switch to a particular gpu
  !> @code
  !! call device_set(1) ! Switch to gpu 1
  !! @endcode
  interface device_set
     module procedure device_set_
  end interface device_set
  !> @}

  !> @{
  !> Get current gpu
  !> @code
  !! integer current = device_get() ! 0 by default
  !! @endcode
  interface device_get
     module procedure device_get_
  end interface device_get
  !> @}

  !> @}

  !> @defgroup time Timing code
  !> @{
  !> Functions useful for timing GPU code

  !> @{
  !> Start timer
  !> @code
  !! double precision elapsed
  !! type(array) A, B
  !! A = randu(1024, 2048)
  !! call timer_start()
  !! B = matmul(A, transpose(A))
  !! elapsed = timer_stop()
  !! @endcode
  interface timer_start
     module procedure timer_start_
  end interface timer_start
  !> @}

  !> @{
  !> Stop timer
  !> @code
  !! double precision elapsed
  !! type(array) A, B
  !! A = randu(1024, 2048)
  !! call timer_start()
  !! B = matmul(A, transpose(A))
  !! elapsed = timer_stop()
  !! @endcode
  interface timer_stop
     module procedure timer_stop_
  end interface timer_stop
  !> @}

  !> @}
  !> @}

  !> @defgroup op Mathematical operations
  !> @{
  !> Arithmetic, relational, logical and other element wise operations


  !> @defgroup arith Arithmetic operations
  !> Basic arithmetic operations
  !> @code
  !! type(array) lhs, rhs, res
  !! real scalar
  !! scalar = 0.5
  !! lhs = constant(1,3, 3)   ! All zeros
  !! rhs = constant(0,3, 3)  ! All ones
  !! res = lhs + rhs    ! All zeros
  !! res = lhs * rhs    ! All ones
  !! res = lhs - scalar ! All 0.5
  !! res = lhs / scalar ! All 2.0
  !! @endcode
  !! @{

  !> @{
  !> Element wise addition
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array, real or integer
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (+)
     module procedure array_plus
     module procedure array_plus_s, array_lplus_s
     module procedure array_plus_d, array_lplus_d
     module procedure array_plus_i, array_lplus_i
  end interface operator (+)
  !> @}

  !> @{
  !> Element wise subtraction
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array, real or integer
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (-)
     module procedure array_minus
     module procedure array_minus_s, array_lminus_s
     module procedure array_minus_d, array_lminus_d
     module procedure array_minus_i, array_lminus_i
     module procedure array_negate
  end interface operator (-)
  !> @}

  !> @{
  !> Element wise multiplication
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array, real or integer
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (*)
     module procedure array_times
     module procedure array_times_s, array_ltimes_s
     module procedure array_times_d, array_ltimes_d
     module procedure array_times_i, array_ltimes_i
  end interface operator (*)
  !> @}

  !> @{
  !> Element wise division
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array, real or integer
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (/)
     module procedure array_div
     module procedure array_div_s, array_ldiv_s
     module procedure array_div_d, array_ldiv_d
     module procedure array_div_i, array_ldiv_i
  end interface operator (/)
  !> @}

  !> @{
  !> Element wise power
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Should be real or integer
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (**)
     module procedure array_pow
     module procedure array_pow_s, array_pow_d, array_pow_i
  end interface operator (**)
  !> @}
  !> @}

  !> @defgroup relate Relational operations (<, >, ==, /=)
  !> Relational operations
  !> @code
  !! type(array) lhs, rhs, res
  !! real scalar
  !! scalar = 0.5
  !! lhs = constant(1,3, 3)    ! All zeros
  !! rhs = constant(0,3, 3)   ! All ones
  !! res = lhs < rhs     ! All false
  !! res = lhs > rhs     ! All true
  !! res = lhs == scalar ! All false
  !! res = lhs /= scalar ! All true
  !! @endcode
  !! @{

  !> @{
  !> Compare greater than
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array or real
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (>)
     module procedure array_ge
     module procedure array_gt_s, array_lgt_s
     module procedure array_gt_d, array_lgt_d
     module procedure array_gt_i, array_lgt_i
  end interface operator (>)
  !> @}

  !> @{
  !> Compare less than
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array or real
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (<)
     module procedure array_le
     module procedure array_lt_s, array_llt_s
     module procedure array_lt_d, array_llt_d
     module procedure array_lt_i, array_llt_i
  end interface operator (<)
  !> @}

  !> @{
  !> Compare greater than or equal to
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array or real
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (>=)
     module procedure array_ge
     module procedure array_ge_s, array_lge_s
     module procedure array_ge_d, array_lge_d
     module procedure array_ge_i, array_lge_i
  end interface operator (>=)
  !> @}

  !> @{
  !> Compare less than or equal to
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array or real
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (<=)
     module procedure array_le
     module procedure array_le_s, array_lle_s
     module procedure array_le_d, array_lle_d
     module procedure array_le_i, array_lle_i
  end interface operator (<=)
  !> @}

  !> @{
  !> Compare equal to
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array or real
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (==)
     module procedure array_eq
     module procedure array_eq_s, array_leq_s
     module procedure array_eq_d, array_leq_d
     module procedure array_eq_i, array_leq_i
  end interface operator (==)
  !> @}

  !> @{
  !> Compare not equal to
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- Can be type array or real
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (/=)
     module procedure array_ne
     module procedure array_ne_s, array_lne_s
     module procedure array_ne_d, array_lne_d
     module procedure array_ne_i, array_lne_i
  end interface operator (/=)
  !> @}
  !> @}

  !> @defgroup logic Logical operations (.and., .or., .not.)
  !> Relational operations
  !> @code
  !! type(array) lhs, rhs, res
  !! real scalar
  !! scalar = 0.5
  !! lhs = randu(3, 3)
  !! rhs = randu(3, 3)
  !! res = lhs < rhs .and. lhs == 0.5
  !! res = rhs > 0.3 .or. rhs <= 0.25
  !! @endcode
  !! @{

  !> @{
  !> Logical AND operator
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- should be type array
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (.and.)
     module procedure array_and
  end interface operator (.and.)
  !> @}

  !> @{
  !> Logical OR operator
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- should be type array
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (.or.)
     module procedure array_or
  end interface operator (.or.)
  !> @}

  !> @{
  !> Logical NOT operator
  !> @param[in] lhs -- Should be type array
  !> @param[in] rhs -- should be type array
  !> @returns output which contains element wise operation of lhs and rhs
  interface operator (.not.)
     module procedure array_not
  end interface operator (.not.)
  !> @}
  !> @}

  !> @defgroup elem Element wise functions (sin, cos, exp, log ..)
  !> @code
  !! type(array) in, out
  !! in = randu(3, 3)     ! Random matrix
  !! out = sin(in)        ! similarily cos(in), tan(in)
  !! out = exp(in)        ! similarily log(in)
  !! out = abs(in - 0.5)  ! absolute values
  !! @endcode
  !! @{

  !> @{
  !> sine of an array
  !> @param[in] in -- Should be type array
  !> @returns output which performs the function element wise
  interface sin
     module procedure array_sin
  end interface sin
  !> @}

  !> @{
  !> co-sine of an array
  !> @param[in] in -- Should be type array
  !> @returns output which performs the function element wise
  interface cos
     module procedure array_cos
  end interface cos
  !> @}

  !> @{
  !> tangent of an array
  interface tan
     module procedure array_tan
  end interface tan
  !> @}

  !> @{
  !> logarithm of an array
  !> @param[in] in -- Should be type array
  !> @returns output which performs the function element wise
  interface log
     module procedure array_log
  end interface log
  !> @}

  !> @{
  !> exponential of an array
  !> @param[in] in -- Should be type array
  !> @returns output which performs the function element wise
  interface exp
     module procedure array_exp
  end interface exp
  !> @}

  !> @{
  !> absolute of an array
  !> @param[in] in -- Should be type array
  !> @returns output which performs the function element wise
  interface abs
     module procedure array_abs
  end interface abs
  !> @}

  !> @}
  !> @}


  !> @defgroup alg Linear Algebra
  !> @{
  !> Linear Algebra routines

  !> @defgroup blas Basic linear algebra (Matrix multiply and dot product)
  !! @{
  !> Matrix multiply, inner and outer products

  !> @{
  !> Matrix multiply
  !> @param[in] A -- type array of size M x K
  !> @param[in] B -- type array of size K x N
  !> @returns C of size M x N which is the matrix product of A, B
  !> @code
  !! type(array) A, B, C
  !! A = randu(20, 25)  ! Random matrix
  !! B = randu(25, 10)  ! Random matrix
  !! C = matmul(A, B)   ! Matrix multiply
  !! @endcode
  interface matmul
     module procedure array_matmul
  end interface matmul
  !> @}

  !> @}


  !> @defgroup dla Factorization: lu, qr, cholesky, singular values
  !> @{
  !> Dense linear algebra: Factorization routines

  !> @{
  !> LU decomposition
  !> Double-precision or complex input requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.
  !> (Double and
  !> @param[out] L -- Optional (Contains lower triangular matrix on exit)
  !> @param[out] U -- Optional (Contains upper triangular matrix on exit)
  !> @param[out] p -- Optional (Contains permutation  matrix on exit, such that L * U = A * p)
  !> @param[in] A - Input matrix. (Contains packed LU matrices, if performed in place)
  !> @code
  !! type(array) A, B
  !! type(array) L, U, p
  !! A = randu(5,5)      ! Generate random matrix
  !! B = A               ! Make a copy of A
  !! call lu(L, U, p, A) ! Out of place LU Decomposition
  !! call lu(B)          ! In place decomposition
  !! !B now has L below the diagonal and U above the diagonal
  !! @endcode
  interface lu
     module procedure array_lu, array_lu_inplace
  end interface lu
  !> @}

  !> @{
  !> Cholesky decomposition
  !> Double-precision or complex input requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.
  !> (Double and
  !> @param[out] R -- Optional (Contains lower triangular matrix on exit, such that A = R*transpose(R))
  !> @param[in] A - Input matrix. (Contains R in the lower triangle, if performed in place)
  !> @code
  !! type(array) A, B
  !! type(array) R
  !! A = randu(5,5)      ! Generate random matrix
  !! B = A               ! Make a copy of A
  !! call cholesky(R, A) ! Out of place cholesky Decomposition
  !! call cholesky(B)    ! In place decomposition
  !! !B now has R below the diagonal
  !! @endcode
  interface cholesky
     module procedure array_cholesky, array_cholesky_inplace
  end interface cholesky
  !> @}

  !> @{
  !> QR decomposition
  !> Double-precision or complex input requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.
  !> (Double and
  !> @param[out] Q -- Contains the orthogonal basis of input matrix
  !> @param[out] R -- Contains a upper upper trianular matrix such that A = Q * R
  !> @param[in]  A -- Input matrix.
  !> @code
  !! type(array) A
  !! type(array) Q, R
  !! A = randu(5,5)   ! Generate random matrix
  !! call qr(Q, R, A)
  !! @endcode
  interface qr
     module procedure array_qr
  end interface qr
  !> @}

  !> @{
  !> singular value decomposition
  !> This function requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.
  !> @param[out] S -- Contains the singular values of input
  !> @param[out] U -- Contains the left unitary matrix
  !> @param[out] V -- Contains the right unitary matrix
  !> @param[in]  A -- Input matrix.
  !> @code
  !! type(array) A
  !! type(array) S, V, U
  !! A = randu(5,5)     ! Generate random matrix
  !! call singular(S, U, V, A)
  !! @endcode
  interface singular
     module procedure array_singular
  end interface singular
  !> @}
  !> @}

  !> @defgroup solve Solving linear systems
  !> @{
  !> Comprehensive options for solving linear systems

  !> @{
  !> Solve a system of equations
  !> Double-precision or complex input requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.
  !> @param[out] R -- Contains the solution to the system of equations
  !> @param[in]  A -- Co-efficient matrix
  !> @param[in]  B -- Observations
  !> @code
  !! type(array) A, B, X0, X
  !! ! Generate random system
  !! A  = randu(5,5)
  !! X0 = randu(5,1)
  !! B = matmul(A, X0)
  !! X = solve(A, B)
  !! @endcode
  interface solve
     module procedure array_solve
  end interface solve
  !> @}
  !> @}

  !> @defgroup linops Other Linear algebra operations: inverse, matrix power, norm, rank
  !> @{
  !> Matrix inverse, power, norm and rank

  !> @{
  !> inverse of a matrix
  !> This function requires <a href=http://accelereyes.com/products/arrayfire_licensing>ArrayFire Pro</a>.
  !> @param[out] R -- Contains the inverse values of input
  !> @param[in]  A -- Input matrix.
  !> @code
  !! type(array) A
  !! type(array) R
  !! A = randu(5,5)      ! Generate random matrix
  !! call inverse( R, A)
  !! @endcode
  interface inverse
     module procedure array_inverse
  end interface inverse
  !> @}

  !> @{
  !> Matrix norm
  !> @param[in] A -- type array of size M x N
  !> @returns B of size N x M which is the matrix norm of A
  !> @code
  !! type(array) A
  !! double precision :: res
  !! real :: p = 1.5
  !! A = randu(10,10) ! Random matrix
  !! res = norm(A)    ! 2-norm
  !! res = norm(A, p) ! p-norm
  !! @endcode
  interface norm
     module procedure array_norm, array_pnorm
  end interface norm
  !> @}
  !> @}
  !> @}

  !> @defgroup data Data Analysis
  !> Routines useful for data analysis
  !> @{

  !> @defgroup sumprod Sum and product
  !> Sum and product of elements in the array
  !> @code
  !! type(array) A, sm, pr
  !! ! Generate random system
  !! A  = randu(5,5)
  !! sm = sum(A)        ! Get sum of elements along columns (default). Same as sum(A,1)
  !! pr = mul(A, 2)     ! Get multiplication of elements along rows
  !! @endcode
  !! @{

  !> @{
  !> Sum of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Sum of the input
  interface sum
     module procedure array_sum
  end interface sum
  !> @}

  !> @{
  !> Product of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Product of input
  interface product
     module procedure array_product
  end interface product
  !> @}
  !> @}

  !> @defgroup minmax Minimum and Maximum
  !> Minimum and maximum values of elements in the array
  !> @code
  !! type(array) A, mn, mx
  !! ! Generate random system
  !! A  = randu(5,5)
  !! mx = max(A)        ! Get maximum of elements along columns (default). Same as max(A,1)
  !! mn = min(A, 2)     ! Get minimum of elements along rows
  !! @endcode
  !! @{

  !> @{
  !> Minimum of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Minimum value of input
  interface min
     module procedure array_min
  end interface min
  !> @}

  !> @{
  !> Maximum of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Maximum value of input
  interface max
     module procedure array_max
  end interface max
  !> @}
  !> @}

  !> @defgroup anyall Test if any / all true
  !> Test if any / all elements in the array are true
  !> @code
  !! type(array) A, an, al
  !! ! Generate random system
  !! A  = randu(5,5)
  !! an = anytrue(A >= 1.0) ! Check if any elment is greater than 1.0
  !! al = alltrue(A >= 0.0) ! Check if all elements are greater than 0.0
  !! @endcode
  !! @{

  !> @{
  !> Find if any element in an array is true, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- true if anything in input is true, else false
  interface anytrue
     module procedure array_anytrue
  end interface anytrue
  !> @}

  !> @{
  !> Find if all elements in an array are true, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- true if everything in input is true, else false
  interface alltrue
     module procedure array_alltrue
  end interface alltrue
  !> @}
  !> @}

  !> @defgroup stat Statistical functions: mean, median, standard deviation, variance
  !> Helpful statistical functions
  !> @code
  !! type(array) A, mn, vr
  !! ! Generate random system
  !! A  = randu(5,5)
  !! an = mean(A) ! Get the mean of A
  !! al = var(A)  ! Get the variance of A
  !! @endcode
  !! @{

  !> @{
  !> Mean of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Mean value of input
  interface mean
     module procedure array_mean
  end interface mean
  !> @}

  !> @{
  !> Standard deviation of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Standard deviation of input
  interface std
     module procedure array_std
  end interface std
  !> @}

  !> @{
  !> Variance of elements in an array, along a given dimension
  !> @param[in]  A -- Input matrix
  !> @param[in]  dim -- Integer (dimension of the operation). Optional. Default: 1
  !> @returns R -- Variance of input
  interface var
     module procedure array_var
  end interface var
  !> @}

  !> @}
  !> @}

contains

  function elements(A) result(num)
    type(array), intent(in) :: A
    integer :: num
    num = product(A%shape)
  end function elements

  function safeidx(d) result(idx)
    integer, dimension(:), intent(in) :: d
    integer, dimension(3) :: idx
    integer, allocatable, dimension(:) :: S
    integer :: f
    integer :: l
    integer :: st

    S = shape(d)

    if (S(1) == 1) then
       f = d(1)
       l = d(1)
       st = 1
    end if

    if (S(1) == 2) then
       f = d(1)
       l = d(2)
       st = 1
    end if

    if (S(1) == 3) then
       f = d(1)
       l = d(2)
       st = d(3)
    end if

    idx = [f-1, l-1, st]

  end function safeidx

  subroutine init_1d(A, S)
    type(array), intent(inout) :: A
    integer, intent(in) :: S(1)
    A%shape(1) = S(1)
    A%shape(2) = 1
    A%shape(3) = 1
    A%shape(4) = 1
    A%rank = 1
  end subroutine init_1d

  subroutine init_2d(A, S)
    type(array), intent(inout) :: A
    integer, intent(in) :: S(2)
    A%shape(1) = S(1)
    A%shape(2) = S(2)
    A%shape(3) = 1
    A%shape(4) = 1
    A%rank = 2
  end subroutine init_2d

  subroutine init_3d(A, S)
    type(array), intent(inout) :: A
    integer, intent(in) :: S(3)
    A%shape(1) = S(1)
    A%shape(2) = S(2)
    A%shape(3) = S(3)
    A%shape(4) = 1
    A%rank = 3
  end subroutine init_3d

  subroutine init_4d(A, S)
    type(array), intent(inout) :: A
    integer, intent(in) :: S(4)
    A%shape(1) = S(1)
    A%shape(2) = S(2)
    A%shape(4) = S(3)
    A%shape(4) = S(4)
    A%rank = 4
  end subroutine init_4d

  subroutine init_eq(L, R)
    type(array), intent(inout) :: L
    type(array), intent(in) :: R
    L%rank  = R%rank
    L%shape = R%shape
  end subroutine init_eq

  function idx_scalar(scalar) result(R)
    integer, intent(in) :: scalar
    type(array) :: R
    integer :: S(1)
    S = [1]
    call init_1d(R, S)
    call af_idx_seq(R%ptr, scalar, scalar, 1, err)
  end function idx_scalar

  function idx_vector(indices) result(R)
    integer, intent(in) :: indices(:)
    type(array) :: R
    call init_1d(R, shape(indices))
    call af_idx_vec(R%ptr, indices, elements(R), err)
  end function idx_vector

  function arr_seq(first, last, step) result(R)
    integer, intent(in), optional :: first
    integer, intent(in), optional :: last
    integer, intent(in), optional :: step
    integer :: f, l, s
    type(array) :: R
    f = 0
    l = 0
    s = 1

    if (present(first)) f = first
    if (present(last )) l = last
    if (present(step )) s = step

    call af_idx_seq(R%ptr, f, l, s, err)
    call init_post(R%ptr, R%shape, R%rank)

  end function arr_seq

  function array_get(in, d1, d2, d3, d4) result(R)
    type(array), intent(in) :: in
    type(array), intent(in) :: d1
    type(array), intent(in), optional :: d2
    integer, dimension(:), intent(in), optional :: d3
    integer, dimension(:), intent(in), optional :: d4
    integer :: dims

    type(array) :: R
    type(C_ptr) :: idx1 = C_NULL_ptr
    type(C_ptr) :: idx2 = C_NULL_ptr
    integer, dimension(3) :: idx3
    integer, dimension(3) :: idx4

    idx1 = d1%ptr
    dims = 1

    if (present(d2)) then
       idx2 = d2%ptr
       dims = 2
    end if

    if (present(d3)) then
       idx3 = safeidx(d3)
       dims = 3
    end if

    if (present(d4)) then
       idx4 = safeidx(d4)
       dims = 4
    end if

    call af_arr_get(R%ptr, in%ptr, idx1, idx2, idx3, idx4, dims, err)
    call init_post(R%ptr, R%shape, R%rank)

  end function array_get

  function array_get2(in, d1, d2, d3) result(R)
    type(array), intent(in) :: in
    type(array), intent(in) :: d1
    integer, dimension(:), intent(in) :: d2
    integer, dimension(:), intent(in), optional :: d3
    integer :: dims

    type(array) :: R
    type(C_ptr) :: idx1 = C_NULL_ptr
    integer, dimension(3) :: idx2
    integer, dimension(3) :: idx3

    idx1 = d1%ptr
    idx2 = safeidx(d2)
    dims = 2

    if (present(d3)) then
       idx3 = safeidx(d3)
       dims = 3
    end if

    call af_arr_get2(R%ptr, in%ptr, idx1, idx2, idx3, dims, err)
    call init_post(R%ptr, R%shape, R%rank)

  end function array_get2

  function array_get_seq(in, d1, d2, d3, d4) result(R)
    type(array), intent(in) :: in
    integer, intent(in) :: d1(:)
    integer, intent(in), optional :: d2(:)
    integer, intent(in), optional :: d3(:)
    integer, intent(in), optional :: d4(:)
    type(array) :: R

    integer, dimension(3) :: idx1
    integer, dimension(3) :: idx2
    integer, dimension(3) :: idx3
    integer, dimension(3) :: idx4
    integer :: dims = 1

    idx1 = safeidx(d1)
    idx2 = safeidx(d1)
    idx3 = safeidx(d1)
    idx4 = safeidx(d1)

    if (present(d2)) then
       idx2 = safeidx(d2)
       dims = 2
    end if

    if (present(d3)) then
       idx3 = safeidx(d3)
       dims = 3
    end if

    if (present(d4)) then
       idx4 = safeidx(d4)
       dims = 4
    end if

    call af_arr_get_seq(R%ptr, in%ptr, idx1, idx2, idx3, idx4, dims, err)
    call init_post(R%ptr, R%shape, R%rank)
  end function array_get_seq

  subroutine array_set(lhs, rhs, d1, d2, d3, d4)
    type(array), intent(in) :: lhs
    type(array), intent(inout) :: rhs
    type(array), intent(in) :: d1
    type(array), intent(in), optional :: d2
    integer, dimension(:), intent(in), optional :: d3
    integer, dimension(:), intent(in), optional :: d4

    type(C_ptr) :: idx1 = C_NULL_ptr
    type(C_ptr) :: idx2 = C_NULL_ptr
    integer, dimension(3) :: idx3
    integer, dimension(3) :: idx4
    integer :: dims

    idx1 = d1%ptr
    dims = 1

    if (present(d2)) then
       idx2 = d2%ptr
       dims = 2
    end if

    if (present(d3)) then
       idx3 = safeidx(d3)
       dims = 3
    end if

    if (present(d4)) then
       idx4 = safeidx(d4)
       dims = 4
    end if

    call af_arr_set(lhs%ptr, rhs%ptr, idx1, idx2, idx3, idx4, dims, err)
  end subroutine array_set

  subroutine array_set2(lhs, rhs, d1, d2, d3)
    type(array), intent(in) :: lhs
    type(array), intent(inout) :: rhs
    type(array), intent(in) :: d1
    integer, dimension(:), intent(in) :: d2
    integer, dimension(:), intent(in), optional :: d3

    type(C_ptr) :: idx1 = C_NULL_ptr
    integer, dimension(3) :: idx2
    integer, dimension(3) :: idx3
    integer :: dims

    idx1 = d1%ptr
    idx2 = safeidx(d2)
    dims = 2

    if (present(d3)) then
       idx3 = safeidx(d3)
       dims = 3
    end if

    call af_arr_set2(lhs%ptr, rhs%ptr, idx1, idx2, idx3, dims, err)
  end subroutine array_set2

  subroutine array_set_seq(R, in, d1, d2, d3, d4)
    type(array), intent(in) :: in
    integer, intent(in) :: d1(:)
    integer, intent(in), optional :: d2(:)
    integer, intent(in), optional :: d3(:)
    integer, intent(in), optional :: d4(:)
    type(array), intent(inout) :: R

    integer, dimension(3) :: idx1
    integer, dimension(3) :: idx2
    integer, dimension(3) :: idx3
    integer, dimension(3) :: idx4
    integer :: dims = 1

    idx1 = safeidx(d1)
    idx2 = safeidx(d1)
    idx3 = safeidx(d1)
    idx4 = safeidx(d1)

    if (present(d2)) then
       idx2 = safeidx(d2)
       dims = 2
    end if

    if (present(d3)) then
       idx3 = safeidx(d3)
       dims = 3
    end if

    if (present(d4)) then
       idx4 = safeidx(d4)
       dims = 4
    end if

    call af_arr_set_seq(R%ptr, in%ptr, idx1, idx2, idx3, idx4, dims, err)
    call init_post(R%ptr, R%shape, R%rank)
  end subroutine array_set_seq

  !> Assigns data to array
  subroutine assign(L, R)
    type(array), intent(inout) :: L
    type(array), intent(in) :: R
    call init_eq(L, R)
    call af_arr_copy(L%ptr, R%ptr, err)
  end subroutine assign

  !> Assigns data to array
  subroutine device1_s(A, B)
    type(array), intent(inout) :: A
    real, intent(in) :: B(:)
    call init_1d(A, shape(B))
    call af_arr_device_s(A%ptr, B, A%shape, err)
  end subroutine device1_s

  !> Assigns data to array
  subroutine device1_d(A, B)
    type(array), intent(inout) :: A
    double precision, intent(in) :: B(:)
    call init_1d(A, shape(B))
    call af_arr_device_d(A%ptr, B, A%shape, err)
  end subroutine device1_d

  !> Assigns data to array
  subroutine device1_c(A, B)
    type(array), intent(inout) :: A
    complex, intent(in) :: B(:)
    call init_1d(A, shape(B))
    call af_arr_device_c(A%ptr, B, A%shape, err)
  end subroutine device1_c

  !> Assigns data to array
  subroutine device1_z(A, B)
    type(array), intent(inout) :: A
    double complex, intent(in) :: B(:)
    call init_1d(A, shape(B))
    call af_arr_device_z(A%ptr, B, A%shape, err)
  end subroutine device1_z

  !> Assigns data to array
  subroutine device2_s(A, B)
    type(array), intent(inout) :: A
    real, intent(in) :: B(:,:)
    call init_2d(A, shape(B))
    call af_arr_device_s(A%ptr, B, A%shape, err)
  end subroutine device2_s

  !> Assigns data to array
  subroutine device2_d(A, B)
    type(array), intent(inout) :: A
    double precision, intent(in) :: B(:,:)
    call init_2d(A, shape(B))
    call af_arr_device_d(A%ptr, B, A%shape, err)
  end subroutine device2_d

  !> Assigns data to array
  subroutine device2_c(A, B)
    type(array), intent(inout) :: A
    complex, intent(in) :: B(:,:)
    call init_2d(A, shape(B))
    call af_arr_device_c(A%ptr, B, A%shape, err)
  end subroutine device2_c

  !> Assigns data to array
  subroutine device2_z(A, B)
    type(array), intent(inout) :: A
    double complex, intent(in) :: B(:,:)
    call init_2d(A, shape(B))
    call af_arr_device_z(A%ptr, B, A%shape, err)
  end subroutine device2_z

  !> Assigns data to array
  subroutine device3_s(A, B)
    type(array), intent(inout) :: A
    real, intent(in) :: B(:,:,:)
    call init_3d(A, shape(B))
    call af_arr_device_s(A%ptr, B, A%shape, err)
  end subroutine device3_s

  !> Assigns data to array
  subroutine device3_d(A, B)
    type(array), intent(inout) :: A
    double precision, intent(in) :: B(:,:,:)
    call init_3d(A, shape(B))
    call af_arr_device_d(A%ptr, B, A%shape, err)
  end subroutine device3_d

  !> Assigns data to array
  subroutine device3_c(A, B)
    type(array), intent(inout) :: A
    complex, intent(in) :: B(:,:,:)
    call init_3d(A, shape(B))
    call af_arr_device_c(A%ptr, B, A%shape, err)
  end subroutine device3_c

  !> Assigns data to array
  subroutine device3_z(A, B)
    type(array), intent(inout) :: A
    double complex, intent(in) :: B(:,:,:)
    call init_3d(A, shape(B))
    call af_arr_device_z(A%ptr, B, A%shape, err)
  end subroutine device3_z

  !> Assigns data to array
  subroutine device4_s(A, B)
    type(array), intent(inout) :: A
    real, intent(in) :: B(:,:,:,:)
    call init_4d(A, shape(B))
    call af_arr_device_s(A%ptr, B, A%shape, err)
  end subroutine device4_s

  !> Assigns data to array
  subroutine device4_d(A, B)
    type(array), intent(inout) :: A
    double precision, intent(in) :: B(:,:,:,:)
    call init_4d(A, shape(B))
    call af_arr_device_d(A%ptr, B, A%shape, err)
  end subroutine device4_d

  !> Assigns data to array
  subroutine device4_c(A, B)
    type(array), intent(inout) :: A
    complex, intent(in) :: B(:,:,:,:)
    call init_4d(A, shape(B))
    call af_arr_device_c(A%ptr, B, A%shape, err)
  end subroutine device4_c

  !> Assigns data to array
  subroutine device4_z(A, B)
    type(array), intent(inout) :: A
    double complex, intent(in) :: B(:,:,:,:)
    call init_4d(A, shape(B))
    call af_arr_device_z(A%ptr, B, A%shape, err)
  end subroutine device4_z

  !> Get the array data back to host
  subroutine host1_s(R, A)
    type(array), intent(in) :: A
    real, intent(inout), dimension(:), allocatable :: R
    allocate(R(A%shape(1)))
    call af_arr_host_s(R, A%ptr, 1, err)
  end subroutine host1_s

  !> Get the array data back to host
  subroutine host1_d(R, A)
    type(array), intent(in) :: A
    double precision, intent(inout), dimension(:), allocatable :: R
    allocate(R(A%shape(1)))
    call af_arr_host_d(R, A%ptr, 1, err)
  end subroutine host1_d

  !> Get the array data back to host
  subroutine host1_c(R, A)
    type(array), intent(in) :: A
    complex, intent(inout), dimension(:), allocatable :: R
    allocate(R(A%shape(1)))
    call af_arr_host_c(R, A%ptr, 1, err)
  end subroutine host1_c

  !> Get the array data back to host
  subroutine host1_z(R, A)
    type(array), intent(in) :: A
    double complex, intent(inout), dimension(:), allocatable :: R
    allocate(R(A%shape(1)))
    call af_arr_host_z(R, A%ptr, 1, err)
  end subroutine host1_z

  !> Get the array data back to host
  subroutine host2_s(R, A)
    type(array), intent(in) :: A
    real, intent(inout), dimension(:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2)))
    call af_arr_host_s(R, A%ptr, 2, err)
  end subroutine host2_s

  !> Get the array data back to host
  subroutine host2_d(R, A)
    type(array), intent(in) :: A
    double precision, intent(inout), dimension(:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2)))
    call af_arr_host_d(R, A%ptr, 2, err)
  end subroutine host2_d

  !> Get the array data back to host
  subroutine host2_c(R, A)
    type(array), intent(in) :: A
    complex, intent(inout), dimension(:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2)))
    call af_arr_host_c(R, A%ptr, 2, err)
  end subroutine host2_c

  !> Get the array data back to host
  subroutine host2_z(R, A)
    type(array), intent(in) :: A
    double complex, intent(inout), dimension(:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2)))
    call af_arr_host_z(R, A%ptr, 2, err)
  end subroutine host2_z

  !> Get the array data back to host
  subroutine host3_s(R, A)
    type(array), intent(in) :: A
    real, intent(inout), dimension(:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3)))
    call af_arr_host_s(R, A%ptr, 3, err)
  end subroutine host3_s

  !> Get the array data back to host
  subroutine host3_d(R, A)
    type(array), intent(in) :: A
    double precision, intent(inout), dimension(:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3)))
    call af_arr_host_d(R, A%ptr, 3, err)
  end subroutine host3_d

  !> Get the array data back to host
  subroutine host3_c(R, A)
    type(array), intent(in) :: A
    complex, intent(inout), dimension(:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3)))
    call af_arr_host_c(R, A%ptr, 3, err)
  end subroutine host3_c

  !> Get the array data back to host
  subroutine host3_z(R, A)
    type(array), intent(in) :: A
    double complex, intent(inout), dimension(:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3)))
    call af_arr_host_z(R, A%ptr, 3, err)
  end subroutine host3_z

  !> Get the array data back to host
  subroutine host4_s(R, A)
    type(array), intent(in) :: A
    real, intent(inout), dimension(:,:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3), A%shape(4)))
    call af_arr_host_s(R, A%ptr, 4, err)
  end subroutine host4_s

  !> Get the array data back to host
  subroutine host4_d(R, A)
    type(array), intent(in) :: A
    double precision, intent(inout), dimension(:,:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3), A%shape(4)))
    call af_arr_host_d(R, A%ptr, 4, err)
  end subroutine host4_d

  !> Get the array data back to host
  subroutine host4_c(R, A)
    type(array), intent(in) :: A
    complex, intent(inout), dimension(:,:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3), A%shape(4)))
    call af_arr_host_c(R, A%ptr, 4, err)
  end subroutine host4_c

  !> Get the array data back to host
  subroutine host4_z(R, A)
    type(array), intent(in) :: A
    double complex, intent(inout), dimension(:,:,:,:), allocatable :: R
    allocate(R(A%shape(1), A%shape(2), A%shape(3), A%shape(4)))
    call af_arr_host_z(R, A%ptr, 4, err)
  end subroutine host4_z

  !> Get the array data back to hostp
  subroutine hostp1_s(R, A)
    type(array), intent(in) :: A
    real, pointer, intent(inout), dimension(:) :: R
    call af_arr_host_s(R, A%ptr, 1, err)
  end subroutine hostp1_s

  !> Get the array data back to hostp
  subroutine hostp1_d(R, A)
    type(array), intent(in) :: A
    double precision, pointer, intent(inout), dimension(:) :: R
    call af_arr_host_d(R, A%ptr, 1, err)
  end subroutine hostp1_d

  !> Get the array data back to hostp
  subroutine hostp1_c(R, A)
    type(array), intent(in) :: A
    complex, pointer, intent(inout), dimension(:) :: R
    call af_arr_host_c(R, A%ptr, 1, err)
  end subroutine hostp1_c

  !> Get the array data back to hostp
  subroutine hostp1_z(R, A)
    type(array), intent(in) :: A
    double complex, pointer, intent(inout), dimension(:) :: R
    call af_arr_host_z(R, A%ptr, 1, err)
  end subroutine hostp1_z

  !> Get the array data back to hostp
  subroutine hostp2_s(R, A)
    type(array), intent(in) :: A
    real, pointer, intent(inout), dimension(:,:) :: R
    call af_arr_host_s(R, A%ptr, 2, err)
  end subroutine hostp2_s

  !> Get the array data back to hostp
  subroutine hostp2_d(R, A)
    type(array), intent(in) :: A
    double precision, pointer, intent(inout), dimension(:,:) :: R
    call af_arr_host_d(R, A%ptr, 2, err)
  end subroutine hostp2_d

  !> Get the array data back to hostp
  subroutine hostp2_c(R, A)
    type(array), intent(in) :: A
    complex, pointer, intent(inout), dimension(:,:) :: R
    call af_arr_host_c(R, A%ptr, 2, err)
  end subroutine hostp2_c

  !> Get the array data back to hostp
  subroutine hostp2_z(R, A)
    type(array), intent(in) :: A
    double complex, pointer, intent(inout), dimension(:,:) :: R
    call af_arr_host_z(R, A%ptr, 2, err)
  end subroutine hostp2_z

  !> Get the array data back to hostp
  subroutine hostp3_s(R, A)
    type(array), intent(in) :: A
    real, pointer, intent(inout), dimension(:,:,:) :: R
    call af_arr_host_s(R, A%ptr, 3, err)
  end subroutine hostp3_s

  !> Get the array data back to hostp
  subroutine hostp3_d(R, A)
    type(array), intent(in) :: A
    double precision, pointer, intent(inout), dimension(:,:,:) :: R
    call af_arr_host_d(R, A%ptr, 3, err)
  end subroutine hostp3_d

  !> Get the array data back to hostp
  subroutine hostp3_c(R, A)
    type(array), intent(in) :: A
    complex, pointer, intent(inout), dimension(:,:,:) :: R
    call af_arr_host_c(R, A%ptr, 3, err)
  end subroutine hostp3_c

  !> Get the array data back to hostp
  subroutine hostp3_z(R, A)
    type(array), intent(in) :: A
    double complex, pointer, intent(inout), dimension(:,:,:) :: R
    call af_arr_host_z(R, A%ptr, 3, err)
  end subroutine hostp3_z

  !> Get the array data back to hostp
  subroutine hostp4_s(R, A)
    type(array), intent(in) :: A
    real, pointer, intent(inout), dimension(:,:,:,:) :: R
    call af_arr_host_s(R, A%ptr, 4, err)
  end subroutine hostp4_s

  !> Get the array data back to hostp
  subroutine hostp4_d(R, A)
    type(array), intent(in) :: A
    double precision, pointer, intent(inout), dimension(:,:,:,:) :: R
    call af_arr_host_d(R, A%ptr, 4, err)
  end subroutine hostp4_d

  !> Get the array data back to hostp
  subroutine hostp4_c(R, A)
    type(array), intent(in) :: A
    complex, pointer, intent(inout), dimension(:,:,:,:) :: R
    call af_arr_host_c(R, A%ptr, 4, err)
  end subroutine hostp4_c

  !> Get the array data back to hostp
  subroutine hostp4_z(R, A)
    type(array), intent(in) :: A
    double complex, pointer, intent(inout), dimension(:,:,:,:) :: R
    call af_arr_host_z(R, A%ptr, 4, err)
  end subroutine hostp4_z

  !> Display array, optionally with a message
  subroutine array_print(A, STR)
    type(array), intent(in) :: A
    character(len=*), intent(in), optional :: STR
    if (present(STR)) write(*,*) STR
    call af_arr_print(A%ptr, err)
  end subroutine array_print

  !> Moddims an input array
  function array_moddims(A, x1, x2, x3, x4) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    integer, intent(in) :: x1
    integer, intent(in), optional :: x2, x3, x4

    R%shape = [x1, 1, 1, 1]
    R%rank = 1
    if (present(x2)) then
       R%shape(2) = x2
       R%rank = 2
    end if
    if (present(x3)) then
       R%shape(3) = x3
       R%rank = 3
    end if
    if (present(x4)) then
       R%shape(4) = x4
       R%rank = 4
    end if

    call af_arr_moddims(R%ptr, A%ptr, R%shape, err)
  end function array_moddims

  !> Flat an input array
  function array_flat(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    R%shape = [1, 1, 1, 1]
    R%shape(1) = elements(A)
    R%rank = 1
    call af_arr_moddims(R%ptr, A%ptr, R%shape, err)
  end function array_flat

  !> Tile an input array
  function array_tile(A, x1, x2, x3, x4) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    integer, intent(in), optional :: x1, x2, x3, x4

    call init_eq(R, A)

    if (present(x1)) then
       R%shape(1) = R%shape(1) * x1
    end if
    if (present(x2)) then
       R%shape(2) = R%shape(2) * x2
       if (x2 > 1 .and. R%rank < 2) R%rank = 2
    end if
    if (present(x3)) then
       R%shape(3) = R%shape(3) * x3
       if (x3 > 1 .and. R%rank < 3) R%rank = 3
    end if
    if (present(x4)) then
       R%shape(4) = R%shape(4) * x4
       if (x4 > 1 .and. R%rank < 4) R%rank = 4
    end if

    call af_arr_tile(R%ptr, A%ptr, R%shape, err)
  end function array_tile

  !> Join two arrays
  function array_join(d, first, second) result(output)
    type(array), intent(in) :: first
    type(array), intent(in) :: second
    integer, intent(in) :: d
    type(array) :: output
    call af_arr_join(d, output%ptr, first%ptr, second%ptr, err)
    call init_post(output%ptr, output%shape, output%rank)
  end function array_join


  !> Generate  uniformly distributed random matrix
  function array_randu(x1, x2, x3, x4, ty) result(R)
    type(array) :: R
    integer, intent(in) :: x1
    integer, intent(in), optional :: x2, x3, x4, ty
    integer :: tt = 1

    R%shape = [x1, 1, 1, 1]
    R%rank = 1
    if (present(x2)) then
       R%shape(2) = x2
       R%rank = 2
    end if
    if (present(x3)) then
       R%shape(3) = x3
       R%rank = 3
    end if
    if (present(x4)) then
       R%shape(4) = x4
       R%rank = 4
    end if

    if (present(ty)) tt = ty

    call af_arr_randu(R%ptr, R%shape, tt, err)
  end function array_randu

  !> Generate  normally distributed random matrix
  function array_randn(x1, x2, x3, x4, ty) result(R)
    type(array) :: R
    integer, intent(in) :: x1
    integer, intent(in), optional :: x2, x3, x4, ty
    integer :: tt = 1

    R%shape = [x1, 1, 1, 1]
    R%rank = 1
    if (present(x2)) then
       R%shape(2) = x2
       R%rank = 2
    end if
    if (present(x3)) then
       R%shape(3) = x3
       R%rank = 3
    end if
    if (present(x4)) then
       R%shape(4) = x4
       R%rank = 4
    end if

    if (present(ty)) tt = ty

    call af_arr_randn(R%ptr, R%shape, tt, err)
  end function array_randn

  !> Generate an array of constant value
  function array_constant(val, x1, x2, x3, x4, ty) result(R)
    type(array) :: R
    integer, intent(in) :: val
    integer, intent(in) :: x1
    integer, intent(in), optional :: x2, x3, x4, ty
    integer :: tt = 1

    R%shape = [x1, 1, 1, 1]
    R%rank = 1
    if (present(x2)) then
       R%shape(2) = x2
       R%rank = 2
    end if
    if (present(x3)) then
       R%shape(3) = x3
       R%rank = 3
    end if
    if (present(x4)) then
       R%shape(4) = x4
       R%rank = 4
    end if

    if (present(ty)) tt = ty

    call af_arr_constant(R%ptr, val, R%shape, tt, err)
  end function array_constant

  !> Generate  an identity matrix
  function array_identity(x1, x2, ty) result(R)
    type(array) :: R
    integer, intent(in) :: x1, x2
    integer, intent(in), optional :: ty
    integer :: tt = 1

    R%shape = [x1, x2, 1, 1]
    R%rank = 2

    if (present(ty)) tt = ty

    call af_arr_identity(R%ptr, R%shape, tt, err)
  end function array_identity

  !> Add two array matrices
  function array_plus(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elplus(R%ptr, A%ptr, B%ptr, err)
  end function array_plus

  !> Subtract two array matrices
  function array_minus(A, B) result(R)
    type(array), intent(in) :: A, B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elminus(R%ptr, A%ptr, B%ptr, err)
  end function array_minus

  !> Multiply two array matrices (element wise)
  function array_times(A, B) result(R)
    type(array), intent(in) :: A, B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_eltimes(R%ptr, A%ptr, B%ptr, err)
  end function array_times

  !> Divide two array matrices (element wise)
  function array_div(A, B) result(R)
    type(array), intent(in) :: A, B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_eldiv(R%ptr, A%ptr, B%ptr, err)
  end function array_div

  !> Element wise power with matrix exponent
  function array_pow(A, B) result(R)
    type(array), intent(in) :: A, B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elpow(R%ptr, A%ptr, B%ptr, err)
  end function array_pow

  !> Negate an array
  function array_negate(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_negate(R%ptr, A%ptr, err)
  end function array_negate

  !> Add scalar to array
  function array_plus_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scplus(R%ptr, A%ptr, B, err)
  end function array_plus_d

  !> Add array to scalar
  function array_lplus_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A + B
  end function array_lplus_d

  !> Add scalar to array
  function array_plus_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A + dble(B)
  end function array_plus_s

  !> Add array to scalar
  function array_lplus_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A + dble(B)
  end function array_lplus_s

  !> Add scalar to array
  function array_plus_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A + dble(B)
  end function array_plus_i

  !> Add array to scalar
  function array_lplus_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A + dble(B)
  end function array_lplus_i

  !> Add scalar to array
  function array_minus_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scminus(R%ptr, A%ptr, B, err)
  end function array_minus_d

  !> Add array to scalar
  function array_lminus_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A - B
  end function array_lminus_d

  !> Add scalar to array
  function array_minus_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A - dble(B)
  end function array_minus_s

  !> Add array to scalar
  function array_lminus_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A - dble(B)
  end function array_lminus_s

  !> Add scalar to array
  function array_minus_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A - dble(B)
  end function array_minus_i

  !> Add array to scalar
  function array_lminus_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A - dble(B)
  end function array_lminus_i

  !> Add scalar to array
  function array_times_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_sctimes(R%ptr, A%ptr, B, err)
  end function array_times_d

  !> Add array to scalar
  function array_ltimes_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A * B
  end function array_ltimes_d

  !> Add scalar to array
  function array_times_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A * dble(B)
  end function array_times_s

  !> Add array to scalar
  function array_ltimes_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A * dble(B)
  end function array_ltimes_s

  !> Add scalar to array
  function array_times_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A * dble(B)
  end function array_times_i

  !> Add array to scalar
  function array_ltimes_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A * dble(B)
  end function array_ltimes_i

  !> Add scalar to array
  function array_div_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scdiv(R%ptr, A%ptr, B, err)
  end function array_div_d

  !> Add array to scalar
  function array_ldiv_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A / B
  end function array_ldiv_d

  !> Add scalar to array
  function array_div_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A / dble(B)
  end function array_div_s

  !> Add array to scalar
  function array_ldiv_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A / dble(B)
  end function array_ldiv_s

  !> Add scalar to array
  function array_div_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A / dble(B)
  end function array_div_i

  !> Add array to scalar
  function array_ldiv_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A / dble(B)
  end function array_ldiv_i

  !> Element wise power with scalar exponent
  function array_pow_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scpow(R%ptr, A%ptr, B, err)
  end function array_pow_d

  !> Element wise power with scalar exponent
  function array_pow_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A ** dble(B)
  end function array_pow_s

  !> Element wise power with scalar exponent
  function array_pow_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A ** dble(B)
  end function array_pow_i

  function array_gt(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elgt(R%ptr, A%ptr, B%ptr, err)
  end function array_gt

  function array_ge(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elge(R%ptr, A%ptr, B%ptr, err)
  end function array_ge

  function array_lt(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_ellt(R%ptr, A%ptr, B%ptr, err)
  end function array_lt

  function array_le(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elle(R%ptr, A%ptr, B%ptr, err)
  end function array_le

  function array_eq(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_eleq(R%ptr, A%ptr, B%ptr, err)
  end function array_eq

  function array_ne(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elne(R%ptr, A%ptr, B%ptr, err)
  end function array_ne

  function array_gt_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scgt(R%ptr, A%ptr, B, err)
  end function array_gt_d

  function array_lgt_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A < B
  end function array_lgt_d

  function array_gt_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    call af_arr_scgt(R%ptr, A%ptr, dble(B), err)
  end function array_gt_s

  function array_lgt_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A < dble(B)
  end function array_lgt_s

  function array_gt_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    call af_arr_scgt(R%ptr, A%ptr, dble(B), err)
  end function array_gt_i

  function array_lgt_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A < dble(B)
  end function array_lgt_i

  function array_lt_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_sclt(R%ptr, A%ptr, B, err)
  end function array_lt_d

  function array_llt_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A > B
  end function array_llt_d

  function array_lt_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    call af_arr_sclt(R%ptr, A%ptr, dble(B), err)
  end function array_lt_s

  function array_llt_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A > dble(B)
  end function array_llt_s

  function array_lt_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    call af_arr_sclt(R%ptr, A%ptr, dble(B), err)
  end function array_lt_i

  function array_llt_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A > dble(B)
  end function array_llt_i

  function array_ge_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scge(R%ptr, A%ptr, B, err)
  end function array_ge_d

  function array_lge_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A <= B
  end function array_lge_d

  function array_ge_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    call af_arr_scge(R%ptr, A%ptr, dble(B), err)
  end function array_ge_s

  function array_lge_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A <= dble(B)
  end function array_lge_s

  function array_ge_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    call af_arr_scge(R%ptr, A%ptr, dble(B), err)
  end function array_ge_i

  function array_lge_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A <= dble(B)
  end function array_lge_i

  function array_le_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scle(R%ptr, A%ptr, B, err)
  end function array_le_d

  function array_lle_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A >= B
  end function array_lle_d

  function array_le_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    call af_arr_scle(R%ptr, A%ptr, dble(B), err)
  end function array_le_s

  function array_lle_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A >= dble(B)
  end function array_lle_s

  function array_le_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    call af_arr_scle(R%ptr, A%ptr, dble(B), err)
  end function array_le_i

  function array_lle_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A >= dble(B)
  end function array_lle_i

  function array_eq_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_sceq(R%ptr, A%ptr, B, err)
  end function array_eq_d

  function array_leq_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A == B
  end function array_leq_d

  function array_eq_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    call af_arr_sceq(R%ptr, A%ptr, dble(B), err)
  end function array_eq_s

  function array_leq_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A == dble(B)
  end function array_leq_s

  function array_eq_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    call af_arr_sceq(R%ptr, A%ptr, dble(B), err)
  end function array_eq_i

  function array_leq_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A == dble(B)
  end function array_leq_i

  function array_ne_d(A, B) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_scne(R%ptr, A%ptr, B, err)
  end function array_ne_d

  function array_lne_d(B, A) result(R)
    type(array), intent(in) :: A
    double precision, intent(in) :: B
    type(array) :: R
    R = A /= B
  end function array_lne_d

  function array_ne_s(A, B) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    call af_arr_scne(R%ptr, A%ptr, dble(B), err)
  end function array_ne_s

  function array_lne_s(B, A) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: B
    type(array) :: R
    R = A /= dble(B)
  end function array_lne_s

  function array_ne_i(A, B) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    call af_arr_scne(R%ptr, A%ptr, dble(B), err)
  end function array_ne_i

  function array_lne_i(B, A) result(R)
    type(array), intent(in) :: A
    integer, intent(in) :: B
    type(array) :: R
    R = A /= dble(B)
  end function array_lne_i

  !> and on two array matrices
  function array_and(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_eland(R%ptr, A%ptr, B%ptr, err)
  end function array_and

  !> or on two array matrices
  function array_or(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    call init_eq(R, A)
    call af_arr_elor(R%ptr, A%ptr, B%ptr, err)
  end function array_or

  !> Not an array
  function array_not(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_not(R%ptr, A%ptr, err)
  end function array_not

  !> sin of array
  function array_sin(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_sin(R%ptr, A%ptr, err)
  end function array_sin

  !> cos of array
  function array_cos(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_cos(R%ptr, A%ptr, err)
  end function array_cos

  !> tan of array
  function array_tan(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_tan(R%ptr, A%ptr, err)
  end function array_tan

  !> log of array
  function array_log(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_log(R%ptr, A%ptr, err)
  end function array_log

  !> absolute of array
  function array_abs(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_abs(R%ptr, A%ptr, err)
  end function array_abs

  !> exponential of array
  function array_exp(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_exp(R%ptr, A%ptr, err)
  end function array_exp

  !> Multiply two array matrices
  function array_matmul(A, B) result(R)
    type(array), intent(in) :: A
    type(array), intent(in) :: B
    type(array) :: R
    R%shape(1) = A%shape(1)
    R%shape(2) = B%shape(2)
    R%shape(3) = B%shape(1)
    R%shape(4) = B%shape(1)
    R%rank = 2
    call af_arr_matmul(R%ptr, A%ptr, B%ptr, err)
  end function array_matmul

  !> Transpose an array
  function array_transpose(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    R%shape(1) = A%shape(2)
    R%shape(2) = A%shape(1)
    R%shape(3) = A%shape(3)
    R%shape(4) = A%shape(4)
    R%rank = 2
    call af_arr_t(R%ptr, A%ptr, err)
  end function array_transpose

  !> Htranspose an array
  function array_htranspose(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    R%shape(1) = A%shape(2)
    R%shape(2) = A%shape(1)
    R%shape(3) = A%shape(3)
    R%shape(4) = A%shape(4)
    R%rank = 2
    call af_arr_h(R%ptr, A%ptr, err)
  end function array_htranspose

  !> Reorder an array
  function array_reorder(A, d1, d2, d3, d4) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d1, d2, d3, d4
    type(array) :: R

    call init_eq(R, A)
    if (present(d1)) R%shape(1) = d1
    if (present(d2)) R%shape(2) = d2
    if (present(d3)) R%shape(3) = d3
    if (present(d4)) R%shape(4) = d4

    call af_arr_t(R%ptr, A%ptr, R%shape, err)
  end function array_reorder

  !> Sort an array
  function array_sort(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_sort(R%ptr, A%ptr, err)
  end function array_sort

  !> Lower triangular matrix of a matrix
  function array_lower(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_lower(R%ptr, A%ptr, err)
  end function array_lower

  !> Upper triangular matrix of a matrix
  function array_upper(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_upper(R%ptr, A%ptr, err)
  end function array_upper

  !> Diag triangular matrix of a matrix
  function array_diag(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_diag(R%ptr, A%ptr, err)
  end function array_diag

  !> Real triangular matrix of a matrix
  function array_real(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_real(R%ptr, A%ptr, err)
  end function array_real

  !> Imag triangular matrix of a matrix
  function array_imag(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_imag(R%ptr, A%ptr, err)
  end function array_imag

  !> Create complex matrix from input matrix
  function array_complex(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_complex(R%ptr, A%ptr, err)
  end function array_complex

  !>  Create complex matrix from two (real, imaginary) matrices
  function array_complex2(Re, Im) result(Cplx)
    type(array), intent(in) :: Re, Im
    type(array) :: Cplx
    call init_eq(Cplx, Re)
    call af_arr_complex2(Cplx%ptr, Re%ptr, Im%ptr, err)
  end function array_complex2

  !> Create conjugate matrix from input matrix
  function array_conjg(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_conjg(R%ptr, A%ptr, err)
  end function array_conjg

  !> Norm of an array
  function array_norm(A) result(R)
    type(array), intent(in) :: A
    double precision :: R
    call af_arr_norm(R, A%ptr, err)
  end function array_norm

  !> Norm of an array
  function array_pnorm(A, p) result(R)
    type(array), intent(in) :: A
    real, intent(in) :: p
    double precision :: R
    call af_arr_pnorm(R, A%ptr, p, err)
  end function array_pnorm

  !> LU decomposition of array
  subroutine array_lu(L, U, p, A)
    type(array), intent(in) :: A
    type(array), intent(inout) :: L
    type(array), intent(inout) :: U
    type(array), intent(inout) :: p
    call init_eq(L, A)
    call init_eq(U, A)
    call init_eq(p, A)
    L%shape(2) = min(A%shape(1), A%shape(2))
    U%shape(1) = L%shape(2)
    p%shape(1) = 1
    p%shape(2) = A%shape(1)
    call af_arr_lu(L%ptr, U%ptr, p%ptr, A%ptr, err)
  end subroutine array_lu

  !> LU decomposition of array
  subroutine array_lu_inplace(A)
    type(array), intent(inout) :: A
    call af_arr_lu_inplace(A%ptr, err)
  end subroutine array_lu_inplace

  !> QR decomposition of array
  subroutine array_qr(Q, R, A)
    type(array), intent(in) :: A
    type(array), intent(inout) :: Q
    type(array), intent(inout) :: R
    call init_eq(Q, A)
    call init_eq(R, A)
    Q%shape(2) = A%shape(1)
    call af_arr_qr(Q%ptr, R%ptr, A%ptr, err)
  end subroutine array_qr

  !> Singular value decomposition of array
  subroutine array_singular(S, U, V, A)
    type(array), intent(in) :: A
    type(array), intent(inout) :: S, U, V
    call init_eq(S, A)
    call init_eq(U, A)
    call init_eq(V, A)
    U%shape(2) = U%shape(1)
    V%shape(1) = V%shape(2)
    call af_arr_singular(S%ptr, U%ptr, V%ptr, A%ptr, err)
  end subroutine array_singular

  !> cholesky decomposition of array
  subroutine array_cholesky(R, A)
    type(array), intent(in) :: A
    type(array), intent(inout) :: R
    call init_eq(R, A)
    call af_arr_cholesky(R%ptr, A%ptr, err)
  end subroutine array_cholesky

  !> cholesky decomposition of array
  subroutine array_cholesky_inplace(A)
    type(array), intent(inout) :: A
    call af_arr_cholesky_inplace(A%ptr, err)
  end subroutine array_cholesky_inplace

  !> Solve a system of equations
  function array_solve(A, B) result(X)
    type(array), intent(in) :: A, B
    type(array) :: X
    call init_eq(X, B)
    X%shape(1) = A%shape(2)
    call af_arr_solve(X%ptr, A%ptr, B%ptr, err)
  end function array_solve

  !> Inverse an array
  function array_inverse(A) result(R)
    type(array), intent(in) :: A
    type(array) :: R
    call init_eq(R, A)
    call af_arr_inverse(R%ptr, A%ptr, err)
  end function array_inverse

  !> Summation of elements in a matrix
  function array_sum (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_sum(R%ptr, A%ptr, dim, err)
  end function array_sum

  !> Product of elements in a matrix
  function array_product (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_product(R%ptr, A%ptr, dim, err)
  end function array_product

  !> Minimum of elements in a matrix
  function array_min (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_min(R%ptr, A%ptr, dim, err)
  end function array_min

  !> Maximum of elements in a matrix
  function array_max (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_max(R%ptr, A%ptr, dim, err)
  end function array_max

  !> Any of elements in a matrix
  function array_anytrue (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_anytrue(R%ptr, A%ptr, dim, err)
  end function array_anytrue

  !> All of elements in a matrix
  function array_alltrue(A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_alltrue(R%ptr, A%ptr, dim, err)
  end function array_alltrue

  !> Mean of elements in a matrix
  function array_mean (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_mean(R%ptr, A%ptr, dim, err)
  end function array_mean

  !> Standard deviation of elements in a matrix
  function array_std (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_stdev(R%ptr, A%ptr, dim, err)
  end function array_std

  !> Variance of elements in a matrix
  function array_var (A, d) result(R)
    type(array), intent(in) :: A
    integer, optional, intent(in) :: d
    type(array) :: R
    integer :: dim = 1
    if (present(d)) dim = d
    call init_eq(R, A)
    R%shape(1) = 1
    call af_arr_var(R%ptr, A%ptr, dim, err)
  end function array_var

  !> Show device info
  subroutine device_info_()
    call af_device_info()
  end subroutine device_info_

  !> Evaluate an expression on device
  subroutine device_eval_(A)
    type(array), intent(in) :: A
    call af_device_eval(A%ptr)
  end subroutine device_eval_

  !> Synchronize on the device
  subroutine device_sync_()
    call af_device_sync()
  end subroutine device_sync_

  !> Get the device number
  function device_get_() result(R)
    integer :: R
    call af_device_get(R)
  end function device_get_

  !> Set a particular device
  subroutine device_set_(R)
    integer :: R
    call af_device_set(R)
  end subroutine device_set_

  !> Show device count
  function device_count_() result(R)
    integer :: R
    call af_device_count(R)
  end function device_count_

  !> Show device start
  subroutine timer_start_()
    call af_timer_start()
  end subroutine timer_start_

  !> Show device stop
  function timer_stop_() result(elapsed)
    double precision :: elapsed
    call af_timer_stop(elapsed)
  end function timer_stop_

end module arrayfire
