program desitter
    implicit none

    ! Define precision
    integer, parameter :: dp = kind(1.0d0)

    ! Constants
    real(dp), parameter :: Lzero = 0.5_dp
    real(dp), parameter :: cf = 1.0_dp / (4.0_dp * 3.141592653589793_dp * Lzero**2)

    ! Variable declarations
    integer :: N, L, R, n_ea, n_dh
    integer :: ind_ea, ind_dh, i, n, l, r
    real(dp) :: ea, dh
    real(dp), allocatable :: EA_array(:), DH_array(:), z_list(:)
    real(dp), allocatable :: fR_array(:,:), fI_array(:,:), F_array(:,:)
    complex(dp), allocatable :: ResultadoIntegral(:,:)

    ! Read command-line arguments
    call get_command_argument(1, N)
    call get_command_argument(2, L)
    call get_command_argument(3, R)
    call get_command_argument(4, n_ea)
    call get_command_argument(5, n_dh)

    ! Allocate arrays based on input arguments
    allocate(EA_array(n_ea))
    allocate(DH_array(n_dh))
    allocate(z_list(100000))
    allocate(fR_array(n_ea, n_dh))
    allocate(fI_array(n_ea, n_dh))
    allocate(F_array(n_ea, n_dh))
    allocate(ResultadoIntegral(n_ea, n_dh))

    ! Initialize EA_array and DH_array
    do i = 1, n_ea
        EA_array(i) = -12.0_dp + (24.0_dp / (n_ea - 1)) * (i - 1)
    end do

    do i = 1, n_dh
        DH_array(i) = 0.175_dp + (0.625_dp / (n_dh - 1)) * (i - 1)
    end do

    ! Initialize arrays to zero
    fR_array = 0.0_dp
    fI_array = 0.0_dp
    F_array = 0.0_dp
    ResultadoIntegral = (0.0_dp, 0.0_dp)

    ! Main loop over EA and DH
    do ind_ea = 1, n_ea
        ea = EA_array(ind_ea)
        do ind_dh = 1, n_dh
            dh = DH_array(ind_dh)

            ! Create z_list for integration
            do i = 1, 100000
                z_list(i) = -dh + (2.0_dp * dh / 99999) * (i - 1)
            end do

            ! Nested loop over n, l, r
            do n = 0, N
                do l = 0, L
                    do r = 0, R
                        if (n == 0 .and. l == 0 .and. r == 0) cycle

                        ! Perform trapezoidal integration
                        fR_array(ind_ea, ind_dh) = integrate(z_list, fReal(z_list, ea, dh, n, l, r))
                        fI_array(ind_ea, ind_dh) = integrate(z_list, fImag(z_list, ea, dh, n, l, r))

                        ResultadoIntegral(ind_ea, ind_dh) = &
                            ResultadoIntegral(ind_ea, ind_dh) + &
                            cmplx(fR_array(ind_ea, ind_dh), fI_array(ind_ea, ind_dh), dp)

                        F_array(ind_ea, ind_dh) = F_array(ind_ea, ind_dh) + &
                            abs(ResultadoIntegral(ind_ea, ind_dh))**2 * cf
                    end do
                end do
            end do

            print *, ind_ea
        end do
    end do

    ! Write results to file
    open(unit=10, file="result.txt")
    do i = 1, n_ea
        write(10, "(100(F20.10,1X))") F_array(i, :)
    end do
    close(10)

    print *, "Results saved to result.txt"

contains

    real(dp) function switch_func(z, DH)
        real(dp), intent(in) :: z, DH
        switch_func = cos(3.141592653589793_dp * z / (2.0_dp * DH))**4
    end function switch_func

    real(dp) function k(n, l, r)
        integer, intent(in) :: n, l, r
        k = sqrt(dble(n * n + l * l + r * r))
    end function k

    real(dp) function nonunit(z, n, l, r)
        real(dp), intent(in) :: z
        integer, intent(in) :: n, l, r
        nonunit = exp(-z) / sqrt(k(n, l, r))
    end function nonunit

    real(dp) function RealDyn(z, EA, n, l, r)
        real(dp), intent(in) :: z, EA
        integer, intent(in) :: n, l, r
        real(dp) :: k_val, exp_val
        k_val = k(n, l, r)
        exp_val = exp(-z)
        RealDyn = cos(z * EA) * cos(2.0_dp * 3.141592653589793_dp * k_val * exp_val / Lzero) + &
                  sin(z * EA) * sin(2.0_dp * 3.141592653589793_dp * k_val * exp_val / Lzero)
    end function RealDyn

    real(dp) function ImDyn(z, EA, n, l, r)
        real(dp), intent(in) :: z, EA
        integer, intent(in) :: n, l, r
        real(dp) :: k_val, exp_val
        k_val = k(n, l, r)
        exp_val = exp(-z)
        ImDyn = cos(z * EA) * sin(2.0_dp * 3.141592653589793_dp * k_val * exp_val / Lzero) - &
                sin(z * EA) * cos(2.0_dp * 3.141592653589793_dp * k_val * exp_val / Lzero)
    end function ImDyn

    real(dp) function fReal(z, EA, DH, n, l, r)
        real(dp), intent(in) :: z(:), EA, DH
        integer, intent(in) :: n, l, r
        integer :: i, sz
        sz = size(z)
        real(dp) :: val(sz)
        do i = 1, sz
            val(i) = switch_func(z(i), DH) * nonunit(z(i), n, l, r) * RealDyn(z(i), EA, n, l, r)
        end do
        fReal = val
    end function fReal

    real(dp) function fImag(z, EA, DH, n, l, r)
        real(dp), intent(in) :: z(:), EA, DH
        integer, intent(in) :: n, l, r
        integer :: i, sz
        sz = size(z)
        real(dp) :: val(sz)
        do i = 1, sz
            val(i) = switch_func(z(i), DH) * nonunit(z(i), n, l, r) * ImDyn(z(i), EA, n, l, r)
        end do
        fImag = val
    end function fImag

    real(dp) function integrate(x, y)
        real(dp), intent(in) :: x(:), y(:)
        integer :: i, sz
        real(dp) :: sum
        sz = size(x)
        sum = 0.0_dp
        do i = 2, sz
            sum = sum + 0.5_dp * (x(i) - x(i - 1)) * (y(i) + y(i - 1))
        end do
        integrate = sum
    end function integrate

end program desitter
