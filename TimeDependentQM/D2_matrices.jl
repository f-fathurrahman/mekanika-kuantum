"""
Build second derivative matrix using 11-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (8*f[i-5] - 125*f[i-4] + 1000*f[i-3] - 6000*f[i-2] + 42000*f[i-1] - 73766*f[i+0] + 
       42000*f[i+1] - 6000*f[i+2] + 1000*f[i+3] - 125*f[i+4] + 8*f[i+5]) / (25200*h**2)
"""
function build_D2_matrix_11pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-5
        mat[i,i]   = -73766.0
        mat[i,i+1] =  42000.0
        mat[i,i+2] =  -6000.0
        mat[i,i+3] =   1000.0
        mat[i,i+4] =   -125.0
        mat[i,i+5] =      8.0
        mat[i+1,i] =  42000.0
        mat[i+2,i] =  -6000.0
        mat[i+3,i] =   1000.0
        mat[i+4,i] =   -125.0
        mat[i+5,i] =      8.0
    end

    mat[N-4,N-4] = -73766.0
    #
    mat[N-4,N-3] =  42000.0
    mat[N-4,N-2] =  -6000.0
    mat[N-4,N-1] =   1000.0
    mat[N-4,N  ] =   -125.0
    #
    mat[N-3,N-4] = 42000.0
    mat[N-2,N-4] = -6000.0
    mat[N-1,N-4] =  1000.0
    mat[N  ,N-4] =  -125.0


    mat[N-3,N-3] = -73766.0
    #
    mat[N-3,N-2] = 42000.0
    mat[N-3,N-1] = -6000.0
    mat[N-3,N  ] =  1000.0
    #
    mat[N-2,N-3] = 42000.0
    mat[N-1,N-3] = -6000.0
    mat[N  ,N-3] =  1000.0


    mat[N-2,N-2] = -73766.0
    #
    mat[N-2,N-1] =  42000.0
    mat[N-2,N  ] =  -6000.0
    #
    mat[N-1,N-2] =  42000.0
    mat[N  ,N-2] =  -6000.0


    mat[N-1,N-1] = -73766.0
    mat[N-1,N]   = 42000.0
    mat[N,N-1]   = 42000.0


    mat[N,N] = -73766.0
    
    return mat/(25200*h^2)
end

"""
Build second derivative matrix using 3-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points
"""
function build_D2_matrix_3pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-1
        mat[i,i] = -2.0
        mat[i,i+1] = 1.0
        mat[i+1,i] = 1.0
    end
    mat[N,N] = -2.0
    return mat/h^2
end

"""
Build second derivative matrix using 5-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (-1*f[i-2] + 16*f[i-1]- 30*f[i+0] + 16*f[i+1] - 1*f[i+2] )/(12*h^2)
"""
function build_D2_matrix_5pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-2
        mat[i,i] = -30.0
        mat[i,i+1] = 16.0
        mat[i,i+2] = -1.0
        mat[i+1,i] = 16.0
        mat[i+2,i] = -1.0
    end

    mat[N-1,N-1] = -30.0
    mat[N-1,N] = 16.0
    mat[N,N-1] = 16.0

    mat[N,N] = -30.0

    return mat/(12*h^2)
end


"""
Build second derivative matrix using 7-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (2*f[i-3] - 27*f[i-2] + 270*f[i-1] - 490*f[i+0] + 270*f[i+1] - 27*f[i+2] + 2*f[i+3] )/(180*h^2)
"""
function build_D2_matrix_7pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-3
        mat[i,i] = -490.0
        mat[i,i+1] = 270.0
        mat[i,i+2] = -27.0
        mat[i,i+3] = 2.0
        mat[i+1,i] = 270.0
        mat[i+2,i] = -27.0
        mat[i+3,i] = 2.0
    end

    mat[N-2,N-2] = -490.0
    #
    mat[N-2,N-1] = 270.0
    mat[N-2,N] = -27.0
    #
    mat[N-1,N-2] = 270.0
    mat[N,N-2] = -27.0

    mat[N-1,N-1] = -490.0
    mat[N-1,N] = 270.0
    mat[N,N-1] = 270.0

    mat[N,N] = -490.0
    return mat/(180.0*h^2)
end

"""
Build second derivative matrix using 9-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (-9*f[i-4] + 128*f[i-3] - 1008*f[i-2] + 8064*f[i-1] - 14350*f[i+0] + 
        8064*f[i+1]-1008*f[i+2]+128*f[i+3]-9*f[i+4])/(5040*h^2)
"""
function build_D2_matrix_9pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-4
        mat[i,i]   = -14350.0
        mat[i,i+1] =   8064.0
        mat[i,i+2] =  -1008.0
        mat[i,i+3] =    128.0
        mat[i,i+4] =     -9.0
        mat[i+1,i] =   8064.0
        mat[i+2,i] =  -1008.0
        mat[i+3,i] =    128.0
        mat[i+4,i] =     -9.0
    end

    mat[N-3,N-3] = -14350.0
    #
    mat[N-3,N-2] =  8064.0
    mat[N-3,N-1] = -1008.0
    mat[N-3,N]   =   128.0
    #
    mat[N-2,N-3] =  8064.0
    mat[N-1,N-3] = -1008.0
    mat[N,N-3]   =   128.0

    mat[N-2,N-2] = -14350.0
    mat[N-2,N-1] =  8064.0
    mat[N-2,N]   = -1008.0
    mat[N-1,N-2] =  8064.0
    mat[N,N-2]   = -1008.0

    mat[N-1,N-1] = -14350.0
    mat[N-1,N]   = 8064.0
    mat[N,N-1]   = 8064.0

    mat[N,N] = -14350.0
    
    return mat/(5040*h^2)
end

"""
Build second derivative matrix using 11-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (8*f[i-5] - 125*f[i-4] + 1000*f[i-3] - 6000*f[i-2] + 42000*f[i-1] - 73766*f[i+0] + 
       42000*f[i+1] - 6000*f[i+2] + 1000*f[i+3] - 125*f[i+4] + 8*f[i+5]) / (25200*h**2)

Periodic boundary condition.
"""
function build_D2_matrix_p_11pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-5
        mat[i,i]   = -73766.0
        mat[i,i+1] =  42000.0
        mat[i,i+2] =  -6000.0
        mat[i,i+3] =   1000.0
        mat[i,i+4] =   -125.0
        mat[i,i+5] =      8.0
        mat[i+1,i] =  42000.0
        mat[i+2,i] =  -6000.0
        mat[i+3,i] =   1000.0
        mat[i+4,i] =   -125.0
        mat[i+5,i] =      8.0
    end

    mat[N-4,N-4] = -73766.0
    #
    mat[N-4,N-3] =  42000.0
    mat[N-4,N-2] =  -6000.0
    mat[N-4,N-1] =   1000.0
    mat[N-4,N  ] =   -125.0
    #
    mat[N-3,N-4] = 42000.0
    mat[N-2,N-4] = -6000.0
    mat[N-1,N-4] =  1000.0
    mat[N  ,N-4] =  -125.0


    mat[N-3,N-3] = -73766.0
    #
    mat[N-3,N-2] = 42000.0
    mat[N-3,N-1] = -6000.0
    mat[N-3,N  ] =  1000.0
    #
    mat[N-2,N-3] = 42000.0
    mat[N-1,N-3] = -6000.0
    mat[N  ,N-3] =  1000.0


    mat[N-2,N-2] = -73766.0
    #
    mat[N-2,N-1] =  42000.0
    mat[N-2,N  ] =  -6000.0
    #
    mat[N-1,N-2] =  42000.0
    mat[N  ,N-2] =  -6000.0


    mat[N-1,N-1] = -73766.0
    mat[N-1,N]   = 42000.0
    mat[N,N-1]   = 42000.0


    mat[N,N] = -73766.0

    mat[1,N  ] = 42000.0
    mat[1,N-1] = -6000.0
    mat[1,N-2] =  1000.0
    mat[1,N-3] =  -125.0
    mat[1,N-4] =     8.0
    #
    mat[N  ,1] = 42000.0
    mat[N-1,1] = -6000.0
    mat[N-2,1] =  1000.0
    mat[N-3,1] =  -125.0
    mat[N-4,1] =     8.0

    mat[2,N  ] = -6000.0
    mat[2,N-1] =  1000.0
    mat[2,N-2] =  -125.0
    mat[2,N-3] =     8.0
    #
    mat[N  ,2] = -6000.0
    mat[N-1,2] =  1000.0
    mat[N-2,2] =  -125.0
    mat[N-3,2] =     8.0

    mat[3,N  ] =  1000.0
    mat[3,N-1] =  -125.0
    mat[3,N-2] =     8.0
    #
    mat[N  ,3] =  1000.0
    mat[N-1,3] =  -125.0
    mat[N-2,3] =     8.0

    mat[4,N  ] =  -125.0
    mat[4,N-1] =     8.0
    #
    mat[N  ,4] =  -125.0
    mat[N-1,4] =     8.0

    mat[5,N] = 8.0
    mat[N,5] = 8.0

    return mat/(25200*h^2)
end

"""
Build second derivative matrix using 3-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

Periodic boundary condition.
"""
function build_D2_matrix_p_3pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-1
        mat[i,i] = -2.0
        mat[i,i+1] = 1.0
        mat[i+1,i] = 1.0
    end
    
    mat[N,N] = -2.0

    mat[1,N] = 1.0
    mat[N,1] = 1.0
    
    return mat/h^2
end



"""
Build second derivative matrix using 5-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (-1*f[i-2] + 16*f[i-1]- 30*f[i+0] + 16*f[i+1] - 1*f[i+2] )/(12*h^2)

Periodic boundary condition.
"""
function build_D2_matrix_p_5pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-2
        mat[i,i] = -30.0
        mat[i,i+1] = 16.0
        mat[i,i+2] = -1.0
        mat[i+1,i] = 16.0
        mat[i+2,i] = -1.0
    end

    mat[N,N] = -30.0

    mat[N-1,N-1] = -30.0
    #
    mat[N-1,N] = 16.0
    mat[N,N-1] = 16.0
    #
    mat[1,N] = 16.0
    mat[1,N-1] = -1.0
    mat[N,1] = 16.0
    mat[N-1,1] = -1.0
    #
    mat[2,N] = -1.0
    mat[N,2] = -1.0

    return mat/(12*h^2)
end


"""
Build second derivative matrix using 7-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (2*f[i-3] - 27*f[i-2] + 270*f[i-1] - 490*f[i+0] + 270*f[i+1] - 27*f[i+2] + 2*f[i+3] )/(180*h^2)

Periodic boundary condition
"""
function build_D2_matrix_p_7pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-3
        mat[i,i] = -490.0
        mat[i,i+1] = 270.0
        mat[i,i+2] = -27.0
        mat[i,i+3] = 2.0
        mat[i+1,i] = 270.0
        mat[i+2,i] = -27.0
        mat[i+3,i] = 2.0
    end

    mat[N,N] = -490.0

    mat[N-2,N-2] = -490.0
    #
    mat[N-2,N-1] = 270.0
    mat[N-2,N] = -27.0
    #
    mat[N-1,N-2] = 270.0
    mat[N,N-2] = -27.0

    mat[N-1,N-1] = -490.0
    mat[N-1,N] = 270.0
    mat[N,N-1] = 270.0

    mat[1,N  ] = 270.0
    mat[1,N-1] = -27.0
    mat[1,N-2] = 2.0
    #
    mat[N  ,1] = 270.0
    mat[N-1,1] = -27.0
    mat[N-2,1] = 2.0

    mat[2,N  ] = -27.0
    mat[2,N-1] = 2.0
    #
    mat[N  ,2] = -27.0
    mat[N-1,2] = 2.0

    mat[3,N] = 2.0
    mat[N,3] = 2.0

    return mat/(180.0*h^2)
end



"""
Build second derivative matrix using 9-points centered
finite difference approximation.

# Arguments
- `N::Int64`: number of grid points
- `h::Float64`: spacing between grid points

f_xx = (-9*f[i-4] + 128*f[i-3] - 1008*f[i-2] + 8064*f[i-1] - 14350*f[i+0] + 
        8064*f[i+1]-1008*f[i+2]+128*f[i+3]-9*f[i+4])/(5040*h^2)

Periodic boundary condition
"""
function build_D2_matrix_p_9pt( N::Int64, h::Float64 )
    mat = zeros(Float64,N,N)
    for i = 1:N-4
        mat[i,i]   = -14350.0
        mat[i,i+1] =   8064.0
        mat[i,i+2] =  -1008.0
        mat[i,i+3] =    128.0
        mat[i,i+4] =     -9.0
        mat[i+1,i] =   8064.0
        mat[i+2,i] =  -1008.0
        mat[i+3,i] =    128.0
        mat[i+4,i] =     -9.0
    end

    mat[N-3,N-3] = -14350.0
    #
    mat[N-3,N-2] =  8064.0
    mat[N-3,N-1] = -1008.0
    mat[N-3,N]   =   128.0
    #
    mat[N-2,N-3] =  8064.0
    mat[N-1,N-3] = -1008.0
    mat[N,N-3]   =   128.0

    mat[N-2,N-2] = -14350.0
    mat[N-2,N-1] =  8064.0
    mat[N-2,N]   = -1008.0
    mat[N-1,N-2] =  8064.0
    mat[N,N-2]   = -1008.0

    mat[N-1,N-1] = -14350.0
    mat[N-1,N]   = 8064.0
    mat[N,N-1]   = 8064.0

    mat[N,N] = -14350.0

    mat[1,N]   =  8064.0
    mat[1,N-1] = -1008.0
    mat[1,N-2] =   128.0
    mat[1,N-3] =    -9.0
    #
    mat[N  ,1] =  8064.0
    mat[N-1,1] = -1008.0
    mat[N-2,1] =   128.0
    mat[N-3,1] =    -9.0

    mat[2,N  ] = -1008.0
    mat[2,N-1] =   128.0
    mat[2,N-2] =    -9.0
    #
    mat[N  ,2] = -1008.0
    mat[N-1,2] =   128.0
    mat[N-2,2] =    -9.0

    mat[3,N  ] = 128.0
    mat[3,N-1] =  -9.0
    #
    mat[N  ,3] = 128.0
    mat[N-1,3] =  -9.0

    mat[4,N] = -9.0
    mat[N,4] = -9.0

    return mat/(5040*h^2)
end

