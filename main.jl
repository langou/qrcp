using LinearAlgebra
using Printf
using Random
rng = MersenneTwister()

function larfg!(x)
   """Computes the Householder transformation for input vector x"""

#  (1) I am not sure documentation works
#
#  (2) we might want to break this routine so that this is not a vector but first
#  component and then the rest

   σ = dot(x[2:end],x[2:end])

   if σ == 0
       τ = 0
       return τ
   end

   β = sqrt(x[1]^2 + σ)

   if x[1] > 0
       x[1] = x[1] + β
       r = -β
   else
       x[1] = x[1] - β
       r = +β
   end

   x[2:end] = x[2:end] / x[1]

   τ = 2.0 / (1.0 + σ / x[1]^2)

   x[1] = r

   return τ

end

function mylarfg!(x1,x2,τ)
   """Computes the Householder transformation for input vector x"""

#  (1) I am not sure documentation works
#
#  (*) maybe change x1 by alpha and x2 by x
#
#  (2) we might want to break this routine so that this is not a vector but first
#  component and then the rest

   σ = dot(x2[1:end],x2[1:end])

   if σ == 0
       τ[1] = 0
       return ;
   end

   β = sqrt(x1[1]^2 + σ)

   if x1[1] > 0
       x1[1] = x1[1] + β
       r = -β
   else
       x1[1] = x1[1] - β
       r = +β
   end

   x2[1:end] = x2[1:end] / x1[1]

   τ[1] = 2.0 / (1.0 + σ / x1[1]^2)

   x1[1] = r

   return ;

end

function geqr2!( A, τ, w )

# TODO: trick with the τ

   m = size(A,1)
   n = size(A,2)

   kend = (m > n ? n : m-1)

   for k=1:kend

       A11 = view( A, k:k, k:k );
       A21 = view( A, k+1:m, k:k );
       A12 = view( A, k:k, k+1:n );
       A22 = view( A, k+1:m, k+1:n );
       w_  = view( w, 1:1, k+1:n );
       τ_1 = view( τ, k:k );

#      τ[k] = larfg!(A_1)
       mylarfg!(A11, A21, τ_1)

#      this is same as larf
       w_ = A12
       w_ += A21'*A22
       A12 .-=  τ_1 * w_
       A22 .-=  A21 * ( τ_1 * w_ )

   end

end

function org2r_v1!( A, τ, w, Q )

   m = size(A,1)
   n = size(A,2)

   kend = (m > n ? n : m-1)

   Q_ = view( Q, 1:m, 1:n );
   Q_ .= zeros(Float64,m,n)
   for k=m:-1:1
      Q_[k,k] = 1.0;
   end 

   for k=kend:-1:1

       A_1 = view( A, k:m, k:k );
       Q_ = view( Q, k:m, k:n );
       w_ = view( w, 1:1, k:n );

       a1 = A_1[1,1];
       A_1[1,1] = 1.0;

       w_ .= A_1' * Q_
       Q_ .-= A_1 * ( τ[k] * w_ )

       A_1[1,1] = a1;

   end

end


function org2r_v2!( A, τ, w, Q )

   m = size(A,1)
   n = size(A,2)

   if ( m > n )

       A_2 = view( A, n+1:m, n:n );
       Q_1 = view( Q, n:n, n:n );
       Q_2 = view( Q, n+1:m, n:n );

       Q_1 .= 1.0 - τ[n]
       Q_2 .= - A_2 * τ[n]

   else

       Q[m,m] = 1.0

   end

   for k=min(m,n)-1:-1:1

       A_2 = view( A, k+1:m, k:k );
       Q_11 = view( Q, k:k, k:k );
       Q_12 = view( Q, k:k, k+1:n );
       Q_21 = view( Q, k+1:m, k:k );
       Q_22 = view( Q, k+1:m, k+1:n );
       w_2 = view( w, 1:1, k+1:n );

       w_2 .= A_2' * Q_22
       w_2 .*= τ[k]

       Q_12 .= - w_2
       Q_22 .-= A_2 * w_2
       Q_11 .= 1.0 - τ[k]
       Q_21 .= - A_2 * τ[k]

   end

end

function org2r!( A, τ, w )

# TODO: trick with the τ
# TODO: full Q (only does thin Q as of now) - get sizes from τ and A

   m = size(A,1)
   n = size(A,2)

   if ( m > n )

       A_2 = view( A, n+1:m, n:n );
       A_1 = view( A, n:n, n:n );
       A_2 = view( A, n+1:m, n:n );

       A_1 .= 1.0 - τ[n]
       A_2 .= - A_2 * τ[n]

   else

       A[m,m] = 1.0

   end

   for k=min(m,n)-1:-1:1

       A_11 = view( A, k:k, k:k );
       A_12 = view( A, k:k, k+1:n );
       A_21 = view( A, k+1:m, k:k );
       A_22 = view( A, k+1:m, k+1:n );
       w_ = view( w, 1:1, k+1:n );

       w_ .= A_21' * A_22
       w_ .*= τ[k]

       A_12 .= - w_
       A_22 .-= A_21 * w_
       A_11 .= 1.0 - τ[k]
       A_21 .= - A_21 * τ[k]

   end

end


# m = 10; n =  4; A = randn(rng,Float64,m,n)
# m =  5; n =  5; A = randn(rng,Float64,m,n)
  m =  6; n = 10; A = randn(rng,Float64,m,n)

if (m >= n)

   Q1 = copy(A)
   R1 = zeros(Float64,n,n)
   τ1 = zeros(Float64,n)
   w = zeros(Float64,1,n)

   geqr2!(Q1,τ1,w)
   R1 = UpperTriangular(Q1[1:n,1:n])
   org2r!(Q1,τ1,w)

   @printf("|| I - Q1ᴴ Q1 ||₁             = %6.2e\n",
        opnorm( I(n) - Q1' * Q1, 1 ) )
   @printf("|| A - Q1 * R1 ||₁ / || A ||₁ = %6.2e\n",
        opnorm( A - Q1 * UpperTriangular(R1), 1) / opnorm(A,1))

   Q2 = copy(A)
   R2 = zeros(Float64,n,n)
   τ2 = zeros(Float64,n)

   LAPACK.geqrf!(Q2,τ2)
   R2 = UpperTriangular(Q2[1:n,1:n])
   LAPACK.orgqr!(Q2,τ2)

   @printf("|| I - Q2ᴴ Q2 ||₁             = %6.2e\n",
           opnorm( I(n) - Q2' * Q2, 1 ) )
   @printf("|| A - Q2 * R2 ||₁ / || A ||₁ = %6.2e\n",
        opnorm( A - Q2 * UpperTriangular(R2), 1) / opnorm(A,1))

else 

   R1 = copy(A)
   Q1 = zeros(Float64,m,m)
   τ1 = zeros(Float64,m)
   w = zeros(Float64,1,n)

   geqr2!(R1,τ1,w)
   Q1 .= R1[1:m,1:m]
   R1 .= triu(R1[1:m,1:n])
   org2r!(Q1,τ1,w)
#  LAPACK.orgqr!(Q1,τ1)

   @printf("|| I - Q1ᴴ Q1 ||₁             = %6.2e\n",
        opnorm( I(m) - Q1' * Q1, 1 ) )
   @printf("|| A - Q1 * R1 ||₁ / || A ||₁ = %6.2e\n",
        opnorm( A - Q1 * R1, 1) / opnorm(A,1))

   R2 = copy(A)
   Q2 = zeros(Float64,m,m)
   τ2 = zeros(Float64,m)

   LAPACK.geqrf!(R2,τ2)
   Q2 .= R2[1:m,1:m]
   R2 .= triu(R2[1:m,1:n])
   LAPACK.orgqr!(Q2,τ2)

   @printf("|| I - Q2ᴴ Q2 ||₁             = %6.2e\n",
           opnorm( I(m) - Q2' * Q2, 1 ) )
   @printf("|| A - Q2 * R2 ||₁ / || A ||₁ = %6.2e\n",
        opnorm( A - Q2 * R2, 1) / opnorm(A,1))

end

;
