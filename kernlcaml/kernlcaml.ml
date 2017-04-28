(* Copyright (c) 2017 by Kernelyze LLC                                       *)
(* Author: Thomas A. Knox                                                    *)
(*                                                                           *)
(* This program is free software: you can redistribute it and/or modify      *)
(* it under the terms of the GNU Affero General Public License as            *)
(* published by the Free Software Foundation, either version 3 of the        *)
(* License, or (at your option) any later version.                           *)
(*                                                                           *)
(* This program is distributed in the hope that it will be useful,           *)
(* but WITHOUT ANY WARRANTY; without even the implied warranty of            *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *)
(* GNU Affero General Public License for more details.                       *)
(*                                                                           *)
(* You should have received a copy of the GNU Affero General Public License  *)
(* along with this program.  If not, see <http://www.gnu.org/licenses/>.     *)

open Bigarray

type vec = (float, float64_elt, fortran_layout) Array1.t
type mat = (float, float64_elt, fortran_layout) Array2.t

(* zeroin is the simple interface to a Brent zero solver.  It takes: 
    a one-variable function
    a lower bound on the interval to search
    an upper bound on the interval to search
    a numerical tolerance

    It returns a value X for which f(X) is nearly zero.
*)
external zeroin : (float -> float) -> float -> float -> float -> float = "caml_zeroin"
(* brent_zero is a richer interface to a Brent zero solver.  It takes:
    a one-variable function
    a lower bound on the interval to search
    an upper bound on the interval to search
    a numerical tolerance
    a maximum number of iterations to perform

    It returns a tuple:
    ( the X-coordinate of the zero 
      the value of the function at the supposed zero
      the number of iterations performed
      an error code )
*)
external brent_zero : (float -> float) -> float -> float -> float -> int -> (float * float * int * int) = "caml_brent_zero"
(* zeroin is the simple interface to a Brent minimizer.  It takes: 
    a one-variable function
    a lower bound on the interval to search
    an upper bound on the interval to search
    a numerical tolerance

    It returns a tuple:
    ( value X for which f(X) is nearly minimal
      the supposed minimal value f(X)
      an error code )
*)
external f_min : (float -> float) -> float -> float -> float -> (float * float * int) = "caml_f_min"
(* brent_min is a richer interface to a Brent minimizer.  It takes:
    a one-variable function
    a lower bound on the interval to search
    an upper bound on the interval to search
    a numerical tolerance
    a maximum number of iterations to perform

    It returns a tuple:
    ( the X for which f(X) is nearly minimal
      the value of the function at this X
      the number of iterations performed
      an error code )
*)
external brent_min : (float -> float) -> float -> float -> float -> int -> (float * float * int * int) = "caml_brent_min"
(* find_all_zeros seeks zeros on each interval defined by consecutive points of a grid.
   It takes:
    a one-variable function
    a grid
    a numerical tolerance

    It returns a vector of X such that f(X) is nearly zero, as determined by the tolerance.
*)
external find_all_zeros : (float -> float) -> vec -> float -> vec = "caml_find_all_zeros"
(* find_rel_optima finds relative optima on each interval defined by consecutive points of a grid.
   It takes:
    a one-variable function
    a grid
    a numerical tolerance 
    
    It returns a 2-column matrix: the first column gives the location of each relative
    optimum and the second gives the corresponding negative absolute value of the objective function. *)
external find_rel_optima : (float -> float) -> vec -> float -> mat = "caml_find_rel_optima"
(* cheb_pts takes an integer, which must be positive, and returns the Chebyshev points
   array of that size (as used extensively in approximation) *)
external cheb_pts : int -> vec = "caml_cheb_pts"
(* matmul takes conforming matrices A and B and returns A * B *)
external matmul : mat -> mat -> mat = "caml_matmul"
(* linear_solve takes a matrix A and a matrix (with the same number of rows as A, but any 
   number of columns) B and computes X such that A X = B *)
external linear_solve : mat -> mat -> mat = "caml_linear_solve"
(* mp_inverse takes a matrix and a numerical tolerance that determines which of
   the matrix's singular values are regarded as negligible.  It returns the
   Moore-Penrose pseudo-inverse of the matrix. *)
external mp_inverse : mat -> float -> mat = "caml_mp_inverse"
(* black_formula implements option pricing under a Black (1976) model; it takes:
    a flag: if true, the option is a call, else it is a put
    a strike price (must be positive)
    a forward price (must be positive)
    a volatility (in log terms, so a number like 0.10 for 10%, since this is a Black model)
    a discount factor

    It returns an option price.
*)
external black_formula : bool -> float -> float -> float -> float -> float = "caml_black_formula"
(* normopt_formula implements option pricing under a normal model, a la Bachelier; it takes:
    a flag: if true, the option is a call, else it is a put
    a strike price (may be negative)
    a forward price (may be negative)
    a volatility (in absolute terms, as this is a normal model not a Black model)
    a discount factor

    It returns an option price.
*)
external normopt_formula : bool -> float -> float -> float -> float -> float = "caml_normopt_formula"
(* borsuk_lower_bound takes:
    a function of two variables
    a rho vector
    a flag (is this optimization over X?  yes, if flag is true; otherwise, it is over Y)
    a numerical tolerance
    a maximum number of iterations
    
    It returns a tuple:
    ( Array: 3 columns, coefficients, then nodes, then errors at nodes
      Discrepancy (convergence measure, should be near zero)
      Number of iters
      Error code )
    *)
external borsuk_lower_bound : 
    (float -> float -> float) -> vec -> bool -> float -> int -> vec -> (mat * float * int * int) 
    = "caml_borsuk_lower_bound_bytecode" "caml_borsuk_lower_bound_native"
(* worst_rho_and_gamma takes:
    a function of two variables
    a lower bound on X
    an upper bound on X
    a lower bound on Y
    an upper bound on Y
    a number of terms
    a numerical tolerance
    a maximum number of iterations
    
    It returns a tuple:
    ( Rho array: 7 columns, rho, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes
      Gamma array: similar to rho array but with gamma quantities
      Discrepancy (convergence measure, should be near zero)
      Number of iters
      Kernel discrepancy (another convergency measure)
      Number of kernel iters
      Error code
      A matrix for rho
      A matrix for gamma ) *)
external worst_rho_and_gamma :
    (float -> float -> float) -> float -> float -> float -> float -> int -> float -> int -> 
        (mat * mat * float * int * float * int * int * mat * mat)
    = "caml_worst_rho_and_gamma_bytecode" "caml_worst_rho_and_gamma_native"
(* kernel_integral_eval takes:
    a function of two variables
    a flag: integral (really a weighted sum) is over X if true, else is over Y
    an array of integral points (the coordinates at which to evalue the variable summed over)
    an array of integral weights
    an array of evaluation points (for the variable not summed over)

    It returns an array of results, one result for each evaluation point.
*)
external kernel_integral_eval :
    (float -> float -> float) -> bool -> vec -> vec -> vec -> vec = "caml_kernel_integral_eval"
(* taylor_rankn_params takes:
    a function of two variables
    a rank (the rank of the desired Taylor approximation, equal to one plus the highest power of X or Y appearing)
    an X coordinate to center the Taylor series on
    a Y coordinate to center the Taylor series on
    a finite-difference step size; this should be larger for higher-order approximations to avoid numerical issues

    It returns a tuple: two matrix factors of the matrix of partial derivatives (appropriately
    scaled by factorials) and the error code.
*)
external taylor_rankn_params:
    (float -> float -> float) -> int -> float -> float -> float -> (mat * mat * int) 
    = "caml_taylor_rankn_params"
(* taylor_rankn_eval takes:
    a rank (the rank of the Taylor approximation, equal to one plus the highest power of X or Y appearing)
    the X coordinate that is the center of the Taylor series
    the Y coordinate that is the center of the Taylor series
    the V matrix factor of the Taylor series matrix (factor of scaled partials matrix)
    the W matrix factor of the Taylor series matrix (factor of scaled partials matrix)
    an array of X coordinates of the evaluation points
    an array of Y coordinates of the evaluation points

    It returns a tuple: a vector of evaluation results and the error code.
*)
external taylor_rankn_eval:
    int -> float -> float -> mat -> mat -> vec -> vec -> (vec * int) 
    = "caml_taylor_rankn_eval_bytecode" "caml_taylor_rankn_eval_native"
(* taylor_rankn_integral_eval takes:
    a rank (the rank of the Taylor approximation, equal to one plus the highest power of X or Y appearing)
    the X coordinate that is the center of the Taylor series
    the Y coordinate that is the center of the Taylor series
    the V matrix factor of the Taylor series matrix (factor of scaled partials matrix)
    the W matrix factor of the Taylor series matrix (factor of scaled partials matrix)
    a flag: the integral (really a weighted sum) is over X if true, else it is over Y
    an array of integral points (values of the coordinate being summed over)
    an array of integral (really sum) weights
    an array of evaluation points (values of the coordinate not being summed over)

    It returns a tuple: a vector of evaluation results and the error code.
*)
external taylor_rankn_integral_eval:
    int -> float -> float -> mat -> mat -> bool -> vec -> vec -> vec -> (vec * int) 
    = "caml_taylor_rankn_integral_eval_bytecode" "caml_taylor_rankn_integral_eval_native"
(* taylor_rankn_integral_coeffs takes:
    a rank (the rank of the Taylor approximation, equal to one plus the highest power of X or Y appearing)
    the X coordinate that is the center of the Taylor series
    the Y coordinate that is the center of the Taylor series
    the V matrix factor of the Taylor series matrix (factor of scaled partials matrix)
    the W matrix factor of the Taylor series matrix (factor of scaled partials matrix)
    a flag: the integral (really a weighted sum) is over X if true, else it is over Y
    an array of integral points (values of the coordinate being summed over)
    an array of integral (really sum) weights

    It returns a tuple: a vector of coefficients and the error code.
*)
external taylor_rankn_integral_coeffs:
    int -> float -> float -> mat -> mat -> bool -> vec -> vec -> (vec * int) 
    = "caml_taylor_rankn_integral_coeffs_bytecode" "caml_taylor_rankn_integral_coeffs_native"
(* eigen_rankn_params computes a Nystrom-method truncated singular function series; it takes:
    a function of two variables
    a rank (the rank of the desired approximation)
    an array of X coordinate points for quadrature
    an array of Y coordinate points for quadrature
    an array of X coordinate weights for quadrature (must be positive)
    an array of Y coordinate weights for quadrature (must be positive)
    
    It returns a tuple: two matrices used in approximating singular functions, 
    the array of singular values 1 through rank, and the error code.
*)
external eigen_rankn_params:
    (float -> float -> float) -> int -> vec -> vec -> vec -> vec -> (mat * mat * vec * int) 
    = "caml_eigen_rankn_params_bytecode" "caml_eigen_rankn_params_native"
(* eigen_rankn_eval evaluates a Nystrom-method truncated singular function series; it takes:
    a function of two variables
    a rank (the rank of the approximation)
    an array of X coordinate points for quadrature
    an array of Y coordinate points for quadrature
    an matrix of X function weights from eigen_rankn_params
    an matrix of Y function weights from eigen_rankn_params
    an array of singular values (the first "rank" of them)
    an array of X coordinate evaluation points
    an array of Y coordinate evaluation points
    
    It returns a tuple: an array of eval results and the error code.
*)
external eigen_rankn_eval:
    (float -> float -> float) -> int -> vec -> vec -> mat -> mat -> vec -> vec -> vec -> (vec * int) 
    = "caml_eigen_rankn_eval_bytecode" "caml_eigen_rankn_eval_native"
(* eigen_rankn_integral_eval takes:
    a function of two variables
    a rank (the rank of the approximation)
    an array of X coordinate points for quadrature
    an array of Y coordinate points for quadrature
    an matrix of X function weights from eigen_rankn_params
    an matrix of Y function weights from eigen_rankn_params
    an array of singular values (the first "rank" of them)
    a flag: the integral (really a weighted sum) is over X if true, else it is over Y
    an array of integral points (values of the coordinate being summed over)
    an array of integral (really sum) weights
    an array of evaluation points (values of the coordinate not being summed over)

    It returns a tuple: a vector of evaluation results and the error code.
*)
external eigen_rankn_integral_eval:
    (float -> float -> float) -> int -> vec -> vec -> mat -> mat -> vec -> bool -> vec -> vec -> vec -> (vec * int) 
    = "caml_eigen_rankn_integral_eval_bytecode" "caml_eigen_rankn_integral_eval_native"
(* eigen_rankn_integral_coeffs takes:
    a function of two variables
    a rank (the rank of the approximation)
    an array of X coordinate points for quadrature
    an array of Y coordinate points for quadrature
    an matrix of X function weights from eigen_rankn_params
    an matrix of Y function weights from eigen_rankn_params
    an array of singular values (the first "rank" of them)
    a flag: the integral (really a weighted sum) is over X if true, else it is over Y
    an array of integral points (values of the coordinate being summed over)
    an array of integral (really sum) weights

    It returns a tuple: a vector of coefficients and the error code.
*)
external eigen_rankn_integral_coeffs:
    (float -> float -> float) -> int -> vec -> vec -> mat -> mat -> vec -> bool -> vec -> vec -> (vec * int) 
    = "caml_eigen_rankn_integral_coeffs_bytecode" "caml_eigen_rankn_integral_coeffs_native"
(* num_opt_rankn_params computes a numerically-optimal approximation under uniform error; it takes:
    a function of two variables
    a lower bound on X
    an upper bound on X
    a lower bound on Y
    an upper bound on Y
    a rank (of the desired approximation)
    a numerical tolerance
    a maximum number of iterations
    
    It returns a tuple:
    ( Rho array
      Gamma array
      V matrix
      W matrix
      The Borsuk lower bound on approximation error (a scalar)
      Error code ) *)
external num_opt_rankn_params:
    (float -> float -> float) -> float -> float -> float -> float -> int -> float -> int -> 
    (vec * vec * mat * mat * float * int) 
    = "caml_num_opt_rankn_params_bytecode" "caml_num_opt_rankn_params_native"
(* num_opt_rankn_eval evaluates a numerically-optimal approximation under uniform error; it takes:
    a function of two variables
    Rho array
    Gamma array
    V matrix
    W matrix
    an array of X coordinate evaluation points
    an array of Y coordinate evaluation points
    
    It returns a tuple: an array of eval results and the error code.
*)
external num_opt_rankn_eval:
    (float -> float -> float) -> vec -> vec -> mat -> mat -> vec -> vec -> (vec * int) 
    = "caml_num_opt_rankn_eval_bytecode" "caml_num_opt_rankn_eval_native"
(* num_opt_rankn_integral_eval takes:
    a function of two variables
    Rho array
    Gamma array
    V matrix
    W matrix
    a flag: the integral (really a weighted sum) is over X if true, else it is over Y
    an array of integral points (values of the coordinate being summed over)
    an array of integral (really sum) weights
    an array of evaluation points (values of the coordinate not being summed over)

    It returns a tuple: a vector of evaluation results and the error code.
*)
external num_opt_rankn_integral_eval:
    (float -> float -> float) -> vec -> vec -> mat -> mat -> bool -> vec -> vec -> vec -> (vec * int) 
    = "caml_num_opt_rankn_integral_eval_bytecode" "caml_num_opt_rankn_integral_eval_native"
(* num_opt_rankn_integral_coeffs takes:
    a function of two variables
    Rho array
    Gamma array
    V matrix
    W matrix
    a flag: the integral (really a weighted sum) is over X if true, else it is over Y
    an array of integral points (values of the coordinate being summed over)
    an array of integral (really sum) weights

    It returns a tuple: a vector of coefficients and the error code.
*)
external num_opt_rankn_integral_coeffs:
    (float -> float -> float) -> vec -> vec -> mat -> mat -> bool -> vec -> vec -> (vec * int) 
    = "caml_num_opt_rankn_integral_coeffs_bytecode" "caml_num_opt_rankn_integral_coeffs_native"

(* Transposition: should consider putting this on the Fortran side *)
let transpose m = 
    let mt = Array2.create float64 fortran_layout (Array2.dim1 m) (Array2.dim2 m) in
     let _ = 
     (  for i = 1 to (Array2.dim1 m) do
            for j = 1 to (Array2.dim2 m) do
                mt.{i,j} <- m.{j,i}
            done;
        done ) in
     mt

(* Get an empty vector *)
let create_vec n = 
    Array1.create float64 fortran_layout n

(* Get an empty matrix *)
let create_mat m n =
    Array2.create float64 fortran_layout m n

(* Obtain a constant vector *)
let constant_vec n fill_with = 
    let cv = Array1.create float64 fortran_layout n in
    let _ = (
        for i = 1 to n do
            cv.{i} <- fill_with
        done) in
    cv

(* Obtain a constant matrix *)
let constant_mat m n fill_with = 
    let cm = Array2.create float64 fortran_layout m n in
    let _ = (
        for i = 1 to m do
            for j = 1 to n do
                cm.{i, j} <- fill_with
            done;
        done) in
    cm

(* Get an equal-spaced grid *)
let equal_spaced_grid n lb ub =
    let g = Array1.create float64 fortran_layout n in
    let _ = ( 
        for i = 1 to n do 
            g.{i} <- (lb +. (ub -. lb) *. ((float_of_int (i - 1)) /. (float_of_int (n - 1))))
        done ) in
    g

(* Obtain a matrix that is suitable for the X-argument in building a surface *)
let x_grid_for_surf m n xlb xub =
    let xg = Array2.create float64 fortran_layout m n in
    let _ = (
        for i = 1 to m do
            for j = 1 to n do
                xg.{i, j} <- (xlb +. (xub -. xlb) *. ((float_of_int (i - 1)) /. (float_of_int (m - 1))))
            done;
        done) in
    xg

(* Obtain a matrix that is suitable for the Y-argument in building a surface *)
let y_grid_for_surf m n ylb yub =
    let yg = Array2.create float64 fortran_layout m n in
    let _ = (
        for i = 1 to m do
            for j = 1 to n do
                yg.{i, j} <- (ylb +. (yub -. ylb) *. ((float_of_int (j - 1)) /. (float_of_int (n - 1))))
            done;
        done) in
    yg

(* One sometimes desires Chebyshev points translated and scaled to an interval different from [-1, 1] *)
let cheb_pts_trans n a b = 
    let cpts = cheb_pts n in
    let _ = (
        for i = 1 to n do
            cpts.{i} <- ( (0.5 *. (b +. a)) +. (0.5 *. (b -. a)) *. cpts.{i} )
        done ) in
    cpts

(* Some crude printing routines for 1-d and 2-d Bigarrays *)
let print_vec v =
    for i = 1 to (Array1.dim v) do
        Printf.printf "%f  %!" v.{i}
    done;
    Printf.printf "\n%!"

let print_mat m = 
    for i = 1 to (Array2.dim1 m) do
        for j = 1 to (Array2.dim2 m) do
            Printf.printf "%f  %!" m.{i,j}
        done;
        Printf.printf "\n%!";
    done