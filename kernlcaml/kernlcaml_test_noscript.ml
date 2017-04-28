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
open Kernlcaml

let test_zeros_optima () = 
    let testgrid = equal_spaced_grid 21 (-10.) 10. in
    let mysin x = sin x in
    let testzeros = find_all_zeros mysin testgrid 0.00000001 in
    let sizeof_zeros = Array1.dim testzeros in
    let _ = (for i = 1 to sizeof_zeros do
        Printf.printf "Zero %d is %f\n" i testzeros.{i}
    done) in
    let testopts = find_rel_optima mysin testzeros 0.000000001 in
    let sizeof_opts = Array2.dim1 testopts in
    let _ = (for i = 1 to sizeof_opts do
        Printf.printf "Rel optimum %d is %f with neg abs objective value %f \n" i testopts.{i,1} testopts.{i,2}
    done) in
    ()

let test_linear_algebra () = 
    let amat = Array2.create float64 fortran_layout 2 2 in
    let _ = (Printf.printf "The A matrix is:\n";
            for i = 1 to 2 do
                for j = 1 to 2 do
                    Array2.set amat i j (sin ( (float_of_int i) -. (float_of_int j) ));
                done;
            done;
            print_mat amat; ) in
    let bmat = Array2.create float64 fortran_layout 2 1 in
    let _ = (Printf.printf "The B matrix is:\n";
            for i = 1 to 2 do
                for j = 1 to 1 do
                    Array2.set bmat i j (cos ( (float_of_int i) -. (float_of_int j) ));
                done;
            done;
            print_mat bmat; ) in
    let cmat = matmul amat bmat in
    let _ = (Printf.printf "The product C = A * B is:\n";
            print_mat cmat; ) in
    let xmat = linear_solve amat cmat in
    let _ = (Printf.printf "The solution X to A * X = C is:\n";
            print_mat xmat; ) in
    let _ = (Printf.printf "Check the solution; this should be zero:\n";
            for i = 1 to 2 do
                for j = 1 to 1 do
                    Printf.printf "%f  " (bmat.{i,j} -. xmat.{i,j})
                done;
                Printf.printf "\n";
            done) in
    ()

let test_borsuk_lb () =
    let grid = cheb_pts 40 in
    let rho = cheb_pts 4 in
    let myexp x y = exp (x *. y) in
    let myinv x y = 1.0 /. (1.0 +. (x *. y) /. 2.0) in
    let mygauss x y = exp ( (y -. x) *. (x -. y) ) in
    let (expres, expdiscr, expniter, experr) = borsuk_lower_bound myexp rho true 1e-14 100 grid in
    let (invres, invdiscr, invniter, inverr) = borsuk_lower_bound myinv rho true 1e-14 100 grid in
    let (gaussres, gaussdiscr, gaussniter, gausserr) = borsuk_lower_bound mygauss rho true 1e-14 100 grid in
    let _ = (
        Printf.printf "Borsuk lower bounds test.\n";
        Printf.printf "Exponential kernel result:\n";
        Printf.printf "Discrepancy: %f Number of iters: %d Error code: %d\n" expdiscr expniter experr;
        Printf.printf "First column is coefficients, second is nodes, third is errors at nodes.\n";
        print_mat expres; 
        Printf.printf "1 / ( 1 + xy/2 ) kernel result:\n";
        Printf.printf "Discrepancy: %f Number of iters: %d Error code: %d\n" invdiscr invniter inverr;
        Printf.printf "First column is coefficients, second is nodes, third is errors at nodes.\n";
        print_mat invres; 
        Printf.printf "Gaussian kernel result:\n";
        Printf.printf "Discrepancy: %f Number of iters: %d Error code: %d\n" gaussdiscr gaussniter gausserr;
        Printf.printf "First column is coefficients, second is nodes, third is errors at nodes.\n";
        print_mat gaussres) in
    ()

let test_worst_rho_gamma() =
    let myexp x y = exp (x *. y) in
    let myinv x y = 1.0 /. (1.0 +. (x *. y) /. 2.0) in
    let mygauss x y = exp ( (y -. x) *. (x -. y) ) in
    let (erhoarr, egammaarr, ediscr, eniter, ekdiscr, eknitr, eerrcode, eamat_rho, eamat_gamma) 
        = worst_rho_and_gamma myexp (-1.) 1.5 (-1.2) 1. 4 1e-15 200 in
    let (lrhoarr, lgammaarr, ldiscr, lniter, lkdiscr, lknitr, lerrcode, lamat_rho, lamat_gamma) 
        = worst_rho_and_gamma myinv (-1.) 1.5 (-1.2) 1. 4 1e-15 200 in
    let (grhoarr, ggammaarr, gdiscr, gniter, gkdiscr, gknitr, gerrcode, gamat_rho, gamat_gamma) 
        = worst_rho_and_gamma mygauss (-1.) 1.5 (-1.2) 1. 4 1e-15 200 in
    let _ = (
        Printf.printf "Worst rho and gamma test.\n";
        Printf.printf "Exponential kernel result:\n";
        Printf.printf "Discrepancy: %f Number of iters: %d Error code: %d\n" ediscr eniter eerrcode;
        Printf.printf "Rho, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes:\n";
        print_mat erhoarr;
        Printf.printf "Gamma, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes:\n";
        print_mat egammaarr; 
        Printf.printf "1 / ( 1 + xy/2 ) kernel result:\n";
        Printf.printf "Discrepancy: %f Number of iters: %d Error code: %d\n" ldiscr lniter lerrcode;
        Printf.printf "Rho, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes:\n";
        print_mat lrhoarr;
        Printf.printf "Gamma, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes:\n";
        print_mat lgammaarr; 
        Printf.printf "Gaussian kernel result:\n";
        Printf.printf "Discrepancy: %f Number of iters: %d Error code: %d\n" gdiscr gniter gerrcode;
        Printf.printf "Rho, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes:\n";
        print_mat grhoarr;
        Printf.printf "Gamma, coeffs, nodes, errors at nodes, Borsuk coeffs, Borsuk nodes, Borsuk errors at nodes:\n";
        print_mat ggammaarr;) in
    ()

let test_taylor_ker f rank xcenter ycenter fin_diff xeval yeval overx intpts intwts intevalpts =
    let (vmat, wmat, errcode) = taylor_rankn_params f rank xcenter ycenter fin_diff in
    let (evalres, evalerr) = taylor_rankn_eval rank xcenter ycenter vmat wmat xeval yeval in
    let (intevalres, intevalerr) = taylor_rankn_integral_eval rank xcenter ycenter vmat wmat overx intpts intwts intevalpts in
    let (intcoeffs, coeffserr) = taylor_rankn_integral_coeffs rank xcenter ycenter vmat wmat overx intpts intwts in
    let _ = (
        Printf.printf "The V matrix:\n";
        print_mat vmat;
        Printf.printf "The W matrix:\n";
        print_mat wmat;
        Printf.printf "Error code in Taylor parameter calculation:\n";
        Printf.printf "%d\n" errcode;
        Printf.printf "Evaluation:\n";
        Printf.printf "The X eval points:\n";
        print_vec xeval;
        Printf.printf "The Y eval points:\n";
        print_vec yeval;
        Printf.printf "The eval results:\n";
        print_vec evalres;
        Printf.printf "Error code in evaluation:\n";
        Printf.printf "%d\n" evalerr;
        Printf.printf "Integral (really sum) evaluation:\n";
        Printf.printf "Is integral over X?\n";
        if overx then
            Printf.printf "Yes\n"
        else
            Printf.printf "No\n";
        Printf.printf "The points over which the sum is taken:\n";
        print_vec intpts;
        Printf.printf "The weights of the sum:\n";
        print_vec intwts;
        Printf.printf "The points at which the integral is evaluated:\n";
        print_vec intevalpts;
        Printf.printf "The results of the integral evaluation:\n";
        print_vec intevalres;
        Printf.printf "Error code in integral evaluation:\n";
        Printf.printf "%d\n" intevalerr;
        Printf.printf "The coefficients of the rank-n integral:\n";
        print_vec intcoeffs;
        Printf.printf "Error code in integral coefficient calculation:\n";
        Printf.printf "%d\n" coeffserr;
    ) in
    ()

let test_taylor () =
    let myexp x y = exp (x *. y) in
    let myinv x y = 1.0 /. (1.0 +. (x *. y) /. 2.0) in
    let mygauss x y = exp ( (y -. x) *. (x -. y) ) in
    let xeval = cheb_pts 10 in
    let yeval = cheb_pts 10 in
    let overx = true in
    let intpts = cheb_pts 20 in
    let intwts = constant_vec (Array1.dim intpts) (1.0 /. (float_of_int (Array1.dim intpts) ) ) in
    let intevalpts = cheb_pts 7 in
    let _ = (
        Printf.printf "Taylor series approximation test.\n";
        Printf.printf "Exponential kernel result:\n";
        test_taylor_ker myexp 3 0.0 0.0 0.01 xeval yeval overx intpts intwts intevalpts;
        Printf.printf "1 / (1 + x * y / 2) kernel result:\n";
        test_taylor_ker myinv 3 0.0 0.0 0.01 xeval yeval overx intpts intwts intevalpts;
        Printf.printf "Gaussian kernel result:\n";
        test_taylor_ker mygauss 3 0.0 0.0 0.01 xeval yeval overx intpts intwts intevalpts;
        ) in
    ()

let test_eigen_ker f rank xpts ypts xwts ywts xeval yeval overx intpts intwts intevalpts =
    let (xfuncwts, yfuncwts, truncsvals, errcode) = eigen_rankn_params f rank xpts ypts xwts ywts in
    let (evalres, evalerr) = eigen_rankn_eval f rank xpts ypts xfuncwts yfuncwts truncsvals xeval yeval in
    let (intevalres, intevalerr) = eigen_rankn_integral_eval f rank xpts ypts xfuncwts yfuncwts truncsvals overx intpts intwts intevalpts in
    let (intcoeffs, coeffserr) = eigen_rankn_integral_coeffs f rank xpts ypts xfuncwts yfuncwts truncsvals overx intpts intwts in
    let _ = (
        Printf.printf "The X function weights matrix:\n";
        print_mat xfuncwts;
        Printf.printf "The Y function weights matrix:\n";
        print_mat yfuncwts;
        Printf.printf "The truncated array of singular values:\n";
        print_vec truncsvals;
        Printf.printf "Error code in singular function parameter calculation:\n";
        Printf.printf "%d\n" errcode;
        Printf.printf "Evaluation:\n";
        Printf.printf "The X eval points:\n";
        print_vec xeval;
        Printf.printf "The Y eval points:\n";
        print_vec yeval;
        Printf.printf "The eval results:\n";
        print_vec evalres;
        Printf.printf "Error code in evaluation:\n";
        Printf.printf "%d\n" evalerr;
        Printf.printf "Integral (really sum) evaluation:\n";
        Printf.printf "Is integral over X?\n";
        if overx then
            Printf.printf "Yes\n"
        else
            Printf.printf "No\n";
        Printf.printf "The points over which the sum is taken:\n";
        print_vec intpts;
        Printf.printf "The weights of the sum:\n";
        print_vec intwts;
        Printf.printf "The points at which the integral is evaluated:\n";
        print_vec intevalpts;
        Printf.printf "The results of the integral evaluation:\n";
        print_vec intevalres;
        Printf.printf "Error code in integral evaluation:\n";
        Printf.printf "%d\n" intevalerr;
        Printf.printf "The coefficients of the rank-n integral:\n";
        print_vec intcoeffs;
        Printf.printf "Error code in integral coefficient calculation:\n";
        Printf.printf "%d\n" coeffserr;
    ) in
    ()

let test_eigen () =
    let myexp x y = exp (x *. y) in
    let myinv x y = 1.0 /. (1.0 +. (x *. y) /. 2.0) in
    let mygauss x y = exp ( (y -. x) *. (x -. y) ) in
    let xpts = cheb_pts 40 in
    let ypts = xpts in
    (* To use singular function approximations effectively, the sum of the weights should
       be the measure of the interval involved for both X and Y. *)
    let xwts = constant_vec (Array1.dim xpts) (2.0 /. (float_of_int (Array1.dim xpts) )) in
    let ywts = xwts in
    let xeval = cheb_pts 10 in
    let yeval = cheb_pts 10 in
    let overx = true in
    let intpts = cheb_pts 20 in
    let intwts = constant_vec (Array1.dim intpts) (1.0 /. (float_of_int (Array1.dim intpts) ) ) in
    let intevalpts = cheb_pts 7 in
    let _ = (
        Printf.printf "Singular function series approximation test.\n";
        Printf.printf "Exponential kernel result:\n";
        test_eigen_ker myexp 3 xpts ypts xwts ywts xeval yeval overx intpts intwts intevalpts;
        Printf.printf "1 / (1 + x * y / 2) kernel result:\n";
        test_eigen_ker myinv 3 xpts ypts xwts ywts xeval yeval overx intpts intwts intevalpts;
        Printf.printf "Gaussian kernel result:\n";
        test_eigen_ker mygauss 3 xpts ypts xwts ywts xeval yeval overx intpts intwts intevalpts;
        ) in
    ()

let test_num_opt_ker f xlb xub ylb yub rank toler maxiter xeval yeval overx intpts intwts intevalpts =
    let (rho, gamma, vmat, wmat, borsuk_lb, errcode) = num_opt_rankn_params f xlb xub ylb yub rank toler maxiter in
    let (evalres, evalerr) = num_opt_rankn_eval f rho gamma vmat wmat xeval yeval in
    let (intevalres, intevalerr) = num_opt_rankn_integral_eval f rho gamma vmat wmat overx intpts intwts intevalpts in
    let (intcoeffs, coeffserr) = num_opt_rankn_integral_coeffs f rho gamma vmat wmat overx intpts intwts in
    let _ = (
        Printf.printf "The rho vector:\n";
        print_vec rho;
        Printf.printf "The gamma vector:\n";
        print_vec gamma;
        Printf.printf "The V matrix:\n";
        print_mat vmat;
        Printf.printf "The W matrix:\n";
        print_mat wmat;
        Printf.printf "The Borsuk lower bound:\n";
        Printf.printf "%f\n" borsuk_lb;
        Printf.printf "Error code in numerically optimal parameter calculation:\n";
        Printf.printf "%d\n" errcode;
        Printf.printf "Evaluation:\n";
        Printf.printf "The X eval points:\n";
        print_vec xeval;
        Printf.printf "The Y eval points:\n";
        print_vec yeval;
        Printf.printf "The eval results:\n";
        print_vec evalres;
        Printf.printf "Error code in evaluation:\n";
        Printf.printf "%d\n" evalerr;
        Printf.printf "Integral (really sum) evaluation:\n";
        Printf.printf "Is integral over X?\n";
        if overx then
            Printf.printf "Yes\n"
        else
            Printf.printf "No\n";
        Printf.printf "The points over which the sum is taken:\n";
        print_vec intpts;
        Printf.printf "The weights of the sum:\n";
        print_vec intwts;
        Printf.printf "The points at which the integral is evaluated:\n";
        print_vec intevalpts;
        Printf.printf "The results of the integral evaluation:\n";
        print_vec intevalres;
        Printf.printf "Error code in integral evaluation:\n";
        Printf.printf "%d\n" intevalerr;
        Printf.printf "The coefficients of the rank-n integral:\n";
        print_vec intcoeffs;
        Printf.printf "Error code in integral coefficient calculation:\n";
        Printf.printf "%d\n" coeffserr;
    ) in
    ()

let test_num_opt () =
    let myexp x y = exp (x *. y) in
    let myinv x y = 1.0 /. (1.0 +. (x *. y) /. 2.0) in
    let mygauss x y = exp ( (y -. x) *. (x -. y) ) in
    let xlb = -1. in
    let xub = 1. in
    let ylb = -1. in
    let yub = 1. in
    let xeval = cheb_pts 10 in
    let yeval = cheb_pts 10 in
    let overx = true in
    let intpts = cheb_pts 20 in
    let intwts = constant_vec (Array1.dim intpts) (1.0 /. (float_of_int (Array1.dim intpts) ) ) in
    let intevalpts = cheb_pts 7 in
    let _ = (
        Printf.printf "Numerically optimal approximation test.\n";
        Printf.printf "Exponential kernel result:\n";
        test_num_opt_ker myexp xlb xub ylb yub 3 1e-15 100 xeval yeval overx intpts intwts intevalpts;
        Printf.printf "1 / (1 + x * y / 2) kernel result:\n";
        test_num_opt_ker myinv xlb xub ylb yub 3 1e-15 100 xeval yeval overx intpts intwts intevalpts;
        Printf.printf "Gaussian kernel result:\n";
        test_num_opt_ker mygauss xlb xub ylb yub 3 1e-15 100 xeval yeval overx intpts intwts intevalpts;
        ) in
    ()

let test_asymm () =
    let blackopt x y = black_formula true x y 0.20 1.0 in
    let xlb = 0.5 in
    let xub = 1. in
    let ylb = 0.9 in
    let yub = 1.2 in
    let xpts = cheb_pts_trans 50 xlb xub in
    let ypts = cheb_pts_trans 50 ylb yub in
    (* To use singular function approximations effectively, the sum of the weights should
       be the measure of the interval involved for both X and Y. *)
    let xwts = constant_vec 50 ( (xub -. xlb) *. (1. /. 50.) ) in
    let ywts = constant_vec 50 ( (yub -. ylb) *. (1. /. 50.) ) in
    let xeval = equal_spaced_grid 10 xlb xub in
    let yeval = equal_spaced_grid 10 ylb yub in
    let overx = true in
    let intpts = equal_spaced_grid 30 xlb xub in
    let intwts = constant_vec (Array1.dim intpts) (1.0 /. (float_of_int (Array1.dim intpts) ) ) in
    let intevalpts = equal_spaced_grid 9 ylb yub in
    let _ = (
        Printf.printf "Tests of an asymmetric kernel (the Black call formula as a function of strike (X) and forward (Y)) on an asymmetric domain:\n";
        Printf.printf "Tests of a Taylor series approximation:\n";
        test_taylor_ker blackopt 3 0.75 1.05 0.01 xeval yeval overx intpts intwts intevalpts;
        Printf.printf "Tests of a singular function approximation:\n";
        test_eigen_ker blackopt 3 xpts ypts xwts ywts xeval yeval overx intpts intwts intevalpts;
        Printf.printf "Tests of a numerically-optimal approximation:\n";
        test_num_opt_ker blackopt xlb xub ylb yub 3 1e-15 100 xeval yeval overx intpts intwts intevalpts;
        ) in
    ()

let () = 
    test_zeros_optima ();
    test_linear_algebra ();
    test_borsuk_lb ();
    test_worst_rho_gamma ();
    test_taylor ();
    test_eigen ();
    test_num_opt ();
    test_asymm ();

