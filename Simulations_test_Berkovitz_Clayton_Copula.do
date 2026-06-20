local simulations = 1000

// Chose sample size
local N = 500 

// Chose the distribution of X
local X_distribution = "exponential"   
//local X_distribution = "log_normal"

local beta_0 = 2.5
local beta = 3

local sigma_eps = 1
local thetas "0.001 0.25 0.5 1 2 5"
local n_thetas = wordcount("`thetas'")
local total = `n_thetas' * `simulations'
set matsize `total'
matrix results = J(`total', 4, .)
local k = 0

foreach theta of local thetas {


forvalues j= 1/`simulations'{
clear
set seed `=`j'+32'
set obs 1000

local i = `j' + `k' * `simulations'
matrix results[`i',1] = `theta'
gen V = rgamma(1/`theta', 1)


gen E1 = -ln(runiform())
gen E2 = -ln(runiform())


gen U_x = (1 + E1/V)^(-1/`theta')
gen U_eps = (1 + E2/V)^(-1/`theta')


sum U_x U_eps


gen x_star  = invnormal(U_x)   
gen eps_star = invnormal(U_eps)  
gen eps = `sigma_eps' * eps_star

matrix results[`i',4] = `theta'

if ("`X_distribution'" == "exponential") {
    generate x = -ln(1 - U_x)
}
else if ("`X_distribution'" == "log_normal") {
    generate x = exp(x_star)
}

generate y = `beta_0' + `beta' * x + eps

// Generate Û_x empirically
sort x
generate U_x_hat = _n / (_N + 1)

// If you want to use a kernel
/*qui summarize x, detail
local sigma = r(sd)
local iqr = r(p75) - r(p25)
local h_factor = min(`sigma', `iqr'/1.34)
local h = 0.9 * `h_factor' * `=_N'^(-1/5)
kdensity x, nograph at(x) generate(f_x_hat) bwidth(`=`h'/2') kernel(gaussian)
generate trapeces = 0.5 * (x - x[_n -1]) * (f_x_hat + f_x_hat[_n - 1])
generate U_x_hat = sum(trapeces)
drop if U_x_hat <=0
drop if U_x_hat >=1
drop trapeces*/


generate x_star_hat = invnorm(U_x_hat)
generate x_star_hat_sq = x_star_hat^2
generate x_star_hat_cb = x_star_hat^3

qui regress y x x_star_hat // x_star_hat_sq x_star_hat_cb
matrix coef = e(b)
local beta_hat =  coef[1,1]
local lambda_hat = coef[1,2]
local beta_0_hat =  coef[1,3]


generate eps_hat = y - `beta_0_hat' - `beta_hat' * x
qui summarize eps_hat
local sigma_eps_hat = sqrt(r(Var) * (r(N)-1) / r(N))
local rho_hat = `lambda_hat' / `sigma_eps_hat'


// Generate Û_eps empirically
sort eps_hat
generate U_eps_hat = _n / (_N + 1)

// If you want to use a kernel.
/*replace eps_hat = atan(eps_hat)
qui summarize eps_hat, detail
local sigma = r(sd)
local iqr = r(p75) - r(p25)
local h_factor = min(`sigma', `iqr'/1.34)
local h = 0.9 * `h_factor' * `=_N'^(-1/5)
kdensity eps_hat, nograph at(eps_hat) generate(f_eps_hat) bwidth(`=`h'/2') kernel(gaussian)
generate trapeces = 0.5 * (eps_hat - eps_hat[_n -1]) * (f_eps_hat + f_eps_hat[_n - 1])
generate U_eps_hat = sum(trapeces)
drop if U_eps_hat <=0
drop if U_eps_hat >=1
drop trapeces*/

generate eps_star_hat = invnorm(U_eps_hat)
qui regress eps_star_hat x_star_hat, noconstant
qui predict v, residual
qui summarize v
local s_v = sqrt(r(Var) * (r(N)-1) / r(N))
generate z = v / `s_v'
drop v

qui regress z x_star_hat x_star_hat_sq x_star_hat_cb, noconstant 
local ll_unrestricted = e(ll)
local n = e(N)

display `ll_unrestricted'

gen z_sq = z^2
qui summarize z_sq
local sum_z_sq = r(sum) 
local ll_restricted = -`n'/2 * log(2*_pi) - 0.5 * `sum_z_sq'

display `ll_restricted'

local LR = -2 * (`ll_restricted' - `ll_unrestricted')
display "LR statistic = " `LR'
display "p-value asintótico (chi2_1) = " 1 - chi2(1, `LR')


matrix results[`i',2] = `LR'
matrix results[`i',3] = 1 - chi2(2, `LR')
matrix results[`i',4] = 1 - chi2(1, `LR')
}
local k = `k'+1
}



clear 

svmat results
rename results1 theta
rename results2 test_Copula_Gaussiana
rename results3 p_score_GL_2
rename results4 p_score_GL_1



generate reject_GL_2 = (p_score_GL_2 < 0.05)
generate reject_GL_1 = (p_score_GL_1 < 0.05)

preserve
collapse (mean) mean_reject=reject_GL_1, by(theta)
export excel using "resultados_Clayton_Copula_`X_distribution'_`N'.xlsx", firstrow(variables) replace
restore
