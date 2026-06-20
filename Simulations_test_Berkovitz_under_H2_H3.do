
local simulations = 1000

// Chose sample size
local N = 500 

// Chose the distribution of X
local X_distribution = "exponential"   
//local X_distribution = "log_normal"


local beta_0 = 2.5
local beta = 3

local sigma_eps = 1
local rho = 0.9

local rhos "0 0.2 0.4 0.6 0.8"
local n_rhos = wordcount("`rhos'")
local total = `n_rhos' * `simulations'
set matsize `total'
matrix results = J(`total', 4, .)
local k = 0

foreach rho of local rhos {

forvalues j= 1/`simulations'{

clear 
set seed `=`j'+32'
local i = `j' + `k' * `simulations'
set obs `N'
matrix results[`i',1] = `rho'
generate U_x = runiform()
generate x_star = invnorm(U_x)
generate eps_star = `rho' * x_star + sqrt(1 - `rho' ^2) * rnormal()


if ("`X_distribution'" == "exponential") {
    generate x = -ln(1 - U_x)
}
else if ("`X_distribution'" == "log_normal") {
    generate x = exp(x_star)
}


generate y = `beta_0' + `beta' * x + eps

// Generate Û_x empirically
sort x
generate U_x_hat = _n / (_N+1)

// If you want to use a kernel
/*
qui summarize x, detail
local sigma = r(sd)
local iqr = r(p75) - r(p25)
local h_factor = min(`sigma', `iqr'/1.34)
local h = 0.9 * `h_factor' * `=_N'^(-1/5)
kdensity xx, nograph at(x) generate(f_x_hat) bwidth(`=`h'/2') // kernel(gaussian)
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
local sigma_eps_hat = sqrt(r(Var)) * ((r(N)-1) / r(N))
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
kdensity eps_hat, nograph at(eps_hat) generate(f_eps_hat) bwidth(`=`h'/4') kernel(gaussian)
generate trapeces = 0.5 * (eps_hat - eps_hat[_n -1]) * (f_eps_hat + f_eps_hat[_n - 1])
generate U_eps_hat = sum(trapeces)
drop if U_eps_hat <=0
drop if U_eps_hat >=1
drop trapeces*/

generate eps_star_hat = invnorm(U_eps_hat)

qui regress eps_star_hat x_star_hat
qui predict v, residual
qui summarize v
local s_v = sqrt(r(Var) * (r(N)-1) / r(N))
generate z = v / `s_v'
drop v


qui regress z x_star_hat x_star_hat_sq x_star_hat_cb, noconstant 

local ll_unrestricted = e(ll)
local n = e(N)

gen z_sq = z^2
qui summarize z_sq
local sum_z_sq = r(sum) 
local ll_restricted = -`n'/2 * log(2*_pi) - 0.5 * `sum_z_sq'

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
rename results1 rho
rename results2 test_Copula_Gaussiana
rename results3 p_score_GL_2
rename results4 p_score_GL_1



generate reject_GL_2 = (p_score_GL_2 < 0.05)
generate reject_GL_1 = (p_score_GL_1 < 0.05)

preserve
collapse (mean) mean_reject=reject_GL_1, by(rho)
export excel using "resultados_Gaussian_Copula_`X_distribution'_`N'.xlsx", firstrow(variables) replace
restore
