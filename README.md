tdistfit
========

Some tools for fitting t-distributions to data.

# EM methods

These methods find the maximum likelihood parameters using the [expectation-maximization algorithm](http://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm). 
Since I am fitting these distributions primarily to calculate entropy I am using covergence of entropy as a stopping criteria for the EM algorithm (rather than the full likelihood) but it is easy to change this if it is not suitable for your purposes.

* `fitt` : fits a multivariate t-distribution using ECME algorithm [^1]
* `fitt_fixnu` : fits a t-distribution with d.o.f. (nu) specified.
* `fitt_commonnu` : fit t-distributions to grouped data, with d.o.f. (nu) common across groups
* `fitt_commonsnu` : fit t-distributions to grouped, with covariance (S) and d.o.f. (nu) common across groups

# Approximate methods

These use a closed form approximation the ML estimate which is faster to compute. However, they didn't work well for me - with the data I was using I sometimes got negative values for terms which should be non-negative (although it seemed to work OK with generated t-distributed samples). 

* `fitt_approx` : fits using the approximate method of Aeschliman et al. [^2]

[^1]: C Liu and D B Rubin, (1995) "ML estimation of the t distribution using EM and its extensions, ECM and ECME", Statistica Sinica, [5, pp19-39](http://www3.stat.sinica.edu.tw/statistica/oldpdf/A5n12.pdf)

[^2]: C Aeschlimna, J Park and KA Cak, "A Novel Parameter Estimation Algorithm for the Multivariate t-Distribution and its Application to Computer Vision" [ECCV 2010](http://link.springer.com/chapter/10.1007%2F978-3-642-15552-9_43)

# License

This project is licensed under the GNU General Public License. For the exact terms please see the [LICENSE file](https://github.com/robince/tdistfit/blob/master/LICENSE).

vim: set ft=markdown:
