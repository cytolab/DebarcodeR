asinh_trans = function(cofactor) {
  transform =   function(x,...) {
    asinh(x/cofactor)/log(10) + log(cofactor/2)/log(10)
  }
  inverse = function(x,...) {
    sinh(x*log(10) - log(cofactor/2))*cofactor
  }
  breaks = function(x) {
    minor.ticks <- (rep(2:9, 6))*(rep(10^(0:5), each = 8))
    major.ticks <- 10^(2:6)
    breaks <-(c(0, major.ticks,- major.ticks))
    breaks[order(breaks)]
  }
  scales::trans_new("asinh", transform, inverse, breaks)
}

