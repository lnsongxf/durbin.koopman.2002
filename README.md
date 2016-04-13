# durbin.koopman.2002
Forward filtering backward smoothing algorithm of Durbin and Koopman (2002)

This respository contains functions for the Kalman filter algorithm proposed in Durbin and Koopman (2002).

It requires the <code>mvtnorm</code> package to draw from multivariate Normal Distributions.

<code>dk.alg2</code> generates a disturbances for the original time series and calls the Kalman filter and smoother.

<code>ffbs.dk</code> runs the Kalman filter and Kalman smoother.
