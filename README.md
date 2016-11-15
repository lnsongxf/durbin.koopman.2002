# Forward Filtering Backward Smoothing Algorithm from Durbin and Koopman (2002)

This respository contains functions for the Kalman filter algorithm proposed in Durbin and Koopman (2002).

It requires the <code>mvtnorm</code> package to draw from a multivariate Normal distribution.

<code>dk.alg2</code> generates disturbances for the original time series and calls the Kalman filter and smoother.

<code>ffbs.dk</code> runs the Kalman filter and Kalman smoother.

## Literature
Durbin, J., Koopman, S., 2002. A simple and efficient simulation smoother for state
space time series analysis. Biometrika 89, 603–616.

## Copyright
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
