Stochastic Oscillator Workshop

First developed for a workshop at Colorado School of Mines as a satellite to the SIAM 
Applied Dynamical Systems meeting in Denver, Colorado, in 2025.

This folder contains four examples: a noisy Stuart-Landau oscillator, a noisy FitzHugh-Nagumo oscillator, a noisy Morris-Lecar oscillator, 
and a noisy heteroclinic oscillator.

Each script uses the general purpose functions TimeSeries.m and Qfunction.m functions to:
  1. generate time-series data in original coordinates
  2. approximate the Q-function using finite-differences
  3. generate time-series data in Q-function coordinates
  4. generate power spectra in original and Q-function coordinates
  5. visualize the results

The scripts may be modified to handle any 2D stochastic oscillator

To run any of the examples, download the TimeSeries.m, Qfunction.m, and PowerSpectrum.m functions from the 'Dependencies' folder, 
and then run any of the examples Main_XXX.m.

Note that the Qfunction.m function may take several minutes to run if very fine discretizations are used.

---

General purpose functions:

- TimeSeries.m
- Qfunction.m
- PowerSpectrum.m

TimeSeries.m generates time-series data for a noisy stochastic oscillator using the Euler-Maruyama method
Qfunction.m discretizes the SKO and generates the Q-function for a 2D stochastic oscillator
PowerSpectrum.m generates power spectra for a noisy stochastic oscillator in original and Q-function coordinates

Each script contains detailed documentation and example input for the interested user.

These functions are used to generate time-series data in Q-function coordinates for a given 2D stochastic oscillator (see the 'Examples' folder).


