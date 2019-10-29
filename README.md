# Cosmic variance cookbook: Python version

The scripts here calculate the root galaxy cosmic variance for a given survey, based on the paper "A Cosmic Variance Cookbook": https://arxiv.org/abs/1001.1737, by Moster+2011. I find it useful for quickly and easily putting error bars on galaxy counts within a given volume.

### Prerequisites

Should work with either Python 2 or 3. Also requires numpy to work.

### Example usage

Example usage calculating an array of root cosmic variances given an array of stellar of mass. The "survey" flag is for a given survey in table 3 in the Moster+2011 paper. Note: my code does not interpolate the stellar mass fits, so the code currently just rounds to the nearest log_m_stellar and uses its fit values.

```
mean_z = 1.25
delta_z = 0.5
log_m_stellar = np.array([10., 10.2, 10.4])
get_cosmic_variance_array(mean_z=mean_z, delta_z=delta_z, log_m_stellar=log_m_stellar, survey="COSMOS")
```

