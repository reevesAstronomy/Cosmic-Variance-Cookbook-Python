import numpy as np

# Cosmic variance calculations
def calc_COSMOS_cosmic_variance(mean_z, delta_z, log_m_stellar, survey="COSMOS"):
    """
    Calculates the cosmic variance in the COSMOS field, based on the 'recipe' in:
    "A COSMIC VARIANCE COOKBOOK" by Moster+2011. See section '3.4. Cookbook for Cosmic Variance'.
    
    Although this function only works for COSMOS, Moster+2011's cookbook also has details
    for how to calculate cosmic variance for UDF, GOODS, GEMS, and EGS.
    
    Inputs
    ------
    mean_z : mean redshift being considered
    delta_z : redshift bin size
    log_m_stellar: log10(mean stellar mass in your bin); make sure to round to the nearest 0.25 or 0.75!
    survey: string indicating which survey field's fitting parameters to use; from Table 3
    
    Outputs
    -------
    delta_gg : root cosmic variance for galaxies (with the above specified inputs)
    
    """
    
    # Step 1: Choose survey field in Table 1 (values below are from Table 3)
    if survey=="UDF":
        sigma_a, sigma_b, beta = 0.251, 0.364, 0.358
    elif survey=="GOODS":
        sigma_a, sigma_b, beta = 0.261, 0.854, 0.684
    elif survey=="GEMS":
        sigma_a, sigma_b, beta = 0.161, 0.520, 0.729
    elif survey=="EGS":
        sigma_a, sigma_b, beta = 0.128, 0.383, 0.673
    elif survey=="COSMOS":
        print("COSMOS")
        sigma_a, sigma_b, beta = 0.069, 0.234, 0.834
    else:
        print("Survey field must be specified (to determine which row of table 3 to use)")

    # Step 2: Choose mean redshift mean_z and redshift bin size delta_z
    # Already specified as inputs to this function

    # Step 3: Choose stellar mass range/value
    # Already specified as inputs to this function; would be nice to add interpolation...
    
    # Step 4a: Calculate dark matter root cosmic variance dm(mean_z, delta_z=0.2) using
    # equation (10) and parameters for your survey from table 3.
    sigma_DM_meanz_deltaz02 = sigma_a / (mean_z**beta + sigma_b) # Equation (10)
    
    # Step 4b: Calculate the galaxy bias b_Mstar_meanz using equation (13) and the
    # parameters for your stellar mass bin or threshold as given in Table 4
    # Table 4: log_mg, b_0, b_1, b_2 (note: uncertainties in b_0, b_1, b_2 neglected here)
    dict_galaxy_bias_fit_parameters = {
    8.75: [0.062, 2.59, 1.025],
    9.25: [0.074, 2.58, 1.039],
    9.75: [0.042, 3.17, 1.147],
    10.25: [0.053, 3.07, 1.225],
    10.75: [0.069, 3.19, 1.269],
    11.25: [0.173, 2.89, 1.438],
    '>8.5': [0.063, 2.62, 1.104],
    '>9.0': [0.085, 2.50, 1.098],
    '>9.5': [0.058, 2.96, 1.192],
    '>10.0': [0.072, 2.90, 1.257],
    '>10.5': [0.093, 3.02, 1.332],
    '>11.0': [0.185, 2.86, 1.448]
    }

    if log_m_stellar not in dict_galaxy_bias_fit_parameters:
        # find closest key
        print("WARNING: input 'log_m_stellar' not found in 'dict_galaxy_bias_fit_parameters'")
        dict_galaxy_bias_fit_parameters_floats = {
        8.75: [0.062, 2.59, 1.025],
        9.25: [0.074, 2.58, 1.039],
        9.75: [0.042, 3.17, 1.147],
        10.25: [0.053, 3.07, 1.225],
        10.75: [0.069, 3.19, 1.269],
        11.25: [0.173, 2.89, 1.438]
        }
        log_m_stellar = min(dict_galaxy_bias_fit_parameters_floats, key=lambda x:abs(x-log_m_stellar+0.001))
        
    gal_bias_params_ = dict_galaxy_bias_fit_parameters[log_m_stellar]
    b_0, b_1, b_2 = gal_bias_params_[0], gal_bias_params_[1], gal_bias_params_[2]
    
    b_Mstellar_meanz = b_0 * (mean_z + 1.)**b_1 + b_2 # Equation 13
    
    # Step 4c: Compute the root cosmic variance for galaxies
    delta_gg = b_Mstellar_meanz * sigma_DM_meanz_deltaz02 * np.sqrt(0.2/delta_z)
    
    return delta_gg
    
def get_cosmic_variance_array(mean_z, delta_z, mass_array, survey="COSMOS"):
    """
    Uses "calc_COSMOS_cosmic_variance" to get the cosmic variance for each log(stellar mass) value
    in the input mass_array.

    Inputs
    ------
    mean_z : average redshift of the redshift cube
    delta_z : redshift cut/thickness
    mass_array : numpy array of log(stellar masses)
    survey: string indicating which survey field's fitting parameters to use; from Table 3

    Outputs
    -------
    cosmic_var_array : a numpy array of root cosmic variances for galaxies
    """
    
    cosmic_var_array = np.zeros(len(mass_array))

    for i, mass in enumerate(mass_array):
        mod = round(round(mass, 1) % 0.5, 2)
        if mod <= 0.0001:
            mass -= 0.25
            mass = round(mass, 2)
        else:
            mod_ = round(round(mass, 1) % 0.25, 2)
            mass_rounded = (mass//0.25)*0.25
            if round(round(mass_rounded, 2) % 0.5, 2) <= 0.0001: mass_rounded+=0.25
            mass = mass_rounded

        if mass > 11.25: 
            mass = '>11.0'

        var_cosmic_ = calc_COSMOS_cosmic_variance(mean_z, delta_z, log_m_stellar=mass, survey=survey)
        cosmic_var_array[i] = var_cosmic_
        
    return cosmic_var_array
