# Evidence for a HURP/EB free mixed-nucleotide zone in kinetochore-microtubules
Implementation of a minimal computational model of HURP dynamics on K-fibres, as described in Castrogiovanni, Inchingolo et al (2021) "Evidence for a HURP/EB free mixed-nucleotide zone in kinetochore-microtubules" https://doi.org/10.1101/2021.07.23.453504

The model is a partial differential equation model implementated in MATLAB.

## Developer

Jonathan U. Harrison (jonathan.u.harrison@warwick.ac.uk),
                Zeeman Institute,
		Mathematics Institute,
		University of Warwick

## Citation Information

This code is provided as supplementary information to the paper

Castrogiovanni, Inchingolo, Harrison, Dudka, Sen, Burroughs, McAinsh, Meraldi (2021) "Evidence for a HURP/EB free mixed-nucleotide zone in kinetochore-microtubules" https://doi.org/10.1101/2021.07.23.453504

## Licensing
This source code is licensed under the MIT License.

## Requirements
The following is required to run this software:

1. [MATLAB](https://uk.mathworks.com/)            version (2018b or later)

## Usage
Navigate to the main folder of the repository. No install should be required. 
On a Unix type operating system, inference can be run on the observed data set given in the `.csv` file `in_vivo_hurp_switch_data_9_cells.csv` via:
	`nohup matlab -r "inference_observed_data; exit" > out.log < /dev/null &`
This should save MCMC (Markov Chain Monte Carlo) output as a `.mat` file and produce plots summarising the MCMC output. Simulations from the fitted model should also be produced. 
Expected run time to perform the MCMC sampling: approximately 4-12 hours.
