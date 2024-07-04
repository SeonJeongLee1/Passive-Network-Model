# Supporting Data for: "Enhancing simulation feasibility for multi-layer 2D MoS2 RRAM devices: reliability performance learnings from passive network model"

This repository contains the MATLAB code and data supporting the paper titled "Enhancing simulation feasibility for multi-layer 2D MoS2 RRAM devices: reliability performance learnings from passive network model". The code provides input parameters and resistance images for simulating the I-V characteristics of individual devices.

## Authors

**Seonjeong Lee<sup>a†</sup>, Yifu Huang<sup>b†</sup>, Yao-Feng Chang<sup>c</sup>, Seungjae Baik<sup>d</sup>, Jack C. Lee<sup>*b</sup>, and Minsuk Koo<sup>*e</sup>

Affiliations:
- <sup>†</sup> School of Electrical and Computer Engineering, University of Seoul, Seoul 02504, South Korea
- <sup>b</sup> Department of Department of Electrical and Computer Engineering, University of Texas at Austin, 10100 Burnet Road, 78758, Austin TX, USA
- <sup>c</sup> Department of Intel Corporation, 2501 NE Century Road, 97124, Hillsboro OR, USA
- <sup>d</sup> Department of Semiconductor Research and Development Center, Samsung Electronics, Hwaseong-si 18448, South Korea
- <sup>e</sup> Department of Department of Computer Science and Engineering, Incheon National University, Incheon 22012, South Korea 

Contact: koo@inu.ac.kr

<sup>†</sup> Seonjeong Lee and Yifu Huang contributed equally.

## File Structure

### MATLAB Code

- **PNM_parameters.m:** Contains input parameters set for simulating the I-V characteristics of individual devices.

#### Key Parameters

- `V2 = -1.0`: voff (Reset voltage)
- `V = 3.5`: von (Set voltage)
- `m = 10`: Number of nodes 
- `n = 40`: Number of nodes
- `num_splits = m - 1`: Parameter to distribute defects from top to bottom in the Passive Network Model (PNM)
- `percent_defect = 20`: Top defect percentage
- `percent_defect3 = 1`: Bottom defect percentage
- `fail1 = 0`: Probability of failure during the reset process
- `fail2 = 0`: Probability of failure during the set process
- `Rv_h_B`: Rv
- `Rh_h_B`: Rh
- `Rd_h_B = sqrt(2) * Rv_h_B`: Rd

### Functions

- **variation_Voltage.m:** Generates variations in `Vset` and `Vreset` based on standard deviations.

- **Resistance Functions:**
  - **Schottky Emission Model:**
    - `R_h`: High resistance
    - `Rn_h`: Individual resistance with Schottky emission model applied
    - `Rh_h_A`: Average high resistance (horizontal)
    - `Rv_h_A`: Average high resistance (vertical)
    - `Rd_h_A`: Average high resistance (diagonal)
    - `layer_Rh_s_A`, `layer_Rv_s_A`, `layer_Rd_s_A`: Adjusted standard deviations for high resistance
  - **Ohmic Behavior Model:**
    - `R_l`: Low resistance
    - `Rh_l_B`: Average low resistance (horizontal)
    - `Rv_l_B`: Average low resistance (vertical)
    - `Rd_l_B`: Average low resistance (diagonal)
    - `layer_Rh_ss_B`, `layer_Rv_ss_B`, `layer_Rd_ss_B`: Adjusted standard deviations for low resistance

## Usage

- Ensure all dependencies and functions (like `variation_Voltage`) are included in your working directory.
- Adjust parameters according to your specific requirements.
- Run the script `PNM_parameters.m` to set up the input parameters and visualize the results.

## Example

```matlab
clear;
% Set voltage and resistance parameters
V = 3.5;
V2 = -1.0;
step = 0.01;
m = 10;
n = 40;

% Define voltage variation function
[Vset_h, Vreset_h] = variation_Voltage(Vset, Vreset, Vset_change, Vreset_change, (n-1)*(m-2));

% Define resistance functions
Rh_h_A = @(VV) Rn_h(VV);

% Assign resistances and create defect distribution
for k = 1:m-1
    ...
end

% Create and display resistance image
imagesc(A_image_R, clims);
colorbar;
colormap(gray);
