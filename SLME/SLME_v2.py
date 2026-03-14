import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from scipy.optimize import minimize

os.chdir(os.path.split(os.path.realpath(__file__))[0])
plt.rcParams['font.sans-serif'] = ['Arial']


def calculate_SLME(material_eV_for_absorbance_data,
                   material_absorbance_data,
                   material_direct_allowed_gap,
                   material_indirect_gap,
                   thickness=50e-6,
                   T=293.15,
                   absorbance_in_inverse_centimeters=True,
                   cut_off_absorbance_below_direct_allowed_gap=True):
    # Absorbance unit: m^-1
    if absorbance_in_inverse_centimeters:
        material_absorbance_data = material_absorbance_data * 100.0

    try:
        solar_spectra_wavelength, solar_spectra_irradiance = np.loadtxt(
            "am1.5G.dat", usecols=[0, 1], unpack=True, skiprows=2
        )
    except OSError:
        print("Could not locate am1.5G.dat in the current directory.")
        sys.exit()

    solar_spectra_wavelength_meters = solar_spectra_wavelength * 1e-9

    c = 299792458.0
    h = 6.62607004081e-34
    h_eV = 4.135667516e-15
    k = 1.3806485279e-23
    k_eV = 8.617330350e-5
    e = 1.602176620898e-19

    delta = material_direct_allowed_gap - material_indirect_gap
    fr = np.exp(-delta / (k_eV * T))

    solar_spectra_photon_flux = solar_spectra_irradiance * (solar_spectra_wavelength_meters / (h * c))
    P_in = simpson(solar_spectra_irradiance, x=solar_spectra_wavelength)

    blackbody_irradiance = (
        2.0 * np.pi * h * c**2 / (solar_spectra_wavelength_meters**5)
    ) * (
        1.0 / (np.exp(h * c / (solar_spectra_wavelength_meters * k * T)) - 1.0)
    )

    blackbody_photon_flux = blackbody_irradiance * (solar_spectra_wavelength_meters / (h * c))

    material_wavelength_for_absorbance_data = (
        (c * h_eV) / (material_eV_for_absorbance_data + 1e-8)
    ) * 1e9

    material_absorbance_data_function = interp1d(
        material_wavelength_for_absorbance_data,
        material_absorbance_data,
        kind='cubic',
        fill_value=(material_absorbance_data[0], material_absorbance_data[-1]),
        bounds_error=False
    )

    material_interpolated_absorbance = np.zeros(len(solar_spectra_wavelength_meters))
    direct_gap_wavelength_nm = ((c * h_eV) / material_direct_allowed_gap) * 1e9

    for i in range(len(solar_spectra_wavelength_meters)):
        if solar_spectra_wavelength[i] < direct_gap_wavelength_nm or not cut_off_absorbance_below_direct_allowed_gap:
            material_interpolated_absorbance[i] = material_absorbance_data_function(solar_spectra_wavelength[i])

    absorbed_by_wavelength = 1.0 - np.exp(-2.0 * material_interpolated_absorbance * thickness)

    J_0_r = e * np.pi * simpson(
        blackbody_photon_flux * absorbed_by_wavelength,
        x=solar_spectra_wavelength_meters
    )
    J_0 = J_0_r / fr

    J_sc = e * simpson(
        solar_spectra_photon_flux * absorbed_by_wavelength,
        x=solar_spectra_wavelength
    )

    def J(V):
        return J_sc - J_0 * (np.exp(e * V / (k * T)) - 1.0)

    def power(V):
        return J(V) * V

    def neg_power(V):
        return -power(V)

    results = minimize(neg_power, x0=[0.0001])
    V_Pmax = results.x.item()
    P_m = power(V_Pmax)

    efficiency = P_m / P_in
    return efficiency


if __name__ == '__main__':
    file_name = "ABSORPTION.dat"
    material_direct_allowed_gap = 1.5
    material_indirect_gap = 1.5

    data = np.loadtxt(file_name, comments="#")
    material_eV_for_absorbance_data = data[:, 0]
    material_absorbance_data = (data[:, 1] + data[:, 2] + data[:, 3]) / 3.0

    SLME = calculate_SLME(
        material_eV_for_absorbance_data=material_eV_for_absorbance_data,
        material_absorbance_data=material_absorbance_data,
        material_direct_allowed_gap=material_direct_allowed_gap,
        material_indirect_gap=material_indirect_gap,
        thickness=50e-6,
        T=293.15
    )

    print("File:", file_name)
    print(f"Standard SLME: {SLME * 100:.2f}%")

    thickness_array = np.logspace(-8, -4, 100)
    SLME_list = []

    for thickness in thickness_array:
        slme_val = calculate_SLME(
            material_eV_for_absorbance_data=material_eV_for_absorbance_data,
            material_absorbance_data=material_absorbance_data,
            material_direct_allowed_gap=material_direct_allowed_gap,
            material_indirect_gap=material_indirect_gap,
            thickness=thickness,
            T=293.15
        )
        SLME_list.append(slme_val * 100.0)

    SLME_array = np.array(SLME_list)
    x_plot = thickness_array * 1e6

    plt.figure()
    plt.semilogx(
        x_plot,
        SLME_array,
        color='tab:blue',
        linewidth=2
    )
    plt.fill_between(
        x_plot,
        SLME_array,
        0,
        color='tab:blue',
        alpha=0.25
    )

    plt.xlabel('Thickness ($\\mu$m)', fontsize=20)
    plt.ylabel('SLME (%)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.minorticks_on()

    ax = plt.gca()
    ax.set_xlim(0.01, ax.get_xlim()[1])
    ax.set_ylim(0, ax.get_ylim()[1])

    plt.tight_layout()
    plt.savefig('SLME_thickness.jpg', dpi=600)
    plt.show()