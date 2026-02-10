# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
# ]
# ///
import marimo

__generated_with = "0.19.7"
app = marimo.App(
    width="full",
    app_title="ARSENAL - Formation of stellar populations",
    layout_file="layouts/plots_marimo.grid.json",
)


@app.cell
def _(mo):
    mo.md(r"""
    # Case 1: Formation of stellar populations
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Initialize a basic single stellar population (SSP).
    Choose:
    * its total initial mass
    * its metals
    * the initial mass function (type, mass range and slope)

    TODO: BSP and CSP
    """)
    return


@app.cell
def _(mo):
    # initial mass of the population - solar masses
    text_mass = mo.ui.text(placeholder="Mass in solar masses", label="Initial mass [MSun]", value = '1e6')
    # metallicity of the population - log10(Fe/H)
    text_metallicity = mo.ui.text(placeholder="Metallicity", label="Metallicity [log10(Fe/H)]", value = '0.012')
    # parameters of the IMF - lowest/highest masses and slope
    text_imf_min = mo.ui.text(placeholder="Lowest mass IMF", label="Lowest mass of the IMF [MSun]", value = '0.08')
    text_imf_max = mo.ui.text(placeholder="Highest mass IMF", label="Highest mass of the IMF [MSun]", value = '100')
    text_imf_slope = mo.ui.slider(start = 2.0, stop = 2.5, step = 0.1, value = 2.3, label = "Slope of the IMF")
    return (
        text_imf_max,
        text_imf_min,
        text_imf_slope,
        text_mass,
        text_metallicity,
    )


@app.cell
def _(mo, text_mass):
    mo.hstack([text_mass, mo.md(f"Initial mass of {text_mass.value} MSun")])
    return


@app.cell
def _(mo, text_metallicity):
    mo.hstack([text_metallicity, mo.md(f"log10(Fe/H) = {text_metallicity.value}")])
    return


@app.cell
def _(mo, text_imf_min):
    mo.hstack([text_imf_min, mo.md(f"Lowest mass of the IMF is {text_imf_min.value} MSun")])
    return


@app.cell
def _(mo, text_imf_max):
    mo.hstack([text_imf_max, mo.md(f"Highest mass of the IMF is {text_imf_max.value} MSun")])
    return


@app.cell
def _(mo, text_imf_slope):
    mo.hstack([text_imf_slope, mo.md(f"Slope of the IMF is {text_imf_slope.value}")])
    return


@app.cell
def _(text_mass, u):
    initial_mass = u.Quantity(text_mass.value, unit = "Msun")
    return (initial_mass,)


@app.cell
def _(text_metallicity, u):
    initial_metals = u.Quantity(text_metallicity.value)
    return (initial_metals,)


@app.cell
def _(mo):
    mo.md(r"""
    Create the instance of the IMF and the formation context
    """)
    return


@app.cell
def _(
    ag,
    initial_mass,
    initial_metals,
    text_imf_max,
    text_imf_min,
    text_imf_slope,
    u,
):
    # create the instance of the IMF
    imf_min_mass = u.Quantity(text_imf_min.value, unit = u.Msun)
    imf_max_mass = u.Quantity(text_imf_max.value, unit = u.Msun)
    imf = ag.dist_funcs.imf.Salpeter(imf_min_mass, imf_max_mass,
                                     alpha = u.Quantity(text_imf_slope.value))
    form_context = ag.FormationContext(imf=imf, mass=initial_mass, metals=initial_metals)

    return form_context, imf, imf_max_mass, imf_min_mass


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Check the distribution of masses in comparison to the underlying PDF from which it is sampled.
    """)
    return


@app.cell
def _(form_context, imf, imf_max_mass, imf_min_mass, np, plt):
    mbins = np.logspace(np.log10(imf_min_mass.value), np.log10(imf_max_mass.value), 50)

    plt.hist(form_context.generate_population().value, bins=mbins,density=True,histtype='step',color=plt.cm.Blues(0.5),label="Sampled IMF")
    plt.plot(mbins, imf._pdf(mbins), color=plt.cm.Reds(0.5), label="IMF PDF")

    plt.yscale("log")
    plt.xscale("log")
    plt.xlim(imf_min_mass.value, imf_max_mass.value)
    plt.xlabel(r"$M \, [M_{\odot}]$")
    plt.ylabel(r"$dN/d\log M$")
    plt.legend()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Plot the specifc SN rate for the population at a sample of different times
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Construct the isochrone for the population, first each star in the discrete population, and then the underlying isochrone that is created to interpolate the isochrone to. Output timing for each to give an idea of how long these things take. As you can see, a good fraction of the time goes in to simply constructing the isochrone.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Create a bunch of different stellar population instances with different choices around discrete versus continuous populations and isochrone-based or track-based interpolation. This should take about 8 seconds to initialize.
    """)
    return


@app.cell
def _(mo):
    run_button = mo.ui.run_button()
    run_button
    return (run_button,)


@app.cell
def _(mo, run_button):
    mo.stop(not run_button.value, mo.md("Click 👆 to run this cell"))
    mo.md("You clicked the button! 🎉")
    return


@app.cell
def _():
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Test timing for calculating the luminosity of the population at a single time.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Compute a history of the population's bolometric luminosity for a list of times. This should take about 40 seconds.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Compare the different results:
    """)
    return


@app.cell
def _():
    return


@app.cell
def _(mo):
    mo.md(r"""
    ### Import the relevant modules
    """)
    return


@app.cell
def _():
    import arsenal_gear as ag
    import marimo as mo
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u
    import astropy

    from time import time

    from scipy.integrate import trapezoid as trapz
    return ag, mo, np, plt, u


if __name__ == "__main__":
    app.run()
