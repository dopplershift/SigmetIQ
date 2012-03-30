import numpy as N
import pylab as P
from matplotlib import cm
from matplotlib.colors import Normalize
from netCDF4 import Dataset
from numpy.fft import fft, fftfreq, fftshift, ifft
from scipy.signal import correlate

from inspection import DataCursor

def fliplr(x):
  return x[::-1]

def periodogram(time_series):
  Z = fliplr(fftshift(fft(time_series)))
  return N.abs(Z)*N.abs(Z)/time_series.size

def bt_spectrum(time_series):
  R = xcorr(time_series, scale='biased')
  return N.abs(fliplr(fftshift(fft(R))))

def xcorr(x, y=None, maxlag=None, scale=None):
    if y is None:
        y = x.conj()

    M = max(x.shape[-1], y.shape[-1])
    if maxlag is None:
        maxlag = M - 1

    R = N.empty(x.shape[:-1]+(2*maxlag+1,), dtype=N.complex64)
    for l in xrange(-maxlag, maxlag+1):
        if l == 0:
            R[...,l+maxlag] = (x*y).sum(axis=-1)
        elif l < 0:
            R[...,l+maxlag] = (x[...,-l:]*y[...,:l]).sum(axis=-1)
        else:
            R[...,l+maxlag] = (x[...,:-l]*y[...,l:]).sum(axis=-1)

#    R = correlate(x, y)
#    M = int(R.size/2.) + 1.0

#    if not maxlag is None:
#        R = R[M-maxlag-1:M+maxlag]

    if scale == 'biased':
        R = R / M
    elif scale == 'unbiased':
        ns = N.arange(M-maxlag, M)
        R = R / N.r_[ns, M, fliplr(ns)]
    elif not scale is None:
        raise ValueError("Scale must be 'biased', 'unbiased', or None")

    return R

def spect_moments(spectrum, Va):
    M = spectrum.shape[-1]

    #Estimate power by taking the sum of the spectrum
    P = sum(spectrum)

    #Find the peak of the spectrum (max Fourier coefficient) -- center calculations
    #about this peak to remove bias due to aliasing
    km = spectrum.argmax(axis=-1)

    #Estimate the velocity by the first moment of the spectrum and the spectrum
    #width from the central second moment
    k = N.arange(-M/2, M/2) + km
    ind = k%M
    Vr = (2 * Va / M) * ((k * spectrum[ind]).sum(axis=-1)/P - M/2)
    sig = (2 * Va / M) * N.sqrt(((k - km)**2 * spectrum[ind]).sum(axis=-1)/P)

    return P, Vr, sig

def covar_ts(H_series, V_series, scale='unbiased', maxlag=2):
    Rhh = xcorr(H_series, maxlag=maxlag, scale=scale)
    Rvv = xcorr(V_series, maxlag=maxlag, scale=scale)
    Rhv1 = xcorr(H_series, V_series.conj(), maxlag=maxlag, scale=scale)
    Rhv2 = xcorr(V_series.conj(), H_series, maxlag=maxlag, scale=scale)

    return Rhh, Rvv, Rhv1, Rhv2

def auto_doppler(time_series, Va):
    maxlag = 2
    R = xcorr(time_series, maxlag=maxlag, scale='biased')
    return doppler_covar(R, Va, maxlag)

def doppler_covar(R, Va, maxlag=2):
    #Power is autocorrelation at lag 0
    P = N.abs(R[..., maxlag])

    #Velocity is angle (phase) of lag 1
    Vr = N.angle(R[..., maxlag + 1]) * -Va / N.pi

    #Spectrum width comes from lags 1 and 2
    log_rat = N.log(abs(R[..., maxlag + 1]/R[..., maxlag + 2]))
    sig = 2 * Va / (N.sqrt(6) * N.pi) * N.sqrt(abs(log_rat))

    return P, Vr, sig

def auto_dual_pol(H_series, V_series):
    maxlag = 1
    Rhh = xcorr(H_series, maxlag=maxlag, scale='unbiased')
    Rvv = xcorr(V_series, maxlag=maxlag, scale='unbiased')
    Rhv1 = xcorr(H_series, V_series.conj(), maxlag=maxlag, scale='unbiased')
    Rhv2 = xcorr(V_series.conj(), H_series, maxlag=maxlag, scale='unbiased')

    return dual_pol_covar(Rhh, Rvv, Rhv1, Rhv2, maxlag)

def dual_pol_covar(Rhh, Rvv, Rhv1, Rhv2, maxlag=2):
    Zdr = N.abs(Rhh[..., maxlag + 1])/N.abs(Rvv[..., maxlag + 1])

    rho_hv = (N.abs(Rhv1[..., maxlag + 1]) + N.abs(Rhv2[..., maxlag + 1])) / (2
        * N.sqrt(N.abs(Rhh[..., maxlag + 1] * Rvv[..., maxlag + 1])))

    phi_dp = N.angle(Rhv1[..., maxlag])

    return Zdr, rho_hv, phi_dp

def auto_all(Rhh, Rvv, Rhv1, Rhv2, Va, maxlag=2):
    PH, VrH, sigH = doppler_covar(Rhh, Va, maxlag)
    PV, VrV, sigV = doppler_covar(Rvv, Va, maxlag)

    Zdr, rho_hv, phi_dp = dual_pol_covar(Rhh, Rvv, Rhv1, Rhv2, maxlag)

    return PH, VrH, sigH, PV, VrV, sigV, Zdr, rho_hv, phi_dp

def four_panel_ppi(az, rng, **moments):
    if len(moments) != 4:
        raise ValueError('4 moments must be passed as keywords')

    deg_to_rad = N.pi / 180.

    #Convert the az and range values to x,y grid
    x = N.sin(az[:,None] * deg_to_rad) * rng[None,]
    y = N.cos(az[:,None] * deg_to_rad) * rng[None,]

    #Set up colormap to draw masked values as white
    cmap = cm.gist_ncar
    cmap.set_bad('white', 1.0)

    panel = 1
    fig = P.figure()
    axes = list()
    d_list = list()
    for name, (data, limits, units) in moments.iteritems():
        d_list.append(data)
        if not axes:
            axes.append(P.subplot(2, 2, panel))
            ax1 = axes[0]
        else:
            axes.append(P.subplot(2, 2, panel, sharex=ax1, sharey=ax1))

        if limits:
            norm = Normalize(*limits)
        else:
            norm = None

        mesh = axes[-1].pcolormesh(x, y, data, shading='flat', cmap=cmap,
            norm=norm)
        cbar = fig.colorbar(mesh, ax=axes[-1])
        cbar.set_label(units)
        P.axis('equal')
        P.xlim(-300, 300)
        P.ylim(-300, 300)
        P.title(name.replace('_', ' '))
        panel += 1

    rngm = (rng[0:-1] + rng[1:]) / 2.
    azm = (az[0:-1] + az[1:]) / 2.
    diffs = az[0:-1] - az[1:]
    azm[diffs>180.] += 180.
    xm = N.sin(azm[:,None] * deg_to_rad) * rngm[None,]
    ym = N.cos(azm[:,None] * deg_to_rad) * rngm[None,]

    try:
        xm = N.ma.MaskedArray(xm, mask=data.mask)
        ym = N.ma.MaskedArray(ym, mask=data.mask)
    except AttributeError:
        pass

    cursor = DataCursor(xm, ym, d_list, axes)

    return fig

def comp_plot(rng, phi, z, az):
    fig = P.figure()
    line1 = P.plot(rng, phi * 180.0/N.pi, 'b')
    P.ylabel('Differential Phase (deg)')
    P.ylim(-180.0, 180.0)
    ax2 = P.twinx()
    line2 = ax2.plot(rng, 10.0*N.log10(z), 'g')
    P.xlabel('Range (km)')
    P.ylabel('Reflectivity Factor (dBZ)')
    P.xlim(0, 300)
    P.ylim(0, 60)
    P.title('Data Along Radial (%.1f deg)' % az)
    P.legend((line1, line2), ('Differential Phase', 'Reflectivity Factor'),
        loc='upper center')
    return fig


if __name__ == '__main__':
    import sys

    MA = N.ma.MaskedArray
    m_to_km = 1 / 1000.
    ppi_start = 685
    pulse_per_rad = 17
    rad_per_deg = N.pi / 180.0
    oun_cal = pow(10, -2.1729)

    #save = 'savefig' in sys.argv
    save = False

    if len(sys.argv) > 1:
        filename = sys.argv[1]
        print 'Filename:', filename
    else:
        filename = 'KOUNTS.20070629.021635.941.ELEVEN.1.H+V.300.nc'
    nc = Dataset(filename, 'r')

    az = nc.variables['Azimuth'][ppi_start:]
    el = nc.variables['Elevation'][ppi_start:]
    ranges = nc.variables['Range'][:]

    Ts = nc.variables['PRT'][ppi_start:]
    pw = nc.variables['PulseWidth'].getValue()
    lam = nc.variables['Wavelength'].getValue()
    noise_pwr = pow(10., nc.variables['NoisePowerH'].getValue() / 10.)
    dattim = nc.variables['Time'][0]

    ih = nc.variables['I_Horizontal']
    iq_h = N.empty((ih.shape[0] - ppi_start, ih.shape[1]), dtype=N.complex64)
    iq_h = ih[ppi_start:] - 1.0j * nc.variables['Q_Horizontal'][ppi_start:]
    num_rad = iq_h.shape[0] / pulse_per_rad
    iq_h.resize((num_rad, pulse_per_rad, iq_h.shape[-1]))
    iq_h = iq_h.swapaxes(1,2)

    iv = nc.variables['I_Vertical']
    iq_v = N.empty((iv.shape[0] - ppi_start, iv.shape[1]), dtype=N.complex64)
    iq_v = iv[ppi_start:] - 1.0j * nc.variables['Q_Vertical'][ppi_start:]
    iq_v.resize((num_rad, pulse_per_rad, iq_v.shape[-1]))
    iq_v = iq_v.swapaxes(1,2)

    nc.close()

    gate_len = ranges[1] - ranges[0]
    rngb = N.r_[0., (ranges[0:-1] + ranges[1:]) / 2., ranges[-1] + gate_len / 2]

    az.resize((num_rad, pulse_per_rad))
    azb = (az[0:-1, pulse_per_rad - 1] + az[1:, 0]) / 2.
    diffs = az[0:-1, pulse_per_rad - 1] - az[1:, 0]
    azb[diffs>180.] += 180.
    azb = N.r_[az[0,0], azb, az[num_rad -1 , pulse_per_rad - 1]]

    el.resize((num_rad, pulse_per_rad))
    elb = (el[0:-1, pulse_per_rad - 1] + el[1:, 0]) / 2.
    elb = N.r_[el[0,0], elb, el[num_rad -1 , pulse_per_rad - 1]]

    Va = lam/(4*Ts[0])

    Rhh, Rvv, Rhv1, Rhv2 = covar_ts(iq_h, iq_v)

    (pwrH, velH, spwH, pwrV, velV, spwV, zdr, rho_hv, phi_dp) = auto_all(Rhh,
        Rvv, Rhv1, Rhv2, Va)

    mask = (pwrH < noise_pwr) | (pwrV < noise_pwr)

    refH = pwrH * ranges**2 * oun_cal
    refV = pwrV * ranges**2 * oun_cal

    fig = four_panel_ppi(azb, rngb * m_to_km,
        Power=(MA(10.0*N.log10(pwrH), mask=mask), None, 'dBm'),
        Reflectivity=(MA(10.0*N.log10(refH), mask=mask), None, 'dBZ'),
        Velocity=(MA(velH, mask=mask), None, 'm/s'),
        Spectrum_Width=(MA(spwH, mask=mask), None, 'm/s'))
    if save:
        fig.savefig('horiz_moments.png', dpi=300)
    fig.axes[0].axis([0,275,0,275])
    fig.canvas.draw()
    if save:
        fig.savefig('horiz_moments_zoom.png', dpi=300)

    fig = four_panel_ppi(azb, rngb * m_to_km,
        Power=(MA(10.0*N.log10(pwrV), mask=mask), None, 'dBm'),
        Reflectivity=(MA(10.0*N.log10(refV), mask=mask), None, 'dBZ'),
        Velocity=(MA(velV, mask=mask), None, 'm/s'),
        Spectrum_Width=(MA(spwV, mask=mask), None, 'm/s'))
    if save:
        fig.savefig('vert_moments.png', dpi=300)
    fig.axes[0].axis([0,275,0,275])
    fig.canvas.draw()
    if save:
        fig.savefig('vert_moments_zoom.png', dpi=300)

    fig = four_panel_ppi(azb, rngb * m_to_km,
        ReflectivityH=(MA(10.0*N.log10(refH), mask=mask), None, 'dBZ'),
        Differential_Reflectivity=(MA(10.0*N.log10(zdr),mask=mask),(-5,10),'dB'),
        Correlation_Coefficient=(MA(rho_hv, mask=mask), (0.3, 1.1), ''),
        Differential_Phase=(MA(phi_dp, mask=mask), (0, N.pi), 'radians'))
    if save:
        fig.savefig('polar_moments.png', dpi=300)
    fig.axes[0].axis([0,275,0,275])
    fig.canvas.draw()
    if save:
        fig.savefig('polar_moments_zoom.png', dpi=300)

    rad = 316
    xl = m_to_km * ranges[[0,-1]]*N.sin(az[rad, 9] * rad_per_deg)
    yl = m_to_km * ranges[[0,-1]]*N.cos(az[rad, 9] * rad_per_deg)
    for i in xrange(0, len(fig.axes), 2):
        fig.axes[i].plot(xl, yl, 'k', linewidth=2, zorder=3)
    fig.axes[0].axis([0,275,0,275])
    fig.canvas.draw()

    if save:
        fig.savefig('polar_moments_line.png', dpi=300)

    fig = comp_plot(ranges * m_to_km, phi_dp[rad], refH[rad], az[rad, 9])
    if save:
        fig.savefig('radial_plot.png', dpi=300)

    phi_dp_pulse = N.angle(iq_h * iq_v.conj())
    phi_dp_min = phi_dp_pulse[rad].min(axis=-1) / rad_per_deg
    phi_dp_max = phi_dp_pulse[rad].max(axis=-1) / rad_per_deg
    phi_dp_avg = phi_dp_pulse[rad].mean(axis=-1) / rad_per_deg
    phi_dp_rng = (phi_dp_avg - phi_dp_min, phi_dp_max - phi_dp_avg)

    fig = P.figure()
    P.errorbar(ranges * m_to_km, phi_dp[rad] / rad_per_deg,
        yerr=phi_dp_rng, fmt='b-', ecolor='g')
    P.xlabel('Range (km)')
    P.ylabel('Differential Phase (deg)')
    P.title('Spread of Differential Phase')
    P.xlim(0,300)
    P.ylim(-360,360)
    if save:
        fig.savefig('radial_plot_range.png', dpi=300)

    if not save:
        P.show()
