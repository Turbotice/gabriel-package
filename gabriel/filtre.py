import numpy as np

def filtre_subpixellaire(indice, k, TF):
    Delta = (k[1]-k[0]) * ( np.abs(TF[indice+1])-np.abs(TF[indice-1])) / ( 2* (2*np.abs(TF[indice])-np.abs(TF[indice+1]) -np.abs(TF[indice-1]) ))
    pike = k[indice]+Delta
    return pike

def bl_filt(t, y, fexc, half_width):
    """
    Simple Blackman filter.
    
    """
    dt = t[1]-t[0]
    hwpts = int(round(half_width * np.abs(1/fexc /dt)))
    nf = hwpts * 2 + 1
    x = np.linspace(-1, 1, nf, endpoint=True)
    x = x[1:-1]   # chop off the useless endpoints with zero weight
    w = 0.42 + 0.5 * np.cos(x * np.pi) + 0.08 * np.cos(x * 2 * np.pi)
    # The so-called ``Blackman Window'' is the specific case for which $ \alpha_0 = 0.42$ $ \alpha_1 = 0.5$ , and $ \alpha_2 = 0.08$
    ytop = np.convolve(y, w, mode='same')
    ybot = np.convolve(np.ones_like(y), w, mode='same')
    
    return ytop / ybot

def filtrage_gld_v2022(t,x, fexc, hwidth = 2):
    # # Exponentielle complexe pour la demodulation:
    c = np.exp(-1j * 2 * np.pi * t * fexc)
    product = x * c
    rotary = x.dtype.kind == 'c'  # complex input
     
    # filter half-width number of points

    demod = bl_filt(t, product, fexc, hwidth)

    if not rotary:   
    #     # The factor of 2 below comes from fact that the
    #     # mean value of a squared unit sinusoid is 0.5.
        demod *= 2
        reconstructed = (demod * np.conj(c))
    if not rotary:
        reconstructed = reconstructed.real
    if np.sign(1/fexc) < 0:
        demod = np.conj(demod)
    #   This is to make the phase increase in time
    #   for both positive and negative demod frequency
    #   when the frequency of the signal exceeds the
    #   frequency of the demodulation.    
    #   reconstructed = (demod * np.conj(c))
    #     # reconstructed = reconstructed.real
    return reconstructed