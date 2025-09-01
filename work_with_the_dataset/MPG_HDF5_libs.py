#
#
# Hans-Rainer Kloeckner
#
# MPIfR 2022
#
#
# This is part of the example scripts to get you a head start using a HDF5 files
# of the MPG telescope 
#

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

import h5py
import numpy as np
import numpy.ma as ma
from functools import reduce



def h5printstructure(item, leading = ''):
    """
    print the structure of the hdf5 file

    https://stackoverflow.com/questions/34330283/how-to-differentiate-between-hdf5-datasets-and-groups-with-h5py
    """
    for key in item:
        if isinstance(item[key], h5py.Dataset):
            print(leading + key + ': ' + str(item[key].shape))            
        else:
            print(leading + key)
            h5printstructure(item[key], leading + ' -')


def h5keypath(h5file, leading = '',kplist = []):
    """
    compiles a complete key path list of the hdf5 file
    """
    for key in h5file:
        if isinstance(h5file[key], h5py.Dataset):
            kplist.append(leading + key)
        else:
            h5keypath(h5file[key], leading + key +'/',kplist)

    return(kplist)



def findkeys(h5file,keys=[],exactmatch=False):
    """
    search hdf5 file for combination of keys to extract the
    key path within the data structure
    exactmatch == True 
    will provide the keypath of only an exact keywords matches
    e.g. exactmatch == True and keys = [timestamp] will NOT include the 
                            path of: monitor/acu/actual_timestamp but
                            exactmatch == False will
    """

    keyplist       = np.array(h5keypath(h5file, leading = '',kplist = []),dtype=str)

    select_on_keys = np.ones(len(keyplist)).astype(dtype=bool)

    for k in keys:
        if exactmatch:
                sellogicforkeys = np.zeros(len(keyplist))
                for h in range(len(keyplist)):
                        klsp = keyplist[h].replace('/',' ').split()
                        for m in klsp:
                                if m == k:
                                        sellogicforkeys[h] = 1
     
        else:
                lookforkeys     = np.char.find(keyplist,k)
                sellogicforkeys = lookforkeys > 0


        select_on_keys  = np.logical_and(select_on_keys,sellogicforkeys.astype(dtype=bool))

    return list(keyplist[select_on_keys])

def index_of(val, in_list):
    try:
        return in_list.index(val)
    except ValueError:
        return None 

def return_equals(*arrays):
    """
    Compares multiple matches and provides
    the indexes of the matches for each input array
    #
    https://stackoverflow.com/questions/30082052/most-efficient-way-to-implement-numpy-in1d-for-muliple-arrays
    """
    matched = reduce(np.intersect1d, arrays)
    indx_per_array = np.array([np.where(np.in1d(array, matched))[0] for array in arrays])

    return indx_per_array
 

def return_equal_data(hdf5data,keys):
    """
    return index of the data array 
    """
    org_indx            = []
    check_data_to_match = []
    for k in keys:
        data = np.array(hdf5data[k])
        org_indx.append(np.zeros(data.shape[0]))
        check_data_to_match.append(data)

    indx_keys = return_equals(*check_data_to_match)

    equal_array_masks = {}
    for k,t in enumerate(keys):
 
        org_indx[k][indx_keys[k]] = 1

        equal_array_masks[t] = {}
        equal_array_masks[t]['selmask'] = indx_keys[k]
        equal_array_masks[t]['mask']    = org_indx[k]

    return equal_array_masks


def get_good_data(data,sel_data,sigma,stats_type='mean',do_info=False):
    """
    """
 
    if stats_type == 'mean':
        stats         = [np.mean(data[sel_data.astype(bool)]),np.std(data[sel_data.astype(bool)])]

    if stats_type == 'kdemean':
        stats         = kdemean(data[sel_data.astype(bool)],accucary=1000)

    # do a selection on pure statistics
    sel  = np.logical_and(data > stats[0] - sigma * stats[1], data < stats[0] + sigma*stats[1])

    if do_info:
        print('data ',len(data),' markes as bad ',len(data)-np.count_nonzero(sel))

    mask      = np.zeros(len(sel_data))
    mask[sel] = 1

    return sel_data * mask




def hrk_fftconvolve2d(x, y, mode="full"):
    """
    just rebuild the scipy convolve2d function to understand its workings

    x and y must be real 2-d numpy arrays.

    mode must be "full" or "valid".

    inspired by:
    https://stackoverflow.com/questions/22325211/valid-and-full-convolution-using-fft2-in-python
    https://stackoverflow.com/questions/50453981/implement-2d-convolution-using-fft
    """
    from numpy.fft import fft2, ifft2

    # this is needed to cope with the dimension of the kernel
    x = np.pad(x,((y.shape[0]-1,0),(y.shape[1]-1,0)),mode='symmetric')

    x_shape = np.array(x.shape)
    y_shape = np.array(y.shape)
    z_shape = x_shape + y_shape - 1
    z = ifft2(fft2(x, z_shape) * fft2(y, z_shape)).real

    if mode == "valid":
        # To compute a valid shape, either np.all(x_shape >= y_shape) or
        # np.all(y_shape >= x_shape).
        valid_shape = x_shape - y_shape + 1
        if np.any(valid_shape < 1):
            valid_shape = y_shape - x_shape + 1
            if np.any(valid_shape < 1):
                raise ValueError("empty result for valid shape")
        start = (z_shape - valid_shape) // 2
        end = start + valid_shape
        z = z[start[0]:end[0], start[1]:end[1]]

    return z





def fg_2d_data(spectrum,spectrum_mask,smooth_kernel=[],sigma=3,dopower='',use_masked_array=[]):
    """
    RFI flagging via 2d-convolution
    """

    # convolution f with g can also be seen as the Multiplication of their Fourier components:
    #
    # Fourier pairs (convolution =*, multiplication = x):    f * g    = F(f) x F(g)
    #                                                        F(f x g) = F(f) * F(g)
    # Mask_Array = Array x Mask
    # 
    # Convolving a masked array with a filter:  Mask_Array * Filter = F(Mask_Array) x F(Filter) = F(Array x Mask) x F(Filter) = (F(Array) * F(Mask)) x F(Filter) 
    #
    # So essentially to do a convolution of a masked array you need to first convolve the unmasked array with the mask and multiply it with the
    # Fouriertransform of the Filter.


    from scipy.signal import convolve2d
    from copy import deepcopy

    total_flags = []

    if len(use_masked_array) == 0:
        use_masked_array = np.ones(len(smooth_kernel))

    # the convention of the masked array for the radio data is different that
    # the one for masked array in numpy
    # Keep in mind: data that is bad is marked as True (or 1) in NUMPY
    #
    new_smo_mask = np.invert(deepcopy(spectrum_mask).astype(bool))  

    # run over each filter
    for smk,usem in zip(smooth_kernel,use_masked_array):

        # 
        new_data_mask_ma_array  = new_smo_mask.astype(bool)

        # switch to either work on the fg data or clean data
        if usem == 1:
            fullmask_data = spectrum * np.invert(new_data_mask_ma_array).astype(int)
        else:
            fullmask_data = spectrum

        # convolution with filter
        conv_image              = convolve2d(fullmask_data,smk,mode='same',boundary='symm')

        if dopower == 'abs':
            smoim               = np.absolute(conv_image)
        elif dopower == 'ang':
            smoim               = np.angle(conv_image)
        else:
            smoim               = conv_image

        smoim_ma                = ma.masked_array(smoim,mask=new_data_mask_ma_array,fill_value=np.nan)


        # Here the statistics 
        # improvement via kde needed
        #
        stats     = [smoim_ma.mean(),smoim_ma.std()] 
        #stats     = kdemean(smoim_ma.compressed())

        # do a selection on pure statistics --> here select bad data 
        sel           = np.logical_or(smoim_ma < stats[0] - sigma * stats[1], smoim_ma > stats[0] + sigma*stats[1])
        new_mask      = np.zeros(spectrum_mask.shape)
        new_mask[sel] = 1

        # combine the masks
        new_smo_mask = np.logical_or(new_mask.astype(bool),new_data_mask_ma_array.astype(bool))

        # estimate the total flags in the mask
        total_flags.append(np.count_nonzero(new_smo_mask))

    return np.invert(new_smo_mask).astype(int), total_flags 


def fg_2d_data_spwd(spectrum_data,spectrum_mask,smooth_kernel=[],sigma=3,dopower='',use_masked_array=[],spwd=2**0,nodc=True):
    """
    RFI flagging via 2d-convolution

    this just organises fg_2d_data 
    """

    spectrum_mask_org = spectrum_mask 

    dc_offset      = 0
    if nodc:
        spectrum_data_no_dc = spectrum_data[:,1:]
        spectrum_mask_no_dc = spectrum_mask[:,1:]
    
        spectrum_data  = spectrum_data_no_dc
        spectrum_mask  = spectrum_mask_no_dc
        dc_offset      = -1
    

    # generates an index array
    channel_idx        = np.arange(spectrum_mask.shape[1])
    channel_idx_spwd   = convert_1d_spwd(channel_idx,spwd)

    # convert the data into spwd
    spectrum_spwd      = convert_2d_spwd(spectrum_data,spwd)
    spectrum_mask_spwd = convert_2d_spwd(spectrum_mask,spwd)

    dims               = spectrum_spwd.shape

    # Flag each spectral window
    #
    fg_info_spwd = []
    for sw in range(dims[1]):

        data_spwd = spectrum_spwd[:,sw,:]
        mask_spwd = spectrum_mask_spwd[:,sw,:]

        idx_spwd  = channel_idx_spwd[sw]
        
        mask, fginfo   = fg_2d_data(data_spwd,mask_spwd,smooth_kernel,sigma,dopower,use_masked_array)
        mask_zero_indx = np.argwhere(mask == 0)
        
        # info for book keeping
        fg_info_spwd.append(fginfo)

        # write the spwd flags into the original mask
        for k in mask_zero_indx:
            
            spectrum_mask_org[k[0],channel_idx_spwd[sw,k[1]]+dc_offset] = 0
            #print(k,channel_idx_spwd[sw,k[1]]+dc_offset)

    # Flag DC term
    if nodc:
        spectrum_mask_org[:,0] = 0

    return spectrum_mask_org, fg_info_spwd 



def convert_1d_spwd(data,spwd=2**0,axis=0):
    """
    reshape the list of data
    """
    return data.reshape(int(data.shape[axis]/spwd),spwd)

def convert_2d_spwd(data,spwd=2**0):
    """
    reshape the array of data
    """
    return data.reshape(data.shape[0],int(data.shape[1]/spwd),spwd)


def complete_fg_mask(mask,sigma=3,axis=0,complete_boundary=9):
    """
    the idea is to check the number of FG in the mask
    and see if a channels range should be flaged
    """
    max_fg  = mask.shape[axis]
    fg_sum  = mask.sum(axis=axis)

    fg_axis    = np.arange(len(fg_sum))
    data_stats = kdemean(fg_sum,accucary=int(max_fg))
   
    if data_stats[1] != 0:

        sel = fg_sum <= data_stats[0] - data_stats[1] * sigma

        complete_fgs = fg_axis[sel]

        if axis == 0:
            for fgc in range(len(complete_fgs)-1):
                if complete_fgs[fgc+1] - complete_fgs[fgc] < complete_boundary: 
                    mask[:,complete_fgs[fgc]:complete_fgs[fgc+1]] = 0
                else:
                    mask[:,fgc] = 0

        else:
            for fgc in complete_fgs:
                if complete_fgs[fgc+1] - complete_fgs[fgc] < complete_boundary: 
                    mask[complete_fgs[fgc]:complete_fgs[fgc+1],:] = 0
                else:
                    mask[fgc,:] = 0

    return mask

def fg_data(spec,spec_mask,smooth_type='hanning',smooth_kernel=3,sigma=3,do_info=False):
    """
    basic flagging applies scipy filter function to the data 
    """
    from scipy.signal import wiener,hanning,convolve,hamming,gaussian,medfilt


    new_data_mask_ma_array  = np.invert(spec_mask.astype(bool))
    fullmask_data           = ma.masked_array(spec,mask=new_data_mask_ma_array,fill_value=np.nan)


    if smooth_type == 'hanning':
        sm_kernel = hanning(smooth_kernel)
        sm_data   = convolve(fullmask_data,sm_kernel,mode='same') / sum(sm_kernel)

    if smooth_type == 'hamming':
        sm_kernel = hamming(smooth_kernel)
        sm_data   = convolve(fullmask_data,sm_kernel,mode='same') / sum(sm_kernel)

    if smooth_type == 'gaussian':
        sm_kernel = gaussian(smooth_kernel,smooth_kernel)
        sm_data   = convolve(fullmask_data,sm_kernel,mode='same') / sum(sm_kernel)

    if smooth_type == 'median':
         sm_data = medfilt(fullmask_data,smooth_kernel)

    if smooth_type == 'wiener':
         sm_data = wiener(fullmask_data,smooth_kernel)


    # here divide the original data with the smoothed one
    n_data    = spec/sm_data
   
    # check for outliers
    #
    data_mask = get_good_data(n_data,spec_mask,sigma,stats_type='kdemean',do_info=do_info)

    return data_mask


def kdemean(x,accucary=1000,doinfo=False):
    """
     use the Kernel Density Estimation (KDE) to determine the mean
    
    (http://jpktd.blogspot.com/2009/03/using-gaussian-kernel-density.html )
    """
    from scipy.stats import gaussian_kde
    from numpy import linspace,min,max,std,mean
    from math import sqrt,log
    
    if mean(x) == std(x):
            if doinfo:
                    print('kde mean = std')
            return(mean(x),std(x))

    if std(x) == 0:
            if doinfo:
                    print('std = 0 ')
            return(mean(x),std(x))

    max_range = max(np.abs([min(x),max(x)]))

    # create instance of gaussian_kde class
    gk     = gaussian_kde(x)

    vra    = linspace(-1*max_range,max_range,accucary)
    vraval = gk.evaluate(vra)

    # get the maximum
    #
    x_value_of_maximum = vra[np.argmax(vraval)]

    # Devide data
    difit = vraval / max(vraval)
    #
    # and select values from 0.5
    sel = difit >= 0.4999

    idx_half_power = list(difit).index(min(difit[sel]))

    if idx_half_power >= accucary -1:
        return(mean(x),std(x))

    delta_accuracy = max([abs(vra[idx_half_power-1] - vra[idx_half_power]),\
                              abs(vra[idx_half_power+1] - vra[idx_half_power])])

    fwhm = abs(x_value_of_maximum - vra[idx_half_power])


    # factor 2 is because only one side is evaluated
    sigma = 2*fwhm/(2*sqrt(2*log(2)))

    # safety net
    # is the KDE is not doing a good job
    #
    if sigma > std(x):
        return(mean(x),std(x))

    return(x_value_of_maximum,abs(sigma)+delta_accuracy)


def mask_merger(orgmask,applymask):
    """
    merge two masks even with different size 
    """
    if orgmask.shape == applymask.shape:
        new_mask = orgmask * applymask
    elif orgmask.shape[0] == applymask.shape[0]:
        new_mask = np.multiply(orgmask,applymask[:,np.newaxis])
    elif orgmask.shape[1] == applymask.shape[0]:
        new_mask = np.transpose(np.multiply(np.transpose(orgmask),applymask[:,np.newaxis]))
    else:
        print('CAUTION MASK COULD NOT BE APPLIED MASK IS ORIGINAL MASK')
        new_mask = orgmask
    return new_mask


def combine_masks(orgmask,listofmask):
    """
    combine all mask into a final on
    """
    new_mask = orgmask
    for k in range(len(listofmask)):
        new_mask = mask_merger(new_mask,listofmask[k])

    return new_mask


def dchannel_axis(obsfile,data_key):
    """
    Generate a delta channel frequency axis  
    """
    frequency_axis = obsfile[data_key+'frequency']
    
    d_channel      = np.diff(frequency_axis[1:])

    d_channel = np.insert(d_channel, 0, 0)        # to cope with the DC term 
    d_channel = np.append(d_channel, d_channel[-1])

    return d_channel



def squash_mask(full_mask,axis,keyname='SQMASK'):
    """
    generate the mask and the selection mask for a specified axis
    axis = 1 would generate time masks
    axis = 0 would generate channel masks
    """

    jmask = {}
    jmask[keyname] = {}
  
    # use numpy masked array environment   (CAUTION: 1 or True indicates BAD DATA)
    # 
    new_data_mask_ma_array  = np.invert(full_mask.astype(bool))
    #fullmask_data           = ma.masked_array(full_mask,mask=new_data_mask_ma_array,fill_value=np.nan)
    fullmask_data           = ma.masked_array(full_mask,mask=new_data_mask_ma_array,fill_value=1)

    # sum over the specified axis and 
    full_sq_mask        = fullmask_data.sum(axis=axis)
    selection_sq_mask   = full_sq_mask > 0

    blank_full_sq_mask  = np.zeros(full_sq_mask.shape)
    blank_full_sq_mask[selection_sq_mask] = 1
    
    jmask[keyname]['mask']    = blank_full_sq_mask
    jmask[keyname]['masksum'] = full_sq_mask
    jmask[keyname]['selmask'] = selection_sq_mask

    return jmask


def integrate_data(full_data,full_mask,integrate_axis,sum_or_mean='sum',fill=-1):
    """
    interate the power spectrum
    axis = 1 would generate time masks
    axis = 0 would generate channel masks
    """
  
    # use numpy masked array environment   (CAUTION: 1 or True indicates BAD DATA)
    # 
    new_data_mask_ma_array  = np.invert(full_mask.astype(bool))
    fullmask_data           = ma.masked_array(full_data,mask=new_data_mask_ma_array,fill_value=np.nan)

    if sum_or_mean == 'sum':
        if fill == -1:
            # sum over the specified axis and 
            full_sq_data        = fullmask_data.sum(axis=integrate_axis)
        else:
            full_sq_data        = fullmask_data.filled(fill).sum(axis=integrate_axis)

    else:
        # mean over the specified axis and 
        full_sq_data        = fullmask_data.mean(axis=integrate_axis)

    return full_sq_data

    


def calc_power(obsfile,data_key,final_mask,final_mask_time,final_mask_channels):
    """
    calculate the power over the entire bandwidth
    """

    #
    # Determine the selected bandwidth
    #
    delta_freq_per_channel      = np.diff(obsfile[data_key+'frequency'][final_mask_channels[data_key]['selmask']])[0]
    number_of_selected_channels = obsfile[data_key+'frequency'][final_mask_channels[data_key]['selmask']].shape[0] 

    selected_frequency_width    = delta_freq_per_channel * number_of_selected_channels 

    #
    # Determine the power of the measure 
    #
    summed_spectrum    = squash_data(obsfile[data_key+'spectrum'],final_mask[data_key],axis=1,sum_or_mean='sum')
    integration_time   = obsfile[data_key+'integration_time'][:].flatten()

    power = (summed_spectrum / integration_time ) * selected_frequency_width 

    return power


def calc_power_spectrum(obsfile,data_key):
    """
    calculate the power spectrum
    counts/s * Hz
    """

    #
    # Determine the channel width 
    #
    mask_DC_chan           = np.ones(len(obsfile[data_key+'frequency'][:]))
    mask_DC_chan[0]        = 0
    sel_exclude_DC_chan    = mask_DC_chan > 0
    #
    delta_freq_per_channel = np.diff(obsfile[data_key+'frequency'][sel_exclude_DC_chan])[0]
    #
    # Determine the power 
    #
    integration_time      = obsfile[data_key+'integration_time'][:].flatten()
    #
    #
    power_spec            = np.divide(obsfile[data_key+'spectrum'],integration_time[:,np.newaxis]) * delta_freq_per_channel 

    return power_spec


def old_calc_power_spectrum(obsfile,data_key,final_mask,sum_mean=''):
    """
    calculate the power spectrum
    """

    #
    # Determine the channel width 
    #
    mask_DC_chan        = np.ones(len(obsfile[data_key+'frequency'][:]))
    mask_DC_chan[0]     = 0
    sel_exclude_DC_chan = mask_DC_chan > 0
    #
    delta_freq_per_channel      = np.diff(obsfile[data_key+'frequency'][sel_exclude_DC_chan])[0]

    #
    # Determine the power 
    #
    integration_time      = obsfile[data_key+'integration_time'][:].flatten()
    #
    power_spec            = np.divide(obsfile[data_key+'spectrum'],integration_time[:,np.newaxis]) * delta_freq_per_channel 


    if sum_mean == 'mean' or sum_mean == 'sum':
        #
        # Take the mean value of the spectrum 
        #
        power_spectrum   = squash_data(power_spec,final_mask[data_key],axis=0,sum_or_mean=sum_mean)

    else:

        # use numpy masked array environment   (CAUTION: 1 or True indicates BAD DATA)
        # 
        new_data_mask_ma_array  = np.invert(final_mask[data_key].astype(bool))
        fullmask_data           = ma.masked_array(power_spec,mask=new_data_mask_ma_array,fill_value=np.nan)

        power_spectrum          =  fullmask_data 


    return power_spectrum




def derivate_data(obsfile,key,toderivate,dertype='abs',sigma=3,accur=1000):
    #
    if dertype == 'abs':
        d_dt = list( np.abs((np.diff(obsfile[key+'/'+toderivate][:].flatten(),prepend=0)) / np.diff(obsfile[key+'/timestamp'][:].flatten(),prepend=0)))
    else:
        d_dt = list( (np.diff(obsfile[key+'/'+toderivate][:].flatten(),prepend=0)) / np.diff(obsfile[key+'/timestamp'][:].flatten(),prepend=0))


    stats  = kdemean(d_dt,accucary=accur)

    cuts1  = stats[0] - sigma * stats[1]
    cuts2  = stats[0] + sigma * stats[1]

    do_plot_hist = False
    if do_plot_hist:
        import matplotlib.pyplot as plt
        print(toderivate,cuts1,cuts2)
        fig, ax = plt.subplots()
        plt.title(toderivate)
        x,bins,patch = ax.hist(d_dt,density=True,bins=accur)
        ax.plot([cuts1,cuts1],[0,max(bins)],'-,','b')
        ax.plot([cuts2,cuts2],[0,max(bins)],'-,','b')
        plt.show()

        sys.exit(-1)

    seldata = np.logical_and(d_dt > cuts1,d_dt < cuts2)
    
    d_dt_mask  = np.zeros(len(seldata))
    d_dt_mask[seldata] = 1

    d_dt_mask[0]       = 0  # erase first entry, due to np.diff

    return d_dt, d_dt_mask, seldata, stats, cuts1, cuts2

def obs_direction(obsfile,key,accur=1000):
    """
    determine the scan direction of an observation
    ra or dec
    """

    tocheck   = ['ra','dec']
    dra_data  = np.diff(obsfile[key+'/'+tocheck[0]][:].flatten()) / np.diff(obsfile[key+'/timestamp'][:].flatten())
    ddec_data = np.diff(obsfile[key+'/'+tocheck[1]][:].flatten()) / np.diff(obsfile[key+'/timestamp'][:].flatten())

    ra_bins,ra_bvalues   = np.histogram(dra_data,accur)
    dec_bins,dec_bvalues = np.histogram(ddec_data,accur)

    otf_t     = [np.abs(ra_bvalues[np.argmax(ra_bins)]),np.abs(dec_bvalues[np.argmax(dec_bins)])]
    
    return tocheck[np.argmax(otf_t)]



def fov_fwhm(obsfreq,diameter,type='fov',outunit='deg'):
    """
    return the fwhm (full width at halve maximum) or fov (field of view)
    input: obsfreq [Hz], dimater in [m]
    Synthesis imaging 1999; P.J. Napier page 41
    """
    from astropy import constants as const
    from math import pi


    if type == 'fov':
        # position of first null is 1.22 
        # fov is twice this value
        value = 2 * 1.22 * (const.c.value / obsfreq) / diameter 
    else:
        # half power beam width
        value = 1.02 * (const.c.value / obsfreq) / diameter 

    # convert radians into degree
    #
    if outunit == 'arcmin':
        factor = 180./pi * 60
    elif outunit == 'arcsec':
        factor = 180./pi * 3600
    else:
        factor = 180./pi

    return value * factor

def sigma2fwhm(sigma):
    """
    convert sigma of a Gaussian to FWHM 

    for a Gaussian funktion, sigma is defined as distance from zero (so to get the full 
    beam need to multiply by 2)

    0.5 = e^(- (WHM**2 )/ (2 sigma**2))
    ln(0.5) = - (WHM**2 )/ (2 sigma**2))
    sqrt(2 * ln(2)) = (WHM/sigma)
    with FWHM = 2 * WHM
    2 * sqrt(2 * ln(2)) = FWHM/sigma
    sqrt( 8 * ln(2)) = 2.35482 = FWHM/sigma
    """
    from math import sqrt,log

    return sqrt( 8 * log(2) ) * sigma

def beaminpixel(Bmaj,Bmin,ddeginpixel):
    """
    Beam area = pi * r**2

    r = sigma

    For a Gaussian function, sigma is defined as distance from zero (so to get the full 
    beam need to multiply by 2)

    = 2 * (pi * sigma**2)      # the motivation of the 2 is related to the flux density integration 
                               # and the ben; the total integration of a Gauss to the FWHM and to 
                               # infinity differs by a factor of 2  
    = 2 * pi * sigma**2 
    = 2 * pi * ( WHM/(sqrt( 2 * ln(2))) )**2 
    = (pi / ln(2))  * WHM**2

    BMaj = 2 * WHM

    = pi / ( 4 * ln(2))  * Bmaj**2
    = 1.13309 * (BMaj * BMin) 

    [arcsec**2] to convert this into pixels divide this by the cellsize squared (arcsec/pixel)**2

    # convert the integrated flux density in an image int a gaussian function

    image is [beam / pixel] to obtain the integrated flux density, sum up all the pixel multiply by the beaminpixel and divide to the 
    number of pixels
    """
    from math import pi, log

    pixelsperbeam = bmaj/(ddegpix) * bmin/(ddegpix) * (pi / ( 4 * log(2) ))

    return pixelperbeam

def filenamecounter(fname,extention='.png'):
    """
    just provide a filename taken exsisting files with the
    same name in the directory into accout
    """
    import os

    filename = fname+extention

    if os.path.isfile(fname+extention):

        counter = 0
        filename = fname+'_'+'{}'+extention
        while os.path.isfile(filename.format(counter)):
            counter += 1
        filename = filename.format(counter)

    return filename
    
def get_obs_info(timestamp_keys,info_idx=0,splittype='/'):
    """
    
    """

    obs_informa = []
    for k in timestamp_keys:
        obs_informa.append(k.split(splittype)[info_idx])
        
    obs_info  = np.unique(obs_informa)

    return obs_info


def source_plant_separation(source_ra,source_dec,planet,time,obs_pos):
    """
        
    """
    from astropy.coordinates import get_body,EarthLocation,SkyCoord
    from astropy.time import Time
    import astropy.units as u


    array_central_pos = EarthLocation(obs_pos[0]*u.deg,obs_pos[1]*u.deg,obs_pos[2]*u.m)
    timestamp         = Time(time,scale='utc',format='unix')
    planet_pos        = get_body(planet,timestamp,location=array_central_pos, ephemeris=None)
    planet_sky_pos    = SkyCoord(ra=planet_pos.ra,dec=planet_pos.dec,frame='gcrs',location=array_central_pos,obstime=timestamp)
    so_coord          = SkyCoord(source_ra,source_dec, unit='deg',frame='icrs')
    separation        = so_coord.separation(planet_sky_pos).deg

    return separation
