#
#
# Hans-Rainer Kloeckner
#
# MPIfR 2024
#
#
# This is an example to get you a head start using a HDF5 files of SKAMPI  
#
# ---------------------------------------------------------------------------------------------
# Generate some spectrum, waterfall and sky postoin plot
# ---------------------------------------------------------------------------------------------
#
#
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------
#
import h5py
import numpy as np
import numpy.ma as ma
import sys
import copy
from optparse import OptionParser
#
from MPG_HDF5_libs import *
#
from astropy.time import Time
from astropy import units as u
#
import matplotlib as plt
import matplotlib.cm as cm
#
from collections import OrderedDict
#

def main():

    # argument parsing
    #
    print('\n== Example how to access and plot SKAMPI Data == \n')

    # ----
    parser       = new_argument_parser()
    (opts, args) = parser.parse_args()

    if opts.datafile == None or opts.help:
        parser.print_help()
        sys.exit()


    # set the parameters
    #
    data_file                 = opts.datafile
    do_use_mask               = opts.use_mask
    doplot_final_spec         = opts.doplot_final_spec
    doplot_final_waterfall    = opts.doplot_final_full_data
    doplot_flag_sky_positions = opts.doplot_sky_positions
    pltsave                   = opts.pltsave


    # ###################################################
    #
    # Some plotting parameter 
    #
    import matplotlib as plt
    import matplotlib.cm as cm
    #
    #
    cmap = copy.copy(cm.cubehelix)
    cmap.set_bad(color='black')
    #
    im_size  = (8.27, 11.69)[::-1]  # A4 landscape
    #im_size  = (8.27, 11.69)       # A4 portrait
    #im_size  = (8, 4)       # A4 portrait
    #
    plt.rcParams['figure.figsize'] = im_size
    #
    #
    DPI = 150
    #
    # ###################################################



    # load the hdf5 file
    #
    #
    obsfile = h5py.File(data_file)

    print('\n=== Investigate masking of the data in ',data_file)


    # Example to optain timestamp data of the file
    #
    #
    # Note that the time stamps of the individual units have different units
    # e.g. timestamps from the acu, the weather station, or the backend 

    find_key = ['timestamp']

    print('\n\t=== extract ',find_key,' ===\n')

    data_key = findkeys(obsfile,keys=find_key,exactmatch=True)

    for d in data_key:
        print('\tkey  : ',d)
        print('\tshape: ',obsfile[d].shape)
        #print('\tfirst data entry: ',obsfile[d][:],'\n')
    print('\t++++++++++++++++++++++++++++\n')



    # get the data keys of the timestamping in the scans and spectra
    #
    timestamp_keys       = findkeys(obsfile,keys=['scan','timestamp'],exactmatch=True)
    #
    spectrum_keys        = findkeys(obsfile,keys=['scan','spectrum'],exactmatch=True)
    # #########  
    #


    # ---------------------------------------------------------------------------------------------
    # Generate a mask to synchronise the data 
    # ---------------------------------------------------------------------------------------------
    #
    print('\n\n   === Syncronise PFB data in time ===\n')
    #
    #
    # the funktion returns the individual mask for each key to assure 
    # the same time stamping of the data
    #
    time_sync_mask       = return_equal_data(obsfile,timestamp_keys)
    #
    # ---------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------
    # Generate a python dic with either the mask of the file or a blank mask
    # ---------------------------------------------------------------------------------------------
    #
    #
    data_mask = {}
    #
    if do_use_mask:

        for s in spectrum_keys:
            mask  = obsfile[s.replace('spectrum','mask')][:].astype(int)
            data_mask[s.replace('/spectrum','')] = mask

    else:

        #
        # generate a dummy for final combination 
        #
        for s in spectrum_keys:
            # generate a blank channel mask
            #
            # get the data
            mask  = obsfile[s.replace('spectrum','mask')][:] 

            # generates a blank mask
            blank_full_mask = np.ones(mask.shape)

            # note that the full dataset includes the DC component of the FFT at index 0 
            # will masked the DC term
            #
            blank_full_mask[:,0] = 0

            # the last channel is always quit high
            #
            blank_full_mask[:,mask.shape[1]-1] = 0

            data_mask[s.replace('/spectrum','')] = blank_full_mask

    #
    # ---------------------------------------------------------------------------------------------


    #
    # Combine the masks into a final one
    #
    final_mask = {}
    #
    for d in timestamp_keys:
            print('\tsync full 2-d mask in time for : ',d.replace('timestamp',''))

            # combine the sync time mask and the full mask
            #
            sync_mask = combine_masks(data_mask[d.replace('/timestamp','')],[np.array(time_sync_mask[d]['mask'])])

            # new set of masks
            #
            final_mask[d.replace('timestamp','')]      = sync_mask

    #
    # ---------------------------------------------------------------------------------------------

  

    # ---------------------------------------------------------------------------------------------
    # Plot the spectrum of the data set
    # ---------------------------------------------------------------------------------------------
    #

    if doplot_final_spec:

        print('\n   === Generate 1d Spectrum plots === \n')

        # Here does the plotting of the data
        #
        import matplotlib.pyplot as plt
        import matplotlib
        #
        for d in timestamp_keys:
            print('\tgenerate plot for : ',d.replace('timestamp',''))

            spectrum_data  = obsfile[d.replace('timestamp','')+'spectrum'][:] 

            sq_masked      = squash_mask(final_mask[d.replace('timestamp','')],0,keyname=d.replace('timestamp',''))
            sq_masked_sum  = sq_masked[d.replace('timestamp','')]['masksum']
            sq_masked_sel  = sq_masked[d.replace('timestamp','')]['selmask']
            channels       = np.arange(len(sq_masked_sum))
            freqs          = obsfile[d.replace('timestamp','')+'frequency'][:]

            spectrum_sum   = integrate_data(spectrum_data,final_mask[d.replace('timestamp','')],integrate_axis=0,sum_or_mean='mean',fill=-1)


            # print the spectrum
            fig, ax = plt.subplots()
            plt.title(d.replace('timestamp',''))
            #ax.scatter(channels[sq_masked_sel],spectrum_sum[sq_masked_sel])
            ax.scatter(freqs[sq_masked_sel],spectrum_sum[sq_masked_sel])
            ax.set_xlabel('frequency')
            ax.set_ylabel('mean spectrum [power, arbitrary units]')
            if pltsave:
                plt_fname = data_file.replace('.hdf5','').replace('.HDF5','')+'_'+d.replace('timestamp','').replace('/','_')+'SPEC'
                plt_fname = filenamecounter(plt_fname,extention='.png')
                fig.savefig(plt_fname,dpi=DPI)
            else:
                plt.show()
        plt.clf()


    # ---------------------------------------------------------------------------------------------
    # Plot the waterfall spectrum of the data set
    # ---------------------------------------------------------------------------------------------
    #
    if doplot_final_waterfall:

        print('\n   === Generates waterfall  plots === \n')

        # Here does the plotting of the data
        #
        import matplotlib.pyplot as plt
        import matplotlib
        #
        for d in timestamp_keys:
            print('\tgenerate plot for : ',d.replace('timestamp',''))

            spectrum_data  = obsfile[d.replace('timestamp','')+'spectrum'][:] 

            # print the waterfall plot
            #
            merged_edge_mask_invert = np.invert(final_mask[d.replace('timestamp','')].astype(bool))
            fullmask_data           = ma.masked_array(spectrum_data,mask=merged_edge_mask_invert,fill_value=np.nan)
            stats                   = [fullmask_data.mean(),fullmask_data.std()]

            fig, ax = plt.subplots()
            plt.title(d.replace('timestamp',''))
            wfplt = ax.imshow(fullmask_data,interpolation='nearest',origin='lower',vmin=stats[0]-3*stats[1],vmax=stats[0]+3*stats[1])
            ax.set_xlabel('channels')
            ax.set_ylabel('time')

            if pltsave:
                plt_fname = data_file.replace('.hdf5','').replace('.HDF5','')+'_'+d.replace('timestamp','').replace('/','_')+'WFPLT'
                plt_fname = filenamecounter(plt_fname,extention='.png')
                fig.savefig(plt_fname,dpi=DPI)
            else:
                plt.show()

        plt.clf()


    #
    #
    # ---------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------
    # Plot the sky and flagged positions
    # ---------------------------------------------------------------------------------------------
    #
    if doplot_flag_sky_positions:

        print('\n   === Generates sky coverage plots === \n')

        # Here does the plotting of the data
        #
        import matplotlib.pyplot as plt
        import matplotlib
        #
        for d in timestamp_keys:
            print('\tgenerate plot for : ',d.replace('timestamp',''))
            

            # determine the spectrum average to be used for colouring 
            #
            spectrum_data  = obsfile[d.replace('timestamp','')+'spectrum'][:] 

            #
            merged_edge_mask_invert = np.invert(final_mask[d.replace('timestamp','')].astype(bool))
            fullmask_data           = ma.masked_array(spectrum_data,mask=merged_edge_mask_invert,fill_value=np.nan)
            stats                   = [fullmask_data.mean(axis=1),fullmask_data.std(axis=1)]

            # get the coordinates 
            #
            ra = obsfile[d.replace('timestamp','ra')][:]
            dec = obsfile[d.replace('timestamp','dec')][:]

            # for better plotting results we sort the data
            #
            sort_it  = np.argsort(stats[0])
            ra_sort  = ra[sort_it]
            dec_sort = dec[sort_it]
            power    = stats[0][sort_it]

            # scatter plot
            #
            fig, ax = plt.subplots()
            plt.title(d.replace('timestamp',''))

            scplt = ax.scatter(ra_sort,dec_sort,c=power)
            ax.set_xlabel('ra')
            ax.set_ylabel('dec')
            fig.colorbar(scplt,ax=ax)

            if pltsave:
                plt_fname = data_file.replace('.hdf5','').replace('.HDF5','')+'_'+d.replace('timestamp','').replace('/','_')+'SKYPLT'
                plt_fname = filenamecounter(plt_fname,extention='.png')
                fig.savefig(plt_fname,dpi=DPI)
            else:
                plt.show()

        plt.clf()

    #
    #
    # ---------------------------------------------------------------------------------------------




def new_argument_parser():

    #
    # some input for better playing around with the example
    #
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)


    parser.add_option('--DATA_FILE', dest='datafile', type=str,
                      help='DATA - HDF5 file of the Prototyp')

    parser.add_option('--DOPLOT_FINAL_SPEC', dest='doplot_final_spec', action='store_true',
                      default=False,help='Plot the final spectrum after Flagging')

    parser.add_option('--DOPLOT_FINAL_WATERFALL', dest='doplot_final_full_data', action='store_true',
                      default=False,help='Plot the final waterfall after Flagging')

    parser.add_option('--DOPLOT_ON_SKY', dest='doplot_sky_positions', action='store_true',
                      default=False,help='Plot the final flags in ra and dec')

    parser.add_option('--DOSAVEPLOT', dest='pltsave', action='store_true',
                      default=False,help='Save the plots as figures')

    parser.add_option('--DOUSEMASK', dest='use_mask', action='store_true',
                      default=False,help='Use the default mask to plot the data')

    parser.add_option('--HELP', dest='help', action='store_true',
                      default=False,help='Show info on input')

    return parser


if __name__ == "__main__":
    main()
