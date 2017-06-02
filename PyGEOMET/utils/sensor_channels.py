#This function returns a list of the available channels
# for the selected sensor
#Also creates a view list with wavelength information
#Note CRTM separates visible and IR spectrums for
# all sensors
#Input:  CRTM sensor
#Output:  channel - input to CRTM
#         names - what the user sees 
def getSensorChannels(sensor):
    
    #GOES-12 Imager (IR Channels)
    if (sensor == 'imgr_g12' or sensor == 'imgr_g13' or
        sensor == 'imgr_g14' or sensor == 'imgr_g15'):
        channels = ['2','3','4','6']
        names = ['2 (3.9 $\mu m$)', '3 (6.5 $\mu m$)',
                 '4 (10.7 $\mu m$)', '6 (13.3 $\mu m$)']
        
    #GOES-16 ABI (IR Channels)    
    elif (sensor == 'abi_gr'):
        channels = ['7','8','9','10',
                    '11','12','13','14',
                    '15','16']
        names = ['7 (3.9 $\mu m$)', '8 (6.2 $\mu m$)',
                 '9 (6.9 $\mu m$)', '10 (7.3 $\mu m$)',
                 '11 (8.4 $\mu m$)', '12 (9.6 $\mu m$)',
                 '13 (10.3 $\mu m$)', '14 (11.2 $\mu m$)',
                 '15 (12.3 $\mu m$)', '16 (13.3 $\mu m$)']        
        
    return channels, names    