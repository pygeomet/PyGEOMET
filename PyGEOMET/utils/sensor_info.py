#This function returns a list of the available channels
# for the selected name and gets the CRTM name name
#Also creates a view list with wavelength information
#Note CRTM separates visible and IR spectrums for
# all names
#Input :  Sensor name - what user sees
#Output:  CRTM name name
#         channel - input to CRTM
#         names - what the user sees 
def getSensorInfo(sname):
    
    #GOES 12-15 Imager (IR Channels)
    if (sname == 'GOES-15 Imager (IR)' or sname == 'GOES-14 Imager (IR)' or
        sname == 'GOES-13 Imager (IR)' or sname == 'GOES-12 Imager (IR)'):
        #Get CRTM sensor name
        if (sname == 'GOES-15 Imager (IR)'):
            sensor = 'imgr_g15'
        elif (sname == 'GOES-14 Imager (IR)'): 
            sensor = 'imgr_g14'
        elif (sname == 'GOES-13 Imager (IR)'):
            sensor = 'imgr_g13'
        elif (sname == 'GOES-12 Imager (IR)'):
            sensor = 'imgr_g12'
        
        channels = ['2','3','4','6']
        cnames = ['2 (3.9 '+u'\u03BC'+'m)', '3 (6.5 '+u'\u03BC'+'m)',
                  '4 (10.7 '+u'\u03BC'+'m)', '6 (13.3 '+u'\u03BC'+'m)']

    #GOES 12-15 (Visible Channel)
    elif (sname == 'GOES-15 Imager (Vis)' or sname == 'GOES-14 Imager (Vis)' or
          sname == 'GOES-13 Imager (Vis)' or sname == 'GOES-12 Imager (Vis)'):
        #Get CRTM sensor name
        if (sname == 'GOES-15 Imager (Vis)'):
            sensor = 'v.imgr_g15'
        elif (sname == 'GOES-14 Imager (Vis)'):
            sensor = 'v.imgr_g14'
        elif (sname == 'GOES-13 Imager (Vis)'):
            sensor = 'v.imgr_g13'
        elif (sname == 'GOES-12 Imager (Vis)'):
            sensor = 'v.imgr_g12'

        channels = ['1']
        cnames = ['1 (0.65 '+u'\u03BC'+'m)']
 
    #GOES-16 ABI (IR Channels)    
    elif (sname == 'GOES-16 ABI (IR)'):
        sensor = 'abi_gr'
        channels = ['7','8','9','10',
                    '11','12','13','14',
                    '15','16']
        cnames = ['7 (3.9 '+u'\u03BC'+'m)', '8 (6.2 '+u'\u03BC'+'m)',
                  '9 (6.9 '+u'\u03BC'+'m)', '10 (7.3 '+u'\u03BC'+'m)',
                  '11 (8.4 '+u'\u03BC'+'m)', '12 (9.6 '+u'\u03BC'+'m)',
                  '13 (10.3 '+u'\u03BC'+'m)', '14 (11.2 '+u'\u03BC'+'m)',
                  '15 (12.3 '+u'\u03BC'+'m)', '16 (13.3 '+u'\u03BC'+'m)']        

    #GOES-16 ABI (Visible Channels)
    elif (sname == 'GOES-16 ABI (Vis)'):
        sensor = 'v.abi_gr'
        channels = ['1','2','3','4','5','6']
        cnames = ['1 (0.47 '+u'\u03BC'+'m)', '2 (0.64 '+u'\u03BC'+'m)',
                  '3 (0.86 '+u'\u03BC'+'m)', '4 (1.37 '+u'\u03BC'+'m)',
                  '5 (1.6 '+u'\u03BC'+'m)', '6 (2.2 '+u'\u03BC'+'m)',]
        
    return sensor, channels, cnames    
