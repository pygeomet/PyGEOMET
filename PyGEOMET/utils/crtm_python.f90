SUBROUTINE CRTM(press, theta, qv, qcloud, qice, qrain, qsnow, qgraupel, &
                qhail, ncloud, nice, nrain, nsnow, ngraupel, nhail, lai, u10, &
                v10, seaice, snowh, coszen, vegfrac, ptop, tsk, ivegtyp, xland, &
                landuse, mp_physics, lat, lon, sensor, channel, coeff_path, &
                request_var, crtm_out, ii, jj, kk) 

use crtm_module

implicit none

!Define Variables
!Input Variables
integer, intent(in) :: ii, jj, kk
real, intent(in) :: press(kk,jj,ii), theta(kk,jj,ii), qv(kk,jj,ii)
real, intent(in) :: qcloud(kk,jj,ii), qice(kk,jj,ii), qrain(kk,jj,ii)
real, intent(in) :: qsnow(kk,jj,ii), qgraupel(kk,jj,ii), qhail(kk,jj,ii)
real, intent(in) :: ncloud(kk,jj,ii), nice(kk,jj,ii), nrain(kk,jj,ii)
real, intent(in) :: nsnow(kk,jj,ii), ngraupel(kk,jj,ii), nhail(kk,jj,ii)
real, intent(in) :: lai(jj,ii), u10(jj,ii)
real, intent(in) :: v10(jj,ii), seaice(jj,ii), snowh(jj,ii)
real, intent(in) :: coszen(jj,ii), vegfrac(jj,ii), tsk(jj,ii), lon(jj,ii)
real, intent(in) :: ptop, ivegtyp(jj,ii), xland(jj,ii), lat(jj,ii)
character(len=10), intent(in) :: landuse, sensor
character(len=250), intent(in) :: coeff_path
character(len=50), intent(in) :: request_var
integer, intent(in) :: mp_physics, channel
real, intent(inout) :: crtm_out(jj,ii)

!Profile variables - CRTM
integer :: n_channels   
integer :: n_clouds
integer :: n_layers
integer, parameter :: n_sensors = 1
integer, parameter :: n_profiles = 1
integer, parameter :: n_absorbers = 2 !H20 and O3
integer, parameter :: n_aerosols = 0

!Processing parameters - CRTM
character(len=20) :: sensor_id(1)
type(crtm_channelinfo_type) :: chinfo(1)
type(crtm_geometry_type)    :: geo(1)
type(crtm_options_type)     :: opt(1)

!Forward declarations
type(crtm_atmosphere_type) :: atm(1)
type(crtm_surface_type)    :: sfc(1)
type(crtm_rtsolution_type),  allocatable :: rts(:,:) !multiple channels    

! Local variables
integer :: usgs_to_npoess(24), igbp_to_npoess(20)
integer :: vegtyp
real :: o3, veg_coverage
real :: lvl_pressure(kk+1), lay_pressure(kk), temp(kk)
!real :: rad(jj,ii), tb(jj,ii)
real, allocatable :: r_eff(:)
integer :: alloc_stat
integer :: err_stat
integer :: i, j, k, n, jday
real :: sat_lat, sat_lon, sat_zenith
real :: hour, sun_zenith, sun_azimuth
real :: land_coverage, water_coverage, snow_coverage, ice_coverage
real, parameter :: rd = 287.05
real, parameter :: cp = 1005.
real, parameter :: r2d = 180./3.141592654
real :: mass_lay

!<<<<<<<<<<<<<<< Start Program >>>>>>>>>>>>>>>!
!Land use conversions
!USGS to NPOESS
usgs_to_npoess = (/15,  1,  5, 11,  6,  6,  6,  7, 13, &
                    6,  8,  9,  8,  9, 12,  1, 18, 18, &
                    5, 10, 10, 10, 10,  1 /)

!IGBP (MODIS)
igbp_to_npoess = (/ 9,  8,  9,  8, 12,  7, 19, 17, 17, &
                    7, 17,  2, 15,  2,  1,  1,  1, 10, &
                   10, 10  /)

sensor_id(1) = sensor

!Initializae CRTM
!Load default NPOESS coefficient file
!Convert input to land classification to NPOESS
err_stat = crtm_init(sensor_id, chinfo, &
                     Load_CloudCoeff = .TRUE., &
                     Load_AerosolCoeff = .TRUE., &
                     IRlandCoeff_File = 'NPOESS.IRland.EmisCoeff.bin', &
                     IRwaterCoeff_File = 'WuSmith.IRwater.EmisCoeff.bin', &
                     MWwaterCoeff_File = 'FASTEM4.MWwater.EmisCoeff.bin', &
                     VISlandCoeff_File = 'NPOESS.VISland.EmisCoeff.bin', &
                     File_Path = trim(coeff_path), &
                     Quiet = .True.)
!Check error
if (err_stat /= SUCCESS) then
   write(6,*) 'Failed to initialize crtm = ', err_stat
endif

!Specify channel subsets for IR bands (channel 3&4)
!GOES-13,14,15 - 3 = water vapor (6.7 um), 4 = IR (10.7 um)
err_stat = crtm_channelinfo_subset( chinfo(1), &
                                    Channel_Subset = (/channel/) )
if ( err_stat /= SUCCESS ) then
   write(6,*) 'Failed to set channels = ', err_stat
endif

! specify numbers of cloud species
! Thompson==8, Ferrier==5,95, WSM6==6, Lin==2
if (mp_physics == 99) then ! Zhao Scheme
   n_clouds = 2 ! GFS uses Zhao scheme
else if( mp_physics == 5 .or. mp_physics == 85 .or. mp_physics == 95) then
   ! change to 6 cloud types because microwave is sensitive to density
   n_clouds = 6  
else if( mp_physics == 8 .or. mp_physics == 6 .or. &
         mp_physics == 2 .or. mp_physics == 28 .or. &
         mp_physics == 10) then
   n_clouds = 5
end if

! Set Ozone concentration
o3 = 0

!Set number of layers based on model input
n_layers = kk

!Allocate more CRTM variables
allocate (r_eff(n_clouds))

! Allocate the forward atmosphere structures
CALL crtm_atmosphere_create( atm , &
                             n_layers , &
                             n_absorbers, &
                             n_clouds , &
                             n_aerosols )
! Check they were created successfully
if ( ANY(.not. crtm_atmosphere_associated( atm ))) then
    write(6,*) 'Failed to create atmosphere for forward model'
endif

!Define atmosphere - CRTM level 1 is TOA
atm(1)%n_layers = n_layers
!Define absorbers
atm(1)%absorber_id(1) = H2O_ID
atm(1)%absorber_id(2) = O3_ID
!And the units
atm(1)%absorber_units(1) = MASS_MIXING_RATIO_UNITS
atm(1)%absorber_units(2) = VOLUME_MIXING_RATIO_UNITS

!Define Clouds based on microphysics
if (mp_physics == 99) then
    !Cloud Types
    atm(1)%cloud(1)%n_layers = n_layers
    atm(1)%cloud(1)%Type = WATER_CLOUD
    atm(1)%cloud(2)%n_layers = n_layers
    atm(1)%cloud(2)%Type = ICE_CLOUD
else if (mp_physics == 5 .or. mp_physics == 85 .or. &
         mp_physics == 95) then
    !Cloud Types
    atm(1)%cloud(1)%n_layers = n_layers
    atm(1)%cloud(1)%Type = WATER_CLOUD
    atm(1)%cloud(2)%n_layers = n_layers
    atm(1)%cloud(2)%Type = ICE_CLOUD
    atm(1)%cloud(3)%n_layers = n_layers
    atm(1)%cloud(3)%Type = RAIN_CLOUD
    atm(1)%cloud(4)%n_layers = n_layers
    atm(1)%cloud(4)%Type = SNOW_CLOUD
    atm(1)%cloud(5)%n_layers = n_layers
    atm(1)%cloud(5)%Type = GRAUPEL_CLOUD
    atm(1)%cloud(6)%n_layers = n_layers
    atm(1)%cloud(6)%Type = HAIL_CLOUD
else if (mp_physics == 8 .or. mp_physics == 28 .or. &
         mp_physics == 6 .or. mp_physics == 2 .or. &
         mp_physics == 10) then
    !Cloud Types
    atm(1)%cloud(1)%n_layers = n_layers
    atm(1)%cloud(1)%Type = WATER_CLOUD
    atm(1)%cloud(2)%n_layers = n_layers
    atm(1)%cloud(2)%Type = ICE_CLOUD
    atm(1)%cloud(3)%n_layers = n_layers
    atm(1)%cloud(3)%Type = RAIN_CLOUD
    atm(1)%cloud(4)%n_layers = n_layers
    atm(1)%cloud(4)%Type = SNOW_CLOUD
    atm(1)%cloud(5)%n_layers = n_layers
    atm(1)%cloud(5)%Type = GRAUPEL_CLOUD
end if
!Define first level pressure - TOA in hPA
atm(1)%Level_Pressure(0) = TOA_PRESSURE

!Loop around Grid
iloop: do i = 1, ii
  jloop: do j = 1, jj

    !Create pressure an temperature arrays - easier to complete in own loop
    do k = 1, n_layers
        lay_pressure(k) = press(k,j,i) 
        temp(k) = theta(k,j,i) * (lay_pressure(k)/1000.)**(rd/cp)
        !Will define level 1 later
        if (k > 1) then
          lvl_pressure(k) = (lay_pressure(k) + lay_pressure(k-1)) * 0.5 
        endif
        !Make sure pressure is decreasing
        if (k > 1 .and. lay_pressure(k) > lay_pressure(k-1)) then
          !write(6,*) "non-monotonic layer pressure", i, j, lay_pressure(k), &
          !           lay_pressure(k-1)
          lay_pressure(k-1) = lay_pressure(k) * 1.0001
        endif
        if (k > 2 .and. lvl_pressure(k) > lvl_pressure(k-1)) then
          !write(6,*) "non-monotonic level pressure", i, j, lvl_pressure(k), &
          !           lvl_pressure(k-1)
          lvl_pressure(k-1) = lvl_pressure(k) * 1.0001
        endif

    enddo
    !Set level pressure top and bottom
    lvl_pressure(n_layers+1) = ptop/100. !convert to hPa
    lvl_pressure(1) = lay_pressure(1) + (lvl_pressure(2) - lay_pressure(2)) 

    !Get the number of channels to process for current sensor
    n_channels = crtm_channelinfo_n_channels( chinfo(1))

    ! Allocate channel-dependent arrays
    allocate( rts(n_channels, n_profiles) , &
              stat = alloc_stat )
    !Check error
    if (alloc_stat /= 0) then
       write(6,*) 'Failed to allocate channel dependent &
                   arrays with stat = ', &
                   alloc_stat, 'sensor# = ', n
       call flush(6)
    endif

    ! Allocate the surface structures
    CALL crtm_surface_create( sfc , &
                              chinfo(1)%n_channels)
    ! Check they were created successfully
    if ( ANY(.not. crtm_surface_associated( sfc ))) then
        write(6,*) 'Failed to create surface for forward model'
        call flush(6)
        stop
    endif

    !Allocate the RTSolution structure
    CALL crtm_rtsolution_create( rts , &
                                 n_layers )
    ! Check they were created successfully
    if ( ANY(.not. crtm_rtsolution_associated( rts )) ) then
        write(6,*) 'Failed to create RT solution'
        stop
    end if

    !Define satellite location
    !GOES East (13)
    if (sensor_id(1)=='imgr_g13' .or. sensor_id(1)=='v.imgr_g13') then
       sat_lat = 0.0
       sat_lon = -75.0
    !Stand-by location
    else if (sensor_id(1)=='imgr_g14' .or. sensor_id(1)=='v.imgr_g14') then
       sat_lat = 0.0
       sat_lon = -105.0
    !GOES West (15)
    else if (sensor_id(1)=='imgr_g15' .or. sensor_id(1)=='v.imgr_g15') then
       sat_lat = 0.0
       sat_lon = -135.0
    !Stand-by location
    else if (sensor_id(1)=='abi_gr' .or. sensor_id(1)=='v.abi_gr') then
       sat_lat = 0.0
       sat_lon = -105.0
    endif

    !Calculate the solar zeinth angle and sat zenith angle
    !jday = 135
    !hour = 17 !17:45:39 UTC
    !CALL zenith(lat(i,j),lon(i,j),sat_lat,sat_lon,sat_zenith)
    !CALL solar_zenith(lati(i,j),lon(i,j),jday,hour,sun_zenith)
    !write(6,*) 'Viewing angle = ', sat_zenith
    !write(6,*) 'Solar zenith angle = ', sun_zenith

    !Set satellite viewing angle
    !geo%sensor_zenith_angle = sat_zenith
    !geo%sensor_scan_angle = 0 ! nadir
    !geo%source_zenith_angle = sun_zenith
    geo%sensor_zenith_angle = 0.
    geo%sensor_scan_angle = 0. ! nadir
    geo%source_zenith_angle = acos(coszen(j,i)) * r2d

    !Define surface sensor data - only necessary for microwave calculations
    sfc(1)%sensordata%n_channels = chinfo(1)%n_channels
    sfc(1)%sensordata%sensor_id = chinfo(1)%sensor_id
    sfc(1)%sensordata%WMO_sensor_id = chinfo(1)%WMO_sensor_id
    sfc(1)%sensordata%WMO_Satellite_id = chinfo(1)%WMO_Satellite_id
    !CRTM needs each sensor to have a WMO ID, GOES-16 doesn't have one
    ! therefore, set it
    if (sensor_id(1) == 'abi_gr' .or. sensor_id(1) == 'v.abi_gr') then 
        sfc(1)%sensordata%WMO_sensor_id = 615
        sfc(1)%sensordata%WMO_Satellite_id = 257 
    endif
    sfc(1)%sensordata%sensor_channel = chinfo(1)%sensor_channel

    !write(6,*) " "
    !write(6,*) "Channels", chinfo(n)%n_channels 
    !write(6,*) "ID", chinfo(n)%sensor_id
    !write(6,*) "WMO sensor id", chinfo(n)%WMO_sensor_id 
    !write(6,*) "WMO sat id", chinfo(n)%WMO_Satellite_id
    !write(6,*) "sensor channel", chinfo(n)%sensor_channel 

    !Only continue calculation if viewing angle is less than
    ! the maximum scaning angle
    if (geo(1)%sensor_zenith_angle <= MAX_SENSOR_SCAN_ANGLE) then

       !Determine surface coverage percentages
       if (xland(j,i) > 1.5) then !Water
         veg_coverage = 0.0
         !Check if surface is covered by snow
         if (snowh(j,i) > 0.1) then !snow covered
           snow_coverage = 1.0
           land_coverage = 0.0
           ice_coverage = 0.0
           water_coverage = 0.0
         else !water
           snow_coverage = 0.0 
           land_coverage = 0.0
           ice_coverage = 0.0
           water_coverage = 1.0
         endif
       else ! land or sea ice
         veg_coverage = vegfrac(j,i)/100. !convert to fraction
         !Check for sea ice 
         if (seaice(j,i) > 0) then !Ice and/or snow
           !Check if surface is covered by snow
           if (snowh(j,i) > 0.1) then !snow covered
             snow_coverage = 1.0
             land_coverage = 0.0
             ice_coverage = 0.0
             water_coverage = 0.0
           else !Ice
             snow_coverage = 0.0
             land_coverage = 0.0
             ice_coverage = 1.0
             water_coverage = 0.0
           endif
         else !land
           !Check if surface is covered by snow
           if (snowh(j,i) > 0.1) then !snow covered
             snow_coverage = 1.0
             land_coverage = 0.0
             ice_coverage = 0.0
             water_coverage = 0.0
           else !Ice
             snow_coverage = 0.0
             land_coverage = 1.0
             ice_coverage = 0.0
             water_coverage = 0.0
           endif
         endif
       endif
       
       !Assign surface type
       sfc(1)%Land_Coverage = land_coverage
       sfc(1)%Water_Coverage = water_coverage
       sfc(1)%Snow_Coverage = snow_coverage
       sfc(1)%Ice_Coverage = ice_coverage

       !Assign a land classification
       if (trim(landuse) .eq. 'USGS') then
           vegtyp = usgs_to_npoess(min(max(1.,ivegtyp(j,i)),24.))
       else !MODIS for now
           vegtyp = igbp_to_npoess(min(max(1.,ivegtyp(j,i)),20.))
       endif      
 
       sfc(1)%Land_Type = vegtyp

       !Assign more surface parameters
       sfc(1)%Wind_Speed = sqrt(u10(j,i)*u10(j,i)+v10(j,i)*v10(j,i))
       sfc(1)%Land_Temperature = tsk(j,i)
       sfc(1)%Snow_temperature = min(tsk(j,i),273.15)
       sfc(1)%Water_temperature = max(tsk(j,i),273.15)
       sfc(1)%Ice_temperature = min(tsk(j,i),273.15)
       sfc(1)%Soil_Moisture_Content = 0.05 !soil_moist(i,j,1) !Get from WRF in g cm-3
       sfc(1)%Vegetation_Fraction = veg_coverage 
       sfc(1)%Soil_Temperature = 283. !soil_temp(i,j,1) !Default
       sfc(1)%Snow_Depth = snowh(j,i)*1000. !Get from WRF in mm
       sfc(1)%LAI = lai(j,i) !From WRF in m2/m2

       !Fill atm arrays - CRTM is in reverse order
       kloop: do k = 1, n_layers
         atm(1)%Level_Pressure(k) = lvl_pressure(n_layers+1-k)
         atm(1)%Pressure(k) = lay_pressure(n_layers+1-k)
         atm(1)%Temperature(k) = temp(n_layers+1-k)
         atm(1)%Absorber(k,1) = max(0.,qv(n_layers+1-k,j,i)*1000.) !H2O g/kg
         atm(1)%Absorber(k,2) = o3
         !if(i==ii.and.j==jj)print*,' '
         !if(i==ii.and.j==jj)print*,'lvl press=',atm(1)%Level_Pressure(k)
         !if(i==ii.and.j==jj)print*,'press=',atm(1)%Pressure(k)
         !if(i==ii.and.j==jj)print*,'q=',atm(1)%Absorber(k,1)
         !if(i==ii.and.j==jj)print*,'q=',qv(n_layers+1-k,j,i)*1000.
         !if(i==ii.and.j==jj)print*,'o3=',atm(1)%Absorber(k,2)
         !if(i==ii.and.j==jj)print*,'o3=',o3
         !if(i==ii.and.j==jj)print*,' '
           
         !Cloud properties - Effective Radius [um] 
         ! and Water Content [kg/m2]
         !Need water content in kg_water/m2 --> need kg_air/m2 factor
         mass_lay = ((atm(1)%Level_Pressure(k)- &
                      atm(1)%Level_Pressure(k-1))*100.)/9.81
         !Microphysics dependent
         if (mp_physics == 99) then  
             atm(1)%Cloud(1)%Water_Content(k) = 0.0
             atm(1)%Cloud(2)%Water_Content(k) = 0.0
             atm(1)%Cloud(1)%Effective_Radius(k) = 0.0
             atm(1)%Cloud(2)%Effective_Radius(k) = 0.0
         else if (mp_physics == 5 .or. mp_physics == 85 .or. &
                  mp_physics == 95) then
             atm(1)%Cloud(1)%Water_Content(k) = 0.0
             atm(1)%Cloud(2)%Water_Content(k) = 0.0
             atm(1)%Cloud(3)%Water_Content(k) = 0.0
             atm(1)%Cloud(4)%Water_Content(k) = 0.0
             atm(1)%Cloud(5)%Water_Content(k) = 0.0
             atm(1)%Cloud(6)%Water_Content(k) = 0.0       
             atm(1)%Cloud(1)%Effective_Radius(k) = 0.0
             atm(1)%Cloud(2)%Effective_Radius(k) = 0.0
             atm(1)%Cloud(3)%Effective_Radius(k) = 0.0
             atm(1)%Cloud(4)%Effective_Radius(k) = 0.0
             atm(1)%Cloud(5)%Effective_Radius(k) = 0.0
             atm(1)%Cloud(6)%Effective_Radius(k) = 0.0
         else if (mp_physics == 8 .or. mp_physics == 28 .or. &
                  mp_physics == 6 .or. mp_physics == 2 .or. mp_physics == 10) then
             atm(1)%Cloud(1)%Water_Content(k) = max(0.,qcloud(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(2)%Water_Content(k) = max(0.,qice(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(3)%Water_Content(k) = max(0.,qrain(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(4)%Water_Content(k) = max(0.,qsnow(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(5)%Water_Content(k) = max(0.,qgraupel(n_layers+1-k,j,i)*mass_lay)
             !Get effective radius of particles
             CALL eff_radius(qcloud(n_layers+1-k,j,i),qice(n_layers+1-k,j,i), &
                             qrain(n_layers+1-k,j,i),qsnow(n_layers+1-k,j,i), &
                             qgraupel(n_layers+1-k,j,i),qv(n_layers+1-k,j,i), &
                             ncloud(n_layers+1-k,j,i),nice(n_layers+1-k,j,i), &
                             nrain(n_layers+1-k,j,i),nsnow(n_layers+1-k,j,i), &
                             ngraupel(n_layers+1-k,j,i),temp(n_layers+1-k), &
                             lay_pressure(n_layers+1-k), &
                             mp_physics,r_eff(:))
             !if (i == ii .and. j == jj) then 
             !    write(6,*), r_eff(1), r_eff(2), r_eff(3), r_eff(4), r_eff(5)
             !    write(6,*), qcloud(n_layers+1-k,j,i)*mass_lay
             !    write(6,*), qice(n_layers+1-k,j,i)*mass_lay
             !    write(6,*), qrain(n_layers+1-k,j,i)*mass_lay
             !    write(6,*), qsnow(n_layers+1-k,j,i)*mass_lay
             !    write(6,*), qgraupel(n_layers+1-k,j,i)*mass_lay
             !endif
             atm(1)%Cloud(1)%Effective_Radius(k) = r_eff(1)
             atm(1)%Cloud(2)%Effective_Radius(k) = r_eff(2)
             atm(1)%Cloud(3)%Effective_Radius(k) = r_eff(3)
             atm(1)%Cloud(4)%Effective_Radius(k) = r_eff(4)
             atm(1)%Cloud(5)%Effective_Radius(k) = r_eff(5)
         endif
         
   
       enddo kloop

       !if (j == jj .and. i == ii) write(6,*) chinfo(1)%Sensor_Channel 
       !Call radiative transfer code
       err_stat = CRTM_Forward( atm         , & ! Input
                                sfc         , & ! Input
                                geo         , & ! Input
                                chinfo(1:1) , & ! Input (n:n) makes it an
                                                ! array which is required
                                rts)!         , & ! Output
                                !Options=opt ) ! Optional input
       if ( err_stat /= SUCCESS ) then
           write(6,*) 'Failed in the forward model'
       endif
         
       !Set the brightness temperature and radiance values from the output
       if (trim(request_var) == 'Radiance') then
         crtm_out(j,i) = rts(1,1)%Radiance
       else if (trim(request_var) == 'Brightness Temperature') then
         crtm_out(j,i) = rts(1,1)%Brightness_Temperature
       else if (trim(request_var) == 'Up Radiance') then
         crtm_out(j,i) = rts(1,1)%Up_Radiance
       else if (trim(request_var) == 'Down Radiance') then
         crtm_out(j,i) = rts(1,1)%Down_Radiance
       else if (trim(request_var) == 'Down Solar Radiance') then
         crtm_out(j,i) = rts(1,1)%Down_Solar_Radiance
       endif
       !rad(i,j,n) = rts(n,1)%Radiance

       ! Deallocate channel-dependent arrays
       deallocate( rts, STAT = alloc_stat )
       if ( alloc_stat /= 0 ) then
           write(6,*) 'Failed to deallocate'
       endif  
    
    else
       write(6,*) "Location out of satellite range :", lat, lon
    endif

   
  enddo jloop
enddo iloop

!Deallocate Arrays
CALL crtm_atmosphere_destroy(atm)
CALL crtm_surface_destroy(sfc)
!CALL crtm_rtsolution_destroy(rts)
err_stat = crtm_destroy(chinfo)
if (err_stat /= success) &
   write(6,*),'ERROR*** crtm_destroy error_status=',err_stat
deallocate(r_eff)

end subroutine crtm

subroutine zenith(lat,lon,satlat,satlon,vzenith)

   !This subroutine calculates the satellite viewing angle
   ! using a spherical approximation
   !Source: T. Soler and D. W. Eisemann (1994), Determination of look angles to 
   !        Geostationary Communication Satellites.

   implicit none

   real, intent(in)  :: lat,lon,satlat,satlon
   real, intent(out) :: vzenith
   real, parameter :: pi = 3.141592654
   real, parameter :: d2r = pi/180.
   real, parameter :: r2d = 180./pi
   real, parameter :: r_earth = 6370. !km
   real, parameter :: r_sat = 42164.  !km - geostationary only
   real :: beta, d
   
   !Angle between earth center and earth location to earth center
   ! and satellite location relative to the earth
   beta = acos(cos((lat-satlat)*d2r)*cos((lon-satlon)*d2r))
   !Distance from earth location to satellite
   d = r_sat * (1. + (r_earth/r_sat)**2 - 2.*(r_earth/r_sat)*cos(beta))**(1./2.)
   !Calculate viewing angle
   vzenith = asin((r_sat/d)*sin(beta))*r2d

end subroutine zenith

subroutine solar_zenith(lat,lon,jday,hour,sun_zenith)

   !This subroutine calculates the solar zenith angle
   ! for a given location at a given time of year
   !Crude calculation to estimate solar time

   implicit none

   integer, intent(in) :: jday
   real, intent(in)  :: lat,lon,hour
   real, intent(out) :: sun_zenith
   real, parameter :: pi = 3.141592654
   real, parameter :: d2r = pi/180.
   real, parameter :: r2d = 180./pi
   integer, parameter :: vernal_equinox = 80. !Spring equinox Julian Day
   real, parameter :: earth_tilt = 23.45
   real :: dec_angle, lambda, solar_time, h

   !Calculate the solar declination angle
   lambda = 360.*(jday-vernal_equinox)/365.25
   dec_angle = sin(earth_tilt*d2r)*sin(lambda*d2r)   
   !Calculate hour angle - 0 at 12 pm local time
   ! 360/24 degrees per hour 
   !Convert from UTC to local solar hour
   solar_time = hour + (lon*(24./360.)) 
   !Hour Angle
   h = 2*pi*(solar_time-12.)/24.
   !Calculate solar zenith angle
   sun_zenith = acos( sin(lat*d2r)*sin(dec_angle) + &
                      cos(lat*d2r)*cos(dec_angle)*cos(h) ) * r2d
   !write(6,*) lambda, dec_angle, solar_time, h, hour, lat*d2r, lon*d2r, sun_zenith 
   !write(6,*) sin(lat*d2r)*sin(dec_angle), cos(lat*d2r)*cos(dec_angle)*cos(h)
   !write(6,*) acos( sin(lat*d2r)*sin(dec_angle) + cos(lat*d2r)*cos(dec_angle)*cos(h) ) 

end subroutine solar_zenith

subroutine eff_radius(qcloud,qice,qrain,qsnow,qgraupel,qv,ncloud,nice,nrain,&
                      nsnow,ngraupel,t,p,mp_physics,re)
   
   !This subroutine calculates the effective
   ! radius of the drop distribution
   ! re = third moment / second moment
   !The calculation is microphysics dependent    

   implicit none

   real, intent(in) :: qcloud, qice, qrain, qsnow, qgraupel, qv
   real, intent(in) :: ncloud, nice, nrain, nsnow, ngraupel
   real, intent(in) :: t, p
   integer, intent(in) :: mp_physics 
   real, intent(inout) :: re(5) !1 = cloud, 2 = ice, 3 = rain, 4 = snow, 5 = graupel

   !Lin microphysics -- Lin et al. (1983), J. Climate and Applied Met.
   !  updates from Rutledge and Hobbs (1984), JAM (Graupel number and density)
   real, parameter :: lin_n0r = 8.e6, lin_n0s = 3.e6, lin_n0g = 4.e6  ![m^-4]
   real, parameter :: lin_rhor = 1000., lin_rhos = 100.   ![kg/m^3] 
   real, parameter :: lin_rhoi = 100., lin_rhog = 400.    ![kg/m^3]
   real, parameter :: lin_ic = 3.e8 ![m^-3] UPP
   real, parameter :: lin_cc = 1.e9 ![m^-3] 

   !Thompson microphysics -- Thompson et al. (2008), MWR.
   !  Some values/equations were adapted from the module_mp_thompson.F code in 
   !    WRF
   real :: wgamma
   real, parameter :: thomp_rhor = 1000., thomp_rhos = 100. ![kg/m^3]
   real, parameter :: thomp_rhog = 500.,  thomp_rhoi = 890. ![kg/m^3]
   real, parameter :: thomp_nt_c = 100.e6 ![m^-3] 
   !Minimum microphysical values from Thompson code
   real, parameter :: R1 = 1.e-12, R2 = 1.e-6 
   !Gamma function ratios from Thompson code
   real, dimension(15), parameter:: gamma_ratio = (/24,60,120,210,336,   &
                   504,720,990,1320,1716,2184,2730,3360,4080,4896/)
   real, parameter :: mu_i = 0.0, mu_r = 0.0, mu_g = 0.0, mu_s = 0.6357
   real, parameter :: bm_r = 3.0, bm_s = 2.0, bm_g = 3.0, bm_i = 3.0
   !Max and min graupel concentration
   real, parameter :: g_min = 1.E4, g_max = 3.e6 
   real :: ci(2), cr(3), crg(3), cig(2), cg(3), cgg(3)
   real :: rc, ri, ni, rr, nr, rg, ng, ng_min
   real :: tc, smob, smob2, smoc, rs, loga, a, b, cs
   integer :: mu_c
   !..For snow moments conversions (from Field et al. 2005)
   real, dimension(10), parameter :: &
      thomp_sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
                    0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
   real, dimension(10), parameter :: &
      thomp_sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
                    0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

   !Morrison microphysics -- Morrison et al. (2009), MWR
   !  Some equations were adapted from the module_mp_morr_two_moment.F code in
   !    WRF
   real, parameter :: morr_rhor = 997., morr_rhoi = 500. ![kg/m^3]
   real, parameter :: morr_rhog = 400., morr_rhos = 100. ![kg/m^3] !hail = 900.
   real, parameter :: morr_qsmall = 1.e-14 !Smallest allowed hydrometeor mixing ratio
   real, parameter :: morr_re_default = 25. ![um] effective radius if q < qsmall
   real, parameter :: morr_nc = 250.e6 ![m^-3] constant from WRF code
   !Slope parameter min and max from Morrison WRF code
   real, parameter :: morr_lamrMax = 1./20.e-6, morr_lamrMin = 1./2800e-6
   real, parameter :: morr_lamiMax = 1./1.e-6, morr_lamiMin = 1./(2.*125.e-6*100.e-6)
   real, parameter :: morr_lamsMax = 1./10.e-6, morr_lamsMin = 1./2000e-6
   real, parameter :: morr_lamgMax = 1./20.e-6, morr_lamgMin = 1./2000e-6
   real :: morr_lamcMin, morr_lamcMax
   real :: morr_mu_c

   !Global (somewhat) varaibles
   real :: nc, lamc, lami, lams, lamr, lamg
   !Minimum mixing ratios needed to compute effective radius
   real, parameter :: min_qc = 1.e-8, min_qr = 1.e-8, min_qi = 1.e-8
   real, parameter :: min_qg = 1.e-8, min_qs = 1.e-8
   real, parameter :: m2um = 1.e6 !meters to micrometers
   real, parameter :: rd = 287.05 !dry air gas constant
   real, parameter :: pi = 3.14159265359
   real :: rho_air
 
   !Start Calculations
   !Initialize 
   re(:) = 0.0

   !Calculate air density from virtual temperature
   rho_air = p/(rd * t * (1. + 0.608 * qv))

   !Thompson
   if (mp_physics == 8) then
         
     !Calcuate cloud effective radius
     rc = MAX(R1, qcloud*rho_air)
     nc = MAX(R2, thomp_nt_c)
     if (rc .le. R1 .or. nc .le. R2) then
       re(1) = 2.49E-6 !Low value from Thompson code
     else
       !Shape parameter
       if (nc .lt. 100) then
           mu_c = 15
       elseif (nc .gt. 1.e10) then
           mu_c = 2
       else
           mu_c = MIN(15, NINT(1.e9/nc) + 2)
       endif
       !Slope parameter
       lamc = (nc*pi*thomp_rhor*gamma_ratio(mu_c)/rc*6.)**(1./3.)
       re(1) = m2um * MAX(2.51E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_c)/lamc), 50.E-6))  
     endif 

     !Calculate ice effective radius
     if (qice .lt. min_qi) then
       re(2) = 0.0
     else
       ci(1) = mu_i + 1.
       ci(2) = bm_i + mu_i + 1.
       cig(1) = wgamma(ci(1))
       cig(2) = wgamma(ci(2))
       ri = MAX(R1, qice*rho_air)
       ni = MAX(R2, nice*rho_air) !Double moment ice - nice is in units [kg^-1]
       if (ri .le. R1 .or. ni .le. R2) then
         re(2) = 4.99E-6 !Low value from Thompson code
       else
         !Slope parameter
         lami = (pi*thomp_rhoi*cig(2)*ni/ri*6*cig(1))**(1./bm_i)
         re(2) = m2um * MAX(5.01E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_i)/lami), 125.E-6)) 
       endif 
     endif

     !Calculate rain effective radius
     if (qrain .lt. min_qr) then
       re(3) = 0.0
     else
       cr(1) = bm_r + 1.
       cr(2) = mu_r + 1.
       cr(3) = bm_r + mu_r + 1.
       crg(1) = wgamma(cr(1))
       crg(2) = wgamma(cr(2))  
       crg(3) = wgamma(cr(3))   
       rr = MAX(R1, qrain*rho_air)
       nr = MAX(R2, nrain*rho_air) !Double moment rain - nrain is in units [kg^-1]
       lamr = (pi*thomp_rhor*crg(3)*nr/rr*crg(2))**(1./bm_r)
       re(3) = m2um * MAX(25.01E-6, MIN(SNGL(0.5D0 * DBLE(3.+mu_r)/lamr), 999.E-6))
     endif

     !Calculate snow effective radius
     cs = bm_s + 1.
     rs = MAX(R1, qsnow*rho_air)
     if (rs .le. R1) then
       re(4) = 9.99E-6 !Low value from Thompson code
     else
       tc = MIN(-0.1, t-273.15)
       smob = rs/0.069
       if (bm_s .gt. (2.0-1.e-3) .and. bm_s .lt. (2.0+1.e-3)) then
         smob2 = smob
       else
         loga = thomp_sa(1) + thomp_sa(2)*tc + thomp_sa(3)*bm_s &
                + thomp_sa(4)*tc*bm_s + thomp_sa(5)*tc*tc &
                + thomp_sa(6)*bm_s*bm_s + thomp_sa(7)*tc*tc*bm_s &
                + thomp_sa(8)*tc*bm_s*bm_s + thomp_sa(9)*tc*tc*tc &
                + thomp_sa(10)*bm_s*bm_s*bm_s
         a = 10.0**loga
         b = thomp_sb(1) + thomp_sb(2)*tc + thomp_sb(3)*bm_s &
                + thomp_sb(4)*tc*bm_s + thomp_sb(5)*tc*tc &
                + thomp_sb(6)*bm_s*bm_s + thomp_sb(7)*tc*tc*bm_s &
                + thomp_sb(8)*tc*bm_s*bm_s + thomp_sb(9)*tc*tc*tc &
                + thomp_sb(10)*bm_s*bm_s*bm_s 
         smob2 = (smob/a)**(1./b)  
       endif
       !..Calculate bm_s+1 (th) moment.  Useful for diameter calcs.
       loga = thomp_sa(1) + thomp_sa(2)*tc + thomp_sa(3)*cs &
              + thomp_sa(4)*tc*cs + thomp_sa(5)*tc*tc &
              + thomp_sa(6)*cs*cs + thomp_sa(7)*tc*tc*cs &
              + thomp_sa(8)*tc*cs*cs + thomp_sa(9)*tc*tc*tc &
              + thomp_sa(10)*cs*cs*cs
       a = 10.0**loga
       b = thomp_sb(1) + thomp_sb(2)*tc + thomp_sb(3)*cs &
              + thomp_sb(4)*tc*cs + thomp_sb(5)*tc*tc &
              + thomp_sb(6)*cs*cs + thomp_sb(7)*tc*tc*cs &
              + thomp_sb(8)*tc*cs*cs + thomp_sb(9)*tc*tc*tc &
              + thomp_sb(10)*cs*cs*cs 
       smoc = a * smob2**b
       re(4) = m2um * MAX(10.E-6, MIN(0.5*(smoc/smob), 999.E-6))
     endif
     
     !Calculate graupel effective radius
     if (qgraupel .lt. min_qg) then
       re(5) = 0.0
     else
       cg(1) = bm_g + 1.
       cg(2) = mu_g + 1.
       cg(3) = bm_g + mu_g + 1.
       cgg(1) = wgamma(cg(1))
       cgg(2) = wgamma(cg(2))
       cgg(3) = wgamma(cg(3)) 
       rg = qgraupel    
       ng = MAX(DBLE(g_min), MIN(200./qgraupel,DBLE(g_max)))
       ng_min = MIN(ng, g_max)
       ng = ng_min
       lamg = ((ng*pi*thomp_rhog*cgg(1)/6.*rg*rho_air)**(1./cg(1))) &
              * (cgg(3)/cgg(2)*cgg(1))**(1./bm_g)
       re(5) = m2um * (3./2.) * (1./lamg)
     endif 

   else if (mp_physics == 10) then !Morrison double moment 

     !Calcuate cloud effective radius
     if (qcloud .lt. morr_qsmall) then
       re(1) = morr_re_default
     else
       nc = (morr_nc/rho_air)
       !Determine the spectral shape parameter for cloud droplets (morr_mu_c)
       !Uses the Martin et al. (1994) formula (adapted from Morrison WRF code)
       morr_mu_c = 0.0005714*(nc/1.e6*rho_air) + 0.2714
       morr_mu_c = (1./morr_mu_c*morr_mu_c) - 1.
       morr_mu_c = MAX(morr_mu_c,2.)
       morr_mu_c = MIN(morr_mu_c,10.)

       lamc = (pi*morr_rhor*wgamma(morr_mu_c+4.)*nc/6.*qcloud*wgamma(morr_mu_c+1.))**(1./3.) 
       !Check if calculated slope parameter is within MIN/MAX values from Morrison
       ! WRF code
       morr_lamcMin = (morr_mu_c+1)/60.e-6
       morr_lamcMax = (morr_mu_c+1)/1.e-6
       if (lamc .lt. morr_lamcMin) lamc = morr_lamcMin
       if (lamc .gt. morr_lamcMax) lamc = morr_lamcMax
       re(1) = m2um * (1./2.) * (wgamma(morr_mu_c+4.)/wgamma(morr_mu_c+3.)) / lamc
     endif

     !Calcuate ice effective radius
     if (qice .lt. morr_qsmall) then
       re(2) = morr_re_default
     else
       lami = (wgamma(4.)*morr_rhoi*pi*nice/6.*qice)**(1./3.)
       !Check if calculated slope parameter is within MIN/MAX values from Morrison
       ! WRF code
       if (lami .lt. morr_lamiMin) lami = morr_lamiMin
       if (lami .gt. morr_lamiMax) lami = morr_lamiMax     
       re(2) = m2um * (3./2.) / lami
     endif

     !Calcuate rain effective radius
     if (qrain .lt. morr_qsmall) then
       re(3) = morr_re_default
     else
       lamr = (morr_rhor*pi*nice/qrain)**(1./3.)
       !Check if calculated slope parameter is within MIN/MAX values from Morrison
       ! WRF code
       if (lamr .lt. morr_lamrMin) lamr = morr_lamrMin
       if (lamr .gt. morr_lamrMax) lamr = morr_lamrMax
       re(3) = m2um * (3./2.) / lamr
     endif

     !Calcuate snow effective radius
     if (qsnow .lt. morr_qsmall) then
       re(4) = morr_re_default
     else
       lams = (wgamma(4.)*morr_rhos*pi*nsnow/6.*qsnow)**(1./3.)
       !Check if calculated slope parameter is within MIN/MAX values from Morrison
       ! WRF code
       if (lams .lt. morr_lamsMin) lams = morr_lamsMin
       if (lams .gt. morr_lamsMax) lams = morr_lamsMax
       re(4) = m2um * (3./2.) / lams
     endif

     !Calcuate graupel effective radius
     if (qgraupel .lt. morr_qsmall) then
       re(5) = morr_re_default
     else
       lamg = (wgamma(4.)*morr_rhog*pi*ngraupel/6.*qgraupel)**(1./3.)
       !Check if calculated slope parameter is within MIN/MAX values from Morrison
       ! WRF code
       if (lamg .lt. morr_lamgMin) lamg = morr_lamgMin
       if (lamg .gt. morr_lamgMax) lamg = morr_lamgMax
       re(5) = m2um * (3./2.) / lamg
     endif

   else !From Lin but default if mp scheme has not yet been implemented

     !Calcuate cloud effective radius - assume constant number concentration
     re(1) = m2um * (3./2.) * ((rho_air * qcloud) / (pi * lin_rhor * lin_cc)) ** (1./3.)

     !Calcuate ice effective radius - assume constant number concentration
     re(2) = m2um * (3./2.) * ((rho_air * qice) / (pi * lin_rhoi * lin_ic)) ** (1./3.) 

     !Calcuate rain effective radius - gamma distribution
     re(3) = m2um * (3./2.) * ((rho_air * qrain) / (pi * lin_rhor * lin_n0r)) ** (1./4.)

     !Calcuate snow effective radius - gamma distribution
     re(4) = m2um * (3./2.) * ((rho_air * qsnow) / (pi * lin_rhos * lin_n0s)) ** (1./4.)

     !Calcuate graupel effective radius - gamma distribution
     re(5) = m2um * (3./2.) * ((rho_air * qgraupel) / (pi * lin_rhog * lin_n0g)) ** (1./4.)

   endif

end subroutine eff_radius

!GAMMA function = EXP(ln(GAMMA))
real function wgamma(x)

implicit none
real :: gammln
real, intent(in) :: x

wgamma = EXP(gammln(x))

end function wgamma

!Numerical Recipes (Fortran Version), Press et al. (1989), pg. 157
!Returns the value ln(gamma(xx)) for xx > 0
real function gammln(xx)

implicit none
real, intent(in) :: xx
real*8, parameter :: stp = 2.50662827465
real*8, dimension(6), parameter :: cof = (/76.18009173,&
                       -86.50532033, 24.01409822,&
                       -1.231739516, .120858003e-2,&
                       -.536382e-5/)

real*8 :: ser, tmp, x, y
integer :: j

x = xx
y = x
tmp = x + 5.5D0
tmp = (x + 0.5D0) * LOG(tmp) - tmp
ser = 1.000000000D0
do j = 1, 6
  y = y + 1.D0
  ser = ser + cof(j)/y
enddo
gammln = tmp + LOG(stp * ser/x)

end function gammln


