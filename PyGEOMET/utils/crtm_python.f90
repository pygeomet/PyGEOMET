SUBROUTINE CRTM(press, theta, qv, qcloud, qice, qrain, qsnow, qgraupel, &
                qhail, lai, u10, v10, seaice, snowh, coszen, &
                vegfrac, ptop, tsk, ivegtyp, xland, landuse, mp_physics, &
                lat, lon, sensor, channel, coeff_path, request_var, crtm_out, ii, jj, kk) 

use crtm_module

implicit none

!Define Variables
!Input Variables
integer, intent(in) :: ii, jj, kk
real, intent(in) :: press(kk,jj,ii), theta(kk,jj,ii), qv(kk,jj,ii)
real, intent(in) :: qcloud(kk,jj,ii), qice(kk,jj,ii), qrain(kk,jj,ii)
real, intent(in) :: qsnow(kk,jj,ii), qgraupel(kk,jj,ii), qhail(kk,jj,ii)
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
         mp_physics == 2 .or. mp_physics == 28)then
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
         mp_physics == 6 .or. mp_physics == 2) then
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
                  mp_physics == 6 .or. mp_physics == 2) then
             atm(1)%Cloud(1)%Water_Content(k) = max(0.,qcloud(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(2)%Water_Content(k) = max(0.,qice(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(3)%Water_Content(k) = max(0.,qrain(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(4)%Water_Content(k) = max(0.,qsnow(n_layers+1-k,j,i)*mass_lay)
             atm(1)%Cloud(5)%Water_Content(k) = max(0.,qgraupel(n_layers+1-k,j,i)*mass_lay)
             !Get effective radius of particles
             CALL eff_radius(qcloud(n_layers+1-k,j,i),qice(n_layers+1-k,j,i), &
                             qrain(n_layers+1-k,j,i),qsnow(n_layers+1-k,j,i), &
                             qgraupel(n_layers+1-k,j,i),qv(n_layers+1-k,j,i), &
                             temp(n_layers+1-k), &
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

subroutine eff_radius(qcloud,qice,qrain,qsnow,qgraupel,qv,t,p,mp_physics,re)
   
   !This subroutine calculates the effective
   ! radius of the drop distribution
   ! re = third moment / second moment
   !The calculation is microphysics dependent    

   implicit none

   real, intent(in) :: qcloud, qice, qrain, qsnow, qgraupel, qv
   real, intent(in) :: t, p
   integer, intent(in) :: mp_physics 
   real, intent(inout) :: re(5) !1 = cloud, 2 = ice, 3 = rain, 4 = snow, 5 = graupel

   !Lin microphysics -- Lin et al. (1983), J. Climate and Applied Met.
   !  updates from Rutledge and Hobbs (1984), JAM (Graupel number and density)
   real, parameter :: lin_n0r = 8.e6, lin_n0s = 3.e6, lin_n0g = 4.e6  ![m^-4]
   real, parameter :: lin_rhor = 1000., lin_rhos = 100.   ![kg/m^3] 
   real, parameter :: lin_rhoi = 100., lin_rhog = 400.    ![kg/m^3]
   real, parameter :: lin_ric = 3.e8 ![m^-3] from UPP? 

   !Global varaibles
   real :: nc, ni, n0r, n0s, n0g, rhog, rhos, rhor, rhoi
   real, parameter :: m2um = 1.e6 !meters to micrometers
   real, parameter :: rd = 287.05 !dry air gas constant
   real, parameter :: pi = 3.14159
   real :: rho_air
 
   !Start Calculations
   !Initialize 
   re(:) = 0.0
 
   !Transfer variables -- changing physics
   n0r = lin_n0r
   n0s = lin_n0s
   n0g = lin_n0g
   nc = lin_ric
   ni = lin_ric
   rhor = lin_rhor
   rhos = lin_rhos
   rhog = lin_rhog
   rhoi = lin_rhoi

   !Calculate air density from virtual temperature
   rho_air = p/(rd * t * (1. + 0.608 * qv)) 

   !Cacluate cloud effective radius - assume constant number concentration
   re(1) = m2um * (3./2.) * ((rho_air * qcloud) / (pi * rhor * nc)) ** (1./3.)

   !Cacluate ice effective radius - assume constant number concentration
   re(2) = m2um * (3./2.) * ((rho_air * qice) / (pi * rhoi * ni)) ** (1./3.) 

   !Cacluate rain effective radius - gamma distribution
   re(3) = m2um * (3./2.) * ((rho_air * qrain) / (pi * rhor * n0r)) ** (1./4.)

   !Cacluate snow effective radius - gamma distribution
   re(4) = m2um * (3./2.) * ((rho_air * qsnow) / (pi * rhos * n0s)) ** (1./4.)

   !Cacluate graupel effective radius - gamma distribution
   re(5) = m2um * (3./2.) * ((rho_air * qgraupel) / (pi * rhog * n0g)) ** (1./4.)

end subroutine eff_radius

