!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute AIMACI_mosaic subroutine for chem(:,:,:,:). 
! * Author: Zihan Xia
! * Email: xiazh@mail.ustc.edu.cn
! Last update: zhxia 10/2024
!
! Contact:
! Chun Zhao, PhD
! Professor 
! University of Science and Technology of China 
! Hefei, Anhui, 230026 
! Email: ChunZhao@ustc.edu.cn
!
!
! Additional information:
! aemol.ustc.edu.cn: resources
!
!************************************************************************
 module module_AIMACI_mosaic 
    use iamas_kind_types
    use iamas_atmchem_vars

 interface
 subroutine AIMACI_mosaic(c_cin, c_cout, batch_size) bind(c, name='AIMACI_mosaic')
    use iso_c_binding
    implicit none
    real(c_double), intent(in) :: c_cin(67,*)
    real(c_double), intent(in) :: c_cout(61,*)
    integer(c_int), value :: batch_size
 end subroutine AIMACI_mosaic
 end interface

 contains
 subroutine AIMACI_mosaic_driver(chem,                                   &
                         rho_phy,t_phy,p_phy,qv,                         &
                         its,ite,kts,kte)
 IMPLICIT NONE

 INTEGER, INTENT(IN) :: its,ite,kts,kte
 REAL, DIMENSION(1:num_chem, kts:kte, its:ite), INTENT(INOUT) :: chem
 REAL, DIMENSION(kts:kte, its:ite), INTENT(IN) :: rho_phy, t_phy, p_phy, qv

!local vars
 real(kind=8), target, allocatable :: CIN(:,:), COUT(:,:)
 real(kind=8), pointer :: c_cin(:,:), c_cout(:,:)
 integer :: batch_size

 integer:: i,k,h

 batch_size = (ite-its+1)*(kte-kts+1) 
 allocate(CIN(67,batch_size))
 allocate(COUT(61,batch_size))
 c_cin => CIN
 c_cout => COUT 

#ifdef USE_DEBUG

 do i = 1,batch_size

 cin(:,i) = (/0.9130055278985415, 0.98679555078927217,                                 &
             0.89512404628847231, 0.53955349910009487, 0.013564209646975574,           &
             -0.0061151145760432483, -6.17416698176994e-09, 4.3720252415208419e-17,    &
             7.2348435940121578e-13, 0.00047873826829603701, 0.00054082790573697036,   &
             7.7108448428227266e-17, 3.6732327187566789e-17, 9.8019542296001613e-05,   &
             -1.3928814487380533e-17, 0.015689379979161828, 1.704639080652799e-06,     &
             7.9461290717567299e-15, 1.7338698166621759e-06, 2.1588122897846039e-07,   &
             1.3415363597879342e-20, 5.1765038100417284e-07, -3.04608222414725e-15,    &
             7.84434485227309e-15, 3.3360000357640226e-15, 5.2912639747924341e-17,     &
             2.0434839405992614e-17, 0.0016552785086368402, 6.1158489466888711e-18,    &
             0.020761958007549972, 5.211399609130818e-07, 2.5750617419098483e-17,      &
             4.7185091217411965e-07, 7.493447893788471e-08, 1.174966130178965e-20,     &
             1.4727622328547549e-05, 1.7706048623090447e-17, 2.3989759349299358e-17,   &
             2.7088062122971482e-15, 2.9660890883409095e-16, 1.0025705976653051e-16,   &
             0.011179486958879171, 3.6647531163620201e-17, 0.02645980678061104,        &
             3.9012261745329891e-07, 1.1897798919670682e-18, 3.5033174917783063e-07,   &
             5.5964006329074814e-08, 5.6640827551761724e-20, 0.00023080447337814902,   &
             7.8706100753508257e-19, 1.0919223723240609e-18, 1.0232663336518518e-14,   &
             1.3931449011783143e-15, 4.4246606666098907e-17, 0.01208783958715392,      &
             8.9209076323768915e-17, 0.028755775402165999, 1.275010394342248e-09,      &
             1.5424520544298527e-19, 4.206189202406803e-07, 6.6919651382817994e-08,    &
             1.199388019396806e-20, 0.00056521483612312091, 9.8017766340973698e-20,    &
             1.3545784126528938e-19, 1.522618596265804e-15/)

 enddo  

#else

    CIN(1,:)= (reshape(t_phy(kts:kte,its:ite),[batch_size]) - 177.03013610839844 ) / (317.97467041015625 - 177.03013610839844)
    CIN(2,:)= (reshape(p_phy(kts:kte,its:ite),[batch_size]) - 5530.49462890625 ) / (105555.1796875 - 5530.49462890625)
    CIN(3,:)= (reshape(rho_phy(kts:kte,its:ite),[batch_size]) - 0.0875762552022934 ) / (1.556276798248291 - 0.0875762552022934)
    CIN(4,:)= (reshape(qv(kts:kte,its:ite),[batch_size]) - 0.0 ) / (0.027006208896636963 - 0.0)
    CIN(5,:)= (reshape(chem(p_h2so4,kts:kte,its:ite),[batch_size]) - 9.952121592105012e-12 ) / (0.0015894059324637055 - 9.952121592105012e-12)
    CIN(6,:)= (reshape(chem(p_hno3,kts:kte,its:ite),[batch_size]) - 4.493148964002103e-09 ) / (0.051831234246492386 - 4.493148964002103e-09)
    CIN(7,:)= (reshape(chem(p_hcl,kts:kte,its:ite),[batch_size]) - 9.50994528652292e-17 ) / (0.0039164284244179726 - 9.50994528652292e-17)
    CIN(8,:)= (reshape(chem(p_nh3,kts:kte,its:ite),[batch_size]) - 9.661394115747498e-17 ) / (0.20336245000362396 - 9.661394115747498e-17)
    CIN(9,:)= (reshape(chem(p_msa,kts:kte,its:ite),[batch_size]) - 9.50994528652292e-17 ) / (8.695050928508863e-05 - 9.50994528652292e-17)
    CIN(10,:)= (reshape(chem(p_so2,kts:kte,its:ite),[batch_size]) - 3.534977679464646e-07 ) / (0.15556906163692474 - 3.534977679464646e-07)
    CIN(11,:)= (reshape(chem(p_no2,kts:kte,its:ite),[batch_size]) - 9.251951382793777e-08 ) / (0.11995994299650192 - 9.251951382793777e-08)
    CIN(12,:)= (reshape(chem(p_so4_a01,kts:kte,its:ite),[batch_size]) - 6.241787029087858e-19 ) / (5.2291460037231445 - 6.241787029087858e-19)
    CIN(13,:)= (reshape(chem(p_no3_a01,kts:kte,its:ite),[batch_size]) - 6.241787029087858e-19 ) / (16.214080810546875 - 6.241787029087858e-19)
    CIN(14,:)= ((reshape(chem(p_cl_a01,kts:kte,its:ite),[batch_size])) - 6.241787029087858e-19 ) / (0.6946127414703369 - 6.241787029087858e-19) !p_cl_a01
    CIN(15,:)= (reshape(chem(p_nh4_a01,kts:kte,its:ite),[batch_size]) - 6.241787029087858e-19 ) / (5.178518772125244 - 6.241787029087858e-19)
    CIN(16,:)= ((reshape(chem(p_na_a01,kts:kte,its:ite),[batch_size])) - 6.241787029087858e-19 ) / (0.04817196726799011 - 6.241787029087858e-19) !p_na_a01
    CIN(17,:)= (reshape(chem(p_oin_a01,kts:kte,its:ite),[batch_size]) - 7.895500985364591e-19 ) / (9.459769248962402 - 7.895500985364591e-19)
    CIN(18,:)= ((reshape(chem(p_dust_a01,kts:kte,its:ite),[batch_size])) - 6.241787029087858e-19 ) / (0.14824660122394562 - 6.241787029087858e-19) !p_dust_a01
    CIN(19,:)= (reshape(chem(p_oc_a01,kts:kte,its:ite),[batch_size]) - 8.599956946336188e-19 ) / (23.83402442932129 - 8.599956946336188e-19)
    CIN(20,:)= (reshape(chem(p_bc_a01,kts:kte,its:ite),[batch_size]) - 8.599956946336188e-19 ) / (3.659121036529541 - 8.599956946336188e-19)
    CIN(21,:)= (reshape(chem(p_water_a01,kts:kte,its:ite),[batch_size]) - 9.714934538947851e-17 ) / (167.65695190429688 - 9.714934538947851e-17)
    CIN(22,:)= (reshape(chem(p_num_a01,kts:kte,its:ite),[batch_size]) - 1.9796543023317188e-16 ) / (110981087232.0 - 1.9796543023317188e-16)
    CIN(23,:)= ((reshape(chem(p_ca_a01,kts:kte,its:ite),[batch_size])) - 6.241787029087858e-19 ) / (0.0005989763885736465 - 6.241787029087858e-19) !p_ca_a01
    CIN(24,:)= ((reshape(chem(p_co3_a01,kts:kte,its:ite),[batch_size])) - 6.241787029087858e-19 ) / (0.0005442982073873281 - 6.241787029087858e-19) !p_co3_a01
    CIN(25,:)= (reshape(chem(p_msa_a01,kts:kte,its:ite),[batch_size]) - 6.241787029087858e-19 ) / (0.24560633301734924 - 6.241787029087858e-19)
    CIN(26,:)= (reshape(chem(p_so4_a02,kts:kte,its:ite),[batch_size]) - 6.241692420305297e-19 ) / (13.656200408935547 - 6.241692420305297e-19)
    CIN(27,:)= (reshape(chem(p_no3_a02,kts:kte,its:ite),[batch_size]) - 6.241692420305297e-19 ) / (56.314762115478516 - 6.241692420305297e-19)
    CIN(28,:)= (reshape(chem(p_cl_a02,kts:kte,its:ite),[batch_size]) - 6.241692420305297e-19 ) / (2.8376691341400146 - 6.241692420305297e-19) !p_cl_a02
    CIN(29,:)= (reshape(chem(p_nh4_a02,kts:kte,its:ite),[batch_size]) - 5.901367051846254e-19 ) / (18.279104232788086 - 5.901367051846254e-19)
    CIN(30,:)= (reshape(chem(p_na_a02,kts:kte,its:ite),[batch_size]) - 6.241692420305297e-19 ) / (1.8478686809539795 - 6.241692420305297e-19) !p_na_a02
    CIN(31,:)= (reshape(chem(p_oin_a02,kts:kte,its:ite),[batch_size]) - 8.599450815198882e-19 ) / (30.915788650512695 - 8.599450815198882e-19)
    CIN(32,:)= (reshape(chem(p_dust_a02,kts:kte,its:ite),[batch_size]) - 4.670115430903549e-19 ) / (425.0188903808594 - 4.670115430903549e-19) !p_dust_a02
    CIN(33,:)= (reshape(chem(p_oc_a02,kts:kte,its:ite),[batch_size]) - 8.599450815198882e-19 ) / (122.57513427734375 - 8.599450815198882e-19)
    CIN(34,:)= (reshape(chem(p_bc_a02,kts:kte,its:ite),[batch_size]) - 8.599450815198882e-19 ) / (16.664186477661133 - 8.599450815198882e-19)
    CIN(35,:)= (reshape(chem(p_water_a02,kts:kte,its:ite),[batch_size]) - 9.752534860872061e-17 ) / (393.0467834472656 - 9.752534860872061e-17)
    CIN(36,:)= (reshape(chem(p_num_a02,kts:kte,its:ite),[batch_size]) - 1.9505069721744122e-16 ) / (10640852992.0 - 1.9505069721744122e-16)
    CIN(37,:)= (reshape(chem(p_ca_a02,kts:kte,its:ite),[batch_size]) - 1.6839498633579893e-19 ) / (1.7172482013702393 - 1.6839498633579893e-19) !p_ca_a02
    CIN(38,:)= (reshape(chem(p_co3_a02,kts:kte,its:ite),[batch_size]) - 7.416211212701199e-20 ) / (2.552623748779297 - 7.416211212701199e-20) !p_co3_a02
    CIN(39,:)= (reshape(chem(p_msa_a02,kts:kte,its:ite),[batch_size]) - 6.241692420305297e-19 ) / (0.05827217549085617 - 6.241692420305297e-19)
    CIN(40,:)= (reshape(chem(p_so4_a03,kts:kte,its:ite),[batch_size]) - 6.239740274059672e-19 ) / (6.48736047744751 - 6.239740274059672e-19)
    CIN(41,:)= (reshape(chem(p_no3_a03,kts:kte,its:ite),[batch_size]) - 6.239740274059672e-19 ) / (5.950143814086914 - 6.239740274059672e-19)
    CIN(42,:)= (reshape(chem(p_cl_a03,kts:kte,its:ite),[batch_size]) - 6.239740274059672e-19 ) / (24.391042709350586 - 6.239740274059672e-19) !p_cl_a03
    CIN(43,:)= (reshape(chem(p_nh4_a03,kts:kte,its:ite),[batch_size]) - 5.951132305448975e-19 ) / (1.4401476383209229 - 5.951132305448975e-19)
    CIN(44,:)= (reshape(chem(p_na_a03,kts:kte,its:ite),[batch_size]) - 6.239740274059672e-19 ) / (15.824102401733398 - 6.239740274059672e-19) !p_na_a03
    CIN(45,:)= (reshape(chem(p_oin_a03,kts:kte,its:ite),[batch_size]) - 7.897418493422065e-19 ) / (6.093714237213135 - 7.897418493422065e-19)
    CIN(46,:)= (reshape(chem(p_dust_a03,kts:kte,its:ite),[batch_size]) - 4.283751257501335e-19 ) / (5207.19482421875 - 4.283751257501335e-19) !p_dust_a03
    CIN(47,:)= (reshape(chem(p_oc_a03,kts:kte,its:ite),[batch_size]) - 8.59332605975081e-19 ) / (25.758081436157227 - 8.59332605975081e-19)
    CIN(48,:)= (reshape(chem(p_bc_a03,kts:kte,its:ite),[batch_size]) - 8.59332605975081e-19 ) / (3.4673867225646973 - 8.59332605975081e-19)
    CIN(49,:)= (reshape(chem(p_water_a03,kts:kte,its:ite),[batch_size]) - 9.260651578745158e-17 ) / (891.1803588867188 - 9.260651578745158e-17)
    CIN(50,:)= (reshape(chem(p_num_a03,kts:kte,its:ite),[batch_size]) - 4.8548950686411305e-17 ) / (2032314496.0 - 4.8548950686411305e-17)
    CIN(51,:)= (reshape(chem(p_ca_a03,kts:kte,its:ite),[batch_size]) - 4.283751257501335e-19 ) / (21.039169311523438 - 4.283751257501335e-19) !p_ca_a03
    CIN(52,:)= (reshape(chem(p_co3_a03,kts:kte,its:ite),[batch_size]) - 2.458277382931029e-19 ) / (31.504108428955078 - 2.458277382931029e-19) !p_co3_a03
    CIN(53,:)= (reshape(chem(p_msa_a03,kts:kte,its:ite),[batch_size]) - 6.239740274059672e-19 ) / (0.05864517018198967 - 6.239740274059672e-19)
    CIN(54,:)= (reshape(chem(p_so4_a04,kts:kte,its:ite),[batch_size]) - 7.130845921338681e-19 ) / (10.924201965332031 - 7.130845921338681e-19)
    CIN(55,:)= (reshape(chem(p_no3_a04,kts:kte,its:ite),[batch_size]) - 7.130845921338681e-19 ) / (9.010531425476074 - 7.130845921338681e-19)
    CIN(56,:)= (reshape(chem(p_cl_a04,kts:kte,its:ite),[batch_size]) - 7.931677212466714e-19 ) / (151.50611877441406 - 7.931677212466714e-19) !p_cl_a04
    CIN(57,:)= (reshape(chem(p_nh4_a04,kts:kte,its:ite),[batch_size]) - 7.130845921338681e-19 ) / (1.491528868675232 - 7.130845921338681e-19)
    CIN(58,:)= (reshape(chem(p_na_a04,kts:kte,its:ite),[batch_size]) - 7.931677212466714e-19 ) / (98.25638580322266 - 7.931677212466714e-19) !p_na_a04
    CIN(59,:)= (reshape(chem(p_oin_a04,kts:kte,its:ite),[batch_size]) - 7.931677212466714e-19 ) / (33.4320068359375 - 7.931677212466714e-19)
    CIN(60,:)= (reshape(chem(p_dust_a04,kts:kte,its:ite),[batch_size]) - 4.820276629427597e-19 ) / (31028.361328125 - 4.820276629427597e-19) !p_dust_a04
    CIN(61,:)= (reshape(chem(p_oc_a04,kts:kte,its:ite),[batch_size]) - 7.130845921338681e-19 ) / (0.3992343842983246 - 7.130845921338681e-19)
    CIN(62,:)= (reshape(chem(p_bc_a04,kts:kte,its:ite),[batch_size]) - 7.130845921338681e-19 ) / (0.16755656898021698 - 7.130845921338681e-19)
    CIN(63,:)= (reshape(chem(p_water_a04,kts:kte,its:ite),[batch_size]) - 7.753711805897284e-17 ) / (5999.90771484375 - 7.753711805897284e-17)
    CIN(64,:)= (reshape(chem(p_num_a04,kts:kte,its:ite),[batch_size]) - 1.0074338145811425e-18 ) / (176208208.0 - 1.0074338145811425e-18)
    CIN(65,:)= (reshape(chem(p_ca_a04,kts:kte,its:ite),[batch_size]) - 4.820276629427597e-19 ) / (125.36711120605469 - 4.820276629427597e-19) !p_ca_a04
    CIN(66,:)= (reshape(chem(p_co3_a04,kts:kte,its:ite),[batch_size]) - 4.820276629427597e-19 ) / (188.03033447265625 - 4.820276629427597e-19) !p_co3_a04
    CIN(67,:)= (reshape(chem(p_msa_a04,kts:kte,its:ite),[batch_size]) - 7.130845921338681e-19 ) / (6.535012245178223 - 7.130845921338681e-19)

#endif

    call AIMACI_mosaic(c_cin, c_cout, batch_size)
    !where (COUT<0.0_8)
    !       COUT = 0.0_8
    !end where

    chem(p_h2so4,kts:kte,its:ite)= reshape(COUT(1,:),[(kte-kts+1),(ite-its+1)]) * (0.00011774363520089537 - 0.0) + 0.0
    chem(p_hno3,kts:kte,its:ite)= reshape(COUT(2,:),[(kte-kts+1),(ite-its+1)]) * (0.03544870764017105 - 4.0900768015494915e-12) + 4.0900768015494915e-12
    chem(p_hcl,kts:kte,its:ite)= reshape(COUT(3,:),[(kte-kts+1),(ite-its+1)]) * (0.004823804367333651 - 0.0) + 0.0
    chem(p_nh3,kts:kte,its:ite)= reshape(COUT(4,:),[(kte-kts+1),(ite-its+1)]) * (0.20322853326797485 - 0.0) + 0.0
    chem(p_msa,kts:kte,its:ite)= reshape(COUT(5,:),[(kte-kts+1),(ite-its+1)]) * (8.109405462164432e-05 - 0.0) + 0.0
    chem(p_so4_a01,kts:kte,its:ite)= reshape(COUT(6,:),[(kte-kts+1),(ite-its+1)]) * (6.801321029663086 - 0.0) + 0.0
    chem(p_no3_a01,kts:kte,its:ite)= reshape(COUT(7,:),[(kte-kts+1),(ite-its+1)]) * (19.573030471801758 - 0.0) + 0.0
    chem(p_cl_a01,kts:kte,its:ite)= reshape(COUT(8,:),[(kte-kts+1),(ite-its+1)]) * (0.8861722946166992 - 0.0) + 0.0
    chem(p_nh4_a01,kts:kte,its:ite)= reshape(COUT(9,:),[(kte-kts+1),(ite-its+1)]) * (6.290355682373047 - 0.0) + 0.0
    chem(p_na_a01,kts:kte,its:ite)= reshape(COUT(10,:),[(kte-kts+1),(ite-its+1)]) * (0.04805875942111015 - 0.0) + 0.0
    !chem(p_seas_a01,kts:kte,its:ite)= (reshape(COUT(8,:),[(kte-kts+1),(ite-its+1)]) * ((0.8861722946166992 - 0.0) + 0.0) + (reshape(COUT(10,:),[(kte-kts+1),(ite-its+1)]) * 0.04805875942111015 - 0.0) + 0.0)
    chem(p_oin_a01,kts:kte,its:ite)= reshape(COUT(11,:),[(kte-kts+1),(ite-its+1)]) * (8.074164390563965 - 0.0) + 0.0
    chem(p_dust_a01,kts:kte,its:ite)= reshape(COUT(12,:),[(kte-kts+1),(ite-its+1)]) * (0.2568328380584717 - 0.0) + 0.0
    chem(p_oc_a01,kts:kte,its:ite)= reshape(COUT(13,:),[(kte-kts+1),(ite-its+1)]) * (16.133548736572266 - 0.0) + 0.0
    chem(p_bc_a01,kts:kte,its:ite)= reshape(COUT(14,:),[(kte-kts+1),(ite-its+1)]) * (2.5001580715179443 - 0.0) + 0.0
    chem(p_water_a01,kts:kte,its:ite)= reshape(COUT(15,:),[(kte-kts+1),(ite-its+1)]) * (252.06390380859375 - 0.0) + 0.0
    chem(p_num_a01,kts:kte,its:ite)= reshape(COUT(16,:),[(kte-kts+1),(ite-its+1)]) * (75314577408.0 - 0.0) + 0.0
    chem(p_ca_a01,kts:kte,its:ite)= reshape(COUT(17,:),[(kte-kts+1),(ite-its+1)]) * (0.001037708600051701 - 0.0) + 0.0 
    chem(p_co3_a01,kts:kte,its:ite)= reshape(COUT(18,:),[(kte-kts+1),(ite-its+1)]) * (9.026831685332581e-05 - 0.0) + 0.0
    chem(p_msa_a01,kts:kte,its:ite)= reshape(COUT(19,:),[(kte-kts+1),(ite-its+1)]) * (0.25581640005111694 - 0.0) + 0.0
    chem(p_so4_a02,kts:kte,its:ite)= reshape(COUT(20,:),[(kte-kts+1),(ite-its+1)]) * (15.489835739135742 - 0.0) + 0.0
    chem(p_no3_a02,kts:kte,its:ite)= reshape(COUT(21,:),[(kte-kts+1),(ite-its+1)]) * (83.33792114257812 - 0.0) + 0.0
    chem(p_cl_a02,kts:kte,its:ite)= reshape(COUT(22,:),[(kte-kts+1),(ite-its+1)]) * (3.357278823852539 - 0.0) + 0.0
    chem(p_nh4_a02,kts:kte,its:ite)= reshape(COUT(23,:),[(kte-kts+1),(ite-its+1)]) * (26.826784133911133 - 0.0) + 0.0
    chem(p_na_a02,kts:kte,its:ite)= reshape(COUT(24,:),[(kte-kts+1),(ite-its+1)]) * (1.8474076986312866 - 0.0) + 0.0
    !chem(p_seas_a02,kts:kte,its:ite)= (reshape(COUT(22,:),[(kte-kts+1),(ite-its+1)]) * (3.357278823852539 - 0.0) + 0.0) + (reshape(COUT(24,:),[(kte-kts+1),(ite-its+1)]) * (1.8474076986312866 - 0.0) + 0.0)
    chem(p_oin_a02,kts:kte,its:ite)= reshape(COUT(25,:),[(kte-kts+1),(ite-its+1)]) * (31.728206634521484 - 0.0) + 0.0
    chem(p_dust_a02,kts:kte,its:ite)= reshape(COUT(26,:),[(kte-kts+1),(ite-its+1)]) * (411.77838134765625 - 0.0) + 0.0
    chem(p_oc_a02,kts:kte,its:ite)= reshape(COUT(27,:),[(kte-kts+1),(ite-its+1)]) * (125.0230712890625 - 0.0) + 0.0
    chem(p_bc_a02,kts:kte,its:ite)= reshape(COUT(28,:),[(kte-kts+1),(ite-its+1)]) * (17.196720123291016 - 0.0) + 0.0
    chem(p_water_a02,kts:kte,its:ite)= reshape(COUT(29,:),[(kte-kts+1),(ite-its+1)]) * (523.324951171875 - 0.0) + 0.0
    chem(p_num_a02,kts:kte,its:ite)= reshape(COUT(30,:),[(kte-kts+1),(ite-its+1)]) * (11442270208.0 - 0.0) + 0.0
    chem(p_ca_a02,kts:kte,its:ite)= reshape(COUT(31,:),[(kte-kts+1),(ite-its+1)]) * (1.663751244544983 - 0.0) + 0.0
    chem(p_co3_a02,kts:kte,its:ite)= reshape(COUT(32,:),[(kte-kts+1),(ite-its+1)]) * (2.458937406539917 - 0.0) + 0.0
    chem(p_msa_a02,kts:kte,its:ite)= reshape(COUT(33,:),[(kte-kts+1),(ite-its+1)]) * (0.059654369950294495 - 0.0) + 0.0
    chem(p_so4_a03,kts:kte,its:ite)= reshape(COUT(34,:),[(kte-kts+1),(ite-its+1)]) * (6.807916164398193 - 0.0) + 0.0
    chem(p_no3_a03,kts:kte,its:ite)= reshape(COUT(35,:),[(kte-kts+1),(ite-its+1)]) * (7.617635250091553 - 0.0) + 0.0
    chem(p_cl_a03,kts:kte,its:ite)= reshape(COUT(36,:),[(kte-kts+1),(ite-its+1)]) * (24.390901565551758 - 0.0) + 0.0
    chem(p_nh4_a03,kts:kte,its:ite)= reshape(COUT(37,:),[(kte-kts+1),(ite-its+1)]) * (2.319885730743408 - 0.0) + 0.0
    chem(p_na_a03,kts:kte,its:ite)= reshape(COUT(38,:),[(kte-kts+1),(ite-its+1)]) * (15.824056625366211 - 0.0) + 0.0
    !chem(p_seas_a03,kts:kte,its:ite)= (reshape(COUT(36,:),[(kte-kts+1),(ite-its+1)]) * (24.390901565551758 - 0.0) + 0.0) + (reshape(COUT(38,:),[(kte-kts+1),(ite-its+1)]) * (15.824056625366211 - 0.0) + 0.0)
    chem(p_oin_a03,kts:kte,its:ite)= reshape(COUT(39,:),[(kte-kts+1),(ite-its+1)]) * (7.16127347946167 - 0.0) + 0.0
    chem(p_dust_a03,kts:kte,its:ite)= reshape(COUT(40,:),[(kte-kts+1),(ite-its+1)]) * (5197.6474609375 - 0.0) + 0.0
    chem(p_oc_a03,kts:kte,its:ite)= reshape(COUT(41,:),[(kte-kts+1),(ite-its+1)]) * (30.05556869506836 - 0.0) + 0.0
    chem(p_bc_a03,kts:kte,its:ite)= reshape(COUT(42,:),[(kte-kts+1),(ite-its+1)]) * (4.055928707122803 - 0.0) + 0.0
    chem(p_water_a03,kts:kte,its:ite)= reshape(COUT(43,:),[(kte-kts+1),(ite-its+1)]) * (1060.1162109375 - 0.0) + 0.0
    chem(p_num_a03,kts:kte,its:ite)= reshape(COUT(44,:),[(kte-kts+1),(ite-its+1)]) * (2028592640.0 - 0.0) + 0.0
    chem(p_ca_a03,kts:kte,its:ite)= reshape(COUT(45,:),[(kte-kts+1),(ite-its+1)]) * (21.00059700012207 - 0.0) + 0.0
    chem(p_co3_a03,kts:kte,its:ite)= reshape(COUT(46,:),[(kte-kts+1),(ite-its+1)]) * (31.383079528808594 - 0.0) + 0.0
    chem(p_msa_a03,kts:kte,its:ite)= reshape(COUT(47,:),[(kte-kts+1),(ite-its+1)]) * (0.06505321711301804 - 0.0) + 0.0
    chem(p_so4_a04,kts:kte,its:ite)= reshape(COUT(48,:),[(kte-kts+1),(ite-its+1)]) * (10.924202919006348 - 0.0) + 0.0
    chem(p_no3_a04,kts:kte,its:ite)= reshape(COUT(49,:),[(kte-kts+1),(ite-its+1)]) * (10.169109344482422 - 0.0) + 0.0
    chem(p_cl_a04,kts:kte,its:ite)= reshape(COUT(50,:),[(kte-kts+1),(ite-its+1)]) * (151.50689697265625 - 0.0) + 0.0
    chem(p_nh4_a04,kts:kte,its:ite)= reshape(COUT(51,:),[(kte-kts+1),(ite-its+1)]) * (1.6307159662246704 - 0.0) + 0.0
    chem(p_na_a04,kts:kte,its:ite)= reshape(COUT(52,:),[(kte-kts+1),(ite-its+1)]) * (98.25704956054688 - 0.0) + 0.0
    !chem(p_seas_a04,kts:kte,its:ite)= (reshape(COUT(50,:),[(kte-kts+1),(ite-its+1)]) * (151.50689697265625 - 0.0) + 0.0) + (reshape(COUT(52,:),[(kte-kts+1),(ite-its+1)]) * (98.25704956054688 - 0.0) + 0.0)
    chem(p_oin_a04,kts:kte,its:ite)= reshape(COUT(53,:),[(kte-kts+1),(ite-its+1)]) * (33.43324279785156 - 0.0) + 0.0
    chem(p_dust_a04,kts:kte,its:ite)= reshape(COUT(54,:),[(kte-kts+1),(ite-its+1)]) * (31043.701171875 - 0.0) + 0.0
    chem(p_oc_a04,kts:kte,its:ite)= reshape(COUT(55,:),[(kte-kts+1),(ite-its+1)]) * (0.3992343842983246 - 0.0) + 0.0
    chem(p_bc_a04,kts:kte,its:ite)= reshape(COUT(56,:),[(kte-kts+1),(ite-its+1)]) * (0.16755656898021698 - 0.0) + 0.0
    chem(p_water_a04,kts:kte,its:ite)= reshape(COUT(57,:),[(kte-kts+1),(ite-its+1)]) * (6920.966796875 - 0.0) + 0.0
    chem(p_num_a04,kts:kte,its:ite)= reshape(COUT(58,:),[(kte-kts+1),(ite-its+1)]) *  (176295360.0 - 0.0) + 0.0
    chem(p_ca_a04,kts:kte,its:ite)= reshape(COUT(59,:),[(kte-kts+1),(ite-its+1)]) * (125.42910766601562 - 0.0) + 0.0
    chem(p_co3_a04,kts:kte,its:ite)= reshape(COUT(60,:),[(kte-kts+1),(ite-its+1)]) * (187.75692749023438 - 0.0) + 0.0
    chem(p_msa_a04,kts:kte,its:ite)= reshape(COUT(61,:),[(kte-kts+1),(ite-its+1)]) * (6.535012245178223 - 0.0) + 0.0

    deallocate(CIN)
    deallocate(COUT)
 
end subroutine AIMACI_mosaic_driver
end module module_AIMACI_mosaic
