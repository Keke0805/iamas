!************************************************************************
! This computer software was prepared by University of Science and Technology of China (USTC) 
! NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY, 
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Module to compute AIPC_saprc99 subroutine for chem(:,:,:,:). 
! * Author: Zihan Xia, ChenJin
! * Email: xiazh@mail.ustc.edu.cn
! Last update: zhxia 09/2024
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
 module module_AIPC_saprc99 
    use iamas_kind_types
    use iamas_atmchem_vars

 interface
 subroutine AIPC_saprc99(c_cin, c_cout, batch_size) bind(c, name='AIPC_saprc99')
    use iso_c_binding
    implicit none
    real(c_double), intent(in) :: c_cin(91,*)
    real(c_double), intent(in) :: c_cout(74,*)
    integer(c_int), value :: batch_size
 end subroutine AIPC_saprc99
 end interface

 contains
 subroutine AIPC_saprc99_driver(chem,                                    &
                         rho_phy,t_phy,p_phy,qv,coszen,                  &
                         ph, is_trop,                                    &
                         its,ite,kts,kte)
 IMPLICIT NONE

 INTEGER, INTENT(IN) :: its,ite,kts,kte
 REAL, DIMENSION(1:num_chem, kts:kte, its:ite), INTENT(INOUT) :: chem
 REAL, DIMENSION(kts:kte, its:ite), INTENT(IN) :: rho_phy, t_phy, p_phy, qv
 REAL, DIMENSION(its:ite), INTENT(IN) :: coszen
 REAL, DIMENSION(1:num_phr, kts:kte, its:ite), INTENT(IN) :: ph
 LOGICAL, DIMENSION(kts:kte, its:ite), INTENT(IN) :: is_trop


!local vars
 real(kind=RKIND), target, allocatable :: CIN(:,:), COUT(:,:)
 real(kind=RKIND), pointer :: c_cin(:,:), c_cout(:,:)
 integer :: batch_size
 real(kind=RKIND), allocatable :: conv(:,:), oconv(:,:), conv_w(:,:), temp_chem(:,:)
 
 ! A few constants 
 REAL(kind=RKIND), PARAMETER  ::  navgdro = 6.022e23   ! molecules/mol
 REAL(kind=RKIND), PARAMETER  ::   &   !molecular weights
      mwh = 1.0079,  mwo = 15.9994!,  mwair = 28.97
 REAL(kind=RKIND), PARAMETER :: mwh2o = 2*mwh + mwo
 !dens2con air
 REAL(kind=RKIND), PARAMETER ::          dens2con_a = 1.e-3     &! kg/m3 -> g/cm3
                              * (1./28.97)          &! -> mole/cm3
                              * navgdro              ! -> molec/cm3
 !dens2con water
 REAL(kind=RKIND), PARAMETER ::          dens2con_w = 1.e-3     &! kg/m3 -> g/cm3
                               * (1./mwh2o)         &! -> mole/cm3
                               * navgdro             ! -> molec/cm3
 integer:: i

 batch_size = (ite-its+1)*(kte-kts+1) 
 allocate(CIN(91,batch_size))
 allocate(COUT(74,batch_size))
 c_cin => CIN
 c_cout => COUT 
 allocate(conv(kte-kts+1,ite-its+1))
 allocate(oconv(kte-kts+1,ite-its+1))
 allocate(conv_w(kte-kts+1,ite-its+1))
 allocate(temp_chem(kte-kts+1,ite-its+1))
 conv= (1.E-6_8 * dens2con_a) / rho_phy(kts:kte,its:ite)
 oconv = 1.E0_8 / conv
 conv_w = dens2con_w / rho_phy(kts:kte,its:ite)

 CIN(1,:)  = (reshape(t_phy(kts:kte,its:ite),[batch_size]) - 174.655517578125) / (318.0460510253906 - 174.655517578125)
 CIN(2,:)  = (reshape(rho_phy(kts:kte,its:ite),[batch_size]) - 0.08756708353757858) / (1.5723837614059448 - 0.08756708353757858)
 CIN(3,:)  = (reshape(p_phy(kts:kte,its:ite)*conv_w,[batch_size])- 1.675948092166633e+24) / (3.063913686305045e+24 - 1.675948092166633e+24)
 CIN(4,:)  = (reshape(qv(kts:kte,its:ite)*conv_w,[batch_size]) - 0.0) / (7.703409361640489e+17 - 0.0) 
 CIN(5,:)  = (reshape(SPREAD(coszen(its:ite),1,kte-kts+1)*conv_w,[batch_size]) - -3.3939950106170937e+20) / (3.4916256618203474e+20 - -3.3939950106170937e+20)
 CIN(6,:)  = (reshape(ph(p_ph_o31d,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.0002218189329141751 - 0.0)
 CIN(7,:)  = (reshape(ph(p_ph_o33p,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.0016130889998748899 - 0.0)
 CIN(8,:)  = (reshape(ph(p_ph_no2,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.03153165802359581 - 0.0)
 CIN(9,:)  = (reshape(ph(p_ph_no3o2,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.08123581856489182 - 0.0)
 CIN(10,:) = (reshape(ph(p_ph_no3o,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.6305895447731018 - 0.0)
 CIN(11,:) = (reshape(ph(p_ph_hno2,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.0069686793722212315 - 0.0)
 CIN(12,:) = (reshape(ph(p_ph_hno3,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(3.1894467156234896e-06 - 0.0)
 CIN(13,:) = (reshape(ph(p_ph_hno4,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(2.225709795311559e-05 - 0.0)
 CIN(14,:) = (reshape(ph(p_ph_h2o2,kts:kte,its:ite)/60.0,[batch_size])  - 0.0) /(2.8756699975929223e-05 - 0.0)
 CIN(15,:) = (reshape(ph(p_ph_ch2or,kts:kte,its:ite)/60.0,[batch_size])  - 0.0) /(0.00013185555872041732 - 0.0)
 CIN(16,:) = (reshape(ph(p_ph_ch2om,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(0.00018208922119811177 - 0.0)
 CIN(17,:) = (reshape(ph(p_ph_ch3o2h,kts:kte,its:ite)/60.0,[batch_size]) - 0.0) /(2.2320618882076815e-05 - 0.0)
 CIN(18,:) = (reshape(chem(p_o3,kts:kte,its:ite)*conv,[batch_size]) - 8748936192.0) / (44287148425216.0 - 8748936192.0)
 CIN(19,:) = (reshape(chem(p_h2o2,kts:kte,its:ite)*conv,[batch_size]) - 0.002086453139781952) / (193582366720.0 - 0.002086453139781952)
 CIN(20,:) = (reshape(chem(p_no,kts:kte,its:ite)*conv,[batch_size]) - 0.001390912919305265) / (930006237184.0 - 0.001390912919305265)
 CIN(21,:) = (reshape(chem(p_no2,kts:kte,its:ite)*conv,[batch_size]) - 151955.859375) / (869180047360.0 - 151955.859375)
 CIN(22,:) = (reshape(chem(p_no3,kts:kte,its:ite)*conv,[batch_size]) - 35.502296447753906) / (5815655424.0 - 35.502296447753906)
 CIN(23,:) = (reshape(chem(p_n2o5,kts:kte,its:ite)*conv,[batch_size]) - 0.0014518782263621688) / (45500178432.0 - 0.0014518782263621688)
 CIN(24,:) = (reshape(chem(p_hono,kts:kte,its:ite)*conv,[batch_size]) - 0.0014587603509426117) / (18956605440.0 - 0.0014587603509426117)
 CIN(25,:) = (reshape(chem(p_hno3,kts:kte,its:ite)*conv,[batch_size]) - 86134.1171875) / (236281266176.0 - 86134.1171875)
 CIN(26,:) = (reshape(chem(p_hno4,kts:kte,its:ite)*conv,[batch_size]) - 993.3319702148438) / (16048924672.0 - 993.3319702148438)
 CIN(27,:) = (reshape(chem(p_so2,kts:kte,its:ite)*conv,[batch_size]) - 37398.2265625) / (1130485645312.0 - 37398.2265625)
 CIN(28,:) = (reshape(chem(p_h2so4,kts:kte,its:ite)*conv,[batch_size]) - 23441.302734375) / (95366389760.0 - 23441.302734375)
 CIN(29,:) = (reshape(chem(p_co,kts:kte,its:ite)*conv,[batch_size]) - 1059878010880.0) / (16950667247616.0 - 1059878010880.0)
 CIN(30,:) = (reshape(chem(p_hcho,kts:kte,its:ite)*conv,[batch_size]) - 14953136.0) / (124107382784.0 - 14953136.0)
 CIN(31,:) = (reshape(chem(p_ccho,kts:kte,its:ite)*conv,[batch_size]) - 23925002.0) / (113933901824.0- 23925002.0)
 CIN(32,:) = (reshape(chem(p_rcho,kts:kte,its:ite)*conv,[batch_size]) - 5534125.0) / (84019200000.0 - 5534125.0)
 CIN(33,:) = (reshape(chem(p_acet,kts:kte,its:ite)*conv,[batch_size]) - 18555.388671875) / (42330890240.0 - 18555.388671875)
 CIN(34,:) = (reshape(chem(p_mek,kts:kte,its:ite)*conv,[batch_size])- 29510.95703125) / (34007363584.0 - 29510.95703125)
 CIN(35,:) = (reshape(chem(p_hcooh,kts:kte,its:ite)*conv,[batch_size]) - 10948.896484375) / (82452439040.0 - 10948.896484375)
 CIN(36,:) = (reshape(chem(p_meoh,kts:kte,its:ite)*conv,[batch_size]) - 19156.693359375) / (133054611456.0 - 19156.693359375)
 CIN(37,:) = (reshape(chem(p_etoh,kts:kte,its:ite)*conv,[batch_size]) - 0.0013220069231465459) / (2.369162082672119 - 0.0013220069231465459)
 CIN(38,:) = (reshape(chem(p_cco_oh,kts:kte,its:ite)*conv,[batch_size]) - 0.13837814331054688) / (127888531456.0 - 0.13837814331054688)
 CIN(39,:) = (reshape(chem(p_rco_oh,kts:kte,its:ite)*conv,[batch_size]) - 0.13837814331054688) / (179384156160.0 - 0.13837814331054688)
 CIN(40,:) = (reshape(chem(p_gly,kts:kte,its:ite)*conv,[batch_size]) - 0.0015090556116774678) / (6388580352.0 - 0.0015090556116774678)
 CIN(41,:) = (reshape(chem(p_mgly,kts:kte,its:ite)*conv,[batch_size]) - 0.0015090556116774678) / (14115335168.0 - 0.0015090556116774678)
 CIN(42,:) = (reshape(chem(p_bacl,kts:kte,its:ite)*conv,[batch_size]) - 0.0013805871130898595) / (7488878080.0 - 0.0013805871130898595)
 CIN(43,:) = (reshape(chem(p_cres,kts:kte,its:ite)*conv,[batch_size]) - 0.0013780617155134678) / (4438989312.0 - 0.0013780617155134678)
 CIN(44,:) = (reshape(chem(p_bald,kts:kte,its:ite)*conv,[batch_size]) - 0.001638737739995122) / (6668132864.0 - 0.001638737739995122)
 CIN(45,:) = (reshape(chem(p_isoprod,kts:kte,its:ite)*conv,[batch_size]) - 0.0015981380129233003) / (23915241472.0 - 0.0015981380129233003)
 CIN(46,:) = (reshape(chem(p_methacro,kts:kte,its:ite)*conv,[batch_size]) - 0.0016600314993411303) / (18255792128.0 - 0.0016600314993411303)
 CIN(47,:) = (reshape(chem(p_mvk,kts:kte,its:ite)*conv,[batch_size]) - 0.0016589591978117824) / (28600266752.0 - 0.0016589591978117824)
 CIN(48,:) = (reshape(chem(p_prod2,kts:kte,its:ite)*conv,[batch_size]) - 5828.4052734375) / (23783424000.0 - 5828.4052734375)
 CIN(49,:) = (reshape(chem(p_dcb1,kts:kte,its:ite)*conv,[batch_size]) - 0.0013459211913868785) / (13960432640.0 - 0.0013459211913868785)
 CIN(50,:) = (reshape(chem(p_dcb2,kts:kte,its:ite)*conv,[batch_size]) - 0.0013414109125733376) / (2707870976.0 - 0.0013414109125733376)
 CIN(51,:) = (reshape(chem(p_dcb3,kts:kte,its:ite)*conv,[batch_size]) - 0.0013414109125733376) / (2045967744.0 - 0.0013414109125733376)
 CIN(52,:) = (reshape(chem(p_ethene,kts:kte,its:ite)*conv,[batch_size]) - 0.0019865958020091057) / (163621634048.0 - 0.0019865958020091057)
 CIN(53,:) = (reshape(chem(p_isoprene,kts:kte,its:ite)*conv,[batch_size]) - 0.001379182212986052) / (132673085440.0 - 0.001379182212986052)
 CIN(54,:) = (reshape(chem(p_c2h6,kts:kte,its:ite)*conv,[batch_size]) - 0.0014128655893728137) / (131573071872.0 - 0.0014128655893728137)
 CIN(55,:) = (reshape(chem(p_c3h8,kts:kte,its:ite)*conv,[batch_size]) - 187.5325469970703) / (178172166144.0 - 187.5325469970703)
 CIN(56,:) = (reshape(chem(p_c2h2,kts:kte,its:ite)*conv,[batch_size]) - 0.0013885501539334655) / (74210066432.0 - 0.0013885501539334655)
 CIN(57,:) = (reshape(chem(p_c3h6,kts:kte,its:ite)*conv,[batch_size]) - 0.0013780617155134678) / (5616316416.0 - 0.0013780617155134678)
 CIN(58,:) = (reshape(chem(p_alk3,kts:kte,its:ite)*conv,[batch_size]) - 11175690.0) / (783694364672.0 - 11175690.0)
 CIN(59,:) = (reshape(chem(p_alk4,kts:kte,its:ite)*conv,[batch_size]) - 0.0053211054764688015) / (482662678528.0 - 0.0053211054764688015)
 CIN(60,:) = (reshape(chem(p_alk5,kts:kte,its:ite)*conv,[batch_size]) - 0.00196873489767313) / (198980730880.0 - 0.00196873489767313)
 CIN(61,:) = (reshape(chem(p_aro1,kts:kte,its:ite)*conv,[batch_size]) - 0.0014093111967667937) / (189882204160.0 - 0.0014093111967667937)
 CIN(62,:) = (reshape(chem(p_aro2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014398486819118261) / (106865385472.0 - 0.0014398486819118261)
 CIN(63,:) = (reshape(chem(p_ole1,kts:kte,its:ite)*conv,[batch_size]) - 0.001668077427893877) / (130782576640.0 - 0.001668077427893877)
 CIN(64,:) = (reshape(chem(p_ole2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014852883759886026) / (119988920320.0 - 0.0014852883759886026)
 CIN(65,:) = (reshape(chem(p_terp,kts:kte,its:ite)*conv,[batch_size]) - 0.0013780617155134678) / (1173483136.0 - 0.0013780617155134678)
 CIN(66,:) = (reshape(chem(p_rno3,kts:kte,its:ite)*conv,[batch_size]) - 297952.96875) / (23537358848.0 - 297952.96875)
 CIN(67,:) = (reshape(chem(p_nphe,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (2908993536.0 - 0.0014160138089209795)
 CIN(68,:) = (reshape(chem(p_phen,kts:kte,its:ite)*conv,[batch_size]) - 0.0014188283821567893) / (532085760.0 - 0.0014188283821567893)
 CIN(69,:) = (reshape(chem(p_pan,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (77684858880.0 - 0.0014160138089209795)
 CIN(70,:) = (reshape(chem(p_pan2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (49266462720.0 - 0.0014160138089209795)
 CIN(71,:) = (reshape(chem(p_pbzn,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (1599492224.0 - 0.0014160138089209795)
 CIN(72,:) = (reshape(chem(p_ma_pan,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (2675750656.0 - 0.0014160138089209795)
 CIN(73,:) = (reshape(chem(p_co2,kts:kte,its:ite)*conv,[batch_size]) - 123146183114752.0) / (9.118790236032205e+16 - 123146183114752.0)
 CIN(74,:) = (reshape(chem(p_cco_ooh,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (114752839680.0 - 0.0014160138089209795)
 CIN(75,:) = (reshape(chem(p_rco_o2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (329765088.0 - 0.0014160138089209795)
 CIN(76,:) = (reshape(chem(p_rco_ooh,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (131599777792.0 - 0.0014160138089209795)
 CIN(77,:) = (reshape(chem(p_xn,kts:kte,its:ite)*conv,[batch_size]) - 0.0014131629141047597) / (20129126400.0 - 0.0014131629141047597)
 CIN(78,:) = (reshape(chem(p_xc,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (86337347584.0 - 0.0014160138089209795)
 CIN(79,:) = (reshape(chem(p_ho,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (576040704.0 - 0.0014160138089209795)
 CIN(80,:) = (reshape(chem(p_ho2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (6171183104.0 - 0.0014160138089209795)
 CIN(81,:) = (reshape(chem(p_c_o2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (4586074624.0 - 0.0014160138089209795)
 CIN(82,:) = (reshape(chem(p_cooh,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (137433899008.0 - 0.0014160138089209795)
 CIN(83,:) = (reshape(chem(p_rooh,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (84788510720.0 - 0.0014160138089209795)
 CIN(84,:) = (reshape(chem(p_ro2_r,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (9205670912.0 - 0.0014160138089209795)
 CIN(85,:) = (reshape(chem(p_r2o2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (5634583040.0 - 0.0014160138089209795)
 CIN(86,:) = (reshape(chem(p_ro2_n,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (890188992.0 - 0.0014160138089209795)
 CIN(87,:) = (reshape(chem(p_cco_o2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (509194368.0 - 0.0014160138089209795)
 CIN(88,:) = (reshape(chem(p_bzco_o2,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (6886416.5 - 0.0014160138089209795) 
 CIN(89,:) = (reshape(chem(p_ma_rco3,kts:kte,its:ite)*conv,[batch_size]) - 0.0014160138089209795) / (56238272.0 - 0.0014160138089209795)
 CIN(90,:) = (reshape(chem(p_sesq,kts:kte,its:ite)*conv,[batch_size]) - 0.0013220069231465459) / (0.023738397285342216 - 0.0013220069231465459)
 CIN(91,:) = (reshape(chem(p_ch4,kts:kte,its:ite)*conv,[batch_size]) - 16766035034112.0) / (413048112152576.0 - 16766035034112.0)

 where (cin > 1.0)
     cin = 1.0
 elsewhere (cin < 0.0)
     cin = 0.0
 end where

 call AIPC_saprc99(c_cin, c_cout, batch_size)
 !where (COUT<0.0_8)
 !       COUT = 0.0_8
 !end where 

 temp_chem = (reshape(COUT(1,:),[(kte-kts+1),(ite-its+1)]) * (44289987969024.0 - 0.001104875234887004) + 0.001104875234887004) * oconv!*real(grid%chemdt*60,8)+chem(p_o3)
 where (is_trop)
    chem(p_o3,kts:kte,its:ite) = temp_chem
 end where 
 
 temp_chem = (reshape(COUT(2,:),[(kte-kts+1),(ite-its+1)]) * (193371668480.0 - 13.832404136657715) + 13.832404136657715) * oconv!*real(grid%chemdt*60,8)+chem(p_h2o2)
 where (is_trop)
    chem(p_h2o2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(3,:),[(kte-kts+1),(ite-its+1)]) * (529722212352.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_no)
 where (is_trop)
    chem(p_no,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(4,:),[(kte-kts+1),(ite-its+1)]) * (1028682350592.0 - 151625.171875) + 151625.171875) * oconv!*real(grid%chemdt*60,8)+chem(p_no2)
 where (is_trop)
    chem(p_no2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(5,:),[(kte-kts+1),(ite-its+1)]) * (6179588096.0 - 9.100873947143555) + 9.100873947143555) * oconv!*real(grid%chemdt*60,8)+chem(p_no3)
 where (is_trop)
    chem(p_no3,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(6,:),[(kte-kts+1),(ite-its+1)]) * (47564476416.0 - 0.2034233808517456) + 0.2034233808517456) * oconv!*real(grid%chemdt*60,8)+chem(p_n2o5)
 where (is_trop)
    chem(p_n2o5,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(7,:),[(kte-kts+1),(ite-its+1)]) * (35072090112.0 - 0.001647678087465465) + 0.001647678087465465) * oconv!*real(grid%chemdt*60,8)+chem(p_hono)
 where (is_trop)
    chem(p_hono,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(8,:),[(kte-kts+1),(ite-its+1)]) * (235947950080.0 - 130814.234375) + 130814.234375) * oconv!*real(grid%chemdt*60,8)+chem(p_hno3)
 where (is_trop)
    chem(p_hno3,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(9,:),[(kte-kts+1),(ite-its+1)]) * (23849406464.0 - 905.9033203125) + 905.9033203125) * oconv!*real(grid%chemdt*60,8)+chem(p_hno4)
 where (is_trop)
    chem(p_hno4,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(10,:),[(kte-kts+1),(ite-its+1)]) * (1117987799040.0 - 37388.453125) +37388.453125) * oconv!*real(grid%chemdt*60,8)+chem(p_so2)
 where (is_trop)
    chem(p_so2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(11,:),[(kte-kts+1),(ite-its+1)]) * (116486684672.0 - 23452.83984375) + 23452.83984375) * oconv!*real(grid%chemdt*60,8)+chem(p_h2so4)
 where (is_trop)
    chem(p_h2so4,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(12,:),[(kte-kts+1),(ite-its+1)]) * (16950055927808.0 - 1059909337088.0) + 1059909337088.0) * oconv!*real(grid%chemdt*60,8)+chem(p_co)
 where (is_trop)
    chem(p_co,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(13,:),[(kte-kts+1),(ite-its+1)]) * (124839337984.0 - 17033054.0) + 17033054.0) * oconv !*real(grid%chemdt*60,8)+chem(p_hcho)
 where (is_trop)
    chem(p_hcho,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(14,:),[(kte-kts+1),(ite-its+1)]) * (113960927232.0 - 23691516.0) + 23691516.0) * oconv!*real(grid%chemdt*60,8)+chem(p_ccho)
 where (is_trop)
    chem(p_ccho,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(15,:),[(kte-kts+1),(ite-its+1)]) * (84089585664.0 - 5869394.5) + 5869394.5) * oconv!*real(grid%chemdt*60,8)+chem(p_rcho)
 where (is_trop)
    chem(p_rcho,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(16,:),[(kte-kts+1),(ite-its+1)]) * (57383546880.0 - 77005.7734375) + 77005.7734375) * oconv!*real(grid%chemdt*60,8)+chem(p_acet)
 where (is_trop)
    chem(p_acet,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(17,:),[(kte-kts+1),(ite-its+1)]) * (39518765056.0 - 174973.90625) + 174973.90625) * oconv!*real(grid%chemdt*60,8)+chem(p_mek)
 where (is_trop)
    chem(p_mek,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(18,:),[(kte-kts+1),(ite-its+1)]) * (84062806016.0 - 62976.43359375) + 62976.43359375) * oconv!*real(grid%chemdt*60,8)+chem(p_hcooh)
 where (is_trop)
    chem(p_hcooh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(19,:),[(kte-kts+1),(ite-its+1)]) * (132863983616.0 - 61178.20703125) + 61178.20703125) * oconv!*real(grid%chemdt*60,8)+chem(p_meoh)
 where (is_trop)
    chem(p_meoh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem  = (reshape(COUT(20,:),[(kte-kts+1),(ite-its+1)]) * (2.369162082672119 - 0.0013220069231465459) + 0.0013220069231465459) * oconv!*real(grid%chemdt*60,8)+chem(p_etoh)
 where (is_trop)
    chem(p_etoh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(21,:),[(kte-kts+1),(ite-its+1)]) * (127888670720.0 - 8973963.0) + 8973963.0) * oconv!*real(grid%chemdt*60,8)+chem(p_cco_oh)
 where (is_trop)
    chem(p_cco_oh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(22,:),[(kte-kts+1),(ite-its+1)]) * (179385040896.0 - 4367208.0) + 4367208.0) * oconv!*real(grid%chemdt*60,8)+chem(p_rco_oh)
 where (is_trop)
    chem(p_rco_oh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(23,:),[(kte-kts+1),(ite-its+1)]) * (9547779072.0 - 0.0019587166607379913) + 0.0019587166607379913) * oconv!*real(grid%chemdt*60,8)+chem(p_gly)
 where (is_trop)
    chem(p_gly,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(24,:),[(kte-kts+1),(ite-its+1)]) * (18614052864.0 - 0.0019222002010792494) + 0.0019222002010792494) * oconv!*real(grid%chemdt*60,8)+chem(p_mgly)
 where (is_trop)
    chem(p_mgly,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(25,:),[(kte-kts+1),(ite-its+1)]) * (8430501888.0 - 0.0014508608728647232) + 0.0014508608728647232) * oconv!*real(grid%chemdt*60,8)+chem(p_bacl)
 where (is_trop)
    chem(p_bacl,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(26,:),[(kte-kts+1),(ite-its+1)]) * (7712128512.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_cres)
 where (is_trop)
    chem(p_cres,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(27,:),[(kte-kts+1),(ite-its+1)]) * (7649209856.0 - 0.0005344006349332631) + 0.0005344006349332631) * oconv!*real(grid%chemdt*60,8)+chem(p_bald)
 where (is_trop)
    chem(p_bald,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(28,:),[(kte-kts+1),(ite-its+1)]) * (23722289152.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_isoprod)
 where (is_trop)
    chem(p_isoprod,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(29,:),[(kte-kts+1),(ite-its+1)]) * (18127364096.0 - 1.9919511032640003e-05) + 1.9919511032640003e-05) * oconv!*real(grid%chemdt*60,8)+chem(p_methacro)
 where (is_trop)
    chem(p_methacro,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(30,:),[(kte-kts+1),(ite-its+1)]) * (28457889792.0 - 0.00023820197384338826) + 0.00023820197384338826) * oconv!*real(grid%chemdt*60,8)+chem(p_mvk)
 where (is_trop)
    chem(p_mvk,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(31,:),[(kte-kts+1),(ite-its+1)]) * (28906188800.0 - 21682.693359375) + 21682.693359375) * oconv!*real(grid%chemdt*60,8)+chem(p_prod2)
 where (is_trop)
    chem(p_prod2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(32,:),[(kte-kts+1),(ite-its+1)]) * (23042674688.0 - 0.00024165434297174215) + 0.00024165434297174215) * oconv!*real(grid%chemdt*60,8)+chem(p_dcb1)
 where (is_trop)
    chem(p_dcb1,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(33,:),[(kte-kts+1),(ite-its+1)]) * (4469666304.0 - 6.992577255005017e-05) + 6.992577255005017e-05) * oconv!*real(grid%chemdt*60,8)+chem(p_dcb2)
 where (is_trop)
    chem(p_dcb2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(34,:),[(kte-kts+1),(ite-its+1)]) * (3367498496.0 - 5.382955350796692e-05) + 5.382955350796692e-05) * oconv !*real(grid%chemdt*60,8)+chem(p_dcb3)
 where (is_trop)
    chem(p_dcb3,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(35,:),[(kte-kts+1),(ite-its+1)]) * (160583221248.0 - 0.0012187260435894132) + 0.0012187260435894132) * oconv!*real(grid%chemdt*60,8)+chem(p_ethene)
 where (is_trop)
    chem(p_ethene,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(36,:),[(kte-kts+1),(ite-its+1)]) * (128428212224.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_isoprene)
 where (is_trop)
    chem(p_isoprene,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(37,:),[(kte-kts+1),(ite-its+1)]) * (131534331904.0 - 0.001412823679856956) + 0.001412823679856956) * oconv!*real(grid%chemdt*60,8)+chem(p_c2h6)
 where (is_trop)
    chem(p_c2h6,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(38,:),[(kte-kts+1),(ite-its+1)]) * (178172166144.0 - 187.5325469970703) + 187.5325469970703) * oconv!*real(grid%chemdt*60,8)+chem(p_c3h8)
 where (is_trop)
    chem(p_c3h8,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(39,:),[(kte-kts+1),(ite-its+1)]) * (43726483456.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_c2h2)
 where (is_trop)
    chem(p_c2h2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(40,:),[(kte-kts+1),(ite-its+1)]) * (2034101504.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_c3h6)
 where (is_trop)
    chem(p_c3h6,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(41,:),[(kte-kts+1),(ite-its+1)]) * (780088705024.0 - 11171813.0) + 11171813.0) * oconv!*real(grid%chemdt*60,8)+chem(p_alk3)
 where (is_trop)
    chem(p_alk3,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(42,:),[(kte-kts+1),(ite-its+1)]) * (477917970432.0 - 0.0052979872561991215) + 0.0052979872561991215) * oconv!*real(grid%chemdt*60,8)+chem(p_alk4)
 where (is_trop)
    chem(p_alk4,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(43,:),[(kte-kts+1),(ite-its+1)]) * (192240009216.0 - 0.0010533719323575497) + 0.0010533719323575497) * oconv!*real(grid%chemdt*60,8)+chem(p_alk5)
 where (is_trop)
    chem(p_alk5,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(44,:),[(kte-kts+1),(ite-its+1)]) * (174464827392.0 - 0.0012793251080438495) + 0.0012793251080438495) * oconv!*real(grid%chemdt*60,8)+chem(p_aro1)
 where (is_trop)
    chem(p_aro1,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(45,:),[(kte-kts+1),(ite-its+1)]) * (100872265728.0 - 6.892842066008598e-05) + 6.892842066008598e-05) * oconv!*real(grid%chemdt*60,8)+chem(p_aro2)
 where (is_trop)
    chem(p_aro2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(46,:),[(kte-kts+1),(ite-its+1)]) * (129196818432.0 - 1.3761081390839536e-05) + 1.3761081390839536e-05) * oconv!*real(grid%chemdt*60,8)+chem(p_ole1)
 where (is_trop)
    chem(p_ole1,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(47,:),[(kte-kts+1),(ite-its+1)]) * (109858570240.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_ole2)
 where (is_trop)
    chem(p_ole2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(48,:),[(kte-kts+1),(ite-its+1)]) * (749003072.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_terp)
 where (is_trop)
    chem(p_terp,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(49,:),[(kte-kts+1),(ite-its+1)]) * (30166493184.0 - 293461.84375) + 293461.84375) * oconv!*real(grid%chemdt*60,8)+chem(p_rno3)
 where (is_trop)
    chem(p_rno3,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(50,:),[(kte-kts+1),(ite-its+1)]) * (4327222784.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_nphe)
 where (is_trop)
    chem(p_nphe,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(51,:),[(kte-kts+1),(ite-its+1)]) * (458057344.0 - 0.00010952115553664044) + 0.00010952115553664044) * oconv!*real(grid%chemdt*60,8)+chem(p_phen)
 where (is_trop)
    chem(p_phen,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(52,:),[(kte-kts+1),(ite-its+1)]) * (78866006016.0 - 8154.09375) + 8154.09375) * oconv!*real(grid%chemdt*60,8)+chem(p_pan)
 where (is_trop)
    chem(p_pan,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(53,:),[(kte-kts+1),(ite-its+1)]) * (49892229120.0 - 16568.75) + 16568.75) * oconv!*real(grid%chemdt*60,8)+chem(p_pan2)
 where (is_trop)
    chem(p_pan2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(54,:),[(kte-kts+1),(ite-its+1)]) * (1606644736.0 - 0.0008262745686806738) + 0.0008262745686806738) * oconv!*real(grid%chemdt*60,8)+chem(p_pbzn)
 where (is_trop)
    chem(p_pbzn,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(55,:),[(kte-kts+1),(ite-its+1)]) * (2676982784.0 - 5.008662265026942e-05) + 5.008662265026942e-05) * oconv!*real(grid%chemdt*60,8)+chem(p_ma_pan)
 where (is_trop)
    chem(p_ma_pan,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(56,:),[(kte-kts+1),(ite-its+1)]) * (9.096986405057331e+16 - 122836249214976.0) + 122836249214976.0) * oconv!*real(grid%chemdt*60,8)+chem(p_co2)
 where (is_trop)
    chem(p_co2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(57,:),[(kte-kts+1),(ite-its+1)]) * (114754846720.0 - 27936.412109375) + 27936.412109375) * oconv !*real(grid%chemdt*60,8)+chem(p_cco_ooh)
 where (is_trop)
    chem(p_cco_ooh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(58,:),[(kte-kts+1),(ite-its+1)]) * (299636672.0 - 0.0308065228164196) + 0.0308065228164196) * oconv!*real(grid%chemdt*60,8)+chem(p_rco_o2)
 where (is_trop)
    chem(p_rco_o2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(59,:),[(kte-kts+1),(ite-its+1)]) * (131599785984.0 - 7347.54833984375) + 7347.54833984375) * oconv!*real(grid%chemdt*60,8)+chem(p_rco_ooh)
 where (is_trop)
    chem(p_rco_ooh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(60,:),[(kte-kts+1),(ite-its+1)]) * (20174114816.0 - 3.1975600719451904) + 3.1975600719451904) * oconv!*real(grid%chemdt*60,8)+chem(p_xn)
 where (is_trop)
    chem(p_xn,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(61,:),[(kte-kts+1),(ite-its+1)]) * (86741893120.0 - 189.1082000732422) + 189.1082000732422) * oconv!*real(grid%chemdt*60,8)+chem(p_xc)
 where (is_trop)
    chem(p_xc,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(62,:),[(kte-kts+1),(ite-its+1)]) * (586730368.0 - 0.024739215150475502) + 0.024739215150475502) * oconv!*real(grid%chemdt*60,8)+chem(p_ho)
 where (is_trop)
    chem(p_ho,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(63,:),[(kte-kts+1),(ite-its+1)]) * (6175001600.0 - 4.339116096496582) + 4.339116096496582) * oconv!*real(grid%chemdt*60,8)+chem(p_ho2)
 where (is_trop)
    chem(p_ho2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(64,:),[(kte-kts+1),(ite-its+1)]) * (4586921472.0 - 253.14813232421875) + 253.14813232421875)* oconv!*real(grid%chemdt*60,8)+chem(p_c_o2)
 where (is_trop)
    chem(p_c_o2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(65,:),[(kte-kts+1),(ite-its+1)]) * (137498804224.0 - 904504.5) + 904504.5) * oconv!*real(grid%chemdt*60,8)+chem(p_cooh)
 where (is_trop)
    chem(p_cooh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(66,:),[(kte-kts+1),(ite-its+1)]) * (85366456320.0 - 10115715.0) + 10115715.0) * oconv!*real(grid%chemdt*60,8)+chem(p_rooh)
 where (is_trop)
    chem(p_rooh,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(67,:),[(kte-kts+1),(ite-its+1)]) * (9000088576.0 - 0.019492199644446373) + 0.019492199644446373)* oconv!*real(grid%chemdt*60,8)+chem(p_ro2_r)
 where (is_trop)
    chem(p_ro2_r,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(68,:),[(kte-kts+1),(ite-its+1)]) * (5399159296.0 - 0.021797887980937958) + 0.021797887980937958) * oconv!*real(grid%chemdt*60,8)+chem(p_r2o2)
 where (is_trop)
    chem(p_r2o2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(69,:),[(kte-kts+1),(ite-its+1)]) * (858729280.0 - 0.014648937620222569) + 0.014648937620222569) * oconv!*real(grid%chemdt*60,8)+chem(p_ro2_n)
 where (is_trop)
    chem(p_ro2_n,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(70,:),[(kte-kts+1),(ite-its+1)]) * (510029440.0 - 0.12708553671836853) + 0.12708553671836853) * oconv!*real(grid%chemdt*60,8)+chem(p_cco_o2)
 where (is_trop)
    chem(p_cco_o2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(71,:),[(kte-kts+1),(ite-its+1)]) * (5249698.5 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_bzco_o2)
 where (is_trop)
    chem(p_bzco_o2,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(72,:),[(kte-kts+1),(ite-its+1)]) * (56520888.0 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_ma_rco3)
 where (is_trop)
    chem(p_ma_rco3,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(73,:),[(kte-kts+1),(ite-its+1)]) * (0.022717637941241264 - 0.0) + 0.0) * oconv!*real(grid%chemdt*60,8)+chem(p_sesq)
 where (is_trop)
    chem(p_sesq,kts:kte,its:ite) = temp_chem
 end where

 temp_chem = (reshape(COUT(74,:),[(kte-kts+1),(ite-its+1)]) * (413047776608256.0 - 16766021402624.0) + 16766021402624.0) * oconv!*real(grid%chemdt*60,8)+chem(p_ch4)
 where (is_trop)
    chem(p_ch4,kts:kte,its:ite) = temp_chem
 end where

    deallocate(CIN)
    deallocate(COUT)
    deallocate(conv)
    deallocate(oconv)
    deallocate(conv_w)
    deallocate(temp_chem)

end subroutine AIPC_saprc99_driver
end module module_AIPC_saprc99
