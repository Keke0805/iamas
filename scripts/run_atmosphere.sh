#!/bin/bash
#SBATCH -p batch
#SBATCH -J u60_2019
#SBATCH -o test_iAMAS.log
#SBATCH -e test_iAMAS.err
#SBATCH -N 4 -n 112
#SBATCH --ntasks-per-node=28

#########################################################################
# 计算节点运行此脚本执行的命令为：sbatch run_atmosphere.sh
# Jiawang Feng, jw525000@mail.ustc.edu.cn;
#
# History:
# 1.  该脚本不支持restart，如果需要, 建议手动
# 2.  加入对blkrun的支持.
# 3.  优化文本输出.
# 4.  同时监控log.atmosphere, 并输出完成进度.
# 5.  qkyang: 实现可以估算每次分块模拟计算还需多少时间(eta time)，2024年05月17日
# 6.  qkyang: added ${codeDir}URBPARM.TBL to support urban model，2024年9月4日
# 7.  qkyang: 增加支持人为热排放和LCZ数据，2025年3月25日
# 8.  zhxia: 开发photolysis模块. 
# 9.  zhxia: 开发emission模块.
# 10. zhxia: 开发化学初始化模块.
# 11. zhxia: 实现读取其他模式输出结果、数据集作为iAMAS初始场.
# 12. zhxia: 开发气相化学模拟框架，实现AI气相化学方案在线运行.
# 13. zhxia: 开发人为源气溶胶化学模拟框架，实现AI气溶胶化学方案在线运行.
# 14. qkyang: 增加支持Shin-Hong边界层方案，2025年7月16日
# 15. zhxia: 增加平流层臭氧linoz方案.
# 16. jwfeng: 增加deep-version dynamic的参数设置，2025年11月26日
# 17. zhxia: 开发KPP框架，实现SAPRC99方案的在线运行.
# 18. xyruan: 分开配置二次无机(SIA)和二次有机(SOA)气溶胶化学方案选项，2026年1月7日
# 19. ybxue: 增加可选择的沙尘源函数功能支持(erod和ssm)，2026年4月15日
# 20. xyruan: 增加CMU云化学方案的运行支持，2026年6月12日
# 21. zhxia: 实现定制化init_final文件中包含的变量
#########################################################################
simStartTime=$(date +%s)
echo "Simulation Start at $(date "+%Y-%m-%d %H:%M:%S")"
#########################################################################
# 以下为可修改参数
# 修改时等号左右不允许出现空格!

debugFile=/dev/null

# blk run说明:
# 0. 运行surface以及第一个init final, 如果第一个init需要chem_input则需要在./output中存在对应时刻的iamas_out
# 1. 运行${blkDays}天
# 2. 将生成的iamas_out以及iamas_diag文件移动到./output中
# 3. 将上个run blk最后时刻提前${spinUpDays}天之前的iamas_out链接为${chemInputFileName}作为下一个blk的chem_input
# 4. 基于上述${chemInputFileName}进行带chem_input的init final, 然后开始循环
# 5. 如果doChemInput=false则不进行上述与chem_input相关的操作

# 重要设置:
# ------------------重要设置开始-------------------------
doBlkRun=false #时间是否切开分块模拟
blkDays=5 # 每次模拟时长
spinUpDays=1 # spin-up时长，输出结果时长则为blkDays-spinUpDays
doChemInput=0 # 1 for do, 0 for don't do
chemOpt=no_chemistry # 'no_chemistry','bin04aer_nature_aq','bin04aer_gas_aq'
doEmisUpdate=false #zhxia add this to control whether the emission_update file is generated
photOpt=no_phot # 'no_phot','photfastj'
gasOPT=no_gaschem # 'no_gaschem','AIPC_saprc99', 'KPP_saprc99'
straGasOPT=no_stra_gaschem # 'no_stra_gaschem','linoz'
CldchemOPT=no_cloudchem # 'no_cloudchem','CMU_cloudchem'
aerosolInorgOPT=no_sia # 'no_sia','AIMACI_mosaic','simple_SIA', 'MOSAIC_mech'
aerosolOrgOPT=no_soa # 'no_soa','simple_SOA'

# 模式起止时间, 重要!

startYear=2018
startMonth=12
startDay=27
startHour=0

endYear=2019
endMonth=2
endDay=1
endHour=0

# 使用的code, 重要!
codeDir=/path_to_iAMAS_source_code/ # code存放的路径, 无需重复编译
buildDir=/path_to_iAMAS_build/
# 包含分辨率/网格信息, 重要!
gridGraphDir=/data/Data_Share/iAMAS/inputdata/meshes/60km/grid_graph_namelist_files/ # 必须包含decompFile, 有init_static.nc,graph.info.part.xx
inputdataDir=/data/Data_Share/iAMAS/inputdata # Intermediate文件存放的路径(例如ERA5制作的气象场), 需要链接为./inputdata
mpTableDir=/lustre/jwfeng/work/MPAS/standard_run_dir/ # MP_THOMPSON tables存放的路径, 无需重复生成

# dt等, 重要!
dt=300
len_disp=60000
apvm_upwinding=0.5 # 在最大分辨率达到4km时设置为0.0
# ------------------重要设置结束-------------------------
staticFileName=init_static.nc #default: init_static.nc
emissionFileName=emission_update.nc #zhxia
surfaceFileName=surface_update.nc
initFinalFileName=init_final.nc
chemInputFileName=chem_input.nc
decompPrefix=graph.info.part.  # decomposition file prefix
# emission设置,zhxia
emis_fgInterval=3600                 # 读取Intermediate文件的时间间隔, 单位: s
emissionInterval=3600                # 生成emission数据的时间间隔, 单位: s
anthro_emis_Prefix=inputdata/emissions/HTAPv32_global_0.1x0.1_SAPRC/emis_anthro # 读取的中间文件prefix
use_anthro_emis=false
biomass_emis_Prefix=inputdata/emissions/FINNv25_global_0.1x0.1_SAPRC/emis_bioms # 读取的中间文件prefix
use_bioms_emis=false
emis_interp_type=variable #fine2coarse, variable, coarse2fine
# surface设置
sfc_fgInterval=21600               # 读取Intermediate文件的时间间隔, 单位: s
surfaceInterval=21600          # 生成surface文件的时间间隔, 单位: s
surfacePrefix=inputdata/metdata/ERA5 # 读取的中间文件prefix
# init final设置
nVertLevels=55
nSoilLevels=4
nFgLevels=38 #38 for ERA5 42 for FNL
nFgSoilLevels=4
nFgGasLevels=49 #zhxia
metPrefix=inputdata/metdata/ERA5
sourceOpt=erod #ybxue: erod or ssm
#sourcePrefix=inputdata/emissions/emissions_source_ssm_500m.dat #ybxue: use SSM file when sourceOpt=ssm
sourcePrefix=inputdata/emissions/emissions_source_dust_dms.dat #ybxue: use EROD file when sourceOpt=erod
AH_emis_Prefix=inputdata/emissions/emissions_anthropogenic_heating_ZJU_2019_1km.dat #qkyang
fix_gasinit_file=inputdata/chemdata/saprc99/iAMAS_Global_gasinit_saprc99.dat #zhxia
use_fix_gasinit=false #zhxia
fix_ozone_init_file=inputdata/chemdata/ozone_TOST/iAMAS_Global_ozone_init_2019_DJF_TOST.dat #zhxia
use_fix_ozone_init=false #zhxia
useSpechumd=true
zTop=30000

# For Deep-version dynamics: jwfeng
epssm_zb=30000 # jwfeng: epssm开始增加的海拔高度
epssm_zt=80000 # jwfeng: epssm不再变化的海拔高度
doKineticVisc=false # deep-version需打开
doRayleighDampU=false # deep-version需打开
nSubSteps=2 # deep-version时需适当增加

# atmosphere设置
# 由于configs无法预先覆盖完全，所以将configs放在数组中，如需增加可以在数组末尾添加,同样的,等号左右不允许出现空格
# 以下没有的config为atmosphere的默认config

# atmosphere &nhyd_model 设置
# Notice: config_start_time and config_run_duration is re-set for every run blocks so not listed below
nhyd_model_configs=(
    config_dt="${dt}"
    config_photdt=1800
    config_chemdt=3600
    config_number_of_sub_steps=${nSubSteps}
    config_len_disp="${len_disp}"
    config_apvm_upwinding="${apvm_upwinding}"
    config_include_kinetic_visc=${doKineticVisc}
    config_del4u_div_factor=1
    config_epssm_zb=${epssm_zb}
    config_epssm_zt=${epssm_zt}
)
damping=(
    config_rayleigh_damp_u=${doRayleighDampU}
)
# atmosphere &physics 设置
physics_configs=(
    config_radt_cld_scheme=\'cld_fraction\'
    config_gwdo_scheme=\'bl_ysu_gwdo\'
    config_o3climatology=true
    isconvection_tstar=false #qkyang: true/false to turn on/off convection Tstar model [ChengYu et al. (2021) PRF model]
    config_convection_scheme=\'cu_grell_freitas\' # 'suite','cu_kain_fritsch','cu_kain_fritsch_sca','cu_tiedtke','cu_ntiedtke','cu_grell_freitas','off'
    config_microp_scheme=\'mp_thompson\' # 'suite','mp_wsm6','mp_thompson','mp_morrison','mp_kessler','off'
    config_radt_lw_scheme=\'rrtmg_lw\' # 'suite','rrtmg_lw','cam_lw','off'
    config_radt_sw_scheme=\'rrtmg_sw\' # `suite',`rrtmg_sw',`cam_sw',`off'
    config_lsm_scheme=\'noah\' #qkyang: =\'noah\' to support urban model
    config_urban_scheme=0 #qkyang: =0,1,2,3 Non-urban, Urban, BEP, BEP_BEM #only works with YSU, MYJ and BouLac PBL schemes
    # Using anthropogenic heating emissions requires AH_emis_Prefix and config_urban_scheme=1
    is_ah_emis=false #qkyang: =true/false turn on/off anthropogenic heating emissions
    config_use_lcz=false #qkyang added to use LCZ data with 11 land cover types for urban model
    # num_urban_layers=400 #qkyang set 400/5400 to support BEP/BEP_BEM urban model
    config_physics_suite=\'convection_permitting\' #\'convection_permitting\' #\'mesoscale_reference\' 
    # config_pbl_scheme=\'bl_mynn\' #qkyang: At present, the compatibility issues between bl_ysu and BEP and BEP_BEM have not been resolved (running will result in an error).
    # config_pbl_scheme=\'bl_ysu\' #\'bl_mynn\'
    # config_sfclayer_scheme=\'sf_monin_obukhov\' #\'sf_mynn\'
    config_sst_update=true
    config_aer_ra_feedback=false
    config_build_turb_n=false #qkyang added to support building turbulence model in YSU model
    config_build_turb_d=false #qkyang added to support daytime building turbulence model in YSU model
    config_topo_wind=0 #qkyang added to support topographic wind model, 0:off, 1:on, 2:on with vertical wind
    config_slope_rad=0 #qkyang added to support 1: enable slope-dependent shortwave correction at surface; 0: disable"
    config_topo_shading=0 #qkyang add to support 1: enable neighboring-cell topographic shading for shortwave direct beam; 0: disable"
)
AHOPTION=0 #qkyang:=1表示使用默认的查表人为热，非is_ah_emis的人为热清单，
# atmosphere &chemistry
chemistry_configs=(
    config_chem_opt="\'${chemOpt}\'"
    config_dust_factor=1.0
    config_seas_factor=1.0
    config_use_aer_IN=false
    do_seas_sst_adjust=true
    config_ctrans_opt=\'general\' # 'general','kfeta','no_ctrans', kfeta must be with cu_kain_fritsch/cu_kain_fritsch_sca
    config_cu_scav_on=true
    config_unsoluble_IN_size=0.5
    config_soluble_IN_size=0.1
    config_drydep_on=true  # dry deposition and mix activation
    config_wetscav_on=true # large-scale wet scavenging
    config_dust_opt=\"gocart\"              # dust emission scheme
    config_source_opt="'${sourceOpt}'"      # dust source function selection
    config_emis_update="${doEmisUpdate}"
    config_use_anthro_emis="${use_anthro_emis}"
    config_use_bioms_emis="${use_bioms_emis}"
    config_use_fix_gasinit="${use_fix_gasinit}"
    config_phot_opt="'${photOpt}'"
    config_gaschem_opt="'${gasOPT}'"
    config_stra_gaschem_opt="'${straGasOPT}'"
    config_cloudchem_opt="'${CldchemOPT}'"
    config_aerchem_sia_opt="'${aerosolInorgOPT}'"
    config_aerchem_soa_opt="'${aerosolOrgOPT}'"
)
# atmosphere &printout
printout_configs=(
    config_print_debug_info=false
    config_print_global_minmax_vel=true
    config_print_detailed_minmax_vel=false
)
# atmosphere &io
io_configs=(
    config_pio_num_iotasks=0 # 做I/O的核数，0代表所有核都参与I/O
    config_pio_stride=1 # 参与I/O的核之间的间隔，与config_pio_num_iotasks的乘积要等于nProcs
)

# Please do not modify the code below this line!
#########################################################################

if [[ $SLURM_NPROCS -eq "" ]]; then
    echo "This script only support running in SLURM environment! Use 'sbatch' to submit it!"
    exit 1
else
    nProcs=$SLURM_NPROCS
fi
echo "Using $nProcs processors"

#########################################################################

# 输出整行的字符,并将需要输出的提示信息放在中间
# $1为行宽度, $2为字符(e.g. "-"), $3为提示信息
function print_whole_line() {
    local delim totalWidth info infoLen preLen postLen
    if [ ! "${2}" ]; then
        delim='-'
    else
        delim=${2}
    fi
    if [ ! "${3}" ]; then
        info=''
        infoLen=0
    else
        info=$3
        infoLen=${#info}
    fi
    totalWidth=$1
    if [[ ${infoLen} -ge ${totalWidth} ]]; then
        echo -e "${info}"
    else
        preLen=$(((totalWidth-infoLen)/2))
        postLen=$((totalWidth-infoLen-preLen))
        echo -n $(yes $delim | sed $preLen'q' | tr -d '\n')
        echo -n -e "${info}"
        echo    $(yes $delim | sed $postLen'q'| tr -d '\n')
    fi
}
infoLength=60
function print_error_info() {
    print_whole_line ${infoLength} "-" "Fatal Error"
    echo -e "\e[41;37;5m${1}\e[0m" #输出红底白字闪烁字符
    print_whole_line ${infoLength} "-" # "\e[31mFatal Error\e[0m"
}

# 检查该路径下的必须文件, 以及参数设置是否合理

# 检查时间格式是否正确
startTimeStr=$(printf "%04d-%02d-%02d %02d:00:00" ${startYear} ${startMonth} ${startDay} ${startHour})
endTimeStr=$(printf "%04d-%02d-%02d %02d:00:00" ${endYear} ${endMonth} ${endDay} ${endHour})
date -d "${startTimeStr}" > ${debugFile}
if [[ $? -ne 0 ]]; then
    print_error_info "Invalid start time: ${startTimeStr}"
    exit 1
fi
date -d "${endTimeStr}" > ${debugFile}
if [[ $? -ne 0 ]]; then
    print_error_info "Invalid end time: ${startTimeStr}"
    exit 1
fi

# 检查开始时间是否早于结束时间
startTimeInt=$(printf "%04d%02d%02d%02d" ${startYear} ${startMonth} ${startDay} ${startHour})
endTimeInt=$(printf "%04d%02d%02d%02d" ${endYear} ${endMonth} ${endDay} ${endHour})
if [[ ${startTimeInt} -ge ${endTimeInt} ]]; then
    print_error_info "Invalid time period: ${startTimeStr} - ${endTimeStr}"
    exit 1
fi

# copy&link essential files
linkFiles=(
    "${gridGraphDir}${staticFileName}"
    "${gridGraphDir}${decompPrefix}${nProcs}"
    "${inputdataDir}"
    "${codeDir}src/atmosphere/physics/schemes/files/GENPARM.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/LANDUSE.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/OZONE_DAT.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/OZONE_LAT.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/OZONE_PLEV.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/RRTMG_LW_DATA"
    "${codeDir}src/atmosphere/physics/schemes/files/RRTMG_LW_DATA.DBL"
    "${codeDir}src/atmosphere/physics/schemes/files/RRTMG_SW_DATA"
    "${codeDir}src/atmosphere/physics/schemes/files/RRTMG_SW_DATA.DBL"
    "${codeDir}src/atmosphere/physics/schemes/files/SOILPARM.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/VEGPARM.TBL"
    "${codeDir}src/atmosphere/physics/schemes/files/URBPARM.TBL" #added by qkyang
    "${codeDir}src/atmosphere/physics/schemes/files/URBPARM_LCZ.TBL" #added by qkyang
    "${mpTableDir}MP_THOMPSON_QIautQS_DATA.DBL"
    "${mpTableDir}MP_THOMPSON_QRacrQG_DATA.DBL"
    "${mpTableDir}MP_THOMPSON_QRacrQS_DATA.DBL"
    "${mpTableDir}MP_THOMPSON_freezeH2O_DATA.DBL"
  )

# 不要copy stream_list.atmosphere.*, 需要预先设置在当前文件夹中
copyFiles=(
    "${buildDir}bin/atmosphere"
    "${buildDir}bin/init_atmosphere"
    "${buildDir}iAMAS/core_atmosphere/namelist.atmosphere"
    "${buildDir}iAMAS/core_init_atmosphere/namelist.init_atmosphere"
    "${buildDir}iAMAS/core_init_atmosphere/streams.init_atmosphere"
    "${buildDir}iAMAS/core_atmosphere/streams.atmosphere"
    "${codeDir}src/atmosphere/chemistry/schemes/AIPC_saprc99.pt"
    "${codeDir}src/atmosphere/chemistry/schemes/AIMACI_mosaic.pt"
    "${codeDir}src/atmosphere/chemistry/schemes/Linoz_March2007.dat"
  )
for fileName in "${linkFiles[@]}" ; do
    if [[ ! -e ${fileName} ]]; then
        print_error_info "${fileName} does not exist"
        exit 1
    else
        if [[ "${fileName}" == "${inputdataDir}" ]]; then
            rm -rf ./inputdata #zhxia add this for solving the error when run this script multiple times in the same dir.
            ln -sf "${fileName}" ./inputdata
        else
            ln -sf "${fileName}" ./
        fi
    fi
done

# zhxia add this for linking emis_anthro to generate emission_update.
# Due to permission, we can't generate thses files in the user dir. You need to login in datauser and run link_generate.sh in /data/Data_Share/iAMAS/inputdata/emissions/HTAPv2_global_0.1x0.1_SAPRC

# added by qkyang to support anthropogenic heating emissions>>>>>>>>>>>>>>
# 将数组元素转换为多行文本
physics_configs_text="${physics_configs[*]}"
physics_configs_text="${physics_configs_text// /$'\n'}"
# 使用 grep 和 awk 提取is_ah_emis值
is_ah_emis=$(echo "$physics_configs_text" | grep '^is_ah_emis=' | awk -F '=' '{print $2}' | sed 's/ #.*//')
config_use_lcz_phys=$(echo "$physics_configs_text" | grep '^config_use_lcz=' | awk -F '=' '{print $2}' | sed 's/ #.*//')
if [ "$config_use_lcz_phys" = true ]; then
    URBPARM_TBL="URBPARM_LCZ.TBL"
else
    config_use_lcz_phys=false
    URBPARM_TBL="URBPARM.TBL"
fi
if [ ! -f $URBPARM_TBL ]; then # 检查URBPARM.TBL是否存在
    echo "Error: File URBPARM.TBL is not exist!"
    exit 1
else
    if [ "$is_ah_emis" = true ]; then
        # 如果 is_ah_emis 为 true，将 AHOPTION: 后面的内容替换为 1
        sed -i 's/^AHOPTION: [0-9]/AHOPTION: 1/' $URBPARM_TBL
    elif [[ "$is_ah_emis" = false || -z "$is_ah_emis" ]]; then
        # 如果 is_ah_emis 为 false 或者为空，将 AHOPTION: 后面的内容替换为 0
        is_ah_emis=false
        if [[ "$AHOPTION" = 1 ]]; then
            echo "Attention!: AHOPTION=1 in URBPARM.TBL, but is_ah_emis=false"
            sed -i 's/^AHOPTION: [0-9]/AHOPTION: 1/' $URBPARM_TBL
        else
            sed -i 's/^AHOPTION: [0-9]/AHOPTION: 0/' $URBPARM_TBL
        fi
    else
        echo "Error: is_ah_emis="$is_ah_emis", which must be true or false"
        exit 1
    fi
fi
# added by qkyang to support anthropogenic heating emissions<<<<<<<<<<<<<

#Added by qkyang to support slope and shading effects>>>>>>>>>>>>>>>>>
config_slope_rad=$(echo "$physics_configs_text" | grep '^config_slope_rad=' | awk -F '=' '{print $2}' | sed 's/ #.*//')
config_topo_shading=$(echo "$physics_configs_text" | grep '^config_topo_shading=' | awk -F '=' '{print $2}' | sed 's/ #.*//')
#Added by qkyang to support slope and shading effects<<<<<<<<<<<<<<<<

# added by qkyang to support building turbulence model>>>>>>>>>>>>>>
config_build_turb_d=$(echo "$physics_configs_text" | grep '^config_build_turb_d=' | awk -F '=' '{print $2}' | sed 's/ #.*//')
config_build_turb_n=$(echo "$physics_configs_text" | grep '^config_build_turb_n=' | awk -F '=' '{print $2}' | sed 's/ #.*//')
# added by qkyang to support building turbulence model<<<<<<<<<<<<<

for fileName in "${copyFiles[@]}" ; do
    if [[ ! -e ${fileName} ]]; then
        print_error_info "${fileName} does not exist!"
        exit 1
    else
        cp "${fileName}" ./
    fi
done
# 检查文件是否存在
essentialFiles=(
    "${staticFileName}"
    "${decompPrefix}${nProcs}"
    "${surfacePrefix}:$(printf "%04d-%02d-%02d_%02d" ${startYear} ${startMonth} ${startDay} ${startHour})"
    "${surfacePrefix}:$(printf "%04d-%02d-%02d_%02d" ${endYear} ${endMonth} ${endDay} ${endHour})"
    atmosphere
    init_atmosphere
    GENPARM.TBL
    LANDUSE.TBL
    MP_THOMPSON_QIautQS_DATA.DBL
    MP_THOMPSON_QRacrQG_DATA.DBL
    MP_THOMPSON_QRacrQS_DATA.DBL
    MP_THOMPSON_freezeH2O_DATA.DBL
    OZONE_DAT.TBL
    OZONE_LAT.TBL
    OZONE_PLEV.TBL
    RRTMG_LW_DATA
    RRTMG_LW_DATA.DBL
    RRTMG_SW_DATA
    RRTMG_SW_DATA.DBL
    SOILPARM.TBL
    VEGPARM.TBL
    URBPARM.TBL
    URBPARM_LCZ.TBL
    namelist.atmosphere
    namelist.init_atmosphere
    stream_list.atmosphere.diagnostics
    stream_list.atmosphere.output
    stream_list.atmosphere.surface
    stream_list.atmosphere.emission
    stream_list.init_atmosphere.init_final
    streams.init_atmosphere
    streams.atmosphere)
for fileName in ${essentialFiles[*]} ; do
    if [[ ! -e ${fileName} ]]; then
        print_error_info "No ${fileName} in current dir. Must include it!"
        exit 1
    fi
done

# 检查dt等重要设置是否正确
if [[ ${dt} -ge $((len_disp/1000*6)) ]]; then
    print_error_info "config_dt=${dt} too long for config_len_disp=${len_disp}!"
    exit 1
fi
if [[ ${len_disp} -lt 4001 && "${apvm_upwinding}" != "0.0" ]]; then
    print_error_info "config_apvm_upwinding should be 0.0 if config_len_disp<=4000!"
    exit 1
fi

# 检查结束
#########################################################################
startTimeStr=`printf "%04d-%02d-%02d_%02d:00:00" ${startYear} ${startMonth} ${startDay} ${startHour}`
endTimeStr=`printf "%04d-%02d-%02d_%02d:00:00" ${endYear} ${endMonth} ${endDay} ${startHour}`

#  set namelist.atmosphere
function set_namelist() {
    # $1为config名字, $2为对应的config设置内容, $3为namelist的文件名
    local configName configValue namelistFileName i
    configName=($1) # 将传递进来的字符串分解为数组
    configValue=($2)
    namelistFileName=($3)
    # 输出修改的行但是先不修改文件
    changedLine=`sed -n "/${configName}/ c\ ${configName} = ${configValue}" $namelistFileName`
    if [[ $changedLine == "" ]]; then
        print_error_info "No ${configName} in ${namelistFileName}. Must include it!"
        exit 1
    else
        echo "=>$changedLine"
    fi
    # 修改文件
    sed -i "/${configName}/ c\ ${configName} = ${configValue}" $namelistFileName
}

# zhxia
function check_emission_file() {
    # 接受参数$1为emissionfile的文件名
    # 检查是否存在符合起止时间条件的emission file, 存在返回1, 不存在返回0
    # 由于xtime一定是等差数列, 所以只要startTime/startTime+interval/endTime-interval/endTime在xtime数组中即可
    # 注意不能用return, 否则在return非0的数时会导致shell认为出错
    local xtime s1 e1
    s1=`date -d "${startTimeStr/_/ } ${emissionInterval} seconds"     "+%Y-%m-%d_%H:%M:%S"`
    e1=`date -d "${endTimeStr/_/ } ${emissionInterval} seconds ago"   "+%Y-%m-%d_%H:%M:%S"`
    if [ -e ${emissionFileName} ]; then
        xtime=`ncdump -v xtime ${1} | sed -n "/xtime =/, /}/ p"`
        if [[ "${xtime}" =~ "${startTimeStr}" && "${xtime}" =~ "${endTimeStr}" && \
              "${xtime}" =~ "${s1}" && \
              "${xtime}" =~ "${e1}" ]]; then
            echo 1
        else
            echo "Emission file exists but the time period doesn't match ${startTimeStr} - ${endTimeStr}"
        fi
    else
        echo "Emission file does not exist"
    fi
}

# zhxia
function run_emission() {
    local emissionStat
    emissionStat=`check_emission_file ${emissionFileName}`
    if [[ "${emissionStat}" == 1 ]]; then
        echo "=> \"${emissionFileName}\" exists, pass"
    else
        echo "=> ${emissionStat}"
        set_namelist config_init_case                10                         ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_start_time               \'${startTimeStr}\'        ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_stop_time                \'${endTimeStr}\'          ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_anthro_emis_prefix       \'${anthro_emis_Prefix}\'  ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_use_anthro_emis          ${use_anthro_emis}         ./namelist.init_atmosphere > ${debugFile} 
        set_namelist config_bioms_emis_prefix        \'${biomass_emis_Prefix}\' ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_use_bioms_emis           ${use_bioms_emis}          ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_emis_interp_type         \'${emis_interp_type}\'    ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_fg_interval              ${emis_fgInterval}         ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_static_interp            false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_vertical_grid            false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_met_interp               false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_input_sst                false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_input_chem               false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_block_decomp_file_prefix \'${decompPrefix}\'        ./namelist.init_atmosphere > ${debugFile}

        # Edit streams.init_atmosphere
        sed -i "/stream name=\"input\"/, /filename_template/ s/filename_template.*$/filename_template=\"${staticFileName}\"/" ./streams.init_atmosphere
        sed -i "/stream name=\"input\"/ a\                  io_type=\"netcdf\"" ./streams.init_atmosphere
        sed -i "/stream name=\"emission\"/, /filename_template/ s/filename_template.*$/filename_template=\"${emissionFileName}\"/" ./streams.init_atmosphere
        sed -i "/stream name=\"emission\"/, /output_interval/ s/output_interval.*$/output_interval=\"${emissionInterval}\" \/>/" ./streams.init_atmosphere
        sed -i "/stream name=\"emission\"/ a\                  clobber_mode=\"replace_files\"" ./streams.init_atmosphere
        sed -i "/stream name=\"emission\"/ a\                  io_type=\"pnetcdf,cdf5\"" ./streams.init_atmosphere
        sed -i '/stream name="lbc"/,/\/stream/ d' ./streams.init_atmosphere # 运行emission时删除stream "lbc"
        sed -i "/stream name=\"chem_input\"/, /\/>/ d" ./streams.init_atmosphere # 运行emission时删除stream "chem_input"
        sed -i '/stream name="output"/,/\/stream/ s|stream_list.init_atmosphere.output|stream_list.init_atmosphere.init_final|' streams.init_atmosphere # zhxia: 通过stream_list.init_atmosphere.init_final控制init_final.nc中所包含的变量
        #sed -i '/stream name="output"/,/\/stream/ d' ./streams.init_atmosphere # 运行emission时删除stream "output"
        sed -i "/^\s*$/ d" ./streams.init_atmosphere # 删除空行
        echo "=> emission start time is: ${startTimeStr}"
        echo "=> emission   end time is: ${endTimeStr}"
        # 开始运行emission
        mpirun -n ${nProcs} ./init_atmosphere
        # 备份emission设置及输出
        mkdir -p ./log/emission
        cp ./namelist.init_atmosphere ./log/emission/namelist.init_atmosphere_emission
        cp ./streams.init_atmosphere  ./log/emission/streams.init_atmosphere_emission
        cp ./log.init_atmosphere*     ./log/emission/
        # 检查运行结果是否正确
        emissionStat=`check_emission_file ${emissionFileName}`
        if [[ ! "${emissionStat}" == 1 ]]; then
            print_error_info "Emission run failed. Please check log.init_atmosphere.*"
            exit 1
        fi
    fi
}


function check_surface_file() {
    # 接受参数$1为surfacefile的文件名
    # 检查是否存在符合起止时间条件的surface file, 存在返回1, 不存在返回0
    # 由于xtime一定是等差数列, 所以只要startTime/startTime+interval/endTime-interval/endTime在xtime数组中即可
    # 注意不能用return, 否则在return非0的数时会导致shell认为出错
    local xtime s1 e1
    s1=`date -d "${startTimeStr/_/ } ${surfaceInterval} seconds"     "+%Y-%m-%d_%H:%M:%S"`
    e1=`date -d "${endTimeStr/_/ } ${surfaceInterval} seconds ago"   "+%Y-%m-%d_%H:%M:%S"`
    if [ -e ${surfaceFileName} ]; then
        xtime=`ncdump -v xtime ${1} | sed -n "/xtime =/, /}/ p"`
        if [[ "${xtime}" =~ "${startTimeStr}" && "${xtime}" =~ "${endTimeStr}" && \
              "${xtime}" =~ "${s1}" && \
              "${xtime}" =~ "${e1}" ]]; then
            echo 1
        else
            echo "Surface file exists but the time period doesn't match ${startTimeStr} - ${endTimeStr}"
        fi
    else
        echo "Surface file does not exist"
    fi
}

function run_surface() {
    local surfaceStat
    # TODO:
    # 检查gridGraphDir中是否存在符合起止时间条件的surface file: 还是算了吧，没必要
    surfaceStat=`check_surface_file ${surfaceFileName}`
    if [[ "${surfaceStat}" == 1 ]]; then
        echo "=> \"${surfaceFileName}\" exists, pass"
    else
        echo "=> ${surfaceStat}"
        set_namelist config_init_case                8                    ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_start_time               \'${startTimeStr}\'  ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_stop_time                \'${endTimeStr}\'    ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_sfc_prefix               \'${surfacePrefix}\' ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_fg_interval              ${sfc_fgInterval}        ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_static_interp            false                ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_vertical_grid            false                ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_met_interp               false                ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_input_sst                true                 ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_input_chem               false                ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_block_decomp_file_prefix \'${decompPrefix}\'  ./namelist.init_atmosphere > ${debugFile}

        # Edit streams.init_atmosphere
        sed -i "/stream name=\"input\"/, /filename_template/ s/filename_template.*$/filename_template=\"${staticFileName}\"/" ./streams.init_atmosphere
        sed -i "/stream name=\"input\"/ a\                  io_type=\"netcdf\"" ./streams.init_atmosphere
        sed -i "/stream name=\"surface\"/, /filename_template/ s/filename_template.*$/filename_template=\"${surfaceFileName}\"/" ./streams.init_atmosphere
        sed -i "/stream name=\"surface\"/, /output_interval/ s/output_interval.*$/output_interval=\"${surfaceInterval}\" \/>/" ./streams.init_atmosphere
        sed -i "/stream name=\"surface\"/ a\                  clobber_mode=\"replace_files\"" ./streams.init_atmosphere
        sed -i "/stream name=\"surface\"/ a\                  io_type=\"pnetcdf,cdf5\"" ./streams.init_atmosphere
        sed -i '/stream name="lbc"/,/\/stream/ d' ./streams.init_atmosphere # 运行surface时删除stream "lbc"
        sed -i "/stream name=\"chem_input\"/, /\/>/ d" ./streams.init_atmosphere # 运行surface时删除stream "chem_input"
        sed -i '/stream name="output"/,/\/stream/ s|stream_list.init_atmosphere.output|stream_list.init_atmosphere.init_final|' streams.init_atmosphere # zhxia: 通过stream_list.init_atmosphere.init_final控制init_final.nc中所包含的变量
        # sed -i '/stream name="output"/,/\/stream/ d' ./streams.init_atmosphere # 运行surface时删除stream "output"
        sed -i "/^\s*$/ d" ./streams.init_atmosphere # 删除空行
        echo "=> surface start time is: ${startTimeStr}"
        echo "=> surface   end time is: ${endTimeStr}"
        # 开始运行surface
        mpirun -n ${nProcs} ./init_atmosphere
        # 备份surface设置及输出
        mkdir -p ./log/surface
        cp ./namelist.init_atmosphere ./log/surface/namelist.init_atmosphere_surface
        cp ./streams.init_atmosphere  ./log/surface/streams.init_atmosphere_surface
        cp ./log.init_atmosphere*     ./log/surface/
        # 检查运行结果是否正确
        surfaceStat=`check_surface_file ${surfaceFileName}`
        if [[ ! "${surfaceStat}" == 1 ]]; then
            print_error_info "Surface run failed. Please check log.init_atmosphere.*"
            exit 1
        fi
    fi
}

# 检查init_final.nc的xtime是否正确,该函数需要传递startTime, 以及是否进行chem_input
function check_init_final_file() {
    # $1 is startTimeStr for this round of init_final, $2 is flag for chem_input, 1 for do chem_input, 0 for don't do.
    # 检查是否存在符合起止时间条件的init_final file, 存在返回1, 不存在返回0
    # 注意不能用return, 否则在return非0的数时会导致shell认为出错
    local sTime dustValues
    if [ -e ${initFinalFileName} ]; then
        sTime=`ncdump -v xtime ${initFinalFileName} | grep -A 1 "xtime =" | grep -v xtime`
        if [[ "${sTime}" =~ "${1}" ]]; then
            dustValues=`ncdump -v dust_a01 ${initFinalFileName} | sed -n "/dust_a01 =/, /}/ p"`
            if [[ ${2} -eq 1 ]]; then # do chem_input时要确保dustValues中存在非0的数, 且不能全是1e-16,暂时判断时存在小数点即可
                if [[ "${dustValues}" =~ \. ]]; then # 此处可以用\.或者"." 不能用"\."
                    echo 1
                else
                    echo "Init final file exists but the value of chem is wrong while doing chem_input"
                fi
            else # 不做chem_input时要确保dustValues中全是1e-16,存在小数点或者全是0即判断失败
                if [[ "${dustValues}" =~ \. ]]; then # 由于全是0的情况与全是1e-16结果没区别
                    echo "Init final file exists but the value of chem is wrong while NOT doing chem_input"
                else
                    echo 1
                fi
            fi
        else
            echo "Init final file exists but the xtime is not matching ${1}"
        fi

    else
        echo "Init final file does not exist"
    fi
}

# 生成init_final.nc, 由于init_final.nc需要重复生成, 所以该函数需要传递startTimeStr, 注意格式!!!
function run_init_final() {
    local logFile initFinalStat
    # $1为init_final的开始时间字符串, $2为是否进行chem_input, 1 for yes, 0 for no
    initFinalStat=`check_init_final_file ${1} ${2}`
    if [[ "${initFinalStat}" == 1 ]]; then
        echo "=> \"${initFinalFileName}\" exists, pass"
    else
        echo "=> ${initFinalStat}"
		cp ${buildDir}iAMAS/core_init_atmosphere/namelist.init_atmosphere ./
		cp ${buildDir}iAMAS/core_init_atmosphere/streams.init_atmosphere ./
        set_namelist config_init_case                7                          ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_theta_adv_order          3                          ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_start_time               \'${1}\'                   ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_nvertlevels              ${nVertLevels}             ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_nsoillevels              ${nSoilLevels}             ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_nfglevels                ${nFgLevels}               ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_nfgsoillevels            ${nFgSoilLevels}           ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_nfggaslevels             ${nFgGasLevels}            ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_met_prefix               \'${metPrefix}\'           ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_source_opt               \'${sourceOpt}\'           ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_source_file              \'${sourcePrefix}\'        ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_AH_emis_file             \'${AH_emis_Prefix}\'      ./namelist.init_atmosphere > ${debugFile}
        set_namelist is_ah_emis                      ${is_ah_emis}              ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_slope_rad                ${config_slope_rad}        ./namelist.init_atmosphere > ${debugFile} #qkyang
        set_namelist config_topo_shading             ${config_topo_shading}     ./namelist.init_atmosphere > ${debugFile} #qkyang
        set_namelist config_build_turb_d             ${config_build_turb_d}     ./namelist.init_atmosphere > ${debugFile} #qkyang
        set_namelist config_build_turb_n             ${config_build_turb_n}     ./namelist.init_atmosphere > ${debugFile} #qkyang
        set_namelist config_fix_gasinit_file         \'${fix_gasinit_file}\'    ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_use_fix_gasinit          ${use_fix_gasinit}         ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_fix_ozone_init_file      \'${fix_ozone_init_file}\' ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_use_fix_ozone_init       ${use_fix_ozone_init}      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_native_gwd_static        false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_use_spechumd             ${useSpechumd}             ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_ztop                     ${zTop}                    ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_nsmterrain               2                          ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_smooth_surfaces          true                       ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_static_interp            false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_vertical_grid            true                       ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_met_interp               true                       ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_input_sst                false                      ./namelist.init_atmosphere > ${debugFile}
        set_namelist config_block_decomp_file_prefix \'${decompPrefix}\'        ./namelist.init_atmosphere > ${debugFile}
        #qkyang: change config_specified_zeta_levels to specific_zw.txt>>>>>>>>>>>>>
        if [ -e "./specific_zw.txt" ]; then #如果specific_zw.txt文件存在，就将config_specified_zeta_levels设置为specific_zw.txt
            echo "Attention: specific_zw.txt exists, set config_specified_zeta_levels=specific_zw.txt!"
            set_namelist config_specified_zeta_levels \'specific_zw.txt\' ./namelist.init_atmosphere > ${debugFile}
        fi
        #qkyang: change config_specified_zeta_levels to specific_zw.txt<<<<<<<<<<<<<
        #qkyang: make sure landuse name is right+++++++++++
        landuse_data=$(ncdump -h ${gridGraphDir}${staticFileName} | sed -n 's/.*:config_landuse_data =\([^;]*\).*/\1/p')
        set_namelist config_landuse_data             ${landuse_data}     ./namelist.init_atmosphere > ${debugFile}
        #qkyang: make sure landuse name is right-----------
        #qkyang: make sure LCZ switch is right>>>>>>>>>>>>
        config_use_lcz_static=$(ncdump -h ${gridGraphDir}${staticFileName} | sed -n 's/.*:config_use_lcz =\([^;]*\).*/\1/p')
        config_use_lcz_static=$(echo "$config_use_lcz_static" | awk '{gsub(/[ '\''"]/, ""); print}') #同时去除空格和引号。
        if [[ $config_use_lcz_static == 'YES' ]]; then
            set_namelist config_use_lcz true ./namelist.init_atmosphere > ${debugFile}
            if [ "$config_use_lcz_phys" = false ]; then
                echo "Error: config_use_lcz=YES in init_static.nc, then config_use_lcz must be true in physics_configs"
                exit 1
            fi
        else
            set_namelist config_use_lcz false ./namelist.init_atmosphere > ${debugFile}
            if [ "$config_use_lcz_phys" = true ]; then
                echo "Error: config_use_lcz not = YES in init_static.nc, then config_use_lcz must be false in physics_configs"
                exit 1                
            fi
        fi
        #qkyang: make sure LCZ switch is right<<<<<<<<<<<

        if [[ ${2} -eq 0 ]]; then
            set_namelist config_input_chem false ./namelist.init_atmosphere
            sed -i "/stream name=\"chem_input\"/, /\/>/ d" ./streams.init_atmosphere # 删除stream "chem_input"
        else
            set_namelist config_input_chem true  ./namelist.init_atmosphere
            sed -i "/stream name=\"chem_input\"/ a\                  io_type=\"pnetcdf,cdf5\"" ./streams.init_atmosphere
            sed -i "/stream name=\"chem_input\"/, /filename_template/ s/filename_template.*$/filename_template=\"${chemInputFileName}\"/" ./streams.init_atmosphere
        fi

        # Edit streams.init_atmosphere
        sed -i "/stream name=\"input\"/, /filename_template/ s/filename_template.*$/filename_template=\"${staticFileName}\"/" ./streams.init_atmosphere
        sed -i "/stream name=\"input\"/ a\                  io_type=\"netcdf\"" ./streams.init_atmosphere
        sed -i "/stream name=\"output\"/, /filename_template/ s/filename_template.*$/filename_template=\"${initFinalFileName}\"/" ./streams.init_atmosphere
        sed -i "/stream name=\"output\"/ a\                  clobber_mode=\"replace_files\"" ./streams.init_atmosphere
        sed -i "/stream name=\"output\"/ a\                  io_type=\"pnetcdf,cdf5\"" ./streams.init_atmosphere
        sed -i '/stream name="output"/,/\/stream/ s|stream_list.init_atmosphere.output|stream_list.init_atmosphere.init_final|' streams.init_atmosphere # zhxia: 通过stream_list.init_atmosphere.init_final控制init_final.nc中所包含的变量
        sed -i "/stream name=\"emission\"/, /\/>/ d" ./streams.init_atmosphere # 运行init final时删除stream "emission"
        sed -i '/stream name="lbc"/,/\/stream/ d' ./streams.init_atmosphere # 运行init final时删除stream "lbc"
        sed -i "/stream name=\"surface\"/, /\/>/ d" ./streams.init_atmosphere # 运行init final时删除stream "surface"
        sed -i "/^\s*$/ d" ./streams.init_atmosphere # 删除空行
        echo "=> init final time is: ${1}"
        # 开始运行init final
        # added by qkyang to support anthropogenic heating emissions>>>>>>>>>>>>>>
        if [ "$is_ah_emis" = true ]; then
            #因为人为热是高分辨率数据，所以限制并行计算核心数为nProcs_init个，避免内存不足
            nProcs_init=24
            linkGraph_init="${gridGraphDir}graph.info.part.${nProcs_init}" #链接10核计算的文件路径
            # 检查源文件是否存在
            if [ -e "$linkGraph_init" ]; then
                # 创建软链接
                ln -sf "$linkGraph_init" ./
            else
                echo "is_ah_emis=true requires $linkGraph_init"
                echo "Error: $linkGraph_init is not found"
                echo "(1) let: is_ah_emis=false"
                echo "(2) or: cd ${gridGraphDir} && gpmetis ./graph.info $nProcs_init"
                exit 1
            fi
            mpirun -n ${nProcs_init} ./init_atmosphere # ${nProcs} ./init_atmosphere
        else
            mpirun -n ${nProcs} ./init_atmosphere   #default
        fi
        # added by qkyang to support anthropogenic heating emissions<<<<<<<<<<<<< 
        # 备份init final设置及输出
        mkdir -p ./log/init_final
        cp ./namelist.init_atmosphere ./log/init_final/namelist.init_atmosphere_init_final_${1}
        cp ./streams.init_atmosphere  ./log/init_final/streams.init_atmosphere_init_final_${1}
        for logFile in $(ls ./log.init_atmosphere*); do
            cp ${logFile} ./log/init_final/${logFile}_${1}
        done
        # 检查运行结果是否正确
        initFinalStat=`check_init_final_file ${1} ${2}`
        if [[ ! "${initFinalStat}" == 1 ]]; then
            print_error_info "Init final run failed. Please check log.init_atmosphere.*"
            exit 1
        fi
    fi
}

function set_namelist_atmosphere() {
    # 首次运行atmosphere时进行的设置
    local configName configValue config
    set_namelist config_block_decomp_file_prefix \'${decompPrefix}\' ./namelist.atmosphere
    for config in ${nhyd_model_configs[*]}; do
        configName=${config%%=*} # 取第一个=前面的内容
        configValue=${config#*=} # 取第一个=后面的内容
        if [[ `sed -n "/${configName}/ p" ./namelist.atmosphere` == "" ]]; then
            echo "=> ${configName} = ${configValue}"
            sed -i "/&nhyd_model/ a\ ${configName} = ${configValue}" ./namelist.atmosphere
        else
            set_namelist ${configName} ${configValue} ./namelist.atmosphere
        fi
    done
    for config in ${physics_configs[*]}; do
        configName=${config%%=*} # 取第一个=前面的内容
        configValue=${config#*=} # 取第一个=后面的内容
        if [[ `sed -n "/${configName}/ p" ./namelist.atmosphere` == "" ]]; then
            echo "=> ${configName} = ${configValue}"
            sed -i "/&physics/ a\ ${configName} = ${configValue}" ./namelist.atmosphere
        else
            set_namelist ${configName} ${configValue} ./namelist.atmosphere
        fi
    done
    for config in ${damping[*]}; do
        configName=${config%%=*} # 取第一个=前面的内容
        configValue=${config#*=} # 取第一个=后面的内容
        if [[ `sed -n "/${configName}/ p" ./namelist.atmosphere` == "" ]]; then
            echo "=> ${configName} = ${configValue}"
            sed -i "/&damping/ a\ ${configName} = ${configValue}" ./namelist.atmosphere
        else
            set_namelist ${configName} ${configValue} ./namelist.atmosphere
        fi
    done
    for config in ${chemistry_configs[*]}; do
        configName=${config%%=*} # 取第一个=前面的内容
        configValue=${config#*=} # 取第一个=后面的内容
        if [[ `sed -n "/${configName}/ p" ./namelist.atmosphere` == "" ]]; then
            echo "=> ${configName} = ${configValue}"
            sed -i "/&chemistry/ a\ ${configName} = ${configValue}" ./namelist.atmosphere
        else
            set_namelist ${configName} ${configValue} ./namelist.atmosphere
        fi
    done
    for config in ${printout_configs[*]}; do
        configName=${config%%=*} # 取第一个=前面的内容
        configValue=${config#*=} # 取第一个=后面的内容
        if [[ `sed -n "/${configName}/ p" ./namelist.atmosphere` == "" ]]; then
            echo "=> ${configName} = ${configValue}"
            sed -i "/&printout/ a\ ${configName} = ${configValue}" ./namelist.atmosphere
        else
            set_namelist ${configName} ${configValue} ./namelist.atmosphere
        fi
    done
    for config in ${io_configs[*]}; do
        configName=${config%%=*} # 取第一个=前面的内容
        configValue=${config#*=} # 取第一个=后面的内容
        if [[ `sed -n "/${configName}/ p" ./namelist.atmosphere` == "" ]]; then
            echo "=> ${configName} = ${configValue}"
            sed -i "/&io/ a\ ${configName} = ${configValue}" ./namelist.atmosphere
        else
            set_namelist ${configName} ${configValue} ./namelist.atmosphere
        fi
    done
}

function set_streams_atmosphere() {
    local tempSurfaceInterval
    local tempEmissionInterval
    # 首次运行atmosphere时进行的设置
    # Edit streams.atmosphere
    sed -i "/stream name=\"input\"/, /filename_template/ s/filename_template.*$/filename_template=\"${initFinalFileName}\"/" ./streams.atmosphere
    sed -i "/stream name=\"input\"/ a\                  io_type=\"pnetcdf,cdf5\"" ./streams.atmosphere

    sed -i "/stream name=\"surface\"/, /filename_template/ s/filename_template.*$/filename_template=\"${surfaceFileName}\"/" ./streams.atmosphere
    sed -i "/stream name=\"surface\"/ a\        io_type=\"pnetcdf,cdf5\"" ./streams.atmosphere
    tempSurfaceInterval=`date -d "20200101 21600 seconds"  "+%H:%M:%S"` # 将seconds转换为小时
    sed -i "/stream name=\"surface\"/, /input_interval/ s/input_interval.*$/input_interval=\"${tempSurfaceInterval}\" >/" ./streams.atmosphere
    sed -i "/stream name=\"emission\"/, /filename_template/ s/filename_template.*$/filename_template=\"${emissionFileName}\"/" ./streams.atmosphere
    sed -i "/stream name=\"emission\"/ a\        io_type=\"pnetcdf,cdf5\"" ./streams.atmosphere
    tempEmissionInterval=`date -d "20200101 3600 seconds"  "+%H:%M:%S"` # 将seconds转换为小时
    sed -i "/stream name=\"emission\"/, /input_interval/ s/input_interval.*$/input_interval=\"${tempEmissionInterval}\" >/" ./streams.atmosphere
    sed -i "/stream name=\"restart\"/, /filename_template/ s/filename_template.*$/filename_template=\"iamas_rst.\$Y-\$M-\$D_\$h.\$m.\$s.nc\"/" ./streams.atmosphere
    sed -i "/stream name=\"restart\"/ a\                  clobber_mode=\"replace_files\"" ./streams.atmosphere
    sed -i "/stream name=\"restart\"/ a\                  io_type=\"pnetcdf,cdf5\"" ./streams.atmosphere
    sed -i "/stream name=\"restart\"/, /output_interval/ s/output_interval.*$/output_interval=\"10_00:00:00\" \/>/" ./streams.atmosphere

    sed -i "/stream name=\"output\"/, /filename_template/ s/filename_template.*$/filename_template=\"iamas_out.\$Y-\$M-\$D_\$h.\$m.\$s.nc\"/" ./streams.atmosphere
    sed -i "/stream name=\"output\"/ a\        clobber_mode=\"replace_files\"" ./streams.atmosphere
    sed -i "/stream name=\"output\"/ a\        io_type=\"pnetcdf,cdf5\"" ./streams.atmosphere
    sed -i "/stream name=\"output\"/, /output_interval/ s/output_interval.*$/output_interval=\"01:00:00\" >/" ./streams.atmosphere

    sed -i "/stream name=\"diagnostics\"/, /filename_template/ s/filename_template.*$/filename_template=\"iamas_diag.\$Y-\$M-\$D_\$h.\$m.\$s.nc\"/" ./streams.atmosphere
    sed -i "/stream name=\"diagnostics\"/ a\        clobber_mode=\"replace_files\"" ./streams.atmosphere
    sed -i "/stream name=\"diagnostics\"/ a\        io_type=\"pnetcdf,cdf5\"" ./streams.atmosphere
    sed -i "/stream name=\"diagnostics\"/, /output_interval/ s/output_interval.*$/output_interval=\"01:00:00\" >/" ./streams.atmosphere

    sed -i "/stream name=\"iau\"/, /\/>/ d" ./streams.atmosphere # 删除stream "iau"
    sed -i "/stream name=\"lbc_in\"/, /\/>/ d" ./streams.atmosphere # 删除stream "lbc_in"
    sed -i "/^\s*$/ d" ./streams.atmosphere # 删除空行
}

print_whole_line ${infoLength} "-"
echo "Be SURE to set up stream_list.atmosphere.* first!"
print_whole_line ${infoLength} "-"
# 首先准备surface file
echo ""

# Run emission, zhxia
if [[ "$doEmisUpdate" == true ]]; then
    currentTime=`date +%s`
    print_whole_line ${infoLength} "-" "Start running emission"
    run_emission
    costTime=$((`date +%s`-currentTime))
    echo "=> cost ${costTime} s"
    print_whole_line ${infoLength} "-"
    echo ""
fi

# Run surface
currentTime=`date +%s`
print_whole_line ${infoLength} "-" "Start running surface"
run_surface
costTime=$((`date +%s`-currentTime))
echo "=> cost ${costTime} s"
print_whole_line ${infoLength} "-"
echo ""

# 设置第一次运行的atmosphere
print_whole_line ${infoLength} "-" "Setting up atmosphere"
echo "=> setting namelist.atmosphere"
set_namelist_atmosphere
echo "=> setting streams.atmosphere"
set_streams_atmosphere
print_whole_line ${infoLength} "-"
# 备份第一次运行的log
mkdir -p ./log/atmosphere
cp ./streams.atmosphere ./log/atmosphere/streams.atmosphere

mkdir -p output

if [[ ${doBlkRun} == "true" ]]; then
    # blk run说明:
    # 0. 运行surface以及第一个init final, 如果第一个init需要chem_input则需要在./output中存在对应时刻的iamas_out
    # 1. 运行${blkDays}天
    # 2. 将生成的iamas_out以及iamas_diag文件移动到./output中
    # 3. 将上个run blk最后时刻提前${spinUpDays}天之前的iamas_out链接为${chemInputFileName}作为下一个blk的chem_input
    # 4. 基于上述${chemInputFileName}进行带chem_input的init final, 然后开始循环
    # 5. 如果doChemInput=false则不进行上述与chem_input相关的操作
    runFlag=yes
    blkStartTime=${startTimeStr}
    while [[ $runFlag == "yes" ]]
    do
        # 准备新一轮的init_final.nc
        echo ""
        print_whole_line ${infoLength} "-" "Initializing init final"
        currentTime=`date +%s`
        if [[ -e ./output/iamas_out.${blkStartTime//:/.}.nc ]]; then
            # 准备chem_input文件
            rm -f ./chem_input.nc
            ln -sf ./output/iamas_out.${blkStartTime//:/.}.nc ./chem_input.nc
            run_init_final ${blkStartTime} ${doChemInput}
        else
            run_init_final ${blkStartTime} 0
        fi
        costTime=$((`date +%s`-currentTime))
        echo "=> cost ${costTime} s"
        print_whole_line ${infoLength} "-"

        tempBlkEndTime=`date -d "${blkStartTime/_/ } ${blkDays} days" "+%Y-%m-%d_%H:%M:%S"`
        if [[ `date +%s -d "${tempBlkEndTime/_/ }"` -ge `date +%s  -d "${endTimeStr/_/ }"` ]]; then
            blkEndTime=${endTimeStr}
            runFlag=no
        else
            blkEndTime=${tempBlkEndTime}
        fi
        runDurationDays=$((($(date +%s -d "${blkEndTime/_/ }" ) - $(date +%s -d "${blkStartTime/_/ }"))/86400))
        runDurationModD=$((($(date +%s -d "${blkEndTime/_/ }" ) - $(date +%s -d "${blkStartTime/_/ }"))%86400))
        runDuration=${runDurationDays}_`date -d "20201231 ${runDurationModD} seconds" "+%H:%M:%S"`
        set_namelist config_start_time ${blkStartTime}  ./namelist.atmosphere > ${debugFile}
        set_namelist config_run_duration ${runDuration} ./namelist.atmosphere > ${debugFile}
        if [[ -e ./output/iamas_out.${blkStartTime//:/.}.nc ]]; then
            if [[ ${doChemInput} == 1 ]]; then
                set_namelist config_input_chem true ./namelist.atmosphere > ${debugFile}
            else
                set_namelist config_input_chem false ./namelist.atmosphere > ${debugFile}
            fi
        else
            set_namelist config_input_chem false ./namelist.atmosphere > ${debugFile}
        fi
        cp ./namelist.atmosphere ./log/atmosphere/namelist.atmosphere_${blkStartTime}-${blkEndTime}
        echo ""
        print_whole_line ${infoLength} "-" "Running block for atmosphere"
        if [[ ${startHour} -gt 0 ]]; then
            echo "=> Block start/end : `date -d "${blkStartTime/_/ }" "+%Y-%m-%d_%H:%M:%S"` / `date -d "${blkEndTime/_/ }" "+%Y-%m-%d_%H:%M:%S"`"
        else
            echo "=> Block start/end : `date -d "${blkStartTime/_/ }" "+%Y-%m-%d"` / `date -d "${blkEndTime/_/ }" "+%Y-%m-%d"`"
        fi
        echo "=> Block duration  : ${runDuration} days"
        currentTime=`date +%s`
        mpirun -n ${nProcs} ./atmosphere &
        mpi_pid=$!
        currRunDate="Initializing ..."
        > log.atmosphere.0000.out #删除原先的文件 revised by qkyang
        printf "=> Progress: %s\r" "${currRunDate}"
        #added by qkyang++++++++++++
        printf "\nqkyang: The estimated eta time was not accurate in the first about 5%% records!\n"
        SystemTimeOld=`date +%s` #读取计算机目前的时间戳（默认值，以防未能读取到第一步模拟时间戳）
        #added by qkyang-----------
        while [[ $(kill -0 ${mpi_pid} > /dev/null 2>&1;echo $?) == "0" ]]; do # 判断进程是否存在
            tmpRunDate=${currRunDate}
            currRunDate="$(sed -n -e "/Begin timestep/p" ./log.atmosphere.0000.out | tail -n 1 | awk '{print $3}')"
            if [[ "${currRunDate}" != "" && "${currRunDate}" != "${tmpRunDate}" ]]; then
                complete_sec=$(($(date +%s -d "${currRunDate/_/ }") - $(date +%s -d "${blkStartTime/_/ }")))
                total_sec=$(($(date +%s -d "${blkEndTime/_/ }") - $(date +%s -d "${blkStartTime/_/ }")))
                percentage=$(echo "scale=5;${complete_sec} / ${total_sec} * 100.0"|bc)
                #replaced by qkyang+++++++++++++++++
                if [ $complete_sec == 0 ];then
                    SystemTimeOld=`date +%s` #读取模拟第一步目前的时间戳
                    sleep 1 #等待1秒，避免计算的模拟耗时时间太短
                    continue
                fi
                SystemTimeNew=`date +%s` #读取计算机目前的时间戳
                CostTime=$(($SystemTimeNew-$SystemTimeOld))     #得到模拟耗费时间
                CostTimeHour=$[$CostTime/3600]                #还需模拟耗时小时
                let Intermediate1=($CostTime-3600*$CostTimeHour)   #中间量
                let CostTimeMinute=(Intermediate1/60) #还需模拟耗时分钟 #
                RestTime=$(($CostTime*${total_sec}/${complete_sec}-$CostTime)) #
                let RestTimeHour=($RestTime/3600)                #还需模拟耗时小时 #
                let Intermediate2=($RestTime-3600*$RestTimeHour)   #中间量 #
                RestTimeMinute=$((Intermediate2/60)) #还需模拟耗时分钟 #
                printf "Clock: %s | Progress: %s [%6.2f%%] | Elapsed: %s | Eta: %s \n"\
                "$(date "+%Y-%m-%d_%H:%M:%S")" "${currRunDate}" "${percentage}" \
                "${CostTimeHour} h ${CostTimeMinute} m" \
                "${RestTimeHour} h ${RestTimeMinute} m"
                #replaced by qkyang-----------------
            fi
            sleep 1
        done
        echo ""
        wait
        if [[ ! -e ./iamas_diag.${blkEndTime//:/.}.nc ]]; then
            print_error_info "Run atmosphere failed. Please check!"
            exit 1
        fi
        # 删除blkStartTime开始${spinUpDays}天内的iamas_out*和iamas_diag*
        echo "=> Deleting iamas files in spin-up days ..."
        for iamas_diag in `ls iamas_diag*.nc`; do
            tempTimeStr=${iamas_diag#*.} # 取第一个.后面的内容
            tempTimeStr=${tempTimeStr%.*} # 取最后一个.前面的内容
            tempTimeStr=${tempTimeStr//./:} # 输出文件中格式为iamas_diag.2010-01-01_00.00.00.nc, 转换为timeStr格式
            if [[ `date +%s -d "${tempTimeStr/_/ }"` -lt `date +%s -d "${blkStartTime/_/ } ${spinUpDays} days"` ]]; then
                rm ${iamas_diag}
            fi
        done
        for iamas_out in `ls iamas_out*.nc`; do
            tempTimeStr=${iamas_out#*.} # 取第一个.后面的内容
            tempTimeStr=${tempTimeStr%.*} # 取最后一个.前面的内容
            tempTimeStr=${tempTimeStr//./:} # 输出文件中格式为iamas_diag.2010-01-01_00.00.00.nc, 转换为timeStr格式
            if [[ `date +%s -d "${tempTimeStr/_/ }"` -lt `date +%s -d "${blkStartTime/_/ } ${spinUpDays} days"` ]]; then
                rm ${iamas_out}
            fi
        done
        # 移动至./output/
        echo "=> Moving iamas files to ./output ..."
        mv ./init_final.nc ./iamas_out*.nc ./iamas_diag*.nc ./output/
        cp ./log.atmosphere.0000.out ./log/atmosphere/log.atmosphere.0000.out_${startTimeStr}-${endTimeStr}
        # 设置下一轮blk开始时间
        blkStartTime=`date -d "${blkEndTime/_/ } ${spinUpDays} days ago" "+%Y-%m-%d_%H:%M:%S"`
        costTime=$((`date +%s`-currentTime))
        echo "=>  cost ${costTime} s"
        print_whole_line ${infoLength} "-"
    done
elif [[ ${doBlkRun} == "false" ]]; then
    # 准备init_final.nc
    echo ""
    print_whole_line ${infoLength} "-" "Initializing init final"
    currentTime=`date +%s`
    if [[ -e ./output/iamas_out.${startTimeStr//:/.}.nc ]]; then
        # 准备chem_input文件
        rm -f ./chem_input.nc
        ln -sf ./output/iamas_out.${startTimeStr//:/.}.nc ./chem_input.nc
        run_init_final ${startTimeStr} ${doChemInput}
    elif [[ -e ./chem_input.nc ]]; then
        run_init_final ${startTimeStr} ${doChemInput}
    else
        run_init_final ${startTimeStr} 0
    fi
    costTime=$((`date +%s`-currentTime))
    echo "=> cost ${costTime} s"
    print_whole_line ${infoLength} "-"

    runDurationDays=$((($(date +%s -d "${endTimeStr/_/ }" ) - $(date +%s -d "${startTimeStr/_/ }"))/86400))
    runDurationModD=$((($(date +%s -d "${endTimeStr/_/ }" ) - $(date +%s -d "${startTimeStr/_/ }"))%86400))
    runDuration=${runDurationDays}_`date -d "20201231 ${runDurationModD} seconds" "+%H:%M:%S"`
    set_namelist config_start_time ${startTimeStr}  ./namelist.atmosphere > ${debugFile}
    set_namelist config_run_duration ${runDuration} ./namelist.atmosphere > ${debugFile}
    if [[ -e ./output/iamas_out.${blkStartTime//:/.}.nc ]]; then
        if [[ ${doChemInput} == 1 ]]; then
            set_namelist config_input_chem true ./namelist.atmosphere > ${debugFile}
        else
            set_namelist config_input_chem false ./namelist.atmosphere > ${debugFile}
        fi
    else
        set_namelist config_input_chem false ./namelist.atmosphere > ${debugFile}
    fi
    cp ./namelist.atmosphere ./log/atmosphere/namelist.atmosphere_${startTimeStr}-${endTimeStr}
    echo ""
    print_whole_line ${infoLength} "-" "Running continue for atmosphere"
    echo "=> Block start/end : ${startTimeStr} / ${endTimeStr}"
    echo "=> Block duration  : ${runDuration} days"

    currentTime=`date +%s`
    mpirun -n ${nProcs} ./atmosphere &
    mpi_pid=$!
    currRunDate="Initializing ..."
    > log.atmosphere.0000.out #删除原先的文件 revised by qkyang
    printf "=> Progress: %s\r" "${currRunDate}"
    #added by qkyang++++++++++++
    printf "\nqkyang: The estimated eta time was not accurate in the first about 5%% records!\n"
    SystemTimeOld=`date +%s` #读取计算机目前的时间戳（默认值，以防未能读取到第一步模拟时间戳） 
    #added by qkyang-----------
    while [[ $(kill -0 ${mpi_pid} > /dev/null 2>&1;echo $?) == "0" ]]; do # 判断进程是否存在
        tmpRunDate=${currRunDate}
        currRunDate="$(sed -n -e "/Begin timestep/p" ./log.atmosphere.0000.out | tail -n 1 | awk '{print $3}')"
        if [[ "${currRunDate}" != "" && "${currRunDate}" != "${tmpRunDate}" ]]; then
            complete_sec=$(($(date +%s -d "${currRunDate/_/ }") + 600 - $(date +%s -d "${startTimeStr/_/ }")))
            total_sec=$(($(date +%s -d "${endTimeStr/_/ }") - $(date +%s -d "${startTimeStr/_/ }")))
            percentage=$(echo "scale=5;${complete_sec} / ${total_sec} * 100.0"|bc)
            #replaced by qkyang+++++++++++++++++
            if [ $complete_sec == 0 ];then
                SystemTimeOld=`date +%s` #读取模拟第一步目前的时间戳
                sleep 1 #等待1秒，避免计算的模拟耗时时间太短
                continue
            fi
            SystemTimeNew=`date +%s` #读取计算机目前的时间戳
            CostTime=$(($SystemTimeNew-$SystemTimeOld))     #得到模拟耗费时间
            CostTimeHour=$[$CostTime/3600]                #还需模拟耗时小时
            let Intermediate1=($CostTime-3600*$CostTimeHour)   #中间量
            let CostTimeMinute=(Intermediate1/60) #还需模拟耗时分钟 #
            RestTime=$(($CostTime*${total_sec}/${complete_sec}-$CostTime)) #
            let RestTimeHour=($RestTime/3600)                #还需模拟耗时小时 #
            let Intermediate2=($RestTime-3600*$RestTimeHour)   #中间量 #
            RestTimeMinute=$((Intermediate2/60)) #还需模拟耗时分钟 #
            printf "Clock: %s | Progress: %s [%6.2f%%] | Elapsed: %s | Eta: %s \n"\
            "$(date "+%Y-%m-%d_%H:%M:%S")" "${currRunDate}" "${percentage}" \
            "${CostTimeHour} h ${CostTimeMinute} m" \
            "${RestTimeHour} h ${RestTimeMinute} m"
            #replaced by qkyang-----------------
        fi
        sleep 1
    done
    echo ""
    wait
    if [[ ! -e ./iamas_diag.${endTimeStr//:/.}.nc ]]; then
        print_error_info "Run atmosphere failed. Please check."
        exit 1
    fi
    # 删除${spinUpDays}天内的iamas_out*和iamas_diag*
    echo "=> Deleting iamas files in spin-up days ..."
    for iamas_diag in `ls iamas_diag*.nc`; do
        tempTimeStr=${iamas_diag#*.} # 取第一个.后面的内容
        tempTimeStr=${tempTimeStr%.*} # 取最后一个.前面的内容
        tempTimeStr=${tempTimeStr//./:} # 输出文件中格式为iamas_diag.2010-01-01_00.00.00.nc, 转换为timeStr格式
        if [[ `date +%s -d "${tempTimeStr/_/ }"` -lt `date +%s -d "${startTimeStr/_/ } ${spinUpDays} days"` ]]; then
            rm ${iamas_diag}
        fi
    done
    for iamas_out in `ls iamas_out*.nc`; do
        tempTimeStr=${iamas_out#*.} # 取第一个.后面的内容
        tempTimeStr=${tempTimeStr%.*} # 取最后一个.前面的内容
        tempTimeStr=${tempTimeStr//./:} # 输出文件中格式为iamas_diag.2010-01-01_00.00.00.nc, 转换为timeStr格式
        if [[ `date +%s -d "${tempTimeStr/_/ }"` -lt `date +%s -d "${startTimeStr/_/ } ${spinUpDays} days"` ]]; then
            rm ${iamas_out}
        fi
    done
    # 移动至./output/
    echo "=> Moving iamas files to ./output ..."
    mv ./init_final.nc ./iamas_out*.nc ./iamas_diag*.nc ./output/
    cp ./log.atmosphere.0000.out ./log/atmosphere/log.atmosphere.0000.out_${startTimeStr}-${endTimeStr}
    costTime=$((`date +%s`-currentTime))
    echo "=> cost ${costTime} s"
    print_whole_line ${infoLength} "-"

else
    print_error_info "\${doBlkRun} must be true/false. Please check."
    exit 1
fi

costTime=$((`date +%s` - simStartTime))
simDurationDays=$((costTime/86400))
simDurationModD=$((costTime%86400))
simDuration=${simDurationDays}_`date -d "20201231 ${simDurationModD} seconds" "+%H:%M:%S"`
echo "Simulation Complete at `date "+%Y-%m-%d %H:%M:%S"`"
echo "Total cost: ${simDuration}"
cp -f ./run_atmosphere.sh ./output/
exit 0

#附件代码(不会执行，非正式运行所需，自行复制到终端执行) added by qkyang+++++++++

#(1)移动模拟结果到storage
#先cd 到当前目录，例如cd /lustre/qkyang/work/Run_iamas/cases/U120km_Urban1_USGS_201903
OutputStorageDir="/storage/aemolcore02/qkyang/iAMAS_Simulation_Results" #added by qkyang, 全部模拟结束后output模拟结果保存目录（会自动根据当前sh脚本路径名创建保存的文件夹名）
TargetStorageDir=$OutputStorageDir/$(echo $(eval pwd) | sed 's/\//+/g')"+output"
echo "Start move outputs to "$TargetStorageDir" at ""`date "+%Y-%m-%d %H:%M:%S"`"
mkdir $TargetStorageDir
mv ./output/init_final.nc ./output/iamas_out*.nc ./output/iamas_diag*.nc $TargetStorageDir
echo "Finished move outputs to $TargetStorageDir at `date "+%Y-%m-%d %H:%M:%S"`"

#(2)清空run_atmosphere.sh所有当前目录生成结果 "删除文件操作，谨慎执行!!!"
#首先进入指定目录，例如：cd /lustre/qkyang/work/Run_iamas/cases/V16km_201903_MODIS_Urban1_SmallZRTest
#find . -maxdepth 1 ! -name '.' ! -name 'run_atmosphere.sh' ! -name 'stream_list.atmosphere.*' -exec rm -rfv {} +

#附件代码(不会执行，非正式运行所需，自行复制到终端执行) added by qkyang----------

