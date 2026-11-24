# 当前代码版本

`main` 保存最近完成 12 小时模拟的 SW64 源码组合；`baseline` 保存原始 iAMAS v2.6.5 源码。

| 版本 | 分支 | 固定标签 | 内容 |
|---|---|---|---|
| 原始源码基线 | `baseline` | `v2.6.5-original-baseline` | 原始发布包，未加入 SW64 兼容与优化修改 |
| 最近 12 小时测试版 | `main` | `sw64-12h-9737442` | SW64 兼容、precomm 优化、vertmx Kokkos、模块细分计时 |

[查看两版代码差异](https://github.com/Keke0805/iamas/compare/baseline...main) · [版本与来源说明](docs/versions/版本说明.md) · [12 小时热点表](docs/12h-9737442/hotspot-table.md)

## 当前版的主要修改

- precomm 共享计算：合并散度与 limiter 过程、减少中间数组反复读写、每 4 个物种复用邻接与几何信息、对部分独立垂直循环使用 SIMD。
- vertmx：增加 Fortran 桥接与 Kokkos C++ 批量计算路径，保留原 Fortran 路径及校验模式。
- 物理、化学、输送和微物理过程：补充分阶段耗时统计。
- 保留 SW64 编译兼容修改和服务器源码中的历史备份文件。

## 最近一次长测

- job：`9737442`；2026-11-23 提交，2026-11-24 完成。
- 360 MPE；55 个垂直层；480 步，每步 90 秒，模拟 12 小时。
- precomm mode 3、vertmx mode 1；没有垂直填充。
- 总运行 max：**10565.3035 秒，约 2 小时 56 分 5 秒**。
- 当前 Kokkos 后端仍为 **Serial/MPE**；这些代码不能表述为已完成 CPE 加速。

源码同时包含多个模式。复现这次运行所用的选择是：

```bash
export IAMAS_PRECOMM_MODE=3
export IAMAS_VERTMX_MODE=1
export OMP_NUM_THREADS=1
export SWCH_USE_HW_COLL=0
```

环境变量只是运行选择；编译依赖和提交资源配置仍需按服务器环境准备。

## 代码入口

- [precomm 模式与优化实现](src/atmosphere/dynamics/iamas_atm_time_integration_scalars_mono_precomm.cpp)：`precomm_simd` 为 mode 3。
- [vertmx C++ 内核](src/atmosphere/chemistry/schemes/module_vertmx_kokkos.cpp)。
- [vertmx Fortran 桥接](src/atmosphere/chemistry/schemes/module_vertmx_kokkos.F)。
- [时间积分与输送调用](src/atmosphere/dynamics/iamas_atm_time_integration.F)。
- [干沉降调用与细分计时](src/atmosphere/chemistry/iamas_atmchem_driver_drydep.F)。

## 来源与复现范围

当前源码按“SW64 precomm 基础树 → vertmx 修改 → 模块计时修改”的顺序组合，对应服务器采用增量对象覆盖链接的测试版本。源文件校验和见 [12h-source-manifest.json](docs/versions/12h-source-manifest.json)；原始基线见 [baseline-source-manifest.json](docs/versions/baseline-source-manifest.json)。

[服务器构建脚本与配置记录](docs/server-build-evidence/)保留用于追溯，其中的绝对路径反映原测试环境，复用时需要调整。本次上传未进行全量重新编译，不将合并后的源码目录声称为已经重新构建验证的发行包。

性能表中的“优化前”是已完成 SW64 适配的运行版本，**不是** `baseline` 分支未经适配的原始发布源码；二者不能混为同一个性能基线。

原项目介绍和安装说明保留于 [README.md](README.md)。
