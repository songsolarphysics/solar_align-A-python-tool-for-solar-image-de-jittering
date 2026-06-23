# solar_align

[English](README.md) | **[中文]**

> 一个用于太阳图像时序数据去抖动 (de-jittering) 的 Python 工具包。
> 移植自蔡钟哲 (J. Chae) 教授的 IDL 程序 `CH_ALIGNOFFSET` 与 `CH_SHIFT_SUB`，
> 在非循环平移、多尺寸输入兼容、偏移量汇总等方面做了现代化改造。
> 已在地基 NVST/Hα 与太空 Solar Orbiter EUI/HRI 数据上实测。

![python](https://img.shields.io/badge/python-%E2%89%A53.8-blue)
![license](https://img.shields.io/badge/license-MIT-green)
![status](https://img.shields.io/badge/status-tested%20on%20NVST%20%2B%20EUI/HRI-success)

---

## 这个工具在做什么

给定一个目录里的若干 FITS 帧——它们由于大气湍流残余、镜筒振动、卫星平台指向抖动等
原因，相邻帧之间存在亚像素到几像素级别的随机偏移——`solar_align` 会用 FFT 互相关
加亚像素精修的方法，把所有帧对齐到一个共同参考帧，然后输出：

- 对齐后的 FITS 文件（文件名不变、形状不变、不循环回绕）
- 每帧偏移量汇总文件（`offsets_summary.csv` 和 `.npz`）

算法核心（FFT 互相关 + 3×3 抛物线亚像素拟合 + 高斯窗 + 多轮迭代精修）完全沿用
Chae 教授 1999 年前后编写的 IDL 实现——这套算法在太阳物理社区里跑了快三十年。

## 主要特性

- **亚像素精度**：单轮约 0.05–0.15 像素，双轮迭代后约 0.03 像素
- **非循环平移**：越界像素填 `NaN`（或用户指定值），不绕回
- **单次插值保证**：多轮迭代时把累积偏移量应用到**原始数据**上一次性平移，避免反复插值导致的模糊
- **自动处理多尺寸输入**：以最小公共尺寸做 center crop，含安全阈值
- **保留 sunpy 元数据**：用 sunpy.map 读写，WCS 与观测时间信息完整保留
- **纯 Python**：无需编译，pip 一行装好

## 安装

```bash
pip install scipy astropy pandas sunpy
```

然后克隆本仓库，或者把 `solar_align.py` 直接放到你的脚本/notebook 旁边。

## 五分钟上手

```python
from solar_align import align_fits_directory

df = align_fits_directory(
    input_dir  = '/path/to/raw_fits',
    output_dir = '/path/to/aligned_fits',
    pattern    = '*.fits',
)
```

运行后会有：

```
/path/to/aligned_fits/
    ├── frame_001.fits          ← 对齐后，文件名不变，形状不变
    ├── frame_002.fits
    ├── ...
    ├── offsets_summary.csv     ← 偏移量汇总，人读
    └── offsets_summary.npz     ← 偏移量汇总，程序读
```

`input_dir` 里按字母序排第一个文件被当作固定参考帧，所有后续帧都对齐到它。

### 检查结果

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('/path/to/aligned_fits/offsets_summary.csv')
fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(df['x_offset'], '.-', label='dx')
ax.plot(df['y_offset'], '.-', label='dy')
ax.set_xlabel('Frame'); ax.set_ylabel('Offset (px)')
ax.legend()
```

CSV 各列含义：

| 列名 | 含义 |
|---|---|
| `filename` | 输出文件名 |
| `time` | 从 header 读取的观测时间 |
| `x_offset`, `y_offset` | 测得的偏移量（像素） |
| `correlation` | 对齐后的归一化互相关系数（1.0 = 完美对齐） |
| `n_iter` | 触发了几轮精修 |

## 仓库结构

```
solar_align/
├── solar_align.py          ← 核心模块
├── fits_to_png.py          ← 把对齐后的 FITS 渲染成 PNG，方便检查
├── make_compare_gif.py     ← 把对齐前后两组 PNG 拼成对比动画 GIF
├── example_usage.py        ← 完整调用示例
├── test_solar_align.py     ← 单元测试 + 集成测试（基于合成数据）
├── README.md               ← 英文版
└── README_zh.md            ← 本文件
```

## API 速览

| 函数 | 作用 |
|---|---|
| `align_fits_directory(input_dir, output_dir, ...)` | 对整个目录批量对齐 |
| `ch_alignoffset(image, reference)` | 测量两幅图像之间的偏移 |
| `ch_shift_sub(image, x0, y0)` | 亚像素平移（非循环） |
| `scan_fits_shapes(input_dir)` | 扫描目录里所有 FITS 的 (ny, nx) 分布 |

完整签名和参数说明见源代码 docstring。

## 已测试过的数据

| 数据集 | 波段 | 仪器位置 | 备注 |
|---|---|---|---|
| NVST | Hα 6562.8 Å | 地基（中国） | 持续性 seeing 抖动 |
| Solar Orbiter EUI/HRI | EUV 174 Å | 太空 | 1–3 秒高 cadence，平台抖动 |
| 合成数据 | — | — | 量化恢复精度用 |

合成数据上的恢复精度：单轮 `≤ 0.15 px`，双轮 `≤ 0.05 px`（30 帧序列，输入 jitter rms ~1 px）。

## 何时它会失灵——失败模式

固定参考 + FFT 互相关的策略有几个已知边界：

- **序列里太阳真的在演化**（相关系数缓慢下降）→ 分段对齐
- **抖动太大**（>50 像素）→ 抛物线拟合不稳定，先做一次粗对齐
- **低信噪比**（cor < 0.7）→ 先做空间平滑或时间平均
- **视场内有显著动态特征**（如耀斑）→ 用 mask（功能开发中）或限制对齐子区域
- **参考帧不好**（seeing 极差、云层）→ 手动挑选优质参考帧

完整的失败模式分析与缓解策略，见我的[公众号文章](#)<sup>※</sup>。

<sup>※</sup> 链接地址：你需要把这里替换成实际公众号文章的 URL

## 把偏移量当作物理信号

`offsets_summary.csv` 里那两列偏移量本身就是物理数据。常见的"二次分析"用法：

- **抖动功率谱**：把 `x_offset` 做 FFT，识别周期性扰动源
- **长期漂移趋势**：低通滤波后判断是热漂移还是机械漂移
- **多仪器联合**：两台仪器同期的偏移序列做互相关，追溯共同扰动事件
- **配准质量加权**：`correlation` 列可以直接当作时序分析的权重

## 致敬

算法核心——FFT 互相关 + 3×3 抛物线亚像素拟合 + 高斯窗 + 多轮迭代精修——
**完全来自蔡钟哲教授 (J. Chae, 首尔大学) 的 IDL 程序**
`CH_ALIGNOFFSET.pro` 和 `CH_SHIFT_SUB.pro`。这两段代码最早写于 1999 年前后，
被打包进 FISS 数据处理流水线，至今仍在太阳物理社区里被广泛使用。本 Python
版本是对这套算法的重新实现和现代化包装，方法学的全部功劳归 Chae 教授。

如果你在已发表的工作中使用了本工具，请同时致谢 Python 移植版本和原始 IDL 程序。

## 许可证

MIT License，详见 `LICENSE`。

## 联系方式

欢迎通过 GitHub Issues 提交 bug 报告、功能建议或 pull request。

作者：谭宋 (Song Tan)，Leibniz Institute for Astrophysics Potsdam (AIP)
联系方式：stan@aip.de
