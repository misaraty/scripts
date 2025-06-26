import os
import matplotlib.pyplot as plt
import matplotlib as mpl
os.chdir(os.path.split(os.path.realpath(__file__))[0])
# mpl.rcParams['font.family'] = 'Arial'

# 从 epw2.out 中提取最大本征值数据
with open('epw2.out', 'r') as f:
    lines = f.readlines()

data = []
capture = False
for i, line in enumerate(lines):
    if "Max. eigenvalue" in line:
        capture = True
        start = i + 1
    if capture and len(data) < 21 and i >= start:
        try:
            parts = line.strip().split()
            if len(parts) >= 2:
                T = float(parts[0])
                eig = float(parts[1])
                data.append((T, eig))
        except:
            pass

# 保存数据文件
with open('data_max_eigenvalue.dat', 'w') as f:
    for T, eig in data:
        f.write(f"{T:.2f} {eig:.6f}\n")

# 分离温度与本征值
temps, eigs = zip(*data)

# 开始绘图
plt.figure()

# 样式设置
plt.plot(temps, eigs, marker='o', markersize=5, linewidth=2, label='Max. eigenvalue', color='tab:blue')
plt.axhline(1.0, color='black', linestyle='--', linewidth=1)  # 1.0 本征值水平线
plt.axvline(4.61, color='gray', linestyle='--', linewidth=1)  # Tc 指示线
plt.text(4.3, 1.2, 'T$_c$ = 4.61 K', fontsize=16)
# plt.annotate('T$_c$ = 4.61 K',
             # xy=(4.61, 1.0),
             # xytext=(4.0, 1.8),
             # fontsize=18,
             # arrowprops=dict(arrowstyle='->', color='black', lw=1.5))

# 坐标轴和标签
plt.xlabel('Temperature (K)', fontsize=16)
plt.ylabel('Max. eigenvalue', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(0, 6)
plt.ylim(0, 5)
# plt.legend(fontsize=16, loc='upper right')

plt.tight_layout()
plt.savefig('fig.jpg', dpi=300)
# plt.show()
plt.close()
