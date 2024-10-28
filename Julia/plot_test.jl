using PyPlot

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.sans-serif"] = "Arial"
# plt.rcParams['font.sans-serif'] = ['Arial']
fontsize1=16
fontsize2=18

x = range(0; stop=2*pi, length=1000)
y = sin.(3 * x + 4 * cos.(2 * x));
plot(x, y, color="tab:blue", label=L"$V_\mathregular{Pb}^{0}$ + e$^{-}$ + h$^{+}$")
plt.xticks(fontsize=fontsize1)
plt.yticks(fontsize=fontsize1)
plt.xlabel(L"$Q$ (amu$^{1/2}$ $\mathregular{\AA}$)", fontsize=fontsize2)
plt.ylabel("Energy (eV)", fontsize=fontsize2)
# plt.axvline(0, color="tab:red", linestyle="dashed") #axvline axhline
# plt.axvline(3, color="tab:red", linestyle="dashed") #axvline axhline
annotate("", xy=(0,0), xytext=(3,0), xycoords="data", arrowprops=Dict("arrowstyle"=>"<->", "linestyle"=>"--", "linewidth"=>1, "shrinkA"=>0, "shrinkB"=>0, "color"=>"tab:red"))
plt.axvspan(xmin=0, xmax=3, facecolor="tab:red", alpha=0.2)
grid("on")
xlim(-1, 7)
ylim(-1.1, 1.1)
# "upper left', "upper right', "lower left', "lower right"
legend(loc="upper right", fontsize=fontsize1)
tight_layout()
savefig("/Users/lenovo/Desktop/test.jpg", dpi=600)
show()





