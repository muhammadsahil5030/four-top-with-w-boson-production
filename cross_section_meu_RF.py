import matplotlib.pyplot as plt
import numpy as np

# LHAPDF sets and corresponding cross section values
lhapdf_sets = ["NNPDF30-nlo-as-0115", "NNPDF30-nlo-as-0117", "NNPDF30-nlo-as-0118", 
               "NNPDF30-nlo-as-0119", "NNPDF30-nlo-as-0121"]
cross_sections = np.array([5.96, 6.44, 6.62, 6.85, 7.26])

# Systematic errors (low and high)
sys_errors_low = np.array([0.020, 0.025, 0.024, 0.026, 0.025])
sys_errors_high = np.array([0.020, 0.025, 0.024, 0.026, 0.025])

# Statistical errors (low and high)
stat_errors_low = np.array([1.99, 1.76, 1.82, 1.88, 2.11])
stat_errors_high = np.array([1.35, 1.53, 1.58, 1.66, 1.72])

plt.figure(figsize=(8, 6))

for i, (cs, stat_low, stat_high, sys_low, sys_high) in enumerate(zip(cross_sections,
                                                                     stat_errors_low, stat_errors_high,
                                                                     sys_errors_low, sys_errors_high)):
    plt.errorbar(cs, i, xerr=[[stat_low], [stat_high]], fmt='o', color='red', ecolor='red', elinewidth=1, capsize=0, label='$\mu_{r,f}$ uncertainty' if i == 0 else "")
    plt.errorbar(cs, i, xerr=[[sys_low], [sys_high]], fmt='o', color='black', label='Statistical uncertainty' if i == 0 else "")


# Adding custom text in the upper left corner
plt.text(1.3, 4.4, r"$t\bar{t}t\bar{t}W$ Signal", fontsize=11, verticalalignment='top')
plt.text(14.0, 3.4, r"$\sigma = 6.62_{-2.6}^{+2.4}$ ab", fontsize=11, verticalalignment='top')
plt.text(14.0, 4.4, r"$\sqrt{s}=$ 13TeV", fontsize=11, verticalalignment='top')
plt.text(14.0, 3.0, r"$\mu_{r,f} = 692$", fontsize=12, verticalalignment='top')
plt.text(14.0, 2.6, r"$\mu_{r,f} = 346$", fontsize=12, verticalalignment='top')
plt.text(14.0, 2.2, r"$mu_{r,f} = 173$", fontsize=12, verticalalignment='top')

# Decrease font size and tilt y-axis labels
plt.yticks(range(len(lhapdf_sets)), lhapdf_sets, rotation=45, ha='right', fontsize=7)  # Adjust fontsize as needed

plt.xticks(np.arange(1.0, 18.0, 1))
#plt.xticks(rotation=0)

x_start = 4.0  
x_end = 9.0    
plt.axvspan(x_start, x_end, color='green', alpha=0.3)
plt.xlabel("$\sigma(ab)$", fontsize=15)
plt.xlim(1.0, 18.0)
plt.title("NLO Cross-Section", fontsize = 17)
plt.legend()
#plt.grid(True)
plt.savefig("/home/msahil/work/analysis/fourtop_wminus/presentation_1/presentation_images/cross_section_plots.png")
plt.show()

