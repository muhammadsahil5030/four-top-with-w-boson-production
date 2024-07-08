import numpy as np
import matplotlib.pyplot as plt

num_var = [3, 6, 9, 12, 15, 18]

test_auc = [0.536, 0.546, 0.690, 0.862, 0.869, 0.904]
train_auc = [0.624, 0.685, 0.761, 0.899, 0.908, 0.932]

plt.ylim(0.50, 0.95)
plt.plot(num_var, test_auc, label='Test AUC')
plt.plot(num_var, train_auc, label= 'Train AUC')
plt.title("Variables Significane")
plt.xlabel("Number of variables")
plt.ylabel("AUC TTTTW")
plt.legend()
plt.savefig("bdt_variables_codes/tmva_bdt_variables_plots/variables_significane.png")
plt.show()
