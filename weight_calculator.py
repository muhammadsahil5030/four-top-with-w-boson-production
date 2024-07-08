#this code is used to simply calculate the weight of a signal and backgrounds
# luminosity = events_number/crossx

import math

luminosity = 300				#inverse fb
events_number = 10000

#--------------------------------- luminosity ----------------------------------

signal_crossx = 0.006628	#fb
lumi_signal = events_number/signal_crossx	#inverse fb

ttbar_crossx = 456000		#fb
lumi_ttbar = events_number/ttbar_crossx	#inverse fb

ttbarH_crossx = 156.0		#fb
lumi_ttbarH = events_number/ttbarH_crossx	#inverse fb

ttbarW_crossx = 33.74		#fb
lumi_ttbarW = events_number/ttbarW_crossx	#inverse fb

ttbarZ_crossx = 128.8		#fb
lumi_ttbarZ = events_number/ttbarZ_crossx	#inverse fb

tripletopbbar_crossx = 0.02054	#fb
lumi_tripletopbbar = events_number/tripletopbbar_crossx	#inverse fb

tripletopwm_crossx = 0.1926	#fb
lumi_tripletopwm = events_number/tripletopwm_crossx	#inverse fb

tripletopwp_crossx = 0.1917	#fb
lumi_tripletopwp = events_number/tripletopwp_crossx	#inverse fb	

#------------------calculating weight of signal and background---------------------

signal_weight = luminosity/lumi_signal
ttbar_weight = luminosity/lumi_ttbar
ttbarH_weight = luminosity/lumi_ttbarH
ttbarW_weight = luminosity/lumi_ttbarW
ttbarZ_weight = luminosity/lumi_ttbarZ
tripletopbbar_weight = luminosity/lumi_tripletopbbar
tripletopwm_weight = luminosity/lumi_tripletopwm
tripletopwp_weight = luminosity/lumi_tripletopwp

"""
signal_weight = (300*0.0004482)/10000
ttbar_weight = (300*202000)/10000
ttbarH_weight = (300*156)/10000
ttbarW_weight = (300*33.74)/10000
ttbarZ_weight = (300*128.8)/10000
tripletop_weight = (300*0.02054)/10000
"""
print("luminosity of signal = ", lumi_signal, ",", "	signal weight = ",signal_weight)
print("luminosity of ttbar background = ", lumi_ttbar, ",", "  ttbar weight = ",ttbar_weight)
print("luminosity of ttbarH background = ", lumi_ttbarH, ",", "  ttbarH weight = ",ttbarH_weight)
print("luminosity of ttbarW background = ", lumi_ttbarW, ",", "  ttbarW weight = ",ttbarW_weight)
print("luminosity of ttbarZ background = ", lumi_ttbarZ, ",", "  ttbarZ weight = ",ttbarZ_weight)
print("luminosity of tripletopbbar background = ", lumi_tripletopbbar, ",", "  tripletop weight = ",tripletopbbar_weight)
print("luminosity of tripletopwm background = ", lumi_tripletopwm, ",", "  tripletop weight = ",tripletopwm_weight)
print("luminosity of tripletopwp background = ", lumi_tripletopwp, ",", "  tripletop weight = ",tripletopwp_weight)

#----------------------------------------END---------------------------------------
