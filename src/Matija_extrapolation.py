import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv('extrapolation_examples.csv')


xaxis =  1/ df['L']
xaxis1 = xaxis[0:7]
t = df['t']
yaxis3= df['E_p']
yaxis3_err = df['E_p_err']
yaxis3 = yaxis3[0:7]
yaxis3_err = yaxis3_err[0:7]


yaxis = df['E_fit'] + df['factor'] / df['L']
yaxis_error = df['E_fit_err'] + df['factor_err'] / df['L']
yaxis1 = yaxis[0:7]

yaxis_error1 = yaxis_error[0:7]

fig, ax = plt.subplots()
ax.errorbar(xaxis1, yaxis3, yerr=yaxis3_err, fmt='o')
ax.plot(xaxis1, yaxis1, linestyle = '--')


ax.set_xlabel(r'$1/L$', fontsize=14)
ax.set_ylabel(r'$E_P/U$', fontsize=14)

plt.show()