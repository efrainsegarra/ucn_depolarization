import numpy as np
import matplotlib.pyplot as plt

# Reduced data from Efrain (email 15/09/2024)
# Rexolite ring G10 scan RUN ??
G10 = [-300, -250, -200, -150, -100, -50, -25, -10, -5, -2, 0, 2, 5, 10, 25, 50, 100, 150, 200, 250, 300]
aBOT10 = [0.469, 0.545, 0.631, 0.715, 0.795, 0.859, 0.882, 0.895, 0.906, 0.907, 0.902, 0.897, 0.915, 0.903, 0.888, 0.858, 0.801, 0.721, 0.645, 0.568, 0.473]
aTOP10 = [0.573, 0.650, 0.721, 0.761, 0.797, 0.814, 0.812, 0.814, 0.816, 0.817, 0.812, 0.803, 0.797, 0.813, 0.810, 0.804, 0.786, 0.758, 0.713, 0.657, 0.586]
# scaling to a(G=0) = 1
aBOT10 = [a/0.905 for a in aBOT10]
aTOP10 = [a/0.817 for a in aTOP10]

# the C.M. offset for a UCN energy group as a function of energy in neV
H = 12.
def z_offset(E):
  E /= 1.02
  if E < H:
    return 2./5.*E - H/2.
  else:
    X = (1.-H/E)**(3./2.)
    return E/(1.-X)*(2./5.-X+3./5.*X**(5./3.)) - H/2.


gamma = 29e-6*2*np.pi
T = 180
a = 3.3
b = 4.7

def alpha(G, E_MIN, E_MAX):
  E_bin = np.arange(E_MIN,E_MAX+0.1, 0.1)
  z_bin = [z_offset(e) for e in E_bin]
  W_bin = [pow(e-E_MIN,a)*pow(E_MAX-e,b) for e in E_bin]
  mean_z = sum([Z*W for Z,W in zip(z_bin,W_bin)])/sum(W_bin)
  print(mean_z)
  s = 0
  for z,W in zip(z_bin,W_bin):
    s += W*np.cos(gamma*(z-mean_z)*G*T)
  return s/sum(W_bin)

fig, ax = plt.subplots(figsize=(9, 7))
G_fine = range(-320,320)
ax.plot(G_fine,[alpha(G, E_MIN=0, E_MAX=95) for G in G_fine],linewidth=1, color='b')
#ax.plot(G_fine,[alpha(G, 0, E_MAX=125) for G in G_fine],linewidth=1, color='k')
#ax.plot(G_fine,[alpha(G, 0, E_MAX=140) for G in G_fine],linewidth=1, color='r')
ax.plot(G10,aBOT10, 'o', color='m', markersize = 9)

ax.set_xlabel('vertical gradient $G_{10}$ (pT/cm)', fontsize=22)
ax.set_ylabel('$\\alpha(T = 180 s)/ \\alpha_0$ - BOT', fontsize=22)
ax.grid(True)
plt.ylim(0.4, 1)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.tight_layout()
plt.savefig('verticalGradientScan.pdf')

plt.show()

