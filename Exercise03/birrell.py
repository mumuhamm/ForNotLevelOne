import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jvp, gamma
import corner

omega_values = np.linspace(0.1, 0.5, 100)
zeta_values  = np.linspace(0.1, 20, 100)
a = 2.5
alpha_k_values = np.zeros((len(omega_values), len(zeta_values)), dtype=complex)
beta_k_values = np.zeros((len(omega_values), len(zeta_values)), dtype=complex)
for i, omega in enumerate(omega_values):
  for j, zeta in enumerate(zeta_values):
    nu = -2j*omega/a 
    nu_abs=np.abs(nu)
    J_nu = jv(nu_abs, zeta)
    J_prime_nu = jvp(nu_abs,zeta, 1 )
    factor = (zeta / 2) ** (1 - 2 * nu)
    alpha_k = (2*np.pi *  factor * gamma(1 + nu) * J_prime_nu * J_nu )/(gamma(1 - nu))
    #print("the omega, and zeta and alpha", omega, zeta, alpha_k)
    alpha_k_values[i, j] = alpha_k
    beta_k = 1 - np.pi * zeta * J_prime_nu * J_nu
    #print("the omega, and zeta and alpha", omega, zeta, beta_k)
    beta_k_values[i, j] = beta_k

# Flatten the arrays for plotting
omega_flat = omega_values[:, None].repeat(len(zeta_values), axis=1).flatten()
alpha_k_real = alpha_k_values.real.flatten()
alpha_k_imag = alpha_k_values.imag.flatten()
beta_k_real = beta_k_values.real.flatten()

# Plot alpha_k (real and imaginary parts) vs omega
plt.figure(figsize=(6, 6))
plt.plot(omega_flat, alpha_k_real, label='Re($\\alpha_k$)', linewidth=1)
plt.plot(omega_flat, alpha_k_imag, label='Im($\\alpha_k$)', linewidth=1)
plt.xlabel(r'$\omega$', fontsize=24)
plt.ylabel(r'$\alpha_k$', fontsize=24)
plt.legend(fontsize=24)
plt.text(0.13, 1.02, r'$\mathbf{FUW}$ $\mathit{Private}$', transform=plt.gca().transAxes, fontsize=12, ha='center', bbox=dict(facecolor='none', edgecolor='none', boxstyle='square'))
plt.tight_layout()
plt.savefig('alphabeta_alpha.pdf')
plt.savefig('alphabeta_alpha.png', dpi=300)
plt.show()

# Plot beta_k (real and imaginary parts) vs omega
plt.figure(figsize=(6,6))
plt.plot(omega_flat, beta_k_real, label='Re($\\beta_k$)', linewidth=1)
plt.xlabel(r'$\omega$', fontsize=24)
plt.ylabel(r'$\beta_k$', fontsize=24)
plt.legend(fontsize=24)
plt.text(0.13, 1.02, 
         r'$\mathbf{FUW}$ $\mathit{Private}$', 
         transform=plt.gca().transAxes, 
         fontsize=12, ha='center', 
         bbox=dict(facecolor='none', 
                   edgecolor='none', boxstyle='square'))
plt.tight_layout()
plt.savefig('alphabeta_beta.pdf')
plt.savefig('alphabeta_beta.png', dpi=300)
plt.show()




# Calculate the probabilities using magnitudes of alpha_k and beta_k
alpha_k_magnitude = np.abs(alpha_k_values)
beta_k_magnitude = np.abs(beta_k_values)
probabilities = alpha_k_magnitude**2 + beta_k_magnitude**2

# Plot likelihood contour
plt.figure(figsize=(8, 6))
X, Y = np.meshgrid(omega_values, zeta_values)
contour = plt.contourf(X, Y, probabilities.T, cmap='viridis')
plt.colorbar(contour, label='Probability')
plt.xlabel(r'$\omega$')
plt.ylabel(r'$\zeta$')
plt.text(0.13, 1.02, r'$\mathbf{FUW}$ $\mathit{Private}$', transform=plt.gca().transAxes, fontsize=12, ha='center')
plt.savefig('alphabetacontour.pdf')
plt.savefig('alphabetacontour.png', dpi=300)
plt.show()



# Create corner plot
data = np.column_stack((alpha_k_magnitude.flatten(), beta_k_magnitude.flatten()))
labels = [r'$|\alpha_k|$', r'$|\beta_k|$']
fig = corner.corner(data, labels=labels, show_titles=True, color='blue', label_kwargs=dict(fontsize=14))
plt.text(-0.72, 1.87, r'$\mathbf{FUW}$ $\mathit{Private}$', transform=plt.gca().transAxes, fontsize=12, ha='center')
plt.savefig('alphabetacorner.pdf')
plt.savefig('alphabetacorner.png', dpi=300)
plt.show()
