
#Parçacık Difüzyonu Simülasyonu

import numpy as np
import matplotlib.pyplot as plt
import pyfiglet

banner = pyfiglet.figlet_format("Parcacik Difuzyonu")
print(banner)

#Değerleri Girme Kısmı

D = float(input("Difüzyon katsayısını (m^2/s) giriniz :"))
L = float(input("Ortam uzunluğunu (m) giriniz:"))
Nx = int(input("Uzay adım sayısını giriniz:"))
Nt = int(input("Zaman adım sayısını giriniz:"))
dt = float(input("Zaman adımını giriniz:"))

#Uzay Adımı Hesaplaması

dx = L / Nx

#Courant-Friedrichs-Lewy Condition, Stabilite Kriteri

alpha = D * dt / dx**2
if alpha > 0.5:
    raise ValueError("Verdiğiniz değerler simülasyonu kararsız kılacağı için hata verdim.")

#Belirli Konsantrasyon Yüksekliği

C = np.zeros(Nx)
C[Nx//2] = 1.0

#Zaman adımına göre difüzyonu hesaplama kısmı
for t in range(Nt):
    C_farkli = C.copy()
    for i in range(1, Nx - 1):
        C_farkli[i] = C[i] + alpha * (C[i+1] - 2*C[i] + C[i-1])
    C = C_farkli.copy()

#Konsantrasyon Gradyanı Hesaplama Kısmı ( Fick'in Birinci Yasası )

dCdx = np.gradient(C, dx)
J = -D * dCdx  # =Difüzyon Akısı



